!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_sm1_psi_tpw (ik, lda, n, m, psi, spsi)
  !----------------------------------------------------------------------------
  !
  ! This subroutine applies the S^{-1} matrix to m wavefunctions psi
  ! and puts the results in spsi.
  ! See Eq.(13) in B. Walker and R. Gebauer, JCP 127, 164106 (2007).
  !
  ! INPUT: ik            k point under consideration
  !        lda           leading dimension of arrays psi, spsi
  !        n             true dimension of psi, spsi
  !        m             number of states psi
  !        psi           the wavefunction to which the S^{-1} 
  !                      matrix is applied
  ! OUTPUT: spsi = S^{-1}*psi
  !
  ! Original routine written by R. Gebauer
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  ! Simplified and generalized to relativistic case by Andrea Dal Corso (2018)
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only
  USE klist,            ONLY : xk, igk_k, ngk
  USE qpoint,           ONLY : ikks, ikqs, nksq
  USE uspp,             ONLY : okvan, vkb, nkb, qq_nt
  USE uspp_param,       ONLY : nh, upf
  USE ions_base,        ONLY : ityp, nat, ntyp=>nsp
  USE becmod,           ONLY : bec_type, becp, calbec
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : intra_bgrp_comm
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE uspp_init,        ONLY : init_us_2
  !
  IMPLICIT NONE
  INTEGER, INTENT(in)      :: lda, n, m, ik
  COMPLEX(DP), INTENT(in)  :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(out) :: spsi(lda*npol,m)
  !
  CALL start_clock( 'lr_sm1_psi_tpw' )
  !
  IF ( gamma_only ) THEN
     CALL sm1_psi_gamma()
  ELSEIF (noncolin) THEN
     CALL sm1_psi_nc()
  ELSE
     CALL sm1_psi_k()
  ENDIF
  !
  CALL stop_clock( 'lr_sm1_psi_tpw' )
  !
  RETURN
  !
CONTAINS
  ! 
    !-----------------------------------------------------------------------
  SUBROUTINE sm1_psi_gamma()
    !-----------------------------------------------------------------------
    !
    ! gamma_only version
    !
    ! Note: 1) vkb must be computed before calling this routine
    !       2) the array bbg must be deallocated somewhere
    !          outside of this routine.
    !
    USE becmod,   ONLY : bec_type,becp,calbec
    USE lrus,     ONLY : bbg
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ibnd
    ! counters
    REAL(DP), ALLOCATABLE :: ps(:,:)
    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda * npol * m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    CALL calbec(n,vkb,psi,becp,m)
    !
    ! Use the array ps as a workspace
    ALLOCATE(ps(nkb,m))
    ps(:,:) = 0.D0
    !
    ! Now let us apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's DO this in 2 steps.
    !
    ! Step 1 : ps = lambda * <beta|psi>
    ! Here, lambda = bbg, and <beta|psi> = becp%r
    !
    call DGEMM( 'N','N',nkb,m,nkb,1.d0,bbg(:,:),nkb,becp%r,nkb,0.d0,ps,nkb)
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta>
    !
    call DGEMM('N','N',2*n,m,nkb,1.d0,vkb,2*lda,ps,nkb,1.d0,spsi,2*lda)
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_gamma

    !-----------------------------------------------------------------------
  SUBROUTINE sm1_psi_k()
    !-----------------------------------------------------------------------
    !
    ! k-points version
    ! Note: the array bbk must be deallocated somewhere
    ! outside of this routine.
    !
    USE lr_lanczos,      ONLY : bbk
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: na, nt, ibnd, ii, jkb
    INTEGER :: ik1, & ! dummy index for k points
               ikk, & ! index of the point k
               ikq, & ! index of the point k+q
               npwq   ! number of the plane-waves at point k+q
    COMPLEX(DP), ALLOCATABLE :: ps(:,:)
    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda*m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    ! Now set up the indices ikk and ikq such that they
    ! correspond to the points k and k+q using the index ik,
    ! which was passed to this routine as an input.
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    !
    IF (n.NE.ngk(ikq)) CALL errore( 'sm1_psiq_k', &
                     & 'Mismatch in the number of plane waves', 1 )
    !
    ! Calculate beta-functions vkb for a given k+q point.
    !
    CALL init_us_2 (n, igk_k(1,ikq), xk(1,ikq), vkb)
    !
    ! Compute the product of the beta-functions vkb with the functions psi
    ! at point k+q, and put the result in becp%k.
    ! becp%k(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    CALL calbec(n, vkb, psi, becp, m)
    !
    ! Use ps as a work space.
    !
    ALLOCATE(ps(nkb,m))
    ps(:,:) = (0.d0,0.d0)
    !
    ! Apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's DO this in 2 steps.
    !
    ! Step 1 : calculate the product 
    ! ps = lambda * <beta|psi>  
    !
    CALL ZGEMM( 'N', 'N', nkb, m, nkb, (1.D0, 0.D0), bbk(1,1,ik), &
                          nkb, becp%k, nkb, (1.D0, 0.D0), ps, nkb )
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta> 
    !
    CALL ZGEMM( 'N', 'N', n, m, nkb, (1.D0, 0.D0), vkb, &
                 lda, ps, nkb, (1.D0, 0.D0), spsi, lda )
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
END SUBROUTINE sm1_psi_k

    !-----------------------------------------------------------------------
SUBROUTINE sm1_psi_nc()
    !-----------------------------------------------------------------------
    !
    ! Noncollinear case
    !
    USE uspp,       ONLY : qq_so
    USE noncollin_module,   ONLY : lspinorb
    USE lr_lanczos, ONLY : bbnc
    USE uspp_init,  ONLY : init_us_2
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ibnd, ii, jkb, ipol, jpol, iis, jkbs
    INTEGER :: ik1, & ! dummy index for k points
               ikk, & ! index of the point k
               ikq, & ! index of the point k+q
               npwq   ! number of the plane-waves at point k+q
    COMPLEX(DP), ALLOCATABLE :: ps(:,:)
    !
    ! Initialize spsi : spsi = psi
    !
    CALL ZCOPY( lda*npol*m, psi, 1, spsi, 1 )
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    ! Now set up the indices ikk and ikq such that they
    ! correspond to the points k and k+q using the index ik,
    ! which was passed to this routine as an input.
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    !
    IF (n/=ngk(ikq)) CALL errore( 'sm1_psiq_nc', &
                     & 'Mismatch in the number of plane waves', 1 )
    !
    ! Calculate beta-functions vkb for a given k+q point.
    !
    CALL init_us_2 (n, igk_k(1,ikq), xk(1,ikq), vkb)
    !
    ! Compute the product of the beta-functions vkb with the functions psi
    ! at point k+q, and put the result in becp%k.
    ! becp%k(ikb,jbnd) = \sum_G vkb^*(ikb,G) psi(G,jbnd) = <beta|psi>
    !
    CALL calbec(n, vkb, psi, becp, m)
    !
    ! Use ps as a work space.
    !
    ALLOCATE(ps(nkb*npol,m))
    ps(:,:) = (0.d0,0.d0)
    !    
    ! Apply the operator S^{-1}, given by Eq.(13), to the functions psi.
    ! Let's DO this in 2 steps.
    !
    ! Step 1 : calculate the product 
    ! ps = lambda * <beta|psi>  
    !
    CALL ZGEMM( 'N', 'N', nkb*npol, m, nkb*npol, (1.D0, 0.D0), bbnc(1,1,ik), &
                  nkb*npol, becp%nc, nkb*npol, (1.D0, 0.D0), ps, nkb*npol )
    !
    ! Step 2 : |spsi> = S^{-1} * |psi> = |psi> + ps * |beta> 
    !
    CALL ZGEMM( 'N', 'N', n, m*npol, nkb, (1.D0, 0.D0), vkb, &
                     & lda, ps, nkb, (1.D0, 0.D0), spsi, lda )
    !
    DEALLOCATE(ps)
    !
    RETURN
    !
END SUBROUTINE sm1_psi_nc

END SUBROUTINE lr_sm1_psi_tpw

!-----------------------------------------------------------------------
SUBROUTINE lr_sm1_initialize()
!-----------------------------------------------------------------------
!
!   This routine initializes the coefficients of S^-1 and saves them in
!   bbg, bbk or bbnc
!
USE kinds,            ONLY : DP
USE control_flags,    ONLY : gamma_only
USE lrus,             ONLY : bbg
USE lr_lanczos,       ONLY : bbk, bbnc
USE klist,            ONLY : xk, ngk, igk_k
USE qpoint,           ONLY : nksq, ikks, ikqs
USE wvfct,            ONLY : npwx
USE becmod,           ONLY : calbec
USE uspp,             ONLY : vkb, nkb, qq_nt, qq_so
USE uspp_param,       ONLY : nh, upf
USE ions_base,        ONLY : ityp,nat,ntyp=>nsp
USE mp,               ONLY : mp_sum
USE mp_global,        ONLY : intra_bgrp_comm
USE noncollin_module, ONLY : noncolin, npol
USE matrix_inversion, ONLY : invmat
USE uspp_init,        ONLY : init_us_2

IMPLICIT NONE
!
INTEGER :: ik1, ikk, ikq, npw, npwq, na, nt, ikb, jkb, ijkb0, ih, jh, ii, &
           ipol, jpol, kpol, ijs, iis, ikbs, jkbs, iks, kjs
REAL(DP),    ALLOCATABLE :: psr(:,:)
COMPLEX(DP), ALLOCATABLE :: ps(:,:), bbnc_nc(:,:)

IF (gamma_only) THEN
   bbg = 0.0d0
   ALLOCATE(psr(nkb,nkb))
ELSE
   IF (noncolin) THEN
      ALLOCATE(bbnc_nc(nkb,nkb))
      bbnc_nc = (0.0d0,0.0d0)
   ELSE   
      bbk = (0.0d0,0.0d0)
   ENDIF
   ALLOCATE(ps(nkb*npol,nkb*npol))
ENDIF

DO ik1 = 1, nksq
   !
   ikk  = ikks(ik1)
   ikq  = ikqs(ik1)
   npw  = ngk(ikk)
   npwq = ngk(ikq)
   !
   ! Calculate beta-functions vkb for a given k+q point.
   !
   CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkb)
   !
   ! Calculate the coefficients B_ij defined by Eq.(15).
   ! B_ij = <beta(i)|beta(j)>, where beta(i) = vkb(i).
   !
   IF (gamma_only) THEN
      CALL calbec (npw, vkb, vkb, bbg, nkb)
   ELSEIF (noncolin) THEN
      CALL zgemm('C','N',nkb,nkb,npwq,(1.d0,0.d0),vkb, &
               npwx,vkb,npwx,(0.d0,0.d0),bbnc_nc,nkb)
      CALL mp_sum(bbnc_nc, intra_bgrp_comm)
   ELSE
      CALL zgemm('C','N',nkb,nkb,npwq,(1.d0,0.d0),vkb, &
               npwx,vkb,npwx,(0.d0,0.d0),bbk(1,1,ik1),nkb)
      CALL mp_sum(bbk(:,:,ik1), intra_bgrp_comm)
   ENDIF
   ! 
   ps(:,:) = (0.d0,0.d0)
   !
   ! Calculate the product of q_nm and B_ij of Eq.(16).
   !
   ijkb0 = 0
   DO nt=1,ntyp
      IF (upf(nt)%tvanp) THEN
         DO na=1,nat
            IF (ityp(na)==nt) THEN
               DO ii=1,nkb
                  DO jh=1,nh(nt)
                     !
                     jkb=ijkb0 + jh
                     !
                     DO ih=1,nh(nt)
                        !
                        ikb = ijkb0 + ih
                        !
                        IF (gamma_only) THEN
                           psr(ikb,ii) = psr(ikb,ii) + bbg(jkb,ii)   &
                                                     * qq_nt(ih,jh,nt) 
                        ELSEIF(noncolin) THEN
                           ijs=0
                           DO ipol=1, npol
                              ikbs = ikb + nkb * ( ipol - 1 )
                              DO jpol=1, npol
                                 iis = ii + nkb * ( jpol - 1 )
                                 ijs = ijs + 1
                                 ps(ikbs,iis) = ps(ikbs,iis) + &
                                       bbnc_nc(jkb,ii)*qq_so(ih,jh,ijs,nt)
                              ENDDO
                           ENDDO
                        ELSE
                           ps(ikb,ii) = ps(ikb,ii) + bbk(jkb,ii,ik1) &
                                                   * qq_nt(ih,jh,nt)
                        ENDIF
                           !
                     ENDDO
                  ENDDO
               ENDDO
               ijkb0 = ijkb0+nh(nt)
            ENDIF
         ENDDO
      ELSE
         DO na = 1, nat
            IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
         ENDDO
      ENDIF
   ENDDO
   !
   ! Add an identity to q_nm * B_ij [see Eq.(16)].
   ! ps = (1 + q*B)
   ! 
   IF (gamma_only) THEN
      DO ii=1,nkb
         psr(ii,ii) = psr(ii,ii) + 1.d0
      ENDDO
   ELSE
      DO ii=1,nkb*npol
         ps(ii,ii) = ps(ii,ii) + (1.d0,0.d0)
      ENDDO
   ENDIF
   !
   ! Invert matrix: (1 + q*B)^{-1}
   ! 
   IF (gamma_only) THEN
      CALL invmat( nkb, psr )
   ELSE
      CALL invmat( nkb*npol, ps )
   ENDIF
   !
   ! Use the array bbk as a work space in order to save the memory.
   !
   IF (gamma_only) THEN
      bbg(:,:)=0.0_DP
   ELSEIF (noncolin) THEN
      bbnc(:,:,ik1) = (0.d0,0.d0)
   ELSE
      bbk(:,:,ik1) = (0.d0,0.d0)
   ENDIF
   !
   ! Finally, let us calculate lambda_nm = -(1+q*B)^{-1} * q
   ! Let us use the array bbk and put there the values of the lambda
   ! coefficients.
   !
   ijkb0 = 0
   DO nt=1,ntyp
      IF (upf(nt)%tvanp) THEN
         DO na=1,nat
            IF (ityp(na)==nt) THEN
               DO ii=1,nkb*npol
                  DO jh=1,nh(nt)
                     !
                     jkb=ijkb0 + jh
                     !
                     DO ih=1,nh(nt)
                        !
                        ikb = ijkb0 + ih
                        !
                        IF (gamma_only) THEN
                           bbg(ii,jkb) = bbg(ii,jkb) &
                                       - ps(ii,ikb) * qq_nt(ih,jh,nt)
 
                        ELSEIF (noncolin) THEN
                           kjs = 0
                           DO kpol=1,npol
                              ikbs = ikb + nkb * (kpol-1)
                              DO jpol=1,npol
                                 jkbs=jkb + nkb * (jpol-1)
                                 kjs=kjs+1
                                 bbnc(ii,jkbs,ik1) = &
                                         bbnc(ii,jkbs,ik1) - &
                                         ps(ii,ikbs)*qq_so(ih,jh,kjs,nt)
                              ENDDO
                           ENDDO
                        ELSE
                           bbk(ii,jkb,ik1) = bbk(ii,jkb,ik1) - &
                                        ps(ii,ikb) * qq_nt(ih,jh,nt)
                        ENDIF 
                        !
                     ENDDO
                  ENDDO
               ENDDO
               ijkb0 = ijkb0+nh(nt)
            ENDIF
         ENDDO
      ELSE
         DO na = 1, nat
            IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
         ENDDO
      ENDIF
   ENDDO
   !
ENDDO ! loop on k points
!
IF (gamma_only) THEN
   DEALLOCATE(psr)
ELSE
   IF (noncolin) DEALLOCATE(bbnc_nc)
   DEALLOCATE(ps)
ENDIF

RETURN
END SUBROUTINE lr_sm1_initialize
