! 
! Copyright (C) 2012-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE add_dkmds_tpw(ik,uact,jpol,becp2,alphadk,bedp,alphapp,weight,zstar)
  !--------=========-------------------------------------------------------
  !
  ! This subroutine calculates a part of the contribution to zstar 
  ! that does not depend on the perturbed wavefunctions. The other
  ! part is calculated in compute_drhous.
  ! It assumes that the variable dpqq has been set. In the noncollinear
  ! and spin_orbit case the variable dpqq_so must be set. 
  !

  USE kinds,     ONLY : DP
  USE constants, ONLY : eps12
  USE uspp,      ONLY : nkb, qq_nt, qq_so
  USE cell_base, ONLY : at
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module, ONLY : noncolin, npol, lspinorb
  USE uspp_param, ONLY: nh
  USE phus,       ONLY: alphap
  USE lrus,      ONLY : becp1, dpqq, dpqq_so
  USE qpoint,    ONLY : ikks
  USE control_lr, ONLY : nbnd_occ
  USE becmod,    ONLY : bec_type

  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum


  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ik, jpol
  REAL(DP), INTENT(IN) :: weight
  COMPLEX(DP), INTENT(INOUT) :: zstar
  COMPLEX(DP), INTENT(IN) :: uact (3 * nat)
  TYPE(bec_type), INTENT(IN) :: alphadk(3), becp2, bedp, alphapp(3)

  INTEGER :: ipol, ijkb0, nt, na, ih, jh, ikb, jkb, ibnd, mu, is, &
             js, ijs, startb, lastb, nbnd_eff, ikk

  COMPLEX(DP), ALLOCATABLE :: ps1(:,:,:), ps2(:,:),  ps3(:,:,:), ps4(:,:) 
  COMPLEX(DP), ALLOCATABLE :: ps1_nc(:,:,:,:), ps2_nc(:,:,:),  &
                              ps3_nc(:,:,:,:), ps4_nc(:,:,:) 

  COMPLEX(DP), EXTERNAL :: zdotc

  CALL start_clock('add_dkmds')
  ikk=ikks(ik)
  nbnd_eff=nbnd_occ(ikk)
  IF (nkb.GT.0) THEN 
     IF (noncolin) THEN
        ALLOCATE (ps1_nc(nkb,npol,nbnd_eff,3))
        ALLOCATE (ps2_nc(nkb,npol,nbnd_eff))
        ALLOCATE (ps3_nc(nkb,npol,nbnd_eff,3))
        ALLOCATE (ps4_nc(nkb,npol,nbnd_eff))
     ELSE
        ALLOCATE (ps1(nkb,nbnd_eff,3))
        ALLOCATE (ps2(nkb,nbnd_eff))
        ALLOCATE (ps3(nkb,nbnd_eff,3))
        ALLOCATE (ps4(nkb,nbnd_eff))
     END IF
  END IF
!
!  Parallelize on the bands
!
  CALL divide (intra_bgrp_comm, nbnd_eff, startb, lastb)

  IF (noncolin) THEN
     ps1_nc = (0.d0, 0.d0)
     ps2_nc = (0.d0, 0.d0)
     ps3_nc = (0.d0, 0.d0)
     ps4_nc = (0.d0, 0.d0)
  ELSE
     ps1 = (0.d0, 0.d0)
     ps2 = (0.d0, 0.d0)
     ps3 = (0.d0, 0.d0)
     ps4 = (0.d0, 0.d0)
  ENDIF

  ijkb0 = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (ityp(na).EQ.nt) THEN
           mu = 3 * (na - 1)
           IF ( ABS (uact (mu + 1) ) + &
                ABS (uact (mu + 2) ) + &
                ABS (uact (mu + 3) ) > eps12) then
              DO ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ipol = 1, 3 
                       DO ibnd=startb, lastb
                          IF (noncolin) then 
                             IF (lspinorb) then
                                ijs=0
                                DO is=1,npol
                                   DO js=1,npol
                                      ijs=ijs+1
                                      ps1_nc(ikb,is,ibnd,ipol)=  &
                                                ps1_nc(ikb,is,ibnd,ipol) &
                                         + qq_so(ih,jh,ijs,nt)*   &
                                              becp1(ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                      ps2_nc(ikb,is,ibnd)=  &
                                                ps2_nc(ikb,is,ibnd) &
                                         + qq_so(ih,jh,ijs,nt)* &
                                           alphap(ipol,ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                      ps3_nc(ikb,is,ibnd,ipol)=  &
                                                ps3_nc(ikb,is,ibnd,ipol) &
                                         + dpqq_so(ih,jh,ijs,jpol,nt)* &
                                              becp1(ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                      ps4_nc(ikb,is,ibnd)=  &
                                                ps4_nc(ikb,is,ibnd) &
                                         + dpqq_so(ih,jh,ijs,jpol,nt)* &
                                           alphap(ipol,ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                   ENDDO
                                ENDDO
                             ELSE
                                DO is=1,npol
                                   ps1_nc (ikb,is,ibnd,ipol)= &
                                        ps1_nc(ikb,is,ibnd,ipol)+  &
                                      qq_nt (ih, jh, nt) *                &
                                      becp1(ik)%nc(jkb,is,ibnd)*  &
                                              uact (mu + ipol)
                                   ps2_nc(ikb,is,ibnd)=  &
                                             ps2_nc(ikb,is,ibnd) &
                                     + qq_nt(ih,jh,nt)* &
                                        alphap(ipol,ik)%nc(jkb,is,ibnd)*  &
                                          uact (mu + ipol)
                                   ps3_nc(ikb,is,ibnd,ipol)=  &
                                            ps3_nc(ikb,is,ibnd,ipol) &
                                         + dpqq(ih,jh,jpol,nt)* &
                                              becp1(ik)%nc(jkb,is,ibnd)*  &
                                              uact (mu + ipol)
                                   ps4_nc(ikb,is,ibnd)=  &
                                            ps4_nc(ikb,is,ibnd) &
                                      + dpqq(ih,jh,jpol,nt)* &
                                           alphap(ipol,ik)%nc(jkb,is,ibnd)*  &
                                           uact (mu + ipol)
                                ENDDO
                             ENDIF
                          ELSE
                             ps1 (ikb,ibnd,ipol)= &
                                ps1(ikb,ibnd,ipol)+  &
                                  qq_nt (ih, jh, nt) *         &
                                  becp1(ik)%k(jkb,ibnd) * uact (mu + ipol)
                             ps2 (ikb,ibnd)=  &
                                 ps2(ikb,ibnd) + &
                                  qq_nt(ih,jh,nt)* &
                                  alphap(ipol,ik)%k(jkb,ibnd)*uact (mu + ipol)
                             ps3(ikb,ibnd,ipol)=  &
                                 ps3(ikb,ibnd,ipol) &
                                  + dpqq(ih,jh,jpol,nt)* &
                                    becp1(ik)%k(jkb,ibnd)* uact (mu + ipol)
                             ps4(ikb,ibnd)=  &
                                    ps4(ikb,ibnd) &
                                 + dpqq(ih,jh,jpol,nt)* &
                                   alphap(ipol,ik)%k(jkb,ibnd)* uact (mu + ipol)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
           ijkb0=ijkb0+nh(nt)
        ENDIF
     ENDDO
  ENDDO
  IF (noncolin) THEN
     CALL mp_sum ( ps1_nc, intra_bgrp_comm )
     CALL mp_sum ( ps2_nc, intra_bgrp_comm )
     CALL mp_sum ( ps3_nc, intra_bgrp_comm )
     CALL mp_sum ( ps4_nc, intra_bgrp_comm )
  ELSE
     CALL mp_sum ( ps1, intra_bgrp_comm )
     CALL mp_sum ( ps2, intra_bgrp_comm )
     CALL mp_sum ( ps3, intra_bgrp_comm )
     CALL mp_sum ( ps4, intra_bgrp_comm )
  ENDIF 
  !
  !  Finally multiply by the becp and alphap terms
  !
  IF ( nkb > 0 ) THEN
     IF (noncolin) THEN
        zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, bedp%nc, 1, ps2_nc, 1)
        zstar=zstar+weight*zdotc(nkb*nbnd_eff*npol, becp2%nc, 1, ps2_nc, 1)
        zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, becp1(ik)%nc, 1, ps4_nc, 1)
        DO ipol=1,3
           zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, alphapp(ipol)%nc, 1, ps1_nc(1,1,1,ipol), 1)
           zstar=zstar+weight*zdotc(nkb*nbnd_eff*npol, alphadk(ipol)%nc, 1, ps1_nc(1,1,1,ipol), 1)
           zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, alphap(ipol,ik)%nc, 1, ps3_nc(1,1,1,ipol), 1)

        ENDDO
     ELSE
        zstar=zstar-weight*zdotc(nkb*nbnd_eff,bedp%k,1,ps2,1)
        zstar=zstar+weight*zdotc(nkb*nbnd_eff,becp2%k,1,ps2,1)
        zstar=zstar-weight*zdotc(nkb*nbnd_eff,becp1(ik)%k,1,ps4,1)
        DO ipol=1,3
           zstar=zstar-weight*zdotc(nkb*nbnd_eff,alphapp(ipol)%k,1,ps1(1,1,ipol),1)
           zstar=zstar+weight*zdotc(nkb*nbnd_eff,alphadk(ipol)%k,1,ps1(1,1,ipol),1)
           zstar=zstar-weight*zdotc(nkb*nbnd_eff,alphap(ipol,ik)%k,1,ps3(1,1,ipol),1)
        END DO
     END IF
  END IF

  IF (noncolin) THEN
     IF (ALLOCATED(ps1_nc)) deallocate(ps1_nc)
     IF (ALLOCATED(ps2_nc)) deallocate(ps2_nc)
     IF (ALLOCATED(ps3_nc)) deallocate(ps3_nc)
     IF (ALLOCATED(ps4_nc)) deallocate(ps4_nc)
  ELSE
     IF (ALLOCATED(ps1)) deallocate(ps1)
     IF (ALLOCATED(ps2)) deallocate(ps2)
     IF (ALLOCATED(ps3)) deallocate(ps3)
     IF (ALLOCATED(ps4)) deallocate(ps4)
  ENDIF

  CALL stop_clock('add_dkmds')
  RETURN

END SUBROUTINE add_dkmds_tpw
