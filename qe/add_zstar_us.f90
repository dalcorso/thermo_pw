!
! Copyright (C) 2012-2017 Andrea Dal Corso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
SUBROUTINE add_zstar_us_tpw()
!--------------------------------------------------------
  ! add the contribution of all modes to the effective charges. 
  ! This term is common to both Born effective charges: those calculated 
  ! as dF/dE and those calculated as dP/du.
  ! 
  ! This subroutine is for the USPP and PAW case
  ! it writes the two arrays zstareu0, zstarue0
  !

  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, ityp, atm
  USE cell_base, ONLY : tpiba
  USE klist,     ONLY : xk, wk, ngk, igk_k
  USE gvect,     ONLY : g
  USE wvfct,     ONLY : npw, npwx, nbnd
  USE becmod,    ONLY : calbec, bec_type, allocate_bec_type, &
                        deallocate_bec_type, becscal
  USE noncollin_module, ONLY : noncolin, npol
  USE wavefunctions, ONLY : evc
  USE uspp,      ONLY : vkb, nkb, okvan
  USE qpoint,    ONLY : nksq, npwq, ikks
  USE efield_mod, ONLY : zstarue0, zstareu0, zstarue0_rec
  USE control_ph, ONLY : zue, zeu, done_start_zstar
  USE control_lr, ONLY : rec_code_read
  USE ph_restart, ONLY : ph_writefile
  USE eqv,       ONLY : dpsi
  USE modes,     ONLY : u, nirr, npert
  USE buffers,   ONLY : get_buffer
  USE units_ph,  ONLY : lrcom, iucom
  USE units_lr,  ONLY : iuwfc, lrwfc
  USE uspp_init, ONLY : init_us_2
  USE partial,   ONLY : done_irr, comp_irr
  USE io_global, ONLY : stdout

  USE mp_global, ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,        ONLY : mp_sum

  IMPLICIT NONE

  INTEGER :: ik, ig, ipol, jpol, nrec, mode, ipert, imode0, npe, irr
  INTEGER :: na, ierr, ikk

  REAL(DP) :: weight
  COMPLEX(DP) :: tpibai, tpibac

  COMPLEX(DP), ALLOCATABLE :: dvkb(:,:,:)
  TYPE(bec_type) :: becp2, alphadk(3), bedp, alphapp(3)
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)


  IF (.NOT. comp_irr(0) .OR. done_irr(0) ) RETURN
  IF (.NOT. (zeu.OR.zue).OR. done_start_zstar ) RETURN
  IF (rec_code_read > -30 ) RETURN

  CALL start_clock('add_zstar_us')
  IF (.NOT. okvan) THEN
     zstarue0=(0.0_DP,0.0_DP)
     zstareu0=(0.0_DP,0.0_DP)
!
!   for norm conserving PP we skip this routine, but save the null zstareu0 
!   in any case in the tensors file
!
     GOTO 100
  ENDIF

   tpibac=CMPLX(tpiba,0.0_DP)
   tpibai=CMPLX(0.0_DP, tpiba)

   ALLOCATE(aux1(npwx*npol, nbnd))
   CALL allocate_bec_type (nkb,nbnd,becp2)
   CALL allocate_bec_type (nkb,nbnd,bedp)
   DO ipol=1,3
      CALL allocate_bec_type (nkb,nbnd,alphadk(ipol))
      CALL allocate_bec_type (nkb,nbnd,alphapp(ipol))
   ENDDO

  ALLOCATE (dvkb(npwx,nkb,3))
  DO ik = 1, nksq
     ikk=ikks(ik)
     npwq = ngk(ikk)
     npw=ngk(ikk)
     weight = wk (ikk)
     IF (nksq>1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
     CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)
     !
     ! Calculates  | d/dk beta >
     !
     CALL dvkb3_tpw(ik,dvkb)
     !
     !   alphadk = <-i d/dk d/du beta|psi>
     !   becp2 = < -i d/dk beta | psi>
     !
     DO jpol = 1, 3
        CALL calbec (npw, dvkb(:,:,jpol), evc, becp2)
        CALL becscal((0.0_DP,-1.0_DP), becp2, nkb, nbnd)

        DO ipol = 1, 3
           DO ig = 1, npw
              aux1 (ig, :) = evc(ig,:) *  & 
                ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
           ENDDO
           IF (noncolin) THEN
              DO ig = 1, npw
                 aux1 (ig+npwx, :) = evc(ig+npwx,:) * & 
                ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
              ENDDO
           ENDIF
           CALL calbec(npw, dvkb(:,:,jpol), aux1, alphadk(ipol))
           CALL becscal(tpibac, alphadk(ipol), nkb, nbnd)
        END DO
        ! Then we read  P_c r_alpha |psi> and store it in dpsi
        nrec = (jpol - 1) * nksq + ik
        CALL get_buffer (dpsi, lrcom, iucom, nrec)
        !
        ! products of the beta functions with dpsi 
        !
        CALL calbec (npw, vkb, dpsi, bedp)
        !
        ! products of the d/du beta functions with dpsi 
        !
        DO ipol = 1, 3
           aux1=(0.d0,0.d0)
           DO ig = 1, npw
              aux1 (ig, :) = dpsi(ig,:) *           &
                   ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
           ENDDO
           IF (noncolin) THEN
              DO ig = 1, npw
                 aux1 (ig+npwx, :) = dpsi(ig+npwx,:) *           &
                      ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
              ENDDO
           ENDIF
           CALL calbec ( npw, vkb, aux1, alphapp(ipol) )
           CALL becscal(tpibai, alphapp(ipol), nkb, nbnd)
        ENDDO
!        DO ipol=1,3
!           IF (noncolin) THEN
!              CALL ZSCAL(nkb*npol*nbnd, tpibai, alphapp(ipol)%nc, 1 )
!           ELSE
!              CALL ZSCAL(nkb*npol*nbnd, tpibai, alphapp(ipol)%k, 1 )
!           ENDIF
!        END DO
        !
        !  We these quantities we can calculate for all modes
        !  <psi| dS/du P_c r|psi> + <psi|(dK(r)/du - dS/du)r|psi>
        !
        imode0 = 0
        DO irr = 1, nirr
           npe=npert(irr)
           DO ipert = 1, npe
              mode = imode0 + ipert
              CALL add_dkmds_tpw(ik,u(1,mode),jpol,becp2,alphadk,bedp,alphapp, &
                             weight, zstarue0(mode,jpol) )
           ENDDO
           imode0 = imode0 + npe
        ENDDO
     ENDDO
  ENDDO

  CALL deallocate_bec_type(becp2)
  CALL deallocate_bec_type(bedp)
  DO ipol=1,3
     CALL deallocate_bec_type (alphadk(ipol))
     CALL deallocate_bec_type (alphapp(ipol))
  ENDDO

  DEALLOCATE(aux1)
  DEALLOCATE(dvkb)
!
!  for recovering reasons we need to keep in memory the 
!  effective charges collected among processors
!
   CALL mp_sum ( zstarue0, inter_pool_comm )
!  WRITE(6,*) ' term Z^{(2a}'
!  DO na = 1, nat
!     WRITE( stdout, '(10x," atom ",i6, a6)') na, atm(ityp(na))
!     WRITE( stdout, '(6x,"Ex  (",6f10.5," )")')  (zstarue0 (3*(na-1)+jpol,1), &
!            jpol = 1, 3)
!     WRITE( stdout, '(6x,"Ey  (",6f10.5," )")')  (zstarue0 (3*(na-1)+jpol,2), &
!            jpol = 1, 3)
!     WRITE( stdout, '(6x,"Ez  (",6f10.5," )")')  (zstarue0 (3*(na-1)+jpol,3), &
!            jpol = 1, 3)
!  ENDDO
 
!
!  The initialization term is the same with both methods to compute
!  the effective charges. 
!
  IF (zeu) THEN
     DO ipol=1,3
        DO mode=1,3*nat
           zstareu0(ipol,mode)=zstarue0(mode,ipol)
        ENDDO
     ENDDO
  ENDIF

100 CONTINUE

  done_start_zstar = .TRUE.
  CALL ph_writefile('tensors',0,0,ierr)

  CALL stop_clock('add_zstar_us')

  RETURN
END SUBROUTINE add_zstar_us_tpw

