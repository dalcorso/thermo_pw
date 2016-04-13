!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine polarizc ( iu )
  !-----------------------------------------------------------------------
  !
  !      calculates the frequency dependent polarizability
  !

  USE io_global,    ONLY : stdout
  USE io_files,     ONLY : iunigk
  USE constants,    ONLY : fpi, rytoev
  USE cell_base,    ONLY : at, bg
  USE klist,        ONLY : wk, nkstot
  USE symme,        ONLY : symmatrix, crys_to_cart
  USE wvfct,        ONLY : npw, npwx, igk
  USE kinds,        ONLY : DP
  USE control_lr,   ONLY : nbnd_occ
  USE lsda_mod,     ONLY : lsda
  USE units_ph,     ONLY : lrdwf, iudwf, lrebar, iuebar
  USE buffers,      ONLY : get_buffer
  USE freq_ph,      ONLY : done_iu, comp_iu, fiu
  USE optical,      ONLY : iu1dwf, lr1dwf, current_w, epsilonc, fru
  USE eqv,          ONLY : dpsi, dvpsi
  USE qpoint,       ONLY : nksq
  USE ph_restart,   ONLY : ph_writefile
  USE cell_base,    ONLY : omega
  USE noncollin_module, ONLY : noncolin
  USE mp_pools,     ONLY : inter_pool_comm
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iu
  !
  ! local variables
  !
  integer :: ibnd, ipol, jpol, nrec, ik, ierr
  ! counter on polarizations
  ! counter on records
  ! counter on k points
  real(kind=DP) :: w, weight, repsilon(3,3), iepsilon(3,3)
  COMPLEX(kind=DP) :: cepsilon(3,3), alpha(3,3)

  complex(kind=DP), EXTERNAL :: zdotc

  call start_clock ('polariz')
  cepsilon(:,:) = (0.0_DP, 0.0_DP)
  if (nksq > 1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (nksq > 1) read (iunigk) npw, igk
     IF (ABS(CMPLX(fru(iu),fiu(iu)))> 1.D-7) THEN
!
!   two wavefunctions at (+w and -w) are added in this case
!
        weight = wk (ik)
     ELSE
!
!   only wavefunction is added in this case
!
        weight = wk(ik) * 2.0_DP
     END IF
     w = fpi * weight / omega
     do ipol = 1, 3
        nrec = (ipol - 1) * nksq + ik
        call get_buffer (dvpsi, lrebar, iuebar, nrec)
        do jpol = 1, 3
           nrec = (jpol - 1) * nksq + ik
           call get_buffer(dpsi, lrdwf, iudwf, nrec)
           do ibnd = 1, nbnd_occ (ik)
              !
              !  this is  <DeltaV*psi(E)|DeltaPsi(E, w)>
              !
              cepsilon(ipol,jpol)=cepsilon(ipol,jpol)-2.d0*w* &
                   zdotc (npw, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1) 
              IF (noncolin) &
                 cepsilon(ipol,jpol)=cepsilon(ipol,jpol)-2.d0*w* &
                   zdotc (npw, dvpsi (1+npwx, ibnd), 1, dpsi (1+npwx, ibnd), 1) 
           enddo
           IF (ABS(CMPLX(fru(iu),fiu(iu))) >= 1.D-7) THEN
!
!   Add the term at -w
!
              call get_buffer(dpsi, lr1dwf, iu1dwf, nrec)
              do ibnd = 1, nbnd_occ (ik)
                 !
                 !  this is  <DeltaV*psi(E)|DeltaPsi(E,-w)>
                 !
                 cepsilon(ipol,jpol)=cepsilon(ipol,jpol)-2.d0*w* &
                     zdotc (npw, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1) 
                 IF (noncolin) &
                 cepsilon(ipol,jpol)=cepsilon(ipol,jpol)-2.d0*w* &
                     zdotc (npw,dvpsi(1+npwx, ibnd),1,dpsi(1+npwx, ibnd),1) 
              enddo
           END IF
        enddo
     enddo
  enddo
  call mp_sum ( cepsilon, intra_bgrp_comm )
  call mp_sum ( cepsilon, inter_pool_comm )
  !
  !      symmetrize
  !
  !       WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
  !       WRITE( stdout,'(10x,"(",3f15.5," )")') ((repsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  repsilon=DREAL(cepsilon)
  iepsilon=DIMAG(cepsilon)
  call crys_to_cart ( repsilon )
  call symmatrix ( repsilon )
  call crys_to_cart ( iepsilon )
  call symmatrix ( iepsilon )
  !
  !    pass to cartesian axis
  !
  !      WRITE( stdout,'(/,10x,"Symmetrized in cartesian axis ",/)')
  !      WRITE( stdout,'(10x,"(",3f15.5," )")') ((repsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  !      WRITE( stdout,'(10x,"(",3f15.5," )")') ((iepsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     repsilon (ipol, ipol) = repsilon (ipol, ipol) + 1.0_DP
  enddo
  cepsilon=CMPLX(repsilon, iepsilon)
  epsilonc(:,:,iu)=cepsilon(:,:)
  !
  !  and print the result
  !
  WRITE( stdout, '(/,10x,"Dielectric constant at &
                     &frequency",f9.4," +",f9.4," i Ry")') current_w 
  WRITE( stdout, '(42x,f9.4," +",f9.4," i eV")') current_w * rytoev
  WRITE( stdout, '(/,10x,"Real part ",f5.2,/)') 
  WRITE( stdout, '(10x,"(",3f18.9," )")') ((repsilon(ipol,jpol), ipol=1,3), jpol=1,3)
  WRITE( stdout, '(/,10x,"Imaginary part ",f5.2,/)') 
  WRITE( stdout, '(10x,"(",3f18.9," )")') ((iepsilon(ipol,jpol), ipol=1,3), jpol=1,3)

  IF (nkstot==1 .OR. (nkstot==2.AND.lsda)) CALL write_polarizc(cepsilon, iu)
!  done_iu(iu)=.TRUE.
!  call ph_writefile('polarization',0,iu,ierr)
  !
  call stop_clock ('polariz')

  return
end subroutine polarizc

  SUBROUTINE write_polarizc(cepsilon, iu)
!
!  This routine write on output the
!
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE cell_base,  ONLY : omega
  USE constants,  ONLY : fpi, BOHR_RADIUS_ANGS
  USE optical,    ONLY : polarc
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iu
  INTEGER :: ipol, jpol
  COMPLEX(DP), INTENT(IN) :: cepsilon(3,3)
  COMPLEX(DP) :: alpha(3,3)
  !
  ! compute the polarization
  !
  alpha=(0.0_DP, 0.0_DP)
  do ipol = 1, 3
     do jpol = 1, 3
        IF (ABS(cepsilon(ipol,jpol)) > 1.D-4) &
        alpha (ipol, jpol)=(3.d0*omega/fpi)*(cepsilon (ipol, jpol) - 1.0_DP )/ &
                                            (cepsilon (ipol, jpol) + 2.0_DP )
     enddo
  enddo
  polarc(:,:,iu)=alpha(:,:)

  WRITE(stdout,'(/,5x,"Polarizability (a.u.)^3",20x,"Polarizability (A^3)")')
  WRITE(stdout,'(5x," Cartesian axis  Real part",/ )')
  WRITE(stdout,'(3f10.4,5x,3f14.4)')((DREAL(polarc(ipol,jpol,iu)), jpol=1,3), &
         (DREAL(polarc(ipol,jpol,iu))*BOHR_RADIUS_ANGS**3, jpol=1,3), ipol=1,3)

  WRITE(stdout,'(/,5x," Cartesian axis  Imaginary part",/ )')
  WRITE(stdout,'(3f10.4,5x,3f14.4)')((DIMAG(polarc(ipol,jpol,iu)), jpol=1,3), &
         (DIMAG(polarc(ipol,jpol,iu))*BOHR_RADIUS_ANGS**3, jpol=1,3), ipol=1,3)
  RETURN
  END SUBROUTINE write_polarizc

