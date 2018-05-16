!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE dv_of_drho_clf

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE dv_of_drho_nlf (dvscf)
  !-----------------------------------------------------------------------
  !
  !  This routine computes the long range hartree potential that
  !  has to be included in EELS calculation even in the lnoloc case.
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : e2, fpi
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE gvect,             ONLY : nl, gstart
  USE cell_base,         ONLY : alat, tpiba2, omega
  USE noncollin_module,  ONLY : nspin_lsda, nspin_mag
  USE qpoint,            ONLY : xq

  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin_mag)
  ! input:  response charge density
  ! output: response Hartree potential
  
  INTEGER :: is
  ! counter on spin polarizations
  REAL(DP) :: qg2
  ! qg2: the modulus of (q+G)^2
  COMPLEX(DP), ALLOCATABLE :: dvaux(:,:)

  CALL start_clock ('dv_of_drho')
  !
  ALLOCATE (dvaux( dfftp%nnr, nspin_mag))
  dvaux (:,:) = (0.d0, 0.d0)
  !
  ! Copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  if (nspin_mag == 2) then
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
  end if
  !
  CALL fwfft ('Dense', dvscf(:,1), dfftp)
  !
  ! Response Hartree potential
  !
  DO is = 1, nspin_lsda
     IF (gstart==2) THEN
        qg2 = xq(1)**2 + xq(2)**2 + xq(3)**2
        IF (qg2 > 1.d-8) THEN
           dvaux(nl(1),is) = dvaux(nl(1),is) + &
                                 & e2 * fpi * dvscf(nl(1),1) / (tpiba2 * qg2)
        ENDIF
     ENDIF
     CALL invfft ('Dense', dvaux (:, is), dfftp)
  ENDDO
  !
  dvscf (:,:) = dvaux (:,:)
  !
  DEALLOCATE (dvaux)
  !
  CALL stop_clock ('dv_of_drho')
  !
  RETURN
  !
END SUBROUTINE dv_of_drho_nlf

END MODULE dv_of_drho_clf
