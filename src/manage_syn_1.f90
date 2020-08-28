!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_syn_1(nwork)
!
!  This subroutine manages the synchronous work after the first set of
!  possibly asynchronous pw.x calculations.
!

USE control_thermo,   ONLY : lev_syn_1, lconv_ke_test, lconv_nk_test,  &
                             lectqha

USE control_elastic_constants, ONLY : ngeom

IMPLICIT NONE
INTEGER :: nwork
  !
  !  In the kinetic energy test write the results
  !
  IF (lconv_ke_test) THEN
     CALL write_e_ke()
     CALL plot_e_ke()
  ENDIF
  !
  ! In the k-point test write the results
  !
  IF (lconv_nk_test) THEN
     CALL write_e_nk()
     CALL plot_e_nk()
  ENDIF
!
!  In a Murnaghan equation calculation determine the lattice constant,
!  bulk modulus and its pressure derivative and write the results.
!  Otherwise interpolate the energy with a quadratic or quartic polynomial.
!
  IF (lev_syn_1) CALL manage_energy_minimum(nwork)
!
!  When computing the elastic constants as a function of temperature
!  here we have the energy for each strained geometry and can compute
!  the T=0 K elastic constants for all the reference geometries.
!
  IF (lectqha) CALL manage_elastic_cons(nwork,ngeom)

  RETURN
END SUBROUTINE manage_syn_1

