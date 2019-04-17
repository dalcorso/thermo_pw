!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE interpolate_cv(vmin_t, celldm_t, ph_cv, cv_t)
!
!  This subroutine receives the isochoric heat capacity at several 
!  geometries, the equilibrium volume (or celldm) as a function 
!  of temperature, and interpolates the isochoric heat capacity at the 
!  equilibrium volume (or celldm) at each temperature.
!  The routine is parallelized over temperatures and must be called
!  by all processors.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE thermo_mod,     ONLY : celldm_geo, omega_geo, no_ph, tot_ngeo
USE control_mur,    ONLY : lmurn
USE temperature,    ONLY : ntemp
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic, &
                           quadratic_var
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic, & 
                           quartic_var
USE cubic_surfaces, ONLY : fit_multi_cubic, evaluate_fit_cubic, cubic_var
USE control_quartic_energy, ONLY : lsolve, poly_degree_cv
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: vmin_t(ntemp), ph_cv(ntemp, tot_ngeo), &
                         celldm_t(6,ntemp)
REAL(DP), INTENT(OUT) :: cv_t(ntemp)

INTEGER :: itemp, igeo, degree, nvar, ndata, startt, lastt
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:), coeff(:)

degree=crystal_parameters(ibrav)
IF (lmurn) degree=1         ! only the volume is variable in this case
IF (poly_degree_cv==4) THEN
   nvar=quartic_var(degree)
ELSEIF (poly_degree_cv==3) THEN
   nvar=cubic_var(degree)
ELSE
   nvar=quadratic_var(degree)
ENDIF
!
ndata = compute_nwork_ph(no_ph,tot_ngeo)

ALLOCATE(x(degree,ndata))
ALLOCATE(x_pos_min(degree))
ALLOCATE(f(ndata))
ALLOCATE(coeff(nvar))
!
!  collect the geometrical data of the geometries for which phonon dispersions
!  have been calculated
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   ndata=ndata+1
   IF (lmurn) THEN
      x(1,ndata)=omega_geo(igeo)
   ELSE
      CALL compress_celldm(celldm_geo(1,igeo), x(1,ndata), degree, ibrav)
   ENDIF
ENDDO
!
!  Divide temperatures among processors
!
CALL divide(world_comm, ntemp, startt, lastt)
cv_t=0.0_DP
DO itemp=startt,lastt
!
!  collect the computed data at this temperatire
!
   ndata=0
   DO igeo=1, tot_ngeo
      IF (no_ph(igeo)) CYCLE
      ndata=ndata+1
      f(ndata)= ph_cv(itemp,igeo)
   ENDDO
!
!  compute the geometry corresponding to this temperature.
!
   IF (lmurn) THEN
      x_pos_min(1)=vmin_t(itemp)
   ELSE
      CALL compress_celldm(celldm_t(1,itemp), x_pos_min, degree, ibrav)
   ENDIF

   IF (poly_degree_cv==4) THEN
!
!   compute the coefficients of the quartic polynomial
!
      CALL fit_multi_quartic(ndata, degree, nvar, lsolve, x, f, coeff)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_quartic(degree, nvar, x_pos_min, cv_t(itemp), coeff)
   ELSEIF( poly_degree_cv==3) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL fit_multi_cubic(ndata, degree, nvar, lsolve, x, f, coeff)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_cubic(degree, nvar, x_pos_min, cv_t(itemp), coeff)
   ELSE
!
!   compute the coefficients of the quadratic polynomial
!
      CALL fit_multi_quadratic(ndata, degree, nvar, x, f, coeff)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_quadratic(degree, nvar, x_pos_min, cv_t(itemp), coeff)

   ENDIF
ENDDO
!
!  collect all temperatures on all processors
!
CALL mp_sum(cv_t, world_comm)

DEALLOCATE(coeff)
DEALLOCATE(f)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_cv
