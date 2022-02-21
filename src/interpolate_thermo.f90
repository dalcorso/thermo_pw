!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_thermo(vmin_t, celldm_t, ph_thermo, thermo_t)
!-----------------------------------------------------------------------
!
!  This subroutine receives the desired thermodynamical quantity at 
!  several geometries (specific heat, entropy or vibrational energy), 
!  the equilibrium volume (or celldm) as a function of temperature, 
!  and interpolates it at the equilibrium volume (or celldm) at each 
!  temperature. The routine is parallelized over temperatures and must 
!  be called by all processors.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE thermo_mod,     ONLY : celldm_geo, omega_geo, no_ph, tot_ngeo
USE control_mur,    ONLY : lmurn
USE temperature,    ONLY : ntemp
USE linear_surfaces, ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_quartic_energy, ONLY : poly_degree_thermo, lsolve
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE polynomial,     ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: vmin_t(ntemp), ph_thermo(ntemp, tot_ngeo), &
                         celldm_t(6,ntemp)

REAL(DP), INTENT(OUT) :: thermo_t(ntemp)

INTEGER :: itemp, igeo, nvar, ndata, startt, lastt
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:)
TYPE(poly1) :: p1
TYPE(poly2) :: p2
TYPE(poly3) :: p3
TYPE(poly4) :: p4

nvar=crystal_parameters(ibrav)
IF (lmurn) nvar=1         ! only the volume is variable in this case
!
ndata = compute_nwork_ph(no_ph,tot_ngeo)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
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
      CALL compress_celldm(celldm_geo(1,igeo), x(1,ndata), nvar, ibrav)
   ENDIF
ENDDO
!
!  Divide temperatures among processors
!
CALL divide(world_comm, ntemp, startt, lastt)
thermo_t=0.0_DP
DO itemp=startt,lastt
!
!  collect the computed data at this temperature
!
   ndata=0
   DO igeo=1, tot_ngeo
      IF (no_ph(igeo)) CYCLE
      ndata=ndata+1
      f(ndata)= ph_thermo(itemp,igeo)
   ENDDO
!
!  compute the geometry corresponding to this temperature.
!
   IF (lmurn) THEN
      x_pos_min(1)=vmin_t(itemp)
   ELSE
      CALL compress_celldm(celldm_t(1,itemp), x_pos_min, nvar, ibrav)
   ENDIF

   IF (poly_degree_thermo == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
      CALL init_poly(nvar,p4)
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4)
!
!  and evaluate the polynomial
!

      CALL evaluate_fit_quartic(nvar, x_pos_min, thermo_t(itemp), p4)
      CALL clean_poly(p4)
   ELSEIF (poly_degree_thermo == 3) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p3)
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, f, p3)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_cubic(nvar, x_pos_min, thermo_t(itemp), p3)
      CALL clean_poly(p3)
   ELSEIF (poly_degree_thermo == 2) THEN
!
!   compute the coefficients of the quadratic polynomial
!
      CALL init_poly(nvar,p2)
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f, p2)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_quadratic(nvar, x_pos_min, thermo_t(itemp), p2)
      CALL clean_poly(p2)
   ELSEIF (poly_degree_thermo == 1) THEN
!
!   compute the coefficients of the linear polynomial
!
      CALL init_poly(nvar,p1)
      CALL fit_multi_linear(ndata, nvar, lsolve, x, f, p1)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_linear(nvar, x_pos_min, thermo_t(itemp), p1)
      CALL clean_poly(p1)

   ENDIF
ENDDO
!
!  collect all temperatures on all processors
!
CALL mp_sum(thermo_t, world_comm)

DEALLOCATE(f)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_thermo

!-----------------------------------------------------------------------
SUBROUTINE interpolate_thermo_p(vmin_p, ph_thermo, thermo_p, itemp)
!-----------------------------------------------------------------------
!
!  This subroutine receives the desired thermodynamical quantity at 
!  several geometries (specific heat, entropy or vibrational energy), 
!  the equilibrium volume as a function of pressure, 
!  and interpolates it at the equilibrium volume at each 
!  pressure. The routine is parallelized over pressures and must 
!  be called by all processors.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE thermo_mod,     ONLY : celldm_geo, omega_geo, no_ph, tot_ngeo
USE control_mur,    ONLY : lmurn
USE temperature,    ONLY : ntemp
USE control_pressure, ONLY : npress
USE linear_surfaces, ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_quartic_energy, ONLY : poly_degree_thermo, lsolve
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE polynomial,     ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
REAL(DP), INTENT(IN)  :: vmin_p(npress), ph_thermo(ntemp, tot_ngeo)

REAL(DP), INTENT(OUT) :: thermo_p(npress)

INTEGER :: igeo, ipress, nvar, ndata, startp, lastp
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:)
TYPE(poly1) :: p1
TYPE(poly2) :: p2
TYPE(poly3) :: p3
TYPE(poly4) :: p4

nvar=crystal_parameters(ibrav)
IF (lmurn) nvar=1         ! only the volume is variable in this case
!
ndata = compute_nwork_ph(no_ph,tot_ngeo)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
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
      CALL errore('intepolate_thermo_p','Not implemented',1)
   ENDIF
ENDDO
!
!  Divide pressures among processors
!
CALL divide(world_comm, npress, startp, lastp)
thermo_p=0.0_DP
DO ipress=startp,lastp
!
!  collect the computed data at this temperature
!
   ndata=0
   DO igeo=1, tot_ngeo
      IF (no_ph(igeo)) CYCLE
      ndata=ndata+1
      f(ndata)= ph_thermo(itemp,igeo)
   ENDDO
!
!  compute the geometry corresponding to this temperature.
!
   IF (lmurn) THEN
      x_pos_min(1)=vmin_p(ipress)
   ELSE
   ENDIF

   IF (poly_degree_thermo == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
      CALL init_poly(nvar,p4)
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4)
!
!  and evaluate the polynomial
!

      CALL evaluate_fit_quartic(nvar, x_pos_min, thermo_p(ipress), p4)
      CALL clean_poly(p4)
   ELSEIF (poly_degree_thermo == 3) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p3)
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, f, p3)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_cubic(nvar, x_pos_min, thermo_p(ipress), p3)
      CALL clean_poly(p3)
   ELSEIF (poly_degree_thermo == 2) THEN
!
!   compute the coefficients of the quadratic polynomial
!
      CALL init_poly(nvar,p2)
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f, p2)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_quadratic(nvar, x_pos_min, thermo_p(ipress), p2)
      CALL clean_poly(p2)
   ELSEIF (poly_degree_thermo == 1) THEN
!
!   compute the coefficients of the linear polynomial
!
      CALL init_poly(nvar,p1)
      CALL fit_multi_linear(ndata, nvar, lsolve, x, f, p1)
!
!  and evaluate the polynomial
!
      CALL evaluate_fit_linear(nvar, x_pos_min, thermo_p(ipress), p1)
      CALL clean_poly(p1)

   ENDIF
ENDDO
!
!  collect all temperatures on all processors
!
CALL mp_sum(thermo_p, world_comm)

DEALLOCATE(f)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_thermo_p

! Copyright (C) 2019 Cristiano Malica

!-----------------------------------------------------------------------
SUBROUTINE interpolate_e0(vmin_t, celldm_t, ph_e0, e0)
!-----------------------------------------------------------------------
!
!  This subroutine receives the zero point energy (ZPE) computed  
!  at several geometries in ph_e0.
!  Then it interpolates the ZPE at the equilibrium volume at T=0 K 
!  (or, more precisely, at the minimum temperature considered 
!  in the calculation) and puts it into e0. 
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE phdos_module,   ONLY : zero_point_energy
USE ph_freq_module, ONLY : zero_point_energy_ph
USE thermo_mod,     ONLY : celldm_geo, omega_geo, no_ph, tot_ngeo
USE control_mur,    ONLY : lmurn
USE temperature,    ONLY : ntemp
USE linear_surfaces, ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_quartic_energy, ONLY : poly_degree_thermo, lsolve
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE polynomial,     ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: vmin_t(ntemp), celldm_t(6,ntemp), ph_e0(tot_ngeo)
REAL(DP), INTENT(OUT) :: e0
INTEGER :: itemp, igeo, nvar, ndata, startt, lastt
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:)
TYPE(poly1) :: p1
TYPE(poly2) :: p2
TYPE(poly3) :: p3
TYPE(poly4) :: p4

nvar=crystal_parameters(ibrav)
IF (lmurn) nvar=1         ! only the volume is variable in this case
                          ! also in the anisotropic case

ndata = compute_nwork_ph(no_ph,tot_ngeo)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
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
      CALL compress_celldm(celldm_geo(1,igeo), x(1,ndata), nvar, ibrav)
   ENDIF
ENDDO

e0=0.0_DP
!
!  collect the computed data at each geometry
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   ndata=ndata+1
   f(ndata)=ph_e0(igeo)
ENDDO
!
!  compute the geometry corresponding to the zero temperature
!
IF (lmurn) THEN
   x_pos_min(1)=vmin_t(1)
ELSE
   CALL compress_celldm(celldm_t(1,1), x_pos_min, nvar, ibrav)
ENDIF

IF (poly_degree_thermo == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
   CALL init_poly(nvar,p4)
   CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4)
!
!  and evaluate the polynomial
!
   CALL evaluate_fit_quartic(nvar, x_pos_min, e0, p4)
   CALL clean_poly(p4)

ELSEIF (poly_degree_thermo == 3) THEN
!
!  compute the coefficients of the cubic polynomial
!
   CALL init_poly(nvar,p3)
   CALL fit_multi_cubic(ndata, nvar, lsolve, x, f, p3)
!
!  and evaluate the polynomial
!
   CALL evaluate_fit_cubic(nvar, x_pos_min, e0, p3)
   CALL clean_poly(p3)
   
ELSEIF (poly_degree_thermo == 2) THEN
!
!  compute the coefficients of the quadratic polynomial
!
   CALL init_poly(nvar,p2)
   CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f, p2)
!
!  and evaluate the polynomial
!
   CALL evaluate_fit_quadratic(nvar, x_pos_min, e0, p2)
   CALL clean_poly(p2)

ELSEIF (poly_degree_thermo == 1) THEN
!
!  compute the coefficients of the linear polynomial
!
   CALL init_poly(nvar,p1)
   CALL fit_multi_linear(ndata, nvar, lsolve, x, f, p1)
!
!  and evaluate the polynomial
!
   CALL evaluate_fit_linear(nvar, x_pos_min, e0, p1)
   CALL clean_poly(p1)
ENDIF

DEALLOCATE(f)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_e0
