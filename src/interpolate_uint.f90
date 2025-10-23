!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_uint(celldm_t, uint_geo_eos, uint_t, tau_t)
!-----------------------------------------------------------------------
!
!  This subroutine receives the internal parameters that at each
!  set of external parameters and interpolates them at the celldm_t.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : celldm_geo_eos, omega_geo_eos, no_ph_eos, &
                           tot_ngeo_eos
USE temperature,    ONLY : ntemp
USE control_atomic_pos, ONLY : nint_var, iconstr_internal, linterpolate_tau
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
REAL(DP), INTENT(IN)  :: uint_geo_eos(nint_var, tot_ngeo_eos), &
                         celldm_t(6,ntemp)

REAL(DP), INTENT(OUT) :: uint_t(nint_var, ntemp), tau_t(3,nat,ntemp)

INTEGER :: itemp, igeo, nvar, ndata, startt, lastt, ivar
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:)
TYPE(poly1) :: p1
TYPE(poly2) :: p2
TYPE(poly3) :: p3
TYPE(poly4) :: p4

IF (.NOT.linterpolate_tau) RETURN
nvar=crystal_parameters(ibrav)
!
ndata = compute_nwork_ph(no_ph_eos,tot_ngeo_eos)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
!
!  collect the geometrical data of the geometries for which phonon dispersions
!  have been calculated
!
ndata=0
DO igeo=1, tot_ngeo_eos
   IF (no_ph_eos(igeo)) CYCLE
   ndata=ndata+1
   CALL compress_celldm(celldm_geo_eos(1,igeo), x(1,ndata), nvar, ibrav)
ENDDO
!
!  Divide temperatures among processors
!
CALL divide(world_comm, ntemp, startt, lastt)
uint_t=0.0_DP
DO ivar=1, nint_var
!
!  collect the computed data
!
   ndata=0
   DO igeo=1, tot_ngeo_eos
      IF (no_ph_eos(igeo)) CYCLE
      ndata=ndata+1
      f(ndata)= uint_geo_eos(ivar,igeo)
   ENDDO
   DO itemp=startt,lastt
!
!  compute the geometry corresponding to this temperature.
!
      CALL compress_celldm(celldm_t(1,itemp), x_pos_min, nvar, ibrav)

      IF (poly_degree_thermo == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
         CALL init_poly(nvar,p4)
         CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4)
!
!  and evaluate the polynomial
!
         CALL evaluate_fit_quartic(nvar, x_pos_min, uint_t(ivar,itemp), p4)
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
         CALL evaluate_fit_cubic(nvar, x_pos_min, uint_t(ivar,itemp), p3)
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
         CALL evaluate_fit_quadratic(nvar, x_pos_min, uint_t(ivar,itemp), &
                                                                          p2)
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
         CALL evaluate_fit_linear(nvar, x_pos_min, uint_t(ivar,itemp), p1)
         CALL clean_poly(p1)

      ENDIF
   ENDDO
ENDDO
!
!  collect all temperatures on all processors
!
CALL mp_sum(uint_t, world_comm)

DO itemp=1, ntemp
   CALL internal_to_tau(celldm_t(1,itemp), tau_t(1,1,itemp),             &
                    uint_t(1,itemp), nat, nint_var, iconstr_internal, 1)
ENDDO

DEALLOCATE(f)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_uint
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_uint_zptt(celldm_p, uint_geo_eos, uint_ptt, tau_ptt) 
                                                            
!-----------------------------------------------------------------------
!
!  This subroutine receives the internal parameters that at each
!  set of external parameters obtained from the minimization of the 
!  energy and interpolates them at the celldm_p.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : celldm_geo_eos, omega_geo_eos, no_ph_eos, &
                           tot_ngeo_eos
USE temperature,    ONLY : ntemp
USE control_pressure,  ONLY : npress
USE control_atomic_pos, ONLY : nint_var, iconstr_internal, linterpolate_tau
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
REAL(DP), INTENT(IN)  :: uint_geo_eos(nint_var, tot_ngeo_eos), &
                         celldm_p(6,npress)

REAL(DP), INTENT(OUT) :: uint_ptt(nint_var, npress), tau_ptt(3,nat,npress)

INTEGER :: ipress, igeo, nvar, ndata, startp, lastp, ivar
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:)
TYPE(poly1) :: p1
TYPE(poly2) :: p2
TYPE(poly3) :: p3
TYPE(poly4) :: p4

IF (.NOT.linterpolate_tau) RETURN
nvar=crystal_parameters(ibrav)
!
ndata = compute_nwork_ph(no_ph_eos,tot_ngeo_eos)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
!
!  collect the geometrical data of the geometries for which phonon dispersions
!  have been calculated
!
ndata=0
DO igeo=1, tot_ngeo_eos
   IF (no_ph_eos(igeo)) CYCLE
   ndata=ndata+1
   CALL compress_celldm(celldm_geo_eos(1,igeo), x(1,ndata), nvar, ibrav)
ENDDO
!
!  Divide pressures among processors
!
CALL divide(world_comm, npress, startp, lastp)
uint_ptt=0.0_DP
DO ivar=1, nint_var
!
!  collect the computed data at this temperature
!
   ndata=0
   DO igeo=1, tot_ngeo_eos
      IF (no_ph_eos(igeo)) CYCLE
      ndata=ndata+1
      f(ndata)= uint_geo_eos(ivar,igeo)
   ENDDO
   DO ipress=startp,lastp
!
!  compute the geometry corresponding to this temperature.
!
      CALL compress_celldm(celldm_p(1,ipress), x_pos_min, nvar, ibrav)

      IF (poly_degree_thermo == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
         CALL init_poly(nvar,p4)
         CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4)
!
!  and evaluate the polynomial
!
         CALL evaluate_fit_quartic(nvar, x_pos_min, uint_ptt(ivar,ipress), p4)
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
         CALL evaluate_fit_cubic(nvar, x_pos_min, uint_ptt(ivar,ipress), p3)
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
         CALL evaluate_fit_quadratic(nvar, x_pos_min, uint_ptt(ivar,ipress), &
                                                                          p2)
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
         CALL evaluate_fit_linear(nvar, x_pos_min, uint_ptt(ivar,ipress), p1)
         CALL clean_poly(p1)
      ENDIF
   ENDDO
ENDDO
!
!  collect all temperatures on all processors
!
CALL mp_sum(uint_ptt, world_comm)

DO ipress=1, npress
   CALL internal_to_tau(celldm_p(1,ipress), tau_ptt(1,1,ipress),             &
                    uint_ptt(1,ipress), nat, nint_var, iconstr_internal, 1)
ENDDO

DEALLOCATE(f)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_uint_zptt
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_uint_t(celldm_t, uint_geo_eos_t, uint_t, tau_t)
!-----------------------------------------------------------------------
!
!  This subroutine receives the internal parameters that at each
!  set of external parameters minimize the free energy and interpolates
!  them at the celldm_t.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : celldm_geo_eos, omega_geo_eos, no_ph_eos, &
                           tot_ngeo_eos
USE temperature,    ONLY : ntemp
USE control_atomic_pos, ONLY : nint_var, iconstr_internal, linternal_thermo
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
REAL(DP), INTENT(IN)  :: uint_geo_eos_t(nint_var, ntemp, tot_ngeo_eos), &
                         celldm_t(6,ntemp)

REAL(DP), INTENT(OUT) :: uint_t(nint_var, ntemp), tau_t(3,nat,ntemp)

INTEGER :: itemp, igeo, nvar, ndata, startt, lastt, ivar
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:)
TYPE(poly1) :: p1
TYPE(poly2) :: p2
TYPE(poly3) :: p3
TYPE(poly4) :: p4

IF (.NOT.linternal_thermo) RETURN
nvar=crystal_parameters(ibrav)
!
ndata = compute_nwork_ph(no_ph_eos,tot_ngeo_eos)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
!
!  collect the geometrical data of the geometries for which phonon dispersions
!  have been calculated
!
ndata=0
DO igeo=1, tot_ngeo_eos
   IF (no_ph_eos(igeo)) CYCLE
   ndata=ndata+1
   CALL compress_celldm(celldm_geo_eos(1,igeo), x(1,ndata), nvar, ibrav)
ENDDO
!
!  Divide temperatures among processors
!
CALL divide(world_comm, ntemp, startt, lastt)
uint_t=0.0_DP
DO ivar=1, nint_var
   DO itemp=startt,lastt
!
!  collect the computed data at this temperature
!
      ndata=0
      DO igeo=1, tot_ngeo_eos
         IF (no_ph_eos(igeo)) CYCLE
         ndata=ndata+1
         f(ndata)= uint_geo_eos_t(ivar,itemp,igeo)
      ENDDO
!
!  compute the geometry corresponding to this temperature.
!
      CALL compress_celldm(celldm_t(1,itemp), x_pos_min, nvar, ibrav)

      IF (poly_degree_thermo == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
         CALL init_poly(nvar,p4)
         CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4)
!
!  and evaluate the polynomial
!
         CALL evaluate_fit_quartic(nvar, x_pos_min, uint_t(ivar,itemp), p4)
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
         CALL evaluate_fit_cubic(nvar, x_pos_min, uint_t(ivar,itemp), p3)
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
         CALL evaluate_fit_quadratic(nvar, x_pos_min, uint_t(ivar,itemp), &
                                                                          p2)
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
         CALL evaluate_fit_linear(nvar, x_pos_min, uint_t(ivar,itemp), p1)
         CALL clean_poly(p1)

      ENDIF
   ENDDO
ENDDO
!
!  collect all temperatures on all processors
!
CALL mp_sum(uint_t, world_comm)

DO itemp=1, ntemp
   CALL internal_to_tau(celldm_t(1,itemp), tau_t(1,1,itemp),             &
                    uint_t(1,itemp), nat, nint_var, iconstr_internal, 1)
ENDDO

DEALLOCATE(f)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_uint_t
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_uint_ptt(celldm_p, uint_geo_eos_t, uint_ptt, tau_ptt, &
                                uint_ptt_p1, tau_ptt_p1, uint_ptt_m1,    &
                                tau_ptt_m1, itemp)
!-----------------------------------------------------------------------
!
!  This subroutine receives the internal parameters that at each
!  set of external parameters minimize the free energy and interpolates
!  them at the celldm_p.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : celldm_geo_eos, omega_geo_eos, no_ph_eos, &
                           tot_ngeo_eos
USE temperature,    ONLY : ntemp
USE control_pressure,  ONLY : npress
USE control_atomic_pos, ONLY : nint_var, iconstr_internal, linternal_thermo
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
INTEGER, INTENT(IN)   :: itemp
REAL(DP), INTENT(IN)  :: uint_geo_eos_t(nint_var, ntemp, tot_ngeo_eos), &
                         celldm_p(6,npress)

REAL(DP), INTENT(INOUT) :: uint_ptt(nint_var, npress), tau_ptt(3,nat,npress)
REAL(DP), INTENT(INOUT) :: uint_ptt_p1(nint_var, npress), &
                                               tau_ptt_p1(3,nat,npress)
REAL(DP), INTENT(INOUT) :: uint_ptt_m1(nint_var, npress), &
                                               tau_ptt_m1(3,nat,npress)

INTEGER :: ipress, igeo, nvar, ndata, startp, lastp, ivar
INTEGER :: compute_nwork_ph
REAL(DP), ALLOCATABLE :: x(:,:), x_pos_min(:), f(:), f_p1(:), f_m1(:)
TYPE(poly1) :: p1, p1_p1, p1_m1
TYPE(poly2) :: p2, p2_p1, p2_m1
TYPE(poly3) :: p3, p3_p1, p3_m1
TYPE(poly4) :: p4, p4_p1, p4_m1

IF (.NOT.linternal_thermo) RETURN
nvar=crystal_parameters(ibrav)
!
ndata = compute_nwork_ph(no_ph_eos,tot_ngeo_eos)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
ALLOCATE(f_p1(ndata))
ALLOCATE(f_m1(ndata))
!
!  collect the geometrical data of the geometries for which phonon dispersions
!  have been calculated
!
ndata=0
DO igeo=1, tot_ngeo_eos
   IF (no_ph_eos(igeo)) CYCLE
   ndata=ndata+1
   CALL compress_celldm(celldm_geo_eos(1,igeo), x(1,ndata), nvar, ibrav)
ENDDO
!
!  Divide pressures among processors
!
CALL divide(world_comm, npress, startp, lastp)
uint_ptt=0.0_DP
uint_ptt_p1=0.0_DP
uint_ptt_m1=0.0_DP
DO ivar=1, nint_var
!
!  collect the computed data at this temperature
!
   ndata=0
   DO igeo=1, tot_ngeo_eos
      IF (no_ph_eos(igeo)) CYCLE
      ndata=ndata+1
      f(ndata)= uint_geo_eos_t(ivar,itemp,igeo)
      f_p1(ndata)=uint_geo_eos_t(ivar, itemp+1, igeo)
      f_m1(ndata)=uint_geo_eos_t(ivar, itemp-1, igeo)
   ENDDO
   DO ipress=startp,lastp
!
!  compute the geometry corresponding to this pressure.
!
      CALL compress_celldm(celldm_p(1,ipress), x_pos_min, nvar, ibrav)

      IF (poly_degree_thermo == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
         CALL init_poly(nvar,p4)
         CALL init_poly(nvar,p4_p1)
         CALL init_poly(nvar,p4_m1)
         CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, p4)
         CALL fit_multi_quartic(ndata, nvar, lsolve, x, f_p1, p4_p1)
         CALL fit_multi_quartic(ndata, nvar, lsolve, x, f_m1, p4_m1)
!
!  and evaluate the polynomial
!
         CALL evaluate_fit_quartic(nvar, x_pos_min, uint_ptt(ivar,ipress), p4)
         CALL evaluate_fit_quartic(nvar, x_pos_min, uint_ptt_p1(ivar,ipress),&
                                                                        p4_p1)
         CALL evaluate_fit_quartic(nvar, x_pos_min, uint_ptt_m1(ivar,ipress), &
                                                                        p4_m1)
         CALL clean_poly(p4)
         CALL clean_poly(p4_p1)
         CALL clean_poly(p4_m1)
      ELSEIF (poly_degree_thermo == 3) THEN
!
!   compute the coefficients of the cubic polynomial
!
         CALL init_poly(nvar,p3)
         CALL init_poly(nvar,p3_p1)
         CALL init_poly(nvar,p3_m1)
         CALL fit_multi_cubic(ndata, nvar, lsolve, x, f, p3)
         CALL fit_multi_cubic(ndata, nvar, lsolve, x, f_p1, p3_p1)
         CALL fit_multi_cubic(ndata, nvar, lsolve, x, f_m1, p3_m1)
!
!  and evaluate the polynomial
!
         CALL evaluate_fit_cubic(nvar, x_pos_min, uint_ptt(ivar,ipress), p3)
         CALL evaluate_fit_cubic(nvar, x_pos_min, uint_ptt_p1(ivar,ipress), &
                                                                p3_p1)
         CALL evaluate_fit_cubic(nvar, x_pos_min, uint_ptt_m1(ivar,ipress), &
                                                                p3_m1)
         CALL clean_poly(p3)
         CALL clean_poly(p3_p1)
         CALL clean_poly(p3_m1)
      ELSEIF (poly_degree_thermo == 2) THEN
!
!   compute the coefficients of the quadratic polynomial
!
         CALL init_poly(nvar,p2)
         CALL init_poly(nvar,p2_p1)
         CALL init_poly(nvar,p2_m1)
         CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f, p2)
         CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f_p1, p2_p1)
         CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f_m1, p2_m1)
!
!  and evaluate the polynomial
!
         CALL evaluate_fit_quadratic(nvar, x_pos_min,uint_ptt(ivar,ipress), &
                                                                          p2)
         CALL evaluate_fit_quadratic(nvar, x_pos_min,uint_ptt_p1(ivar,ipress),&
                                                                       p2_p1)
         CALL evaluate_fit_quadratic(nvar, x_pos_min,uint_ptt_m1(ivar,ipress),&
                                                                       p2_m1)
         CALL clean_poly(p2)
         CALL clean_poly(p2_p1)
         CALL clean_poly(p2_m1)
      ELSEIF (poly_degree_thermo == 1) THEN
!
!   compute the coefficients of the linear polynomial
!
         CALL init_poly(nvar,p1)
         CALL init_poly(nvar,p1_p1)
         CALL init_poly(nvar,p1_m1)
         CALL fit_multi_linear(ndata, nvar, lsolve, x, f, p1)
         CALL fit_multi_linear(ndata, nvar, lsolve, x, f_p1, p1_p1)
         CALL fit_multi_linear(ndata, nvar, lsolve, x, f_m1, p1_m1)
!
!  and evaluate the polynomial
!
         CALL evaluate_fit_linear(nvar, x_pos_min, uint_ptt(ivar,ipress), p1)
         CALL evaluate_fit_linear(nvar, x_pos_min, uint_ptt_p1(ivar,ipress), &
                                                             p1_p1)
         CALL evaluate_fit_linear(nvar, x_pos_min, uint_ptt_m1(ivar,ipress), &
                                                             p1_m1)
         CALL clean_poly(p1)
         CALL clean_poly(p1_p1)
         CALL clean_poly(p1_m1)
      ENDIF
   ENDDO
ENDDO
!
!  collect all temperatures on all processors
!
CALL mp_sum(uint_ptt, world_comm)
CALL mp_sum(uint_ptt_p1, world_comm)
CALL mp_sum(uint_ptt_m1, world_comm)

DO ipress=1, npress
   CALL internal_to_tau(celldm_p(1,ipress), tau_ptt(1,1,ipress),             &
                 uint_ptt(1,ipress), nat, nint_var, iconstr_internal, 1)
   CALL internal_to_tau(celldm_p(1,ipress), tau_ptt_p1(1,1,ipress),          &
                 uint_ptt_p1(1,ipress), nat, nint_var, iconstr_internal, 1)
   CALL internal_to_tau(celldm_p(1,ipress), tau_ptt_m1(1,1,ipress),          &
                 uint_ptt_m1(1,ipress), nat, nint_var, iconstr_internal, 1)
ENDDO

DEALLOCATE(f)
DEALLOCATE(f_p1)
DEALLOCATE(f_m1)
DEALLOCATE(x_pos_min)
DEALLOCATE(x)
!
RETURN
END SUBROUTINE interpolate_uint_ptt
!
