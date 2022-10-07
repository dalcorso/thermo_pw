!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit()
  !-----------------------------------------------------------------------
  !
  !   This routine receives the total energy for several values of 
  !   celldm and fits them with a quadratic polynomial of dimension 
  !   equal to the number of indipendent parameters in celldm. 
  !   It finds the minimum of the quadratic function.
  !   If lquartic is true it fits the data with a quartic polynomial 
  !   and starting from the quadratic minimum it finds the minimum of
  !   the quartic polynomial.
  !   
  !   the output of this routine is celldm0(:) and emin the energy at
  !   the mimimum.
  !
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : ibrav
  USE control_mur,  ONLY : emin
  USE equilibrium_conf, ONLY : celldm0
  USE thermo_mod,   ONLY : celldm_geo, omega_geo, energy_geo
  USE control_pressure, ONLY : pressure, pressure_kb
  USE control_quadratic_energy, ONLY : hessian_v, hessian_e, x_pos_min, &
                                       p2, nvar
  USE control_quartic_energy, ONLY : p4, x_min_4, lquartic, lsolve,     &
                                     hessian4_v, hessian4_e
  USE lattices,     ONLY : expand_celldm, crystal_parameters
  USE quadratic_surfaces, ONLY : fit_multi_quadratic,  &
                          find_quadratic_extremum,     &
                          write_quadratic_hessian,     &
                          print_quadratic_polynomial,  &
                          summarize_fitting_data,      &
                          compare_quadratic_fit,       &
                          introduce_quadratic_fit, print_chisq_quadratic
  USE quartic_surfaces, ONLY : fit_multi_quartic, &
                          find_quartic_extremum, print_quartic_polynomial, &
                          print_chisq_quartic, introduce_quartic_fit,      &
                          write_quartic_hessian
  USE polynomial,   ONLY : init_poly
  USE vector_mod,       ONLY : write_vector
  USE io_global,    ONLY : stdout
  IMPLICIT NONE

  INTEGER  :: tncoeff4, ndata, ncoeff
  REAL(DP), ALLOCATABLE :: x(:,:), y(:), f(:)
  REAL(DP) :: ymin, ymin4
  INTEGER  :: compute_nwork
  !
  ! Only the first image does the calculation
  !
  celldm0(:)=0.0_DP
  emin=0.0_DP
  !
  WRITE(stdout,'(/,5x,71("-"))')
  !
  nvar=crystal_parameters(ibrav)
  !
  ndata = compute_nwork()

  IF (pressure>0.0_DP) THEN
     WRITE(stdout,'(/,5x,"Fitting the enthalpy with a function of the crystal &
                                         &parameters.")') 
     WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  ELSE
     WRITE(stdout,'(/,5x,"Fitting the energy with a function of the crystal &
                                         &parameters.")') 
  ENDIF


  ALLOCATE(x(nvar,ndata))
  ALLOCATE(x_pos_min(nvar))
  ALLOCATE(hessian_e(nvar))
  ALLOCATE(hessian_v(nvar, nvar))
  ALLOCATE(f(ndata))

  CALL init_poly(nvar,p2)
  ncoeff=1+nvar+p2%ncoeff2
  CALL introduce_quadratic_fit(nvar, ncoeff, ndata)

  f(:)=energy_geo(:) + pressure * omega_geo(:)

  CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)
  !
!  CALL summarize_fitting_data(nvar, ndata, x, f)

  CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,p2)
 
  CALL print_quadratic_polynomial(nvar, p2)

  IF (pressure>0.0_DP) THEN
     WRITE(stdout,'(7x,"Enthalpy (1)      Fitted enthalpy (2)   &
                                                    &DeltaE (1)-(2)")') 
  ELSE
     WRITE(stdout,'(7x,"Energy (1)      Fitted energy (2)   &
                                                    &DeltaE (1)-(2)")') 
  ENDIF
  CALL compare_quadratic_fit(ndata, nvar, x, f, p2)
  CALL print_chisq_quadratic(ndata, nvar, x, f, p2)

  CALL find_quadratic_extremum(nvar,x_pos_min,ymin,p2)
  !
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')

  CALL write_vector(nvar,x_pos_min)
  CALL print_energy(ymin)
  !
  CALL write_quadratic_hessian(nvar,p2,hessian_v,hessian_e)
  !
  CALL expand_celldm(celldm0, x_pos_min, nvar, ibrav)
  emin=ymin
  !
  IF (lquartic) THEN
     ALLOCATE(x_min_4(nvar))
     ALLOCATE(hessian4_e(nvar))
     ALLOCATE(hessian4_v(nvar, nvar))
     CALL init_poly(nvar,p4)
     tncoeff4=1+nvar+p4%ncoeff2+p4%ncoeff3+p4%ncoeff4
     CALL introduce_quartic_fit(nvar,tncoeff4,ndata)

     CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,p4)
     !
     CALL print_quartic_polynomial(nvar,p4) 
!    WRITE(stdout,'(/,7x,"Energy (1)    Fitted energy (2)   DeltaE (1)-(2)")') 
     CALL print_chisq_quartic(ndata, nvar, x, f, p4)
!
!   searching the minimum starting from the minimum of the quadratic
!
     x_min_4=x_pos_min
     CALL find_quartic_extremum(nvar,x_min_4,ymin4,p4)

     WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     CALL write_vector(nvar,x_min_4)
     CALL print_energy(ymin4)
     CALL write_quartic_hessian(nvar, x_min_4, p4, hessian4_v, hessian4_e)
!
     CALL expand_celldm(celldm0, x_min_4, nvar, ibrav)
     emin=ymin4
  ENDIF

  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
END SUBROUTINE quadratic_fit

!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_t(itemp, celldm_t, free_e_min_t, ph_free_ener, &
                                  ndatatot )
  !-----------------------------------------------------------------------
  !
  !   This routine receives the free energy at temperature itemp
  !   for several values of celldm and fits it with a quadratic function 
  !   of dimension equal to the number of indipendent parameters in celldm. 
  !   Then it adds the quadratic polynomial that fits the energy
  !   (or enthalphy) and finds the minimum.
  !   If lquartic is true it interpolates the free energy with a polynomial 
  !   of degree poly_degree_ph, adds the quartic polynomial that fits 
  !   the energy (or enthalpy) and finds its minimum starting from
  !   the minimum of the quadratic polynomial.
  !
  !   The output of this routine is celldm_t at the given temperature itemp
  !   and the free energy at the minimum free_e_min_t
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE thermo_mod,  ONLY : celldm_geo, no_ph
  USE control_quadratic_energy, ONLY : nvar, enthalpy_p2 => p2
  USE control_quartic_energy, ONLY :  p4, lquartic, poly_degree_ph, lsolve
  USE lattices,    ONLY : compress_celldm, expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, &
                      find_two_quadratic_extremum,    &
                      print_quadratic_polynomial,     &
                      summarize_fitting_data,         &
                      introduce_quadratic_fit, print_chisq_quadratic

  USE cubic_surfaces, ONLY : fit_multi_cubic, print_chisq_cubic
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      fit_multi_quartic, find_two_quartic_extremum,     &
                      find_quartic_cubic_extremum, print_chisq_quartic, &
                      find_quartic_linear_extremum

  USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  INTEGER  :: itemp
  INTEGER  :: ndata, ndatatot 
  REAL(DP) :: ph_free_ener(ndatatot), celldm_t(6), free_e_min_t
  REAL(DP), ALLOCATABLE :: x(:,:), f(:), x_pos_min(:)
  TYPE(poly1) :: pt1
  TYPE(poly2) :: pt2
  TYPE(poly3) :: pt3             
  TYPE(poly4) :: pt4             
  REAL(DP) :: ymin
  INTEGER  :: idata, ncoeff
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ndata = compute_nwork_ph(no_ph,ndatatot)

  ALLOCATE(x(nvar,ndata))
  ALLOCATE(x_pos_min(nvar))
  ALLOCATE(f(ndata))
  CALL init_poly(nvar,pt2)

  ncoeff=1+nvar+pt2%ncoeff2
  IF (MOD(itemp-1,50)==0) &
     CALL introduce_quadratic_fit(nvar, ncoeff, ndata)

  ndata=0
  DO idata=1,ndatatot
     IF (no_ph(idata)) CYCLE
     ndata=ndata+1
     CALL compress_celldm(celldm_geo(1,idata), x(1,ndata), nvar, ibrav)
     f(ndata)=ph_free_ener(idata)
  END DO
  !
  !CALL summarize_fitting_data(nvar, ndata, x, f)
  !
  CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f, pt2)

  CALL print_quadratic_polynomial(nvar, pt2)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  CALL print_chisq_quadratic(ndata, nvar, x, f, pt2)

  CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, enthalpy_p2, pt2)
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
  CALL write_vector(nvar,x_pos_min)
  CALL print_genergy(ymin)

  IF (lquartic) THEN
     IF (poly_degree_ph==4) THEN
        WRITE(stdout,'(/,5x, "Fit improved with a fourth order polynomial")') 
        CALL init_poly(nvar,pt4)
        CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, pt4)
        CALL print_chisq_quartic(ndata, nvar, x, f, pt4)
        CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, p4, pt4) 
        CALL clean_poly(pt4)
        WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     ELSEIF (poly_degree_ph==3) THEN
        WRITE(stdout,'(/,5x, "Fit improved with a third order polynomial")') 
        CALL init_poly(nvar,pt3)
        CALL fit_multi_cubic(ndata, nvar, lsolve, x, f, pt3)
        CALL print_chisq_cubic(ndata, nvar, x, f, pt3)
        CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, p4, pt3)
        CALL clean_poly(pt3)
        WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
     ELSEIF (poly_degree_ph==1) THEN
        WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
        CALL init_poly(nvar,pt1)
        CALL fit_multi_linear(ndata, nvar, lsolve, x, f, pt1)
        CALL print_chisq_linear(ndata, nvar, x, f, pt1)
        CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, p4, pt1)
        CALL clean_poly(pt1)
        WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
     ELSE
        WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
        CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, p4, pt2)
     ENDIF
     CALL write_vector(nvar,x_pos_min)
     CALL print_genergy(ymin)
  ENDIF

  free_e_min_t=ymin
  CALL expand_celldm(celldm_t, x_pos_min, nvar, ibrav)

  CALL clean_poly(pt2)
  DEALLOCATE(f)
  DEALLOCATE(x_pos_min)
  DEALLOCATE(x)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t

!-------------------------------------------------------------------
SUBROUTINE set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)
!-------------------------------------------------------------------
!
!  this rouotine receives and array of values of celldm celldm_geo(6,ndata)
!  and transform it in a compact array x(nvar,ndata), where nvar depends
!  on the Bravais lattice
!
USE kinds, ONLY : DP
USE lattices, ONLY : compress_celldm

IMPLICIT NONE
INTEGER, INTENT(IN)     :: ibrav, nvar, ndata
REAL(DP), INTENT(IN)    :: celldm_geo(6,ndata)
REAL(DP), INTENT(INOUT) :: x(nvar,ndata)

INTEGER :: idata

DO idata=1,ndata
   CALL compress_celldm(celldm_geo(1,idata), x(1,idata), nvar, ibrav)
ENDDO

RETURN
END SUBROUTINE set_x_from_celldm

!-------------------------------------------------------------------
SUBROUTINE print_genergy(ymin)
!-------------------------------------------------------------------
USE kinds, ONLY : DP
USE control_pressure, ONLY : pressure
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) :: ymin

IF (pressure > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Gibbs energy at the extremum:",f22.12," Ry")') ymin
ELSE
   WRITE(stdout,'(5x,"Free energy at the extremum:",f22.12," Ry")') ymin
ENDIF

RETURN
END SUBROUTINE print_genergy

!-------------------------------------------------------------------
SUBROUTINE print_energy(ymin)
!-------------------------------------------------------------------

USE kinds, ONLY : DP
USE control_pressure, ONLY : pressure
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) :: ymin

IF (pressure > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Enthalpy at the extremum:",f21.12," Ry")') ymin
ELSE
   WRITE(stdout,'(5x,"Energy at the extremum:",f21.12," Ry")') ymin
ENDIF

RETURN
END SUBROUTINE print_energy
