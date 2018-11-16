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
                               coeff, degree
  USE control_quartic_energy, ONLY : nvar4, coeff4, x_min_4, lquartic, lsolve
  USE lattices,     ONLY : expand_celldm, crystal_parameters
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_fit_extremum,  &
                          write_fit_hessian, print_quadratic_polynomial,  &
                          summarize_fitting_data, write_vector,           &
                          introduce_quadratic_fit, print_chisq_quadratic, &
                          quadratic_var
  USE quartic_surfaces, ONLY : fit_multi_quartic, quartic_var, &
                          find_quartic_extremum, print_quartic_polynomial, &
                          print_chisq_quartic, introduce_quartic_fit
  USE io_global,    ONLY : stdout
  IMPLICIT NONE

  INTEGER  :: nvar, ndata
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
  degree=crystal_parameters(ibrav)
  nvar=quadratic_var(degree)
  !
  ndata = compute_nwork()

  IF (pressure>0.0_DP) THEN
     WRITE(stdout,'(/,5x,"Enthalpy")') 
     WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  ELSE
     WRITE(stdout,'(/,5x,"Energy")') 
  ENDIF

  CALL introduce_quadratic_fit(degree, nvar, ndata)

  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(hessian_e(degree))
  ALLOCATE(hessian_v(degree, degree))
  ALLOCATE(f(ndata))
  ALLOCATE(coeff(nvar))

  f(:)=energy_geo(:) + pressure * omega_geo(:)

  CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_geo)
  !
  CALL summarize_fitting_data(degree, ndata, x, f)

  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)
  
  CALL print_quadratic_polynomial(degree, nvar, coeff)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  CALL print_chisq_quadratic(ndata, degree, nvar, x, f, coeff)

  CALL find_fit_extremum(degree,nvar,x_pos_min,ymin,coeff)
  !
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')

  CALL write_vector(degree,x_pos_min)
  CALL print_energy(ymin)
  !
  CALL write_fit_hessian(degree,nvar,coeff,hessian_v,hessian_e)
  !
  CALL expand_celldm(celldm0, x_pos_min, degree, ibrav)
  emin=ymin
  !
  IF (lquartic) THEN

     nvar4=quartic_var(degree)
     CALL introduce_quartic_fit(degree, nvar4, ndata)

     ALLOCATE(x_min_4(degree))
     ALLOCATE(coeff4(nvar4))

     CALL fit_multi_quartic(ndata,degree,nvar4,lsolve,x,f,coeff4)
     !
     CALL print_quartic_polynomial(degree, nvar4, coeff4)
!    WRITE(stdout,'(/,7x,"Energy (1)    Fitted energy (2)   DeltaE (1)-(2)")') 
     CALL print_chisq_quartic(ndata, degree, nvar, x, f, coeff4)
!
!   searching the minimum starting from the minimum of the quadratic
!
     x_min_4=x_pos_min
     CALL find_quartic_extremum(degree,nvar4,x_min_4,ymin4,coeff4)

     WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     CALL write_vector(degree,x_min_4)
     CALL print_energy(ymin4)
!
     CALL expand_celldm(celldm0, x_min_4, degree, ibrav)
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
  !   This routine receives the total free energy for several values of 
  !   celldm and fits them with a quadratic function of dimension 
  !   equal to the number of indipendent parameters in celldm. 
  !   Then it adds the quadratic polynomium that fits the energy
  !   (or enthalphy) and finds the minimum.
  !   If lquartic and lquartic_ph are both true it interpolates the
  !   data with a quartic polynomium, adds the quartic polynomium that
  !   fits the energy (or enthalphy) and finds its minimum starting from
  !   the minimum of the quadratic polynomium.
  !
  !   The output of this routine is celldm_t at the given temperature itemp
  !   and the free energy at the minimum free_e_min_t
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE thermo_mod,  ONLY : celldm_geo, no_ph
  USE control_quadratic_energy, ONLY : degree, nvar, enthalpy_coeff => coeff
  USE control_quartic_energy, ONLY :  nvar4, coeff4, lquartic, lquartic_ph, &
                                      lsolve
  USE lattices,    ONLY : compress_celldm, expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout

  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_two_fit_extremum, &
                      print_quadratic_polynomial, quadratic_var, &
                      summarize_fitting_data, write_vector, &
                      introduce_quadratic_fit, print_chisq_quadratic

  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum, &
                      fit_multi_quartic, find_two_quartic_extremum,  &
                      print_chisq_two_quartic, print_chisq_quartic_quadratic

  IMPLICIT NONE
  INTEGER  :: itemp
  INTEGER  :: ndata, ndatatot
  REAL(DP) :: ph_free_ener(ndatatot), celldm_t(6), free_e_min_t
  REAL(DP), ALLOCATABLE :: x(:,:), f(:), coeff(:), x_pos_min(:), coefft4(:)
  REAL(DP) :: ymin
  INTEGER  :: idata
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  degree=crystal_parameters(ibrav)
  nvar=quadratic_var(degree)
  !
  ndata = compute_nwork_ph(no_ph,ndatatot)

  IF (MOD(itemp-1,50)==0) &
     CALL introduce_quadratic_fit(degree, nvar, ndata)

  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(f(ndata))
  ALLOCATE(coeff(nvar))
  ALLOCATE(coefft4(nvar4))

  ndata=0
  DO idata=1,ndatatot
     IF (.NOT. no_ph(idata)) THEN
        ndata=ndata+1
        CALL compress_celldm(celldm_geo(1,idata), x(1,ndata), degree, ibrav)
        f(ndata)=ph_free_ener(idata)
     END IF
  END DO
  !
  !CALL summarize_fitting_data(degree, ndata, x, f)
  !
  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)

  CALL print_quadratic_polynomial(degree, nvar, coeff)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  CALL print_chisq_quadratic(ndata, degree, nvar, x, f, coeff)

  CALL find_two_fit_extremum(degree,nvar,x_pos_min,ymin,enthalpy_coeff,coeff)
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
  CALL write_vector(degree,x_pos_min)
  CALL print_genergy(ymin)

  IF (lquartic) THEN
     IF (lquartic_ph) THEN
        WRITE(stdout,'(/,5x, "Fit improved with a fourth order polynomial")') 
        CALL fit_multi_quartic(ndata,degree,nvar4,lsolve,x,f,coefft4)
        CALL print_chisq_two_quartic(ndata, degree, nvar4, x, f, &
                                                           coeff4, coefft4)
        CALL find_two_quartic_extremum(degree,nvar4,x_pos_min,&
                                                       ymin,coeff4,coefft4)
     ELSE
        WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
        CALL print_chisq_quartic_quadratic(ndata, degree, nvar4, nvar, &
                                                x, f, coeff4, coeff)
        CALL find_quartic_quadratic_extremum(degree,nvar4,nvar,x_pos_min,&
                                                           ymin,coeff4,coeff)
     ENDIF
     WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     CALL write_vector(degree,x_pos_min)
     CALL print_genergy(ymin)
  ENDIF

  free_e_min_t=ymin
  CALL expand_celldm(celldm_t, x_pos_min, degree, ibrav)

  DEALLOCATE(x_pos_min)
  DEALLOCATE(coefft4)
  DEALLOCATE(coeff)
  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t

INTEGER FUNCTION compute_nwork_ph(no_ph,ndatatot)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndatatot
LOGICAL, INTENT(IN) :: no_ph(ndatatot)

INTEGER :: idata, counter_ndata

counter_ndata=0
DO idata=1,ndatatot
   IF (.NOT. no_ph(idata)) counter_ndata=counter_ndata+1
ENDDO
compute_nwork_ph=counter_ndata

RETURN
END FUNCTION compute_nwork_ph

SUBROUTINE set_x_from_celldm(ibrav, degree, ndata, x, celldm_geo)
!
!  this rouotine receives and array of values of celldm celldm_geo(6,ndata)
!  and transform it in a compact array x(degree,ndata), where degree depends
!  on the Bravais lattice
!
USE kinds, ONLY : DP
USE lattices, ONLY : compress_celldm

IMPLICIT NONE
INTEGER, INTENT(IN)     :: ibrav, degree, ndata
REAL(DP), INTENT(IN)    :: celldm_geo(6,ndata)
REAL(DP), INTENT(INOUT) :: x(degree,ndata)

INTEGER :: idata

DO idata=1,ndata
   CALL compress_celldm(celldm_geo(1,idata), x(1,idata), degree, ibrav)
ENDDO

RETURN
END SUBROUTINE set_x_from_celldm

SUBROUTINE print_genergy(ymin)
USE kinds, ONLY : DP
USE control_pressure, ONLY : pressure
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) :: ymin

IF (pressure > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Gibbs energy at the extremum",f22.12)') ymin
ELSE
   WRITE(stdout,'(5x,"Free energy at the extremum",f22.12)') ymin
ENDIF

RETURN
END SUBROUTINE print_genergy

SUBROUTINE print_energy(ymin)

USE kinds, ONLY : DP
USE control_pressure, ONLY : pressure
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) :: ymin

IF (pressure > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Enthalpy at the extremum",f18.12)') ymin
ELSE
   WRITE(stdout,'(5x,"Energy at the extremum",f18.12)') ymin
ENDIF

RETURN
END SUBROUTINE print_energy
