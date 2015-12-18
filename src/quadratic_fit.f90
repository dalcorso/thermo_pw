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
  !   celldm and fits them with a quadratic function of dimension 
  !   equal to the number of indipendent parameters in celldm. 
  !   It finds also the minimum of the quadratic function.
  !   
  !   the output of this routine is celldm0(:)
  !
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : ibrav
  USE control_mur, ONLY : celldm0, emin
  USE mp_images, ONLY : my_image_id, root_image
  USE thermo_mod, ONLY : celldm_geo, energy_geo, omega_geo
  USE control_pressure, ONLY : pressure, pressure_kb
  USE control_quadratic_energy, ONLY : hessian_v, hessian_e, x_pos_min, &
                               coeff, degree
  USE io_global, ONLY : stdout
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_fit_extremum, &
                                 write_fit_hessian, evaluate_fit_quadratic, &
                                 print_quadratic_polynomial, &
                                 summarize_fitting_data, write_vector, &
                                 introduce_quadratic_fit
  IMPLICIT NONE
  INTEGER :: nvar, ndata
  REAL(DP), ALLOCATABLE :: x(:,:), y(:), f(:)
  REAL(DP) :: ymin, chisq, aux
  INTEGER :: idata
  INTEGER :: compute_nwork
  !
  ! Only the first image does the calculation
  !
  celldm0(:)=0.0_DP
  IF (my_image_id /= root_image) RETURN
  !
  CALL compute_degree(ibrav,degree,nvar)
  !
  ndata = compute_nwork()

  IF (pressure>0.0_DP) THEN
     WRITE(stdout,'(/,5x,"Enthalpy")') 
     WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  ELSE
     WRITE(stdout,'(/,5x,"Energy")') 
  ENDIF

  CALL introduce_quadratic_fit(degree, nvar, ndata)

  
  IF ( ALLOCATED(x_pos_min) ) DEALLOCATE(x_pos_min)
  IF ( ALLOCATED(hessian_e) ) DEALLOCATE(hessian_e)
  IF ( ALLOCATED(hessian_v) ) DEALLOCATE(hessian_v)
  IF ( ALLOCATED(coeff) ) DEALLOCATE(coeff)

  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(hessian_e(degree))
  ALLOCATE(hessian_v(degree, degree))
  ALLOCATE(f(ndata))
  ALLOCATE(coeff(nvar))

  DO idata=1,ndata
     f(idata)=energy_geo(idata) + pressure * omega_geo(idata)
  END DO

  CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_geo)
  !
  CALL summarize_fitting_data(degree, ndata, x, f)

  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)
  
  CALL print_quadratic_polynomial(degree, nvar, coeff)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  chisq=0.0_DP
  DO idata=1,ndata
     CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
     chisq = chisq + (aux - f(idata))**2
  ENDDO

  CALL find_fit_extremum(degree,nvar,x_pos_min,ymin,coeff)

  WRITE(stdout,'(/,5x,"Extremum found at:")')

  CALL write_vector(degree,x_pos_min)

  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x,"Enthalpy at the extremum",f18.12)') ymin
  ELSE
     WRITE(stdout,'(5x,"Energy at the extremum",f18.12)') ymin
  END IF
  WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

  CALL write_fit_hessian(degree,nvar,coeff,hessian_v,hessian_e)

  CALL set_celldm_from_xmin(ibrav, degree, x_pos_min, celldm0)

  emin=ymin

  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
END SUBROUTINE quadratic_fit

!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_t(itemp)
  !-----------------------------------------------------------------------
  !
  !   This routine receives the total free energy for several values of 
  !   celldm and fits them with a quadratic function of dimension 
  !   equal to the number of indipendent parameters in celldm. 
  !
  !   The output of this routine is celldm_t at the given temperature itemp
  !
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : ibrav
  USE mp_images, ONLY : my_image_id, root_image
  USE thermo_mod, ONLY : celldm_geo, energy_geo, no_ph
  USE control_quadratic_energy, ONLY : degree, nvar, coeff_t, &
                                       enthalpy_coeff => coeff 
  USE temperature, ONLY : temp
  USE control_pressure, ONLY : pressure, pressure_kb
  USE thermodynamics, ONLY : ph_free_ener
  USE anharmonic, ONLY : celldm_t
  USE io_global, ONLY : stdout
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_two_fit_extremum, &
                                 write_fit_hessian, evaluate_fit_quadratic,  &
                                 print_quadratic_polynomial, &
                                 summarize_fitting_data, write_vector, &
                                 introduce_quadratic_fit
  IMPLICIT NONE
  INTEGER :: itemp
  INTEGER :: ndata, ndatatot
  REAL(DP), ALLOCATABLE :: x(:,:), f(:), coeff(:), x_pos_min(:), &
            celldm_data(:,:)
  REAL(DP) :: ymin, chisq, aux
  INTEGER :: idata
  INTEGER :: compute_nwork, compute_nwork_ph
  !
  ! Only the first image does the calculation
  !
  IF (my_image_id /= root_image) RETURN
  !
  CALL compute_degree(ibrav, degree, nvar)
  !
  ndatatot= compute_nwork()
  ndata = compute_nwork_ph(no_ph,ndatatot)

  WRITE(stdout,'(/,5x,70("-"))')
  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x, "Gibbs energy from phdos, at T= ", f12.6)') temp(itemp)
     WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  ELSE
     WRITE(stdout,'(5x, "Helmholtz free energy from phdos, at T= ", f12.6)') &
                                                                 temp(itemp)
  ENDIF

  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(f(ndata))
  ALLOCATE(coeff(nvar))
  ALLOCATE(celldm_data(6, ndata))

  IF (MOD(itemp-1,50)==0) &
     CALL introduce_quadratic_fit(degree, nvar, ndata)

  ndata=0
  DO idata=1,ndatatot
     IF (.NOT.no_ph(idata)) THEN
        ndata=ndata+1
        celldm_data(:,ndata)=celldm_geo(:,idata)
        f(ndata)=ph_free_ener(itemp,idata)
     END IF
  END DO

  CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_data)
  !
  !CALL summarize_fitting_data(degree, ndata, x, f)
  !
  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)

  CALL print_quadratic_polynomial(degree, nvar, coeff)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  chisq=0.0_DP
  DO idata=1,ndata
     CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
     chisq = chisq + (aux - f(idata))**2
  ENDDO

  CALL find_two_fit_extremum(degree,nvar,x_pos_min,ymin,enthalpy_coeff,coeff)

  WRITE(stdout,'(/,5x,"Extremum found at:")')

  CALL write_vector(degree,x_pos_min)

  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x,"Gibbs energy at the extremum",f18.12)') ymin
  ELSE
     WRITE(stdout,'(5x,"Free energy at the extremum",f18.12)') ymin
  END IF
  WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

  CALL set_celldm_from_xmin(ibrav, degree, x_pos_min, celldm_t(1,itemp))

  coeff_t(1:nvar,itemp) = coeff(1:nvar)

  DEALLOCATE(x_pos_min)
  DEALLOCATE(coeff)
  DEALLOCATE(celldm_data)
  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_t_ph(itemp)
  !-----------------------------------------------------------------------
  !
  !   This routine receives the total free energy for several values of 
  !   celldm and fits them with a quadratic function of dimension 
  !   equal to the number of indipendent parameters in celldm. 
  !
  !   The output of this routine is celldmf_t at the given temperature itemp
  !
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : ibrav
  USE mp_images, ONLY : my_image_id, root_image
  USE thermo_mod, ONLY : celldm_geo, energy_geo, no_ph
  USE control_quadratic_energy, ONLY : degree, nvar, coeff_t, &
                         enthalpy_coeff => coeff
  USE temperature, ONLY : temp
  USE control_pressure, ONLY : pressure, pressure_kb
  USE ph_freq_thermodynamics, ONLY : phf_free_ener
  USE ph_freq_anharmonic, ONLY : celldmf_t
  USE io_global, ONLY : stdout
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_two_fit_extremum, &
                                 write_fit_hessian, evaluate_fit_quadratic,  &
                                 print_quadratic_polynomial, &
                                 summarize_fitting_data, write_vector, &
                                 introduce_quadratic_fit
  IMPLICIT NONE
  INTEGER :: itemp
  INTEGER :: ndata, ndatatot
  REAL(DP), ALLOCATABLE :: x(:,:), f(:), coeff(:), x_pos_min(:), &
                           celldm_data(:,:)
  REAL(DP) :: ymin, chisq, aux
  INTEGER :: idata
  INTEGER :: compute_nwork, compute_nwork_ph
  !
  ! Only the first image does the calculation
  !
  IF (my_image_id /= root_image) RETURN
  !
  CALL compute_degree(ibrav, degree, nvar)
  !
  ndatatot= compute_nwork()
  ndata = compute_nwork_ph(no_ph,ndatatot)

  WRITE(stdout,'(/,5x,70("+"))')
  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x, "Gibbs energy from integration, at T= ", f12.6)') &
                                                                   temp(itemp)
     WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  ELSE
     WRITE(stdout,'(5x, "Helmholtz Free energy from integration, at T= ", &
                                                      &f12.6)') temp(itemp)
  ENDIF

  IF (MOD(itemp-1,50)==0) &
     CALL introduce_quadratic_fit(degree, nvar, ndata)

  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(f(ndata))
  ALLOCATE(coeff(nvar))
  ALLOCATE(celldm_data(6,ndata))

  ndata=0
  DO idata=1,ndatatot
     IF (.NOT. no_ph(idata)) THEN
        ndata=ndata+1
        celldm_data(:,ndata)=celldm_geo(:,idata)
        f(ndata)=phf_free_ener(itemp,idata)
     END IF
  END DO

  CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_data)
  !
  !CALL summarize_fitting_data(degree, ndata, x, f)
  !
  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)

  CALL print_quadratic_polynomial(degree, nvar, coeff)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  chisq=0.0_DP
  DO idata=1,ndata
     CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
     chisq = chisq + (aux - f(idata))**2
  ENDDO

  CALL find_two_fit_extremum(degree,nvar,x_pos_min,ymin,enthalpy_coeff,coeff)

  WRITE(stdout,'(/,5x,"Extremum found at:")')
  CALL write_vector(degree,x_pos_min)

  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x,"Gibbs energy at the extremum",f18.12)') ymin
  ELSE
     WRITE(stdout,'(5x,"Free energy at the extremum",f18.12)') ymin
  END IF
  WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

  CALL set_celldm_from_xmin(ibrav, degree, x_pos_min, celldmf_t(1,itemp))

  DEALLOCATE(x_pos_min)
  DEALLOCATE(coeff)
  DEALLOCATE(celldm_data)
  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t_ph

SUBROUTINE compute_degree(ibrav, degree, nvar)
!
!  number of coefficients of the quadratic equation
!  degrees       nvar
!  1               3, 
!  2               6, 
!  3              10, 
!  4              15, 
!  5              21, 
!  6              28
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
INTEGER, INTENT(OUT) :: degree, nvar

 SELECT CASE (ibrav)
    CASE(1,2,3)
       degree=1
    CASE(4,5,6,7)
       degree=2
    CASE(8,9,91,10,11)
       degree=3
    CASE(12,-12,13,-13)
       degree=4
    CASE DEFAULT
       degree=6
 END SELECT
 nvar = (1 + degree)*(degree + 2) / 2
 RETURN
END SUBROUTINE compute_degree

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
USE kinds, ONLY : DP

IMPLICIT NONE

INTEGER, INTENT(IN) :: ibrav, degree, ndata
REAL(DP), INTENT(IN) :: celldm_geo(6,ndata)
REAL(DP), INTENT(INOUT) :: x(degree,ndata)

INTEGER :: idata, pdata

  SELECT CASE (ibrav)
     CASE(1,2,3) 
       DO idata=1,ndata
          x(1,idata)=celldm_geo(1,idata)
       ENDDO 
     CASE(4,5,6,7)
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           IF (ibrav==5) THEN
              x(2,idata)=ACOS(celldm_geo(4,idata))
           ELSE
              x(2,idata)=celldm_geo(3,idata)
           ENDIF
        END DO 
     CASE(8,9,91,10,11)
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           x(2,idata)=celldm_geo(2,idata)
           x(3,idata)=celldm_geo(3,idata)
        ENDDO
     CASE(12,-12,13,-13) 
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           x(2,idata)=celldm_geo(2,idata)
           x(3,idata)=celldm_geo(3,idata)
           IF (ibrav>0) THEN

!   c unique
!
              x(4,idata)=ACOS(celldm_geo(4,idata))
           ELSE
!
!   b unique
!
              x(4,idata)=ACOS(celldm_geo(5,idata))
           ENDIF
        ENDDO
     CASE DEFAULT
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           x(2,idata)=celldm_geo(2,idata)
           x(3,idata)=celldm_geo(3,idata)
           x(4,idata)=ACOS(celldm_geo(4,idata))
           x(5,idata)=ACOS(celldm_geo(5,idata))
           x(6,idata)=ACOS(celldm_geo(6,idata))
        ENDDO
  END SELECT

RETURN
END SUBROUTINE set_x_from_celldm

SUBROUTINE set_celldm_from_xmin(ibrav, degree, x, celldm)
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, degree
REAL(DP), INTENT(IN) :: x(degree)
REAL(DP), INTENT(INOUT) :: celldm(6)

  celldm=0.0_DP
  SELECT CASE (ibrav)
     CASE(1,2,3) 
        celldm(1)=x(1)
     CASE(4,5,6,7)
        celldm(1)=x(1)
        IF (ibrav==5) THEN
           celldm(4)=COS(x(2))
        ELSE
           celldm(3)= x(2)
        ENDIF
     CASE(8,9,91,10,11)
        celldm(1)=x(1)
        celldm(2)=x(2)
        celldm(3)=x(3)
     CASE(12,-12,13,-13) 
        celldm(1)=x(1)
        celldm(2)=x(2)
        celldm(3)=x(3)
        IF (ibrav>0) THEN
!
!   c unique
!
           celldm(4)=COS(x(4))
        ELSE
!
!   b unique
!
           celldm(5)=COS(x(4))
        ENDIF
     CASE DEFAULT
        celldm(1)=x(1)
        celldm(2)=x(2)
        celldm(3)=x(3)
        celldm(4)=COS(x(4))
        celldm(5)=COS(x(5))
        celldm(6)=COS(x(6))
  END SELECT

RETURN
END SUBROUTINE set_celldm_from_xmin
