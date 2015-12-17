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
                                 summarize_fitting_data
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
  !  number of coefficients of the quadratic equation
  !  degrees       nvar
  !  1               3, 
  !  2               6, 
  !  3              10, 
  !  4              15, 
  !  5              21, 
  !  6              28
  !
  ndata = compute_nwork()

  WRITE(stdout,'(/,5x,"Fitting the energy with a quadratic function")') 
  WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  WRITE(stdout,'(5x,"Number of degrees of freedom: ",i5)')  degree
  WRITE(stdout,'(5x,"Coefficients of the quadratic:",i5)')  nvar
  WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata
  
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

  SELECT CASE (ibrav)
     CASE(1,2,3) 
       DO idata=1,ndata
          x(1,idata)=celldm_geo(1,idata)
          f(idata)=energy_geo(idata) + pressure * omega_geo(idata) 
       ENDDO 
     CASE(4,5,6,7)
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           IF (ibrav==5) THEN
              x(2,idata)=ACOS(celldm_geo(4,idata))
           ELSE
              x(2,idata)=celldm_geo(3,idata)
           ENDIF
           f(idata)=energy_geo(idata) + pressure * omega_geo(idata)
        END DO 
     CASE(8,9,91,10,11)
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           x(2,idata)=celldm_geo(2,idata)
           x(3,idata)=celldm_geo(3,idata)
           f(idata)=energy_geo(idata) + pressure * omega_geo(idata)
        ENDDO
     CASE(12,-12,13,-13) 
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           x(2,idata)=celldm_geo(2,idata)
           x(3,idata)=celldm_geo(3,idata)
           IF (ibrav>0) THEN
!
!   c unique
!
              x(4,idata)=ACOS(celldm_geo(4,idata))
           ELSE
!
!   b unique
!
              x(4,idata)=ACOS(celldm_geo(5,idata))
           ENDIF
           f(idata)=energy_geo(idata) + pressure * omega_geo(idata)
        ENDDO
     CASE DEFAULT
        DO idata=1,ndata
           x(1,idata)=celldm_geo(1,idata)
           x(2,idata)=celldm_geo(2,idata)
           x(3,idata)=celldm_geo(3,idata)
           x(4,idata)=ACOS(celldm_geo(4,idata))
           x(5,idata)=ACOS(celldm_geo(5,idata))
           x(6,idata)=ACOS(celldm_geo(6,idata))
           f(idata)=energy_geo(idata) + pressure * omega_geo(idata)
        ENDDO
  END SELECT
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
  WRITE(stdout,'(23x,"x1=",f12.6)') x_pos_min(1)
  IF (degree>1) WRITE(stdout,'(23x,"x2=",f12.6)') x_pos_min(2)
  IF (degree>2) WRITE(stdout,'(23x,"x3=",f12.6)') x_pos_min(3)
  IF (degree>3) WRITE(stdout,'(23x,"x4=",f12.6)') x_pos_min(4)
  IF (degree>4) WRITE(stdout,'(23x,"x5=",f12.6)') x_pos_min(5)
  IF (degree>5) WRITE(stdout,'(23x,"x6=",f12.6)') x_pos_min(6)
  WRITE(stdout,'(5x,"Energy at the extremum",f18.12)') ymin
  WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

  CALL write_fit_hessian(degree,nvar,coeff,hessian_v,hessian_e)

  SELECT CASE (ibrav)
     CASE(1,2,3) 
        celldm0(1)=x_pos_min(1)
     CASE(4,5,6,7)
        celldm0(1)=x_pos_min(1)
        IF (ibrav==5) THEN
           celldm0(4)=COS(x_pos_min(2))
        ELSE
           celldm0(3)= x_pos_min(2)
        ENDIF
     CASE(8,9,91,10,11)
        celldm0(1)=x_pos_min(1)
        celldm0(2)=x_pos_min(2)
        celldm0(3)=x_pos_min(3)
     CASE(12,-12,13,-13) 
        celldm0(1)=x_pos_min(1)
        celldm0(2)=x_pos_min(2)
        celldm0(3)=x_pos_min(3)
        IF (ibrav>0) THEN
!
!   c unique
!
           celldm0(4)=COS(x_pos_min(4))
        ELSE
!
!   b unique
!
           celldm0(5)=COS(x_pos_min(4))
        ENDIF
     CASE DEFAULT
        celldm0(1)=x_pos_min(1)
        celldm0(2)=x_pos_min(2)
        celldm0(3)=x_pos_min(3)
        celldm0(4)=COS(x_pos_min(4))
        celldm0(5)=COS(x_pos_min(5))
        celldm0(6)=COS(x_pos_min(6))
  END SELECT
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
  USE thermo_mod, ONLY : celldm_geo, energy_geo, omega_geo, no_ph
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
                                 summarize_fitting_data
  IMPLICIT NONE
  INTEGER :: itemp
  INTEGER :: ndata, ndatatot
  REAL(DP), ALLOCATABLE :: x(:,:), y(:), f(:), coeff(:), x_pos_min(:)
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
  !  number of coefficients of the quadratic equation
  !  degrees       nvar
  !  1               3, 
  !  2               6, 
  !  3              10, 
  !  4              15, 
  !  5              21, 
  !  6              28
  !
  ndatatot= compute_nwork()
  ndata = compute_nwork_ph(no_ph,ndatatot)

  WRITE(stdout,'(5x, "Phdos Free energy, at T= ", f12.6)') temp(itemp)
  WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  WRITE(stdout,'(/,5x,"Fitting the phdos free energy with a quadratic &
                      &function")') 
  WRITE(stdout,'(5x,"Number of degrees of freedom: ",i5)')  degree
  WRITE(stdout,'(5x,"Coefficients of the quadratic:",i5)')  nvar
  WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata
  
  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(f(ndata))
  ALLOCATE(coeff(nvar))

  ndata=0
  SELECT CASE (ibrav)
     CASE(1,2,3) 
       DO idata=1,ndatatot
          IF (.NOT. no_ph(idata)) THEN
             ndata=ndata+1
             x(1,ndata)=celldm_geo(1,idata)
             f(ndata)=ph_free_ener(itemp,idata) 
          ENDIF
       ENDDO 
     CASE(4,5,6,7)
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              IF (ibrav==5) THEN
                 x(2,ndata)=ACOS(celldm_geo(4,idata))
              ELSE
                 x(2,ndata)=celldm_geo(3,idata)
              ENDIF
              f(ndata)=ph_free_ener(itemp,idata) 
           ENDIF
        END DO 
     CASE(8,9,91,10,11)
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              x(2,ndata)=celldm_geo(2,idata)
              x(3,ndata)=celldm_geo(3,idata)
              f(ndata)= ph_free_ener(itemp,idata) 
           END IF
        ENDDO
     CASE(12,-12,13,-13) 
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              x(2,ndata)=celldm_geo(2,idata)
              x(3,ndata)=celldm_geo(3,idata)
              IF (ibrav>0) THEN
!
!   c unique
!
                 x(4,ndata)=ACOS(celldm_geo(4,idata))
              ELSE
!
!   b unique
!
                 x(4,ndata)=ACOS(celldm_geo(5,idata))
              ENDIF
              f(ndata)=ph_free_ener(itemp,idata)
           ENDIF
        ENDDO
     CASE DEFAULT
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              x(2,ndata)=celldm_geo(2,idata)
              x(3,ndata)=celldm_geo(3,idata)
              x(4,ndata)=ACOS(celldm_geo(4,idata))
              x(5,ndata)=ACOS(celldm_geo(5,idata))
              x(6,ndata)=ACOS(celldm_geo(6,idata))
              f(ndata)=ph_free_ener(itemp,idata) 
           END IF
        ENDDO
  END SELECT
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
  WRITE(stdout,'(23x,"x1=",f16.9)') x_pos_min(1)
  IF (degree>1) WRITE(stdout,'(23x,"x2=",f16.9)') x_pos_min(2)
  IF (degree>2) WRITE(stdout,'(23x,"x3=",f16.9)') x_pos_min(3)
  IF (degree>3) WRITE(stdout,'(23x,"x4=",f16.9)') x_pos_min(4)
  IF (degree>4) WRITE(stdout,'(23x,"x5=",f16.9)') x_pos_min(5)
  IF (degree>5) WRITE(stdout,'(23x,"x6=",f16.9)') x_pos_min(6)
  WRITE(stdout,'(5x,"Energy at the extremum",f18.12)') ymin
  WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

  celldm_t(:,itemp)=0.0_DP
  SELECT CASE (ibrav)
     CASE(1,2,3) 
        celldm_t(1,itemp)=x_pos_min(1)
     CASE(4,5,6,7)
        celldm_t(1,itemp)=x_pos_min(1)
        IF (ibrav==5) THEN
           celldm_t(4,itemp)=COS(x_pos_min(2))
        ELSE
           celldm_t(3,itemp)= x_pos_min(2)
        ENDIF
     CASE(8,9,91,10,11)
        celldm_t(1,itemp)=x_pos_min(1)
        celldm_t(2,itemp)=x_pos_min(2)
        celldm_t(3,itemp)=x_pos_min(3)
     CASE(12,-12,13,-13) 
        celldm_t(1,itemp)=x_pos_min(1)
        celldm_t(2,itemp)=x_pos_min(2)
        celldm_t(3,itemp)=x_pos_min(3)
        IF (ibrav>0) THEN
!
!   c unique
!
           celldm_t(4,itemp)=COS(x_pos_min(4))
        ELSE
!
!   b unique
!
           celldm_t(5,itemp)=COS(x_pos_min(4))
        ENDIF
     CASE DEFAULT
        celldm_t(1,itemp)=x_pos_min(1)
        celldm_t(2,itemp)=x_pos_min(2)
        celldm_t(3,itemp)=x_pos_min(3)
        celldm_t(4,itemp)=COS(x_pos_min(4))
        celldm_t(5,itemp)=COS(x_pos_min(5))
        celldm_t(6,itemp)=COS(x_pos_min(6))
  END SELECT
  coeff_t(1:nvar,itemp) = coeff(1:nvar)

  DEALLOCATE(x_pos_min)
  DEALLOCATE(coeff)
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
  USE thermo_mod, ONLY : celldm_geo, energy_geo, omega_geo, no_ph
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
                                 summarize_fitting_data
  IMPLICIT NONE
  INTEGER :: itemp
  INTEGER :: ndata, ndatatot
  REAL(DP), ALLOCATABLE :: x(:,:), y(:), f(:), coeff(:), x_pos_min(:)
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
  !  number of coefficients of the quadratic equation
  !  degrees       nvar
  !  1               3, 
  !  2               6, 
  !  3              10, 
  !  4              15, 
  !  5              21, 
  !  6              28
  !
  ndatatot= compute_nwork()
  ndata = compute_nwork_ph(no_ph,ndatatot)

  WRITE(stdout,'(5x, "Phdos Free energy, at T= ", f12.6)') temp(itemp)
  WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  WRITE(stdout,'(/,5x,"Fitting the phdos free energy with a quadratic &
                      &function")') 
  WRITE(stdout,'(5x,"Number of degrees of freedom: ",i5)')  degree
  WRITE(stdout,'(5x,"Coefficients of the quadratic:",i5)')  nvar
  WRITE(stdout,'(5x,"Number of fitting data:",7x,i5,/)')  ndata
  
  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(f(ndata))
  ALLOCATE(coeff(nvar))

  ndata=0
  SELECT CASE (ibrav)
     CASE(1,2,3) 
       DO idata=1,ndatatot
          IF (.NOT. no_ph(idata)) THEN
             ndata=ndata+1
             x(1,ndata)=celldm_geo(1,idata)
             f(ndata)=phf_free_ener(itemp,idata) 
          ENDIF
       ENDDO 
     CASE(4,5,6,7)
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              IF (ibrav==5) THEN
                 x(2,ndata)=ACOS(celldm_geo(4,idata))
              ELSE
                 x(2,ndata)=celldm_geo(3,idata)
              ENDIF
              f(ndata)=phf_free_ener(itemp,idata) 
           ENDIF
         END DO 
     CASE(8,9,91,10,11)
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              x(2,ndata)=celldm_geo(2,idata)
              x(3,ndata)=celldm_geo(3,idata)
              f(ndata)= phf_free_ener(itemp,idata) 
           ENDIF
        ENDDO
     CASE(12,-12,13,-13) 
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              x(2,ndata)=celldm_geo(2,idata)
              x(3,ndata)=celldm_geo(3,idata)
              IF (ibrav>0) THEN
!
!   c unique
!
                 x(4,ndata)=ACOS(celldm_geo(4,idata))
              ELSE
!
!   b unique
!
                 x(4,ndata)=ACOS(celldm_geo(5,idata))
              ENDIF
              f(ndata)= phf_free_ener(itemp,idata) 
           ENDIF
        ENDDO
     CASE DEFAULT
        DO idata=1,ndatatot
           IF (.NOT. no_ph(idata)) THEN
              ndata=ndata+1
              x(1,ndata)=celldm_geo(1,idata)
              x(2,ndata)=celldm_geo(2,idata)
              x(3,ndata)=celldm_geo(3,idata)
              x(4,ndata)=ACOS(celldm_geo(4,idata))
              x(5,ndata)=ACOS(celldm_geo(5,idata))
              x(6,ndata)=ACOS(celldm_geo(6,idata))
              f(ndata)= phf_free_ener(itemp,idata) 
           ENDIF
        ENDDO
  END SELECT
  !
  !CALL summarize_fitting_data(degree, ndata, x, f)
  !
  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)

  CALL print_quadratic_polynomial(degree, nvar, coeff)

  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  chisq=0.0_DP
  DO idata=1,ndata
     CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
     chisq = chisq + (aux - f(idata))**2
  ENDDO

  CALL find_two_fit_extremum(degree,nvar,x_pos_min,ymin,enthalpy_coeff,coeff)

  WRITE(stdout,'(/,5x,"Extremum found at:")')
  WRITE(stdout,'(23x,"x1=",f16.9)') x_pos_min(1)
  IF (degree>1) WRITE(stdout,'(23x,"x2=",f16.9)') x_pos_min(2)
  IF (degree>2) WRITE(stdout,'(23x,"x3=",f16.9)') x_pos_min(3)
  IF (degree>3) WRITE(stdout,'(23x,"x4=",f16.9)') x_pos_min(4)
  IF (degree>4) WRITE(stdout,'(23x,"x5=",f16.9)') x_pos_min(5)
  IF (degree>5) WRITE(stdout,'(23x,"x6=",f16.9)') x_pos_min(6)
  WRITE(stdout,'(5x,"Energy at the extremum",f18.12)') ymin
  WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

  celldmf_t(:,itemp)=0.0_DP
  SELECT CASE (ibrav)
     CASE(1,2,3) 
        celldmf_t(1,itemp)=x_pos_min(1)
     CASE(4,5,6,7)
        celldmf_t(1,itemp)=x_pos_min(1)
        IF (ibrav==5) THEN
           celldmf_t(4,itemp)=COS(x_pos_min(2))
        ELSE
           celldmf_t(3,itemp)= x_pos_min(2)
        ENDIF
     CASE(8,9,91,10,11)
        celldmf_t(1,itemp)=x_pos_min(1)
        celldmf_t(2,itemp)=x_pos_min(2)
        celldmf_t(3,itemp)=x_pos_min(3)
     CASE(12,-12,13,-13) 
        celldmf_t(1,itemp)=x_pos_min(1)
        celldmf_t(2,itemp)=x_pos_min(2)
        celldmf_t(3,itemp)=x_pos_min(3)
        IF (ibrav>0) THEN
!
!   c unique
!
           celldmf_t(4,itemp)=COS(x_pos_min(4))
        ELSE
!
!   b unique
!
           celldmf_t(5,itemp)=COS(x_pos_min(4))
        ENDIF
     CASE DEFAULT
        celldmf_t(1,itemp)=x_pos_min(1)
        celldmf_t(2,itemp)=x_pos_min(2)
        celldmf_t(3,itemp)=x_pos_min(3)
        celldmf_t(4,itemp)=COS(x_pos_min(4))
        celldmf_t(5,itemp)=COS(x_pos_min(5))
        celldmf_t(6,itemp)=COS(x_pos_min(6))
  END SELECT

  DEALLOCATE(x_pos_min)
  DEALLOCATE(coeff)
  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t_ph

SUBROUTINE compute_degree(ibrav, degree, nvar)
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
