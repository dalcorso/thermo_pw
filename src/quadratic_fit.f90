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
  !   It finds also the minimum of the quadratic function.
  !   
  !   the output of this routine is celldm0(:)
  !
  !
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : ibrav
  USE control_mur,  ONLY : emin
  USE equilibrium_conf, ONLY : celldm0
  USE mp_images,    ONLY : my_image_id, root_image
  USE thermo_mod,   ONLY : celldm_geo, energy_geo, omega_geo
  USE control_pressure, ONLY : pressure, pressure_kb
  USE control_quadratic_energy, ONLY : hessian_v, hessian_e, x_pos_min, &
                               coeff, degree
  USE control_quartic_energy, ONLY : nvar4, coeff4, x_min_4, lquartic, lsolve
  USE io_global,    ONLY : stdout
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_fit_extremum, &
                                 write_fit_hessian, evaluate_fit_quadratic, &
                                 print_quadratic_polynomial, &
                                 summarize_fitting_data, write_vector, &
                                 introduce_quadratic_fit
  USE quartic_surfaces, ONLY : fit_multi_quartic, compute_quartic_var, &
                               find_quartic_extremum, evaluate_fit_quartic, &
                               print_quartic_polynomial, introduce_quartic_fit
  IMPLICIT NONE

  INTEGER  :: nvar, ndata
  REAL(DP), ALLOCATABLE :: x(:,:), y(:), f(:)
  REAL(DP) :: ymin, chisq, aux, ymin4
  INTEGER  :: idata
  INTEGER  :: compute_nwork
  !
  ! Only the first image does the calculation
  !
  celldm0(:)=0.0_DP
  IF (my_image_id /= root_image) RETURN
  !
  WRITE(stdout,'(/,5x,71("-"))')
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
  WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

  CALL find_fit_extremum(degree,nvar,x_pos_min,ymin,coeff)
  !
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')

  CALL write_vector(degree,x_pos_min)
  !
  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x,"Enthalpy at the extremum",f18.12)') ymin
  ELSE
     WRITE(stdout,'(5x,"Energy at the extremum",f18.12)') ymin
  END IF
  !
  CALL write_fit_hessian(degree,nvar,coeff,hessian_v,hessian_e)
  !
  CALL set_celldm_from_xmin(ibrav, degree, x_pos_min, celldm0)
  !
  emin=ymin
  !
  IF (lquartic) THEN

     nvar4=compute_quartic_var(degree)
     CALL introduce_quartic_fit(degree, nvar4, ndata)

     ALLOCATE(x_min_4(degree))
     ALLOCATE(coeff4(nvar4))

     CALL fit_multi_quartic(ndata,degree,nvar4,lsolve,x,f,coeff4)
     !
     CALL print_quartic_polynomial(degree, nvar4, coeff4)
!    WRITE(stdout,'(/,7x,"Energy (1)    Fitted energy (2)   DeltaE (1)-(2)")') 
     chisq=0.0_DP
     DO idata=1,ndata
        CALL evaluate_fit_quartic(degree,nvar4,x(1,idata),aux,coeff4)
!       WRITE(stdout,'(2f10.4,2f19.12,e19.12)') x(1,idata), x(2,idata), f(idata), &
!                                                   aux, f(idata)-aux
        chisq = chisq + (aux - f(idata))**2
     ENDDO
     WRITE(stdout,'(/,5x,"chi square=",e18.5,/)') chisq
!
!   searching the minimum starting from the minimum of the quadratic
!
     x_min_4=x_pos_min
     CALL find_quartic_extremum(degree,nvar4,x_min_4,ymin4,coeff4)

     WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     CALL write_vector(degree,x_min_4)
     IF (pressure > 0.0_DP) THEN
        WRITE(stdout,'(5x,"Enthalpy at the extremum",f18.12)') ymin4
     ELSE
        WRITE(stdout,'(5x,"Energy at the extremum",f18.12)') ymin4
     END IF
     CALL set_celldm_from_xmin(ibrav, degree, x_min_4, celldm0)
     emin=ymin4
  ENDIF

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
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE mp_images,   ONLY : my_image_id, root_image
  USE thermo_mod,  ONLY : celldm_geo, energy_geo, no_ph, omega_geo
  USE control_quadratic_energy, ONLY : degree, nvar, coeff_t, &
                                       enthalpy_coeff => coeff 
  USE control_quartic_energy, ONLY :  nvar4, coeff4, lquartic, lquartic_ph, &
                                      lsolve
  USE temperature, ONLY : temp
  USE control_pressure, ONLY : pressure, pressure_kb
  USE thermodynamics, ONLY : ph_free_ener
  USE anharmonic,  ONLY : celldm_t
  USE io_global,   ONLY : stdout
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_two_fit_extremum, &
                          write_fit_hessian, evaluate_fit_quadratic,  &
                          print_quadratic_polynomial, &
                          summarize_fitting_data, write_vector, &
                          introduce_quadratic_fit
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum, &
                          evaluate_quartic_quadratic, fit_multi_quartic, &
                          find_two_quartic_extremum, evaluate_two_quartic
  IMPLICIT NONE
  INTEGER :: itemp
  INTEGER :: ndata, ndatatot
  REAL(DP), ALLOCATABLE :: x(:,:), f(:), coeff(:), x_pos_min(:), &
            celldm_data(:,:), fun(:), coefft4(:)
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

  IF (MOD(itemp-1,50)==0) &
     CALL introduce_quadratic_fit(degree, nvar, ndata)

  ALLOCATE(x(degree,ndata))
  ALLOCATE(x_pos_min(degree))
  ALLOCATE(f(ndata))
  ALLOCATE(fun(ndata))
  ALLOCATE(coeff(nvar))
  ALLOCATE(coefft4(nvar4))
  ALLOCATE(celldm_data(6, ndata))

  ndata=0
  DO idata=1,ndatatot
     IF (.NOT.no_ph(idata)) THEN
        ndata=ndata+1
        celldm_data(:,ndata)=celldm_geo(:,idata)
        f(ndata)=ph_free_ener(itemp,idata)
        fun(ndata)=f(ndata)+energy_geo(idata)+ pressure * omega_geo(idata)
     END IF
  END DO

  CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_data)
  !
  !CALL summarize_fitting_data(degree, ndata, x, f)
  !
  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)


!  CALL print_quadratic_polynomial(degree, nvar, coeff)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  chisq=0.0_DP
  DO idata=1,ndata
     CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
     chisq = chisq + (aux - f(idata))**2
  ENDDO
  WRITE(stdout,'(5x,"chi square=",e18.5)') chisq

  CALL find_two_fit_extremum(degree,nvar,x_pos_min,ymin,enthalpy_coeff,coeff)
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
  CALL write_vector(degree,x_pos_min)
  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x,"Gibbs energy at the extremum",f22.12)') ymin
  ELSE
     WRITE(stdout,'(5x,"Free energy at the extremum",f22.12)') ymin
  END IF

  IF (lquartic) THEN
     IF (lquartic_ph) THEN
        WRITE(stdout,'(/,5x, "Fit improved with a fourth order polynomial")') 
        CALL fit_multi_quartic(ndata,degree,nvar4,lsolve,x,f,coefft4)
     ELSE
        WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
     ENDIF
     chisq=0.0_DP
     DO idata=1,ndata
        IF (lquartic_ph) THEN
           CALL evaluate_two_quartic(degree,nvar4,x(1,idata),aux,&
                                                         coeff4,coefft4)
        ELSE
           CALL evaluate_quartic_quadratic(degree,nvar4,nvar,x(1,idata),aux,&
                                                         coeff4,coeff)
        END IF
!       WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
        chisq = chisq + (aux - fun(idata))**2
     ENDDO
     WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

     IF (lquartic_ph) THEN
        CALL find_two_quartic_extremum(degree,nvar4,x_pos_min,&
                                                        ymin,coeff4,coefft4)
     ELSE
        CALL find_quartic_quadratic_extremum(degree,nvar4,nvar,x_pos_min,&
                                                            ymin,coeff4,coeff)
     ENDIF
     WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     CALL write_vector(degree,x_pos_min)
     IF (pressure > 0.0_DP) THEN
        WRITE(stdout,'(5x,"Gibbs energy at the extremum",f22.12)') ymin
     ELSE
        WRITE(stdout,'(5x,"Free energy at the extremum",f22.12)') ymin
     END IF
  END IF

  CALL set_celldm_from_xmin(ibrav, degree, x_pos_min, celldm_t(1,itemp))

  coeff_t(1:nvar,itemp) = coeff(1:nvar)

  DEALLOCATE(x_pos_min)
  DEALLOCATE(coeff)
  DEALLOCATE(coefft4)
  DEALLOCATE(celldm_data)
  DEALLOCATE(fun)
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
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE mp_images,   ONLY : my_image_id, root_image
  USE thermo_mod,  ONLY : celldm_geo, energy_geo, no_ph, omega_geo
  USE control_quadratic_energy, ONLY : degree, nvar, coeff_t, &
                         enthalpy_coeff => coeff
  USE temperature, ONLY : temp
  USE control_quartic_energy, ONLY :  nvar4, coeff4, lquartic, lquartic_ph, &
                                      lsolve
  USE control_pressure, ONLY : pressure, pressure_kb
  USE ph_freq_thermodynamics, ONLY : phf_free_ener
  USE ph_freq_anharmonic, ONLY : celldmf_t
  USE io_global,   ONLY : stdout
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_two_fit_extremum, &
                          write_fit_hessian, evaluate_fit_quadratic,  &
                          print_quadratic_polynomial, &
                          summarize_fitting_data, write_vector, &
                          introduce_quadratic_fit
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum, &
                          evaluate_quartic_quadratic, fit_multi_quartic, &
                          find_two_quartic_extremum, evaluate_two_quartic

  IMPLICIT NONE
  INTEGER :: itemp
  INTEGER :: ndata, ndatatot
  REAL(DP), ALLOCATABLE :: x(:,:), f(:), coeff(:), x_pos_min(:), &
                           celldm_data(:,:), fun(:), coefft4(:)
  REAL(DP) :: ymin, chisq, aux
  INTEGER  :: idata
  INTEGER  :: compute_nwork, compute_nwork_ph
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
  ALLOCATE(fun(ndata))
  ALLOCATE(coeff(nvar))
  ALLOCATE(coefft4(nvar4))
  ALLOCATE(celldm_data(6,ndata))

  ndata=0
  DO idata=1,ndatatot
     IF (.NOT. no_ph(idata)) THEN
        ndata=ndata+1
        celldm_data(:,ndata)=celldm_geo(:,idata)
        f(ndata)=phf_free_ener(itemp,idata)
        fun(ndata)=f(ndata)+energy_geo(idata)+ pressure * omega_geo(idata)
     END IF
  END DO

  CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_data)
  !
  !CALL summarize_fitting_data(degree, ndata, x, f)
  !
  CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)

!  CALL print_quadratic_polynomial(degree, nvar, coeff)

!  WRITE(stdout,'(/,7x,"Energy (1)      Fitted energy (2)   DeltaE (1)-(2)")') 
  chisq=0.0_DP
  DO idata=1,ndata
     CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,coeff)
!     WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
     chisq = chisq + (aux - f(idata))**2
  ENDDO
  WRITE(stdout,'(5x,"chi square=",e18.5)') chisq

  CALL find_two_fit_extremum(degree,nvar,x_pos_min,ymin,enthalpy_coeff,coeff)
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
  CALL write_vector(degree,x_pos_min)
  IF (pressure > 0.0_DP) THEN
     WRITE(stdout,'(5x,"Gibbs energy at the extremum",f22.12)') ymin
  ELSE
     WRITE(stdout,'(5x,"Free energy at the extremum",f22.12)') ymin
  END IF

  IF (lquartic) THEN
     IF (lquartic_ph) THEN
        WRITE(stdout,'(/,5x, "Fit improved with a fourth order polynomial")') 
        CALL fit_multi_quartic(ndata,degree,nvar4,lsolve,x,f,coefft4)
     ELSE
        WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
     ENDIF
     chisq=0.0_DP
     DO idata=1,ndata
        IF (lquartic_ph) THEN
           CALL evaluate_two_quartic(degree,nvar4,x(1,idata),aux,coeff4,coefft4)
        ELSE
           CALL evaluate_quartic_quadratic(degree,nvar4,nvar,x(1,idata),aux,&
                                                     coeff4,coeff)
        ENDIF
!       WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
        chisq = chisq + (aux - fun(idata))**2
     ENDDO
     WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq


     IF (lquartic_ph) THEN
        CALL find_two_quartic_extremum(degree,nvar4,x_pos_min,&
                                                        ymin,coeff4,coefft4)
     ELSE
        CALL find_quartic_quadratic_extremum(degree,nvar4,nvar,x_pos_min,&
                                                            ymin,coeff4,coeff)
     ENDIF
     WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     CALL write_vector(degree,x_pos_min)
     IF (pressure > 0.0_DP) THEN
        WRITE(stdout,'(5x,"Gibbs energy at the extremum",f22.12)') ymin
     ELSE
        WRITE(stdout,'(5x,"Free energy at the extremum",f22.12)') ymin
     END IF
  END IF

  CALL set_celldm_from_xmin(ibrav, degree, x_pos_min, celldmf_t(1,itemp))

  DEALLOCATE(x_pos_min)
  DEALLOCATE(coefft4)
  DEALLOCATE(coeff)
  DEALLOCATE(celldm_data)
  DEALLOCATE(fun)
  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t_ph

SUBROUTINE compute_degree(ibrav, degree, nvar)
!
!  For each Bravais lattice this routine gives the number of independent
!  crystallographic parameters and the number of variables of a 
!  quadratic polynomial of these variables
!
!  degrees       nvar
!  1               3, 
!  2               6, 
!  3              10, 
!  4              15, 
!  5              21, 
!  6              28
!
USE lattices, ONLY : crystal_parameters
IMPLICIT NONE
INTEGER, INTENT(IN)  :: ibrav
INTEGER, INTENT(OUT) :: degree, nvar

degree = crystal_parameters(ibrav)
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
!
!  this rouotine receives and array of values of celldm celldm_geo(6,ndata)
!  and transform it in a compact array x(degree,ndata), where degree depends
!  on the Bravais lattice
!
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, degree, ndata
REAL(DP), INTENT(IN) :: celldm_geo(6,ndata)
REAL(DP), INTENT(INOUT) :: x(degree,ndata)

INTEGER :: idata

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
      ENDDO 
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
!
!  This routine receives a set of crystallographic parameters in the
!  array x(degree) and transform it in the celldm array
!
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER,  INTENT(IN)    :: ibrav, degree
REAL(DP), INTENT(IN)    :: x(degree)
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
