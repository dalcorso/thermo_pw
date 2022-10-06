!
! Copyright (C) 2014-2021 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE write_e_omega()
!-------------------------------------------------------------------------
!
! This routine receives as input the energies at the computed geometries
! and gives as output the energy (or the enthalpy at the given pressure
! if this is non zero) and the pressure as a function of the volume. 
! In the anisotropic case it gives also the crystal parameters as a 
! function of the pressure. 
!
! Note that the extrapolation can be done only in a limited range of
! pressures with respect to the pressure of the computed geometries.
!
! In output we write also the volume omega0 extrapolated at zero pressure.
! It can be used to plot V/V0 as a function of pressure. Note however 
! that the most accurate V0 is calculated when the computed geometries 
! are about it.
!
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE data_files,       ONLY : flevdat
USE cell_base,        ONLY : ibrav
USE thermo_mod,       ONLY : omega_geo, celldm_geo, energy_geo
USE control_mur,      ONLY : omegap0
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_thermo,   ONLY : lgeo_to_file
USE control_pressure, ONLY : pressure, pressure_kb, pmin, pmax, deltap, &
                             npress, press
USE control_quartic_energy, ONLY : lquartic, lsolve
USE geometry_file,      ONLY : write_geometry_output
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum, &
                             print_chisq_quadratic
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             print_quartic_polynomial, print_chisq_quartic
USE polynomial,       ONLY : poly2, poly4, init_poly, clean_poly
USE lattices,         ONLY : compress_celldm, expand_celldm, crystal_parameters
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode, stdout

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1
INTEGER  :: i, iu_mur, ipress, idata, nvar, ndata
INTEGER  :: find_free_unit, compute_nwork
REAL(DP) :: ymin, ymin4
REAL(DP) :: compute_omega_geo
REAL(DP), ALLOCATABLE :: f(:), x(:,:), x_pos_min(:), x_min_4(:), e(:), &
                         omega(:), celldmp(:,:) 
TYPE(poly2) :: p2
TYPE(poly4) :: p4

IF (my_image_id /= root_image) RETURN
!
!  The name of the output files that will contain the volume, energy, pressure
!
filename="energy_files/"//TRIM(flevdat)//'_mur'
CALL add_pressure(filename)
!
!  and file with pressure celldm 
!
filename1="energy_files/"//TRIM(flevdat)//'_mur_celldm'
CALL add_pressure(filename1)
ndata=compute_nwork()

WRITE(stdout,'(/,2x,76("*"))')
WRITE(stdout,'(5x,"Computing the volume as a function of pressure")')
WRITE(stdout,'(5x,"From p_min=",f10.5," kbar to p_max=",f10.5," kbar")') &
       pmin*ry_kbar, pmax*ry_kbar    
WRITE(stdout,'(5x,"Computing",i5," pressures, deltap=",f10.5," kbar")') &
                                                     npress, deltap*ry_kbar
WRITE(stdout,'(/,2x,76("*"))')


nvar=crystal_parameters(ibrav)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
ALLOCATE(e(npress))
ALLOCATE(omega(npress))
ALLOCATE(celldmp(6,npress))
CALL init_poly(nvar,p2)

IF (lquartic) ALLOCATE(x_min_4(nvar))
DO idata=1, ndata
   CALL compress_celldm(celldm_geo(1,idata), x(1,idata), nvar, ibrav)
END DO

DO ipress=1, npress
!
!  fit the enthalpy with a quadratic polynomial and find the minimum
!
   DO idata=1, ndata
     f(idata)=energy_geo(idata) + press(ipress) * omega_geo(idata) / ry_kbar
   END DO
   CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,p2)
!   CALL print_chisq_quadratic(ndata, nvar, x, f, p2)
   CALL find_quadratic_extremum(nvar,x_pos_min,ymin,p2)
   IF (lquartic) THEN
!
!   fit the enthalpy with a quartic polynomial and find the minimum
!
      CALL init_poly(nvar,p4)
      CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,p4)
!      CALL print_quartic_polynomial(nvar,p4)
!      CALL print_chisq_quartic(ndata, nvar, x, f, p4)
      x_min_4=x_pos_min
      CALL find_quartic_extremum(nvar,x_min_4,ymin4,p4)
      CALL expand_celldm(celldmp(1,ipress), x_min_4, nvar, ibrav)
!
!   find the volume that corresponds to the minimum geometry at this pressure
!
      omega(ipress)=compute_omega_geo(ibrav,celldmp(1,ipress))
      e(ipress)=ymin4 - press(ipress) * omega(ipress) / ry_kbar
      CALL clean_poly(p4)
   ELSE
      CALL expand_celldm(celldmp(1,ipress), x_pos_min, nvar, ibrav)
      omega(ipress)=compute_omega_geo(ibrav,celldmp(1,ipress))
      e(ipress) = ymin - press(ipress) * omega(ipress) / ry_kbar
   ENDIF
ENDDO
CALL find_omega0(press/ry_kbar,omega,npress,omegap0)

IF (lgeo_to_file) CALL write_geometry_output(npress, press, celldmp)

IF (vmin_input == 0.0_DP) vmin_input=omega(npress) * 0.98_DP
IF (vmax_input == 0.0_DP) vmax_input=omega(1) * 1.02_DP

IF (ionode) THEN
!
!  Print the energy, enthalpy, and pressure as a function of the volume
!
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (pressure_kb /= 0.0_DP) THEN
      WRITE(iu_mur,'( "#",2x,"omega (a.u.)**3",6x,"enthalpy (Ry)",6x,  &
                            & "enthalpy(p) (Ry)",5x,"pressure (kbar)")')
   ELSE
      WRITE(iu_mur,'( "#",3x,"omega (a.u.)**3",7x,"energy (Ry)",8x,&
                            &"enthalpy(p)(Ry)",4x,"pressure (kbar)")')
   END IF
   DO ipress=1,npress
      WRITE(iu_mur,'(f18.10,3f20.10)') omega(ipress), e(ipress)+ &
            pressure*omega(ipress), e(ipress)+&
                         press(ipress)*omega(ipress)/ry_kbar, press(ipress) 
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
!
!  print the crystal parameter as a function of pressure
!
   OPEN(UNIT=iu_mur, FILE=TRIM(filename1), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_mur,'( "# pressure (kbar)",3x,"celldm(1)",7x,"celldm(2)",7x,&
              &"celldm(3)",7x,"celldm(4)",7x,"celldm(5)",7x,"celldm(6)" )')
   DO ipress=1,npress
      WRITE(iu_mur,'(7f15.8)') press(ipress), celldmp(1, ipress), &
               celldmp(2, ipress), celldmp(3, ipress), celldmp(4, ipress), &
               celldmp(5, ipress), celldmp(6, ipress)
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
END IF

DEALLOCATE(x)
DEALLOCATE(x_pos_min)
DEALLOCATE(f)
DEALLOCATE(e)
DEALLOCATE(omega)
DEALLOCATE(celldmp)
CALL clean_poly(p2)

IF (lquartic) DEALLOCATE(x_min_4)

RETURN
END SUBROUTINE write_e_omega

!-------------------------------------------------------------------------
SUBROUTINE write_e_omega_t(itemp, phf, ndatatot)
!-------------------------------------------------------------------------
!
! This routine receives as input the energies at the computed geometries
! and gives as output the energy (or the enthalpy at the given pressure
! if this is non zero) and the pressure as a function of the volume. 
! In the anisotropic case it gives also the crystal parameters as a 
! function of the pressure. 
!
! Note that the extrapolation can be done only in a limited range of
! pressures with respect to the pressure of the computed geometries.
!
! In output we write also the volume omega0 extrapolated at zero pressure.
! It can be used to plot V/V0 as a function of pressure. Note however 
! that the most accurate V0 is calculated when the computed geometries 
! are about it.
!
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE data_files,       ONLY : flevdat
USE cell_base,        ONLY : ibrav
USE thermo_mod,       ONLY : omega_geo, celldm_geo, energy_geo, no_ph
USE control_vol,      ONLY : vmin_input, vmax_input
USE control_mur,      ONLY : omegap0
USE temperature,      ONLY : temp
USE control_pressure, ONLY : pressure, pressure_kb, pmin, pmax, npress, &
                             press, deltap
USE control_quartic_energy, ONLY : lquartic, lsolve
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum, &
                             print_chisq_quadratic
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             print_quartic_polynomial, print_chisq_quartic
USE polynomial,       ONLY : poly2, poly4, init_poly, clean_poly
USE lattices,         ONLY : compress_celldm, expand_celldm, crystal_parameters
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode, stdout

IMPLICIT NONE

INTEGER  :: itemp
INTEGER  :: ndatatot
REAL(DP) :: phf(ndatatot)

CHARACTER(LEN=256) :: filename, filename1
CHARACTER(LEN=6) :: int_to_char
INTEGER  :: i, iu_mur, ipress, idata, nvar, ndata
INTEGER  :: find_free_unit, compute_nwork
REAL(DP) :: ymin, ymin4
REAL(DP) :: compute_omega_geo
REAL(DP), ALLOCATABLE :: f(:), x(:,:), x_pos_min(:), x_min_4(:), e(:), &
                         omega(:), celldmp(:,:), f1(:), ome(:)
TYPE(poly2) :: p2
TYPE(poly4) :: p4

IF (my_image_id /= root_image) RETURN
IF (itemp<=0) RETURN
!
!  The name of the output files that will contain the volume, energy, pressure
!
filename="anhar_files/"//TRIM(flevdat)//'_mur.'//TRIM(int_to_char(itemp))
CALL add_pressure(filename)

filename1="anhar_files/"//TRIM(flevdat)//'_mur_celldm.'//&
                                                TRIM(int_to_char(itemp))
CALL add_pressure(filename1)

ndata=compute_nwork()

WRITE(stdout,'(/,2x,76("*"))')
WRITE(stdout,'(5x,"Computing the volume as a function of pressure at &
                            &T=",f8.2," K")') temp(itemp)
WRITE(stdout,'(5x,"From p_min=",f10.5," kbar to p_max=",f10.5," kbar")') &
                                              pmin*ry_kbar, pmax*ry_kbar    
WRITE(stdout,'(5x,"Computing",i5," pressures, deltap=",f10.5," kbar")') &
                                                     npress, deltap*ry_kbar
WRITE(stdout,'(/,2x,76("*"))')


nvar=crystal_parameters(ibrav)

ALLOCATE(x(nvar,ndata))
ALLOCATE(x_pos_min(nvar))
ALLOCATE(f(ndata))
ALLOCATE(f1(ndata))
ALLOCATE(e(npress))
ALLOCATE(ome(ndata))
ALLOCATE(omega(npress))
ALLOCATE(celldmp(6,npress))
CALL init_poly(nvar,p2)

IF (lquartic) ALLOCATE(x_min_4(nvar))
ndata=0
DO idata=1,ndatatot
   IF (no_ph(idata)) CYCLE
   ndata=ndata+1
   CALL compress_celldm(celldm_geo(1,idata), x(1,ndata), nvar, ibrav)
   f1(ndata)=energy_geo(idata)+phf(idata)
   ome(ndata)=omega_geo(idata)
END DO

DO ipress=1, npress
!
!  fit the enthalpy with a quadratic polynomial and find the minimum
!
   DO idata=1, ndata
     f(idata)=f1(idata) + press(ipress) * ome(idata)/ry_kbar
   END DO
   CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,p2)
!   CALL print_chisq_quadratic(ndata, nvar, x, f, p2)
   CALL find_quadratic_extremum(nvar,x_pos_min,ymin,p2)
   IF (lquartic) THEN
!
!   fit the enthalpy with a quartic polynomial and find the minimum
!
      CALL init_poly(nvar,p4)
      CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,p4)
!      CALL print_quartic_polynomial(nvar,p4)
!      CALL print_chisq_quartic(ndata, nvar, x, f, p4)
      x_min_4=x_pos_min
      CALL find_quartic_extremum(nvar,x_min_4,ymin4,p4)
      CALL expand_celldm(celldmp(1,ipress), x_min_4, nvar, ibrav)
!
!   find the volume that corresponds to the minimum geometry at this pressure
!
      omega(ipress)=compute_omega_geo(ibrav,celldmp(1,ipress))
      e(ipress)=ymin4 - press(ipress) * omega(ipress)/ry_kbar
      CALL clean_poly(p4)
   ELSE
      CALL expand_celldm(celldmp(1,ipress), x_pos_min, nvar, ibrav)
      omega(ipress)=compute_omega_geo(ibrav,celldmp(1,ipress))
      e(ipress) = ymin - press(ipress) * omega(ipress)/ry_kbar
   ENDIF
ENDDO

IF (vmin_input == 0.0_DP) vmin_input=omega(npress) * 0.98_DP
IF (vmax_input == 0.0_DP) vmax_input=omega(1) * 1.02_DP

IF (ionode) THEN
!
!  Print the energy, enthalpy, and pressure as a function of the volume
!
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   CALL write_mur_start_line(itemp, iu_mur)
   DO ipress=1,npress
      WRITE(iu_mur,'(f18.10,3f20.10)') omega(ipress), e(ipress)+ &
            pressure*omega(ipress), e(ipress)+press(ipress)*omega(ipress), &
            press(ipress) * ry_kbar
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
!
!  print the crystal parameter as a function of pressure
!
   OPEN(UNIT=iu_mur, FILE=TRIM(filename1), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_mur,'( "# pressure (kbar)",3x,"celldm(1)",7x,"celldm(2)",7x,&
              &"celldm(3)",7x,"celldm(4)",7x,"celldm(5)",7x,"celldm(6) T=",&
                       &f8.2," K" )') temp(itemp)
   DO ipress=1,npress
      WRITE(iu_mur,'(7f15.8)') press(ipress), celldmp(1, ipress), &
               celldmp(2, ipress), celldmp(3, ipress), celldmp(4, ipress), &
               celldmp(5, ipress), celldmp(6, ipress)
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
END IF

DEALLOCATE(x)
DEALLOCATE(x_pos_min)
DEALLOCATE(f)
DEALLOCATE(f1)
DEALLOCATE(e)
DEALLOCATE(omega)
DEALLOCATE(ome)
DEALLOCATE(celldmp)
CALL clean_poly(p2)

IF (lquartic) DEALLOCATE(x_min_4)

RETURN
END SUBROUTINE write_e_omega_t

!-----------------------------------------------------------------------
SUBROUTINE find_omega0(press,omega,npress,omega0)
!-----------------------------------------------------------------------
!
!  Given the volume at many pressures, find the volume at zero pressure,
!  linearly interpolating among the two pressures that contain the zero.
!
USE kinds,            ONLY : DP
IMPLICIT NONE
INTEGER  :: npress
REAL(DP) :: press(npress), omega(npress)
REAL(DP) :: omega0

INTEGER :: ipress, ipress0

ipress0=0
DO ipress=2, npress
   IF (press(ipress) <= 0.0_DP) ipress0=ipress
ENDDO
IF (ipress0==0.OR.ipress0==npress) THEN
!
!  In this case the zero pressure is out of range. Set omega0 zero so
!  the plot is not done.
!
   omega0=0.0_DP
   RETURN
ENDIF

omega0=omega(ipress0) - (omega(ipress0+1) - omega(ipress0)) * press(ipress0) &
                      / (press(ipress0+1) - press(ipress0))

RETURN
END SUBROUTINE find_omega0
