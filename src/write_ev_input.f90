!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_ev_input(file_dat)
  !-----------------------------------------------------------------------
  !
  !  This routine receives the summary of the cell volumes and of
  !  the energy at each volume and writes the input file 
  !  for the ev_sub subroutine.
  !
  !
  USE kinds, ONLY : DP
  USE mp_images, ONLY : my_image_id, root_image
  USE thermo_mod, ONLY : omega_geo, energy_geo, ngeo
  USE control_pressure, ONLY : pressure
  USE io_global, ONLY : ionode
  IMPLICIT NONE
  CHARACTER(LEN=256) :: file_dat
  INTEGER :: iu_ev, igeom
  !
  IF (my_image_id /= root_image) RETURN
  !
  CALL write_ev_driver(file_dat) 
  !
  IF (ionode) THEN
     iu_ev=2
     OPEN(UNIT=iu_ev, FILE=TRIM(file_dat), STATUS='UNKNOWN', FORM='FORMATTED')
     DO igeom=1,ngeo(1)
        WRITE(iu_ev,'(2e30.15)') omega_geo(igeom), energy_geo(igeom) + &
                                        pressure * omega_geo(igeom)
     ENDDO
     !
     CLOSE(iu_ev)
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE write_ev_input

SUBROUTINE write_ev_driver(file_dat)

USE kinds, ONLY : DP
USE io_global, ONLY : ionode
USE control_pressure, ONLY : pressure, pressure_kb
IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: file_dat
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char
INTEGER :: iu_ev
!
IF (ionode) THEN
   iu_ev=2
   filename='input_ev'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
     
   OPEN(UNIT=iu_ev, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_ev,'("au")')
   WRITE(iu_ev,'("hex")')
   WRITE(iu_ev,'("4")')
   WRITE(iu_ev,'(a)') TRIM(file_dat)
   WRITE(iu_ev,'(a)') TRIM(file_dat)//'.ev.out'
   CLOSE(iu_ev)
ENDIF
!
RETURN
END SUBROUTINE write_ev_driver

SUBROUTINE do_ev()
!
!  This subroutine computes the equilibrium volume and bulk modulus
!
USE kinds, ONLY : DP
USE control_mur, ONLY : vmin, b0, b01, emin
USE control_pressure, ONLY : pressure, pressure_kb
USE data_files, ONLY : flevdat
USE io_global, ONLY : meta_ionode_id, stdout
USE mp_world, ONLY : world_comm
USE mp, ONLY : mp_bcast

IMPLICIT NONE
CHARACTER(LEN=256) :: file_dat, filename
CHARACTER(LEN=8) :: float_to_char

  file_dat=TRIM(flevdat) 
  filename='input_ev'
  IF (pressure /= 0.0_DP) THEN
      file_dat=TRIM(file_dat)//'.'//TRIM(float_to_char(pressure_kb,1))
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
  ENDIF

  CALL write_ev_input(file_dat)
  CALL ev_sub(vmin, b0, b01, emin, filename )

  CALL mp_bcast(vmin, meta_ionode_id, world_comm)
  CALL mp_bcast(b0, meta_ionode_id, world_comm)
  CALL mp_bcast(b01, meta_ionode_id, world_comm)
  CALL mp_bcast(emin, meta_ionode_id, world_comm)
  !
  RETURN
END SUBROUTINE do_ev

SUBROUTINE do_ev_t(itemp)
!
!  This subroutine compute the equilibrium volume and bulk modulus
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo, no_ph
USE control_mur,    ONLY : vmin, b0, b01
USE thermodynamics, ONLY : ph_free_ener
USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, free_e_min_t
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure, pressure_kb
USE data_files,     ONLY : flevdat
USE quadratic_surfaces, ONLY : polifit
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER, PARAMETER :: m1=5
INTEGER :: idata, ndata, i1
REAL(DP) :: a(m1), x(ngeo(1)), y(ngeo(1)), aux, aux1


  IF (my_image_id /= root_image) RETURN

!  WRITE(stdout,*) 
  ndata=0
  DO idata=1,ngeo(1)
     IF (.NOT. no_ph(idata)) THEN
        ndata=ndata+1
        x(ndata)=omega_geo(idata)
        y(ndata)=ph_free_ener(itemp,idata)
!       WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
     ENDIF
  ENDDO
  CALL polifit(x, y, ndata, a, m1)

  CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, a, m1, vm)
  aux = (vmin / vm)**b01 * b0
  DO i1=3,m1
     aux=aux+ (i1-2.0_DP) * (i1 - 1.0_DP) * a(i1) * vm ** (i1-2) * ry_kbar
  ENDDO
  aux1= b0 * b01 / aux * ( vmin / vm )** b01
  DO i1=3,m1
     aux1=aux1- (i1-2.0_DP)**2*(i1 - 1.0_DP) * a(i1) * vm ** (i1-2) * ry_kbar &
                                 / aux
  ENDDO

  vmin_t(itemp)=vm
  b0_t(itemp)=aux
  b01_t(itemp)=aux1

!  free_e_min_t(itemp)=free_e
  !
  CALL compute_celldm_geo(vmin_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("-"))')
  WRITE(stdout,'(5x, "free energy from phonon dos, at T= ", f12.6)') temp(itemp)
  IF (pressure /= 0.0_DP) &
     WRITE(stdout, '(5x,"pressure = ",f15.6," kbar")') pressure_kb
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
                                 &" a.u.")') celldm_(1)
  WRITE(stdout,'(5x, "The bulk modulus is ",24x,f12.3,"  kbar")')  b0_t(itemp)
  WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",&
                               &f9.3)')  b01_t(itemp)
  WRITE(stdout,'(2x,76("-"),/)')

  RETURN
END SUBROUTINE do_ev_t

SUBROUTINE do_ev_t_ph(itemp) 
!
!  This subroutine compute the equilibrium volume and bulk modulus
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo, no_ph
USE constants,      ONLY : ry_kbar
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE ph_freq_anharmonic,     ONLY : vminf_t, b0f_t, b01f_t, free_e_minf_t
USE control_mur,    ONLY : vmin, b0, b01
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure, pressure_kb
USE quadratic_surfaces, ONLY : polifit
USE data_files,     ONLY : flevdat
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER, PARAMETER :: m1=5
INTEGER :: idata, ndata, i1
REAL(DP) :: a(m1), x(ngeo(1)), y(ngeo(1)), aux, aux1


  IF (my_image_id /= root_image) RETURN

!  WRITE(stdout,*) 
  ndata=0
  DO idata=1,ngeo(1)
     IF (.NOT. no_ph(idata)) THEN
        ndata=ndata+1
        x(ndata)=omega_geo(idata)
        y(ndata)=phf_free_ener(itemp,idata) 
!       WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
     ENDIF
  ENDDO
  CALL polifit(x, y, ndata, a, m1)

  CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, a, m1, vm)
  aux = (vmin / vm)**b01 * b0
  DO i1=3,m1
     aux=aux+ (i1-2.0_DP) * (i1 - 1.0_DP) * a(i1) * vm ** (i1-2) * ry_kbar
  ENDDO
  aux1= b0 * b01 / aux * ( vmin / vm )** b01
  DO i1=3,m1
     aux1=aux1- (i1-2.0_DP)**2*(i1 - 1.0_DP) * a(i1) * vm ** (i1-2) * ry_kbar &
                                 / aux
  ENDDO

  vminf_t(itemp)=vm
  b0f_t(itemp)=aux
  b01f_t(itemp)=aux1

!  free_e_minf_t(itemp)=free_e
  !
  CALL compute_celldm_geo(vminf_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x, "free energy from phonon frequencies at T=",f12.4)') &
                                                               temp(itemp)
  IF (pressure /= 0.0_DP) &
     WRITE(stdout, '(5x,"pressure = ",f15.6," kbar")') pressure_kb
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",16x,f12.4,&
                                 &" a.u.")') celldm_(1)
  WRITE(stdout,'(5x, "The bulk modulus is ",31x,f12.3,"  kbar")')  b0f_t(itemp)
  WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",5x,&
                               &f11.3)')  b01f_t(itemp)

  WRITE(stdout,'(2x,76("+"),/)')

  RETURN
END SUBROUTINE do_ev_t_ph

SUBROUTINE find_min_mur_pol(v0, b0, b01, a, m1, vm)
!
!  This routine minimizes a function equal to the Murnaghan equation
!  with parameters v0, b0, and b01 and a polynomial of degree m1-1 of
!  the form a(1) + a(2) * v + a(3) * v**2 +... a(m1) * v**(m1-1)
!  NB: b0 must be in atomic units not kbar.
!
USE kinds, ONLY : DP 
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, a(m1)
REAL(DP), INTENT(OUT) :: vm
REAL(DP), PARAMETER :: tol=1.D-12
INTEGER, PARAMETER :: maxstep=100
REAL(DP) :: fx, fxp, v1, v1old
REAL(DP) :: compute_fun, compute_fun_deriv
INTEGER :: istep

v1=v0*0.8_DP
v1old=v1
DO istep=1,maxstep
   fx = compute_fun(v1, v0, b0, b01, a, m1)
   fxp = compute_fun_deriv(v1, v0, b0, b01, a, m1)
!
!  Newton method
!
   v1 = v1 - fx/fxp
!   WRITE(stdout,'(5x,"Step", i4, " V1=", f20.12, " f= ", f20.12)') istep, v1, fx
   IF (ABS(v1-v1old) < tol .OR. ABS(fx) < tol ) GOTO 100
   v1old=v1
ENDDO
CALL errore('find_min_mur_pol','minimum not found',1)
100 CONTINUE
vm=v1
!WRITE(stdout,'("Vmin", 3f20.12)') vm

RETURN
END SUBROUTINE find_min_mur_pol

FUNCTION compute_fun(v, v0, b0, b01, a, m1)
USE kinds, ONLY : DP 
IMPLICIT NONE
REAL(DP) :: compute_fun
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v, v0, b0, b01, a(m1)
REAL(DP) :: aux
INTEGER :: i1

aux=b0 * ( 1.0_DP  - (v0 / v ) **b01 ) / b01 + a(2)

DO i1=3, m1
   aux = aux + (i1-1.0_DP) * a(i1) * v**(i1-2)
ENDDO

compute_fun=aux
RETURN
END FUNCTION compute_fun

FUNCTION compute_fun_deriv(v, v0, b0, b01, a, m1)
USE kinds, ONLY : DP 
IMPLICIT NONE
REAL(DP) :: compute_fun_deriv
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v, v0, b0, b01, a(m1)
REAL(DP) :: aux
INTEGER :: i1

aux=b0 * (v0 / v )**( b01+1.0_DP)/ v0 + a(3) * 2.0_DP
DO i1=4, m1
   aux = aux + (i1-1.0_DP) * (i1-2.0_DP) * a(i1) * v**(i1-3)
ENDDO
compute_fun_deriv=aux

RETURN
END FUNCTION compute_fun_deriv
