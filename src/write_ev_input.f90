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
  USE mp_images,        ONLY : my_image_id, root_image
  USE thermo_mod,       ONLY : omega_geo, energy_geo, ngeo
  USE control_pressure, ONLY : pressure
  USE io_global,        ONLY : ionode

  IMPLICIT NONE
  CHARACTER(LEN=256) :: file_dat
  CHARACTER(LEN=256) :: filedata
  INTEGER :: iu_ev, igeom
  INTEGER :: find_free_unit
  !
  IF (my_image_id /= root_image) RETURN
  !
  filedata=TRIM(file_dat)
  CALL write_ev_driver(filedata) 
  !
  IF (ionode) THEN
     iu_ev=find_free_unit()
     OPEN(UNIT=iu_ev, FILE=TRIM(filedata), STATUS='UNKNOWN', FORM='FORMATTED')
     DO igeom=1,ngeo(1)
        WRITE(iu_ev,'(2e30.15)') omega_geo(igeom), energy_geo(igeom) + &
                                        pressure * omega_geo(igeom)
     ENDDO
     !
     CLOSE(UNIT=iu_ev,STATUS='keep')
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE write_ev_input

!-----------------------------------------------------------------------
SUBROUTINE write_ev_driver(file_dat)
!-----------------------------------------------------------------------

USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: file_dat
CHARACTER(LEN=256) :: filename
INTEGER :: iu_ev
INTEGER :: find_free_unit
!
IF (ionode) THEN
   iu_ev=find_free_unit()
   filename='energy_files/input_ev'
   CALL add_pressure(filename)
     
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

!-----------------------------------------------------------------------
SUBROUTINE do_ev()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!
USE control_mur, ONLY : vmin, b0, b01, emin
USE data_files,  ONLY : flevdat
USE io_global,   ONLY : meta_ionode_id, stdout
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_bcast

IMPLICIT NONE
CHARACTER(LEN=256) :: file_dat, filename

  file_dat="energy_files/"//TRIM(flevdat) 
  CALL add_pressure(file_dat)
  filename='energy_files/input_ev'
  CALL add_pressure(filename)

  CALL write_ev_input(file_dat)
  CALL ev_sub(vmin, b0, b01, emin, filename )

  CALL mp_bcast(vmin, meta_ionode_id, world_comm)
  CALL mp_bcast(b0, meta_ionode_id, world_comm)
  CALL mp_bcast(b01, meta_ionode_id, world_comm)
  CALL mp_bcast(emin, meta_ionode_id, world_comm)
  !
  RETURN
END SUBROUTINE do_ev

!-----------------------------------------------------------------------
SUBROUTINE do_ev_t(itemp)
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!  at a given temperature. It uses the free energy computed by phonon dos
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo, no_ph
USE control_mur,    ONLY : vmin, b0, b01, emin
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_eldos,  ONLY : lel_free_energy
USE thermodynamics, ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, free_e_min_t, a_t
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE data_files,     ONLY : flevdat
USE polyfit_mod,    ONLY : polyfit
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER  :: m1
INTEGER :: idata, ndata, i1
REAL(DP) :: a(poly_degree_ph+1), x(ngeo(1)), y(ngeo(1)), aux, aux1
REAL(DP) :: compute_mur_fun


  IF (my_image_id /= root_image) RETURN

  m1=poly_degree_ph+1
!  WRITE(stdout,*) 
  ndata=0
  DO idata=1,ngeo(1)
     IF (no_ph(idata)) CYCLE
     ndata=ndata+1
     x(ndata)=omega_geo(idata)
     y(ndata)=ph_free_ener(itemp,idata)
     IF (lel_free_energy) y(ndata)=y(ndata)+el_free_ener(itemp,idata)
!    WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
  ENDDO
  CALL polyfit(x, y, ndata, a, poly_degree_ph)
  CALL save_free_energy_on_file(itemp) 

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
  free_e_min_t(itemp)=emin + compute_mur_fun(vm, vmin, b0/ry_kbar, b01, a, m1)
  a_t(1:m1,itemp)=a(1:m1)
  !
  CALL compute_celldm_geo(vmin_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("-"))')
  WRITE(stdout,'(5x, "free energy from phonon dos, at T= ", f12.6)') &
                                                                temp(itemp)
  CALL summarize_mur(celldm_(1), b0_t(itemp), b01_t(itemp), &
                                                     free_e_min_t(itemp))

  WRITE(stdout,'(2x,76("-"),/)')

  RETURN
END SUBROUTINE do_ev_t

!-----------------------------------------------------------------------
SUBROUTINE do_ev_t_ph(itemp) 
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!  at a given temperature. It uses the free energy computed by the
!  sum over the phonon frequencies
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo, no_ph
USE constants,      ONLY : ry_kbar
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics,      ONLY : el_free_ener
USE ph_freq_anharmonic,     ONLY : vminf_t, b0f_t, b01f_t, free_e_minf_t
USE control_mur,    ONLY : emin, vmin, b0, b01
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_eldos,  ONLY : lel_free_energy
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE polyfit_mod,    ONLY : polyfit
USE data_files,     ONLY : flevdat
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER :: m1
INTEGER :: idata, ndata, i1
REAL(DP) :: a(poly_degree_ph+1), x(ngeo(1)), y(ngeo(1)), aux, aux1
REAL(DP) :: compute_mur_fun


  IF (my_image_id /= root_image) RETURN

  m1=poly_degree_ph+1
!  WRITE(stdout,*) 
  ndata=0
  DO idata=1,ngeo(1)
     IF (no_ph(idata)) CYCLE
     ndata=ndata+1
     x(ndata)=omega_geo(idata)
     y(ndata)=phf_free_ener(itemp,idata) 
     IF (lel_free_energy) y(ndata)=y(ndata)+el_free_ener(itemp,idata)
!    WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
  ENDDO
  CALL polyfit(x, y, ndata, a, poly_degree_ph)

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
  free_e_minf_t(itemp)=emin + compute_mur_fun(vm, vmin, b0/ry_kbar, b01, a, m1)
  !
  CALL compute_celldm_geo(vminf_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x, "free energy from phonon frequencies at T=",f12.4)') &
                                                               temp(itemp)
  CALL summarize_mur(celldm_(1), b0f_t(itemp), b01f_t(itemp), &
                                                      free_e_minf_t(itemp))

  WRITE(stdout,'(2x,76("+"),/)')

  RETURN
END SUBROUTINE do_ev_t_ph

!-----------------------------------------------------------------------
SUBROUTINE do_ev_t_el(itemp)
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!  at a given temperature. It uses the free energy computed by phonon dos
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo
USE control_mur,    ONLY : vmin, b0, b01, emin
USE control_quartic_energy, ONLY : poly_degree_ph
USE el_thermodynamics, ONLY : el_free_ener
USE el_anharmonic,  ONLY : vmine_t, b0e_t, b01e_t, free_e_mine_t
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE data_files,     ONLY : flevdat
USE polyfit_mod,    ONLY : polyfit
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER  :: m1
INTEGER :: idata, ndata, i1
REAL(DP) :: a(poly_degree_ph+1), x(ngeo(1)), y(ngeo(1)), aux, aux1
REAL(DP) :: compute_mur_fun

  IF (my_image_id /= root_image) RETURN

  m1=poly_degree_ph+1
  WRITE(stdout,*) ngeo(1), central_geo
  ndata=0
  DO idata=1,ngeo(1)
     ndata=ndata+1
     x(ndata)=omega_geo(idata)
     y(ndata)=el_free_ener(itemp,idata)
     WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
  ENDDO
  WRITE(stdout,*)
  CALL polyfit(x, y, ndata, a, poly_degree_ph)

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

  vmine_t(itemp)=vm
  b0e_t(itemp)=aux
  b01e_t(itemp)=aux1
  free_e_mine_t(itemp)=emin + compute_mur_fun(vm, vmin, b0/ry_kbar, b01, a, m1)
  !
  CALL compute_celldm_geo(vmine_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("-"))')
  WRITE(stdout,'(5x, "free energy from electron dos, at T= ", f12.6)') &
                                                                temp(itemp)
  IF (pressure_kb /= 0.0_DP) &
     WRITE(stdout, '(5x,"pressure = ",f15.6," kbar")') pressure_kb
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
                                 &" a.u.")') celldm_(1)
  WRITE(stdout,'(5x, "The bulk modulus is ",24x,f12.3,"  kbar")')  b0e_t(itemp)
  WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",&
                               &f9.3)')  b01e_t(itemp)
  IF (pressure_kb /= 0.0_DP) THEN
     WRITE(stdout,'(5x,"The Gibbs energy at the minimum is    ",&
                                    6x,f20.9," Ry")') free_e_mine_t(itemp)
  ELSE
     WRITE(stdout,'(5x,"The free energy at the minimum is",6x,f20.9," Ry")') &
                                                 free_e_mine_t(itemp)
  END IF

  WRITE(stdout,'(2x,76("-"),/)')

  RETURN
END SUBROUTINE do_ev_t_el

!-----------------------------------------------------------------------
SUBROUTINE find_min_mur_pol(v0, b0, b01, a, m1, vm)
!-----------------------------------------------------------------------
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

v1=v0
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

!-----------------------------------------------------------------------
FUNCTION compute_fun(v, v0, b0, b01, a, m1)
!-----------------------------------------------------------------------
!
!  This function computes the first derivative of the Murnaghan equation
!  plus a polynomial of degree m1-1 with respect to the volume
!
USE kinds, ONLY : DP 
IMPLICIT NONE
REAL(DP) :: compute_fun
INTEGER,  INTENT(IN) :: m1
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
!
!-----------------------------------------------------------------------
FUNCTION compute_fun_deriv(v, v0, b0, b01, a, m1)
!-----------------------------------------------------------------------
!
!  This function computes the second derivative of the Murnaghan equation
!  plus a polynomial of degree m1-1 with respect to the volume
!
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

!-----------------------------------------------------------------------
FUNCTION compute_mur_fun(v, v0, b0, b01, a, m1)
!-----------------------------------------------------------------------
!
!   This function computes the sum of a Murnaghan equation and
!   a polynomial. Compute fun computes the derivative with respect
!   to the volume of this function.
!
USE kinds, ONLY : DP 
IMPLICIT NONE
REAL(DP) :: compute_mur_fun
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v, v0, b0, b01, a(m1)
REAL(DP) :: aux
INTEGER :: i1

aux = b0 * v0 * (( v0 / v )**(b01-1.0_DP) / (b01-1.0_DP) + (v/v0) )/b01 &
           - v0 * b0 / (b01-1.0_DP) + a(1) 
DO i1=2, m1
   aux = aux +  a(i1) * v**(i1-1)
ENDDO
compute_mur_fun=aux

RETURN
END FUNCTION compute_mur_fun

!-------------------------------------------------------------------------
SUBROUTINE summarize_mur(celldm, b0, b01, free_e_min)
!-------------------------------------------------------------------------
!
USE kinds,     ONLY : DP
USE control_pressure, ONLY : pressure_kb
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) ::  celldm, b0, b01, free_e_min

IF (pressure_kb /= 0.0_DP) &
   WRITE(stdout, '(5x,"pressure = ",f15.6," kbar")') pressure_kb
WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
                                 &" a.u.")') celldm
WRITE(stdout,'(5x, "The bulk modulus is ",24x,f12.3,"  kbar")')  b0
WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",&
                               &f9.3)')  b01
IF (pressure_kb /= 0.0_DP) THEN
   WRITE(stdout,'(5x,"The Gibbs energy at the minimum is    ",&
                                    6x,f20.9," Ry")') free_e_min
ELSE
   WRITE(stdout,'(5x,"The free energy at the minimum is",6x,f20.9," Ry")') &
                                                 free_e_min
END IF

RETURN
END SUBROUTINE summarize_mur

SUBROUTINE save_free_energy_on_file(itemp)
!
!  This subroutines saves on file the vibrational free energy and
!  (if calculated) the electronic free energy every temp_nstep temperatures.
!
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : omega_geo, ngeo, no_ph
USE temperature,       ONLY : ntemp, temp_nstep
USE thermodynamics,    ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE control_eldos,     ONLY : lel_free_energy
USE data_files,        ONLY : flevdat
USE io_global,         ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filedata
CHARACTER(LEN=6) :: int_to_char
INTEGER, INTENT(IN) :: itemp
INTEGER :: idata, iu_free
INTEGER :: find_free_unit

IF (temp_nstep>ntemp) RETURN
IF (MOD(itemp-1,temp_nstep)>0) RETURN

filedata="anhar_files/"//TRIM(flevdat)//"_free."//TRIM(int_to_char(itemp))
CALL add_pressure(filedata)

IF (ionode) THEN
   iu_free=find_free_unit()
   OPEN(UNIT=iu_free, FILE=TRIM(filedata), STATUS='UNKNOWN', FORM='FORMATTED')
   DO idata=1,ngeo(1)
      IF (no_ph(idata)) CYCLE
      IF (lel_free_energy) THEN
         WRITE(iu_free,'(3f25.14)') omega_geo(idata), &
                        ph_free_ener(itemp,idata), el_free_ener(itemp,idata)
      ELSE
         WRITE(iu_free,'(2f25.14)') omega_geo(idata), &
                                    ph_free_ener(itemp,idata)
      ENDIF
   ENDDO
   CLOSE(UNIT=iu_free, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE save_free_energy_on_file
