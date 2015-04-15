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
  !  This routine receives the summary of the lattice constants and of
  !  the energy at each lattice constant and writes the input file 
  !  for the ev_sub subroutine.
  !
  !
  USE kinds, ONLY : DP
  USE mp_images, ONLY : my_image_id, root_image, nproc_image
  USE mp_world,  ONLY : world_comm
  USE thermo_mod, ONLY : omega_geo, energy_geo, ngeo
  USE io_global, ONLY : ionode
  IMPLICIT NONE
  CHARACTER(LEN=256) :: file_dat
  INTEGER :: iu_ev, igeom
  !
  !  First collect the results from all images
  !
  IF (my_image_id /= root_image) RETURN
  !
  CALL write_ev_driver(file_dat) 
  !
  IF (ionode) THEN
     iu_ev=2
     OPEN(UNIT=iu_ev, FILE=TRIM(file_dat), STATUS='UNKNOWN', FORM='FORMATTED')
     DO igeom=1,ngeo
        WRITE(iu_ev,'(2e30.15)') omega_geo(igeom), energy_geo(igeom)
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
USE io_global, ONLY : ionode
IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: file_dat
INTEGER :: iu_ev
!
IF (ionode) THEN
   iu_ev=2
   OPEN(UNIT=iu_ev, FILE='input_ev', STATUS='UNKNOWN', FORM='FORMATTED')
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
!  This subroutine compute the equilibrium volume and bulk modulus
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : omega_geo, alat_geo, ngeo
USE control_mur, ONLY : vmin, b0, b01, emin
USE control_thermo, ONLY : flevdat
USE io_global, ONLY : meta_ionode_id, stdout
USE mp_world, ONLY : world_comm
USE mp, ONLY : mp_bcast

IMPLICIT NONE
CHARACTER(LEN=256) :: file_dat
REAL(DP) :: compute_alat_geo

  file_dat=TRIM(flevdat) 
  CALL write_ev_input(file_dat)
  CALL ev_sub(vmin, b0, b01, emin)
  CALL mp_bcast(vmin, meta_ionode_id, world_comm)
  CALL mp_bcast(b0, meta_ionode_id, world_comm)
  CALL mp_bcast(b01, meta_ionode_id, world_comm)
  CALL mp_bcast(emin, meta_ionode_id, world_comm)
  !
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",f12.4," a.u.")') &
               compute_alat_geo(vmin, alat_geo(ngeo/2+1), omega_geo(ngeo/2+1))
  WRITE(stdout,'(5x, "The bulk modulus is ",15x,f12.2," kbar")')  b0
  WRITE(stdout,'(2x,76("+"),/)')

END SUBROUTINE do_ev

SUBROUTINE do_ev_t(itemp)
!
!  This subroutine compute the equilibrium volume and bulk modulus
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, alat_geo, omega_geo, energy_geo
USE thermodynamics, ONLY : ph_free_ener
USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, free_e_min_t
USE temperature,    ONLY : ntemp, temp
USE control_thermo, ONLY : flevdat
USE io_global,      ONLY : meta_ionode_id, stdout,  meta_ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
CHARACTER(LEN=256) :: file_dat
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: compute_alat_geo
REAL(DP) :: free_e, vm, b0, b01

  IF (my_image_id /= root_image) RETURN

  file_dat=TRIM(flevdat)//TRIM(int_to_char(itemp))
  CALL write_ev_driver(file_dat)

  IF (meta_ionode) THEN
     iu_ev=2
     OPEN(UNIT=iu_ev, FILE=TRIM(file_dat), STATUS='UNKNOWN', FORM='FORMATTED')
     DO igeom=1,ngeo
        WRITE(iu_ev,'(2f20.10)') omega_geo(igeom), energy_geo(igeom) + &
                                            ph_free_ener(itemp,igeom)
     ENDDO
     !
     CLOSE(iu_ev)
     !
  END IF

  CALL ev_sub(vm, b0, b01, free_e)
  vmin_t(itemp)=vm
  b0_t(itemp)=b0
  b01_t(itemp)=b01
  free_e_min_t(itemp)=free_e
  !
  WRITE(stdout,'(/,2x,76("-"))')
  WRITE(stdout,'(5x, "phdos free energy, at T= ", f12.6)') temp(itemp)
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",f12.6," a.u.")') &
      compute_alat_geo(vmin_t(itemp), alat_geo(ngeo/2+1), omega_geo(ngeo/2+1))
  WRITE(stdout,'(5x, "The bulk modulus is ",15x,f12.4," kbar")')  b0_t(itemp)
  WRITE(stdout,'(5x, "The bulk modulus derivative is ",15x,f12.4)')  b01_t(itemp)
  WRITE(stdout,'(2x,76("-"),/)')

END SUBROUTINE do_ev_t

SUBROUTINE do_ev_t_ph(itemp)
!
!  This subroutine compute the equilibrium volume and bulk modulus
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, alat_geo, omega_geo, energy_geo
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE ph_freq_anharmonic,     ONLY : vminf_t, b0f_t, b01f_t, free_e_minf_t
USE temperature,    ONLY : ntemp, temp
USE control_thermo, ONLY : flevdat
USE io_global,      ONLY : meta_ionode_id, stdout,  meta_ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
CHARACTER(LEN=256) :: file_dat
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: free_e, vm, b0, b01
REAL(DP) :: compute_alat_geo

  IF (my_image_id /= root_image) RETURN

  file_dat=TRIM(flevdat)//'_ph'//TRIM(int_to_char(itemp))
  CALL write_ev_driver(file_dat)

  IF (meta_ionode) THEN
     iu_ev=2
     OPEN(UNIT=iu_ev, FILE=TRIM(file_dat), STATUS='UNKNOWN', FORM='FORMATTED')
     DO igeom=1,ngeo
        WRITE(iu_ev,'(2f20.10)') omega_geo(igeom), energy_geo(igeom) + &
                                            phf_free_ener(itemp,igeom)
     ENDDO
     !
     CLOSE(iu_ev)
     !
  END IF

  CALL ev_sub(vm, b0, b01, free_e)
  vminf_t(itemp)=vm
  b0f_t(itemp)=b0
  b01f_t(itemp)=b01
  free_e_minf_t(itemp)=free_e
  !
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x, "ph_freq free energy at T=",f12.4)') temp(itemp)
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",f12.6," a.u.")') &
      compute_alat_geo(vminf_t(itemp), alat_geo(ngeo/2+1), omega_geo(ngeo/2+1))
  WRITE(stdout,'(5x, "The bulk modulus is ",15x,f12.4," kbar")')  b0f_t(itemp)
  WRITE(stdout,'(5x, "The bulk modulus derivative is ",15x,f12.4)')  b01f_t(itemp)
  WRITE(stdout,'(2x,76("+"),/)')

END SUBROUTINE do_ev_t_ph

