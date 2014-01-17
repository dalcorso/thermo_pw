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
  USE thermo_mod, ONLY : alat_geo, energy_geo
  USE thermodynamics, ONLY : ngeo
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
        WRITE(iu_ev,'(2f20.10)') alat_geo(igeom), energy_geo(igeom)
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
   WRITE(iu_ev,'("fcc")')
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
USE thermo_mod, ONLY : vmin, b0, emin
USE control_thermo, ONLY : flevdat
USE io_global, ONLY : meta_ionode_id, stdout
USE mp_world, ONLY : world_comm
USE mp, ONLY : mp_bcast

IMPLICIT NONE
CHARACTER(LEN=256) :: file_dat

  file_dat=TRIM(flevdat)
  CALL write_ev_input(file_dat)
  CALL ev_sub(vmin, b0, emin)
  CALL mp_bcast(vmin, meta_ionode_id, world_comm)
  CALL mp_bcast(b0, meta_ionode_id, world_comm)
  CALL mp_bcast(emin, meta_ionode_id, world_comm)
  !
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",f12.4," a.u.")') &
                                           (vmin * 4.0_DP)**(1.0_DP/3.0_DP)
  WRITE(stdout,'(5x, "The bulk modulus is ",15x,f12.2," kbar")')  b0
  WRITE(stdout,'(2x,76("+"),/)')

END SUBROUTINE do_ev

SUBROUTINE do_ev_t(itemp)
!
!  This subroutine compute the equilibrium volume and bulk modulus
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : alat_geo, energy_geo, vmin_t, b0_t, free_e_min_t
USE control_thermo, ONLY : flevdat
USE io_global, ONLY : meta_ionode_id, stdout,  meta_ionode
USE mp_images, ONLY : my_image_id, root_image
USE mp_world, ONLY : world_comm
USE mp, ONLY : mp_bcast
USE thermodynamics, ONLY : ntemp, temp, ph_free_ener, ngeo

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
CHARACTER(LEN=256) :: file_dat
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: free_e, vm, b0

  IF (my_image_id /= root_image) RETURN

  IF (.NOT.ALLOCATED(vmin_t)) ALLOCATE(vmin_t(ntemp))
  IF (.NOT.ALLOCATED(b0_t)) ALLOCATE(b0_t(ntemp))
  IF (.NOT.ALLOCATED(free_e_min_t)) ALLOCATE(free_e_min_t(ntemp))
  
  file_dat=TRIM(flevdat)//TRIM(int_to_char(itemp))
  CALL write_ev_driver(file_dat)

  IF (meta_ionode) THEN
     iu_ev=2
     OPEN(UNIT=iu_ev, FILE=TRIM(file_dat), STATUS='UNKNOWN', FORM='FORMATTED')
     DO igeom=1,ngeo
        WRITE(iu_ev,'(2f20.10)') alat_geo(igeom), energy_geo(igeom) + &
                                            ph_free_ener(itemp,igeom)
     ENDDO
     !
     CLOSE(iu_ev)
     !
  END IF

  CALL ev_sub(vm, b0, free_e)
  vmin_t(itemp)=vm
  b0_t(itemp)=b0
  free_e_min_t(itemp)=free_e

!  CALL mp_bcast(vmin_t(itemp), meta_ionode_id, world_comm)
!  CALL mp_bcast(b0_t(itemp), meta_ionode_id, world_comm)
!  CALL mp_bcast(free_e_min_t(itemp), meta_ionode_id, world_comm)
  !
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x, "The equilibrium lattice constant at T=",f12.4, &
       & " is ",f12.4," a.u.")') temp(itemp,1), &
                                 (vmin_t(itemp)*4.0_DP)**(1.0_DP/3.0_DP)
  WRITE(stdout,'(5x, "The bulk modulus is ",15x,f12.2," kbar")')  b0_t(itemp)
  WRITE(stdout,'(2x,76("+"),/)')

END SUBROUTINE do_ev_t

