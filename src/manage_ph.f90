!
! Copyright (C) 2017-2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE manage_ph()
!--------------------------------------------------------------------------
!
!  This routine controls the calculation of the phonon dispersions 
!  and harmonic thermodynamic quantities when geometries are made 
!  in sequence
!
! thermo_pw variables
!
  USE thermo_mod,       ONLY : no_ph, max_geometries, start_geometry, &
                               last_geometry, dynmat_on_file
  USE control_thermo,   ONLY : lq2r
  USE control_phrun,    ONLY : auxdyn
  USE initial_conf,     ONLY : collect_info_save
  USE collect_info,     ONLY : destroy_collect_info_type
!
!  variables of pw.x or of the phonon
!
  USE io_global,        ONLY : stdout
  USE output,           ONLY : fildyn
  USE control_ph,       ONLY : always_run
!
!  parallel control
!
  USE mp_asyn,          ONLY : stop_signal_activated
  USE mp_world,         ONLY : world_comm
  USE mp,               ONLY : mp_barrier

IMPLICIT NONE

INTEGER  :: part, nwork, igeom, ph_geometries
LOGICAL  :: do_ph, check_dyn_file_exists
CHARACTER(LEN=80) :: message

always_run=.TRUE.
CALL start_clock( 'PHONON' )
ALLOCATE(collect_info_save(1))
!
!  check which dynamical matrices are already on file. Needs fildyn,
!  so it works only with after_disp or if one sets fildyn in the
!  thermo_control file.
!
CALL check_dynmat_all_geo_on_file()
!
!  main loop on the geometries, one after the other
!
ph_geometries=0
DO igeom=start_geometry,last_geometry
!
!  first check all the possibilities of not computing this geometry
!
   IF (no_ph(igeom).OR.stop_signal_activated) CYCLE  ! skip this geometry
   IF (dynmat_on_file(igeom)) GOTO 50  ! skip but compute the thermodynamic
   IF (ph_geometries+1 > max_geometries) THEN
      WRITE(stdout,'(5x,"The code stops because max_geometries is",&
                            &i4)') max_geometries
      GOTO 1000   ! skip all the other geometries
   ENDIF
   !
   ! The phonons of this geometry are calculated
   !
   WRITE(message,'(5x,"Computing geometry ", i5)') igeom
   CALL decorated_write(message)
   !
   ! This command selects where the phonon reads the pw.x output
   ! and allows to run the phonon at different geometries
   !
   CALL set_outdir_name(igeom)
   !
   ! Read the pw.x output and the phonon input
   !
   CALL thermo_ph_readin()
   CALL save_ph_variables()
   !
   ! If all the dynamical matrices for this geometry are on file, skip
   ! the calculation. We recheck here because now we know fildyn.
   !
   CALL set_fildyn_name(igeom)
   IF (check_dyn_file_exists(fildyn)) GOTO 100  ! clean pw and ph variables
   !                                            ! and compute thermodynamic
   !
   ! Increase the counter of the geometries calculated in this run
   !
   ph_geometries=ph_geometries+1
   !
   ! Checking the status of the calculation and if necessary initialize
   ! the q mesh and all the representations. Then save the information
   ! on the collect_info_save structure
   !
   CALL check_initial_status_tpw(auxdyn,0)
   !
   CALL collect_the_info(collect_info_save(1))
   !
   ! Now initialize the asynchronous work
   !
   part=2
   CALL initialize_thermo_work(nwork, part)
   !
   ! Asynchronous work starts here. No communication is allowed except 
   ! through the master/slave mechanism
   !
   CALL run_thermo_asynchronously(nwork, part, igeom)
   !
   IF (stop_signal_activated) GOTO 100  ! clean pw.x and ph.x variables
   !  
   ! Return to synchronous work. Collect the work of all images and
   ! write the dynamical matrix on file. The collection is divided
   ! among images
   !
   CALL mp_barrier(world_comm)
   CALL divide_collection_work(igeom)
   CALL manage_collection(auxdyn, igeom)
   CALL clean_collection_work()
   !
100 CONTINUE
   !
   ! finally we close the ph calculation.
   !
   CALL close_ph_geometry(.TRUE.)
   CALL destroy_collect_info_type(collect_info_save(1))
50 CONTINUE
   !
   ! Now we set all the files names with the number of this geometry and 
   ! compute and write the harmonic thermodynamic properties
   !
   IF (.NOT.stop_signal_activated) THEN
      CALL set_files_names(igeom)
      auxdyn=fildyn
      IF (lq2r) CALL manage_ph_postproc(igeom)
   ENDIF
ENDDO  ! loop on geometries
CALL restore_files_names()
1000 CONTINUE

RETURN
END SUBROUTINE manage_ph
