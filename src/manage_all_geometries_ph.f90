!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE manage_all_geometries_ph()
!--------------------------------------------------------------------------
!
!  This routine controls the calculation of the phonon dispersions 
!  when geometries are made in parallel. It first reads the input files
!  of all the phonon calculations and computes all the tasks to do.
!  Then it runs asynchronously all the tasks and finally collects the 
!  results.
!
  USE input_parameters, ONLY : outdir

  USE thermo_mod,       ONLY : what, ibrav_geo, celldm_geo, no_ph,      &
                               start_geometry, last_geometry, phgeo_on_file
  USE control_ph,       ONLY : always_run, ldisp, low_directory_check
  USE cell_base,        ONLY : ibrav, celldm

  USE mp_asyn,          ONLY : with_asyn_images, stop_signal_activated
  USE mp_images,        ONLY : nimage, my_image_id
  USE mp,               ONLY : mp_barrier
  USE mp_world,         ONLY : world_comm
  USE output,           ONLY : fildyn
  USE control_thermo,   ONLY : outdir_thermo, after_disp, set_internal_path, &
                               lq2r
  USE io_global,        ONLY : stdout

IMPLICIT NONE

INTEGER  :: part, nwork, igeom, iaux
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: std, ldcs, after_disp_save, something_todo

CHARACTER (LEN=256) :: auxdyn=' '

always_run=.TRUE.
CALL start_clock( 'PHONON' )
CALL check_phgeo_on_file()
IF (after_disp) THEN
   ldisp=.TRUE.
ELSE
!
!  Initialize the work of all the geometries and analyze what is on disk
!  This routine must be called by all processors
!
   CALL check_geo_initial_status(something_todo)
!
!  Initialize the asynchronous work
!
   IF (something_todo) THEN
      auxdyn=fildyn
      part=2
      CALL initialize_thermo_work(nwork, part, iaux)
!
!  Asynchronous work starts here. No communication is
!  allowed except though the master/slave mechanism
!  
      CALL run_thermo_asynchronously(nwork, part, 1, auxdyn)

      CALL deallocate_asyn()
      IF (stop_signal_activated) GOTO 100
   ENDIF
ENDIF
!
!  Now all calculations are done, we collect the results. 
!  Each image acts independently. The total number of collection tasks
!  are divided between images. The processors must be resynchronized here
!  otherwise some partial dynamical matrix could be missing.
!
CALL mp_barrier(world_comm)
IF (.NOT.(after_disp.AND.(what=='mur_lc_t'.OR. &
                               what=='elastic_constants_t'))) THEN
   DO igeom=start_geometry, last_geometry
      IF (no_ph(igeom)) CYCLE
      IF (phgeo_on_file(igeom)) CYCLE
      CALL check_stc_g(igeom, nimage, my_image_id, std)
      IF (.NOT.std) CYCLE
      WRITE(stdout,'(/,5x,40("%"))') 
      WRITE(stdout,'(5x,"Collecting geometry ", i5)') igeom
      WRITE(stdout,'(5x,40("%"),/)') 
      outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(igeom))//'/'
      !
      ! ... reads the phonon input
      !
      CALL initialize_geometry_and_ph(.TRUE., igeom, auxdyn)
      ldcs=low_directory_check
      low_directory_check=.TRUE.
      CALL check_initial_geometry(auxdyn)
      low_directory_check=ldcs

      CALL manage_collection(auxdyn, igeom)

      CALL close_ph_geometry(.TRUE.)
   ENDDO
   CALL restore_files_names()
ENDIF
!
!  resynchronize all processors, otherwise some dynamical matrix could be
!  missing. The calculation of the thermodynamical properties is
!  parallelized over all processors, so the following routines must
!  be called by all images.
!
CALL mp_barrier(world_comm)

after_disp_save=after_disp
DO igeom=start_geometry, last_geometry
   IF (no_ph(igeom)) CYCLE
   after_disp=after_disp_save
   IF (phgeo_on_file(igeom)) after_disp=.TRUE.
   WRITE(stdout,'(/,5x,40("%"))') 
   WRITE(stdout,'(5x,"Computing thermodynamic properties", i5)') igeom
   WRITE(stdout,'(5x,40("%"),/)') 
   outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(igeom))//'/'
   !  
   !  The geometry must be reset here
   !
   ibrav=ibrav_geo(igeom)
   celldm(:)=celldm_geo(:,igeom)
   IF (set_internal_path) CALL set_bz_path()
   IF (after_disp.AND.igeom==start_geometry) CALL initialize_file_names()

   CALL set_files_names(igeom)
   auxdyn=fildyn
!
!  Compute the dispersions and the thermodynamic properties
!
   IF (lq2r) CALL manage_ph_dispersions(auxdyn, igeom)

   CALL restore_files_names()
ENDDO
100 CONTINUE
after_disp=after_disp_save
CALL restore_files_names()

RETURN
END SUBROUTINE manage_all_geometries_ph
