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
                               max_geometries, start_geometry, last_geometry
  USE control_ph,       ONLY : always_run, ldisp, trans, low_directory_check
  USE freq_ph,          ONLY : fpol
  USE control_lr,       ONLY : lgamma
  USE cell_base,        ONLY : ibrav, celldm

  USE mp_asyn,          ONLY : with_asyn_images, stop_signal_activated
  USE mp_images,        ONLY : nimage, my_image_id
  USE mp,               ONLY : mp_barrier
  USE mp_world,         ONLY : world_comm
  USE output,           ONLY : fildyn
  USE control_thermo,   ONLY : outdir_thermo, after_disp, set_internal_path, &
                               lq2r, lectqha
  USE io_global,        ONLY : stdout

IMPLICIT NONE

INTEGER  :: part, nwork, igeom, ph_geometries, iaux
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: std, ldcs

CHARACTER (LEN=256) :: auxdyn=' '

ph_geometries=0
always_run=.TRUE.
CALL start_clock( 'PHONON' )
IF (after_disp) THEN
   ldisp=.TRUE.
ELSE
!
!  Initialize the work of all the geometries and analyze what is on disk
!  This routine must be called by all processors
!
   CALL check_geo_initial_status()
!
!  Initialize the asynchronous work
!
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
      IF (trans) THEN
         IF (with_asyn_images) CALL collect_everything(auxdyn, igeom)
      ELSEIF (fpol) THEN
         IF (lgamma) THEN
            CALL collect_all_epsilon()
            CALL plot_epsilon_omega_opt()
         ELSE
            CALL collect_all_chi()
            CALL plot_epsilon_omega_q()
         ENDIF
      ENDIF
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

DO igeom=start_geometry, last_geometry
   IF (no_ph(igeom)) CYCLE
   WRITE(stdout,'(/,5x,40("%"))') 
   WRITE(stdout,'(5x,"Computing thermodynamic properties", i5)') igeom
   WRITE(stdout,'(5x,40("%"),/)') 
   outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(igeom))//'/'
   !
   ! ... reads the phonon input
   !
   IF (after_disp.AND.(what=='mur_lc_t'.OR.&
                                what=='elastic_constants_t')) THEN
!
!  The geometry is read by thermo_ph_readin from the output files of pw.x,
!  except in the case where after_disp=.TRUE.. In this case we have to
!  set it here.
!
      ibrav=ibrav_geo(igeom)
      celldm(:)=celldm_geo(:,igeom)
      IF (set_internal_path) CALL set_bz_path()
      IF (igeom==1) CALL initialize_file_names()
   ELSE
      IF (lectqha.AND.set_internal_path) CALL set_bz_path()
   ENDIF
   CALL set_files_names(igeom)
   auxdyn=fildyn
!
!  Set the BZ path for the present geometry
!
   CALL set_paths_disp()
!
!  And compute the thermodynamic properties
!
   IF (lq2r) CALL manage_ph_dispersions(auxdyn, igeom)

ENDDO
100 CONTINUE
CALL restore_files_names()

RETURN
END SUBROUTINE manage_all_geometries_ph

SUBROUTINE initialize_geometry_and_ph(recover_, igeom, auxdyn)
!
!  This routine substitutes phq_readin. It reads the output of pw.x and
!  initializes the same variables initialized by phq_readin, but does not 
!  read the ph_control input file.
!  It must be called by all processors of an image.
!
USE input_parameters, ONLY : outdir
USE control_ph,       ONLY : tmp_dir_ph, tmp_dir_phq, rec_code_read, recover
USE ph_restart,       ONLY : ph_readfile
USE output,           ONLY : fildyn
USE save_ph,          ONLY : save_ph_input_variables, tmp_dir_save
USE io_files,         ONLY : tmp_dir, check_tempdir
USE mp_images,        ONLY : my_image_id
USE mp_pools,         ONLY : kunit
USE ions_base,        ONLY : nat

IMPLICIT NONE
LOGICAL, INTENT(IN) :: recover_
INTEGER, INTENT(IN) :: igeom

CHARACTER (LEN=256) :: auxdyn
LOGICAL :: exst, parallelfs
INTEGER :: ierr, ierr1
CHARACTER(LEN=256), EXTERNAL :: trimcheck
CHARACTER(LEN=6) :: int_to_char

kunit=1
tmp_dir=trimcheck(outdir)
tmp_dir_save=tmp_dir
tmp_dir_ph= TRIM (tmp_dir) // '_ph' // TRIM(int_to_char(my_image_id)) //'/'
CALL check_tempdir ( tmp_dir_ph, exst, parallelfs )
tmp_dir_phq=tmp_dir_ph

CALL read_file()
CALL allocate_part(nat)
CALL allocate_ph_tpw()
CALL save_ph_input_variables()
CALL set_files_names(igeom)
auxdyn=fildyn
tmp_dir=tmp_dir_save

rec_code_read=-1000
IF (recover_) THEN
   recover=.TRUE.
   CALL ph_readfile('init', 0, 0, ierr)
   CALL ph_readfile('status_ph', 0, 0, ierr1)
ENDIF

RETURN
END SUBROUTINE initialize_geometry_and_ph
