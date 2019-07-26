!
! Copyright (C) 2017 Andrea Dal Corso
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
!  when geometries are made in sequence
!
  USE input_parameters, ONLY : outdir

  USE ions_base,        ONLY : nat
  USE thermo_mod,       ONLY : what, ibrav_geo, celldm_geo, no_ph,      &
                               max_geometries, start_geometry, last_geometry
  USE control_ph,       ONLY : with_ext_images, always_run, ldisp, trans, &
                               rec_code_read, xmldyn
  USE initial_conf,     ONLY : collect_info_save
  USE disp,             ONLY : nqs, comp_iq, done_iq
  USE grid_irr_iq,      ONLY : comp_irr_iq, done_irr_iq, irr_iq
  USE ph_restart,       ONLY : destroy_status_run
  USE freq_ph,          ONLY : fpol
  USE save_ph,          ONLY : clean_input_variables
  USE control_lr,       ONLY : lgamma
  USE cell_base,        ONLY : ibrav, celldm


  USE mp_asyn,          ONLY : with_asyn_images, stop_signal_activated

  USE io_global,        ONLY : meta_ionode, stdout, meta_ionode_id
  USE output,           ONLY : fildyn
  USE control_thermo,   ONLY : outdir_thermo, after_disp, set_internal_path, &
                               lq2r, lectqha
  USE control_qe,       ONLY : use_ph_images
  USE collect_info,     ONLY : comm_collect_info, init_collect_info, &
                               destroy_collect_info_type
  USE mp_images,        ONLY : nimage, my_image_id, inter_image_comm
  USE mp_world,         ONLY : world_comm
  USE mp,               ONLY : mp_sum, mp_barrier

IMPLICIT NONE

INTEGER  :: part, nwork, igeom, exit_status, ph_geometries, iaux
LOGICAL  :: check_dyn_file_exists, do_ph, has_xml, fninit, after_disp_save
CHARACTER(LEN=6) :: int_to_char

CHARACTER (LEN=256) :: auxdyn

with_ext_images=with_asyn_images
ph_geometries=0
always_run=.TRUE.
CALL start_clock( 'PHONON' )
IF (use_ph_images) ALLOCATE(collect_info_save(1))
fninit=.FALSE.
after_disp_save=after_disp
DO igeom=start_geometry,last_geometry
   IF (no_ph(igeom).OR.stop_signal_activated) CYCLE
   after_disp=after_disp_save
   auxdyn=TRIM(fildyn)//'.g'//TRIM(int_to_char(igeom))//'.'
   IF (auxdyn(1:18)/='dynamical_matrices') &
          auxdyn='dynamical_matrices/'//TRIM(auxdyn)
   IF (check_dyn_file_exists(auxdyn)) after_disp=.TRUE.
   auxdyn=' '
   WRITE(stdout,'(/,5x,40("%"))') 
   WRITE(stdout,'(5x,"Computing geometry ", i5)') igeom
   WRITE(stdout,'(5x,40("%"),/)') 
   outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(igeom))//'/'
   !
   ! ... reads the phonon input
   !
   IF (after_disp) THEN
      ldisp=.TRUE.
   ELSE
      CALL thermo_ph_readin()
      CALL save_ph_variables()
   ENDIF
   IF (.NOT.fninit) THEN
      CALL initialize_file_names()
      fninit=.TRUE.
   ENDIF
   CALL set_files_names(igeom)

   auxdyn=fildyn
   IF (after_disp) THEN
!
!  The geometry is read by thermo_ph_readin from the output files of pw.x,
!  except in the case where after_disp=.TRUE.. In this case we have to
!  set it here.
!
      IF (what=='mur_lc_t'.OR.what=='elastic_constants_t') THEN
         ibrav=ibrav_geo(igeom)
         celldm(:)=celldm_geo(:,igeom)
      ENDIF
      IF (set_internal_path) CALL set_bz_path()
      do_ph=.FALSE.
   ELSE
      do_ph=.TRUE.
      IF (trans) do_ph=.NOT. check_dyn_file_exists(auxdyn) 
      IF (lectqha.AND.set_internal_path) CALL set_bz_path()
   ENDIF
!
!  Set the BZ path for the present geometry
!
   CALL set_paths_disp()

   IF ( do_ph ) THEN
      ph_geometries=ph_geometries+1
      IF (ph_geometries > max_geometries) THEN
         WRITE(stdout,'(5x,"The code stops because max_geometries is",&
                               &i4)') max_geometries
         GOTO 1000
      ENDIF
   !
   ! ... Checking the status of the calculation and if necessary initialize
   ! ... the q mesh and all the representations
   !
      CALL check_initial_status_tpw(auxdyn)
!
!  With images the recover is done only at the level of complete tasks
!  so rec_code_read is reinitialized after the writing of the modes.
!
      IF (with_asyn_images) rec_code_read=-40
      !
      IF (use_ph_images) THEN
         CALL init_collect_info(collect_info_save(1), nqs, nat, &
                nimage, my_image_id+1, comp_irr_iq, done_irr_iq, comp_iq, &
                done_iq, irr_iq)

         CALL comm_collect_info(collect_info_save(1), inter_image_comm)
      ENDIF
      !
      part=2
      CALL initialize_thermo_work(nwork, part, iaux)
      !
      !  Asyncronous work starts again. No communication is
      !  allowed except though the master workers mechanism
      !
      CALL run_thermo_asynchronously(nwork, part, igeom, auxdyn)
      IF (stop_signal_activated) GOTO 100
      !  
      !   return to syncronous work. Collect the work of all images and
      !   writes the dynamical matrix
      !
      CALL mp_barrier(world_comm)
      IF (trans) THEN
         IF (with_asyn_images) CALL collect_everything(auxdyn, igeom)
      ELSEIF (fpol) THEN
         IF (lgamma) THEN
            IF (nimage>1) CALL collect_all_epsilon()
            CALL plot_epsilon_omega_opt()
         ELSE
            IF (nimage>1) CALL collect_all_chi()
            CALL plot_epsilon_omega_q()
         ENDIF
      ENDIF
      !
   END IF

   IF (lq2r) CALL manage_ph_dispersions(auxdyn, igeom)

100 CALL deallocate_asyn()
   IF (.NOT. after_disp) THEN
      CALL clean_pw(.TRUE.)
      CALL close_phq(.FALSE.)
      CALL clean_input_variables()
      CALL destroy_status_run()
      CALL deallocate_ph_tpw()
      CALL deallocate_part()
      IF (use_ph_images) CALL destroy_collect_info_type(collect_info_save(1))
   ELSE
      CALL restore_files_names()
   ENDIF
ENDDO
after_disp=after_disp_save
IF (.NOT.fninit) THEN
   CALL initialize_file_names()
   fninit=.TRUE.
ENDIF
CALL restore_files_names()
1000 CONTINUE

RETURN
END SUBROUTINE manage_ph

