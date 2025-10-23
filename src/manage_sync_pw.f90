!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------
SUBROUTINE manage_sync_pw()
!-------------------------------------------------
!
!  This driver runs a self consistent calculation of pw.x 
!  possibly followed by a pw.x band structure calculation if required 
!  by the task. Only the first image does the calculation. The
!  relaxed coordinates are sent to all the other images.
!
  USE kinds,            ONLY : DP
  USE control_thermo,   ONLY : lbands_syn_1
  USE ions_base,        ONLY : tau
  USE cell_base,        ONLY : at, omega, celldm
  USE control_atomic_pos, ONLY : linternal_thermo
  USE input_parameters, ONLY : outdir

  USE io_files,         ONLY : tmp_dir, wfc_dir, check_tempdir
  USE control_pwrun,    ONLY : do_punch
  USE control_2d_bands, ONLY : only_bands_plot
  USE control_xrdp,     ONLY : lxrdp
  USE io_global,        ONLY : stdout, meta_ionode_id
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE mp_images,        ONLY : nimage, my_image_id, root_image
  USE mp_asyn,          ONLY : with_asyn_images

IMPLICIT NONE

INTEGER  :: exit_status
LOGICAL  :: run
CHARACTER (LEN=80) :: message

   with_asyn_images=.FALSE.
   CALL set_outdir_name(1)
   CALL set_tmp_dir(outdir)

   IF (my_image_id==root_image) THEN
!
!   do the self consistent calculation at the new lattice constant
!
      do_punch=.TRUE.
      IF (.NOT.only_bands_plot) THEN
         WRITE(message,'(5x,"Doing a self-consistent calculation", i5)') 
         CALL decorated1_write(message)
         CALL check_existence(0,1,run)
         IF (run) THEN
!
!    with linternal_thermo relaxation was disabled. Here we enable it
!
            IF (linternal_thermo) CALL initialize_relaxation()
            CALL do_pwscf(exit_status, .TRUE.)
            CALL save_existence(0,1)
         END IF

         IF (lxrdp) CALL manage_xrdp('.scf')

      ENDIF
   ENDIF
   IF (lbands_syn_1) CALL manage_bands()
   CALL mp_bcast(tau, meta_ionode_id, world_comm)
   CALL mp_bcast(celldm, meta_ionode_id, world_comm)
   CALL mp_bcast(at, meta_ionode_id, world_comm)
   CALL mp_bcast(omega, meta_ionode_id, world_comm)
   CALL set_equilibrium_conf(celldm, tau, at, omega)
   with_asyn_images=(nimage>1)
   !
   RETURN
END SUBROUTINE manage_sync_pw
