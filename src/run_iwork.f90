!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------
SUBROUTINE run_iwork(iwork, part, iq, irr, igeom)
!----------------------------------------------------------------
!
! This subroutine runs pw.x or ph.x for the iwork of part part.
! It also checks if the pw.x calculation has already been done
! and copy the output of pw.x in the variables of thermo_pw.
!
USE control_thermo,       ONLY : lpwscf, lpwband, lstress, lphonon, lberry, &
                                 lef, geometry, all_geometries_together
USE control_qe,           ONLY : use_ph_images
USE control_phrun,        ONLY : auxdyn
USE thermo_mod,           ONLY : energy_geo, ef_geo, iwho
USE elastic_constants,    ONLY : sigma_geo
USE piezoelectric_tensor, ONLY : polar_geo, tot_b_phase, nppl
USE ener,                 ONLY : etot, ef
USE force_mod,            ONLY : sigma
USE freq_ph,              ONLY : fpol
USE io_global,            ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN)             :: iwork, part, iq, irr, igeom
INTEGER                         :: igeom1
LOGICAL                         :: run, exit_status
!
! Initialize run for the band calculation
!
run=.TRUE.
!
!  Start with a message of what is computed
!
IF (lphonon(iwork).AND..NOT.fpol) THEN
   igeom1=igeom
   IF (all_geometries_together) igeom1=geometry(iwork)
   IF (use_ph_images) THEN
      WRITE(stdout,'(5x,"Doing the work of image", i5, " of geometry", i5 )')&
                           iq, igeom1
      CALL print_image_work()
   ELSE
      WRITE(stdout,'(5x,"Doing q-point",&
            &i5," irrep", i5, " of geometry", i5 )') iq, irr, igeom1
   ENDIF
ENDIF
WRITE(stdout,'(2x,76("+"),/)')
!
!  And now do the real calculation. Only one of the three option 
!  must be true and the input has been already set by set_thermo_work_todo.
!  First pw.x
!
IF (lpwscf(iwork).OR.lpwband(iwork)) THEN
   IF (lpwscf(iwork)) CALL check_existence(iwork,part,run)
   IF (run) THEN
      CALL do_pwscf(exit_status, lpwscf(iwork))
      IF (lpwscf(iwork)) energy_geo(iwork)=etot
      IF (lef(iwork)) ef_geo(iwork)=ef
      IF (lstress(iwork)) sigma_geo(:,:,iwork)=sigma(:,:)
      IF (lberry(iwork)) CALL do_berry(exit_status, polar_geo(1,iwork), &
                                 tot_b_phase(1,iwork), nppl)
      IF (lpwscf(iwork)) CALL save_existence(iwork,part)
      IF (lpwscf(iwork)) CALL save_geometry(iwork,part,iwho)
   ENDIF
ENDIF
!
!  pw.x for polarization calculation
!
!
!  A phonon calculation
!
IF (lphonon(iwork)) THEN
   CALL do_phonon_tpw(auxdyn) 
   CALL collect_grid_files_tpw()
   IF (all_geometries_together) CALL close_ph_geometry(.TRUE.)
ENDIF

RETURN
END SUBROUTINE run_iwork

