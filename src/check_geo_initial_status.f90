!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE check_geo_initial_status(something_todo)
!--------------------------------------------------------------------------
!
!  This routine analyses in sequence all the xml files written by
!  pw.x, generates a mesh of q vectors and a set of displacement patterns 
!  for each geometry if they are not already on file, or if recover
!  is false. Moreover it collects the informations on the calculations
!  that need to be done in each geometry and the calculations that have
!  been already done. It does this calling the routine 
!  check_initial_geometry in each geometry (which is a shortened version
!  of check_initial_status) and collecting all the results in
!  collect_info_save.
!  It is not parallelized over images so it must be called by all processors. 
!  Some communications occurs via the world_comm, since
!  all images must have the same set of displacement patterns and of q vectors.
!
!  If use_ph_images is .TRUE. it uses the phonon image_q_irr_tpw routine 
!  to compute the set of q and irreducible representations that must be 
!  calculated by each image.
!
!  The output of this routine are the data on the collect_info_type
!  See the collect_info.f90 
!
!  collect_info_save
!
USE thermo_mod,       ONLY : no_ph, start_geometry, last_geometry, &
                             tot_ngeo, phgeo_on_file
USE initial_conf,     ONLY : collect_info_save

USE control_qe,       ONLY : use_ph_images
USE control_ph,       ONLY : recover

USE mp,               ONLY : mp_barrier
USE mp_world,         ONLY : world_comm

IMPLICIT NONE

LOGICAL, INTENT(OUT)  :: something_todo

INTEGER             :: igeom
CHARACTER (LEN=256) :: auxdyn=' '
LOGICAL             :: fninit

ALLOCATE(collect_info_save(tot_ngeo))
!
!  loop on all the geometries calculated in this run
!
fninit=.FALSE.
something_todo=.FALSE.
DO igeom=1, tot_ngeo
   collect_info_save(igeom)%nqs=0
ENDDO
!
!  loop on the geometries
!
DO igeom=start_geometry,last_geometry
   IF (no_ph(igeom)) CYCLE
   IF (phgeo_on_file(igeom)) CYCLE
   something_todo=.TRUE.
   CALL set_outdir_name(igeom)
   !
   !    reads the pw.x output, the phonon input and generate 
   !    the grid of q points and modes. 
   !    For the fist geometry we read also ph_control for the info 
   !    on the q point mesh
   !
   IF (.NOT.fninit) THEN
      CALL thermo_ph_readin()
      CALL save_ph_variables()
      fninit=.TRUE.
   ELSE
      CALL fast_phq_readin(recover, igeom)
   ENDIF
   CALL set_fildyn_name(igeom)
   CALL check_initial_geometry(auxdyn,0)
   !
   ! collect the info on what has to be calculated and what already exists.
   !
   CALL collect_the_info(collect_info_save(igeom))
   !
   CALL close_ph_geometry(.FALSE.)
ENDDO
!
!  here there is a syncronization of all processors.
!
CALL mp_barrier(world_comm)
CALL restore_files_names()

RETURN
END SUBROUTINE check_geo_initial_status
