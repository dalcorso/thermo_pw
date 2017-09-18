!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE check_geo_initial_status()
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
!
USE input_parameters, ONLY : outdir

USE thermo_mod,       ONLY : no_ph, start_geometry, last_geometry, &
                             tot_ngeo
USE ions_base,        ONLY : nat
USE disp,             ONLY : nqs, comp_iq, done_iq
USE grid_irr_iq,      ONLY : comp_irr_iq, done_irr_iq, irr_iq

USE output,           ONLY : fildyn
USE control_thermo,   ONLY : outdir_thermo
USE control_qe,       ONLY : use_ph_images
USE control_ph,       ONLY : recover

USE initial_conf,     ONLY : collect_info_save
USE collect_info,     ONLY : comm_collect_info, init_collect_info

USE mp_images,        ONLY : inter_image_comm, nimage, my_image_id
USE mp,               ONLY : mp_barrier
USE mp_world,         ONLY : world_comm


IMPLICIT NONE

INTEGER  :: igeom
INTEGER  :: nima, pos
CHARACTER(LEN=6) :: int_to_char
CHARACTER (LEN=256) :: auxdyn=' '

ALLOCATE(collect_info_save(tot_ngeo))
!
!  loop on all the geometries calculated in this run
!
DO igeom=start_geometry,last_geometry
   IF (no_ph(igeom)) CYCLE
   outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(igeom))//'/'
   !
   ! ... reads the phonon input and generate q grid and modes
   !
   IF (igeom==start_geometry) THEN
      CALL thermo_ph_readin()
      CALL save_ph_variables()
      CALL initialize_file_names()
      CALL set_files_names(igeom)

      auxdyn=fildyn
      CALL check_initial_geometry(auxdyn)
   ELSE
      CALL initialize_geometry_and_ph(recover, igeom, auxdyn)
      CALL check_initial_geometry(auxdyn)
   ENDIF
   !
   ! collect the info on what has to be calculated and what already exists.
   !
   IF (use_ph_images) THEN
      pos=my_image_id+1
      nima=nimage
   ELSE
      pos=1
      nima=1
   ENDIF
   CALL init_collect_info(collect_info_save(igeom), nqs, nat, nima, pos, &
                         comp_irr_iq, done_irr_iq, comp_iq, done_iq, irr_iq)

   CALL close_ph_geometry(.FALSE.)
ENDDO
!
!  If the phonon images are used, all images must have the same information 
!  on what must be calculated by each image, so the collect_info structure 
!  is shared between all images for each geometry. 
!  Otherwise here there is a syncronization of all processors.
!
IF (use_ph_images) THEN 
   DO igeom=start_geometry, last_geometry
      CALL comm_collect_info(collect_info_save(igeom), inter_image_comm)
   ENDDO
ELSE
   CALL mp_barrier(world_comm)
ENDIF

CALL restore_files_names()

RETURN
END SUBROUTINE check_geo_initial_status

SUBROUTINE save_ph_variables()

USE control_ph,   ONLY : epsil, zeu, zue
USE initial_conf, ONLY : epsil_save, zeu_save, zue_save

IMPLICIT NONE

epsil_save=epsil
zeu_save=zeu
zue_save=zue

RETURN
END SUBROUTINE save_ph_variables
