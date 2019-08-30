!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE run_thermo_asynchronously(nwork, part, igeom, auxdyn)
  USE kinds,           ONLY : DP
  USE mp,              ONLY : mp_bcast
  USE mp_world,        ONLY : world_comm, nproc
  USE io_global,       ONLY : ionode, ionode_id, meta_ionode_id, stdout
  USE mp_images,       ONLY : nimage, root_image, my_image_id, & 
                              intra_image_comm
  USE mp_asyn,         ONLY : asyn_master_init, asyn_worker_init, &
                              asyn_close, asyn_master, asyn_worker, &
                              asyn_master_work, with_asyn_images, &
                              asyn_stop, stop_signal_activated
  USE thermo_mod,      ONLY : energy_geo
  USE control_thermo,  ONLY : lpwscf, lstress, lphonon, lberry, &
                              all_geometries_together, max_seconds_tpw, &
                              geometry
  USE control_qe,      ONLY : use_ph_images
  USE elastic_constants, ONLY : sigma_geo
  USE piezoelectric_tensor, ONLY : polar_geo, nppl
  USE ener,            ONLY : etot
  USE force_mod,       ONLY : sigma
  USE freq_ph,         ONLY : fpol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nwork, part, igeom
  CHARACTER (LEN=256), INTENT(IN) :: auxdyn

  INTEGER :: iq, irr, igeom1
  INTEGER, ALLOCATABLE :: proc_num(:)
  INTEGER :: proc_per_image, iwork, image, exit_status
  LOGICAL :: all_done_asyn, run
  CHARACTER (LEN=256) :: auxdyn_loc
  REAL(DP) :: get_clock
  !
  IF ( nwork == 0 ) RETURN
  auxdyn_loc=auxdyn
  !
  !  Two calculations are now possible. If with_asyn_images is .TRUE. each
  !  image does a different task, and the first image is the master that 
  !  keeps into account how much work there is to do and who does what.
  !  Otherwise the works are done one after the other in the traditional
  !  fashion
  !
  with_asyn_images = ( nimage > 1 )
  IF ( with_asyn_images .AND. nwork > 1) THEN
     IF (my_image_id == root_image) THEN
  !
  !  This is the master and it has to rule all the other images
  !
        all_done_asyn=.FALSE.
        ALLOCATE(proc_num(0:nimage-1))
        proc_num(0)=0 
        proc_per_image = nproc / nimage
        DO image=1,nimage-1
           proc_num(image) = proc_num(image-1) + proc_per_image
        ENDDO
  !
  !    and initialize the asynchronous communication
  !
        IF (ionode) CALL asyn_master_init(nimage, nwork, proc_num, world_comm)
  !
  !    enter into a loop of listening, answering, and working
  !
        iwork=1
        DO WHILE ( .NOT. all_done_asyn )
!
!  The root processor of this image acts as the master. 
!  See if some worker is ready to work and give it something
!  to do.
!
           IF (ionode) CALL asyn_master(all_done_asyn)
           CALL mp_bcast(all_done_asyn, ionode_id, intra_image_comm) 
!
!   If the maximum cpu time is elapsed try to stop everything
!
           IF (ionode) THEN
              IF (get_clock('THERMO_PW')> max_seconds_tpw.AND. &
                 .NOT. stop_signal_activated) CALL asyn_stop()
           ENDIF
           IF (iwork > 0) THEN
              all_done_asyn=.FALSE.
!
!          Now also the master can do something. Ask for some work.
!
              IF (ionode) CALL asyn_master_work(iwork)
              CALL mp_bcast(iwork, ionode_id, intra_image_comm) 
!
!   And do the work
!
              IF (iwork>0) THEN
                 CALL set_thermo_work_todo(iwork, part, iq, irr, auxdyn_loc)
                 IF (lpwscf(iwork)) THEN
                    WRITE(stdout,'(/,2x,76("+"))')
                    WRITE(stdout,'(5x,"I am the master and now I do geometry",&
                                       & i5," /", i5)') iwork, nwork
                    WRITE(stdout,'(2x,76("+"),/)')
                 ELSE IF (lphonon(iwork)) THEN
                    WRITE(stdout,'(/,2x,76("+"))')
                    IF (fpol) THEN
                       WRITE(stdout,'(5x,"I am the master and now I do &
                            &frequency", i5," /", i5)') iwork, nwork
                    ELSE
                       IF (all_geometries_together) THEN
                          igeom1=geometry(iwork)
                       ELSE
                          igeom1=igeom
                       ENDIF
                       IF (use_ph_images) THEN
                          WRITE(stdout,'(5x,"Work",i5,"/ ",i5)') iwork, nwork 
                          WRITE(stdout,'(5x,"I am the master and now I &
                       &do the work of image",i5," of geometry",i5)') iq, igeom1
                          CALL print_image_work()
                       ELSE
                          WRITE(stdout,'(5x,"Work",i5,"/ ",i5)') iwork, nwork 
                          WRITE(stdout,'(5x,"I am the master and now I &
                                     &do point",i5," irrep", i5, " of  &
                                      geometry", i5 )') iq, irr, igeom1
                       ENDIF
                    ENDIF
                    WRITE(stdout,'(2x,76("+"),/)')
                 END IF
                 IF (lpwscf(iwork)) THEN
                    CALL check_existence(iwork,part,igeom,run)
                    IF (run) THEN
                       CALL do_pwscf(exit_status, .TRUE.)
                       energy_geo(iwork)=etot
                       IF (lstress(iwork)) THEN
                          sigma_geo(:,:,iwork)=sigma(:,:)
                       ENDIF
                       CALL save_existence(iwork,part,igeom)
                    ENDIF
                 ENDIF
                 IF (lberry(iwork)) CALL do_berry(exit_status, &
                                       polar_geo(1,iwork),nppl)
                 IF (lphonon(iwork)) THEN
                    CALL do_phonon_tpw(auxdyn_loc) 
                    CALL collect_grid_files_tpw()
                    IF (all_geometries_together) CALL close_ph_geometry(.TRUE.)
                 ENDIF
              ENDIF
           ENDIF
           !
        ENDDO
        !
        DEALLOCATE(proc_num)
        IF (ionode) CALL asyn_close()
     ELSE
!
!  This is a worker and asks the master what to do.
!  First initializes the worker stuff. It declares the identity of the master 
!  and which communicator to use
!
        IF (ionode) CALL asyn_worker_init(meta_ionode_id, world_comm)
        iwork=1
        DO WHILE (iwork > 0)
!
!       The root_processor of each image asks the master for some work
!       to do and sends the info to all processors of its image
!       This is a blocking request, all the image blocks here until the
!       master answers.
!
           IF (ionode) CALL asyn_worker(iwork)
           CALL mp_bcast(iwork, root_image, intra_image_comm) 
!
!       and then do the work
!
           IF (iwork>0) THEN
              CALL set_thermo_work_todo(iwork, part, iq, irr, auxdyn_loc)
              IF (lpwscf(iwork)) THEN
                 WRITE(stdout,'(/,2x,76("+"))')
                 WRITE(stdout,'(5x,"I am image ", i5, " and now I do &
                           &geometry", i5," /",i5)') my_image_id, iwork, nwork
                 WRITE(stdout,'(2x,76("+"),/)')
              ELSE IF (lphonon(iwork)) THEN
                 WRITE(stdout,'(/,2x,76("+"))')
                 IF (all_geometries_together) THEN
                    igeom1=geometry(iwork)
                 ELSE
                    igeom1=igeom
                 ENDIF
                 IF (fpol) THEN
                    WRITE(stdout,'(5x,"I am image ",i5," and now I do &
                     &frequency", i5, " /", i5 )') my_image_id, iwork, nwork
                 ELSE
                    IF (use_ph_images) THEN
                       WRITE(stdout,'(5x,"Work",i5,"/ ",i5)') iwork, nwork 
                       WRITE(stdout,'(5x,"I am image ",i5," and now I do &
                          the work of image", i5, " of geometry", i5 )') &
                           my_image_id, iq, igeom1
                       CALL print_image_work()
                    ELSE
                       WRITE(stdout,'(5x,"Work",i5,"/ ",i5)') iwork, nwork 
                       WRITE(stdout,'(5x,"I am image ",i5," and now I do point", &
                     i5," irrep", i5, " of geometry", i5 )') my_image_id, iq, &
                                                          irr, igeom1 
                    END IF
                 END IF
                 WRITE(stdout,'(2x,76("+"),/)')
              END IF
 
              IF (lpwscf(iwork)) THEN
                 WRITE(stdout,'(/,2x,76("+"))')
                 CALL check_existence(iwork,part,igeom,run)
                 IF (run) THEN
                    CALL do_pwscf(exit_status, .TRUE.)
                    energy_geo(iwork)=etot
                    IF (lstress(iwork)) THEN
                       sigma_geo(:,:,iwork)=sigma(:,:)
                    ENDIF
                    CALL save_existence(iwork,part,igeom)
                 ENDIF
                 WRITE(stdout,'(2x,76("+"),/)')
              END IF
              IF (lberry(iwork)) CALL do_berry(exit_status, &
                                       polar_geo(1,iwork), nppl)
              IF (lphonon(iwork)) THEN
                 CALL do_phonon_tpw(auxdyn_loc) 
                 CALL collect_grid_files_tpw()
                 IF (all_geometries_together) CALL close_ph_geometry(.TRUE.)
              END IF
           END IF
        END DO
     END IF
  ELSE
!
!  This is the standard case. Asynchronous images are not used. There is
!  only the master that does all the works one after the other.
!
     IF (my_image_id == root_image) THEN
        DO iwork = 1, nwork
           CALL set_thermo_work_todo(iwork, part, iq, irr, auxdyn_loc)
           WRITE(stdout,'(/,2x,76("+"))')

           IF (lpwscf(iwork)) THEN
              WRITE(stdout,'(5x,"Doing geometry", i5)') iwork
           ELSEIF (lphonon(iwork)) THEN
              IF (all_geometries_together) THEN
                 igeom1=geometry(iwork)
              ELSE
                 igeom1=igeom
              ENDIF
              IF (fpol) THEN
                 WRITE(stdout,'(5x,"Doing frequency", i5)') iwork
              ELSE
                 IF (use_ph_images) THEN
                    WRITE(stdout,'(5x,"Doing all calculation&
                          & of geometry", i5 )') igeom1
                 ELSE
                    WRITE(stdout,'(5x,"Now I do point", &
                     i5," irrep", i5, " of geometry", i5 )') iq, irr, igeom1
                 ENDIF
              ENDIF
           ENDIF
           WRITE(stdout,'(2x,76("+"),/)')
            
           IF (lpwscf(iwork)) THEN
              CALL check_existence(iwork,part,igeom,run)
              IF (run) THEN
                 CALL do_pwscf(exit_status, .TRUE.)
                 energy_geo(iwork)=etot
                 IF (lstress(iwork)) THEN
                    sigma_geo(:,:,iwork)=sigma(:,:)
                 ENDIF
                 CALL save_existence(iwork,part,igeom)
              ENDIF
           END IF
           IF (lberry(iwork)) CALL do_berry(exit_status, &
                                            polar_geo(1,iwork), nppl)
           IF (lphonon(iwork)) THEN
              CALL do_phonon_tpw(auxdyn_loc)
              IF (all_geometries_together) CALL close_ph_geometry(.TRUE.)
           ENDIF
        END DO
     END IF
  END IF
  CALL mp_bcast(stop_signal_activated, meta_ionode_id, world_comm)
!
RETURN
END SUBROUTINE run_thermo_asynchronously
