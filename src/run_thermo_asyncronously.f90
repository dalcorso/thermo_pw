!
! Copyright (C) 2013-2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------
SUBROUTINE run_thermo_asynchronously(nwork, part, igeom)
!-------------------------------------------------------------
!
!  This routine is the main driver of the asynchronous calculation.
!  It uses the routines of mp_asyn to coordinate the asynchronous
!  run of the nwork tasks. For each task it calls the routine 
!  set_thermo_work_todo that sets the input and run_iwork that
!  call pw.x or ph.x.
!
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
  USE control_thermo,  ONLY : max_seconds_tpw
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nwork, part, igeom

  INTEGER :: iq, irr
  INTEGER, ALLOCATABLE :: proc_num(:)
  INTEGER :: proc_per_image, iwork, image
  LOGICAL :: all_done_asyn
  REAL(DP) :: get_clock
  !
  IF ( nwork == 0 ) RETURN
  !
  !  Two calculations are now possible. If with_asyn_images is .TRUE. each
  !  image does a different task, and the first image is the master that 
  !  keeps into account how much work there is to do and who does what.
  !  Otherwise the works are done one after the other in the traditional
  !  fashion
  !
  IF ( with_asyn_images .AND. nwork > 1) THEN
     IF (my_image_id == root_image) THEN
        !
        !  This is the master and it has to rule all the other images
        !  It has to call asyn_master_init and give the number of workers
        !  the number of tasks and the number of one processor
        !  per image, proc_num, which is the ionode of each image
        !  used to communicate among images
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
              !  Now also the master can do something. Ask for some work.
              !
              IF (ionode) CALL asyn_master_work(iwork)
              CALL mp_bcast(iwork, ionode_id, intra_image_comm) 
              !
              !   And do the work
              !
              IF (iwork>0) THEN
                 CALL set_thermo_work_todo(iwork, part, iq, irr)
                 WRITE(stdout,'(/,2x,76("+"))')
                 WRITE(stdout,'(5x,"I am the master doing work",&
                                       &i5," / ", i5)') iwork, nwork
                 CALL run_iwork(iwork, part, iq, irr, igeom)
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
        !  First initializes the worker stuff. It declares the identity 
        !  of the master and which communicator to use
        !
        IF (ionode) CALL asyn_worker_init(meta_ionode_id, world_comm)
        iwork=1
        DO WHILE (iwork > 0)
           !
           !  The root_processor of each image asks the master for some work
           !  to do and sends the info to all processors of its image
           !  This is a blocking request, all the image blocks here until the
           !  master answers.
           !
           IF (ionode) CALL asyn_worker(iwork)
           CALL mp_bcast(iwork, root_image, intra_image_comm) 
           !
           !  and then do the work
           !
           IF (iwork>0) THEN
              CALL set_thermo_work_todo(iwork, part, iq, irr)
              WRITE(stdout,'(/,2x,76("+"))')
              WRITE(stdout,'(5x,"I am image ",i5," doing work",&
                                   &i5," / ", i5)') my_image_id, iwork, nwork
              CALL run_iwork(iwork, part, iq, irr, igeom)
           END IF
        END DO
     END IF
  ELSE
     !
     !  This is the case without images. There is only one image,
     !  the master, that does all the works one after the other.
     !  This part is called also if there is only one work. In this
     !  case the master does it.
     !
     IF (my_image_id == root_image) THEN
        DO iwork = 1, nwork
           CALL set_thermo_work_todo(iwork, part, iq, irr)
           WRITE(stdout,'(/,2x,76("+"))')
           WRITE(stdout,'(5x,"Doing work",i5," / ", i5)') iwork, nwork
           CALL run_iwork(iwork, part, iq, irr, igeom) 
        END DO
     END IF
  END IF
  CALL mp_bcast(stop_signal_activated, meta_ionode_id, world_comm)
  CALL deallocate_asyn()
  !
RETURN
END SUBROUTINE run_thermo_asynchronously
