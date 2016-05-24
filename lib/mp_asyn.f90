!
! Copyright (C) 2013 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE mp_asyn
!
!  This module gives to the calling program the functionality of  
!  master and workers. The master initializes the variables of the module 
!  saying how many processors there are (nproc),
!  how much work (nwork) has to be done, which mpi communicator is used 
!  between master and workers, and an array proc_num(0:nworkers) that 
!  for the master (0) and for all the workers (1:nworkers) gives the 
!  number of the processor within the communicator. 
!  The typical use of this module is inter images asynchronous
!  communication. The images can divide the total work into nwork
!  pieces and carry it out asynchronously. The images should not 
!  communicate with each other: only intra-image communications are allowed
!  because in general images run asynchronously. They can however 
!  syncronize the work using this module, that tells to the representative 
!  of each image which is the task that the image has to carry out. 
!
!  The master calls once asyn_init, and as often as possible asyn_master.
!  The master calls also asyn_close when the work is finished.
!  The workers call asyn_worker when they need some work to do and
!  wait until the master answers. The master can call asyn_master_work 
!  to ask for same work to do. 
! 
!  If the code is running after a serial compilation the routines of this
!  module should not be called.
!
IMPLICIT NONE
#ifdef __MPI
INCLUDE 'mpif.h'
#endif
SAVE 
PRIVATE

INTEGER, PARAMETER :: READY=0   ! ready message
INTEGER, PARAMETER :: NO_WORK=-1! end of work message

INTEGER :: nwork                ! number of jobs to do
INTEGER :: nworkers             ! number of workers
INTEGER :: master               ! the master processor
INTEGER :: asyn_comm            ! the communicator for the asynchronous work
INTEGER :: tag=1                ! tag used in each message

INTEGER, ALLOCATABLE :: proc_num(:) ! gives the processor number of 
                                ! the master and of the workers inside the
                                ! communicator
INTEGER, ALLOCATABLE :: req(:)  ! request number in irecv 
INTEGER, ALLOCATABLE :: buf(:)  ! buffer to contain the messages
INTEGER, ALLOCATABLE :: sent(:) ! the master use this array to save the
                                ! processor that is doing a job
INTEGER, ALLOCATABLE :: doing(:) ! the master use this array to know
                                 ! what a worker is doing
LOGICAL, ALLOCATABLE :: done(:) ! master use this array to known which
                                ! jobs have been already done
LOGICAL, ALLOCATABLE :: done_proc(:) ! master use this array to know
                   ! to which processor it sent the NO_WORK message
LOGICAL :: with_asyn_images=.FALSE. ! The calling program must set this
                                ! variable to true to use this module

PUBLIC asyn_master_init, asyn_worker_init, asyn_close, asyn_master, &
       asyn_worker, asyn_master_work, with_asyn_images

CONTAINS

SUBROUTINE initialize_master(nproc_, nwork_, proc_num_, comm)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nproc_, nwork_, comm
INTEGER, INTENT(IN) :: proc_num_(0:nproc_-1)

nworkers=nproc_-1
master=proc_num_(0)
nwork=nwork_
asyn_comm=comm
ALLOCATE(done(nwork))
ALLOCATE(sent(nwork))
ALLOCATE(doing(0:nworkers))
done=.FALSE.
sent=-1
doing=NO_WORK
tag=1
!
! In the case there is only an image, only the master does some work
! and no communication occurs
!
IF (nworkers==0) RETURN

ALLOCATE(proc_num(nworkers))
ALLOCATE(done_proc(nworkers))
ALLOCATE(buf(nworkers))
ALLOCATE(req(nworkers))
done_proc=.FALSE.
proc_num=proc_num_(1:nworkers)

RETURN
END SUBROUTINE initialize_master

SUBROUTINE asyn_master_init(nproc_, nwork_, proc_num_, comm)
!
!  This routine initializes the asynchronous work. It is called by the master
!  when it knows:
!  nwork_ : number of works to do_. 
!  nproc_ : is the number of processors master+workers. 
!  proc_num_(0:nproc_-1) : mapping between master+workers and physical
!                          processors. proc_num_ is the number of the
!                          master or of the workers within the communicator
!  comm                  : the communicator
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nproc_, nwork_, comm
INTEGER, INTENT(IN) :: proc_num_(0:nproc_-1)
INTEGER :: iproc                ! counter on workers
INTEGER :: ierr                 ! error variable


CALL initialize_master(nproc_, nwork_, proc_num_, comm)

IF (nworkers==0) RETURN
!
! during initialization the master listens to all the workers for the
! READY message, without blocking
!
#ifdef __MPI
DO iproc=1,nworkers
   CALL mpi_irecv(buf(iproc),1,MPI_INTEGER,proc_num(iproc),tag,asyn_comm,&
                  req(iproc),ierr)
ENDDO
#endif

RETURN
END SUBROUTINE asyn_master_init

SUBROUTINE asyn_worker_init(master_, comm)
IMPLICIT NONE
INTEGER, INTENT(IN) :: master_ ! the master node
INTEGER, INTENT(IN) :: comm ! the comunicator

master=master_
asyn_comm=comm
tag=1

RETURN
END SUBROUTINE asyn_worker_init

SUBROUTINE asyn_master(all_done) 
!
!  This routine is called by the master. It checks if some worker is
!  available to do work and possibly send the work to do to the worker
!
IMPLICIT NONE
LOGICAL, INTENT(OUT) :: all_done    ! when this variable becomes true
                                    ! the master must stop calling the
                                    ! routine because all workers received
                                    ! the NO_WORK message.
LOGICAL :: work_finished            ! if this becomes .true. the work is 
                                    ! finished
#ifdef __MPI
INTEGER :: status_(MPI_STATUS_SIZE) ! status of the probe function
#endif
INTEGER :: iproc, iwork             ! counters
INTEGER :: ierr                     ! error variable
LOGICAL :: exst                     ! if true the worker is willing to work

#ifdef __MPI
DO iproc=1,nworkers
!
!  If this processor has already received the end of work message 
!  check the next
!
   IF (done_proc(iproc)) CYCLE
!
!   check if a message arrived from proc_num(iproc), without blocking
!
   CALL mpi_test(req(iproc),exst,status_,ierr)
!
!   if no message arrived from proc_num(iproc) check the next worker
!
   IF (.NOT.exst) CYCLE
!
!  here the processor iproc has answered, mark his work as done
!
   IF (doing(iproc)>0) done(doing(iproc))=.TRUE.
!
!   if proc_num(iproc) is ready to work, send it the work to do if
!   there is some work available, otherwise keep the worker idle waiting
!   for the work to do
!
   IF (choose_next(iwork, work_finished)) THEN
!
!   here the master blocks, it must be sure that proc_num(iproc) received the
!   message. At this point the work is considered as sent and the previous
!   job sent to iproc as done
!
      CALL mpi_send(iwork,1,MPI_INTEGER,proc_num(iproc),tag,asyn_comm, ierr)

      sent(iwork)=iproc
      doing(iproc)=iwork
!
!   listen again from processor proc_num(iproc), without blocking
!
      CALL mpi_irecv(buf(iproc),1,MPI_INTEGER,proc_num(iproc),tag, &
                         asyn_comm,req(iproc),ierr)
   ENDIF
!
!    if there is no more work to do tell iproc to exit
!
   IF (work_finished) THEN
      CALL mpi_send(NO_WORK,1,MPI_INTEGER,proc_num(iproc),tag, &
                                              asyn_comm,ierr)
      doing(iproc)=NO_WORK
      done_proc(iproc)=.TRUE.
   ENDIF
ENDDO
!
!  check if all processors have received the NO_WORK message
!
all_done=.TRUE.
DO iproc=1,nworkers
   all_done=all_done.AND.done_proc(iproc)
ENDDO
#endif

RETURN
END SUBROUTINE asyn_master

SUBROUTINE asyn_master_work(master_work) 
!
!  This subroutine is called by the master if it is ready to do some
!  work. In this case the routine gives it some job or -1 if there
!  is no more work to do 
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: master_work
INTEGER :: iwork
LOGICAL :: all_done, there_is_work, work_finished
!
!  When this routine is called the master has done his work
!
IF (doing(0)>0) done(doing(0))=.TRUE.
!
!  Check if there is some work available
!
there_is_work=choose_next(iwork, work_finished)
!
! If there is no work to do and the work is not finished the master 
! must loop here checking if the situation changes.
!
DO WHILE (.NOT.(there_is_work .OR. work_finished)) 
   CALL asyn_master(all_done)
   there_is_work=choose_next(iwork, work_finished)
END DO
!
! There is some work to do or the work is finished
!
IF (there_is_work) THEN
   sent(iwork)=0
   doing(0)=iwork
   master_work=iwork  
ELSE IF (work_finished) THEN
   doing(0)=NO_WORK
   master_work=-1
ENDIF

RETURN
END SUBROUTINE asyn_master_work

SUBROUTINE asyn_worker(worker_work) 
!
!  This subroutine is called by the worker if it is ready to do some
!  work. In this case the routine gives it some work or -1 if there
!  is no more work to do 
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: worker_work ! worker work
#ifdef __MPI
INTEGER :: status_(MPI_STATUS_SIZE) ! status of the receive function
INTEGER :: ierr                     ! error variable
INTEGER :: iwork                    ! auxiliary
!
!  the worker first tells the master that it is ready to work
!  and blocks waiting for the master to listen
!
CALL mpi_send(READY,1,MPI_INTEGER,master,tag,asyn_comm,ierr)
!
!  when the code arrives here the master has recognized that 
!  this worker wants something to do and send the job to do
!
CALL mpi_recv(iwork,1,MPI_INTEGER,master,tag,asyn_comm,status_,ierr)
worker_work=iwork
#endif
RETURN
END SUBROUTINE asyn_worker

LOGICAL FUNCTION choose_next(iwork, work_finished)
!
!  This function chooses the next work to do.
!  There are three possible outputs:
!  a) All works have been done. In this case work_finished becomes .TRUE.
!  b) There is some work to do. In this case the function returns .TRUE.
!     and iwork gives the work to do.
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: iwork
LOGICAL, INTENT(OUT) :: work_finished
INTEGER :: iw, ip
!
!  case a, all works have been done
!
choose_next = .FALSE.
iwork=0
work_finished=.TRUE.
DO iw=1, nwork
   work_finished = work_finished .AND. (sent(iw) /= -1)
ENDDO
IF (work_finished) RETURN
!
!  First check that for a given work there is no pending dependence
!
DO iw=1,nwork
   IF ( sent(iw) == -1 .AND. iwork==0 ) iwork=iw
END DO
choose_next = (iwork>0)

RETURN
END FUNCTION choose_next

SUBROUTINE asyn_close()
!
! deallocate everything and close
!
IMPLICIT NONE

IF (nworkers > 0) THEN
   DEALLOCATE(req)
   DEALLOCATE(buf)
   DEALLOCATE(done_proc)
   DEALLOCATE(proc_num)
END IF
DEALLOCATE(sent)
DEALLOCATE(done)
DEALLOCATE(doing)

RETURN
END SUBROUTINE asyn_close

END MODULE mp_asyn
