! Copyright (C) 2016-2019 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE write_ph_freq(igeom)
  !-----------------------------------------------------------------------
  !
  !  This routine computes the frequencies needed to make
  !  the intergrals over the Brillouin zone for the geometry igeom
  !  and saves them in ph_freq_save.
  !  When with_eigen=.FALSE., before diagonalizing the dynamical matrices, 
  !  the routine checks if the frequencies are on file. In that case 
  !  it reads the frequencies and exit. If the frequencies are not on file
  !  the routine computes them and saves them on file. 
  !  When with_eigen=.TRUE. the dynamical matrices are always diagonalized 
  !  and the frequencies and eigenvectors are not saved on disk.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE ions_base,     ONLY : nat, tau, ityp
  USE cell_base,     ONLY : ibrav, at, bg

  USE control_dosq,  ONLY : nq1_d, nq2_d, nq3_d
  USE data_files,    ONLY : fldosfrq
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE control_thermo, ONLY : with_eigen
  USE ph_freq_module, ONLY : init_ph_freq
  USE matdyn_mod,     ONLY : matdyn_interp
  USE mp_world,       ONLY : world_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: igeom

  CHARACTER(LEN=256) :: filename
  INTEGER :: iq, nq, nqx, startq, lastq
  REAL(DP), ALLOCATABLE :: q(:,:), wq(:), dos_q(:,:)
  LOGICAL :: file_exist
  LOGICAL :: check_file_exists
  !
  filename='phdisp_files/'//TRIM(fldosfrq)
  file_exist=check_file_exists(filename)

  IF (.NOT.with_eigen .AND. file_exist) THEN
!
!  If the eigenvectors are not needed and the frequencies are already on disk, 
!  read them and exit
!
     WRITE(stdout,'(/,2x,76("-"))')
     WRITE(stdout,'(5x,"Frequencies for BZ integrations read from file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filename)
     WRITE(stdout,'(2x,76("-"),/)')

     CALL readd_ph_freq(ph_freq_save(igeom),filename)

     RETURN
  ELSE
!
!   Otherwise computes the q points, the weights and the frequencies.
!   First generates the q points
!
     nqx = nq1_d * nq2_d * nq3_d
     ALLOCATE(q(3,nqx))
     ALLOCATE(wq(nqx))
     CALL gen_qpoints_tpw (ibrav, at, bg, nat, tau, ityp, &
          nq1_d, nq2_d, nq3_d, nqx, nq, q, wq)
!
!  To avoid to keep allocated too much memory we copy the q points in 
!  dos_q and their weight in ph_freq_save (note that this quantity is
!  distributed among all processors). 
!
     ALLOCATE ( dos_q(3, nq) )
     dos_q(1:3,1:nq)=q(:,1:nq)

     CALL divide(world_comm, nq, startq, lastq)
     CALL init_ph_freq(ph_freq_save(igeom), nat, nq1_d, nq2_d, nq3_d, &
                                         startq, lastq, nq, with_eigen)
     ph_freq_save(igeom)%wg(:)=wq(startq:lastq)
     DEALLOCATE(q)
     DEALLOCATE(wq)
!
!   here computes the frequencies
!
     IF (with_eigen) THEN
        CALL matdyn_interp(nq, dos_q, ph_freq_save(igeom)%nu(:,:),  & 
                          startq, lastq, ph_freq_save(igeom)%displa(:,:,:))
     ELSE
        CALL matdyn_interp(nq, dos_q, ph_freq_save(igeom)%nu(:,:), startq, &
                                                                   lastq)
     ENDIF
!
!  The q vectors used to compute the vibrational density of states are not 
!  available outside this routine
!
     DEALLOCATE(dos_q)
  END IF

  IF (.NOT.with_eigen) THEN
!
!   If with_eigen=.TRUE. the eigenvalues and eigenvectors are not saved
!   on disk, but recalculated every time this routine is called.
!
     WRITE(stdout,'(/,2x,76("+"))')
     WRITE(stdout,'(5x,"Writing frequencies for BZ integration on file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filename)
     WRITE(stdout,'(2x,76("+"),/)')

     CALL cwrite_ph_freq(ph_freq_save(igeom),filename)
  END IF
  !
  RETURN
END SUBROUTINE write_ph_freq

!---------------------------------------------------------------------
SUBROUTINE readd_ph_freq(ph_freq_data,filename)
!---------------------------------------------------------------------
!
!   This routine reads and distributes the data of a ph_freq_type.
!   The data are assumed to be in a single file and are distributed
!   to all the available processors. 
!
USE ions_base,      ONLY : nat
USE control_thermo, ONLY : with_eigen
USE ph_freq_module, ONLY : read_ph_freq_data, init_ph_freq, destroy_ph_freq, &
                           ph_freq_type
USE control_dosq,   ONLY : nq1_d, nq2_d, nq3_d
USE io_global,      ONLY : meta_ionode_id, meta_ionode
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_bcast
IMPLICIT NONE

TYPE(ph_freq_type) :: ph_freq_data
CHARACTER(LEN=256) :: filename

TYPE(ph_freq_type) :: buffer
INTEGER :: nq, startq, lastq

IF (meta_ionode) CALL read_ph_freq_data(buffer,filename)

CALL mp_bcast(buffer%nq, meta_ionode_id, world_comm)
nq=buffer%nq
IF (.NOT.meta_ionode) THEN
   ALLOCATE(buffer%nu(3*nat,nq))
   ALLOCATE(buffer%wg(nq))
ENDIF
CALL mp_bcast(buffer%nu, meta_ionode_id, world_comm)
CALL mp_bcast(buffer%wg, meta_ionode_id, world_comm)
CALL divide(world_comm, buffer%nq, startq, lastq)

CALL init_ph_freq(ph_freq_data, nat, nq1_d, nq2_d, nq3_d, startq, lastq, &
                                                          nq, with_eigen)
ph_freq_data%nu(:,:)=buffer%nu(:,startq:lastq)
ph_freq_data%wg(:)=buffer%wg(startq:lastq)

IF (with_eigen) THEN
   IF (.NOT.meta_ionode) ALLOCATE(buffer%displa(3*nat,3*nat,nq))
   CALL mp_bcast(buffer%displa, meta_ionode_id, world_comm)
   ph_freq_data%displa(:,:,:)=buffer%displa(:,:,startq:lastq)
ENDIF

CALL destroy_ph_freq(buffer)

RETURN
END SUBROUTINE readd_ph_freq

!---------------------------------------------------------------------
SUBROUTINE cwrite_ph_freq(ph_freq_data, filename)
!---------------------------------------------------------------------
!
!   This routine collects the data of a ph_freq_type and writes them on file.
!   The data are written only by the meta_ionode. 
!   All processors of the world communicator are assumed to have a 
!   piece of data
!
USE kinds,          ONLY : DP
USE ions_base,      ONLY : nat
USE control_thermo, ONLY : with_eigen
USE ph_freq_module, ONLY : write_ph_freq_data, init_ph_freq, destroy_ph_freq, &
                           ph_freq_type
USE control_dosq,   ONLY : nq1_d, nq2_d, nq3_d
USE io_global,      ONLY : meta_ionode
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum
IMPLICIT NONE

TYPE(ph_freq_type) :: ph_freq_data
CHARACTER(LEN=256) :: filename

TYPE(ph_freq_type) :: buffer
INTEGER :: nq, startq, lastq
!
!  Initialize a buffer with the correct dinension
!
nq=ph_freq_data%nq
CALL init_ph_freq(buffer, nat, nq1_d, nq2_d, nq3_d, 1, nq, nq, with_eigen)
!
!  find which q belong to this processor
!
CALL divide(world_comm, nq, startq, lastq)
!
!  copy them into buffer
!
buffer%nu=0.0_DP
buffer%wg=0.0_DP
buffer%nu(:,startq:lastq)=ph_freq_data%nu
buffer%wg(startq:lastq)=ph_freq_data%wg
IF (with_eigen) buffer%displa(:,:,startq:lastq)=ph_freq_data%displa(:,:,:)
!
!  and collect buffer on all processors
!
CALL mp_sum(buffer%nu,world_comm)
CALL mp_sum(buffer%wg,world_comm)
!
!  meta_ionode saves the collected data on file
!
IF (meta_ionode) CALL write_ph_freq_data(buffer,filename)
!
CALL destroy_ph_freq(buffer)

RETURN
END SUBROUTINE cwrite_ph_freq
