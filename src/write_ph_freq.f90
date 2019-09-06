! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE write_ph_freq(igeom)
  !-----------------------------------------------------------------------
  !
  !  This routine allocates space for saving the frequences of the
  !  igeom geometry, computes and save them in this variable.
  !  The behaviour of the routine depends on
  !  the flag with_eigen. When with_eigen is .FALSE. before diagonalizing
  !  the dynamical matrices the routine checks if the frequencies are on
  !  file. In that case it reads the frequencies and exit. If with_eigen
  !  is .TRUE. the dynamical matrices are always diagonalized because the
  !  eigenvector are not saved on disk.
  !  The variables allocated by this routine:
  !  freq_save, z_save, dos_q, dos_wq they must be deallocated outside the
  !  routine after they have been used. 
  !  Presently they are deallocated by this routine or by deallocate_q2r
  !  at the end of the run.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE ions_base,     ONLY : nat, tau, ityp, amass
  USE cell_base,     ONLY : ibrav, at, bg, celldm
  USE phonon_save,   ONLY : freq_save, z_save
  USE thermo_mod,    ONLY : tot_ngeo
  USE constants,     ONLY : amu_ry
  USE control_dosq,  ONLY : nq1_d, nq2_d, nq3_d, dos_q, dos_wq, dos_nqs
  USE data_files,    ONLY : fldosfrq
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE control_thermo, ONLY : with_eigen
  USE ph_freq_module, ONLY : init_ph_freq, read_ph_freq_data, &
                             write_ph_freq_data, ph_freq_type, destroy_ph_freq
  USE mp,             ONLY : mp_sum, mp_bcast
  USE mp_world,       ONLY : mpime, world_comm, nproc
  USE io_global,      ONLY : meta_ionode, meta_ionode_id
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: igeom

  CHARACTER(LEN=256) :: filename

  INTEGER :: iq, na, nta, ipol, jmode, nq, nqx, iq_eff, ntetra, &
             startq, lastq, iproc, imode
  INTEGER, ALLOCATABLE :: nq_all(:)
  INTEGER, ALLOCATABLE :: start_proc(:)
  INTEGER, ALLOCATABLE :: last_proc(:)
  REAL(DP) :: masst, unorm(3*nat)
  REAL(DP), ALLOCATABLE :: q(:,:), wq(:)
  LOGICAL  :: file_exist
  TYPE(ph_freq_type) :: buffer
  CHARACTER(LEN=6) :: int_to_char
  LOGICAL :: check_file_exists
  !
  filename='phdisp_files/'//TRIM(fldosfrq)
  file_exist=check_file_exists(filename)
  !
  ! Allocate space for the q points and the frequencies
  !
  IF (ALLOCATED(dos_q)) DEALLOCATE(dos_q)
  IF (ALLOCATED(dos_wq)) DEALLOCATE(dos_wq)
  ntetra = 6 * nq1_d * nq2_d * nq3_d
  nqx = nq1_d * nq2_d * nq3_d

  IF (.NOT.with_eigen .AND. file_exist) THEN
!
! if eigenvectors are not needed and the frequencies are already on disk, 
! read them and exit
!
     WRITE(stdout,'(/,2x,76("-"))')
     WRITE(stdout,'(5x,"Frequencies for BZ integrations read from file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filename)
     WRITE(stdout,'(2x,76("-"),/)')

     IF (meta_ionode) CALL read_ph_freq_data(buffer,filename)
     CALL mp_bcast(buffer%nq, meta_ionode_id, world_comm)
     nq=buffer%nq
     IF (.NOT.meta_ionode) THEN
        ALLOCATE(buffer%nu(3*nat,nq))
        ALLOCATE(buffer%wg(nq))
     ENDIF
     CALL mp_bcast(buffer%wg, meta_ionode_id, world_comm)
     CALL mp_bcast(buffer%nu, meta_ionode_id, world_comm)
     CALL divide(world_comm, buffer%nq, startq, lastq)
     CALL init_ph_freq(ph_freq_save(igeom), nat, nq1_d, nq2_d, nq3_d, &
                                         startq, lastq, nq, with_eigen)

     ALLOCATE ( freq_save(3*nat, startq:lastq ))
     ALLOCATE(dos_wq(startq:lastq))
     freq_save(:,startq:lastq) = buffer%nu(:,startq:lastq)
     dos_wq(startq:lastq) = buffer%wg(startq:lastq)
     dos_nqs=nq
     iq_eff=0
     DO iq=startq,lastq
        iq_eff=iq_eff+1
        ph_freq_save(igeom)%nu(:,iq_eff)=buffer%nu(:,iq)
        ph_freq_save(igeom)%wg(iq_eff)=buffer%wg(iq)
     ENDDO
     CALL destroy_ph_freq(buffer) 
     RETURN
  ELSE
!
!   otherwise recompute the q points, the weights and the frequencies
!
     ALLOCATE(q(3,nqx))
     ALLOCATE(wq(nqx))
     CALL gen_qpoints_tpw (ibrav, at, bg, nat, tau, ityp, nq1_d, nq2_d, nq3_d, &
          ntetra, nqx, nq, q, wq)
     dos_nqs=nq
     CALL divide(world_comm, nq, startq, lastq)
     ALLOCATE ( dos_q(3, nq) )
     ALLOCATE ( dos_wq(startq:lastq) )
     ALLOCATE ( freq_save(3*nat, startq:lastq) )
     IF (with_eigen) ALLOCATE ( z_save(3*nat, 3*nat, startq:lastq) )
!
!  This avoid to keep allocated too much memory
!
     dos_q(1:3,1:nq)=q(:,1:nq)
     dos_wq(startq:lastq)=wq(startq:lastq)
     DEALLOCATE(q)
     DEALLOCATE(wq)
     CALL matdyn_interp(dos_nqs, dos_q, startq, lastq, with_eigen)
  END IF
!
!   initialize space to save the frequencies
!
  CALL init_ph_freq(ph_freq_save(igeom), nat, nq1_d, nq2_d, nq3_d, &
                                         startq, lastq, nq, with_eigen)
!
!   save the frequencies
!
  iq_eff=0
  DO iq=startq, lastq
     iq_eff=iq_eff+1
     ph_freq_save(igeom)%wg(iq_eff)=dos_wq(iq)
     ph_freq_save(igeom)%nu(:,iq_eff)=freq_save(:,iq)
  ENDDO

  IF (with_eigen) THEN
!
!  The eigenvectors are not saved on disk. They are recalculated each time.
!
     iq_eff=0
     DO iq=startq, lastq
        iq_eff=iq_eff+1
        unorm=0.0_DP
        DO na = 1,nat
           nta = ityp(na)
           masst=SQRT(amu_ry*amass(nta))
           DO ipol = 1,3
              jmode=(na-1)*3+ipol
              ph_freq_save(igeom)%displa(jmode,:,iq_eff)=&
                                              z_save(jmode,:,iq)*masst
              unorm(:)=unorm(:)+ABS(ph_freq_save(igeom)&
                                        %displa(jmode,:,iq_eff))**2
           END DO
        END DO
        DO imode=1,3*nat
           ph_freq_save(igeom)%displa(:,imode,iq_eff)=&
              ph_freq_save(igeom)%displa(:,imode,iq_eff)/SQRT(unorm(imode))
        ENDDO
        z_save(:,:,iq)= ph_freq_save(igeom)%displa(:,:,iq_eff)             
     END DO
  ELSE
     ALLOCATE(nq_all(nproc))
     ALLOCATE(start_proc(nproc))
     ALLOCATE(last_proc(nproc))
     nq_all=0
     nq_all(mpime+1) = ph_freq_save(igeom)%nq_eff
     CALL mp_sum(nq_all,world_comm)
     nq=0
     DO iproc=1, nproc
        start_proc(iproc)=nq+1
        nq=nq+nq_all(iproc)
        last_proc(iproc)=nq
     ENDDO
     
     CALL init_ph_freq(buffer, nat, nq1_d, nq2_d, nq3_d, 1, nq, &
                                                         nq, with_eigen) 
     
     buffer%nu=0.0_DP
     buffer%wg=0.0_DP
     buffer%nu(:,start_proc(mpime+1):last_proc(mpime+1))=&
                             ph_freq_save(igeom)%nu
     buffer%wg(start_proc(mpime+1):last_proc(mpime+1))=&
                             ph_freq_save(igeom)%wg

     CALL mp_sum(buffer%nu,world_comm)
     CALL mp_sum(buffer%wg,world_comm)

     DEALLOCATE(nq_all)
     DEALLOCATE(start_proc)
     DEALLOCATE(last_proc)
     
     WRITE(stdout,'(/,2x,76("+"))')
     WRITE(stdout,'(5x,"Writing frequencies for BZ integration on file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filename)
     WRITE(stdout,'(2x,76("+"),/)')
     IF (meta_ionode) CALL write_ph_freq_data(buffer,filename)
     CALL destroy_ph_freq(buffer)
  END IF
  !
  RETURN
END SUBROUTINE write_ph_freq
!
