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
  !  igeom geometry and save them in this variable.
  !  The behaviour of the routine depends on
  !  the flag with_eigen. When with_eigen is .FALSE. before diagonalizing
  !  the dynamical matrices the routine checks if the frequencies are on
  !  file. In that case it reads the frequencies and exit. If with_eigen
  !  is .TRUE. the dynamical matrices are always diagonalized because the
  !  eigenvector are not saved on disk.
  !  The variables allocated by this routine:
  !  freq_save, z_save, dos_q, dos_wq they must be deallocated outside the
  !  routine
  !  separately after they have been used. A call to clean_ifc_variables
  !  cleanup these variables, among other things.
  !
  USE kinds,         ONLY : DP
  USE mp_images,     ONLY : my_image_id, root_image
  USE io_global,     ONLY : ionode, stdout
  USE ions_base,     ONLY : nat, tau, ityp, amass
  USE cell_base,     ONLY : ibrav, at, bg, celldm
  USE phonon_save,   ONLY : freq_save, z_save
  USE thermo_mod,    ONLY : tot_ngeo
  USE constants,     ONLY : amu_ry
  USE control_dosq,  ONLY : nq1_d, nq2_d, nq3_d, dos_q, dos_wq, dos_nqs
  USE ifc,           ONLY : atm, zeu, m_loc
  USE data_files,    ONLY : fldosfrq
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE control_thermo, ONLY : with_eigen
  USE ph_freq_module, ONLY : init_ph_freq, read_ph_freq_data, &
                             write_ph_freq_data
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: igeom

  CHARACTER(LEN=256) :: filename

  INTEGER :: iq, na, nta, ipol, jmode, nq, nqx, ntetra
  REAL(DP) :: masst
  LOGICAL  :: file_exist

  CHARACTER(LEN=6) :: int_to_char
  LOGICAL :: check_file_exists
  !
  filename='phdisp_files/'//TRIM(fldosfrq)
  file_exist=check_file_exists(filename)
  IF (.NOT.ALLOCATED(ph_freq_save)) ALLOCATE(ph_freq_save(tot_ngeo))

  IF ( my_image_id /= root_image ) RETURN
  !
  ! Allocate space for the q points and the frequencies
  !
  IF (ALLOCATED(dos_q)) DEALLOCATE(dos_q)
  IF (ALLOCATED(dos_wq)) DEALLOCATE(dos_wq)
  ntetra = 6 * nq1_d * nq2_d * nq3_d
  nqx = nq1_d * nq2_d * nq3_d
  ALLOCATE ( dos_q(3,nqx) )
  ALLOCATE ( dos_wq(nqx) )
  ALLOCATE ( freq_save(3*nat, nqx) )
  IF (with_eigen) ALLOCATE ( z_save(3*nat, 3*nat, nqx) )

  IF (.NOT.with_eigen .AND. file_exist) THEN
!
! if eigenvectors are not needed and the frequencies are already on disk, 
! read them and exit
!
     WRITE(stdout,'(/,2x,76("-"))')
     WRITE(stdout,'(5x,"Frequencies for BZ integrations read from file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filename)
     WRITE(stdout,'(2x,76("-"),/)')

     CALL read_ph_freq_data(ph_freq_save(igeom),filename)
     nq=ph_freq_save(igeom)%nq
     freq_save(:,1:nq) = ph_freq_save(igeom)%nu(:,1:nq)
     dos_wq(1:nq) = ph_freq_save(igeom)%wg(1:nq)
     dos_nqs=nq
     RETURN
  ELSE
!
!   otherwise recompute the q points, the weights and the frequencies
!
     CALL gen_qpoints (ibrav, at, bg, nat, tau, ityp, nq1_d, nq2_d, nq3_d, &
          ntetra, nqx, nq, dos_q, dos_wq)
     dos_nqs=nq
     CALL matdyn_interp(dos_nqs, dos_q, with_eigen)
  END IF
!
!   initialize space to save the frequencies
!
  CALL init_ph_freq(ph_freq_save(igeom), nat, nq1_d, nq2_d, nq3_d, nq, &
                                                              with_eigen)
!
!   save the frequencies
!
  DO iq=1, nq
     ph_freq_save(igeom)%wg(iq)=dos_wq(iq)
     ph_freq_save(igeom)%nu(:,iq)=freq_save(:,iq)
  ENDDO

  IF (with_eigen) THEN
!
!  The eigenvectors are not saved on disk. They are recalculated each time.
!
     DO iq=1, nq
        DO na = 1,nat
           nta = ityp(na)
           masst=SQRT(amu_ry*amass(nta))
           DO ipol = 1,3
              jmode=(na-1)*3+ipol
              ph_freq_save(igeom)%displa(jmode,:,iq)=z_save(jmode,:,iq)*masst
           END DO
        END DO
     END DO
  ELSE
     WRITE(stdout,'(/,2x,76("+"))')
     WRITE(stdout,'(5x,"Writing frequencies for BZ integration on file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filename)
     WRITE(stdout,'(2x,76("+"),/)')
     CALL write_ph_freq_data(ph_freq_save(igeom),filename)
  END IF
  !
  RETURN
END SUBROUTINE write_ph_freq
!
!-----------------------------------------------------------------------
SUBROUTINE gen_qpoints (ibrav, at_, bg_, nat, tau, ityp, nk1, nk2, nk3, &
     ntetra, nqx, nq, q, wq)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : set_sym_bl, find_sym, s, irt, nsym, &
                         nrot, t_rev, time_reversal,  sname, &
                         allfrac, remove_sym
  USE initial_conf, ONLY : nr1_save, nr2_save, nr3_save
  !
  IMPLICIT NONE
  ! input
  INTEGER :: ibrav, nat, nk1, nk2, nk3, ntetra, ityp(*)
  REAL(DP) :: at_(3,3), bg_(3,3), tau(3,nat)
  ! output
  INTEGER :: nqx, nq, tetra(4,ntetra)
  REAL(DP) :: q(3,nqx), wq(nqx)
  ! local
  REAL(DP) :: xqq(3), mdum(3,nat)
  LOGICAL :: magnetic_sym=.FALSE., skip_equivalence=.FALSE.
  !
  time_reversal = .true.
  t_rev(:) = 0
  xqq (:) =0.d0
  at = at_
  bg = bg_
  CALL set_sym_bl ( )
  !
  CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, bg, nqx, &
                           0,0,0, nk1,nk2,nk3, nq, q, wq)
  !
  CALL find_sym ( nat, tau, ityp, .NOT.time_reversal, mdum )
  IF ( .NOT. allfrac ) CALL remove_sym ( nr1_save, nr2_save, nr3_save )
  !
  CALL irreducible_BZ (nrot, s, nsym, time_reversal, magnetic_sym, &
                       at, bg, nqx, nq, q, wq, t_rev)
  !
!  IF (ntetra /= 6 * nk1 * nk2 * nk3) &
!       CALL errore ('gen_qpoints','inconsistent ntetra',1)
  !
!  write(stdout,*) 'tetrahedra'
!  CALL tetrahedra (nsym, s, time_reversal, t_rev, at, bg, nqx, 0, 0, 0, &
!       nk1, nk2, nk3, nq, q, ntetra, tetra)
  !
  RETURN
END SUBROUTINE gen_qpoints
!
