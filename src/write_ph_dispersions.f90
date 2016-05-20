! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE write_ph_dispersions()
  !-----------------------------------------------------------------------
  !
  !  This routine computes the phonon frequencies along the path
  !  of q points defined by thermo_pw and writes the result on file.
  !  The symmetry of the modes is analyzed and the frequencies and
  !  their symmetry are written in two files with a format
  !  that can be interpreted by the plotband_sub routine.
  !  The routine uses the matdyn_interp routine to actually interpolate
  !  the dynamical matrices on the q point path.
  !  If the file with the frequencies is already found on disk the
  !  routine exits.
  !
  USE kinds,      ONLY : DP
  USE mp_images,  ONLY : my_image_id, root_image
  USE io_global,  ONLY : ionode, stdout
  USE constants,  ONLY : ry_to_cmm1
  USE ions_base,  ONLY : nat, tau, ityp, nsp, amass
  USE symm_base,  ONLY : set_sym
  USE phonon_save, ONLY : freq_save, z_save
  USE thermo_sym, ONLY : code_group_save
  USE rap_point_group,  ONLY : code_group
  USE control_pwrun, ONLY : nr1_save, nr2_save, nr3_save
  USE control_paths, ONLY : disp_q, disp_wq, disp_nqs, high_sym_path
  USE control_ph,    ONLY : xmldyn
  USE ifc,           ONLY : m_loc, atm, zeu, has_zstar
  USE noncollin_module, ONLY : nspin_mag
  USE data_files,     ONLY : flfrq, flvec
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: filefrq, filename, filevec
  INTEGER :: nqs, nta, ipol, ios, code_group_old, n, i, iq, nq, iout
  LOGICAL :: lo_to_split
  CHARACTER(LEN=15), ALLOCATABLE :: name_rap_mode(:)
  REAL(DP) :: ps, qh, dq(3), q1(3), q2(3), modq1, modq2, dqmod, dqmod_save
  REAL(DP), ALLOCATABLE :: w2(:,:)
  INTEGER, ALLOCATABLE :: num_rap_mode(:,:), qcode_group(:), aux_ind(:)
  LOGICAL, ALLOCATABLE :: high_sym(:)
  LOGICAL :: check_file_exists
  !
  filefrq="phdisp_files/"//TRIM(flfrq)
  IF (check_file_exists(filefrq)) THEN
     WRITE(stdout,'(/,2x,76("-"))')
     WRITE(stdout,'(5x,"Frequencies for dispersions are on file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filefrq)
     WRITE(stdout,'(2x,76("-"),/)')
     RETURN
  ENDIF

  IF ( my_image_id /= root_image ) RETURN

  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x,"Interpolating the dynamical matrices")')
  WRITE(stdout,'(5x,"Frequencies written on file ")')
  WRITE(stdout,'(5x,a)') TRIM(filefrq)
  WRITE(stdout,'(2x,76("+"),/)')

  ALLOCATE(freq_save(3*nat, disp_nqs))
  ALLOCATE(z_save(3*nat, 3*nat, disp_nqs))
  ALLOCATE(w2(3*nat, disp_nqs))
!
!  we always need the eigenvectors to make the symmetry analysis
!
  CALL matdyn_interp(.TRUE.)

  iout=0
  IF (flvec/=' ') THEN
     iout=4
     filevec="phdisp_files/"//flvec
     IF (ionode) OPEN(UNIT=iout, FILE=TRIM(filevec), STATUS='unknown', &
                                                       FORM='formatted')
  END IF
  DO iq=1, disp_nqs
     w2(:,iq) = SIGN((freq_save(:,iq)/ry_to_cmm1)**2, freq_save(:,iq) )
     IF (iout/=0.AND.ionode) CALL writemodes(nat,disp_q(1,iq),w2(1,iq), &
                                                       z_save(1,1,iq),iout)
  END DO
  IF (iout/=0.AND.ionode) CLOSE(unit=iout)
  
  nq=disp_nqs
  ALLOCATE ( num_rap_mode(3*nat,nq) )
  ALLOCATE ( high_sym(nq) )
  ALLOCATE ( qcode_group(nq) )
  ALLOCATE ( name_rap_mode(3*nat) )

  IF (xmldyn) CALL set_sym(nat, tau, ityp, nspin_mag, m_loc, nr1_save, &
                                                      nr2_save, nr3_save )
  num_rap_mode=-1
!
!  Initialize high_sym
!
  high_sym(1:nq)=high_sym_path(1:nq)
!
! Now at all q points for which we can perform the symmetry analysis check
! the bands. If there some symmetry change is detected change also 
! high_symmetry
!
  DO n=1, nq
     lo_to_split=.FALSE.
     IF (n==1) THEN
        high_sym(n)=.TRUE.
        IF (nq>1) THEN
           dq(:) = disp_q(:,2) - disp_q(:,1)
        ELSE
           dq(:) = 0.0_DP
        END IF
        dqmod_save = sqrt( dq(1)**2 + dq(2)**2 + dq(3)**2 )
        code_group_old=0
     ENDIF
  !
  ! Cannot use the small group of \Gamma to analize the symmetry
  ! of the mode if there is an electric field.
  !
     qh = SQRT(disp_q(1,n)**2+disp_q(2,n)**2+disp_q(3,n)**2)
     IF (qh < 1.d-9 .AND. has_zstar) lo_to_split=.TRUE.
     IF (xmldyn.AND..NOT.lo_to_split) THEN
        CALL find_representations_mode_q(nat,nsp,disp_q(:,n), &
                    w2(:,n),z_save(:,:,n),tau,ityp,amass,name_rap_mode, &
                    num_rap_mode(:,n), nspin_mag)
        qcode_group(n)=code_group
        IF (n==1) THEN
           code_group_old=code_group
        ELSE
           dq(:) = disp_q(:,n) - disp_q(:,n-1)
           dqmod= sqrt( dq(1)**2 + dq(2)**2 + dq(3)**2 )
           IF (dqmod < 1.D-6) THEN
              !
              !   In this case is_high_sym does not change because the point
              !   is the same
              high_sym(n)=high_sym(n-1)
              !
           ELSE IF (dqmod < 5.0_DP * dqmod_save) THEN
!
!    In this case the two points are considered close
!
              IF (.NOT. high_sym(n-1)) &
                 high_sym(n) = code_group /= code_group_old .OR. high_sym(n)

              dqmod_save= MAX(dqmod_save * 0.5_DP, dqmod)
           ELSE
              high_sym(n)=.TRUE.
           ENDIF
           code_group_old=code_group
        ENDIF
!     write(6,'(2i5, 3f15.5,l5)') n, qcode_group(n), q(:,n), high_sym(n)
     ELSEIF (lo_to_split) THEN
!
!  At gamma the group is the point group of the solid
!
        qcode_group(n)=code_group_save
     ENDIF
  END DO
  !
  IF (flfrq.NE.' '.AND.ionode) THEN
     OPEN (unit=2,file=filefrq ,status='unknown',form='formatted')
     WRITE(2, '(" &plot nbnd=",i6,", nks=",i6," /")') 3*nat, nq
     DO n=1, nq
        WRITE(2, '(10x,3f10.6)')  disp_q(1,n), disp_q(2,n), disp_q(3,n)
        WRITE(2,'(6f10.4)') (freq_save(i,n), i=1,3*nat)
     END DO
     CLOSE(unit=2)
  END IF
  !
  !  If the force constants are in the xml format we write also
  !  the file with the representations of each mode
  !
  IF (flfrq.NE.' '.AND.xmldyn) THEN
     ALLOCATE(aux_ind(nq))
     aux_ind=0
     CALL find_aux_ind_xk(disp_q(1,1), disp_q(1,2), aux_ind(2))
     CALL find_aux_ind_xk(disp_q(1,nq), disp_q(1,nq-1), aux_ind(nq-1))
     DO n=2,nq-1
!        write(6,'(3f15.5,3l5)') q(:,n), high_sym(n-1), high_sym(n), high_sym(n+1)
        IF (high_sym(n).AND..NOT.high_sym(n+1)) &
           CALL find_aux_ind_xk(disp_q(1,n), disp_q(1,n+1), aux_ind(n+1))
        IF (high_sym(n).AND..NOT.high_sym(n-1)) &
           CALL find_aux_ind_xk(disp_q(1,n), disp_q(1,n-1), aux_ind(n-1))
     ENDDO

     IF (ionode) THEN
        filename=TRIM(filefrq)//'.rap'
        OPEN (UNIT=2,FILE=TRIM(filename),STATUS='unknown',FORM='formatted')
        WRITE(2, '(" &plot_rap nbnd_rap=",i6,", nks_rap=",i6," /")') 3*nat, nq
        DO n=1, nq
           WRITE(2,'(10x,3f10.6,l6,2i6)') disp_q(1,n), disp_q(2,n), &
                                          disp_q(3,n), high_sym(n),&
                                          qcode_group(n), aux_ind(n) 
           WRITE(2,'(6i10)') (num_rap_mode(i,n), i=1,3*nat)
        END DO
        CLOSE(unit=2)
     ENDIF

     DEALLOCATE(aux_ind)
  ENDIF
  !
  DEALLOCATE (w2) 
  DEALLOCATE (freq_save) 
  DEALLOCATE (z_save) 
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (name_rap_mode)
  DEALLOCATE (qcode_group)
  DEALLOCATE (high_sym)
  disp_nqs=0
  !
  RETURN
END SUBROUTINE write_ph_dispersions
!
SUBROUTINE find_representations_mode_q ( nat, ntyp, xq, w2, u, tau, ityp, &
                  amass, name_rap_mode, num_rap_mode, nspin_mag )

  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : find_sym, s, sr, ftau, irt, nsym, &
                         nrot, t_rev, time_reversal, sname, copy_sym, &
                         s_axis_to_cart
  USE rap_point_group,  ONLY : code_group, gname

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, ntyp, nspin_mag
  REAL(DP), INTENT(IN) :: xq(3), amass(ntyp), tau(3,nat)
  REAL(DP), INTENT(IN) :: w2(3*nat)
  INTEGER, INTENT(IN) :: ityp(nat)
  COMPLEX(DP), INTENT(IN) :: u(3*nat,3*nat)
  CHARACTER(15), INTENT(OUT) :: name_rap_mode(3*nat)
  INTEGER, INTENT(OUT) :: num_rap_mode(3*nat)
  REAL(DP) :: gi (3, 48), gimq (3), sr_is(3,3,48), rtau(3,48,nat)
  INTEGER :: irotmq, nsymq, nsym_is, isym, i, ierr
  LOGICAL :: minus_q, search_sym, sym(48), magnetic_sym
!
!  find the small group of q
!
  time_reversal=.TRUE.
  IF (.NOT.time_reversal) minus_q=.FALSE.

  sym(1:nsym)=.true.
  call smallg_q (xq, 0, at, bg, nsym, s, ftau, sym, minus_q)
  nsymq=copy_sym(nsym,sym )
  call s_axis_to_cart ()
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
!
!  if the small group of q is non symmorphic,
!  search the symmetries only if there are no G such that Sq -> q+G
!
  search_sym=.TRUE.
  IF ( ANY ( ftau(:,1:nsymq) /= 0 ) ) THEN
     DO isym=1,nsymq
        search_sym=( search_sym.and.(abs(gi(1,isym))<1.d-8).and.  &
                                    (abs(gi(2,isym))<1.d-8).and.  &
                                    (abs(gi(3,isym))<1.d-8) )
     END DO
  END IF
!
!  Set the representations tables of the small group of q and
!  find the mode symmetry
!
  IF (search_sym) THEN
     magnetic_sym=(nspin_mag==4)
     CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)
     sym (1:nsym) = .TRUE.
     CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
     CALL find_mode_sym_new (u, w2, tau, nat, nsymq, sr, irt, xq,    &
             rtau, amass, ntyp, ityp, 1, .FALSE., .FALSE., num_rap_mode, ierr)

     CALL print_mode_sym(w2, num_rap_mode, .FALSE.)

  ELSE
     CALL find_group(nsymq,sr,gname,code_group)
  ENDIF
  RETURN
END SUBROUTINE find_representations_mode_q
