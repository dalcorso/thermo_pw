!
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
  USE constants,  ONLY : pi
  USE mp_images,  ONLY : my_image_id, root_image
  USE io_global,  ONLY : ionode, stdout
  USE proj_rap_point_group, ONLY : code_groupq_ext, qptype, lqproj, qgauge
  USE lattices,   ONLY : same_star
  USE point_group,ONLY : nsym_group
  USE constants,  ONLY : ry_to_cmm1
  USE ions_base,  ONLY : nat, tau, ityp, nsp, amass
  USE symm_base,  ONLY : set_sym, nsym, s, allfrac, remove_sym
  USE fft_base,   ONLY : dfftp
  USE cell_base,  ONLY : at
  USE thermo_sym, ONLY : code_group_save
  USE rap_point_group,  ONLY : code_group
  USE initial_conf, ONLY : nr1_save, nr2_save, nr3_save
  USE control_paths, ONLY : disp_q, disp_nqs, high_sym_path, nrap_plot, &
                            rap_plot, dkmod_save
  USE control_ph,    ONLY : xmldyn, search_sym
  USE control_lr,    ONLY : lgamma
  USE matdyn_mod,    ONLY : matdyn_interp
  USE lr_symm_base,  ONLY : rtau
  USE ifc,           ONLY : m_loc, has_zstar
  USE io_bands,      ONLY : write_bands, write_representations
  USE noncollin_module, ONLY : nspin_mag
  USE data_files,     ONLY : flfrq, flvec
  USE ph_symmetry,    ONLY : initialize_gcode_old
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: filefrq, filename, filevec
  INTEGER :: nqs, nta, ipol, ios, code_group_old, n, i, iq, nq, iout
  LOGICAL :: lo_to_split
  CHARACTER(LEN=15), ALLOCATABLE :: name_rap_mode(:)
  REAL(DP) :: ps, qh, dq(3), q1(3), q2(3), modq1, modq2, dqmod
  REAL(DP), ALLOCATABLE :: w2(:,:)
  INTEGER, ALLOCATABLE :: num_rap_mode(:,:), qcode_group(:), aux_ind(:), &
                          qcode_group_ext(:), ptypeq(:,:), lprojq(:)
  REAL(DP), ALLOCATABLE :: gaugeq(:,:), freq_save(:,:)
  COMPLEX(DP), ALLOCATABLE :: z_save(:,:,:)

  INTEGER :: qcode_old, isym
  LOGICAL, ALLOCATABLE :: high_sym(:), same_next(:)
  LOGICAL :: check_file_exists
  !
!
!  Set the BZ path for the present geometry
!
  CALL set_paths_disp()

  filefrq="phdisp_files/"//TRIM(flfrq)
  IF (check_file_exists(filefrq)) THEN
     WRITE(stdout,'(/,2x,76("-"))')
     WRITE(stdout,'(5x,"Frequencies for dispersions are on file ")') 
     WRITE(stdout,'(5x,a)') TRIM(filefrq)
     WRITE(stdout,'(2x,76("-"),/)')
     RETURN
  ENDIF

  IF (disp_nqs==0) RETURN

  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x,"Interpolating the dynamical matrices for dispersion")')
  WRITE(stdout,'(5x,"Frequencies written on file ")')
  WRITE(stdout,'(5x,a)') TRIM(filefrq)
  WRITE(stdout,'(2x,76("+"),/)')

  ALLOCATE(freq_save(3*nat, disp_nqs))
  ALLOCATE(z_save(3*nat, 3*nat, disp_nqs))
  ALLOCATE(w2(3*nat, disp_nqs))
!
!  we always need the eigenvectors to make the symmetry analysis
!
  CALL matdyn_interp(disp_nqs, disp_q, freq_save, 1, disp_nqs, z_save)

  IF ( my_image_id /= root_image ) THEN
     DEALLOCATE (freq_save) 
     DEALLOCATE (z_save) 
     DEALLOCATE (w2) 
     RETURN
  ENDIF

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
  ALLOCATE ( rtau(3,48,nat) )
  ALLOCATE ( num_rap_mode(3*nat,nq) )
  ALLOCATE ( high_sym(nq) )
  ALLOCATE ( qcode_group(nq) )
  ALLOCATE ( name_rap_mode(3*nat) )
  ALLOCATE ( qcode_group_ext(nq) )
  ALLOCATE ( lprojq(nq) )
  ALLOCATE ( ptypeq(3,nq) )
  ALLOCATE ( gaugeq(48,nq) )

  IF (xmldyn) THEN
     CALL set_sym(nat, tau, ityp, nspin_mag, m_loc)
     IF ( .NOT. allfrac ) CALL remove_sym ( dfftp%nr1, dfftp%nr2, dfftp%nr3 )
  ENDIF

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
  IF (nq > 0.AND.search_sym) WRITE(stdout,'(/,5x,70("*"))')
  qcode_old=0
  DO n=1, nq
     lo_to_split=.FALSE.
     IF (n==1) THEN
        high_sym(n)=.TRUE.
        IF (nq>1) THEN
           dq(:) = disp_q(:,2) - disp_q(:,1)
        ELSE
           dq(:) = 0.0_DP
        END IF
        code_group_old=0
        CALL initialize_gcode_old(0)
     ENDIF
  !
  ! Cannot use the small group of \Gamma to analize the symmetry
  ! of the mode if there is an electric field.
  !
     qh = SQRT(disp_q(1,n)**2+disp_q(2,n)**2+disp_q(3,n)**2)
     IF (qh < 1.d-9) THEN
        lgamma=.TRUE.
        IF (has_zstar) lo_to_split=.TRUE.
     ELSE
        lgamma=.FALSE.
     ENDIF
     IF (search_sym) WRITE(stdout, '(/,20x,"q=(",2(f10.5,","),f10.5,"  )")') &
                                                                 disp_q(:,n)
     IF (xmldyn.AND..NOT.lo_to_split) THEN
        IF (n>1) qcode_old=qcode_group(n-1)
        CALL find_representations_mode_q(nat,nsp,disp_q(:,n), &
                    w2(:,n),z_save(:,:,n),tau,ityp,amass, &
                    num_rap_mode(:,n), nspin_mag, qcode_old)
        qcode_group(n)=code_group
        qcode_group_ext(n)=code_groupq_ext
        ptypeq(:,n)=qptype(:)
        lprojq(n)=lqproj
        gaugeq(:,n)=qgauge(:)
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
           ELSE IF (dqmod < 5.0_DP * dkmod_save) THEN
!
!    In this case the two points are considered close
!
              IF (.NOT. high_sym(n-1)) &
                 high_sym(n) = code_group /= code_group_old .OR. high_sym(n)

           ELSE
              high_sym(n)=.TRUE.
           ENDIF
           code_group_old=code_group
        ENDIF
        IF (search_sym) WRITE(stdout,'(/,5x,70("*"))')
!     WRITE(stdout,'(2i5, 3f15.5,l5)') n, qcode_group(n), q(:,n), high_sym(n)
     ELSEIF (lo_to_split) THEN
!
!  At gamma the group is the point group of the solid
!
        qcode_group(n)=code_group_save
        qcode_group_ext(n)=0
        ptypeq(:,n)=1
        lprojq(n)=0
        gaugeq(:,n)=0.0_DP

        IF (search_sym) THEN
           WRITE(stdout, '(/,5x,"Mode symmetry analysis not available &
                                                      &for this point")') 
           WRITE(stdout, '(/,5x,70("*"))')
        ENDIF
     ENDIF
  END DO
  !
  IF (flfrq.NE.' ') CALL write_bands(nq, 3*nat, disp_q, freq_save, 1.0_DP, &
                                                          filefrq)
  !
  !  If the force constants are in the xml format we write also
  !  the file with the representations of each mode
  !
  IF (flfrq.NE.' '.AND.xmldyn.AND.search_sym) THEN
     ALLOCATE(aux_ind(nq))
     ALLOCATE(same_next(nq))
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

     DO n=1, nq
        IF (n==nq) THEN
           same_next(n)=.FALSE.
        ELSE
           same_next(n)= same_star(nsym, s, disp_q(1,n), disp_q(1,n+1), at)
        ENDIF
     ENDDO

     filename=TRIM(filefrq)//'.rap'
     CALL write_representations(nq, 3*nat, disp_q, num_rap_mode, high_sym,  &
                       qcode_group, aux_ind, qcode_group_ext, ptypeq, lprojq, &
                       same_next, gaugeq, filename, 0, nq)

     DEALLOCATE(same_next)
     DEALLOCATE(aux_ind)
  ENDIF
  !
  DEALLOCATE (w2) 
  DEALLOCATE (freq_save) 
  DEALLOCATE (z_save) 
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (name_rap_mode)
  DEALLOCATE (high_sym)
  DEALLOCATE (qcode_group)
  DEALLOCATE (ptypeq)
  DEALLOCATE (qcode_group_ext)
  DEALLOCATE (lprojq)
  DEALLOCATE (rtau)
  !
  RETURN
END SUBROUTINE write_ph_dispersions
!
SUBROUTINE find_representations_mode_q ( nat, ntyp, xq, w2, u, tau, ityp, &
                  amass, num_rap_mode, nspin_mag, qcode_old )

  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at, bg
  USE symm_base,    ONLY : s, irt, nsym, time_reversal, copy_sym, &
                           s_axis_to_cart, inverse_s
  USE lr_symm_base, ONLY : gi, nsymq, rtau
  USE control_ph,   ONLY : search_sym
  USE ph_symmetry,  ONLY : manage_ph_symmetry

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nat, ntyp, nspin_mag, qcode_old
  REAL(DP), INTENT(IN) :: xq(3), amass(ntyp), tau(3,nat)
  REAL(DP), INTENT(IN) :: w2(3*nat)
  INTEGER, INTENT(IN)  :: ityp(nat)
  COMPLEX(DP), INTENT(IN) :: u(3*nat,3*nat)
  INTEGER, INTENT(OUT) :: num_rap_mode(3*nat)

  REAL(DP) :: gimq (3)
  INTEGER :: irotmq
  LOGICAL :: minus_q, sym(48), magnetic_sym
  LOGICAL :: symmorphic_or_nzb
!
!  This routine assumes that u contains the eigenvectors of the dynamical
!  matrix
!
!
!  find the small group of q and set the quantities needed to
!  symmetrize a phonon mode
!
  IF (.NOT.search_sym) RETURN
  time_reversal=(nspin_mag/=4)
  minus_q=.TRUE.
  IF (nspin_mag/=4) minus_q=.FALSE.

  sym(1:nsym)=.true.
  call smallg_q_tpw (xq, 0, at, bg, nsym, s, sym, minus_q)
  nsymq=copy_sym(nsym,sym)
  CALL s_axis_to_cart ()
  CALL set_giq_tpw (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  CALL inverse_s()
  CALL sgam_lr (at, bg, nsymq, s, irt, tau, rtau, nat)
!
!  if the small group of q is non symmorphic,
!  search the symmetries only if there are no G such that Sq -> q+G
!
  CALL manage_ph_symmetry(u, w2, num_rap_mode, xq, .TRUE., -1)

  RETURN
END SUBROUTINE find_representations_mode_q
