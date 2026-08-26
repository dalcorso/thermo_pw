!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine new_solve_eq(iu, flag)
  !-----------------------------------------------------------------------
  !! This routine is a driver for the solution of the linear system which
  !! defines the change of the wavefunction due to an electric field.
  !! It performs the following tasks:
  !! a) computes the bare potential term times \(|\psi\rangle \);
  !! b) adds to it the screening term \(\Delta V_\text{SCF}|psi\rangle\).
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) compute also the SCF part
  !!    of the response Hubbard potential;
  !! c) applies \(P_c^+\) (orthogonalization to valence states);
  !! d) calls \(\texttt{cgsolve_all}\) to solve the linear system;
  !! e) computes \(\Delta \rho\), \(\Delta V_\text{SCF}|psi\rangle\) and
  !!    symmetrizes them;
  !! f) if \(\text{lda_plus_u}=\text{TRUE}\) compute also the response
  !!    occupation matrices dnsscf.
  !! Step b, c, d are done in \(\text{sternheimer_kernel}\).
  !
  USE kinds,                 ONLY : DP
  USE constants,             ONLY : e2, fpi, rytoev
  USE ions_base,             ONLY : nat
  USE cell_base,             ONLY : bg, tpiba2
  USE io_global,             ONLY : ionode
  USE io_files,              ONLY : diropn
  USE mp,                    ONLY : mp_sum
  USE klist,                 ONLY : ltetra, lgauss, xk, ngk, igk_k
  USE gvecs,                 ONLY : doublegrid
  USE gvect,                 ONLY : gg
  USE optical,               ONLY : current_w, fru, chirr, chirz, chizr, &
                                    chizz, epsm1
  USE freq_ph,               ONLY : fiu
  USE fft_base,              ONLY : dfftp, dffts
  USE lsda_mod,              ONLY : lsda, current_spin, isk
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : save_buffer, get_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : vkb
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : nspin_mag, noncolin, domag
  USE paw_variables,         ONLY : okpaw
  USE ldaU,                  ONLY : lda_plus_u
  USE qpoint,                ONLY : xq
  USE units_ph,              ONLY : lrdrho, iudrho, lrbar, iubar
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE output,                ONLY : fildrho
  USE control_ph,            ONLY : ext_recover, zeu
  USE recover_mod,           ONLY : write_rec, read_rec
  USE qpoint,                ONLY : nksq, ikks, ikqs
  USE control_flags,         ONLY : use_gpu
  USE control_lr,            ONLY : lgamma, convt, rec_code_read, rec_code, where_rec
  USE uspp_init,             ONLY : init_us_2
  USE dfpt_type,             ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data
  USE dfpt_kernels_tpw,      ONLY : dfpt_kernel_tpw
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE eqv,                   ONLY : dvpsi, evq, dpsi
  USE fft_interfaces,        ONLY : fwfft
  USE io_global,             ONLY : stdout
  !
  IMPLICIT NONE
  !
  EXTERNAL :: stop_smoothly_ph

  INTEGER, INTENT(IN) :: iu   ! frequency to compute in the list
  INTEGER, INTENT(IN) :: flag ! if 1 compute the charge-charge and
                              ! charge magnetization responses
                              ! if 2 and lsda computes the magnetization
                              ! magnetization response
  !
  LOGICAL :: exst
  !!
  INTEGER :: ikk, npw, iter0, ipert, ik, is, npert, nrec, npwq, ikq
  !! counters
  INTEGER :: nnrs
  !! dffts%nnr, defined used for acc
  REAL(DP) :: dr2
  !! self-consistency error
  TYPE(dfpt_data_type) :: dfpt_data
  !! Data that describes linear response quantities
  COMPLEX(DP), ALLOCATABLE :: drhoscf_aux(:,:,:)
  COMPLEX(DP) :: w
  REAL(DP) :: xqmod2
  LOGICAL :: finite_frequency
  !
  CALL start_clock ('solve_eq')
  npert=1
  !
  CALL allocate_dfpt_data(dfpt_data, npert)
  !
  nnrs = dffts%nnr
  !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:nnrs, 1:nspin_mag, 1:npert))

  w=CMPLX(fru(iu),fiu(iu))
  finite_frequency=(ABS(w) > 1.D-7)
  !
  if (rec_code_read == -20.AND.ext_recover) then
     ! restarting in Electric field calculation
     CALL read_rec(dr2, iter0, dfpt_data)
  else if (rec_code_read > -20 .AND. rec_code_read <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
     dr2 = 0.d0
  endif
  !
  IF (rec_code_read > -20) convt=.TRUE.
  !
  if (convt) go to 155
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  IF ((lgauss.AND..NOT.finite_frequency)) &
          CALL errore ('solve_eq', 'insert a finite frequency', 1)
  !
  ! Compute P_c^+ x psi for all polarization and k points and store in buffer
  !
  DO ik = 1, nksq
     DO ipert = 1, npert
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq = ngk(ikq)
        IF (lsda) current_spin = isk(ikk)
        nrec = (ipert - 1) * nksq + ik
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        IF (nksq > 1) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikk)
           !$acc update device(evc)
           CALL get_buffer(evq, lrwfc, iuwfc, ikq)
           !$acc update device(evq)
        ENDIF
        !
        CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, use_gpu)
        !$acc update host(vkb)
        !
        ! computes P_c^+ x psi_kpoint, written to buffer iuebar.
        !
        !
        !  At the first iteration dvbare_q*psi_kpoint is calculated
        !  and written to file
        !
        CALL dveqpsi_us (ik)
        !
        !  with flag=2 the perturbation is a magnetic field along z
        !
        IF (lsda.AND.current_spin==2.AND.flag==2) dvpsi=-dvpsi
        CALL save_buffer (dvpsi, lrbar, iubar, nrec)
        !
     ENDDO ! ipert
  ENDDO ! ik
  !
  ! Set records for restart
  !
  rec_code = -20  ! Electric field
  where_rec = 'solve_e...'
  !
  !   Solve DFPT fixed-point equation
  !
  CALL dfpt_kernel_tpw('EELS', npert, iter0, lrbar, iubar, dr2,  &
              dfpt_data, 1, 0, write_rec_callback = write_rec, &
              stop_callback = stop_smoothly_ph, w_freq=w)
  !
  IF ( fildrho /= ' ') THEN
     IF ( ionode ) THEN
        INQUIRE (UNIT = iudrho, OPENED = exst)
        IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
        CALL diropn (iudrho, TRIM(fildrho)//'.Ew', lrdrho, exst)
     END IF
     DO ipert=1,npert
        CALL davcio_drho(dfpt_data%drhop(1,1,ipert),lrdrho, iudrho,ipert,+1)
     END DO
  END IF
155 CONTINUE
!
!  compute here the susceptibility and the inverse of the dielectric
!  constant
!
!  CALL compute_susceptibility(drhoscfout)

  DO is=1,nspin_mag
     CALL fwfft ('Rho', dfpt_data%drhop(:,is,1), dfftp)
  END DO
  IF (flag==1) THEN
     chirr(iu)=(0.0_DP,0.0_DP)
     chizr(iu)=(0.0_DP,0.0_DP)
     epsm1(iu)=(0.0_DP,0.0_DP)
  ELSE
     chirz(iu)=(0.0_DP,0.0_DP)
     chizz(iu)=(0.0_DP,0.0_DP)
  ENDIF
  xqmod2=(xq(1)**2+xq(2)**2+xq(3)**2)*tpiba2
  IF (ABS(gg(1))<1.d-8) THEN
     IF (flag==1) THEN
        chirr(iu) = dfpt_data%drhop(dfftp%nl(1),1,1) 
        IF (lsda) chirr(iu) = chirr(iu) + dfpt_data%drhop(dfftp%nl(1),2,1)
        epsm1(iu) = CMPLX(1.0_DP,0.0_DP)+ chirr(iu)*fpi*e2/xqmod2
        IF (lsda) chizr(iu) = dfpt_data%drhop(dfftp%nl(1),1,1) - &
                              dfpt_data%drhop(dfftp%nl(1),2,1)
     ELSE IF (lsda) THEN
        chizz(iu)=dfpt_data%drhop(dfftp%nl(1),1,1)-dfpt_data%drhop(dfftp%nl(1),2,1)
        chirz(iu)=dfpt_data%drhop(dfftp%nl(1),1,1)+dfpt_data%drhop(dfftp%nl(1),2,1)
     END IF
  END IF

  IF (flag==1) THEN
     CALL mp_sum(epsm1(iu),intra_bgrp_comm)        
     CALL mp_sum(chirr(iu),intra_bgrp_comm)        
     CALL mp_sum(chizr(iu),intra_bgrp_comm)        
  ELSE
     CALL mp_sum(chizz(iu),intra_bgrp_comm)        
     CALL mp_sum(chirz(iu),intra_bgrp_comm)        
  END IF

  IF (flag==1) THEN
     WRITE(stdout, '(/,6x,"Inverse dielectric constant at &
                        &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
     WRITE(stdout, '(46x,f9.4," +",f9.4," i eV")') current_w * rytoev
     WRITE(stdout,'(/,6x,"epsilon^-1(q,w) =",2f15.6)') epsm1(iu)

     WRITE( stdout, '(/,5x,"Charge-charge susceptibility:")') 

     WRITE(stdout,'(/,6x,"chirr(q,w) =",2f15.6)') chirr(iu)
     IF (lsda) THEN
        WRITE(stdout,'(/,6x,"m_z-charge susceptibility:")')
        WRITE(stdout,'(/,6x,"chizr(q,w) =",2f15.6)') chizr(iu)
     ENDIF
  ELSE IF (lsda) THEN
     WRITE( stdout, '(/,6x,"m_z - m_z susceptibility at &
                     &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
     WRITE( stdout, '(43x,f9.4," +",f9.4," i eV")') current_w * rytoev
     WRITE(stdout,'(/,6x,"chizz(q,w) =",2f15.6)') chizz(iu)
     WRITE(stdout,'(/,6x,"chirz(q,w) =",2f15.6)') chirz(iu)
  END IF
!
  !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)
  CALL deallocate_dfpt_data(dfpt_data)
  !
  call stop_clock ('solve_eq')
  return
end subroutine new_solve_eq
