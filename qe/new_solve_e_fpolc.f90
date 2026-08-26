!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine new_solve_e_fpolc(iu)
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
  USE ions_base,             ONLY : nat
  USE cell_base,             ONLY : bg
  USE io_global,             ONLY : ionode
  USE io_files,              ONLY : diropn
  USE mp,                    ONLY : mp_sum
  USE klist,                 ONLY : ltetra, lgauss, xk, ngk, igk_k
  USE gvecs,                 ONLY : doublegrid
  USE gvect,                 ONLY : gg
  USE optical,               ONLY : current_w, fru
  USE freq_ph,               ONLY : fiu
  USE fft_base,              ONLY : dfftp, dffts
  USE lsda_mod,              ONLY : lsda, current_spin, isk
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : vkb
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : nspin_mag, noncolin, domag
  USE paw_variables,         ONLY : okpaw
  USE ldaU,                  ONLY : lda_plus_u
  USE units_ph,              ONLY : lrdrho, iudrho, lrebar, iuebar
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE magnetic_charges,      ONLY : alpha_me
  USE output,                ONLY : fildrho
  USE control_ph,            ONLY : ext_recover, zeu
  USE recover_mod,           ONLY : write_rec, read_rec
  USE qpoint,                ONLY : nksq, ikks
  USE control_lr,            ONLY : lgamma, convt, rec_code_read, rec_code, where_rec
  USE uspp_init,             ONLY : init_us_2
  USE dfpt_type,             ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data
  USE dfpt_kernels_tpw,      ONLY : dfpt_kernel_tpw
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE eqv,                   ONLY : dvpsi
  USE fft_interfaces,        ONLY : fwfft
  !
  IMPLICIT NONE
  !
  EXTERNAL :: stop_smoothly_ph

  INTEGER, INTENT(IN) :: iu   ! frequency to compute in the list
  !
  LOGICAL :: exst
  !!
  COMPLEX(DP) :: drhoscfs (dffts%nnr, nspin_mag, 3)
  !!
  INTEGER :: ikk, npw, iter0, ipol, ik, is, icart
  !! counters
  INTEGER :: nnr
  !! dffts%nnr, defined used for acc
  REAL(DP) :: dr2
  !! self-consistency error
  TYPE(dfpt_data_type) :: dfpt_data
  !! Data that describes linear response quantities
  COMPLEX(DP), ALLOCATABLE :: drhoscf_aux(:,:,:)
  COMPLEX(DP) :: alpha_work(3,3), w
  !
  call start_clock ('solve_e')
  !
  !  This routine is task group aware
  !
  CALL allocate_dfpt_data(dfpt_data, 3)
  !
  nnr = dffts%nnr
  !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:nnr, 1:nspin_mag, 1:3))

  w=CMPLX(fru(iu),fiu(iu))
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
  if ( (lgauss .or. ltetra) .or..not.lgamma) call errore ('solve_e', &
       'called in the wrong case', 1)
  !
  ! Compute P_c^+ x psi for all polarization and k points and store in buffer
  !
  DO ik = 1, nksq
     DO ipol = 1, 3
        ikk = ikks(ik)
        npw = ngk(ikk)
        IF (lsda) current_spin = isk(ikk)
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        IF (nksq > 1) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikk)
           !$acc update device(evc)
        ENDIF
        !
        CALL init_us_2(npw, igk_k(1, ikk), xk(1, ikk), vkb, .true.)
        !$acc update host(vkb)
        !
        ! computes P_c^+ x psi_kpoint, written to buffer iuebar.
        !
        CALL dvpsi_e(ik, ipol)
        !
     ENDDO ! ipol
  ENDDO ! ik
  !
  ! Set records for restart
  !
  rec_code = -20  ! Electric field
  where_rec = 'solve_e...'
  !
  !   Solve DFPT fixed-point equation
  !
  CALL dfpt_kernel_tpw('EPSILC', 3, iter0, lrebar, iuebar, dr2,  &
              dfpt_data, 1, 0, write_rec_callback = write_rec, &
              stop_callback = stop_smoothly_ph, w_freq=w)
  !
  IF ( fildrho /= ' ') THEN
     IF ( ionode ) THEN
        INQUIRE (UNIT = iudrho, OPENED = exst)
        IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
        CALL diropn (iudrho, TRIM(fildrho)//'.Ew', lrdrho, exst)
     END IF
     DO ipol=1,3
        CALL davcio_drho(dfpt_data%drhop(1,1,ipol),lrdrho, iudrho,ipol,+1)
     END DO
  END IF
155 CONTINUE
  !
  ! The dfpt_data are deallocated but the smooth charge is used in the
  ! zstar. We copy it here
  !
  drhoscfs=dfpt_data%drhos
  !
  !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)
  CALL deallocate_dfpt_data(dfpt_data)
  !
  call stop_clock ('solve_e')
  return
end subroutine new_solve_e_fpolc
