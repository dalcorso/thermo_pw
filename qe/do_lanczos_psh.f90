!
! Copyright (C) 2017-2018 Andrea Dal Corso 
! This routine has been obtained by modifyng the routine solve_e and
! its generalizations solve_e_fpol and solve_eq,
! by substituting the call to the conjugate gradient with 
! a Lanczos step: application of (H-eS) and possibly P_c^+ dV_Hxc.
! The application of S^-1 has also been included where needed.
!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE do_lanczos_psh()
  !-----------------------------------------------------------------------
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field.
  !    The electric field perturbation can be x, y, and z (when q=0) or e^{iqr}
  !    when q is different form gamma.
  !    It performs the following tasks:
  !     a) computes the bare potential term  V_ext | psi > at the first
  !        iteration.
  !     b) applies P_c^+ (orthogonalization to valence states) and possibly
  !        S^-1.
  !     c) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !     d) computes (H-Se) dpsi
  !     e) adds to it the screening term P_c^+ Delta V_{SCF} | psi > if needed
  !     f) apply S^-1 if needed
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE cell_base,             ONLY : at
  USE io_global,             ONLY : stdout, ionode
  USE klist,                 ONLY : xk, wk, ngk, igk_k, lgauss, ltetra
  USE qpoint,                ONLY : nksq, ikks, ikqs
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE fft_helper_subroutines, ONLY : fftx_ntgrp
  USE fft_interfaces,        ONLY : fft_interpolate
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                 ONLY : nbnd, npwx,  et
  USE wavefunctions,         ONLY : evc
  USE eqv,                   ONLY : dpsi, dvpsi, evq
  USE becmod,                ONLY : becp, calbec
  USE scf,                   ONLY : rho
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : nhm
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_symmetry,          ONLY : paw_desymmetrize
  USE paw_add_symmetry,      ONLY : paw_deqsymmetrize
  USE control_ph,            ONLY : lnoloc, ext_recover, recover
  USE control_lr,            ONLY : alpha_pv, nbnd_occ, lgamma, lgamma_gamma, &
                                    rec_code_read, rec_code, convt
  USE lrus,                  ONLY : int3, int3_paw
  USE dv_of_drho_lr,         ONLY : dv_of_drho
  USE dv_of_drho_clf,        ONLY : dv_of_drho_nlf
  USE lr_global,             ONLY : rpert, evc0, evq0, sevq0, d0psi, d0psi2
  USE lr_lanczos,            ONLY : lanczos_steps, evc1, sevc1, evc1_new,    &
                                    lanczos_restart_step
  USE uspp_init,             ONLY : init_us_2
  USE recover_mod,           ONLY : read_rec, write_rec
  USE units_lr,              ONLY : lrwfc, iuwfc
  USE buffers,               ONLY : get_buffer
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm, ntask_groups
  USE mp,                    ONLY : mp_sum

  IMPLICIT NONE

  COMPLEX(DP) , ALLOCATABLE, TARGET :: &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  COMPLEX(DP) , POINTER ::             &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  COMPLEX(DP) , ALLOCATABLE ::         &
                   dbecsum(:,:,:,:),   & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   aux1 (:,:),         &   ! auxiliary space to apply potential
                   tg_dv(:,:),         &   ! task group variables
                   tg_psic(:,:)          

  INTEGER :: kter, iter, iter0, ipol, jpol, ibnd, ik, ikp, ikk, ikq, is, &
             npw, npwq, incr, v_siz, ig
  ! counters or indices
  LOGICAL :: lmet

  REAL(DP) :: weight, alpha_pv0, anorm, dr2
  ! weight of k points and store alpha_pv
  REAL(DP) :: tcpu, get_clock
  LOGICAL :: exst

  CALL start_clock ('do_lanczos')
  alpha_pv0=alpha_pv
  convt=.FALSE.
  lmet = (lgauss .OR. ltetra)
  IF (lmet.AND.lgamma) THEN
     CALL errore('do_lanczos','only q/=0 allowed for metals',1)
  ELSE
     IF (lmet) alpha_pv=0.0_DP
  ENDIF


  ALLOCATE (dvscfin( dfftp%nnr, nspin_mag, rpert))
  IF (doublegrid) THEN
     ALLOCATE (dvscfins(dffts%nnr, nspin_mag, rpert))
  ELSE
     dvscfins => dvscfin
  ENDIF
  ALLOCATE (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, rpert))
  IF (noncolin) ALLOCATE (dbecsum_nc (nhm, nhm, nat, nspin, rpert))

  ALLOCATE (aux1(dffts%nnr,npol))
  !
  !  This routine is task group aware
  !
  incr=1
  IF ( dffts%has_task_groups ) THEN
     !
     v_siz =  dffts%nnr_tg
     ALLOCATE( tg_dv   ( v_siz, nspin_mag ) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = fftx_ntgrp(dffts)
     !
  ENDIF
  IF (recover) THEN
     IF (rec_code_read == -20.AND.ext_recover) &
        CALL read_rec(dr2, iter0, rpert, dvscfin, dvscfins)
     iter0=0
     CALL lr_restart_tpw (iter0, recover)
     iter0=iter0+1
  ELSE
     iter0=1
     dr2=0.0_DP
  ENDIF
  !
  !   The outside loop is over the lanczos steps
  !
  DO kter = iter0, lanczos_steps+1

!     write(6,*) 'kter', kter
     iter = kter
!
!   A loop over k points
!
     DO ik = 1, nksq
        ikk=ikks(ik)
        ikq=ikqs(ik)
        npw=ngk(ikk)
        npwq=ngk(ikq)
        IF (lsda) current_spin = isk (ikk)
        !
        CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
        !
        IF (iter==iter0) THEN
           !
           ! at the first iteration reads unperturbed wavefuctions psi_k 
           ! in G_space, for all bands. When q /= 0 reads also evq 
           !
           IF (nksq>1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
           evc0(:,:,ik)=evc(:,:)
           IF (.NOT.lgamma .AND.nksq>1) THEN
              CALL get_buffer (evq, lrwfc, iuwfc, ikq)
              evq0(:,:,ik)=evq(:,:)
           ENDIF
           IF (okvan) THEN
              CALL calbec (npwq, vkb, evq, becp, nbnd)
              CALL s_psi (npwx, npwq, nbnd, evq, sevq0(1,1,ik))
           ELSE
              sevq0(:,:,ik)=evq(:,:)
           ENDIF
        ELSE
           !
           !  At the following iterations only copy the correct wavefunction
           !
           evc(:,:)=evc0(:,:,ik)
           IF (.NOT.lgamma .AND.nksq.gt.1) evq(:,:)=evq0(:,:,ik)
        ENDIF
        !
        IF (iter>1) CALL g2_kin(ikq)
        !
        !  And a loop over the perturbations
        !
        DO ipol = 1, rpert
           ikp = ik + nksq * (ipol - 1)  
           
           IF (iter == iter0) THEN
              !
              ! At the first iteration compute the right hand side of the
              ! linear systems and the starting vectors.
              ! computes P_c^+ V_ext u_kv into dvpsi array.
              ! At q=0 V_ext is r_1, r_2 or r_3 (r in crystal coordinates)
              ! at q/=0 it is e^{iqr}. In the US case there is also an
              ! augmentation term computed here.
              !
              IF (lgamma) THEN
                 CALL dvpsi_e(ik,ipol)
              ELSE
                 CALL dveqpsi_us(ik)
                 d0psi2(:,:,ik)=dvpsi(:,:)
              ENDIF

              !
              ! Orthogonalize dvpsi to valence states: Apply P_c^+ and change
              ! sign.
              !
              CALL orthogonalize(dvpsi, evq, ikk, ikq, sevq0(:,:,ik), npwq, &
                                                                      .TRUE.)
              !
              !  save here P_c^+ V_ext u_kv needed to compute the 
              !  susceptibility
              !
              d0psi(:,:,ik,ipol)=-dvpsi(:,:)
              !
              !  In the US case apply S^-1
              !
              IF (okvan) THEN
                 dpsi(:,:) = (0.0d0, 0.0d0)
                 CALL lr_sm1_initialize()
                 CALL lr_sm1_psi_tpw(ik,npwx,npwq,nbnd_occ(ikk),dvpsi,dpsi)
                 dvpsi(:,:) = dpsi(:,:)
              ENDIF

              !
              ! At the first iteration evc1 is not known and is set 
              ! to S^-1 P_c^+ V_ext u_kv
              !
              IF (iter==1) evc1(:,:,ikp,1)=-dvpsi(:,:)
              !
           ENDIF
           IF (iter>1) THEN
              !
              !  Here we compute (H-eS) psi^1, at the iteration iter-1 
              !  Since ch_psi_all apply also alpha_pv Q, we set alpha_pv
              !  to zero.
              !
              alpha_pv=0.0_DP
              CALL ch_psi_all (npwq, evc1(:,:,ikp,1), sevc1(:,:,ikp,1), &
                                    et(:,ikk), ik, nbnd_occ(ikk))
              IF (.NOT.lmet) alpha_pv=alpha_pv0
!
!  We are still computing the iteration iter-1 and we have to skip the
!  interaction if is odd
!
              IF (MOD(iter-1,2)/=0) GOTO 111
              !
              ! Now apply P_c^+ dV_Hxc. Note that since this is the
              ! iter-1 step, made at the iter step ich and inch must be
              ! reversed.
              ! First apply dV_Hxc.
              !
              IF ( dffts%has_task_groups ) THEN
                 IF (noncolin) THEN
                    CALL tg_cgather( dffts, dvscfins(:,1,ipol), tg_dv(:,1))
                    IF (domag) THEN
                       DO jpol=2,4
                          CALL tg_cgather( dffts, dvscfins(:,jpol,ipol), &
                                                          tg_dv(:,jpol))
                       ENDDO
                    ENDIF
                 ELSE
                    CALL tg_cgather( dffts, dvscfins(:,current_spin,ipol), &
                                                             tg_dv(:,1))
                 ENDIF
              ENDIF
              dvpsi=(0.0_DP,0.0_DP)
              DO ibnd = 1, nbnd_occ (ikk), incr
                 IF ( dffts%has_task_groups ) THEN
                    CALL cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, &
                                         nbnd_occ (ikk) )
                    CALL apply_dpot(v_siz, tg_psic, tg_dv, 1)
                    CALL cft_wave_tg (ik, dvpsi, tg_psic, -1, v_siz, ibnd, &
                                         nbnd_occ (ikk))
                 ELSE
                    CALL cft_wave (ik, evc (1, ibnd), aux1, +1)
                    CALL apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipol), &
                                                             current_spin)
                    CALL cft_wave (ik, dvpsi (1, ibnd), aux1, -1)
                 ENDIF
              ENDDO
              !
              !  Add the US contribution if any
              !
              CALL adddvscf(ipol,ik)
              !
              ! Apply -P_c^+
              !
              CALL orthogonalize(dvpsi, evq, ikk, ikq, &
                                     & sevq0(:,:,ik), npwq, .TRUE.)
              !
              !  the orthogonalize routine changes the sign of dvpsi, so here
              !  we subtract.
              !
              sevc1(:,:,ikp,1) = sevc1(:,:,ikp,1) - dvpsi(:,:)
              !
              ! Ultrasoft case: apply the S^{-1} operator.
              ! evc1_new = S^{-1} * sevc1_new
              ! If not ultrasoft: evc1_new = sevc1_new
              !
111           CONTINUE
              IF (okvan) THEN
                 CALL lr_sm1_psi_tpw (ik,npwx,npwq,nbnd_occ(ikk), &
                                sevc1(1,1,ikp,1), evc1_new(1,1,ikp,1))
              ELSE
                 evc1_new(:,:,ikp,1)=sevc1(:,:,ikp,1)
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     !
     !  here we do the lanczos step of iteration iter-1
     !
     IF (iter > 1) CALL psh_lanczos_step(iter-1,.FALSE.)
     !
     !  and save the results to restart each lanczos_restart_step
     !
     IF ((lanczos_restart_step>0 .AND. iter /= 1 .AND. &
         ( mod(iter-1,lanczos_restart_step)==0)).OR.(iter==lanczos_steps+1)) &
                                CALL lanczos_write_restart_tpw(iter-1)
     !
     !  The independent particle approximation corresponds to set to zero 
     !  dV_Hxc
     !  The interaction is not computed at the odd iteraction. Now the
     !  iteration if iter.
     !
     IF ((lnoloc.AND.lgamma).OR.MOD(iter,2)/=0) THEN
        dvscfins=(0.d0,0.d0)
        IF (okvan) int3=(0.0_DP,0.0_DP)
        CYCLE
     ENDIF
     !
     !  compute the charge density with the updated solution,
     !  ready for the next iteration.
     !
     dvscfin(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)
     IF (noncolin) dbecsum_nc=(0.d0,0.d0)
     !
     !   Another loop on k points
     !
     DO ik = 1, nksq
        ikk=ikks(ik)
        ikq=ikqs(ik)
        npw=ngk(ikk)
        npwq = ngk(ikq)
        IF (lsda) current_spin = isk (ikk)
        !
        !  In the US case we need vkb to compute dbecsum
        !
        IF (okvan) CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
        !
        !  and on perturbations
        !
        DO ipol = 1, rpert
           ikp = ik + nksq * (ipol - 1)  
           !
           evc(:,:)=evc0(:,:,ik)
           !
           !  This is the charge at iteration iter, so we use these 
           !  perturbed wavefunctions
           !
           dpsi(:,:)=evc1(:,:,ikp,1)
           !
           weight=wk(ikk)
           IF (noncolin) THEN
              CALL incdrhoscf_nc(dvscfin(1,1,ipol), weight, ik,         &
                                 dbecsum_nc(1,1,1,1,ipol), dpsi, 1.0_DP)
           ELSE
              CALL incdrhoscf(dvscfin(1,current_spin,ipol), weight, ik, &
                                 dbecsum(1,1,current_spin,ipol), dpsi)
           ENDIF
        ENDDO   ! on perturbations
     ENDDO      ! on k points
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions
     !  coming from each slice of bands
     !
     IF (noncolin) THEN
        CALL mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        CALL mp_sum ( dbecsum, intra_bgrp_comm )
     END IF
     !
     !  Interpolates the smooth part of the induced charge in the thick
     !  grid if needed
     !
     IF (doublegrid) THEN
        DO ipol=1,rpert
           DO is=1,nspin_mag
              CALL fft_interpolate (dffts, dvscfin(:,is,ipol), dfftp, &
                                                      dvscfin(:,is,ipol))
           ENDDO
        ENDDO
     ENDIF
     !
     !  Rotate dbecsum_nc in the spin-orbit case
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, rpert)
     !
     !  And add the augmentation part of the induced charge    
     !
     IF (lgamma) THEN
        CALL addusddense (dvscfin, dbecsum)
     ELSE
        CALL addusddenseq (dvscfin, dbecsum)
     ENDIF
     !
     !  Collect the contribution of all pools
     !
     CALL mp_sum ( dvscfin, inter_pool_comm )
     IF (okpaw) CALL mp_sum ( dbecsum, inter_pool_comm )
     !
     !   dvscfin contains the (unsymmetrized) linear charge response
     !   for the rpert perturbations (3 in the lgamma case and one in the
     !   q /= 0 case) - symmetrize it
     !
     IF (.NOT.lgamma_gamma) THEN
        IF (lgamma) THEN
           CALL symmetrize_drho(dvscfin, dbecsum, 0, 3, 2 )
!           CALL psyme (dvscfin)
!           IF ( noncolin.and.domag ) CALL psym_dmage(dvscfin)
        ELSE
           CALL symmetrize_drho(dvscfin, dbecsum, 0, 3, 3 )
!           CALL psymeq (dvscfin)
!           IF ( noncolin.AND.domag ) CALL psym_dmageq(dvscfin)
        ENDIF
     ENDIF
     !
     !   calculate the corresponding linear potential response
     !
     DO ipol=1,rpert
        IF (lnoloc) THEN
           CALL dv_of_drho_nlf (dvscfin (1, 1, ipol))
        ELSE
           CALL dv_of_drho (dvscfin (1, 1, ipol))
        ENDIF 
     ENDDO
     !
     !  And interpolate the potential on the smooth grid if needed
     !
     IF (doublegrid) THEN
        DO ipol = 1, rpert
           DO is=1,nspin_mag
              CALL fft_interpolate (dfftp, dvscfin(:,is,ipol), dffts, &
                                           dvscfins(:,is,ipol))
           ENDDO
        ENDDO
     ENDIF
!
!    In the PAW case computes the change of the D coefficients, after
!    symmetrization of the dbecsum terms
!
     IF (okpaw) THEN
        IF (noncolin.AND.domag) THEN
!           CALL PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
        ELSE
           dbecsum=2.0_DP * dbecsum
           IF (lgamma) THEN
              IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
           ELSE
              CALL PAW_deqsymmetrize(dbecsum)
           ENDIF
           CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,rpert)
        ENDIF
     ENDIF
!
!   In the US and PAW case computes the integral of dV_Hxc and the 
!   augmentation function. This quantity is needed in adddvscf.
!
     CALL newdq(dvscfin,rpert)

     tcpu = get_clock ('PHONON')
     WRITE( stdout, '(/,5x," iter # ",i8," total cpu time :",f8.1,&
                                  &" secs ")') iter, tcpu
     !
     FLUSH( stdout )
     IF ((lanczos_restart_step>0 .AND. iter /= 1 .AND. &
         (mod(iter-1,lanczos_restart_step)==0)).OR.iter==lanczos_steps+1) THEN
        rec_code=-20
        CALL write_rec('solve_e...', 0, dr2, iter, .FALSE., rpert, dvscfin)
     ENDIF
     !
  ENDDO  ! Lanczos iterations

  DEALLOCATE (aux1)
  DEALLOCATE (dbecsum)
  IF (doublegrid) DEALLOCATE (dvscfins)
  DEALLOCATE (dvscfin)
  IF (noncolin) DEALLOCATE(dbecsum_nc)
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE( tg_dv  )
     DEALLOCATE( tg_psic)
  ENDIF

  alpha_pv=alpha_pv0
  CALL stop_clock ('do_lanczos')

  RETURN
END SUBROUTINE do_lanczos_psh
