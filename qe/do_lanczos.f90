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
SUBROUTINE do_lanczos()
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
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,              ONLY : domag
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag
  USE wvfct,                 ONLY : nbnd, npwx,  et
  USE wavefunctions_module,  ONLY : evc
  USE eqv,                   ONLY : dpsi, dvpsi, evq
  USE becmod,                ONLY : becp, calbec
  USE scf,                   ONLY : rho
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : nhm
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_symmetry,          ONLY : paw_desymmetrize
  USE paw_add_symmetry,      ONLY : paw_deqsymmetrize
  USE control_ph,            ONLY : lnoloc,  lgamma_gamma
  USE control_lr,            ONLY : alpha_pv, nbnd_occ, lgamma
  USE lrus,                  ONLY : int3, int3_paw
  USE dv_of_drho_lr,         ONLY : dv_of_drho
  USE lr_lanczos,            ONLY : rpert, lanczos_steps, evc0, evq0, sevc0, &
                                    evc1, sevc1, evc1_new, d0psi, d0psi2
  USE units_ph,              ONLY : lrwfc, iuwfc
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
             ich, inch, npw, npwq, incr, v_siz, ig
  ! counters or indices
  LOGICAL :: lmet

  REAL(DP) :: weight, alpha_pv0, anorm
  ! weight of k points and store alpha_pv
  REAL(DP) :: tcpu, get_clock

  CALL start_clock ('do_lanczos')

  alpha_pv0=alpha_pv
  lmet = (lgauss .OR. ltetra)
  IF (lmet.AND.lgamma) THEN
     CALL errore('do_lanczos','only q/=0 allowed for metals',1)
  ELSE
     IF (lmet) alpha_pv=0.0_DP
  ENDIF

  iter0=0

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
  IF ( dffts%have_task_groups ) THEN
     !
     v_siz =  dffts%nnr_tg
     ALLOCATE( tg_dv   ( v_siz, nspin_mag ) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = fftx_ntgrp(dffts)
     !
  ENDIF
  !
  !   The outside loop is over the lanczos steps
  !
  DO kter = 1, lanczos_steps+1

!     write(6,*) 'kter', kter
     iter = kter + iter0
!
!  ich and inch selects the component that has (inch) the Hxc potential 
!  depending on the iterations they must be switched
!
     IF ( MOD(iter,2)==0) THEN
        ich=2
        inch=1
     ELSE
        ich=1
        inch=2
     ENDIF
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
        IF (iter==1) THEN
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
              CALL s_psi (npwx, npwq, nbnd, evq, sevc0(1,1,ik))
           ELSE
              sevc0(:,:,ik)=evq(:,:)
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
           
           IF (iter == 1) THEN
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
                 anorm = DSQRT(at(1,ipol)**2 + at(2,ipol)**2 + at(3,ipol)**2)
                 dvpsi(:,:) = dvpsi(:,:) / anorm
              ELSE
                 CALL dveqpsi_us(ik)
                 d0psi2(:,:,ik)=dvpsi(:,:)
              ENDIF
              !
              ! Orthogonalize dvpsi to valence states: Apply P_c^+ and change
              ! sign.
              !
              CALL orthogonalize(dvpsi, evq, ikk, ikq, sevc0(1,1,ik), npwq, &
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
              evc1(:,:,ikp,1)=-dvpsi(:,:)
              evc1(:,:,ikp,2)=-dvpsi(:,:)
              !
           ELSE
              !
              !  Here we compute (H-eS) psi^1, at the iteration iter-1 
              !  Since ch_psi_all apply also alpha_pv Q, we set alpha_pv
              !  to zero.
              !
              alpha_pv=0.0_DP
              CALL ch_psi_all (npwq, evc1(:,:,ikp,ich), sevc1(:,:,ikp,ich), &
                                    et(:,ikk), ik, nbnd_occ(ikk))
              CALL ch_psi_all (npwq, evc1(:,:,ikp,inch), sevc1(:,:,ikp,inch), &
                                    et(:,ikk), ik, nbnd_occ(ikk))
              IF (.NOT.lmet) alpha_pv=alpha_pv0
              !
              ! Now apply P_c^+ dV_Hxc. Note that since this is the
              ! iter-1 step, made at the iter step ich and inch must be
              ! reversed.
              ! First apply dV_Hxc.
              !
              IF ( dffts%have_task_groups ) THEN
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
                 IF ( dffts%have_task_groups ) THEN
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
                                     & sevc0(:,:,ik), npwq, .TRUE.)
              !
              !  the orthogonalize routine changes the sign of dvpsi, so here
              !  we subtract.
              !
              sevc1(:,:,ikp,ich) = sevc1(:,:,ikp,ich) - dvpsi(:,:)
              !
              ! Ultrasoft case: apply the S^{-1} operator.
              ! evc1_new = S^{-1} * sevc1_new
              ! If not ultrasoft: evc1_new = sevc1_new
              !
              IF (okvan) THEN
                 CALL lr_sm1_psi_tpw (ik,npwx,npwq,nbnd_occ(ikk), &
                                sevc1(1,1,ikp,ich), evc1_new(1,1,ikp,ich))
                 CALL lr_sm1_psi_tpw (ik,npwx,npwq,nbnd_occ(ikk), &
                                sevc1(1,1,ikp,inch), evc1_new(1,1,ikp,inch))
              ELSE
                 evc1_new(:,:,ikp,ich)=sevc1(:,:,ikp,ich)
                 evc1_new(:,:,ikp,inch)=sevc1(:,:,ikp,inch)
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     !
     !  here we do the lanczos step of iteration iter-1
     !
     IF (iter > 1) CALL nh_lanczos_step(iter-1,.FALSE.)
     IF (iter==lanczos_steps+1) EXIT
     !
     !  The independent particle approximation corresponds to set to zero 
     !  dV_Hxc
     !
     IF (lnoloc) THEN
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
           dpsi(:,:)=evc1(:,:,ikp,inch)
           !
           weight=wk(ikk)
           IF (noncolin) THEN
              CALL incdrhoscf_nc(dvscfin(1,1,ipol), weight, ik,         &
                                 dbecsum_nc(1,1,1,1,ipol), dpsi)
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
              CALL cinterpolate (dvscfin(1,is,ipol), dvscfin(1,is,ipol), 1)
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
           CALL psyme (dvscfin)
           IF ( noncolin.and.domag ) CALL psym_dmage(dvscfin)
        ELSE
           CALL psymeq (dvscfin)
           IF ( noncolin.AND.domag ) CALL psym_dmageq(dvscfin)
        ENDIF
     ENDIF
     !
     !   calculate the corresponding linear potential response
     !
     DO ipol=1,rpert
        CALL dv_of_drho (dvscfin (1, 1, ipol), .false.)
     ENDDO
     !
     !  And interpolate the potential on the smooth grid if needed
     !
     IF (doublegrid) THEN
        DO ipol = 1, rpert
           DO is=1,nspin_mag
              CALL cinterpolate (dvscfin(1,is,ipol),dvscfins(1,is,ipol),-1)
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
     !
  ENDDO  ! Lanczos iterations

  DEALLOCATE (aux1)
  DEALLOCATE (dbecsum)
  IF (doublegrid) DEALLOCATE (dvscfins)
  DEALLOCATE (dvscfin)
  IF (noncolin) DEALLOCATE(dbecsum_nc)
  IF ( dffts%have_task_groups ) THEN
     DEALLOCATE( tg_dv  )
     DEALLOCATE( tg_psic)
  ENDIF
  alpha_pv=alpha_pv0

  CALL stop_clock ('do_lanczos')

  RETURN
END SUBROUTINE do_lanczos
