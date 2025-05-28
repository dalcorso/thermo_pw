!
! Copyright (C) 2018 Andrea Dal Corso 
! This routine has been obtained by modifyng the routine solve_e 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE do_cg_e(drhoscfs)
  !-----------------------------------------------------------------------
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field.
  !    The electric field perturbation can be x, y, and z (when q=0) or e^{iqr}
  !    when q is different form gamma.
  !    It performs the following tasks:
  !     a) computes the bare potential term  V_ext | psi > at the first
  !        iteration.
  !     b) applies P_c^+ (orthogonalization to valence states) 
  !     c) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !     d) computes (H-Se) dpsi
  !     e) adds to it the screening term P_c^+ Delta V_{SCF} | psi > if needed
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode
  USE klist,                 ONLY : xk, wk, ngk, igk_k, nkstot
  USE qpoint,                ONLY : nksq, ikks, ikqs
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE fft_helper_subroutines, ONLY : fftx_ntgrp
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                 ONLY : nbnd, npwx, g2kin, et
  USE wavefunctions,         ONLY : evc
  USE eqv,                   ONLY : dpsi, dvpsi, evq
  USE becmod,                ONLY : becp, calbec
  USE scf,                   ONLY : rho, v_of_0
  USE uspp,                  ONLY : okvan, vkb, nlcc_any
  USE uspp_param,            ONLY : nhm
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_symmetry,          ONLY : paw_desymmetrize
  USE paw_add_symmetry,      ONLY : paw_deqsymmetrize
  USE control_ph,            ONLY : tr2_ph, convt, lnoloc, lgamma_gamma, zeu, &
                                    niter_ph
  USE control_lr,            ONLY : alpha_pv, nbnd_occ, lgamma
  USE lrus,                  ONLY : int3, int3_paw
  USE dv_of_drho_lr,         ONLY : dv_of_drho
  USE lr_global,             ONLY : rpert, evc0, evq0, sevq0, &
                                    d0psi, d0psi2
  USE lr_cg,                 ONLY : evc1, res, pres, dir, dir_new, prec_vec
  USE units_lr,              ONLY : lrwfc, iuwfc, lrdwf, iudwf
  USE uspp_init,             ONLY : init_us_2
  USE buffers,               ONLY : get_buffer, save_buffer
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm, ntask_groups
  USE mp,                    ONLY : mp_sum
  USE fft_interfaces,        ONLY : fft_interpolate


  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: drhoscfs(dffts%nnr, nspin_mag, rpert)

  COMPLEX(DP) , ALLOCATABLE, TARGET :: &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  COMPLEX(DP) , ALLOCATABLE ::         &
                   dbecsum(:,:,:,:),   & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   aux1 (:,:),         &   ! auxiliary space to apply potential
                   tg_dv(:,:),         &   ! task group variables
                   tg_psic(:,:)          
  COMPLEX(DP), ALLOCATABLE ::          &
                   drhoscf (:,:,:)     ! The change of the scf charge, thick
                                       ! mesh
  COMPLEX(DP), POINTER ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  REAL(DP), ALLOCATABLE :: h_dia(:,:), & ! used to build the preconditioning
                           s_dia(:,:)


  INTEGER :: kter, iter, iter0, ipol, jpol, ibnd, ik, ikp, ikk, ikq, is, &
             npw, npwq, incr, v_siz, ig
  ! counters or indices

  REAL(DP) :: weight, thresh, dr2, aa
  ! weight of k points and store alpha_pv
  REAL(DP) :: tcpu, get_clock

  CALL start_clock ('do_cg_e')

  IF (.NOT.lgamma) CALL errore('do_cg_e','Called in the wrong case',1)

  thresh=tr2_ph*rpert*nbnd*nkstot
  dr2=0.0_DP
  convt=.FALSE.
  iter0=0

  ALLOCATE (drhoscf(dfftp%nnr, nspin_mag, rpert))
  ALLOCATE (dvscfin(dfftp%nnr, nspin_mag, rpert))

  IF (doublegrid) THEN
     ALLOCATE (dvscfins(dffts%nnr, nspin_mag, rpert))
  ELSE
     dvscfins => dvscfin
  ENDIF

  ALLOCATE (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, rpert))
  IF (noncolin) ALLOCATE (dbecsum_nc (nhm, nhm, nat, nspin, rpert))

  ALLOCATE (aux1(dffts%nnr,npol))
  ALLOCATE (h_dia(npwx,npol))
  ALLOCATE (s_dia(npwx,npol))

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
  !
  !   The outside loop is over the conjugate gradient steps
  !
  DO kter = 1, niter_ph

!     WRITE(stdout,*) 'kter', kter
     iter = kter + iter0
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
        CALL g2_kin(ikq)
        !
        IF (iter==1) THEN
           !
           ! at the first iteration reads unperturbed wavefuctions psi_k 
           ! in G_space, for all bands. When q /= 0 reads also evq 
           !
           IF (nksq>1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
           evc0(:,:,ik)=evc(:,:)
           IF (okvan) THEN
              CALL calbec (npwq, vkb, evq, becp, nbnd)
              CALL s_psi (npwx, npwq, nbnd, evq, sevq0(1,1,ik))
           ELSE
              sevq0(:,:,ik)=evq(:,:)
           ENDIF
!
!   At the first iteration compute the preconditioning
!
           h_dia=0.0_DP
           s_dia=0.0_DP
           CALL usnldiag( npwq, h_dia, s_dia )

           DO ibnd = 1, nbnd_occ (ikk)
              DO ig = 1, npwq
                 aa=g2kin(ig)+v_of_0+h_dia(ig,1)-et(ibnd,ikk)*s_dia(ig,1)
                 prec_vec(ig,ibnd,ik)= 1.0_DP/ MAX(1.0_DP, aa)
                 IF (noncolin) &
                     prec_vec(ig+npwx,ibnd,ik)= 1.0_DP/ MAX(1.0_DP, aa)
              ENDDO
           ENDDO
        ELSE
           !
           !  At the following iterations only copy the correct wavefunction
           !
           evc(:,:)=evc0(:,:,ik)
        ENDIF
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
              CALL dvpsi_e(ik,ipol)
              !
              ! Orthogonalize dvpsi to valence states: Apply P_c^+ and change
              ! sign.
              !
              CALL orthogonalize(dvpsi, evq0(1,1,ik), ikk, ikq, sevq0(1,1,ik), npwq, &
                                                                      .TRUE.)
              !
              !  save here P_c^+ V_ext u_kv needed to compute the 
              !  susceptibility
              !
              d0psi(:,:,ik,ipol)=dvpsi(:,:)
              !
              ! At the first iteration evc1 is not known and is set to zero
              ! the residual to -P_c^+ V_ext u_kv
              !
              evc1(:,:,ikp,1)=(0.0_DP, 0.0_DP)
              res(:,:,ikp,1)=dvpsi(:,:)
              !
           ELSE
              !
              !  Here we compute (H-eS) psi^1, at the iteration iter-1 
              !  Since ch_psi_all apply also alpha_pv Q, we set alpha_pv
              !  to zero.
              !
              CALL ch_psi_all (npwq, dir(:,:,ikp,1), dir_new(:,:,ikp,1), &
                                    et(:,ikk), ik, nbnd_occ(ikk))
              !
              ! Now apply P_c^+ dV_Hxc. 
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
              CALL orthogonalize(dvpsi, evq0(1,1,ik), ikk, ikq, &
                                     & sevq0(1,1,ik), npwq, .TRUE.)
              !
              !  the orthogonalize routine changes the sign of dvpsi, so here
              !  we subtract.
              !
              dir_new(:,:,ikp,1) = dir_new(:,:,ikp,1) - dvpsi(:,:)
           ENDIF
        ENDDO
     ENDDO
     !
     !  here we do the cojugate-gradient step. At the first iteration
     !  we only precondition the residual vector.
     !
     IF (iter > 1) THEN
        CALL h_pcg_step(convt,thresh,dr2)
     ELSE
        CALL lr_prec(res, pres)
        dir=pres
     ENDIF
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
     drhoscfs(:,:,:)=(0.d0,0.d0)
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
        evc(:,:)=evc0(:,:,ik)
        !
        !  and on perturbations
        !
        DO ipol = 1, rpert
           ikp = ik + nksq * (ipol - 1)  
           !
           !  At self consistence we compute the charge density with
           !  the solution of the linear system, during the iterations
           !  with the vector to which we apply the matrix

           !
           IF (convt) THEN
              dpsi(:,:)=evc1(:,:,ikp,1)
           ELSE
              dpsi(:,:)=dir(:,:,ikp,1)
           ENDIF
           !
           weight=wk(ikk)
           IF (noncolin) THEN
              CALL incdrhoscf_nc(drhoscfs(1,1,ipol), weight, ik,         &
                                 dbecsum_nc(1,1,1,1,ipol), dpsi, 1.0_DP)
           ELSE
              CALL incdrhoscf(drhoscfs(1,current_spin,ipol), weight, ik, &
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
              CALL fft_interpolate (dffts, drhoscfs(:,is,ipol), dfftp, &
                                                drhoscf(:,is,ipol) )
           ENDDO
        ENDDO
     ELSE
        CALL ZCOPY (nspin_mag*dfftp%nnr*rpert, drhoscfs, 1, drhoscf, 1)
     ENDIF
     !
     !  Rotate dbecsum_nc in the spin-orbit case
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, rpert)
     !
     !  And add the augmentation part of the induced charge    
     !
     CALL addusddense (drhoscf, dbecsum)
     !
     !  Collect the contribution of all pools. At self-consistence 
     !  the uncollected charge is needed
     !
     CALL mp_sum ( drhoscf, inter_pool_comm )
     IF (okpaw) CALL mp_sum ( dbecsum, inter_pool_comm )
     !
     !   dvscfin contains the (unsymmetrized) linear charge response
     !   for the rpert perturbations (3 in the lgamma case and one in the
     !   q /= 0 case) - symmetrize it
     !
     IF (.NOT.lgamma_gamma) THEN
        CALL psymdvscf (drhoscf)
        IF ( noncolin.and.domag ) CALL psym_dmage(drhoscf)
     ENDIF
     !
     !   calculate the corresponding linear potential response
     !
     dvscfin = drhoscf
     DO ipol=1,rpert
        CALL dv_of_drho (dvscfin (1, 1, ipol))
     ENDDO
     !
     !  And interpolate the potential on the smooth grid if needed
     !
     IF (doublegrid) THEN
        DO ipol = 1, rpert
           DO is=1,nspin_mag
              CALL fft_interpolate (dfftp, dvscfin(:,is,ipol),dffts, &
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
           IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
           CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,rpert)
        ENDIF
     ENDIF
!
!   In the US and PAW case computes the integral of dV_Hxc and the 
!   augmentation function. This quantity is needed in adddvscf.
!
     CALL newdq(dvscfin,rpert)

     tcpu = get_clock ('PHONON')
     WRITE( stdout, '(/,5x," iter # ",i8," total cpu time :",f8.1," secs ")') &
                                                                    iter, tcpu
     WRITE( stdout, "(5x,' thresh=',es10.3,15x,' |ddv_scf|^2 = ',es10.3 )") &
                                               thresh, dr2

     FLUSH( stdout )
     !
     IF (convt) EXIT
     !
  ENDDO  ! cg iterations
!
!   now save the solution on disk
!
  IF (convt) THEN
     DO ik=1, nksq
        DO ipol = 1, rpert
           ikp = ik + nksq * (ipol - 1)  
           CALL save_buffer(evc1(1,1,ikp,1), lrdwf, iudwf, ikp)
        ENDDO
     ENDDO
     IF (zeu) THEN
        CALL addnlcc_zstar_eu_us_tpw( drhoscf )
        CALL zstar_eu_us_tpw (dvscfin)
     ENDIF
  ENDIF

  DEALLOCATE (h_dia)
  DEALLOCATE (s_dia)
  DEALLOCATE (aux1)
  DEALLOCATE (dbecsum)
  DEALLOCATE (dvscfin)
  DEALLOCATE (drhoscf)
  IF (doublegrid) DEALLOCATE (dvscfins)
  IF (noncolin) DEALLOCATE (dbecsum_nc)
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE( tg_dv )
     DEALLOCATE( tg_psic )
  ENDIF

  CALL stop_clock ('do_cg_e')

  RETURN
END SUBROUTINE do_cg_e
