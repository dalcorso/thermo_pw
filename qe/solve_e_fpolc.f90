!
! Copyright (C) 2016 Andrea Dal Corso 
! This routine has been obtained by modifyng the routine solve_e_fpol
! of the QE distribution to deal with the calculations at w and -w. 
!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_e_fpolc(iu)
  !-----------------------------------------------------------------------
  !
  !    This routine generalizes to finite complex frequencies the routine
  !    solve_e of the Quantum ESPRESSO distribution. 
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field.
  !    It performs the following tasks:
  !     a) computes the bare potential term  x | psi >
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : diropn
  USE klist,                 ONLY : lgauss, xk, wk
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                 ONLY : nbnd, npwx, g2kin,  et
  USE klist,                 ONLY : ngk, igk_k
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer, save_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag, domag
  USE scf,                   ONLY : rho, v_of_0
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_symmetry,          ONLY : paw_desymmetrize
  USE eqv,                   ONLY : dpsi, dvpsi
  USE units_ph,              ONLY : lrdrho, iudrho
  USE units_lr,              ONLY : lrdwf, iudwf, lrwfc, iuwfc
  USE output,                ONLY : fildrho
  USE control_flags,         ONLY : use_gpu
  USE control_ph,            ONLY : ext_recover, rec_code, &
                                    lnoloc, convt, tr2_ph, nmix_ph, &
                                    alpha_mix, lgamma_gamma, niter_ph, &
                                    flmixdpot, rec_code_read
  USE control_lr,            ONLY : alpha_pv, nbnd_occ, lgamma
  USE lrus,                  ONLY : int3_paw
  USE qpoint,                ONLY : nksq
  USE recover_mod,           ONLY : read_rec, write_rec
  USE optical,               ONLY : current_w, fru, iu1dwf, lr1dwf
  USE freq_ph,               ONLY : fiu
  USE linear_solvers,        ONLY : ccg_many_vectors
  USE dv_of_drho_lr,         ONLY : dv_of_drho
  USE mp_asyn,               ONLY : asyn_master, with_asyn_images
  USE mp_images,             ONLY : my_image_id, root_image
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE fft_helper_subroutines, ONLY : fftx_ntgrp
  USE fft_interfaces,         ONLY : fft_interpolate
  USE uspp_init,            ONLY : init_us_2
  USE apply_dpot_mod,        ONLY : apply_dpot_bands

  implicit none

  INTEGER, INTENT(IN) :: iu
  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  complex(DP), allocatable :: h_diag (:,:)
  complex(DP), allocatable :: h_diag1 (:,:)
  real(DP), allocatable :: h_diagr (:,:)
  real(DP), allocatable :: h_dia (:,:), s_dia(:,:)
  ! h_diag: diagonal part of the Hamiltonian

  complex(DP) , allocatable, target ::      &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  complex(DP) , pointer ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  complex(DP) , allocatable ::   &
                   dpsi1(:,:),   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   mixin(:), mixout(:), &  ! auxiliary for paw mixing
                   aux1 (:,:),  ps (:,:), &
                   tg_dv(:,:), &
                   tg_psic(:,:), aux2(:,:), dvpsi1(:,:)

  logical :: conv_root, exst, all_done_asyn
  ! conv_root: true if linear system is converged

  integer :: kter, iter0, ipol, ibnd, iter, lter, ik, ig, is, nrec, ndim, ios, &
             nmix_ph_eff, nnr, nnrs
  ! counters
  integer :: ltaver, lintercall, incr, jpol, v_siz, npw, npwq

  real(DP) :: tcpu, get_clock
  ! timing variables

  COMPLEX(DP) :: w, aa  !frequency
  REAL(DP) :: aar, weight, alpha_pv0
  LOGICAL :: ldpsi1

  EXTERNAL ch_psi_all, cg_psi
  EXTERNAL apply_ac, ccg_psi_tpw, scal_prod
  COMPLEX(DP) :: scal_prod

  call start_clock ('solve_e')
  !
  !  This routine is task group aware
  !
  w=CMPLX(fru(iu),fiu(iu))
  ldpsi1=ABS(w)>1.D-7
  alpha_pv0=alpha_pv
  alpha_pv=alpha_pv0 + REAL(w)
!
!  the default nmix_ph is slightly too small for this routine.
!
  nmix_ph_eff=max(nmix_ph,8)

  allocate (dvscfin( dfftp%nnr, nspin_mag, 3))
  nnr=dfftp%nnr
  if (doublegrid) then
     allocate (dvscfins(dffts%nnr, nspin_mag, 3))
     nnrs=dffts%nnr
  else
     dvscfins => dvscfin
     nnrs=nnr
  endif
  !$acc enter data create(dvscfins(1:nnrs, 1:nspin_mag, 1:3))
  allocate (dvscfout(dfftp%nnr, nspin_mag, 3))
  IF (okpaw) THEN
     ALLOCATE (mixin(dfftp%nnr*nspin_mag*3+(nhm*(nhm+1)*nat*nspin_mag*3)/2) )
     ALLOCATE (mixout(dfftp%nnr*nspin_mag*3+(nhm*(nhm+1)*nat*nspin_mag*3)/2) )
  ENDIF
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 3))
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin, 3))
  IF (ldpsi1) THEN
     allocate (dpsi1(npwx*npol, nbnd))
     allocate (dvpsi1(npwx*npol, nbnd))
     allocate (h_diag(npwx*npol, nbnd))
     allocate (h_diag1(npwx*npol, nbnd))
     allocate (h_dia(npwx,npol))
     allocate (s_dia(npwx,npol))
  ELSE
     allocate (h_diagr(npwx*npol, nbnd))
  ENDIF
  allocate (aux1(dffts%nnr,npol))
  allocate (aux2(npwx*npol, nbnd))
  !$acc enter data create (aux2(1:npwx*npol, 1:nbnd) )
  IF (okpaw) mixin=(0.0_DP,0.0_DP)


  if (rec_code_read == -20.AND.ext_recover) then
     ! restarting in Electric field calculation
     IF (okpaw) THEN
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins, dvscfout, dbecsum)
        CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixin, dvscfin, dbecsum, ndim, -1 )
     ELSE
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins)
     ENDIF
  else if (rec_code_read > -20 .AND. rec_code_read <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
  endif
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
  IF ( ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
  end if
  IF (rec_code_read > -20) convt=.TRUE.
  !
  if (convt) go to 155
  !
  if (lgauss) call errore ('solve_e_fpolc', 'For metals a q is needed', 1)
  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph

!     write(6,*) 'kter', kter
     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)
     IF (noncolin) dbecsum_nc=(0.d0,0.d0)

     do ik = 1, nksq
        if (lsda) current_spin = isk (ik)
        npw=ngk(ik)
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ik)
        npwq = npw
        call init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb, use_gpu)
        !$acc update host(vkb)
        !
        ! compute the kinetic energy
        !
        CALL g2_kin(ik)
        IF (ldpsi1) THEN
           h_diag=(0.0_DP,0.0_DP)
           h_diag1=(0.0_DP,0.0_DP)
           h_dia=0.0_DP
           s_dia=0.0_DP
           CALL usnldiag( npw, h_dia, s_dia )
           DO ibnd = 1, nbnd_occ (ik)
              !
              !  we precondition with the inverse of the diagonal elements 
              !  of the matrix
              !
              DO ig = 1, npw
                 aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                    (et(ibnd,ik)+w)*s_dia(ig,1) 
                 IF (ABS(aa)<1._DP) aa=(1.0_DP,0.0_DP)
                 h_diag(ig,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa 
                 aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                    (et(ibnd,ik)-w)*s_dia(ig,1) 
                 IF (ABS(aa)<1._DP) aa=(1.0_DP,0.0_DP)
                 h_diag1(ig,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa 
              END DO
              IF (noncolin) THEN
                 DO ig = 1, npw
                    aa=g2kin(ig)+v_of_0+h_dia(ig,2)- &
                       (et(ibnd,ik)+w)*s_dia(ig,2) 
                    IF (ABS(aa)<1._DP) aa=(1.0_DP,0.0_DP)
                    h_diag(ig+npwx,ibnd)=CMPLX(1.d0, 0.d0,kind=DP)/aa 
                    aa=g2kin(ig)+v_of_0+h_dia(ig,2)- &
                       (et(ibnd,ik)-w)*s_dia(ig,2) 
                    IF (ABS(aa)<1._DP) aa=(1.0_DP,0.0_DP)
                    h_diag1(ig+npwx,ibnd)=CMPLX(1.d0, 0.d0,kind=DP)/aa 
                 END DO
              ENDIF
           !
           enddo
        ELSE
           CALL h_prec (ik, evc, h_diagr)
        ENDIF
        !
        do ipol = 1, 3
           !
           ! computes/reads P_c^+ x psi_kpoint into dvpsi array
           !
           call dvpsi_e (ik, ipol)
           !
           if (iter > 1) then
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
!              IF ( dffts%has_task_groups ) THEN
!                 IF (noncolin) THEN
!                    CALL tg_cgather( dffts, dvscfins(:,1,ipol), &
!                                                                tg_dv(:,1))
!                    IF (domag) THEN
!                       DO jpol=2,4
!                          CALL tg_cgather( dffts, dvscfins(:,jpol,ipol), &
!                                                             tg_dv(:,jpol))
!                       ENDDO
!                    ENDIF
!                 ELSE
!                    CALL tg_cgather( dffts, dvscfins(:,current_spin,ipol), &
!                                                             tg_dv(:,1))
!                 ENDIF
!              ENDIF
!              aux2=(0.0_DP,0.0_DP)
!              do ibnd = 1, nbnd_occ (ik), incr
!                 IF ( dffts%has_task_groups ) THEN
!                    call cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, &
!                                      nbnd_occ (ik) )
!                    call apply_dpot(v_siz, tg_psic, tg_dv, 1)
!                    call cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, &
!                                      nbnd_occ (ik))
!                 ELSE
!                    call cft_wave (ik, evc (1, ibnd), aux1, +1)
!                    call apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipol), current_spin)
!                    call cft_wave (ik, aux2 (1, ibnd), aux1, -1)
!                 ENDIF
!              enddo
              CALL apply_dpot_bands(ik, nbnd_occ(ik), &
                                dvscfins(:, :, ipol), evc, aux2)

              dvpsi=dvpsi+aux2
              !
              call adddvscf(ipol,ik)
              !
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
           !
           CALL orthogonalize(dvpsi, evc, ik, ik, dpsi, npwq, .false.)
           !
           !  dvpsi is saved because the ccg_many_vectors routine corrupts it
           !
           IF (ldpsi1) dvpsi1(:,:)=dvpsi(:,:)
           !
           if (iter == 1) then
              !
              !  At the first iteration dpsi and dvscfin are set to zero,
              !
              dpsi(:,:)=(0.d0,0.d0)
              IF (ldpsi1) dpsi1(:,:)=(0.d0,0.d0)
              dvscfin(:,:,:)=(0.d0,0.d0)
              !
              ! starting threshold for the iterative solution of the linear
              ! system
              !
              thresh = 1.d-2
              if (lnoloc) thresh = 1.d-5
           else
              ! starting value for  delta_psi is read from iudwf
              !
              nrec = (ipol - 1) * nksq + ik
              call get_buffer (dpsi, lrdwf, iudwf, nrec)
              IF (ldpsi1) call get_buffer (dpsi1, lr1dwf, iu1dwf, nrec)
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (0.1d0 * sqrt (dr2), 1.0d-2)
           endif
           !
           ! iterative solution of the linear system (H-e-w)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           conv_root = .true.

           current_w=w
           IF (ldpsi1) THEN
!
!  Complex or imaginary frequency. Use bicojugate gradient.
!
              CALL ccg_many_vectors (apply_ac,ccg_psi_tpw,scal_prod,dvpsi, &
                dpsi, h_diag, npwx*npol, npw, thresh, ik, lter, conv_root, &
                                   anorm, nbnd_occ(ik))
           ELSE
!
!    zero frequency. The standard QE solver
!
             CALL errore('solve_e_fpolc','not programmed',1)
             call cgsolve_all (ch_psi_all,cg_psi,et(1,ik),dvpsi,dpsi, &
               h_diagr,npwx,npw,thresh,ik,lter,conv_root,anorm,nbnd_occ(ik), &
                                                                     npol)
           END IF

!           WRITE(6,*) 'w+', lter
           ltaver = ltaver + lter
           lintercall = lintercall + 1
           if (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                &         ' solve_e_fpolc: root not converged ',es10.3)") ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           nrec = (ipol - 1) * nksq + ik
           call save_buffer(dpsi, lrdwf, iudwf, nrec)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight=wk(ik)
           IF ( with_asyn_images.AND.my_image_id==root_image.AND.ionode ) &
                              CALL asyn_master(all_done_asyn)
           IF (ldpsi1) THEN
!
!    complex frequency, two wavefunctions must be computed
!
              weight=wk(ik)/2.0_DP
!
!   In this case compute also the wavefunction at frequency -w.
!
              current_w=-w
              CALL ccg_many_vectors (apply_ac,ccg_psi_tpw,scal_prod,dvpsi1, &
                             dpsi1, h_diag1, npwx*npol, npw, thresh, ik, &
                             lter, conv_root, anorm, nbnd_occ(ik))
!              WRITE(6,*) 'w-', lter
              ltaver = ltaver + lter
              lintercall = lintercall + 1
              if (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                &         ' solve_e: root not converged ',es10.3)") ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
              nrec = (ipol - 1) * nksq + ik
              call save_buffer(dpsi1, lr1dwf, iu1dwf, nrec)
!
!    add the two wavefunctions before computing the contribution to the charge
!
              CALL DAXPY(npwx*nbnd_occ(ik)*npol*2, 1.0_DP, dpsi1, 1, dpsi, 1)
           END IF
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           IF (noncolin) THEN
              call incdrhoscf_nc(dvscfout(1,1,ipol), weight, ik, &
                                 dbecsum_nc(1,1,1,1,ipol), dpsi, 1.0_DP)
           ELSE
              call incdrhoscf (dvscfout(1,current_spin,ipol), weight, &
                         ik, dbecsum(1,1,current_spin,ipol), dpsi)
           END IF
           IF ( with_asyn_images.AND.my_image_id==root_image.AND.ionode ) &
                              CALL asyn_master(all_done_asyn)
        enddo   ! on polarizations
     enddo      ! on k points
     current_w=w

     IF (lnoloc) THEN
        dvscfout=(0.d0,0.d0)
        convt=.TRUE.
        GOTO 1001
     ENDIF
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions
     !  coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        call mp_sum ( dbecsum, intra_bgrp_comm )
     END IF

     if (doublegrid) then
        do is=1,nspin_mag
           do ipol=1,3
              call fft_interpolate (dffts, dvscfout(:,is,ipol), dfftp, &
                                           dvscfout(:,is,ipol))
           enddo
        enddo
     endif
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 3)
     !
     call addusddense (dvscfout, dbecsum)
     !
     !   dvscfout contains the (unsymmetrized) linear charge response
     !   for the three polarizations - symmetrize it
     !
     call mp_sum ( dvscfout, inter_pool_comm )
     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
     if (.not.lgamma_gamma) then
        IF (ldpsi1) THEN
           CALL symmetrize_drho(dvscfout, dbecsum, 0, 3, 4)
        ELSE
           CALL symmetrize_drho(dvscfout, dbecsum, 0, 3, 2)
        ENDIF
     endif
     IF ( with_asyn_images.AND.my_image_id==root_image.AND.ionode ) &
                              CALL asyn_master(all_done_asyn)
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     do ipol=1,3
        if (fildrho.ne.' ') call davcio_drho(dvscfout(1,1,ipol),lrdrho, &
             iudrho,ipol,+1)
        call dv_of_drho (dvscfout (1, 1, ipol))
     enddo
     !
     !   mix the new potential with the old
     !
     IF (okpaw) THEN
     !
     !  In this case we mix also dbecsum
     !
        call setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixout, dvscfout, dbecsum, ndim, -1 )
        call mix_potential_tpw (2*3*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, 3*tr2_ph/npol, iter, &
                          flmixdpot, convt)
        call setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )
     ELSE
        call mix_potential_tpw (2*3*dfftp%nnr*nspin_mag, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2, 3 * tr2_ph / npol, iter, flmixdpot, convt)
     ENDIF
 
     if (doublegrid) then
        do is=1,nspin_mag
           do ipol = 1, 3
              call fft_interpolate (dfftp, dvscfin(:,is,ipol), &
                                    dffts, dvscfins(:,is,ipol))
           enddo
        enddo
     endif

     IF (okpaw) THEN
        IF (noncolin.AND.domag) THEN
!           call PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
        ELSE
!
!    The presence of c.c. in the formula gives a factor 2.0
!
           dbecsum=2.0_DP * dbecsum
           IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
           call PAW_dpotential(dbecsum,rho%bec,int3_paw,3)
        ENDIF
     ENDIF

     call newdq(dvscfin,3)

1001 CONTINUE

     CALL mp_sum(ltaver,inter_pool_comm)
     CALL mp_sum(lintercall,inter_pool_comm)
     averlt = DBLE (ltaver) / DBLE (lintercall)

     tcpu = get_clock ('PHONON')
     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / 3
     WRITE( stdout, "(5x,' thresh=',es10.3, ' alpha_mix = ',f6.3, &
          &      ' |ddv_scf|^2 = ',es10.3 )") thresh, alpha_mix (kter), dr2
     !
     FLUSH( stdout )
     !
     ! rec_code: state of the calculation
     ! rec_code=-20 Electric Field
     !
     rec_code=-20
     IF (okpaw) THEN
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 3, dvscfin, &
                                                       dvscfout, dbecsum)
     ELSE
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 3, dvscfin)
     ENDIF

     if (check_stop_now()) call stop_smoothly_ph (.false.)

     if (convt) goto 155

  enddo
155 continue

  deallocate (aux1)
  IF (ldpsi1) THEN
     DEALLOCATE (dpsi1)
     DEALLOCATE (dvpsi1)
     DEALLOCATE (h_diag)
     DEALLOCATE (h_diag1)
     DEALLOCATE(h_dia)
     DEALLOCATE(s_dia)
  ELSE
     DEALLOCATE (h_diagr)
  ENDIF
  DEALLOCATE (dbecsum)
  DEALLOCATE (dvscfout)

  IF (okpaw) THEN
     DEALLOCATE(mixin)
     DEALLOCATE(mixout)
  ENDIF
  !$acc exit data delete(dvscfins)
  IF (doublegrid) DEALLOCATE (dvscfins)
  DEALLOCATE (dvscfin)
  IF (noncolin) DEALLOCATE(dbecsum_nc)
  !$acc exit data delete(aux2)
  DEALLOCATE(aux2)
  IF ( dffts%has_task_groups ) THEN
     !
     DEALLOCATE( tg_dv  )
     DEALLOCATE( tg_psic)
     !
  ENDIF

  alpha_pv=alpha_pv0

  CALL stop_clock ('solve_e')
  RETURN
END SUBROUTINE solve_e_fpolc
