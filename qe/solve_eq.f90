!
! Copyright (C) 2016 Andrea Dal Corso 
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_eq(iu, flag)
  !-----------------------------------------------------------------------
  !
  !    This routine generalizes to finite complex frequencies and 
  !    finite q vectors the routine solve_e of the Quantum ESPRESSO 
  !    distribution. 
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field 
  !    of finite wavevector q and complex frequency omega.
  !    It performs the following tasks:
  !     a) computes the bare potential term  e^{iqr} | psi >
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system at zero
  !        frequency or ccg_many_vectors
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !
  USE kinds,                 ONLY : DP
  USE constants,             ONLY : e2, fpi, rytoev
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : prefix, diropn
  USE cell_base,             ONLY : tpiba2
  USE fft_interfaces,        ONLY : fwfft
  USE fft_interfaces,        ONLY : fft_interpolate
  USE klist,                 ONLY : lgauss, xk, wk
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                 ONLY : nbnd, npwx, g2kin,  et
  USE klist,                 ONLY : ngk, igk_k
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer, save_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : okvan, nkb, vkb
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag, domag
  USE scf,                   ONLY : rho, v_of_0
  USE gvect,                 ONLY : gg
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_add_symmetry,      ONLY : paw_deqsymmetrize

  USE eqv,                   ONLY : dpsi, dvpsi, evq
  USE units_lr,              ONLY : lrdwf, iudwf, lrwfc, iuwfc
  USE units_ph,              ONLY : lrdrho, iudrho, lrbar, iubar
  USE output,                ONLY : fildrho
  USE control_flags,         ONLY : use_gpu
  USE control_ph,            ONLY : ext_recover
  USE control_lr,            ONLY : lgamma, alpha_pv, nbnd_occ, rec_code, &
                                    alpha_mix, lgamma_gamma, niter_ph, &
                                    flmixdpot, rec_code_read, convt, tr2_ph, &
                                    rec_code, lnoloc
  USE dfpt_type,             ONLY : dfpt_data_type, allocate_dfpt_data,   &
                                    deallocate_dfpt_data, dfpt_dvscfp_to_dvscfs
  USE lrus,                  ONLY : int3_paw
  USE qpoint,                ONLY : xq, nksq, ikks, ikqs
  USE recover_mod,           ONLY : read_rec, write_rec

  USE optical,               ONLY : current_w, fru, iu1dwf, lr1dwf, chirr, &
                                    chirz, chizr, chizz, epsm1
  USE freq_ph,               ONLY : fiu
  USE linear_solvers,        ONLY : ccg_many_vectors
  USE dv_of_drho_lr,         ONLY : dv_of_drho
  USE dv_of_drho_clf,        ONLY : dv_of_drho_nlf
  USE mp_asyn,               ONLY : asyn_master, with_asyn_images
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp_images,             ONLY : root_image, my_image_id
  USE mp,                    ONLY : mp_sum
  USE fft_helper_subroutines, ONLY : fftx_ntgrp
  USE uspp_init,            ONLY : init_us_2
  USE apply_dpot_mod,       ONLY : apply_dpot_bands
  USE incdrhoscf_mod,       ONLY : incdrhoscf, incdrhoscf_nc

  implicit none

  INTEGER, INTENT(IN) :: iu
  INTEGER, INTENT(IN) :: flag   ! if 1 compute the charge-charge and
                                ! charge magnetization responses
                                ! if 2 and lsda computes the magnetization
                                ! magnetization response
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

  TYPE(dfpt_data_type) :: dfpt_data
  !! Data that describes linear response quantities
 
  complex(DP) , allocatable ::   &
                   dpsi1(:,:),   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   mixin(:), mixout(:), &  ! auxiliary for paw mixing
                   aux1 (:,:),  ps (:,:), &
                   tg_dv(:,:), &
                   tg_psic(:,:), aux2(:,:), dvpsi1(:,:)

  complex(DP), EXTERNAL :: zdotc      ! the scalar product function

  logical :: conv_root, exst, all_done_asyn, time_reversed
  ! conv_root: true if linear system is converged

  integer :: kter, iter0, ipol, ibnd, iter, lter, ik, ikk, ikq, &
             ig, is, nrec, ndim, npw, npwq, nnr, nnrs, ios
  ! counters
  integer :: ltaver, lintercall, incr, jpol, v_siz
  real(DP) :: xqmod2, alpha_pv0

  real(DP) :: tcpu, get_clock
  ! timing variables

  COMPLEX(DP) :: w  !frequency
  REAL(DP) :: aa, weight
  LOGICAL :: ldpsi1

  external ch_psi_all, cg_psi
  external apply_ac, ccg_psi_tpw, scal_prod
  COMPLEX(DP) :: scal_prod


  call start_clock ('solve_eq')
  !
  !  This routine is task group aware
  !
  w=CMPLX(fru(iu),fiu(iu))
  ldpsi1=ABS(w)>1.D-7
  alpha_pv0=alpha_pv
  alpha_pv=alpha_pv0 + REAL(w)


  nnr=dfftp%nnr
  nnrs=dffts%nnr

  CALL allocate_dfpt_data(dfpt_data, 1)
  !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:nnrs, 1:nspin_mag, 1:1))
  allocate (dvscfout(dfftp%nnr, nspin_mag, 1))
  IF (okpaw) THEN
     ALLOCATE (mixin(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
     ALLOCATE (mixout(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
  ENDIF
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin, 1))
  IF (ldpsi1) THEN
     allocate (dpsi1(npwx*npol,nbnd))
     allocate (dvpsi1(npwx*npol,nbnd))
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
        CALL read_rec(dr2, iter0, dfpt_data)
        CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixin, dfpt_data%dvscfp, dfpt_data%dbecsum, ndim, -1 )
     ELSE
        CALL read_rec(dr2, iter0, dfpt_data)
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
  if ((lgauss.and..not.ldpsi1)) &
          call errore ('solve_eq', 'insert a finite frequency', 1)
  !
  !   The outside loop is over the iterations
  !
  time_reversed=.FALSE.
  do kter = 1, niter_ph

!     write(6,*) 'kter', kter
     FLUSH( stdout )
     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dfpt_data%drhos(:,:,:)=(0.d0,0.d0)
     dfpt_data%dbecsum(:,:,:,:)=(0.d0,0.d0)
     IF (noncolin) dbecsum_nc=(0.d0,0.d0)

     do ik = 1, nksq
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq = ngk(ikq)
        if (lsda) current_spin = isk (ikk)
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        IF (nkb>0) THEN
           call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb, use_gpu)
           !$acc update host(vkb)
        ENDIF
        !
        ! reads unperturbed wavefuctions psi(k) and psi(k+q)
        !
        if (nksq.gt.1) then
           call get_buffer (evc, lrwfc, iuwfc, ikk)
           call get_buffer (evq, lrwfc, iuwfc, ikq)
        endif
        !
        ! compute the kinetic energy
        !
        CALL g2_kin(ikq)

        IF (ldpsi1) THEN
           h_diag=(0.0_DP,0.0_DP)
           h_diag1=(0.0_DP,0.0_DP)
           h_dia=0.0_DP
           s_dia=0.0_DP
           CALL usnldiag( npwq, npol, h_dia, s_dia )

           DO ibnd = 1, nbnd_occ (ikk)
              !
              DO ig = 1, npwq
                 aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                    (et(ibnd,ikk)+w)*s_dia(ig,1) 
                 IF (ABS(aa)<1.0_DP) aa=1.0_DP
                 h_diag(ig,ibnd)=CMPLX(1.0d0, 0.d0,kind=DP) / aa
                 aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                    (et(ibnd,ikk)-w)*s_dia(ig,1) 
                 IF (ABS(aa)<1.0_DP) aa=1.0_DP
                 h_diag1(ig,ibnd)=CMPLX(1.0d0, 0.d0,kind=DP) / aa
              END DO
              !
              IF (noncolin) THEN
                 do ig = 1, npwq
                    aa=g2kin(ig)+v_of_0+h_dia(ig,2)- &
                       (et(ibnd,ikk)+w)*s_dia(ig,2)
                    IF (ABS(aa)<1.0_DP) aa=1.0_DP
                    h_diag(ig+npwx,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa
                    aa=g2kin(ig)+v_of_0+h_dia(ig,2)- &
                       (et(ibnd,ikk)-w)*s_dia(ig,2)
                    IF (ABS(aa)<1.0_DP) aa=1.0_DP
                    h_diag1(ig+npwx,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa
                 enddo
              END IF
            enddo
         ELSE
            CALL h_prec (ik, evc, h_diagr)
            do ibnd = 1, nbnd_occ (ikk)
               !
               do ig = 1, npwq
                  aa=1.0_DP / h_diagr(ig,ibnd)-et(ibnd,ikk)-REAL(w,KIND=DP)
                  h_diagr(ig,ibnd)=1.d0 /max(1.0d0,aa)
               end do
               IF (noncolin) THEN
                  do ig = 1, npwq
                     h_diagr(ig+npwx,ibnd)= h_diagr(ig,ibnd)
                  enddo
               END IF
            enddo
         ENDIF
           !
        !
        do ipol = 1, 1
           nrec = (ipol - 1) * nksq + ik
           !
           if (iter > 1) then
              !
              ! After the first iteration dvbare_q*psi_kpoint is read from file
              !
              call get_buffer (dvpsi, lrbar, iubar, nrec)
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
              CALL apply_dpot_bands(ik, nbnd_occ(ikk), &
                                dfpt_data%dvscfs(:, :, ipol), evc, aux2)
              dvpsi=dvpsi+aux2
              !
              call adddvscf(ipol,ik,time_reversed)
              !
           else
              !
              !  At the first iteration dvbare_q*psi_kpoint is calculated
              !  and written to file
              !
              CALL dveqpsi_us (ik)
              !
              !  with flag=2 the perturbation is a magnetic field along z
              !
              IF (lsda.AND.current_spin==2.AND.flag==2) dvpsi=-dvpsi
              call save_buffer (dvpsi, lrbar, iubar, nrec)
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
           !
           IF (ldpsi1) THEN
              dvpsi1=dvpsi
              CALL orthogonalize_omega(dvpsi1, evq, ikk, ikq, dpsi, npwq, -w)
           ENDIF
           CALL orthogonalize_omega(dvpsi, evq, ikk, ikq, dpsi, npwq, w)
           !
           if (iter == 1) then
              !
              !  At the first iteration dpsi and dvscfin are set to zero,
              !

              dpsi(:,:)=(0.d0,0.d0)
              IF (ldpsi1) dpsi1(:,:)=(0.d0,0.d0)
              dfpt_data%dvscfp(:,:,:)=(0.d0,0.d0)
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
              thresh = min (0.1d0 * sqrt (dr2), thresh)
           endif
           !
           ! iterative solution of the linear system (H-e)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           conv_root = .true.

           current_w=w
           IF (ldpsi1) THEN
!
!  Complex or imaginary frequency. Use bicojugate gradient.
!
              CALL ccg_many_vectors (apply_ac,ccg_psi_tpw,scal_prod,dvpsi, &
                dpsi, h_diag, npwx*npol, npwq, thresh, ik, lter, conv_root, &
                                   anorm, nbnd_occ(ikk))
           ELSE
!
!    zero frequency. The standard QE solver
!
              CALL cgsolve_all (ch_psi_all,cg_psi,et(1,ikk),dvpsi,dpsi, &
               h_diagr,npwx,npwq,thresh,ik,lter,conv_root,anorm,&
                                                          nbnd_occ(ikk),npol)
           END IF

           ltaver = ltaver + lter
           lintercall = lintercall + 1
           if (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                &         ' solve_e: root not converged ',es10.3)") ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           nrec = (ipol - 1) * nksq + ik
           call save_buffer(dpsi, lrdwf, iudwf, nrec)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight=wk(ikk)
           IF (ldpsi1) THEN
!
!    complex frequency, two wavefunctions must be computed
!
              weight=wk(ikk)/2.0_DP
!
!   In this case compute also the wavefunction at frequency -w.
!
              current_w=-w
              CALL ccg_many_vectors (apply_ac,ccg_psi_tpw,scal_prod,dvpsi1, &
                             dpsi1, h_diag1, npwx*npol, npwq, thresh, ik, &
                             lter, conv_root, anorm, nbnd_occ(ikk))
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
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
              CALL DAXPY(npwx*nbnd_occ(ikk)*npol*2, 1.0_DP, dpsi1, 1, dpsi, 1)
           END IF
           IF (noncolin) THEN
              CALL incdrhoscf_nc(dfpt_data%drhos(1,1,ipol), weight, ik, &
                                 dbecsum_nc(1,1,1,1,ipol), dpsi, 1.0_DP)
           ELSE
              CALL incdrhoscf (dfpt_data%drhos(1,current_spin,ipol), weight, &
                         ik, dfpt_data%dbecsum(1,1,current_spin,ipol), dpsi)
           END IF
        ENDDO   ! on polarizations
        IF ( with_asyn_images.AND.my_image_id==root_image.AND.ionode ) &
                           CALL asyn_master(all_done_asyn)
     ENDDO      ! on k points
     current_w=w
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions
     !  coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        call mp_sum ( dfpt_data%dbecsum, intra_bgrp_comm )
     END IF

     if (doublegrid) then
        do is=1,nspin_mag
           call fft_interpolate (dffts, dfpt_data%drhos(:,is,1), dfftp, &
                                                      dvscfout(:,is,1))
        enddo
     else
       CALL zcopy (nspin_mag*dfftp%nnr, dfpt_data%drhos, 1, &
                                                      dvscfout, 1)
     endif
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, &
                                              dfpt_data%dbecsum, 1)
     !
     call addusddenseq (dvscfout, dfpt_data%dbecsum)
     !
     !   dvscfout contains the (unsymmetrized) linear charge response
     !   for the three polarizations - symmetrize it
     !
     call mp_sum ( dvscfout, inter_pool_comm )
     
     IF (okpaw) call mp_sum ( dfpt_data%dbecsum, inter_pool_comm )

     CALL symmetrize_drho(dvscfout, dfpt_data%dbecsum, 0, 3, 3)

     dfpt_data%drhop(:,:,1)=dvscfout(:,:,1)
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     IF (lnoloc) THEN
        CALL dv_of_drho_nlf (dvscfout (1, 1, 1))
     ELSE
        CALL dv_of_drho (dvscfout (1, 1, 1))
     ENDIF
     !
     !   mix the new potential with the old
     !
     IF (okpaw) THEN
     !
     !  In this case we mix also dbecsum
     !
        call setmixout(dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag)/2, &
                    mixout, dvscfout, dfpt_data%dbecsum, ndim, -1 )
        CALL mix_potential_tpw (2*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, tr2_ph/npol, iter, flmixdpot, &
                         convt)
        call setmixout(dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag)/2, &
                       mixin, dfpt_data%dvscfp, dfpt_data%dbecsum, ndim, 1 )
     ELSE
        CALL mix_potential_tpw(2*dfftp%nnr*nspin_mag, dvscfout, &
                                             dfpt_data%dvscfp,  &
             alpha_mix (kter), dr2,  tr2_ph / npol, iter, flmixdpot, convt)
     ENDIF

     CALL dfpt_dvscfp_to_dvscfs(dfpt_data)

     IF (okpaw) THEN
        IF (noncolin.AND.domag) THEN
!           call PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
        ELSE
!
!    The presence of c.c. in the formula gives a factor 2.0
!
           dfpt_data%dbecsum=2.0_DP * dfpt_data%dbecsum
           IF (.NOT. lgamma_gamma) CALL PAW_deqsymmetrize(dfpt_data%dbecsum)
           call PAW_dpotential(dfpt_data%dbecsum,rho%bec,int3_paw,1)
        ENDIF
     ENDIF

     call newdq(dfpt_data%dvscfp,1)

1001 CONTINUE

     CALL mp_sum(ltaver,inter_pool_comm)
     CALL mp_sum(lintercall,inter_pool_comm)

     averlt = DBLE (ltaver) / DBLE (lintercall)

     tcpu = get_clock ('PHONON')
     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
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
        CALL write_rec('solve_e...', 0, dr2, iter, convt, dfpt_data)
     ELSE
        CALL write_rec('solve_e...', 0, dr2, iter, convt, dfpt_data)
     ENDIF

     if (check_stop_now()) call stop_smoothly_ph (.false.)

     if (convt) goto 155

  enddo
155 continue
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

  deallocate (aux1)
  IF (ldpsi1) THEN
     deallocate (dpsi1)
     deallocate (dvpsi1)
     deallocate (h_diag)
     deallocate (h_diag1)
     deallocate (h_dia)
     deallocate (s_dia)
  ELSE
     deallocate (h_diagr)
  ENDIF
  deallocate (dvscfout)
  IF (okpaw) THEN
     DEALLOCATE(mixin)
     DEALLOCATE(mixout)
  ENDIF
  if (noncolin) deallocate(dbecsum_nc)
  !$acc exit data delete(aux2)
  deallocate(aux2)
  IF ( dffts%has_task_groups ) THEN
     !
     DEALLOCATE( tg_dv  )
     DEALLOCATE( tg_psic)
     !
  ENDIF
  alpha_pv=alpha_pv0

  !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)
  CALL deallocate_dfpt_data(dfpt_data)


  call stop_clock ('solve_eq')
  return
end subroutine solve_eq
