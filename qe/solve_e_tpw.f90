!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e_tpw(drhoscf)
  !-----------------------------------------------------------------------
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field.
  !    It performs the following tasks:
  !     a) computes the bare potential term  x | psi >
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !        If lda_plus_u=.true. compute also the SCF part
  !        of the response Hubbard potential.
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !     f) If lda_plus_u=.true. compute also the response occupation
  !        matrices dnsscf
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : prefix, diropn
  USE klist,                 ONLY : ltetra, lgauss, xk, wk, ngk, igk_k
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                 ONLY : nbnd, npwx, g2kin, et
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer, save_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : okvan, nkb, vkb, nlcc_any
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag, domag
  USE scf,                   ONLY : rho
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_symmetry,          ONLY : paw_desymmetrize

  USE units_ph,              ONLY : lrdrho,  iudrho
  USE units_lr,              ONLY : lrdwf, iudwf, lrwfc, iuwfc
  USE output,                ONLY : fildrho
  USE control_flags,         ONLY : use_gpu
  USE control_ph,            ONLY : ext_recover, rec_code, &
                                    lnoloc, convt, tr2_ph, nmix_ph, zeu, &
                                    alpha_mix, lgamma_gamma, niter_ph, &
                                    flmixdpot, rec_code_read
  USE recover_mod,           ONLY : read_rec, write_rec

  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE fft_helper_subroutines, ONLY : fftx_ntgrp

  USE lrus,                  ONLY : int3_paw
  USE qpoint,                ONLY : nksq, ikks, ikqs
  USE eqv,                   ONLY : dpsi, dvpsi
  USE control_lr,            ONLY : nbnd_occ, lgamma
  USE dv_of_drho_lr
  USE fft_interfaces,        ONLY : fft_interpolate, fwfft
  USE ldaU,                  ONLY : lda_plus_u
  USE magnetic_charges,      ONLY : alpha_me
  USE cell_base,             ONLY : bg
  USE gvect,                 ONLY : gg
  USE uspp_init,             ONLY : init_us_2
  USE apply_dpot_mod,        ONLY : apply_dpot_bands

  implicit none

  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP), allocatable :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian

  complex(DP) , allocatable, target ::      &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  complex(DP) , pointer ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  complex(DP) :: drhoscf (dffts%nnr, nspin_mag, 3)
  complex(DP), allocatable:: drhoscf_aux(:,:,:), alpha_work(:,:)
  complex(DP), pointer ::  drhoscfh (:,:,:) ! change of the scf potential (output)
  complex(DP), allocatable, target :: &
                   dvscfout (:,:,:) ! change of the scf potential (output)
  complex(DP) , allocatable ::   &
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   mixin(:), mixout(:), &  ! auxiliary for paw mixing
                   ps (:,:), &
                   tg_dv(:,:), &
                   tg_psic(:,:), aux2(:,:)

  complex(DP), EXTERNAL :: zdotc      ! the scalar product function

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: npw, npwq
  integer :: kter, iter0, ipol, ibnd, iter, lter, ik, ig, is, nrec, ndim, ios
  ! counters
  integer :: ltaver, lintercall, incr, jpol, v_siz, nnr, nnrs
  integer :: ikk, ikq
  integer :: icart, jcart

  real(DP) :: tcpu, get_clock
  ! timing variables

  external ch_psi_all, cg_psi

  call start_clock ('solve_e')
  !
  !  This routine is task group aware
  !
  IF (noncolin.AND.domag) then 
     allocate (drhoscf_aux(dfftp%nnr, nspin_mag, 3))
     allocate (alpha_work(3,3))
     drhoscf_aux = (0.0_DP, 0.0_DP)
     alpha_work = (0.0_DP, 0.0_DP)
  END IF
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
  IF (nlcc_any.AND.zeu) THEN
     ALLOCATE(drhoscfh(dfftp%nnr, nspin_mag, 3))
  ELSE
    drhoscfh => dvscfout
  ENDIF
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 3))
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin, 3))
  allocate (h_diag(npwx*npol, nbnd))
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
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if ( (lgauss .or. ltetra) .or..not.lgamma) call errore ('solve_e', &
       'called in the wrong case', 1)

  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph

!     write(6,*) 'kter', kter
     FLUSH( stdout )
     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     drhoscf(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)
     IF (noncolin) dbecsum_nc=(0.d0,0.d0)
     !
     ! DFPT+U: at each iteration calculate dnsscf,
     ! i.e. the scf variation of the occupation matrix ns.
     !
     IF (lda_plus_u .AND. (iter.NE.1)) &
        CALL dnsq_scf (3, .false., 0, 1, .false.)
     !
     do ik = 1, nksq
        ikk=ikks(ik)
        ikq=ikqs(ik)
        npw = ngk(ikk)
        npwq= npw     ! q=0 always in this routine
        if (lsda) current_spin = isk (ikk)
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ikk)
        !
        ! compute beta functions and kinetic energy for k-point ik
        ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
        !
        IF ( nkb>0 ) THEN
           CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb, use_gpu)
           !$acc update host(vkb)
        ENDIF
        CALL g2_kin(ikk)
        !
        ! compute preconditioning matrix h_diag used by cgsolve_all
        !
        CALL h_prec (ik, evc, h_diag)
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
              CALL apply_dpot_bands(ik, nbnd_occ(ikk), &
                                dvscfins(:, :, ipol), evc, aux2)
              !
              dvpsi=dvpsi+aux2
              !
              call adddvscf(ipol,ik)
              !
              ! DFPT+U: add to dvpsi the scf part of the response
              ! Hubbard potential dV_hub
              !
              if (lda_plus_u) call adddvhubscf (ipol, ik)
              !
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
           !
           CALL orthogonalize(dvpsi, evc, ikk, ikk, dpsi, npwq, .false.)
           !
           if (iter == 1) then
              !
              !  At the first iteration dpsi and dvscfin are set to zero,
              !
              dpsi(:,:)=(0.d0,0.d0)
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
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (0.1d0 * sqrt (dr2), 1.0d-2)
           endif
           !
           ! iterative solution of the linear system (H-e)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
           !

           conv_root = .true.

           call cgsolve_all (ch_psi_all,cg_psi,et(1,ikk),dvpsi,dpsi, &
              h_diag,npwx,npw,thresh,ik,lter,conv_root,anorm,&
              nbnd_occ(ikk),npol)

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
           IF (noncolin) THEN
              call incdrhoscf_nc(drhoscf(1,1,ipol),wk(ikk),ik, &
                                 dbecsum_nc(1,1,1,1,ipol), dpsi, 1.0_DP)
           ELSE
              call incdrhoscf (drhoscf(1,current_spin,ipol), wk(ikk), &
                            ik, dbecsum(1,1,current_spin,ipol), dpsi)
           ENDIF
        enddo   ! on polarizations
     enddo      ! on k points
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
              call fft_interpolate (dffts, drhoscf(:,is,ipol), &
                                    dfftp, drhoscfh(:,is,ipol))
           enddo
        enddo
     else
        CALL zcopy (nspin_mag*dfftp%nnr*3, drhoscf, 1, drhoscfh, 1)
     endif
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 3)
     !
     call addusddense (drhoscfh, dbecsum)
     !
     !   dvscfout contains the (unsymmetrized) linear charge response
     !   for the three polarizations - symmetrize it
     !
     call mp_sum ( drhoscfh, inter_pool_comm )
     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
     CALL symmetrize_drho(drhoscfh, dbecsum, 0, 3, 2)
!     if (.not.lgamma_gamma) then
!        call psyme (drhoscfh)
!        IF ( noncolin.and.domag ) CALL psym_dmage(drhoscfh)
!     endif
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     do ipol=1,3
        if (fildrho.ne.' ') call davcio_drho(drhoscfh(1,1,ipol),lrdrho, &
             iudrho,ipol,+1)
        IF (lnoloc) then
           dvscfout(:,:,ipol)=(0.d0,0.d0)
        ELSE
          IF (nlcc_any.and.zeu) call zcopy (dfftp%nnr*nspin_mag,&
                         drhoscfh(1,1,ipol),1,dvscfout(1,1,ipol),1)
           call dv_of_drho (dvscfout (1, 1, ipol))
        ENDIF
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
        call mix_potential (2*3*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, 3*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
        call setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )

     ELSE
        call mix_potential (2*3*dfftp%nnr*nspin_mag, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2, 3 * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
     ENDIF
     if (doublegrid) then
        do is=1,nspin_mag
           do ipol = 1, 3
              call fft_interpolate (dfftp, dvscfin(:,is,ipol), dffts, &
                                           dvscfins(:,is,ipol))
           enddo
        enddo
     endif

     IF (okpaw) THEN
!
!    The presence of c.c. in the formula gives a factor 2.0
!
        dbecsum=2.0_DP * dbecsum
        IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
        call PAW_dpotential(dbecsum,rho%bec,int3_paw,3)
     ENDIF

     call newdq(dvscfin,3)

     CALL mp_sum(ltaver, inter_pool_comm)
     CALL mp_sum(lintercall, inter_pool_comm)
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
IF (convt) THEN
!
!    Add the contribution of this change of the self-consistent potential
!    to the effective charges
!
   IF (zeu) THEN
      CALL addnlcc_zstar_eu_us_tpw( drhoscfh )
      CALL zstar_eu_us_tpw (dvscfin)
   ENDIF
ENDIF
!
!   In the noncolinear magnetic case, we compute the integral of the response 
!   of the magnetization to the electric field, which allows 
!   to compute the frozen-ion spin component of the magnetoelectric tensor.
!
IF (noncolin.AND.domag) then 
   alpha_me(:,:) = (0.0_DP, 0.0_DP)
   DO ipol = 1, 3
      drhoscf_aux(:,:,:) = drhoscfh(:,:,:)
      DO is=2,nspin_mag
         CALL fwfft ('Rho', drhoscf_aux(:,is,ipol), dfftp)
         IF (ABS(gg(1)).LT.1.d-8) THEN 
            alpha_me(ipol,is-1)= drhoscf_aux(dfftp%nl(1),is,ipol)
         ENDIF
      ENDDO
   ENDDO
   CALL mp_sum(alpha_me,intra_bgrp_comm)
!
!   and we bring to cartesian coordinates the components 
!   of the magnetoelectric tensor.
!
   DO icart=1,3
      alpha_work(:,icart) = alpha_me(1,icart) * bg(:,1) + &
                            alpha_me(2,icart) * bg(:,2) + &
                            alpha_me(3,icart) * bg(:,3)
   ENDDO

   alpha_me = alpha_work
ENDIF

  deallocate (h_diag)
  deallocate (dbecsum)
  deallocate (dvscfout)
  IF (nlcc_any.AND.zeu) DEALLOCATE(drhoscfh)
  IF (okpaw) THEN
     DEALLOCATE(mixin)
     DEALLOCATE(mixout)
  ENDIF
  !$acc exit data delete(dvscfins)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  if (noncolin) deallocate(dbecsum_nc)
  !$acc exit data delete(aux2)
  deallocate(aux2)
  IF ( dffts%has_task_groups ) THEN
     !
     DEALLOCATE( tg_dv  )
     DEALLOCATE( tg_psic)
     !
  ENDIF
  IF (allocated(drhoscf_aux)) deallocate(drhoscf_aux)
  IF (allocated(alpha_work)) deallocate(alpha_work)

  call stop_clock ('solve_e')
  return
end subroutine solve_e_tpw
