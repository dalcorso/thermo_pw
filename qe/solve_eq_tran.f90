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
subroutine solve_eq_tran(iu, flag)
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
  USE klist,                 ONLY : lgauss, xk, wk, ngk, igk_k
  USE gvect,                 ONLY : g
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts, dtgs
  USE fft_parallel,          ONLY : tg_cgather
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,              ONLY : domag
  USE wvfct,                 ONLY : nbnd, npwx, g2kin,  et
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer, save_buffer
  USE wavefunctions_module,  ONLY : evc
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : upf, nhm
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag, nspin_lsda
  USE scf,                   ONLY : rho, v_of_0
  USE gvect,                 ONLY : gg, nl
  USE paw_variables,         ONLY : okpaw
  USE paw_add_onecenter,     ONLY : paw_deqtranpotential
  USE paw_add_symmetry,      ONLY : paw_deqsymmetrize

  USE eqv,                   ONLY : dpsi, dvpsi, evq
  USE units_ph,              ONLY : lrdwf, iudwf, lrwfc, iuwfc, lrdrho, &
                                    iudrho, lrbar, iubar
  USE output,                ONLY : fildrho
  USE control_ph,            ONLY : ext_recover, rec_code, &
                                    lnoloc, convt, tr2_ph, nmix_ph, &
                                    alpha_mix, lgamma_gamma, niter_ph, &
                                    flmixdpot, rec_code_read
  USE control_lr,            ONLY : lgamma, alpha_pv, nbnd_occ
  USE dv_of_drho_lr,         ONLY : dv_of_drho
  USE lrus,                  ONLY : int3_paw
  USE qpoint,                ONLY : xq, nksq, ikks, ikqs
  USE recover_mod,           ONLY : read_rec, write_rec

  USE optical,               ONLY : current_w, fru, iu1dwf, lr1dwf, &
                                    chipm, chimp, chixx, chixy
  USE freq_ph,               ONLY : fiu
  USE linear_solvers,        ONLY : ccg_many_vectors
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm, ntask_groups
  USE mp,                    ONLY : mp_sum

  implicit none

  INTEGER, INTENT(IN) :: iu
  INTEGER, INTENT(IN) :: flag  ! if 1 computes the \chi+- response
                               ! if 2 computes the \chi-+ response and
                               ! all the transverse tensor \chi_xx=\chi_yy and 
                               ! \chi_xy=-chi_yx if \chi+- has been
                               ! computed before
                               !

  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  complex(DP), allocatable :: h_diag (:,:)
  complex(DP), allocatable :: h_diag1 (:,:)
  real(DP), allocatable :: h_dia (:,:), s_dia(:,:)
  real(DP), allocatable :: h_diagr (:,:)
  ! h_diag: diagonal part of the Hamiltonian

  complex(DP) , allocatable, target ::      &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  complex(DP) , pointer ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  complex(DP) , allocatable ::   &
                   dpsi1(:,:),   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   drhoscfout (:), & ! change of the scf charge (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   mixin(:), mixout(:), &  ! auxiliary for paw mixing
                   aux1 (:,:),  ps (:,:), &
                   tg_dv(:,:), &
                   tg_psic(:,:), aux2(:,:), dvpsi1(:,:)


  complex(DP), EXTERNAL :: zdotc      ! the scalar product function

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: kter, iter0, ipol, ibnd, iter, lter, ik, ikk, ikq, &
             ig, is, nrec, ndim, ios, spin_psi0, nmix_ph_eff, npw, npwq
  integer, allocatable :: save_ikqs(:)
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


  call start_clock ('solve_eq_tran')
  !
  !  This routine is task group aware
  !
  w=CMPLX(fru(iu),fiu(iu))
  ldpsi1=ABS(w)>1.D-7
  nmix_ph_eff=max(nmix_ph,8)
  alpha_pv0=alpha_pv
  alpha_pv=alpha_pv0 + REAL(w)

  allocate (dvscfin( dfftp%nnr, nspin_mag, 1))
  if (doublegrid) then
     allocate (dvscfins(dffts%nnr, nspin_mag, 1))
  else
     dvscfins => dvscfin
  endif
  allocate (dvpsi1(npwx,nbnd))
  allocate (dvscfout(dfftp%nnr, nspin_mag, 1))
  allocate (drhoscfout(dfftp%nnr))
  IF (okpaw) THEN
     ALLOCATE (mixin(dfftp%nnr+(nhm*(nhm+1)*nat)/2) )
     ALLOCATE (mixout(dfftp%nnr+(nhm*(nhm+1)*nat)/2) )
  ENDIF
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 1))
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin, 1))
  IF (ldpsi1) THEN
     allocate (dpsi1(npwx*npol,nbnd))
     allocate (h_diag(npwx*npol, nbnd))
     allocate (h_diag1(npwx*npol, nbnd))
     allocate (h_dia(npwx, npol))
     allocate (s_dia(npwx, npol))
  ELSE
     allocate (h_diagr(npwx*npol, nbnd))
  ENDIF
  allocate (aux1(dffts%nnr,npol))
  allocate( save_ikqs(nksq) )
  save_ikqs(:)=ikqs(:)

  allocate (aux2(npwx*npol, nbnd))
  IF (okpaw) mixin=(0.0_DP,0.0_DP)

  if (rec_code_read == -20.AND.ext_recover) then
     ! restarting in Electric field calculation
     IF (okpaw) THEN
        CALL read_rec(dr2, iter0, 1, dvscfin, dvscfins, dvscfout, dbecsum)
        CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixin, dvscfin, dbecsum, ndim, -1 )
     ELSE
        CALL read_rec(dr2, iter0, 1, dvscfin, dvscfins)
     ENDIF
  else if (rec_code_read > -20 .AND. rec_code_read <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
  endif
  incr=1
  IF ( dtgs%have_task_groups ) THEN
     !
     v_siz =  dtgs%tg_nnr * dtgs%nogrp
     ALLOCATE( tg_dv   ( v_siz, nspin_mag ) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = dtgs%nogrp
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
  do kter = 1, niter_ph

!     write(6,*) 'kter', kter
     FLUSH (stdout)
     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)
     IF (noncolin) dbecsum_nc=(0.d0,0.d0)

     do ik = 1, nksq
        ikk = ikks(ik)
        npw=ngk(ikk)
        spin_psi0=isk(ikk)
!
!   when psi0 has spin up we need to apply the Hamiltonian for spin down
!   and viceversa
!
        current_spin=spin_psi0+1
        IF (current_spin==3) current_spin=1
!
! takes the k+q point with opposite spin.
! this should work also with pools
!
        IF (spin_psi0==1) THEN
           ikq = save_ikqs(ik) + 2
        ELSE
           ikq = save_ikqs(ik) - 2
        ENDIF
        ikqs(ik)=ikq
        IF (spin_psi0==isk(ikq)) &
           CALL errore('solve_eq_tran','problem with spin',1)
        npwq=ngk(ikq)
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
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
           CALL usnldiag( npwq, h_dia, s_dia )

           do ibnd = 1, nbnd_occ (ikk)
              !
              do ig = 1, npwq
                 aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                    (et(ibnd,ikk)+w)*s_dia(ig,1)
                 h_diag(ig,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa
                 aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                    (et(ibnd,ikk)-w)*s_dia(ig,1)
                 h_diag1(ig,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa
              end do
            end do
         ELSE
            CALL h_prec (ik, evc, h_diagr) 
            h_diagr=0.0_DP
            do ibnd = 1, nbnd_occ (ikk)
               !
               do ig = 1, npwq
                  aa=1.0_DP/h_diagr(ig,ibnd)-et(ibnd,ikk)-REAL(w,KIND=DP)
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
              IF( dtgs%have_task_groups ) THEN
                 CALL tg_cgather( dffts, dtgs, dvscfins(:,flag,ipol), &
                                                             tg_dv(:,flag))
              ENDIF

              aux2=(0.0_DP,0.0_DP)
              do ibnd = 1, nbnd_occ (ikk), incr
                 IF ( dtgs%have_task_groups ) THEN
                    call cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, &
                                      nbnd_occ (ikk) )
                    call apply_dpot(v_siz, tg_psic, tg_dv, 1)
                    call cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, &
                                   nbnd_occ (ikk))
                 ELSE
                    call cft_wave (ik, evc (1, ibnd), aux1, +1)
                    call apply_dpot(dffts%nnr, aux1, &
                                    dvscfins(1,1,ipol), 1)
                    call cft_wave (ik, aux2 (1, ibnd), aux1, -1)
                 ENDIF
              enddo
              dvpsi=dvpsi+aux2
              call adddvscf_tran(spin_psi0, ik)
              !
           else
              !
              !  At the first iteration dvbare_q*psi_kpoint is calculated
              !  and written to file. The external perturbation is assumed
              !  to be along x so b_ext+ = b_ext-
              !
              CALL dveqpsi_us (ik)
              call save_buffer (dvpsi, lrbar, iubar, nrec)
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
           !

           IF (ldpsi1) THEN
              dvpsi1=dvpsi
              CALL orthogonalize_omega(dvpsi1, evq, ikk, ikq, dpsi, npwq, -w)
           END IF
           CALL orthogonalize_omega(dvpsi, evq, ikk, ikq, dpsi, npwq, w)
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
           IF (ldpsi1) THEN
!
!    complex frequency, two wavefunctions must be computed
!
              weight=wk(ikk)/2.0_DP
           ELSE
!
!    zero frequency. Only one pertubed wavefunction is needed
!
              weight=wk(ikk)
           END IF
           
           IF (spin_psi0==flag) &
              call incdrhoscf (dvscfout(1,1,ipol), weight, ik, &
                                      dbecsum(1,1,1,ipol), dpsi)

           IF (ldpsi1) THEN
!
!   In this case compute also the wavefunction at frequency -w.
!
              current_w=-w
              CALL ccg_many_vectors (apply_ac,ccg_psi_tpw,scal_prod,dvpsi1, &
                             dpsi1, h_diag, npwx*npol, npwq, thresh, ik, &
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
              IF (spin_psi0 /= flag) &
                 call incdrhoscf (dvscfout(1,1,ipol), weight, &
                                          ik, dbecsum(1,1,1,ipol), dpsi1)
           END IF
        enddo   ! on polarizations
     enddo      ! on k points
     current_w=w

     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions
     !  coming from each slice of bands
     !
     call mp_sum ( dbecsum, intra_bgrp_comm )

     if (doublegrid) call cinterpolate (dvscfout, dvscfout, 1)
     !
     nspin_mag=1
     call addusddenseq (dvscfout, dbecsum)
     nspin_mag=2
     !
     !   dvscfout contains the (unsymmetrized) linear charge response
     !   for the three polarizations - symmetrize it
     !
     call mp_sum ( dvscfout, inter_pool_comm )
     
     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
!
!   the y component of the magnetization is obtained multiplying by -i
!
     dvscfout(:,1,1) = 2.0_DP*dvscfout(:,1,1)
     CALL psymeq(dvscfout)
     drhoscfout(:) = dvscfout(:,1,1)
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     IF (lnoloc) THEN
        dvscfout=(0.d0,0.d0)
        convt=.TRUE.
        GOTO 1001
     ELSE
        call dv_of_drho_tran (dvscfout)
     ENDIF
     !
     !   mix the new potential with the old
     !
     IF (okpaw) THEN
     !
     !  In this case we mix also dbecsum
     !
        call setmixout(dfftp%nnr,(nhm*(nhm+1)*nat)/2, &
                    mixout, dvscfout, dbecsum, ndim, -1 )
        call mix_potential (2*dfftp%nnr+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, tr2_ph/npol, iter, &
                         nmix_ph_eff, flmixdpot, convt)
        call setmixout(dfftp%nnr,(nhm*(nhm+1)*nat)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )
     ELSE
        call mix_potential (2*dfftp%nnr, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2,  tr2_ph / npol, iter, nmix_ph_eff, flmixdpot, convt)
     ENDIF

     if (doublegrid) call cinterpolate (dvscfin,dvscfins,-1)

     IF (okpaw) THEN
        IF (noncolin.AND.domag) THEN
!           call PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
        ELSE
!
!    The presence of c.c. in the formula gives a factor 2.0
!
           dbecsum=2.0_DP * dbecsum
           IF (.NOT. lgamma_gamma) CALL PAW_deqsymmetrize(dbecsum)
           call PAW_deqtranpotential(dbecsum,rho%bec,int3_paw)
        ENDIF
     ENDIF

     call newdq(dvscfin,1)

1001 CONTINUE

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
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 1, dvscfin, &
                                                       dvscfout, dbecsum)
     ELSE
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 1, dvscfin)
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

  CALL fwfft ('Dense', drhoscfout, dfftp)

  IF (flag==1) chipm(iu)=(0.0_DP,0.0_DP)

  IF (flag==2) THEN
     chimp(iu)=(0.0_DP,0.0_DP)
     chixx(iu)=(0.0_DP,0.0_DP)
     chixy(iu)=(0.0_DP,0.0_DP)
  END IF

  IF (ABS(gg(1))<1.d-8) THEN
     IF (flag==1) THEN
        chipm(iu)=drhoscfout(nl(1))  
     ELSE
        chimp(iu)=drhoscfout(nl(1))  
        IF (ABS(chipm(iu))>1.D-9) THEN
           chixx(iu)=(chipm(iu)+chimp(iu))*0.25_DP
           chixy(iu)=(chipm(iu)-chimp(iu))*0.25_DP *(0.0_DP,-1.0_DP)
        ENDIF
     ENDIF
  ENDIF

  IF (flag==1) THEN
     CALL mp_sum(chipm(iu),intra_bgrp_comm)        
  ELSE
     CALL mp_sum(chimp(iu),intra_bgrp_comm)        
     CALL mp_sum(chixx(iu),intra_bgrp_comm)        
     CALL mp_sum(chixy(iu),intra_bgrp_comm)        
  ENDIF

  WRITE( stdout, '(/,6x,"Transverse susceptibility at &
                     &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)

  WRITE( stdout, '(44x,f9.4," +",f9.4," i eV")') current_w * rytoev

  IF (flag==1 .OR. ABS(chipm(iu))>1.D-9) &
      WRITE(stdout,'(/,6x,"chi+-(q,w) =", 2f15.6)') chipm(iu)

  IF (flag==2) THEN
     WRITE(stdout,'(/,6x,"chi-+(q,w) =", 2f15.6)') chimp(iu)
     WRITE(stdout,'(/,6x,"chixx(q,w) =", 2f15.6)') chixx(iu)
     WRITE(stdout,'(/,6x,"chixy(q,w) =", 2f15.6)') chixy(iu)
  END IF

  ikqs(:)=save_ikqs(:)
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
  deallocate (dbecsum)
  deallocate (dvscfout)
  IF (okpaw) THEN
     DEALLOCATE(mixin)
     DEALLOCATE(mixout)
  ENDIF
  deallocate (drhoscfout)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  if (noncolin) deallocate(dbecsum_nc)
  deallocate (aux2)
  deallocate (save_ikqs)
  IF ( dtgs%have_task_groups ) THEN
     !
     DEALLOCATE( tg_dv  )
     DEALLOCATE( tg_psic)
     !
  ENDIF
  alpha_pv=alpha_pv0

  call stop_clock ('solve_eq_tran')
  return
end subroutine solve_eq_tran
