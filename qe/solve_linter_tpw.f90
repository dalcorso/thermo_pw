!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE solve_linter_tpw (irr, imode0, npe, drhoscf)
  !-----------------------------------------------------------------------
  !
  !    Driver routine for the solution of the linear system which
  !    defines the change of the wavefunction due to a lattice distorsion
  !    It performs the following tasks:
  !     a) computes the bare potential term Delta V | psi >
  !        and an additional term in the case of US pseudopotentials
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, diropn
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions,        ONLY : evc
  USE cell_base,            ONLY : at, omega
  USE klist,                ONLY : ltetra, lgauss, xk, wk, ngk, igk_k
  USE gvecs,                ONLY : doublegrid
  USE fft_base,             ONLY : dfftp, dffts
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, npwx
  USE scf,                  ONLY : rho, vrs
  USE uspp,                 ONLY : okvan, nkb, vkb, deeq_nc
  USE uspp_param,           ONLY : nhm
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, domag
  USE paw_variables,        ONLY : okpaw
  USE paw_onecenter,        ONLY : paw_dpotential
  USE buffers,              ONLY : save_buffer, get_buffer
  USE control_flags,        ONLY : use_gpu
  USE control_ph,           ONLY : rec_code, niter_ph, convt, &
                                   alpha_mix, rec_code_read, &
                                   where_rec, ext_recover
  USE el_phon,              ONLY : elph
  USE modes,                ONLY : u
  USE uspp,                 ONLY : nlcc_any
  USE units_ph,             ONLY : iudrho, lrdrho, iubar, lrbar, &
                                   iudvscf, iuint3paw, lint3paw
  USE units_lr,             ONLY : iuwfc, lrwfc, iudwf, lrdwf
  USE output,               ONLY : fildrho, fildvscf
  USE phus,                 ONLY : becsumort, alphap, int1_nc
  USE recover_mod,          ONLY : write_rec
  USE recover_mod_tpw,      ONLY : read_rec_tpw
  ! used to write fildrho:
  USE dfile_autoname,       ONLY : dfile_name
  USE save_ph,              ONLY : tmp_dir_save
  ! used oly to write the restart file
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp
  USE mp_asyn,              ONLY : asyn_master, with_asyn_images
  USE mp_images,            ONLY : my_image_id, root_image
  USE mp,                   ONLY : mp_sum
  USE efermi_shift_tpw,     ONLY : ef_shift_tpw, ef_shift_paw_tpw,  def
  USE lrus,         ONLY : int3_paw, becp1
  USE eqv,          ONLY : dvpsi, dpsi, evq
  USE qpoint,       ONLY : xq, nksq, ikks, ikqs
  USE qpoint_aux,   ONLY : ikmks, ikmkmqs, becpt, alphapt
  USE control_lr,   ONLY : lgamma, alpha_pv
  USE nc_mag_aux,   ONLY : int1_nc_save, deeq_nc_save, int3_save
  USE dv_of_drho_lr, ONLY : dv_of_drho
  USE fft_interfaces, ONLY : fft_interpolate, fwfft
  USE ldaU,         ONLY : lda_plus_u
  USE magnetic_charges,     ONLY : mag_charge_mode
  USE gvect,                ONLY : gg
  USE uspp_init,            ONLY : init_us_2

  IMPLICIT NONE

  INTEGER :: irr, npe, imode0
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes

  COMPLEX(DP) :: drhoscf (dfftp%nnr, nspin_mag, npe)
  ! output: the change of the scf charge

  REAL(DP) , ALLOCATABLE :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian
  REAL(DP) :: thresh, averlt, dr2, rsign
  ! thresh: convergence threshold
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  ! rsign : sign or the term in the magnetization
  REAL(DP) :: dos_ef, weight, aux_avg (2)
  ! Misc variables for metals
  ! dos_ef: density of states at Ef

  COMPLEX(DP), ALLOCATABLE, TARGET :: dvscfin(:,:,:)
  ! change of the scf potential
  COMPLEX(DP), POINTER :: dvscfins (:,:,:)
  ! change of the scf potential (smooth part only)
  COMPLEX(DP), ALLOCATABLE :: drhoscfh (:,:,:), dvscfout (:,:,:), & 
       drhoscf_aux(:,:,:)
  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  COMPLEX(DP), ALLOCATABLE :: ldos (:,:), ldoss (:,:), mixin(:),  &
       dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:,:), &
       drhoc(:), dbecsum_aux (:,:,:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  REAL(DP), ALLOCATABLE :: becsum1(:,:,:)

  LOGICAL ::             &
             exst,       & ! used to open the recover file
             lmetq0,     & ! true if xq=(0,0,0) in a metal
             all_done_asyn

  INTEGER :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             ltaver,     & ! average counter
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ndim,       & ! dimension of dbecsum
             is,         & ! counter on spin polarizations
             nrec,       & ! the record number for dvpsi and dpsi
             ipol,       & ! counter on polarization
             isolv,      & ! counter on linear systems
             nsolv,      & ! number of linear systems
             ikmk,       & ! index of mk
             ikmkmq,     & ! index of mk-mq
             npw,        & ! number of plane waves at k  
             npwq,       & ! number of plane waves at k+q
             nnr, nnrs,  & ! number of FFT mesh point
             mode          ! mode index

  INTEGER  :: iq_dummy
  REAL(DP) :: tcpu, get_clock ! timing variables
  CHARACTER(LEN=256) :: filename

  !
  IF (rec_code_read > 20 ) RETURN

  CALL start_clock ('solve_linter')
!
!  This routine is task group aware
!
  nsolv=1
  IF (noncolin.AND.domag) THEN 
     nsolv=2
     ALLOCATE (drhoscf_aux(dfftp%nnr, nspin_mag, npe))
  ENDIF

  ALLOCATE (dvscfin ( dfftp%nnr , nspin_mag , npe))
  nnr=dfftp%nnr
  IF (doublegrid) THEN
     allocate (dvscfins (dffts%nnr , nspin_mag , npe))
     nnrs=dffts%nnr
  ELSE
     dvscfins => dvscfin
     nnrs=nnr
  ENDIF
  !$acc enter data create(dvscfins(1:nnr, 1:nspin_mag, 1:npe))
  ALLOCATE (drhoscfh ( dfftp%nnr, nspin_mag , npe))
  ALLOCATE (dvscfout ( dfftp%nnr, nspin_mag , npe))
  ALLOCATE (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin_mag , npe))
  IF (okpaw) THEN
     ALLOCATE (mixin(dfftp%nnr*nspin_mag*npe+&
                                    (nhm*(nhm+1)*nat*nspin_mag*npe)/2) )
     mixin=(0.0_DP,0.0_DP)
  ELSE
     ALLOCATE(mixin(1))
  ENDIF
  IF (noncolin) ALLOCATE (dbecsum_nc (nhm,nhm, nat , nspin , npe, nsolv))
  ALLOCATE (h_diag ( npwx*npol, nbnd))
  ALLOCATE (drhoc(dfftp%nnr))
  IF (noncolin.AND.domag.AND.okvan) THEN
     ALLOCATE (int3_save( nhm, nhm, nat, nspin_mag, npe, 2))
     ALLOCATE (dbecsum_aux ( (nhm * (nhm + 1))/2 , nat , nspin_mag , npe))
  ENDIF

  !
  IF (rec_code_read == 10.AND.ext_recover) THEN
     ! restart from Phonon calculation
     IF (okpaw) THEN
        CALL read_rec_tpw(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh, &
                                                              dbecsum)
        IF (convt) THEN
           CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
        ELSE
           CALL setmixout(npe*dfftp%nnr*nspin_mag,&
           (nhm*(nhm+1)*nat*nspin_mag*npe)/2,mixin,dvscfin,dbecsum,ndim,-1)
        ENDIF
     ELSE
        CALL read_rec_tpw(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh)
     ENDIF
     rec_code=0
  ELSE
     iter0 = 0
     convt =.FALSE.
     where_rec='no_recover'
  ENDIF

  IF (ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     filename = dfile_name(xq, at, fildrho, TRIM(tmp_dir_save)//prefix, & 
                                          generate=.true., index_q=iq_dummy)
     CALL diropn (iudrho, filename, lrdrho, exst)
  ENDIF

  IF (convt) GOTO 155
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
  IF (lmetq0) THEN
     ALLOCATE ( ldos ( dfftp%nnr  , nspin_mag) )
     ALLOCATE ( ldoss( dffts%nnr , nspin_mag) )
     ALLOCATE (becsum1 ( (nhm * (nhm + 1))/2 , nat , nspin_mag))
     CALL localdos ( ldos , ldoss , becsum1, dos_ef )
     IF (.NOT.okpaw) DEALLOCATE(becsum1)
  endif
  !
  ! In this case it has recovered after computing the contribution
  ! to the dynamical matrix. This is a new iteration that has to
  ! start from the beginning.
  !
  IF (iter0==-1000) iter0=0
  !
  !   The outside loop is over the iterations
  !
  DO kter = 1, niter_ph
     iter = kter + iter0

     ltaver = 0

     lintercall = 0
     drhoscf(:,:,:) = (0.d0, 0.d0)
     dbecsum(:,:,:,:) = (0.d0, 0.d0)
     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     !
     ! DFPT+U: at each ph iteration calculate dnsscf,
     ! i.e. the scf variation of the occupation matrix ns.
     !
     IF (lda_plus_u .AND. (iter.NE.1)) &
        CALL dnsq_scf (npe, lmetq0, imode0, irr, .true.)
     !
     DO ik = 1, nksq
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq)
        !
        IF (lsda) current_spin = isk (ikk)
        !
        ! compute beta functions and kinetic energy for k-point ikq
        ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
        !
        IF ( nkb > 0 .AND.(iter==1.OR.alpha_pv>0.0_DP)) THEN
           CALL init_us_2( npwq, igk_k(1,ikq), xk(1,ikq), vkb, use_gpu)
           !$acc update host(vkb)
        ENDIF
        CALL g2_kin (ikq) 
        !
        ! Start the loop on the two linear systems, one at B and one at
        ! -B
        !
        DO isolv=1,nsolv
           IF (isolv==2) THEN
              ikmk = ikmks(ik)
              ikmkmq = ikmkmqs(ik)
              rsign=-1.0_DP
           ELSE
              ikmk=ikk
              ikmkmq=ikq
              rsign=1.0_DP
           ENDIF
           !
           ! read unperturbed wavefunctions psi(k) and psi(k+q)
           !
           IF (nksq.GT.1.OR.nsolv==2) THEN
              IF (lgamma) THEN
                 call get_buffer (evc, lrwfc, iuwfc, ikmk)
              ELSE
                 CALL get_buffer (evc, lrwfc, iuwfc, ikmk)
                 CALL get_buffer (evq, lrwfc, iuwfc, ikmkmq)
              ENDIF
           ENDIF
           !
           ! compute preconditioning matrix h_diag used by cgsolve_all
           !
           CALL h_prec (ik, evq, h_diag)
           !
           DO ipert = 1, npe
              mode = imode0 + ipert
              nrec = (ipert - 1) * nksq + ik + (isolv-1) * npe * nksq
              !
              !  and now adds the contribution of the self consistent term
              !
              IF (where_rec =='solve_lint'.OR.iter>1) THEN
                 !
                 ! After the first iteration dvbare_q*psi_kpoint is read i
                 ! from file
                 !
                 CALL get_buffer (dvpsi, lrbar, iubar, nrec)
                 !
                 ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
                 !
                 CALL add_dvscf_rhs( dvscfins, isolv, ipert, ik, npe )

              ELSE
                 !
                 !  At the first iteration dvbare_q*psi_kpoint is calculated
                 !  and written to file
                 !
                 IF (isolv==1) THEN
                    CALL dvqpsi_us_tpw (ik, u (1, mode),.false., becp1, &
                                                                    alphap )
                    ! DFPT+U: At the first ph iteration the bare perturbed 
                    ! Hubbard potential dvbare_hub_q * psi_kpoint 
                    ! is calculated and added to dvpsi.
                    !
                    IF (lda_plus_u) CALL dvqhub_barepsi_us (ik, u(1,mode))
                    !
                 ELSE
                    IF (okvan) THEN
                       deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,2)
                       int1_nc(:,:,:,:,:)=int1_nc_save(:,:,:,:,:,2)
                    ENDIF
                    CALL dvqpsi_us_tpw (ik, u (1, mode),.false., becpt, &
                                                                 alphapt)
                    IF (okvan) THEN
                       deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,1)
                       int1_nc(:,:,:,:,:)=int1_nc_save(:,:,:,:,:,1)
                    ENDIF
                 ENDIF
                 CALL save_buffer (dvpsi, lrbar, iubar, nrec)
              ENDIF
              !
              ! Ortogonalize dvpsi to valence states: ps = <evq|dvpsi>
              ! Apply -P_c^+.
              !
              CALL orthogonalize(dvpsi, evq, ikmk, ikmkmq, dpsi, npwq, .false.)
              !
              IF (where_rec=='solve_lint'.OR.iter > 1) THEN
                 !
                 ! starting value for delta_psi is read from iudwf
                 !
                 CALL get_buffer( dpsi, lrdwf, iudwf, nrec)
                 !
                 ! threshold for iterative solution of the linear system
                 !
                 thresh = MIN (1.d-1 * SQRT (dr2), 1.d-2)
              ELSE
                 !
                 !  At the first iteration dpsi and dvscfin are set to zero
                 !
                 dpsi(:,:) = (0.d0, 0.d0)
                 dvscfin (:, :, ipert) = (0.d0, 0.d0)
                 !
                 ! starting threshold for iterative solution of the linear system
                 !
                 thresh = 1.0d-2
              ENDIF
              !
              ! iterative solution of the linear system (H-eS)*dpsi=dvpsi,
              ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
              !
              CALL solve_linear_system(dvpsi, dpsi, h_diag, thresh, ik, &
                                                            isolv, lter)
              ltaver = ltaver + lter
              lintercall = lintercall + 1
              !
              ! writes delta_psi on iunit iudwf, k=kpoint,
              !
              !               if (nksq.gt.1 .or. npert(irr).gt.1)
              CALL save_buffer (dpsi, lrdwf, iudwf, nrec)
              !
              ! calculates dvscf, sum over k => dvscf_q_ipert
              !
              weight = wk (ikk)
              IF (nsolv==2) weight=weight/2.0_DP
              IF (noncolin) THEN
                 CALL incdrhoscf_nc(drhoscf(1,1,ipert),weight,ik, &
                                 dbecsum_nc(1,1,1,1,ipert,isolv), dpsi, rsign)
              ELSE
                 CALL incdrhoscf (drhoscf(1,current_spin,ipert), weight, ik, &
                            dbecsum(1,1,current_spin,ipert), dpsi)
              END IF
              ! on perturbations
           ENDDO
           ! on isolv
        ENDDO
        !  on k points
     ENDDO  
     !
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        call mp_sum ( dbecsum, intra_bgrp_comm )
     ENDIF

!     CALL compute_augmented_drho(drhoscf, dbecsum, drhoscfh, npe)

     IF (doublegrid) THEN
        DO is = 1, nspin_mag
           DO ipert = 1, npe
              CALL fft_interpolate (dffts, drhoscf(:,is,ipert), &
                                               dfftp, drhoscfh(:,is,ipert))
           ENDDO
        ENDDO
     ELSE
        CALL zcopy (npe*nspin_mag*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     ENDIF
     !
     !  In the noncolinear, spin-orbit case rotate dbecsum
     !
     IF (noncolin.and.okvan) THEN
        CALL set_dbecsum_nc(dbecsum_nc, dbecsum, npe)
        IF (nsolv==2) THEN
           dbecsum_aux=(0.0_DP,0.0_DP)
           CALL set_dbecsum_nc(dbecsum_nc(1,1,1,1,1,2), dbecsum_aux, npe)
           dbecsum(:,:,1,:)=dbecsum(:,:,1,:)+dbecsum_aux(:,:,1,:)
           dbecsum(:,:,2:4,:)=dbecsum(:,:,2:4,:)-dbecsum_aux(:,:,2:4,:)
        ENDIF
     ENDIF
     !
     !    Now we compute for all perturbations the total charge and potential
     !
     CALL addusddens (drhoscfh, dbecsum, imode0, npe, 0)
     !
     !   Reduce the delta rho across pools
     !
     CALL mp_sum ( drhoscf, inter_pool_comm )
     CALL mp_sum ( drhoscfh, inter_pool_comm )
     IF (okpaw) CALL mp_sum ( dbecsum, inter_pool_comm )
     !
     ! q=0 in metallic case deserve special care (e_Fermi can shift)
     !
     IF (okpaw) THEN
        IF (lmetq0) &
           CALL ef_shift_paw_tpw (drhoscfh, dbecsum, ldos, ldoss, becsum1, &
                                                  dos_ef, irr, npe, .false.)
        DO ipert=1,npe
           dbecsum(:,:,:,ipert)=2.0_DP *dbecsum(:,:,:,ipert) &
                               +becsumort(:,:,:,imode0+ipert)
        ENDDO
     ELSE
        IF (lmetq0) call ef_shift_tpw(drhoscfh,ldos,ldoss,dos_ef,irr,npe,.false.)
     ENDIF
     !
     !   After the loop over the perturbations we have the linear change
     !   in the charge density for each mode of this representation.
     !   Here we symmetrize them ...
     !
     CALL symmetrize_drho(drhoscfh, dbecsum, irr, npe, 1)
     !
     !   ... save them on disk and
     !   compute the corresponding change in scf potential
     !
     DO ipert = 1, npe
        IF (fildrho.NE.' ') then 
           CALL davcio_drho (drhoscfh(1,1,ipert), lrdrho, iudrho, &
                                                      imode0+ipert, +1)
!           close(iudrho)
        ENDIF
        
        CALL zcopy (dfftp%nnr*nspin_mag,drhoscfh(1,1,ipert),1,&
                                                      dvscfout(1,1,ipert),1)
        !
        ! Compute the response of the core charge density
        ! IT: Should the condition "imode0+ipert > 0" be removed?
        !
        IF (imode0+ipert > 0) THEN
           CALL addcore (imode0+ipert, drhoc)
        ELSE
           drhoc(:) = (0.0_DP,0.0_DP) 
        ENDIF
        !
        ! Compute the response HXC potential
        CALL dv_of_drho (dvscfout(1,1,ipert), .true., drhoc)
     ENDDO
     !
     !   And we mix with the old potential
     !
     CALL manage_mixing(dvscfout, dvscfin, dbecsum, mixin, npe, iter, kter, &
                                                             dr2, convt ) 
     !
     !  put the fermi energy shift if needed
     !
     IF (lmetq0.AND.convt) THEN
        IF (okpaw) THEN
           call ef_shift_paw_tpw (drhoscf, dbecsum, ldos, ldoss, becsum1, &
                                                  dos_ef, irr, npe, .true.)
        ELSE
            call ef_shift_tpw (drhoscf, ldos, ldoss, dos_ef, irr, npe, .true.)
        ENDIF
     ENDIF
     !
     ! check that convergence has been reached on ALL processors in this image
     !
     CALL check_all_convt(convt)

     IF (doublegrid) THEN
        DO ipert = 1, npe
           DO is = 1, nspin_mag
              CALL fft_interpolate (dfftp, dvscfin(:,is,ipert), dffts, &
                                                        dvscfins(:,is,ipert))
           ENDDO
        ENDDO
     ENDIF
!
!   calculate here the change of the D1-~D1 coefficients due to the phonon
!   perturbation in the PAW case and the int3 integrals in the US/PAW case
!
     CALL compute_int3_coeff(dvscfin, dbecsum, npe)

     CALL mp_sum ( ltaver, inter_pool_comm )
     CALL mp_sum ( lintercall, inter_pool_comm )
     averlt = DBLE (ltaver) / lintercall

     tcpu = get_clock ('PHONON')

     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / npe
     WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
     !
     !    Here we save the information for recovering the run from this poin
     !
     FLUSH( stdout )
     !
     rec_code=10
     IF (okpaw) THEN
        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
                                               dvscfin, drhoscfh, dbecsum)
     ELSE
        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
                                               dvscfin, drhoscfh)
     ENDIF

     IF ( check_stop_now() ) CALL stop_smoothly_ph (.false.)
     IF ( with_asyn_images.AND.my_image_id==root_image.AND.ionode ) &
                           CALL asyn_master(all_done_asyn)
     IF (convt) GOTO 155
  ENDDO
155 iter0=0
  !
  !    Here we compute the magnetic charge associated with the
  !    modes of the irreducible representation.
  !
  IF (noncolin.AND.domag.AND.lgamma) THEN 
     DO ipert = 1, npe
        drhoscf_aux(:,:,:) = drhoscfh(:,:,:)
        DO is=2,nspin_mag
           CALL fwfft ('Rho', drhoscf_aux(:,is,ipert), dfftp)
           IF (ABS(gg(1)).LT.1.d-8) THEN 
              mag_charge_mode(imode0+ipert,is-1)= &
                              omega*drhoscf_aux(dfftp%nl(1),is,ipert)
           END IF
        END DO
        CALL mp_sum(mag_charge_mode(imode0+ipert,1:3),intra_bgrp_comm)
     ENDDO
  END IF
  !
  !    A part of the dynamical matrix requires the integral of
  !    the self consistent change of the potential and the variation of
  !    the charge due to the displacement of the atoms.
  !    We compute it here.
  !
  IF (convt) THEN
     CALL drhodvus (irr, imode0, dvscfin, npe)
     IF (fildvscf.NE.' ') THEN
        DO ipert = 1, npe
           IF (lmetq0) then
                dvscfin(:,:,ipert) = dvscfin(:,:,ipert)-def(ipert)
                IF (doublegrid) dvscfins(:,:,ipert) = dvscfins(:,:,ipert) &
                                                                 -def(ipert)
           ENDIF
           CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
                         imode0 + ipert, +1 )
           IF (okpaw.AND.me_bgrp==0) CALL davcio( int3_paw(:,:,:,:,ipert), &
                                     lint3paw, iuint3paw, imode0+ipert, + 1 )
        ENDDO
        IF (elph) CALL elphel (irr, npe, imode0, dvscfins)
     ENDIF
  ENDIF

  IF (convt.AND.nlcc_any) CALL addnlcc (imode0, drhoscfh, npe)

  IF (ALLOCATED(ldoss)) DEALLOCATE (ldoss)
  IF (ALLOCATED(ldos)) DEALLOCATE (ldos)
  DEALLOCATE (h_diag)
  DEALLOCATE (dbecsum)
  IF (ALLOCATED(becsum1)) DEALLOCATE (becsum1)
  DEALLOCATE (mixin)
  IF (noncolin) DEALLOCATE (dbecsum_nc)
  DEALLOCATE (dvscfout)
  DEALLOCATE (drhoscfh)
  !$acc exit data delete(dvscfins)
  IF (doublegrid) DEALLOCATE (dvscfins)
  DEALLOCATE (dvscfin)
  DEALLOCATE (drhoc)
  IF (noncolin.AND.domag.AND.okvan) THEN
     DEALLOCATE (int3_save)
     DEALLOCATE (dbecsum_aux)
  ENDIF
  IF (ALLOCATED(drhoscf_aux)) DEALLOCATE(drhoscf_aux)

  CALL stop_clock ('solve_linter')
  RETURN
END SUBROUTINE solve_linter_tpw

