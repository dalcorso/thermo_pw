!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
MODULE response_kernels_tpw
CONTAINS
SUBROUTINE sternheimer_kernel_tpw(first_iter, npert, lrdvpsi, &
         iudvpsi, thresh, dvscfins, all_conv, avg_iter, drhoout, dbecsum,    &
         dbecsum_nc, exclude_hubbard, w_freq)
   !----------------------------------------------------------------------------
   !! Compute the density response to the perturbation dV = dV_bare + dV_ind by the
   !! non-interacting susceptibility. Solve Sternheimer equation
   !! (H[B] - e + w) * dpsi = dvpsi = -P_c+ * (dV_bare + dV_ind[B]) * psi.
   !! If w is not zero or time_reversal is not true (because B is nonzero)
   !! computes also
   !! (H - e - w) * dpsi = dvpsi = -R_c+ * (dV_bare + dV_ind[-B]) * psi.
   !!
   !! where R_c is the projector with the time reversed wavefunctions computed
   !! at minus B.
   !! Perfoms the following tasks:
   !!  a) reads the bare potential term Delta V | psi > from buffer iudvpsi.
   !!  b) adds to it the screening term Delta V_{SCF} | psi >.
   !!     If lda_plus_u=.true. compute also the SCF part
   !!     of the response Hubbard potential.
   !!  c) applies P_c^+ (orthogonalization to valence states).
   !!  d) calls cgsolve_all to solve the linear system.
   !!  e) returns the Delta rho, and if lda_plus_u=.true. also return dbecsum
   !!
   !! dV_bare * psi is read from buffer iudvpsi, so they must be already calculated.
   !! dV_ind is given by input dvscfins, and dV_ind * psi is calculated in apply_dpot_bands.
   !!
   !! For USPPs, adddvscf is called, so relevant arrays must be already calculated.
   !! For DFT+U, adddvhubscf is called, so relevant arrays must be already calculated.
   !!
   !! Results are added to drhoout, dbecsum, dbecsum_nc, so these output arrays should
   !! be initialized before calling this subroutine.
   !!
   !! Input:
   !!    - first_iter: true if the first iteration, where dvscfins = 0
   !!    - npert: number of perturbations
   !!    - lrdvpsi: record length for the buffer storing dV_bare * psi
   !!    - iudvpsi: unit for the buffer storing dV_bare * psi
   !!    - thresh: threshold for solving Sternheimer equation
   !!    - dvscfins: dV_ind calculated in the previous iteration
   !!    - exclude_hubbard: If TRUE, do not add the response of the Hubbard potential.
   !!                       Used in hp.x (Default: FALSE)
   !!
   !! Output:
   !!    - avg_iter: average number of iterations for the linear equation solver
   !!    - drhoout: induced charge density (dffts, without augmentation term)
   !!    - dbecsum: becsum with dpsi
   !!    - dbecsum_nc: becsum with dpsi. Optional, used if noncolin is true.
   !!
   !----------------------------------------------------------------------------
   USE kinds,                 ONLY : DP
   USE io_global,             ONLY : stdout
   USE mp,                    ONLY : mp_sum
   USE mp_pools,              ONLY : inter_pool_comm
   USE buffers,               ONLY : get_buffer, save_buffer
   USE fft_base,              ONLY : dffts
   USE ions_base,             ONLY : nat
   USE klist,                 ONLY : xk, wk, ngk, igk_k
   USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                 ONLY : nbnd, npwx, et, g2kin
   USE scf,                   ONLY : v_of_0
   USE wavefunctions,         ONLY : evc
   USE noncollin_module,      ONLY : noncolin, domag, npol, nspin_mag
   USE uspp,                  ONLY : vkb
   USE uspp_param,            ONLY : nhm
   USE uspp_init,             ONLY : init_us_2
   USE ldaU,                  ONLY : lda_plus_u
   USE units_lr,              ONLY : iuwfc, lrwfc, lrdwf, iudwf
   USE control_lr,            ONLY : alpha_pv, nbnd_occ, lgamma
   USE optical,               ONLY : current_w
   USE qpoint,                ONLY : nksq, ikks, ikqs
   USE qpoint_aux,            ONLY : ikmks, ikmkmqs, becpt
   USE eqv,                   ONLY : dpsi, dvpsi, evq
   USE apply_dpot_mod,        ONLY : apply_dpot_bands
   USE lr_nc_mag,             ONLY : lr_apply_time_reversal
   USE incdrhoscf_mod,        ONLY : incdrhoscf, incdrhoscf_nc
   USE linear_solvers,        ONLY : ccg_many_vectors
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: first_iter
   !! true if the first iteration.
   LOGICAL, INTENT(IN), OPTIONAL :: exclude_hubbard
   !! true if ignoring the Hubbard response term
   INTEGER, INTENT(IN) :: npert
   !! number of perturbations
   INTEGER, INTENT(IN) :: lrdvpsi
   !! record length for the buffer storing dV_bare * psi
   INTEGER, INTENT(IN) :: iudvpsi
   !! unit for the buffer storing dV_bare * psi
   REAL(DP), INTENT(IN) :: thresh
   !! threshold for solving the linear equation
   LOGICAL, INTENT(OUT) :: all_conv
   !! True if converged at all k points and perturbations
   REAL(DP), INTENT(OUT) :: avg_iter
   !! average number of iterations for the linear equation solver
   COMPLEX(DP), INTENT(INOUT) :: dvscfins(dffts%nnr, nspin_mag, npert)
   !! dV_ind calculated in the previous iteration
   COMPLEX(DP), INTENT(INOUT) :: drhoout(dffts%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   COMPLEX(DP), INTENT(IN), OPTIONAL :: w_freq  ! complex frequency
   !
   LOGICAL :: conv_root
   !! true if linear system is converged
   LOGICAL :: exclude_hubbard_
   !! Local variable to set the default of exclude_hubbard to false
   LOGICAL :: time_reversed
   !! true if solving for time reversed wave functions
   INTEGER :: ikk, ikq, npw, npwq, ipert, num_iter, ik, nrec, ikmk, ikmkmq, &
              ibnd, ig, nrect
   !! counters
   INTEGER :: tot_num_iter
   !! total number of iterations in cgsolve_all
   INTEGER :: tot_cg_calls
   !! total number of cgsolve_all calls
   INTEGER :: nsolv
   !! number of linear systems solved each iteration
   INTEGER :: isolv
   !! counter on linear systems solved each iteration
   REAL(DP) :: anorm
   !! the norm of the error of the conjugate gradient solution
   REAL(DP) :: rsign
   !! sign of the term in the magnetization
   REAL(DP) :: weight
   !! weight of the k point
   REAL(DP) :: alpha_pv0
   !! save alpha_pv
   REAL(DP), ALLOCATABLE :: h_diag(:, :)
   !! diagonal part of the Hamiltonian, used for preconditioning
   COMPLEX(DP), ALLOCATABLE :: h_diagc(:, :)
   !! diagonal part of the Hamiltonian, used for preconditioning complex form
   COMPLEX(DP) , ALLOCATABLE :: aux2(:, :)
   !! temporary storage used in apply_dpot_bands
   LOGICAL :: time_reversal
   !! if .TRUE. time reversal is available
   LOGICAL :: magnetic_field
   !! .TRUE. if there is a noncolinear magnetization
   COMPLEX(DP) :: w
   !! frequency
   EXTERNAL ch_psi_all, cg_psi
   !! functions passed to cgsolve_all (static case)
   EXTERNAL apply_ac, ccg_psi_tpw, scal_prod
   !! functions passed to ccgsolve_all (dynamic case, TPW routines)
   EXTERNAL ch_psi_all_complex, ccg_psi
   !! functions passed to ccgsolve_all (dynamic case, QE routine)
   COMPLEX(DP) :: scal_prod
   !! functions passed to cgsolve_many_vectors
   !
   ! Initialization
   !
   CALL start_clock("sth_kernel")
   !
   exclude_hubbard_ = .FALSE.
   IF (PRESENT(exclude_hubbard)) exclude_hubbard_ = exclude_hubbard
   !
   magnetic_field=(noncolin.AND.domag)
   time_reversal=.NOT. magnetic_field
   !
   !  set the frequency
   !
   w=(0.0,0.0)
   IF (PRESENT(w_freq)) w=w_freq
   
   nsolv=1
   IF (.NOT.time_reversal.OR.ABS(w)>0.0_DP)  nsolv=2

   alpha_pv0=alpha_pv
   alpha_pv=alpha_pv0 + REAL(w)
   !
   IF (ABS(w)>0.0_DP) THEN
      ALLOCATE(h_diagc(npwx*npol, nbnd))
   ELSE
      ALLOCATE(h_diag(npwx*npol, nbnd))
      h_diag = 0.d0
   ENDIF
   ALLOCATE(aux2(npwx*npol, nbnd))
   aux2 = (0.d0, 0.d0)
   !
   !$acc enter data create(aux2(1:npwx*npol, 1:nbnd))
   !

   all_conv = .TRUE.
   tot_num_iter = 0
   tot_cg_calls = 0
   !
   DO ik = 1, nksq
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)

      IF (lsda) current_spin = isk(ikk)
      !
      ! Set time-reversed k and k+q points
      !
      DO isolv=1, nsolv
         IF (magnetic_field.AND.isolv==2) THEN
            time_reversed=.TRUE.
            ikmk = ikmks(ik)
            ikmkmq = ikmkmqs(ik)
            rsign = -1.0_DP
         !
         !  change the sign of the magnetic field if required
         !
            CALL lr_apply_time_reversal(first_iter, 2, dvscfins)
         ELSE
            time_reversed=.FALSE.
            ikmk = ikk
            ikmkmq = ikq
            rsign = 1.0_DP
         ENDIF
         !
         !
         ! reads unperturbed wavefunctions psi_k in G_space, for all bands
         ! if q=0, evq is a pointer to evc
         !
         IF (nksq > 1 .OR. nsolv==2) THEN
            IF (lgamma) THEN
               !civn: in this case evq is a pointer to evc
               CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
               !$acc update device(evc)
            ELSE
               !civn: in this case evq is allocated separately and needs to be updated on device
               CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
               !$acc update device(evc)
               CALL get_buffer(evq, lrwfc, iuwfc, ikmkmq)
               !$acc update device(evq)
            ENDIF
         ENDIF
         !
         ! compute beta functions and kinetic energy for k-point ik
         ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
         !
         CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
         !$acc update host(vkb)
         !
         ! compute the kinetic energy g2kin: (k+q+G)^2
         !
         CALL g2_kin(ikq)
         !
         ! compute preconditioning matrix h_diag used by (c)cgsolve_all 
         !
         IF (ABS(w)>0.0_DP) THEN
            !
            ! finite frequency
            !
            IF (isolv==1) THEN
               CALL h_prec_freq(ik, w, h_diagc)
            ELSE
               CALL h_prec_freq(ik, -w, h_diagc)
            ENDIF    
         ELSE
            !
            !  static case
            !
            CALL h_prec(ik, evq, h_diag)
            !
         ENDIF
         !
         DO ipert = 1, npert
            !
            ! read P_c^+ x psi_kpoint into dvpsi.
            !
            nrec = (ipert - 1) * nksq + ik 
            nrect = nrec
            IF (isolv==2) nrect = nrect + npert * nksq
            !
            IF (magnetic_field) THEN
               CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrect)
            ELSE
               CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrec)
            ENDIF
            !
            IF (.NOT. first_iter) THEN
               !
               ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
               ! dvscf_q from previous iteration (mix_potential)
               !
               CALL apply_dpot_bands(ik, nbnd_occ(ikk), dvscfins(:, :, ipert),&
                                                                   evc, aux2)

               dvpsi = dvpsi + aux2
               !
               !  In the case of US pseudopotentials there is an additional
               !  self-consistent term which comes from the dependence of D on
               !  V_{eff} on the bare change of the potential
               !
               CALL adddvscf(ipert, ik, time_reversed)
               !
               ! DFPT+U: add to dvpsi the scf part of the response
               ! Hubbard potential dV_hub
               !
               IF (lda_plus_u .AND. (.NOT. exclude_hubbard_)) &
                         CALL adddvhubscf(ipert, ik)
               !
            ENDIF ! .NOT. first_iter
            !
            ! Orthogonalize dvpsi to valence states
            !
            IF (ABS(w)>0.0_DP) THEN
               IF (isolv==1) THEN
                  CALL orthogonalize_omega(dvpsi, evq, ikmk, ikmkmq, dpsi, &
                                                                  npwq, w)
               ELSE
                  CALL orthogonalize_omega(dvpsi, evq, ikmk, ikmkmq, dpsi, &
                                                                  npwq, -w)
               ENDIF
            ELSE
               CALL orthogonalize(dvpsi, evq, ikmk, ikmkmq, dpsi, npwq, &
                                                                    .FALSE.)
            ENDIF
            !
            ! Initial guess for dpsi
            !
            IF (first_iter) THEN
               !
               !  At the first iteration dpsi is set to zero
               !
               dpsi(:, :) = (0.d0,0.d0)
            ELSE
               !
               ! starting value for delta_psi is read from iudwf
               !
               CALL get_buffer(dpsi, lrdwf, iudwf, nrect)
            ENDIF
            !
            ! iterative solution of the linear system (H-e)*dpsi=dvpsi
            ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
            !
            conv_root = .TRUE.
            !
            ! TODO: should nbnd_occ(ikk) be nbnd_occ(ikmk)?
            IF (ABS(w)>0.0_DP) THEN
               current_w=w
               IF (isolv==2) current_w=-w
               !
               !  Complex or imaginary frequency. Use bicojugate gradient.
               !
               IF (.TRUE.) THEN
                  CALL ccgsolve_all(ch_psi_all_complex, ccg_psi, et(1, ikmk), &
                          dvpsi, dpsi, h_diagc, npwx, npwq, thresh, ik, &
                          num_iter, conv_root, anorm, nbnd_occ(ikk), npol, & 
                          current_w)
               ELSE
                  CALL ccg_many_vectors (apply_ac,ccg_psi_tpw,scal_prod, &
                    dvpsi, dpsi, h_diagc, npwx*npol, npwq, thresh, ik,   &
                    num_iter, conv_root, anorm, nbnd_occ(ikk))
               ENDIF

            ELSE
               CALL cgsolve_all(ch_psi_all, cg_psi, et(1, ikmk), dvpsi, dpsi, &
                             h_diag, npwx, npwq, thresh, ik, num_iter, &
                             conv_root, anorm, nbnd_occ(ikk), npol)
            ENDIF
            weight=wk(ikk)
!            IF (nsolv==2) weight=weight/2.0_DP
! 
!    this division by two is made by sterheimer_postproc
!
            !
            tot_num_iter = tot_num_iter + num_iter
            tot_cg_calls = tot_cg_calls + 1
            !
            IF (.NOT. conv_root) THEN
               all_conv = .FALSE.
               WRITE( stdout, "(5x, 'kpoint', i4, ' sternheimer_kernel: &
               &root not converged, thresh < ', es10.3)") ik, anorm
            ENDIF
            !
            ! writes delta_psi on iunit iudwf, k=kpoint,
            !
            CALL save_buffer(dpsi, lrdwf, iudwf, nrect)
            !
            ! calculates dvscf, sum over k => dvscf_q_ipert
            !
            IF (noncolin) THEN
               CALL incdrhoscf_nc(drhoout(1,1,ipert), weight, ik, &
                               dbecsum_nc(1,1,1,1,ipert), dpsi, rsign)
            ELSE
               CALL incdrhoscf(drhoout(1,current_spin,ipert), weight, &
                            ik, dbecsum(1,1,current_spin,ipert), dpsi)
            ENDIF
         ENDDO ! ipert
         !
         !  reset the original magnetic field if it was changed
         !
         IF (time_reversed) CALL lr_apply_time_reversal(first_iter, 1, &
                                                               dvscfins)
      ENDDO ! isolv
   ENDDO ! ik
   current_w=w
   !
   CALL mp_sum(tot_num_iter, inter_pool_comm)
   CALL mp_sum(tot_cg_calls, inter_pool_comm)
   avg_iter = REAL(tot_num_iter, DP) / REAL(tot_cg_calls, DP)
   !
   !$acc exit data delete(aux2)
   !
   DEALLOCATE(aux2)
   IF (ABS(w)>0.0_DP) THEN
      DEALLOCATE(h_diagc)
   ELSE
      DEALLOCATE(h_diag)
   ENDIF
   alpha_pv=alpha_pv0
   !
   CALL stop_clock("sth_kernel")
   !
!----------------------------------------------------------------------------
END SUBROUTINE sternheimer_kernel_tpw
!------------------------------------------------------------------------------
!
END MODULE response_kernels_tpw
