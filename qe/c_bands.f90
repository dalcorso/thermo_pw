!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! Copyright (C) 2017 Dal Corso Andrea
!
!  This version of the routine should be called only by the phonon code.  
!  It allows the computation of the bands in the k points outside 
!  the irreducible wedge that corresponds to the full point group 
!  by rotating the wavefunctions in the irreducible wedge. 
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE c_bands_nscf_tpw( )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine for Hamiltonian diagonalization routines
  ! ... specialized to non-self-consistent calculations (no electric field)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer
  USE basis,                ONLY : starting_wfc
  USE fft_base,             ONLY : dfftp
  USE symm_base,            ONLY : s, sr, ftau, invs, t_rev
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g, nl, gg, ngm
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_flags,        ONLY : ethr, restart, isolve, io_level, iverbosity
  USE save_ph,              ONLY : tmp_dir_save
  USE io_files,             ONLY : tmp_dir, prefix
  USE buffers,              ONLY : open_buffer, close_buffer
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions_module, ONLY : evc
  USE control_lr,           ONLY : lgamma
  USE mp_asyn,              ONLY : asyn_master, with_asyn_images
  USE mp_images,            ONLY : my_image_id, root_image
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum, mp_min
  USE io_global,            ONLY : ionode
  USE check_stop,           ONLY : check_stop_now
  USE band_computation,     ONLY : diago_bands, isym_bands, ik_origin, nks0
  USE noncollin_module,     ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  REAL(DP) :: avg_iter, ethr_
  ! average number of H*psi products
  INTEGER :: ik_, ik, ik_eff, ibnd, nkdum, npw, ios, ipol, ind1, ind2, gk(3)
  INTEGER :: ishift, ik_diago, ik_sym
  ! ik_: k-point already done in a previous run
  ! ik : counter on k points
  LOGICAL :: all_done_asyn
  !
  REAL(DP), EXTERNAL :: get_clock
  COMPLEX(DP), ALLOCATABLE :: psic(:,:,:), evcr(:,:)
  COMPLEX(DP) :: d_spin(2,2)
  INTEGER :: has_e, ig, iuawfc, lrawfc
  LOGICAL :: exst_mem, exst
  !
  CALL start_clock( 'c_bands' )
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands(ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  !
  DO ik = 1, ik_
     CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
  END DO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSE IF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")')
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve)
  END IF
  IF (tmp_dir /= tmp_dir_save) THEN
     iuawfc = 20
     lrawfc = nbnd * npwx * npol
     CALL open_buffer (iuawfc, 'wfc', lrawfc, io_level, exst_mem, exst, &
                                                         tmp_dir_save)
     IF (.NOT.exst.AND..NOT.exst_mem) THEN
        CALL errore ('c_bands_tpw', 'file '//trim(prefix)//'.wfc not found', 1)
     END IF
  ENDIF
!
!  find the minimum number of k points diagonalized in all pools
!
  nkdum=0
  DO ik=1,nks
     IF (ik <= ik_) CYCLE
     IF (diago_bands(ik)) nkdum=nkdum+1
  ENDDO
  CALL mp_min(nkdum,inter_pool_comm)
  !
  ! ... For each k point (except those already calculated if restarting)
  ! ... diagonalizes the hamiltonian
  !
  ik_diago=0
  ik_sym=0
  DO ishift=0,1
     IF ( ik_>0 .AND. MOD(ik_,2)==0 .AND. ishift==0) CYCLE
!
!  first compute the k vectors (odd index) and then the k+q
!
     k_loop: DO ik = ik_+1+ishift, nks, 2
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
        IF (diago_bands(ik)) THEN
           ik_diago=ik_diago+1
           current_k = ik
           IF ( lsda ) current_spin = isk(ik)
           call g2_kin( ik )
           ! 
           ! ... More stuff needed by the hamiltonian: nonlocal projectors
           !
           IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
           !
           ! ... Needed for LDA+U
           !
           IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
               CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )
           !
           ! ... calculate starting  wavefunctions
           !
           IF ( iverbosity > 0 ) WRITE( stdout, 9001 ) ik
           !
           IF ( TRIM(starting_wfc) == 'file' ) THEN
              !
              CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
              !
           ELSE
              !
              CALL init_wfc ( ik )
              !
           END IF
           !
           ! ... diagonalization of bands for k-point ik
           !
           call diag_bands ( 1, ik, avg_iter )
           !
           ! ... save wave-functions (unless disabled in input)
           !
           IF ( io_level > -1 ) CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
           !
           ! ... beware: with pools, if the number of k-points on different
           ! ... pools differs, make sure that all processors are still in
           ! ... the loop on k-points before checking for stop condition
           !
!           nkdum  = kunit * ( nkstot / kunit / npool )
           IF (ik_diago .le. nkdum) THEN
              !
              ! ... stop requested by user: save restart information,
              ! ... save wavefunctions to file
              !
              IF (check_stop_now()) THEN
                 CALL save_in_cbands(ik, ethr, avg_iter, et )
                 RETURN
              END IF
           ENDIF
           !
           ! report about timing
           !
           IF ( iverbosity > 0 ) THEN
              WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
              FLUSH( stdout )
           ENDIF
           IF ( with_asyn_images.AND.my_image_id==root_image.AND.ionode ) &
                        CALL asyn_master(all_done_asyn)
           !
        ELSE
           IF (ik_origin(ik)/=ik) THEN
              ik_sym=ik_sym+1
              CALL get_buffer ( evc, nwordwfc, iunwfc, ik_origin(ik) )
              ALLOCATE(psic(dfftp%nnr, npol, nbnd))
              ALLOCATE(evcr(npwx*npol, nbnd))
              npw=ngk(ik_origin(ik))
              psic=(0.0_DP,0.0_DP)
              DO ipol=1, npol
                 ind1=1+(ipol-1)*npwx
                 ind2=npw+(ipol-1)*npwx
                 DO ibnd=1,nbnd
                    psic(nl(igk_k(1:npw,ik_origin(ik))),ipol,ibnd) = &
                                                       evc(ind1:ind2,ibnd)
                    CALL invfft ('Dense', psic(:,ipol,ibnd), dfftp)
                 ENDDO
              ENDDO
              CALL compute_gk(xk(1,ik), xk(1,ik_origin(ik)), &
                 s(1,1,invs(isym_bands(ik))), t_rev(invs(isym_bands(ik))), gk)
              IF (noncolin) THEN
                 has_e=1
                 CALL find_u(sr(1,1,isym_bands(ik)),d_spin) 
              ENDIF
              CALL rotate_all_psi_tpw(ik,psic,evcr,s(1,1,invs(isym_bands(ik))),&
                     ftau(1,invs(isym_bands(ik))), d_spin, has_e, gk)
              CALL save_buffer ( evcr, nwordwfc, iunwfc, ik )
              et(1:nbnd,ik)=et(1:nbnd,ik_origin(ik))
              DEALLOCATE(evcr)
              DEALLOCATE(psic)
           ELSE
             CALL errore('c_bands','Problem the code should not arrive here',1)
           ENDIF
        ENDIF
     END DO k_loop
  ENDDO
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  CALL mp_sum( ik_diago, inter_pool_comm )
  CALL mp_sum( ik_sym, inter_pool_comm )
  avg_iter = avg_iter / ik_diago
  !
  WRITE( stdout, '(/,5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)' ) &
       ethr, avg_iter

  WRITE( stdout, '(/,5X,"Total points",I5)' ) nkstot
  WRITE( stdout, '(5X,"Diagonalized points",I5,3X,"Symmetrized points", I5)' )&
                                              ik_diago, ik_sym
  IF (tmp_dir /= tmp_dir_save) CALL close_buffer(iuawfc,'keep')
  !
  CALL stop_clock( 'c_bands' )
  !
  RETURN
  !
  ! formats
  !
9001 FORMAT(/'     Computing kpt #: ',I5 )
9000 FORMAT( '     total cpu time spent up to now is ',F10.1,' secs' )
  !
END SUBROUTINE c_bands_nscf_tpw

SUBROUTINE compute_gk(xk, xk_orig, s, t_rev, gk)

USE kinds, ONLY : DP
USE cell_base, ONLY : at, bg
IMPLICIT NONE

REAL(DP) :: xk(3), xk_orig(3)
INTEGER :: s(3,3), t_rev, gk(3)

REAL(DP) :: xkc(3), xko(3), xks(3)
INTEGER  :: kpol

xkc=xk
xko=xk_orig
CALL cryst_to_cart(1,xkc,at,-1)
CALL cryst_to_cart(1,xko,at,-1)

DO kpol=1,3
   xks(kpol)=s(kpol,1)*xko(1)+s(kpol,2)*xko(2)+s(kpol,3)*xko(3)
END DO
IF (t_rev==1) xks = - xks

gk = NINT(xks - xkc)

RETURN
END SUBROUTINE compute_gk
