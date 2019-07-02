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
  USE buffers,              ONLY : get_buffer, save_buffer, open_buffer, &
                                   close_buffer
  USE basis,                ONLY : starting_wfc
  USE fft_base,             ONLY : dfftp
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g, gg, ngm
  USE fft_interfaces,       ONLY : invfft
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_flags,        ONLY : ethr, restart, isolve, io_level, iverbosity
  USE save_ph,              ONLY : tmp_dir_save
  USE io_files,             ONLY : tmp_dir, prefix
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions,        ONLY : evc
  USE control_lr,           ONLY : lgamma
  USE mp_asyn,              ONLY : asyn_master, with_asyn_images
  USE mp_images,            ONLY : my_image_id, root_image, intra_image_comm
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm, my_pool_id
  USE mp_orthopools,        ONLY : intra_orthopool_comm, mp_start_orthopools,&
                                   mp_stop_orthopools, me_orthopool
  USE mp,                   ONLY : mp_sum, mp_min, mp_get
  USE io_global,            ONLY : ionode
  USE check_stop,           ONLY : check_stop_now
  USE band_computation,     ONLY : diago_bands, ik_origin, sym_for_diago
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : domag
  !
  IMPLICIT NONE
  !
  REAL(DP) :: avg_iter, ethr_, mem_to_alloc, max_mem, t_trev
  ! average number of H*psi products
  INTEGER :: ik_, ik, ik_eff, ibnd, nkdum, npw, ios, ipol, ind1, ind2, gk(3)

  ! ik_: k-point already done in a previous run
  ! ik : counter on k points
  INTEGER :: ik_diago, ik_sym, ikk, divide, queue
  LOGICAL :: all_done_asyn
  !
  REAL(DP), EXTERNAL :: get_clock
  COMPLEX(DP), ALLOCATABLE :: psic(:,:,:), evcr(:,:,:), evcbuffer(:,:,:,:), &
                              evcrecv(:,:,:,:)
  REAL(DP), ALLOCATABLE :: aux_xk(:,:), aux_et(:,:), mem(:)
  INTEGER, ALLOCATABLE :: working_pool(:), ikd(:), need(:,:), index_send(:,:), &
                          index_recv(:,:), ticket(:), nkdiag_loc(:), nkrecv_loc(:)
  INTEGER, EXTERNAL :: local_kpoint_index, global_kpoint_index
  INTEGER :: ig, iks, ik1, ikg, ikdiag, ikdiag_loc, nkdiag, ipool, iuawfc, & 
             lrawfc
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
     ikk=global_kpoint_index(nkstot,ik)
     IF (diago_bands(ikk)) nkdum=nkdum+1
  ENDDO
  CALL mp_min(nkdum,inter_pool_comm)
  !
  ! ... For each k point (except those already calculated if restarting)
  ! ... diagonalizes the hamiltonian
  !
  ik_diago=0
  ik_sym=0
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     ikk=global_kpoint_index(nkstot,ik)
     IF (diago_bands(ikk)) THEN
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
        !  In the noncolinear magnetic case we have k, k+q, -k -k-q and
        !  to the last two wavefunctions we must apply t_rev.
        !  When lgamma is true we have only k and -k
        !
        IF (noncolin.AND.domag) THEN
           IF (lgamma.AND. MOD(ik,2)==0) THEN
              CALL start_clock( 't_rev' )
              CALL apply_trev(evc, ik, ik-1)
              CALL stop_clock( 't_rev' )
           ELSEIF (.NOT.lgamma.AND.(MOD(ik,4)==3.OR.MOD(ik,4)==0)) THEN
              CALL start_clock( 't_rev' )
              CALL apply_trev(evc, ik, ik-2)
              CALL stop_clock( 't_rev' )
           ENDIF
        ENDIF
        !
        ! ... save wave-functions (unless disabled in input)
        !
        IF ( io_level > -1 ) CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
        !
        !
        IF (ik_diago .LE. nkdum) THEN
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
     ENDIF
  ENDDO k_loop

  IF (sym_for_diago) THEN
     CALL mp_start_orthopools(intra_image_comm)
     !
     !   et and xk collected on all pools
     !
     ALLOCATE(aux_et(nbnd,nkstot))
     ALLOCATE(aux_xk(3,nkstot))
     CALL poolcollect(    3, nks, xk, nkstot, aux_xk)
     CALL poolcollect( nbnd, nks, et, nkstot, aux_et)
!
!   count how many points have been diagonalized and prepare the indices
!   to transfer them among pools.
!
     nkdiag=0
     DO ik=1, nkstot
        IF (diago_bands(ik)) nkdiag=nkdiag+1
     ENDDO
!
!   This part sets: 
!   need, for each pool and each k points is 1 if the pool needs it, or 0. 
!   (all pools have the same info, ipool from 1 to npool). 
!   index_recv, each pool has its own and is the position on the receive
!   vector of the k point ikdiag
!   index_send, each pool has it own and is the position on the send vector
!   of the k point ikdiag
!   working_pool is the pool that has diagonalized (from 1 to npool) ikdiag
!   ikd is the global index of the k point in the nkdiag list (all pools
!   have the same info).
!   ticket is used to divide the sending and receiving array into
!   smaller pieces of a maximum size of 1 GByte
!
!
     ALLOCATE(need(nkdiag,npool))     
     ALLOCATE(working_pool(nkdiag))
     ALLOCATE(ikd(nkdiag))
     ALLOCATE(ticket(nkdiag))

     max_mem=1._DP   ! in Gbytes

     DO divide=1, nkdiag

        working_pool=0
        need=0
        ticket=0

        ALLOCATE(index_send(nkdiag,divide))             
        ALLOCATE(index_recv(nkdiag,divide))     
        ALLOCATE(nkdiag_loc(divide))     
        ALLOCATE(nkrecv_loc(divide))     
        ALLOCATE(mem(divide))

        index_send=0
        index_recv=0
        nkdiag_loc=0
        nkrecv_loc=0

        DO ikdiag=1, nkdiag
           ticket(ikdiag)=INT(dfloat((ikdiag-1))/nkdiag*divide+1)
        END DO

        DO queue=1, divide
           nkdiag=0
           DO ik=1, nkstot
              IF (diago_bands(ik)) THEN
                 nkdiag=nkdiag+1
                 ikd(nkdiag)=ik
                 ikk=local_kpoint_index(nkstot,ik)
                 IF (ikk /=-1) THEN
                    IF (ticket(nkdiag)==queue) THEN
                       nkdiag_loc(queue)=nkdiag_loc(queue)+1
                       index_send(nkdiag,queue)=nkdiag_loc(queue)
                       working_pool(nkdiag)=my_pool_id+1
                    END IF
                 ELSE
                    nks_loop: DO ikk=1,nks
                       IF (ticket(nkdiag)==queue) THEN
                          ikg=global_kpoint_index(nkstot,ikk)
                          IF (.NOT.diago_bands(ikg).AND.ik_origin(ikg)==ik) THEN
                             nkrecv_loc(queue)=nkrecv_loc(queue)+1
                             need(nkdiag, me_orthopool+1)=1
                             index_recv(nkdiag,queue)=nkrecv_loc(queue)
                             EXIT nks_loop
                          ENDIF
                       END IF
                    ENDDO nks_loop
                 ENDIF
              ENDIF
           ENDDO
           mem(queue)=dfftp%nnr*npol*nbnd*&
                      (nkdiag_loc(queue)+nkrecv_loc(queue))*16.0_DP/1.D9
        END DO

        mem_to_alloc=maxval(mem)

        IF (mem_to_alloc<=max_mem) THEN 
           EXIT
        ELSE
           DEALLOCATE(index_recv)     
           DEALLOCATE(index_send)     
           DEALLOCATE(nkdiag_loc)     
           DEALLOCATE(nkrecv_loc)     
        END IF

        DEALLOCATE(mem)
     END DO

     CALL mp_sum(need, intra_orthopool_comm)
     CALL mp_sum(working_pool, intra_orthopool_comm)

     ALLOCATE(psic(dfftp%nnr, npol, nbnd))
     ALLOCATE(evcr(dfftp%nnr, npol, nbnd))

     DO queue=1,divide
        IF (npool>1) THEN
           IF (nkdiag_loc(queue)==0) nkdiag_loc(queue)=1
           IF (nkrecv_loc(queue)==0) nkrecv_loc(queue)=1
           ALLOCATE(evcbuffer(dfftp%nnr,npol,nbnd,nkdiag_loc(queue)))
           ALLOCATE(evcrecv(dfftp%nnr,npol,nbnd,nkrecv_loc(queue)))
           evcbuffer=(0.0_DP, 0.0_DP)
           evcrecv=(0.0_DP, 0.0_DP)
        ENDIF
!
!   Here brings the diagonalized wavefunction in real space and put them 
!   in the evcbuffer if there are pools, otherwise rotate and save them
!
        ikdiag_loc=0
        DO ikdiag=1, nkdiag
           IF (ticket(ikdiag)/=queue) CYCLE
           ikk=local_kpoint_index(nkstot, ikd(ikdiag))
           IF (ikk /=-1) THEN
              ikdiag_loc=ikdiag_loc+1
              CALL get_buffer ( evc, nwordwfc, iunwfc, ikk )
              npw=ngk(ikk)
              psic=(0.0_DP,0.0_DP)
              DO ipol=1, npol
                 ind1=1+(ipol-1)*npwx
                 ind2=npw+(ipol-1)*npwx
                 DO ibnd=1,nbnd
                    psic(dfftp%nl(igk_k(1:npw,ikk)),ipol,ibnd) = &
                         evc(ind1:ind2,ibnd)
                    CALL invfft ('Rho', psic(:,ipol,ibnd), dfftp)
                 ENDDO
              ENDDO
              IF (npool==1) THEN
                 DO ik1 = 1, nks
                    ik=global_kpoint_index(nkstot, ik1)
                    IF (diago_bands(ik)) CYCLE
                    IF (ik_origin(ik)/=ikd(ikdiag)) CYCLE
                    ik_sym=ik_sym+1
                    CALL rotate_and_save_psic(psic, evcr, aux_xk, ik, ik1, &
                         ik_origin(ik))
                    et(1:nbnd,ik1)=aux_et(1:nbnd,ik_origin(ik))
                 ENDDO
              ELSE
                 evcbuffer(:,:,:,ikdiag_loc)=psic(:,:,:)
              ENDIF
           ENDIF
        ENDDO
        
        IF (npool>1) THEN
!
!   pools receive the wavefunction that they need to obtain the
!   rotated ones from the pools that have computed them
!

           DO ikdiag=1,nkdiag
              IF (ticket(ikdiag)/=queue) CYCLE
              DO ipool=1, npool
                 IF (need(ikdiag,ipool)==1) THEN
                    IF (me_orthopool==(working_pool(ikdiag)-1)) THEN
!
!    I am the pool that calculated the wavefunction. Send it to the
!    pool that needs it
!
                       DO ibnd=1,nbnd
                          CALL mp_get(evcbuffer(:,:,ibnd,index_send(ikdiag,queue)), &
                               evcbuffer(:,:,ibnd,index_send(ikdiag,queue)), &
                               me_orthopool, ipool-1, working_pool(ikdiag)-1, &
                               ibnd, intra_orthopool_comm)
                       ENDDO
                    ELSEIF (me_orthopool==(ipool-1)) THEN
!
!   I am the pool that needs the ikdiag wavefunctions so I receive them 
!   here
!
                       DO ibnd=1,nbnd
                          CALL mp_get(evcrecv(:,:,ibnd,index_recv(ikdiag,queue)), &
                               evcrecv(:,:,ibnd,index_recv(ikdiag,queue)), &
                               me_orthopool, ipool-1, working_pool(ikdiag)-1, &
                               ibnd, intra_orthopool_comm)
                       ENDDO
                    ELSE
!
!  I am not interested in ikdiag, but call mp_get due to the presence of
!  a barrier. I do not communicate here. 
!
                       DO ibnd=1,nbnd
                          CALL mp_get(evcrecv(:,:,ibnd,1), evcrecv(:,:,ibnd,1), &
                               me_orthopool, ipool-1, working_pool(ikdiag)-1, &
                               ibnd, intra_orthopool_comm)
                       ENDDO
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
!
!   finally use the received wavefunctions to obtain the rotated ones
!
           DO ikdiag=1, nkdiag
              IF(ticket(ikdiag)/=queue) CYCLE
              DO ikk = 1, nks
                 ik=global_kpoint_index(nkstot, ikk)
                 IF (diago_bands(ik)) CYCLE
                 IF (ik_origin(ik)/=ikd(ikdiag)) CYCLE
                 ik_sym=ik_sym+1
                 IF (need(ikdiag,me_orthopool+1)==1) THEN
                    psic(:,:,:)=evcrecv(:,:,:,index_recv(ikdiag,queue))
                 ELSE
                    psic(:,:,:)=evcbuffer(:,:,:,index_send(ikdiag,queue))
                 ENDIF
                 CALL rotate_and_save_psic(psic, evcr, aux_xk, ik, ikk, &
                      ik_origin(ik))
                 et(1:nbnd,ikk)=aux_et(1:nbnd,ik_origin(ik))
              ENDDO
           ENDDO
        ENDIF
        IF (npool>1) THEN
           DEALLOCATE(evcbuffer)
           DEALLOCATE(evcrecv)
        ENDIF
     END DO

     CALL mp_stop_orthopools()
     !
     DEALLOCATE(evcr)
     DEALLOCATE(psic)
     DEALLOCATE(ikd)
     DEALLOCATE(aux_et)
     DEALLOCATE(aux_xk)
     DEALLOCATE(working_pool)
     DEALLOCATE(index_recv)     
     DEALLOCATE(index_send)
     DEALLOCATE(nkdiag_loc)     
     DEALLOCATE(nkrecv_loc)          
     DEALLOCATE(need)
     DEALLOCATE(ticket)
     !
  ENDIF
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
 
  CALL print_clock( 't_rev' )
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
!
!   Note that rotate_all_psi applies only the symmetry without 
!   time reversal. So instead of gk we have to give here -gk.
!   Here we cannot ignore completely time reversal because only TSk=k+G.
!   So to apply completely the symmetry one has to use apply_trev on the
!   output of rotate_all_psi calculated with -G.
!
IF (t_rev==1) gk=-gk

RETURN
END SUBROUTINE compute_gk

SUBROUTINE rotate_and_save_psic(psic, evcr, aux_xk, ik, ikk, iko) 
!
!  Input variables: psic with the wavefunction to rotate in real space
!  evcr : where the rotated function is written in reciprocal space
!  aux_xk all the list of k points for all pools
!  ik global index of the current k point
!  iko global index of the k point that rotated gives the current k point
!  ikk local index of the current k point

USE kinds,      ONLY : DP
USE symm_base,  ONLY : s, sr, ftau, invs, t_rev
USE fft_base,   ONLY : dfftp
USE klist,      ONLY : xk, nkstot, ngk, igk_k
USE wvfct,      ONLY : nbnd, npwx
USE wavefunctions, ONLY : evc
USE fft_interfaces,       ONLY : fwfft
USE io_files,          ONLY : iunwfc, nwordwfc
USE band_computation,  ONLY : isym_bands
USE noncollin_module,  ONLY : noncolin, npol
USE control_lr,        ONLY : lgamma
USE spin_orb,          ONLY : domag
USE buffers,           ONLY : save_buffer

IMPLICIT NONE

INTEGER :: ik, iko, ikk
REAL(DP) :: aux_xk(3,nkstot)
COMPLEX(DP) :: psic(dfftp%nnr, npol, nbnd), evcr(dfftp%nnr, npol, nbnd)
COMPLEX(DP) :: d_spin(2,2)
INTEGER :: has_e, ibnd, ipol, ig, iks, npw, gk(3)

CALL compute_gk(xk(1,ikk), aux_xk(1,iko), &
          s(1,1,invs(isym_bands(ik))), t_rev(invs(isym_bands(ik))), gk)

IF (noncolin) THEN
   has_e=1
   CALL find_u(sr(1,1,isym_bands(ik)),d_spin) 
ENDIF

iks=0
IF (noncolin.AND.domag) THEN
   IF (t_rev(isym_bands(ik))==0) THEN
      !
      ! In this case the symmetry has no time reversal, but
      ! we need the time reversed wavefunctions for -k and -k-q
      !
      CALL rotate_all_psi_r_tpw(psic,evcr, s(1,1,invs(isym_bands(ik))), &
                           ftau(1,invs(isym_bands(ik))), d_spin, has_e, gk)
      IF (lgamma.AND. MOD(ik,2)==0) THEN
         CALL start_clock( 't_rev' )
         CALL apply_trev_r(evcr)
         CALL stop_clock( 't_rev' )
         iks=-1
      ELSEIF (.NOT.lgamma.AND.(MOD(ik,4)==3.OR.MOD(ik,4)==0)) THEN
         CALL start_clock( 't_rev' )
         CALL apply_trev_r(evcr)
         CALL stop_clock( 't_rev' )
         iks=-2
      ENDIF
   ELSE
   !
   ! In this case the symmetry has time reversal, and
   ! we need the time reversed wavefunctions for -k and -k-q
   ! so for these case we do not apply any time reversal
   ! but must use the k indeces of k and k+q
   ! for obtaining k+q we must apply time reversal because
   ! the symmetry has it
   !
      CALL rotate_all_psi_r_tpw(psic,evcr, s(1,1,invs(isym_bands(ik))),&
                            ftau(1,invs(isym_bands(ik))), d_spin, has_e, gk)
      IF (lgamma.AND.MOD(ik,2)==0) THEN
         iks=-1
      ELSEIF (lgamma.AND.MOD(ik,2)==1) THEN
         CALL apply_trev_r(evcr)
      ELSEIF (.NOT.lgamma.AND.(MOD(ik,4)==1.OR.MOD(ik,4)==2)) THEN
         CALL apply_trev_r(evcr)
      ELSEIF (.NOT.lgamma.AND.(MOD(ik,4)==3.OR.MOD(ik,4)==0)) THEN
         iks=-2
      ENDIF
   ENDIF
ELSE
   CALL rotate_all_psi_r_tpw(psic, evcr, s(1,1,invs(isym_bands(ik))),&
                   ftau(1,invs(isym_bands(ik))), d_spin, has_e, gk)
ENDIF

DO ibnd=1,nbnd
   DO ipol=1,npol
      CALL fwfft ('Rho', evcr(:,ipol,ibnd), dfftp)
   ENDDO
ENDDO

npw=ngk(ikk+iks)
evc=(0.0_DP,0.0_DP)
DO ibnd=1,nbnd
   DO ig=1, npw
      evc(ig, ibnd)=evcr(dfftp%nl(igk_k(ig,ikk+iks)),1,ibnd)
      IF (npol==2) &
         evc(ig+npwx, ibnd)=evcr(dfftp%nl(igk_k(ig,ikk+iks)),2,ibnd)
   ENDDO
ENDDO
CALL save_buffer ( evc, nwordwfc, iunwfc, ikk )

RETURN
END SUBROUTINE rotate_and_save_psic
