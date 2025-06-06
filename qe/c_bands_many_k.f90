! Copyright (C) 2023 Dal Corso Andrea
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE c_bands_many_k(iter)
!--------------------------------------------------------------------------
!
#if defined(__CUDA)
USE cudafor
#endif
USE kinds,            ONLY : DP
USE io_files,         ONLY : iunwfc, nwordwfc
USE control_qe,       ONLY : many_k
USE wvfct,            ONLY : et, nbnd, nbndx, npwx, g2kin, btype, current_k
USE klist,            ONLY : nks, ngk, igk_k, xk, nkstot
USE many_k_mod,       ONLY : evck_d, g2kink_d, vkbk_d, h_diagk_d, s_diagk_d, &
                             nkblocks, nksbx, nksb, startkb, current_ikb, &
                             initialize_fft_factors, deallocate_fft_factors, &
                             evck, alloc_many_k, allocate_becps_many_k, &
                             deallocate_becps_many_k, &
                             initialize_device_variables
USE scf,              ONLY : v_of_0, vrs
USE uspp,             ONLY : nkb, vkb, okvan
USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type
USE fft_base,         ONLY : dffts
USE wavefunctions,    ONLY : evc
USE g_psi_mod,        ONLY : h_diag, s_diag
USE control_flags,    ONLY : ethr, use_gpu, lscf, use_gpu
USE gcscf_module,     ONLY : lgcscf
USE buffers,          ONLY : get_buffer, save_buffer
USE noncollin_module, ONLY : npol
USE io_global,        ONLY : stdout
USE uspp_init,        ONLY : init_us_2
USE uspp_param,       ONLY : lmaxkb, nhm
USE lsda_mod,         ONLY : lsda, current_spin, isk
USE mp,               ONLY : mp_sum
USE mp_pools,         ONLY : inter_pool_comm
USE mp_bands,         ONLY : intra_bgrp_comm

IMPLICIT NONE
#if defined(__CUDA)
#include<init_us_2_interf.f90>
#include<hdiag_interf.f90>
#include<ylm2_interf.f90>
#endif
INTEGER, INTENT(IN)  :: iter
INTEGER, ALLOCATABLE :: nhpsi(:), dav_iter(:), notcnv(:), ikt(:)
INTEGER  :: ik, ik1, ikb, i, j, ipol, ntry, ierr
REAL(DP) :: avg_iter, time, scnds
LOGICAL  :: lrot
EXTERNAL h_psii, s_psii, g_psii
EXTERNAL h_psii_gpu, s_psii_gpu, g_psii_dev
#if defined(__CUDA)
REAL(DP), DEVICE, ALLOCATABLE :: et_d(:,:)
REAL(DP), DEVICE, ALLOCATABLE :: ylm_d(:,:,:)
REAL(DP), DEVICE, ALLOCATABLE :: vkb1_d(:,:,:)
INTEGER, DEVICE, ALLOCATABLE :: ikt_d(:)
#endif

CALL start_clock( 'c_bands' )

ALLOCATE( nhpsi( nksbx ))
ALLOCATE( dav_iter( nksbx ))
ALLOCATE( notcnv( nksbx ))
ALLOCATE( ikt( nksbx ))
#if defined(__CUDA)
   ALLOCATE( et_d( nbnd, nks ))
   ALLOCATE(ikt_d( nksbx ))
   et_d(1:nbnd,1:nks)=et(1:nbnd,1:nks)
#endif

CALL allocate_becps_many_k(1,1)
CALL initialize_fft_factors(1,1)
CALL initialize_device_variables()



lrot=(iter==1)
avg_iter=0.0_DP
DO ikb=1,nkblocks
   current_ikb=ikb
   IF (.NOT.alloc_many_k.OR.nkblocks>1) THEN
      DO ik = startkb(ikb)+1, startkb(ikb)+nksb(ikb)
         ik1=ik-startkb(ikb)
         ikt(ik1)=ik
      ENDDO
      alloc_many_k=.TRUE.
      DO ik = startkb(ikb)+1, startkb(ikb)+nksb(ikb)
         ik1=ik-startkb(ikb)
         IF ( nks > 1 .AND. lscf) THEN
            CALL get_buffer ( evck(1,nbnd*(ik1-1)+1), &
                  nwordwfc, iunwfc, ik )
         ELSEIF (.NOT.lscf) THEN
           !
           current_k = ik
           !
           IF ( lsda ) current_spin = isk(ik)
           !
           CALL g2_kin( ik )
           !
           ! ... More stuff needed by the hamiltonian: nonlocal projectors
           !
           IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), &
                                                    xk(1,ik), vkb, .true.)
           CALL init_wfc ( ik )

           evck(1:npwx*npol,nbnd*(ik1-1)+1:nbnd*ik1)=evc(1:npwx*npol,1:nbnd)
           !
         ELSEIF (nks==1) THEN
           evck(1:npwx*npol,nbnd*(ik1-1)+1:nbnd*ik1)=evc(1:npwx*npol,1:nbnd)
         ENDIF
      ENDDO
      evck_d=evck
#if defined(__CUDA)
      ikt_d=ikt
      ALLOCATE(ylm_d((lmaxkb + 1) **2, npwx, nksb(ikb)))
      ALLOCATE(vkb1_d(nhm, npwx, nksb(ikb)))
      IF ( nkb > 0) CALL ylm2_dev<<<dim3(nksb(ikb),npwx,1),&
               dim3(1,1,1)>>>(ylm_d, lmaxkb, npwx, nksb(ikb), ikt_d )
      ierr=cudaDeviceSynchronize()
      IF ( nkb > 0) CALL init_us_2_kernel<<<dim3(nksb(ikb),npwx,1),&
                     dim3(1,1,1)>>>(vkbk_d, ylm_d, vkb1_d, nhm, lmaxkb, &
                              nkb, npwx, nksb(ikb), ikt_d)
      ierr=cudaDeviceSynchronize()
      DEALLOCATE(vkb1_d)
      DEALLOCATE(ylm_d)
      CALL hdiag_kernel<<<dim3(nksb(ikb),npwx,1),dim3(1,1,1)>>>&
             (g2kink_d, h_diagk_d, s_diagk_d, npwx, v_of_0, npol, nkb, &
             nksb(ikb), ikb )
      ierr=cudaDeviceSynchronize()
#else
      DO ik = startkb(ikb)+1, startkb(ikb)+nksb(ikb)
         ik1=ik-startkb(ikb)
         CALL g2_kin( ik )
         DO i=1,ngk(ik)
            g2kink_d(i,ik1)=g2kin(i)
         ENDDO
         IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), &
                                                xk(1,ik), vkb, .TRUE. )
         DO i=1,npwx
            vkbk_d(i,nkb*(ik1-1)+1:nkb*ik1)=vkb(i,1:nkb)
         ENDDO
         DO j=1, ngk(ik)
            DO ipol=1, npol
               h_diagk_d(j, ipol, ik1) = g2kink_d(j,ik1) + v_of_0
            ENDDO
         END DO
         !
         CALL usnldiag( ngk(ik), h_diagk_d(:,:,ik1), s_diagk_d(:,:,ik1))
      ENDDO
#endif
   ENDIF
!
   ntry=0

   david_loop: DO
      IF (use_gpu) THEN
#if defined(__CUDA)
         CALL cegterg_vk ( h_psii_gpu, s_psii_gpu, okvan, g_psii_dev, &
         ngk(startkb(ikb)+1), npwx, nbnd, nbndx, npol, evck_d, ethr, &
         et_d(:,startkb(ikb)+1), btype(:,startkb(ikb)+1), notcnv, lrot, &
         dav_iter, nhpsi, nksb(ikb), nkb)
         et(:,startkb(ikb)+1:startkb(ikb)+nksb(ikb))= &
         et_d(:,startkb(ikb)+1:startkb(ikb)+nksb(ikb))
#endif
      ELSE
         CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm )
         CALL cegterg_vk ( h_psii, s_psii, okvan, g_psii, &
              ngk(startkb(ikb)+1), npwx, nbnd, nbndx, npol, evck_d, ethr, &
              et(:,startkb(ikb)+1), btype(:,startkb(ikb)+1), notcnv, lrot, &
              dav_iter, nhpsi, nksb(ikb), nkb)
         CALL deallocate_bec_type( becp )
      ENDIF
      ntry = ntry + 1
      !
      ! ... exit condition
      !
      IF ( test_exit_cond(nksb(ikb)) ) EXIT david_loop
   ENDDO david_loop

   evck(1:npwx*npol,1:nbnd*nksb(ikb))=evck_d(1:npwx*npol,1:nbnd*nksb(ikb))
   DO ik = startkb(ikb)+1, startkb(ikb)+nksb(ikb)
      ik1=ik-startkb(ikb)
      IF ( nks > 1 ) THEN
         CALL save_buffer (evck(1,nbnd*(ik1-1)+1), nwordwfc, iunwfc, ik )
      ELSE
         evc(1:npwx*npol,1:nbnd)=evck(1:npwx*npol,nbnd*(ik1-1)+1:nbnd*ik1)
      ENDIF
   ENDDO

   DO ik=1,nksb(ikb)
      avg_iter=avg_iter+dav_iter(ik)
   ENDDO
ENDDO
!$acc update device(et)

CALL mp_sum( avg_iter, inter_pool_comm )
avg_iter = avg_iter / nkstot

CALL print_clock( 'cegterg' )
CALL print_clock( 'cegterg:init' )
CALL print_clock( 'cegterg:diag' )
CALL print_clock( 'cegterg:update' )
CALL print_clock( 'cegterg:overlap' )
CALL print_clock( 'cegterg:last' )
CALL print_clock( 'h_psi_dev' )
!
WRITE( stdout, &
     '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
     ethr, avg_iter
!
CALL stop_clock( 'c_bands' ); !write (*,*) 'stop c_bands' ; FLUSH(6)

CALL deallocate_fft_factors()
CALL deallocate_becps_many_k()

DEALLOCATE( notcnv )
DEALLOCATE( dav_iter )
DEALLOCATE( nhpsi )
DEALLOCATE( ikt )
#if defined(__CUDA)
   DEALLOCATE( ikt_d )
   DEALLOCATE( et_d )
#endif

RETURN

CONTAINS
   !-----------------------------------------------------------------------
  FUNCTION test_exit_cond(nksb)
    !-----------------------------------------------------------------------
    !! This logical function is .TRUE. when iterative diagonalization
    !! is converged.
    !
    IMPLICIT NONE
    !
    LOGICAL :: test_exit_cond
    INTEGER :: ik, notconv, nksb
    !
    notconv=0
    DO ik=1,nksb
       notconv=notconv+notcnv(ik)
    ENDDO
    IF ( lscf .AND. lgcscf ) THEN
       !
       ! ... tight condition for GC-SCF
       !
       test_exit_cond = .NOT. ( ( ntry <= 8 ) .AND. ( notconv > 0 ) )
       !
    ELSE
       !
       test_exit_cond = .NOT. ( ( ntry <= 5 ) .AND. &
            ( ( .NOT. lscf .AND. ( notconv > 0 ) ) .OR. &
            (       lscf .AND. ( notconv > 5 ) ) ) )
       !
    END IF
    !
  END FUNCTION test_exit_cond


END SUBROUTINE c_bands_many_k
