!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! Copyright (C) 2023 Andrea Dal Corso (generalization to many k) 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
#if ! defined (__CUDA)
#define cublasDGEMM DGEMM
#endif
!
!----------------------------------------------------------------------------
SUBROUTINE s_psii_gpu( lda, n, m, psi_d, spsi_d, ik )
  !--------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  !! \(\textit{Wrapper routine}\): performs bgrp parallelization on 
  !! non-distributed bands if suitable and required, calls old S\psi
  !! routine s_psi_ . See comments in h_psi.f90 about band 
  !! parallelization.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE xc_lib,           ONLY : exx_is_active
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_allgather, mp_size, &
                               mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  INTEGER, INTENT(IN) :: ik
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi_d(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi_d(lda*npol,m)
  !! S matrix dot wavefunctions psi
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, spsi_d
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
  !
  CALL start_clock_gpu( 's_psi_bgrp' )
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     ! use band parallelization here
     ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
     CALL divide_all( inter_bgrp_comm,m,m_start,m_end, recv_counts,displs )
     CALL mp_type_create_column_section( spsi_d(1,1), 0, lda*npol, lda*npol, column_type )
     !
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL s_psii__gpu( lda, n, m_end-m_start+1, psi_d(1,m_start), spsi_d(1,m_start), ik )
     CALL mp_allgather(spsi_d, column_type, recv_counts, displs, inter_bgrp_comm )
     !
     CALL mp_type_free( column_type )
     DEALLOCATE( recv_counts )
     DEALLOCATE( displs )
  ELSE
     ! don't use band parallelization here
     CALL s_psii__gpu( lda, n, m, psi_d, spsi_d, ik )
  ENDIF
  !
  CALL stop_clock_gpu( 's_psi_bgrp' )
#endif
  !
  RETURN
  !
END SUBROUTINE s_psii_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE s_psii__gpu( lda, n, m, psi_d, spsi_d, ik )
  !----------------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  USE kinds,            ONLY : DP
  USE uspp,             ONLY : nkb, okvan, ofsbeta, vkb
  USE uspp_param,       ONLY : upf, nh, nhm
  USE many_k_mod,       ONLY : becpk_d, pssk_d, vkbk_d
  USE ions_base,        ONLY : nat, nsp, ityp
  USE control_flags,    ONLY : gamma_only 
  USE noncollin_module, ONLY : npol, noncolin, lspinorb
  USE realus,           ONLY : real_space, invfft_orbital_gamma,     &
                               fwfft_orbital_gamma, calbec_rs_gamma, &
                               s_psir_gamma, invfft_orbital_k,       &
                               fwfft_orbital_k, calbec_rs_k, s_psir_k
  USE wavefunctions,    ONLY : psic
  USE fft_base,         ONLY : dffts
#if defined (__CUDA)
  USE cublas
  USE device_memcpy_m,  ONLY : dev_memcpy
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  INTEGER, INTENT(IN) :: ik
  !! the k point
  COMPLEX(DP), INTENT(IN) :: psi_d(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi_d(lda*npol,m)
  !! S matrix dot wavefunctions psi
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, spsi_d
  !
  COMPLEX(DP), PINNED, ALLOCATABLE :: psi_host(:,:)
  COMPLEX(DP), PINNED, ALLOCATABLE ::spsi_host(:,:)
  !
  INTEGER :: ibnd
  !
  LOGICAL :: need_host_copy
  !
  ! ... initialize  spsi
  !
  !CALL dev_memcpy( spsi_d , psi_d )
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
  !
  need_host_copy = real_space
  IF (need_host_copy) THEN
      ALLOCATE(psi_host(lda*npol,m), spsi_host(lda*npol,m))
      psi_host  = psi_d
      spsi_host = spsi_d
  END IF
  !
  CALL start_clock_gpu( 's_psi' )  
  !
  ! ... The product with the beta functions
  !
  IF ( gamma_only ) THEN
     !
     IF ( real_space ) THEN
        !
        DO ibnd = 1, m, 2
!SdG: the becp are already computed ! no need to invfft psi to real space.
!           CALL invfft_orbital_gamma( psi_host, ibnd, m ) 
!SdG: we just need to clean psic in real space ...
           CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
!SdG: ... before computing the us-only contribution ...
           CALL s_psir_gamma( ibnd, m )
!SdG: ... and add it to spsi (already containing psi).
           CALL fwfft_orbital_gamma( spsi_host, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        spsi_d = spsi_host
        !
     ELSE
        !
!        CALL s_psi_gamma_gpu()
        !
     ENDIF
     !
  ELSEIF ( noncolin ) THEN
     !
     CALL s_psi_nc_gpu()
     !
  ELSE 
     !
     IF ( real_space ) THEN
        !
        DO ibnd = 1, m
!SdG: the becp are already computed ! no need to invfft psi to real space.
!           CALL invfft_orbital_k( psi, ibnd, m )
!SdG: we just need to clean psic in real space ...
           CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
!SdG: ... before computing the us-only contribution ...
           CALL s_psir_k( ibnd, m )
!SdG: ... and add it to spsi (already containing psi).
           CALL fwfft_orbital_k( spsi_host, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        spsi_d = spsi_host
        !
     ELSE
        !
        CALL s_psi_k_gpu()
        !
     ENDIF    
     !
  ENDIF    
  !
  CALL stop_clock_gpu( 's_psi' )
#else
  COMPLEX(DP), ALLOCATABLE :: psi_host(:,:)
  COMPLEX(DP), ALLOCATABLE ::spsi_host(:,:)
#endif
  !
  RETURN
  !
  CONTAINS
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_k_gpu()
       !-----------------------------------------------------------------------
       !! k-points version of \(\textrm{s_psi}\) routine.
       !
       USE device_fbuff_m,   ONLY : dev_buf
       USE uspp,             ONLY : qq_at
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
         ! counters
       COMPLEX(DP), POINTER :: qqc_d(:,:,:)
#if defined(__CUDA)
       attributes(DEVICE) ::  qqc_d
       !
!       pssk_d(1:nkb,1:m,ik) = ( 0.D0, 0.D0 )
       !
       ! qq is real:  copy it into a complex variable to perform
       ! a zgemm - simple but sub-optimal solution
       !
       ! here we need to use qq_at (device) instead of qq_nt_d otherwise real space augmentation brakes!
       !  qq_nt_d would be much faster and works for calculations without real space augmentation
       CALL dev_buf%lock_buffer( qqc_d, (/ nhm, nhm, nat/), ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_k_gpu ', ' cannot allocate buffer (qqc_d) ', ABS(ierr) )

       !$acc parallel loop collapse(3) present(qq_at)
       DO na = 1, nat
          DO jh = 1, nhm
             DO ih = 1, nhm
                qqc_d(ih,jh, na) = CMPLX ( qq_at(ih,jh, na), 0.0_dp, KIND=dp )
             END DO
          END DO
       END DO

       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN

                   CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                     qqc_d(1,1,na), nhm, becpk_d(ofsbeta(na)+1,1,1,ik), nkb, &
                     (0.0_dp,0.0_dp), pssk_d(ofsbeta(na)+1,1,1,ik), nkb )
                   !
                END IF
             END DO
          END IF
       END DO

       CALL dev_buf%release_buffer(qqc_d, ierr)
#endif
       !
       RETURN
       !
     END SUBROUTINE s_psi_k_gpu     
     !
     !
     !-----------------------------------------------------------------------
      SUBROUTINE s_psi_nc_gpu ( )
     !-----------------------------------------------------------------------
       !
       !! k-points noncolinear/spinorbit version of \(\textrm{s_psi}\) routine.
       !
       USE device_fbuff_m,   ONLY : dev_buf
       USE uspp,             ONLY : qq_at, qq_so
       !
       IMPLICIT NONE
       !
       !    here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ipol, ierr
       ! counters
       COMPLEX(DP), POINTER :: ps_d(:,:,:)
       COMPLEX(DP), POINTER :: qqc_d(:,:,:)
#if defined(__CUDA)
       attributes(DEVICE) :: ps_d, qqc_d
       ! the product vkb and psi
       !
       ! sync if needed
       !
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, npol, m /), ierr)
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_nc_gpu ', ' cannot allocate buffer (ps_d) ', ABS(ierr) )

       ps_d(1:nkb,1:npol,1:m) = (0.D0,0.D0)
       !
       IF ( .NOT. lspinorb ) THEN
          CALL dev_buf%lock_buffer( qqc_d, (/ nhm, nhm, nat /), ierr )
          IF( ierr /= 0 .and. ierr /= -1 ) &
             CALL errore( ' s_psi_nc_gpu ', ' cannot allocate buffer (qqc_d) ', ABS(ierr) )
          ! Possibly convert only what's needed??
          !$acc parallel loop collapse(3) present(qq_at)
          DO na = 1, nat
             DO jh = 1, nhm
                DO ih = 1, nhm
                   qqc_d(ih, jh, na) = CMPLX ( qq_at(ih,jh, na), 0.0_dp, KIND=dp )
                END DO
             END DO
          END DO
       END IF
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             !
             IF ( .NOT. lspinorb ) THEN
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      DO ipol=1,npol
                         CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                              qqc_d(1,1, na), nhm, &
                              becpk_d(ofsbeta(na)+1,ipol,1,ik), nkb*npol, &
                              (0.0_dp,0.0_dp), ps_d(ofsbeta(na)+1,ipol,1), nkb*npol )
                       END DO
                    END IF
                END DO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      !$acc host_data use_device(qq_so)
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,1,nt), nhm, becpk_d(ofsbeta(na)+1,1,1,ik) , nkb*npol, &
                           (0.0_dp,0.0_dp), ps_d(ofsbeta(na)+1,1,1), nkb*npol )
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,2,nt), nhm, becpk_d(ofsbeta(na)+1,2,1,ik), nkb*npol, &
                           (1.0_dp,0.0_dp), ps_d(ofsbeta(na)+1,1,1), nkb*npol )
                      !
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,3,nt), nhm, becpk_d(ofsbeta(na)+1,1,1,ik), nkb*npol, &
                           (0.0_dp,0.0_dp), ps_d(ofsbeta(na)+1,2,1), nkb*npol )
                      CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           qq_so(1,1,4,nt), nhm, becpk_d(ofsbeta(na)+1,2,1,ik), nkb*npol, &
                           (1.0_dp,0.0_dp), ps_d(ofsbeta(na)+1,2,1), nkb*npol )
                      !$acc end host_data
                    END IF
                END DO
             END IF
          END IF
       END DO
       IF ( .NOT. lspinorb ) CALL dev_buf%release_buffer(qqc_d, ierr)

!!$acc data present(vkbk_d(:,:))
!!$acc host_data use_device(vkbk_d)
       call ZGEMM ('N', 'N', n, m*npol, nkb, (1.d0, 0.d0), &
          vkbk_d(1,nkb*(ik-1)+1), lda, ps_d, nkb, (1.d0, 0.d0), &
          spsi_d(1,1), lda)
!!$acc end host_data
!!$acc end data

       CALL dev_buf%release_buffer(ps_d, ierr)
#endif

       RETURN

    END SUBROUTINE s_psi_nc_gpu

END SUBROUTINE s_psii__gpu

