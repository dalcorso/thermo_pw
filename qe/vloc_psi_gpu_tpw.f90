!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! Workaround for missing interface
#if ! defined (__CUDA)
#define tg_gather_gpu tg_gather
#endif
!-----------------------------------------------------------------------
SUBROUTINE vloc_psii_k_gpu( lda, n, m, psi_d, v_d, hpsi_d, ik )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points. GPU double.
  !
  !   fft to real space
  !   product with the potential v on the smooth grid
  !   back to reciprocal space
  !   addition to the hpsi
  !
  USE parallel_include
  USE kinds,         ONLY : DP
  USE wvfct,         ONLY : current_k
  USE klist,         ONLY : igk_k
  USE mp_bands,      ONLY : me_bgrp
  USE control_flags, ONLY : many_fft
  USE many_k_mod,    ONLY : startkb, current_ikb
  USE fft_base,      ONLY : dffts
  USE fft_wave
  USE fft_helper_subroutines
#if defined(__CUDA)
  USE device_fbuff_m,    ONLY : dev_buf
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m, ik
  COMPLEX(DP), INTENT(IN) :: psi_d(lda,m)
  COMPLEX(DP), INTENT(INOUT):: hpsi_d(lda,m)
  REAL(DP), INTENT(IN) :: v_d(dffts%nnr)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, v_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: ibnd, ebnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic(:), vpsi(:,:)
  ! ... Task Groups
  LOGICAL :: use_tg
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:), tg_vpsi(:,:)
  REAL(DP), POINTER :: tg_v_d(:)
  REAL(DP) :: scnds, time, tot_time_fft, tot_time_vloc
  COMPLEX(DP) :: bsum
  REAL(DP) :: asum
#if defined(__CUDA)
  attributes(DEVICE) :: tg_v_d
  !
  INTEGER :: dffts_nnr, idx, group_size, hm_vec(3)
  INTEGER :: ierr, brange
  !
  CALL start_clock_gpu ('vloc_psi')
  use_tg = dffts%has_task_groups
  !
  ALLOCATE( psi(lda,m) )
  !$acc data create( psi )
  !$acc kernels
  psi = psi_d
  !$acc end kernels
  !
  incr = many_fft
  !
  IF( use_tg ) THEN
     CALL start_clock_gpu ('vloc_psi:tg_gather')
     dffts_nnr =  dffts%nnr_tg
     incr = fftx_ntgrp(dffts)
     CALL dev_buf%lock_buffer( tg_v_d, dffts_nnr, ierr )
     ALLOCATE( tg_psic(dffts_nnr), tg_vpsi(dffts_nnr,incr) )
     CALL tg_gather_gpu( dffts, v_d, tg_v_d )
     CALL stop_clock_gpu ('vloc_psi:tg_gather')
  ELSE
     dffts_nnr = dffts%nnr
     ALLOCATE( psic(dffts_nnr*incr), vpsi(dffts_nnr,incr) )
  ENDIF
  !
  IF( use_tg ) THEN
     !
     CALL tg_get_nnr( dffts, right_nnr )
     !
     !$acc data create(tg_psic,tg_vpsi)
     DO ibnd = 1, m, fftx_ntgrp(dffts)
        !
        CALL tgwave_g2r( psi(1:n,ibnd:m), tg_psic, dffts, n, igk_k(:,ik) )
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        !$acc parallel loop
        DO j = 1, dffts%nr1x*dffts%nr2x* right_nr3
           tg_psic(j) = tg_psic(j) * tg_v_d(j)
        ENDDO
        !
        brange = m-ibnd+1
        !
        CALL tgwave_r2g( tg_psic, tg_vpsi(:,1:brange), dffts, n, igk_k(:,ik) )
        !
        DO idx = 1, fftx_ntgrp(dffts)
           IF ( idx+ibnd-1 <= m ) THEN
              !$acc parallel loop
              DO j = 1, n
                 hpsi_d(j,ibnd+idx-1) = hpsi_d(j,ibnd+idx-1) + tg_vpsi(j,idx-1)
              ENDDO
           ENDIF
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ELSEIF (many_fft > 1) THEN
     !
     !$acc data create(psic,vpsi)
     DO ibnd = 1, m, incr
        !
        group_size = MIN(many_fft,m-(ibnd-1))
        hm_vec(1)=group_size ; hm_vec(2)=n ; hm_vec(3)=group_size
        ebnd = ibnd+group_size-1
        !
        call cpu_time(time)
        CALL wave_g2r( psi(:,ibnd:ebnd), psic, dffts, &
                 igk=igk_k(:,ik+startkb(current_ikb)), howmany_set=hm_vec )
        call cpu_time(scnds)
        time=scnds-time
        tot_time_fft=tot_time_fft+time
        !
        call cpu_time(time)
        CALL start_clock('hpsi:psic')

        !$acc parallel loop collapse(2)
        DO idx = 0, group_size-1
           DO j = 1, dffts_nnr
              psic(idx*dffts_nnr+j) = psic(idx*dffts_nnr+j) * v_d(j)
           ENDDO
        ENDDO
        CALL stop_clock('hpsi:psic')
        call cpu_time(scnds)
        time=scnds-time
        tot_time_vloc=tot_time_vloc+time
        !
        call cpu_time(time)
        CALL wave_r2g( psic, vpsi, dffts, &
                 igk=igk_k(:,ik+startkb(current_ikb)), howmany_set=hm_vec )
        call cpu_time(scnds)
        time=scnds-time
        tot_time_fft=tot_time_fft+time
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, group_size-1
           DO j = 1, n
              hpsi_d(j,ibnd+idx) = hpsi_d(j,ibnd+idx) + vpsi(j,idx+1)
           ENDDO
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ELSE
     !
     !$acc data create(psic,vpsi)
     DO ibnd = 1, m
        !
        CALL wave_g2r( psi(1:n,ibnd:ibnd), psic, dffts, &
                                      igk=igk_k(:,ik+startkb(current_ikb)) )
        !
        CALL start_clock('hpsi:psic')
        !$acc parallel loop
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v_d(j)
        ENDDO
        CALL stop_clock('hpsi:psic')
        !
        CALL wave_r2g( psic, vpsi(1:n,:), dffts, &
                                   igk=igk_k(:,ik+startkb(current_ikb)) )
        !
        !$acc parallel loop
        DO i = 1, n
           hpsi_d(i,ibnd) = hpsi_d(i,ibnd) + vpsi(i,1)
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psi )
  !
  IF( use_tg ) THEN
     DEALLOCATE( tg_psic, tg_vpsi )
     CALL dev_buf%release_buffer( tg_v_d, ierr )
  ELSE
     DEALLOCATE( psic, vpsi )
  ENDIF
  !
  CALL stop_clock_gpu( 'vloc_psi' )
#endif
  !
99 format ( 20 ('(',2f12.9,')') )
  !
  RETURN
  !
END SUBROUTINE vloc_psii_k_gpu
!
