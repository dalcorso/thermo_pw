!
! Copyright (C) 2003-2013 PWSCF group
! Copyright (C) 2023 Andrea Dal Corso (generalization to many k)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psik_k( lda, n, m, psi, v, hpsi, ik )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points:
  !
  !! * fft to real space;
  !! * product with the potential v on the smooth grid;
  !! * back to reciprocal space;
  !! * addition to the hpsi.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE many_k_mod,             ONLY : startkb, current_ikb
  USE mp_bands,               ONLY : me_bgrp
  USE fft_base,               ONLY : dffts
  USE fft_wave
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, tg_get_group_nr3
  USE wavefunctions,          ONLY : psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  INTEGER, INTENT(IN) :: ik
  !! the k point
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  !! Hamiltonian dot psi
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !! the total pot. in real space (smooth grid) for current spin
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr
  INTEGER :: i, iin, right_nnr, right_nr3, right_inc
  COMPLEX(DP), ALLOCATABLE :: vpsi(:,:)
  ! ... chunking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
  ! ... Task Groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:), tg_vpsi(:,:)
  INTEGER :: idx, brange, dffts_nnr
  !
  CALL start_clock( 'vloc_psi' )
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     CALL start_clock( 'vloc_psi:tg_gather' )
     dffts_nnr = dffts%nnr_tg
     incr = fftx_ntgrp(dffts)
     ALLOCATE( tg_v(dffts_nnr) )
     ALLOCATE( tg_psic(dffts_nnr), tg_vpsi(lda,incr) )
     CALL tg_gather( dffts, v, tg_v )
     CALL stop_clock( 'vloc_psi:tg_gather' )
  ELSE
     dffts_nnr = dffts%nnr
     ALLOCATE( vpsi(lda,1) )
  ENDIF
  !
  IF ( use_tg ) THEN
     !
     CALL tg_get_nnr( dffts, right_nnr )
     !
     ! ... compute the number of chuncks
     numblock = (n+blocksize-1)/blocksize
     !
     DO ibnd = 1, m, fftx_ntgrp(dffts)
        !
        CALL tgwave_g2r( psi(:,ibnd:m), tg_psic, dffts, n, igk_k(:,ik) )
        !
!        write (6,*) 'wfc R '
!        write (6,99) (tg_psic(i), i=1,400)
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        !$omp parallel do
        DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
           tg_psic(j) = tg_psic(j) * tg_v(j)
        ENDDO
        !$omp end parallel do
        !
!        write (6,*) 'v psi R ' 
!        write (6,99) (tg_psic(i), i=1,400)
        !
        brange = m-ibnd+1
        !
        CALL tgwave_r2g( tg_psic, tg_vpsi(:,1:brange), dffts, n, igk_k(:,ik) )
        !
        !$omp parallel do collapse(2)
        DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
           DO j = 1, numblock
              DO iin = (j-1)*blocksize+1, MIN(j*blocksize,n)
                 hpsi(iin,ibnd+idx) = hpsi(iin,ibnd+idx) + tg_vpsi(iin,idx+1)
              ENDDO
           ENDDO
        ENDDO
        !$omp end parallel do
        !
     ENDDO
     !
  ELSE
     !
     DO ibnd = 1, m
        !
        CALL wave_g2r( psi(1:n,ibnd:ibnd), psic, dffts, &
             igk=igk_k(:,ik+startkb(current_ikb)) )
        !
!        write (6,*) 'wfc R '
!        write (6,99) (psic(i), i=1,400)
        !
        !$omp parallel do
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v(j)
        ENDDO
        !$omp end parallel do
        !
!        write (6,*) 'v psi R '
!        write (6,99) (psic(i), i=1,400)
        !
        CALL wave_r2g( psic(1:dffts_nnr), vpsi(1:n,:), dffts, &
                     igk=igk_k(:,ik+startkb(current_ikb)) )
        !
        !$omp parallel do
        DO i = 1, n
           hpsi(i,ibnd) = hpsi(i,ibnd) + vpsi(i,1)
        ENDDO
        !$omp end parallel do
        !
!        write (6,*) 'v psi G ', ibnd
!        write (6,99) (psic(i), i=1,400)
        !
     ENDDO
     !
  ENDIF
  !
  IF ( use_tg ) THEN
     DEALLOCATE( tg_psic, tg_vpsi )
     DEALLOCATE( tg_v )
  ELSE
     DEALLOCATE( vpsi )
  ENDIF
  !
  CALL stop_clock( 'vloc_psi' )
  !
99 format ( 20 ('(',2f12.9,')') )
  !
  RETURN
  !
END SUBROUTINE vloc_psik_k
!
