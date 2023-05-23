!
! Copyright (C) 2001-2003 PWSCF group
! Copyright (C) 2023 Andrea Dal Corso (generalization to many k)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define TEST_NEW_PRECONDITIONING
!
!-----------------------------------------------------------------------
SUBROUTINE g_psii( lda, n, m, npol, psi, e, ik )
  !-----------------------------------------------------------------------
  !! This routine computes an estimate of the inverse Hamiltonian
  !! and applies it to m wavefunctions.
  !
  USE kinds
  USE many_k_mod,     ONLY : h_diagk_d, s_diagk_d
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: lda
  !! input: the leading dimension of psi
  INTEGER,intent(in) :: n
  !! input: the real dimension of psi
  INTEGER,intent(in) :: m
  !! input: the number of coordinates of psi
  INTEGER, intent(in) :: npol
  !! input: the number of bands
  INTEGER, INTENT(IN) :: ik
  !! input: the k point
  COMPLEX(DP) :: psi(lda, npol, m)
  !! inp/out: the psi vector
  REAL(DP), intent(in) :: e(m)
  !! input: the eigenvectors
  !
  !  ... local variables
  !
  INTEGER :: ipol
  ! counter of coordinates of psi
  REAL(DP), PARAMETER :: eps = 1.0d-4
  ! a small number
  REAL(DP) :: x, scala, denm
  !
  INTEGER :: k, i
  ! counter on psi functions
  ! counter on G vectors
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: iblock, numblock
  ! chunking parameters
  !
#if ! defined (__CUDA)
  CALL start_clock( 'g_psi' )
  !
  ! compute the number of chuncks
  numblock  = (n+blocksize-1)/blocksize
  !
#ifdef TEST_NEW_PRECONDITIONING
  scala = 1.d0
  !$omp parallel do collapse(3) private(x, denm)
  DO k = 1, m
     DO ipol=1, npol
        DO iblock = 1, numblock
           DO i = (iblock-1)*blocksize+1, MIN( iblock*blocksize, n )
              x = (h_diagk_d(i,ipol,ik) - e(k)*s_diagk_d(i,ipol,ik))*scala
              denm = 0.5_dp*(1.d0+x+SQRT(1.d0+(x-1)*(x-1.d0)))/scala
              psi (i, ipol, k) = psi (i, ipol, k) / denm
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$omp end parallel do
#else
  !$omp parallel do collapse(3) private(denm)
  DO ipol=1,npol
     DO k = 1, m
        DO iblock = 1, numblock
           DO i = (iblock-1)*blocksize+1, MIN( iblock*blocksize, n )
              denm = h_diagk_d(i,ipol,ik) - e(k) * s_diagk_d(i,ipol,ik)
              !
              ! denm = g2+v(g=0) - e(k)
              !
                 IF (ABS(denm) < eps) denm = SIGN( eps, denm )
              !
              ! denm = sign( max( abs(denm),eps ), denm )
              !
              psi(i, ipol, k) = psi(i, ipol, k) / denm
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$omp end parallel do
#endif
  !
  CALL stop_clock( 'g_psi' )
#endif
  !
  RETURN
  !
END SUBROUTINE g_psii
