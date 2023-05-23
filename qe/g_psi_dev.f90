!
! Copyright (C) 2001-2003 PWSCF group
! Copyright (C) 2023 Andrea Dal Corso (global version)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define TEST_NEW_PRECONDITIONING
!
#if defined(__CUDA)
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE g_psii_dev( lda, outk, npw, nveck, nb1k, stx, &
                                          npol, psi, e, nvecx, nset )
  !-----------------------------------------------------------------------

  !! This routine computes an estimate of the inverse Hamiltonian
  !! and applies it to nveck wavefunctions for each set. 
  !  The wavefunctions start from nb1k for each set.
  !  The array outk allows to skip completely some set,
  !  while npw is the number of plane waves per set.
  !  Psi contains the wavefunctions that are separated by nvecx 
  !  positions per set. stx gives the starting position in this
  !  array of each set. The routine runs on GPU and is parallelized
  !  on the G vectors on the bands and on the sets.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : h_diagk_d, s_diagk_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  !! input: the number of sets
  !! input: the dimension of each sets in psi
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in psi 
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction.
  COMPLEX(DP), INTENT(INOUT), DEVICE :: psi(lda, npol, nvecx * nset)
  !! inp/out: the psi vector
  REAL(DP), INTENT(IN), DEVICE :: e(nvecx,nset)
  !! input: the current eigenvalues
  !
  !  ... local variables
  !
  INTEGER :: ipol
  ! counter of coordinates of psi
  REAL(DP) :: x, scala, denm
  !
  INTEGER :: k0, k1, ik, i, k, n, m, nb1, stx_
  !
  ! counter on psi functions
  ! counter on G vectors
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  n = npw(ik)
  nb1=nb1k(ik)
  stx_=stx(ik)
  m=nveck(ik) 

  i=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i>n) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>m) RETURN
  k = k0 + stx_ + nb1 - 1
  k1= k0 + nb1 - 1
#ifdef TEST_NEW_PRECONDITIONING
  scala = 1.d0
  DO ipol=1, npol
     x = (h_diagk_d(i,ipol,ik) - e(k1,ik)*s_diagk_d(i,ipol,ik))*scala
     denm = 0.5_dp*(1.d0+x+SQRT(1.d0+(x-1)*(x-1.d0)))/scala
     psi (i, ipol, k) = psi (i, ipol, k) / denm
  ENDDO
#endif
  !
  RETURN
  !
END SUBROUTINE g_psii_dev
!
#endif
