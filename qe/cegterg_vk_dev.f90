!
! Copyright (C) 2023- Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  Parts of the cegterg routines moved to GPU 
!
#if defined(__CUDA)
!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_init(outk_d, nbasek_d, kdimk_d, &
                    stx_d, hpsi, psi, hc, sc, kdmx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: kdmx, nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: kdimk_d(nk)
INTEGER, DEVICE :: stx_d(nk)
COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)

INTEGER :: ik, nbase, kdim, nb1, notcnv, stx_, i, j, k, i0, j0
COMPLEX(DP) :: asum   

ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN

nbase=nbasek_d(ik)
kdim=kdimk_d(ik)
stx_=stx_d(ik)
!
i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x 
IF (i0>nbase) RETURN
j0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y 
IF (j0>nbase) RETURN
i= i0 + stx_
j= j0 + stx_
asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(psi(k,i)) * hpsi(k,j)
ENDDO
hc(i0,j0,ik)= asum

asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(psi(k,i)) * psi(k,j)
ENDDO
sc(i0,j0,ik)= asum
RETURN
END SUBROUTINE cegterg_init

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_init_us (outk_d, nbasek_d, kdimk_d, &
                      stx_d, hpsi, spsi, psi, hc, sc, kdmx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: kdmx, nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: kdimk_d(nk)
INTEGER, DEVICE :: stx_d(nk)
COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: spsi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)

INTEGER :: ik, nbase, kdim, nb1, notcnv, stx_, i, j, k, i0, j0
COMPLEX(DP) :: asum   

ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN

nbase=nbasek_d(ik)
kdim=kdimk_d(ik)
stx_=stx_d(ik)
!
i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x 
IF (i0>nbase) RETURN
j0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y 
IF (j0>nbase) RETURN
i= i0 + stx_
j= j0 + stx_
asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(psi(k,i)) * hpsi(k,j)
ENDDO
hc(i0,j0,ik)= asum

asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(psi(k,i)) * spsi(k,j)
ENDDO
sc(i0,j0,ik)= asum
RETURN
END SUBROUTINE cegterg_init_us

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_overlap(outk_d, nbasek_d, kdimk_d, &
              notcnvk_d, nb1k_d, stx_d, hpsi, psi, hc, sc, kdmx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: kdmx, nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: kdimk_d(nk)
INTEGER, DEVICE :: notcnvk_d(nk)
INTEGER, DEVICE :: nb1k_d(nk)
INTEGER, DEVICE :: stx_d(nk)
COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)

INTEGER :: ik, nbase, kdim, nb1, notcnv, stx_, i, j, k, i0, j0, i1
COMPLEX(DP) :: asum   

ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN

nbase=nbasek_d(ik)
kdim=kdimk_d(ik)
nb1=nb1k_d(ik)
notcnv=notcnvk_d(ik)
stx_=stx_d(ik)
!
i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x 
IF (i0>notcnv) RETURN
j0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y 
IF (j0>(nbase+notcnv)) RETURN
i= i0 + stx_+ nb1 - 1
j= j0 + stx_
i1= i0 + nb1 - 1  
asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(hpsi(k,i)) * psi(k,j)
ENDDO
hc(i1,j0,ik)= asum

asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(psi(k,i)) * psi(k,j)
ENDDO
sc(i1,j0,ik)= asum

RETURN
END SUBROUTINE cegterg_overlap

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_overlap_us (outk_d, nbasek_d, kdimk_d, &
          notcnvk_d, nb1k_d, stx_d, hpsi, spsi, psi, hc, sc, kdmx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: kdmx, nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: kdimk_d(nk)
INTEGER, DEVICE :: notcnvk_d(nk)
INTEGER, DEVICE :: nb1k_d(nk)
INTEGER, DEVICE :: stx_d(nk)
COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: spsi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)

INTEGER :: ik, nbase, kdim, nb1, notcnv, stx_, i, j, k, i0, j0, i1
COMPLEX(DP) :: asum   

ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN

nbase=nbasek_d(ik)
kdim=kdimk_d(ik)
nb1=nb1k_d(ik)
notcnv=notcnvk_d(ik)
stx_=stx_d(ik)
!
i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x 
IF (i0>notcnv) RETURN
j0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y 
IF (j0>(nbase+notcnv)) RETURN
i= i0 + stx_+ nb1 - 1
j= j0 + stx_
i1= i0 + nb1 - 1  
asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(hpsi(k,i)) * psi(k,j)
ENDDO
hc(i1,j0,ik)= asum

asum=(0.0_DP,0.0_DP)
DO k=1, kdim
   asum = asum + CONJG(spsi(k,i)) * psi(k,j)
ENDDO
sc(i1,j0,ik)= asum
RETURN
END SUBROUTINE cegterg_overlap_us

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE compute_dot_ew (outk_d, npw, nbasek_d, &
                      notcnvk_d, stx_d, psi, ew, npol, npwx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: npwx, nvecx, npol, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk), npw(nk), notcnvk_d(nk), stx_d(nk)
COMPLEX(DP), DEVICE :: psi(npwx*npol,nvecx*nk)
REAL(DP), DEVICE :: ew(nvecx,nk)

INTEGER :: ik, nbase, notcnv, stx_, i, i0, nbn, np
COMPLEX(DP) :: asum   

ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN

np=npw(ik)
nbase=nbasek_d(ik)
notcnv=notcnvk_d(ik)
stx_=stx_d(ik)
!
i0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y 
IF (i0>notcnv) RETURN

nbn = stx_ + nbase + i0
!
asum=(0.0_DP,0.0_DP)
DO i=1, np
   asum=asum+CONJG(psi(i,nbn))*psi(i,nbn)
ENDDO
ew(i0,ik) = asum
IF (npol>1) THEN
   asum=(0.0_DP,0.0_DP)
   DO i=1, np
      asum=asum+CONJG(psi(npwx+i,nbn))*psi(npwx+i,nbn)
   ENDDO
   ew(i0,ik) = ew(i0,ik)+asum
ENDIF

RETURN
END SUBROUTINE compute_dot_ew

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_upd0(outk_d, nbasek_d, st_d, aux_d, &
                    conv_d, nb1k_d, vc, ew, e, nvecx, nvec, kter, nk) 
!-------------------------------------------------------------------------
USE util_param, ONLY : DP
IMPLICIT NONE
INTEGER, VALUE :: nk, nvecx, nvec, kter
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: st_d(nk)
INTEGER, DEVICE :: aux_d(nk)
INTEGER, DEVICE :: nb1k_d(nk)
LOGICAL, DEVICE :: conv_d(nvec*nk)
COMPLEX(DP), DEVICE :: vc(nvecx,nvecx,nk)
REAL(DP), DEVICE :: ew(nvecx,nk), e(nvec*nk)

INTEGER :: ik, np, n, i, nbase, st_

ik=(BlockIdx%x-1) * BlockDim%x + ThreadIdx%x
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN
!
nbase=nbasek_d(ik)
st_=st_d(ik)
!
aux_d(ik) = kter !; write(*,*) kter, notcnv, conv
!
np = 0
!
DO n = 1, nvec
   !
   IF ( .NOT. conv_d(st_+n) ) THEN
      !
      ! ... this root not yet converged ... 
      !
      np = np + 1
      !
      ! ... reorder eigenvectors so that coefficients for unconverged
      ! ... roots come first. This allows to use quick matrix-matrix 
      ! ... multiplications to set a new basis vector (see below)
      !
      IF ( np /= n ) THEN
         DO i = 1, nvecx
            vc(i,np,ik) = vc(i,n,ik)
         END DO 
      END IF 
      !
      ! ... for use in g_psi
      !
      ew(nbase+np,ik) = e(st_+n)
      !
   END IF
   !
END DO
!
nb1k_d(ik) = nbase + 1

RETURN
END SUBROUTINE cegterg_upd0

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_upd3(outk_d, nbasek_d, kdimk_d, &
        notcnvk_d, nb1k_d, stx_d, hpsi, psi, vc, kdmx, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: kdmx, nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: kdimk_d(nk)
INTEGER, DEVICE :: notcnvk_d(nk)
INTEGER, DEVICE :: nb1k_d(nk)
INTEGER, DEVICE :: stx_d(nk)
COMPLEX(DP), DEVICE :: hpsi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: psi(kdmx,nvecx*nk)
COMPLEX(DP), DEVICE :: vc(nvecx,nvecx,nk)

INTEGER :: ik, nbase, kdim, nb1, notcnv, stx_, i, j, k, i0, j0, k0
COMPLEX(DP) :: asum   

ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN

nbase=nbasek_d(ik)
kdim=kdimk_d(ik)
nb1=nb1k_d(ik)
notcnv=notcnvk_d(ik)
stx_=stx_d(ik)
!
i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x 
IF (i0>kdim) RETURN
j0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y 
IF (j0>notcnv) RETURN
j= j0 + stx_+nb1 - 1  
asum=(0.0_DP,0.0_DP)
DO k=1, nbase
   asum = asum + hpsi(i0,k+stx_) * vc(k,j0,ik)
ENDDO
psi(i0,j)= psi(i0,j)+asum

RETURN
END SUBROUTINE cegterg_upd3

!-------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE cegterg_herm(outk_d, nbasek_d, nb1k_d, &
                              hc, sc, nvecx, nk)
!-------------------------------------------------------------------------
USE cudafor
USE util_param, ONLY : DP
IMPLICIT NONE

INTEGER, VALUE :: nvecx, nk
LOGICAL, DEVICE :: outk_d(nk)
INTEGER, DEVICE :: nbasek_d(nk)
INTEGER, DEVICE :: nb1k_d(nk)
COMPLEX(DP), DEVICE :: hc(nvecx,nvecx,nk)
COMPLEX(DP), DEVICE :: sc(nvecx,nvecx,nk)

INTEGER :: ik, nbase, nb1, n, m

ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ik > nk) RETURN
IF (outk_d(ik)) RETURN

nbase=nbasek_d(ik)
nb1=nb1k_d(ik)
!
n=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x 
IF (n>nbase) RETURN

IF ( n>=nb1 ) THEN
   hc(n,n,ik) = CMPLX( REAL( hc(n,n,ik) ), 0.D0 ,KIND=DP)
   sc(n,n,ik) = CMPLX( REAL( sc(n,n,ik) ), 0.D0 ,KIND=DP)
ENDIF

DO m = MAX(n+1,nb1), nbase
   !
   hc(n,m,ik) = CONJG( hc(m,n,ik) )
   sc(n,m,ik) = CONJG( sc(m,n,ik) )
   !
END DO

RETURN
END SUBROUTINE cegterg_herm

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE copy_psi_gpu( lda, outk, enter, kdim, stx, st, &
                              npol, psi_d, evc, nvec, nvecx, nk )
  !-----------------------------------------------------------------------
  !
  ! This routine copies evc on psi for all k points. Note that it copies
  ! a fixed number of bands, nvec, for each k point.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nk, nvec, nvecx
  !! input: the number of k point
  !! input: the number of bands per k
  !! input: the maximum number of bands per k
  LOGICAL, INTENT(IN), DEVICE :: outk(nk), enter(nk)
  !! input: variables that allow to skip some k points
  INTEGER, INTENT(IN), DEVICE :: kdim(nk)
  !! input: the real dimension of psi
  INTEGER, INTENT(IN), DEVICE :: stx(nk), st(nk)
  !! input: the number of coordinates of psi
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of polarizations
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx * nk)
  COMPLEX(DP), DEVICE :: evc(lda*npol, nvec * nk)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, i0, ik, k, k1, n, stx_, st_
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nk) RETURN
  IF (outk(ik).OR..NOT.enter(ik)) RETURN
  n = kdim(ik)
  st_=st(ik)
  stx_=stx(ik)

  i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i0 > n) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>nvec) RETURN
  k = stx_ + k0 
  k1= st_ + k0 

  psi_d(i0,k)=evc(i0,k1)
  !
  RETURN
  !
END SUBROUTINE copy_psi_gpu
!
#endif
