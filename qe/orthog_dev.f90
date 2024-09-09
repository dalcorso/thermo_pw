!
! Copyright (C) 2024 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
SUBROUTINE orthogonalize_dev(st_d, outk_d, kdimk_d, npwk_d, nveck_d, &
                   nb1k_d, ikblk_d, nbndk_d, ikb, nk, npe, nsolv,    &
                   dvpsik_d, evqk_d, sevqk_d, ortho_ps_d, npol, npwx,  &
                   nbnd, nksbx_ph)

USE cudafor
USE util_param,     ONLY : DP
IMPLICIT NONE
#include<orthog_interf.f90>
INTEGER, INTENT(IN) :: npwx, npol, ikb, nk, npe, nsolv, nbnd, nksbx_ph
!! input: the leading dimension of psi
!! input: the number of components of each wavefunctions
!! input: the current block of k vectors
!! input: the number of k vectors in the block
!! input: the number of perturbations
!! input: the number of linear systems solved
!! input: the number of bands in each k
!! input: the maximum number of k in each block
INTEGER, INTENT(IN), DEVICE :: st_d(nk*npe*nsolv)
!! input: start of each set of functions
LOGICAL, INTENT(IN), DEVICE :: outk_d(nk*npe*nsolv)
!! input: when .TRUE. the set is not calculated
INTEGER, INTENT(IN), DEVICE :: kdimk_d(nk*nsolv), npwk_d(nk*nsolv)
!! input: kdim is npwx*npol in the noncollinear case, or equal to npwk_d
!!        in the collinear case
!! input: the number of plane waves for each set
INTEGER, INTENT(IN), DEVICE :: ikblk_d(nk*nsolv)
!! input: index of k in the current block of k points
INTEGER, INTENT(IN), DEVICE :: nveck_d(nk*nsolv)
!! input: the number of vectors to compute.
INTEGER, INTENT(IN), DEVICE :: nb1k_d(nk*npe*nsolv)
!! input: the starting point of any set
INTEGER, INTENT(IN), DEVICE :: nbndk_d(nk*npe*nsolv)
!! input: the number of bands to compute for each set
COMPLEX(DP), INTENT(IN), DEVICE :: evqk_d(npwx*npol, nbnd*nk*nsolv)
!! inp: the psi vector
COMPLEX(DP), INTENT(INOUT), DEVICE :: sevqk_d(npwx*npol, nbnd*nk*nsolv)
!! out: the spsi vector
COMPLEX(DP), INTENT(INOUT), DEVICE :: dvpsik_d(npwx*npol, nbnd*nk*npe*nsolv)
!! the functions to orthogonalize
COMPLEX(DP), INTENT(INOUT), DEVICE :: ortho_ps_d(nksbx_ph*nbnd*nsolv, &
                                               nbnd*nksbx_ph*npe*nsolv)
!
!  Local variables
!
INTEGER :: ierr

CALL start_clock('ortho_dev')
CALL orthog_ps<<<dim3(nk*npe*nsolv,nbnd,nbnd),&
        dim3(1,1,1)>>>( st_d, nbndk_d, ikb, nk, npe, nsolv, &
                   dvpsik_d, evqk_d, ortho_ps_d, npol, npwx, nbnd, nksbx_ph )
ierr=cudaDeviceSynchronize()
CALL stop_clock('ortho_dev')
CALL start_clock('ortho_spsi')
CALL s_psik_dev(npwx, outk_d, kdimk_d, npwk_d, nveck_d, nb1k_d, &
        st_d, st_d, ikblk_d, npol, evqk_d, sevqk_d,  &
        nbnd, nbnd, nk*nsolv)
ierr=cudaDeviceSynchronize()
CALL stop_clock('ortho_spsi')
CALL start_clock('ortho_last')
CALL orthog_last<<<dim3(nk*npe*nsolv,nbnd,npwx/32+1),&
     dim3(1,1,32)>>>( st_d, nbndk_d, ikb, nk, npe, nsolv, &
              dvpsik_d, sevqk_d, ortho_ps_d, npol, npwx, nbnd, nksbx_ph )
ierr=cudaDeviceSynchronize()
CALL stop_clock('ortho_last')

RETURN
END SUBROUTINE orthogonalize_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE orthog_ps( st, nbndk, ikb, nk, npe, nsolv, &
                     dvpsik_d, evqk_d, ortho_ps_d, npol, npwx, nbnd, nksbx_ph )
!-----------------------------------------------------------------------
  !
  !  This routine computes the scalar product of the vector evqk_d and
  !  dvpsik_d.
  !  Each thread computes one k point, one perturbation, and one band
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d
  USE wvfct_gpum,     ONLY : et => et_d
  USE many_k_mod,     ONLY : lgauss => lgauss_d, ltetra => ltetra_d,    &
                             degauss=> degauss_d, ngauss => ngauss_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d, ikmks=>ikmks_d,    &
                             ikmkmqs=>ikmkmqs_d, startkb_ph => startkb_ph_d,&
                             nbnd_occ => nbnd_occ_d, alpha_pv => alpha_pv_d,&
                             ef => ef_d
  !
  IMPLICIT NONE
#include<wgauss_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, npwx, nbnd, ikb, nksbx_ph
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of plane waves
  !! input: the number of bands
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: the number of bands to compute for each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  COMPLEX(DP), DEVICE, INTENT(IN) :: evqk_d(npwx*npol, nbnd*nk*nsolv)
  !! inp/out: the bloch wavefunctions
  COMPLEX(DP), DEVICE, INTENT(IN) :: dvpsik_d(npwx*npol, nbnd*nk*npe*nsolv)
  !! the functions to orthogonalize
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: ortho_ps_d(nksbx_ph*nbnd*nsolv, &
                                                 nbnd*nksbx_ph*npe*nsolv)
  !! the scalar products

  INTEGER :: isp, isolv, ipert, ik1, ik, st_, mb, ikmk, ikmkmq, ikk, ikq, &
             k2, k1, npwq, nbnd_eff, id, ibnd, jbnd, ig, ikwf
  REAL(DP) :: wg1, w0g, wgp, deltae, theta, wwg
  COMPLEX(DP) :: aux
  !
  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe*nsolv) RETURN
  isp=(id-1)/nk+1
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1
  ik1=MOD(id-1,nk)+1
  st_=st(id)
  mb=nbndk(id)
  ik=ik1+startkb_ph(ikb)
  IF (isolv==2) THEN
     ikmk = ikmks(ik)
     ikmkmq = ikmkmqs(ik)
  ELSE
     ikmk=ikks(ik)
     ikmkmq=ikqs(ik)
  ENDIF

  ikk=ikmk
  ikq=ikmkmq

  ibnd=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (ibnd>nbnd_occ(ikk)) RETURN
  k2 = ibnd + st_

  nbnd_eff=nbnd_occ(ikq)
  npwq=ngk(ikq)
  IF (ltetra .OR. lgauss) nbnd_eff=nbnd

  jbnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (jbnd>nbnd_eff) RETURN

  ikwf=ik1 + (isolv-1)*nk
  k1 = nbnd*(ikwf-1) + jbnd
!
!  make the scalar product between the Bloch wavefunctions and 
!  the functions to orthogonalize.
!  In the metallic case compute also the weights in the projector
!

  IF (lgauss) THEN
     CALL wgauss_dev ((ef-et(ibnd,ikk)) / degauss, ngauss, wg1)
     CALL w0gauss_dev((ef-et(ibnd,ikk)) / degauss, ngauss, w0g) 
     w0g = w0g / degauss
     CALL wgauss_dev ( (ef - et (jbnd, ikq) ) / degauss, ngauss, wgp)
     deltae = et (jbnd, ikq) - et (ibnd, ikk)
     CALL wgauss_dev (deltae / degauss, 0, theta)
     wwg = wg1 * (1.d0 - theta) + wgp * theta
     IF (jbnd <= nbnd_occ (ikq) ) THEN
        IF (ABS (deltae) > 1.0d-5) THEN
           wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
        ELSE
           !
           !  if the two energies are too close takes the limit
           !  of the 0/0 ratio
           !
           wwg = wwg - alpha_pv * theta * w0g
        ENDIF
     ENDIF
  ELSE
     wwg=1.0_DP
  ENDIF


  aux=(0.0_DP,0.0_DP)
  DO ig = 1, npwq
     aux=aux + CONJG(evqk_d(ig,k1)) * dvpsik_d(ig,k2)
  ENDDO
  IF (npol==2) THEN
     DO ig = 1, npwq
        aux= aux + CONJG(evqk_d(ig+npwx,k1)) * dvpsik_d(ig+npwx,k2)
     ENDDO
  ENDIF
  ortho_ps_d(k1, k2) = wwg*aux

RETURN
END SUBROUTINE orthog_ps

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE orthog_last( st, nbndk, ikb, nk, npe, nsolv, &
                   dvpsik_d, sevqk_d, ortho_ps_d, npol, npwx, nbnd, nksbx_ph )
!-----------------------------------------------------------------------
  !
  !  Each thread computes one k point, one perturbation, and one band
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d
  USE wvfct_gpum,     ONLY : et => et_d
  USE many_k_mod,     ONLY : lgauss => lgauss_d, ltetra => ltetra_d,    &
                             degauss=> degauss_d, ngauss => ngauss_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d, ikmks=>ikmks_d,    &
                             ikmkmqs=>ikmkmqs_d, startkb_ph => startkb_ph_d,&
                             nbnd_occ => nbnd_occ_d, ef => ef_d
  !
  IMPLICIT NONE
#include<wgauss_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, npwx, nbnd, ikb, nksbx_ph
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of plane waves
  !! input: the number of bands
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: the number of bands to compute for each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  COMPLEX(DP), DEVICE, INTENT(IN) :: sevqk_d(npwx*npol, nbnd*nk*nsolv)
  !! inp/out: the bloch wavefunctions multiplied by S
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: dvpsik_d(npwx*npol, nbnd*nk*npe*nsolv)
  !! the functions to orthogonalize
  COMPLEX(DP), DEVICE, INTENT(IN) :: ortho_ps_d(nksbx_ph*nbnd*nsolv, &
                                                 nbnd*nksbx_ph*npe*nsolv)
  !
  !  Local variables
  !
  INTEGER :: isp, isolv, ipert, ik1, ik, st_, mb, ikmk, ikmkmq, ikk, ikq, &
             k2, k1, npwq, nbnd_eff, id, ibnd, jbnd, ig, ikwf

  REAL(DP) :: wg1
  !
  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe*nsolv) RETURN
  isp=(id-1)/nk+1
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1
  ik1=MOD(id-1,nk)+1
  st_=st(id)
  mb=nbndk(id)
  ik=ik1+startkb_ph(ikb)
  IF (isolv==2) THEN
     ikmk = ikmks(ik)
     ikmkmq = ikmkmqs(ik)
  ELSE
     ikmk=ikks(ik)
     ikmkmq=ikqs(ik)
  ENDIF

  ikk=ikmk
  ikq=ikmkmq

  ibnd=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (ibnd>nbnd_occ(ikk)) RETURN
  k2 = ibnd + st_

  nbnd_eff=nbnd_occ(ikq)
  npwq=ngk(ikq)
  IF (ltetra .OR. lgauss) nbnd_eff=nbnd

  ig=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ig>npwq) RETURN

  ikwf=ik1 + (isolv-1)*nk
!
!  dvpsik_d(ig,k2) has already the vector to orthogonalize.
!  Here we multiply by theta in the metallic case and ortho_ps contains
!  ortho_ps = beta * <evq|dvpsi>
!
  wg1=1.0_DP
  IF (lgauss) CALL wgauss_dev ((ef-et(ibnd,ikk)) / degauss, ngauss, wg1)
!  DO ig=1, npwq
     dvpsik_d(ig,k2)=-wg1*dvpsik_d(ig,k2)
!  ENDDO
  IF (npol==2) THEN
!     DO ig=1, npwq
        dvpsik_d(ig+npwx,k2)=-wg1*dvpsik_d(ig+npwx,k2)
!     ENDDO
  ENDIF

  DO jbnd=1, nbnd_eff
     k1 = nbnd*(ikwf-1) + jbnd
!     DO ig = 1, npwq
        dvpsik_d(ig,k2)=dvpsik_d(ig,k2) + sevqk_d(ig,k1) * ortho_ps_d(k1, k2)
!     ENDDO
     IF (npol==2) THEN
!        DO ig = 1, npwq
           dvpsik_d(ig+npwx,k2)=dvpsik_d(ig+npwx,k2) + sevqk_d(ig+npwx,k1) &
                                                     * ortho_ps_d(k1, k2)
!        ENDDO
     ENDIF
  ENDDO

RETURN
END SUBROUTINE orthog_last
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE wgauss_dev (x, n, wgauss)
  !-----------------------------------------------------------------------
  !! This function computes the approximate theta function for the
  !! given order n, at the point x:
  !
  !! * \( n \geq 0 \): Methfessel-Paxton case. See PRB 40, 3616 (1989).
  !! * \( n=-1 \): cold smearing (Marzari-Vanderbilt-DeVita-Payne,
  !!   see PRL 82, 3296 (1999)):
  !!   $$ \frac{1}{2} \text{erf}\(x-\frac{1}{\sqrt(2)}\) + \frac{1}{\sqrt{2\pi}} \exp
  !!   {-\(x-\frac{1}{sqrt{2}}\)^2} + 1/2 $$
  !! * \( n=-99 \): Fermi-Dirac case:
  !!   $$ \frac{1.0}{1.0+\exp{-x}} $$
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  implicit none
  real(DP), DEVICE :: wgauss
  !! output: the value of the function
  real(DP), VALUE :: x
  !! input: the argument of the function
  integer, VALUE :: n
  !! input: the order of the function
  !
  ! ... local variables
  !
  real(DP) :: a, hp, arg, hd, xp
  ! the coefficient a_n
  ! the hermitean function
  ! the argument of the exponential
  ! the hermitean function
  ! auxiliary variable (cold smearing)
  integer :: i, ni
  ! counter on the n indices
  ! counter on 2n
  real(DP), parameter :: maxarg = 200.d0
  ! maximum value for the argument of the exponential
  
  ! Fermi-Dirac smearing
  if (n.eq. - 99) then
     if (x.lt. - maxarg) then
        wgauss = 0.d0
     elseif (x.gt.maxarg) then
        wgauss = 1.d0
     else
        wgauss = 1.0d0 / (1.0d0 + exp ( - x) )
     endif
     return

  endif
  ! Cold smearing
  if (n.eq. - 1) then
     xp = x - 1.0d0 / sqrt (2.0d0)
     arg = min (maxarg, xp**2)
     wgauss = 0.5d0 * erf(xp) + 1.0d0 / sqrt (2.0d0 * pi) * exp ( - &
          arg) + 0.5d0
     return

  endif
  ! Methfessel-Paxton and plain gaussian cases 
  arg = -x 
  IF (arg .LT. sqrt(maxarg)) THEN 
     wgauss = 0.5_DP * ERFC( arg)
  ELSE 
     wgauss = 0._DP
  END IF 
  if (n.eq.0) return
  hd = 0.d0
  arg = min (maxarg, x**2)
  hp = exp ( - arg)
  ni = 0
  a = 1.d0 / sqrt (pi)
  do i = 1, n
     hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
     ni = ni + 1
     a = - a / (DBLE (i) * 4.0d0)
     wgauss = wgauss - a * hd
     hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
     ni = ni + 1
  enddo
  return
end subroutine wgauss_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE w0gauss_dev (x, n, w0gauss)
  !-----------------------------------------------------------------------
  !! The derivative of wgauss,  an approximation to the delta function:
  !
  !! * (n>=0): derivative of the corresponding Methfessel-Paxton \(\text{wgauss}\)
  !! * (n=-1 ): derivative of cold smearing:
  !!            $$ \frac{1}{\sqrt{\pi}}\text{exp}(-(x-\frac{1}{\sqrt{2}})^2)(2-\sqrt{2}x) $$
  !! * (n=-99): derivative of Fermi-Dirac function: \(0.5/(1.0+\text{cosh}(x))\)
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : sqrtpm1
  implicit none
  real(DP), DEVICE :: w0gauss
  !! output: the value of the function
  real(DP), VALUE :: x
  !! input: the point where to compute the function
  integer, VALUE :: n
  !! input: the order of the smearing function
  !
  ! ... local variables
  !
  real(DP) :: a, arg, hp, hd
  ! the coefficients a_n
  ! the argument of the exponential
  ! the hermite function
  ! the hermite function

  integer :: i, ni
  ! counter on n values
  ! counter on 2n values

  ! Fermi-Dirac smearing

  if (n.eq. - 99) then
     if (abs (x) .le.36.0) then
        w0gauss = 1.0d0 / (2.0d0 + exp ( - x) + exp ( + x) )
        ! in order to avoid problems for large values of x in the e
     else
        w0gauss = 0.d0
     endif
     return

  endif
  ! cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
  if (n.eq. - 1) then
     arg = min (200.d0, (x - 1.0d0 / sqrt (2.0d0) ) **2)
     w0gauss = sqrtpm1 * exp ( - arg) * (2.0d0 - sqrt ( 2.0d0) * x)
     return

  endif

  if (n.gt.10 .or. n.lt.0) PRINT *, 'wgauss_dev n wrong' 

  ! Methfessel-Paxton
  arg = min (200.d0, x**2)
  w0gauss = exp ( - arg) * sqrtpm1
  if (n.eq.0) return
  hd = 0.0d0
  hp = exp ( - arg)
  ni = 0
  a = sqrtpm1
  do i = 1, n
     hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
     ni = ni + 1
     a = - a / (DBLE (i) * 4.0d0)
     hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
     ni = ni + 1
     w0gauss = w0gauss + a * hp
  enddo
  return
END SUBROUTINE w0gauss_dev

#endif
