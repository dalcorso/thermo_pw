!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! Copyright (C) 2023 Andrea Dal Corso for the device routine 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
#if defined(__CUDA)
ATTRIBUTES(DEVICE) SUBROUTINE ylmr2_dev (lmax2, g, gg, ylm)
  !-----------------------------------------------------------------------
  !
  !     Real spherical harmonics ylm(G) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical 
  !     Recipes but avoiding the calculation of factorials that generate 
  !     overflow for lmax > 11
  !     Last modified May 2nd, 2021, by PG
  !     Device version by Andrea Dal Corso (2023)
  !
  USE upf_kinds, ONLY : DP
  USE upf_const, ONLY : pi, fpi
  !
  IMPLICIT NONE
  !
  integer, VALUE, INTENT(IN) :: lmax2
  real(DP), VALUE, INTENT(IN) :: gg 
  real(DP), INTENT(IN) :: g(3)
  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  ! Incorrect results will ensue if gg != g(1)^2 + g(2)^2 +g(3)^2 
  !
  real(DP), DEVICE, INTENT(OUT) :: ylm (lmax2)
  !
  ! local variables
  !
  real(DP), parameter :: eps = 1.0d-9
  real(DP) :: cost , sent, phi 
  real(DP) :: c, gmod
  integer :: lmax, l, m, lm, lm1, lm2
  !
  if (lmax2 < 1) return
  do lmax = 0, 25
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  STOP 1
10 continue
  !
  if (lmax == 0) then
     ylm(1) =  SQRT (1.d0 / fpi)
     return
  end if
  !
  gmod = SQRT (gg)
  if (gmod < eps) then
     cost = 0.d0
  else
     cost = g(3)/gmod
  endif
  sent = SQRT(MAX(0.0_dp,1.0_dp-cost*cost))
  !
  !  cost = cos(theta), sent = sin(theta), with theta = polar angle
  !
  !  The Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m),
  !  where P(:,l,m) are the Legendre Polynomials (0 <= m <= l),
  !  and are currently stored into Ylm(:,lm) where lm = l**2+1+2*m
  !  (one might also store them into an array Q(l,m) for each ig)
  !
  ylm (1) = 1.d0
  ylm (2) = cost
  ylm (4) =-sent/sqrt(2.d0)
  do l = 2, lmax
     !
     !  recursion on l for Q(:,l,m)
     !
     do m = 0, l - 2
        lm = (l  )**2 + 1 + 2*m
        lm1= (l-1)**2 + 1 + 2*m
        lm2= (l-2)**2 + 1 + 2*m
        ylm(lm) = cost*(2*l-1)/sqrt(DBLE(l*l-m*m))*ylm(lm1) &
                - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*ylm(lm2)
     end do
     lm = (l  )**2 + 1 + 2*l
     lm1= (l  )**2 + 1 + 2*(l-1)
     lm2= (l-1)**2 + 1 + 2*(l-1)
     ylm(lm1) = cost * sqrt(DBLE(2*l-1)) * ylm(lm2)
     ylm(lm ) = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent*ylm(lm2) 
  enddo
  !
  ! now add cos(phi), sin(phi), and other factors to get the true Ylm
  ! beware the arc tan, it is defined modulo pi
  !
  if (g(1) > eps) then
     phi  = atan( g(2)/g(1) )
  else if (g(1) < -eps) then
     phi  = atan( g(2)/g(1) ) + pi
  else
     phi  = sign( pi/2.d0,g(2) )
  end if
  lm = 1
  ylm(1) = ylm(1) / sqrt(fpi)
  !
  do l = 1, lmax
     !
     c = sqrt (DBLE(2*l+1) / fpi)
     !
     ! Y_lm, m = 0
     !
     lm = lm + 1
     ylm(lm) = c * ylm(lm)
     !
     do m = 1, l
        !
        ! Y_lm, m > 0
        !
        lm = lm + 2
        ylm(lm-1) = c * sqrt(2.d0) * ylm(lm) * cos (m*phi)
        !
        ! Y_lm, m < 0
        !
        ylm(lm  ) = c * sqrt(2.d0) * ylm(lm) * sin (m*phi)
        !
     end do
  end do
  !
  return
end subroutine ylmr2_dev
#endif
