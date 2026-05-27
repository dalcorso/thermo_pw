!
! Copyright (C) 2026 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE spherical_mod
!---------------------------------------------------------------------------
!
!   This module contains the support routines for dealing with
!   spherical coordinates. Presently it contains only two routines:
!   cartesian_to_spherical passes from cartesian to spherical coordinates 
!   spherical_to_cartesian passes from spherical to cartesian coordinates
!   angles are given in radiants
!

  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC cartesian_to_spherical, spherical_to_cartesian

CONTAINS
!
!---------------------------------------------------------------------------
SUBROUTINE cartesian_to_spherical(x, x_spher)
!---------------------------------------------------------------------------
!
USE constants, ONLY : eps12
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: x(3)
REAL(DP), INTENT(OUT) :: x_spher(3)

REAL(DP) :: rmod, theta, phi, pi

pi=4.0_DP * ATAN(1.0_DP)
rmod=SQRT(x(1)**2+x(2)**2+x(3)**2)
IF (rmod>1.D-8*eps12) THEN
   theta=ACOS(x(3)/rmod)
   phi=ATAN2(x(2),x(1))
ELSE
!
!   null vector
!
   theta=0.0_DP
   phi=0.0_DP
ENDIF
IF (phi<0.0_DP) phi=phi+2.0_DP*pi
x_spher(1) = rmod
x_spher(2) = theta
x_spher(3) = phi

RETURN
END SUBROUTINE cartesian_to_spherical

!---------------------------------------------------------------------------
SUBROUTINE spherical_to_cartesian(x, x_spher)
!---------------------------------------------------------------------------
!
USE constants, ONLY : eps12
!
IMPLICIT NONE
REAL(DP), INTENT(OUT) :: x(3)
REAL(DP), INTENT(IN) :: x_spher(3)

REAL(DP) :: rmod, theta, phi, pi, cost, sint, cosp, sinp

rmod  = x_spher(1)
theta = x_spher(2)
phi   = x_spher(3)

cost=COS(theta)
sint=SIN(theta)
cosp=COS(phi)
sinp=SIN(phi)
x(1)=rmod*sint*cosp
x(2)=rmod*sint*sinp
x(3)=rmod*cost

RETURN
END SUBROUTINE spherical_to_cartesian

END MODULE spherical_mod
