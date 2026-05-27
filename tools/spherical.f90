!
! Copyright (C) 2026 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
PROGRAM spherical_coordinates
!--------------------------------------------------------------------
!
!  This program receives in input a vector in 3d in cartesian
!  coordinates and gives as output it polar coordinates.
!  Angles are given in radiant and in degrees.
!  0 <= theta < 180, 0<= phi < 360
!  It can also do the reverse. Given the polar coordinates of a vector
!  it gives its cartesian coordinates.
!
USE kinds, ONLY : DP
USE spherical_mod, ONLY : cartesian_to_spherical, spherical_to_cartesian
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end

USE io_global,        ONLY : stdin, stdout, meta_ionode

IMPLICIT NONE
CHARACTER(LEN=9) :: code='SPHERICAL'

REAL(DP) :: x(3), x_spher(3), rmod, theta, phi, pi
INTEGER  :: iwhat, ideg

CALL mp_startup ( start_images=.TRUE. )
CALL environment_start ( code )

pi=4.0_DP * ATAN(1.0_DP)
IF (meta_ionode) THEN
   WRITE(stdout,'(/,5x,"What do you want to do?")')
   WRITE(stdout,'(5x,"1) From cartesian to spherical")')
   WRITE(stdout,'(5x,"2) From spherical to cartesian")')
   READ(stdin,*) iwhat
   IF (iwhat==1) THEN
      WRITE(stdout,'(/,5x,"Give the cartesian coordinates x,y,z")')
      READ(stdin,*) x(1), x(2), x(3)
      WRITE(stdout,'(3f20.12)') x(1), x(2), x(3)
      CALL cartesian_to_spherical(x, x_spher) 
      WRITE(stdout,'(5x,"The spherical coordinates r, theta, phi")')
      WRITE(stdout,'(3f20.12)') x_spher(1), x_spher(2), x_spher(3)
      WRITE(stdout,'(3f20.12)') x_spher(1), x_spher(2)*180.0_DP/pi, &
                                x_spher(3)*180.0_DP/pi
   ELSEIF (iwhat==2) THEN
      WRITE(stdout,'(/,5x,"(1) Radiant or (2) degree?")')
      READ(stdin,*) ideg
      WRITE(stdout,'(/,5x,"Give the spherical coordinates r, theta, phi")')
      READ(stdin,*) rmod, theta, phi
      WRITE(stdout,'(3f20.12)') rmod,theta,phi
      IF (ideg==2) THEN
         theta=theta/180._DP * pi
         phi=phi/180._DP * pi
      ENDIF
      x_spher(1)=rmod
      x_spher(2)=theta
      x_spher(3)=phi
      CALL spherical_to_cartesian(x, x_spher)
      WRITE(stdout,'(5x,"The cartesian coordinates are x,y,z")')
      WRITE(stdout,'(3f20.12)') x(1), x(2), x(3)
   ENDIF
ENDIF

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM spherical_coordinates
