!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE elastic_constants
!
!   this module contains the support routines for the calculation
!   of the elastic constant
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL(DP) :: el_con(6,6)   ! The elastic constant
  REAL(DP), ALLOCATABLE :: internal_strain(:,:,:)  ! 6,3,nat internal strain
 
  REAL(DP), ALLOCATABLE :: sigma_geo(:,:,:) ! The stress tensor computed 
                                            ! for each strain
  REAL(DP), ALLOCATABLE :: epsilon_geo(:,:,:) ! The strain tensor for each
                                            ! geometry
  REAL(DP), ALLOCATABLE :: epsilon_voigt(:,:) ! the strain tensor as a 6D array
                                            ! for each geometry

  PUBLIC trans_epsilon, el_con, sigma_geo, epsilon_geo, apply_strain, &
         epsilon_voigt, compute_elastic_constants, print_elastic_constants

CONTAINS
!
SUBROUTINE trans_epsilon(eps_voigt, eps_tensor, flag)
!
!  This routine transforms the strain tensor from a one dimensional
!  vector with 6 components, to a 3 x 3 symmetric tensor and viceversa
!
!  flag
!   1      6   --> 3 x 3
!  -1    3 x 3 --> 6 
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP) :: eps_voigt(6), eps_tensor(3,3)
INTEGER :: flag
INTEGER :: i, j

IF ( flag == 1 ) THEN
   eps_tensor(1,1) = eps_voigt(1)
   eps_tensor(1,2) = eps_voigt(6) * 0.5_DP
   eps_tensor(1,3) = eps_voigt(5) * 0.5_DP
   eps_tensor(2,2) = eps_voigt(2) 
   eps_tensor(2,3) = eps_voigt(4) * 0.5_DP
   eps_tensor(3,3) = eps_voigt(3) 
   DO i = 1, 3
      DO j = 1, i-1
         eps_tensor(i,j) = eps_tensor(j,i)
      END DO
   END DO
ELSEIF ( flag == -1 ) THEN
  eps_voigt(1) = eps_tensor(1,1)
  eps_voigt(2) = eps_tensor(2,2)
  eps_voigt(3) = eps_tensor(3,3)
  eps_voigt(4) = eps_tensor(2,3) * 2.0_DP
  eps_voigt(5) = eps_tensor(1,3) * 2.0_DP
  eps_voigt(6) = eps_tensor(1,2) * 2.0_DP
ELSE
   CALL errore('trans_epsilon', 'unknown flag', 1)
ENDIF

RETURN
END SUBROUTINE trans_epsilon

SUBROUTINE print_elastic_constants(what)
!
!  This routine writes on output the elastic constants
!
IMPLICIT NONE
CHARACTER(LEN=*) :: what
INTEGER :: i, j

IF (what(1:3)=='fi_') THEN
   WRITE(stdout, '(/,5x,"Frozen ions")')
ELSE
   WRITE(stdout, *)
ENDIF
WRITE(stdout,'(5x,"Elastic constants C_ij (kbar) ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)

DO i=1,6
   WRITE(stdout,'(i5, 6f12.5)') i, (el_con(i,j), j=1,6)
ENDDO
WRITE(stdout,'(/)')

RETURN
END SUBROUTINE print_elastic_constants

SUBROUTINE compute_elastic_constants(sigma_geo, epsil_geo, nwork, ngeo, ibrav)
!
!  This routine computes the elastic constants by fitting the stress-strain
!  relation with a second order polynomial. This is calculated
!  on the basis of the Bravais lattice and Laue type.
!
USE constants, ONLY : ry_kbar
!USE symm_aux, ONLY : ilaue
IMPLICIT NONE
REAL(DP), INTENT(IN) :: sigma_geo(3,3,nwork), epsil_geo(3,3,nwork)
INTEGER, INTENT(IN) :: nwork, ngeo, ibrav
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo)
INTEGER :: i, j, igeo

el_con=0.0_DP
IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
!
!  These relationships are true for cubic solids
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,igeo)
      y(igeo)=sigma_geo(3,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(1,1) = -alpha(2)
   el_con(2,2) = el_con(1,1)
   el_con(3,3) = el_con(1,1)
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,igeo)
      y(igeo)=sigma_geo(1,1,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(1,2) = -alpha(2)
   el_con(1,3) = el_con(1,2)
   el_con(2,3) = el_con(1,2)
   DO i=1,3
      DO j=1,i-1
         el_con(i,j)=el_con(j,i)
      ENDDO
   ENDDO
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,ngeo+igeo)
      y(igeo)=sigma_geo(2,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(4,4)=-alpha(2) * 0.5_DP
   el_con(5,5)=el_con(4,4)
   el_con(6,6)=el_con(4,4)
ENDIF
el_con = el_con * ry_kbar

RETURN
END SUBROUTINE compute_elastic_constants

SUBROUTINE apply_strain(a, b, epsil)
!
!  This routine receives as input a vector a and a strain tensor \epsil and
!  gives as output a vector b = a + \epsil a 
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: a(3), epsil(3,3)
REAL(DP), INTENT(OUT) :: b(3)
INTEGER :: i

b = a 
DO i=1,3
   b(:)=b(:) + epsil(:,i) * a(i)
ENDDO

RETURN
END SUBROUTINE apply_strain

END MODULE elastic_constants
