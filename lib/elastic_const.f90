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

SUBROUTINE print_elastic_constants(frozen_ions)
!
!  This routine writes on output the elastic constants
!
IMPLICIT NONE
LOGICAL, INTENT(IN) :: frozen_ions
INTEGER :: i, j

IF (frozen_ions) THEN
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
!  cubic case
!
!  c_11 = c_22 = c_33
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,igeo)
      y(igeo)=sigma_geo(3,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(1,1) = -alpha(2)
   el_con(2,2) = el_con(1,1)
   el_con(3,3) = el_con(1,1)
!
! c_12 = c_13 = c_23
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,igeo)
      y(igeo)=sigma_geo(1,1,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(1,2) = -alpha(2)
   el_con(2,1) = el_con(1,2)
   el_con(1,3) = el_con(1,2)
   el_con(3,1) = el_con(1,3)
   el_con(2,3) = el_con(1,2)
   el_con(3,2) = el_con(2,3)
!
! c_44 = c_55 = c_66
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,ngeo+igeo)
      y(igeo)=sigma_geo(2,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(4,4)=-alpha(2) * 0.5_DP
   el_con(5,5)=el_con(4,4)
   el_con(6,6)=el_con(4,4)
!
ELSEIF (ibrav==4) THEN
!
!  hexagonal case
!
!  c_11 = c_22
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,1,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(1,1) = -alpha(2)
   el_con(2,2) = el_con(1,1)
!
!  c_12
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(2,2,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(1,2) = -alpha(2)
   el_con(2,1) = el_con(1,2)
!
!  c_13
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(3,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(1,3) = -alpha(2)
   el_con(3,1) = el_con(1,3)
   el_con(2,3) = el_con(1,3)
   el_con(3,2) = el_con(2,3)
!
!  c_33
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,ngeo+igeo)
      y(igeo)=sigma_geo(3,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(3,3)=-alpha(2) 
!
!  c_44
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(2,3,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 ) 
   el_con(4,4)=-alpha(2) * 0.5_DP
   el_con(5,5)=el_con(4,4) 
   el_con(6,6)=0.5_DP*(el_con(1,1)-el_con(1,2))
ELSEIF (ibrav==5) THEN
!
!  trigonal case. We consider only the lowest symmetry class and calculate
!  c_15 in all cases
!
!
!  c_11 = c_22
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,1,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,1) = -alpha(2)
   el_con(2,2) = el_con(1,1)
!
!  c_12 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(2,2,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,2) = -alpha(2)
   el_con(2,1) = el_con(1,2)
!
!  c_13 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(3,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,3) = -alpha(2)
   el_con(3,1) = el_con(1,3)
!
!  c_14 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(2,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,4) = -alpha(2)
   el_con(4,1) = el_con(1,4)
   el_con(2,4) = -el_con(1,4)
   el_con(4,2) = el_con(2,4)
   el_con(5,6) = el_con(1,4)
   el_con(6,5) = el_con(5,6)
!
!  c_15 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,5) = -alpha(2)
   el_con(5,1) = el_con(1,5)
   el_con(2,5) = -el_con(1,5)
   el_con(5,2) = el_con(2,5)
   el_con(4,6) = el_con(2,5)
   el_con(6,4) = el_con(4,6)
!
!  c_33 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,ngeo+igeo)
      y(igeo)=sigma_geo(3,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,3) = -alpha(2)
!
!  c_44 = c_55
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(2,3,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(4,4) = -alpha(2) * 0.5_DP
   el_con(5,5) = el_con(4,4)

   el_con(6,6) = 0.5_DP * ( el_con(1,1) - el_con(1,2) )

ELSEIF (ibrav==6 .OR. ibrav==7) THEN
!
!  tetragonal case. Do not distinguish the two different classes. Computes
!  all elements, in some cases c_16=-c_26 will be zero, but we compute it
!  in any case.
!
!  c_11 = c_22
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,1,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,1) = -alpha(2)
   el_con(2,2) = el_con(1,1)
!
!  c_12 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(2,2,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,2) = -alpha(2)
   el_con(2,1) = el_con(1,2)
!
! c_13
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(3,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,1) = -alpha(2)
   el_con(1,3) = el_con(3,1)
   el_con(2,3) = el_con(3,1)
   el_con(3,2) = el_con(3,1)
!
! c_16, c_26
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,2,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,6) = -alpha(2)
   el_con(6,1) = el_con(1,6)
   el_con(2,6) = - el_con(1,6)
   el_con(6,2) = el_con(2,6)
!
! c_33
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,ngeo+igeo)
      y(igeo)=sigma_geo(3,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,3) = -alpha(2)
!
! c_44 = c_55
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(2,3,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(4,4) = -alpha(2) * 0.5_DP
   el_con(5,5) = el_con(4,4)
!
! c_66 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,2,3*ngeo+igeo)
      y(igeo)=sigma_geo(1,2,3*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(6,6) = -alpha(2) * 0.5_DP

ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
!
!  orthorombic case
!
!  c_11 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,1,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,1) = -alpha(2)
!
!  c_12
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(2,2,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,2) = -alpha(2)
   el_con(2,1) = el_con(1,2)
!
! c_13
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(3,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,1) = -alpha(2)
   el_con(1,3) = el_con(3,1)
!
! c_22
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,2,ngeo+igeo)
      y(igeo)=sigma_geo(2,2,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(2,2) = -alpha(2)
!
! c_23
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,2,ngeo+igeo)
      y(igeo)=sigma_geo(3,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(2,3) = -alpha(2)
   el_con(3,2) = el_con(2,3)
!
! c_33
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(3,3,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,3) = -alpha(2)
!
! c_44
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,3*ngeo+igeo)
      y(igeo)=sigma_geo(2,3,3*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(4,4) = -alpha(2) * 0.5_DP
!
! c_55
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,3,4*ngeo+igeo)
      y(igeo)=sigma_geo(1,3,4*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(5,5) = -alpha(2) * 0.5_DP
!
! c_66
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,2,5*ngeo+igeo)
      y(igeo)=sigma_geo(1,2,5*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(6,6) = -alpha(2) * 0.5_DP
ELSE
!
!  generic implementation, quite slow but should work with any lattice.
!  Computes all the elements of the elastic constants matrix, requires
!  6 * ngeo_strain self consistent calculations
!
!  c_11 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,1,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,1) = -alpha(2)
!
!  c_12 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(2,2,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,2) = -alpha(2)
   el_con(2,1) = el_con(1,2)
!
!  c_13 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(3,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,3) = -alpha(2)
   el_con(3,1) = el_con(1,3)
!
!  c_14 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(2,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,4) = -alpha(2)
   el_con(4,1) = el_con(1,4)
!
!  c_15 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,3,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,5) = -alpha(2)
   el_con(5,1) = el_con(1,5)
!
!  c_16 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,1,igeo)
      y(igeo)=sigma_geo(1,2,igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(1,6) = -alpha(2)
   el_con(6,1) = el_con(1,6)
!
!  c_22 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,2,ngeo+igeo)
      y(igeo)=sigma_geo(2,2,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(2,2) = -alpha(2)
!
!  c_23 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,2,ngeo+igeo)
      y(igeo)=sigma_geo(3,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(2,3) = -alpha(2)
   el_con(3,2) = el_con(2,3)
!  
!  c_24 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,2,ngeo+igeo)
      y(igeo)=sigma_geo(2,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(2,4) = -alpha(2)
   el_con(4,2) = el_con(2,4)
!  
!  c_25 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,2,ngeo+igeo)
      y(igeo)=sigma_geo(1,3,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(2,5) = -alpha(2)
   el_con(5,2) = el_con(2,5)
!  
!  c_26 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,2,ngeo+igeo)
      y(igeo)=sigma_geo(1,2,ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(2,6) = -alpha(2)
   el_con(6,2) = el_con(2,6)
!  
!  c_33 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(3,3,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,3) = -alpha(2)
!  
!  c_34 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(2,3,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,4) = -alpha(2)
   el_con(4,3) = el_con(3,4)
!  
!  c_35 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(1,3,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,5) = -alpha(2)
   el_con(5,3) = el_con(3,5)
!  
!  c_36 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(3,3,2*ngeo+igeo)
      y(igeo)=sigma_geo(1,2,2*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(3,6) = -alpha(2)
   el_con(6,3) = el_con(3,6)
!  
!  c_44
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,3*ngeo+igeo)
      y(igeo)=sigma_geo(2,3,3*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(4,4) = -alpha(2) * 0.5_DP
!  
!  c_45 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,3*ngeo+igeo)
      y(igeo)=sigma_geo(1,3,3*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(4,5) = -alpha(2) * 0.5_DP
   el_con(5,4) = el_con(4,5)
!  
!  c_46 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(2,3,3*ngeo+igeo)
      y(igeo)=sigma_geo(1,2,3*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(4,6) = -alpha(2) * 0.5_DP
   el_con(6,4) = el_con(4,6)
!  
!  c_55 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,3,4*ngeo+igeo)
      y(igeo)=sigma_geo(1,3,4*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(5,5) = -alpha(2) * 0.5_DP
!  
!  c_56 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,3,4*ngeo+igeo)
      y(igeo)=sigma_geo(1,2,4*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(5,6) = -alpha(2) * 0.5_DP
   el_con(6,5) = el_con(5,6)
!  
!  c_66 
!
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(1,2,5*ngeo+igeo)
      y(igeo)=sigma_geo(1,2,5*ngeo+igeo)
   ENDDO
   CALL polifit( x, y, ngeo, alpha, m1 )
   el_con(6,6) = -alpha(2) * 0.5_DP
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
