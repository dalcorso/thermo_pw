!
! Copyright (C) 2014-2015 Andrea Dal Corso 
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
!   TODO : energy for the tetragonal of class C_4h
!          energy for the monoclinic case
!          rombohedral cell, both energy and stress for S_6 and D_3d
!          Advanced: monoclinic and rombohedral both stress and energy C_2h
!


  USE kinds,     ONLY : DP
  USE constants, ONLY : ry_kbar
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL(DP) :: el_con(6,6)   ! The elastic constant
  REAL(DP) :: el_compliances(6,6)   ! The elastic compliances
 
  REAL(DP), ALLOCATABLE :: sigma_geo(:,:,:) ! The stress tensor computed 
                                            ! for each strain
  REAL(DP), ALLOCATABLE :: epsilon_geo(:,:,:) ! The strain tensor for each
                                            ! geometry
  REAL(DP), ALLOCATABLE :: epsilon_voigt(:,:) ! the strain tensor as a 6D array
                                            ! for each geometry
  REAL(DP) :: press=0.0_DP                  ! estimated pressure


  PUBLIC sigma_geo, epsilon_geo,  epsilon_voigt, &     ! public variables
         el_con, el_compliances, press,   &            ! public variables
         compute_elastic_constants,       &            !
         compute_elastic_constants_ene,   &            ! computing routines
         compute_elastic_compliances,     &            !
         print_elastic_constants,         &            ! public printing routines
         print_elastic_compliances,       &            !
         write_elastic, read_elastic,     &            ! public I/O on file
         macro_elasticity, print_macro_elasticity, &   ! public auxiliary tools
         print_sound_velocities, &                     ! public auxiliary tools
         compute_sound, voigt_index, el_cons_voigt     ! public auxiliary tools

CONTAINS

SUBROUTINE print_elastic_constants(elc, frozen_ions)
!
!  This routine writes on output the elastic constants
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: elc(6,6)
LOGICAL, INTENT(IN) :: frozen_ions
INTEGER :: i, j

WRITE(stdout,'(/,20x,40("-"),/)')
IF (frozen_ions) THEN
   WRITE(stdout, '(/,5x,"Frozen ions")')
ELSE
   WRITE(stdout, *)
ENDIF

IF (press /= 0.0_DP) THEN
   WRITE(stdout,'(5x,"Estimated pressure (kbar) ",f14.5)') press * ry_kbar
   WRITE(stdout,'(5x,"Elastic constants defined from stress-strain",f14.5)')
END IF

WRITE(stdout,'(5x,"Elastic constants C_ij (kbar) ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)

DO i=1,6
   WRITE(stdout,'(i5, 6f12.5)') i, (elc(i,j), j=1,6)
ENDDO

WRITE(stdout,'(/,5x,"1 bar = 10^5 Pa; 10 kbar = 1 GPa; 1 atm = 1.01325 bar;&
             & 1 Pa = 1 N/m^2")')
WRITE(stdout,'(5x,"1 Pa = 10 dyn/cm^2; 1 Mbar = 10^11 Pa")')
WRITE(stdout,'(5x,"1 torr = 1 mm Hg = 1/760 bar = 7.5006 x 10^-3 Pa",/)')

RETURN
END SUBROUTINE print_elastic_constants

SUBROUTINE print_elastic_compliances(els, frozen_ions)
!
!  This routine writes on output the elastic compliances
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: els(6,6)
LOGICAL, INTENT(IN) :: frozen_ions
INTEGER :: i, j

WRITE(stdout,'(/,20x,40("-"),/)')
IF (frozen_ions) THEN
   WRITE(stdout, '(/,5x,"Frozen ions")')
ELSE
   WRITE(stdout, *)
ENDIF
WRITE(stdout,'(5x,"Elastic compliances  S_ij (1/Mbar) ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)

DO i=1,6
   WRITE(stdout,'(i5, 6f12.5)') i, (els(i,j)*1.E3_DP, j=1,6)
ENDDO
WRITE(stdout,'(/,5x,"1/Mbar = 1/10^{11} Pa; 1 Pa = 1 N/m^2")')

RETURN
END SUBROUTINE print_elastic_compliances 

SUBROUTINE write_elastic(filename)
!
!  This routine writes the elastic constants and compliances on file.
!  It must be called after computing the elastic constant
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER :: find_free_unit
INTEGER :: outunit, ios, i, j

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_elastic','ploblem opening output file', ABS(ios))

IF (ionode) THEN
   DO i=1,6
      WRITE(outunit,'(4e20.10)') (el_con(i,j), j=1,6)
   ENDDO
   WRITE(outunit,*)
   DO i=1,6
      WRITE(outunit,'(4e20.10)') (el_compliances(i,j), j=1,6)
   END DO
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_elastic

SUBROUTINE read_elastic(filename, exists)
!
!  This routine reads the elastic constants and compliances from file.
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
LOGICAL, INTENT(OUT) :: exists
INTEGER :: inunit, ios, i, j
INTEGER :: find_free_unit

IF (ionode) THEN
   inunit=find_free_unit()
   OPEN(UNIT=inunit, FILE=TRIM(filename), STATUS='old', FORM='formatted', &
       ERR=100, IOSTAT=ios)
ENDIF

IF (ionode) THEN
   DO i=1,6
      READ(inunit,'(4e20.10)',ERR=100,IOSTAT=ios) (el_con(i,j), j=1,6)
   ENDDO
   READ(inunit,*)
   DO i=1,6
      READ(inunit,'(4e20.10)',ERR=100,IOSTAT=ios) (el_compliances(i,j), j=1,6)
   END DO
   CLOSE(inunit)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
IF (ios /= 0) THEN
   exists=.FALSE.
   RETURN
ENDIF
CALL mp_bcast(el_con,ionode_id,intra_image_comm)
CALL mp_bcast(el_compliances,ionode_id,intra_image_comm)
exists=.TRUE.

RETURN
END SUBROUTINE read_elastic

SUBROUTINE compute_elastic_constants(sigma_geo, epsil_geo, nwork, ngeo_strain, &
                                     ibrav, laue, m1)
!
!  This routine computes the elastic constants by fitting the stress-strain
!  relationship with a polynomial of order m1-1. This is calculated
!  on the basis of the Bravais lattice and Laue type.
!  With reference to the strain indicated in the strain.f90 module the
!  the routine expects to find in sigma_geo the stress that corresponds
!  to the following strain: 
!
!  Laue class         ibrav            strain
!  29,32 T_h, O_h      1               E  F 
!  29,32 T_h, O_h      2,3             E  F3 
!  19,23 C_6h, D_6h    4               C  E  H  
!  18,22 C_4h, D_4h    6,7             C  E  H  G
!  20    D_2h          8,9,10,11       C  D  E  G  H  I
!  25,27 D_3d, S_6     5               C  E  I
!  16    C_2h          12,-12, 13,-13  C  D  E  G  H  I
!  2     C_i           14              C  D  E  G  H  I
!
!  For each strain the routine expects ngeo_strain values of strain for a total
!  of nwork values.
!  Both sigma_geo and epsilon_geo are in cartesian coordinates in the
!  axis of the unstrained lattice.
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: sigma_geo(3,3,nwork), epsil_geo(3,3,nwork)
INTEGER,  INTENT(IN) :: nwork, ngeo_strain, ibrav, laue, m1

INTEGER :: npos

el_con=0.0_DP
SELECT CASE (laue)
   CASE (29,32)
!
!  cubic case (T_h and O_h Laue classes)
!
!  c_11 = c_22 = c_33
!
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,1) = el_con(3,3)
      el_con(2,2) = el_con(3,3)
!
! c_12 = c_13 = c_23
!
      CALL el_cons_ij(1, 3, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(1,3)
      el_con(2,1) = el_con(1,2)
      el_con(3,1) = el_con(1,3)
      el_con(2,3) = el_con(1,2)
      el_con(3,2) = el_con(2,3)
!
! c_44 = c_55 = c_66
!
      npos=ngeo_strain+1
      CALL el_cons_ij(4, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(5,5)=el_con(4,4)
      el_con(6,6)=el_con(4,4)

   CASE (19,23)
!
!  hexagonal case (C_6h and D_6h Laue classes)
!
!  c_11 = c_22
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(2,2) = el_con(1,1)
!
!  c_12
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
!  c_13
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,3) = el_con(3,1)
      el_con(2,3) = el_con(1,3)
      el_con(3,2) = el_con(2,3)
!
!  c_33
!
      npos=ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!  c_44
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(5, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(4,4)=el_con(5,5) 
      el_con(6,6)=0.5_DP*(el_con(1,1)-el_con(1,2))

   CASE (18,22)
!
!  tetragonal case, (C_4h and D_4h Laue classes)
!
!  c_11 = c_22
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(2,2) = el_con(1,1)
!
!  c_12 
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
! c_13
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,3) = el_con(3,1)
      el_con(2,3) = el_con(3,1)
      el_con(3,2) = el_con(2,3)
!
! c_33
!
      npos=ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_44 = c_55
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(5, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(4,4) = el_con(5,5)
!
! c_66 
!
      npos=3*ngeo_strain+1
      CALL el_cons_ij(6, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!   This part is non zero for the C_4h Laue class
!
      IF (laue==18) THEN
!
! c_16, c_26
!
         CALL el_cons_ij(6, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
         el_con(1,6) = el_con(6,1)
         el_con(2,6) = - el_con(1,6)
         el_con(6,2) = el_con(2,6)
     END IF

   CASE (20)
!
!  orthorhombic case, (D_2h Laue class)
!
!  c_11 
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_12
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
! c_13
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,3) = el_con(3,1)
!
! c_22
!
      npos=ngeo_strain+1
      CALL el_cons_ij(2, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_23
!
      CALL el_cons_ij(3, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,3) = el_con(3,2)
!
! c_33
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1, 1, npos), &
                                         sigma_geo(1, 1, npos), m1)
!
! c_66
!
      npos=3*ngeo_strain+1
      CALL el_cons_ij(6, 6, ngeo_strain, epsil_geo(1, 1, npos), &
                                         sigma_geo(1, 1, npos), m1)
!
! c_55
!
      npos=4*ngeo_strain+1
      CALL el_cons_ij(5, 5, ngeo_strain, epsil_geo(1, 1, npos), &
                                         sigma_geo(1, 1, npos), m1)
!
! c_44
!
      npos=5*ngeo_strain+1
      CALL el_cons_ij(4, 4, ngeo_strain, epsil_geo(1, 1, npos), &
                                         sigma_geo(1, 1, npos), m1)
   CASE (25,27)
!
!  trigonal case, (D_3d and S_6 Laue classes)
!
!
!  c_11 = c_22
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(2,2) = el_con(1,1)
!
!  c_12 
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
!  c_13 = c_23
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,3) = el_con(3,1)
      el_con(2,3) = el_con(1,3)
      el_con(3,2) = el_con(2,3)
!
!  c_14 
!
      CALL el_cons_ij(4, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,4) = el_con(4,1)
      el_con(2,4) = -el_con(1,4)
      el_con(4,2) = el_con(2,4)
      el_con(5,6) = el_con(1,4)
      el_con(6,5) = el_con(5,6)
!
!  c_33 
!
      npos=ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!  c_44 = c_55
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(4, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(5,5) = el_con(4,4)

      el_con(6,6) = 0.5_DP * ( el_con(1,1) - el_con(1,2) )
!
!   This part need to be computed only in the S_6 Laue class
!
!  c_15 
!
      IF (laue==27) THEN
         CALL el_cons_ij(5, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
         el_con(1,5) = el_con(5,1)
         el_con(2,5) = -el_con(1,5)
         el_con(5,2) = el_con(2,5)
         el_con(4,6) = el_con(2,5)
         el_con(6,4) = el_con(4,6)
      END IF

   CASE (16)
!
!   monoclinic case, (C_2h Laue class)
!
!  c_11 
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_12
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
! c_13
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,3) = el_con(3,1)
!
! c_22
!
      npos=ngeo_strain+1
      CALL el_cons_ij(2, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_23
!
      CALL el_cons_ij(3, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,3) = el_con(3,2)
!
! c_33
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_66
!
      npos=3*ngeo_strain+1
      CALL el_cons_ij(6, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_55
!
      npos=4*ngeo_strain+1
      CALL el_cons_ij(5, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_44
!
      npos=5*ngeo_strain+1
      CALL el_cons_ij(4, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!  c15
!
      IF (ibrav==-12.OR.ibrav==-13) THEN
!
!  monoclinic case unique axis b
!
!
!  c15
!
         CALL el_cons_ij(5, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
         el_con(1,5) = el_con(5,1)
!
!  c25
!
         npos=ngeo_strain+1
         CALL el_cons_ij(5, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(2,5) = el_con(5,2)
!
!  c35
!
         npos=2*ngeo_strain+1
         CALL el_cons_ij(5, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(3,5) = el_con(5,3)
!
!  c46
!
         npos=3*ngeo_strain+1
         CALL el_cons_ij(4, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(6,4) = el_con(4,6)

      ELSE
!
!  monoclinic case unique axis c
!
!  c16
!
         CALL el_cons_ij(6, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
         el_con(1,6) = el_con(6,1)
!
!  c26
!
         npos=ngeo_strain+1
         CALL el_cons_ij(6, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(2,6) = el_con(6,2)
!
!  c36
!
         npos=2*ngeo_strain+1
         CALL el_cons_ij(6, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)

         el_con(3,6) = el_con(6,3)
!
!  c45
!
         npos=5*ngeo_strain+1
         CALL el_cons_ij(5, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(4,5) = el_con(5,4)

      END IF

   CASE(2) 
!
!  generic implementation, quite slow but should work with any lattice.
!  Computes all the elements of the elastic constants matrix, requires
!  6 * ngeo_strain self consistent calculations
!
!  c_11 
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_12 
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
!  c_13 
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,3) = el_con(3,1)
!
!  c_14 
!
      CALL el_cons_ij(4, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,4) = el_con(4,1)
!
!  c_15 
!
      CALL el_cons_ij(5, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,5) = el_con(5,1)
!
!  c_16 
!
      CALL el_cons_ij(6, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,6) = el_con(6,1)
!
!  c_22 
!
      npos=ngeo_strain+1
      CALL el_cons_ij(2, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!  c_23 
!
      CALL el_cons_ij(3, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,3) = el_con(3,2)
!  
!  c_24 
!
      CALL el_cons_ij(4, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,4) = el_con(4,2)
!  
!  c_25 
!
      CALL el_cons_ij(5, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,5) = el_con(5,2)
!  
!  c_26 
!
      CALL el_cons_ij( 6, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                          sigma_geo(1,1,npos), m1)
      el_con(2,6) = el_con(6,2)
!  
!  c_33 
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_34 
!
      CALL el_cons_ij(4, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(3,4) = el_con(4,3)
!  
!  c_35 
!
      CALL el_cons_ij(5, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(3,5) = el_con(5,3)
!  
!  c_36 
!
      CALL el_cons_ij(6, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(3,6) = el_con(6,3)
!  
!  c_44
!
      npos=5*ngeo_strain+1
      CALL el_cons_ij(4, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_45 
!
      CALL el_cons_ij(5, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(4,5) = el_con(5,4)
!  
!  c_46 
!
      CALL el_cons_ij(6, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(4,6) = el_con(6,4)
!  
!  c_55 
!
      npos=4*ngeo_strain+1
      CALL el_cons_ij(5, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_56 
!
      CALL el_cons_ij(6, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(5,6) = el_con(6,5)
!  
!  c_66 
!
      npos=3*ngeo_strain+1
      CALL el_cons_ij(6, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
CASE DEFAULT
   CALL errore('compute_elastic_constants', 'Unknown Laue class ', 1)
END SELECT
el_con = el_con * ry_kbar

RETURN
END SUBROUTINE compute_elastic_constants

SUBROUTINE compute_elastic_constants_ene(energy_geo, epsil_geo, nwork, &
                                        ngeo_strain, ibrav, laue, omega, m1)
!
!  This routine computes the elastic constants by fitting the total
!  energy-strain relation with a polynomial of order m1-1. 
!  The components of the strain and the elastic constant tensors are
!  set on the basis of the Bravais lattice and Laue class.
!  The pressure is estimated and the calculated elastic constants are
!  corrected to keep into account for the difference between stress-strain
!  coefficients and energy coefficients. The output of this routine are
!  the elastic constants defined from the linear relationship between
!  stress and strain.
!
!  With reference to the strain indicated in the strain.f90 module the
!  the routine expects to find in energy_geo the total energy that corresponds
!  to the following strain: 
!
!  Laue class         ibrav            strain
!  29,32 T_h,  O_h     1               A   E   F 
!  29,32 T_h,  O_h     2,3             A   E   F3 
!  19,23 C_6h, D_6h    4               C   E   B1  A  H  
!  18,22 C_4h, D_4h    6,7             E   C   B   B1 G  H
!  20    D_2h          8,9,10,11       C   D   E   B  B1 B2  G  H  I
!  25,27 D_3d, S_6     5               Not available yet
!  16    C_2h          12,-12, 13,-13  Not available yet
!  2     C_i           14              Not available yet
!
!  For each strain the routine expects the energy for ngeo_strain values of 
!  strain for a total of nwork values.
!  epsilon_geo is in cartesian coordinates in the axis of the unstrained lattice.
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: epsil_geo(3,3,nwork), omega
REAL(DP), INTENT(IN) :: energy_geo(nwork)
INTEGER,  INTENT(IN) :: nwork, ngeo_strain, ibrav, laue, m1
REAL(DP) :: alpha(m1)
REAL(DP) :: b0, a0, aux
INTEGER  :: base_data
CHARACTER(LEN=41), PARAMETER :: FRMT='(/,5x,"Estimated pressure",f12.5," kbar")'

el_con=0.0_DP

WRITE(stdout,'(/,5x,"Fitting ",i3," functions,",i4, &
                  &" data per function, number of data =",i5)') &
                                      nwork/ngeo_strain, ngeo_strain, nwork
SELECT CASE (laue)
   CASE (29,32)
!
!  cubic case, (T_h and O_h Laue classes)
!
!
!   alpha(3) is multiplied by 2 because of the definition of the quadratic
!   interpolating polynomial. The second derivative of the polynomial with
!   respect to x is 2.0 * alpha(3)
!
      CALL el_cons_ij_ene(1, 1, 'C_11+2C_12', ngeo_strain, epsil_geo, &
                                                    energy_geo, alpha, m1)

      press=- alpha(2) / 3.0_DP / omega
      WRITE(stdout, FRMT) press*ry_kbar

      b0 = 1.0_DP / 3.0_DP / omega * ( 2.0_DP * alpha(3) ) +  &
           2.0_DP * press
      WRITE(stdout,'(5x,"Bulk Modulus",f18.5," kbar")') b0/3.0_DP*ry_kbar

      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_11', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
      el_con(2,2) = el_con(1,1)
      el_con(3,3) = el_con(1,1)
!
! c_12 = c_13 = c_23
!
      el_con(1,2) = 0.5_DP* (b0 - el_con(1,1)) 
      el_con(2,1) = el_con(1,2)
      el_con(3,1) = el_con(1,2)
      el_con(1,3) = el_con(1,2)
      el_con(2,3) = el_con(1,2)
      el_con(3,2) = el_con(1,2)
!
! c_44 = c_55 = c_66. Note that we should derive with respect to \epsilon_4,
! while the second derivative of the energy is with respect to \epsilon_23^2,
! hence the factor 1/4
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(2, 3, 'C_44', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(4,4) = 0.25_DP / 3.0_DP / omega * ( 2.0_DP * alpha(3) ) &
                        - 0.5_DP * press

      el_con(5,5)=el_con(4,4)
      el_con(6,6)=el_con(4,4)

   CASE (19,23)
!
!  hexagonal case, (C_6h and D_6h)
!
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
!
!  C_11=C_22
!
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
      el_con(2,2) = el_con(1,1)

      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
!
! C_33
!
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )

      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega / 2.0_DP
      WRITE(stdout,FRMT) press*ry_kbar
!
!  C_13=C_23
!
      aux = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
      el_con(1,3) = (aux - el_con(1,1) - el_con(3,3)) * 0.5_DP + press
      el_con(3,1) = el_con(1,3)
      el_con(2,3) = el_con(1,3)
      el_con(3,2) = el_con(2,3)

      base_data=3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega / 3.0_DP
      WRITE(stdout,FRMT) press*ry_kbar
      aux = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_12
!
!  pressure is added 3 times because in the expression we should
!  use C, while el_con(1,3) is B_13. To bring it to C_13 we must subtract P
!  so adding 2P we get C_12, to print in output B_12 we add another P.
!
      el_con(1,2) = (aux - 2.0_DP * el_con(1,1) - el_con(3,3) &
                         - 4.0_DP * el_con(1,3) ) * 0.5_DP + 3.0_DP * press
      el_con(2,1) = el_con(1,2)

!
!  C_44=C_55
!
      base_data=4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(5,5) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) - 0.5_DP * press
      el_con(4,4) = el_con(5,5)

      el_con(6,6) = (el_con(1,1) - el_con(1,2)) * 0.5_DP

   CASE(18,22)
!
!  tetragonal case, (C_4h and D_4h Laue classes)
!
!  C_33 
!
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_11 = C_22
!
      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
      el_con(2,2) = el_con(1,1)
!
!  C_12
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
      el_con(1,2) = 1.0_DP/omega*( 2.0_DP*alpha(3) )*0.5_DP-el_con(1,1) + &
                                                                     press 
      el_con(2,1) = el_con(1,2)
!
!  C_13=C_23
!
      base_data=3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
      el_con(1,3) = (1.0_DP/omega*(2.0_DP*alpha(3)) - el_con(3,3) &
                                              -el_con(1,1)) * 0.5_DP + press
      el_con(3,1) = el_con(1,3)
      el_con(2,3) = el_con(1,3)
      el_con(3,2) = el_con(2,3)
!
!  C_66
!
      base_data=4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_66', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(6,6) = 0.25_DP/omega*(2.0_DP*alpha(3)) - 0.5_DP * press
!
!  C_44=C_55
!
      base_data = 5 * ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(5,5) = 0.25_DP/omega*(2.0_DP*alpha(3)) - 0.5_DP * press
      el_con(4,4) = el_con(5,5)

   CASE(20)
!
!  Orthorombic case (D_2h Laue class)  
!
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
!
!  C_11
!
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )

      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_22', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
!
!  C_22
!
      el_con(2,2) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )

      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
!
!  C_33
!
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )

      base_data = 3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
      el_con(1,2) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(2,2) ) * 0.5_DP + press
      el_con(2,1)=el_con(1,2)

      base_data = 4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
!
!  C_13
!
      el_con(1,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(3,3) ) * 0.5_DP + press
      el_con(3,1)=el_con(1,3)

      base_data = 5*ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_23', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,FRMT) press*ry_kbar
!
!  C_23
!
      el_con(2,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                   - el_con(3,3) ) * 0.5_DP + press
      el_con(3,2)=el_con(2,3)

      base_data = 6*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_66', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
!
!  C_66
!
      el_con(6,6) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) - 0.5_DP * press

      base_data = 7*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
!
!  C_55
!
      el_con(5,5) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) - 0.5_DP * press

      base_data = 8*ngeo_strain+1
      CALL el_cons_ij_ene(2, 3, 'C_44', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
!
!  C_44
!
      el_con(4,4) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) - 0.5_DP * press

   CASE DEFAULT
      CALL errore('compute_elastic_constants_ene',&
                                   'Case not yet available',1)
END SELECT
el_con = el_con * ry_kbar

RETURN
END SUBROUTINE compute_elastic_constants_ene

SUBROUTINE compute_elastic_compliances(cmn,smn)
!
! This routine receives as input the elastic constants matrix and computes 
! its inverse, the elastic compliance matrix
!
USE kinds, ONLY : DP
USE matrix_inversion, ONLY : invmat
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: cmn(6,6)
REAL(DP), INTENT(INOUT) :: smn(6,6)

CALL invmat(6, cmn, smn)

RETURN
END SUBROUTINE compute_elastic_compliances

SUBROUTINE voigt_index(m, n, mn, flag)
!
!  If flag is .true., this routine receives two indeces 1<= m, n <=3 and
!  gives the voigt index 1<=mn<=6 corresponding to these two indices,
!  If flag is .false. it receive mn and gives as output m and n, m<=n
!
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: m, n, mn
LOGICAL, INTENT(IN) :: flag 
INTEGER :: voigt(3,3), mind(6), nind(6)
DATA voigt / 1, 6, 5, 6, 2, 4, 5, 4, 3 / 
DATA mind  / 1, 2, 3, 2, 1, 1 /
DATA nind  / 1, 2, 3, 3, 3, 2 /

IF (flag) THEN
   IF (m<1.OR.m>3.OR.n<1.OR.n>3) &
      CALL errore('voigt_index','m or n out or range',1)
   mn=voigt(m,n) 
ELSE
   IF (mn<1.OR.mn>6) &
      CALL errore('voigt_index','mn out of range',1)
   m=mind(mn)
   n=nind(mn)
ENDIF

RETURN
END SUBROUTINE voigt_index

SUBROUTINE el_cons_ij(pq, mn, ngeo, epsil_geo, sigma_geo, m1)
USE kinds, ONLY : DP
USE quadratic_surfaces, ONLY : polifit, write_poli

IMPLICIT NONE
INTEGER, INTENT(IN) :: mn, pq, ngeo
REAL(DP), INTENT(IN) :: epsil_geo(3,3,ngeo), sigma_geo(3,3,ngeo)
INTEGER :: igeo, m, n, p, q, mnin, pqin
INTEGER :: m1                  ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo)

WRITE(stdout,'(/,20x,40("-"),/)')
mnin=mn
CALL voigt_index(m,n,mnin,.FALSE.)
pqin=pq
CALL voigt_index(p,q,pqin,.FALSE.)
WRITE(stdout,'(/,5x,"Elastic constant ",2i5)') pq, mn
WRITE(stdout,'(/,10x,"strain",7x,"stress (kbar)")') 
DO igeo=1,ngeo
   x(igeo)=epsil_geo(m,n,igeo)
   y(igeo)=sigma_geo(p,q,igeo)
   WRITE(stdout,'(2f18.10)') x(igeo), y(igeo)*ry_kbar
ENDDO
CALL polifit( x, y, ngeo, alpha, m1 )
CALL write_poli(alpha,m1)
el_con(pq, mn) = -alpha(2)
!
!  The elastic constant tensor relates the stress to the strain in voigt
!  notation. Since e_23 = 0.5 e_4, e_13 = 0.5 e_5, e_12 = 0.5 e_6 we have
!  to divide by 2 the elements of the elastic constant calculated with 
!  off diagonal strain components
!
IF (m /= n) el_con(pq, mn) = el_con(pq, mn) * 0.5_DP

RETURN
END SUBROUTINE el_cons_ij

SUBROUTINE el_cons_ij_ene(m, n, label, ngeo, epsil_geo, energy_geo, alpha, m1)
USE kinds, ONLY : DP
USE quadratic_surfaces, ONLY : polifit, write_poli

IMPLICIT NONE
CHARACTER(LEN=*) :: label
INTEGER,  INTENT(IN)  :: m, n, ngeo, m1
REAL(DP), INTENT(IN)  :: epsil_geo(3,3,ngeo), energy_geo(ngeo)
REAL(DP), INTENT(OUT) :: alpha(m1)          ! the polynomial coefficients
INTEGER  :: igeo
REAL(DP) :: x(ngeo), y(ngeo)

WRITE(stdout,'(5x,a)')  TRIM(label)
WRITE(stdout,'(8x,"strain",7x,"Energy (Ry)")')
DO igeo=1,ngeo
   x(igeo)=epsil_geo(m,n,igeo)
   y(igeo)=energy_geo(igeo)
   WRITE(stdout,'(2f18.10)') x(igeo), y(igeo)
ENDDO
CALL polifit( x, y, ngeo, alpha, m1 )
CALL write_poli(alpha,m1)

RETURN
END SUBROUTINE el_cons_ij_ene

SUBROUTINE el_cons_voigt(elconv, elcon, flag)
!
!  This routine transform an elastic constant tensor in the 6x6 Voigt
!  form into a four index tensor 3x3x3x3 (flag=.false.) or viceversa 
!  (flag=.true.)
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: elcon(3,3,3,3)
REAL(DP), INTENT(INOUT) :: elconv(6,6)
LOGICAL, INTENT(IN) :: flag

INTEGER :: ij, mn, i, j, m, n

IF (flag) THEN
   elconv=0.0_DP
   DO ij=1,6
      CALL voigt_index(i,j,ij,.FALSE.)
      DO mn=1,6
         CALL voigt_index(m,n,mn,.FALSE.)
         elconv(ij,mn) = elcon(i,j,m,n) 
      ENDDO
   ENDDO
ELSE
   elcon=0.0_DP
   DO i=1,3
      DO j=1,3
         CALL voigt_index(i,j,ij,.TRUE.)
         DO m=1,3
            DO n=1,3
               CALL voigt_index(m,n,mn,.TRUE.)
               elcon(i,j,m,n) = elconv(ij,mn)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF

RETURN
END SUBROUTINE el_cons_voigt

SUBROUTINE macro_elasticity( ibrav, cmn, smn, b0v,  &
                             e0v, g0v, nuv, b0r, e0r, g0r, nur )
!
!  This routine collects some relationships that link the elastic constants 
!  to the parameters of the macroscopic elasticity and to
!  the polycristalline averages.
!  
!  It receives as input the elastic constants and compliances and gives
!  as output:
!
!  b0v : The Voigt average of the bulk modulus
!  e0v : The Voigt average Young modulus
!  g0v : The Voigt average Shear modulus
!  nuv : The Voigt average Poisson ratio 
!  b0r : The Reuss average of the bulk modulus
!  e0r : The Reuss average Young modulus
!  g0r : The Reuss average Shear modulus
!  nur : The Reuss average Poisson ratio 
!  
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
REAL(DP), INTENT(IN) :: cmn(6,6), smn(6,6)
REAL(DP), INTENT(OUT) :: b0v, e0v, g0v, nuv, b0r, e0r, g0r, nur
REAL(DP) :: c11v, c12v, c44v


c11v= (3.0_DP / 15.0_DP)*(cmn(1,1) + cmn(2,2) + cmn(3,3)) +   &
      (2.0_DP / 15.0_DP)*(cmn(1,2) + cmn(2,3) + cmn(1,3)) +   &
      (4.0_DP / 15.0_DP)*(cmn(4,4) + cmn(5,5) + cmn(6,6)) 
c12v= (1.0_DP / 15.0_DP)*(cmn(1,1) + cmn(2,2) + cmn(3,3)) +   &
      (4.0_DP / 15.0_DP)*(cmn(1,2) + cmn(2,3) + cmn(1,3)) -   &
      (2.0_DP / 15.0_DP)*(cmn(4,4) + cmn(5,5) + cmn(6,6)) 
c44v= (c11v-c12v)*0.5_DP
b0v=( c11v + 2.0_DP * c12v) / 3.0_DP
e0v = (c11v - c12v)*(c11v+2.0_DP*c12v) / (c11v+c12v)
g0v = c44v
nuv = e0v/ (2.0_DP * g0v) - 1.0_DP
b0r=1.0_DP/(smn(1,1) + smn(2,2) + smn(3,3) + &
    2.0_DP*smn(1,2) + 2.0_DP*smn(1,3)+ 2.0_DP*smn(2,3))
e0r = 15.0_DP /( 3.0_DP *(smn(1,1) + smn(2,2) + smn(3,3)) +  &
                 2.0_DP *(smn(1,2) + smn(2,3) + smn(1,3)) +  &
                         (smn(4,4) + smn(5,5) + smn(6,6)) )
g0r = 15.0_DP /( 4.0_DP *(smn(1,1) + smn(2,2) + smn(3,3)) - &
                 4.0_DP *(smn(1,2) + smn(2,3) + smn(1,3)) + &
                 3.0_DP *(smn(4,4) + smn(5,5) + smn(6,6)) )
nur = e0r/ (2.0_DP * g0r) - 1.0_DP

RETURN
END SUBROUTINE macro_elasticity

SUBROUTINE print_macro_elasticity(ibrav, cmn, smn, macro_el, flag)
!
! If flag is .true. print the macroscopic elastic properties. If it
! is false it only computes them
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: cmn(6,6), smn(6,6)
REAL(DP), INTENT(INOUT) :: macro_el(8)
INTEGER, INTENT(IN) :: ibrav
LOGICAL, INTENT(IN) :: flag
REAL(DP) :: b0v, e0v, g0v, nuv, b0r, e0r, g0r, nur

CALL macro_elasticity( ibrav, cmn, smn, b0v, &
                             e0v, g0v, nuv, b0r, e0r, g0r, nur )

IF (flag) THEN
   WRITE(stdout,'(/,20x,40("-"),/)')

   WRITE(stdout, '(/,5x, "Voigt approximation:")') 

   WRITE(stdout, '(5x, "Bulk modulus  B = ",f12.5," kbar")') b0v
   WRITE(stdout, '(5x, "Young modulus E = ",f12.5," kbar")') e0v
   WRITE(stdout, '(5x, "Shear modulus G = ",f12.5," kbar")') g0v
   WRITE(stdout, '(5x,"Poisson Ratio n = ",f12.5)') nuv

   WRITE(stdout, '(/,5x, "Reuss approximation:")') 
   WRITE(stdout, '(5x, "Bulk modulus  B = ",f12.5," kbar")') b0r
   WRITE(stdout, '(5x, "Young modulus E = ",f12.5," kbar")') e0r
   WRITE(stdout, '(5x, "Shear modulus G = ",f12.5," kbar")') g0r
   WRITE(stdout, '(5x,"Poisson Ratio n = ",f12.5)') nur

   WRITE(stdout, '(/,5x, "Voigt-Reuss-Hill average of the two &
                                                         &approximations:")') 
   WRITE(stdout, '(5x, "Bulk modulus  B = ",f12.5," kbar")') (b0v+b0r)*0.5_DP
   WRITE(stdout, '(5x, "Young modulus E = ",f12.5," kbar")') (e0v+e0r)*0.5_DP
   WRITE(stdout, '(5x, "Shear modulus G = ",f12.5," kbar")') (g0v+g0r)*0.5_DP
   WRITE(stdout, '(5x,"Poisson Ratio n = ",f12.5)') (e0v+e0r)/    &
                                                 (2.d0*(g0v+g0r))-1.0_DP
END IF

macro_el(1)=b0v
macro_el(2)=e0v
macro_el(3)=g0v
macro_el(4)=nuv
macro_el(5)=b0r
macro_el(6)=e0r
macro_el(7)=g0r
macro_el(8)=nur

RETURN
END SUBROUTINE print_macro_elasticity

SUBROUTINE print_sound_velocities(ibrav, cmn, smn, density, vp, vb, vg)
!
!  In input the elastic constants are in kbar, the elastic compliances 
!  in kbar^-1 and the density in Kg/m^3. The sound velocity is printed
!  in m/sec
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
REAL(DP), INTENT(IN) :: cmn(6,6), smn(6,6), density
REAL(DP), INTENT(OUT) :: vp, vb, vg 
REAL(DP) :: b0v, e0v, g0v, nuv, b0r, e0r, g0r, nur

REAL(DP) :: g0, b0

CALL macro_elasticity( ibrav, cmn, smn, b0v, e0v, g0v, nuv, b0r, e0r, g0r, nur )

WRITE(stdout, '(/,5x, "Voigt-Reuss-Hill average; sound velocities:",/)') 

g0 = ( g0r + g0v ) * 0.5_DP
b0 = ( b0r + b0v ) * 0.5_DP

IF (b0 + 4.0_DP * g0 / 3.0_DP > 0.0_DP) THEN
   vp = SQRT( ( b0 + 4.0_DP * g0 / 3.0_DP ) * 1.D8 / density )
   WRITE(stdout, '(5x, "Compressional V_P = ",f12.3," m/s")') vp
ELSE
   vp=0.0_DP
   WRITE(stdout, '(5x, "The system is unstable for compressional deformations")') 
ENDIF
IF (b0 > 0.0_DP) THEN
   vb = SQRT( b0 * 1.D8 / density )
   WRITE(stdout, '(5x, "Bulk          V_B = ",f12.3," m/s")') vb
ELSE
   vb =0.0_DP
   WRITE(stdout, '(5x, "The system is unstable")') 
END IF
IF (g0 > 0.0_DP) THEN
   vg = SQRT( g0 * 1.D8 / density )
   WRITE(stdout, '(5x, "Shear         V_G = ",f12.3," m/s")') vg
ELSE
   vg = 0.0_DP
   WRITE(stdout, '(5x, "The system is unstable for shear deformations")') 
ENDIF


RETURN
END SUBROUTINE print_sound_velocities

SUBROUTINE set_sound_mat(elcon, qvec, soundmat)
!
!  This routine receives the elastic constants in the format C_{ijkl}
!  a direction for the propagation of the sound waves and gives as
!  output the matrix that has to be diagonalized to obtain the sound
!  speed in the input direction
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: elcon(3,3,3,3)   ! the elastic constants
REAL(DP), INTENT(IN) :: qvec(3)           ! the direction of the sound waves
REAL(DP), INTENT(INOUT) :: soundmat(3,3)  ! the matrix that must be 
                                          ! diagonalized to compute the 
                                          ! sound speed
INTEGER :: i,j,k,l

DO i=1,3
   DO l=1,3
      soundmat(i,l)=0.0_DP
      DO j=1,3
         DO k=1,3
            soundmat(i,l) = soundmat(i,l) + elcon(i,j,k,l)*qvec(j)*qvec(k)
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE set_sound_mat

SUBROUTINE compute_sound(elcon, qvec, density, sound_speed, sound_disp)
!
!  This routine receives as input the elastic constants in Voigt notation 
!  C_{ij}, the direction of propagation of the sound, and the
!  density of the solid. It gives as output the sound velocity and the
!  eigenvectors that give the polarization of the sound wave.
!  
!  In input the density is in kg/m^3 and the elastic constants in kbar.
!  On output the speed of sound is in m/s.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(IN) :: elcon(3,3,3,3)
REAL(DP), INTENT(IN) :: qvec(3)
REAL(DP), INTENT(IN) :: density

REAL(DP), INTENT(INOUT) :: sound_speed(3)
REAL(DP), INTENT(INOUT) :: sound_disp(3,3)

REAL(DP) :: soundmat(3,3)
INTEGER :: i
!
!  set up the matrix to diagonalize
!
CALL set_sound_mat(elcon, qvec, soundmat)
!
!  and diagonalize it
!
!WRITE(stdout,*) soundmat(1,1), soundmat(1,2), soundmat(1,3)
!WRITE(stdout,*) soundmat(2,1), soundmat(2,2), soundmat(2,3)
!WRITE(stdout,*) soundmat(3,1), soundmat(3,2), soundmat(3,3)

CALL rdiagh( 3, soundmat, 3, sound_speed, sound_disp )
!
!  the factor 1D8 converts from kbar to N/m^2
!
DO i=1,3
   IF (sound_speed(i) > 0.0_DP) THEN
      sound_speed(i) = SQRT( sound_speed(i) * 1.D8 / density )
   ELSE
      sound_speed(i) = 0.0_DP
   ENDIF
ENDDO

RETURN
END SUBROUTINE compute_sound

END MODULE elastic_constants
