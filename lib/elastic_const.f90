!
! Copyright (C) 2014-2019 Andrea Dal Corso 
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
  REAL(DP) :: strs(3,3)=0.0_DP              ! estimated stress


  PUBLIC sigma_geo, epsilon_geo, epsilon_voigt, &      ! public variables
         el_con, el_compliances, press,   &            ! public variables
         compute_elastic_constants,       &            ! 
         elastic_constants_from_compliances, &         !
         compute_elastic_constants_ene,   &            ! computing routines
         compute_elastic_compliances,     &            !
         print_elastic_constants,         &            ! public printing 
                                                       ! routines
         print_elastic_compliances,       &            !
         print_el_cons_info,              &            ! print the number
                                                       ! of calculations 
         write_elastic, read_elastic,     &            ! public I/O on file
         macro_elasticity, print_macro_elasticity, &   ! public auxiliary tools
         print_sound_velocities, &                     ! public auxiliary tools
         compute_sound, &                              ! public auxiliary tools
         correct_for_stress,  &                        ! corrects the second
                                         ! energy derivatives to give the
                                         ! stress-strain elastic constants
         correct_for_pressure, &         ! corrects for pressure
         expand_el_cons,       &         ! expand the elastic constant read
                                         ! from file in a full tensor
         write_el_cons_on_file,     &    ! write elastic constants on file
         read_el_cons_from_file,    &    ! read elastic constants from file
         write_macro_el_on_file,    &    ! write macro-elasticity variables on file 
         write_sound_on_file             ! write sound velocities on file 


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

SUBROUTINE compute_elastic_constants(sigma_geo, epsil_geo, nwork, ngeo_strain,&
                                     ibrav, laue, m1)
!
!  This routine computes the elastic constants by fitting the stress-strain
!  relationship with a polynomial of degree m1-1. The nonvanishing
!  components of the tensor are calculated on the basis of the Bravais 
!  lattice and Laue type.
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
!  25,27 D_3d, S_6     4, 5            C  E  I
!  16    C_2h          12,-12, 13,-13  C  D  E  G  H  I
!  2     C_i           14              C  D  E  G  H  I
!
!  For each strain the routine expects ngeo_strain values of strain 
!  for a total of nwork values.
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
!  c_21 = c_12
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
!  c_31
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(3,2) = el_con(3,1)
!
!  c_33
!
      npos=ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!  c_13  could be different from c_31 if the cell has a stress
!
      CALL el_cons_ij(1, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,3) = el_con(1,3)
!
!  c_44 = c_55
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
!  c_21 = c_12
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
! c_31
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(3,2) = el_con(3,1)
!
! c_33
!
      npos=ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_13 could be different from c_31 if the cell has a stress
!
      CALL el_cons_ij(1, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,3) = el_con(1,3)
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
         CALL el_cons_ij(1, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(6,1) = el_con(1,6)
         el_con(6,2) = -el_con(1,6)
         el_con(2,6) = el_con(6,2)

     END IF

   CASE (20)
!
!  orthorhombic case, (D_2h Laue class)
!
! c_11 
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
! c_21
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
! c_31
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
! c_22
!
      npos=ngeo_strain+1
      CALL el_cons_ij(2, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_12 could be different from c_21 if the cell has a stress
!
      CALL el_cons_ij(1, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_32
!
      CALL el_cons_ij(3, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_33
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1, 1, npos), &
                                         sigma_geo(1, 1, npos), m1)
!
! c_13 could be different from c_31 if the cell has a stress
!
      CALL el_cons_ij(1, 3, ngeo_strain, epsil_geo(1, 1, npos), &
                                         sigma_geo(1, 1, npos), m1)
!
! c_23 could be different from c_32 if the cell has a stress
!
      CALL el_cons_ij(2, 3, ngeo_strain, epsil_geo(1, 1, npos), & 
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
!  c_21 = c_12 
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(1,2) = el_con(2,1)
!
!  c_31 = c_32
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
      el_con(3,2) = el_con(3,1)
!
!  c_41 = c_14 
!
      IF (ibrav/=4) THEN
         CALL el_cons_ij(4, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
         el_con(1,4) = el_con(4,1)
         el_con(2,4) = -el_con(1,4)
         el_con(4,2) = el_con(2,4)
         el_con(5,6) = el_con(1,4)
         el_con(6,5) = el_con(5,6)
      ENDIF
!
!  c_33 
!
      npos=ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)

!
!  c_13 could be different from c_31 if the cell has a stress
!
      CALL el_cons_ij(1, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
      el_con(2,3)=el_con(1,3)
!
!  c_44 = c_55
!
      npos=2*ngeo_strain+1
      IF (ibrav==4) THEN
         CALL el_cons_ij(5, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(4,4) = el_con(5,5)
      ELSE
         CALL el_cons_ij(4, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
         el_con(5,5) = el_con(4,4)
!
!  This part needs to be computed only in the S_6 Laue class for the
!  trigonal lattice
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
      ENDIF

      el_con(6,6) = 0.5_DP * ( el_con(1,1) - el_con(1,2) )

   CASE (16)
!
!   monoclinic case, (C_2h Laue class)
!
!  c_11 
!
      CALL el_cons_ij(1, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_21
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
! c_31
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
! c_22
!
      npos=ngeo_strain+1
      CALL el_cons_ij(2, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_12 could be different from c_21 if the cell has a stress
!
      CALL el_cons_ij(1, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_32
!
      CALL el_cons_ij(3, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_33
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_13 could be different from c_31 if the cell has a stress
!
      CALL el_cons_ij(1, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
! c_23 could be different from c_32 if the cell has a stress
!
      CALL el_cons_ij(2, 3, ngeo_strain, epsil_geo(1,1,npos), &
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

      IF (ibrav==12.OR.ibrav==13) THEN
!
!  monoclinic case unique axis c
!
!
!  c_61
!
         CALL el_cons_ij(6, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_62
!
         npos=ngeo_strain+1
         CALL el_cons_ij(6, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c63
!
         npos=2*ngeo_strain+1
         CALL el_cons_ij(6, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_54
!
         npos=5*ngeo_strain+1
         CALL el_cons_ij(5, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_16
!
         npos=3*ngeo_strain+1
         CALL el_cons_ij(1, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_26
!
         CALL el_cons_ij(2, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_36
!
         CALL el_cons_ij(3, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_45
!
         npos=4*ngeo_strain+1
         CALL el_cons_ij(4, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)

      ENDIF

      IF (ibrav==-12.OR.ibrav==-13) THEN
!
!  monoclinic case unique axis b
!
!
!   c_15
!
         npos=4*ngeo_strain+1
         CALL el_cons_ij(1, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!   c_25
!
         CALL el_cons_ij(2, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!   c_35
!
         CALL el_cons_ij(3, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_64
!
         npos=5*ngeo_strain+1
         CALL el_cons_ij(6, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_51
!
         CALL el_cons_ij(5, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_52
!
         npos=ngeo_strain+1
         CALL el_cons_ij(5, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_53
!
         npos=2*ngeo_strain+1
         CALL el_cons_ij(5, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
!
!  c_46
!
         npos=3*ngeo_strain+1
         CALL el_cons_ij(4, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                            sigma_geo(1,1,npos), m1)
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
!  c_21 
!
      CALL el_cons_ij(2, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_31 
!
      CALL el_cons_ij(3, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_41 
!
      CALL el_cons_ij(4, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_51 
!
      CALL el_cons_ij(5, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_61 
!
      CALL el_cons_ij(6, 1, ngeo_strain, epsil_geo, sigma_geo, m1)
!
!  c_22 
!
      npos=ngeo_strain+1
      CALL el_cons_ij(2, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!  c_12 
!
      CALL el_cons_ij(1, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!
!  c_32
!
      CALL el_cons_ij(3, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_42 
!
      CALL el_cons_ij(4, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_52 
!
      CALL el_cons_ij(5, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_62 
!
      CALL el_cons_ij(6, 2, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_33 
!
      npos=2*ngeo_strain+1
      CALL el_cons_ij(3, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_13 
!
      CALL el_cons_ij(1, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_23 
!
      CALL el_cons_ij(2, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_43 
!
      CALL el_cons_ij(4, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_53
!
      CALL el_cons_ij(5, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_63 
!
      CALL el_cons_ij(6, 3, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_44
!
      npos=5*ngeo_strain+1
      CALL el_cons_ij(4, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_14
!
      CALL el_cons_ij(1, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_24
!
      CALL el_cons_ij(2, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_34
!
      CALL el_cons_ij(3, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_54 
!
      CALL el_cons_ij(5, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_64 
!
      CALL el_cons_ij(6, 4, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_55 
!
      npos=4*ngeo_strain+1
      CALL el_cons_ij(5, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_15 
!
      CALL el_cons_ij(1, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_25 
!
      CALL el_cons_ij(2, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_35 
!
      CALL el_cons_ij(3, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_45 
!
      CALL el_cons_ij(4, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_65 
!
      CALL el_cons_ij(6, 5, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_66 
!
      npos=3*ngeo_strain+1
      CALL el_cons_ij(6, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_16 
!
      CALL el_cons_ij(1, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_26 
!
      CALL el_cons_ij(2, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_36 
!
      CALL el_cons_ij(3, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_46 
!
      CALL el_cons_ij(4, 6, ngeo_strain, epsil_geo(1,1,npos), &
                                         sigma_geo(1,1,npos), m1)
!  
!  c_56 
!
      CALL el_cons_ij(5, 6, ngeo_strain, epsil_geo(1,1,npos), &
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
!  set on the basis of the Bravais lattice and the Laue class.
!  The stress is estimated and the calculated elastic constants are
!  corrected keeping into account the difference between stress-strain
!  coefficients and second derivatives of the energy with respect to
!  strain. The output of this routine are the stress-strain elastic 
!  constants.
!
!  With reference to the strain indicated in the strain.f90 module the
!  routine expects to find in energy_geo the total energy that 
!  corresponds to the following strain: 
!
!  Laue class         ibrav            strain
!  29,32 T_h,  O_h     1               A   E   F
!  29,32 T_h,  O_h     2,3             A   E   F3
!  19,23 C_6h, D_6h    4               C   E   B1  A   H
!  22    D_4h          6,7             E   C   B   B1  G   H
!  18    C_4h          6,7             E   C   B   B1  G   H   CG
!  20    D_2h          8,9,10,11       C   D   E   B   B1  B2  G   H  I
!  25,27 D_3d,S_6      4               C   E   B1  A   H
!  25,   D_3d,         5               C   E   B1  A   H   CI
!  27,   S_6           5               C   E   B1  A   H   CI  CG
!  16    C_2h          12, 13,         C   D   E   B   B1  B2  G   H  I
!                                      CG  DG  EG  HI
!  16    C_2h          -12,-13,        C   D   E   B   B1  B2  G   H  I
!                                      CH  DH  EH  GI 
!  2     C_i           14              C   D   E   B   B1  B2  G   H   I
!                                      CG  CH  CI  DG  DH  DI  EG  EH  EI
!                                      GH  IH  IG
!
!  For each strain the routine expects the energy for ngeo_strain values of 
!  strain for a total of nwork values.
!  epsilon_geo is in cartesian coordinates in the axis of the unstrained 
!  lattice.
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: epsil_geo(3,3,nwork), omega
REAL(DP), INTENT(IN) :: energy_geo(nwork)
INTEGER,  INTENT(IN) :: nwork, ngeo_strain, ibrav, laue, m1
REAL(DP) :: alpha(m1)
REAL(DP) :: b0, a0, s11, s22, s33, s12, s13, s23, bmat(6,6)
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

      strs=0.0_DP
      strs(1,1)=-press
      strs(2,2)=-press
      strs(3,3)=-press

   CASE (19,23)
!
!  hexagonal case, (C_6h and D_6h)
!
!  C_11=C_22
!
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      s11= alpha(2) / omega
      WRITE(stdout,'("S11=",f15.8," kbar")') s11*ry_kbar
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
      el_con(2,2) = el_con(1,1)
!
! C_33
!
      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s33= alpha(2) / omega
      WRITE(stdout,'("S33=",f15.8," kbar")') s33*ry_kbar
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_13=C_23
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega / 2.0_DP
      WRITE(stdout,'("-(S11+S33)/2",f15.8," kbar",f15.8," kbar")') &
                                 press*ry_kbar, -(s11+s33)*0.5_DP*ry_kbar

      el_con(1,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                                     - el_con(3,3)) * 0.5_DP 
      el_con(3,1) = el_con(1,3)
      el_con(2,3) = el_con(1,3)
      el_con(3,2) = el_con(2,3)
!
!  C_12
!
      base_data=3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega / 3.0_DP
      WRITE(stdout,'("-(2*S11+S33)/3", f15.8," kbar", f15.8, " kbar")') &
                            press*ry_kbar, -(2.0_DP*s11+s33)/3.0_DP*ry_kbar

      el_con(1,2) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) &
                             - 2.0_DP * el_con(1,1) - el_con(3,3) &
                             - 4.0_DP * el_con(1,3) ) * 0.5_DP 
      el_con(2,1) = el_con(1,2)

!
!  C_44=C_55
!
      base_data=4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(5,5) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
      el_con(4,4) = el_con(5,5)

      el_con(6,6) = (el_con(1,1) - el_con(1,2)) * 0.5_DP

      strs=0.0_DP
      strs(1,1)=s11
      strs(2,2)=s11
      strs(3,3)=s33
      CALL correct_for_stress(bmat,el_con,strs)
      el_con=bmat
         
   CASE(25,27)
!
!   Trigonal and hexagonal systems (D_3d and S_6)
!
!  C_11=C_22
!
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      s11= alpha(2) / omega
      WRITE(stdout,'("S11=",f15.8," kbar")') s11*ry_kbar
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
      el_con(2,2) = el_con(1,1)
!
! C_33
!
      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s33= alpha(2) / omega
      WRITE(stdout,'("S33=",f15.8, " kbar")') s33*ry_kbar
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_13=C_23
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega / 2.0_DP
      WRITE(stdout,'("-(S11+S33)/2",f15.8," kbar",f15.8," kbar")') &
                                press*ry_kbar, -(s11+s33)*0.5_DP*ry_kbar

      el_con(1,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                                     - el_con(3,3)) * 0.5_DP 
      el_con(3,1) = el_con(1,3)
      el_con(2,3) = el_con(1,3)
      el_con(3,2) = el_con(2,3)
!
!  C_12
!
      base_data=3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega / 3.0_DP
      WRITE(stdout,'("-(2*S11+S33)/3", f15.8," kbar", f15.8," kbar")') &
                     press*ry_kbar, -(2.0_DP*s11+s33)/3.0_DP*ry_kbar

      el_con(1,2) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - &
                                2.0_DP * el_con(1,1) - el_con(3,3) &
                              - 4.0_DP * el_con(1,3) ) * 0.5_DP 
      el_con(2,1) = el_con(1,2)
!
!  C_44=C_55. The factor 1/4 is due to the definition of epsilon_4
!
      base_data=4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(5,5) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
      el_con(4,4) = el_con(5,5)

      el_con(6,6) = (el_con(1,1) - el_con(1,2)) * 0.5_DP
!
!  C_14 = -C_24 = C_56
!
      IF (ibrav==5) THEN
         base_data=5*ngeo_strain+1
         CALL el_cons_ij_ene(1, 1, 'C_14', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         el_con(1,4) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                                       - el_con(4,4))*0.5_DP
         el_con(4,1) = el_con(1,4)

         el_con(2,4) = -el_con(1,4) 
         el_con(4,2) = el_con(2,4)

         el_con(5,6) = el_con(1,4) 
         el_con(6,5) = el_con(5,6)

         IF (laue==27) THEN
!
!  C_25 = -C_15 = C_46
!
            base_data=6*ngeo_strain+1
            CALL el_cons_ij_ene(2, 2, 'C_25', ngeo_strain, &
                epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

            el_con(2,5) = ( 1.0_DP / omega * (2.0_DP * alpha(3) ) &
                                     - el_con(2,2) - el_con(5,5))*0.5_DP
            el_con(5,2) = el_con(2,5)

            el_con(1,5) = -el_con(2,5) 
            el_con(5,1) = el_con(1,5)

            el_con(4,6) = el_con(2,5) 
            el_con(6,4) = el_con(4,6)
         ENDIF
      ENDIF

      strs=0.0_DP
      strs(1,1)=s11
      strs(2,2)=s11
      strs(3,3)=s33
      CALL correct_for_stress(bmat,el_con,strs)
      el_con=bmat

   CASE(18,22)
!
!  tetragonal case, (C_4h and D_4h Laue classes)
!
!  C_33 
!
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      s33=alpha(2) / omega
      WRITE(stdout,'("S33=",f15.9," kbar")') s33*ry_kbar
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_11 = C_22
!
      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      s11=alpha(2) / omega
      WRITE(stdout,'("S11=",f15.9," kbar")') s11*ry_kbar
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
      el_con(2,2) = el_con(1,1)
!
!  C_12
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega
      WRITE(stdout,'("-S11=", f15.8," kbar", f15.8," kbar")') &
                                         press*ry_kbar, -s11*ry_kbar
      el_con(1,2) = 0.5_DP/omega*( 2.0_DP*alpha(3) ) - el_con(1,1) 
      el_con(2,1) = el_con(1,2)
!
!  C_13=C_23
!
      base_data=3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(2 S11+S33)/2",f15.8," kbar",f15.8," kbar" )') &
                       press*ry_kbar, -(2.0_DP*s11+s33)/2.0_DP*ry_kbar
      el_con(1,3) = (1.0_DP/omega*(2.0_DP*alpha(3)) - el_con(3,3) &
                                              -el_con(1,1)) * 0.5_DP 
      el_con(3,1) = el_con(1,3)
      el_con(2,3) = el_con(1,3)
      el_con(3,2) = el_con(2,3)
!
!  C_66 the factor 1/4 is for the definition of \epsilon_4
!
      base_data=4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_66', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(6,6) = 0.25_DP/omega*(2.0_DP*alpha(3)) 
!
!  C_44=C_55
!
      base_data = 5 * ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      el_con(5,5) = 0.25_DP/omega*(2.0_DP*alpha(3))
      el_con(4,4) = el_con(5,5)
      IF (laue==18) THEN
!
!   C_16
!
         base_data=6*ngeo_strain+1
         CALL el_cons_ij_ene(1, 1, 'C_16', ngeo_strain, &
              epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         el_con(1,6) = (1.0_DP/omega*(2.0_DP*alpha(3)) - el_con(1,1) &
                                                       - el_con(6,6) )* 0.5_DP
         el_con(6,1) = el_con(1,6)
         el_con(2,6) = -el_con(1,6)
         el_con(6,2) = el_con(2,6)
      ENDIF

      strs=0.0_DP
      strs(1,1)=s11
      strs(2,2)=s11
      strs(3,3)=s33
      CALL correct_for_stress(bmat,el_con,strs)
      el_con=bmat

   CASE(20)
!
!  Orthorombic case (D_2h Laue class)  
!
!  C_11
!
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      s11=alpha(2) / omega
      WRITE(stdout,'("S11=",f15.8," kbar")') s11*ry_kbar
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_22
!
      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_22', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s22=alpha(2) / omega
      WRITE(stdout,'("S22=",f15.9," kbar")') s22*ry_kbar
      el_con(2,2) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_33
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s33=alpha(2) / omega
      WRITE(stdout,'("S33= ",f15.9," kbar")') s33*ry_kbar
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_12
!
      base_data = 3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S12)/2",f15.8," kbar",f15.8," kbar")') &
                                   press*ry_kbar, -(s11+s22)/2.0_DP*ry_kbar
      el_con(1,2) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(2,2) ) * 0.5_DP 
      el_con(2,1)=el_con(1,2)
!
!  C_13
!
      base_data = 4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S33)/2",f15.8," kbar",f15.8," kbar")') &
                                  press*ry_kbar, -(s11+s33)/2.0_DP*ry_kbar
      el_con(1,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(3,3) ) * 0.5_DP 
      el_con(3,1)=el_con(1,3)
!
!  C_23
!
      base_data = 5*ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_23', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S33)/2",f15.8," kbar",f15.8," kbar")') &
                               press*ry_kbar, -(s22+s33)/2.0_DP*ry_kbar
      el_con(2,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                   - el_con(3,3) ) * 0.5_DP 
      el_con(3,2)=el_con(2,3)
!
!  C_66. The factor 1/4 is due to the definition of \epsilon_6
!
      base_data = 6*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_66', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      el_con(6,6) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
!
!  C_55
!
      base_data = 7*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      el_con(5,5) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
!
!  C_44
!
      base_data = 8*ngeo_strain+1
      CALL el_cons_ij_ene(2, 3, 'C_44', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      el_con(4,4) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 

      strs=0.0_DP
      strs(1,1)=s11
      strs(2,2)=s22
      strs(3,3)=s33
      CALL correct_for_stress(bmat,el_con,strs)
      el_con=bmat

   CASE(16)
!
!   Monoclinic (C_2h Laue class)
!
!
!  C_11
!
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      s11=alpha(2) / omega
      WRITE(stdout,'("S11=",f15.8," kbar")') s11*ry_kbar
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_22
!
      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_22', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s22=alpha(2) / omega
      WRITE(stdout,'("S22=",f15.8," kbar")') s22*ry_kbar
      el_con(2,2) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_33
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s33=alpha(2) / omega
      WRITE(stdout,'("S33=",f15.8," kbar")') s33*ry_kbar
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_12
!
      base_data = 3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S22)/2",f15.8," kbar",f15.8," kbar")') &
                              press*ry_kbar, -(s11+s22)/2.0_DP*ry_kbar
      el_con(1,2) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(2,2) ) * 0.5_DP 
      el_con(2,1)=el_con(1,2)
!
!  C_13
!
      base_data = 4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S33)/2",f15.8," kbar", f15.8," kbar")') &
                                  press*ry_kbar, -(s11+s33)/2.0_DP*ry_kbar
      el_con(1,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(3,3) ) * 0.5_DP 
      el_con(3,1)=el_con(1,3)
!
!  C_23
!
      base_data = 5*ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_23', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega
      WRITE(stdout,'("-(S22+S33)/2",f15.8," kbar",f15.8," kbar")') &
                              press*ry_kbar, -(s22+s33)/2.0_DP*ry_kbar
      el_con(2,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                   - el_con(3,3) ) * 0.5_DP 
      el_con(3,2)=el_con(2,3)
!
!  C_66
!
      base_data = 6*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_66', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s12=- alpha(2) / omega
      WRITE(stdout,'("S12=",f15.8," kbar")') s12*ry_kbar
      el_con(6,6) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
!
!  C_55
!
      base_data = 7*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s13=- alpha(2) / omega
      WRITE(stdout,'("S13=",f15.8," kbar")') s13*ry_kbar
      el_con(5,5) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
!
!  C_44
!
      base_data = 8*ngeo_strain+1
      CALL el_cons_ij_ene(2, 3, 'C_44', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s23=- alpha(2) / omega
      WRITE(stdout,'("S23=",f15.8," kbar")') s23*ry_kbar
      el_con(4,4) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 

      IF (ibrav>0) THEN
!
!   Monoclinic lattice: unique axis c
!
!
!  C_16 
!
         base_data = 9*ngeo_strain+1
         CALL el_cons_ij_ene(1, 1, 'C_16', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("-(S11+S12)",f15.8," kbar",f15.8," kbar")') &
                   press*ry_kbar, -(s11+s12)*ry_kbar 
         el_con(1,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                               - el_con(6,6)) *0.5_DP
         el_con(6,1) = el_con(1,6)
!
!  C_26
!
         base_data = 10*ngeo_strain+1
         CALL el_cons_ij_ene(2, 2, 'C_26', ngeo_strain, &
                epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("-(S22+S12)",f15.8," kbar",f15.8," kbar")') &
                                           press*ry_kbar, -(s22+s12)*ry_kbar 
         el_con(2,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                                  - el_con(6,6)) *0.5_DP
         el_con(6,2) = el_con(2,6)
!
!  C_36 
!
         base_data = 11*ngeo_strain+1
         CALL el_cons_ij_ene(3, 3, 'C_36', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("-(S33+S12)",f15.8," kbar",f15.8," kbar")') & 
                                          press*ry_kbar, -(s33+s12)*ry_kbar 
         el_con(3,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(3,3) &
                                               - el_con(6,6)) *0.5_DP
         el_con(6,3) = el_con(3,6)
!
!  C_45
!
         base_data = 12*ngeo_strain+1
         CALL el_cons_ij_ene(1, 3, 'C_45', ngeo_strain, &
               epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("(S13+S23)",f15.8," kbar",f15.8," kbar")') &
                                        press*ry_kbar, -(s13+s23)*ry_kbar 
         el_con(4,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(4,4) &
                                                  - el_con(5,5) ) * 0.5_DP
         el_con(5,4) = el_con(4,5)
      ELSE
!
!   Monoclinic case: unique axis b
!
!
!  C_15 
!
         base_data = 9*ngeo_strain+1
         CALL el_cons_ij_ene(1, 1, 'C_15', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("-(S11+S13)",f15.8," kbar",f15.8," kbar")') &
                                         press*ry_kbar, -(s11+s13)*ry_kbar 
         el_con(1,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                               - el_con(5,5)) *0.5_DP
         el_con(5,1) = el_con(1,5)
!
!  C_25  
!
         base_data = 10*ngeo_strain+1
         CALL el_cons_ij_ene(2, 2, 'C_25', ngeo_strain, &
               epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("-(S22+S13)",f15.8," kbar",f15.8," kbar")') &
                                         press*ry_kbar, -(s22+s13)*ry_kbar 
         el_con(2,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                                - el_con(5,5)) *0.5_DP
         el_con(5,2) = el_con(2,5)
!
!  C_35
!
         base_data = 11*ngeo_strain+1
         CALL el_cons_ij_ene(1, 1, 'C_35', ngeo_strain, &
                epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("-(S33+S13)",f15.8," kbar",f15.8," kbar")') &
                                         press*ry_kbar, -(s33+s13)*ry_kbar 
         el_con(3,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(3,3) &
                                                  - el_con(5,5)) *0.5_DP
         el_con(5,3) = el_con(3,5)
!
!  C_46 
!
         base_data = 12*ngeo_strain+1
         CALL el_cons_ij_ene(1, 2, 'C_46', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

         press=- alpha(2) / omega
         WRITE(stdout,'("-(S23+S12)",f15.8," kbar",f15.8," kbar")') &
                                         press*ry_kbar, -(s23+s12)*ry_kbar 
         el_con(4,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(4,4) &
                                                  - el_con(6,6)) *0.5_DP
         el_con(6,4) = el_con(4,6)
      ENDIF

      strs=0.0_DP
      strs(1,1)=s11
      strs(2,2)=s22
      strs(3,3)=s33
      strs(1,2)=s12
      strs(2,1)=s12
      strs(1,3)=s13
      strs(3,1)=s13
      strs(2,3)=s23
      strs(3,2)=s23
      CALL correct_for_stress(bmat,el_con,strs)
      el_con=bmat

   CASE(2)
!
!   Triclinic (C_i Laue class)
!
!
!  C_11
!
      CALL el_cons_ij_ene(1, 1, 'C_11', ngeo_strain, epsil_geo, &
                                                 energy_geo, alpha, m1)

      s11=alpha(2) / omega
      WRITE(stdout,'("S11=",f15.8," kbar")') s11*ry_kbar
      el_con(1,1) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_22
!
      base_data=ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_22', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s22=alpha(2) / omega
      WRITE(stdout,'("S22=",f15.8," kbar")') s22*ry_kbar
      el_con(2,2) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_33
!
      base_data=2*ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_33', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s33=alpha(2) / omega
      WRITE(stdout,'("S33=",f15.8," kbar")') s33*ry_kbar
      el_con(3,3) = 1.0_DP / omega * ( 2.0_DP * alpha(3) )
!
!  C_12
!
      base_data = 3*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_12', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S12)/2",f15.8," kbar",f15.8," kbar")') &
                                  press*ry_kbar, -(s11+s22)/2.0_DP*ry_kbar
      el_con(1,2) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(2,2) ) * 0.5_DP 
      el_con(2,1)=el_con(1,2)
!
!  C_13
!
      base_data = 4*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_13', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S33)/2",f15.8," kbar",f15.8," kbar")') &
                                 press*ry_kbar, -(s11+s33)/2.0_DP*ry_kbar
      el_con(1,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                   - el_con(3,3) ) * 0.5_DP 
      el_con(3,1)=el_con(1,3)
!
!  C_23
!
      base_data = 5*ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_23', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)
      press=- alpha(2) / omega
      WRITE(stdout,'("-(S22+S33)/2",f15.8," kbar",f15.8," kbar")') &
                                  press*ry_kbar, -(s22+s33)/2.0_DP*ry_kbar
      el_con(2,3) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                   - el_con(3,3) ) * 0.5_DP 
      el_con(3,2)=el_con(2,3)
!
!  C_66
!
      base_data = 6*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_66', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s12=- alpha(2) / omega
      WRITE(stdout,'("S12=",f15.8," kbar")') s12*ry_kbar
      el_con(6,6) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
!
!  C_55
!
      base_data = 7*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_55', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s13=- alpha(2) / omega
      WRITE(stdout,'("S13=",f15.8," kbar")') s13*ry_kbar
      el_con(5,5) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
!
!  C_44
!
      base_data = 8*ngeo_strain+1
      CALL el_cons_ij_ene(2, 3, 'C_44', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      s23=- alpha(2) / omega
      WRITE(stdout,'("S12=",f15.8," kbar")') s23*ry_kbar
      el_con(4,4) = 0.25_DP / omega * ( 2.0_DP * alpha(3) ) 
!
!  C_16
!
      base_data = 9*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_16', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S12)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s11+s12)*ry_kbar 
      el_con(1,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                               - el_con(6,6)) *0.5_DP
      el_con(6,1) = el_con(1,6)
!
!  C_15
!
      base_data = 10*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_15', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S13)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s11+s13)*ry_kbar 
      el_con(1,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                               - el_con(5,5)) *0.5_DP
      el_con(5,1) = el_con(1,5)
!
!  C_14
!
      base_data = 11*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_14', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S11+S23)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s11+s23)*ry_kbar 
      el_con(1,4) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(1,1) &
                                               - el_con(4,4) ) * 0.5_DP
      el_con(4,1) = el_con(1,4)
!
!  C_26
!
      base_data = 12*ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_26', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S22+S12)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s22+s12)*ry_kbar 
      el_con(2,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                               - el_con(6,6)) *0.5_DP
      el_con(6,2) = el_con(2,6)
!
!  C_25
!
      base_data = 13*ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_25', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S22+S13)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s22+s13)*ry_kbar 
      el_con(2,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                               - el_con(5,5)) *0.5_DP
      el_con(5,2) = el_con(2,5)
!
!  C_24
!
      base_data = 14*ngeo_strain+1
      CALL el_cons_ij_ene(2, 2, 'C_24', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S22+S23)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s22+s23)*ry_kbar 
      el_con(2,4) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(2,2) &
                                               - el_con(4,4) ) * 0.5_DP
      el_con(4,2) = el_con(2,4)
!
!  C_36
!
      base_data = 15*ngeo_strain+1
      CALL el_cons_ij_ene(3, 3, 'C_36', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S33+S12)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s33+s12)*ry_kbar 
      el_con(3,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(3,3) &
                                               - el_con(6,6)) *0.5_DP
      el_con(6,3) = el_con(3,6)
!
!  C_35
!
      base_data = 16*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_35', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S33+S13)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s33+s13)*ry_kbar 
      el_con(3,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(3,3) &
                                               - el_con(5,5)) *0.5_DP
      el_con(5,3) = el_con(3,5)
!
!  C_34
!
      base_data = 17*ngeo_strain+1
      CALL el_cons_ij_ene(1, 1, 'C_34', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S33+S23)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s33+s23)*ry_kbar
      el_con(3,4) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(3,3) &
                                               - el_con(4,4) ) * 0.5_DP
      el_con(4,3) = el_con(3,4)
!
!  C_56
!
      base_data = 18*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_56', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S12+S23)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s12+s23)*ry_kbar
      el_con(5,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(5,5) &
                                               - el_con(6,6)) *0.5_DP
      el_con(6,5) = el_con(5,6)
!
!  C_46
!
      base_data = 19*ngeo_strain+1
      CALL el_cons_ij_ene(1, 2, 'C_46', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S23+S12)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s23+s12)*ry_kbar
      el_con(4,6) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(4,4) &
                                               - el_con(6,6)) *0.5_DP
      el_con(6,4) = el_con(4,6)
!
!  C_45
!
      base_data = 20*ngeo_strain+1
      CALL el_cons_ij_ene(1, 3, 'C_45', ngeo_strain, &
             epsil_geo(1,1,base_data), energy_geo(base_data), alpha, m1)

      press=- alpha(2) / omega
      WRITE(stdout,'("-(S13+S23)",f15.8," kbar",f15.8," kbar")') &
                                          press*ry_kbar, -(s13+s23)*ry_kbar
      el_con(4,5) = (1.0_DP / omega * ( 2.0_DP * alpha(3) ) - el_con(4,4) &
                                               - el_con(5,5) ) * 0.5_DP
      el_con(5,4) = el_con(4,5)

      strs=0.0_DP
      strs(1,1)=s11
      strs(2,2)=s22
      strs(3,3)=s33
      strs(1,2)=s12
      strs(2,1)=s12
      strs(1,3)=s13
      strs(3,1)=s13
      strs(2,3)=s23
      strs(3,2)=s23
      CALL correct_for_stress(bmat,el_con,strs)
      el_con=bmat

   CASE DEFAULT
      CALL errore('compute_elastic_constants_ene',&
                                   'Case not yet available',1)
END SELECT
press = (strs(1,1) + strs(2,2) + strs(3,3))/ 3.0_DP
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

SUBROUTINE elastic_constants_from_compliances(cmn,smn)
!
! This routine receives as input the elastic compliances matrix smn
! and computes its inverse, the elastic constants matrix
!
USE kinds, ONLY : DP
USE matrix_inversion, ONLY : invmat
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: smn(6,6)
REAL(DP), INTENT(INOUT) :: cmn(6,6)

CALL invmat(6, smn, cmn)

RETURN
END SUBROUTINE elastic_constants_from_compliances

SUBROUTINE el_cons_ij(pq, mn, ngeo, epsil_geo, sigma_geo, m1)
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit, write_poly
USE voigt, ONLY : voigt_extract_indices

IMPLICIT NONE
INTEGER, INTENT(IN) :: mn, pq, ngeo
REAL(DP), INTENT(IN) :: epsil_geo(3,3,ngeo), sigma_geo(3,3,ngeo)
INTEGER :: igeo, m, n, p, q, mnin, pqin
INTEGER :: m1                  ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo)

WRITE(stdout,'(/,20x,40("-"),/)')
mnin=mn
CALL voigt_extract_indices(m,n,mnin)
pqin=pq
CALL voigt_extract_indices(p,q,pqin)
WRITE(stdout,'(/,5x,"Elastic constant ",2i5)') pq, mn
WRITE(stdout,'(/,10x,"strain",7x,"stress (kbar)")') 
DO igeo=1,ngeo
   x(igeo)=epsil_geo(m,n,igeo)
   y(igeo)=sigma_geo(p,q,igeo)
   WRITE(stdout,'(2f18.10)') x(igeo), y(igeo)*ry_kbar
ENDDO
CALL polyfit( x, y, ngeo, alpha, m1-1 )
CALL write_poly(alpha,m1-1)
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
USE polyfit_mod, ONLY : polyfit, write_poly

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
CALL polyfit( x, y, ngeo, alpha, m1-1 )
CALL write_poly(alpha,m1-1)

RETURN
END SUBROUTINE el_cons_ij_ene

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

SUBROUTINE correct_for_stress(bmat, el_cons_, stres)
!
!  This routine receives as input the second derivatives of the 
!  energy with respect to strain in el_cons(6,6) and gives as output
!  the stress-strain elastic constants contained in the matrix 
!  bmat(6,6). It receives also the stress of the unperturbed system
!  in stres(3,3)
!
USE voigt, ONLY : voigt_extract_indices
IMPLICIT NONE

REAL(DP), INTENT(IN) :: el_cons_(6,6)
REAL(DP), INTENT(IN) :: stres(3,3)
REAL(DP), INTENT(INOUT) :: bmat(6,6)

INTEGER :: ij, kl, ijm, klm, i, j, k, l

DO ij=1,6
   ijm=ij
   CALL voigt_extract_indices(i,j,ijm)
   DO kl=1,6
      klm=kl
      CALL voigt_extract_indices(k,l,klm)
      bmat(ij,kl)= el_cons_(ij,kl) - 0.5_DP * (2.0_DP*stres(i,j)*delta(k,l) &
                    -0.5_DP*(stres(j,l)*delta(i,k) +      &
                             stres(j,k)*delta(i,l) +      &
                             stres(i,l)*delta(j,k) +      &
                             stres(i,k)*delta(j,l) ) )
   ENDDO
ENDDO

RETURN
END SUBROUTINE correct_for_stress

SUBROUTINE correct_for_pressure(bmat, el_cons_, pressure)
!
!  This routine receives as input the second derivatives of the 
!  energy with respect to strain in el_cons(6,6) and gives as output
!  the stress-strain elastic constants contained in the matrix 
!  bmat(3,3,3,3). It receives also the pressure of the unperturbed system
!  and assumes that this is the only stress present on the system
!
IMPLICIT NONE

REAL(DP), INTENT(IN) :: el_cons_(6,6)
REAL(DP), INTENT(IN) :: pressure
REAL(DP), INTENT(INOUT) :: bmat(6,6)

REAL(DP) :: stres(3,3)
INTEGER :: i

stres=0.0_DP
DO i=1,3
   stres(i,i)=-pressure
ENDDO

CALL correct_for_stress(bmat, el_cons_, stres)

RETURN
END SUBROUTINE correct_for_pressure

INTEGER FUNCTION delta(i,j)
!
!   delta function between two integers
!
INTEGER, INTENT(IN) :: i, j

delta=0
IF (i==j) delta=1

RETURN
END FUNCTION delta

SUBROUTINE print_el_cons_info(elastic_algorithm, laue, ibrav, ngeo_strain)

IMPLICIT NONE
INTEGER, INTENT(IN) :: laue, ibrav, ngeo_strain
CHARACTER(LEN=*), INTENT(IN) :: elastic_algorithm

INTEGER :: nstrain

IF (elastic_algorithm=='standard'.OR.elastic_algorithm=='advanced') THEN
   SELECT CASE (laue) 
!
!  cubic T_h (m-3), O_h (m-3m)
!
      CASE (29,32)
         WRITE(stdout,'(/,5x,"It requires two strains: e1 and e4")') 
         nstrain=2
      CASE (27,25)
!
!  trigonal S_6 (-3), D_3d (-3m)
!
         WRITE(stdout,'(/,5x,"It requires three strains: e1, e3, and e4.")') 
         nstrain=3
      CASE (19,23)
! 
!  hexagonal C_6h (6/m), D_6h (6/mmm)
!
         WRITE(stdout,'(/,5x,"It requires three strains: e1, e3, and e5.")') 
         nstrain=3
      CASE (18,22)
!
!  tetragonal C_4h (4/m),  tetragonal D_4h (4/mmm)
!
         WRITE(stdout,'(/,5x,"It requires four strains: e1, e3, e4, e6")') 
         nstrain=4
      CASE (16,20,2)
!
!    monoclinic case, class C_2h (2/m), orthorhombic D_2h (mmm) 
!    triclinic case, class C_i(-1)
!
         WRITE(stdout,'(/,5x,"It requires all six strains.")') 
         nstrain=6

   END SELECT
ELSE
   SELECT CASE (laue) 
!
!  cubic T_h (m-3), O_h (m-3m)
!
      CASE (29,32)
         WRITE(stdout,'(/,5x,"It requires three strains.")') 
         nstrain=3
      CASE (25)
         IF (ibrav==4) THEN
            WRITE(stdout,'(/,5x,"It requires five strains.")') 
            nstrain=5
         ELSE
            WRITE(stdout,'(/,5x,"It requires six strains.")') 
            nstrain=6
         ENDIF
      CASE (27)
         IF (ibrav==4) THEN
            WRITE(stdout,'(/,5x,"It requires five strains.")') 
            nstrain=5
         ELSE
            WRITE(stdout,'(/,5x,"It requires seven strains.")') 
            nstrain=7
         ENDIF
      CASE (19,23)
!
!  hexagonal C_6h (6/m), D_6h (6/mmm)
!
         WRITE(stdout,'(/,5x,"It requires five strains.")') 
         nstrain=5
      CASE (22)
!
!  tetragonal D_4h (4/mmm)
!
         WRITE(stdout,'(/,5x,"It requires six strains.")') 
         nstrain=6
      CASE (18)
!
!  tetragonal C_4h (4/m).
!
         WRITE(stdout,'(/,5x,"It requires seven strains.")') 
         nstrain=7
      CASE (20)
!
!  orthorhombic D_2h (mmm) 
!
         WRITE(stdout,'(/,5x,"It requires nine strains.")') 
         nstrain=9
      CASE (16)
!
!    monoclinic case, class C_2h (2/m), orthorhombic D_2h (mmm) 
!
         WRITE(stdout,'(/,5x,"It requires 13 strains.")') 
         nstrain=13
      CASE (2)
!
!    triclinic case, class C_i (-1)
!
         WRITE(stdout,'(/,5x,"It requires 21 strains.")') 
         nstrain=21
   END SELECT
ENDIF

WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                         nstrain*ngeo_strain 
RETURN
END SUBROUTINE print_el_cons_info

SUBROUTINE write_el_cons_on_file(temp, ntemp, ibrav, laue, el_cons_t, b0, &
                                                             filename, iflag)
!
!  iflag=0 writes the elastic constants
!  iflag=1 writes the elastic compliances
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, laue, iflag
REAL(DP), INTENT(IN) :: temp(ntemp), el_cons_t(6,6,ntemp), b0(ntemp)
REAL(DP) :: b0_t(ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_el_cons, ios
INTEGER :: find_free_unit

iu_el_cons=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_el_cons, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_el_cons_on_file','opening elastic constants (T) file',&
                                                             ABS(ios))
!
! Plot the compressibility togheter with the elastic complinances
!

IF (iflag/=0) THEN
   b0_t=1.0_DP/b0
ELSE
   b0_t=b0
ENDIF

IF (meta_ionode) THEN
   SELECT CASE (laue)
      CASE(29,32)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
                  & 13x, "     C_12 ", 13x, "     C_44 ")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
                  & 13x, "     S_12 ", 13x, "     S_44 ")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,4e20.12)') temp(itemp), b0_t(itemp), &
                 el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                 el_cons_t(4,4,itemp)
         ENDDO
      CASE(25)
!
!     D_3d
!            
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
                  & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, &
                       &" C_44 ", 13x, " C_14")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
                  & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, &
                       &" S_44 ", 13x, " S_14")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,7e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(1,4,itemp)
         ENDDO
      CASE(27)
!
!     S_6
!            
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", " C_11 ", 13x, &
                  &" C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44", 13x, &
                  &" C_14", 13x, "C_25" )')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", " S_11 ", 13x, &
                  &" S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44", 13x, &
                  &" S_14", 13x, "S_25" )')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,8e20.12)')  temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(1,4,itemp), &
                  el_cons_t(2,5,itemp)
         ENDDO
      CASE(19,23)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", 13x,&
                  &" C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44" )')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", 13x,&
                  &" S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44" )')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,6e20.12)') temp(itemp),  b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp)
         ENDDO
      CASE(22)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", 13x,&
                  &" C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44", 13x, & 
                  &" C_66 " )')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", 13x,&
                  &" S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44", 13x, & 
                  &" S_66 " )')
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,7e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp)
         ENDDO
      CASE(20)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B", 13x, " C_11 ", 13x,&
                  &" C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", 13x,&
                  &" C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K", 13x, " S_11 ", 13x,&
                  &" S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", 13x,&
                  &" S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,10e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                  el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                  el_cons_t(6,6,itemp)
         ENDDO
      CASE(18)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
                  & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44", &
                  & 13x, " C_66 ", 13x, " C_16 ")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
                  & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44", &
                  & 13x, " S_66 ", 13x, " S_16 ")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,8e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp), &
                  el_cons_t(1,6,itemp)
         ENDDO
      CASE(16)
         IF (ibrav < 0) THEN
            !
            !  b unique
            !
            IF (iflag==0) THEN
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
               & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", &
               & 13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ", &
               & 13x, " C_15 ", 13x, " C_25 ", 13x, " C_35 ", 13x, " C_46 ")')
            ELSE
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
               & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", &
               & 13x, " S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ", &
               & 13x, " S_15 ", 13x, " S_25 ", 13x, " S_35 ", 13x, " S_46 ")')
            ENDIF
            DO itemp=2,ntemp-1
               WRITE(iu_el_cons,'(e16.8,14e20.12)') temp(itemp), b0_t(itemp),&
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,5,itemp), &
                     el_cons_t(2,5,itemp), el_cons_t(3,5,itemp), &
                     el_cons_t(4,6,itemp)
            ENDDO
         ELSE
            !
            !  c unique
            !
            IF (iflag==0) THEN
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ",13x," C_11 ",  &
               & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", &
               & 13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ", &
               & 13x, " C_16 ", 13x, " C_26 ", 13x, " C_36 ", 13x, " C_45 ")')
            ELSE
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ",13x," S_11 ",  &
               & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", &
               & 13x, " S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ", &
               & 13x, " S_16 ", 13x, " S_26 ", 13x, " S_36 ", 13x, " S_45 ")')
            ENDIF
            DO itemp=2,ntemp-1
               WRITE(iu_el_cons,'(e16.8,14e20.12)') temp(itemp), b0_t(itemp),&
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,6,itemp), &
                     el_cons_t(2,6,itemp), el_cons_t(3,6,itemp), &
                     el_cons_t(4,5,itemp)
            ENDDO
         ENDIF
      CASE(2)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
               & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", &
               & 13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ", &
               & 13x, " C_14 ", 13x, " C_15 ", 13x, " C_16 ", 13x, " C_24 ", &
               & 13x, " C_25 ", 13x, " C_26 ", 13x, " C_34 ", 13x, " C_35 ", &
               & 13x, " C_36 ", 13x, " C_45 ", 13x, " C_46 ", 13x, " C_56 ")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
               & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", &
               & 13x, " S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ", &
               & 13x, " S_14 ", 13x, " S_15 ", 13x, " S_16 ", 13x, " S_24 ", &
               & 13x, " S_25 ", 13x, " S_26 ", 13x, " S_34 ", 13x, " S_35 ", &
               & 13x, " S_36 ", 13x, " S_45 ", 13x, " S_46 ", 13x, " S_56 ")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,24e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                  el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                  el_cons_t(6,6,itemp), el_cons_t(1,4,itemp), &
                  el_cons_t(1,5,itemp), el_cons_t(1,6,itemp), &
                  el_cons_t(2,4,itemp), el_cons_t(2,5,itemp), &
                  el_cons_t(2,6,itemp), el_cons_t(3,4,itemp), &
                  el_cons_t(3,5,itemp), el_cons_t(3,6,itemp), &
                  el_cons_t(4,5,itemp), el_cons_t(4,6,itemp), &
                  el_cons_t(5,6,itemp)
         ENDDO
   END SELECT
   CLOSE(iu_el_cons)
ENDIF

RETURN
END SUBROUTINE write_el_cons_on_file

!
! Copyright (C) 2019 Cristiano Malica
!
SUBROUTINE read_el_cons_from_file(temp, ntemp, ibrav, laue, el_cons_t, b0_t, &
                                                             filename)
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast

IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, laue
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: b0_t(ntemp), el_cons_t(6,6,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

REAL(DP) :: rdum

INTEGER :: itemp, iu_el_cons, ios
INTEGER :: find_free_unit

iu_el_cons=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_el_cons, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('read_el_cons_on_file','opening elastic constants (T) file',&
                                                             ABS(ios))

el_cons_t=0.0_DP

IF (meta_ionode) THEN
   SELECT CASE (laue)
      CASE(29,32)
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,4e20.12)') rdum, b0_t(itemp), &
                 el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                 el_cons_t(4,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature', 1) 
         ENDDO
      CASE(25)
!
!     D_3d
!            
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,7e20.12)') rdum, b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(1,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature', 1) 
         ENDDO
      CASE(27)
!
!     S_6
!            
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,8e20.12)') rdum, b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(1,4,itemp), &
                  el_cons_t(2,5,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature', 1) 
         ENDDO
      CASE(19,23)
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,6e20.12)') rdum,  b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature', 1) 
         ENDDO
      CASE(22)
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,7e20.12)') rdum, b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature', 1) 
         ENDDO
      CASE(20)
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,10e20.12)') rdum, b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                  el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                  el_cons_t(6,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature', 1) 
         ENDDO
      CASE(18)
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,8e20.12)') rdum, b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp), &
                  el_cons_t(1,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature', 1) 
         ENDDO
      CASE(16)
         IF (ibrav < 0) THEN
            !
            !  b unique
            !
            READ(iu_el_cons,*) 
            DO itemp=2,ntemp-1
               READ(iu_el_cons,'(e16.8,14e20.12)') rdum, b0_t(itemp),&
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,5,itemp), &
                     el_cons_t(2,5,itemp), el_cons_t(3,5,itemp), &
                     el_cons_t(4,6,itemp)
               IF (ABS(rdum-temp(itemp))>1D-5) &
                  CALL errore('read_el_cons_on_file','uncorrect temperature',&
                                                                         1) 
            ENDDO
         ELSE
            !
            !  c unique
            !
            READ(iu_el_cons,*) 
            DO itemp=2,ntemp-1
               READ(iu_el_cons,'(e16.8,14e20.12)') rdum, b0_t(itemp),&
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,6,itemp), &
                     el_cons_t(2,6,itemp), el_cons_t(3,6,itemp), &
                     el_cons_t(4,5,itemp)
               IF (ABS(rdum-temp(itemp))>1D-5) &
                  CALL errore('read_el_cons_on_file','uncorrect temperature',&
                                                                         1) 
            ENDDO
         ENDIF
      CASE(2)
         READ(iu_el_cons,*) 
         DO itemp=2,ntemp-1
            READ(iu_el_cons,'(e16.8,24e20.12)') rdum, b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                  el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                  el_cons_t(6,6,itemp), el_cons_t(1,4,itemp), &
                  el_cons_t(1,5,itemp), el_cons_t(1,6,itemp), &
                  el_cons_t(2,4,itemp), el_cons_t(2,5,itemp), &
                  el_cons_t(2,6,itemp), el_cons_t(3,4,itemp), &
                  el_cons_t(3,5,itemp), el_cons_t(3,6,itemp), &
                  el_cons_t(4,5,itemp), el_cons_t(4,6,itemp), &
                  el_cons_t(5,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_el_cons_on_file','uncorrect temperature',1) 
         ENDDO
   END SELECT
   CLOSE(iu_el_cons)
   CALL expand_el_cons(el_cons_t, laue, ibrav, ntemp, temp)
ENDIF
CALL mp_bcast(el_cons_t, meta_ionode_id, world_comm)

RETURN
END SUBROUTINE read_el_cons_from_file

SUBROUTINE expand_el_cons(el_cons_t, laue, ibrav, ntemp, temp)

USE kinds,      ONLY : DP

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ibrav, laue, ntemp
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: el_cons_t(6,6,ntemp)
INTEGER :: itemp, i, j

SELECT CASE (laue)

   CASE(29,32)
!
!  cubic T_h (m-3), O_h (m-3m)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(3,3,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(1,3,itemp)=el_cons_t(1,2,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,2,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(6,6,itemp)=el_cons_t(4,4,itemp) 
      END DO

   CASE(25) 
!
!  trigonal D_3d (-3m)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(2,4,itemp)=-el_cons_t(1,4,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(5,6,itemp)=el_cons_t(1,4,itemp)
         el_cons_t(6,6,itemp)=(el_cons_t(1,1,itemp)-&
                         el_cons_t(1,2,itemp))/2.0_DP
      END DO

   CASE (27)
!
!  trigonal S_6 (-3)
!
      DO itemp=2,ntemp-1
         el_cons_t(1,5,itemp)=-el_cons_t(2,5,itemp)
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(2,4,itemp)=-el_cons_t(1,4,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(5,6,itemp)=el_cons_t(1,4,itemp)
         el_cons_t(6,6,itemp)=(el_cons_t(1,1,itemp)-&
                         el_cons_t(1,2,itemp))/2.0_DP
      END DO

   CASE (19,23)
!
!  hexagonal C_6h (6/m), D_6h (6/mmm)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
         el_cons_t(6,6,itemp)=(el_cons_t(1,1,itemp)-&
                         el_cons_t(1,2,itemp))/2.0_DP
      END DO

   CASE(22)
!
!  tetragonal D_4h (4/mmm)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
      END DO

   CASE(20)
!
!  orthorhombic D_2h (mmm)
!
      !There are no other elastic constants 
      !dependent from those read.

   CASE(18)
!
!  tetragonal C_4h (4/m)
!
      DO itemp=2,ntemp-1
         el_cons_t(2,2,itemp)=el_cons_t(1,1,itemp)
         el_cons_t(2,3,itemp)=el_cons_t(1,3,itemp)
         el_cons_t(2,6,itemp)=-el_cons_t(1,6,itemp)
         el_cons_t(5,5,itemp)=el_cons_t(4,4,itemp)
      END DO

   CASE(16)
!
!    monoclinic case, class C_2h (2/m) 
!
!    There are no other elastic constants 
!    dependent from those read in both
!    b-unique and c-unique cases.

   CASE(2)
!
!    triclinic case or generic 
!
!    There are no other elastic constants 
!    dependent from those read.

END SELECT 

DO i=1, 6
   DO j=i+1, 6
      el_cons_t(j,i,:)=el_cons_t(i,j,:)
   END DO
END DO

RETURN
END SUBROUTINE expand_el_cons

SUBROUTINE write_macro_el_on_file(temp, ntemp, macro_el_t, filename)
!
! This routine creates a file with macro-elasticity variables as a function 
! of temperature.
!

USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), macro_el_t(8,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

REAL(DP) :: macro_el_t_aver(4,ntemp) !Reuss-Voigt-Hill variables
INTEGER :: itemp, iu_macro_el, ios
INTEGER :: find_free_unit

iu_macro_el=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_macro_el, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_macro_el_on_file','opening macro elasticity file',&
                                                             ABS(ios))

IF (meta_ionode) THEN
   WRITE(iu_macro_el,'("#",2x,"b0: bulk modulus (kbar), e0: Young modulus (kbar), g0: & 
                                  shear modulus (kbar), nu: Poisson ratio")')
   WRITE(iu_macro_el,'("#",2x,"v: Voigt average, r: Reuss average")')
   WRITE(iu_macro_el,'("#",2x,"T  ", 20x, "b0v ", 13x, "e0v ", 13x, "g0v ",  &
               & 13x, "  nuv", 13x, "   b0r", 13x, "e0r ", 13x, "g0r ", 13x, "nur ")')
   DO itemp=2,ntemp-1
      
      macro_el_t_aver(1,itemp) = (macro_el_t(1,itemp)+macro_el_t(5,itemp))*0.5_DP
      macro_el_t_aver(2,itemp) = (macro_el_t(2,itemp)+macro_el_t(6,itemp))*0.5_DP
      macro_el_t_aver(3,itemp) = (macro_el_t(3,itemp)+macro_el_t(7,itemp))*0.5_DP
      macro_el_t_aver(4,itemp) = (macro_el_t(2,itemp)+macro_el_t(6,itemp))/ &
                          (2.d0*(macro_el_t(3,itemp)+macro_el_t(7,itemp)))-1.0_DP

      WRITE(iu_macro_el,'(e16.8, 8e18.10)') temp(itemp), macro_el_t(1,itemp), &
              macro_el_t(2,itemp), macro_el_t(3,itemp), macro_el_t(4,itemp), &
              macro_el_t(5,itemp), macro_el_t(6,itemp), macro_el_t(7,itemp), &
              macro_el_t(8,itemp)
   ENDDO

   CLOSE(iu_macro_el)

END IF

iu_macro_el=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_macro_el, FILE=TRIM(filename)//'_aver', FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=40, IOSTAT=ios)
40 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_macro_el_on_file','opening macro elasticity file',&
                                                             ABS(ios))

IF (meta_ionode) THEN
   WRITE(iu_macro_el,'("#",2x,"Macro-elasticity variable within the Voigt-Reuss-Hill &
                                                            approximation")')
   WRITE(iu_macro_el,'("#",2x,"b0: bulk modulus (kbar), e0: Young modulus (kbar), g0: & 
                                  shear modulus (kbar), nu: Poisson ratio")')
   WRITE(iu_macro_el,'("#",2x,"v: Voigt average, r: Reuss average")')
   WRITE(iu_macro_el,'("#",2x,"T  ", 20x, "b0 ", 13x, "e0 ", 13x, "g0 ",  &
                                                             & 13x, "  nu")')
   DO itemp=2,ntemp-1

      WRITE(iu_macro_el,'(e16.8, 8e18.10)') temp(itemp), macro_el_t_aver(1,itemp), &
              macro_el_t_aver(2,itemp), macro_el_t_aver(3,itemp), macro_el_t_aver(4,itemp)
   ENDDO

   CLOSE(iu_macro_el)

ENDIF

RETURN
END SUBROUTINE write_macro_el_on_file

SUBROUTINE write_sound_on_file(temp, ntemp, v_t, filename)
!
! This routine creates a file with sound velocities as a function 
! of temperature.
!

USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), v_t(3,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_sound, ios
INTEGER :: find_free_unit

iu_sound=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_sound, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_sound_on_file','opening sound velocities file',&
                                                             ABS(ios))

IF (meta_ionode) THEN
   WRITE(iu_sound,'("#",2x,"V_P: compressional velocity (m/s), V_B: bulk &
                      velocity (m/s), V_G: shear velocity (m/s)")')
   WRITE(iu_sound,'("#",2x,"T  ", 20x, "V_P ", 13x, "V_B ", 13x, "V_G ")')
   DO itemp=2,ntemp-1
      WRITE(iu_sound,'(e16.8, 8e18.10)') temp(itemp), v_t(1,itemp), &
                                                v_t(2,itemp), v_t(3,itemp)
   ENDDO

   CLOSE(iu_sound)
ENDIF

RETURN
END SUBROUTINE write_sound_on_file

END MODULE elastic_constants
