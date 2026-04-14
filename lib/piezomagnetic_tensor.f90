!
! Copyright (C) 2026 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE piezomagnetic_tensor
!---------------------------------------------------------------------------
!
!   this module contains the support routines for the calculation
!   of the piezomagnetic tensor. It should not be necessary since
!   the modulus for piezoelectricity should be sufficient, but
!   the piezoelectric module needs to be cleaned
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL(DP) :: piezom_tensor(3,6)   ! The piezomagnetic tensor h_{alpha,m}

  REAL(DP), ALLOCATABLE :: mag_strain(:,:) ! The magnetization for each strain
                                    ! in units of Bohr magnetons
!
!   Some array to simplify dealing with piezomagnetic tensor
!
  INTEGER, PARAMETER :: pm_elements=18

  CHARACTER(LEN=6) :: pm_names(pm_elements)

  DATA  pm_names /                                                   &
         'h_{11}', 'h_{12}', 'h_{13}', 'h_{14}', 'h_{15}', 'h_{16}', &
         'h_{21}', 'h_{22}', 'h_{23}', 'h_{24}', 'h_{25}', 'h_{26}', &
         'h_{31}', 'h_{32}', 'h_{33}', 'h_{34}', 'h_{35}', 'h_{36}' /


  INTEGER, PARAMETER :: pm_types=24

  INTEGER :: pm_code_group(pm_types)  ! code of the point group for each type
                                      ! This tensor requires magnetic point
                                      ! groups and in input one should give
                                      ! the b_birss code group since the
                                      ! piezomagnetic tensor is odd with 
                                      ! respect to time reversal and of odd
                                      ! rank
                                      ! Please check the magnetic_point_group
                                      ! module.
  DATA  pm_code_group / 30, 28, 26, 24, 21, 21, 17, 15, 14, 13, 13, 12, &
                         11, 10, 9, 8, 7, 6, 5, 4, 4, 3, 3, 1 /

  INTEGER  :: pm_present(pm_elements, pm_types)

  DATA pm_present / &
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 30  T_d
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 28  T
       0,0,0,2,3,0, 0,0,0,0,0,0, 1,0,0,0,0,4, & ! 26  S_4
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,2, & ! 24  D_2d
       1,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 21  D_3h 106 m perp x2
       0,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,0,0, & ! 21  D_3h 107 m perp x1
       1,0,0,0,0,0, 0,2,0,0,0,0, 0,0,0,0,0,0, & ! 17  C_3h
       0,0,0,0,3,0, 0,0,0,0,0,0, 1,0,2,0,0,0, & ! 15  C_6v
       0,0,0,0,3,0, 0,0,0,0,0,0, 1,0,2,0,0,0, & ! 14  C_4v
       0,0,0,0,4,0, 0,1,0,0,0,0, 2,0,3,0,0,0, & ! 13  C_3v 72  m perp x1
       1,0,0,0,4,0, 0,0,0,0,0,0, 2,0,3,0,0,0, & ! 13  C_3v 73  m perp x2
       0,0,0,0,4,0, 0,0,0,5,0,0, 1,2,3,0,0,0, & ! 12  C_2v
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 11  D_6
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 10  D_4
       1,0,0,2,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & !  9  D_3
       0,0,0,1,0,0, 0,0,0,0,2,0, 0,0,0,0,0,3, & !  8  D_2
       0,0,0,3,4,0, 0,0,0,0,0,0, 1,0,2,0,0,0, & !  7  C_6
       0,0,0,3,4,0, 0,0,0,0,0,0, 1,0,2,0,0,0, & !  6  C_4
       1,0,0,5,6,0, 0,2,0,0,0,0, 3,0,4,0,0,0, & !  5  C_3
       0,0,0,4,0,5, 1,2,3,0,6,0, 0,0,0,7,0,8, & !  4  C_2 b-unique
       0,0,0,4,5,0, 0,0,0,6,7,0, 1,2,3,0,0,8, & !  4  C_2 c-unique
       1,2,3,0,4,0, 0,0,0,5,0,6, 7,8,9,0,10,0, & !  3  C_s b-unique
       1,2,3,0,0,4, 5,6,7,0,0,8, 0,0,0,9,10,0, & !  3  C_s c-unique
       1,2,3,4,5,6, 7,8,9,10,11,12, 13,14,15,16,17,18 / !  1  C_1

!
!  Magnetic point groups corresponding to T_d and T
!  T, T_d, 4'32', -4'3m', m-3m'
!
!  Magnetic point groups corresponding to S_4
!  4', \bar 4', 4'/m
!
!  Magnetic point groups corresponding to D_2d
!  4'22, 4'mm', -4'm2, -4'2m', 4'/mmm'
!
!  Magnetic point groups corresponding to D_3h
!  6'22', 6'mm', -6'2m', -6'2'm, 6'/m'mm'
!
!  Magnetic point groups corresponding to C_3h
!  6', -6', 6'/m'
!
!  Magnetic point groups corresponding to C_6v, C_4v
!  42'2', 4m'm', -42'm', 4/mm'm', 62'2', 6m'm', -6m'2', 6/mm'm
!
!  Magnetic point groups corresponding to C_3v
!  32', 3m', -3m'
!
!  Magnetic point groups corresponding to C_2v
!  2'2'2, m'm2', m'm'2, mm'm'
!
!  Magnetic point groups corresponding to D_6, D_4
!  622, 6mm, -62m, 6/mmm, 422, 4mm, -42m, 4/mmm
!
!  Magnetic point groups corresponding to D_3
!  32, 3m, -3m
!
!  Magnetic point groups corresponding to D_2
!  222, mm2, mmm
!
!  Magnetic point groups corresponding to C_6, C_4
!  4, -4, 4/m, 6, -6, 6/m, 
!
!  Magnetic point groups corresponding to C_3
!  3, -3
!
!  Magnetic point groups corresponding to C_2
!  2, m, 2/m
!
!  Magnetic point groups corresponding to C_s
!  m', 2', 2'/m'
!
!  Magnetic point groups corresponding to C_1
!  1, -1
!

  PUBLIC compute_piezom_tensor, mag_strain,               &
         print_piezom_tensor, piezom_tensor,              &
         pm_elements, clean_piezom_tensor,                &
         print_piezom_info, allocate_piezom,              &
         compute_magnetization_equil,                     &
         deallocate_piezom, write_piezom_tensor,          &
         read_piezom_tensor, write_piezom_tensor_on_file, &
         read_piezom_tensor_from_file,                    &
         pm_names, pm_types, pm_present,                  &
         pm_code_group, get_pm_type

CONTAINS
!
!---------------------------------------------------------------------------
SUBROUTINE print_piezom_tensor(piezom_tensor, fact, label)
!---------------------------------------------------------------------------
!
!  This routine writes on output the piezoelectric tensor
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label
REAL(DP), INTENT(IN) :: piezom_tensor(3,6)
REAL(DP) :: fact
INTEGER :: i, j

WRITE(stdout,'(5x, a)') TRIM(label)
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
!
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (piezom_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_piezom_tensor

!-------------------------------------------------------------------------
SUBROUTINE allocate_piezom(nwork)
!-------------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: nwork

ALLOCATE(mag_strain(3,nwork))
mag_strain=0.0_DP
piezom_tensor=0.0_DP

RETURN
END SUBROUTINE allocate_piezom
!
!-------------------------------------------------------------------------
SUBROUTINE deallocate_piezom()
!-------------------------------------------------------------------------
IMPLICIT NONE

IF (ALLOCATED(mag_strain)) DEALLOCATE(mag_strain)

RETURN
END SUBROUTINE deallocate_piezom
!
!-------------------------------------------------------------------------
SUBROUTINE write_piezom_tensor(filename,mag0)
!-------------------------------------------------------------------------
!
!  This routine writes the piezomagnetic tensor on file.
!  It must be called after computing the piezomagnetic tensor.
!  It saves: 
!  the spontaneous magnetization of the unperturbed geometry
!  the piezomagnetic tensor
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: mag0(3)
REAL(DP) :: fact
INTEGER :: find_free_unit
INTEGER :: outunit, ios, i, j

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', &
        FORM='formatted', ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_piezom_tensor','ploblem opening output file', ABS(ios))

IF (ionode) THEN
   WRITE(outunit,'("Spontaneous magnetization (bohr magnetons)")')
   WRITE(outunit,'(3e19.10)') (mag0(i), i=1,3)
   WRITE(outunit,*)
   WRITE(outunit,'("The piezomagnetic tensor (bohr magnetons)")')
   DO i=1,3
      WRITE(outunit,'(6e19.10)') (piezom_tensor(i,j), j=1,6)
   ENDDO
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_piezom_tensor

!-------------------------------------------------------------------------
SUBROUTINE read_piezom_tensor(filename, mag0, exists)
!-------------------------------------------------------------------------
!
!  This routine reads the piezomagnetic tensor from file.
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP) :: mag0(3)
LOGICAL, INTENT(OUT) :: exists
INTEGER :: inunit, ios, i, j
INTEGER :: find_free_unit

IF (ionode) THEN
   inunit=find_free_unit()
   OPEN(UNIT=inunit, FILE=TRIM(filename), STATUS='old', FORM='formatted', &
       ERR=100, IOSTAT=ios)
ENDIF

IF (ionode) THEN
   READ(inunit,*)
   READ(inunit,'(3e19.10)') (mag0(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (piezom_tensor(i,j), j=1,6)
   ENDDO
   CLOSE(inunit)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
IF (ios /= 0) THEN
   exists=.FALSE.
   RETURN
ENDIF
CALL mp_bcast(mag0,ionode_id,intra_image_comm)
CALL mp_bcast(piezom_tensor,ionode_id,intra_image_comm)
exists=.TRUE.

RETURN
END SUBROUTINE read_piezom_tensor
!
!---------------------------------------------------------------------------
SUBROUTINE compute_piezom_tensor(mag_geo, epsil_geo, nwork, &
                                ngeo, ibrav, code_group, code_group_ext)
!---------------------------------------------------------------------------
!
!  This routine computes the piezomagnetic tensor h_{\alpha,m} by fitting the 
!  magnetization strain relation with a second order polynomial. This is 
!  calculated on the basis of the B point group as defined by Birss
!  and obtained from the magnetic point group of the solid.
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: mag_geo(3,nwork), epsil_geo(3,3,nwork)
INTEGER, INTENT(IN) :: ngeo, ibrav, code_group, code_group_ext, nwork
INTEGER :: i, j, igeo, alpha, ind, mn
LOGICAL :: check_group_ibrav

WRITE(stdout,'(/,20x,40("-"),/)')
piezom_tensor=0.0_DP
IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group) 
     CASE(2,16,18,19,20,22,23,25,27,29,32) 
     CASE(3)
!
!  C_s   Monoclinic
!  Magnetic point groups corresponding to C_s
!  m', 2', 2'/m'
!
!        WRITE(stdout,'(5x,"( h11  h12  h13   .   h15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .   h24   .   h26 )")') 
!        WRITE(stdout,'(5x,"( h31  h32  h33   .   h35   .  )")') 
!
        CALL piezom_ij(1, 1, ngeo, epsil_geo, mag_geo )
        CALL piezom_ij(1, 2, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1) )
        CALL piezom_ij(1, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                        mag_geo(1,2*ngeo+1) )
        CALL piezom_ij(2, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                                        mag_geo(1,5*ngeo+1) )
        CALL piezom_ij(3, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                                        mag_geo(1,4*ngeo+1) )


        IF (ibrav==-12.OR.ibrav==-13) THEN

           CALL piezom_ij(3, 1, ngeo, epsil_geo, mag_geo)
           CALL piezom_ij(3, 2, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
           CALL piezom_ij(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                     mag_geo(1,2*ngeo+1))
           CALL piezom_ij(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                     mag_geo(1,3*ngeo+1))
           CALL piezom_ij(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                     mag_geo(1,4*ngeo+1))

        ELSE
!                WRITE(stdout,'(5x,"( h11  h12  h13   .    .   h16 )")') 
!                WRITE(stdout,'(5x,"( h21  h22  h23   .    .   h26 )")') 
!                WRITE(stdout,'(5x,"(  .    .    .   h34  h35   .  )")') 
!
           CALL piezom_ij(2, 1, ngeo, epsil_geo, mag_geo)
           CALL piezom_ij(2, 2, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
           CALL piezom_ij(2, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                            mag_geo(1,2*ngeo+1))
           CALL piezom_ij(3, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                            mag_geo(1,3*ngeo+1))
           CALL piezom_ij(1, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                            mag_geo(1,5*ngeo+1))

        ENDIF
     CASE(4)
!
!  C_2   Monoclinic
!  Magnetic point groups corresponding to C_2
!  2, m, 2/m
!
        CALL piezom_ij(1, 4, ngeo, epsil_geo(1,1,3*ngeo+1),mag_geo(1,3*ngeo+1))
        CALL piezom_ij(2, 5, ngeo, epsil_geo(1,1,4*ngeo+1),mag_geo(1,4*ngeo+1))
        CALL piezom_ij(3, 6, ngeo, epsil_geo(1,1,5*ngeo+1),mag_geo(1,5*ngeo+1))

        IF (ibrav==-12) THEN
!            WRITE(stdout,'(5x,"(  .    .    .   h14   .   h16 )")') 
!            WRITE(stdout,'(5x,"( h21  h22  h23   .   h25   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   h34   .   h36 )")') 

           CALL piezom_ij(2, 1, ngeo, epsil_geo, mag_geo)
           CALL piezom_ij(2, 2, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
           CALL piezom_ij(2, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                     mag_geo(1,2*ngeo+1))
           CALL piezom_ij(3, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                                     mag_geo(1,3*ngeo+1))
           CALL piezom_ij(1, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                                     mag_geo(1,5*ngeo+1))
        ELSE
!            WRITE(stdout,'(5x,"(  .    .    .   h14  h15   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   h24  h25   .  )")') 
!            WRITE(stdout,'(5x,"( h31  h32  h33   .    .   h36 )")') 
!
            CALL piezom_ij(3, 1, ngeo, epsil_geo, mag_geo)
            CALL piezom_ij(3, 2, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
            CALL piezom_ij(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                      mag_geo(1,2*ngeo+1))
            CALL piezom_ij(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                      mag_geo(1,3*ngeo+1))
            CALL piezom_ij(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                      mag_geo(1,4*ngeo+1))

         ENDIF
      CASE(5)
!
!   C_3
!
!        WRITE(stdout,'(5x,"( g11 -g11   .   g14  g15 -g22 )")') 
!        WRITE(stdout,'(5x,"(-g22  g22   .   g15 -g14 -g11 )")') 
!        WRITE(stdout,'(5x,"( g31 -g31  g33   .    .    .  )")') 
!
         CALL piezom_ij(1, 1, ngeo, epsil_geo, mag_geo)
         piezom_tensor(1,2) = -piezom_tensor(1,1)
         piezom_tensor(2,6) = -piezom_tensor(1,1)
         CALL piezom_ij(2, 1, ngeo, epsil_geo, mag_geo)
         piezom_tensor(2,2) = -piezom_tensor(2,1)
         piezom_tensor(1,6) =  piezom_tensor(2,1)
         CALL piezom_ij(3, 1, ngeo, epsil_geo, mag_geo)
         piezom_tensor(3,2) = -piezom_tensor(3,1)
         CALL piezom_ij(3, 3, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
         CALL piezom_ij(1, 4, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                mag_geo(1,2*ngeo+1))
         piezom_tensor(2,5) = -piezom_tensor(1,4)
         CALL piezom_ij(2, 4, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                mag_geo(1,2*ngeo+1))
         piezom_tensor(1,5) = -piezom_tensor(2,4)
      CASE(6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!  Magnetic point groups corresponding to C_6, C_4
!  4, -4, 4/m, 6, -6, 6/m, 
!
!             WRITE(stdout,'(5x,"(  .    .    .   h14  h15   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .   h15 -h14   .  )")') 
!             WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 
!
        CALL piezom_ij(3, 1, ngeo, epsil_geo,mag_geo)
        piezom_tensor(3,2) = piezom_tensor(3,1)
        CALL piezom_ij(3, 3, ngeo, epsil_geo(1,1,ngeo+1),mag_geo(1,ngeo+1))
        CALL piezom_ij(1, 4, ngeo, epsil_geo(1,1,2*ngeo+1),mag_geo(1,2*ngeo+1))
        piezom_tensor(2,5) = -piezom_tensor(1,4)
        CALL piezom_ij(2, 4, ngeo, epsil_geo(1,1,2*ngeo+1),mag_geo(1,2*ngeo+1))
        piezom_tensor(1,5) = piezom_tensor(2,4)

     CASE(8)
!
!  D_2 (222) Orthorombic
!  Magnetic point groups corresponding to D_2
!  222, mm2, mmm
!
!         WRITE(stdout,'(5x,"(  .    .    .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   h25   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   h36 )")') 

        CALL piezom_ij(1, 4, ngeo, epsil_geo, mag_geo)
        CALL piezom_ij(2, 5, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
        CALL piezom_ij(3, 6, ngeo, epsil_geo(1,1,2*ngeo+1),mag_geo(1,2*ngeo+1))

      CASE(9)
!
! D_3  Trigonal 
! Magnetic point groups corresponding to D_3
! 32, 3m, -3m
!
!         WRITE(stdout,'(5x,"( h11 -h11   .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -h14 -h11 )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezom_ij(1, 1, ngeo, epsil_geo, mag_geo)
        piezom_tensor(1,2) = -piezom_tensor(1,1)
        piezom_tensor(2,6) = -piezom_tensor(1,1)
        CALL piezom_ij(1, 4, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
        piezom_tensor(2,5) = -piezom_tensor(1,4)

     CASE(10,11)
!
! D_4  tetragonal, D_6 hexagonal
! Magnetic point groups corresponding to D_6, D_4
! 622, 6mm, -62m, 6/mmm, 422, 4mm, -42m, 4/mmm
!
!         WRITE(stdout,'(/,5x,"(  .    .    .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   h14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezom_ij(1, 4, ngeo, epsil_geo, mag_geo)
        piezom_tensor(2,5) = -piezom_tensor(1,4)

     CASE(12)
!
! C_2v  Orthorombic
! Magnetic point groups corresponding to C_2v
! 2'2'2, m'm2', m'm'2, mm'm'
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   h15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   h24   .    .  )")') 
!         WRITE(stdout,'(5x,"( h31  h32  h33   .    .    .  )")') 

        CALL piezom_ij(3, 1, ngeo, epsil_geo, mag_geo)
        CALL piezom_ij(3, 2, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
        CALL piezom_ij(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1),mag_geo(1,2*ngeo+1))
        CALL piezom_ij(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1),mag_geo(1,3*ngeo+1))
        CALL piezom_ij(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1),mag_geo(1,4*ngeo+1))

     CASE(13)
!
! C_3v  Trigonal. Assuming m perpendicular to x1
!  Magnetic point groups corresponding to C_3v
!  32', 3m', -3m'
!
        IF (code_group_ext==72) THEN
!
!         m perpendicular to x1
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   h15  h21 )")') 
!         WRITE(stdout,'(5x,"( h21 -h21   .   h15   .    .  )")') 
!         WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 
!
           CALL piezom_ij(2, 1, ngeo, epsil_geo, mag_geo)
           piezom_tensor(2,2) = -piezom_tensor(2,1)
           piezom_tensor(1,6) = piezom_tensor(2,1)
        ELSEIF (code_group_ext==73) THEN
!
!         m perpendicular to x2
!
!         WRITE(stdout,'(5x,"( h11 -h11   .    .   h15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   h15   .  -h11 )")') 
!         WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 
!
           CALL piezom_ij(1, 1, ngeo, epsil_geo, mag_geo)
           piezom_tensor(1,2) = -piezom_tensor(1,1)
           piezom_tensor(2,6) = -piezom_tensor(1,1)
        ENDIF

        CALL piezom_ij(3, 1, ngeo, epsil_geo, mag_geo)
        piezom_tensor(3,2) = piezom_tensor(3,1)
        CALL piezom_ij(3, 3, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
        CALL piezom_ij(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                        mag_geo(1,2*ngeo+1))
        piezom_tensor(2,4) = piezom_tensor(1,5)

     CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
! Magnetic point groups corresponding to C_6v, C_4v
! 42'2', 4m'm', -42'm', 4/mm'm', 62'2', 6m'm', -6m'2', 6/mm'm
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   h15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   h15   .    .  )")') 
!         WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 

        CALL piezom_ij(3, 1, ngeo, epsil_geo(1,1,1), mag_geo(1,1))
!
!   The first strain is of the form (epsilon,epsilon,0,0,0,0) to keep the
!   cell tetragonal or hexagonal. Hence the factor 0.5.
!
        piezom_tensor(3,1) = piezom_tensor(3,1) * 0.5_DP
        piezom_tensor(3,2) = piezom_tensor(3,1)
        CALL piezom_ij(3, 3, ngeo, epsil_geo(1,1,ngeo+1),mag_geo(1,ngeo+1))
        CALL piezom_ij(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1),mag_geo(1,2*ngeo+1))
        piezom_tensor(2,4) = piezom_tensor(1,5)

     CASE(17)
!
! C_3h hexagonal
! Magnetic point groups corresponding to C_3h
! 6', -6', 6'/m'
!
!             WRITE(stdout,'(5x,"( h11 -h11   .    .    .   h21 )")') 
!             WRITE(stdout,'(5x,"( h21 -h21   .    .    .  -h11 )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezom_ij(1, 1, ngeo, epsil_geo, mag_geo)
        piezom_tensor(1,2) = -piezom_tensor(1,1)
        piezom_tensor(2,6) = -piezom_tensor(1,1)
        CALL piezom_ij(2, 1, ngeo, epsil_geo, mag_geo)
        piezom_tensor(2,2) = -piezom_tensor(2,1)
        piezom_tensor(1,6) =  piezom_tensor(2,1)

      CASE(21)
!
! D_3h hexagonal
! Magnetic point groups corresponding to D_3h
! 6'22', 6'mm', -6'2m', -6'2'm, 6'/m'mm'
!
        IF (code_group_ext==107) THEN
!
!      mirror perpendicular to x1
!
!             WRITE(stdout,'(5x,"(  .    .    .    .    .   h21 )")') 
!             WRITE(stdout,'(5x,"( h21 -h21   .    .    .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
!
           CALL piezom_ij(2, 1, ngeo, epsil_geo, mag_geo)
           piezom_tensor(2,2) = -piezom_tensor(2,1)
           piezom_tensor(1,6) =  piezom_tensor(1,2)
        ELSEIF (code_group_ext==106) THEN
!
!      mirror perpendicular to x2
!     
!             WRITE(stdout,'(5x,"( g11 -g11   .    .    .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .  -g11 )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
!
           CALL piezom_ij(1, 1, ngeo, epsil_geo, mag_geo)
           piezom_tensor(1,2) = -piezom_tensor(1,1)
           piezom_tensor(2,6) = -piezom_tensor(1,1)
        ENDIF


     CASE(24)
!
! D_2d tetragonal: axis 2 || x1
!
!         WRITE(stdout,'(5x,"(  .    .    .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   h14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   h36 )")') 

        CALL piezom_ij(1, 4, ngeo, epsil_geo, mag_geo)
        piezom_tensor(2,5) = piezom_tensor(1,4)
        CALL piezom_ij(3, 6, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))

     CASE(26)
!
! S_4 tetragonal
! Magnetic point groups corresponding to S_4
! 4', \bar 4', 4'/m
!
!        WRITE(stdout,'(5x,"(  .    .    .   h14  h15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .  -h15  h14   .  )")') 
!        WRITE(stdout,'(5x,"( h31 -h31   .    .    .   h36 )")') 

        CALL piezom_ij(3, 1, ngeo, epsil_geo, mag_geo)
        piezom_tensor(3,2) = -piezom_tensor(3,1)
        CALL piezom_ij(1, 4, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
        piezom_tensor(2,5) = piezom_tensor(1,4)
        CALL piezom_ij(2, 4, ngeo, epsil_geo(1,1,ngeo+1), mag_geo(1,ngeo+1))
        piezom_tensor(1,5) = -piezom_tensor(2,4)
        CALL piezom_ij(3, 6, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                        mag_geo(1,2*ngeo+1))

     CASE(28,30)
!
! T, T_d cubic
!  Magnetic point groups corresponding to T_d and T
!  T, T_d, 4'32', -4'3m', m-3m'
!
!             WRITE(stdout,'(5x,"(  .    .    .   h14   .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .   h14   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .   h14 )")') 
!
        CALL piezom_ij(1, 4, ngeo, epsil_geo, mag_geo)
        piezom_tensor(2,5) = piezom_tensor(1,4)
        piezom_tensor(3,6) = piezom_tensor(1,4)
     CASE(31)
!
!  O All components vanish. We do nothing.
!
     CASE DEFAULT
!
!  C_1 
!
        DO mn=1,6
           ind = ngeo * (mn-1)
           DO alpha=1,3
              CALL piezom_ij(alpha, mn, ngeo, epsil_geo(1,1,ind+1), &
                                             mag_geo(1,ind+1))
           ENDDO
        ENDDO
  END SELECT
ELSE
   DO mn=1,6
      ind = ngeo * (mn-1)
      DO alpha=1,3
         CALL piezom_ij(alpha, mn, ngeo, epsil_geo(1,1,ind+1), &
                                        mag_geo(1,ind+1))
      ENDDO
   ENDDO
ENDIF
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE compute_piezom_tensor
!
!---------------------------------------------------------------------------
SUBROUTINE clean_piezom_tensor(piezo, ibrav, code_group, code_group_ext)
!---------------------------------------------------------------------------
!
! This routine receives a piezoelectric tensor and a point group
! and sets to zero all components that have not been calculated.
! The point group should be the b_birss point group that correspond
! to a given magnetic point group.
!
!
IMPLICIT NONE
REAL(DP), INTENT(INOUT) :: piezo(3,6)
INTEGER, INTENT(IN) ::  ibrav, code_group, code_group_ext
INTEGER :: i, j, igeo, alpha, ind, mn
LOGICAL :: check_group_ibrav

REAL(DP) :: piezom(3,6)

piezom=0.0_DP
IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group) 
     CASE(2,16,18,19,20,22,23,25,27,29,32) 
     CASE(3)
!
!  C_s   Monoclinic
!
        piezom(1,1)=piezo(1,1)
        piezom(1,2)=piezo(1,2)
        piezom(1,3)=piezo(1,3)
        piezom(2,6)=piezo(2,6)
        piezom(3,5)=piezo(3,5)

        IF (ibrav==-12.OR.ibrav==-13) THEN
!
!        WRITE(stdout,'(5x,"( h11  h12  h13   .   h15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .   h24   .   h26 )")') 
!        WRITE(stdout,'(5x,"( h31  h32  h33   .   h35   .  )")') 
!
           piezom(3,1)=piezo(3,1)
           piezom(3,2)=piezo(3,2)
           piezom(3,3)=piezo(3,3)
           piezom(2,4)=piezo(2,4)
           piezom(1,5)=piezo(1,5)

        ELSE
!                WRITE(stdout,'(5x,"( h11  h12  h13   .    .   h16 )")') 
!                WRITE(stdout,'(5x,"( h21  h22  h23   .    .   h26 )")') 
!                WRITE(stdout,'(5x,"(  .    .    .   h34  h35   .  )")') 
!
           piezom(2,1)=piezo(2,1)
           piezom(2,2)=piezo(2,2)
           piezom(2,3)=piezo(2,3)
           piezom(1,6)=piezo(1,6)
           piezom(3,4)=piezo(3,4)

        ENDIF
     CASE(4)
!
!  C_2   Monoclinic
!
        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)
        piezom(3,6)=piezo(3,6)

        IF (ibrav==-12) THEN
!            WRITE(stdout,'(5x,"(  .    .    .   h14   .   h16 )")') 
!            WRITE(stdout,'(5x,"( h21  h22  h23   .   h25   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   h34   .   h36 )")') 

           piezom(2,1)=piezo(2,1)
           piezom(2,2)=piezo(2,2)
           piezom(2,3)=piezo(2,3)
           piezom(1,6)=piezo(1,6)
           piezom(3,4)=piezo(3,4)
        ELSE
!
!            WRITE(stdout,'(5x,"(  .    .    .   h14  h15   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   h24  h25   .  )")') 
!            WRITE(stdout,'(5x,"( h31  h32  h33   .    .   h36 )")') 
!
           piezom(3,1)=piezo(3,1)
           piezom(3,2)=piezo(3,2)
           piezom(3,3)=piezo(3,3)
           piezom(1,5)=piezo(1,5)
           piezom(2,4)=piezo(2,4)
         ENDIF
     CASE(5)
!
!   C_3
!
!            WRITE(stdout,'(5x,"( h11 -h11   .   h14  h15 -h22 )")') 
!            WRITE(stdout,'(5x,"(-h22  h22   .   h15 -h14 -h11 )")') 
!            WRITE(stdout,'(5x,"( h31 -h31  h33   .    .    .  )")') 
!
           piezom(1,1)=piezo(1,1)
           piezom(1,2)=piezo(1,2)
           piezom(1,4)=piezo(1,4)
           piezom(1,5)=piezo(1,5)
           piezom(1,6)=piezo(1,6)
           piezom(2,1)=piezo(2,1)
           piezom(2,2)=piezo(2,2)
           piezom(2,4)=piezo(2,4)
           piezom(2,5)=piezo(2,5)
           piezom(2,6)=piezo(2,6)
           piezom(3,1)=piezo(3,1)
           piezom(3,2)=piezo(3,2)
           piezom(3,3)=piezo(3,3)

      CASE(6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!
!             WRITE(stdout,'(5x,"(  .    .    .   h14  h15   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .   h24 -h14   .  )")') 
!             WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 
!
        piezom(3,1)=piezo(3,1)
        piezom(3,2)=piezo(3,2)
        piezom(3,3)=piezo(3,3)
        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)
        piezom(2,4)=piezo(2,4)
        piezom(1,5)=piezo(1,5)

     CASE(8)
!
!  D_2 (222) Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   h25   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   h36 )")') 

        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)
        piezom(3,6)=piezo(3,6)

      CASE(9)
!
! D_3  Trigonal 
!
!         WRITE(stdout,'(5x,"( h11 -h11   .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -h14 -h11 )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        piezom(1,1)=piezo(1,1)
        piezom(1,2)=piezo(1,2)
        piezom(2,6)=piezo(2,6)
        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)

     CASE(10,11)
!
! D_4  tetragonal, D_6 hexagonal
!
!         WRITE(stdout,'(/,5x,"(  .    .    .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -h14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)

     CASE(12)
!
! C_2v  Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   h15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   h24   .    .  )")') 
!         WRITE(stdout,'(5x,"( h31  h32  h33   .    .    .  )")') 

        piezom(3,1)=piezo(3,1)
        piezom(3,2)=piezo(3,2)
        piezom(3,3)=piezo(3,3)
        piezom(2,4)=piezo(2,4)
        piezom(1,5)=piezo(1,5)

     CASE(13)
!
! C_3v  Trigonal. 
!
       IF (code_group_ext==72) THEN
!
!         m perpendicular to x1
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   h15 -h21 )")') 
!         WRITE(stdout,'(5x,"( h21 -h21   .   h15   .    .  )")') 
!         WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 
!
           piezom(2,1)=piezo(2,1)
           piezom(2,2)=piezo(2,2)
           piezom(1,6)=piezo(1,6)
       ELSEIF (code_group_ext==73) THEN
!
!         m perpendicular to x2
!
!          WRITE(stdout,'(5x,"( h11 -h11   .    .   h15   .  )")') 
!          WRITE(stdout,'(5x,"(  .    .    .   h15   .  -e11 )")') 
!          WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 
           piezom(1,1)=piezo(1,1)
           piezom(1,2)=piezo(1,2)
           piezom(2,6)=piezo(2,6)
        ENDIF
        piezom(3,1)=piezo(3,1)
        piezom(3,2)=piezo(3,2)
        piezom(3,3)=piezo(3,3)
        piezom(1,5)=piezo(1,5)
        piezom(2,4)=piezo(2,4)

     CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   h15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   h15   .    .  )")') 
!         WRITE(stdout,'(5x,"( h31  h31  h33   .    .    .  )")') 

        piezom(3,1)=piezo(3,1)
        piezom(3,2)=piezo(3,2)
        piezom(3,3)=piezo(3,3)
        piezom(1,5)=piezo(1,5)
        piezom(2,4)=piezo(2,4)

     CASE(17)
!
! C_3h hexagonal
!
!             WRITE(stdout,'(5x,"( h11 -h11   .    .    .  -h21 )")') 
!             WRITE(stdout,'(5x,"( h21 -h21   .    .    .   h11 )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        piezom(1,1)=piezo(1,1)
        piezom(1,2)=piezo(1,2)
        piezom(2,6)=piezo(2,6)
        piezom(2,1)=piezo(2,1)
        piezom(2,2)=piezo(2,2)
        piezom(1,6)=piezo(1,6)

      CASE(21)
!
! D_3h hexagonal
!
         IF (code_group_ext==107) THEN
!
!             WRITE(stdout,'(5x,"(  .    .    .    .    .  -h21 )")') 
!             WRITE(stdout,'(5x,"( h21 -h21   .    .    .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
            piezom(2,1)=piezo(2,1)
            piezom(2,2)=piezo(2,2)
            piezom(1,6)=piezo(1,6)
         ELSEIF (code_group_ext==106) THEN
!
!             WRITE(stdout,'(5x,"( h11 -h11   .    .    .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .  -h11 )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
!
            piezom(1,1)=piezo(1,1)
            piezom(1,2)=piezo(1,2)
            piezom(2,6)=piezo(2,6)
         ENDIF

     CASE(24)
!
! D_2d tetragonal: axis 2 || x1
!
!         WRITE(stdout,'(5x,"(  .    .    .   h14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   h14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   h36 )")') 

        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)
        piezom(3,6)=piezo(3,6)

     CASE(26)
!
! S_4 tetragonal
!
!        WRITE(stdout,'(5x,"(  .    .    .   h14  h15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .  -h15  h14   .  )")') 
!        WRITE(stdout,'(5x,"( h31 -h31   .    .    .   h36 )")') 

        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)
        piezom(3,1)=piezo(3,1)
        piezom(3,2)=piezo(3,2)
        piezom(1,5)=piezo(1,5)
        piezom(2,4)=piezo(2,4)
        piezom(3,6)=piezo(3,6)

     CASE(28,30)
!
! T, T_d cubic
!
!             WRITE(stdout,'(5x,"(  .    .    .   h14   .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .   h14   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .   h14 )")') 
!
        piezom(1,4)=piezo(1,4)
        piezom(2,5)=piezo(2,5)
        piezom(3,6)=piezo(3,6)
     CASE(31)
!
!  O cubic group. In this case all elements vanish.
!
     CASE DEFAULT
!
!  C_1 
!
      piezom(:,:)=piezo(:,:)
  END SELECT
ELSE
   piezom(:,:)=piezo(:,:)
ENDIF
piezo(:,:)=piezom(:,:)

RETURN
END SUBROUTINE clean_piezom_tensor

!---------------------------------------------------------------------------
SUBROUTINE piezom_ij(ialpha, mn, ngeo, epsil_geo, mag_geo)
!---------------------------------------------------------------------------
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit
USE voigt, ONLY : voigt_extract_indices

IMPLICIT NONE
INTEGER, INTENT(IN) :: mn, ialpha, ngeo
REAL(DP), INTENT(IN) :: epsil_geo(3,3,ngeo), mag_geo(3,ngeo)
INTEGER :: igeo, m, n, mnin
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo)

mnin=mn
CALL voigt_extract_indices(m,n,mnin)
WRITE(stdout,'("Piezomagnetic_tensor(",i1,",",i1,"): &
                    &strain,    polarization")') ialpha, mn
DO igeo=1,ngeo
   x(igeo)=epsil_geo(m,n,igeo)
   y(igeo)=mag_geo(ialpha,igeo)
   WRITE(stdout,'(20x,f15.10,3x,f15.10)') x(igeo), y(igeo)
ENDDO
CALL polyfit( x, y, ngeo, alpha, m1-1 )
piezom_tensor(ialpha, mn) = alpha(2)
!
!  The piezomagnetic tensor is defined from 
!  M_i = \sum_{jk} h_{ijk} \epsilon_{jk}
!  therefore when we apply an off-diagonal strain we obtain a
!  magnetization
!  M_i = h_{ijk} \epsilon_{jk} + h_{ikj} \epsilon_{jk} 
!      = (h_{ijk}+h_{ikj}) \epsilon_{jk}
!  So dM_i / \epsilon_{jk} is twice the piezomagnetic tensor 
!  h_{ijk} and we have to divide by 2 to get h_{ijk}.
!
!  Similarly, in Voigt notation we have
!  h_{im} = h_{ijk} while \epsilon_{m} = 2 \epsilon_{jk} 
!  for the off diagonal components and \epsilon_{m} = \epsilon{jk} for
!  the diagonal ones so we have still
!  M_i= h_{im} \epsilon_{m} 
!  We have h_{im} = d P_i \over d \epsilon_{m} and for the off diagonal
!  components  d P_i \over 2 d\epsilon_{jk} and since we have computed  
!  d P_i \over d\epsilon_{jk} above we still need to divide by 2.
!
IF (m /= n) piezom_tensor(ialpha, mn) = piezom_tensor(ialpha, mn) * 0.5_DP

RETURN
END SUBROUTINE piezom_ij
!
!----------------------------------------------------------------------------
SUBROUTINE print_piezom_info(code_group,ibrav,ngeo_strain)
!----------------------------------------------------------------------------
!
!  This routine receives the code of the Birss B point group and the
!  code of the Bravais lattice and write on output how many strain are
!  necessary to calculate the piezomagnetic tensor.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, ibrav, ngeo_strain

INTEGER :: nstrain

SELECT CASE (code_group) 
   CASE (2,16,18,19,20,22,23,25,27,29,31,32) 
      nstrain=0
   CASE (5,6,7)
!
!  C_3 trigonal, C_4, tetragonal, C_6 hexagonal
!
      WRITE(stdout,'(/,5x,"It requires three strains: e1, e3, and e4")')
      nstrain=3
   CASE (8)
!
!  D_2 (222) Orthorombic
!
      WRITE(stdout,'(/,5x,"It requires three strains: e4, e5, and e6")')
      nstrain=3
   CASE (9)
!
! D_3  Trigonal 
!
      WRITE(stdout,'(/,5x,"It requires two strains: e1 and e4")')
      nstrain=2
   CASE (10,11,28,30)
!
! D_4  tetragonal, D_6 hexagonal, T, T_d cubic
!
      WRITE(stdout,'(/,5x,"It requires one strain: e4")')
      nstrain=1
   CASE (12)
!
! C_2v  Orthorombic
!
      WRITE(stdout,'(/,5x,"It requires five strains: e1, e2, e3, e4, &
                                                               &and e5 ")')
      nstrain=5
   CASE (13,14,15)
!
! C_3v  Trigonal, C_4v tetragonal, C_6v hexagonal
!
      WRITE(stdout,'(/,5x,"It requires three strain: e1, e3, and e4 ")')
      nstrain=3
   CASE (17,21)
!
! C_3h or D_3h hexagonal
!
      WRITE(stdout,'(/,5x,"It requires one strain: e1 ")')
      nstrain=1
   CASE (24)
!
! D_2d tetragonal: axis 2 || x1
!
      WRITE(stdout,'(/,5x,"It requires two strains: e4 and e6")')
      nstrain=2
   CASE (26)
!
! S_4 tetragonal
!
      WRITE(stdout,'(/,5x,"It requires three strains: e1, e4, and e6")')
      nstrain=3
   CASE DEFAULT ! CASE(1,3,4)
!
!  C_1, C_s and C_2 
!
      WRITE(stdout,'(/,5x,"It requires all six strains")')
      nstrain=6
   END SELECT

   IF(nstrain>0) THEN
     WRITE(stdout,'(5x,"for a total of",i3,&
                                  &" scf calculations")') nstrain*ngeo_strain
   ENDIF
RETURN
END SUBROUTINE print_piezom_info
!
!-----------------------------------------------------------------------------
SUBROUTINE compute_magnetization_equil(mag_geo, epsil_geo, mag0, ngeo, &
                                                               ngeo_strain)
!-----------------------------------------------------------------------------

USE constants, ONLY : electron_si, electronmass_si, h_planck_si
USE polyfit_mod, ONLY : polyfit

IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeo, ngeo_strain
REAL(DP), INTENT(IN) :: epsil_geo(ngeo), mag_geo(3,ngeo)
REAL(DP), INTENT(OUT) :: mag0(3)

INTEGER :: igeo, ipol, istep, nstep, count_igeo
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo,3)
REAL(DP) :: polar_at(3), fact

WRITE(stdout,'(/,20x,40("-"),/)')
nstep=ngeo/ngeo_strain
count_igeo=0
DO istep=1, nstep
   DO igeo=1,ngeo_strain
      count_igeo=count_igeo+1
      x(igeo)=epsil_geo(count_igeo)
      DO ipol=1,3 ! loop on direct lattice vectors
         y(igeo,ipol)=mag_geo(ipol,count_igeo)
      ENDDO
   ENDDO
   DO ipol=1,3
      WRITE(stdout,'("Strain type = ", i2," Magnetization(",i1,")")') istep, &
                                                                     ipol
      DO igeo=1,ngeo_strain
         WRITE(stdout,'(2f15.10)') x(igeo), y(igeo,ipol)
      ENDDO
      CALL polyfit( x, y(1,ipol), ngeo_strain, alpha, m1-1 )
      mag0(ipol)=alpha(1)
   ENDDO
!
!  and print it
!
   WRITE(stdout,'(/,"Magnetic moment of the equilibrium structure &
                             &strain type: ",i4)') istep
   WRITE(stdout,'(5x,"M=(", 2(f10.5,","), f10.5, "   ) mu_B units")') &
                                                             mag0(:)
   fact= h_planck_si * electron_si *0.5_DP / electronmass_si
   WRITE(stdout,'(5x,"M=(", 2(f10.5,","), f10.5, "   ) J/T ",/)') &
                                                         mag0(:) * fact
ENDDO

WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE compute_magnetization_equil
!
!-------------------------------------------------------------------------
SUBROUTINE write_piezom_tensor_on_file(temp, ntemp, ibrav, code_group, &
                         code_group_ext, piezom_tensor_t, filename, iflag)
!-------------------------------------------------------------------------
!
!  iflag=0 writes the piezoelectric tensor as a function of temperature
!  iflag=2 writes the piezoelectric tensor as a function of pressure
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, code_group, code_group_ext, iflag
REAL(DP), INTENT(IN) :: temp(ntemp), piezom_tensor_t(3,6,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_piezo, ios
INTEGER :: find_free_unit
CHARACTER(LEN=7) :: label
!
! If the piezoelectric tensor vanishes return
!
SELECT CASE (code_group)
   CASE(2,16,18,19,20,22,23,25,27,29,31,32)
   RETURN
END SELECT

iu_piezo=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_piezo, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_piezom_on_file','opening piezo tensor (T) file',&
                                                             ABS(ios))
!
!  Choose if to plot as a function of temperature or pressure
!
IF (iflag<2) THEN
   label='  T(K) '
ELSE
   label='p(kbar)'
ENDIF

IF (meta_ionode) THEN

   WRITE(iu_piezo,'("#    piezomagnetic tensor in mu_B units")')
   WRITE(iu_piezo,'("#    multiply by XXX to have it in J/T/m^3")')
   SELECT CASE (code_group)
      CASE(30,28)
!
!   T_d, T
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_14 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,e20.12)') temp(itemp),  &
                piezom_tensor_t(1,4,itemp)
         ENDDO
      CASE(26)
!
!   S_4
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_31 ", &
                   &13x, " h_14 ", 13x, " h_15 ", 13x, " h_36 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,4e20.12)') temp(itemp),  &
                piezom_tensor_t(3,1,itemp),                 &  
                piezom_tensor_t(1,4,itemp),                 &
                piezom_tensor_t(1,5,itemp),                 &
                piezom_tensor_t(3,6,itemp)
         ENDDO
      CASE(24)
!
!   D_2d
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_14 ", &
                   &13x, " h_36 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,2e20.12)') temp(itemp),  &
                piezom_tensor_t(1,4,itemp),                 &  
                piezom_tensor_t(3,6,itemp)                
         ENDDO
      CASE(21)
!
!   D_3h
!
         IF (code_group_ext==107) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_22 ")') label

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,e20.12)') temp(itemp),  &
                               piezom_tensor_t(2,2,itemp)                
            ENDDO
         ELSEIF (code_group_ext==106) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ")') label

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,e20.12)') temp(itemp),  &
                               piezom_tensor_t(1,1,itemp)                
            ENDDO
         ENDIF

      CASE(17)
!
!   C_3h
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ", &
                                     &13x, " h_22 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,2e20.12)') temp(itemp),  &
                                   piezom_tensor_t(1,1,itemp),          &  
                                   piezom_tensor_t(2,2,itemp)                
         ENDDO

      CASE(14,15)
!
!   C_4v and C_6v
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_31 ", &
                  & 13x, "     h_33 ", 13x, "     h_15 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,3e20.12)') temp(itemp),  &
                piezom_tensor_t(3,1,itemp), &
                piezom_tensor_t(3,3,itemp), &
                piezom_tensor_t(1,5,itemp)
         ENDDO
      CASE(13)
!
!   C_3v
!
         IF (code_group_ext==72) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_22 ", &
                      &13x, " h_31 ", 13x, " h_33 ", 13x, " h_15 ")') label

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,4e20.12)') temp(itemp),  &
                      piezom_tensor_t(2,2,itemp),               &  
                      piezom_tensor_t(3,1,itemp),               &  
                      piezom_tensor_t(3,3,itemp),               &  
                      piezom_tensor_t(1,5,itemp)                
            ENDDO
         ELSEIF (code_group_ext==73) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ", &
                      &13x, " h_31 ", 13x, " h_33 ", 13x, " h_15 ")') label

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,4e20.12)') temp(itemp),  &
                      piezom_tensor_t(1,1,itemp),               &  
                      piezom_tensor_t(3,1,itemp),               &  
                      piezom_tensor_t(3,3,itemp),               &  
                      piezom_tensor_t(1,5,itemp)                
            ENDDO
         ENDIF
      CASE(12)
!
!   C_2v
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_31 ", &
                   &13x, " h_32 ", 13x, " h_33 ", 13x, " h_15 ",&
                   &13x, " h_24 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),       &
                          piezom_tensor_t(3,1,itemp),            &  
                          piezom_tensor_t(3,2,itemp),            &  
                          piezom_tensor_t(3,3,itemp),            &  
                          piezom_tensor_t(1,5,itemp),            &  
                          piezom_tensor_t(2,4,itemp)                
         ENDDO
      CASE(11,10)
!
!   D_6, D_4
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_14 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,e20.12)') temp(itemp),  &
                               piezom_tensor_t(1,4,itemp)                
         ENDDO

      CASE(9)
!
!   D_3
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ", 13x, " h_14 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,2e20.12)') temp(itemp),  &
                               piezom_tensor_t(1,1,itemp),  &  
                               piezom_tensor_t(1,4,itemp)
         ENDDO
      CASE(8)
!
!   D_2
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_14 ", &
                   &13x, " h_25 ", 13x, " h_36 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,3e20.12)') temp(itemp),  &
                piezom_tensor_t(1,4,itemp),                &  
                piezom_tensor_t(2,5,itemp),                &  
                piezom_tensor_t(3,6,itemp)  
         ENDDO
      CASE(6,7)
!
!  C_6  C_4
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_31 ", &
                      &13x, " h_33 ", 13x, " h_14 ", 13x, " h_15 ",&
                      &13x, " h_24 ")') label
         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                piezom_tensor_t(3,1,itemp),                 &  
                piezom_tensor_t(3,3,itemp),                 &  
                piezom_tensor_t(1,4,itemp),                 &  
                piezom_tensor_t(1,5,itemp),                 & 
                piezom_tensor_t(2,4,itemp)  
         ENDDO
      CASE(5)
!
!   C_3
!
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ", &
                   &13x, " h_22 ", 13x, " h_31 ", 13x, " h_33 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,6e20.12)') temp(itemp),  &
                piezom_tensor_t(1,1,itemp),                &  
                piezom_tensor_t(2,2,itemp),                &  
                piezom_tensor_t(3,1,itemp),                &  
                piezom_tensor_t(3,3,itemp),                &
                piezom_tensor_t(1,4,itemp),                &
                piezom_tensor_t(1,5,itemp)

         ENDDO
      CASE(4)
!
!   C_2
!
         IF (ibrav==-12.OR.ibrav==-13) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_21 ", &
                 &13x, " h_22 ", 13x, " h_23 ", 13x, " h_14 ", &
                 &13x, " h_16 ", 13x, " h_25 ", 13x, " h_34 ", &
                 &13x, " h_36 ")') label

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,8e20.12)') temp(itemp),  &
                   piezom_tensor_t(2,1,itemp),                &  
                   piezom_tensor_t(2,2,itemp),                &  
                   piezom_tensor_t(2,3,itemp),                &  
                   piezom_tensor_t(1,4,itemp),                &  
                   piezom_tensor_t(1,6,itemp),                &  
                   piezom_tensor_t(2,5,itemp),                &  
                   piezom_tensor_t(3,4,itemp),                &  
                   piezom_tensor_t(3,6,itemp)  
            ENDDO
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_21 ", &
                   &13x, " h_22 ", 13x, " h_23 ", 13x, " h_14 ", &
                   &13x, " h_15 ", 13x, " h_25 ", 13x, " h_24 ", &
                   &13x, " h_36 ")') label
            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,8e20.12)') temp(itemp),  &
                   piezom_tensor_t(2,1,itemp),                &  
                   piezom_tensor_t(2,2,itemp),                &  
                   piezom_tensor_t(2,3,itemp),                &  
                   piezom_tensor_t(1,4,itemp),                &  
                   piezom_tensor_t(1,5,itemp),                &  
                   piezom_tensor_t(2,5,itemp),                &  
                   piezom_tensor_t(2,4,itemp),                &  
                   piezom_tensor_t(3,6,itemp)  
            ENDDO
         ENDIF
      CASE(3)
!
!   C_s
!
         IF (ibrav==-12.OR.ibrav==-13) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ",             &
                   &13x, " h_12 ", 13x, " h_13 ", 13x, " h_15",  &
                   &13x, " h_24 ", 13x, " h_26 ", 13x, " h_31 ", &
                   &13x, " h_32 ", 13x, " h_33 ", 13x, " h_35 ")') label

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,10e20.12)') temp(itemp),  &
                   piezom_tensor_t(1,1,itemp),                &  
                   piezom_tensor_t(1,2,itemp),                &  
                   piezom_tensor_t(1,3,itemp),                &  
                   piezom_tensor_t(1,5,itemp),                &  
                   piezom_tensor_t(2,4,itemp),                &  
                   piezom_tensor_t(2,6,itemp),                &  
                   piezom_tensor_t(3,1,itemp),                &  
                   piezom_tensor_t(3,2,itemp),                &
                   piezom_tensor_t(3,3,itemp),                &
                   piezom_tensor_t(3,5,itemp)  
            ENDDO
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ",             &
                   &13x, " h_12 ", 13x, " h_13 ", 13x, " h_16",  &
                   &13x, " h_21 ", 13x, " h_22 ", 13x, " h_23 ", &
                   &13x, " h_26 ", 13x, " h_34 ", 13x, " h_35 ")') label

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,10e20.12)') temp(itemp),  &
                   piezom_tensor_t(1,1,itemp),                &  
                   piezom_tensor_t(1,2,itemp),                &  
                   piezom_tensor_t(1,3,itemp),                &  
                   piezom_tensor_t(1,6,itemp),                &  
                   piezom_tensor_t(2,1,itemp),                &  
                   piezom_tensor_t(2,2,itemp),                &  
                   piezom_tensor_t(2,3,itemp),                &  
                   piezom_tensor_t(2,6,itemp),                &
                   piezom_tensor_t(3,4,itemp),                &
                   piezom_tensor_t(3,5,itemp)  
            ENDDO
         ENDIF
      CASE DEFAULT
         WRITE(iu_piezo,'("#",5x, a7, 13x, " h_11 ",             &
                   &13x, " h_12 ", 13x, " h_13 ", 13x, " h_14",  &
                   &13x, " h_15 ", 13x, " h_16 ", 13x, " h_21 ", &
                   &13x, " h_22 ", 13x, " h_23 ", 13x, " h_24 ", &
                   &13x, " h_25 ", 13x, " h_26 ", 13x, " h_31 ", &
                   &13x, " h_32 ", 13x, " h_33 ", 13x, " h_34 ", &
                   &13x, " h_35 ", 13x, " h_36 ")') label

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,18e20.12)') temp(itemp), &
                  piezom_tensor_t(1,1,itemp), piezom_tensor_t(1,2,itemp), &
                  piezom_tensor_t(1,3,itemp), piezom_tensor_t(1,4,itemp), &
                  piezom_tensor_t(1,5,itemp), piezom_tensor_t(1,6,itemp), &
                  piezom_tensor_t(2,1,itemp), piezom_tensor_t(2,2,itemp), &
                  piezom_tensor_t(2,3,itemp), piezom_tensor_t(2,4,itemp), &
                  piezom_tensor_t(2,5,itemp), piezom_tensor_t(2,6,itemp), &
                  piezom_tensor_t(3,1,itemp), piezom_tensor_t(3,2,itemp), &
                  piezom_tensor_t(3,3,itemp), piezom_tensor_t(3,4,itemp), &
                  piezom_tensor_t(3,5,itemp), piezom_tensor_t(3,6,itemp)
         ENDDO
   END SELECT
   CLOSE(iu_piezo)
ENDIF

RETURN
END SUBROUTINE write_piezom_tensor_on_file
!
!-------------------------------------------------------------------------
SUBROUTINE read_piezom_tensor_from_file(temp, ntemp, ibrav, code_group, &
                               code_group_ext, piezom_tensor_t, filename)
!-------------------------------------------------------------------------
!
! This routine reads a file with the temperature (or pressure) and
! the inequivalent components of the piezomagnetic tensor.
! It fills the entire tensor.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, code_group, code_group_ext
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: piezom_tensor_t(3,6,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_piezo, ios
INTEGER :: find_free_unit
REAL(DP) :: rdum
CHARACTER(LEN=7) :: label
!
! If the piezomagnetic tensor vanishes return
!
SELECT CASE (code_group)
   CASE(2,16,18,19,20,22,23,25,27,29,31,32)
   RETURN
END SELECT

iu_piezo=find_free_unit()
piezom_tensor_t=0.0_DP
IF (meta_ionode) &
   OPEN(UNIT=iu_piezo, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('read_piezom_from_file','opening piezo tensor (T) file',&
                                                             ABS(ios))
IF (meta_ionode) THEN

   READ(iu_piezo,*)
   READ(iu_piezo,*)
   SELECT CASE (code_group)
      CASE(30,28)
!
!   T_d, T
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(1,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                  'incorrect temperature', 1)
         ENDDO
      CASE(26)
!
!   S_4
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(3,1,itemp),  &  
                piezom_tensor_t(1,4,itemp),  &
                piezom_tensor_t(1,5,itemp),  &
                piezom_tensor_t(3,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                      'incorrect temperature', 1)
         ENDDO
      CASE(24)
!
!   D_2d
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,*) rdum, piezom_tensor_t(1,4,itemp),   &  
                                           piezom_tensor_t(3,6,itemp) 
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                          'incorrect temperature', 1)
         ENDDO
      CASE(21)
!
!   D_3h
!
         READ(iu_piezo,*)

         IF (code_group_ext==107) THEN
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum, piezom_tensor_t(2,2,itemp)                
               IF (ABS(rdum-temp(itemp))>1D-5) &
                  CALL errore('read_piezom_from_file',&
                                     'incorrect temperature', 1)
            ENDDO
         ELSEIF (code_group_ext==106) THEN
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum, piezom_tensor_t(1,1,itemp)                
               IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file',&
                                     'incorrect temperature', 1)
            ENDDO
         ENDIF
      CASE(17)
!
!   C_3h
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(1,1,itemp),  &  
                piezom_tensor_t(2,2,itemp)                
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                  'incorrect temperature', 1)
         ENDDO

      CASE(14,15)
!
!   C_4v and C_6v
!
         READ(iu_piezo,*)
         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(3,1,itemp), &
                 piezom_tensor_t(3,3,itemp), piezom_tensor_t(1,5,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                    'incorrect temperature', 1)
         ENDDO
      CASE(13)
!
!   C_3v
!
         READ(iu_piezo,*)

         IF (code_group_ext==72) THEN
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum, piezom_tensor_t(2,2,itemp),  &  
                   piezom_tensor_t(3,1,itemp),                &  
                   piezom_tensor_t(3,3,itemp),                &  
                   piezom_tensor_t(1,5,itemp)                
               IF (ABS(rdum-temp(itemp))>1D-5) &
                  CALL errore('read_piezom_from_file', &
                                          'incorrect temperature', 1)
            ENDDO
         ELSEIF (code_group_ext==73) THEN
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum, piezom_tensor_t(1,1,itemp),  &  
                   piezom_tensor_t(3,1,itemp),                &  
                   piezom_tensor_t(3,3,itemp),                &  
                   piezom_tensor_t(1,5,itemp)                
               IF (ABS(rdum-temp(itemp))>1D-5) &
                  CALL errore('read_piezom_from_file', &
                                          'incorrect temperature', 1)
            ENDDO
         ENDIF
      CASE(12)
!
!   C_2v
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(3,1,itemp),   &  
                piezom_tensor_t(3,2,itemp),                &  
                piezom_tensor_t(3,3,itemp),                &  
                piezom_tensor_t(1,5,itemp),                &
                piezom_tensor_t(2,4,itemp)                
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                  'incorrect temperature', 1)
         ENDDO
      CASE(11,10)
!
!   D_6, D_4
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(1,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                  'incorrect temperature', 1)
         ENDDO

      CASE(9)
!
!   D_3
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(1,1,itemp),  &  
                piezom_tensor_t(1,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file', &
                                     'incorrect temperature', 1)
         ENDDO
      CASE(8)
!
!   D_2
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(1,4,itemp),  &  
                piezom_tensor_t(2,5,itemp),                &  
                piezom_tensor_t(3,6,itemp)  
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file',&
                                  'incorrect temperature', 1)
         ENDDO
      CASE(6, 7)
!
!   C_6 C_4
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(3,1,itemp),     &  
                piezom_tensor_t(3,3,itemp),                &  
                piezom_tensor_t(1,4,itemp),                &  
                piezom_tensor_t(1,5,itemp)                  
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file',&
                                 'incorrect temperature', 1)
         ENDDO
      CASE(5)
!
!   C_3
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, piezom_tensor_t(1,1,itemp),  &  
                piezom_tensor_t(2,2,itemp),                &  
                piezom_tensor_t(3,1,itemp),                &  
                piezom_tensor_t(3,3,itemp),                &  
                piezom_tensor_t(1,4,itemp),                &
                piezom_tensor_t(1,5,itemp)  
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezom_from_file',&
                                  'incorrect temperature', 1)
         ENDDO
      CASE(4)
!
!   C_2
!
         READ(iu_piezo,*)

         IF (ibrav==-12.OR.ibrav==-13) THEN
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum,                      &
                piezom_tensor_t(2,1,itemp),                &  
                piezom_tensor_t(2,2,itemp),                &  
                piezom_tensor_t(2,3,itemp),                &  
                piezom_tensor_t(1,4,itemp),                &  
                piezom_tensor_t(1,6,itemp),                &  
                piezom_tensor_t(2,5,itemp),                &  
                piezom_tensor_t(3,4,itemp),                &  
                piezom_tensor_t(3,6,itemp)  
                IF (ABS(rdum-temp(itemp))>1D-5) &
                   CALL errore('read_piezom_from_file',&
                                    'incorrect temperature',1)
            ENDDO
         ELSE
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum,                      &
                piezom_tensor_t(2,1,itemp),                &  
                piezom_tensor_t(2,2,itemp),                &  
                piezom_tensor_t(2,3,itemp),                &  
                piezom_tensor_t(1,4,itemp),                &  
                piezom_tensor_t(1,5,itemp),                &  
                piezom_tensor_t(2,5,itemp),                &  
                piezom_tensor_t(2,4,itemp),                &  
                piezom_tensor_t(3,6,itemp)  
                IF (ABS(rdum-temp(itemp))>1D-5) &
                   CALL errore('read_piezom_from_file',&
                                    'incorrect temperature',1)
            ENDDO
         ENDIF
      CASE(3)
!
!   C_s
!
         READ(iu_piezo,*)

         IF (ibrav==-12.OR.ibrav==-13) THEN
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum, piezom_tensor_t(1,1,itemp),  &  
                   piezom_tensor_t(1,2,itemp),                &  
                   piezom_tensor_t(1,3,itemp),                &  
                   piezom_tensor_t(1,5,itemp),                &  
                   piezom_tensor_t(2,4,itemp),                &  
                   piezom_tensor_t(2,6,itemp),                &  
                   piezom_tensor_t(3,1,itemp),                &  
                   piezom_tensor_t(3,2,itemp),                &
                   piezom_tensor_t(3,3,itemp),                &
                   piezom_tensor_t(3,5,itemp)  
                   IF (ABS(rdum-temp(itemp))>1D-5) &
                      CALL errore('read_piezom_from_file',&
                                            'incorrect temperature',1)
            ENDDO
         ELSE
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum, piezom_tensor_t(1,1,itemp),  &  
                   piezom_tensor_t(1,2,itemp),                &  
                   piezom_tensor_t(1,3,itemp),                &  
                   piezom_tensor_t(1,6,itemp),                &  
                   piezom_tensor_t(2,1,itemp),                &  
                   piezom_tensor_t(2,2,itemp),                &  
                   piezom_tensor_t(2,3,itemp),                &  
                   piezom_tensor_t(2,6,itemp),                &
                   piezom_tensor_t(3,4,itemp),                &
                   piezom_tensor_t(3,5,itemp)  
                   IF (ABS(rdum-temp(itemp))>1D-5) &
                      CALL errore('read_piezom_from_file',&
                                            'incorrect temperature',1)
            ENDDO
         ENDIF
      CASE DEFAULT
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, &
                  piezom_tensor_t(1,1,itemp), piezom_tensor_t(1,2,itemp), &
                  piezom_tensor_t(1,3,itemp), piezom_tensor_t(1,4,itemp), &
                  piezom_tensor_t(1,5,itemp), piezom_tensor_t(1,6,itemp), &
                  piezom_tensor_t(2,1,itemp), piezom_tensor_t(2,2,itemp), &
                  piezom_tensor_t(2,3,itemp), piezom_tensor_t(2,4,itemp), &
                  piezom_tensor_t(2,5,itemp), piezom_tensor_t(2,6,itemp), &
                  piezom_tensor_t(3,1,itemp), piezom_tensor_t(3,2,itemp), &
                  piezom_tensor_t(3,3,itemp), piezom_tensor_t(3,4,itemp), &
                  piezom_tensor_t(3,5,itemp), piezom_tensor_t(3,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
                CALL errore('read_piezom_from_file',&
                                              'incorrect temperature',1)
         ENDDO
   END SELECT
   CLOSE(iu_piezo)
   CALL expand_piezom_tensor(piezom_tensor_t, code_group, code_group_ext, &
                                                      ibrav, ntemp, temp)
ENDIF
CALL mp_bcast(piezom_tensor_t, meta_ionode_id, world_comm) 

RETURN
END SUBROUTINE read_piezom_tensor_from_file

!-------------------------------------------------------------------------
SUBROUTINE expand_piezom_tensor(piezom_tensor_t, code_group, code_group_ext, &
                              ibrav, ntemp, temp)
!-------------------------------------------------------------------------
!
! This routine reconstruct the complete piezomagnetic tensor from the
! components read from file
!
USE kinds,      ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, code_group, code_group_ext
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: piezom_tensor_t(3,6,ntemp)

INTEGER :: itemp
!
! If the piezomagnetic tensor vanishes return
!
SELECT CASE (code_group)
   CASE(2,16,18,19,20,22,23,25,27,29,31,32)
   RETURN
END SELECT

SELECT CASE (code_group)
   CASE(30,28)
!
!   T_d, T
!
         DO itemp=2,ntemp-1
            piezom_tensor_t(2,5,itemp) = piezom_tensor_t(1,4,itemp)
            piezom_tensor_t(3,6,itemp) = piezom_tensor_t(1,4,itemp)
         ENDDO
    CASE(26)
!
!   S_4
!
         DO itemp=2,ntemp-1
            piezom_tensor_t(3,2,itemp)=-piezom_tensor_t(3,1,itemp)
            piezom_tensor_t(2,5,itemp)= piezom_tensor_t(1,4,itemp)
            piezom_tensor_t(2,4,itemp)=-piezom_tensor_t(1,5,itemp)
         ENDDO
   CASE(24)
!
!   D_2d
!
         DO itemp=2,ntemp-1
            piezom_tensor_t(2,5,itemp)= piezom_tensor_t(1,4,itemp)
         ENDDO
   CASE(21)
!
!   D_3h
!
         IF (code_group_ext==107) THEN
            DO itemp=2,ntemp-1
               piezom_tensor_t(2,2,itemp)=-piezom_tensor_t(2,1,itemp)
               piezom_tensor_t(1,6,itemp)=-piezom_tensor_t(2,1,itemp)
            ENDDO
         ELSEIF (code_group_ext==106) THEN
            DO itemp=2,ntemp-1
               piezom_tensor_t(1,2,itemp)=-piezom_tensor_t(1,1,itemp)
               piezom_tensor_t(2,6,itemp)=-piezom_tensor_t(1,1,itemp)
            ENDDO
         ENDIF

   CASE(17)
!
!   C_3h
!

         DO itemp=2,ntemp-1
            piezom_tensor_t(1,2,itemp)=-piezom_tensor_t(1,1,itemp)
            piezom_tensor_t(2,1,itemp)=-piezom_tensor_t(2,2,itemp)
            piezom_tensor_t(1,6,itemp)= piezom_tensor_t(2,2,itemp)
            piezom_tensor_t(2,6,itemp)= piezom_tensor_t(1,1,itemp)
         ENDDO
   CASE(14,15)
!
!   C_4v and C_6v
!
         DO itemp=2,ntemp-1
            piezom_tensor_t(3,2,itemp)= piezom_tensor_t(3,1,itemp)
            piezom_tensor_t(2,4,itemp)= piezom_tensor_t(1,5,itemp)
         ENDDO
   CASE(13)
!
!   C_3v
!

         IF (code_group_ext==72) THEN
            DO itemp=2,ntemp-1
               piezom_tensor_t(2,1,itemp)=-piezom_tensor_t(2,2,itemp)
               piezom_tensor_t(1,6,itemp)=-piezom_tensor_t(2,2,itemp)
               piezom_tensor_t(3,2,itemp)= piezom_tensor_t(3,1,itemp)
               piezom_tensor_t(2,4,itemp)= piezom_tensor_t(1,5,itemp)
            ENDDO
         ELSEIF (code_group_ext==73) THEN
            DO itemp=2,ntemp-1
               piezom_tensor_t(1,2,itemp)=-piezom_tensor_t(1,1,itemp)
               piezom_tensor_t(2,6,itemp)=-piezom_tensor_t(1,1,itemp)
               piezom_tensor_t(3,2,itemp)= piezom_tensor_t(3,1,itemp)
               piezom_tensor_t(2,4,itemp)= piezom_tensor_t(1,5,itemp)
            ENDDO
         ENDIF
   CASE(12)
!
!   C_2v ! all elements are independent
!

   CASE(11,10)
!
!   D_6, D_4
!

         DO itemp=2,ntemp-1
            piezom_tensor_t(2,5,itemp)=-piezom_tensor_t(1,4,itemp)
         ENDDO

   CASE(9)
!
!   D_3
!

         DO itemp=2,ntemp-1
            piezom_tensor_t(1,2,itemp)=-piezom_tensor_t(1,1,itemp)
            piezom_tensor_t(2,5,itemp)=-piezom_tensor_t(1,4,itemp)
            piezom_tensor_t(2,6,itemp)=-piezom_tensor_t(1,1,itemp)
         ENDDO
   CASE(8)
!
!   D_2  ! all elements are independent
!
   CASE(6,7)
!
!   C_6 C_4
!
         DO itemp=2,ntemp-1
            piezom_tensor_t(3,2,itemp)= piezom_tensor_t(3,1,itemp)
            piezom_tensor_t(2,5,itemp)=-piezom_tensor_t(1,4,itemp)
            piezom_tensor_t(2,4,itemp)= piezom_tensor_t(1,5,itemp)
         ENDDO
   CASE(5)
!
!   C_3
!
       DO itemp=2,ntemp-1
          piezom_tensor_t(1,2,itemp)=-piezom_tensor_t(1,1,itemp)
          piezom_tensor_t(2,6,itemp)=-piezom_tensor_t(1,1,itemp)
          piezom_tensor_t(2,1,itemp)=-piezom_tensor_t(2,2,itemp)
          piezom_tensor_t(1,6,itemp)=-piezom_tensor_t(2,2,itemp)
          piezom_tensor_t(3,2,itemp)= piezom_tensor_t(3,1,itemp)
          piezom_tensor_t(2,5,itemp)=-piezom_tensor_t(1,4,itemp)
          piezom_tensor_t(2,4,itemp)= piezom_tensor_t(1,5,itemp)
       ENDDO

   CASE DEFAULT ! CASE(1,3,4)
!
!     No symmetry in the default case
!
END SELECT

RETURN
END SUBROUTINE expand_piezom_tensor
!
!-----------------------------------------------------------------------
FUNCTION get_pm_type(code_group, code_group_ext, ibrav)
!-----------------------------------------------------------------------
INTEGER :: get_pm_type
INTEGER, INTENT(IN) :: code_group, code_group_ext, ibrav

INTEGER :: itype, aux_type

aux_type=0
DO itype=1,pm_types
   IF (pm_code_group(itype)==code_group) aux_type=itype
ENDDO
!
!  Special cases
!
IF (code_group==3) THEN
!
! C_s b-unique or c-unique
!
   IF (ibrav==-12.OR.ibrav==-13) THEN
      aux_type=22
   ELSEIF (ibrav==12.OR.ibrav==13) THEN
      aux_type=23
   ENDIF
ENDIF

IF (code_group==4) THEN
!
! C_2 b-unique or c-unique
!
   IF (ibrav==-12.OR.ibrav==-13) THEN
      aux_type=20
   ELSEIF (ibrav==12.OR.ibrav==13) THEN
      aux_type=21
   ENDIF
ENDIF

IF (code_group==13) THEN
!
!  C_3v mirror perpendicular to x1 or to x2
!
   IF (code_group_ext==72) THEN
      aux_type=10
   ELSEIF (code_group_ext==73) THEN
      aux_type=11
   ENDIF
ENDIF

IF (code_group==21) THEN
!
!  D_3h mirror perpendicular to x1 or to x2
!
   IF (code_group_ext==106) THEN
      aux_type=5
   ELSEIF (code_group_ext==107) THEN
      aux_type=6
   ENDIF
ENDIF

IF (aux_type==0) CALL errore('get_pt_type','code_group not available',1)
get_pm_type=aux_type
RETURN
END FUNCTION get_pm_type

END MODULE piezomagnetic_tensor
