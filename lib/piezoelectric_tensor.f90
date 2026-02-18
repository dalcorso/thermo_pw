!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE piezoelectric_tensor
!---------------------------------------------------------------------------
!
!   this module contains the support routines for the calculation
!   of the piezoelectric tensor
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL(DP) :: g_piezo_tensor(3,6)   ! The piezoelectric tensor g_{\alpha,m}
                                    ! improper piezoelectric tensor linking
                                    ! polarization and strain
  REAL(DP) :: eg_piezo_tensor(3,6)  ! The piezoelectric tensor eg_{\alpha,m}
                                    ! proper piezoelectric tensor obtained
                                    ! transforming the improper one
  REAL(DP) :: e_piezo_tensor(3,6)   ! The piezoelectric tensor e_{\alpha,m}
                                    ! proper piezoelectric tensor linking
                                    ! polarization and strain
  REAL(DP) :: d_piezo_tensor(3,6)   ! The piezoelectric tensor d_{\alpha,m}
                                    ! linking polarization and stress
!
!  The following four variables are auxiliary, they can contain another
!  copy of the piezoelectric tensor. Usually the clamped ion term
! 
  REAL(DP) :: g_piezo_tensor_fi(3,6)   ! The piezoelectric tensor g_{\alpha,m}
                                    ! improper piezoelectric tensor linking
                                    ! polarization and strain
  REAL(DP) :: eg_piezo_tensor_fi(3,6)  ! The piezoelectric tensor eg_{\alpha,m}
                                    ! proper piezoelectric tensor obtained
                                    ! transforming the improper one
  REAL(DP) :: e_piezo_tensor_fi(3,6)   ! The piezoelectric tensor e_{\alpha,m}
                                    ! proper piezoelectric tensor linking
                                    ! polarization and strain
  REAL(DP) :: d_piezo_tensor_fi(3,6)   ! The piezoelectric tensor d_{\alpha,m}
                                    ! linking polarization and stress

  REAL(DP), ALLOCATABLE :: polar_strain(:,:) ! The polarization for each strain
                                    ! in units of e bohr/Omega, Omega in bohr^3
  REAL(DP), ALLOCATABLE :: tot_b_phase(:,:) ! Total Berry phase (elec. + ions)
                                    ! in each direction for each strain

  INTEGER :: nppl      ! number of points per line in berry phase calculation

!
!   Some array to simplify dealing with piezoelectric tensor
!
  INTEGER, PARAMETER :: pt_elements=18

  CHARACTER(LEN=6) :: pt_names(pt_elements)

  DATA  pt_names /                                                   &
         'e_{11}', 'e_{12}', 'e_{13}', 'e_{14}', 'e_{15}', 'e_{16}', &
         'e_{21}', 'e_{22}', 'e_{23}', 'e_{24}', 'e_{25}', 'e_{26}', &
         'e_{31}', 'e_{32}', 'e_{33}', 'e_{34}', 'e_{35}', 'e_{36}' /


  CHARACTER(LEN=6) :: ptd_names(pt_elements)

  DATA  ptd_names /                                                  &
         'd_{11}', 'd_{12}', 'd_{13}', 'd_{14}', 'd_{15}', 'd_{16}', &
         'd_{21}', 'd_{22}', 'd_{23}', 'd_{24}', 'd_{25}', 'd_{26}', &
         'd_{31}', 'd_{32}', 'd_{33}', 'd_{34}', 'd_{35}', 'd_{36}' /

  INTEGER, PARAMETER :: pt_types=20

  INTEGER :: pt_code_group(pt_types)  ! code of the point group for each type
  DATA  pt_code_group / 30, 28, 26, 24, 21, 17, 15, 14, 13, 12, 11, 10, &
                         9, 8, 7, 6, 5, 4, 3, 1 /

  INTEGER  :: pt_present(pt_elements, pt_types)

  DATA pt_present / &
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 30  T_d
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 28  T
       0,0,0,2,3,0, 0,0,0,0,0,0, 1,0,0,0,0,4, & ! 26  S_4
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,2, & ! 24  D_2d
       1,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 21  D_3h
       1,0,0,0,0,0, 0,2,0,0,0,0, 0,0,0,0,0,0, & ! 17  C_3h
       0,0,0,0,3,0, 0,0,0,0,0,0, 1,0,2,0,0,0, & ! 15  C_6v
       0,0,0,0,3,0, 0,0,0,0,0,0, 1,0,2,0,0,0, & ! 14  C_4v
       0,0,0,0,4,0, 0,1,0,0,0,0, 2,0,3,0,0,0, & ! 13  C_3v
       0,0,0,0,4,0, 0,0,0,5,0,0, 1,2,3,0,0,0, & ! 12  C_2v
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 11  D_6
       0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & ! 10  D_4
       1,0,0,2,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, & !  9  D_3
       0,0,0,1,0,0, 0,0,0,0,2,0, 0,0,0,0,0,3, & !  8  D_2
       0,0,0,2,3,0, 0,0,0,0,0,0, 1,0,0,0,0,0, & !  7  C_6
       0,0,0,3,4,0, 0,0,0,0,0,0, 1,0,2,0,0,0, & !  6  C_4
       1,0,0,5,6,0, 0,2,0,0,0,0, 3,0,4,0,0,0, & !  5  C_3
       0,0,0,4,0,5, 1,2,3,0,6,0, 0,0,0,7,0,8, & !  4  C_2
       1,2,3,0,4,0, 0,0,0,5,0,6, 7,8,9,0,10,0, & !  3  C_s
       1,2,3,4,5,6, 7,8,9,10,11,12, 13,14,15,16,17,18 / !  1  C_1

  PUBLIC g_piezo_tensor, polar_strain, compute_improper_piezo_tensor, &
         compute_proper_piezo_tensor,                     &
         print_piezo_tensor, e_piezo_tensor,              &
         eg_piezo_tensor, pt_elements,                    &
         compute_d_piezo_tensor, d_piezo_tensor, nppl,    &
         compute_polarization_equil,                      &
         proper_improper_piezo, clean_piezo_tensor,       &
         print_piezo_info, tot_b_phase, allocate_piezo,   &
         deallocate_piezo, write_piezo_tensor,            &
         read_piezo_tensor, write_piezo_tensor_on_file,   &
         read_piezo_tensor_from_file,                     &
         read_piezo_tensor_fi,                            &
         pt_names, ptd_names, pt_types, pt_present,       &
         pt_code_group, get_pt_type, compute_relax_piezo, &
         e_piezo_tensor_fi, eg_piezo_tensor_fi,           &
         d_piezo_tensor_fi

CONTAINS
!
!---------------------------------------------------------------------------
SUBROUTINE print_piezo_tensor(piezo_tensor, fact, label, frozen_ions)
!---------------------------------------------------------------------------
!
!  This routine writes on output the piezoelectric tensor
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: label
LOGICAL, INTENT(IN) :: frozen_ions
REAL(DP), INTENT(IN) :: piezo_tensor(3,6)
REAL(DP) :: fact
INTEGER :: i, j
CHARACTER(LEN=30) :: fi_string

fi_string=''
IF (frozen_ions) fi_string="Frozen ions"

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x, a)') TRIM(label)
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
!
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (piezo_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_piezo_tensor

!-------------------------------------------------------------------------
SUBROUTINE allocate_piezo(nwork)
!-------------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: nwork

ALLOCATE(polar_strain(3,nwork))
ALLOCATE(tot_b_phase(3,nwork))
polar_strain=0.0_DP
tot_b_phase=0.0_DP
g_piezo_tensor=0.0_DP
eg_piezo_tensor=0.0_DP
e_piezo_tensor=0.0_DP
d_piezo_tensor=0.0_DP

RETURN
END SUBROUTINE allocate_piezo
!
!-------------------------------------------------------------------------
SUBROUTINE deallocate_piezo()
!-------------------------------------------------------------------------
IMPLICIT NONE

IF (ALLOCATED(polar_strain)) DEALLOCATE(polar_strain)
IF (ALLOCATED(tot_b_phase)) DEALLOCATE(tot_b_phase)

RETURN
END SUBROUTINE deallocate_piezo
!
!-------------------------------------------------------------------------
SUBROUTINE write_piezo_tensor(filename,polar0)
!-------------------------------------------------------------------------
!
!  This routine writes the piezoelectric tensor on file.
!  It must be called after computing the piezoelectric tensor.
!  It saves: 
!  the spontaneous polarization of the unperturbed geometry
!  the improper piezoelectric tensor
!  the corrected proper piezoelectric tensor
!  the proper piezoelectric tensor
!
USE io_global, ONLY : ionode, ionode_id
USE constants, ONLY : electron_si, bohr_radius_si
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: polar0(3)
REAL(DP) :: fact
INTEGER :: find_free_unit
INTEGER :: outunit, ios, i, j

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_piezo_tensor','ploblem opening output file', ABS(ios))

fact= electron_si / (bohr_radius_si)**2
IF (ionode) THEN
   WRITE(outunit,'("Spontaneous polarization (e/bohr**2)")')
   WRITE(outunit,'(3e19.10)') (polar0(i), i=1,3)
   WRITE(outunit,*)
   WRITE(outunit,'("Improper piezoelectric tensor (e/bohr**2)")')
   DO i=1,3
      WRITE(outunit,'(6e19.10)') (g_piezo_tensor(i,j), j=1,6)
   ENDDO
   WRITE(outunit,*)
   WRITE(outunit,'("Proper corrected piezoelectric tensor (e/bohr**2)")')
   DO i=1,3
      WRITE(outunit,'(6e19.10)') (eg_piezo_tensor(i,j), j=1,6)
   END DO
   WRITE(outunit,*)
   WRITE(outunit,'("Proper piezoelectric tensor (e/bohr**2)")')
   DO i=1,3
      WRITE(outunit,'(6e19.10)') (e_piezo_tensor(i,j), j=1,6)
   END DO
   WRITE(outunit,*)
   WRITE(outunit,'("Proper piezoelectric tensor (C/m**2)")')
   DO i=1,3
      WRITE(outunit,'(6e19.10)') (e_piezo_tensor(i,j)*fact, j=1,6)
   END DO
   WRITE(outunit,*)
   WRITE(outunit,'("Strain piezoelectric tensor (pC/N)")')
   DO i=1,3
      WRITE(outunit,'(6e19.10)') (d_piezo_tensor(i,j)*fact*1.D4, j=1,6)
   END DO
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_piezo_tensor

!-------------------------------------------------------------------------
SUBROUTINE read_piezo_tensor(filename, polar0, exists)
!-------------------------------------------------------------------------
!
!  This routine reads the piezoelectric tensor from file.
!
USE io_global, ONLY : ionode, ionode_id
USE constants, ONLY : electron_si, bohr_radius_si
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP) :: polar0(3)
LOGICAL, INTENT(OUT) :: exists
REAL(DP) :: e_piezo_tensor_cm2(3,6), fact
INTEGER :: inunit, ios, i, j
INTEGER :: find_free_unit

IF (ionode) THEN
   inunit=find_free_unit()
   OPEN(UNIT=inunit, FILE=TRIM(filename), STATUS='old', FORM='formatted', &
       ERR=100, IOSTAT=ios)
ENDIF

fact = electron_si / bohr_radius_si**2
IF (ionode) THEN
   READ(inunit,*)
   READ(inunit,'(3e19.10)') (polar0(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (g_piezo_tensor(i,j), j=1,6)
   ENDDO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (eg_piezo_tensor(i,j), j=1,6)
   END DO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (e_piezo_tensor(i,j), j=1,6)
   END DO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (e_piezo_tensor_cm2(i,j),&
                                                                       j=1,6)
   END DO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (d_piezo_tensor(i,j), j=1,6)
   END DO
   d_piezo_tensor=d_piezo_tensor / fact / 1.D-4 ! bring it back to units 
                                                ! e/(bohr**2 * kbar)
   CLOSE(inunit)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
IF (ios /= 0) THEN
   exists=.FALSE.
   RETURN
ENDIF
CALL mp_bcast(polar0,ionode_id,intra_image_comm)
CALL mp_bcast(g_piezo_tensor,ionode_id,intra_image_comm)
CALL mp_bcast(eg_piezo_tensor,ionode_id,intra_image_comm)
CALL mp_bcast(e_piezo_tensor,ionode_id,intra_image_comm)
CALL mp_bcast(d_piezo_tensor,ionode_id,intra_image_comm)
exists=.TRUE.

RETURN
END SUBROUTINE read_piezo_tensor
!
!-------------------------------------------------------------------------
SUBROUTINE read_piezo_tensor_fi(filename, polar0_fi, exists)
!-------------------------------------------------------------------------
!
!  This routine reads the piezoelectric tensor from file and
!  save it in the auxiliary variable with the fi extension
!
USE io_global, ONLY : ionode, ionode_id
USE constants, ONLY : electron_si, bohr_radius_si
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP) :: polar0_fi(3)
LOGICAL, INTENT(OUT) :: exists
REAL(DP) :: e_piezo_tensor_cm2(3,6), fact
INTEGER :: inunit, ios, i, j
INTEGER :: find_free_unit

IF (ionode) THEN
   inunit=find_free_unit()
   OPEN(UNIT=inunit, FILE=TRIM(filename), STATUS='old', FORM='formatted', &
       ERR=100, IOSTAT=ios)
ENDIF

fact = electron_si / bohr_radius_si**2
IF (ionode) THEN
   READ(inunit,*)
   READ(inunit,'(3e19.10)') (polar0_fi(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (g_piezo_tensor_fi(i,j), j=1,6)
   ENDDO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (eg_piezo_tensor_fi(i,j), j=1,6)
   END DO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (e_piezo_tensor_fi(i,j), j=1,6)
   END DO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (e_piezo_tensor_cm2(i,j),&
                                                                       j=1,6)
   END DO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(6e19.10)',ERR=100,IOSTAT=ios) (d_piezo_tensor_fi(i,j), j=1,6)
   END DO
   d_piezo_tensor_fi=d_piezo_tensor_fi / fact / 1.D-4 ! bring it back to units 
                                                ! e/(bohr**2 * kbar)
   CLOSE(inunit)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
IF (ios /= 0) THEN
   exists=.FALSE.
   RETURN
ENDIF
CALL mp_bcast(polar0_fi,ionode_id,intra_image_comm)
CALL mp_bcast(g_piezo_tensor_fi,ionode_id,intra_image_comm)
CALL mp_bcast(eg_piezo_tensor_fi,ionode_id,intra_image_comm)
CALL mp_bcast(e_piezo_tensor_fi,ionode_id,intra_image_comm)
CALL mp_bcast(d_piezo_tensor_fi,ionode_id,intra_image_comm)
exists=.TRUE.

RETURN
END SUBROUTINE read_piezo_tensor_fi
!---------------------------------------------------------------------------
SUBROUTINE compute_improper_piezo_tensor(polar_geo, epsil_geo, nwork, &
                                ngeo, ibrav, code_group)
!---------------------------------------------------------------------------
!
!  This routine computes the piezoelectric tensor g_{\alpha,m} by fitting the 
!  polarization strain relation with a second order polynomial. This is 
!  calculated on the basis of the solid point group.
!  The polarization enters in units of e/bohr^2
!
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: polar_geo(3,nwork), epsil_geo(3,3,nwork)
INTEGER, INTENT(IN) :: ngeo, ibrav, code_group, nwork
INTEGER :: i, j, igeo, alpha, ind, mn
LOGICAL :: check_group_ibrav

WRITE(stdout,'(/,20x,40("-"),/)')
g_piezo_tensor=0.0_DP
IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group) 
     CASE(2,16,18,19,20,22,23,25,27,29,32) 
     CASE(3)
!
!  C_s   Monoclinic
!
!        WRITE(stdout,'(5x,"( g11  g12  g13   .   d15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .   g24   .   d26 )")') 
!        WRITE(stdout,'(5x,"( g31  g32  g33   .   d35   .  )")') 
!
        CALL piezo_ij(1, 1, ngeo, epsil_geo, polar_geo )
        CALL piezo_ij(1, 2, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1) )
        CALL piezo_ij(1, 3, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1) )
        CALL piezo_ij(2, 6, ngeo, epsil_geo(1,1,5*ngeo+1), polar_geo(1,5*ngeo+1) )

        IF (ibrav==-12) THEN

           CALL piezo_ij(3, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
           CALL piezo_ij(3, 2, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
           CALL piezo_ij(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                     polar_geo(1,2*ngeo+1))
           CALL piezo_ij(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                     polar_geo(1,3*ngeo+1))
           CALL piezo_ij(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                     polar_geo(1,4*ngeo+1))
           CALL piezo_ij(3, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                     polar_geo(1,4*ngeo+1))

        ELSE
!                WRITE(stdout,'(5x,"( d11  d12  d13   .    .   d16 )")') 
!                WRITE(stdout,'(5x,"( d21  d22  d23   .    .   d26 )")') 
!                WRITE(stdout,'(5x,"(  .    .    .   d16  d26   .  )")') 
!
           CALL piezo_ij(2, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
           CALL piezo_ij(2, 2, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
           CALL piezo_ij(2, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                            polar_geo(1,2*ngeo+1))
           CALL piezo_ij(1, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                            polar_geo(1,5*ngeo+1))

        ENDIF
     CASE(4)
!
!  C_2   Monoclinic
!
        CALL piezo_ij(1, 4, ngeo, epsil_geo(1,1,3*ngeo+1), polar_geo(1,3*ngeo+1))
        CALL piezo_ij(2, 5, ngeo, epsil_geo(1,1,4*ngeo+1), polar_geo(1,4*ngeo+1))
        CALL piezo_ij(3, 6, ngeo, epsil_geo(1,1,5*ngeo+1), polar_geo(1,5*ngeo+1))

        IF (ibrav==-12) THEN
!            WRITE(stdout,'(5x,"(  .    .    .   d14   .   d16 )")') 
!            WRITE(stdout,'(5x,"( d21  d22  d23   .   d25   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   d34   .   d36 )")') 

           CALL piezo_ij(2, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
           CALL piezo_ij(2, 2, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
           CALL piezo_ij(2, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                     polar_geo(1,2*ngeo+1))
           CALL piezo_ij(1, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                                     polar_geo(1,5*ngeo+1))
           CALL piezo_ij(3, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                                     polar_geo(1,3*ngeo+1))
        ELSE
!            WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   d24  d25   .  )")') 
!            WRITE(stdout,'(5x,"( d31  d32  d33   .    .   d36 )")') 
!
            CALL piezo_ij(3, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
            CALL piezo_ij(3, 2, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
            CALL piezo_ij(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                      polar_geo(1,2*ngeo+1))
            CALL piezo_ij(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                      polar_geo(1,4*ngeo+1))
            CALL piezo_ij(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                      polar_geo(1,3*ngeo+1))

         ENDIF

      CASE(6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!
!             WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .   d24 -d14   .  )")') 
!             WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 
!
        CALL piezo_ij(3, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(3,2) = g_piezo_tensor(3,1)
        CALL piezo_ij(3, 3, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
        CALL piezo_ij(1, 4, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1))
        g_piezo_tensor(2,5) = -g_piezo_tensor(1,4)
        CALL piezo_ij(2, 4, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1))
        CALL piezo_ij(1, 5, ngeo, epsil_geo(1,1,3*ngeo+1), polar_geo(1,3*ngeo+1))

     CASE(8)
!
!  D_2 (222) Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   d25   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   d36 )")') 

        CALL piezo_ij(1, 4, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        CALL piezo_ij(2, 5, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
        CALL piezo_ij(3, 6, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1))

      CASE(9)
!
! D_3  Trigonal 
!
!         WRITE(stdout,'(5x,"( d11 -d11   .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -d14 2d11 )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezo_ij(1, 1, ngeo, epsil_geo, polar_geo)
        g_piezo_tensor(1,2) = -g_piezo_tensor(1,1)
        g_piezo_tensor(2,6) = 2.0_DP * g_piezo_tensor(1,1)
        CALL piezo_ij(1, 4, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
        g_piezo_tensor(2,5) = -g_piezo_tensor(1,4)

     CASE(10,11)
!
! D_4  tetragonal, D_6 hexagonal
!
!         WRITE(stdout,'(/,5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -d14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezo_ij(1, 4, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(2,5) = -g_piezo_tensor(1,4)

     CASE(12)
!
! C_2v  Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   d24   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d32  d33   .    .    .  )")') 

        CALL piezo_ij(3, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        CALL piezo_ij(3, 2, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
        CALL piezo_ij(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1))
        CALL piezo_ij(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), polar_geo(1,3*ngeo+1))
        CALL piezo_ij(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), polar_geo(1,4*ngeo+1))

     CASE(13)
!
! C_3v  Trigonal. Assuming m perpendicular to x1
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15 -d21 )")') 
!         WRITE(stdout,'(5x,"( d21 -d21   .   d15   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 
!
        CALL piezo_ij(2, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(2,2) = -g_piezo_tensor(2,1)
        g_piezo_tensor(1,6) = g_piezo_tensor(2,2)
        CALL piezo_ij(3, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(3,2) = g_piezo_tensor(3,1)
        CALL piezo_ij(3, 3, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
        CALL piezo_ij(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1))
        g_piezo_tensor(2,4) = g_piezo_tensor(1,5)

     CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   d15   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 

        CALL piezo_ij(3, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(3,1) = g_piezo_tensor(3,1) * 0.5_DP
        g_piezo_tensor(3,2) = g_piezo_tensor(3,1)
        CALL piezo_ij(3, 3, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
        CALL piezo_ij(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1))
        g_piezo_tensor(2,4) = g_piezo_tensor(1,5)

     CASE(17)
!
! C_3h hexagonal
!
!             WRITE(stdout,'(5x,"( d11 -d11   .    .    .  -d12 )")') 
!             WRITE(stdout,'(5x,"( d12 -d12   .    .    .   d11 )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezo_ij(1, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(1,2) = -g_piezo_tensor(1,1)
        g_piezo_tensor(2,6) = g_piezo_tensor(1,1)
        CALL piezo_ij(2, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(2,2) = -g_piezo_tensor(2,1)
        g_piezo_tensor(1,6) = -g_piezo_tensor(2,1)

      CASE(21)
!
! D_3h hexagonal
!
!             WRITE(stdout,'(5x,"(  .    .    .    .    .  -d21 )")') 
!             WRITE(stdout,'(5x,"( d21 -d21   .    .    .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
        CALL piezo_ij(2, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(2,2) = -g_piezo_tensor(2,1)
        g_piezo_tensor(1,6) = -g_piezo_tensor(1,2)

     CASE(24)
!
! D_2d tetragonal: axis 2 || x1
!
!         WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   d14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   d36 )")') 

        CALL piezo_ij(1, 4, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(2,5) = g_piezo_tensor(1,4)
        CALL piezo_ij(3, 6, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))

     CASE(26)
!
! S_4 tetragonal
!
!        WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .  -d15  d14   .  )")') 
!        WRITE(stdout,'(5x,"( d31 -d31   .    .    .   d36 )")') 

        CALL piezo_ij(1, 4, ngeo, epsil_geo(1,1,ngeo+1), polar_geo(1,ngeo+1))
        g_piezo_tensor(2,5) = g_piezo_tensor(1,4)
        CALL piezo_ij(3, 1, ngeo, epsil_geo(1,1,1), polar_geo(1,1))
        g_piezo_tensor(3,2) = -g_piezo_tensor(3,1)
        CALL piezo_ij(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1), polar_geo(1,2*ngeo+1))
        g_piezo_tensor(2,4) = -g_piezo_tensor(1,5)
        CALL piezo_ij(3, 6, ngeo, epsil_geo(1,1,3*ngeo+1), polar_geo(1,3*ngeo+1))

     CASE(28,30)
!
! T, T_d cubic
!
!             WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .   d14   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .   d14 )")') 
!
        CALL piezo_ij(1, 4, ngeo, epsil_geo, polar_geo)
        g_piezo_tensor(2,5) = g_piezo_tensor(1,4)
        g_piezo_tensor(3,6) = g_piezo_tensor(1,4)
     CASE(31)

     CASE DEFAULT
!
!  C_1 
!
        DO mn=1,6
           ind = ngeo * (mn-1)
           DO alpha=1,3
              CALL piezo_ij(alpha, mn, ngeo, epsil_geo(1,1,ind+1), &
                                             polar_geo(1,ind+1))
           ENDDO
        ENDDO
  END SELECT
ELSE
   DO mn=1,6
      ind = ngeo * (mn-1)
      DO alpha=1,3
         CALL piezo_ij(alpha, mn, ngeo, epsil_geo(1,1,ind+1), &
                                        polar_geo(1,ind+1))
      ENDDO
   ENDDO
ENDIF
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE compute_improper_piezo_tensor

!---------------------------------------------------------------------------
SUBROUTINE compute_proper_piezo_tensor(tot_b_phase, epsil_geo, nwork, ngeo, &
                                           ibrav, code_group, at)
!---------------------------------------------------------------------------
!
!  This routine computes the proper piezoelectric tensor e_{\alpha,m} by 
!  fitting the total berry phase strain relation with a second order 
!  polynomial and computing its derivative with respect to strain. 
!  Finally the proper piezoelectric tensor is computed from Eq.24 of
!  D. Vanderbilt Jour. of Phys. and Chemistry of Solids 61, 147 (2000).
!  This is calculated on the basis of the solid point group.
!
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: tot_b_phase(3,nwork), epsil_geo(3,3,nwork), &
                        at(3,3)
INTEGER, INTENT(IN) :: ngeo, ibrav, code_group, nwork
INTEGER :: i, j, igeo, alpha, ind, mn
LOGICAL :: check_group_ibrav

WRITE(stdout,'(/,20x,40("-"),/)')
e_piezo_tensor=0.0_DP
IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group) 
     CASE(2,16,18,19,20,22,23,25,27,29,32) 
     CASE(3)
!
!  C_s   Monoclinic
!
!        WRITE(stdout,'(5x,"( g11  g12  g13   .   d15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .   g24   .   d26 )")') 
!        WRITE(stdout,'(5x,"( g31  g32  g33   .   d35   .  )")') 
!
        CALL piezo_ij_bp(1, 1, ngeo, epsil_geo, tot_b_phase, at )
        CALL piezo_ij_bp(1, 2, ngeo, epsil_geo(1,1,ngeo+1), tot_b_phase(1,ngeo+1), at )
        CALL piezo_ij_bp(1, 3, ngeo, epsil_geo(1,1,2*ngeo+1), tot_b_phase(1,2*ngeo+1), at )
        CALL piezo_ij_bp(2, 6, ngeo, epsil_geo(1,1,5*ngeo+1), tot_b_phase(1,5*ngeo+1), at )

        IF (ibrav==-12) THEN

           CALL piezo_ij_bp(3, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
           CALL piezo_ij_bp(3, 2, ngeo, epsil_geo(1,1,ngeo+1), &
                                                 tot_b_phase(1,ngeo+1),at)
           CALL piezo_ij_bp(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                     tot_b_phase(1,2*ngeo+1), at)
           CALL piezo_ij_bp(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                     tot_b_phase(1,3*ngeo+1), at)
           CALL piezo_ij_bp(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                     tot_b_phase(1,4*ngeo+1), at)
           CALL piezo_ij_bp(3, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                     tot_b_phase(1,4*ngeo+1), at)

        ELSE
!                WRITE(stdout,'(5x,"( d11  d12  d13   .    .   d16 )")') 
!                WRITE(stdout,'(5x,"( d21  d22  d23   .    .   d26 )")') 
!                WRITE(stdout,'(5x,"(  .    .    .   d16  d26   .  )")') 
!
           CALL piezo_ij_bp(2, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
           CALL piezo_ij_bp(2, 2, ngeo, epsil_geo(1,1,ngeo+1), &
                           tot_b_phase(1,ngeo+1), at)
           CALL piezo_ij_bp(2, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                            tot_b_phase(1,2*ngeo+1), at)
           CALL piezo_ij_bp(1, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                            tot_b_phase(1,5*ngeo+1), at)

        ENDIF
     CASE(4)
!
!  C_2   Monoclinic
!
        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                              tot_b_phase(1,3*ngeo+1), at)
        CALL piezo_ij_bp(2, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                              tot_b_phase(1,4*ngeo+1), at)
        CALL piezo_ij_bp(3, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                              tot_b_phase(1,5*ngeo+1), at)

        IF (ibrav==-12) THEN
!            WRITE(stdout,'(5x,"(  .    .    .   d14   .   d16 )")') 
!            WRITE(stdout,'(5x,"( d21  d22  d23   .   d25   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   d34   .   d36 )")') 

           CALL piezo_ij_bp(2, 1, ngeo, epsil_geo(1,1,1), &
                                              tot_b_phase(1,1), at)
           CALL piezo_ij_bp(2, 2, ngeo, epsil_geo(1,1,ngeo+1),   &       
                                              tot_b_phase(1,ngeo+1), at)
           CALL piezo_ij_bp(2, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                              tot_b_phase(1,2*ngeo+1), at)
           CALL piezo_ij_bp(1, 6, ngeo, epsil_geo(1,1,5*ngeo+1), &
                                              tot_b_phase(1,5*ngeo+1), at)
           CALL piezo_ij_bp(3, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                              tot_b_phase(1,3*ngeo+1), at)
        ELSE
!            WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   d24  d25   .  )")') 
!            WRITE(stdout,'(5x,"( d31  d32  d33   .    .   d36 )")') 
!
            CALL piezo_ij_bp(3, 1, ngeo, epsil_geo(1,1,1), &
                                               tot_b_phase(1,1), at)
            CALL piezo_ij_bp(3, 2, ngeo, epsil_geo(1,1,ngeo+1), &
                                               tot_b_phase(1,ngeo+1), at)
            CALL piezo_ij_bp(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                      tot_b_phase(1,2*ngeo+1), at)
            CALL piezo_ij_bp(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                      tot_b_phase(1,4*ngeo+1), at)
            CALL piezo_ij_bp(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                      tot_b_phase(1,3*ngeo+1), at)

         ENDIF

      CASE(6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!
!             WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .   d24 -d14   .  )")') 
!             WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 
!
        CALL piezo_ij_bp(3, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        g_piezo_tensor(3,2) = g_piezo_tensor(3,1)
        CALL piezo_ij_bp(3, 3, ngeo, epsil_geo(1,1,ngeo+1), &
                                                tot_b_phase(1,ngeo+1), at)
        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                tot_b_phase(1,2*ngeo+1), at)
        g_piezo_tensor(2,5) = -g_piezo_tensor(1,4)
        CALL piezo_ij_bp(2, 4, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                tot_b_phase(1,2*ngeo+1), at)
        CALL piezo_ij_bp(1, 5, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                                tot_b_phase(1,3*ngeo+1), at)

     CASE(8)
!
!  D_2 (222) Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   d25   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   d36 )")') 

        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        CALL piezo_ij_bp(2, 5, ngeo, epsil_geo(1,1,ngeo+1), &
                                                   tot_b_phase(1,ngeo+1), at)
        CALL piezo_ij_bp(3, 6, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                 tot_b_phase(1,2*ngeo+1), at)

      CASE(9)
!
! D_3  Trigonal 
!
!         WRITE(stdout,'(5x,"( d11 -d11   .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -d14 2d11 )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezo_ij_bp(1, 1, ngeo, epsil_geo, tot_b_phase, at)
        g_piezo_tensor(1,2) = -g_piezo_tensor(1,1)
        g_piezo_tensor(2,6) = 2.0_DP * g_piezo_tensor(1,1)
        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo(1,1,ngeo+1), &
                                            tot_b_phase(1,ngeo+1), at)
        g_piezo_tensor(2,5) = -g_piezo_tensor(1,4)

     CASE(10,11)
!
! D_4  tetragonal, D_6 hexagonal
!
!         WRITE(stdout,'(/,5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -d14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        g_piezo_tensor(2,5) = -g_piezo_tensor(1,4)

     CASE(12)
!
! C_2v  Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   d24   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d32  d33   .    .    .  )")') 

        CALL piezo_ij_bp(3, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        CALL piezo_ij_bp(3, 2, ngeo, epsil_geo(1,1,ngeo+1), &
                                                tot_b_phase(1,ngeo+1), at)
        CALL piezo_ij_bp(3, 3, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                tot_b_phase(1,2*ngeo+1), at)
        CALL piezo_ij_bp(2, 4, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                                tot_b_phase(1,3*ngeo+1), at)
        CALL piezo_ij_bp(1, 5, ngeo, epsil_geo(1,1,4*ngeo+1), &
                                                tot_b_phase(1,4*ngeo+1), at)

     CASE(13)
!
! C_3v  Trigonal. Assuming m perpendicular to x1
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15 -d21 )")') 
!         WRITE(stdout,'(5x,"( d21 -d21   .   d15   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 
!
        CALL piezo_ij_bp(2, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(2,2) = -e_piezo_tensor(2,1)
        e_piezo_tensor(1,6) = e_piezo_tensor(2,2)
        CALL piezo_ij_bp(3, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(3,2) = e_piezo_tensor(3,1)
        CALL piezo_ij_bp(3, 3, ngeo, epsil_geo(1,1,ngeo+1), &
                                                  tot_b_phase(1,ngeo+1), at)
        CALL piezo_ij_bp(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                  tot_b_phase(1,2*ngeo+1), at)
        e_piezo_tensor(2,4) = e_piezo_tensor(1,5)

     CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   d15   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 

        CALL piezo_ij_bp(3, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(3,1) = e_piezo_tensor(3,1) * 0.5_DP
        e_piezo_tensor(3,2) = e_piezo_tensor(3,1)
        CALL piezo_ij_bp(3, 3, ngeo, epsil_geo(1,1,ngeo+1), &
                                                  tot_b_phase(1,ngeo+1), at)
        CALL piezo_ij_bp(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                  tot_b_phase(1,2*ngeo+1), at)
        e_piezo_tensor(2,4) = e_piezo_tensor(1,5)

     CASE(17)
!
! C_3h hexagonal
!
!             WRITE(stdout,'(5x,"( d11 -d11   .    .    .  -d12 )")') 
!             WRITE(stdout,'(5x,"( d12 -d12   .    .    .   d11 )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        CALL piezo_ij_bp(1, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(1,2) = -e_piezo_tensor(1,1)
        e_piezo_tensor(2,6) = e_piezo_tensor(1,1)
        CALL piezo_ij_bp(2, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(2,2) = -e_piezo_tensor(2,1)
        e_piezo_tensor(1,6) = -e_piezo_tensor(2,1)

      CASE(21)
!
! D_3h hexagonal
!
!             WRITE(stdout,'(5x,"(  .    .    .    .    .  -d21 )")') 
!             WRITE(stdout,'(5x,"( d21 -d21   .    .    .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
        CALL piezo_ij_bp(2, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(2,2) = -e_piezo_tensor(2,1)
        e_piezo_tensor(1,6) = -e_piezo_tensor(1,2)

     CASE(24)
!
! D_2d tetragonal: axis 2 || x1
!
!         WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   d14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   d36 )")') 

        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(2,5) = e_piezo_tensor(1,4)
        CALL piezo_ij_bp(3, 6, ngeo, epsil_geo(1,1,ngeo+1), &
                           tot_b_phase(1,ngeo+1), at)

     CASE(26)
!
! S_4 tetragonal
!
!        WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .  -d15  d14   .  )")') 
!        WRITE(stdout,'(5x,"( d31 -d31   .    .    .   d36 )")') 

        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo(1,1,ngeo+1), &
                                               tot_b_phase(1,ngeo+1), at)
        e_piezo_tensor(2,5) = e_piezo_tensor(1,4)
        CALL piezo_ij_bp(3, 1, ngeo, epsil_geo(1,1,1), tot_b_phase(1,1), at)
        e_piezo_tensor(3,2) = -e_piezo_tensor(3,1)
        CALL piezo_ij_bp(1, 5, ngeo, epsil_geo(1,1,2*ngeo+1), &
                                                   tot_b_phase(1,2*ngeo+1), at)
        e_piezo_tensor(2,4) = -e_piezo_tensor(1,5)
        CALL piezo_ij_bp(3, 6, ngeo, epsil_geo(1,1,3*ngeo+1), &
                                             tot_b_phase(1,3*ngeo+1), at)

     CASE(28,30)
!
! T, T_d cubic
!
!             WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .   d14   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .   d14 )")') 
!
        CALL piezo_ij_bp(1, 4, ngeo, epsil_geo, tot_b_phase, at)
        e_piezo_tensor(2,5) = e_piezo_tensor(1,4)
        e_piezo_tensor(3,6) = e_piezo_tensor(1,4)
     CASE(31)

     CASE DEFAULT
!
!  C_1 
!
        DO mn=1,6
           ind = ngeo * (mn-1)
           DO alpha=1,3
              CALL piezo_ij_bp(alpha, mn, ngeo, epsil_geo(1,1,ind+1), &
                                             tot_b_phase(1,ind+1), at)
           ENDDO
        ENDDO
  END SELECT
ELSE
   DO mn=1,6
      ind = ngeo * (mn-1)
      DO alpha=1,3
         CALL piezo_ij_bp(alpha, mn, ngeo, epsil_geo(1,1,ind+1), &
                                        tot_b_phase(1,ind+1), at)
      ENDDO
   ENDDO
ENDIF
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE compute_proper_piezo_tensor

!---------------------------------------------------------------------------
SUBROUTINE clean_piezo_tensor(piezo, ibrav, code_group)
!---------------------------------------------------------------------------
!
! This routine receives a piezoelectric tensor and a point group
! and sets to zero all components that have not been calculated.
!
!
IMPLICIT NONE
REAL(DP), INTENT(INOUT) :: piezo(3,6)
INTEGER, INTENT(IN) ::  ibrav, code_group
INTEGER :: i, j, igeo, alpha, ind, mn
LOGICAL :: check_group_ibrav

REAL(DP) :: epiezo(3,6)

epiezo=0.0_DP
IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group) 
     CASE(2,16,18,19,20,22,23,25,27,29,32) 
     CASE(3)
!
!  C_s   Monoclinic
!
!        WRITE(stdout,'(5x,"( g11  g12  g13   .   d15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .   g24   .   d26 )")') 
!        WRITE(stdout,'(5x,"( g31  g32  g33   .   d35   .  )")') 
!
        epiezo(1,1)=piezo(1,1)
        epiezo(1,2)=piezo(1,2)
        epiezo(1,3)=piezo(1,3)
        epiezo(2,6)=piezo(2,6)
        IF (ibrav==-12) THEN

           epiezo(3,1)=piezo(3,1)
           epiezo(3,2)=piezo(3,2)
           epiezo(3,3)=piezo(3,3)
           epiezo(2,4)=piezo(2,4)
           epiezo(1,5)=piezo(1,5)
           epiezo(3,5)=piezo(3,5)

        ELSE
!                WRITE(stdout,'(5x,"( d11  d12  d13   .    .   d16 )")') 
!                WRITE(stdout,'(5x,"( d21  d22  d23   .    .   d26 )")') 
!                WRITE(stdout,'(5x,"(  .    .    .   d16  d26   .  )")') 
!
           epiezo(2,1)=piezo(2,1)
           epiezo(2,2)=piezo(2,2)
           epiezo(2,3)=piezo(2,3)
           epiezo(1,6)=piezo(1,6)

        ENDIF
     CASE(4)
!
!  C_2   Monoclinic
!
        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)
        epiezo(3,6)=piezo(3,6)

        IF (ibrav==-12) THEN
!            WRITE(stdout,'(5x,"(  .    .    .   d14   .   d16 )")') 
!            WRITE(stdout,'(5x,"( d21  d22  d23   .   d25   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   d34   .   d36 )")') 

           epiezo(2,1)=piezo(2,1)
           epiezo(2,2)=piezo(2,2)
           epiezo(2,3)=piezo(2,3)
           epiezo(1,6)=piezo(1,6)
           epiezo(3,4)=piezo(3,4)
        ELSE
!            WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!            WRITE(stdout,'(5x,"(  .    .    .   d24  d25   .  )")') 
!            WRITE(stdout,'(5x,"( d31  d32  d33   .    .   d36 )")') 
!
           epiezo(3,1)=piezo(3,1)
           epiezo(3,2)=piezo(3,2)
           epiezo(3,3)=piezo(3,3)
           epiezo(1,5)=piezo(1,5)
           epiezo(2,4)=piezo(2,4)
         ENDIF

      CASE(6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!
!             WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .   d24 -d14   .  )")') 
!             WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 
!
        epiezo(3,1)=piezo(3,1)
        epiezo(3,2)=piezo(3,2)
        epiezo(3,3)=piezo(3,3)
        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)
        epiezo(2,4)=piezo(2,4)
        epiezo(1,5)=piezo(1,5)

     CASE(8)
!
!  D_2 (222) Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   d25   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   d36 )")') 

        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)
        epiezo(3,6)=piezo(3,6)

      CASE(9)
!
! D_3  Trigonal 
!
!         WRITE(stdout,'(5x,"( d11 -d11   .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -d14 2d11 )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        epiezo(1,1)=piezo(1,1)
        epiezo(1,2)=piezo(1,2)
        epiezo(2,6)=piezo(2,6)
        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)

     CASE(10,11)
!
! D_4  tetragonal, D_6 hexagonal
!
!         WRITE(stdout,'(/,5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .  -d14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)

     CASE(12)
!
! C_2v  Orthorombic
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   d24   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d32  d33   .    .    .  )")') 

        epiezo(3,1)=piezo(3,1)
        epiezo(3,2)=piezo(3,2)
        epiezo(3,3)=piezo(3,3)
        epiezo(2,4)=piezo(2,4)
        epiezo(1,5)=piezo(1,5)

     CASE(13)
!
! C_3v  Trigonal. Assuming m perpendicular to x1
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15 -d21 )")') 
!         WRITE(stdout,'(5x,"( d21 -d21   .   d15   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 
!
        epiezo(2,1)=piezo(2,1)
        epiezo(2,2)=piezo(2,2)
        epiezo(1,6)=piezo(1,6)
        epiezo(3,1)=piezo(3,1)
        epiezo(3,2)=piezo(3,2)
        epiezo(3,3)=piezo(3,3)
        epiezo(1,5)=piezo(1,5)
        epiezo(2,4)=piezo(2,4)

     CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   d15   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 

        epiezo(3,1)=piezo(3,1)
        epiezo(3,2)=piezo(3,2)
        epiezo(3,3)=piezo(3,3)
        epiezo(1,5)=piezo(1,5)
        epiezo(2,4)=piezo(2,4)

     CASE(17)
!
! C_3h hexagonal
!
!             WRITE(stdout,'(5x,"( d11 -d11   .    .    .  -d12 )")') 
!             WRITE(stdout,'(5x,"( d12 -d12   .    .    .   d11 )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 

        epiezo(1,1)=piezo(1,1)
        epiezo(1,2)=piezo(1,2)
        epiezo(2,6)=piezo(2,6)
        epiezo(2,1)=piezo(2,1)
        epiezo(2,2)=piezo(2,2)
        epiezo(1,6)=piezo(1,6)

      CASE(21)
!
! D_3h hexagonal
!
!             WRITE(stdout,'(5x,"(  .    .    .    .    .  -d21 )")') 
!             WRITE(stdout,'(5x,"( d21 -d21   .    .    .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
        epiezo(2,1)=piezo(2,1)
        epiezo(2,2)=piezo(2,2)
        epiezo(1,6)=piezo(1,6)

     CASE(24)
!
! D_2d tetragonal: axis 2 || x1
!
!         WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .   d14   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .    .    .   d36 )")') 

        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)
        epiezo(3,6)=piezo(3,6)

     CASE(26)
!
! S_4 tetragonal
!
!        WRITE(stdout,'(5x,"(  .    .    .   d14  d15   .  )")') 
!        WRITE(stdout,'(5x,"(  .    .    .  -d15  d14   .  )")') 
!        WRITE(stdout,'(5x,"( d31 -d31   .    .    .   d36 )")') 

        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)
        epiezo(3,1)=piezo(3,1)
        epiezo(3,2)=piezo(3,2)
        epiezo(1,5)=piezo(1,5)
        epiezo(2,4)=piezo(2,4)
        epiezo(3,6)=piezo(3,6)

     CASE(28,30)
!
! T, T_d cubic
!
!             WRITE(stdout,'(5x,"(  .    .    .   d14   .    .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .   d14   .  )")') 
!             WRITE(stdout,'(5x,"(  .    .    .    .    .   d14 )")') 
!
        epiezo(1,4)=piezo(1,4)
        epiezo(2,5)=piezo(2,5)
        epiezo(3,6)=piezo(3,6)
     CASE(31)

     CASE DEFAULT
!
!  C_1 
!
      epiezo(:,:)=piezo(:,:)
  END SELECT
ELSE
   epiezo(:,:)=piezo(:,:)
ENDIF
piezo(:,:)=epiezo(:,:)

RETURN
END SUBROUTINE clean_piezo_tensor

!---------------------------------------------------------------------------
SUBROUTINE piezo_ij(ialpha, mn, ngeo, epsil_geo, polar_geo)
!---------------------------------------------------------------------------
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit
USE voigt, ONLY : voigt_extract_indices

IMPLICIT NONE
INTEGER, INTENT(IN) :: mn, ialpha, ngeo
REAL(DP), INTENT(IN) :: epsil_geo(3,3,ngeo), polar_geo(3,ngeo)
INTEGER :: igeo, m, n, mnin
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo)

mnin=mn
CALL voigt_extract_indices(m,n,mnin)
WRITE(stdout,'("Piezoelectric_tensor(",i1,",",i1,"): &
                    &strain,    polarization")') ialpha, mn
DO igeo=1,ngeo
   x(igeo)=epsil_geo(m,n,igeo)
   y(igeo)=polar_geo(ialpha,igeo)
   WRITE(stdout,'(20x,f15.10,3x,f15.10)') x(igeo), y(igeo)
ENDDO
CALL polyfit( x, y, ngeo, alpha, m1-1 )
g_piezo_tensor(ialpha, mn) = alpha(2)
!
!  The stress piezoelectric tensor is defined from 
!  P_i = \sum_{jk} e_{ijk} \epsilon_{jk}
!  therefore when we apply an off-diagonal strain we obtain a
!  polarization
!  P_i = e_{ijk} \epsilon_{jk} + e_{ikj} \epsilon_{jk} 
!      = (e_{ijk}+e_{ikj}) \epsilon_{jk}
!  So dP_i / \epsilon_{jk} is twice the piezoelectric tensor 
!  e_{ijk} and we have to divide by 2 to get e_{ijk}.
!
!  Similarly, in Voigt notation we have
!  e_{im} = e_{ijk} while \epsilon_{m} = 2 \epsilon_{jk} 
!  for the off diagonal components and \epsilon_{m} = \epsilon{jk} for
!  the diagonal ones so we have still
!  P_i= e_{im} \epsilon_{m} 
!  We have e_{im} = d P_i \over d \epsilon_{m} and for the off diagonal
!  components  d P_i \over 2 d\epsilon_{jk} and since we have computed  
!  d P_i \over d\epsilon_{jk} above we still need to divide by 2.
!
IF (m /= n) g_piezo_tensor(ialpha, mn) = g_piezo_tensor(ialpha, mn) * 0.5_DP

RETURN
END SUBROUTINE piezo_ij

!---------------------------------------------------------------------------
SUBROUTINE piezo_ij_bp(ialpha, mn, ngeo, epsil_geo, tot_b_phase, at)
!---------------------------------------------------------------------------
!
! On output the proper piezoelectric coefficient must be multiplied by 
! alat and divided by the volume of the umperturbed unit cell to have 
! it in e/bohr**2.
!
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit
USE voigt, ONLY : voigt_extract_indices

IMPLICIT NONE
INTEGER, INTENT(IN) :: mn, ialpha, ngeo
REAL(DP), INTENT(IN) :: epsil_geo(3,3,ngeo), tot_b_phase(3,ngeo), at(3,3)
INTEGER :: igeo, m, n, mnin, ipol
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo)

REAL(DP) :: polar_at(3)

mnin=mn
CALL voigt_extract_indices(m,n,mnin)
DO ipol=1,3 ! loop on the three components in the crystal basis
   WRITE(stdout,'("Piezoelectric_tensor(",i1,",",i1,"): &
             &strain,        Berry phase, direction= ",i1)') ialpha, mn, ipol
   DO igeo=1,ngeo
      x(igeo)=epsil_geo(m,n,igeo)
      y(igeo)=tot_b_phase(ipol,igeo)
      WRITE(stdout,'(20x,f15.10,3x,f15.10)') x(igeo), y(igeo)
   ENDDO
   CALL polyfit( x, y, ngeo, alpha, m1-1 )
   polar_at(ipol)=alpha(2)
ENDDO
e_piezo_tensor(ialpha,mn) = polar_at(1)*at(ialpha,1) + &
                            polar_at(2)*at(ialpha,2) + &
                            polar_at(3)*at(ialpha,3)
!
!  The piezoelectric tensor relates the polarization to the strain in voigt
!  notation. Since e_23 = 0.5 e_4, e_13 = 0.5 e_5, e_12 = 0.5 e_6 we have
!  to divide by 2 the elements of the piezoelectric tensor calculated with 
!  off diagonal strain components
!
IF (m /= n) e_piezo_tensor(ialpha, mn) = e_piezo_tensor(ialpha, mn) * 0.5_DP
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE piezo_ij_bp
!
!----------------------------------------------------------------------------
SUBROUTINE compute_d_piezo_tensor(smn)
!----------------------------------------------------------------------------
!
! This routine computes the direct piezoelectric tensor that links the
! induced polarization to the stress
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: smn(6,6)
INTEGER :: alpha, m, n

d_piezo_tensor=0.0_DP
DO alpha=1,3
   DO m=1,6
      DO n=1,6
         d_piezo_tensor(alpha,m) = d_piezo_tensor(alpha,m) + &
                                   e_piezo_tensor(alpha,n) * smn(n, m)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE compute_d_piezo_tensor

!----------------------------------------------------------------------------
SUBROUTINE print_piezo_info(code_group,ibrav,ngeo_strain)
!----------------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, ibrav, ngeo_strain

INTEGER :: nstrain

SELECT CASE (code_group) 
   CASE (2,16,18,19,20,22,23,25,27,29,32) 
      nstrain=0
   CASE (3)
!
!  C_s   Monoclinic
!
      WRITE(stdout,'(/,5x,"It requires five strains: e1, e2, e3, e4, and e5")')
      nstrain=5
   CASE (4)
!
!  C_2   Monoclinic
!
      WRITE(stdout,'(/,5x,"It requires all six strains")')
      nstrain=6
   CASE (6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!
      WRITE(stdout,'(/,5x,"It requires four strains: e1, e3, e4, and e5")')
      nstrain=4
   CASE (8)
!
!  D_2 (222) Orthorombic
!
      WRITE(stdout,'(/,5x,"It requires two strains: e4, e5, and e6")')
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
! C_3v  Trigonal. Assuming m perpendicular to x1
! C_4v tetragonal, C_6v hexagonal
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
      nstrain=1
   CASE (31)
      nstrain=0
   CASE DEFAULT
!
!  C_1 
!
      WRITE(stdout,'(/,5x,"It requires all six strains")')
      nstrain=6
   END SELECT

   IF(nstrain>0) THEN
     WRITE(stdout,'(5x,"for a total of",i3,&
                                  &" scf calculations")') nstrain*ngeo_strain
     WRITE(stdout,'(5x,"and " i4" Berry phase calculations")') &
                                                        nstrain*ngeo_strain*3
   ENDIF
RETURN
END SUBROUTINE print_piezo_info
!
!-----------------------------------------------------------------------------
SUBROUTINE compute_polarization_equil(polar_geo, epsil_geo, polar0, ngeo, &
                                                               ngeo_strain)
!-----------------------------------------------------------------------------

USE constants, ONLY : electron_si, bohr_radius_si
USE polyfit_mod, ONLY : polyfit

IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeo, ngeo_strain
REAL(DP), INTENT(IN) :: epsil_geo(ngeo), polar_geo(3,ngeo)
REAL(DP), INTENT(OUT) :: polar0(3)

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
         y(igeo,ipol)=polar_geo(ipol,count_igeo)
      ENDDO
   ENDDO
   DO ipol=1,3
      WRITE(stdout,'("Strain type = ", i2," Polarization(",i1,")")') istep, &
                                                                     ipol
      DO igeo=1,ngeo_strain
         WRITE(stdout,'(2f15.10)') x(igeo), y(igeo,ipol)
      ENDDO
      CALL polyfit( x, y(1,ipol), ngeo_strain, alpha, m1-1 )
      polar_at(ipol)=alpha(1)
   ENDDO
!
!  save the polarization to export to the calling routine
!
   polar0(:)=polar_at(:)
!
!  and print it
!
   WRITE(stdout,'(/,"Spontaneous polarization of the equilibrium structure &
                             &strain type: ",i4)') istep
   WRITE(stdout,'(5x,"P=(", 2(f10.5,","), f10.5, "   ) e/(a.u.)^2")') &
                                                             polar_at(:)
   fact= electron_si / (bohr_radius_si)**2
   WRITE(stdout,'(5x,"P=(", 2(f10.5,","), f10.5, "   ) C/m^2",/)') &
                                                         polar_at(:) * fact
ENDDO

WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE compute_polarization_equil

!--------------------------------------------------------------------------
SUBROUTINE proper_improper_piezo(polar, piezo_in, piezo_out, dir)
!--------------------------------------------------------------------------
!
! This routine receives the piezoelectric tensor piezo_in which
! is in the proper or improper form and transforms it in the other
! form using the formula
! \tilde e_ijk = e_ijk + \delta_jk P_i - \delta_ij P_k
! or the inverse
! e_ijk = \tilde e_ijk - \delta_jk P_i + \delta_ij P_k
! where P is the polarization, \tilde e_ijk is the proper piezoelectric
! tensor and e_ijk is the improper one. Piezo_in is the proper piezoelectric
! tensor if dir=1 or the improper one if dir=-1. piezo_out is the
! other tensor. Both tensors in input and output are 3x6 matrices in 
! Voigt notation.
! 
!
USE kinds, ONLY : DP
USE voigt, ONLY : to_voigt3
IMPLICIT NONE

INTEGER, INTENT(IN) :: dir
REAL(DP), INTENT(INOUT)  :: polar(3), piezo_in(3,6)
REAL(DP), INTENT(OUT) ::  piezo_out(3,6)

REAL(DP) :: piezo_in_aux(3,3,3), piezo_voigt_aux(3,6)
INTEGER :: ipol, jpol, kpol, i, m


IF (dir==1) THEN
!
!  from proper to improper
!
   piezo_voigt_aux(:,:)=piezo_in(:,:)
   DO i=1,3
      DO m=3,6
         piezo_voigt_aux(i,m)=piezo_in(i,m)*2.0_DP
      ENDDO
   ENDDO
   CALL to_voigt3(piezo_voigt_aux, piezo_in_aux, 1.0_DP, .FALSE.)
   DO ipol=1,3
      DO jpol=1,3
         DO kpol=1,3
            IF (jpol==kpol) &
               piezo_in_aux(ipol,jpol,kpol) = piezo_in_aux(ipol,jpol,kpol) &
                                                             - polar(ipol)
            IF (ipol==jpol) &
               piezo_in_aux(ipol,jpol,kpol) = piezo_in_aux(ipol,jpol,kpol) &
                                                             + polar(kpol)
         ENDDO
      ENDDO
   ENDDO
   CALL to_voigt3(piezo_out, piezo_in_aux, 1.0_DP, .TRUE.)
   DO i=1,3
      DO m=3,6
         piezo_out(i,m)=piezo_out(i,m)*0.5_DP
      ENDDO
   ENDDO
ELSEIF (dir==-1) THEN
!
!  from improper to proper
!
   piezo_voigt_aux(:,:)=piezo_in(:,:)
   DO i=1,3
      DO m=3,6
         piezo_voigt_aux(i,m)=piezo_in(i,m)*2.0_DP
      ENDDO
   ENDDO
   CALL to_voigt3(piezo_voigt_aux, piezo_in_aux, 1.0_DP, .FALSE.)
   DO ipol=1,3
      DO jpol=1,3
         DO kpol=1,3
            IF (jpol==kpol) &
               piezo_in_aux(ipol,jpol,kpol) = piezo_in_aux(ipol,jpol,kpol) &
                                                               + polar(ipol)
            IF (ipol==jpol) &
               piezo_in_aux(ipol,jpol,kpol) = piezo_in_aux(ipol,jpol,kpol) &
                                                               - polar(kpol)
         ENDDO
      ENDDO
   ENDDO
   CALL to_voigt3(piezo_out, piezo_in_aux, 1.0_DP, .TRUE.)
   DO i=1,3
      DO m=3,6
         piezo_out(i,m)=piezo_out(i,m)*0.5_DP
      ENDDO
   ENDDO
ELSE
   CALL errore('proper_improper_piezo','wrong dir',1)
ENDIF

RETURN
END SUBROUTINE proper_improper_piezo
!
!-------------------------------------------------------------------------
SUBROUTINE write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
                                  e_piezo_tensor_t, filename, iflag, iwhat)
!-------------------------------------------------------------------------
!
!  iflag=0 writes the piezoelectric tensor as a function of temperature
!  iflag=2 writes the piezoelectric tensor as a function of pressure
!  iwhat=1 writes e_piezo_tensor
!  iwhat=2 writes d_piezo_tensor
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, code_group, iflag, iwhat
REAL(DP), INTENT(IN) :: temp(ntemp), e_piezo_tensor_t(3,6,ntemp)
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
   CALL errore('write_piezo_on_file','opening piezo tensor (T) file',&
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

   IF (iwhat==1) THEN
      WRITE(iu_piezo,'("#    stress piezoelectric tensor in e/bohr^2 units")')
      WRITE(iu_piezo,'("#    multiply by 0.57214766E+02 to have it in C/m^2")')
   ELSE
      WRITE(iu_piezo,'("#    strain piezoelectric tensor in e/bohr^2/kbar &
                                                                     &units")')
      WRITE(iu_piezo,'("#    multiply by 0.57214766E+06 to have it in pC/N")')
   ENDIF
   SELECT CASE (code_group)
      CASE(30,28)
!
!   T_d, T
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_14 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_14 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,4,itemp)
         ENDDO
      CASE(26)
!
!   S_4
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_31 ", &
                   &13x, " e_14 ", 13x, " e_15 ", 13x, " e_36 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_31 ", &
                   &13x, " d_14 ", 13x, " d_15 ", 13x," d_36 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(1,4,itemp),                &
                e_piezo_tensor_t(1,5,itemp),                &
                e_piezo_tensor_t(3,6,itemp)
         ENDDO
      CASE(24)
!
!   D_2d
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_14 ", &
                   &13x, " e_36 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_14 ", &
                   &13x, " d_36 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,4,itemp),                &  
                e_piezo_tensor_t(3,6,itemp)                
         ENDDO
      CASE(21)
!
!   D_3h
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_11 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_11 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,1,itemp)                
         ENDDO
      CASE(17)
!
!   C_3h
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_11 ", &
                   &13x, " e_22 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_11 ", &
                   &13x, " d_22 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,1,itemp),                &  
                e_piezo_tensor_t(2,2,itemp)                
         ENDDO

      CASE(14,15)
!
!   C_4v and C_6v
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_31 ", &
                  & 13x, "     e_33 ", 13x, "     e_15 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_31 ", &
                  & 13x, "     d_33 ", 13x, "     d_15 ")') label
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(3,1,itemp), e_piezo_tensor_t(3,3,itemp), &
                e_piezo_tensor_t(1,5,itemp)
         ENDDO
      CASE(13)
!
!   C_3v
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_22 ", &
                   &13x, " e_31 ", 13x, " e_33 ", 13x, " e_15 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_22 ", &
                   &13x, " d_31 ", 13x, " d_33 ", 13x, " d_15 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(2,2,itemp),                &  
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,3,itemp),                &  
                e_piezo_tensor_t(1,5,itemp)                
         ENDDO
      CASE(12)
!
!   C_2v
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_31 ", &
                   &13x, " e_32 ", 13x, " e_33 ", 13x, " e_15 ",&
                   &13x, " e_24 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_31 ", &
                   &13x, " d_32 ", 13x, " d_33 ", 13x, " d_15 ",&
                   &13x, " d_24 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,2,itemp),                &  
                e_piezo_tensor_t(3,3,itemp),                &  
                e_piezo_tensor_t(1,5,itemp),                &  
                e_piezo_tensor_t(2,4,itemp)                
         ENDDO
      CASE(11,10)
!
!   D_6, D_4
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_14 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_14 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,4,itemp)                
         ENDDO

      CASE(9)
!
!   D_3
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_11 ", 13x, " e_14 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_11 ", 13x, " d_14 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,1,itemp),                &  
                e_piezo_tensor_t(1,4,itemp)
         ENDDO
      CASE(8)
!
!   D_2
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_14 ", &
                   &13x, " e_25 ", 13x, " e_36 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_14 ", &
                   &13x, " d_25 ", 13x, " d_36 ")') label
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,4,itemp),                &  
                e_piezo_tensor_t(2,5,itemp),                &  
                e_piezo_tensor_t(3,6,itemp)  
         ENDDO
      CASE(6,7)
!
!  C_6  C_4
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_31 ", &
                   &13x, " e_33 ", 13x, " e_14 ", 13x, " e_15 ",&
                   &13x, " e_24 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_31 ", &
                   &13x, " d_33 ", 13x, " d_14 ", 13x, " d_15 ",&
                   &13x, " d_24 ")') label
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,3,itemp),                &  
                e_piezo_tensor_t(1,4,itemp),                &  
                e_piezo_tensor_t(1,5,itemp),                & 
                e_piezo_tensor_t(2,4,itemp)  
         ENDDO
      CASE(5)
!
!   C_3
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_11 ", &
                   &13x, " e_22 ", 13x, " e_31 ", 13x, " e_33 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_11 ", &
                   &13x, " d_22 ", 13x, " d_31 ", 13x, " d_33 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,1,itemp),                &  
                e_piezo_tensor_t(2,2,itemp),                &  
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,3,itemp)  
         ENDDO
      CASE(4)
!
!   C_2
!
         IF (ibrav==-12) THEN
            IF (iwhat==1) THEN
               WRITE(iu_piezo,'("#",5x, a7, 13x, " e_21 ", &
                   &13x, " e_22 ", 13x, " e_23 ", 13x, " e_14 ", &
                   &13x, " e_16 ", 13x, " e_25 ", 13x, " e_34 ", &
                   &13x, " e_36 ")') label
            ELSE
               WRITE(iu_piezo,'("#",5x, a7, 13x, " d_21 ", &
                   &13x, " d_22 ", 13x, " d_23 ", 13x, " d_14 ", &
                   &13x, " d_16 ", 13x, " d_25 ", 13x, " d_34 ", &
                   &13x, " d_36 ")') label

            ENDIF

            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                   e_piezo_tensor_t(2,1,itemp),                &  
                   e_piezo_tensor_t(2,2,itemp),                &  
                   e_piezo_tensor_t(2,3,itemp),                &  
                   e_piezo_tensor_t(1,4,itemp),                &  
                   e_piezo_tensor_t(1,6,itemp),                &  
                   e_piezo_tensor_t(2,5,itemp),                &  
                   e_piezo_tensor_t(3,4,itemp),                &  
                   e_piezo_tensor_t(3,6,itemp)  
            ENDDO
         ELSE
            IF (iwhat==1) THEN
               WRITE(iu_piezo,'("#",5x, a7, 13x, " e_21 ", &
                   &13x, " e_22 ", 13x, " e_23 ", 13x, " e_14 ", &
                   &13x, " e_15 ", 13x, " e_25 ", 13x, " e_24 ", &
                   &13x, " e_36 ")') label
            ELSE
               WRITE(iu_piezo,'("#",5x, a7, 13x, " d_21 ", &
                   &13x, " d_22 ", 13x, " d_23 ", 13x, " d_14 ", &
                   &13x, " d_15 ", 13x, " d_25 ", 13x, " d_24 ", &
                   &13x, " d_36 ")') label

            ENDIF
            DO itemp=2,ntemp-1
               WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                   e_piezo_tensor_t(2,1,itemp),                &  
                   e_piezo_tensor_t(2,2,itemp),                &  
                   e_piezo_tensor_t(2,3,itemp),                &  
                   e_piezo_tensor_t(1,4,itemp),                &  
                   e_piezo_tensor_t(1,5,itemp),                &  
                   e_piezo_tensor_t(2,5,itemp),                &  
                   e_piezo_tensor_t(2,4,itemp),                &  
                   e_piezo_tensor_t(3,6,itemp)  
            ENDDO
         ENDIF
      CASE(3)
!
!   C_s
!
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_11 ",             &
                   &13x, " e_12 ", 13x, " e_13 ", 13x, " e_15",  &
                   &13x, " e_24 ", 13x, " e_26 ", 13x, " e_31 ", &
                   &13x, " e_32 ", 13x, " e_33 ", 13x, " e_35 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_11 ",             &
                   &13x, " d_12 ", 13x, " d_13 ", 13x, " d_15",  &
                   &13x, " d_24 ", 13x, " d_26 ", 13x, " d_31 ", &
                   &13x, " d_32 ", 13x, " d_33 ", 13x, " d_35 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
                e_piezo_tensor_t(1,1,itemp),                &  
                e_piezo_tensor_t(1,2,itemp),                &  
                e_piezo_tensor_t(1,3,itemp),                &  
                e_piezo_tensor_t(1,5,itemp),                &  
                e_piezo_tensor_t(2,4,itemp),                &  
                e_piezo_tensor_t(2,6,itemp),                &  
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,2,itemp),                &
                e_piezo_tensor_t(3,3,itemp),                &
                e_piezo_tensor_t(3,5,itemp)  
         ENDDO

      CASE DEFAULT
         IF (iwhat==1) THEN
            WRITE(iu_piezo,'("#",5x, a7, 13x, " e_11 ",             &
                   &13x, " e_12 ", 13x, " e_13 ", 13x, " e_14",  &
                   &13x, " e_15 ", 13x, " e_16 ", 13x, " e_21 ", &
                   &13x, " e_22 ", 13x, " e_23 ", 13x, " e_24 ", &
                   &13x, " e_25 ", 13x, " e_26 ", 13x, " e_31 ", &
                   &13x, " e_32 ", 13x, " e_33 ", 13x, " e_34 ", &
                   &13x, " e_35 ", 13x, " e_36 ")') label
         ELSE
            WRITE(iu_piezo,'("#",5x, a7, 13x, " d_11 ",             &
                   &13x, " d_12 ", 13x, " d_13 ", 13x, " d_14",  &
                   &13x, " d_15 ", 13x, " d_16 ", 13x, " d_21 ", &
                   &13x, " d_22 ", 13x, " d_23 ", 13x, " d_24 ", &
                   &13x, " d_25 ", 13x, " d_26 ", 13x, " d_31 ", &
                   &13x, " d_32 ", 13x, " d_33 ", 13x, " d_34 ", &
                   &13x, " d_35 ", 13x, " d_36 ")') label
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,'(e16.8,21e20.12)') temp(itemp), &
                  e_piezo_tensor_t(1,1,itemp), e_piezo_tensor_t(1,2,itemp), &
                  e_piezo_tensor_t(1,3,itemp), e_piezo_tensor_t(1,4,itemp), &
                  e_piezo_tensor_t(1,5,itemp), e_piezo_tensor_t(1,6,itemp), &
                  e_piezo_tensor_t(2,1,itemp), e_piezo_tensor_t(2,2,itemp), &
                  e_piezo_tensor_t(2,3,itemp), e_piezo_tensor_t(2,4,itemp), &
                  e_piezo_tensor_t(2,5,itemp), e_piezo_tensor_t(2,6,itemp), &
                  e_piezo_tensor_t(3,1,itemp), e_piezo_tensor_t(3,2,itemp), &
                  e_piezo_tensor_t(3,3,itemp), e_piezo_tensor_t(3,4,itemp), &
                  e_piezo_tensor_t(3,5,itemp), e_piezo_tensor_t(3,6,itemp)
         ENDDO
   END SELECT
   CLOSE(iu_piezo)
ENDIF

RETURN
END SUBROUTINE write_piezo_tensor_on_file
!
!-------------------------------------------------------------------------
SUBROUTINE read_piezo_tensor_from_file(temp, ntemp, ibrav, code_group, &
                                  e_piezo_tensor_t, filename)
!-------------------------------------------------------------------------
!
! This routine reads a file with the temperature (or pressure) and
! the inequivalent components of the piezoelectric tensor.
! It fills the entire tensor.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, code_group
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: e_piezo_tensor_t(3,6,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_piezo, ios
INTEGER :: find_free_unit
REAL(DP) :: rdum
CHARACTER(LEN=7) :: label
!
! If the piezoelectric tensor vanishes return
!
SELECT CASE (code_group)
   CASE(2,16,18,19,20,22,23,25,27,29,31,32)
   RETURN
END SELECT

iu_piezo=find_free_unit()
e_piezo_tensor_t=0.0_DP
IF (meta_ionode) &
   OPEN(UNIT=iu_piezo, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('read_piezo_from_file','opening piezo tensor (T) file',&
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
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(26)
!
!   S_4
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(3,1,itemp),  &  
                e_piezo_tensor_t(1,4,itemp),  &
                e_piezo_tensor_t(1,5,itemp),  &
                e_piezo_tensor_t(3,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(24)
!
!   D_2d
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            WRITE(iu_piezo,*) rdum, e_piezo_tensor_t(1,4,itemp),   &  
                                           e_piezo_tensor_t(3,6,itemp) 
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(21)
!
!   D_3h
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,1,itemp)                
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(17)
!
!   C_3h
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,1,itemp),  &  
                e_piezo_tensor_t(2,2,itemp)                
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO

      CASE(14,15)
!
!   C_4v and C_6v
!
         READ(iu_piezo,*)
         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(3,1,itemp), &
                 e_piezo_tensor_t(3,3,itemp), e_piezo_tensor_t(1,5,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(13)
!
!   C_3v
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(2,2,itemp),  &  
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,3,itemp),                &  
                e_piezo_tensor_t(1,5,itemp)                
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(12)
!
!   C_2v
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(3,1,itemp),   &  
                e_piezo_tensor_t(3,2,itemp),                &  
                e_piezo_tensor_t(3,3,itemp),                &  
                e_piezo_tensor_t(1,5,itemp),                &
                e_piezo_tensor_t(2,4,itemp)                
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(11,10)
!
!   D_6, D_4
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO

      CASE(9)
!
!   D_3
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,1,itemp),  &  
                e_piezo_tensor_t(1,4,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(8)
!
!   D_2
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,4,itemp),  &  
                e_piezo_tensor_t(2,5,itemp),                &  
                e_piezo_tensor_t(3,6,itemp)  
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(6, 7)
!
!   C_6 C_4
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(3,1,itemp),     &  
                e_piezo_tensor_t(3,3,itemp),                &  
                e_piezo_tensor_t(1,4,itemp),                &  
                e_piezo_tensor_t(1,5,itemp),                &
                e_piezo_tensor_t(2,4,itemp)  
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(5)
!
!   C_3
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,1,itemp),  &  
                e_piezo_tensor_t(2,2,itemp),                &  
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,3,itemp)  
            IF (ABS(rdum-temp(itemp))>1D-5) &
               CALL errore('read_piezo_from_file','incorrect temperature', 1)
         ENDDO
      CASE(4)
!
!   C_2
!
         READ(iu_piezo,*)

         IF (ibrav==-12) THEN
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum,                       &
                e_piezo_tensor_t(2,1,itemp),                &  
                e_piezo_tensor_t(2,2,itemp),                &  
                e_piezo_tensor_t(2,3,itemp),                &  
                e_piezo_tensor_t(1,4,itemp),                &  
                e_piezo_tensor_t(1,6,itemp),                &  
                e_piezo_tensor_t(2,5,itemp),                &  
                e_piezo_tensor_t(3,4,itemp),                &  
                e_piezo_tensor_t(3,6,itemp)  
                IF (ABS(rdum-temp(itemp))>1D-5) &
                   CALL errore('read_piezo_from_file','incorrect temperature',1)
            ENDDO
         ELSE
            DO itemp=2,ntemp-1
               READ(iu_piezo,*) rdum,                       &
                e_piezo_tensor_t(2,1,itemp),                &  
                e_piezo_tensor_t(2,2,itemp),                &  
                e_piezo_tensor_t(2,3,itemp),                &  
                e_piezo_tensor_t(1,4,itemp),                &  
                e_piezo_tensor_t(1,5,itemp),                &  
                e_piezo_tensor_t(2,5,itemp),                &  
                e_piezo_tensor_t(2,4,itemp),                &  
                e_piezo_tensor_t(3,6,itemp)  
                IF (ABS(rdum-temp(itemp))>1D-5) &
                   CALL errore('read_piezo_from_file','incorrect temperature',1)
            ENDDO
         ENDIF

      CASE(3)
!
!   C_s
!
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, e_piezo_tensor_t(1,1,itemp),  &  
                e_piezo_tensor_t(1,2,itemp),                &  
                e_piezo_tensor_t(1,3,itemp),                &  
                e_piezo_tensor_t(1,5,itemp),                &  
                e_piezo_tensor_t(2,4,itemp),                &  
                e_piezo_tensor_t(2,6,itemp),                &  
                e_piezo_tensor_t(3,1,itemp),                &  
                e_piezo_tensor_t(3,2,itemp),                &
                e_piezo_tensor_t(3,3,itemp),                &
                e_piezo_tensor_t(3,5,itemp)  
                IF (ABS(rdum-temp(itemp))>1D-5) &
                   CALL errore('read_piezo_from_file','incorrect temperature',1)
         ENDDO

      CASE DEFAULT
         READ(iu_piezo,*)

         DO itemp=2,ntemp-1
            READ(iu_piezo,*) rdum, &
                  e_piezo_tensor_t(1,1,itemp), e_piezo_tensor_t(1,2,itemp), &
                  e_piezo_tensor_t(1,3,itemp), e_piezo_tensor_t(1,4,itemp), &
                  e_piezo_tensor_t(1,5,itemp), e_piezo_tensor_t(1,6,itemp), &
                  e_piezo_tensor_t(2,1,itemp), e_piezo_tensor_t(2,2,itemp), &
                  e_piezo_tensor_t(2,3,itemp), e_piezo_tensor_t(2,4,itemp), &
                  e_piezo_tensor_t(2,5,itemp), e_piezo_tensor_t(2,6,itemp), &
                  e_piezo_tensor_t(3,1,itemp), e_piezo_tensor_t(3,2,itemp), &
                  e_piezo_tensor_t(3,3,itemp), e_piezo_tensor_t(3,4,itemp), &
                  e_piezo_tensor_t(3,5,itemp), e_piezo_tensor_t(3,6,itemp)
            IF (ABS(rdum-temp(itemp))>1D-5) &
                CALL errore('read_piezo_from_file','incorrect temperature',1)
         ENDDO
   END SELECT
   CLOSE(iu_piezo)
   CALL expand_piezo_tensor(e_piezo_tensor_t, code_group, ibrav, ntemp, temp)
ENDIF
CALL mp_bcast(e_piezo_tensor_t, meta_ionode_id, world_comm) 

RETURN
END SUBROUTINE read_piezo_tensor_from_file

!-------------------------------------------------------------------------
SUBROUTINE expand_piezo_tensor(e_piezo_tensor_t, code_group, ibrav, &
                                                       ntemp, temp)
!-------------------------------------------------------------------------
!
! This routine reconstruct the complete piezoelectric tensor from the
! components read from file
!
USE kinds,      ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, code_group
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: e_piezo_tensor_t(3,6,ntemp)

INTEGER :: itemp
!
! If the piezoelectric tensor vanishes return
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
            e_piezo_tensor_t(2,5,itemp) = e_piezo_tensor_t(1,4,itemp)
            e_piezo_tensor_t(3,6,itemp) = e_piezo_tensor_t(1,4,itemp)
         ENDDO
    CASE(26)
!
!   S_4
!
         DO itemp=2,ntemp-1
            e_piezo_tensor_t(3,2,itemp)=-e_piezo_tensor_t(3,1,itemp)
            e_piezo_tensor_t(2,5,itemp)= e_piezo_tensor_t(1,4,itemp)
            e_piezo_tensor_t(2,4,itemp)=-e_piezo_tensor_t(1,5,itemp)
         ENDDO
   CASE(24)
!
!   D_2d
!
         DO itemp=2,ntemp-1
            e_piezo_tensor_t(2,5,itemp)= e_piezo_tensor_t(1,4,itemp)
         ENDDO
   CASE(21)
!
!   D_3h
!
         DO itemp=2,ntemp-1
            e_piezo_tensor_t(2,2,itemp)=-e_piezo_tensor_t(2,1,itemp)
            e_piezo_tensor_t(1,6,itemp)=-e_piezo_tensor_t(2,1,itemp)
         ENDDO
   CASE(17)
!
!   C_3h
!

         DO itemp=2,ntemp-1
            e_piezo_tensor_t(1,2,itemp)=-e_piezo_tensor_t(1,1,itemp)
            e_piezo_tensor_t(2,1,itemp)=-e_piezo_tensor_t(2,2,itemp)
            e_piezo_tensor_t(1,6,itemp)= e_piezo_tensor_t(2,2,itemp)
            e_piezo_tensor_t(2,6,itemp)= e_piezo_tensor_t(1,1,itemp)
         ENDDO
   CASE(14,15)
!
!   C_4v and C_6v
!
         DO itemp=2,ntemp-1
            e_piezo_tensor_t(3,2,itemp)= e_piezo_tensor_t(3,1,itemp)
            e_piezo_tensor_t(2,4,itemp)= e_piezo_tensor_t(1,5,itemp)
         ENDDO
   CASE(13)
!
!   C_3v
!

         DO itemp=2,ntemp-1
            e_piezo_tensor_t(1,2,itemp)=-e_piezo_tensor_t(2,2,itemp)
            e_piezo_tensor_t(1,6,itemp)= e_piezo_tensor_t(2,2,itemp)
            e_piezo_tensor_t(3,2,itemp)= e_piezo_tensor_t(3,1,itemp)
            e_piezo_tensor_t(2,4,itemp)= e_piezo_tensor_t(1,5,itemp)
         ENDDO
   CASE(12)
!
!   C_2v ! all elements are independent
!

   CASE(11,10)
!
!   D_6, D_4
!

         DO itemp=2,ntemp-1
            e_piezo_tensor_t(2,5,itemp)=-e_piezo_tensor_t(1,4,itemp)
         ENDDO

   CASE(9)
!
!   D_3
!

         DO itemp=2,ntemp-1
            e_piezo_tensor_t(1,2,itemp)=-e_piezo_tensor_t(1,1,itemp)
            e_piezo_tensor_t(2,5,itemp)=-e_piezo_tensor_t(1,4,itemp)
            e_piezo_tensor_t(2,6,itemp)=2.0_DP*e_piezo_tensor_t(1,1,itemp)
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
            e_piezo_tensor_t(3,2,itemp)= e_piezo_tensor_t(3,1,itemp)
            e_piezo_tensor_t(2,5,itemp)=-e_piezo_tensor_t(1,4,itemp)
         ENDDO
   CASE(5)
!
!   C_3
!
!      To be checked it does not exist in previous routines
!
!         DO itemp=2,ntemp-1
!            WRITE(iu_piezo,'(e16.8,5e20.12)') temp(itemp),  &
!                e_piezo_tensor_t(1,1,itemp),                &  
!                e_piezo_tensor_t(2,2,itemp),                &  
!                e_piezo_tensor_t(3,1,itemp),                &  
!                e_piezo_tensor_t(3,3,itemp)  
!         ENDDO

   CASE(4)
!
!   C_2
!
          ! all elements are independent

   CASE(3)
!
!   C_s   ! to be checked
!

   CASE DEFAULT
!     No symmetry in the default case
!
END SELECT

RETURN
END SUBROUTINE expand_piezo_tensor

!-----------------------------------------------------------------------
FUNCTION get_pt_type(code_group)
!-----------------------------------------------------------------------
INTEGER :: get_pt_type
INTEGER, INTENT(IN) :: code_group

INTEGER :: itype, aux_type

aux_type=0
DO itype=1,pt_types
   IF (pt_code_group(itype)==code_group) aux_type=itype
ENDDO
IF (aux_type==0) CALL errore('get_pt_type','code_group not available',1)

get_pt_type=aux_type
RETURN
END FUNCTION get_pt_type

!-------------------------------------------------------------------------
SUBROUTINE compute_relax_piezo(ibrav, code_group, nat, max_nint_var, &
                               nint_var_ec, stypes, piezo_relax,     &
                               zeu_eq, dtau_dint, duint_depsilon)
!-------------------------------------------------------------------------

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, code_group, nat, max_nint_var, stypes
INTEGER, INTENT(IN) :: nint_var_ec(stypes)
REAL(DP), INTENT(OUT) :: piezo_relax(3,6)
REAL(DP), INTENT(IN) :: zeu_eq(3,3,nat)
REAL(DP), INTENT(IN) :: dtau_dint(3,nat,max_nint_var,stypes)
REAL(DP), INTENT(IN) :: duint_depsilon(max_nint_var,stypes)

INTEGER :: itypes, iint, ipol, jpol, i, j, na
REAL(DP) :: piezo_relax_aux(3,stypes)
LOGICAL :: check_group_ibrav


piezo_relax=0.0_DP
piezo_relax_aux=0.0_DP
DO itypes=1, stypes
   DO iint=1, nint_var_ec(itypes)
      DO ipol=1,3
         DO jpol=1,3
            DO na=1,nat
               piezo_relax_aux(ipol,itypes) = piezo_relax_aux(ipol,itypes) +  &
                     zeu_eq(ipol,jpol,na) * dtau_dint(jpol,na,iint,itypes) * &
                                            duint_depsilon(iint,itypes)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group)
     CASE(2,16,18,19,20,22,23,25,27,29,32)
     CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
!
!         WRITE(stdout,'(5x,"(  .    .    .    .   d15   .  )")') 
!         WRITE(stdout,'(5x,"(  .    .    .   d15   .    .  )")') 
!         WRITE(stdout,'(5x,"( d31  d31  d33   .    .    .  )")') 

        IF (stypes /= 3) CALL errore('compute_relax_piezo','problem with &
                         C_4v of C_6v',1)
!
!   The factor 0.5 is due to the fact that in strain type 1 we have both
!   epsilon_xx and epsilon_yy
!   The factor 0.5 in piezo_relax(1,5) is due to the fact that
!   e_5 = 2.0 * e_xz
!
        piezo_relax(3,1) = piezo_relax_aux(3,1) * 0.5_DP 
        piezo_relax(3,2) = piezo_relax(3,1)
        piezo_relax(3,3) = piezo_relax_aux(3,2)
        piezo_relax(1,5) = piezo_relax_aux(1,3) * 0.5_DP
        piezo_relax(2,4) = piezo_relax(1,5)
     CASE DEFAULT
        CALL errore('compute_relax_piezo', 'point group not implemented',1) 
  END SELECT
ELSE
   IF (stypes /= 6) CALL errore('compute_relax_piezo','problem with &
                         generic piezo',1)
   DO i=1,3
      DO j=1,6
         piezo_relax(i,j) = piezo_relax_aux(i,j)
      ENDDO
   ENDDO
ENDIF


RETURN
END SUBROUTINE compute_relax_piezo

END MODULE piezoelectric_tensor
