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
 
  REAL(DP), ALLOCATABLE :: polar_geo(:,:) ! The polarization for each strain
                                    ! in units of e bohr/Omega, Omega in bohr^3
  REAL(DP), ALLOCATABLE :: tot_b_phase(:,:) ! Total Berry phase (elec. + ions)
                                    ! in each direction for each strain

  INTEGER :: nppl      ! number of points per line in berry phase calculation

  PUBLIC g_piezo_tensor, polar_geo,  compute_improper_piezo_tensor, &
         compute_proper_piezo_tensor,                     &
         print_d_piezo_tensor, print_g_piezo_tensor,      &
         print_e_piezo_tensor, e_piezo_tensor,            &
         eg_piezo_tensor, print_eg_piezo_tensor,          &
         compute_d_piezo_tensor, d_piezo_tensor, nppl,    &
         compute_polarization_equil,                      &
         proper_improper_piezo, clean_piezo_tensor,       &
         print_piezo_info, tot_b_phase, allocate_piezo,   &
         deallocate_piezo, write_piezo_tensor,            &
         read_piezo_tensor

CONTAINS
!
!---------------------------------------------------------------------------
SUBROUTINE print_d_piezo_tensor(frozen_ions)
!---------------------------------------------------------------------------
!
!  This routine writes on output the piezoelectric tensor
!
USE kinds, ONLY : DP
USE constants, ONLY : electron_si, bohr_radius_si
USE cell_base, ONLY : alat, omega
IMPLICIT NONE
LOGICAL, INTENT(IN) :: frozen_ions
REAL(DP) :: fact
INTEGER :: i, j
CHARACTER(LEN=30) :: fi_string

fi_string=''
IF (frozen_ions) fi_string="Frozen ions"

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Piezoelectric tensor d_ij [pC/N] ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
!
!  the factor 10000.0 comes from the fact that the compliances were in 1/kbar
!  that is in 1/10^8 Pa, to have pC/N we need to multiply and divide by 10000
!
fact= electron_si / (bohr_radius_si)**2 * 10000.0_DP
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (d_piezo_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_d_piezo_tensor

!---------------------------------------------------------------------------
SUBROUTINE print_g_piezo_tensor(frozen_ions)
!---------------------------------------------------------------------------
!
!  This routine writes on output the piezoelectric tensor
!  The piezoelectric tensor enters in units of e / bohr**2
!
USE kinds, ONLY : DP
USE constants, ONLY : electron_si, bohr_radius_si
IMPLICIT NONE
LOGICAL, INTENT(IN) :: frozen_ions
REAL(DP) :: fact
INTEGER :: i, j
CHARACTER(LEN=30) :: fi_string

fi_string=''
IF (frozen_ions) fi_string="Frozen ions"

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Improper piezoelectric tensor &
                                &gamma_ij [ 10^{-2} e/(a.u.)^2 ]")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= 100.0_DP 
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (g_piezo_tensor(i,j)*fact, j=1,6)
ENDDO

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Improper piezoelectric tensor gamma_ij [C/m^2] ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= electron_si / (bohr_radius_si)**2 
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (g_piezo_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_g_piezo_tensor
!
!---------------------------------------------------------------------------
SUBROUTINE print_e_piezo_tensor(frozen_ions)
!---------------------------------------------------------------------------
!
!  This routine writes on output the proper piezoelectric tensor
!  The piezoelectric tensor enters in units of e / bohr**2
!
USE kinds, ONLY : DP
USE constants, ONLY : electron_si, bohr_radius_si
IMPLICIT NONE
LOGICAL, INTENT(IN) :: frozen_ions
REAL(DP) :: fact
INTEGER :: i, j
CHARACTER(LEN=30) :: fi_string

fi_string=''
IF (frozen_ions) fi_string="Frozen ions"

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Proper piezoelectric tensor &
                                &gamma_ij [ 10^{-2} e/(a.u.)^2 ]")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= 100.0_DP 
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (e_piezo_tensor(i,j)*fact, j=1,6)
ENDDO

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Proper piezoelectric tensor gamma_ij [C/m^2] ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= electron_si / (bohr_radius_si)**2 
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (e_piezo_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_e_piezo_tensor
!
!---------------------------------------------------------------------------
SUBROUTINE print_eg_piezo_tensor(frozen_ions)
!---------------------------------------------------------------------------
!
!  This routine writes on output the proper piezoelectric tensor
!  obtained transforming the improper one with the 
!  polarization of the unperturbed structure.
!  The piezoelectric tensor enters in units of e / bohr**2
!
USE kinds, ONLY : DP
USE constants, ONLY : electron_si, bohr_radius_si
IMPLICIT NONE
LOGICAL, INTENT(IN) :: frozen_ions
REAL(DP) :: fact
INTEGER :: i, j
CHARACTER(LEN=30) :: fi_string

fi_string=''
IF (frozen_ions) fi_string="Frozen ions"
WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Proper transformed piezoelectric tensor gamma_ij [C/m^2] ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= electron_si / (bohr_radius_si)**2 
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (eg_piezo_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_eg_piezo_tensor

!-------------------------------------------------------------------------
SUBROUTINE allocate_piezo(nwork)
!-------------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: nwork

ALLOCATE(polar_geo(3,nwork))
ALLOCATE(tot_b_phase(3,nwork))
polar_geo=0.0_DP
tot_b_phase=0.0_DP

RETURN
END SUBROUTINE allocate_piezo
!
!-------------------------------------------------------------------------
SUBROUTINE deallocate_piezo()
!-------------------------------------------------------------------------
IMPLICIT NONE

IF (ALLOCATED(polar_geo)) DEALLOCATE(polar_geo)
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
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: polar0(3)
INTEGER :: find_free_unit
INTEGER :: outunit, ios, i, j

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_piezo_tensor','ploblem opening output file', ABS(ios))

IF (ionode) THEN
   WRITE(outunit,'("Spontaneous polarization (e/bohr**2)")')
   WRITE(outunit,'(4e20.10)') (polar0(i), i=1,3)
   WRITE(outunit,*)
   WRITE(outunit,'("Improper piezoelectric tensor (e/bohr**2)")')
   DO i=1,3
      WRITE(outunit,'(4e20.10)') (g_piezo_tensor(i,j), j=1,6)
   ENDDO
   WRITE(outunit,*)
   WRITE(outunit,'("Proper corrected piezoelectric tensor (e/bohr**2)")')
   DO i=1,3
      WRITE(outunit,'(4e20.10)') (eg_piezo_tensor(i,j), j=1,6)
   END DO
   WRITE(outunit,*)
   WRITE(outunit,'("Proper piezoelectric tensor (e/bohr**2)")')
   DO i=1,3
      WRITE(outunit,'(4e20.10)') (e_piezo_tensor(i,j), j=1,6)
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
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP) :: polar0(3)
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
   WRITE(inunit,'(4e20.10)') (polar0(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(4e20.10)',ERR=100,IOSTAT=ios) (g_piezo_tensor(i,j), j=1,6)
   ENDDO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(4e20.10)',ERR=100,IOSTAT=ios) (eg_piezo_tensor(i,j), j=1,6)
   END DO
   READ(inunit,*)
   READ(inunit,*)
   DO i=1,3
      READ(inunit,'(4e20.10)',ERR=100,IOSTAT=ios) (e_piezo_tensor(i,j), j=1,6)
   END DO
   CLOSE(inunit)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
IF (ios /= 0) THEN
   exists=.FALSE.
   RETURN
ENDIF
CALL mp_bcast(g_piezo_tensor,ionode_id,intra_image_comm)
CALL mp_bcast(eg_piezo_tensor,ionode_id,intra_image_comm)
CALL mp_bcast(e_piezo_tensor,ionode_id,intra_image_comm)
exists=.TRUE.

RETURN
END SUBROUTINE read_piezo_tensor
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
LOGICAL check_group_ibrav

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
LOGICAL check_group_ibrav

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
LOGICAL check_group_ibrav

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
!  The piezoelectric tensor relates the polarization to the strain in voigt
!  notation. Since e_23 = 0.5 e_4, e_13 = 0.5 e_5, e_12 = 0.5 e_6 we have
!  to divide by 2 the elements of the piezoelectric tensor calculated with 
!  off diagonal strain components
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
   CALL to_voigt3(piezo_voigt_aux, piezo_in_aux, .FALSE.)
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
   CALL to_voigt3(piezo_out, piezo_in_aux, .TRUE.)
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
   CALL to_voigt3(piezo_voigt_aux, piezo_in_aux, .FALSE.)
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
   CALL to_voigt3(piezo_out, piezo_in_aux, .TRUE.)
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

END MODULE piezoelectric_tensor
