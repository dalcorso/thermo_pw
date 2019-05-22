!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE piezoelectric_tensor
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

  REAL(DP) :: d_piezo_tensor(3,6)   ! The piezoelectric tensor d_{\alpha,m}
 
  REAL(DP), ALLOCATABLE :: polar_geo(:,:) ! The polarization for each strain

  INTEGER :: nppl      ! number of points per line in berry phase calculation

  PUBLIC g_piezo_tensor, polar_geo, compute_piezo_tensor, &
         print_d_piezo_tensor, print_g_piezo_tensor,      &
         compute_d_piezo_tensor, d_piezo_tensor, nppl 

CONTAINS
!
SUBROUTINE print_d_piezo_tensor(frozen_ions)
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
fact= alat * electron_si / (bohr_radius_si)**2 / omega * 10000.0_DP
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (d_piezo_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_d_piezo_tensor

SUBROUTINE print_g_piezo_tensor(frozen_ions)
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
WRITE(stdout,'(5x,"Piezoelectric tensor gamma_ij * a^2 / e")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= alat ** 3  / omega
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (g_piezo_tensor(i,j)*fact, j=1,6)
ENDDO

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Piezoelectric tensor gamma_ij [ 10^{-2} e/(a.u.)^2 ]")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= alat * 100.0_DP / omega
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (g_piezo_tensor(i,j)*fact, j=1,6)
ENDDO

WRITE(stdout,'(/,5x,a)') TRIM(fi_string)
WRITE(stdout,'(5x,"Piezoelectric tensor gamma_ij [C/m^2] ")')
WRITE(stdout,'(4x,"i j=",i9,5i12)') (i, i=1,6)
fact= alat * electron_si / (bohr_radius_si)**2 / omega
DO i=1,3
   WRITE(stdout,'(i5, 6f12.5)') i, (g_piezo_tensor(i,j)*fact, j=1,6)
ENDDO
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_g_piezo_tensor

SUBROUTINE compute_piezo_tensor(polar_geo, epsil_geo, nwork, ngeo, &
                                           ibrav, code_group)
!
!  This routine computes the piezoelectric tensor g_{\alpha,m} by fitting the 
!  polarization strain relation with a second order polynomial. This is 
!  calculated on the basis of the solid point group.
!
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: polar_geo(3,nwork), epsil_geo(3,3,nwork)
INTEGER, INTENT(IN) :: ngeo, ibrav, code_group, nwork
INTEGER :: i, j, igeo, alpha, ind, mn
LOGICAL check_group_ibrav

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

RETURN
END SUBROUTINE compute_piezo_tensor

SUBROUTINE piezo_ij(ialpha, mn, ngeo, epsil_geo, polar_geo)
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit
USE voigt, ONLY : voigt_index

IMPLICIT NONE
INTEGER, INTENT(IN) :: mn, ialpha, ngeo
REAL(DP), INTENT(IN) :: epsil_geo(3,3,ngeo), polar_geo(3,ngeo)
INTEGER :: igeo, m, n, mnin
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo)

WRITE(stdout,'(/,20x,40("-"),/)')
mnin=mn
CALL voigt_index(m,n,mnin,.FALSE.)
WRITE(stdout,'("Piezo ",2i5)') ialpha, mn
DO igeo=1,ngeo
   x(igeo)=epsil_geo(m,n,igeo)
   y(igeo)=polar_geo(ialpha,igeo)
   WRITE(stdout,'(2f15.10)') x(igeo), y(igeo)
ENDDO
CALL polyfit( x, y, ngeo, alpha, m1 )
g_piezo_tensor(ialpha, mn) = alpha(2)
WRITE(stdout,'(/,20x,40("-"),/)')
!
!  The piezoelectric tensor relates the polarization to the strain in voigt
!  notation. Since e_23 = 0.5 e_4, e_13 = 0.5 e_5, e_12 = 0.5 e_6 we have
!  to divide by 2 the elements of the piezoelectric tensor calculated with 
!  off diagonal strain components
!
IF (m /= n) g_piezo_tensor(ialpha, mn) = g_piezo_tensor(ialpha, mn) * 0.5_DP

RETURN
END SUBROUTINE piezo_ij

SUBROUTINE compute_d_piezo_tensor(smn)
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
                                   g_piezo_tensor(alpha,n) * smn(n, m)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE compute_d_piezo_tensor

END MODULE piezoelectric_tensor
