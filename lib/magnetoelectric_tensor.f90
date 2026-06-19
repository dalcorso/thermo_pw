!
! Copyright (C) 2026 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE magnetoelectric_tensor
!---------------------------------------------------------------------------
!
!   this module contains the support routines for the calculation
!   of the magnetoelectric tensor. 
!
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

PUBLIC compute_magnetoelectric_tensor, set_magnetic_field_for_mt, &
       write_magnetoelectric_tensor

CONTAINS
!
!--------------------------------------------------------------------
SUBROUTINE compute_magnetoelectric_tensor(code_group, bfield_geo, &
           polar_geo, nwork, magnetoelectric)
!--------------------------------------------------------------------
!
USE kinds, ONLY : DP
USE constants, ONLY : munought_si, electron_si, h_planck_si, &
                      electronmass_si, pi, bohr_radius_si, hartree_si

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, code_group
REAL(DP), INTENT(IN) :: bfield_geo(3,nwork), polar_geo(3,nwork)
REAL(DP), INTENT(OUT) :: magnetoelectric(3,3)
CHARACTER(LEN=11) :: group_name

INTEGER :: ind, magnetoelec_nc, ngeo
REAL(DP) :: fact, mu_b

WRITE(stdout,'(/,5x,"Using point group ",a,":",/)') &
                                           TRIM(group_name(code_group))
magnetoelectric=0.0_DP
SELECT CASE (code_group)
   CASE(2,16,18,19,20,22,23,25,27,29)
     WRITE(stdout,'(/,5x,"Symmetry group has inversion. Axial tensors &
                     &of rank 2 vanish")')
   CASE(1)  
!
!   C_1
!
!      WRITE(stdout,'(/,5x, "( a11  a12  a13 )")')
!      WRITE(stdout,'(  5x, "( a21  a22  a23 )")')
!      WRITE(stdout,'(  5x, "( a31  a32  a11 )")')
      magnetoelec_nc = 3
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      CALL magnetoelectric_ij(2, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(2,1))
      CALL magnetoelectric_ij(3, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(3,1))
      ind=ngeo+1
      CALL magnetoelectric_ij(1, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(1,2))
      CALL magnetoelectric_ij(2, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(2,2))
      CALL magnetoelectric_ij(3, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,2))
      ind=2*ngeo+1
      CALL magnetoelectric_ij(1, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(1,3))
      CALL magnetoelectric_ij(2, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(2,3))
      CALL magnetoelectric_ij(3, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,3))
   CASE(3)  
!
!  C_s
!
!      WRITE(stdout,'(/,5x, "(  .    .   a13 )")')
!      WRITE(stdout,'(  5x, "(  .    .   a23 )")')
!      WRITE(stdout,'(  5x, "( a31  a32   .  )")')
      magnetoelec_nc = 3
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(3, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(3,1))
      ind=ngeo+1
      CALL magnetoelectric_ij(3, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,2))
      ind=2*ngeo+1
      CALL magnetoelectric_ij(1, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(1,3))
      CALL magnetoelectric_ij(2, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(2,3))
   CASE(5,6,7)  
!
!     C_3, C_4, C_6
!
!      WRITE(stdout,'(/,5x, "( a11  a12   .  )")')
!      WRITE(stdout,'(  5x, "(-a12  a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 2
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      CALL magnetoelectric_ij(2, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(2,1))
      magnetoelectric(2,2)=magnetoelectric(1,1)
      magnetoelectric(1,2)=-magnetoelectric(2,1)
      ind=ngeo+1
      CALL magnetoelectric_ij(3, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,3))
   CASE(4)  
!
!  C_4
!
!      WRITE(stdout,'(/,5x, "( a11  a12   .  )")')
!      WRITE(stdout,'(  5x, "( a21  a22   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 3
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      CALL magnetoelectric_ij(2, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(2,1))
      ind=ngeo+1
      CALL magnetoelectric_ij(1, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(1,2))
      CALL magnetoelectric_ij(2, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(2,2))
      ind=2*ngeo+1
      CALL magnetoelectric_ij(3, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,3))
   CASE(8)  
!
!  D_2
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .   a22   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 3
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      ind=ngeo+1
      CALL magnetoelectric_ij(2, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(2,2))
      ind=2*ngeo+1
      CALL magnetoelectric_ij(3, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,3))
   CASE(9,10,11)  
!
!   D_3, D_4, D_6
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .   a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 2
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      magnetoelectric(2,2)=magnetoelectric(1,1)
      ind=ngeo+1
      CALL magnetoelectric_ij(3, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,3))

   CASE(12)  
!
!   C_2v
!
!      WRITE(stdout,'(/,5x, "(  .   a12   .  )")')
!      WRITE(stdout,'(  5x, "( a21   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .    .    .  )")')
      magnetoelec_nc = 2
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(2, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(2,1))
      ind=ngeo+1
      CALL magnetoelectric_ij(1, 2, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(1,2))

   CASE(14,15)  
!
!   C_3v, C_4v, C_6v
!
!      WRITE(stdout,'(/,5x, "(  .   a12   .  )")')
!      WRITE(stdout,'(  5x, "(-a12   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .    .    .  )")')
!
      magnetoelec_nc = 1
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(2, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(2,1))
      magnetoelectric(1,2)=-magnetoelectric(2,1)
   CASE(24)  
!
!   D_2d
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .  -a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .    .  )")')
      magnetoelec_nc = 1
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      magnetoelectric(2,2)=-magnetoelectric(1,1)

   CASE(26)  
!
!   S_4 
!
!      WRITE(stdout,'(/,5x, "( a11  a12   .  )")')
!      WRITE(stdout,'(  5x, "( a12 -a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')

      magnetoelec_nc = 2
      ngeo=nwork/magnetoelec_nc
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      magnetoelectric(2,2)=-magnetoelectric(1,1)
      CALL magnetoelectric_ij(2, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(2,1))
      magnetoelectric(1,2)=magnetoelectric(2,1)
      ind=ngeo+1
      CALL magnetoelectric_ij(3, 3, ngeo, bfield_geo(1,ind), &
                                    polar_geo(1,ind), magnetoelectric(3,3))

   CASE(28,30,31)  
!
!   T, O
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .   a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a11 )")')
      ngeo=nwork
      CALL magnetoelectric_ij(1, 1, ngeo, bfield_geo, polar_geo, &
                                                     magnetoelectric(1,1))
      magnetoelectric(2,2)=magnetoelectric(1,1)
      magnetoelectric(3,3)=magnetoelectric(1,1)
   CASE DEFAULT 
END SELECT
!
!  Bring the units to the SI. 
!  The magnetoelectric tensor is defined as chi=mu_0 d P / d B where P is the
!  polarization in C/m^2 and B is the magnetic field in Testa. 
!  B / mu_0 is in A/m as the magnetic field H.
!  In input the magnetic field is multiplied by the Bohr magneton and 
!  is in Ry while the polarization is in atomic units e/(a.u.)^2.
!  We first calculate mu_b the Bohr magneton in SI. 
!
mu_b = electron_si * h_planck_si / 4.0_DP / pi / electronmass_si
!  
!  Multiplication of the magnetic field by hartree_si/2
!  gives mu_b B in Joule, and dividing by mu_B gives the 
!  the magnetic field in Tesla: 
!
fact = 2.0_DP * mu_b / hartree_si
!
WRITE(stdout,'(/,5x,"Multiply the magnetic field by", 1pe15.8, &
                       &" to convert it to Tesla")') 1.0_DP / fact
!
!  Here we bring the polarization into C/m^2
!
fact=fact * electron_si / (bohr_radius_si)**2
!
! we multiply by mu_0 and for 10^12 so we will have the 
! magnetoelectric tensor in ps/m
!
fact = fact * munought_si * 1.D12
!
magnetoelectric(:,:)=magnetoelectric(:,:) * fact


RETURN
END SUBROUTINE compute_magnetoelectric_tensor

!--------------------------------------------------------------------
SUBROUTINE set_magnetic_field_for_mt(do_field,magnetoelec_nc,code_group)
!--------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group
INTEGER, INTENT(OUT) :: magnetoelec_nc
LOGICAL, INTENT(OUT) :: do_field(3)

magnetoelec_nc=0
do_field=.FALSE.
SELECT CASE (code_group)
   CASE(2,16,18,19,20,22,23,25,27,29)
     WRITE(stdout,'(/,5x,"Symmetry group has inversion. Axial tensors &
                     &of rank 2 vanish")')
   CASE(1)  
!
!   C_1
!
!      WRITE(stdout,'(/,5x, "( a11  a12  a13 )")')
!      WRITE(stdout,'(  5x, "( a21  a22  a23 )")')
!      WRITE(stdout,'(  5x, "( a31  a32  a11 )")')
      magnetoelec_nc = 3
      do_field=.TRUE.
   CASE(3)  
!
!  C_s
!
!      WRITE(stdout,'(/,5x, "(  .    .   a13 )")')
!      WRITE(stdout,'(  5x, "(  .    .   a23 )")')
!      WRITE(stdout,'(  5x, "( a31  a32   .  )")')
      magnetoelec_nc = 3
      do_field=.TRUE.
   CASE(5,6,7)  
!
!     C_3, C_4, C_6
!
!      WRITE(stdout,'(/,5x, "( a11  a12   .  )")')
!      WRITE(stdout,'(  5x, "(-a12  a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 2
      do_field(1)=.TRUE.
      do_field(3)=.TRUE.
   CASE(4)  
!
!  C_4
!
!      WRITE(stdout,'(/,5x, "( a11  a12   .  )")')
!      WRITE(stdout,'(  5x, "( a21  a22   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 3
      do_field=.TRUE.
   CASE(8)  
!
!  D_2
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .   a22   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 3
      do_field=.TRUE.
   CASE(9,10,11)  
!
!   D_3, D_4, D_6
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .   a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')
      magnetoelec_nc = 2
      do_field(1)=.TRUE.
      do_field(3)=.TRUE.
   CASE(12)  
!
!   C_2v
!
!      WRITE(stdout,'(/,5x, "(  .   a12   .  )")')
!      WRITE(stdout,'(  5x, "( a21   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .    .    .  )")')
      magnetoelec_nc = 2
      do_field(1)=.TRUE.
      do_field(2)=.TRUE.
   CASE(14,15)  
!
!   C_3v, C_4v, C_6v
!
!      WRITE(stdout,'(/,5x, "(  .   a12   .  )")')
!      WRITE(stdout,'(  5x, "(-a12   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .    .    .  )")')
!
      magnetoelec_nc = 1
      do_field(1)=.TRUE.
   CASE(24)  
!
!   D_2d
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .  -a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .    .  )")')
      magnetoelec_nc = 1
      do_field(1)=.TRUE.
   CASE(26)  
!
!   S_4 
!
!      WRITE(stdout,'(/,5x, "( a11  a12   .  )")')
!      WRITE(stdout,'(  5x, "( a12 -a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a33 )")')

      magnetoelec_nc = 2
      do_field(1)=.TRUE.
      do_field(3)=.TRUE.
   CASE(28,30,31)  
!
!   T, O
!
!      WRITE(stdout,'(/,5x, "( a11   .    .  )")')
!      WRITE(stdout,'(  5x, "(  .   a11   .  )")')
!      WRITE(stdout,'(  5x, "(  .    .   a11 )")')
      magnetoelec_nc = 1
      do_field(1)=.TRUE.
   CASE DEFAULT 
      magnetoelec_nc = 3
      do_field=.TRUE.
END SELECT

RETURN
END SUBROUTINE set_magnetic_field_for_mt
!
!---------------------------------------------------------------------------
SUBROUTINE magnetoelectric_ij(ipol, jpol, ngeo, bfield_geo, polar_geo, &
                                                        magnetoelectric)
!---------------------------------------------------------------------------
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit
USE constants, ONLY : electron_si, bohr_radius_si

IMPLICIT NONE
INTEGER, INTENT(IN)   :: ipol, jpol, ngeo
REAL(DP), INTENT(IN)  :: bfield_geo(3,ngeo), polar_geo(3,ngeo)
REAL(DP), INTENT(OUT) :: magnetoelectric
INTEGER :: igeo
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ngeo), y(ngeo), fact
!
!  Here we bring the polarization into C/m^2
!
fact=electron_si / (bohr_radius_si)**2

WRITE(stdout,'(5x,"Component  (",i1,",",i1,"):        &
                    & B (Ry)      polarization (C/m^2)")') ipol, jpol
DO igeo=1,ngeo
   x(igeo)=bfield_geo(jpol,igeo) 
   y(igeo)=polar_geo(ipol,igeo) 
   WRITE(stdout,'(24x,f15.10,3x,f15.10)') x(igeo), y(igeo)*fact
ENDDO
CALL polyfit( x, y, ngeo, alpha, m1-1 )
magnetoelectric = alpha(2)
!
RETURN
END SUBROUTINE magnetoelectric_ij
!
!-------------------------------------------------------------------------
SUBROUTINE write_magnetoelectric_tensor(filename, magnetoelectric)
!-------------------------------------------------------------------------
!
!  This routine writes the magnetoelectric tensor on file. 
!
!  The file contains a 3x3 matrix in ps/m.
!  The name of the file that contains these quantities is given as 
!  input.
!
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: magnetoelectric(3,3)

INTEGER :: find_free_unit
INTEGER :: outunit, ios, ipol, jpol

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', &
        FORM='formatted', ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_magnetoelectric_tensor',&
                           'ploblem opening output file', ABS(ios))

IF (ionode) THEN
   WRITE(outunit,'("Magnetoelectric tensor (ps/m)")')
   DO ipol=1,3
      WRITE(outunit,'(3e19.10)') (magnetoelectric(ipol,jpol), jpol=1,3)
   ENDDO
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_magnetoelectric_tensor

END MODULE magnetoelectric_tensor
