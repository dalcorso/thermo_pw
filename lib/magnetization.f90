!
! Copyright (C) 2026 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE magnetization_vector
!---------------------------------------------------------------------------
!
!   this module contains the support routines for the calculation
!   of the magnetization. Presently it contains: 
!   Routines to read and write the magnetization on file;
!   Routines that contains the form of the vector. 
!   Routines to compute the magnetic susceptibility
!   
!   Magnetization is an c-axial vector of odd rank (1) and the
!   point group to use to deal with it must be the B group of Birss.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE

  INTEGER :: mod_tot
!
!   Some array to simplify dealing with pyroelectric tensor
!
  INTEGER, PARAMETER :: m_elements=3
  
  CHARACTER(LEN=6) :: m_names(m_elements)

  DATA  m_names / 'm_{1} ', 'm_{2} ', 'm_{3} ' /

  INTEGER, PARAMETER :: m_types=12

  INTEGER :: m_code_group(m_types)  ! code of the point group for each type
  DATA  m_code_group /  1, 3, 3, 4, 4, 5, 6, 7, 12, 13, 14, 15 /

  INTEGER  :: m_present(m_elements, m_types)

  DATA m_present / &
       1,2,3, & ! 1  C_1
       1,0,3, & ! 3  C_s  ! b unique
       1,2,0, & ! 3  C_s  ! c unique 
       0,2,0, & ! 4  C_2  ! b unique
       0,0,3, & ! 4  C_2  ! c unique
       0,0,3, & ! 5  C_3
       0,0,3, & ! 6  C_4
       0,0,3, & ! 7  C_6
       0,0,3, & ! 12 C_2v
       0,0,3, & ! 13 C_3v
       0,0,3, & ! 14 C_4v
       0,0,3  / ! 15 C_6v

!
!  Magnetic point groups corresponding to C_1 
!  1, -1, 
!
!  Magnetic point groups corresponding to C_s 
!  2', m', 2'/m' 
!
!  Magnetic point groups corresponding to C_2 
!  2, m, 2/m 
!
!  Magnetic point groups corresponding to C_3 
!  3, -3
!
!  Magnetic point groups corresponding to C_4 
!  4, -4, 4/m
!
!  Magnetic point groups corresponding to C_6 
!  6, -6, 6/m
!
!  Magnetic point groups corresponding to C_2v 
!  2'2'2, m'm'2, m'm2'
!
!  Magnetic point groups corresponding to C_2v 
!  2'2'2, m'm'2, m'm2', m'm'm
!
!  Magnetic point groups corresponding to C_3v 
!  32', 3m', -3m'
!
!  Magnetic point groups corresponding to C_4v 
!  42'2', 4m'm', -42'm', 4/mm'm'
!
!  Magnetic point groups corresponding to C_6v 
!  62'2', 6m'm', -62'm', 6/mm'm'
!

  PUBLIC write_magnetization, read_magnetization,    &
         m_names, m_types, m_present, m_code_group,  &
         get_m_type, m_elements, compute_magnetic_susceptibility,  &
         write_magnetic_susceptibility, set_magnetic_field_for_ms, &
         check_magn_susc, compute_magnetic_susceptibility_from_magnetoelectric

CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE write_magnetization(filename,mag)
!-------------------------------------------------------------------------
!
!  This routine writes the magnetization on file.
!  It saves: 
!  the magnetization of the current geometry
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: mag(3)
REAL(DP) :: fact
INTEGER :: find_free_unit
INTEGER :: outunit, i, ios

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_magnetization','ploblem opening output file', ABS(ios))

fact= 1.0
IF (ionode) THEN
   WRITE(outunit,'("Magnetization (cartesian &
                                       &coordinates) (Bohr magnetons)")')
   WRITE(outunit,'(3e20.10)') (mag(i), i=1,3)
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_magnetization
!
!-------------------------------------------------------------------------
SUBROUTINE read_magnetization(filename,mag)
!-------------------------------------------------------------------------
!
!  This routine reads the magnetization from file.
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(OUT) :: mag(3)
INTEGER :: find_free_unit
INTEGER :: inunit, i, ios

IF (ionode) THEN
   inunit=find_free_unit()
   OPEN(UNIT=inunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('read_magnetization','ploblem opening input file', ABS(ios))

IF (ionode) THEN
   READ(inunit, *)
   READ(inunit,'(3e20.10)') (mag(i), i=1,3)
   CLOSE(inunit)
ENDIF
CALL mp_bcast(mag,ionode_id,intra_image_comm)

RETURN
END SUBROUTINE read_magnetization
!
!-----------------------------------------------------------------------
FUNCTION get_m_type(code_group, ibrav)
!-----------------------------------------------------------------------
INTEGER :: get_m_type
INTEGER, INTENT(IN) :: code_group, ibrav

INTEGER :: itype, aux_type

aux_type=0
DO itype=1,m_types
   IF (m_code_group(itype)==code_group) aux_type=itype
ENDDO
IF (aux_type==0) CALL errore('get_py_type','code_group not available',1)
!
!  for code group 3 (C_s) and 4 (C_2) we must distinguish if the 
!  monoclinic lattice is b unique or c unique
!
IF (code_group==3) THEN
   IF (ibrav < 0) THEN 
      aux_type=2
   ELSE
      aux_type=3
   ENDIF
ENDIF
IF (code_group==4) THEN
   IF (ibrav < 0) THEN 
      aux_type=4
   ELSE
      aux_type=5
   ENDIF
ENDIF
get_m_type=aux_type

RETURN
END FUNCTION get_m_type
!
!----------------------------------------------------------------------------
SUBROUTINE compute_magnetic_susceptibility(mag_geo, bfield_geo, nwork, &
                         ibrav, omega, magnetic_susceptibility)
!----------------------------------------------------------------------------
!
!  This routine computes the magnetic susceptibility assuming that
!  the magnetic fields are those set in set_magnetic_field_for_ms
!  It works in all solids.
!
USE kinds,  ONLY : DP
USE constants, ONLY : munought_si, electron_si, h_planck_si, &
                      electronmass_si, pi, bohr_radius_si, hartree_si
IMPLICIT NONE

INTEGER :: nwork, ibrav
REAL(DP), INTENT(IN) :: mag_geo(3,nwork), bfield_geo(3, nwork), omega
REAL(DP), INTENT(OUT) :: magnetic_susceptibility(3,3)

INTEGER :: ipol, jpol, ngeo_b, ind, mag_susc_nc
REAL(DP) :: mu_b, omega_m3, fact

magnetic_susceptibility(:,:)=0.0_DP
SELECT CASE(ibrav) 
   CASE(1,2,3)
      mag_susc_nc = 1
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      magnetic_susceptibility(2,2)= magnetic_susceptibility(1,1)
      magnetic_susceptibility(3,3)= magnetic_susceptibility(1,1)
   CASE(4,5,6,7)
      mag_susc_nc = 2
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      magnetic_susceptibility(2,2)= magnetic_susceptibility(1,1)
      ind = ngeo_b+1 
      CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))
   CASE(8,9,10,11)
      mag_susc_nc = 3
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      ind = ngeo_b+1 
      CALL compute_one_mag_susc(2,2,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(2,2))
      ind = 2*ngeo_b+1 
      CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))
   CASE(12,13)
      mag_susc_nc = 3
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      CALL compute_one_mag_susc(2,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(2,1))
      magnetic_susceptibility(1,2)= magnetic_susceptibility(2,1)
      ind = ngeo_b+1 
      CALL compute_one_mag_susc(2,2,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(2,2))
      ind = 2*ngeo_b+1 
      CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))
   CASE(-12,-13)
      mag_susc_nc = 3
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      CALL compute_one_mag_susc(3,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(3,1))
      magnetic_susceptibility(1,3)= magnetic_susceptibility(3,1)
      ind = ngeo_b+1 
      CALL compute_one_mag_susc(2,2,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(2,2))
      ind = 2*ngeo_b+1 
      CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))
   CASE DEFAULT
      mag_susc_nc = 3
      ngeo_b=nwork/mag_susc_nc
      DO ipol = 1, 3
         DO jpol = 1, 3
            ind = (jpol-1) * ngeo_b + 1 
            CALL compute_one_mag_susc(ipol,jpol,mag_geo(1,ind), &
                     bfield_geo(1,ind), ngeo_b, &
                     magnetic_susceptibility(ipol,jpol))
         ENDDO
      ENDDO
END SELECT
!
!  Bring the units to the SI. 
!  The susceptibility is defined as chi=mu_0 d M / d B where M is the
!  magnetization in A/m and B is the magnetic field in Testa. 
!  mu_0 M is in Tesla so the susceptibility is adimentional.
!  In input the volume is in (a.u.)^3, the magnetic field is multiplied by 
!  the Bohr magneton and is in Ry and the magnetic moment per cell is in units 
!  of Bohr magnetons. mag_susc contains the derivative of the magnetic 
!  moment per cell with respect to the magnetic field multiplied by mu_b 
!  in Ry. 
!  We first calculate mu_b the Bohr magneton in SI. 
!
mu_b = electron_si * h_planck_si / 4.0_DP / pi / electronmass_si
!  
!  Multiplication of the magnetic field by hartree_si/2
!  gives mu_b B in Joule, and dividing by mu_B gives the 
!  the magnetic field in Tesla: 
!
fact = 2.0_DP * mu_b  / hartree_si 
!
WRITE(stdout,'(/,5x,"Multiply the magnetic field by", 1pe15.8, &
                       &" to convert it to Tesla")') 1.0_DP / fact
!
!  omega in cube meters
!
omega_m3= omega * (bohr_radius_si)**3
!
!  Another mu_b transforms the magnetic
!  moment per cell into SI A m^2. Dividing by the volume 
!  omega_m3 which is in m^3, we obtain the magnetization M
!  in A/m. 
!
!  magnetization in A/m from magnetic moment per cell in unit of mu_b
!
fact= mu_b * fact / omega_m3
!
! finally adimensional susceptibility is obtained multiplying by mu_0
!
fact = fact * munought_si 

magnetic_susceptibility(:,:)=magnetic_susceptibility(:,:) * fact

RETURN
END SUBROUTINE compute_magnetic_susceptibility
!
!----------------------------------------------------------------------------
SUBROUTINE compute_magnetic_susceptibility_from_magnetoelectric(mag_geo, &
                         bfield_geo, nwork, ibrav, code_group, omega,    &
                         magnetic_susceptibility)
!----------------------------------------------------------------------------
!
!  This routine computes the magnetic susceptibility assuming that
!  the magnetic fields are those set in set_magnetic_field_for_mt
!  that are used to compute the magnetoelectric tensor. In some cases
!  these magnetic fields are sufficient to compute the magnetic susceptibility
!  and no new calculation is needed. In other cases this routine exit and
!  one has to use what='scf_magnetic_susceptibility'. This happens for all
!  solids with inversion symmetry for instance.
!
!  code group here should be the use used for the magnetoelectric tensor
!  so the b_birss_code_group
!
USE kinds,  ONLY : DP
USE constants, ONLY : munought_si, electron_si, h_planck_si, &
                      electronmass_si, pi, bohr_radius_si, hartree_si
IMPLICIT NONE

INTEGER :: nwork, code_group, ibrav
REAL(DP), INTENT(IN) :: mag_geo(3,nwork), bfield_geo(3, nwork), omega
REAL(DP), INTENT(OUT) :: magnetic_susceptibility(3,3)

INTEGER :: ipol, jpol, ngeo_b, ind, mag_susc_nc
REAL(DP) :: mu_b, omega_m3, fact

magnetic_susceptibility(:,:)=0.0_DP
SELECT CASE(code_group) 
   CASE(2,16,18,19,20,22,23,25,27,29,32,12,13,14,15,17,21,24)
   CASE(26,9,10,11,5,6,7)
!
!  S_4, D_3, D_4, D_6, C_3, C_4, C_6 
!

      mag_susc_nc = 2
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      magnetic_susceptibility(2,2)= magnetic_susceptibility(1,1)
      ind = ngeo_b+1 
      CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))
   CASE(8)
 !
 !  D_2
 !
      mag_susc_nc = 3
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      ind = ngeo_b+1 
      CALL compute_one_mag_susc(2,2,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(2,2))
      ind = 2*ngeo_b+1 
      CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))
   CASE(28,30,31)
 !
 ! Cubic solids
 !
      mag_susc_nc = 1
      ngeo_b=nwork/mag_susc_nc
      CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
      magnetic_susceptibility(2,2)= magnetic_susceptibility(1,1)
      magnetic_susceptibility(3,3)= magnetic_susceptibility(1,1)
   CASE(3,4) 
      IF (ibrav==-12.OR.ibrav==-13) THEN
         mag_susc_nc = 3
         ngeo_b=nwork/mag_susc_nc
         CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
         CALL compute_one_mag_susc(3,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(3,1))
         magnetic_susceptibility(1,3)= magnetic_susceptibility(3,1)
         ind = ngeo_b+1 
         CALL compute_one_mag_susc(2,2,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(2,2))
         ind = 2*ngeo_b+1 
         CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))
      ELSEIF (ibrav==12.OR.ibrav==13) THEN
         mag_susc_nc = 3
         ngeo_b=nwork/mag_susc_nc
         CALL compute_one_mag_susc(1,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(1,1))
         CALL compute_one_mag_susc(2,1,mag_geo,bfield_geo,ngeo_b,&
                               magnetic_susceptibility(2,1))
         magnetic_susceptibility(1,2)= magnetic_susceptibility(2,1)
         ind = ngeo_b+1 
         CALL compute_one_mag_susc(2,2,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(2,2))
         ind = 2*ngeo_b+1 
         CALL compute_one_mag_susc(3,3,mag_geo(1,ind),bfield_geo(1,ind),ngeo_b,&
                               magnetic_susceptibility(3,3))

      ENDIF
   CASE DEFAULT
!
!  1 C_1
!
      mag_susc_nc = 3
      ngeo_b=nwork/mag_susc_nc
      DO ipol = 1, 3
         DO jpol = 1, 3
            ind = (jpol-1) * ngeo_b + 1 
            CALL compute_one_mag_susc(ipol,jpol,mag_geo(1,ind), &
                     bfield_geo(1,ind), ngeo_b, &
                     magnetic_susceptibility(ipol,jpol))
         ENDDO
      ENDDO
END SELECT
!
!  Bring the units to the SI. 
!  The susceptibility is defined as chi=mu_0 d M / d B where M is the
!  magnetization in A/m and B is the magnetic field in Testa. 
!  mu_0 M is in Tesla so the susceptibility is adimentional.
!  In input the volume is in (a.u.)^3, the magnetic field is multiplied by 
!  the Bohr magneton and is in Ry and the magnetic moment per cell is in units 
!  of Bohr magnetons. mag_susc contains the derivative of the magnetic 
!  moment per cell with respect to the magnetic field multiplied by mu_b 
!  in Ry. 
!  We first calculate mu_b the Bohr magneton in SI. 
!
mu_b = electron_si * h_planck_si / 4.0_DP / pi / electronmass_si
!  
!  Multiplication of the magnetic field by hartree_si/2
!  gives mu_b B in Joule, and dividing by mu_B gives the 
!  the magnetic field in Tesla: 
!
fact = 2.0_DP * mu_b  / hartree_si 
!
!  omega in cube meters
!
omega_m3= omega * (bohr_radius_si)**3
!
!  Another mu_b transforms the magnetic
!  moment per cell into SI A m^2. Dividing by the volume 
!  omega_m3 which is in m^3, we obtain the magnetization M
!  in A/m. 
!
!  magnetization in A/m from magnetic moment per cell in unit of mu_b
!
fact= mu_b * fact / omega_m3
!
! finally adimensional susceptibility is obtained multiplying by mu_0
!
fact = fact * munought_si 

magnetic_susceptibility(:,:)=magnetic_susceptibility(:,:) * fact

RETURN
END SUBROUTINE compute_magnetic_susceptibility_from_magnetoelectric
!
!----------------------------------------------------------------------------
LOGICAL FUNCTION check_magn_susc(code_group)
!----------------------------------------------------------------------------
!
!  This function gives .TRUE. if the point group is one of those for 
!  which the magnetoelectric tensor requires the same magnetic fields 
!  of the magnetic susceptibility. For these groups the magnetic 
!  susceptibility is computed when one computes the magnetoelectric tensor.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group

check_magn_susc=.FALSE.
SELECT CASE (code_group)
   CASE(2,16,18,19,20,22,23,25,27,29,32,12,13,14,15,17,21,24)
   CASE(1,3,4,5,6,7,8,9,10,11,26,28,30,31)
       check_magn_susc=.TRUE.
   CASE DEFAULT
END SELECT

RETURN
END FUNCTION check_magn_susc
!----------------------------------------------------------------------------
SUBROUTINE set_magnetic_field_for_ms(do_field,mag_susc_nc,ibrav)
!----------------------------------------------------------------------------
USE kinds,  ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN)  :: ibrav
INTEGER, INTENT(OUT) :: mag_susc_nc
LOGICAL, INTENT(OUT) :: do_field(3)

do_field=.FALSE.
SELECT CASE(ibrav) 
   CASE(1,2,3)
      mag_susc_nc = 1
      do_field(1)=.TRUE.
   CASE(4,5,6,7)
      mag_susc_nc = 2
      do_field(1) =.TRUE.
      do_field(3) =.TRUE.
   CASE DEFAULT
      mag_susc_nc = 3
      do_field=.TRUE.
END SELECT

RETURN
END SUBROUTINE set_magnetic_field_for_ms

!---------------------------------------------------------------------------
SUBROUTINE compute_one_mag_susc(ipol, jpol, mag_geo, bfield_geo, ndata, &
                                mag_susc)
!---------------------------------------------------------------------------
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit

IMPLICIT NONE
INTEGER, INTENT(IN) :: ipol, jpol, ndata
REAL(DP), INTENT(IN) :: mag_geo(3,ndata), bfield_geo(3,ndata)
REAL(DP), INTENT(OUT) :: mag_susc
INTEGER :: igeo
INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
REAL(DP) :: alpha(m1)          ! the polynomial coefficients
REAL(DP) :: x(ndata), y(ndata)


WRITE(stdout,'(/,5x,"Component (", i1,",", i1,")")') ipol, jpol

DO igeo=1,ndata
   x(igeo)=bfield_geo(jpol,igeo)
   y(igeo)=mag_geo(ipol,igeo)
   WRITE(stdout,'(20x,f15.10,3x,f15.10)') x(igeo), y(igeo)
ENDDO
CALL polyfit( x, y, ndata, alpha, m1-1 )

mag_susc=alpha(2)

RETURN
END SUBROUTINE compute_one_mag_susc
!
!-------------------------------------------------------------------------
SUBROUTINE write_magnetic_susceptibility(filename, mag_susc)
!-------------------------------------------------------------------------
!
!  This routine writes the magnetic susceptibility on file.
!
!  The file contains the 3x3 matrix of the magnetic susceptibility. 
!  The name of the file that contains these quantities is given as 
!  input.
!
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: mag_susc(3,3)

INTEGER :: find_free_unit
INTEGER :: outunit, ios, ipol, jpol

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', &
        FORM='formatted', ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_magnetic_susceptibility',&
                           'ploblem opening output file', ABS(ios))

IF (ionode) THEN
   WRITE(outunit,'("Magnetic susceptibility SI (adimensional) ")')
   DO ipol=1,3
      WRITE(outunit,'(3e19.10)') (mag_susc(ipol,jpol), jpol=1,3)
   ENDDO
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_magnetic_susceptibility
!
END MODULE magnetization_vector
