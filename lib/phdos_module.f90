!
! Copyright (C) 2010 Quantum ESPRESSO group. Inspired to the fqha.f90
! routine of QE and to the programs in the QHA directory by 
! Eyvaz Isaev, Department of Physics, Chemistry, and Biophysics (IFM), 
! Linkoping University, Sweden.
! Theoretical Physics Department, Moscow State Institute of Steel and Alloys,
! Russia.
! Materials Theory Group, Institute of Physics and Materials Science, 
! Uppsala University, Sweden.
!

!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE phdos_module
!
!  This module provide methods to read phonon dos files and to calculate 
!  their contribution to the free energy. It defines a type phdos that
!  contains the phonon dos as a function of energy.
!
! USE kinds, ONLY : dp
!
USE kinds, ONLY : DP
IMPLICIT NONE
SAVE
PRIVATE
REAL(DP), PARAMETER :: kb=8.6173324d-5/13.6058d0 ! Boltzmann constant in Ry/K
REAL(DP), PARAMETER :: ry_to_cmm1= 8065.5d0 * 13.6058d0
REAL(DP), PARAMETER :: kb1=1.0d0/8065.5d0/8.6173324d-5 ! inverse Boltzmann 
                                                       ! constant in cm^{-1}/K

TYPE phdos_type
   INTEGER :: number_of_points    ! rhe number of points
   REAL(DP) :: de                 ! interval of the mesh of frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: nu(:)     ! the frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: phdos(:)  ! the phdos (states/ cm-1)
END TYPE phdos_type

PUBLIC :: phdos_type, read_phdos_data, zero_point_energy, free_energy, &
          vib_energy, vib_entropy, specific_heat_cv, &
          integrated_dos, set_phdos, destroy_phdos
          
CONTAINS

SUBROUTINE set_phdos(phdos,ndiv,deltanu)
IMPLICIT NONE
TYPE(phdos_type), INTENT(INOUT) :: phdos
INTEGER, INTENT(IN) :: ndiv
REAL(DP), INTENT(IN) :: deltanu

phdos%number_of_points=ndiv
phdos%de=deltanu
ALLOCATE(phdos%nu(ndiv))
ALLOCATE(phdos%phdos(ndiv))

RETURN
END SUBROUTINE set_phdos

SUBROUTINE destroy_phdos(phdos)
IMPLICIT NONE
TYPE(phdos_type), INTENT(INOUT) :: phdos

IF (ALLOCATED(phdos%nu)) DEALLOCATE(phdos%nu)
IF (ALLOCATED(phdos%phdos)) DEALLOCATE(phdos%phdos)

RETURN
END SUBROUTINE destroy_phdos

SUBROUTINE read_phdos_data(phdos, filename)
!
!  This subroutine reads the phdos from a file. It allocates space,
!  opens and closes the phdos file.
!
TYPE(phdos_type), INTENT(INOUT) :: phdos
CHARACTER(LEN=256), INTENT(IN) :: filename
INTEGER :: iunit, ios
INTEGER, PARAMETER :: ndivx=100000
REAL(DP), ALLOCATABLE :: nu(:), dos(:)
REAL(DP) :: de, de_
INTEGER :: i, ndiv

iunit=65
OPEN(file=TRIM(filename), unit=iunit, status='old', form='formatted',  &
     err=100, iostat=ios)
100 IF (ios /= 0) STOP 'opening file'

ALLOCATE(nu(ndivx))
ALLOCATE(dos(ndivx))
de = 0d0
DO i=1,ndivx
     ! nu(i) = frequencies (cm^{-1}), dos(i) in states/cm^{-1} 
   READ(iunit, *, END=20, ERR=10, IOSTAT=ios) nu(i),dos(i)
   IF ( nu(i) < -1.d0 ) THEN
      STOP ' wrong grid: omega < 0'
   ELSE IF ( nu(i) < 0.d0 ) THEN
      nu(i) = 0.d0
   END IF
   IF ( i > 1 ) THEN
      de = nu(i) - nu(i-1)
      IF ( i > 2 ) THEN
         de_ = nu(i) - nu(i-1)
         IF ( ABS(de - de_) > 1.0d-4 ) STOP ' wrong grid: not uniform'
      END IF
   END IF
   ndiv=i
ENDDO
10 IF (ios /= 0 ) STOP 'problem reading phdos'
20 continue

phdos%number_of_points=ndiv
phdos%de=de
ALLOCATE(phdos%nu(ndiv))
ALLOCATE(phdos%phdos(ndiv))
phdos%nu(:) = nu(1:ndiv)
phdos%phdos(:) = dos(1:ndiv)

DEALLOCATE(nu)
DEALLOCATE(dos)
CLOSE(iunit)

RETURN
END SUBROUTINE read_phdos_data

SUBROUTINE zero_point_energy(phdos, ener)
!
!  This subroutine receives as input a phdos and computes the zero point 
!  energy that corresponds to that phdos. The output energy is in Ry.
!
!USE constants, ONLY : RY_TO_CMM1
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(OUT) :: ener
INTEGER :: ndiv

ndiv=phdos%number_of_points
ener = 0.5_DP * phdos%de*dot_product(phdos%phdos(1:ndiv),phdos%nu(1:ndiv))
! result is in cm^{-1}, bring it to Ry
ener = ener / ry_to_cmm1 

RETURN
END SUBROUTINE zero_point_energy

SUBROUTINE free_energy(phdos, temp, ener)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational free energy at that temperature. ener contains
!  only the vibrational contribution WITHOUT the zero point energy. 
!  
!
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: ener

INTEGER :: ndiv, i
REAL(DP) :: nu, arg, temp1

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   arg= kb1 * nu * temp1
   IF (nu > 0.0_DP) &
      ener = ener + phdos%phdos(i)* kb * temp * LOG( 1.0_DP - EXP( - arg ) )
ENDDO
ener = ener*phdos%de

RETURN
END SUBROUTINE free_energy

SUBROUTINE vib_energy(phdos, temp, ener)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational energy at that temperature. ener contains
!  the energy WITHOUT the zero point energy. 
!  
!
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: ener

INTEGER :: ndiv, i
REAL(DP) :: nu, temp1, arg

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   arg= kb1 * nu * temp1
   IF (nu > 0.d0 .AND. arg < 700._DP) ener = ener + phdos%phdos(i)* nu /  & 
                                           ( EXP( arg ) - 1.0_DP )
ENDDO
ener = ener * phdos%de / ry_to_cmm1

RETURN
END SUBROUTINE vib_energy

SUBROUTINE vib_entropy(phdos, temp, entr)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational entropy at that temperature. 
!  
!
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: entr
REAL(DP) :: ener, free_ener

CALL free_energy(phdos, temp, free_ener)
CALL vib_energy(phdos, temp, ener)

IF (temp > 0.0_DP) THEN
   entr = ( ener - free_ener ) / temp
ELSE
   entr = 0.0_DP
ENDIF

RETURN
END SUBROUTINE vib_entropy

SUBROUTINE specific_heat_cv(phdos, temp, cv)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the constant volume specific heat at that temperature. 
!  The output cv is in Ry / K.
!
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: cv

INTEGER :: ndiv, i
REAL(DP) :: nu, temp1, arg

cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   arg= kb1 * nu * temp1
   IF (nu > 0.d0 .AND. arg < 650._DP) cv = cv + phdos%phdos(i) * EXP(arg) * &
                                        ( arg / ( EXP( arg ) - 1.0_DP )) ** 2 
ENDDO
cv = cv * phdos%de * kb

RETURN
END SUBROUTINE specific_heat_cv

SUBROUTINE integrated_dos(phdos, tot_dos)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational energy at that temperature. ener contains
!  the energy WITHOUT the zero point energy. 
!  
!
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(OUT) :: tot_dos
INTEGER :: ndiv, i
REAL(DP) :: nu

tot_dos=0.0_DP
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   IF (nu > 0.d0) tot_dos = tot_dos + phdos%phdos(i)
ENDDO

tot_dos = tot_dos * phdos%de
RETURN
END SUBROUTINE integrated_dos

END MODULE phdos_module
