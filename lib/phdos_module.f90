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
USE constants, ONLY :  k_boltzmann_ry, ry_to_cmm1
IMPLICIT NONE
SAVE
PRIVATE

REAL(DP), PARAMETER :: kb=k_boltzmann_ry ! Boltzmann constant in Ry/K
REAL(DP), PARAMETER :: kb1=1.0_DP/kb/ry_to_cmm1 ! inverse Boltzmann 
                                                ! constant in cm^{-1}/K

TYPE phdos_type
   INTEGER :: number_of_points    ! rhe number of points
   REAL(DP) :: de                 ! interval of the mesh of frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: nu(:)     ! the frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: phdos(:)  ! the phdos (states/ cm-1)
END TYPE phdos_type

PUBLIC :: phdos_type, read_phdos_data, zero_point_energy, free_energy, &
          vib_energy, vib_entropy, specific_heat_cv, fecv, &
          integrated_dos, set_phdos, destroy_phdos, find_minimum_maximum
          
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
USE mp_images, ONLY : intra_image_comm
USE io_global, ONLY : ionode_id, ionode, stdout
USE mp,        ONLY : mp_bcast
IMPLICIT NONE
TYPE(phdos_type), INTENT(INOUT) :: phdos
CHARACTER(LEN=256), INTENT(IN) :: filename
INTEGER :: iunit, ios
INTEGER, PARAMETER :: ndivx=100000
REAL(DP), ALLOCATABLE :: nu(:), dos(:)
REAL(DP) :: de, de_
INTEGER :: i, ndiv
INTEGER :: find_free_unit

IF (ionode) THEN
   iunit=find_free_unit()
   OPEN(file=TRIM(filename), unit=iunit, status='old', &
       form='formatted', err=100, iostat=ios)
ENDIF
100 CALL mp_bcast(ios, ionode_id, intra_image_comm)
IF (ios /= 0) CALL errore('read_phdos_data', &
                          'opening file'//TRIM(filename), ABS(ios))

ALLOCATE(nu(ndivx))
ALLOCATE(dos(ndivx))
de = 0d0
IF (ionode) THEN
   DO i=1,ndivx
       ! nu(i) = frequencies (cm^{-1}), dos(i) in states/cm^{-1} 
      READ(iunit, *, END=20, ERR=10, IOSTAT=ios) nu(i),dos(i)
      IF ( nu(i) < -1.d0 ) THEN
         write(stdout,*) i, nu(i), dos(i)
         CALL errore('read_phdos_data','negative frequencies',1)
      ELSE IF ( nu(i) < 0.d0 ) THEN
         nu(i) = 0.d0
      END IF
      IF ( i ==2 ) de_ = nu(2) - nu(1)
      IF (i > 2) THEN
         de = nu(i) - nu(i-1)
         IF ( ABS(de - de_) > 1.0d-4 ) &
            CALL errore('read_phdos_data','nonuniform grid',1)
      END IF
      ndiv=i
   ENDDO
10 IF (ios /= 0 ) CALL errore('read_phdos_data', 'problem reading phdos', 1)
20 continue
ENDIF
CALL mp_bcast(ndiv,ionode_id,intra_image_comm)
CALL mp_bcast(de,ionode_id,intra_image_comm)
CALL mp_bcast(nu,ionode_id,intra_image_comm)
CALL mp_bcast(dos,ionode_id,intra_image_comm)

phdos%number_of_points=ndiv
phdos%de=de
ALLOCATE(phdos%nu(ndiv))
ALLOCATE(phdos%phdos(ndiv))
phdos%nu(:) = nu(1:ndiv)
phdos%phdos(:) = dos(1:ndiv)

DEALLOCATE(nu)
DEALLOCATE(dos)
IF (ionode) CLOSE(iunit)

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
REAL(DP) :: nu, arg, temp1, earg

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   arg= kb1 * nu * temp1
   earg = EXP( - arg )
   IF (nu >  0.0_DP) &
      ener = ener + phdos%phdos(i)* kb * temp * LOG( 1.0_DP - earg )
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
REAL(DP) :: nu, temp1, arg, earg

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   arg= kb1 * nu * temp1
   earg = EXP( -arg )
   IF (nu > 0.d0) ener = ener + phdos%phdos(i)* nu * earg/  & 
                                           ( 1.0_DP - earg ) 
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
REAL(DP) :: nu, temp1, arg, earg

cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   arg= kb1 * nu * temp1
   earg = EXP( - arg )
   IF (nu > 0.d0 ) cv = cv + phdos%phdos(i) * earg * &
                                        ( arg / ( 1.0_DP - earg )) ** 2 
ENDDO
cv = cv * phdos%de * kb

RETURN
END SUBROUTINE specific_heat_cv

SUBROUTINE fecv(phdos, temp, free_ener, ener, cv)
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational free energy at that temperature. ener contains
!  only the vibrational contribution WITHOUT the zero point energy. 
!  
!
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: free_ener, ener, cv

INTEGER :: ndiv, i
REAL(DP) :: nu, arg, temp1, earg, g

free_ener=0.0_DP
ener=0.0_DP
cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=phdos%number_of_points
DO i=1,ndiv
   nu=phdos%nu(i)
   g=phdos%phdos(i)
   arg= kb1 * nu * temp1
   earg = EXP( - arg )
   IF (nu > 0.0_DP ) THEN
      free_ener = free_ener + g * kb * temp * LOG( 1.0_DP - earg )
      ener = ener + g * nu * earg / ( 1.0_DP - earg )
      cv = cv + g * earg * ( arg / ( 1.0_DP - earg )) ** 2
   ENDIF
ENDDO
free_ener = free_ener*phdos%de
ener = ener * phdos%de / ry_to_cmm1
cv = cv * phdos%de * kb

RETURN
END SUBROUTINE fecv

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

SUBROUTINE find_minimum_maximum(phdos, freqmin, freqmax)
!
!  find the range of the phdos frequencies
!
IMPLICIT NONE
TYPE(phdos_type), INTENT(IN) :: phdos
REAL(DP), INTENT(OUT) :: freqmin, freqmax
INTEGER :: ndiv

ndiv=phdos%number_of_points
freqmin=phdos%nu(1)
freqmax=phdos%nu(ndiv)

RETURN
END SUBROUTINE find_minimum_maximum

END MODULE phdos_module
