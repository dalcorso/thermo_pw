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

REAL(DP), PARAMETER :: thr_ph=1.D-3   ! a phonon with frequency smaller than
                                      ! this is considered of zero frequency.
REAL(DP), PARAMETER :: thr_taylor=1.D-3   ! When the argument of the functions
                                          ! is smaller than this threshold use
                                          ! a Taylor expansion expression 
                                          ! to avoid the numerical error due
                                          ! to the calculation of 
                                          ! 1.0-(1.0-epsilon) 


TYPE phdos_type
   INTEGER :: number_of_points    ! rhe number of points
   REAL(DP) :: de                 ! interval of the mesh of frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: nu(:)     ! the frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: phdos(:)  ! the phdos (states/ cm-1)
END TYPE phdos_type

TYPE gen_phdos_type
   INTEGER :: number_of_points    ! the number of points
   INTEGER :: nat                 ! the number of atoms
   REAL(DP) :: de                 ! interval of the mesh of frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: nu(:)     ! the frequencies (cm-1)
   REAL(DP), ALLOCATABLE :: phdos(:,:,:)  ! the generalized phdos (states/ cm-1)
END TYPE gen_phdos_type

PUBLIC :: phdos_type, read_phdos_data, zero_point_energy, free_energy, &
          vib_energy, vib_entropy, specific_heat_cv, fecv, &
          phdos_debye_factor, integrated_dos, set_phdos, destroy_phdos, &
          find_minimum_maximum, gen_phdos_type, set_gen_phdos, &
          read_genphdos_data, destroy_gen_phdos
          
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

SUBROUTINE set_gen_phdos(phdos,ndiv,nat,deltanu)
IMPLICIT NONE
TYPE(gen_phdos_type), INTENT(INOUT) :: phdos
INTEGER, INTENT(IN) :: ndiv, nat
REAL(DP), INTENT(IN) :: deltanu

phdos%number_of_points=ndiv
phdos%nat=nat
phdos%de=deltanu
ALLOCATE(phdos%nu(ndiv))
ALLOCATE(phdos%phdos(6,nat,ndiv))

RETURN
END SUBROUTINE set_gen_phdos

SUBROUTINE destroy_phdos(phdos)
IMPLICIT NONE
TYPE(phdos_type), INTENT(INOUT) :: phdos

IF (ALLOCATED(phdos%nu)) DEALLOCATE(phdos%nu)
IF (ALLOCATED(phdos%phdos)) DEALLOCATE(phdos%phdos)

RETURN
END SUBROUTINE destroy_phdos

SUBROUTINE destroy_gen_phdos(phdos)
IMPLICIT NONE
TYPE(gen_phdos_type), INTENT(INOUT) :: phdos

IF (ALLOCATED(phdos%nu)) DEALLOCATE(phdos%nu)
IF (ALLOCATED(phdos%phdos)) DEALLOCATE(phdos%phdos)

RETURN
END SUBROUTINE destroy_gen_phdos

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

SUBROUTINE read_genphdos_data(phdos, nat, filename)
!
!  This subroutine reads the phdos from a file. It allocates space,
!  opens and closes the phdos file.
!
USE mp_images, ONLY : intra_image_comm
USE io_global, ONLY : ionode_id, ionode, stdout
USE mp,        ONLY : mp_bcast

IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
TYPE(gen_phdos_type), INTENT(INOUT) :: phdos
CHARACTER(LEN=256), INTENT(IN) :: filename
INTEGER :: iunit, ios
INTEGER, PARAMETER :: ndivx=100000
REAL(DP), ALLOCATABLE :: nu(:), dos(:,:,:)
REAL(DP) :: de, de_
INTEGER :: i, ijpol, na, ndiv
INTEGER :: find_free_unit
CHARACTER(LEN=256) :: filename_loc
CHARACTER(LEN=6) :: int_to_char

IF (ionode) iunit=find_free_unit()

ALLOCATE(nu(ndivx))
ALLOCATE(dos(6,nat,ndivx))
de = 0d0
DO na=1,nat
   filename_loc=TRIM(filename)//'.'//TRIM(int_to_char(na))
   IF (ionode) OPEN(FILE=TRIM(filename_loc), UNIT=iunit, STATUS='old', &
                       FORM='formatted', ERR=100, IOSTAT=ios)
100 CALL mp_bcast(ios, ionode_id, intra_image_comm)
   IF (ios /= 0) CALL errore('read_phdos_data', &
                          'opening file'//TRIM(filename_loc), ABS(ios))
   IF (ionode) THEN
      DO i=1,ndivx
    ! nu(i) = frequencies (cm^{-1}), dos(i) in states/cm^{-1} 
         READ(iunit, *, END=20, ERR=10, IOSTAT=ios) nu(i), &
                                                (dos(ijpol,na,i),ijpol=1,6)
         ndiv=i
      ENDDO
20    ios=0
   ENDIF
10 CALL mp_bcast(ios, ionode_id, intra_image_comm)
   IF (ios /= 0 ) CALL errore('read_genphdos_data', 'problem reading phdos', 1)
   IF (ionode) CLOSE(iunit)
ENDDO
CALL mp_bcast(ndiv,ionode_id,intra_image_comm)
CALL mp_bcast(nu,ionode_id,intra_image_comm)
CALL mp_bcast(dos,ionode_id,intra_image_comm)

DO i=1,ndiv
   IF ( nu(i) < -1.d0 ) THEN
      WRITE(stdout,*) i, nu(i), dos(ijpol,1,i)
      CALL errore('read_genphdos_data','negative frequencies',1)
   ELSE IF ( nu(i) < 0.d0 ) THEN
      nu(i) = 0.d0
   END IF
   IF ( i ==2 ) THEN
      de_ = nu(2) - nu(1)
   ELSEIF (i > 2) THEN
      de = nu(i) - nu(i-1)
      IF ( ABS(de - de_) > 1.0d-4 ) &
               CALL errore('read_genphdos_data','nonuniform grid',1)
   END IF
ENDDO

phdos%number_of_points=ndiv
phdos%de=de
ALLOCATE(phdos%nu(ndiv))
ALLOCATE(phdos%phdos(6,nat,ndiv))
phdos%nu(:) = nu(1:ndiv)
phdos%phdos(1:6,1:nat,1:ndiv) = dos(1:6,1:nat,1:ndiv)

DEALLOCATE(nu)
DEALLOCATE(dos)

RETURN
END SUBROUTINE read_genphdos_data

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

SUBROUTINE phdos_debye_factor(phdos, temp, b_fact,  &
                                              nat, amass, nsp, ityp)

!
USE constants, ONLY : h_planck_si, c_si, amu_si

IMPLICIT NONE 
TYPE(gen_phdos_type), INTENT(IN) :: phdos
INTEGER :: nsp, nat, ityp(nat)
REAL(DP), INTENT(IN) :: temp, amass(nsp)
REAL(DP), INTENT(INOUT) :: b_fact(3,3,nat)

INTEGER :: ndiv, i, na, ipol, jpol, ijpol
REAL(DP) :: nu, g, temp1, arg, expt, fact, tfact

b_fact = 0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv = phdos%number_of_points

fact = 1.D20 * h_planck_si / c_si / 100.0_DP / amu_si 

DO na=1,nat
   ijpol=0
   DO ipol=1,3
      DO jpol=ipol,3
         ijpol=ijpol+1
         DO i=1,ndiv
            nu = phdos%nu(i)
            g = phdos%phdos(ijpol,na,i)
            arg = kb1 * nu * temp1
            IF (arg > thr_taylor) THEN
               expt = EXP(-arg)
               tfact = (1.0_DP+expt)/(1.0_DP-expt)
               b_fact(ipol,jpol,na) = b_fact(ipol,jpol,na) + tfact * g / nu
            ELSEIF (nu > thr_ph) THEN
               b_fact(ipol,jpol,na) = b_fact(ipol,jpol,na) + (2.0_DP/arg &
                       + arg/6.0_DP - arg**3/360.0_DP) * g / nu
            ENDIF
         ENDDO
         IF (ipol/=jpol) b_fact(jpol,ipol,na)=b_fact(ipol,jpol,na)
      ENDDO
   ENDDO
   b_fact(:,:,na) = b_fact(:,:,na) * fact * phdos%de / amass(ityp(na))
ENDDO

RETURN
END SUBROUTINE phdos_debye_factor


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
