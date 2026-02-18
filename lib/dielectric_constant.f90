!
! Copyright (C) 2025 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE dielectric_constant
!---------------------------------------------------------------------------
!
!   this module contains the support routines for the calculation
!   of the dielectric constant. Presently it contains routines to read
!   and write it on file.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE

  PUBLIC write_epsilon_infty_on_file, write_zeu_on_file,   &
         write_dielectric_properties_to_file,              &
         read_dielectric_properties_from_file,             &
         polar_mode_permittivity_tpw

CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE write_dielectric_properties_to_file(filename, nat, epsilon_infty, &
                                                           zeu)
!-------------------------------------------------------------------------
!
!  This routine writes the dielectric constant and the zeu effective
!  charge on file. These files can be read as the dielectric properties
!  of the equilibrium configurations.
!
!  The file contains the 3x3 matrix of the dielectric constant
!  The 3x3xnat file of the zeu Born effective charges for all atoms
!  The name of the file that contains these quantities is gives as 
!  input.
!
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: epsilon_infty(3,3), zeu(3,3,nat)

INTEGER :: find_free_unit
INTEGER :: outunit, ios, ipol, jpol, na

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', &
        FORM='formatted', ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_dielectric_properties_on_file',&
                           'ploblem opening output file', ABS(ios))

IF (ionode) THEN
   WRITE(outunit,'("Dielectric constant epsilon_infty")')
   DO ipol=1,3
      WRITE(outunit,'(3e19.10)') (epsilon_infty(ipol,jpol), jpol=1,3)
   ENDDO
   WRITE(outunit,*)
   DO na=1,nat
      WRITE(outunit,'("Born effective charge atom ",i5)') na
      DO ipol=1,3
         WRITE(outunit,'(6e19.10)') (zeu(ipol,jpol,na), jpol=1,3)
      ENDDO
   ENDDO
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_dielectric_properties_to_file
!
!-------------------------------------------------------------------------
SUBROUTINE read_dielectric_properties_from_file(filename, nat, epsilon_infty, &
                                                           zeu)
!-------------------------------------------------------------------------
!
!  This routine writes the dielectric constant and the zeu effective
!  charge on file. These files can be read as the dielectric properties
!  of the equilibrium configurations.
!
!  The file contains the 3x3 matrix of the dielectric constant
!  The 3x3xnat file of the zeu Born effective charges for all atoms
!  The name of the file that contains these quantities is gives as 
!  input.
!
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(OUT) :: epsilon_infty(3,3), zeu(3,3,nat)

INTEGER :: find_free_unit
INTEGER :: outunit, ios, ipol, jpol, na

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', &
        FORM='formatted', ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_dielectric_properties_on_file',&
                           'ploblem opening output file', ABS(ios))

IF (ionode) THEN
   READ(outunit,*)
   DO ipol=1,3
      READ(outunit,*) (epsilon_infty(ipol,jpol), jpol=1,3)
   ENDDO
   READ(outunit,*)
   DO na=1,nat
      READ(outunit,*) 
      DO ipol=1,3
         READ(outunit,*) (zeu(ipol,jpol,na), jpol=1,3)
      ENDDO
   ENDDO
   CLOSE(outunit)
ENDIF
CALL mp_bcast(zeu,ionode_id,intra_image_comm)
CALL mp_bcast(epsilon_infty,ionode_id,intra_image_comm)

RETURN
END SUBROUTINE read_dielectric_properties_from_file
!
!-------------------------------------------------------------------------
SUBROUTINE write_epsilon_infty_on_file(temp, ntemp, epsilon_infty_t, &
                                       astring, filename, iflag)
!-------------------------------------------------------------------------
!
!  iflag=0 writes the dielectric constant as a function of temperature
!  iflag=2 writes the dielectric constant as a function of pressure
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, iflag
REAL(DP), INTENT(IN) :: temp(ntemp), epsilon_infty_t(3,3,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: astring
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_epsilon, ios
INTEGER :: find_free_unit
CHARACTER(LEN=7) :: label

iu_epsilon=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_epsilon, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_epsilon_infty_on_file','opening dielectric (T) file',&
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

   WRITE(iu_epsilon,'(a)') astring
   WRITE(iu_epsilon,'("#",5x, a7, 13x, " e_11 ", 13x, " e_12",  &
                      &13x, " e_13", 13x, " e_22", 13x, " e_23",&
                      &13x, " e_33")') label

   DO itemp=2,ntemp-1
      WRITE(iu_epsilon,'(e16.8,6e20.12)') temp(itemp),             &
           epsilon_infty_t(1,1,itemp), epsilon_infty_t(1,2,itemp), &
           epsilon_infty_t(1,3,itemp), epsilon_infty_t(2,2,itemp), &
           epsilon_infty_t(2,3,itemp), epsilon_infty_t(3,3,itemp)
   ENDDO
   CLOSE(iu_epsilon)
ENDIF

RETURN
END SUBROUTINE write_epsilon_infty_on_file
!
!-------------------------------------------------------------------------
SUBROUTINE write_zeu_on_file(temp, ntemp, nat, ibrav, code_group, &
                                  zeu_t, filename, iflag)
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
INTEGER, INTENT(IN) :: ntemp, ibrav, nat, code_group, iflag
REAL(DP), INTENT(IN) :: temp(ntemp), zeu_t(3,3,nat,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_zeu, ios, na
INTEGER :: find_free_unit
CHARACTER(LEN=7) :: label
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename_tot

iu_zeu=find_free_unit()
DO na=1, nat
   filename_tot=TRIM(filename)//".at_"//TRIM(int_to_char(na))
   IF (meta_ionode) &
      OPEN(UNIT=iu_zeu, FILE=TRIM(filename_tot), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30    CALL mp_bcast(ios, meta_ionode_id, world_comm)
      CALL errore('write_zeu_on_file','opening zeu (T) file',&
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
      WRITE(iu_zeu,'("#    Born effective charges atom ",i5)') na
      WRITE(iu_zeu,'("#",5x, a7, 13x, " Z_11 ", 13x, " Z_12",      &
                         &13x, " Z_13", 13x, " Z_21", 13x, " e_22",&
                         &13x, " Z_23", 13x, " Z_31", 13x, " e_32",&
                         &13x, " e_33")') label

      DO itemp=2,ntemp-1
         WRITE(iu_zeu,'(e16.8,9e20.12)') temp(itemp),          &
              zeu_t(1,1,na,itemp), zeu_t(1,2,na,itemp),        &
              zeu_t(1,3,na,itemp), zeu_t(2,1,na,itemp),        &
              zeu_t(2,2,na,itemp), zeu_t(2,3,na,itemp),        &
              zeu_t(3,1,na,itemp), zeu_t(3,2,na,itemp), zeu_t(3,3,na,itemp)
      ENDDO
      CLOSE(iu_zeu)
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_zeu_on_file
!
!  This routine has been imported from QE 
!  Copyright (C) 2001-2016 Quantum ESPRESSO group
!
!  This routine has been simplified mainly to remove writing and
!  to pass eps0 as a output variable that can be used by other codes.
!
!----------------------------------------------------------------------
SUBROUTINE polar_mode_permittivity_tpw( nat, eps_infty, z, zstar, w2, omega, &
                                                        eps_new)
  !----------------------------------------------------------------------

  !
  ! Algorithm from Fennie and Rabe, Phys. Rev. B 68, 184111 (2003)
  !
  USE kinds, ONLY: DP
  USE constants, ONLY : pi, tpi, fpi, eps8, eps12, &
                       ELECTRON_SI, BOHR_RADIUS_SI, AMU_SI, C_SI, &
                       EPSNOUGHT_SI, AMU_RY, RY_TO_THZ, RY_TO_CMM1
  implicit none
  !number of atoms
  integer, intent(in) :: nat
 
  !electronic part of the permittivity
  real(DP), intent(in) :: eps_infty(3,3)
 
  !displacement eigenvectors
  complex(DP), intent(in) :: z(3*nat,3*nat)
 
  !born effective charges
  real(DP), intent(in) :: zstar(3,3,nat)
 
  !square of the phonon frequencies
  real(DP), intent(in) :: w2(3*nat)
 
  !cell volume
  real(DP), intent(in) :: omega
 
  !combined permittivity
  REAL(DP), INTENT(OUT) :: eps_new(3,3)

  !mode index
  integer :: imode
 
  !atom index
  integer :: iat
 
  !atom vector component index
  integer :: iat_component
 
  !Cartesian direction indices
  integer :: i, j
 
  !mode effective charge
  real(DP) :: meffc(3)
 
  !total effective plasma frequency
  real(DP) :: weff_tot, freq
 
  !polar mode contribution to the permittivity
  real(DP) :: deps(3,3)
 
 
  !Conversion factor for plasma frequency from Rydberg atomic units to SI
  real(DP) :: plasma_frequency_si 
  !Conversion factor for permittivity from Rydberg atomic units to SI
  real(DP) :: permittivity_si 

  ! some compiler do not like SQRT in initialization expressions
  plasma_frequency_si = ELECTRON_SI/sqrt(EPSNOUGHT_SI*BOHR_RADIUS_SI**3*AMU_SI)
  permittivity_si = plasma_frequency_si**2 / (fpi * pi)

  eps_new=eps_infty
  
  !Calculate contributions by mode
  DO imode = 1,3*nat
   
     ! Calculate the mode effective charge
     meffc = 0.0d0
     DO i = 1 , 3
        DO iat = 1 , nat
           DO j = 1, 3
              iat_component = 3*(iat-1) + j
 
              ! Equation (3) of Finnie and Rabe
              ! Rydberg units = (e / sqrt(2)) * Bohr
              meffc(i) = meffc(i) + zstar(i,j,iat)*z(iat_component,imode)* &
                                    sqrt(AMU_RY)
 
           END DO
        END DO
     END DO
 
     ! Calculate the polar mode contribution to the permittivity
     deps = 0.0d0
     ! Use only hard modes (frequency^2 > 10^-8 Ry)
     ! eps12 pass from 1/THz to 1/Hz
     !
     IF (w2(imode) > eps8) THEN
        DO i = 1 , 3
           DO j = 1 , 3
    
              ! Equation (2) of Finnie and Rabe
              deps(i,j) = (permittivity_si*eps12**2/omega)*meffc(i)*meffc(j)/&
                          (w2(imode)*RY_TO_THZ**2)
    
           END DO
        END DO
     END IF
    
     ! Add polar mode contribution to the total permittivity
     DO i = 1 , 3
        DO j = 1 , 3
           eps_new(i,j) = eps_new(i,j) + deps(i,j)
        END DO
     END DO
  END DO
  
RETURN 
END SUBROUTINE polar_mode_permittivity_tpw

END MODULE dielectric_constant
