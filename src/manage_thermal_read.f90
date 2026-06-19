!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE manage_thermal_read()
!-----------------------------------------------------------------------
!
!  This routine computes and writes on file the harmonic thermodynamic 
!  quantities calculated using the phonon dos
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : start_geometry, last_geometry
USE temperature,    ONLY : ntemp, temp
USE thermodynamics, ONLY : ph_ener, ph_free_ener, ph_entropy, ph_ce,     &
                           ph_b_fact, ph_e0
USE ph_freq_thermodynamics, ONLY : phf_ener, phf_free_ener, phf_entropy, &
                           phf_ce, phf_b_fact, phf_e0
USE control_thermo, ONLY : with_eigen, ltherm_dos, ltherm_freq

USE data_files,     ONLY : fltherm

IMPLICIT NONE

INTEGER :: igeom
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filetherm
LOGICAL  :: do_read

IF (ltherm_dos) THEN
   DO igeom=start_geometry, last_geometry

      filetherm='therm_files/'//TRIM(fltherm)//'.g'//TRIM(int_to_char(igeom))
      CALL read_thermo(ntemp, temp, ph_ener(1,igeom), ph_free_ener(1,igeom), &
                    ph_entropy(1,igeom), ph_ce(1,igeom), ph_e0(igeom), &
                                                        do_read, filetherm)
      IF (with_eigen.AND.do_read) &
         CALL read_b_factor(ntemp, ph_b_fact(1,1,1,1,igeom),filetherm) 
      IF (.NOT.do_read) CALL errore('manage_thermo_read','problem readin',&
                                            igeom)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO igeom=start_geometry, last_geometry
      filetherm='therm_files/'//TRIM(fltherm)//'.g'//&
                                TRIM(int_to_char(igeom))//'_ph'
      CALL read_thermo(ntemp, temp, phf_ener(1,igeom), phf_free_ener(1,igeom),&
                    phf_entropy(1,igeom), phf_ce(1,igeom), phf_e0(igeom), &
                    do_read, filetherm)
      IF (with_eigen.AND.do_read) &
         CALL read_b_factor(ntemp, phf_b_fact(1,1,1,1,igeom),filetherm) 
      IF (.NOT.do_read) CALL errore('manage_thermo_read','problem readin ph',&
                                            igeom)
   ENDDO
ENDIF

RETURN
END SUBROUTINE manage_thermal_read
