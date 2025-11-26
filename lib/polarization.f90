!
! Copyright (C) 2025 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE polarization_vector
!---------------------------------------------------------------------------
!
!   this module contains the support routines for the calculation
!   of the polarization. Presently it contains only a routine to read
!   and write it on file.
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE

  INTEGER :: mod_tot

  PUBLIC write_polarization, read_polarization, mod_tot

CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE write_polarization(filename,polar,berry_phase)
!-------------------------------------------------------------------------
!
!  This routine writes the Berry phase and the polarization on file.
!  It must be called after computing the polarization
!  It saves: 
!  the spontaneous polarization of the current geometry
!
USE io_global, ONLY : ionode, ionode_id
USE constants, ONLY : electron_si, bohr_radius_si
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: polar(3), berry_phase(3)
REAL(DP) :: fact
INTEGER :: find_free_unit
INTEGER :: outunit, i, ios

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_polarization','ploblem opening output file', ABS(ios))

fact= electron_si / (bohr_radius_si)**2
IF (ionode) THEN
   WRITE(outunit,'("Berry phase calculated along the three &
                                          &reciprocal vectors")')
   WRITE(outunit,'(3e20.10)') (berry_phase(i), i=1,3)
   WRITE(outunit,*)
   WRITE(outunit,'("Spontaneous polarization (cartesian &
                                       &coordinates) (e/bohr**2)")')
   WRITE(outunit,'(3e20.10)') (polar(i), i=1,3)
   WRITE(outunit,*)
   WRITE(outunit,'("Spontaneous polarization (cartesian &
                                       &coordinates) (C/m**2)")')
   WRITE(outunit,'(3e20.10)') (polar(i)*fact, i=1,3)
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_polarization

!-------------------------------------------------------------------------
SUBROUTINE read_polarization(filename,polar,berry_phase)
!-------------------------------------------------------------------------
!
!  This routine writes the Berry phase and the polarization on file.
!  It must be called after computing the polarization
!  It saves: 
!  the spontaneous polarization of the current geometry
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(OUT) :: polar(3), berry_phase(3)
REAL(DP) :: polar_cm2(3)
INTEGER :: find_free_unit
INTEGER :: inunit, i, ios

IF (ionode) THEN
   inunit=find_free_unit()
   OPEN(UNIT=inunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('read_polarization','ploblem opening input file', ABS(ios))

IF (ionode) THEN
   READ(inunit,*)
   READ(inunit,'(3e20.10)') (berry_phase(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   READ(inunit,'(3e20.10)') (polar(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   READ(inunit,'(3e20.10)') (polar_cm2(i), i=1,3)
   CLOSE(inunit)
ENDIF
CALL mp_bcast(berry_phase,ionode_id,intra_image_comm)
CALL mp_bcast(polar,ionode_id,intra_image_comm)

RETURN
END SUBROUTINE read_polarization

END MODULE polarization_vector
