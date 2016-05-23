! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_files_names(igeom)

USE postscript_files,  ONLY : flpstherm, flpsdisp, flpsdos
USE data_files,        ONLY : fltherm, flfrc, flfrq, fldos, flpband, flvec, &
                              fldosfrq
USE control_gnuplot,   ONLY : flgnuplot
USE internal_files_names, ONLY : fildyn_thermo, flfrc_thermo, flfrq_thermo, &
                       fldos_thermo, fltherm_thermo, flpband_thermo,       &
                       flpsdos_thermo, flpstherm_thermo, flgnuplot_thermo, &
                       flpsdisp_thermo, flvec_thermo, fldosfrq_thermo
USE output,            ONLY : fildyn
  !
IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=6) :: int_to_char

IF (igeom==1) fildyn_thermo="dynamical_matrices/"//TRIM(fildyn)
IF (igeom==1) flfrc_thermo=TRIM(flfrc)
IF (igeom==1) flfrq_thermo=TRIM(flfrq)
IF (igeom==1) flvec_thermo=TRIM(flvec)
IF (igeom==1) fldos_thermo=TRIM(fldos)
IF (igeom==1) fldosfrq_thermo=TRIM(fldosfrq)
IF (igeom==1) fltherm_thermo=TRIM(fltherm)
IF (igeom==1) flpband_thermo=TRIM(flpband)

IF (igeom==1) flgnuplot_thermo=TRIM(flgnuplot)

IF (igeom==1) flpsdos_thermo=TRIM(flpsdos)
IF (igeom==1) flpstherm_thermo=TRIM(flpstherm)
IF (igeom==1) flpsdisp_thermo=TRIM(flpsdisp)

fildyn=TRIM(fildyn_thermo)//'.g'//TRIM(int_to_char(igeom))//'.'
flfrc=TRIM(flfrc_thermo)//'.g'//TRIM(int_to_char(igeom))
flfrq=TRIM(flfrq_thermo)//'.g'//TRIM(int_to_char(igeom))
flvec=TRIM(flvec_thermo)//'.g'//TRIM(int_to_char(igeom))
fldos=TRIM(fldos_thermo)//'.g'//TRIM(int_to_char(igeom))
fldosfrq=TRIM(fldosfrq_thermo)//'.g'//TRIM(int_to_char(igeom))
fltherm=TRIM(fltherm_thermo)//'.g'//TRIM(int_to_char(igeom))
flpband=TRIM(flpband_thermo)//'.g'//TRIM(int_to_char(igeom))

flgnuplot=TRIM(flgnuplot_thermo)//'.g'//TRIM(int_to_char(igeom))

flpsdos=TRIM(flpsdos_thermo)//'.g'//TRIM(int_to_char(igeom))
flpstherm=TRIM(flpstherm_thermo)//'.g'//TRIM(int_to_char(igeom))
flpsdisp=TRIM(flpsdisp_thermo)//'.g'//TRIM(int_to_char(igeom))

RETURN
END SUBROUTINE set_files_names

!----------------------------------------------------------------------------
SUBROUTINE restore_files_names()
!----------------------------------------------------------------------------

USE postscript_files,  ONLY : flpstherm, flpsdisp, flpsdos
USE data_files,        ONLY : fltherm, flfrc, flfrq, fldos,  flpband, fldosfrq 
USE internal_files_names, ONLY : fildyn_thermo, flfrc_thermo, flfrq_thermo, &
                        fldos_thermo, fltherm_thermo, flpband_thermo,       &
                        flpsdos_thermo, flpstherm_thermo, flgnuplot_thermo, &
                        flpsdisp_thermo, fldosfrq_thermo
USE control_gnuplot,   ONLY : flgnuplot
USE output,            ONLY : fildyn
  !
IMPLICIT NONE

fildyn = TRIM(fildyn_thermo)
flfrc = TRIM(flfrc_thermo)
flfrq = TRIM(flfrq_thermo)
fldos = TRIM(fldos_thermo)
fldosfrq = TRIM(fldosfrq_thermo)
fltherm = TRIM(fltherm_thermo)
flpband = TRIM(flpband_thermo)

flgnuplot =TRIM(flgnuplot_thermo)

flpsdos = TRIM(flpsdos_thermo)
flpstherm = TRIM(flpstherm_thermo)
flpsdisp = TRIM(flpsdisp_thermo)

RETURN
END SUBROUTINE restore_files_names
