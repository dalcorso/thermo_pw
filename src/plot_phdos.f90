!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE plot_phdos()

USE kinds,      ONLY : DP
USE control_dosq,  ONLY : freqmin, freqmax
USE data_files, ONLY : fldos
USE postscript_files, ONLY : flpsdos

IMPLICIT NONE
CHARACTER(LEN=256) :: filedos

filedos="phdisp_files/"//TRIM(fldos)
CALL simple_plot('_dos', filedos, flpsdos, 'frequency (cm^{-1})', &
                'DOS (states / cm^{-1} / cell)', 'color_red', freqmin, &
                                   freqmax, 0.0_DP, 0.0_DP)
RETURN
END SUBROUTINE plot_phdos
