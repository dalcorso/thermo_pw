!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE plot_form_factors()

USE kinds, ONLY : DP
USE ions_base, ONLY : nsp, atm
USE control_xrdp, ONLY : smin, smax, flpsformf, flformf 

IMPLICIT NONE
INTEGER :: it
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename, filenameps

DO it=1,nsp
   filename=TRIM(flformf)//'.'//TRIM(int_to_char(it))
   filenameps=TRIM(flpsformf)//'.'//TRIM(int_to_char(it))
   CALL simple_plot('_formf', filename, filenameps, 's = sin({/Symbol q}) /&
        &{/Symbol l}    (\305^{-1})', atm(it)//' atomic scattering factor f (s)', &
            'color_red', smin, smax, 0.0_DP, 0.0_DP)
END DO

RETURN
END SUBROUTINE plot_form_factors
