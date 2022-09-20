!
! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE manage_plot_anhar()
!---------------------------------------------------------------------
!
USE kinds,                 ONLY : DP

IMPLICIT NONE

CALL plot_anhar_energy()
CALL plot_anhar_volume() 
CALL plot_anhar_press()
CALL plot_anhar_bulk() 
CALL plot_anhar_dbulk()
CALL plot_anhar_beta()
CALL plot_anhar_heat()
CALL plot_anhar_gamma()
CALL plot_hugoniot()
!
CALL plot_thermo_anhar()
CALL plot_dw_anhar()
CALL plot_t_debye()

RETURN
END SUBROUTINE manage_plot_anhar

