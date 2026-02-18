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
!  This routine call all the routines that plot quasi-anharmonic quantities
!  in the p-V thermodynamics
!
USE control_mur, ONLY : lmurn
IMPLICIT NONE

CALL plot_anhar_energy()
CALL plot_anhar_volume() 
CALL plot_anhar_density() 
CALL plot_anhar_press()
CALL plot_anhar_bulk() 
CALL plot_anhar_dbulk()
CALL plot_anhar_beta()
CALL plot_anhar_heat()
CALL plot_anhar_gamma()
CALL plot_anhar_thermo()
IF (lmurn) THEN
   CALL plot_anhar_dw()
!
   CALL plot_hugoniot()
!
   CALL plot_t_debye()
ENDIF
RETURN
END SUBROUTINE manage_plot_anhar
!
!---------------------------------------------------------------------
SUBROUTINE manage_plot_anhar_anis()
!---------------------------------------------------------------------
!
!  This routine calls all the routines that plot quasi-anharmonic quantities
!  in the stress-strain thermodynamics and then call the previous routine
!  to plot the p-V thermodynamic quantities
!
IMPLICIT NONE

CALL plot_anhar_anis_celldm()
CALL plot_anhar_anis_alpha()
CALL plot_anhar_anis_uint()
CALL plot_anhar_anis_uint_zsisa()
CALL plot_anhar_anis_dw()
CALL plot_thermal_stress()
CALL plot_generalized_gruneisen()
!
!  and here the p-V thermodynamic quantities
!
CALL manage_plot_anhar()

RETURN
END SUBROUTINE manage_plot_anhar_anis
!
!---------------------------------------------------------------------------
SUBROUTINE manage_plot_elastic()
!---------------------------------------------------------------------------
USE thermo_mod,   ONLY : what

IMPLICIT NONE
IF (what=='mur_lc_t') THEN
   CALL plot_elastic_t(0,.TRUE.)
   CALL plot_elastic_t(1,.TRUE.)
   CALL plot_elastic_d_t1(0,.FALSE.)
   CALL plot_elastic_d_t1(1,.FALSE.)
   CALL plot_elastic_pt(0,.TRUE.)
   CALL plot_elastic_pt(1,.TRUE.)
   CALL plot_elastic_d_pt(0,.TRUE.)
   CALL plot_elastic_d_pt(1,.TRUE.)
   CALL plot_elastic_ptt(0,.TRUE., .FALSE.)
   CALL plot_elastic_ptt(1,.TRUE., .FALSE.)
   CALL plot_elastic_d_ptt(0,.TRUE., .FALSE.)
   CALL plot_elastic_d_ptt(1,.TRUE., .FALSE.)
   CALL plot_elastic_ptt(0,.FALSE., .TRUE.)
   CALL plot_elastic_ptt(1,.FALSE., .TRUE.)

   CALL plot_macro_el_new_t()
   CALL plot_sound_t()
   CALL plot_macro_el_new_pt()
   CALL plot_sound_pt()
   CALL plot_macro_el_new_ptt()
   CALL plot_sound_ptt()
ENDIF

RETURN
END SUBROUTINE manage_plot_elastic
!
!---------------------------------------------------------------------------
SUBROUTINE manage_plot_piezo()
!---------------------------------------------------------------------------
USE thermo_mod,   ONLY : what

IMPLICIT NONE
IF (what=='mur_lc_t') THEN
   CALL plot_piezo_t(0)
   CALL plot_piezo_d_t(0)
   CALL plot_piezo_pt()
   CALL plot_piezo_d_pt()
   CALL plot_piezo_ptt()
   CALL plot_piezo_d_ptt()
   CALL plot_pyro_t(0)
ENDIF

RETURN
END SUBROUTINE manage_plot_piezo
