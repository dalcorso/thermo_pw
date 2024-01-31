! Copyright (C) 2024 Dal Corso Andrea
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE clean_pw_tpw( flag)
  !----------------------------------------------------------------------
  !    
  ! This is a temporary routine for QE7.3. I hope in the next versions
  ! it will not be needed.
  ! At the moment we need to clean the atomic tables, qrad tables and
  ! vloc tables that depend on the volume before running a new pw run, 
  ! otherwise we get wrong tables. These tables should not depend on
  ! the volume to be reusable.
  !
  USE rhoat_mod, ONLY : deallocate_tab_rhoat
  USE vloc_mod,  ONLY : deallocate_tab_vloc
  USE qrad_mod,  ONLY : deallocate_tab_qrad
  USE rhoc_mod,  ONLY : deallocate_tab_rhc
  IMPLICIT NONE
  LOGICAL :: flag

  CALL deallocate_tab_vloc()
  CALL deallocate_tab_qrad()
  CALL deallocate_tab_rhoat()
  CALL deallocate_tab_rhc()
  CALL clean_pw(flag)
  !
  RETURN
  !
END SUBROUTINE clean_pw_tpw
