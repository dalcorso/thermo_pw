!
! Copyright (C) 2020 Andrea Dal Corso and Cristiano Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE manage_el_free_energy
!---------------------------------------------------------------------
!
USE thermo_mod,           ONLY : tot_ngeo, start_geometry, last_geometry
USE control_thermo,       ONLY : outdir_thermo, lectqha
USE input_parameters,     ONLY : outdir

IMPLICIT NONE

INTEGER :: igeom
CHARACTER(LEN=6) :: int_to_char

DO igeom=start_geometry, last_geometry
   outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(igeom))//'/'
   CALL set_tmp_dir(outdir)
   CALL set_el_files_names(igeom)
   CALL dos_sub()
   CALL plot_dos()
   CALL write_el_thermo(igeom)
   CALL plot_el_thermo()
ENDDO

IF (.NOT.lectqha) CALL manage_el_anhar()

RETURN
END SUBROUTINE manage_el_free_energy
