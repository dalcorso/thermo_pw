!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE plot_e_ke()
!--------------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside flkeconv
!  
!
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpskeconv
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_set_gfact,           &
                             gnuplot_set_eref,            &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_write_file_mul_data
USE data_files,      ONLY : flkeconv
USE control_conv,    ONLY : ke, nkeden, ncutoffene
USE thermo_mod,      ONLY : energy_geo
USE mp_images,       ONLY : my_image_id, root_image
USE io_global,       ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char
CHARACTER(LEN=12), ALLOCATABLE :: color(:)
LOGICAL :: first_step, last_step
INTEGER :: iden
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_keconv'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpskeconv)//TRIM(flext)
CALL gnuplot_write_header(filename, ke(1), ke(ncutoffene), 0.0_DP, &
                                                   0.0_DP, 1.0_DP, flext ) 
CALL gnuplot_xlabel('Kinetic energy (Ry)',.FALSE.) 
CALL gnuplot_ylabel('Total energy (mRy)',.FALSE.) 
CALL gnuplot_set_eref(energy_geo(ncutoffene),.FALSE.) 
CALL gnuplot_set_gfact(1000._DP,.FALSE.) 

ALLOCATE(color(nkeden))
color='color_green'
color(1)='color_red'
IF (nkeden>1) color(nkeden)='color_blue'
DO iden=1,nkeden
   first_step=(iden==1)
   last_step=(iden==nkeden)
   filename='energy_files/'//TRIM(flkeconv)//TRIM(int_to_char(iden))&
                                                //'/'//TRIM(flkeconv)
   CALL gnuplot_write_file_mul_data(filename,1,2,color(iden),first_step,&
                                          last_step,.FALSE.)
ENDDO
DEALLOCATE(color)
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_e_ke
