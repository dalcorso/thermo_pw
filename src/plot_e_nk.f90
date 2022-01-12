!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE plot_e_nk()
!---------------------------------------------------------------------
!
!  This is a driver to plot the quantities written inside flnkconv
!
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsnkconv
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_set_eref,            &
                             gnuplot_set_gfact,           &
                             gnuplot_write_file_mul_data
USE thermo_mod,       ONLY : energy_geo
USE data_files,       ONLY : flnkconv
USE control_conv,     ONLY : nk_test, nnk, nsigma
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char
CHARACTER(LEN=12), ALLOCATABLE :: color(:)
INTEGER :: isigma
INTEGER :: ierr, system
REAL(DP) :: xmin, xmax
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_nkconv'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsnkconv)//TRIM(flext)
xmin = nk_test(1,1)
xmax = nk_test(1,nnk)
CALL gnuplot_write_header(filename, xmin, xmax, 0.0_DP, 0.0_DP, 1.0_DP, flext ) 

CALL gnuplot_xlabel(' nk ',.FALSE.) 
CALL gnuplot_ylabel('Total energy error (mRy)',.FALSE.) 
CALL gnuplot_set_eref(energy_geo(nnk),.FALSE.) 
CALL gnuplot_set_gfact(1000._DP,.FALSE.) 

ALLOCATE(color(nsigma))
color='color_green'
color(1)='color_red'
IF (nsigma>1) color(nsigma)='color_blue'
DO isigma=1,nsigma
   first_step=(isigma==1)
   last_step=(isigma==nsigma)
   filename='energy_files/'//TRIM(flnkconv)//TRIM(int_to_char(isigma))&
                                                      //'/'//TRIM(flnkconv)
   CALL gnuplot_write_file_mul_data(filename,1,4,color(isigma),first_step, &
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
END SUBROUTINE plot_e_nk
