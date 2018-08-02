!
! Copyright (C) 2015 Andrea Dal Corso
! Copyright (C) 2018 Cristiano Malica for plot_gen_phdos
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE plot_phdos()

USE kinds,            ONLY : DP
USE control_dosq,     ONLY : freqmin, freqmax
USE data_files,       ONLY : fldos
USE control_thermo,   ONLY : with_eigen
USE postscript_files, ONLY : flpsdos
USE control_gnuplot,  ONLY : flext

IMPLICIT NONE
CHARACTER(LEN=256) :: filedos, filepsdos

filedos="phdisp_files/"//TRIM(fldos)
filepsdos=TRIM(flpsdos)//TRIM(flext)
CALL simple_plot('_dos', filedos, filepsdos, 'frequency (cm^{-1})', &
                'DOS (states / cm^{-1} / cell)', 'color_red', freqmin, &
                                   freqmax, 0.0_DP, 0.0_DP)
IF (with_eigen) CALL plot_gen_phdos()

RETURN
END SUBROUTINE plot_phdos
!
SUBROUTINE plot_gen_phdos()
!
!  This is a driver to plot the generalized phdos
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE cell_base,        ONLY : ibrav
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE postscript_files, ONLY : flpsdos
USE control_dosq,     ONLY : freqmin, freqmax
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data
USE data_files,       ONLY : fldos
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, psfilename
INTEGER :: ierr, system, na
CHARACTER(LEN=6) :: int_to_char

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_gphdos'
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpsdos)//'_gphdos'//TRIM(flext)

CALL gnuplot_write_header(psfilename, freqmin, freqmax, 0.0_DP, 0.0_DP, &
                                                        1.0_DP, flext )

CALL gnuplot_xlabel('frequency (cm^{-1})', .FALSE.)

DO na=1,nat

   filename="phdisp_files/"//TRIM(fldos)//'.'//TRIM(int_to_char(na))
!
!   First the diagonal components
!
   CALL gnuplot_ylabel('g-DOS_{ii} (1 / cm^{-1}) (atom '// &
                                        TRIM(int_to_char(na))//')',.FALSE.)

   CALL gnuplot_write_file_mul_data(filename,1,2,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.FALSE., &
                                                   .FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data(filename,1,7,'color_green', .FALSE., &
                                                   .TRUE. ,.FALSE.)
!
!  And then the off diagonal, only for noncubic solids
!
   IF (ibrav/=1.AND.ibrav/=2.AND.ibrav/=3) THEN

      CALL gnuplot_ylabel('g-DOS_{ij} (1 / cm^{-1}) (atom '// &
                                        TRIM(int_to_char(na))//')',.FALSE.)

      CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE., &
                                                   .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.FALSE., &
                                                   .FALSE.,.FALSE.)
      CALL gnuplot_write_file_mul_data(filename,1,6, 'color_green',.FALSE., &
                                                   .TRUE.,.FALSE.)
   ENDIF
ENDDO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_gen_phdos
