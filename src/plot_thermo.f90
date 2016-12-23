!
! Copyright (C) 2014-15 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_thermo()
!
!  This is a driver to plot the quantities written inside fltherm
!  
!
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE postscript_files, ONLY : flpstherm
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE data_files,      ONLY : fltherm
USE temperature,     ONLY : tmin, tmax
USE mp_images,       ONLY : root_image, my_image_id
USE io_global,       ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, filetherm
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_therm'
CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(flpstherm, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(flpstherm, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF
filetherm="therm_files/"//TRIM(fltherm)
filename=TRIM(filetherm)//'_ph'
CALL gnuplot_xlabel('T (K)', .FALSE.) 
CALL gnuplot_ylabel('Vibrational energy (kJ / (N mol))',.FALSE.) 
CALL gnuplot_set_fact(1313.3130_DP, .FALSE.) 

IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,2,'color_red',.TRUE.,&
                                                     .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue', &
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Vibrational free energy (kJ / (N mol))', .FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,3,'color_red',.TRUE.,&
                                                     .NOT.ltherm_freq, .FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',&
                                                .NOT.ltherm_dos,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1313313.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Entropy (J / K / (N mol))',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,4,'color_red',.TRUE., &
                                                   .NOT.ltherm_freq,.FALSE.)
IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.NOT.ltherm_dos,&
                                                             .TRUE.,.FALSE.)

CALL gnuplot_ylabel('Heat capacity C_v (J / K / (N mol))',.FALSE.) 
IF (ltherm_dos) &
   CALL gnuplot_write_file_mul_data(filetherm,1,5,'color_red',.TRUE.,&
                 .NOT.ltherm_freq,.FALSE.)

IF (ltherm_freq) &
   CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.NOT.ltherm_dos,&
                                                           .TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_thermo

SUBROUTINE plot_thermo_debye(igeom)
!
!  This is a driver to plot the quantities written inside fltherm_debye
!
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot
USE postscript_files, ONLY : flpstherm
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE data_files,       ONLY : fltherm
USE temperature,      ONLY : tmin, tmax
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: gnu_filename, filename, psfilename
CHARACTER(LEN=6) :: int_to_char
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_debye.g'//&
                                                   TRIM(int_to_char(igeom))
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpstherm)//'_debye.g'//TRIM(int_to_char(igeom))
IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(psfilename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(psfilename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF
filename='therm_files/'//TRIM(fltherm)//'_debye.g'//TRIM(int_to_char(igeom))
CALL gnuplot_xlabel('T (K)', .FALSE.) 
CALL gnuplot_ylabel('Debye vibrational energy (kJ / (N mol))',.FALSE.) 
CALL gnuplot_set_fact(1313.3130_DP, .FALSE.) 

CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Debye vibrational free energy (kJ / (N mol))', .FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',.TRUE.,.TRUE., .FALSE.)

CALL gnuplot_set_fact(1313313.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Debye entropy (J / K / (N mol))',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Debye heat capacity C_v (J / K / (N mol))',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_thermo_debye

SUBROUTINE plot_el_thermo()
!
!  This is a driver to plot the quantities written inside fltherm_el_thermo
!
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot
USE postscript_files, ONLY : flpseltherm
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE klist,            ONLY : degauss
USE ktetra,           ONLY : ltetra
USE data_files,       ONLY : fleltherm
USE temperature,      ONLY : tmin, tmax
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, psfilename
INTEGER :: system
INTEGER :: ierr

IF ( my_image_id /= root_image ) RETURN

IF (degauss==0.0_DP.AND..NOT.ltetra) RETURN
gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_eltherm'
CALL gnuplot_start(gnu_filename)

psfilename=TRIM(flpseltherm)
IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(psfilename, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(psfilename, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF
filename='therm_files/'//TRIM(fleltherm)
CALL gnuplot_xlabel('T (K)', .FALSE.) 
CALL gnuplot_ylabel('Electron energy (kJ / (N mol))',.FALSE.) 
CALL gnuplot_set_fact(1313.3130_DP, .FALSE.) 

CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Electron free energy (kJ / (N mol))', .FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',.TRUE.,.TRUE., .FALSE.)

CALL gnuplot_set_fact(1313313.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Electron entropy (J / K / (N mol))',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('Electron heat capacity C_v (J / K / (N mol))',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.TRUE.,.TRUE.,.FALSE.)
CALL gnuplot_set_fact(1.0_DP, .FALSE.) 

CALL gnuplot_ylabel('Electron chemical potential {/Symbol m} (Ry)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,1,6,'color_blue',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_el_thermo
