!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_elastic_t()
!
!  This is a driver to plot the elastic constants as a function of
!  temperature
!
USE kinds,            ONLY : DP
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot
USE control_thermo,   ONLY : ltherm_dos, ltherm_freq
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE postscript_files, ONLY : flpsanhar
USE anharmonic,       ONLY : lelastic
USE ph_freq_anharmonic,  ONLY : lelasticf
USE temperature,      ONLY : tmin, tmax
USE control_elastic_constants, ONLY : el_con_ibrav_geo
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, filenameps, filelastic
INTEGER :: system
INTEGER :: ierr, ibrav

IF ( my_image_id /= root_image ) RETURN

ibrav=el_con_ibrav_geo(1)

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_cons"
CALL gnuplot_start(gnu_filename)

filenameps=TRIM(flpsanhar)//".el_cons"
IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, 1.0_DP ) 
ENDIF
filelastic="anhar_files/"//TRIM(flanhar)//".el_cons"
filename=TRIM(filelastic)//"_ph"

CALL gnuplot_xlabel('T (K)', .FALSE.) 
CALL gnuplot_ylabel('C_{11} (kbar)',.FALSE.) 

CALL gnuplot_set_fact(1.0_DP, .FALSE.) 

IF (lelastic) &
   CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
IF (lelasticf) &
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)

CALL gnuplot_ylabel('C_{12} (kbar)',.FALSE.) 
IF (lelastic) &
   CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
IF (lelasticf) &
   CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
ENDIF

IF (ibrav==4) THEN
   CALL gnuplot_ylabel('C_{13} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{33} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,6,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,6,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,7,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,7,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_elastic_t
