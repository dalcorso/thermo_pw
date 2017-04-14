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
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE control_elastic_constants, ONLY : el_con_ibrav_geo
USE data_files,       ONLY : flanhar
USE postscript_files, ONLY : flpsanhar
USE anharmonic,       ONLY : lelastic
USE ph_freq_anharmonic,  ONLY : lelasticf
USE temperature,      ONLY : tmin, tmax
USE thermo_sym,       ONLY : laue
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filename, filenameps, filelastic
INTEGER :: ibrav

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
IF (laue==32.OR.laue==29) THEN
   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
ENDIF

IF (laue==2.OR.laue==18.OR.laue==19.OR.laue==20.OR.laue==22.OR.laue==23) THEN
!
!  tetrahonal, hexagonal or orthorhombic
!
   CALL gnuplot_ylabel('C_{13} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
ENDIF

IF (laue==18.OR.laue==22.OR.laue==19.OR.laue==23) THEN
!
!  tetragonal or hexagonal
!
   CALL gnuplot_ylabel('C_{33} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,6,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,6,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
ENDIF

IF (laue==18.OR.laue==22) THEN
!
!  tetragonal C_4h or D_4h
!
   CALL gnuplot_ylabel('C_{66} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,7,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,7,'color_blue',.NOT.lelastic,&
                                                         .TRUE.,.FALSE.)
ENDIF

IF (laue==18) THEN
!
!  tetragonal C_4h
!
   CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,8,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,8,'color_blue',.NOT.lelastic,&
                                                         .TRUE.,.FALSE.)
ENDIF

IF (laue==2.OR.laue==16.OR.laue==20) THEN
!
!  triclinic C_i, monoclinic C_2h, or orthorhombic D_2h
!
   CALL gnuplot_ylabel('C_{22} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{23} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,6,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,6,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{33} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,7,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,7,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,8,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,8,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{55} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,9,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,9,'color_blue',.NOT.lelastic,&
                                                     .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{66} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,10,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,10,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)
END IF

IF (laue==16) THEN
   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{15} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.) 
   ENDIF

   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,11,'color_red',.TRUE.,&
                                                  .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,11,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)
  
   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('C_{26} (kbar)',.FALSE.) 
   ENDIF

   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,12,'color_red',.TRUE.,&
                                                  .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,12,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)
  
   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{35} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('C_{36} (kbar)',.FALSE.) 
   ENDIF

   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,13,'color_red',.TRUE.,&
                                                  .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,13,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)
  

   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{46} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('C_{45} (kbar)',.FALSE.) 
   ENDIF

   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,14,'color_red',.TRUE.,&
                                                  .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,14,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)
ENDIF

IF (laue==2) THEN
   CALL gnuplot_ylabel('C_{14} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,11,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,11,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{15} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,12,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,12,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,13,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,13,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{24} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,14,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,14,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,15,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,15,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{26} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,16,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,16,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{34} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,17,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,17,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{35} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,18,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,18,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{36} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,19,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,19,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{45} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,20,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,20,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{46} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,21,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,21,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{56} (kbar)',.FALSE.) 
   IF (lelastic) &
      CALL gnuplot_write_file_mul_data(filelastic,1,22,'color_red',.TRUE.,&
                                                     .NOT.lelasticf,.FALSE.)
   IF (lelasticf) &
      CALL gnuplot_write_file_mul_data(filename,1,22,'color_blue', &
                                              .NOT.lelastic,.TRUE.,.FALSE.)

ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_elastic_t
