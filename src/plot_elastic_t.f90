!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_elastic_t(iflag, with_s)
!
!  This is a driver to plot the elastic constants as a function of
!  temperature
!
USE kinds,            ONLY : DP
USE thermo_mod,       ONLY : ibrav_geo
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                             gnuplot_ylabel, &
                             gnuplot_xlabel, &
                             gnuplot_write_file_mul_data, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar
USE postscript_files, ONLY : flpsanhar
USE anharmonic,       ONLY : lelastic
USE ph_freq_anharmonic,  ONLY : lelasticf
USE control_grun,     ONLY : lb0_t
USE temperature,      ONLY : tmin, tmax
USE thermo_sym,       ONLY : laue
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
LOGICAL :: with_s
CHARACTER(LEN=256) :: gnu_filename, filename, filenameps, filelastic, &
                      filelastic_s, filename_s
INTEGER :: ibrav
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lelastic.OR.lelasticf).OR..NOT.lb0_t) RETURN

ibrav=ibrav_geo(1)

IF (iflag==0) THEN
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_cons"
   filenameps=TRIM(flpsanhar)//".el_cons"//TRIM(flext)
   filelastic="anhar_files/"//TRIM(flanhar)//".el_cons"
   filelastic_s="anhar_files/"//TRIM(flanhar)//".el_cons_s"
ELSE
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_comp"
   filenameps=TRIM(flpsanhar)//".el_comp"//TRIM(flext)
   filelastic="anhar_files/"//TRIM(flanhar)//".el_comp"
   filelastic_s="anhar_files/"//TRIM(flanhar)//".el_comp_s"
ENDIF
CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
ENDIF
filename=TRIM(filelastic)//"_ph"
filename_s=TRIM(filelastic_s)//"_ph"

CALL gnuplot_xlabel('T (K)', .FALSE.) 
IF (iflag==0) THEN
   CALL gnuplot_ylabel('C_{11} (kbar)',.FALSE.) 
   CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
ELSE
   CALL gnuplot_ylabel('S_{11} (Mbar^{-1})',.FALSE.) 
   CALL gnuplot_set_fact(1.D3, .FALSE.) 
ENDIF

IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filelastic_s,1,3,'color_green',.FALSE.,&
                                                     .NOT.lelasticf,.FALSE.)
ENDIF

IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filename,1,3,'color_blue',.NOT.lelastic,&
                                                     .NOT.with_s,.FALSE.)
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filename_s,1,3,'color_orange',.FALSE.,&
                                                     .TRUE.,.FALSE.)
ENDIF
IF (iflag==0) THEN
   CALL gnuplot_ylabel('C_{12} (kbar)',.FALSE.) 
ELSE
   CALL gnuplot_ylabel('S_{12} (Mbar^{-1})',.FALSE.) 
ENDIF
IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filelastic_s,1,4,'color_green',.FALSE.,&
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filename,1,4,'color_blue',.NOT.lelastic,&
                                                     .NOT.with_s,.FALSE.)
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filename_s,1,4,'color_orange',.FALSE.,&
                                                     .TRUE.,.FALSE.)
ENDIF
IF (laue==32.OR.laue==29) THEN
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{44} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
      CALL gnuplot_write_file_mul_data(filelastic_s,1,5,'color_green',.FALSE.,&
                                                     .NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue', &
                                     .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,5,'color_orange',&
                                                      .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==2.OR.laue==18.OR.laue==19.OR.laue==20.OR.laue==22.OR.laue==23.OR.&
    laue==25.OR.laue==27) THEN
!
!  tetragonal, hexagonal, trigonal, or orthorhombic
!
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{13} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{13} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,5,'color_green',&
                                         .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue', &
                                  .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,5,'color_orange',&
                                                   .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==18.OR.laue==22.OR.laue==19.OR.laue==23.OR.laue==25.OR.laue==27) THEN
!
!  tetragonal or hexagonal
!
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{33} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{33} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,6,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,6,'color_green', &
                                 .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,6,'color_blue',&
                               .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,6,'color_orange', &
                                                 .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{44} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,7,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,7,'color_green', &
                                         .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,7,'color_blue', &
                                         .NOT.lelastic, .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,7,'color_orange', &
                                       .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==25.OR.laue==27) THEN
!
!  trigonal D_3d or S_6
!
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{14} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{14} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,8,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,8,'color_green', &
                                            .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,8,'color_blue',&
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,8,'color_orange', &
                                            .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==27) THEN
!
!  trigonal S_6
!
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{25} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,9,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,9,'color_green', &
                                            .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,9,'color_blue', &
                                      .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,9,'color_orange', &
                                                  .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==18.OR.laue==22) THEN
!
!  tetragonal C_4h or D_4h
!
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{66} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{66} (Mbar^{-1})',.FALSE.) 
   ENDIF

   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,8,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,8,'color_green',&
                                   .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,8,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,8,'color_orange', &
                                      .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==18) THEN
!
!  tetragonal C_4h
!
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{16} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,9,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,9,'color_green', &
                                     .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,9,'color_blue', &
                                     .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,9,'color_orange', &
                                                 .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==2.OR.laue==16.OR.laue==20) THEN
!
!  triclinic C_i, monoclinic C_2h, or orthorhombic D_2h
!
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{22} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{22} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,6,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,6,'color_green',&
                                          .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,6,'color_blue', &
                                .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,6,'color_orange', &
                                   .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{23} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{23} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,7,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,7,'color_green', &
                                   .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,7,'color_blue', &
                                  .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,7,'color_orange', &
                                      .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{33} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{33} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,8,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,8,'color_green',&
                                              .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,8,'color_blue', &
                                           .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,8,'color_orange', &
                                                      .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{44} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,9,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,9,'color_green', &
                                           .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,9,'color_blue', &
                                   .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,9,'color_orange', &
                               .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{55} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{55} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,10,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,10,'color_green',&
                               .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,10,'color_blue', &
                                 .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,10,'color_orange', &
                                      .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{66} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{66} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,11,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,11,'color_green',&
                                           .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,11,'color_blue', &
                                         .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,11,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
END IF

IF (laue==16) THEN
   IF (ibrav>0) THEN
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{15} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{15} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ELSE
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{16} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ENDIF

   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,12,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,12,'color_green', &
                                              .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,12,'color_blue', &
                                         .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,12,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (ibrav>0) THEN
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{25} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ELSE
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{26} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{26} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ENDIF

   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,13,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,13,'color_green', &
                                          .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,13,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,13,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
  
   IF (ibrav>0) THEN
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{35} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{35} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ELSE
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{36} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{36} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ENDIF

   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,14,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,14,'color_green', &
                                       .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,14,'color_blue', &
                                         .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,14,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF 

   IF (ibrav>0) THEN
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{46} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{46} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ELSE
      IF (iflag==0) THEN
         CALL gnuplot_ylabel('C_{45} (kbar)',.FALSE.) 
      ELSE
         CALL gnuplot_ylabel('S_{45} (Mbar^{-1})',.FALSE.) 
      ENDIF
   ENDIF

   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,15,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,15,'color_green',&
                                         .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,15,'color_blue', &
                                        .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,15,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==2) THEN
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{14} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{14} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,12,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,12,'color_green', &
                                    .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,12,'color_blue', &
                                         .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,12,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{15} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{15} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,13,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,13,'color_green',&
                                      .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,13,'color_blue', &
                                        .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,13,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{16} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,14,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,14,'color_green',&
                                          .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,14,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,14,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{24} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{24} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,15,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,15,'color_green', &
                                       .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,15,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,15,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{25} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,16,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,16,'color_green',&
                                  .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,16,'color_blue', &
                                         .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,16,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)

   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{26} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{26} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,17,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,17,'color_green',&
                               .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,17,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,17,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{34} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{34} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,18,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,18,'color_green',&
                                             .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,18,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,18,'color_blue', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{35} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{35} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,19,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,19,'color_green',&
                                      .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,19,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,19,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{36} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{36} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,20,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,20,'color_green',&
                           .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,20,'color_blue', &
                                         .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,20,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{45} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{45} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,21,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,21,'color_green',&
                                        .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,21,'color_blue', &
                                           .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,21,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF

   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{46} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{46} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,22,'color_red',.TRUE.,&
                                                  .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,22,'color_green',&
                                          .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,22,'color_blue', &
                                          .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,22,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
   IF (iflag==0) THEN
      CALL gnuplot_ylabel('C_{56} (kbar)',.FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S_{56} (Mbar^{-1})',.FALSE.) 
   ENDIF
   IF (lelastic) THEN
      CALL gnuplot_write_file_mul_data(filelastic,1,23,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filelastic_s,1,23,'color_green',&
                                             .FALSE.,.NOT.lelasticf,.FALSE.)
   ENDIF
   IF (lelasticf) THEN
      CALL gnuplot_write_file_mul_data(filename,1,23,'color_blue', &
                                         .NOT.lelastic,.NOT.with_s,.FALSE.)
      IF (with_s) &
         CALL gnuplot_write_file_mul_data(filename_s,1,23,'color_orange', &
                                              .FALSE.,.TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (iflag==0) THEN
   CALL gnuplot_ylabel('Bulk modulus B (kbar)',.FALSE.)
   CALL gnuplot_set_fact(1.0_DP, .FALSE.)
ELSE
   CALL gnuplot_ylabel('Compressibility K (Mbar^{-1})',.FALSE.)
   CALL gnuplot_set_fact(1.D3, .FALSE.)
ENDIF


IF (lelastic) THEN
   CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,&
                                                     .NOT.with_s,.FALSE.)
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filelastic_s,1,2,'color_green',.FALSE.,&
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   CALL gnuplot_write_file_mul_data(filename,1,2,'color_blue',.NOT.lelastic,&
                                                  .NOT.with_s,.FALSE.)
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filename_s,1,2,'color_orange',.FALSE.,&
                                                     .TRUE.,.FALSE.)
ENDIF
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_elastic_t
