!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE plot_elastic_t(iflag, with_s)
!---------------------------------------------------------------------
!
!  This is a driver to plot the elastic constants as a function of
!  temperature
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : ibrav_geo
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar, fl_el_cons
USE postscript_files, ONLY : flpsanhar
USE control_elastic_constants,  ONLY : lelastic, lelasticf, lelastic_p, &
                             el_cons_t_available
USE control_grun,     ONLY : lb0_t
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : tmin, tmax
USE thermo_sym,       ONLY : laue
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
LOGICAL :: with_s
CHARACTER(LEN=256) :: gnu_filename, filenameps, filelastic, filelastic_s
INTEGER :: ibrav
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lelastic.OR.lelasticf.OR.lelastic_p).OR..NOT.lb0_t) RETURN

ibrav=ibrav_geo(1)

IF (MOD(iflag,2)==0) THEN
   IF (iflag==0) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_cons"
      filenameps=TRIM(flpsanhar)//".el_cons"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_cons"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_cons_s"
   ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//".el_cons_p"
      filenameps=TRIM(flpsanhar)//".el_cons_p"//TRIM(flext)
      filelastic="elastic_constants/"//TRIM(fl_el_cons)//".el_cons_p"
      IF (with_s) CALL errore('plot_elastic_t','Called in the wrong case',1)
   ENDIF
ELSE
   IF (iflag==1) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_comp"
      filenameps=TRIM(flpsanhar)//".el_comp"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_comp"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_comp_s"
   ELSE
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//".el_comp_p"
      filenameps=TRIM(flpsanhar)//".el_comp_p"//TRIM(flext)
      filelastic="elastic_constants/"//TRIM(fl_el_cons)//".el_comp_p"
      IF (with_s) CALL errore('plot_elastic_t','Called in the wrong case',1)
   ENDIF
ENDIF
CALL gnuplot_start(gnu_filename)

IF (iflag<2) THEN
   IF (tmin ==1._DP) THEN
      CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ELSE
      CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ENDIF
   CALL gnuplot_xlabel('T (K)', .FALSE.) 
ELSE
   CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, &
                                              0.0_DP, 0.0_DP, 1.0_DP, flext ) 
   CALL gnuplot_xlabel('p (kbar)', .FALSE.) 
ENDIF

IF (MOD(iflag,2)==0) THEN
   CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
   CALL plot_one_elastic_constant(1, 3, 'C_{11} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
ELSE
   CALL gnuplot_set_fact(1.D3, .FALSE.) 
   CALL plot_one_elastic_constant(1, 3, 'S_{11} (Mbar^{-1})', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
ENDIF

IF (MOD(iflag,2)==0) THEN
   CALL plot_one_elastic_constant(1, 4, 'C_{12} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
ELSE
   CALL plot_one_elastic_constant(1, 4, 'S_{12} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
ENDIF

IF (laue==32.OR.laue==29) THEN
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 5, 'C_{44} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 5, 'S_{44} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (laue==2.OR.laue==18.OR.laue==19.OR.laue==20.OR.laue==22.OR.laue==23.OR.&
    laue==25.OR.laue==27) THEN
!
!  tetragonal, hexagonal, trigonal, or orthorhombic
!
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 5, 'C_{13} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 5, 'S_{13} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (laue==18.OR.laue==22.OR.laue==19.OR.laue==23.OR.laue==25.OR.laue==27) THEN
!
!  tetragonal or hexagonal
!
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 6, 'C_{33} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 6, 'S_{33} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 7, 'C_{44} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 7, 'S_{44} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (laue==25.OR.laue==27) THEN
!
!  trigonal D_3d or S_6
!
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 8, 'C_{14} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 8, 'S_{14} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (laue==27) THEN
!
!  trigonal S_6
!
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 9, 'C_{25} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 9, 'S_{25} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (laue==18.OR.laue==22) THEN
!
!  tetragonal C_4h or D_4h
!
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 8, 'C_{66} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 8, 'S_{66} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (laue==18) THEN
!
!  tetragonal C_4h
!
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 9, 'C_{16} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 9, 'S_{16} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (laue==2.OR.laue==16.OR.laue==20) THEN
!
!  triclinic C_i, monoclinic C_2h, or orthorhombic D_2h
!
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 6, 'C_{22} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 6, 'S_{22} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 7, 'C_{23} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 7, 'S_{23} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 8, 'C_{33} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 8, 'S_{33} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 9, 'C_{44} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 9, 'S_{44} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 10, 'C_{55} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 10, 'S_{55} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 11, 'C_{66} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 11, 'S_{66} (Mbar^{-1})', lelastic, &
             lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
END IF

IF (laue==16) THEN
   IF (ibrav>0) THEN
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 12, 'C_{15} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 12, 'S_{15} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ELSE
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 12, 'C_{16} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 12, 'S_{16} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ENDIF

   IF (ibrav>0) THEN
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 13, 'C_{25} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 13, 'S_{25} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ELSE
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 13, 'C_{26} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 13, 'S_{26} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ENDIF
  
   IF (ibrav>0) THEN
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 14, 'C_{35} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 14, 'S_{35} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ELSE
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 14, 'C_{36} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 14, 'S_{36} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ENDIF

   IF (ibrav>0) THEN
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 15, 'C_{46} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 15, 'S_{46} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ELSE
      IF (MOD(iflag,2)==0) THEN
         CALL plot_one_elastic_constant(1, 15, 'C_{45} (kbar)', lelastic,  & 
             lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         CALL plot_one_elastic_constant(1, 15, 'S_{45} (Mbar^{-1})', &
         lelastic, lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ENDIF
ENDIF

IF (laue==2) THEN
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 12, 'C_{14} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 12, 'S_{14} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 13, 'C_{15} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 13, 'S_{15} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 14, 'C_{16} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 14, 'S_{16} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 15, 'C_{24} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 15, 'S_{24} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 16, 'C_{25} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 16, 'S_{25} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 17, 'C_{26} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 17, 'S_{26} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 18, 'C_{34} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 18, 'S_{34} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 19, 'C_{35} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 19, 'S_{35} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 20, 'C_{36} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 20, 'S_{36} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 21, 'C_{45} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 21, 'S_{45} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 22, 'C_{46} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 22, 'S_{46} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF

   IF (MOD(iflag,2)==0) THEN
      CALL plot_one_elastic_constant(1, 23, 'C_{56} (kbar)', lelastic,  & 
           lelasticf, filelastic, filelastic_s, with_s)
   ELSE
      CALL plot_one_elastic_constant(1, 23, 'S_{56} (Mbar^{-1})', lelastic, &
           lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
ENDIF

IF (MOD(iflag,2)==0) THEN
   CALL gnuplot_set_fact(1.0_DP, .FALSE.)
   CALL plot_one_elastic_constant(1, 2, 'Bulk modulus B (kbar)', lelastic,  & 
        lelasticf, filelastic, filelastic_s, with_s)
ELSE
   CALL gnuplot_set_fact(1.D3, .FALSE.)
   CALL plot_one_elastic_constant(1, 2, 'Compressibility K (Mbar^{-1})', &
        lelastic, lelasticf, filelastic, filelastic_s, with_s)
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_elastic_t

!---------------------------------------------------------------------
SUBROUTINE plot_one_elastic_constant(i, j, label, lelastic, lelasticf, &
                                           filelastic, filelastic_s, with_s)
!---------------------------------------------------------------------
USE gnuplot, ONLY : gnuplot_ylabel, &
                    gnuplot_write_file_mul_data
USE color_mod, ONLY : color
USE control_elastic_constants, ONLY : ngeom
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j
LOGICAL, INTENT(IN) :: lelastic, lelasticf, with_s
CHARACTER(LEN=*) :: label
CHARACTER(LEN=256), INTENT(IN) :: filelastic, filelastic_s

INTEGER :: igeom
CHARACTER(LEN=256) :: filename, filename_s
CHARACTER(LEN=6) :: int_to_char

INTEGER :: ic

IF (with_s.AND.ngeom>1) CALL errore('plot_one_elastic_constant',&
                      &'with_s and ngeom>1 are incompatible',1)
filename=TRIM(filelastic)//"_ph"
filename_s=TRIM(filelastic_s)//"_ph"

CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
IF (lelastic) THEN
   IF (ngeom==1) THEN
      CALL gnuplot_write_file_mul_data(filename,i,j,'color_red', .TRUE.,&
                                                      .NOT.with_s,.FALSE.)
   ELSE
      DO igeom=1, ngeom
         filename=TRIM(filelastic)//".g"//TRIM(int_to_char(igeom))
         ic=MOD(igeom-1,8)+1
         CALL gnuplot_write_file_mul_data(filename,i,j,color(ic),&
                   (igeom==1), (igeom==ngeom).AND..NOT.lelasticf,.FALSE.)
      END DO
   ENDIF
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filelastic_s,i,j,'color_green',.FALSE.,&
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   IF (ngeom==1) THEN
      CALL gnuplot_write_file_mul_data(filename,i,j,'color_blue', &
                                        .NOT.lelastic,.NOT.with_s,.FALSE.)
   ELSE
!
!   In this case with_s is false
!
      DO igeom=1, ngeom
         filename=TRIM(filelastic)//".g"//TRIM(int_to_char(igeom))//"_ph"
         ic=MOD(igeom-1,8)+1
         CALL gnuplot_write_file_mul_data(filename,i,j,color(ic), &
             .NOT.lelastic.AND.(igeom==1), (igeom==ngeom),.FALSE.)
      ENDDO
   ENDIF
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filename_s,i,j,'color_orange',.FALSE.,&
                                                     .TRUE.,.FALSE.)
ENDIF
IF (.NOT.(lelastic.OR.lelasticf)) THEN
   CALL gnuplot_write_file_mul_data(filelastic,i,j,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)
ENDIF

RETURN
END SUBROUTINE plot_one_elastic_constant
!
!---------------------------------------------------------------------
SUBROUTINE plot_elastic_t1(iflag, with_s)
!---------------------------------------------------------------------
!
!  This is a driver to plot the elastic constants (iflag=0,2),
!  or the elastic compliance (iflag=1,3) as a function of
!  temperature (iflag=0,1) or of pressure (iflag=2,3).
!  It is supposed to substitute plot_elastic_t.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : ibrav_geo
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE color_mod,        ONLY : color
USE data_files,       ONLY : flanhar, fl_el_cons
USE elastic_constants, ONLY : get_ec_type, ec_present, ect_names, ecm_names
USE postscript_files, ONLY : flpsanhar, flps_el_cons
USE control_elastic_constants,  ONLY : lelastic, lelasticf, lelastic_p, &
                             el_cons_t_available
USE control_grun,     ONLY : lb0_t
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : tmin, tmax
USE thermo_sym,       ONLY : laue
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
LOGICAL :: with_s
CHARACTER(LEN=256) :: gnu_filename, filenameps, filelastic, filelastic_s
INTEGER :: ibrav, ec_type, iec, ic, ic_last
INTEGER :: ierr, system
LOGICAL :: first, last

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lelastic.OR.lelasticf.OR.lelastic_p).OR..NOT.lb0_t) RETURN

ibrav=ibrav_geo(1)

IF (MOD(iflag,2)==0) THEN
   IF (iflag==0) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_cons"
      filenameps=TRIM(flpsanhar)//".el_cons"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_cons"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_cons_s"
   ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//".el_cons_p"
      filenameps=TRIM(flps_el_cons)//".el_cons_p"//TRIM(flext)
      filelastic="elastic_constants/"//TRIM(fl_el_cons)//".el_cons_p"
      IF (with_s) CALL errore('plot_elastic_t','Called in the wrong case',1)
   ENDIF
ELSE
   IF (iflag==1) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_comp"
      filenameps=TRIM(flpsanhar)//".el_comp"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_comp"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_comp_s"
   ELSE
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//".el_comp_p"
      filenameps=TRIM(flps_el_cons)//".el_comp_p"//TRIM(flext)
      filelastic="elastic_constants/"//TRIM(fl_el_cons)//".el_comp_p"
      IF (with_s) CALL errore('plot_elastic_t','Called in the wrong case',1)
   ENDIF
ENDIF
CALL gnuplot_start(gnu_filename)

IF (iflag<2) THEN
   IF (tmin ==1._DP) THEN
      CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ELSE
      CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
   ENDIF
   CALL gnuplot_xlabel('T (K)', .FALSE.) 
ELSE
   CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, &
                                              0.0_DP, 0.0_DP, 1.0_DP, flext ) 
   CALL gnuplot_xlabel('p (kbar)', .FALSE.) 
ENDIF

ec_type=get_ec_type(laue, ibrav)
ic=0
DO iec=1,21
   IF (ec_present(iec, ec_type)>0) THEN 
      IF (MOD(iflag,2)==0) THEN
         ic=ic+1
         CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
         CALL plot_one_elastic_constant(1, ec_present(iec, ec_type)+2,    &
                          ect_names(iec)//' (kbar)', lelastic,           & 
                          lelasticf, filelastic, filelastic_s, with_s)
      ELSE
         ic=ic+1
         CALL gnuplot_set_fact(1.D3, .FALSE.) 
         CALL plot_one_elastic_constant(1, ec_present(iec, ec_type)+2,    &
                         ecm_names(iec)//' (Mbar^{-1})', lelastic,       & 
                     lelasticf, filelastic, filelastic_s, with_s)
      ENDIF
   ENDIF
ENDDO
ic_last=ic

IF (MOD(iflag,2)==0) THEN
   CALL gnuplot_set_fact(1.0_DP, .FALSE.)
   CALL plot_one_elastic_constant(1, 2, 'Bulk modulus B (kbar)', lelastic,  &
        lelasticf, filelastic, filelastic_s, with_s)
ELSE
   CALL gnuplot_set_fact(1.D3, .FALSE.)
   CALL plot_one_elastic_constant(1, 2, 'Compressibility K (Mbar^{-1})', &
        lelastic, lelasticf, filelastic, filelastic_s, with_s)
ENDIF

IF (iflag>1) THEN
!
!  Here make a plot of all elastic constants or compliances in the same plot
!
   ic=0
   IF (MOD(iflag,2)==0) THEN
      CALL gnuplot_ylabel('C (kbar)',.FALSE.)
      CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
   ELSE
      CALL gnuplot_ylabel('S (Mbar^{-1})',.FALSE.)
      CALL gnuplot_set_fact(1.D3, .FALSE.) 
   ENDIF
   first=.TRUE.
   DO iec=1,21
      IF (ec_present(iec, ec_type)>0) THEN 
         ic=ic+1
         last=(ic==ic_last)
         CALL gnuplot_write_file_mul_data(filelastic,1,&
              ec_present(iec, ec_type)+2,color(MOD(ic,8)),first,last,.FALSE.)
         first=.FALSE.
      ENDIF
   ENDDO
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_elastic_t1
