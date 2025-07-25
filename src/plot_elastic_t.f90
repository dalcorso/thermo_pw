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
!  temperature or pressure
!  iflag    0   elastic constants as a function of temperature
!           1   elastic compliances as a function of temperature
!           2   elastic constants as a function of pressure
!           3   elastic compliances as a function of pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact
USE data_files,       ONLY : flanhar, fl_el_cons
USE postscript_files, ONLY : flpsanhar
USE initial_conf,     ONLY : ibrav_save
USE control_elastic_constants,  ONLY : lelastic, lelasticf, lelastic_p
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

ibrav=ibrav_save

IF (MOD(iflag,2)==0) THEN
   IF (iflag==0) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_cons"
      filenameps=TRIM(flpsanhar)//".el_cons"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_cons"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_cons_s"
   ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_cons_p"
      filenameps=TRIM(flpsanhar)//".el_cons_p"//TRIM(flext)
      filelastic="elastic_constants/"//TRIM(fl_el_cons)//".el_cons_p"
      IF (with_s) CALL errore('plot_elastic_t','Called in the wrong case',1)
   ENDIF
ELSE
   IF (iflag==1) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_comp"
      filenameps=TRIM(flpsanhar)//".el_comp"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_comp"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_comp_s"
   ELSE
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_comp_p"
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
   IF (laue==32.OR.laue==29) THEN
      CALL plot_one_elastic_constant(1, 6, 'Shear anisotropy A', lelastic,  & 
        lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
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
USE thermo_mod, ONLY : what
USE gnuplot, ONLY : gnuplot_ylabel, &
                    gnuplot_write_file_mul_data
USE color_mod, ONLY : color
USE control_elastic_constants, ONLY : ngeom, all_geometry_done_geo
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, j
LOGICAL, INTENT(IN) :: lelastic, lelasticf, with_s
CHARACTER(LEN=*) :: label
CHARACTER(LEN=256), INTENT(IN) :: filelastic, filelastic_s

INTEGER :: igeom, first_ngeom, last_ngeom
CHARACTER(LEN=256) :: filename, filename_s, filename_ph, filename_s_ph
CHARACTER(LEN=6) :: int_to_char

INTEGER :: ic

IF (with_s.AND.ngeom>1) CALL errore('plot_one_elastic_constant',&
                      &'with_s and ngeom>1 are incompatible',1)

filename=TRIM(filelastic)
filename_s=TRIM(filelastic_s)
filename_ph=TRIM(filelastic)//'_ph'
filename_s_ph=TRIM(filelastic_s)//'_ph'

IF (lelastic.OR.lelasticf) THEN
   last_ngeom=1
   first_ngeom=0
   DO igeom=1,ngeom
      IF (all_geometry_done_geo(igeom).AND.first_ngeom==0) first_ngeom=igeom
      IF (all_geometry_done_geo(igeom)) last_ngeom=igeom
   ENDDO
ENDIF

CALL gnuplot_ylabel(TRIM(label),.FALSE.) 
IF (lelastic) THEN
   IF (what/='elastic_constants_geo') THEN
      CALL gnuplot_write_file_mul_data(filename,i,j,'color_red', .TRUE.,&
                          .NOT.with_s.AND..NOT.lelasticf,.FALSE.)
   ELSE
      DO igeom=1, ngeom
         IF (.NOT.all_geometry_done_geo(igeom)) CYCLE
         filename=TRIM(filelastic)//".g"//TRIM(int_to_char(igeom))
         ic=MOD(igeom-1,8)+1
         CALL gnuplot_write_file_mul_data(filename,i,j,color(ic),&
                   (igeom==first_ngeom), (igeom==last_ngeom).AND. &
                                             .NOT.lelasticf,.FALSE.)
      END DO
   ENDIF
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filelastic_s,i,j,'color_green',.FALSE.,&
                                                     .NOT.lelasticf,.FALSE.)
ENDIF
IF (lelasticf) THEN
   IF (what/='elastic_constants_geo') THEN
      CALL gnuplot_write_file_mul_data(filename_ph,i,j,'color_blue', &
                                        .NOT.lelastic,.NOT.with_s,.FALSE.)
   ELSE
!
!   In this case with_s is false
!
      DO igeom=1, ngeom
         IF (.NOT.all_geometry_done_geo(igeom)) CYCLE
         filename_ph=TRIM(filelastic)//".g"//TRIM(int_to_char(igeom))//"_ph"
         ic=MOD(igeom-1,8)+1
         CALL gnuplot_write_file_mul_data(filename_ph,i,j,color(ic), &
             .NOT.lelastic.AND.(igeom==first_ngeom), (igeom==last_ngeom),&
             .FALSE.)
      ENDDO
   ENDIF
   IF (with_s) &
      CALL gnuplot_write_file_mul_data(filename_s_ph,i,j,'color_orange', &
                                        .FALSE.,.TRUE.,.FALSE.)
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
USE control_elastic_constants,  ONLY : lelastic, lelasticf, lelastic_p
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
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_cons"
      filenameps=TRIM(flpsanhar)//".el_cons"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_cons"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_cons_s"
   ELSE
!
!  In this case plot the elastic constants as a function of pressure
!
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_cons_p"
      filenameps=TRIM(flps_el_cons)//".el_cons_p"//TRIM(flext)
      filelastic="elastic_constants/"//TRIM(fl_el_cons)//".el_cons_p"
      IF (with_s) CALL errore('plot_elastic_t1','Called in the wrong case',1)
   ENDIF
ELSE
   IF (iflag==1) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_comp"
      filenameps=TRIM(flpsanhar)//".el_comp"//TRIM(flext)
      filelastic="anhar_files/"//TRIM(flanhar)//".el_comp"
      filelastic_s="anhar_files/"//TRIM(flanhar)//".el_comp_s"
   ELSE
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_el_comp_p"
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
   IF (laue==32.OR.laue==29) THEN
      CALL plot_one_elastic_constant(1, 6, 'Shear anisotropy A', lelastic,  &
        lelasticf, filelastic, filelastic_s, with_s)
   ENDIF
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

!---------------------------------------------------------------------
SUBROUTINE plot_elastic_pt(iflag, with_s)
!---------------------------------------------------------------------
!
!  This is a driver to plot the elastic constants (iflag=0),
!  or the elastic compliance (iflag=1) as a function of
!  temperature at a set of pressures.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : ibrav_geo
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE color_mod,        ONLY : color
USE data_files,       ONLY : flanhar, fl_el_cons
USE elastic_constants, ONLY : get_ec_type, ec_present, ect_names, ecm_names
USE postscript_files, ONLY : flpsanhar, flps_el_cons
USE control_elastic_constants,  ONLY : lelastic, lelasticf, lelastic_p
USE control_grun,     ONLY : lb0_t
USE control_pressure, ONLY : pmin, pmax, npress_plot, ipress_plot, press
USE temperature,      ONLY : tmin, tmax
USE color_mod,        ONLY : color
USE thermo_sym,       ONLY : laue
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
LOGICAL :: with_s
CHARACTER(LEN=256) :: gnu_filename, filenameps, filelastic, filelastic_s, &
                      filelastic_ph, filelastic_s_ph, label
INTEGER :: ibrav, ec_type, iec, ic, ic_last, istep, ipressp, ipress
INTEGER :: ierr, system
LOGICAL :: first, last, first_step, last_step

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lelastic.OR.lelasticf.OR.lelastic_p).OR..NOT.lb0_t) RETURN
IF (npress_plot==0) RETURN

ibrav=ibrav_geo(1)

IF (MOD(iflag,2)==0) THEN
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_cons_p"
   filenameps=TRIM(flpsanhar)//".el_cons_p"//TRIM(flext)
ELSE
   gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_comp_p"
   filenameps=TRIM(flpsanhar)//".el_comp_p"//TRIM(flext)
ENDIF
CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                    1.0_DP, flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                    1.0_DP, flext ) 
ENDIF
CALL gnuplot_xlabel('T (K)', .FALSE.) 

ec_type=get_ec_type(laue, ibrav)
DO iec=1,21
   IF (ec_present(iec, ec_type)>0) THEN 
      istep=0
      IF (MOD(iflag,2)==0) THEN
         DO ipressp=1, npress_plot
            first_step=(ipressp==1)
            last_step=(ipressp==npress_plot)
            ipress=ipress_plot(ipressp)
            istep=MOD(istep,8)+1
            filelastic="anhar_files/"//TRIM(flanhar)//'.el_cons_press'
            CALL add_value(filelastic,press(ipress))
            filelastic_s="anhar_files/"//TRIM(flanhar)//'.el_cons_s_press'
            CALL add_value(filelastic_s,press(ipress))
            filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_ph_press'
            CALL add_value(filelastic_ph,press(ipress))
            filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_s_ph_press'
            CALL add_value(filelastic_s_ph,press(ipress))
            IF (first_step) THEN
               CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
               WRITE(label,'(a," (kbar)")') TRIM(ect_names(iec)) 
               CALL gnuplot_ylabel(TRIM(label),.FALSE.)
            ENDIF
            IF (ltherm_dos) THEN
               CALL gnuplot_write_file_mul_data(filelastic,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          first_step,((.NOT.with_s).AND.(last_step.AND.&
                            .NOT.ltherm_freq)),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          .FALSE.,(last_step.AND..NOT.ltherm_freq),.FALSE.)
            ENDIF
            IF (ltherm_freq) THEN
               CALL gnuplot_write_file_mul_data(filelastic_ph,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),&
                          ((.NOT.with_s).AND.last_step),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          .FALSE.,last_step,.FALSE.)
            ENDIF

         ENDDO
      ELSE
         DO ipressp=1, npress_plot
            first_step=(ipressp==1)
            last_step=(ipressp==npress_plot)
            ipress=ipress_plot(ipressp)
            istep=MOD(istep,8)+1
            filelastic="anhar_files/"//TRIM(flanhar)//'.el_comp_press'
            CALL add_value(filelastic,press(ipress))
            filelastic_s="anhar_files/"//TRIM(flanhar)//'.el_comp_s_press'
            CALL add_value(filelastic_s,press(ipress))
            filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_comp_ph_press'
            CALL add_value(filelastic_ph,press(ipress))
            filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.el_comp_s_ph_press'
            CALL add_value(filelastic_s_ph,press(ipress))
            IF (first_step) THEN
               CALL gnuplot_set_fact(1.D3, .FALSE.) 
               WRITE(label,'(a," (Mbar^{-1})")') TRIM(ecm_names(iec)) 
               CALL gnuplot_ylabel(TRIM(label),.FALSE.)
            ENDIF
            IF (ltherm_dos) THEN
               CALL gnuplot_write_file_mul_data(filelastic,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          first_step,((.NOT.with_s).AND.(last_step.AND.&
                             .NOT.ltherm_freq)),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s,1, &
                          ec_present(iec, ec_type)+2,color(istep), &
                          .FALSE.,(last_step.AND..NOT.ltherm_freq),.FALSE.)
            ENDIF
            IF (ltherm_freq) THEN
               CALL gnuplot_write_file_mul_data(filelastic_ph,1,   &
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),        &
                          ((.NOT.with_s).AND.last_step),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          .FALSE.,last_step,.FALSE.)
            ENDIF
         ENDDO
      ENDIF
   ENDIF
ENDDO

istep=0
IF (MOD(iflag,2)==0) THEN
   DO ipressp=1, npress_plot
      first_step=(ipressp==1)
      last_step=(ipressp==npress_plot)
      ipress=ipress_plot(ipressp)
      istep=MOD(istep,8)+1
      filelastic="anhar_files/"//TRIM(flanhar)//'.el_cons_press'
      CALL add_value(filelastic,press(ipress))
      filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_ph_press'
      CALL add_value(filelastic_ph,press(ipress))
      IF (first_step) THEN
         CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
         WRITE(label,'("Bulk modulus B (kbar)")') 
         CALL gnuplot_ylabel(TRIM(label),.FALSE.)
      ENDIF
      IF (ltherm_dos) THEN
         CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep), &
                          first_step,(last_step.AND..NOT.ltherm_freq),.FALSE.)
      ENDIF
      IF (ltherm_freq) THEN
         CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep), &
                           (first_step.AND..NOT.ltherm_dos),last_step,.FALSE.)
      ENDIF
   ENDDO
ELSE
   DO ipressp=1, npress_plot
      first_step=(ipressp==1)
      last_step=(ipressp==npress_plot)
      ipress=ipress_plot(ipressp)
      istep=MOD(istep,8)+1
      filelastic="anhar_files/"//TRIM(flanhar)//'.el_comp_press'
      CALL add_value(filelastic,press(ipress))
      filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_comp_ph_press'
      CALL add_value(filelastic_ph,press(ipress))
      IF (first_step) THEN
         CALL gnuplot_set_fact(1.D3, .FALSE.) 
         WRITE(label,'("Compressibility K (Mbar^{-1})")') 
         CALL gnuplot_ylabel(TRIM(label),.FALSE.)
      ENDIF
      IF (ltherm_dos) THEN
         CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep), &
                          first_step,(last_step.AND..NOT.ltherm_freq),.FALSE.)
      ENDIF
      IF (ltherm_freq) THEN
         CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),last_step,.FALSE.)
      ENDIF
   ENDDO
ENDIF

istep=0
IF (MOD(iflag,2)==0.AND.(laue==32.OR.laue==29)) THEN
   DO ipressp=1, npress_plot
      first_step=(ipressp==1)
      last_step=(ipressp==npress_plot)
      ipress=ipress_plot(ipressp)
      istep=MOD(istep,8)+1
      filelastic="anhar_files/"//TRIM(flanhar)//'.el_cons_press'
      CALL add_value(filelastic,press(ipress))
      filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_ph_press'
      CALL add_value(filelastic_ph,press(ipress))
      IF (first_step) THEN
         CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
         WRITE(label,'("Shear anisotropy A")') 
         CALL gnuplot_ylabel(TRIM(label),.FALSE.)
      ENDIF
      IF (ltherm_dos) THEN
         CALL gnuplot_write_file_mul_data(filelastic,1,6,color(istep), &
                          first_step,(last_step.AND..NOT.ltherm_freq),.FALSE.)
      ENDIF
      IF (ltherm_freq) THEN
         CALL gnuplot_write_file_mul_data(filelastic_ph,1,6,color(istep), &
                           (first_step.AND..NOT.ltherm_dos),last_step,.FALSE.)
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
END SUBROUTINE plot_elastic_pt
!
!---------------------------------------------------------------------
SUBROUTINE plot_elastic_ptt(iflag, with_t, with_s)
!---------------------------------------------------------------------
!
!  This is a driver to plot the elastic constants (iflag=0),
!  or the elastic compliance (iflag=1) as a function of
!  pressure at a set of temperatures.
!  It can print the isothermal elastic constantsi (with_t=.TRUE.), 
!  the isoentropic one (with_s=.TRUE.) 
!  or both (with_t=.TRUE., with_s=.TRUE.). 
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : ibrav_geo
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE color_mod,        ONLY : color
USE data_files,       ONLY : flanhar, fl_el_cons
USE elastic_constants, ONLY : get_ec_type, ec_present, ect_names, ecm_names
USE postscript_files, ONLY : flpsanhar, flps_el_cons
USE control_elastic_constants,  ONLY : lelastic, lelasticf, lelastic_p
USE control_grun,     ONLY : lb0_t
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : temp, ntemp_plot, itemp_plot
USE color_mod,        ONLY : color
USE thermo_sym,       ONLY : laue
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
INTEGER :: iflag
LOGICAL :: with_t, with_s
CHARACTER(LEN=256) :: gnu_filename, filenameps, filelastic, &
                      filelastic_ph, filelastic_s, filelastic_s_ph, label
INTEGER :: ibrav, ec_type, iec, ic, ic_last, istep, itempp, itemp
INTEGER :: ierr, system
LOGICAL :: first, last, first_step, last_step

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lelastic.OR.lelasticf.OR.lelastic_p).OR..NOT.lb0_t) RETURN
IF (ntemp_plot==0) RETURN

ibrav=ibrav_geo(1)

IF (MOD(iflag,2)==0) THEN
   IF (with_t) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_cons_t"
      filenameps=TRIM(flpsanhar)//".el_cons_t"//TRIM(flext)
   ELSEIF (with_s) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_cons_s_t"
      filenameps=TRIM(flpsanhar)//".el_cons_s_t"//TRIM(flext)
   ENDIF
   IF (with_t.AND.with_s) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_cons_st_t"
      filenameps=TRIM(flpsanhar)//".el_cons_st_t"//TRIM(flext)
   ENDIF
ELSE
   IF (with_t) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_comp_t"
      filenameps=TRIM(flpsanhar)//".el_comp_t"//TRIM(flext)
   ELSEIF(with_s) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_comp_s_t"
      filenameps=TRIM(flpsanhar)//".el_comp_s_t"//TRIM(flext)
   ENDIF
   IF (with_t.AND.with_s) THEN
      gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_anhar_el_comp_st_t"
      filenameps=TRIM(flpsanhar)//".el_comp_st_t"//TRIM(flext)
   ENDIF
ENDIF
CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, 0.0_DP, &
                              0.0_DP, 1.0_DP, flext ) 
CALL gnuplot_xlabel('p (kbar)', .FALSE.) 

ec_type=get_ec_type(laue, ibrav)
DO iec=1,21
   IF (ec_present(iec, ec_type)>0) THEN 
      istep=0
      IF (MOD(iflag,2)==0) THEN
         DO itempp=1, ntemp_plot
            first_step=(itempp==1)
            last_step=(itempp==ntemp_plot)
            itemp=itemp_plot(itempp)
            istep=MOD(istep,8)+1
            filelastic="anhar_files/"//TRIM(flanhar)//'.el_cons_temp'
            CALL add_value(filelastic,temp(itemp))
            filelastic_s="anhar_files/"//TRIM(flanhar)//'.el_cons_s_temp'
            CALL add_value(filelastic_s,temp(itemp))
            filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_ph_temp'
            CALL add_value(filelastic_ph,temp(itemp))
            filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_s_ph_temp'
            CALL add_value(filelastic_s_ph,temp(itemp))
            IF (first_step) THEN
               CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
               WRITE(label,'(a," (kbar)")') TRIM(ect_names(iec)) 
               CALL gnuplot_ylabel(TRIM(label),.FALSE.)
            ENDIF
            IF (ltherm_dos) THEN
               IF (with_t) &
                  CALL gnuplot_write_file_mul_data(filelastic,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          first_step,((.NOT.with_s).AND.(last_step.AND.&
                           .NOT.ltherm_freq)),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.with_t),            &
                          (last_step.AND..NOT.ltherm_freq),.FALSE.)
            ENDIF
            IF (ltherm_freq) THEN
               IF (with_t) &
                  CALL gnuplot_write_file_mul_data(filelastic_ph,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),&
                          ((.NOT.with_s).AND.last_step),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.ltherm_dos.AND..NOT.with_t),&
                           last_step,.FALSE.)
            ENDIF
         ENDDO
      ELSE
         DO itempp=1, ntemp_plot
            first_step=(itempp==1)
            last_step=(itempp==ntemp_plot)
            itemp=itemp_plot(itempp)
            istep=MOD(istep,8)+1
            filelastic="anhar_files/"//TRIM(flanhar)//'.el_comp_temp'
            CALL add_value(filelastic,temp(itemp))
            filelastic_s="anhar_files/"//TRIM(flanhar)//'.el_comp_s_temp'
            CALL add_value(filelastic_s,temp(itemp))
            filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_comp_ph_temp'
            CALL add_value(filelastic_ph,temp(itemp))
            filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.el_comp_s_ph_temp'
            CALL add_value(filelastic_s_ph,temp(itemp))
            IF (first_step) THEN
               CALL gnuplot_set_fact(1.D3, .FALSE.) 
               WRITE(label,'(a," (Mbar^{-1})")') TRIM(ecm_names(iec)) 
               CALL gnuplot_ylabel(TRIM(label),.FALSE.)
            ENDIF
            IF (ltherm_dos) THEN
               IF (with_t) &
                  CALL gnuplot_write_file_mul_data(filelastic,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          first_step,((.NOT.with_s).AND.(last_step.AND.&
                                       .NOT.ltherm_freq)),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.with_t),&
                          (last_step.AND..NOT.ltherm_freq),.FALSE.)
            ENDIF
            IF (ltherm_freq) THEN
               IF (with_t) &
                  CALL gnuplot_write_file_mul_data(filelastic_ph,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.ltherm_dos),&
                          ((.NOT.with_s).AND.last_step),.FALSE.)
               IF (with_s) &
                  CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,&
                          ec_present(iec, ec_type)+2,color(istep), &
                          (first_step.AND..NOT.ltherm_dos.AND..NOT.with_t),&
                                        last_step,.FALSE.)
            ENDIF
         ENDDO
      ENDIF
   ENDIF
ENDDO

istep=0
IF (MOD(iflag,2)==0) THEN
   DO itempp=1, ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filelastic="anhar_files/"//TRIM(flanhar)//'.el_cons_temp'
      CALL add_value(filelastic,temp(itemp))
      filelastic_s="anhar_files/"//TRIM(flanhar)//'.el_cons_s_temp'
      CALL add_value(filelastic_s,temp(itemp))
      filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_ph_temp'
      CALL add_value(filelastic_ph,temp(itemp))
      filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_s_ph_temp'
      CALL add_value(filelastic_s_ph,temp(itemp))
      IF (first_step) THEN
         CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
         WRITE(label,'("Bulk modulus B (kbar)")') 
         CALL gnuplot_ylabel(TRIM(label),.FALSE.)
      ENDIF
      IF (ltherm_dos) THEN
         IF (with_t) &
            CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep), &
                          first_step,(last_step.AND..NOT.ltherm_freq      &
                          .AND..NOT.with_s),.FALSE.)
         IF (with_s) &
            CALL gnuplot_write_file_mul_data(filelastic_s,1,2,color(istep), &
                    (first_step.AND..NOT.with_t),                           &
                    (last_step.AND..NOT.ltherm_freq),.FALSE.)
      ENDIF
      IF (ltherm_freq) THEN
         IF (with_t) &
            CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep), &
                (first_step.AND..NOT.ltherm_dos),(last_step.AND..NOT.with_s),&
                .FALSE.)
         IF (with_s) &
            CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,color(istep), &
                (first_step.AND..NOT.ltherm_dos.AND..NOT.with_s),&
                 last_step,.FALSE.)
      ENDIF
   ENDDO
ELSE
   DO itempp=1, ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filelastic="anhar_files/"//TRIM(flanhar)//'.el_comp_temp'
      CALL add_value(filelastic,temp(itemp))
      filelastic_s="anhar_files/"//TRIM(flanhar)//'.el_comp_s_temp'
      CALL add_value(filelastic_s,temp(itemp))
      filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_comp_ph_temp'
      CALL add_value(filelastic_ph,temp(itemp))
      filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.el_comp_s_ph_temp'
      CALL add_value(filelastic_s_ph,temp(itemp))
      IF (first_step) THEN
         CALL gnuplot_set_fact(1.D3, .FALSE.) 
         WRITE(label,'("Compressibility K (Mbar^{-1})")') 
         CALL gnuplot_ylabel(TRIM(label),.FALSE.)
      ENDIF
      IF (ltherm_dos) THEN
         IF (with_t) &
            CALL gnuplot_write_file_mul_data(filelastic,1,2,color(istep), &
                     first_step,(last_step.AND..NOT.ltherm_freq.AND.      &
                     .NOT.with_s),.FALSE.)
         IF (with_s) &
            CALL gnuplot_write_file_mul_data(filelastic_s,1,2,color(istep), &
                     (first_step.AND..NOT.with_t),                        &
                     (last_step.AND..NOT.ltherm_freq),.FALSE.)
      ENDIF
      IF (ltherm_freq) THEN
         IF (with_t) &
            CALL gnuplot_write_file_mul_data(filelastic_ph,1,2,color(istep), &
                (first_step.AND..NOT.ltherm_dos),(last_step.AND..NOT.with_s),&
                                                                .FALSE.)
         IF (with_s) &
            CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,2,color(istep), &
                (first_step.AND..NOT.ltherm_dos.AND..NOT.with_t),&
                 last_step,.FALSE.)
      ENDIF
   ENDDO
ENDIF

istep=0
IF (MOD(iflag,2)==0.AND.(laue==32.OR.laue==29)) THEN
   DO itempp=1, ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filelastic="anhar_files/"//TRIM(flanhar)//'.el_cons_temp'
      CALL add_value(filelastic,temp(itemp))
      filelastic_s="anhar_files/"//TRIM(flanhar)//'.el_cons_s_temp'
      CALL add_value(filelastic_s,temp(itemp))
      filelastic_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_ph_temp'
      CALL add_value(filelastic_ph,temp(itemp))
      filelastic_s_ph="anhar_files/"//TRIM(flanhar)//'.el_cons_s_ph_temp'
      CALL add_value(filelastic_s_ph,temp(itemp))
      IF (first_step) THEN
         CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
         WRITE(label,'("Shear anisotropy A")') 
         CALL gnuplot_ylabel(TRIM(label),.FALSE.)
      ENDIF
      IF (ltherm_dos) THEN
         IF (with_t) &
            CALL gnuplot_write_file_mul_data(filelastic,1,6,color(istep), &
                          first_step,(last_step.AND..NOT.ltherm_freq      &
                          .AND..NOT.with_s),.FALSE.)
         IF (with_s) &
            CALL gnuplot_write_file_mul_data(filelastic_s,1,6,color(istep), &
                    (first_step.AND..NOT.with_t),                           &
                    (last_step.AND..NOT.ltherm_freq),.FALSE.)
      ENDIF
      IF (ltherm_freq) THEN
         IF (with_t) &
            CALL gnuplot_write_file_mul_data(filelastic_ph,1,6,color(istep), &
                (first_step.AND..NOT.ltherm_dos),(last_step.AND..NOT.with_s),&
                .FALSE.)
         IF (with_s) &
            CALL gnuplot_write_file_mul_data(filelastic_s_ph,1,6,color(istep), &
                (first_step.AND..NOT.ltherm_dos.AND..NOT.with_s),&
                 last_step,.FALSE.)
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
END SUBROUTINE plot_elastic_ptt

!---------------------------------------------------------------------
SUBROUTINE plot_macro_elastic_p()
!---------------------------------------------------------------------
!
!  This routine plots the polycristalline averages of the
!  bulk modulus, Young modulus, shear modulus and poisson ration
!  calculated at T = 0 K as a function of pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE data_files,       ONLY : fl_el_cons
USE postscript_files, ONLY : flps_el_cons
USE control_elastic_constants,  ONLY : lelastic_p
USE control_pressure, ONLY : pmin, pmax
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filenameps, filelastic
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.lelastic_p) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_macro_el_p"
filenameps=TRIM(flps_el_cons)//".macro_el_p"//TRIM(flext)
filelastic="elastic_constants/"//TRIM(fl_el_cons)//".macro_el_p_aver"

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, &
                                            0.0_DP, 0.0_DP, 1.0_DP, flext ) 
CALL gnuplot_xlabel('p (kbar)', .FALSE.) 

CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Bulk modulus (B) (kbar)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Young modulus (E) (kbar)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Shear modulus (G) (kbar)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.0_DP, .FALSE.) 
CALL gnuplot_ylabel('Poisson ratio ({/Symbol n})',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_ylabel('B, E, G (kbar)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,&
                                                            .FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_blue',.FALSE.,&
                                                            .FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_green',.FALSE.,&
                                                            .TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_macro_elastic_p
!
!---------------------------------------------------------------------
SUBROUTINE plot_sound_speed_p()
!---------------------------------------------------------------------
!
!  This is a driver to plot the speed of sounds as a 
!  function of pressure at T = 0 K
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_set_fact, gnuplot_ylabel,     &
                             gnuplot_write_file_mul_data
USE data_files,       ONLY : fl_el_cons
USE postscript_files, ONLY : flps_el_cons
USE control_elastic_constants,  ONLY : lelastic_p
USE control_pressure, ONLY : pmin, pmax
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filenameps, filelastic
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.lelastic_p) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_sound_vel_p"
filenameps=TRIM(flps_el_cons)//".sound_vel_p"//TRIM(flext)
filelastic="elastic_constants/"//TRIM(fl_el_cons)//".sound_vel_p"

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, pmin*ry_kbar, pmax*ry_kbar, &
                                            0.0_DP, 0.0_DP, 1.0_DP, flext ) 
CALL gnuplot_xlabel('p (kbar)', .FALSE.) 

CALL gnuplot_set_fact(1.D-3, .FALSE.) 
CALL gnuplot_ylabel(' V_P (km/sec)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.D-3, .FALSE.) 
CALL gnuplot_ylabel('V_B (km/sec)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,.TRUE.,.FALSE.)

CALL gnuplot_set_fact(1.D-3, .FALSE.) 
CALL gnuplot_ylabel('V_S (km/sec)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE.,.TRUE.,.FALSE.)
!
!   All three speeds in the same plot.
!
CALL gnuplot_set_fact(1.D-3, .FALSE.) 
CALL gnuplot_ylabel('V_P, V_B, V_S (km/sec)',.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,&
                                                            .FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_blue',.FALSE.,&
                                                            .FALSE.,.FALSE.)
CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_green',.FALSE.,&
                                                            .TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_sound_speed_p
