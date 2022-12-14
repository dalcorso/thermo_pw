!
! Copyright (C) 2021 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
SUBROUTINE plot_geo_p()
!-------------------------------------------------------------------
!
!  This is a driver to plot the geometry as a function of pressure 
!  It plots: 
!     the celldm parameters as a function of p
!     the volume as a function of p
!     the ratio volume/volume0 as a function of p 
! volume0 is the volume at zero pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsmur
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_set_gfact,           &
                             gnuplot_write_command,       &
                             gnuplot_write_file_mul_data
USE data_files,       ONLY : flevdat
USE control_mur,      ONLY : omegap0
USE control_pressure, ONLY : pmin, pmax
USE initial_conf,     ONLY : ibrav_save
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1, filename2, gnu_filename, label
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_mur_celldm'
CALL add_pressure(gnu_filename)
filename=TRIM(flpsmur)//"_celldm"
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, pmin*ry_kbar, pmax*ry_kbar, &
                     0.0_DP, 0.0_DP, 1.0_DP, flext ) 

CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 

filename1="energy_files/"//TRIM(flevdat)//'_mur_celldm'
CALL add_pressure(filename1)
filename2="energy_files/"//TRIM(flevdat)//'_mur'
CALL add_pressure(filename2)
!
!  plot a
!
CALL gnuplot_ylabel('a (a.u.)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename1,1,2,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
!
!  plot b/a
!
IF (ABS(ibrav_save)>7) THEN
   CALL gnuplot_ylabel('b/a',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,3,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
ENDIF
!
!  plot c/a
!
IF (ibrav_save==4.OR.ABS(ibrav_save)>5) THEN
   CALL gnuplot_ylabel('c/a',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,4,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
ENDIF
!
!  plot cos(alpha)
!
IF (ibrav_save==5.OR.ibrav_save==12.OR.ibrav_save==13.OR.ibrav_save==14) THEN
   CALL gnuplot_ylabel('cos({\Symbol a})',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,5,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
ENDIF
!
!  plot cos(beta)
!
IF (ibrav_save==-12.OR.ibrav_save==-13.OR.ibrav_save==14) THEN
   CALL gnuplot_ylabel('cos({\Symbol b})',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,5,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
ENDIF
!
!  plot cos(gamma)
!
IF (ibrav_save==14) THEN
   CALL gnuplot_ylabel('cos({\Symbol c})',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename1,1,5,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
ENDIF
!
!   plot the volume
!
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename2,4,1,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
!
!  plot V/V_0 as a function of the
!
IF (omegap0>0.0_DP) THEN
   WRITE(label,'("set xrange [0:",f12.5,"]")') pmax*ry_kbar
   CALL gnuplot_write_command(TRIM(label),.FALSE.)
   CALL gnuplot_set_gfact(1.0_DP/omegap0,.FALSE.)
   CALL gnuplot_ylabel('V/V_0',.FALSE.) 
   CALL gnuplot_write_file_mul_data(filename2,4,1,'color_red',.TRUE.,.TRUE.,&
                                                                     .FALSE.)
ENDIF
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_geo_p

!-------------------------------------------------------------------
SUBROUTINE plot_geo_p_t()
!-------------------------------------------------------------------
!
!  This is a driver to plot the geometry as a function of pressure 
!  It plots: 
!     the celldm parameters as a function of p
!     the volume as a function of p
!     the ratio volume/volume0 as a function of p 
! volume0 is the volume at zero pressure
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsmur
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,  &
                             gnuplot_write_header,        &
                             gnuplot_ylabel,              &
                             gnuplot_xlabel,              &
                             gnuplot_set_gfact,           &
                             gnuplot_write_command,       &
                             gnuplot_write_file_mul_data
USE data_files,       ONLY : flevdat, flanhar
USE control_mur,      ONLY : omegap0
USE control_pressure, ONLY : pmin, pmax
USE temperature,      ONLY : temp, ntemp, ntemp_plot, itemp_plot
USE color_mod,        ONLY : color
USE initial_conf,     ONLY : ibrav_save
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, gnu_filename, label
CHARACTER(LEN=6) :: int_to_char
INTEGER :: istep, itemp, itempp, ierr, system
LOGICAL :: first_step, last_step

IF ( my_image_id /= root_image ) RETURN

IF (ntemp_plot==0) RETURN

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//'_mur_celldm_t'
CALL add_pressure(gnu_filename)
filename=TRIM(flpsmur)//"_celldm_t"
CALL add_pressure(filename)
filename=TRIM(filename)//TRIM(flext)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filename, pmin*ry_kbar, pmax*ry_kbar, &
                     0.0_DP, 0.0_DP, 1.0_DP, flext ) 

CALL gnuplot_xlabel('pressure (kbar)',.FALSE.) 

!
!  plot a
!
CALL gnuplot_ylabel('a (a.u.)',.FALSE.) 

istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flevdat)//'_mur_celldm'//&
                                               TRIM(int_to_char(itemp))
   CALL add_pressure(filename)
   CALL gnuplot_write_file_mul_data(filename,1,2,color(istep),first_step, &
                                 last_step, .FALSE.)
ENDDO
!
!  plot b/a
!
IF (ABS(ibrav_save)>7) THEN
   CALL gnuplot_ylabel('b/a',.FALSE.) 
   istep=0
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flevdat)//'_mur_celldm'//&
                                               TRIM(int_to_char(itemp))
      CALL add_pressure(filename)
      CALL gnuplot_write_file_mul_data(filename,1,3,color(istep),first_step,&
                                    last_step,.FALSE.)
   ENDDO
ENDIF
!
!  plot c/a
!
IF (ibrav_save==4.OR.ABS(ibrav_save)>5) THEN
   CALL gnuplot_ylabel('c/a',.FALSE.) 
   istep=0
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flevdat)//'_mur_celldm'//&
                                               TRIM(int_to_char(itemp))
      CALL add_pressure(filename)
      CALL gnuplot_write_file_mul_data(filename,1,4,color(istep), &
                               first_step,last_step,.FALSE.)
   ENDDO
ENDIF
!
!  plot cos(alpha)
!
IF (ibrav_save==5.OR.ibrav_save==12.OR.ibrav_save==13.OR.ibrav_save==14) THEN
   CALL gnuplot_ylabel('cos({\Symbol a})',.FALSE.) 
   istep=0
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flevdat)//'_mur_celldm'//&
                                               TRIM(int_to_char(itemp))
      CALL add_pressure(filename)
      CALL gnuplot_write_file_mul_data(filename,1,5,color(istep),&
                              first_step,last_step,.FALSE.)
   ENDDO
ENDIF
!
!  plot cos(beta)
!
IF (ibrav_save==-12.OR.ibrav_save==-13.OR.ibrav_save==14) THEN
   CALL gnuplot_ylabel('cos({\Symbol b})',.FALSE.) 
   istep=0
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flevdat)//'_mur_celldm'//&
                                               TRIM(int_to_char(itemp))
      CALL add_pressure(filename)
      CALL gnuplot_write_file_mul_data(filename,1,5,color(istep), &
                          first_step,last_step,.FALSE.)
   ENDDO
ENDIF
!
!  plot cos(gamma)
!
IF (ibrav_save==14) THEN
   CALL gnuplot_ylabel('cos({\Symbol c})',.FALSE.) 
   istep=0
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)

      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flevdat)//'_mur_celldm'//&
                                               TRIM(int_to_char(itemp))
      CALL add_pressure(filename)
      CALL gnuplot_write_file_mul_data(filename,1,5,color(istep),&
                                       first_step,last_step,.FALSE.)
   ENDDO
ENDIF
!
!   plot the volume
!
CALL gnuplot_ylabel('Volume ((a.u.)^3)',.FALSE.) 
istep=0
DO itempp=1,ntemp_plot
   first_step=(itempp==1)
   last_step=(itempp==ntemp_plot)
   itemp=itemp_plot(itempp)
   istep=MOD(istep,8)+1
   filename="anhar_files/"//TRIM(flanhar)//'.mur_temp'
   CALL add_value(filename,temp(itemp))
   CALL gnuplot_write_file_mul_data(filename,4,1,color(istep),&
                                          first_step,last_step,.FALSE.)
ENDDO
!
!  plot V/V_0 as a function of the
!
IF (omegap0>0.0_DP) THEN
   WRITE(label,'("set xrange [0:",f12.5,"]")') pmax*ry_kbar
   CALL gnuplot_write_command(TRIM(label),.FALSE.)
   CALL gnuplot_set_gfact(1.0_DP/omegap0,.FALSE.)
   CALL gnuplot_ylabel('V/V_0',.FALSE.) 
   istep=0
   DO itempp=1,ntemp_plot
      first_step=(itempp==1)
      last_step=(itempp==ntemp_plot)
      itemp=itemp_plot(itempp)
      istep=MOD(istep,8)+1
      filename="anhar_files/"//TRIM(flevdat)//'_mur'//&
                                               TRIM(int_to_char(itemp))
      CALL add_pressure(filename)
      CALL gnuplot_write_file_mul_data(filename,4,1,color(istep),&
                        first_step, last_step, .FALSE.)
   ENDDO
ENDIF
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_geo_p_t
