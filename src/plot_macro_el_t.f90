!
! Copyright (C) 2020 Cristiano Malica 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE plot_macro_el_t
!-------------------------------------------------------------------------
!
!  This is a driver to plot the macro-elasticity variables (MEVs) as a function of
!  temperature. It plots both the MEVs computed from isothermal elastic constants (ECs)
!  and from adiabatic ECs in two separated files. It compares the results obtained
!  from Voigt average (color red), Reuss average (color blue) and
!  Voigt-Reuss-Hill average (color green).
!  It also plots the three sound velocities (V_P, V_B, V_G) by comparing the isothermal 
!  (red) with the adiabatic (blue) result in the single plot. 
!  It is used by both quasi-static or quasi-harmonic temperature dependent ECs. 
!  If both lelastic and lelasticf are .TRUE. it uses only the ECs obtained from the 
!  free energy calculated from the phonon dos.
!
USE control_gnuplot,  ONLY : flgnuplot, flext
USE data_files,       ONLY : flanhar
USE postscript_files, ONLY : flpsanhar
USE anharmonic,       ONLY : lelastic
USE ph_freq_anharmonic,  ONLY : lelasticf
USE control_grun,     ONLY : lb0_t
USE mp_images,        ONLY : root_image, my_image_id

IMPLICIT NONE
CHARACTER(LEN=256) :: gnu_filename, filenameps, filelastic

IF ( my_image_id /= root_image ) RETURN
IF (.NOT.(lelastic.OR.lelasticf).OR..NOT.lb0_t) RETURN

IF (lelasticf) filelastic="anhar_files/"//TRIM(flanhar)//".macro_el_ph"
IF (lelastic)  filelastic="anhar_files/"//TRIM(flanhar)//".macro_el"
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_macro_el"
filenameps=TRIM(flpsanhar)//".macro_el"//TRIM(flext)

CALL plot_one_macro_el(filelastic, gnu_filename, filenameps)

IF (lelasticf) filelastic="anhar_files/"//TRIM(flanhar)//".macro_el_s_ph"
IF (lelastic)  filelastic="anhar_files/"//TRIM(flanhar)//".macro_el_s"
gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_macro_el_s"
filenameps=TRIM(flpsanhar)//".macro_el_s"//TRIM(flext)

CALL plot_one_macro_el(filelastic, gnu_filename, filenameps)


CALL plot_sound()

RETURN
END SUBROUTINE plot_macro_el_t

!-------------------------------------------------------------------------
SUBROUTINE plot_one_macro_el(filelastic, gnu_filename, filenameps)
!-------------------------------------------------------------------------

USE kinds,            ONLY : DP
USE temperature,      ONLY : tmin, tmax
USE control_gnuplot,  ONLY : gnuplot_command, lgnuplot, flext
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_ylabel, gnuplot_write_file_mul_data
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256), INTENT(INOUT) :: filelastic, gnu_filename, filenameps
CHARACTER(LEN=256) :: filelastic_aver
INTEGER :: ierr, system

filelastic_aver=TRIM(filelastic)//"_aver"

filelastic=TRIM(filelastic)
filelastic_aver=TRIM(filelastic_aver)
gnu_filename=TRIM(gnu_filename)
filenameps=TRIM(filenameps)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext ) 
ENDIF

CALL gnuplot_xlabel('T (K)', .FALSE.) 

CALL gnuplot_ylabel('Bulk modulus (kbar)',.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,6,'color_blue',.FALSE., &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic_aver,1,2,'color_green',.FALSE., &
                                                              .TRUE.,.FALSE.)

CALL gnuplot_ylabel('Young modulus (kbar)',.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,7,'color_blue',.FALSE., &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic_aver,1,3,'color_green',.FALSE., &
                                                              .TRUE.,.FALSE.)

CALL gnuplot_ylabel('Shear modulus (kbar)',.FALSE.)
        
CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE., &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,8,'color_blue',.FALSE., &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic_aver,1,4,'color_green',.FALSE., &
                                                              .TRUE.,.FALSE.)

CALL gnuplot_ylabel('Poisson ratio',.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,5,'color_red',.TRUE., &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,9,'color_blue',.FALSE., &
                                                              .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic_aver,1,5,'color_green',.FALSE., &
                                                              .TRUE.,.FALSE.)
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_one_macro_el

!-------------------------------------------------------------------------
SUBROUTINE plot_sound()
!-------------------------------------------------------------------------

USE kinds,               ONLY : DP
USE temperature,         ONLY : tmin, tmax
USE data_files,          ONLY : flanhar
USE postscript_files,    ONLY : flpsanhar
USE anharmonic,          ONLY : lelastic
USE ph_freq_anharmonic,  ONLY : lelasticf
USE control_gnuplot,     ONLY : gnuplot_command, flgnuplot, lgnuplot, flext
USE gnuplot,             ONLY : gnuplot_start, gnuplot_end,           &
                                gnuplot_write_header, gnuplot_xlabel, &
                                gnuplot_ylabel, gnuplot_write_file_mul_data
USE io_global,        ONLY : ionode
USE mp_images,        ONLY : root_image, my_image_id

IMPLICIT NONE

CHARACTER(LEN=256) :: filelastic, filelastic_s, gnu_filename, filenameps
INTEGER :: ierr, system

filelastic=TRIM(filelastic)
gnu_filename=TRIM(gnu_filename)
filenameps=TRIM(filenameps)

IF (lelasticf) THEN
   filelastic="anhar_files/"//TRIM(flanhar)//".sound_vel_ph"
   filelastic_s="anhar_files/"//TRIM(flanhar)//".sound_vel_s_ph"
END IF

IF (lelastic) THEN
   filelastic="anhar_files/"//TRIM(flanhar)//".sound_vel"
   filelastic_s="anhar_files/"//TRIM(flanhar)//".sound_vel_s"
END IF

gnu_filename="gnuplot_files/"//TRIM(flgnuplot)//"_sound_vel"
filenameps=TRIM(flpsanhar)//".sound_vel"//TRIM(flext)

CALL gnuplot_start(gnu_filename)

IF (tmin ==1._DP) THEN
   CALL gnuplot_write_header(filenameps, 0.0_DP, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ELSE
   CALL gnuplot_write_header(filenameps, tmin, tmax, 0.0_DP, 0.0_DP, &
                                                       1.0_DP, flext )
ENDIF

CALL gnuplot_xlabel('T (K)', .FALSE.)

CALL gnuplot_ylabel('V_{P} (m/s)',.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,2,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic_s,1,2,'color_blue',.FALSE., &
                                                              .TRUE.,.FALSE.)
CALL gnuplot_ylabel('V_{B} (m/s)',.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,3,'color_red',.TRUE.,  &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic_s,1,3,'color_blue',.FALSE., &
                                                              .TRUE.,.FALSE.)
CALL gnuplot_ylabel('V_{G} (m/s)',.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic,1,4,'color_red',.TRUE., &
                                                             .FALSE.,.FALSE.)

CALL gnuplot_write_file_mul_data(filelastic_s,1,4,'color_blue',.FALSE., &
                                                              .TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
           ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN

END SUBROUTINE plot_sound

