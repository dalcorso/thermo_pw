!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_epsilon_omega_opt()
!
!  This is a driver to plot the quantities inside the epsilon or polariz
!  files. The number of plots depends on the geometry and on the
!  fact that the system is an isolated molecule or not.
!
!  For a solid it makes two plots one for the real and one for the 
!  imaginary part of epsilon. Depending on the Bravais lattice it
!  plots one or more curves (up to three) on the same plot: 
!  ibrav = 1,2,3  it plots only epsilon_xx 
!  ibrav = 4,5,6,7  it plots only epsilon_xx (red) and epsilon_zz (green)
!  ibrav = 8,9,10,11 it plots epsilon_xx (red), epsilon_yy (green), 
!                   epsilon_zz (blue).
!  ibrav = 12, 13, 14 at each frequency it diagonalizes the dielectric
!                    tensor and plots the three eigenvalues as a function
!                    of the frequency. Note that the real and the imaginary
!                    part can have different principal directions and
!                    these might change with the frequency.
!  
!  The code assumes that the system is a molecule if there is only one
!  k point (or two in the lsda case). In this case it plots only
!  the trace of the dielectric constant tensor (real and imaginary).
!  Moreover it plots also the product of the trace of the imaginary part
!  and the frequency, proportional to the average photoabsorption cross
!  section.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : rytoev
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsepsilon
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,   &
                            gnuplot_write_header,          &
                            gnuplot_ylabel,                &
                            gnuplot_xlabel,                &
                            gnuplot_write_command,         &
                            gnuplot_write_file_mul_data
USE lsda_mod,         ONLY : lsda
USE klist,            ONLY : nkstot
USE data_files,       ONLY : flepsilon
USE cell_base,        ONLY : ibrav
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, string
INTEGER :: im
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_epsilon'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsepsilon)//TRIM(flext)
CALL gnuplot_write_header(filename, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                                                    rytoev, flext ) 

CALL gnuplot_xlabel('{/Symbol w}  (eV)',.FALSE.) 

DO im=1,2 
   IF (nkstot==1.OR.(nkstot==2.AND.lsda)) THEN
!
!  molecular case, plot the trace of the dielectric tensor multiplied
!  by the frequency
!      
      IF (im==1) THEN
         CALL gnuplot_ylabel('tr {/Symbol e}_1 ({/Symbol w})',.FALSE.) 
         filename='dynamical_matrices/'//TRIM(flepsilon)//'_re'
         string="plot """//TRIM(filename)//""" u ($1*13.605698066):(($3+$4+$5)&
             &/3.0) w l lw 3 lc rgb color_red"
      ELSE
         CALL gnuplot_ylabel('tr {/Symbol e}_2 ({/Symbol w})',.FALSE.) 
         filename='dynamical_matrices/'//TRIM(flepsilon)//'_im'
         string="plot """//TRIM(filename)//""" u ($1*13.605698066):(($3+$4+$5)&
             &/3.0) w l lw 3 lc rgb color_red"
      ENDIF

      CALL gnuplot_write_command(string,.FALSE.)
      IF (im==2) THEN
         CALL gnuplot_ylabel('Absorption   {/Symbol s}_2 ({/Symbol w}) (a.u.)^2',&
                                                        .FALSE.) 
         filename='dynamical_matrices/polariz_im'
!
!   the factor corresponds to 2*pi / c where c is the speed of light in a.u.
!   the formula assumes that the frequency is given in Ry and the 
!   polarizability in (bohr)^3. The absorption is in (bohr)^2.
!
         string="plot """//TRIM(filename)//""" u ($1*13.605698066):(($3+$4+$5)&
             &*$1*0.045850618/3.0) w l lw 3 lc rgb color_red"
         CALL gnuplot_write_command(string,.FALSE.)
      ENDIF
   ELSE
!
!  solid case, depending on the Bravais lattice plot the dielectric constant
!
      IF (im==1) THEN
         CALL gnuplot_ylabel('{/Symbol e}_1 ({/Symbol w})',.FALSE.) 
         filename='dynamical_matrices/'//TRIM(flepsilon)//'_re'
      ELSE
         CALL gnuplot_ylabel('{/Symbol e}_2 ({/Symbol w})',.FALSE.) 
         filename='dynamical_matrices/'//TRIM(flepsilon)//'_im'
      ENDIF

      IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
         CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)
      ELSE IF (ibrav==4.OR.ibrav==5.OR.ibrav==6.OR.ibrav==7) THEN
         CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE.,&
                                                .FALSE.,.FALSE.)
         CALL gnuplot_write_file_mul_data(filename,1,5,'color_green',.FALSE.,&
                                                .TRUE.,.FALSE.)
      ELSE IF (ibrav==8.OR.ibrav==9.OR.ibrav==10.OR.ibrav==11) THEN
         CALL gnuplot_write_file_mul_data(filename,1,3,'color_red',.TRUE.,&
                                                .FALSE.,.FALSE.)
         CALL gnuplot_write_file_mul_data(filename,1,4,'color_green',.FALSE.,&
                                                .FALSE.,.FALSE.)
         CALL gnuplot_write_file_mul_data(filename,1,5,'color_blue',.FALSE.,&
                                                .TRUE.,.FALSE.)
      ELSE IF (ABS(ibrav)==12.OR.ABS(ibrav)==13.OR.ibrav==14) THEN
!
!  To be programmed
!
      END IF
   END IF
END DO

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_epsilon_omega_opt

SUBROUTINE plot_epsilon_omega_q()
!
!  This is a driver to plot the quantities inside the epsilon files. It makes 
!  four plots. The first two contain the real and the imaginary part of the 
!  inverse of the dielectric constant (q w), the other two the real and 
!  imaginary part of the dielectric constant of (q w).
!
USE kinds,            ONLY : DP
USE constants,        ONLY : rytoev
USE control_gnuplot,  ONLY : flgnuplot, lgnuplot, gnuplot_command, flext
USE postscript_files, ONLY : flpsepsilon
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,   &
                            gnuplot_write_header,          &
                            gnuplot_ylabel,                &
                            gnuplot_xlabel,                &
                            gnuplot_write_command,         &
                            gnuplot_write_file_mul_data_minus, &
                            gnuplot_write_file_mul_data
USE data_files,       ONLY : flepsilon
USE mp_images,        ONLY : my_image_id, root_image
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: gnu_filename, filename, string
INTEGER :: im
INTEGER :: ierr, system

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_files/'//TRIM(flgnuplot)//'_epsilon'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsepsilon)//TRIM(flext)
CALL gnuplot_write_header(filename, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                                                               rytoev, flext ) 

CALL gnuplot_xlabel('{/Symbol w}  (eV)',.FALSE.) 

CALL gnuplot_ylabel('Re 1 / {/Symbol e} (q, {/Symbol w})',.FALSE.) 
filename='dynamical_matrices/'//TRIM(flepsilon)
CALL gnuplot_write_file_mul_data(filename,2,4,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_ylabel('- Im 1 / {/Symbol e} (q, {/Symbol w})',.FALSE.) 
CALL gnuplot_write_file_mul_data_minus(filename,2,5,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_ylabel('{/Symbol e}_1 (q, {/Symbol w})',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,2,6,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_ylabel('{/Symbol e}_2 (q, {/Symbol w})',.FALSE.) 
CALL gnuplot_write_file_mul_data(filename,2,7,'color_red',.TRUE.,&
                                                .TRUE.,.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!IF (lgnuplot.AND.ionode) &
!   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename), WAIT=.FALSE.)

RETURN
END SUBROUTINE plot_epsilon_omega_q
