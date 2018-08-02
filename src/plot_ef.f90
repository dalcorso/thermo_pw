!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE plot_ef(filband, flgnuplot, filenameps)
  !-----------------------------------------------------------------------
  !
  !  this routine writes a gnuplot script to make a two dimensional 
  !  countour plot of the fermi energy. 
  !
USE constants, ONLY : rytoev
USE control_gnuplot, ONLY : lgnuplot, gnuplot_command, flext
USE klist,     ONLY : ltetra, lgauss
USE ener, ONLY : ef
uSE gnuplot,   ONLY : gnuplot_start, gnuplot_end,             &
                      gnuplot_start_2dplot,                   &
                      gnuplot_set_contour, gnuplot_do_2dplot, &
                      gnuplot_xlabel,     &
                      gnuplot_ylabel, gnuplot_close_2dplot_prep
USE efermi_plot, ONLY : n1, n2, kxmin, kxmax, kymin, kymax, has_ef
USE wvfct,       ONLY : nbnd
USE mp_images, ONLY : my_image_id, root_image
USE io_global, ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filband, flgnuplot, filenameps
INTEGER :: ncountours, band_with_ef, ibnd, ierr

INTEGER :: system
CHARACTER(LEN=256) :: tablefile, filename, gnu_filename, flnameps
CHARACTER(LEN=6)   :: int_to_char
CHARACTER(LEN=12)  :: color(8)
CHARACTER(LEN=50)  :: xlabel, ylabel

IF ( my_image_id /= root_image ) RETURN
!
! The plot is done only in the metallic case
!
IF (.NOT.(lgauss.OR.ltetra)) RETURN

gnu_filename=TRIM(flgnuplot)//'_ef'
color(1)='color_red'
color(2)='color_green'
color(3)='color_blue'
color(4)='color_yellow'
color(5)='color_pink'
color(6)='color_cyan'
color(7)='color_orange'
color(8)='color_black'
tablefile='gnuplot_files/table'
CALL gnuplot_start(gnu_filename)

CALL gnuplot_start_2dplot(nbnd, n1, n2)

ncountours=0
DO ibnd=1, nbnd
   IF (has_ef(ibnd)) THEN
      ncountours=ncountours+1
      filename=TRIM(filband) // '.' // TRIM(int_to_char(ibnd))
      CALL gnuplot_set_contour(filename, ef*rytoev, &
                                     color(MOD(ncountours-1,8)+1), tablefile)
   ENDIF
ENDDO

CALL gnuplot_close_2dplot_prep()

xlabel=' '
ylabel=' '

flnameps=TRIM(filenameps)//'_ef'//TRIM(flext)
CALL gnuplot_do_2dplot(flnameps, kxmin, kxmax, kymin, kymax, xlabel, &
                                                  ylabel, tablefile, flext)
CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
     ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!  IF (lgnuplot.AND.ionode) &
!     CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
!                                       //TRIM(gnu_filename),
!                                       WAIT=.FALSE.)
RETURN
END SUBROUTINE plot_ef
