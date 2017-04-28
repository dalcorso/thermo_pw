!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE initialize_plot_dispersion(kx, nks, xscale, ymin, ymax, eref,    &
     print_eref, nlines, start_point, last_point, label_disp_q,             &
     point_group_path, projective, lprojk, ylabel, gnu_filename, filenameps)
!
!   This routine initializes the gnuplot script and write in the script
!   the commands to plot the frame where the bands are plotted. It
!   writes the kpoints labels and when an enhace_plot is required it
!   writes the names of the point groups in the external points of 
!   each panel and in the center of each panel. 
!   The routine receives as input:
!   nks : the number of k points
!   kx(nks)  : the x coordinate of the plot for each k point
!   xscale,  : the scale of the x coordinates 
!   ymin, ymax : the minimum and maximum y of the plot
!   eref     : the position where to plot a horizontal line
!   print_eref : if true set a horizontal line at 0.0
!   nlines   : the number of panels in which the plot is divided
!   start_point(nlines) ! the stanting point of each panel
!   last_point(nlines)  ! the last point of each panel
!   label_disp_q(nqaux) ! for each k point label in the plot its position
!                         in the list of k point. The labels are passed
!                         using the control path module  
!   point_group_path(nlines,3)  ! for each line the name of the point
!                       ! group of the first point (in position 1) of the
!                       ! last point (in position 2) and of the line
!                       ! (in position 3). 
!   projective(nlines)  ! For each line it says if that line use projective
!                       ! representations
!   lprojk(nks)         ! For each point it says if that point use projective
!                       ! representation (1), switched projective 
!                       ! representation (2), or if it is a zone border but
!                       ! use standard representations (3).
!   ylabel              ! the string to write on y
!   gnu_filename        ! the name of gnuplot script
!   filenameps          ! the name of the postscript file.
!

USE kinds,           ONLY : DP
USE control_paths,   ONLY : nqaux, letter_path, q_in_band_form
USE control_2d_bands, ONLY : nkz
USE control_bands,   ONLY : enhance_plot
USE gnuplot,         ONLY : gnuplot_start, gnuplot_write_header,           &
                            gnuplot_ylabel, gnuplot_xlabel,                &
                            gnuplot_write_vertical_line,                   &
                            gnuplot_write_horizontal_segment,              &
                            gnuplot_put_label_yl, gnuplot_write_label,     &
                            gnuplot_write_label_yl, gnuplot_write_command, &
                            gnuplot_set_eref, gnuplot_print_objects,       &
                            gnuplot_unset_border, gnuplot_unset_xticks,    &
                            gnuplot_rectangle_yl  

USE io_global,     ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: gnu_filename, filenameps
CHARACTER(LEN=*), INTENT(IN) :: ylabel
INTEGER, INTENT(IN) :: nks, nlines
INTEGER, INTENT(IN) :: start_point(nlines), last_point(nlines), &
                       projective(nlines)
INTEGER, INTENT(IN) :: lprojk(nks)
INTEGER, INTENT(IN) :: label_disp_q(nqaux)

REAL(DP), INTENT(IN) :: kx(nks)
REAL(DP), INTENT(IN) :: ymin, ymax, eref, xscale 
CHARACTER(LEN=12), INTENT(IN) :: point_group_path(nlines,3)
LOGICAL, INTENT(IN) :: print_eref

INTEGER :: ilines, n, m, nks_, center_point, color_map(nks), colored_point(nks)
REAL(DP) :: x(2) 
CHARACTER(LEN=30) :: xlabel, color(5)
CHARACTER(LEN=256) :: command
!
!   header and commands specific for the bands
!
nks_=nks/nkz
CALL gnuplot_start(gnu_filename)
CALL gnuplot_write_header(filenameps, kx(1), kx(nks_), ymin, ymax, xscale ) 
CALL gnuplot_write_command('band_lw=2',.FALSE.)
!
!  shift above the maximum or the baseline of the point groups names
!
command="lshift=(ymax - ymin)/35."
CALL gnuplot_write_command(TRIM(command), .FALSE.)
!
!  shift below the minimum or the k points labels
!
command="shift=-(ymax - ymin)/40."
CALL gnuplot_write_command(TRIM(command), .FALSE.)
!
!  The x label is not set in a band plot. Use the following if you want
!  to compare with experiments and use A^-1
!
xlabel="k ({\305}^{-1})"
CALL gnuplot_xlabel(xlabel,.TRUE.) 
CALL gnuplot_ylabel(ylabel,.FALSE.) 
CALL gnuplot_set_eref(eref,.FALSE.)   
CALL gnuplot_unset_xticks(.FALSE.) 
CALL gnuplot_unset_border(.FALSE.) 
color(1)='color_black'
color(2)='color_orange'
color(3)='color_yellow'
color(4)='color_pink'
color(5)='color_gray'
!
!   For each point sets its color depending if the representation is
!   projective or standard
!
color_map=1
IF (enhance_plot) THEN
   DO n=1,nks
      IF (lprojk(n)==1) color_map(n)=2
   ENDDO
ENDIF
!
! First plot a vertical line in all the initial and final points of each panel
! Print also the horizontal lines of the panels and the Fermi level
! if requested in input. Note that the first and last line of the plot
! are always black, only the label will be colored if needed
!
colored_point=0
DO ilines=1, nlines
   IF (colored_point(start_point(ilines)) /= -1) THEN
      IF (ilines /= 1) THEN
         CALL gnuplot_write_vertical_line(kx(start_point(ilines)), 2, &
                   'front', color(color_map(start_point(ilines))), .FALSE.)
      ELSE
         CALL gnuplot_write_vertical_line(kx(start_point(ilines)), 2, &
                   'front', color(color_map(1)), .FALSE.)
      ENDIF
      colored_point(start_point(ilines))=-1 
   ENDIF
   IF (colored_point(last_point(ilines)) /= -1) THEN
      IF (ilines/=nlines) THEN
         CALL gnuplot_write_vertical_line(kx(last_point(ilines)), 2, &
                   'front', color(color_map(last_point(ilines))), .FALSE.)
      ELSE
         CALL gnuplot_write_vertical_line(kx(last_point(ilines)), 2, &
                   'front', color(1), .FALSE.)
      ENDIF
      colored_point(last_point(ilines))=-1 
   ENDIF
   x(1)=kx(start_point(ilines))
   x(2)=kx(last_point(ilines))
   CALL gnuplot_write_horizontal_segment(x, 'ymin', 2, 'front', &
                                                             color(1), .FALSE.)
   CALL gnuplot_write_horizontal_segment(x, 'ymax', 2, 'front', &
                                                             color(1), .FALSE.)
   IF (print_eref) &
      CALL gnuplot_write_horizontal_segment(x, '0.0', 2, 'front', color(1), &
                                                                   .FALSE.)
ENDDO
!
! If q_in_band_form some additional label might need a vertical line. Put
! it if it is not a starting or last point of a panel.
!
IF (q_in_band_form) THEN
   DO n=1, nqaux
      IF (n/=1.AND.n/=nqaux.AND.colored_point(label_disp_q(n)) /= -1 ) &
      CALL gnuplot_write_vertical_line(kx(label_disp_q(n)), 1, 'front', &
                        color(color_map(label_disp_q(n))), .FALSE.)
   ENDDO
END IF

IF (enhance_plot) THEN
!
!  with an enhanced plot, write the point group of the initial point
!  and of the line. If the last point of a panel does not coincide with 
!  the first point of the next panel write also the final point group
!
   DO ilines = 1, nlines
      CALL gnuplot_put_label_yl(kx(start_point(ilines)), ' ymax + lshift ', &
           ilines, point_group_path(ilines,1), &
           color(color_map(start_point(ilines))), .FALSE., 'center')

      IF ( ilines==nlines .OR. ABS(kx(last_point(ilines))- &
                                          kx(start_point(ilines+1)))>1.D-6 ) &
         CALL gnuplot_put_label_yl(kx(last_point(ilines)), ' ymax + lshift ', &
              nlines+ilines, point_group_path(ilines,2),  &
              color(color_map(last_point(ilines))), .FALSE., 'center')
      center_point=(start_point(ilines)+last_point(ilines))/2
      CALL gnuplot_put_label_yl(kx(center_point), ' ymin + lshift ',   &
          2*nlines+ilines, point_group_path(ilines,3), 'color_black',  &
                                                        .FALSE., 'center')
   ENDDO
!
!  Now put a rectancle in all the panels that correspond to zone
!  border lines. Reset the color map to the colors of the rectangles
!
   color_map=0
   DO ilines=1,nlines
      x(1)=kx(start_point(ilines))
      x(2)=kx(last_point(ilines))
      color_map(start_point(ilines):last_point(ilines))=projective(ilines)+2
      IF (color_map(start_point(ilines)) > 2) &
         CALL gnuplot_rectangle_yl(x, '1.0', &
                               color(color_map(start_point(ilines))))
   ENDDO
!
!   finally check if there is some part of the path that is at
!   zone border but has not been colored. In this case the panel must 
!   become gray
!
!   avoid that a rectangle starts at the last point of a line to not color 
!   the gaps in the plot
!
   DO ilines=1,nlines
      IF (color_map(last_point(ilines))==2) color_map(last_point(ilines))=5
   ENDDO
!
!  then check all the points that are at zone border. If they are not
!  colored, they might be the start of a gray rectangle.
!
   DO n=1,nks_-1
      IF (color_map(n)==2.AND.lprojk(n)/=0.AND.lprojk(n+1)==2) THEN
         x(1)=kx(n)
         x(2)=0.0_DP
         color_map(n)=5
         DO m=n+1,nks_
            IF (lprojk(m)/=3.AND.x(2)==0.0_DP)  THEN
               IF (lprojk(m)==0) THEN
                  x(2)=kx(m-1)
               ELSE
                  x(2)=kx(m)
               ENDIF
               CALL gnuplot_rectangle_yl(x, '1.0', 'color_gray')
            ELSEIF (x(2)==0.0_DP) THEN
               color_map(m)=5
            ENDIF
         ENDDO
      ENDIF
   ENDDO    
ENDIF
!
!  Finally the k point labels are written in black
!
DO n=1, nqaux
   IF (letter_path(n) /='') &
         CALL gnuplot_write_label_yl(kx(label_disp_q(n)), &
                          ' ymin + shift ', letter_path(n),.FALSE.)
END DO

RETURN
END SUBROUTINE initialize_plot_dispersion

SUBROUTINE plot_dispersion(nlines, nrap, has_points, dorap, with_lines, fileout)
!
!   This routine writes a part of the gnuplot script with the commands
!   to plot the band structure in a series of panels. It is supposed
!   that the information on the bands (x and y coordinates)
!   are in the file fileout.#line 
!   (when there is one representation or representations are not used)
!   or in the files fileout.#line.#rap when there is more than one 
!   representation. The routine receives the followin info:
!   nlines : the number of panels
!   nrap(nlines) : the number of representations per panel
!   dorap(12,nlines) : for each line and each representation this
!                      variable must be .true. to plot the bands of that
!                      representation
!   has_points(12,nlines) : for each line and each representation this
!                      variable is .true. if there bands of that representation
!   with_lines(nlines) : if .FALSE. the routine plots the bands with points 
!                instead of countinous lines for that line.
!
!
USE point_group,   ONLY : color_rap
USE gnuplot,       ONLY : gnuplot_write_file_data,  &
                          gnuplot_write_file_mul_point, &
                          gnuplot_write_command

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: fileout
INTEGER, INTENT(IN) :: nlines
INTEGER, INTENT(IN) :: nrap(nlines)
LOGICAL, INTENT(IN) :: has_points(12,nlines), dorap(12,nlines)
LOGICAL, INTENT(IN) :: with_lines(nlines)

INTEGER :: ilines, irap
INTEGER :: ir, first_rap, last_rap, first_line, last_line
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6), EXTERNAL :: int_to_char

IF (ANY(.NOT.with_lines(:))) &
   CALL gnuplot_write_command('point_size=0.3', .FALSE.)

first_rap=0
first_line=1
last_rap=0
last_line=1
DO ilines=1, nlines
   DO irap=1, nrap(ilines)
      IF (dorap(irap,ilines).AND.has_points(irap,ilines).AND.&
                                                        first_rap==0 ) THEN
         first_rap=irap 
         first_line=ilines
      ENDIF
      IF (has_points(irap,ilines).AND.dorap(irap,ilines)) THEN
         last_rap=irap 
         last_line=ilines
      ENDIF
   END DO
END DO

DO ilines = 1, nlines
   DO irap=1, nrap(ilines)
      IF (nrap(ilines)==1) THEN
         filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines))
      ELSE
         filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines)) &
                               //  "." // TRIM(int_to_char(irap))
      ENDIF
      IF (has_points(irap,ilines).AND.dorap(irap,ilines)) THEN
         IF (first_line==ilines .AND. irap==first_rap .AND. &
              last_line==ilines .AND. irap==last_rap ) THEN 
! 
!   Special case, only one line and one representation
!
            IF (with_lines(ilines)) THEN
               CALL gnuplot_write_file_data(filename,'band_lw',&
                             color_rap(irap),.TRUE.,.TRUE.,.FALSE.)
            ELSE
               CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                                      'color_red', .TRUE., .TRUE., .FALSE.)
            ENDIF
         ELSEIF (first_line==ilines .AND. irap==first_rap) THEN
            IF (with_lines(ilines)) THEN
               CALL gnuplot_write_file_data(filename,'band_lw',&
                             color_rap(irap),.TRUE.,.FALSE.,.FALSE.)
            ELSE
               CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                                      'color_red', .TRUE., .FALSE., .FALSE.)
            ENDIF
         ELSEIF (last_line==ilines .AND. irap==last_rap) THEN
            IF (with_lines(ilines)) THEN
               CALL gnuplot_write_file_data(filename,'band_lw',&
                           color_rap(irap),.FALSE.,.TRUE.,.FALSE.)
            ELSE
               CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                                      'color_red', .FALSE., .TRUE., .FALSE.)
            ENDIF
         ELSE
            IF (with_lines(ilines)) THEN
               CALL gnuplot_write_file_data(filename,'band_lw',&
                       color_rap(irap),.FALSE.,.FALSE.,.FALSE.)
            ELSE
               CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                               'color_red', .FALSE., .FALSE., .FALSE.)
            ENDIF
         ENDIF
      ENDIF  
   ENDDO
ENDDO

RETURN
END SUBROUTINE plot_dispersion

