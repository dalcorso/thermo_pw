!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_dispersion(kx, e, xk, nks, nbnd, emin, emax, eref, nlines, &
     nrap, rap, gcodek, e_rap, nbnd_rapk, start_rapk, has_points, &
     start_point, last_point, nrap_plot, rap_plot, label_disp_q, icode, exist_rap, &
     igeom, fileout)

USE kinds,           ONLY : DP
USE constants,       ONLY : bohr_radius_si
USE ions_base,       ONLY : nat
USE cell_base,       ONLY : tpiba
USE klist,           ONLY : degauss
USE control_paths,   ONLY : nqaux, letter_path, q_in_band_form
USE postscript_files, ONLY : flpsband, flpsdisp, flpsgrun
USE control_gnuplot, ONLY : flgnuplot, gnuplot_command, lgnuplot
USE control_2d_bands, ONLY : nkz, lprojpbs, identify_sur, sym_divide, &
                          force_bands
USE point_group,   ONLY : color_rap
USE gnuplot,       ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                          gnuplot_write_file_data, gnuplot_ylabel, &
                          gnuplot_write_vertical_line, gnuplot_write_label, &
                          gnuplot_write_horizontal_line, &
                          gnuplot_write_label_yl, gnuplot_write_command, &
                          gnuplot_set_eref, gnuplot_unset_xticks,   &
                          gnuplot_write_file_mul_point,             &
                          gnuplot_xlabel,  &
                          gnuplot_print_objects, gnuplot_write_command
USE io_global,     ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: fileout
INTEGER, INTENT(IN) :: nks, nbnd, nlines, icode, igeom
INTEGER, INTENT(IN) :: nrap(nlines), gcodek(nks), rap(nbnd,nks)
INTEGER, INTENT(IN) :: start_point(nks), last_point(nks)
INTEGER, INTENT(IN) :: nbnd_rapk(12,nks), start_rapk(12,nks)
INTEGER, INTENT(IN) :: nrap_plot(nks), rap_plot(12,nks)
INTEGER, INTENT(IN) :: label_disp_q(nqaux)
REAL(DP), INTENT(IN) :: kx(nks), e(nbnd,nks), xk(3,nks)
REAL(DP), INTENT(IN) :: e_rap(nbnd,nks)
REAL(DP), INTENT(IN) :: emin, emax, eref 
LOGICAL, INTENT(IN) :: has_points(nlines,12)
LOGICAL, INTENT(IN) :: exist_rap

INTEGER :: ilines, irap, n, nks_, ierr
INTEGER :: nrapp, ir, first_rap, last_rap, first_line, last_line
INTEGER :: system
REAL(DP) :: shift, xscale
LOGICAL, ALLOCATABLE :: dorap(:,:)
LOGICAL :: with_lines
CHARACTER(LEN=256) :: gnu_filename, filename, command
CHARACTER(LEN=30) :: xlabel
CHARACTER(LEN=6), EXTERNAL :: int_to_char

IF (icode==1) THEN
   gnu_filename=TRIM(flgnuplot)//'_band'
ELSEIF (icode==2) THEN
   gnu_filename=TRIM(flgnuplot)//'_disp'
ELSEIF (icode==3) THEN
   gnu_filename=TRIM(flgnuplot)//'_grun'
ELSEIF (icode==4) THEN
   gnu_filename=TRIM(flgnuplot)//'_grun_freq'
ENDIF
gnu_filename="gnuplot_files/"//TRIM(gnu_filename)
!
!  The Gruneisen parameters of dispersions that do not have representations
!  must be plotted with points
!

CALL gnuplot_start(gnu_filename)

with_lines=.TRUE.
IF (icode==1) THEN
   filename=TRIM(flpsband)
ELSEIF (icode==2) THEN
   filename=TRIM(flpsdisp)
ELSEIF (icode==3) THEN
   filename=TRIM(flpsgrun)
   with_lines=.FALSE.
ELSEIF (icode==4) THEN
   filename=TRIM(flpsgrun)//'_freq'
ENDIF
IF (igeom > 1) filename=filename//TRIM(int_to_char(igeom))

nks_=nks / nkz 
xscale=tpiba / bohr_radius_si / 1.D10
CALL gnuplot_write_header(filename, kx(1), kx(nks_), emin, emax, xscale ) 
xlabel="k ({\305}^{-1})"
CALL gnuplot_xlabel(xlabel, .TRUE.) 
CALL gnuplot_unset_xticks(.FALSE.) 
IF (icode==1) THEN
   CALL gnuplot_ylabel('Energy (eV)',.FALSE.) 
   IF (degauss > 0.0_DP) CALL gnuplot_write_horizontal_line(0.0_DP, 2, &
                                         'front', 'color_black', .FALSE.)
ELSEIF (icode==2.OR.icode==4) THEN
   CALL gnuplot_ylabel('Frequency (cm^{-1})',.FALSE.) 
ELSEIF (icode==3) THEN
   CALL gnuplot_ylabel('Mode-Gr\374neisen parameters  {/Symbol g}_{/Symbol n}&
                 &({/Helvetica-Bold q})',.FALSE.) 
   CALL gnuplot_write_command('point_size=0.3', .FALSE.)
ENDIF
DO ilines = 2, nlines
   CALL gnuplot_write_vertical_line(kx(start_point(ilines)), 2, 'front', &
                      'color_black', .FALSE.)
END DO
CALL gnuplot_set_eref(eref,.FALSE.)   
!
!  The letters are below the minimum of this quantity
!
command="shift=-(ymax - ymin)/40."
CALL gnuplot_write_command(TRIM(command), .FALSE.)

IF (q_in_band_form) THEN
   DO n=1, nqaux
      IF (n /= 1 .AND. n /= nqaux ) &
         CALL gnuplot_write_vertical_line(kx(label_disp_q(n)), 1, 'front', &
                                       'color_black', .FALSE.)
      CALL gnuplot_write_label_yl(kx(label_disp_q(n)), ' ymin + shift ', &
             letter_path(n),.FALSE.)
   ENDDO
ELSE
   DO n=1, nqaux
      IF (letter_path(n) /='') &
         CALL gnuplot_write_label_yl(kx(label_disp_q(n)), &
                          ' ymin + shift ', letter_path(n),.FALSE.)
   END DO
END IF

CALL gnuplot_write_command('band_lw=2',.FALSE.)

IF (lprojpbs) CALL proj_band_structure(kx, e, nks, nbnd, emin, emax, eref, &
                  e_rap, nrap, nbnd_rapk, start_rapk, nlines, start_point, &
                  last_point, nrap_plot, rap_plot )

IF (identify_sur) CALL plot_surface_states(nbnd, nks, nlines, kx, e_rap, &
                            emin, emax, eref, nrap, nbnd_rapk, start_rapk, &
                            start_point, last_point, nrap_plot, rap_plot )

IF ( nkz == 1 .OR. ( nkz > 1 .AND. .NOT. lprojpbs) .OR. force_bands ) THEN
!
! In this case the projected band structure have been read from input 
! or it is not available because this is a standard calculation, we can
! plot the single bands
! 
!
   IF (.NOT.exist_rap) THEN
      filename=TRIM(fileout)
      IF (with_lines) THEN
         CALL gnuplot_write_file_data(filename,'band_lw','color_red',&
                                           .TRUE.,.TRUE.,.FALSE.)
      ELSE
         CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                                         'color_red', .TRUE., .TRUE., .FALSE.)
      ENDIF
   ELSE
      ALLOCATE(dorap(nlines,12))
      first_rap=0
      first_line=1
      last_rap=0
      last_line=1
      DO ilines=1, nlines
         DO irap=1, nrap(ilines)
            dorap(ilines,irap)=.TRUE.
            IF (sym_divide) THEN
               dorap(ilines,irap)=.FALSE.
               nrapp= nrap_plot(start_point(ilines)+1)
               DO ir=1,nrapp
                  dorap(ilines,irap) = dorap(ilines,irap) &
                       .OR.(rap_plot(ir,start_point(ilines)+1)==irap) 
               END DO
               IF (nrapp==0) dorap(ilines,irap)=.TRUE.
            END IF
            IF (dorap(ilines,irap).AND.has_points(ilines,irap) &
                            .AND.first_rap==0 ) THEN
               first_rap=irap 
               first_line=ilines
            ENDIF
            IF (has_points(ilines,irap).AND.dorap(ilines,irap)) THEN
               last_rap=irap 
               last_line=ilines
            ENDIF
         END DO
      END DO

      DO ilines = 1, nlines
         IF (nrap(ilines)==1) THEN
            filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines))
            IF (has_points(ilines,1)) THEN
               IF (first_line==ilines .AND. first_rap==1 ) THEN
                  IF (with_lines) THEN
                     CALL gnuplot_write_file_data(filename,'band_lw', &
                                       'color_red', .TRUE.,.FALSE.,.FALSE.)
                  ELSE
                     CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                                         'color_red', .TRUE., .FALSE., .FALSE.)
                  ENDIF
               ELSEIF (last_line==ilines.AND.last_rap==1) THEN
                 IF (with_lines) THEN
                     CALL gnuplot_write_file_data(filename,'band_lw',&
                                         'color_red',.FALSE.,.TRUE.,.FALSE.)
                  ELSE
                     CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                                         'color_red', .FALSE., .TRUE., .FALSE.)
                  ENDIF
               ELSE
                 IF (with_lines) THEN
                     CALL gnuplot_write_file_data(filename,'band_lw',&
                                  'color_red', .FALSE.,.FALSE.,.FALSE.)
                  ELSE
                     CALL gnuplot_write_file_mul_point(filename, 1, 2, &
                                         'color_red', .FALSE., .FALSE., .FALSE.)
                  ENDIF
               END IF
            END IF
         ELSE
            DO irap=1, nrap(ilines)
               filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines)) &
                                    //  "." // TRIM(int_to_char(irap))
               IF (has_points(ilines,irap).AND.dorap(ilines,irap)) THEN
                  IF (first_line==ilines .AND. irap==first_rap) THEN
                     CALL gnuplot_write_file_data(filename,'band_lw',&
                                    color_rap(irap),.TRUE.,.FALSE.,.FALSE.)
                  ELSEIF (last_line==ilines .AND. irap==last_rap) THEN
                     CALL gnuplot_write_file_data(filename,'band_lw',&
                              color_rap(irap),.FALSE.,.TRUE.,.FALSE.)
                  ELSE
                     CALL gnuplot_write_file_data(filename,'band_lw',&
                              color_rap(irap),.FALSE.,.FALSE.,.FALSE.)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF  
      ENDDO
      DEALLOCATE(dorap)
   ENDIF
ELSE
   CALL gnuplot_print_objects()
END IF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_dispersion

