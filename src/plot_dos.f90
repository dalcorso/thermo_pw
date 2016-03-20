!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_dos( )
!
USE kinds, ONLY : DP
USE constants, ONLY : rytoev
USE gnuplot,   ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                      gnuplot_xlabel, gnuplot_ylabel, gnuplot_write_command, &
                      gnuplot_write_file_mul_data
USE control_gnuplot, ONLY : gnuplot_command, lgnuplot, flgnuplot
USE control_dos,     ONLY : save_ndos
USE io_files,   ONLY : prefix
USE data_files, ONLY : fldos
USE postscript_files, ONLY : flpsdos
USE ener,             ONLY : ef
USE klist,     ONLY : degauss
USE ktetra,    ONLY : ltetra
USE lsda_mod,  ONLY : nspin
USE io_global, ONLY : ionode
USE mp_images, ONLY : my_image_id, root_image

IMPLICIT NONE

INTEGER  :: ierr

CHARACTER(LEN=256) :: gnu_filename, filename, fildos, ylabel, xlabel
REAL(DP), ALLOCATABLE :: e(:), dos(:), int_dos(:)
INTEGER :: n
INTEGER :: iu_dos
REAL(DP) :: ymax, ymin, ymin1, ymax1

IF ( my_image_id /= root_image ) RETURN

ALLOCATE(e(save_ndos))
ALLOCATE(dos(save_ndos))
ALLOCATE(int_dos(save_ndos))

IF ( ionode ) THEN
   IF ( fldos == ' ' ) THEN
      fildos = trim(prefix)//'.dos'
   ELSE
      fildos=TRIM(fldos)
   ENDIF

   iu_dos=2
   OPEN (unit=iu_dos, file=fildos, status='unknown', form='formatted')

   READ(iu_dos,*)
   DO n=1, save_ndos
      READ(iu_dos,*) e(n), dos(n), int_dos(n)
   ENDDO
   CLOSE(iu_dos)
END IF

ymin=1.D10
ymax=0.0_DP
DO n=1,save_ndos
   IF (dos(n) > ymax) ymax=dos(n)
   IF (dos(n) < ymin) ymin=dos(n)
END DO
ymax=ymax*1.1_DP

gnu_filename = TRIM(flgnuplot)//'_dos'
ymin1=1.D10
ymax1=0.0_DP
DO n=1,save_ndos
   IF (int_dos(n) > ymax1) ymax1=int_dos(n)
   IF (int_dos(n) < ymin1) ymin1=int_dos(n)
END DO
ymax1=ymax1*1.1_DP

gnu_filename = TRIM(flgnuplot)//'_dos'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpsdos)

CALL gnuplot_write_header(filename, e(1), e(save_ndos), ymin, ymax, 1.0_DP )

xlabel='Energy (eV)'
IF (nspin==2) THEN
   ylabel='dos (states / (eV spin cell))'
ELSE
   ylabel='dos (states / (eV cell) )'
ENDIF

CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
CALL gnuplot_xlabel(TRIM(xlabel), .FALSE.)
CALL gnuplot_write_command('plot_width=2',.FALSE.)

IF (degauss>0.0_DP.OR.ltetra) THEN
   WRITE(ylabel,'("set arrow from ",f13.5,",",f13.5," to ",f13.5,",",&
                 &f13.5," nohead lw 2")') ef*rytoev,ymin,ef*rytoev,ymax
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   WRITE(ylabel,'("set label ""E_F"" at",f13.5,",",f13.5)') ef*rytoev*1.05_DP, &
                                                                  ymax*0.9_DP
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
END IF

CALL gnuplot_write_file_mul_data(fldos,1,2,'color_red',.TRUE.,.TRUE.,.FALSE.)
IF (nspin==2) THEN
   CALL gnuplot_write_file_mul_data(fldos,1,3,'color_green',.TRUE.,.TRUE.,.FALSE.)
   ylabel='Integrated dos (states/cell)'
   CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
   WRITE(ylabel,'("set yrange[",f13.6,":",f13.6,"]")') ymin1, ymax1
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   CALL gnuplot_write_command('unset arrow',.FALSE.)
   CALL gnuplot_write_command('unset label',.FALSE.)
   CALL gnuplot_write_file_mul_data(fldos,1,4,'color_blue',.TRUE.,.TRUE.,.FALSE.)
ELSE
   ylabel='Integrated dos (states/cell)'
   CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
   WRITE(ylabel,'("set yrange[",f13.6,":",f13.6,"]")') ymin1, ymax1
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   CALL gnuplot_write_command('unset arrow',.FALSE.)
   CALL gnuplot_write_command('unset label',.FALSE.)
   CALL gnuplot_write_file_mul_data(fldos,1,3,'color_blue',.TRUE.,.TRUE.,.FALSE.)
ENDIF

CALL gnuplot_end()

DEALLOCATE(e)
DEALLOCATE(dos)
DEALLOCATE(int_dos)

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_dos
