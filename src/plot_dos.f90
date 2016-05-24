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
                      gnuplot_write_file_mul_data, &
                      gnuplot_write_file_mul_data_minus
USE control_gnuplot, ONLY : gnuplot_command, lgnuplot, flgnuplot
USE control_eldos,   ONLY : save_ndos
USE control_bands,   ONLY : emin_input
USE io_files,   ONLY : prefix
USE data_files, ONLY : fleldos
USE postscript_files, ONLY : flpseldos
USE ener,             ONLY : ef
USE klist,     ONLY : degauss
USE ktetra,    ONLY : ltetra
USE lsda_mod,  ONLY : nspin
USE io_global, ONLY : ionode
USE mp_images, ONLY : my_image_id, root_image

IMPLICIT NONE

INTEGER  :: ierr

CHARACTER(LEN=256) :: gnu_filename, filename, fileeldos, ylabel, xlabel
REAL(DP), ALLOCATABLE :: e(:), dos(:), ddos(:), int_dos(:)
INTEGER :: n
INTEGER :: iu_dos
REAL(DP) :: ymax, ymin, ymin1, ymax1, e1

IF ( my_image_id /= root_image ) RETURN

ALLOCATE(e(save_ndos))
ALLOCATE(dos(save_ndos))
IF (nspin==2) ALLOCATE(ddos(save_ndos))
ALLOCATE(int_dos(save_ndos))

fileeldos='therm_files/'//TRIM(fleldos)
IF ( ionode ) THEN
   iu_dos=2
   OPEN (unit=iu_dos, file=TRIM(fileeldos), status='unknown', form='formatted')

   READ(iu_dos,*)
   IF (nspin==2) THEN
      DO n=1, save_ndos
         READ(iu_dos,*) e(n), dos(n), ddos(n), int_dos(n)
      ENDDO
   ELSE
      DO n=1, save_ndos
         READ(iu_dos,*) e(n), dos(n), int_dos(n)
      ENDDO
   END IF
   CLOSE(iu_dos)
END IF

ymin=1.D10
ymax=0.0_DP
DO n=1,save_ndos
   IF (dos(n) > ymax) ymax=dos(n)
   IF (dos(n) < ymin) ymin=dos(n)
   IF (nspin==2) THEN
      IF (ddos(n) > ymax) ymax=ddos(n)
      IF (ddos(n) < ymin) ymin=ddos(n)
   END IF 
END DO
ymax=ymax*1.1_DP

gnu_filename = 'gnuplot_files/'//TRIM(flgnuplot)//'_dos'
ymin1=1.D10
ymax1=0.0_DP
DO n=1,save_ndos
   IF (int_dos(n) > ymax1) ymax1=int_dos(n)
   IF (int_dos(n) < ymin1) ymin1=int_dos(n)
END DO
ymax1=ymax1*1.1_DP

gnu_filename = 'gnuplot_files/'//TRIM(flgnuplot)//'_eldos'
CALL gnuplot_start(gnu_filename)

filename=TRIM(flpseldos)

xlabel='Energy (eV)'
e1=e(1)
IF (emin_input /= 0.0_DP) e1=emin_input
IF (nspin==2) THEN
   CALL gnuplot_write_header(filename, e1, e(save_ndos), -ymax, ymax, 1.0_DP )
   ylabel='dos (states / (spin  eV  cell) )'
ELSE
   CALL gnuplot_write_header(filename, e1, e(save_ndos), ymin, ymax, 1.0_DP )
   ylabel='dos (states / (eV cell))'
ENDIF

CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
CALL gnuplot_xlabel(TRIM(xlabel), .FALSE.)
CALL gnuplot_write_command('plot_width=2',.FALSE.)

IF (degauss>0.0_DP.OR.ltetra) THEN
   IF (nspin==2) THEN
      WRITE(ylabel,'("set arrow from ",f13.5,",-ymax to ",f13.5,&
                      &", ymax nohead lw 2")') ef*rytoev,ef*rytoev
   ELSE
      WRITE(ylabel,'("set arrow from ",f13.5,",",f13.5," to ",f13.5,",",&
                 &f13.5," nohead lw 2")') ef*rytoev,ymin,ef*rytoev,ymax
   END IF
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   WRITE(ylabel,'("set label ""E_F"" at",f13.5,",",f13.5)') &
                                    ef*rytoev*1.05_DP, ymax*0.92_DP
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   IF (nspin==2) THEN 
      WRITE(ylabel,'("set arrow from xmin,0.0 to xmax, 0.0 nohead lw 2")')
      CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   END IF
END IF

IF (nspin==2) THEN
   CALL gnuplot_write_file_mul_data(fileeldos,1,2,'color_red',.TRUE.,.FALSE.,.FALSE.)
   CALL gnuplot_write_file_mul_data_minus(fileeldos,1,3,'color_blue',.FALSE., &
                                                  .TRUE., .FALSE.)
   ylabel='Integrated dos (states / cell)'
   CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
   WRITE(ylabel,'("set yrange[",f13.6,":",f13.6,"]")') ymin1, ymax1
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   CALL gnuplot_write_command('unset arrow',.FALSE.)
   CALL gnuplot_write_command('unset label',.FALSE.)
   CALL gnuplot_write_file_mul_data(fileeldos,1,4,'color_blue',.TRUE.,.TRUE.,&
                                                                      .FALSE.)
ELSE
   CALL gnuplot_write_file_mul_data(fileeldos,1,2,'color_red',.TRUE.,.TRUE.,&
                                                                      .FALSE.)
   ylabel='Integrated dos (states / cell)'
   CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
   WRITE(ylabel,'("set yrange[",f13.6,":",f13.6,"]")') ymin1, ymax1
   CALL gnuplot_write_command(TRIM(ylabel),.FALSE.)
   CALL gnuplot_write_command('unset arrow',.FALSE.)
   CALL gnuplot_write_command('unset label',.FALSE.)
   CALL gnuplot_write_file_mul_data(fileeldos,1,3,'color_blue',.TRUE.,.TRUE.,&
                                                                      .FALSE.)
ENDIF

CALL gnuplot_end()

DEALLOCATE(e)
DEALLOCATE(dos)
IF (nspin==2) DEALLOCATE(ddos)
DEALLOCATE(int_dos)

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_dos
