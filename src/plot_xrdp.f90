!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_xrdp(exten)
!
!  This is a driver to plot the X-ray diffraction intensities as a 
!  function of the scattering angle 2 theta
!  
!
USE kinds,           ONLY : DP
USE control_gnuplot, ONLY : flgnuplot, lgnuplot, gnuplot_command
USE control_xrdp,    ONLY : flxrdp, flpsxrdp
USE gnuplot,         ONLY : gnuplot_start, gnuplot_end,  &
                            gnuplot_write_header,        &
                            gnuplot_ylabel,              &
                            gnuplot_xlabel,              &
                            gnuplot_write_command,       &
                            gnuplot_line
USE control_pressure, ONLY : pressure, pressure_kb
USE mp_images,       ONLY : my_image_id, root_image, intra_image_comm
USE mp,              ONLY : mp_bcast
USE io_global,       ONLY : ionode, ionode_id

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: exten

CHARACTER(LEN=256) :: filename, gnu_filename
CHARACTER(LEN=8) :: float_to_char
INTEGER :: system
INTEGER :: ierr, ios
INTEGER, PARAMETER :: data_max=10000
INTEGER :: miller(3,data_max), multiplicity(data_max)
REAL(DP) :: theta2(data_max), dist(data_max), intensity(data_max), x(2), y(2)
INTEGER :: ndata, idata

IF ( my_image_id /= root_image ) RETURN

filename=TRIM(flxrdp)//TRIM(exten)
IF (pressure /= 0.0_DP) &
   filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

IF (ionode) OPEN (UNIT=47, FILE=TRIM(filename), STATUS='old', &
                                    FORM='formatted',ERR=100,IOSTAT=ios)
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
CALL errore('plot_xrdp','problem opening the data file',ios)

IF (ionode) THEN
   READ(47,*)
   READ(47,*)
   ndata=0
   DO idata=1,data_max
      READ(47,'(3i5,5x,i4,f16.4,f18.8,f12.2)',END=200) miller(:,idata), &
                  multiplicity(idata), theta2(idata), dist(idata), & 
                  intensity(idata)
      ndata=ndata+1
   ENDDO
200 CONTINUE
   CLOSE (47)
ENDIF
CALL mp_bcast(ndata,ionode_id,intra_image_comm)
CALL mp_bcast(miller(3,1:ndata),ionode_id,intra_image_comm)
CALL mp_bcast(multiplicity(1:ndata),ionode_id,intra_image_comm)
CALL mp_bcast(theta2(1:ndata),ionode_id,intra_image_comm)
CALL mp_bcast(dist(1:ndata),ionode_id,intra_image_comm)
CALL mp_bcast(intensity(1:ndata),ionode_id,intra_image_comm)

gnu_filename=TRIM(flgnuplot)//TRIM(exten)//'_xrdp'
IF (pressure /= 0.0_DP) gnu_filename=TRIM(gnu_filename)//'.'//&
                                     TRIM(float_to_char(pressure_kb,1))
CALL gnuplot_start(gnu_filename)
filename=TRIM(flpsxrdp)//TRIM(exten)
IF (pressure /= 0.0_DP) &
   filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
CALL gnuplot_write_header(TRIM(filename), 0.0_DP, 180.0_DP, 0.0_DP, &
                                                  125.0_DP, 1.0_DP) 

CALL gnuplot_xlabel('2 {/Symbol q} (degrees)',.FALSE.) 
CALL gnuplot_ylabel('Intensity (arbitrary units)',.FALSE.) 
CALL gnuplot_write_command('set xtics 25',.FALSE.)
CALL gnuplot_write_command('set ytics 25',.FALSE.)

DO idata=1,ndata
   x(1)=theta2(idata)
   y(1)=0.0_DP
   x(2)=x(1)
   y(2)=intensity(idata)
   CALL gnuplot_line(x, y, '4', 'front', 'color_blue')
END DO

CALL gnuplot_write_command('plot x+1000',.FALSE.)

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_xrdp
