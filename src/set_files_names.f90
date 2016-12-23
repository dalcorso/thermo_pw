! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_files_names(igeom)

USE postscript_files,  ONLY : flpstherm, flpsdisp, flpsdos
USE data_files,        ONLY : fltherm, flfrc, flfrq, fldos, flpband, flvec, &
                              fldosfrq
USE control_gnuplot,   ONLY : flgnuplot
USE internal_files_names, ONLY : fildyn_thermo, flfrc_thermo, flfrq_thermo, &
                       fldos_thermo, fltherm_thermo, flpband_thermo,       &
                       flpsdos_thermo, flpstherm_thermo, flgnuplot_thermo, &
                       flpsdisp_thermo, flvec_thermo, fldosfrq_thermo
USE output,            ONLY : fildyn
  !
IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=6) :: int_to_char

IF (igeom==1) fildyn_thermo="dynamical_matrices/"//TRIM(fildyn)
IF (igeom==1) flfrc_thermo=TRIM(flfrc)
IF (igeom==1) flfrq_thermo=TRIM(flfrq)
IF (igeom==1) flvec_thermo=TRIM(flvec)
IF (igeom==1) fldos_thermo=TRIM(fldos)
IF (igeom==1) fldosfrq_thermo=TRIM(fldosfrq)
IF (igeom==1) fltherm_thermo=TRIM(fltherm)
IF (igeom==1) flpband_thermo=TRIM(flpband)

IF (igeom==1) flgnuplot_thermo=TRIM(flgnuplot)

IF (igeom==1) flpsdos_thermo=TRIM(flpsdos)
IF (igeom==1) flpstherm_thermo=TRIM(flpstherm)
IF (igeom==1) flpsdisp_thermo=TRIM(flpsdisp)

fildyn=TRIM(fildyn_thermo)//'.g'//TRIM(int_to_char(igeom))//'.'
flfrc=TRIM(flfrc_thermo)//'.g'//TRIM(int_to_char(igeom))
flfrq=TRIM(flfrq_thermo)//'.g'//TRIM(int_to_char(igeom))
flvec=TRIM(flvec_thermo)//'.g'//TRIM(int_to_char(igeom))
fldos=TRIM(fldos_thermo)//'.g'//TRIM(int_to_char(igeom))
fldosfrq=TRIM(fldosfrq_thermo)//'.g'//TRIM(int_to_char(igeom))
fltherm=TRIM(fltherm_thermo)//'.g'//TRIM(int_to_char(igeom))
flpband=TRIM(flpband_thermo)//'.g'//TRIM(int_to_char(igeom))

flgnuplot=TRIM(flgnuplot_thermo)//'.g'//TRIM(int_to_char(igeom))

flpsdos=TRIM(flpsdos_thermo)//'.g'//TRIM(int_to_char(igeom))
flpstherm=TRIM(flpstherm_thermo)//'.g'//TRIM(int_to_char(igeom))
flpsdisp=TRIM(flpsdisp_thermo)//'.g'//TRIM(int_to_char(igeom))

RETURN
END SUBROUTINE set_files_names

!----------------------------------------------------------------------------
SUBROUTINE restore_files_names()
!----------------------------------------------------------------------------

USE postscript_files,  ONLY : flpstherm, flpsdisp, flpsdos
USE data_files,        ONLY : fltherm, flfrc, flfrq, fldos,  flpband, fldosfrq 
USE internal_files_names, ONLY : fildyn_thermo, flfrc_thermo, flfrq_thermo, &
                        fldos_thermo, fltherm_thermo, flpband_thermo,       &
                        flpsdos_thermo, flpstherm_thermo, flgnuplot_thermo, &
                        flpsdisp_thermo, fldosfrq_thermo
USE control_gnuplot,   ONLY : flgnuplot
USE output,            ONLY : fildyn
  !
IMPLICIT NONE

fildyn = TRIM(fildyn_thermo)
flfrc = TRIM(flfrc_thermo)
flfrq = TRIM(flfrq_thermo)
fldos = TRIM(fldos_thermo)
fldosfrq = TRIM(fldosfrq_thermo)
fltherm = TRIM(fltherm_thermo)
flpband = TRIM(flpband_thermo)

flgnuplot =TRIM(flgnuplot_thermo)

flpsdos = TRIM(flpsdos_thermo)
flpstherm = TRIM(flpstherm_thermo)
flpsdisp = TRIM(flpsdisp_thermo)

RETURN
END SUBROUTINE restore_files_names

SUBROUTINE set_files_for_plot(icode, file_disp, filedata, filerap, fileout, &
                                     gnu_filename, filenameps)
!
!   This routine receives as input a code of what we want to plot
!   icode  1  band structure
!   icode  2  phonon
!   icode  3  gruneisen parameters
!   icode  4  phonon interpolated from the fit of the frequencies
!
!   and sets the names of the files:
!   filedata  the file that contains the data to plot
!   fileout   the file where the plotting routine writes the data in
!             a format readable by gnuplot
!   gnu_filename the name of the file that will contain the gnuplot script
!   filenameps the name of the file that will contain the postscript picture
!   The names of the files are built using the part of the filename 
!   contained in the thermo_pw modules and are provided on the dummy veriables
!

USE thermo_mod,       ONLY : tot_ngeo, central_geo
USE control_thermo,   ONLY : spin_component
USE lsda_mod,         ONLY : nspin

USE data_files,       ONLY : flpgrun, flpband, filband, flfrq, flgrun
USE postscript_files, ONLY : flpsband, flpsdisp, flpsgrun
USE control_gnuplot,  ONLY : flgnuplot

USE io_global,        ONLY : stdout

IMPLICIT NONE

INTEGER, INTENT(IN) :: icode
CHARACTER(LEN=256), INTENT(IN) :: file_disp
CHARACTER(LEN=256), INTENT(OUT) :: filedata, filerap, fileout, gnu_filename, &
                                   filenameps 
CHARACTER(LEN=6) :: int_to_char
!
!  first the file with the data 
!
  IF (icode==1) THEN
     filedata = "band_files/"//TRIM(filband)
     IF (nspin==2) &
        filedata = "band_files/"//TRIM(filband)// &
                              '.'//TRIM(int_to_char(spin_component))
  ELSEIF (icode==2) THEN
     filedata = "phdisp_files/"//TRIM(flfrq)
  ELSEIF (icode==3) THEN
     filedata = "anhar_files/"//TRIM(flgrun)
  ELSEIF (icode==4) THEN
     filedata = "anhar_files/"//TRIM(flgrun)//'_freq'
  END IF
!
!  the the name of the files were the symmetry information is found
!
  IF (icode==1.OR.icode==2) THEN
     filerap=TRIM(filedata)//".rap"
  ELSEIF (icode==3.OR.icode==4) THEN
     filerap = "phdisp_files/"//TRIM(file_disp)//'.g'//&
                                        TRIM(int_to_char(central_geo))//".rap"
  ENDIF
!
!  then the name of the file where the data are written
!
  IF (icode==1.OR.icode==2) THEN
     IF (flpband == ' ') THEN
        fileout=' '
     ELSEIF (icode==1) THEN
        fileout="band_files/"//TRIM(flpband)
        IF (nspin==2) &
           fileout="band_files/"//TRIM(flpband)//"."//&
                                  TRIM(int_to_char(spin_component))
     ELSEIF (icode==2) THEN
        fileout="phdisp_files/"//TRIM(flpband)
     ENDIF
  ELSEIF (icode==3) THEN
     IF (flpgrun == ' ') THEN
        fileout=' '
     ELSE
        fileout="anhar_files/"//TRIM(flpgrun)
     ENDIF
  ELSEIF (icode==4) THEN
     IF (flpgrun == ' ') THEN
        fileout=' '
     ELSE
        fileout="anhar_files/"//TRIM(flpgrun)//'_freq'
     ENDIF
  ENDIF
!
!  then the name of the file with the gnuplot script
!
  IF (icode==1) THEN
     gnu_filename=TRIM(flgnuplot)//'_band'
     IF (nspin==2) gnu_filename=TRIM(flgnuplot)//'_band.'&
                   //TRIM(int_to_char(spin_component))
  ELSEIF (icode==2) THEN
     gnu_filename=TRIM(flgnuplot)//'_disp'
  ELSEIF (icode==3) THEN
     gnu_filename=TRIM(flgnuplot)//'_grun'
  ELSEIF (icode==4) THEN
     gnu_filename=TRIM(flgnuplot)//'_grun_freq'
  ENDIF
!
!  And the name of the postscript file
!
  IF (icode==1) THEN
     filenameps=TRIM(flpsband)
     IF (nspin==2.AND.spin_component==1) filenameps=TRIM(flpsband)//'_up'
     IF (nspin==2.AND.spin_component==2) filenameps=TRIM(flpsband)//'_down'
  ELSEIF (icode==2) THEN
     filenameps=TRIM(flpsdisp)
  ELSEIF (icode==3) THEN
     filenameps=TRIM(flpsgrun)
  ELSEIF (icode==4) THEN
     filenameps=TRIM(flpsgrun)//'_freq'
  ENDIF

  gnu_filename="gnuplot_files/"//TRIM(gnu_filename)

  IF (icode==1) THEN
     WRITE(stdout,'(5x,"Writing bands in gnuplot format to file ",a)') &
                                                            TRIM(fileout)
  ELSEIF (icode==2) THEN
     WRITE(stdout,'(5x,"Writing phonons in gnuplot format to file ",a)') &
                                                            TRIM(fileout)
  ELSEIF (icode==3) THEN
     WRITE(stdout,'(5x,"Writing Gruneisen parameters in gnuplot format &
                                            &to file ",a)') TRIM(fileout)
  ELSEIF (icode==4) THEN
     WRITE(stdout,'(5x,"Writing interpolated phonons in gnuplot format  &
                                            &to file ",a)') TRIM(fileout)
  ENDIF


RETURN
END SUBROUTINE set_files_for_plot

