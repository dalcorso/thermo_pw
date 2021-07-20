!
! Copyright (C) 2020 Cristiano Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
PROGRAM pdec
!------------------------------------------------------------------------

! This program reads the temperature dependent elastic constants (ECs)
! files for each pressure in the anhar_files directory and writes in 
! the directory pressure_files for each temperature a file with the 
! pressure dependent ECs and plots them.
! The user have to provide the filename in input (only the part without 
! the value of the pressure).
! The output file names (data files and plots) are identified by the 
! input filename plus the temperature.
! The Laue class index (required in input) follows the same covention 
! used in thermo_pw routines (see for instance print_el_cons_shape 
! subroutine in lib/nye.f90 file).  

USE kinds,   ONLY : DP
USE mp_global,    ONLY : mp_startup, mp_global_end
USE environment,  ONLY : environment_start, environment_end
USE clib_wrappers,  ONLY : f_mkdir_safe
USE elastic_constants, ONLY : read_el_cons_from_file
USE nye,     ONLY : print_el_cons_shape
USE polyfit_mod, ONLY : polyfit, compute_poly
USE io_global,    ONLY : stdout

IMPLICIT NONE

CHARACTER(LEN=16)   :: code='pdec'

INTEGER :: ibrav, laue, ntemp, npress, ntemp_press
INTEGER :: ios, iupress, lines
INTEGER :: len_path, ipress, var, itemp, jtemp, ipol, jpol
REAL(DP), ALLOCATABLE :: p(:), temp(:), temp_press(:)
REAL(DP), ALLOCATABLE :: el_cons_t_press(:,:,:,:), b0_t_press(:,:)
REAL(DP), ALLOCATABLE :: el_cons_tinput_p(:,:,:,:), b0_tinput_p(:,:)
INTEGER(DP), ALLOCATABLE :: ind(:)
REAL(DP) :: tmin, tmax, deltat
REAL(DP) :: a(2)
INTEGER  :: system, ierr, stdin
CHARACTER (LEN=256) :: input_filename, filename_path, partial_filename, str_tmp, rdum_str
CHARACTER (LEN=256) :: file_press, filepressure, filelastic, filefirstpressure
CHARACTER(LEN=8) :: real_to_char

LOGICAL :: exst, pres, midl

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

ios = f_mkdir_safe( 'pressure_files' )
stdin=5
WRITE(stdout,*)
WRITE(stdout,'(5x,"Insert the part of the elastic constants file name before the pressure:")')
READ(stdin,*) input_filename

WRITE(stdout,'(5x,"Bravais index of your solid ? ")')
READ(stdin,*) ibrav

WRITE(stdout,'(5x,"Laue class of your solid ? ")')
READ(stdin,*) laue

WRITE(stdout,'(5x,"The shape of elastic constant tensor of your system is:")')
CALL print_el_cons_shape(laue,ibrav)

filename_path = './anhar_files/'//TRIM(input_filename)

INQUIRE(FILE=TRIM(filename_path),EXIST=exst)

ierr=system('ls '//TRIM(filename_path)//'.* > ec_files_pressure.txt')

file_press = 'ec_files_pressure.txt'

CALL lines_number(file_press, npress)

IF (exst) npress=npress+1

WRITE(stdout,*)
WRITE(stdout,'(5x,"The number of pressures considered is",I6)') npress

ALLOCATE(p(npress))
ALLOCATE(ind(npress)) 

p=0.0_DP
ind=0

OPEN (UNIT=iupress, FILE=TRIM(file_press), ERR=50, IOSTAT=ios)

50 CALL errore('pdec','opening file', ABS(ios))

DO ipress=1,npress
   IF (exst.AND.ipress==1) THEN 
      p(ipress)=0.0_DP
      CYCLE
   ENDIF 

   READ(iupress, '(A)', IOSTAT=ios) partial_filename
   str_tmp = partial_filename(LEN_TRIM(filename_path)+2:)
   READ(str_tmp,*) p(ipress) 
END DO

CLOSE(UNIT=iupress)

!Sort the array of pressures (required in the case negative pressures
!are present)

CALL hpsort(npress, p, ind)

WRITE(stdout,'(5x,"The minimum pressure considered is",F6.1," kbar")') p(1)
WRITE(stdout,'(5x,"The maximum pressure considered is",F6.1," kbar")') p(npress)
WRITE(stdout,*)

filefirstpressure=filename_path//TRIM(real_to_char(p(1)))

CALL lines_number(filefirstpressure, ntemp)
WRITE(stdout,'(5x,"Number of temperature",I6)') ntemp

ntemp=ntemp+1

OPEN (UNIT=iupress, FILE=TRIM(filefirstpressure), ERR=30, IOSTAT=ios)

30 CALL errore('pdec','opening elastic constants (T) file', ABS(ios)) 

READ(iupress,*)
READ(iupress, '(e16.8, A)', IOSTAT=ios) tmin, rdum_str
READ(iupress, '(e16.8, A)', IOSTAT=ios) deltat, rdum_str 

CLOSE(UNIT=iupress)

deltat=deltat-tmin
tmin=tmin-deltat
tmax=tmin+(ntemp-1)*deltat

ALLOCATE( temp(ntemp) )
ALLOCATE( b0_t_press(ntemp,npress) )
ALLOCATE( el_cons_t_press(6,6,ntemp,npress) )

temp(ntemp)=0.0_DP

DO itemp=1, ntemp
   temp(itemp)= tmin + (itemp-1)*deltat
   !WRITE(stdout,'(5x,"Temperature ",f15.8)') temp(itemp)
ENDDO

WRITE(stdout,*)
WRITE(stdout,'(5x,"Number of temperatures for the elastic constants vs pressure")')
WRITE(stdout,'(5x,"plots (select a value from 1 to 10):")')
READ(stdin,*) ntemp_press

IF (ntemp_press<1 .OR. ntemp_press>10) THEN
   WRITE(stdout,'(5x,"Incorrect number of temperatures.")')
   STOP 1
ENDIF

ALLOCATE(temp_press(ntemp_press))

DO itemp=1, ntemp_press
   WRITE(stdout,'(5x,"Temperature",I2," (K):")') itemp
   READ(stdin,*) temp_press(itemp)
   IF (temp_press(itemp)<tmin.OR.temp_press(itemp)>tmax) THEN
      WRITE(stdout,'(5x,"Invalid temperature")')
      STOP 1
   END IF
END DO

WRITE(stdout,'(5x,"Temperature",F6.1,":")') temp_press

el_cons_t_press=0.0_DP
b0_t_press=0.0_DP

DO ipress=1, npress
   filepressure=TRIM(filename_path)//'.'//real_to_char(p(ipress))
   IF (filepressure==TRIM(filename_path)//'.0.0') &
                                 filepressure=TRIM(filename_path)
   CALL read_el_cons_from_file(temp, ntemp, ibrav, laue, &
     el_cons_t_press(:,:,:,ipress), b0_t_press(:,ipress), filepressure)
ENDDO

ALLOCATE( el_cons_tinput_p(6,6,ntemp_press,npress) )
ALLOCATE( b0_tinput_p(ntemp_press,npress) )

el_cons_tinput_p=0.0_DP
b0_tinput_p=0.0_DP

DO itemp=1, ntemp_press
   DO jtemp=1, ntemp-1
      pres=(temp_press(itemp)==temp(jtemp))
      midl=(temp_press(itemp)>temp(jtemp).AND.temp_press(itemp)<temp(jtemp+1))
      IF (.NOT.pres.AND..NOT.midl) THEN
         CYCLE
      ELSE IF (pres) THEN
         el_cons_tinput_p(:,:,itemp,:)=el_cons_t_press(:,:,jtemp,:)
         b0_tinput_p(itemp,:)=b0_t_press(jtemp,:)
      ELSE IF (midl) THEN
         DO ipress=1, npress
            CALL polyfit(temp(jtemp:jtemp+1),b0_t_press(jtemp:jtemp+1,ipress),2,a,1) 
            CALL compute_poly(temp_press(itemp),1,a,b0_tinput_p(itemp,ipress))
            DO ipol=1, 6
               DO jpol=1, 6
                  IF (el_cons_t_press(ipol,jpol,jtemp,ipress)==0) CYCLE
                  CALL polyfit(temp(jtemp:jtemp+1), &
                               el_cons_t_press(ipol,jpol,jtemp:jtemp+1,ipress),2,a,1)
                  CALL compute_poly(temp_press(itemp),1,a, &
                               el_cons_tinput_p(ipol,jpol,itemp,ipress))
               END DO
            END DO
         END DO
      END IF
   END DO
END DO

! Write the files of the pressure dependent elastic constants at each input
! temperature

DO itemp = 1, ntemp_press
   filelastic='pressure_files/'//TRIM(input_filename)//'.'//real_to_char(temp_press(itemp))
   CALL write_el_cons_p(p, npress, ibrav, laue, el_cons_tinput_p(:,:,itemp,:), &
                        b0_tinput_p(itemp,:), filelastic)
ENDDO

! For each input temperature plot the pressure dependent elastic constants

DO itemp = 1, ntemp_press
   filelastic=TRIM(input_filename)//'.'//real_to_char(temp_press(itemp))
   CALL plot_elastic_p(ibrav, laue, filelastic)
ENDDO

CALL environment_end(code)
CALL mp_global_end ()

END PROGRAM pdec

!------------------------------------------------------------------------
SUBROUTINE lines_number(filename, lines)
!------------------------------------------------------------------------

!Count the number of lines written in a file

USE kinds, ONLY : DP

IMPLICIT NONE
CHARACTER (LEN=256) :: filename
INTEGER :: var, iupress
INTEGER, INTENT(OUT) :: lines

OPEN(UNIT=iupress, FILE=TRIM(filename))
lines = 0
var = 0

DO WHILE (var==0)
   !WRITE(*,*) 'Il valore di var e', var
   READ(iupress, *, iostat=var)
   IF (var/=0) EXIT
   lines = lines + 1
END DO

!WRITE(*,*) 'The number of lines is', lines

CLOSE(UNIT=iupress)

RETURN
END SUBROUTINE lines_number

!------------------------------------------------------------------------
FUNCTION real_to_char(y)
!------------------------------------------------------------------------

USE kinds, ONLY : DP

IMPLICIT NONE
REAL(DP) :: y
CHARACTER(LEN=8) :: real_to_char

WRITE(real_to_char,'(F6.1)') y

real_to_char = ADJUSTL(real_to_char)

RETURN
END FUNCTION real_to_char

!------------------------------------------------------------------------
SUBROUTINE write_el_cons_p(p, npress, ibrav, laue, el_cons_p, b0_p, filename)
!------------------------------------------------------------------------
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress, ibrav, laue
REAL(DP), INTENT(IN) :: p(npress), el_cons_p(6,6,npress), b0_p(npress)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: iu_el_cons, ios, ipress
INTEGER :: find_free_unit

iu_el_cons=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_el_cons, FILE=TRIM(filename), FORM='formatted', &
                                   STATUS='UNKNOWN', ERR=30,IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_el_cons_p','opening elastic constants (P) file',&
                                                             ABS(ios))
IF (meta_ionode) THEN
   SELECT CASE (laue)
      CASE(29,32)
         WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", 13x, " C_11 ", &
                            13x, "     C_12 ", 13x, "   C_44")')
         DO ipress=1,npress
            WRITE(iu_el_cons,'(f6.1,4e20.12)') p(ipress), b0_p(ipress),&
                 el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                 el_cons_p(4,4,ipress)
         ENDDO

      CASE(25)
!
!     D_3d
!            
         WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", 13x, " C_11 ", &
                   13x, " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, &
                                             " C_44 ", 13x, " C_14")')
         DO ipress=1,npress
            WRITE(iu_el_cons,'(e16.8,7e20.12)') p(ipress), b0_p(ipress),&
                  el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                  el_cons_p(1,3,ipress), el_cons_p(3,3,ipress), &
                  el_cons_p(4,4,ipress), el_cons_p(1,4,ipress)
         ENDDO

     CASE(27)
!
!     S_6
!            
         WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", " C_11 ", 13x, &
                   " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44",13x, &
                   " C_14", 13x, "C_25" )')
         DO ipress=1,npress
            WRITE(iu_el_cons,'(e16.8,8e20.12)')  p(ipress), b0_p(ipress), &
                  el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                  el_cons_p(1,3,ipress), el_cons_p(3,3,ipress), &
                  el_cons_p(4,4,ipress), el_cons_p(1,4,ipress), &
                  el_cons_p(2,5,ipress)
         ENDDO

     CASE(19,23)
         WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", 13x, " C_11 ",13x,&
                  " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44")')
         DO ipress=1,npress
            WRITE(iu_el_cons,'(e16.8,6e20.12)') p(ipress), b0_p(ipress), &
                  el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                  el_cons_p(1,3,ipress), el_cons_p(3,3,ipress), &
                  el_cons_p(4,4,ipress)
         ENDDO

     CASE(22)
         WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", 13x, " C_11 ",13x,&
                  " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44",13x,& 
                  " C_66 " )')
         DO ipress=1,npress
            WRITE(iu_el_cons,'(e16.8,7e20.12)') p(ipress), b0_p(ipress),&
                  el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                  el_cons_p(1,3,ipress), el_cons_p(3,3,ipress), &
                  el_cons_p(4,4,ipress), el_cons_p(6,6,ipress)
         ENDDO

     CASE(20)
         WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B", 13x, " C_11 ", 13x,&
                   " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23", 13x,&
                   " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66")')
         DO ipress=1,npress
            WRITE(iu_el_cons,'(e16.8,10e20.12)') p(ipress), b0_p(ipress), &
                  el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                  el_cons_p(1,3,ipress), el_cons_p(2,2,ipress), &
                  el_cons_p(2,3,ipress), el_cons_p(3,3,ipress), &
                  el_cons_p(4,4,ipress), el_cons_p(5,5,ipress), &
                  el_cons_p(6,6,ipress)
         ENDDO

     CASE(18)
         WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", 13x, " C_11 ", &
                  13x, " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x,"C_44", &
                  13x, " C_66 ", 13x, " C_16 ")')
         DO ipress=1,npress
            WRITE(iu_el_cons,'(e16.8,8e20.12)') p(ipress), b0_p(ipress),&
                  el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                  el_cons_p(1,3,ipress), el_cons_p(3,3,ipress), &
                  el_cons_p(4,4,ipress), el_cons_p(6,6,ipress), &
                  el_cons_p(1,6,ipress)
         ENDDO

     CASE(16)
         IF (ibrav < 0) THEN
            !
            !  b unique
            !
            WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", 13x, " C_11", &
              &13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23", &
              &13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66", &
              &13x, " C_15 ", 13x, " C_25 ", 13x, " C_35 ", 13x, " C_46")')
            DO ipress=1,npress
               WRITE(iu_el_cons,'(e16.8,14e20.12)') p(ipress),b0_p(ipress),&
                     el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                     el_cons_p(1,3,ipress), el_cons_p(2,2,ipress), &
                     el_cons_p(2,3,ipress), el_cons_p(3,3,ipress), &
                     el_cons_p(4,4,ipress), el_cons_p(5,5,ipress), &
                     el_cons_p(6,6,ipress), el_cons_p(1,5,ipress), &
                     el_cons_p(2,5,ipress), el_cons_p(3,5,ipress), &
                     el_cons_p(4,6,ipress)
            ENDDO
         ELSE
            !
            !  c unique
            !
            WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ",13x," C_11 ",  &
             &13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23", &
             &13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66", &
             &13x, " C_16 ", 13x, " C_26 ", 13x, " C_36 ", 13x, " C_45")')
            DO ipress=1,npress
               WRITE(iu_el_cons,'(e16.8,14e20.12)') p(ipress), b0_p(ipress),&
                     el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                     el_cons_p(1,3,ipress), el_cons_p(2,2,ipress), &
                     el_cons_p(2,3,ipress), el_cons_p(3,3,ipress), &
                     el_cons_p(4,4,ipress), el_cons_p(5,5,ipress), &
                     el_cons_p(6,6,ipress), el_cons_p(1,6,ipress), &
                     el_cons_p(2,6,ipress), el_cons_p(3,6,ipress), &
                     el_cons_p(4,5,ipress)
            ENDDO
         ENDIF

     CASE(2)
        WRITE(iu_el_cons,'("#",2x,"P  ", 10x, " B ", 13x, " C_11 ", &
               &13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, "C_23 ", &
               &13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, "C_66 ", &
               &13x, " C_14 ", 13x, " C_15 ", 13x, " C_16 ", 13x, "C_24 ", &
               &13x, " C_25 ", 13x, " C_26 ", 13x, " C_34 ", 13x, "C_35 ", &
               &13x, " C_36 ", 13x, " C_45 ", 13x, " C_46 ", 13x, "C_56 ")')
        DO ipress=1,npress
            WRITE(iu_el_cons,'(e16.8,24e20.12)') p(ipress),b0_p(ipress), &
                  el_cons_p(1,1,ipress), el_cons_p(1,2,ipress), &
                  el_cons_p(1,3,ipress), el_cons_p(2,2,ipress), &
                  el_cons_p(2,3,ipress), el_cons_p(3,3,ipress), &
                  el_cons_p(4,4,ipress), el_cons_p(5,5,ipress), &
                  el_cons_p(6,6,ipress), el_cons_p(1,4,ipress), &
                  el_cons_p(1,5,ipress), el_cons_p(1,6,ipress), &
                  el_cons_p(2,4,ipress), el_cons_p(2,5,ipress), &
                  el_cons_p(2,6,ipress), el_cons_p(3,4,ipress), &
                  el_cons_p(3,5,ipress), el_cons_p(3,6,ipress), &
                  el_cons_p(4,5,ipress), el_cons_p(4,6,ipress), &
                  el_cons_p(5,6,ipress)
        ENDDO

   END SELECT
   CLOSE(iu_el_cons)
ENDIF

RETURN
END SUBROUTINE write_el_cons_p

!------------------------------------------------------------------------
SUBROUTINE plot_elastic_p(ibrav, laue, filename)
!------------------------------------------------------------------------
!
!  This is a driver to plot the elastic constants as a function of
!  temperature
!
USE kinds,            ONLY : DP
USE gnuplot,          ONLY : gnuplot_start, gnuplot_end,           &
                             gnuplot_write_header, gnuplot_xlabel, &
                             gnuplot_ylabel, gnuplot_set_fact, &
                             gnuplot_write_file_mul_point
USE io_global, ONLY : ionode

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, laue
CHARACTER (LEN=256), INTENT(IN) :: filename
CHARACTER(LEN=256) :: gnuplot_command, gnu_filename, filenameps, filelastic
INTEGER :: ierr, system

gnuplot_command='gnuplot'

gnu_filename="pressure_files/gnuplot_tmp."//TRIM(filename)
filenameps="pressure_files/"//TRIM(filename)//'.ps'
filelastic="pressure_files/"//TRIM(filename)

CALL gnuplot_start(gnu_filename)

CALL gnuplot_write_header(filenameps, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                                                        1.0_DP, '.ps' )
CALL gnuplot_xlabel('P (kbar)', .FALSE.)

CALL gnuplot_ylabel('Bulk modulus (kbar)',.FALSE.)
CALL gnuplot_write_file_mul_point(filelastic,1,2,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)

CALL gnuplot_ylabel('C_{11} (kbar)',.FALSE.)
CALL gnuplot_write_file_mul_point(filelastic,1,3,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)

CALL gnuplot_ylabel('C_{12} (kbar)',.FALSE.)
CALL gnuplot_write_file_mul_point(filelastic,1,4,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)

IF (laue==32.OR.laue==29) THEN
   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,5,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)
ENDIF

IF (laue==2.OR.laue==18.OR.laue==19.OR.laue==20.OR.laue==22.OR.laue==23.OR.&
    laue==25.OR.laue==27) THEN
!
!  tetragonal, hexagonal, trigonal, or orthorhombic
!
   CALL gnuplot_ylabel('C_{13} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,5,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)
ENDIF

IF (laue==18.OR.laue==22.OR.laue==19.OR.laue==23.OR.laue==25.OR.laue==27) THEN
!
!  tetragonal or hexagonal
!
   CALL gnuplot_ylabel('C_{33} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,6,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,7,'color_red',.TRUE.,&
                                                     .TRUE.,.FALSE.)
ENDIF

IF (laue==25.OR.laue==27) THEN
!
!  trigonal D_3d or S_6
!
   CALL gnuplot_ylabel('C_{14} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,8,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)
ENDIF

IF (laue==27) THEN
!
!  trigonal S_6
!
   CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,9,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)
ENDIF

IF (laue==18.OR.laue==22) THEN
!
!  tetragonal C_4h or D_4h
!
   CALL gnuplot_ylabel('C_{66} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,9,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)
ENDIF

IF (laue==18) THEN
!
!  tetragonal C_4h
!
   CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,9,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)
ENDIF

IF (laue==2.OR.laue==16.OR.laue==20) THEN
!
!  triclinic C_i, monoclinic C_2h, or orthorhombic D_2h
!

   CALL gnuplot_ylabel('C_{22} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,6,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{23} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,7,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{33} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,8,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{44} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,9,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{55} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,10,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)

   CALL gnuplot_ylabel('C_{66} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,11,'color_red',.TRUE.,&
                                                    .TRUE.,.FALSE.)
END IF

IF (laue==16) THEN
   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{15} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,12,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ELSE
      CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,12,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ENDIF

   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,13,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ELSE
      CALL gnuplot_ylabel('C_{26} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,13,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ENDIF

   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{35} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,14,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ELSE
      CALL gnuplot_ylabel('C_{35} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,14,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ENDIF

   IF (ibrav>0) THEN
      CALL gnuplot_ylabel('C_{46} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,15,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ELSE
      CALL gnuplot_ylabel('C_{45} (kbar)',.FALSE.)
      CALL gnuplot_write_file_mul_point(filelastic,1,15,'color_red',.TRUE.,&
                                                             .TRUE.,.FALSE.)
   ENDIF
ENDIF

IF (laue==2) THEN
   CALL gnuplot_ylabel('C_{14} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,12,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{15} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,13,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{16} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,14,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{24} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,15,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{25} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,16,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{26} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,17,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{34} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,18,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{35} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,19,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{45} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,21,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{46} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,22,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
   CALL gnuplot_ylabel('C_{56} (kbar)',.FALSE.)
   CALL gnuplot_write_file_mul_point(filelastic,1,23,'color_red',.TRUE.,&
                                                          .TRUE.,.FALSE.)
ENDIF

CALL gnuplot_end()

IF (ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_elastic_p
