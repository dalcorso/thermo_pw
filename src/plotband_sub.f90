!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plotband_sub(icode,igeom,file_disp)
  !
  ! reads data files produced by "bands_sub", produces
  ! * data file ready for plotting with gnuplot, xmgr or the like
  ! Important notice:
  ! - k-points processed by bands.x should be along a continuous path
  !
  ! icode allows for fine control of the output. Presently the 
  ! implemented options are:
  !
  ! 1   band structure plot
  ! 2   phonon dispersion plot
  ! 3   gruneisen parameters plot         
  ! 
  USE kinds, ONLY : DP
  USE control_thermo, ONLY : filband, flfrq, flgrun
  USE control_bands, ONLY : emin_input, emax_input, flpband
  USE control_grun,  ONLY : flpgrun, grunmin_input, grunmax_input
  USE control_paths, ONLY : letter_path, label_list, label_disp_q, npk_label, &
                            nqaux
  USE thermo_mod,    ONLY : ngeo
  USE constants,     ONLY : rytoev
  USE ener,          ONLY : ef
  USE klist,         ONLY : degauss, nelec
  USE ifc,           ONLY : freqmin_input, freqmax_input
  USE mp,            ONLY : mp_bcast
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_images,     ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: file_disp
  INTEGER, INTENT(IN) :: icode, igeom
  REAL(DP), ALLOCATABLE :: e(:,:), k(:,:), e_in(:), kx(:)
  REAL(DP), ALLOCATABLE :: e_rap(:,:), k_rap(:,:)
  REAL(DP) :: k1(3), k2(3), ps
  REAL(DP) :: emin, emax, eps=1.d-4
  REAL(DP) :: mine, dxmod, dxmod_save, eref, modk1, modk2
  INTEGER, ALLOCATABLE :: nbnd_rapk(:), rap(:,:)
  INTEGER, ALLOCATABLE :: npoints(:)
  INTEGER, ALLOCATABLE :: point(:), nrap(:)
  INTEGER :: nks = 0, nbnd = 0, nlines
  INTEGER :: nks_rap = 0, nbnd_rap = 0
  INTEGER :: ilines, irap, ibnd, ipoint, jnow, ios, i, j, n, ierr
  LOGICAL, ALLOCATABLE :: high_symmetry(:), is_in_range(:), &
                          is_in_range_rap(:), todo(:,:), has_points(:,:)
  LOGICAL :: exist_rap
  CHARACTER(LEN=256) :: filename, filedata, fileout
  CHARACTER(LEN=6), EXTERNAL :: int_to_char

  NAMELIST /plot/ nks, nbnd
  NAMELIST /plot_rap/ nks_rap, nbnd_rap

  IF ( my_image_id /= root_image ) RETURN

  IF (icode==1.OR.icode==2) THEN
    IF (flpband == ' ') RETURN
    fileout=TRIM(flpband)
  ELSEIF (icode==3) THEN
    IF (flpgrun == ' ') RETURN
    fileout=TRIM(flpgrun)
  ENDIF
  IF (icode==1) filedata = TRIM(filband)
  IF (icode==2) filedata = TRIM(flfrq)
  IF (icode==3) filedata = TRIM(flgrun)

  IF (ionode) &
     OPEN(UNIT=1,FILE=TRIM(filedata),FORM='formatted',STATUS='OLD',ERR=10,&
                                                  IOSTAT=ios)
10 CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('plotband_sub','opening band file',ABS(ios))

  IF (ionode) READ (1, plot, IOSTAT=ios)
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('plotband_sub','reading plot namelist',ABS(ios))
  CALL mp_bcast(nks, ionode_id, intra_image_comm)
  CALL mp_bcast(nbnd, ionode_id, intra_image_comm)
  !
  IF (nks <= 0 .or. nbnd <= 0) THEN
     CALL errore('plotband_sub','reading plot namelist',ABS(ios))
  ELSE
     WRITE(stdout, '("Reading ",i4," bands at ",i6," k-points")') nbnd, nks
  ENDIF

  IF (icode==3) THEN
     filename = TRIM(file_disp)//'.g'//TRIM(int_to_char(ngeo/2+1))//".rap"
  ELSE
     filename=TRIM(filedata)//".rap"
  ENDIF
  exist_rap=.TRUE.
  IF (ionode) OPEN(UNIT=21, FILE=TRIM(filename), FORM='formatted', &
                   STATUS='old', ERR=100, IOSTAT=ios)
100 CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) exist_rap=.FALSE.

  IF (exist_rap) THEN
     IF (ionode) READ (21, plot_rap, ERR=110, IOSTAT=ios)
110  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     IF (ios == 0 ) THEN
        CALL mp_bcast(nks_rap, ionode_id, intra_image_comm)
        CALL mp_bcast(nbnd_rap, ionode_id, intra_image_comm)
        IF (nks_rap/=nks.or.nbnd_rap/=nbnd) THEN
           WRITE(stdout,'("file with representations not compatible &
                          & with bands")')
           exist_rap=.FALSE.
        ENDIF
     ELSE
        WRITE(6,'("Problem reading representation file")')
        exist_rap=.FALSE.
     ENDIF
  ENDIF
  !
  ALLOCATE (e(nbnd,nks))
  ALLOCATE (k(3,nks)) 
  ALLOCATE (e_in(nks)) 
  ALLOCATE (kx(nks)) 
  ALLOCATE (npoints(nks)) 
  ALLOCATE (high_symmetry(nks))
  ALLOCATE (is_in_range(nbnd))
  ALLOCATE (point(nks))
  ALLOCATE(nrap(nks))

  IF (exist_rap) THEN
     ALLOCATE(nbnd_rapk(nks))
     ALLOCATE(e_rap(nbnd,nks))
     ALLOCATE(rap(nbnd,nks))
     ALLOCATE(k_rap(3,nks))
     ALLOCATE(todo(nbnd,2))
     ALLOCATE(is_in_range_rap(nbnd))
  ENDIF

  high_symmetry=.FALSE.

  IF (ionode) THEN
     ierr=0
     DO n=1,nks
        READ(1,*,end=220,err=220) ( k(i,n), i=1,3 )
        READ(1,*,end=220,err=220) (e(i,n),i=1,nbnd)
        IF (exist_rap) THEN
           READ(21,*,end=221,err=221) (k_rap(i,n),i=1,3), high_symmetry(n)
           READ(21,*,end=221,err=221) (rap(i,n),i=1,nbnd)
           IF (abs(k(1,n)-k_rap(1,n))+abs(k(2,n)-k_rap(2,n))+  &
               abs(k(3,n)-k_rap(3,n))  > eps ) THEN
               WRITE(stdout,'("Incompatible k points in rap file")')
               exist_rap=.false.
           ENDIF
        ENDIF
     ENDDO
     CLOSE(unit=1)
     IF (exist_rap) CLOSE(unit=21)
     GOTO 222
220   ierr=1
     GOTO 222
221   ierr=2
  ENDIF
222   CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF (ierr==1) THEN
!
!   This is a fatal error, stop the run
!
     WRITE(stdout, '("Error reading k-point # ",i4)') n
     CALL errore('plotband_sub','problem reading data',1)
  ENDIF
  IF (ierr==2) THEN
!
!   This is not a fatal error, continue the run with representations
!
     DEALLOCATE(nbnd_rapk)
     DEALLOCATE(e_rap)
     DEALLOCATE(rap)
     DEALLOCATE(k_rap)
     DEALLOCATE(todo)
     DEALLOCATE(is_in_range_rap)
     CLOSE(UNIT=21,STATUS='KEEP')
     exist_rap=.FALSE.
  ENDIF 
  CALL mp_bcast(k, ionode_id, intra_image_comm)
  CALL mp_bcast(e, ionode_id, intra_image_comm)
  CALL mp_bcast(exist_rap, ionode_id, intra_image_comm)
  IF (exist_rap) THEN
     CALL mp_bcast(k_rap, ionode_id, intra_image_comm)
     CALL mp_bcast(rap, ionode_id, intra_image_comm)
     CALL mp_bcast(high_symmetry, ionode_id, intra_image_comm)
  ENDIF
!
!  Now find the high symmetry points in addition to those already identified
!  in the representation file
!
  DO n=1,nks
     IF (n==1 .OR. n==nks) THEN
        high_symmetry(n) = .true.
     ELSE
        k1(:) = k(:,n) - k(:,n-1)
        k2(:) = k(:,n+1) - k(:,n)
        modk1=sqrt( k1(1)*k1(1) + k1(2)*k1(2) + k1(3)*k1(3) )
        modk2=sqrt( k2(1)*k2(1) + k2(2)*k2(2) + k2(3)*k2(3) )
        IF (modk1 <1.d-6 .OR. modk2 < 1.d-6) CYCLE
        ps = ( k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3) ) / &
             modk1 / modk2 
        high_symmetry(n) = (ABS(ps-1.d0) >1.0d-4).OR.high_symmetry(n)
!
!  The gamma point is a high symmetry point
!
        IF (k(1,n)**2+k(2,n)**2+k(3,n)**2 < 1.0d-9) high_symmetry(n)=.true.
!
!   save the typical lenght of dk
!
        IF (n==2) dxmod_save = sqrt( k1(1)**2 + k1(2)**2 + k1(3)**2)

     ENDIF
  ENDDO

  kx(1) = 0.d0
  DO n=2,nks
     dxmod=sqrt ( (k(1,n)-k(1,n-1))**2 + &
                  (k(2,n)-k(2,n-1))**2 + &
                  (k(3,n)-k(3,n-1))**2 )
     IF (dxmod > 5*dxmod_save) THEN
!
!   A big jump in dxmod is a sign that the point k(:,n) and k(:,n-1)
!   are quite distant and belong to two different lines. We put them on
!   the same point in the graph 
!
        kx(n)=kx(n-1)
     ELSEIF (dxmod > 1.d-5) THEN
!
!  This is the usual case. The two points k(:,n) and k(:,n-1) are in the
!  same path.
!
        kx(n) = kx(n-1) + dxmod
        dxmod_save = dxmod
     ELSE
!
!  This is the case in which dxmod is almost zero. The two points coincide
!  in the graph, but we do not save dxmod.
!
        kx(n) = kx(n-1) +  dxmod

     ENDIF
  ENDDO

  emin=1.d10
  emax=-1.d10
  DO n=1,nks
     DO i=1,nbnd
        emin = min(emin, e(i,n))
        emax = max(emax, e(i,n))
     ENDDO
  ENDDO

  IF (icode==1) THEN
!
!   bands can be shifted
!
     eref=-1d20
     IF (degauss > 0.0_DP) THEN
        eref=ef * rytoev
     ELSE
        DO ibnd=1, NINT(nelec/2)
           IF (e(ibnd,1) > eref) eref=e(ibnd,1)
        END DO
     END IF
     IF (emin_input /= 0.0_DP) emin=emin_input + eref
     IF (emax_input /= 0.0_DP) emax=emax_input + eref
  ELSE IF (icode==2) THEN
!
!   no shift for phonon, but take the maximum energy slightly above the 
!   dispersion and an integer value
!
     eref=0.0_DP
     emax=NINT(emax*1.05_DP)
     emin=NINT(emin)
     IF (freqmin_input /= 0.0_DP ) emin = freqmin_input
     IF (freqmax_input /= 0.0_DP ) emax = freqmax_input
  ELSEIF (icode==3) THEN
!
!  no shift for gruneisen parameters, and slightly enlarge the maximum 
!  and minimum values, or set the input values.
!
     eref=0.0_DP
     emax=emax*1.05_DP
     IF (emin < 0.0_DP) THEN
        emin=emin*1.05_DP
     ELSE
        emin=emin*0.95_DP
     ENDIF
     IF (grunmin_input /= 0.0_DP ) emin = grunmin_input
     IF (grunmax_input /= 0.0_DP ) emax = grunmax_input
  ELSE
     CALL errore('plotband_sub','Problem with icode',1)
  ENDIF
!
!  Since the minimum and maximum energies are given in input we can
!  sign the bands that are completely outside this range.
!
  is_in_range = .false.
  DO i=1,nbnd
     is_in_range(i) = any (e(i,1:nks) >= emin .and. e(i,1:nks) <= emax)
  ENDDO
!
!  Now we compute how many paths there are: nlines
!  The first point of this path: point(iline)
!  How many points are in each path: npoints(iline)
!
  DO n=1,nks
     IF (high_symmetry(n)) THEN
        IF (n==1) THEN
!
!   first point. Initialize the number of lines, and the number of point
!   and say that this line start at the first point
!
           nlines=1
           npoints(1)=1
           point(1)=1
        ELSEIF (n==nks) THEN
!
!    Last point. Here we save the last point of this line, but
!    do not increase the number of lines
!
           npoints(nlines) = npoints(nlines)+1
           point(nlines+1)=n
        ELSE
!
!   Middle line. The current line has one more point, and there is a new
!   line that has to be initialized. It has one point and its first point
!   is the current k.
!
           npoints(nlines) = npoints(nlines)+1
           nlines=nlines+1
           npoints(nlines) = 1
           point(nlines)=n
        ENDIF
!        IF (n==1) THEN
!           WRITE( stdout,'("high-symmetry point: ",3f7.4,&
!                         &"   x coordinate   0.0000")') (k(i,n),i=1,3)
!        ELSE
!           WRITE( stdout,'("high-symmetry point: ",3f7.4,&
!                         &"   x coordinate",f9.4)') (k(i,n),i=1,3), kx(n)
!        ENDIF
     ELSE
!
!   This k is not an high symmetry line so we just increase the number of
!   points of this line.
!
        npoints(nlines) = npoints(nlines)+1
     ENDIF
  ENDDO
  !
  ALLOCATE(has_points(nlines,12))
  has_points=.FALSE.
  IF (.NOT.exist_rap) THEN
!
!  Here the symmetry analysis has not been done. So simply save the bands
!  on output. The odd one from left to right, the even one from right to
!  left.
!
     IF (ionode) OPEN (UNIT=2, FILE=TRIM(fileout), FORM='formatted',&
                 STATUS='unknown', IOSTAT=ios)
     CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('plotband_sub','opening file'//TRIM(flpband),ABS(ios))
     ! draw bands
     DO i=1,nbnd
        IF (is_in_range(i)) THEN
           IF ( mod(i,2) /= 0) THEN
              IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e(i,n),n=1,nks)
           ELSE
              IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e(i,n),n=nks,1,-1)
           ENDIF
        ENDIF
     ENDDO
     IF (ionode) CLOSE (UNIT = 2, STATUS='KEEP')
  ELSE
!
!   In this case we write a diffent file for each line and for each
!   representation. Each file contains the bands of that representation.
!   The file is called filename.#line.#rap
!
!   First determine for each line how many representations are there
!   in each line
!
     DO ilines=1,nlines
        nrap(ilines)=0
        DO ipoint=1,npoints(ilines)-2
           n=point(ilines) + ipoint
           DO ibnd=1,nbnd
              nrap(ilines)=max(nrap(ilines),rap(ibnd,n))
           ENDDO
        ENDDO
        IF (nrap(ilines) > 12) CALL errore("plotband_sub",&
                                           "Too many representations",1)
     ENDDO
!
!   Then, for each line and for each representation along that line
!
     DO ilines=1,nlines
        IF (nrap(ilines)==0) THEN
!
!   Along this line the symmetry decomposition has not been done.
!   Plot all the bands as in the standard case
!
           filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines))

           IF (ionode) OPEN(UNIT=2, FILE=TRIM(filename), FORM='formatted',&
               STATUS='unknown', IOSTAT=ios)
           CALL mp_bcast(ios, ionode_id, intra_image_comm)
           ! draw bands
           DO i=1,nbnd
              IF (is_in_range(i)) THEN
                 IF ( mod(i,2) /= 0) THEN
                    IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e(i,n), &
                                            n=point(ilines), point(ilines+1))
                 ELSE
                    IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e(i,n), &
                                   n=point(ilines+1), point(ilines),-1 )
                 ENDIF
                 has_points(ilines,1)=.TRUE.
              ENDIF
           ENDDO
           IF (ionode) CLOSE (unit = 2)
        ENDIF

        todo=.true.
        DO irap=1, nrap(ilines)
!
!     open a file
!
           filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines)) &
                                   //  "." // TRIM(int_to_char(irap))
           IF (ionode) OPEN (UNIT=2, FILE=TRIM(filename), FORM='formatted', &
                       STATUS='unknown', IOSTAT=ios)
           CALL mp_bcast(ios, ionode_id, intra_image_comm)
           CALL errore("plotband_sub","opening file"//TRIM(filename), &
                                                       ABS(ios)) 
!  For each k point along this line selects only the bands which belong
!  to the irap representation
           nbnd_rapk=100000
           DO n=point(ilines)+1, point(ilines+1)-1
              nbnd_rapk(n)=0
              DO i=1,nbnd
                 IF (rap(i,n)==irap) THEN
                    nbnd_rapk(n) = nbnd_rapk(n) + 1
                    e_rap(nbnd_rapk(n),n)=e(i,n)
                 ENDIF
              ENDDO
           ENDDO
!
!   on the two high symmetry points the representation is different. So for each
!   band choose the closest eigenvalue available.
!
           DO i=1,nbnd_rapk(point(ilines)+1)
              mine=1.e8
              DO j=1,nbnd
                 IF (abs(e_rap(i,point(ilines)+1)-e(j,point(ilines)))<mine &
                                                        .and. todo(j,1)) THEN
                    e_rap(i,point(ilines))=e(j,point(ilines))
                    mine=abs( e_rap(i,point(ilines)+1)-e(j,point(ilines)))
                    jnow=j
                 ENDIF
              ENDDO
              todo(jnow,1)=.false.
           ENDDO
           DO i=1,nbnd_rapk(point(ilines+1)-1)
              mine=1.e8
              DO j=1,nbnd
                 IF (abs(e_rap(i,point(ilines+1)-1)- &
                          e(j,point(ilines+1)))<mine .and. todo(j,2)) THEN
                    e_rap(i,point(ilines+1))=e(j,point(ilines+1))
                    mine=abs(e_rap(i,point(ilines+1)-1)-e(j,point(ilines+1)) )
                    jnow=j
                 ENDIF
              ENDDO
              todo(jnow,2)=.false.
           ENDDO
           is_in_range_rap=.false.
           DO i=1,minval(nbnd_rapk)
              is_in_range_rap(i) = any (e_rap(i,point(ilines):point(ilines+1))&
                    >= emin .and. e(i,point(ilines):point(ilines+1)) <= emax)
           ENDDO
           DO i=1,minval(nbnd_rapk)
              IF (is_in_range_rap(i)) THEN
                 IF ( mod(i,2) /= 0) THEN
                    IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e_rap(i,n), &
                                        n=point(ilines),point(ilines+1))
                 ELSE
                    IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e_rap(i,n), &
                                       n=point(ilines+1),point(ilines),-1)
                 ENDIF
                 has_points(ilines,irap)=.TRUE.
              ENDIF
           ENDDO
           IF (minval(nbnd_rapk)==0) THEN
              IF (ionode) CLOSE (UNIT = 2, STATUS='DELETE')
           ELSE
              IF (ionode) CLOSE (UNIT = 2, STATUS='KEEP')
           ENDIF
        ENDDO
     ENDDO
  ENDIF

  WRITE(stdout,'("bands in xmgr format written to file ",a)') flpband
  !
  IF (emin_input /= 0.0_DP) THEN
     emin=emin_input
  ELSE
     emin=emin-eref
  END IF

  IF (emax_input /= 0.0_DP) THEN
     emax=emax_input 
  ELSE
     emax=emax-eref
  END IF

  CALL write_gnuplot_file(kx, e, nks, nbnd, emin, emax, eref, nlines, nrap, &
                  has_points, point, icode, exist_rap, igeom, fileout)


  DEALLOCATE(point)
  DEALLOCATE(is_in_range)
  DEALLOCATE(high_symmetry)
  DEALLOCATE(npoints) 
  DEALLOCATE(kx) 
  DEALLOCATE(e_in) 
  DEALLOCATE(k) 
  DEALLOCATE(e)
  DEALLOCATE(has_points)
  DEALLOCATE(nrap)

  IF (exist_rap) THEN
     DEALLOCATE(is_in_range_rap)
     DEALLOCATE(todo)
     DEALLOCATE(k_rap)
     DEALLOCATE(rap)
     DEALLOCATE(e_rap)
     DEALLOCATE(nbnd_rapk)
  ENDIF

  RETURN
END SUBROUTINE plotband_sub

SUBROUTINE write_gnuplot_file(kx, e, nks, nbnd, emin, emax, eref, nlines, &
                  nrap, has_points, point, icode, exist_rap, igeom, fileout)
USE kinds,           ONLY : DP
USE klist,           ONLY : degauss
USE control_paths,   ONLY : nqaux, label_disp_q, letter_path, &
                            q_in_band_form
USE control_gnuplot, ONLY : flgnuplot, flpsband, flpsdisp, &
                            flpsgrun, gnuplot_command, lgnuplot
USE gnuplot,       ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                          gnuplot_write_file_data, gnuplot_ylabel, &
                          gnuplot_write_vertical_line, gnuplot_write_label, &
                          gnuplot_write_horizontal_line, &
                          gnuplot_write_label_yl, gnuplot_write_command, &
                          gnuplot_set_eref, gnuplot_unset_xticks
USE io_global,     ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: fileout
INTEGER, INTENT(IN) :: nks, nbnd, nlines, icode, igeom
INTEGER, INTENT(IN) :: point(nks), nrap(nks)
REAL(DP), INTENT(IN) :: kx(nks), e(nbnd,nks)
REAL(DP), INTENT(IN) :: emin, emax, eref 
LOGICAL, INTENT(IN) :: has_points(nlines,12)
LOGICAL, INTENT(IN) :: exist_rap

INTEGER :: ilines, irap, n, ierr
INTEGER :: system
REAL(DP) :: shift
CHARACTER(LEN=256) :: gnu_filename, filename, command
CHARACTER(LEN=30) :: colore(12)
CHARACTER(LEN=6), EXTERNAL :: int_to_char

IF (icode==1) THEN
   gnu_filename=TRIM(flgnuplot)//'_band'
ELSEIF (icode==2) THEN
   gnu_filename=TRIM(flgnuplot)//'_disp'
ELSEIF (icode==3) THEN
   gnu_filename=TRIM(flgnuplot)//'_grun'
ENDIF

CALL gnuplot_start(gnu_filename)

IF (icode==1) THEN
   filename=TRIM(flpsband)
ELSEIF (icode==2) THEN
   filename=TRIM(flpsdisp)
ELSEIF (icode==3) THEN
   filename=TRIM(flpsgrun)
ENDIF
IF (igeom > 1) filename=filename//TRIM(int_to_char(igeom))

CALL gnuplot_write_header(filename, kx(1), kx(nks), emin, emax ) 
CALL gnuplot_unset_xticks(.FALSE.) 
IF (icode==1) THEN
   CALL gnuplot_ylabel('Energy (eV)',.FALSE.) 
   IF (degauss > 0.0_DP) CALL gnuplot_write_horizontal_line(0.0_DP, 2, &
                                         'front', 'black', .FALSE.)
ELSEIF (icode==2) THEN
   CALL gnuplot_ylabel('Frequency (cm^{-1})',.FALSE.) 
ELSEIF (icode==3) THEN
   CALL gnuplot_ylabel('{/Symbol g}_{/Symbol n}({/Helvetica-Bold q})',.FALSE.) 
ENDIF
DO ilines = 2, nlines
   CALL gnuplot_write_vertical_line(kx(point(ilines)), 2, 'front', 'black', &
                                     .FALSE.)
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
                                       'black', .FALSE.)
   CALL gnuplot_write_label_yl(kx(label_disp_q(n)), ' ymin + shift ', letter_path(n),.FALSE.)
   ENDDO
ELSE
   DO n=1, nqaux
      IF (letter_path(n) /='') &
         CALL gnuplot_write_label_yl(kx(label_disp_q(n)), &
                          ' ymin + shift ', letter_path(n),.FALSE.)
   END DO
END IF

colore(1)='red'
colore(2)='green'
colore(3)='blue'
colore(4)='cyan'
colore(5)='magenta'
colore(6)='yellow'
colore(7)='pink'
colore(8)='black'
colore(10)='grey'
colore(11)='light-blue'
colore(12)='orange'

IF (.NOT.exist_rap) THEN
   filename=TRIM(fileout)
   CALL gnuplot_write_file_data(filename,'red',.TRUE.,.TRUE.,.FALSE.)
ELSE
   DO ilines = 1, nlines
      IF (nrap(ilines)==0) THEN
         filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines))
         IF (has_points(ilines,1)) THEN
            IF (ilines==1) THEN
               CALL gnuplot_write_file_data(filename,'red',.TRUE.,.FALSE.,&
                                                     .FALSE.)
            ELSEIF (ilines==nlines) THEN
               CALL gnuplot_write_file_data(filename,'red',.FALSE.,.TRUE.,&
                                                     .FALSE.)
            ELSE
               CALL gnuplot_write_file_data(filename,'red',.FALSE.,.FALSE.,&
                                                     .FALSE.)
            END IF
         END IF
      ELSE
         DO irap=1, nrap(ilines)
            filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines)) &
                                 //  "." // TRIM(int_to_char(irap))
            IF (has_points(ilines,irap)) THEN
               IF (ilines==1.AND.irap==1) THEN
                  CALL gnuplot_write_file_data(filename,colore(irap),&
                                              .TRUE.,.FALSE.,.FALSE.)
               ELSEIF (ilines==nlines.AND.irap==nrap(nlines)) THEN
                  CALL gnuplot_write_file_data(filename,colore(irap),&
                                              .FALSE.,.TRUE.,.FALSE.)
               ELSE
                  CALL gnuplot_write_file_data(filename,colore(irap),&
                                              .FALSE.,.FALSE.,.FALSE.)
               ENDIF
            ENDIF
         ENDDO
      ENDIF  
   ENDDO
ENDIF

CALL gnuplot_end()

IF (lgnuplot.AND.ionode) &
   ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

RETURN
END SUBROUTINE write_gnuplot_file
