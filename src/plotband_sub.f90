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
  ! 4   phonon dispersion plot of interpolated frequencies
  !
  ! 
  USE kinds, ONLY : DP
  USE control_bands, ONLY : emin_input, emax_input 
  USE control_grun,  ONLY : grunmin_input, grunmax_input
  USE data_files,    ONLY : flpgrun, flpband, filband, flfrq, flgrun
  USE control_paths, ONLY : letter_path, label_list, label_disp_q, npk_label, &
                            nqaux, nrap_plot, rap_plot
  USE control_2d_bands, ONLY : nkz, aux_ind_sur, identify_sur, lprojpbs, &
                               sym_divide, lsurface_state, lsurface_state_rap
  USE thermo_mod,    ONLY : tot_ngeo
  USE thermo_sym,    ONLY : code_group_save
  USE constants,     ONLY : rytoev
  USE point_group,   ONLY : convert_rap, has_sigma_h
  USE ions_base,     ONLY : nat
  USE ener,          ONLY : ef
  USE klist,         ONLY : degauss, nelec
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,      ONLY : lspinorb
  USE ifc,           ONLY : freqmin_input, freqmax_input
  USE mp,            ONLY : mp_bcast
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_images,     ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: file_disp
  INTEGER, INTENT(IN) :: icode, igeom
  REAL(DP), ALLOCATABLE :: e(:,:), k(:,:), kx(:)
  REAL(DP), ALLOCATABLE :: e_rap(:,:), k_rap(:,:), k_eff(:,:), e_eff(:,:)
  REAL(DP) :: k1(3), k2(3), ps
  REAL(DP) :: emin, emax, eps=1.d-4
  REAL(DP) :: dxmod, dxmod_save, eref, modk1, modk2
  INTEGER, ALLOCATABLE :: nbnd_rapk(:,:), rap(:,:), rap_eff(:,:), gcodek(:), &
                          aux_ind(:), gcodek_eff(:), aux_ind_eff(:)
  INTEGER, ALLOCATABLE :: start_rapk(:,:)
  INTEGER, ALLOCATABLE :: start_point(:), last_point(:), nrap(:)
  INTEGER, ALLOCATABLE :: start_point_eff(:), last_point_eff(:), rapin(:)
  INTEGER, ALLOCATABLE :: nrap_plot_eff(:), rap_plot_eff(:,:), nbnd_count(:)
  INTEGER, ALLOCATABLE :: label_disp_q_eff(:)
  INTEGER :: nks = 0, nbnd = 0, nlines, nks_
  INTEGER :: code_group_line, nband_occ
  INTEGER :: nks_rap = 0, nbnd_rap = 0
  INTEGER :: ilines, irap, ibnd, ipoint, jnow, ios, i, j, n, ik, ikz, &
             ike, ik2, spe, lpe, nbc, iq, tot_points, ishift, ierr
  INTEGER :: start_shift, last_shift, central_geo
  LOGICAL, ALLOCATABLE :: high_symmetry(:), is_in_range(:), &
                          is_in_range_rap(:), has_points(:,:), &
                          lsurface_state_eff(:,:)
  LOGICAL :: exist_rap, type1
  CHARACTER(LEN=256) :: filename, filedata, fileout
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  CHARACTER(LEN=11) :: group_name

  NAMELIST /plot/ nks, nbnd
  NAMELIST /plot_rap/ nks_rap, nbnd_rap

  IF ( my_image_id /= root_image ) RETURN


  IF (icode==1.OR.icode==2) THEN
    IF (flpband == ' ') RETURN
    fileout=TRIM(flpband)
  ELSEIF (icode==3) THEN
    IF (flpgrun == ' ') RETURN
    fileout=TRIM(flpgrun)
  ELSEIF (icode==4) THEN
    IF (flpgrun == ' ') RETURN
    fileout=TRIM(flpgrun)//'_freq'
  ENDIF
  IF (icode==1) filedata = TRIM(filband)
  IF (icode==2) filedata = TRIM(flfrq)
  IF (icode==3) filedata = TRIM(flgrun)
  IF (icode==4) filedata = TRIM(flgrun)//'_freq'


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
     WRITE(stdout, '(/,5x,"Reading ",i4," bands at ",i6," k-points")') nbnd, nks
  ENDIF

  IF (icode==3.OR.icode==4) THEN
     central_geo=tot_ngeo/2
     IF (MOD(tot_ngeo,2)==1) central_geo=central_geo+1
     filename = TRIM(file_disp)//'.g'//TRIM(int_to_char(central_geo))//".rap"
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
        WRITE(stdout,'("Problem reading representation file")')
        exist_rap=.FALSE.
     ENDIF
  ENDIF
  !
  ALLOCATE (e(nbnd,nks))
  ALLOCATE (k(3,nks)) 
  ALLOCATE (high_symmetry(nks))
  ALLOCATE (rap(nbnd,nks))
  ALLOCATE (gcodek(nks))
  ALLOCATE (aux_ind(nks))
  IF (exist_rap) THEN
     ALLOCATE(k_rap(3,nks))
  ENDIF

  high_symmetry=.FALSE.
  rap=-1

  IF (ionode) THEN
     ierr=0
     DO n=1,nks
        READ(1,*,end=220,err=220) (k(i,n), i=1,3 )
        READ(1,*,end=220,err=220) (e(i,n),i=1,nbnd)
        IF (exist_rap) THEN
           READ(21,*,end=221,err=221) (k_rap(i,n),i=1,3), high_symmetry(n), &
                                      gcodek(n), aux_ind(n)
           READ(21,*,end=221,err=221) (rap(i,n),i=1,nbnd)
           IF (abs(k(1,n)-k_rap(1,n))+abs(k(2,n)-k_rap(2,n))+  &
               abs(k(3,n)-k_rap(3,n))  > eps .AND. icode /=3 &
                                                  .AND. icode /=4 ) THEN
               WRITE(stdout,'("Incompatible k points in rap file")')
               CLOSE(unit=21)
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
!   This is not a fatal error, continue the run without representations
!
     DEALLOCATE(k_rap)
     CLOSE(UNIT=21,STATUS='KEEP')
     exist_rap=.FALSE.
  ENDIF 
  CALL mp_bcast(k, ionode_id, intra_image_comm)
  CALL mp_bcast(e, ionode_id, intra_image_comm)
  CALL mp_bcast(exist_rap, ionode_id, intra_image_comm)
  IF (exist_rap) THEN
     CALL mp_bcast(k_rap, ionode_id, intra_image_comm)
     CALL mp_bcast(rap, ionode_id, intra_image_comm)
     CALL mp_bcast(gcodek, ionode_id, intra_image_comm)
     CALL mp_bcast(aux_ind, ionode_id, intra_image_comm)
     CALL mp_bcast(high_symmetry, ionode_id, intra_image_comm)
  ENDIF
  IF (identify_sur) CALL identify_surface_states(nat,nbnd,nks,e,rap)
  WRITE(stdout,'(5x,"Starting the generation of the plot",/)') 
!
!  at this point we have:
!  k, e   for all k and bands
!  If the representations exist
!  rap    for all k and bands
!  gcodek for all k
!  aux_ind for all k
!  If this is a surface band calculation we have also
!  lsurface_states for all k and band tell us if this is a surface state
!

!
!   In a projected band structure calculation, we change the representations
!   of the points groups of higher symmetry that might occur for particular 
!   values of k_z and bring the representations to those of the subgroup
!   of the surface. It is assumed that the second k_z if nzk is odd, or the 
!   first k_z if nkz is even has no more symmetry than the surface
!     
  ALLOCATE (rapin(nbnd))
  nks_= nks / nkz
  IF (exist_rap) THEN
      IF (nkz > 1 .AND. lprojpbs .AND. sym_divide) THEN
         type1=has_sigma_h(code_group_save)
         IF (type1) THEN
            ishift=nks_
         ELSE
            ishift=0
         ENDIF
         DO ikz=1,nkz
            DO ik = 1, nks_
               ike = ik + nks_ * ( ikz - 1 )
               ik2 = ik + ishift
               IF (gcodek(ike) /= gcodek(ik2)) THEN
                  rapin(:)=rap(:,ike)
                  CALL convert_rap(nbnd,rapin,rap(1,ike),&
                       gcodek(ike), gcodek(ik2), aux_ind_sur(ik,ikz),lspinorb)
                  gcodek(ike)=gcodek(ik2)
                  aux_ind(ike) = aux_ind(ik2)
!
!   a point must be high symmetry in all planes.
!
                  high_symmetry(ike)=high_symmetry(ik2)
               ENDIF
            END DO
         END DO
      END IF
  END IF
!
!  Now find the high symmetry points in addition to those already identified
!  in the representation file. This is necessary because the representation
!  file might be missing. In this case we identify here the high
!  symmetry points
!
  DO n=1,nks
     IF (n==1 .OR. n==nks) THEN
        high_symmetry(n) = .TRUE.
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
!
! At this point we should count how many lines do we have and
! set the first and the last point of each line
!
  ALLOCATE (start_point(nks))
  ALLOCATE (last_point(nks))
  DO n=1,nks
     IF (high_symmetry(n)) THEN
        IF (n==1) THEN
!
!   first point. Initialize the number of lines
!   and say that this line start at the first point
!
           nlines=1
           start_point(1)=1
        ELSEIF (n==nks) THEN
!
!    Last point. Here we save the last point of this line, but
!    do not increase the number of lines
!
           last_point(nlines)=n
        ELSEIF (MOD(n,nks_)==0) THEN
!
!   In this case do nothing. This is the last point of a graph, but we have
!   not to increase the number of lines, this is done by the next point
!   which is the first point of the graph
!
           last_point(nlines)=n
        ELSE
!
!   Middle line. In general this is both the last point of the current 
!   line and the start point of the next line. This is not the case 
!   when the previous point was also a high symmetry point
!
           IF (.NOT. high_symmetry(n-1)) THEN
              last_point(nlines)=n
           ENDIF
           IF (.NOT. high_symmetry(n+1)) THEN
              nlines=nlines+1
              start_point(nlines)=n
           ENDIF
        ENDIF
     ENDIF
  ENDDO
  tot_points=0
  DO ilines=1,nlines
     tot_points=tot_points+last_point(ilines)-start_point(ilines)+1
!     write(stdout,*) 'ilines, start_point, last_point', ilines,  &
!                 start_point(ilines), last_point(ilines)
  ENDDO
  WRITE(stdout,'(5x,"Number of lines:", i4, " Total number of points:",i8)') &
                                        nlines, tot_points

  ALLOCATE (kx(tot_points)) 
  ALLOCATE (k_eff(3,tot_points)) 
  ALLOCATE (e_eff(nbnd,tot_points)) 
  ALLOCATE (gcodek_eff(tot_points))
  ALLOCATE (rap_eff(nbnd,tot_points))
  ALLOCATE (aux_ind_eff(tot_points))
  ALLOCATE (start_point_eff(nlines))
  ALLOCATE (last_point_eff(nlines))
  ALLOCATE (nrap_plot_eff(tot_points))
  ALLOCATE (rap_plot_eff(12,tot_points))
  ALLOCATE (nbnd_count(tot_points))
  ALLOCATE (is_in_range(nbnd))
  ALLOCATE (nrap(nlines))
  ALLOCATE (nbnd_rapk(12,tot_points))
  ALLOCATE (start_rapk(12,tot_points))
  ALLOCATE (has_points(nlines,12))
  ALLOCATE (label_disp_q_eff(nqaux))
  ALLOCATE (e_rap(nbnd,tot_points))
  IF (identify_sur) THEN
     ALLOCATE(lsurface_state_eff(nbnd,tot_points))
     ALLOCATE(lsurface_state_rap(nbnd,tot_points))
  ENDIF

  IF (exist_rap) THEN
     ALLOCATE(is_in_range_rap(nbnd))
  ENDIF
!
!   Create a new list of points and eigenvalue so that there is no
!   more overlapping points that belong to two lines
!
  ik=0
  DO ilines=1,nlines
     DO n=start_point(ilines), last_point(ilines)
        ik=ik+1
        IF (ik>tot_points) CALL errore('plotband_sub','problem with points',1)
        k_eff(:,ik)=k(:,n)
        e_eff(:,ik)=e(:,n)
        rap_eff(:,ik)=rap(:,n)
        aux_ind_eff(ik)=aux_ind(n)
        gcodek_eff(ik)=gcodek(n)
        IF (identify_sur) lsurface_state_eff(:,ik)=lsurface_state(:,n)
        nrap_plot_eff(ik)=nrap_plot(n)
        rap_plot_eff(:,ik)=rap_plot(:,n)
!
!  nrap_plot and rap_plot are wrong on the first and last point that might be
!  in common to more than one line. The effective value are set acconding to the
!  values in the same line and are correct. Set also the effective starting and
!  last point
!
        IF (n==start_point(ilines)) THEN
           start_point_eff(ilines)=ik
           nrap_plot_eff(ik)=nrap_plot(n+1)
           rap_plot_eff(:,ik)=rap_plot(:,n+1)
        ENDIF
        IF (n==last_point(ilines)) THEN
           last_point_eff(ilines)=ik
           nrap_plot_eff(ik)=nrap_plot(n-1)
           rap_plot_eff(:,ik)=rap_plot(:,n-1)
        ENDIF
        DO iq=1, nqaux
           IF (label_disp_q(iq)==n) label_disp_q_eff(iq)=ik
        ENDDO
     END DO
  END DO
  IF (ik /= tot_points) CALL errore('plotband_sub','Missing points',1)
!
!  now compute the x coordinate on the plot, but use the effective k
!  points
!
  kx(1) = 0.d0
  nks_=tot_points / nkz
  DO n=2,tot_points
     dxmod=sqrt ( (k_eff(1,n)-k_eff(1,n-1))**2 + &
                  (k_eff(2,n)-k_eff(2,n-1))**2 + &
                  (k_eff(3,n)-k_eff(3,n-1))**2 )
     IF (MOD(n-1,nks_)==0) THEN
!
!   This is a PBS calculation and now we switch to the next plot
!
         kx(n)=0.0_DP
     ELSEIF (dxmod > 5*dxmod_save) THEN
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

  IF (nkz>1) WRITE(stdout,'(/,5x,"Projected band structure calculation, &
                      &nkz=",i5,/)') nkz
  DO ilines=1,nlines
     WRITE(stdout,'(5x,"Line ",i5," starts at",f13.6," ends at ",f13.6)') &
           ilines, kx(start_point_eff(ilines)), kx(last_point_eff(ilines))
!     WRITE(stdout,'(5x,"Line ",i5," starts at",i5," ends at ",i5)') &
!           ilines, start_point_eff(ilines), last_point_eff(ilines)
  ENDDO
!
!  shift the bands if necessary
!
  emin=1.d10
  emax=-1.d10
  DO n=1, tot_points
     DO i=1,nbnd
        emin = min(emin, e_eff(i,n))
        emax = max(emax, e_eff(i,n))
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
        nband_occ=NINT(nelec/2)
        IF (noncolin) nband_occ=nband_occ*2        
        DO ibnd=1, nband_occ
           IF (e_eff(ibnd,1) > eref) eref=e_eff(ibnd,1)
        END DO
     END IF
     IF (emin_input /= 0.0_DP) emin=emin_input + eref
     IF (emax_input /= 0.0_DP) emax=emax_input + eref
  ELSE IF (icode==2.OR.icode==4) THEN
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
     is_in_range(i) = ANY ( e_eff(i,1:tot_points) >= emin .AND. &
                            e_eff(i,1:tot_points) <= emax )
  ENDDO
  !
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
              IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e_eff(i,n), &
                                                               n=1,tot_points)
           ELSE
              IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e_eff(i,n),&
                                                              n=tot_points,1,-1)
           ENDIF
        ENDIF
     ENDDO
     IF (ionode) CLOSE (UNIT = 2, STATUS='KEEP')
     DO n=1,tot_points
        e_rap(:,n)=e_eff(:,n)
        nbnd_rapk(1,n)=nbnd
        start_rapk(1,n)=1
        IF (identify_sur) lsurface_state_rap(:,n)=lsurface_state_eff(:,n)
     ENDDO
  ELSE
!
!   In this case we write a different file for each line and for each
!   representation. Each file contains the bands of that representation.
!   The file is called filename.#line.#rap
!
!
!   Check that all representations along a line refer to the same point group. 
!   Usually  this is not the case in the first and last point of a line.
!   In this case use the group-subgroup compatibility tables
!   to transform the representation of the higher symmetry group 
!   into the low symmetry group.
!

    nbnd_rapk=10000
    nbnd_count=0
    start_rapk=0
    WRITE(stdout,'(/,5x,"Representations per line:")')
    DO ilines=1,nlines
       spe=start_point_eff(ilines)
       lpe=last_point_eff(ilines) 
       code_group_line=gcodek_eff(spe+1)
       DO n=spe,lpe
          IF (gcodek_eff(n) /= code_group_line) THEN
!
!   it is the point before or after the one that change symmetry 
!   that carry the information about the group subgroup relationship
!   (aux_ind)
!
             IF (n==spe) THEN
                rapin(:)=rap_eff(:,n)
                CALL convert_rap(nbnd,rapin,rap_eff(1,n),gcodek_eff(n), &
                                  code_group_line, aux_ind_eff(n+1),lspinorb)
             ELSEIF (n==lpe) THEN
                rapin(:)=rap_eff(:,n)
                CALL convert_rap(nbnd,rapin,rap_eff(1,n),gcodek_eff(n), &
                                  code_group_line, aux_ind_eff(n-1),lspinorb)
             ELSE
                CALL errore('plotband_sub','unexpected change of symmetry',1)
             ENDIF
          END IF
       END DO
!
!  check if the first and last point have representations or not
!
       start_shift=0 
       last_shift=0 
       DO i=1,nbnd
          IF (rap_eff(i,spe) <= 0) start_shift=1
          IF (rap_eff(i,lpe) <= 0) last_shift=1
       ENDDO
!
!   If the point at line border has not the representations 
!   we copy those of the closest point. This is equivalent to make all 
!   avoided crossing between the line border and the point after of before.
!
       IF (start_shift==1) rap_eff(:,spe)=rap_eff(:,spe+1)
       IF (last_shift==1) rap_eff(:,lpe)=rap_eff(:,lpe-1)
!
!   Determine for each line which is the representation with the highest
!   number 
!
       nrap(ilines)=1
       DO n=spe, lpe
          DO ibnd=1,nbnd
             nrap(ilines)=MAX(nrap(ilines),rap_eff(ibnd,n))
          ENDDO
       ENDDO
       IF (nrap(ilines) > 12) CALL errore("plotband_sub",&
                                           "Too many representations",1)
       WRITE(stdout,'(5x,"line ",i7, " nrap",i7,a11)') ilines, nrap(ilines), &
                                         TRIM(group_name(code_group_line))

       IF (nrap(ilines)==1) THEN
!
!   Along this line the symmetry decomposition has not been done or
!   the group has only one representation.
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
                   IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e_eff(i,n), & 
                                                            n=spe, lpe)
                ELSE
                   IF (ionode) WRITE (2,'(2f10.4)') (kx(n), e_eff(i,n), &
                                                            n=lpe, spe,-1 )
                ENDIF
                has_points(ilines,1)=.TRUE.
             ENDIF
          ENDDO
          IF (ionode) CLOSE (unit = 2)
          DO n=spe,lpe
             e_rap(:,n)=e_eff(:,n)
             nbnd_rapk(1,n)=nbnd
             start_rapk(1,n)=1
             IF (identify_sur) lsurface_state_rap(:,n)=lsurface_state_eff(:,n)
          ENDDO
       ELSE
          DO irap=1, nrap(ilines)
!
!     along this line there are several representations
!     open a file for each representation
!
             filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines)) &
                                   //  "." // TRIM(int_to_char(irap))
             IF (ionode) OPEN (UNIT=2, FILE=TRIM(filename), FORM='formatted', &
                       STATUS='unknown', IOSTAT=ios)
             CALL mp_bcast(ios, ionode_id, intra_image_comm)
             CALL errore("plotband_sub","opening file"//TRIM(filename), &
                                                       ABS(ios)) 
!
!  For each k point along this line selects only the bands which belong
!  to the irap representation
!
             DO n = spe, lpe 
                nbnd_rapk(irap, n)=0
                DO i=1,nbnd
                   IF (rap_eff(i,n)==irap) THEN
                      nbnd_rapk(irap,n) = nbnd_rapk(irap,n) + 1
                      nbnd_count(n)=nbnd_count(n)+1
                      nbc=nbnd_count(n)
                      IF (nbnd_rapk(irap,n)==1) THEN
                         start_rapk(irap,n)=nbc
                      END IF
                      e_rap(nbc,n)=e_eff(i,n)
                      IF (identify_sur) &
                         lsurface_state_rap(nbc,n)=lsurface_state_eff(i,n)
                   ENDIF
                ENDDO
             ENDDO
!
!   find the bands that are within the range of energies
!
             is_in_range_rap=.FALSE.
             DO i=1,MINVAL(nbnd_rapk(irap,spe:lpe))
                DO n=spe,lpe
                   is_in_range_rap(i) = is_in_range_rap(i) .OR. &
                         (e_rap(start_rapk(irap,n)+i-1,n)  >= emin .AND. &
                          e_rap(start_rapk(irap,n)+i-1,n) <= emax ) 
                ENDDO
             ENDDO
             DO i=1,MINVAL(nbnd_rapk(irap,spe:lpe))
                IF (is_in_range_rap(i)) THEN
                   IF ( mod(i,2) /= 0) THEN
                      IF (ionode) WRITE (2,'(2f10.4)') (kx(n), &
                         e_rap(start_rapk(irap,n)+i-1,n), n=spe,lpe)
                   ELSE
                      IF (ionode) WRITE (2,'(2f10.4)') (kx(n), &
                         e_rap(start_rapk(irap,n)+i-1,n), n=lpe,spe,-1)
                   END IF
                   has_points(ilines,irap)=.TRUE.
                END IF
             END DO
             IF (MINVAL(nbnd_rapk(irap,spe:lpe))==0) THEN
                IF (ionode) CLOSE (UNIT = 2, STATUS='DELETE')
             ELSE
                IF (ionode) CLOSE (UNIT = 2, STATUS='KEEP')
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  ENDIF

  IF (icode==1) THEN
     WRITE(stdout,'(5x,"Bands in gnuplot format written to file ",a)') flpband
  ELSEIF (icode==2) THEN
     WRITE(stdout,'(5x,"Phonons in gnuplot format written to file ",a)') flpband
  ELSEIF (icode==3) THEN
     WRITE(stdout,'(5x,"Gruneisen parameters in gnuplot format written &
                                            &to file ",a)') flpband
  ELSEIF (icode==4) THEN
     WRITE(stdout,'(5x,"Interpolated phonons in gnuplot format written &
                                            &to file ",a)') flpband
  ENDIF
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

  CALL write_gnuplot_file(kx, e_eff, k_eff, tot_points, nbnd, emin, emax, &
                          eref, nlines, nrap, rap_eff, gcodek_eff, &
                          e_rap, nbnd_rapk, start_rapk, &
                          has_points, start_point_eff, &
                          last_point_eff, nrap_plot_eff, rap_plot_eff, &
                          label_disp_q_eff, &
                          icode, exist_rap, igeom, fileout)

  DEALLOCATE(start_point)
  DEALLOCATE(last_point)
  DEALLOCATE(is_in_range)
  DEALLOCATE(high_symmetry)
  DEALLOCATE(kx) 
  DEALLOCATE(k) 
  DEALLOCATE(e)
  DEALLOCATE(has_points)
  DEALLOCATE(nrap)
  DEALLOCATE(gcodek)
  DEALLOCATE(rap)
  DEALLOCATE(aux_ind)
  DEALLOCATE(rap_eff)
  DEALLOCATE(k_eff)
  DEALLOCATE(e_eff)
  DEALLOCATE(gcodek_eff)
  DEALLOCATE(aux_ind_eff)
  DEALLOCATE(start_point_eff)
  DEALLOCATE(last_point_eff)
  DEALLOCATE(nrap_plot_eff)
  DEALLOCATE(rap_plot_eff)
  DEALLOCATE(nbnd_rapk)
  DEALLOCATE(nbnd_count)
  DEALLOCATE(start_rapk)
  DEALLOCATE(label_disp_q_eff)
  DEALLOCATE(rapin)
  DEALLOCATE(e_rap)
  IF (identify_sur) THEN
     DEALLOCATE(lsurface_state_eff)
     DEALLOCATE(lsurface_state_rap)
  ENDIF

  IF (exist_rap) THEN
     DEALLOCATE(is_in_range_rap)
     DEALLOCATE(k_rap)
  ENDIF


  RETURN
END SUBROUTINE plotband_sub

SUBROUTINE write_gnuplot_file(kx, e, xk, nks, nbnd, emin, emax, eref, nlines, &
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
END SUBROUTINE write_gnuplot_file

SUBROUTINE identify_surface_states(nat, nbnd, nkstot, e, rap)
!
!  This routine searches, among the bands the surface states using the
!  information given in input:
!  sur_layers is the number of surface atoms,
!  sur_thr is the percentage (from 0 to 1) of charge that a state must
!          have on the surface atoms to be considered as a surface state
!  averag  is the projection on each layer of each state
!
USE kinds,            ONLY : DP
USE control_2d_bands, ONLY : averag, vacuum, nlayers, sur_thr, sur_layers, &
                             sym_divide, surface1, surface2, lsurface_state, &
                             subtract_vacuum
USE io_global,        ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: nat, nbnd, nkstot
REAL(DP),INTENT(IN) :: e(nbnd,nkstot)
INTEGER, INTENT(IN) :: rap(nbnd,nkstot)

REAL(DP) :: suggested_sur_thr, maxsum
REAL(DP), ALLOCATABLE :: sumna(:,:)
LOGICAL, ALLOCATABLE :: plot(:)
INTEGER :: na, ibnd, ik, iter, ilayers, npoints, surface_layer

CALL read_state_densities()
ALLOCATE(lsurface_state(nbnd,nkstot))
ALLOCATE(sumna(nbnd,nkstot))
ALLOCATE(plot(nlayers))
!
!  surface1 and surface2 are the surface layers
!  
!
! assume always two equal surfaces
!
IF (sur_layers * 2 > nat ) &
   CALL errore('identify_surface_states','too many surface layers',1)
plot=.FALSE.
DO na=1, sur_layers
   plot(surface1-na+1)=.TRUE.
   plot(surface2+na-1)=.TRUE.
ENDDO
WRITE(stdout,'(/,5x, "Identifing surface states using charge &
                                &density per layer" )') 
WRITE(stdout,'(5x, "with nlayers=",i5," layers per surface",/)') &
         sur_layers
DO ilayers=1,nlayers
   IF (plot(ilayers)) WRITE(stdout,'(5x, "Surface layer", i5)') ilayers
ENDDO

WRITE(stdout,'(/,5x,"Layers are sets of FFT mesh planes perpendicular to the z")')
WRITE(stdout,'(5x, "direction. The first layer contains the origin,")')
WRITE(stdout,'(5x, "the other layers continue with positive z up to &
                   &alat*celldm(3)")')
                               
!
!  plot the bands that have more than plot_thr percentage of the charge
!  on the selected atoms
!
sumna=0.0_DP
maxsum=-1.D20
DO ik=1,nkstot
   DO ibnd=1,nbnd
      DO ilayers=1,nlayers
         IF (plot(ilayers)) sumna(ibnd,ik)=sumna(ibnd,ik)+averag(ilayers,1,ibnd, ik)
      ENDDO
      IF (subtract_vacuum) sumna(ibnd,ik)=sumna(ibnd,ik)-vacuum(1,ibnd,ik)
      IF (sumna(ibnd,ik)>maxsum) maxsum=sumna(ibnd,ik)
   ENDDO
ENDDO

DO iter=1, 6
   suggested_sur_thr = MIN(0.6_DP, MAX(0.2_DP, maxsum-0.15_DP))
   npoints = 0
   DO ik=1,nkstot
      DO ibnd=1,nbnd
         IF (sumna(ibnd,ik)> suggested_sur_thr) npoints=npoints+1
      ENDDO
   ENDDO
   IF (npoints > 10) EXIT
   IF (suggested_sur_thr == 0.2_DP) EXIT
END DO
!
! the suggested sur_thr is not larger than 0.6 and not smaller than 0.2
!
WRITE(stdout, '(/,5x,"Maximum density on the chosen layers", f15.3)') maxsum 
WRITE(stdout, '(5x,"Number of layers", i5)') nlayers

IF (sur_thr == 0.0_DP) THEN
   sur_thr=suggested_sur_thr
   WRITE(stdout,'(5x,"Using sur_thr =",f15.3)') sur_thr
ELSE
   WRITE(stdout,'(5x,"Suggested sur_thr =",f15.3)') suggested_sur_thr
   WRITE(stdout,'(5x,"Using sur_thr =",f15.3,/)') sur_thr
ENDIF
WRITE(stdout,'(25x,30("-"),/)') 

lsurface_state=.FALSE.
WRITE(stdout,'(5x,"Searching surface states for ",i6," k-points and ",&
                            &i6, " bands:")') nkstot, nbnd
IF (subtract_vacuum) &
WRITE(stdout,'(5x,"Vacuum charge has been subtracted")')
WRITE(stdout,'(5x,"ik,    ibnd,   charge on surface layers vacuum charge &
                                                 & energy rap")')
DO ik=1,nkstot
   DO ibnd=1,nbnd
      IF (sumna(ibnd,ik)> sur_thr) THEN
         WRITE(stdout,'(5x,2i8,3f13.7,i5)') ik, ibnd, sumna(ibnd,ik), &
                                              vacuum(1,ibnd,ik), e(ibnd,ik), &
                                              rap(ibnd,ik)
         lsurface_state(ibnd,ik)=.TRUE.
      ENDIF
   ENDDO
ENDDO

DEALLOCATE(plot)
DEALLOCATE(sumna)

RETURN
END SUBROUTINE identify_surface_states

SUBROUTINE plot_surface_states(nbnd, nks, nlines, kx, e_rap, emin, emax, eref, &
                  nrap, nbnd_rapk, start_rapk, start_point, last_point, &
                  nrap_plot, rap_plot )
!
!  This routine writes on the gnuplot scripts the command to
!  plot the surface states. The surface states must have been already
!  identified along each line 
!
USE kinds,            ONLY : DP
USE control_2d_bands, ONLY : sym_divide, lsurface_state_rap
USE control_bands,    ONLY : lsym
USE gnuplot,          ONLY : gnuplot_line, gnuplot_circle, gnuplot_write_command

IMPLICIT NONE
INTEGER, INTENT(IN)  :: nbnd, nks, nlines
INTEGER, INTENT(IN)  :: nrap(nlines), nbnd_rapk(12,nks), start_rapk(12,nks)
INTEGER, INTENT(IN)  :: start_point(nlines), last_point(nlines)
INTEGER, INTENT(IN)  :: nrap_plot(nks), rap_plot(12,nks)
REAL(DP), INTENT(IN) :: kx(nks), e_rap(nbnd, nks)
REAL(DP), INTENT(IN) :: emin, emax, eref
REAL(DP) :: x(2), y(2), ys(2), delta
INTEGER :: ilines, ibnd, ibnd1, jbnd, ik, irap, ir, nrapp
LOGICAL :: dorap

CALL gnuplot_write_command('surface_state_width=2',.FALSE.)
CALL gnuplot_write_command('surface_state_color="blue"',.FALSE.)
CALL gnuplot_write_command('surface_state_radius=0.002*xscale',.FALSE.)

DO ilines=1, nlines
   DO ik=start_point(ilines),last_point(ilines)
      IF (lsym) THEN
         DO irap=1,nrap(ilines)
            dorap=.TRUE.
            IF (sym_divide) THEN
               dorap=.FALSE.
               nrapp= nrap_plot(ik)
               DO ir=1,nrapp
                  dorap=dorap.OR.(rap_plot(ir,ik)==irap)
               END DO
               IF (nrapp==0) dorap=.TRUE.
            END IF
            IF (dorap) THEN
               DO jbnd=1, nbnd_rapk(irap,ik)
                  ibnd=start_rapk(irap,ik)+jbnd-1
                  IF (lsurface_state_rap(ibnd,ik)) THEN
                     x(1)=kx(ik)
                     ys(1)=e_rap(ibnd,ik) - eref
                     y(1)=MIN(emax, ys(1))
                     y(1)=MAX(emin, y(1))
                     IF (.NOT.(y(1)==emax .OR. y(1)==emin) ) THEN
                        CALL gnuplot_circle(x(1),y(1),'surface_state_radius',&
                                       '2', &
                                       'surface_state_color')
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO ibnd=1, nbnd
            IF (lsurface_state_rap(ibnd,ik)) THEN
               x(1)=kx(ik)
               ys(1)=e_rap(ibnd,ik) - eref
               y(1)=MIN(emax, ys(1))
               y(1)=MAX(emin, y(1))
               IF (.NOT.(y(1)==emax .OR. y(1)==emin) ) THEN
                   CALL gnuplot_circle(x(1),y(1),'surface_state_radius',&
                                       '2', &
                                       'surface_state_color')
               ENDIF
            ENDIF
         ENDDO
      END IF
   ENDDO
ENDDO

RETURN
END SUBROUTINE plot_surface_states
