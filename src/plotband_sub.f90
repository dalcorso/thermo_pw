!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! Copyright (C) 2016 Andrea Dal Corso
!                    A major rewriting and clean up of the routine.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plotband_sub(icode, filedata, filerap, fileout, &
                                      gnu_filename, filenameps )
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
  USE constants, ONLY : bohr_radius_si, pi
  USE control_bands, ONLY : emin_input, emax_input 
  USE control_grun,  ONLY : grunmin_input, grunmax_input
  USE control_dosq,  ONLY : freqmin_input, freqmax_input
  USE control_paths, ONLY : label_disp_q, nqaux, high_sym_path, disp_nqs, &
                            nrap_plot, rap_plot, long_path
  USE control_2d_bands, ONLY : nkz, identify_sur, lprojpbs, &
                               sym_divide, lsurface_state, lsurface_state_rap, &
                               force_bands
  USE cell_base,     ONLY : celldm, tpiba
  USE point_group,   ONLY : convert_rap_new, group_name_schoenflies
  USE control_gnuplot, ONLY : gnuplot_command, lgnuplot
  USE gnuplot,         ONLY : gnuplot_end, gnuplot_print_objects
  USE ions_base,     ONLY : nat
  USE spin_orb,      ONLY : lspinorb
  USE thermo_sym,    ONLY : code_group_save
  USE io_bands,      ONLY : read_parameters, read_representations, read_bands
  USE mp,            ONLY : mp_bcast
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_images,     ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: icode
  CHARACTER(LEN=256), INTENT(IN) :: filedata, filerap, fileout, &
                                               gnu_filename, filenameps
!
!  path variables
!
  INTEGER :: nks = 0, nbnd = 0,  nks_rap = 0, nbnd_rap = 0, nks_         
  INTEGER, ALLOCATABLE :: rap(:,:), gcodek(:), aux_ind(:), gcodek_ext(:), &
                          ptypek(:,:), lprojk(:)                   
  REAL(DP), ALLOCATABLE :: e(:,:), k(:,:), k_rap(:,:), gaugek(:,:) 
  LOGICAL, ALLOCATABLE :: high_symmetry(:), same_next(:)           
!
! effective path variables
!
  INTEGER :: tot_points, tot_points_
  REAL(DP), ALLOCATABLE :: kx(:), k_eff(:,:), e_eff(:,:), gaugek_eff(:,:)
  INTEGER, ALLOCATABLE  :: rap_eff(:,:), gcodek_eff(:), aux_ind_eff(:), &
                           gcodek_ext_eff(:), ptypek_eff(:,:), lprojk_eff(:), &
                           nrap_plot_eff(:), rap_plot_eff(:,:)
  LOGICAL, ALLOCATABLE  :: same_next_eff(:), lsurface_state_eff(:,:)
  CHARACTER(LEN=12), ALLOCATABLE :: point_group_path(:,:)
!
!  lines control
!
  INTEGER :: nlines
  INTEGER, ALLOCATABLE :: start_point(:), last_point(:), nrap(:), &
                          nbnd_rapk_min(:,:), start_point_eff(:), &
                          last_point_eff(:), label_disp_q_eff(:), &
                          projective(:)
  LOGICAL, ALLOCATABLE :: is_in_range(:,:), is_in_range_rap(:,:,:), &
                          has_points(:,:), dorap(:,:), with_lines(:)
!
! bands ordered with representation
!
  INTEGER, ALLOCATABLE :: start_rapk(:,:), nbnd_rapk(:,:) 
  REAL(DP), ALLOCATABLE :: e_rap(:,:)
!
!  Auxiliary
! 
  REAL(DP), PARAMETER :: eps=1.d-4

  REAL(DP) :: emin, emax, ymax, ymin, dxmod, dxmod_save, eref, dgap, &
              factor_dx, sizeb, sizec, xscale

  INTEGER, ALLOCATABLE :: rapin(:),  nbnd_count(:)

  INTEGER :: code_group_line, code_group_ext_line, ilines, irap, ibnd, &
             i, n, ik, spe, lpe, nbc, iq, ishift, ir, count0, &
             nrapp, cpe, start_shift, last_shift, ncentral, nlines_, &
             iunout, ios

  LOGICAL :: exist_rap, type1, lso, print_eref, norap

  CHARACTER(LEN=256) :: filename, ylabel
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  CHARACTER(LEN=11) :: group_name

  INTEGER :: find_free_unit
  INTEGER :: ierr, system

  IF ( my_image_id /= root_image ) RETURN

  lso=.FALSE.
  print_eref=.FALSE.
  IF (icode==1.AND.lspinorb) lso=.TRUE.

  CALL read_parameters(nks, nbnd, filedata)
!
!   if there are less than four points one cannot do a plot so we exit
!
  IF (nks<4) RETURN

  ALLOCATE (e(nbnd,nks))
  ALLOCATE (k(3,nks)) 
  ALLOCATE (k_rap(3,nks))

  ALLOCATE (high_symmetry(nks))
  ALLOCATE (rap(nbnd,nks))
  ALLOCATE (gcodek(nks))
  ALLOCATE (aux_ind(nks))
  ALLOCATE (gcodek_ext(nks))
  ALLOCATE (lprojk(nks))
  ALLOCATE (gaugek(48,nks))
  ALLOCATE (ptypek(3,nks))
  ALLOCATE (same_next(nks))

  CALL read_bands(nks, nbnd, k, e, filedata)

  CALL read_representations(nks, nbnd, k_rap, rap, high_symmetry,      &
                    gcodek, aux_ind, gcodek_ext, ptypek, lprojk, &
                    same_next, gaugek, exist_rap, filerap)
  DO n=1, nks
     IF (ABS(k(1,n)-k_rap(1,n))+ABS(k(2,n)-k_rap(2,n))+  &
         ABS(k(3,n)-k_rap(3,n)) > eps.AND.icode /=3.AND.icode/=4 ) THEN
         WRITE(stdout,'(/,5x,"Incompatible k points in rap file")')
         exist_rap=.FALSE.
     END IF
  ENDDO
!
!   If we do not have a representation file we set default values 
!   for the representations
!
  IF (.NOT. exist_rap) THEN
     high_symmetry(1:nks)=high_sym_path(1:nks)
     rap=-1
     gcodek=0
     aux_ind=0
     gcodek_ext=0
     lprojk=0
     gaugek=0.0_DP
     ptypek=1
     same_next=.FALSE.
  END IF

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
!  We start by converting the representations for projected band structure
!  calculations.
!
   IF (nkz > 1 .AND. lprojpbs .AND. sym_divide .AND. exist_rap) &
      CALL convert_rap_surface(nbnd, nks, nkz, high_symmetry, gcodek, aux_ind, &
                          gcodek_ext, ptypek, rap, gaugek, lprojk)
!
! At this point we should count how many lines we have and
! set the first and the last point of each line
!
  ALLOCATE (start_point(nks))
  ALLOCATE (last_point(nks))

  nks_=nks / nkz
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
           IF (.NOT. high_symmetry(n-1)) last_point(nlines)=n
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
  WRITE(stdout,'(5x,"Number of lines:", i4, " Total number of points:",i8,/)') &
                                        nlines, tot_points
!
!   First the variables that characterize the lines
!
  ALLOCATE (point_group_path(nlines,3))
  ALLOCATE (dorap(12,nlines))
  ALLOCATE (is_in_range(nbnd,nlines))
  ALLOCATE (has_points(12,nlines))
  ALLOCATE (with_lines(nlines))
  ALLOCATE (projective(nlines))
  ALLOCATE (nrap(nlines))
  ALLOCATE (nbnd_rapk_min(12,nlines))
  ALLOCATE (is_in_range_rap(nbnd,12,nlines))
!
!  then we have to define a new path, that we call effective path,
!  made of tot_points where each starting and last point of each line 
!  must be separated, so that it can have a different point group, label 
!  and set of representations. Actually the starting and last point of 
!  each line must be classified with the same representations of the line.
!  The effective path has the same variables as the real path with the
!  exception of high_symmetry which is no more necessary. Instead we
!  add the representation information nrap_plot_eff and rap_plot_eff
!  that say in each point which representations have to be plotted and
!  label_disp_q_eff which says where the labels go in the effective path. 
!
  ALLOCATE (k_eff(3,tot_points)) 
  ALLOCATE (e_eff(nbnd,tot_points)) 

  ALLOCATE (rap_eff(nbnd,tot_points))
  ALLOCATE (gcodek_eff(tot_points))
  ALLOCATE (aux_ind_eff(tot_points))
  ALLOCATE (gcodek_ext_eff(tot_points))
  ALLOCATE (ptypek_eff(3,tot_points))
  ALLOCATE (gaugek_eff(48,tot_points))
  ALLOCATE (lprojk_eff(tot_points))
  ALLOCATE (same_next_eff(tot_points))

  ALLOCATE (nrap_plot_eff(tot_points))
  ALLOCATE (rap_plot_eff(12,tot_points))
  ALLOCATE (label_disp_q_eff(nqaux))
!
!  the lines must be characterized also on the effective path
!
  ALLOCATE (start_point_eff(nlines))
  ALLOCATE (last_point_eff(nlines))
!
! and if there are representations the bands must be reordered. See below
!
  ALLOCATE (nbnd_rapk(12,tot_points))
  ALLOCATE (start_rapk(12,tot_points))
  ALLOCATE (e_rap(nbnd,tot_points))
!
!  Finally we need a few additional working variables
!
  ALLOCATE (nbnd_count(tot_points))
  ALLOCATE (rapin(nbnd))

!
!  and the x coordinate in the effective path
!
  ALLOCATE (kx(tot_points)) 
!
!  some variables that identify surface states
!
  IF (identify_sur) THEN
     ALLOCATE(lsurface_state_eff(nbnd,tot_points))
     ALLOCATE(lsurface_state_rap(nbnd,tot_points))
  ENDIF
!
! Now build the effective path and copy inside the relevant variables
!
  ik=0
  label_disp_q_eff=0
  DO ilines=1,nlines
     ncentral=(start_point(ilines)+last_point(ilines))/2
     DO n=start_point(ilines), last_point(ilines)
        ik=ik+1
        IF (ik>tot_points) CALL errore('plotband_sub','problem with points',1)

        k_eff(:,ik)=k(:,n)
        e_eff(:,ik)=e(:,n)

        rap_eff(:,ik)=rap(:,n)
        gcodek_eff(ik)=gcodek(n)
        aux_ind_eff(ik)=aux_ind(n)
        gcodek_ext_eff(ik)=gcodek_ext(n)
        ptypek_eff(:,ik)=ptypek(:,n)
        gaugek_eff(:,ik)=gaugek(:,n)
        lprojk_eff(ik)=lprojk(n)
        same_next_eff(ik)=same_next(n)

        nrap_plot_eff(ik)=nrap_plot(n)
        rap_plot_eff(:,ik)=rap_plot(:,n)

        IF (identify_sur) lsurface_state_eff(:,ik)=lsurface_state(:,n)
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
!
!  search if this n has a label and put the label on ik
!
        DO iq=1, nqaux
           IF (label_disp_q(iq)==n) label_disp_q_eff(iq)=ik
        ENDDO
     END DO
  END DO

  IF (ik /= tot_points) CALL errore('plotband_sub','Missing points',1)
!
!   Now set the limits of the y coordinate of the plot
!
!
! Find the maximum and minimum of the plot 
!
  emin=1.d10
  emax=-1.d10
  DO n=1, tot_points
     DO ibnd=1,nbnd
        emin = MIN(emin, e_eff(ibnd,n))
        emax = MAX(emax, e_eff(ibnd,n))
     ENDDO
  ENDDO
!
!  if necessary shift the bands and adjust the maximum and minimum
!
  IF (icode==1) THEN
!
!   bands can be shifted
!
     CALL compute_eref_band(e_eff, nbnd, eref, print_eref)

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
!  Now for each line set the variables relevant for the line
!  These are:
!
!  point_group_path(nlines,3)   ! the label of the point group, initial
!                               ! final and middle
!  nrap(nlines)                 ! the maximum number of the representation
!                               ! in each line 
!  dorap(12,nlines)             ! For each representation if it has to be
!                               ! plotted
!  is_in_range(nbnd,nlines)     ! for each band if it is between emin and emax
!  is_in_range_rap(nbnd,12,nlines) ! for each band in a given rap if it is 
!                               ! between emin and emax
!  has_points(12,nlines)        ! for each representation if there are bands
!                               ! of that representations
!  with_lines(nlines)           ! if .TRUE. join the bands with lines
!  projective(nlines)           ! if 1,2,3 the line is projective or zone border
!
!  nbnd_rapk_min(12,nlines)     ! the number of bands for each rap. The
!                               ! minimum values among the points of a line.
!
!
  WRITE(stdout,'(/,5x,"Representations per line:")')

  nbnd_count=0
  start_rapk=0
  is_in_range = .FALSE.
  is_in_range_rap=.FALSE.
  has_points=.FALSE.

  DO ilines=1,nlines
     spe=start_point_eff(ilines)
     lpe=last_point_eff(ilines) 
     cpe=(spe+lpe)/2.0_DP
!
!   defaults if there is not the representation file
!
     with_lines(ilines)=.TRUE.
!
!   gruneisen parameters are plotted as points if there is no representation
!
     IF (icode==3.AND..NOT.exist_rap) with_lines(ilines)=.FALSE.

     projective(ilines)=0
     dorap(:,ilines)=.TRUE.
     point_group_path(ilines,1:3)=' '
     nrap(ilines)=1
!
!   by default all the bands in a single representation
!
     DO n=spe,lpe
        e_rap(:,n)=e_eff(:,n)
        nbnd_rapk(1,n)=nbnd
        start_rapk(1,n)=1
        IF (identify_sur) lsurface_state_rap(:,n)=lsurface_state_eff(:,n)
     ENDDO
!
!    Now check if some band is in the plotting range  
!
     DO ibnd=1,nbnd
        is_in_range(ibnd,ilines) = ANY ( e_eff(ibnd,spe:lpe) >= emin .AND. &
                                         e_eff(ibnd,spe:lpe) <= emax )
     ENDDO

     IF (exist_rap) THEN
        code_group_line=gcodek_eff(spe+1)
        code_group_ext_line=gcodek_ext_eff(spe+1)
!
!   set the point group name of the first, last and central points of each
!   line
!
        CALL group_name_schoenflies(gcodek_eff(spe),point_group_path(ilines,1))
        CALL group_name_schoenflies(gcodek_eff(lpe),point_group_path(ilines,2))
        CALL group_name_schoenflies(gcodek_eff(cpe),point_group_path(ilines,3))
!
!   Then find if the line has projective representations and if 
!   the representations have switched
!   projective
!      0     not projective representations
!      1     projective not switched
!      2     projective switched
!      3     zone border point but standard representation
!
        IF (lprojk_eff(cpe)==1) THEN
           IF (icode==1) THEN
              IF (lspinorb.AND.ptypek_eff(1,cpe)==1.OR. &
                    (.NOT.lspinorb.AND.ptypek_eff(1,cpe)==-1)) THEN
                 projective(ilines)=2
              ELSE
                 projective(ilines)=1
              ENDIF
           ELSE
              IF (ptypek_eff(1,cpe)==-1) THEN
                 projective(ilines)=2
              ELSE
                 projective(ilines)=1
              ENDIF
           ENDIF
        ELSEIF (lprojk_eff(cpe)==2) THEN
!
!   zone border point but standard symmetry analysis
!
           projective(ilines)=3
        ENDIF
!
!  Now convert the representations of the border points to those of the
!  central point group
!
        DO n=spe,lpe
           IF (gcodek_eff(n) /= code_group_line) THEN
!
!   it is the point before or after the one that change symmetry 
!   that carry the information about the group subgroup relationship
!   (aux_ind)
!
              IF (n==spe) THEN
                 rapin(:)=rap_eff(:,n)
!                WRITE(6,*) 'first', k_eff(1,n), k_eff(2,n), k_eff(3,n)
                 CALL convert_rap_new(nbnd,rapin,rap_eff(1,n),       &
                             gcodek_ext_eff(n), code_group_ext_line, &
                             aux_ind_eff(n+1), ptypek_eff(1,n),      &
                             ptypek_eff(1,n+1), gaugek_eff(1,n),     &
                             gaugek_eff(1,n+1))
              ELSEIF (n==lpe) THEN
                 rapin(:)=rap_eff(:,n)
!                WRITE(6,*) 'last', k_eff(1,n), k_eff(2,n), k_eff(3,n)
                 CALL convert_rap_new(nbnd,rapin,rap_eff(1,n),       &
                             gcodek_ext_eff(n), code_group_ext_line, &
                             aux_ind_eff(n-1), ptypek_eff(1,n),      &
                             ptypek_eff(1,n-1), gaugek_eff(1,n),     &
                             gaugek_eff(1,n-1))
              ELSE
                 CALL errore('plotband_sub','unexpected change of symmetry',1)
              ENDIF
           ENDIF
        ENDDO
!
!    check if the initial and final points have representations
!
        start_shift=0 
        last_shift=0 
        DO ibnd=1,nbnd
           IF (rap_eff(ibnd,spe) <= 0) start_shift=1
           IF (rap_eff(ibnd,lpe) <= 0) last_shift=1
        ENDDO
!
!   if not we copy those of the closest point. This is equivalent to make all 
!   avoided crossing between the line border and the point after of before.
!
        IF (start_shift==1) rap_eff(:,spe)=rap_eff(:,spe+1)
        IF (last_shift==1) rap_eff(:,lpe)=rap_eff(:,lpe-1)
!
!   Determine for each line which is the representation with the highest
!   number 
!
        nrap(ilines)=1
        count0=0
        norap=.FALSE.
        DO n=spe, lpe
           DO ibnd=1,nbnd
              nrap(ilines)=MAX(nrap(ilines),rap_eff(ibnd,n))
              norap=(norap.OR.(rap_eff(ibnd,n)< 1))
              IF (rap_eff(ibnd,n)<1) count0=count0+1
           ENDDO
        ENDDO
!
!   if along this line too many bands have no symmetry classification, 
!   disable the plot of the representations and put dots in the plot along
!   this line
!
        IF (count0/(lpe-spe) > 4) THEN
           nrap(ilines)=1
           with_lines(ilines)=.FALSE.
        ENDIF

        IF (nrap(ilines) > 12) CALL errore("plotband_sub",&
                                           "Too many representations",1)
        WRITE(stdout,'(5x, "Line ", i4, " maximum repres. number", i4, &
                     " point group ", i4, 2x, a11)') ilines, nrap(ilines),    &
                       code_group_ext_line, TRIM(group_name(code_group_line))

        IF (icode==3.AND.norap) with_lines(ilines)=.FALSE.
!
!   here set dorap, for each line which representations have to be plotted.
!   This is controlled by nrap_plot_eff and rap_plot_eff
!
        nrapp= nrap_plot_eff(spe)
        DO irap=1, nrap(ilines)
           IF (sym_divide) THEN
              dorap(irap,ilines)=(nrapp==0)
              DO ir=1,nrapp
                 dorap(irap,ilines) = dorap(irap,ilines) &
                       .OR.(rap_plot_eff(ir,spe)==irap)
              ENDDO
           ENDIF
        ENDDO
!
!   Finally manage the representations. We count how many bands
!   there are for each representation and reorder the bands so that there
!   are first those of representation 1, then those of representation 2 etc.
!   nbnd_rapk(irap,n) says how many bands of each representation per k point
!   start_rapk(irap,n) where the bands of representation irap starts
!                      in the reordered band list
!   
        
        IF (nrap(ilines)/=1) THEN
!
!  For each k point along this line selects only the bands which belong
!  to the irap representation. 
!
           DO irap=1,nrap(ilines)
              DO n = spe, lpe 
                 nbnd_rapk(irap, n)=0
                 DO ibnd=1,nbnd
                    IF (rap_eff(ibnd,n)==irap) THEN
                       nbnd_rapk(irap,n) = nbnd_rapk(irap,n) + 1
                       nbnd_count(n)=nbnd_count(n)+1
                       nbc=nbnd_count(n)
                       IF (nbnd_rapk(irap,n)==1) start_rapk(irap,n)=nbc
                       e_rap(nbc,n)=e_eff(ibnd,n)
                       IF (identify_sur) &
                          lsurface_state_rap(nbc,n)=lsurface_state_eff(ibnd,n)
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
!
!   For a given representation we take all the bands that are present 
!   in all k point along a line. Note that some part of these bands could
!   be outside the emin-emax range if this is provided from input.
!
           DO irap=1,nrap(ilines)
              nbnd_rapk_min(irap,ilines)=MINVAL(nbnd_rapk(irap,spe:lpe)) 
           ENDDO

           DO irap=1,nrap(ilines)
              DO n=spe,lpe
                 DO i=1,nbnd_rapk_min(irap,ilines)
                    ibnd=start_rapk(irap,n)+i-1
                    is_in_range_rap(i,irap,ilines)=&
                                            is_in_range_rap(i,irap,ilines) &
                     .OR.(e_rap(ibnd,n) >= emin .AND. e_rap(ibnd,n) <= emax ) 
                 ENDDO
              ENDDO
           ENDDO
!
!    check if the representations have some point
!
           DO irap=1,nrap(ilines)
              DO ibnd=1,nbnd_rapk_min(irap,ilines)
                 IF (is_in_range_rap(ibnd,irap,ilines)) &
                                      has_points(irap,ilines)=.TRUE.
              ENDDO
           ENDDO
        ELSE
           DO ibnd=1,nbnd
              IF (is_in_range(ibnd,ilines)) has_points(1,ilines)=.TRUE.
           ENDDO
        ENDIF
     ELSE
        WRITE(stdout,'(5x, "Line ", i7, " Representations not available")')
        DO ibnd=1,nbnd
           IF (is_in_range(ibnd,ilines)) has_points(1,ilines)=.TRUE.
        ENDDO
     ENDIF
  ENDDO
!
!  In the pbs case we cannot use the symmetry information of the first and
!  last layer if sym_divide has not transformed them
!
  IF (nkz > 1 .AND. lprojpbs .AND. (.NOT.sym_divide) .AND. exist_rap ) &
     CALL convert_proj_surface(nlines, tot_points, point_group_path,   &
                               projective, lprojk_eff)
!
!   here try to estimate the size of the path
!
  IF (celldm(2)>0.0_DP) THEN
     sizeb=celldm(2)
  ELSE
     sizeb=1.0_DP
  ENDIF

  IF (celldm(3)>0.0_DP) THEN
     sizec=celldm(3)
  ELSE
     sizec=1.0_DP
  ENDIF
  factor_dx = MAX(6.0_DP, 2.0_DP * sizeb, 2.0_DP * sizec, &
                                       2.0_DP/sizeb, 2.0_DP/sizec)

  dxmod_save = SQRT( (k(1,2)-k(1,1))**2 +  &
                     (k(2,2)-k(2,1))**2 +  &
                     (k(3,2)-k(3,1))**2 )
!
!  and set a spacing for panels that are not contiguous.
!
  dgap=15.0_DP * dxmod_save
  IF (.NOT.long_path) dgap=0.0_DP
!
!  now compute the x coordinate on the plot, but use the effective k
!  points
!
  kx(1) = 0.d0
  nks_=tot_points/nkz
  DO n=2,tot_points
     dxmod=sqrt ( (k_eff(1,n)-k_eff(1,n-1))**2 + &
                  (k_eff(2,n)-k_eff(2,n-1))**2 + &
                  (k_eff(3,n)-k_eff(3,n-1))**2 )
     IF (MOD(n-1,nks_)==0) THEN
!
!   This is a PBS calculation and now we switch to the next plot
!
         kx(n)=0.0_DP
     ELSEIF (dxmod > factor_dx*dxmod_save) THEN
!
!   A big jump in dxmod is a sign that the point k(:,n) and k(:,n-1)
!   are quite distant and belong to two different lines. We put them on
!   the same point in the graph if the two k point belong to the same star
!   up to a reciprocal lattice vector, or we introduce a gap.
!
        IF (same_next_eff(n-1)) THEN
           kx(n)=kx(n-1) 
        ELSE
           kx(n)=kx(n-1) + dgap
        ENDIF
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
  WRITE(stdout,'(/,5x,"x coordinate on the plot: ",/)')
  DO ilines=1,nlines
     WRITE(stdout,'(5x,"Line ",i5," starts at",f13.6," ends at ",f13.6)') &
           ilines, kx(start_point_eff(ilines)), kx(last_point_eff(ilines))
  ENDDO
  !
  IF (ionode) iunout=find_free_unit()
!
!   In this case we write a different file for each line and for each
!   representation. Each file contains the bands of that representation.
!   The file is called filename.#line.#rap
!
  DO ilines=1,nlines
     spe=start_point_eff(ilines)
     lpe=last_point_eff(ilines) 
     IF (nrap(ilines)==1) THEN
!
!   Along this line the symmetry decomposition has not been done or
!   the group has only one representation.
!   Plot all the bands as in the standard case
!
        IF (has_points(1,ilines)) THEN
           filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines))

           IF (ionode) OPEN(UNIT=iunout, FILE=TRIM(filename), &
                            FORM='formatted', STATUS='unknown', IOSTAT=ios)
           CALL mp_bcast(ios, ionode_id, intra_image_comm)
           IF (ionode) THEN
              DO ibnd=1,nbnd
                 IF (is_in_range(ibnd,ilines)) THEN
                    WRITE (iunout,'(2f14.7)') (kx(n), e_eff(ibnd,n), &
                                                                  n=spe, lpe)
                    WRITE(iunout,*)
                 ENDIF
              ENDDO
              CLOSE (UNIT=iunout, STATUS='KEEP')
           ENDIF
        ENDIF
     ELSE
        DO irap=1, nrap(ilines)
!
!     along this line there are several representations
!     open a file for each representation
!
           IF (has_points(irap,ilines)) THEN
              filename=TRIM(fileout) // "." // TRIM(int_to_char(ilines)) &
                                    //  "." // TRIM(int_to_char(irap))
              IF (ionode) OPEN (UNIT=iunout, FILE=TRIM(filename), &
                             FORM='formatted', STATUS='unknown', IOSTAT=ios)
              CALL mp_bcast(ios, ionode_id, intra_image_comm)
              CALL errore("plotband_sub","opening file"//TRIM(filename), &
                                                       ABS(ios)) 
              IF (ionode) THEN
                 DO i=1,nbnd_rapk_min(irap,ilines)
                    IF (is_in_range_rap(i,irap,ilines)) THEN
                       WRITE (iunout,'(2f14.7)') (kx(n), &
                             e_rap(start_rapk(irap,n)+i-1,n), n=spe,lpe)
                       WRITE(iunout,*) 
                    ENDIF
                 ENDDO
                 CLOSE (UNIT = iunout, STATUS='KEEP')
              ENDIF
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !
  IF (emin_input /= 0.0_DP) THEN
     ymin=emin_input
  ELSE
     ymin=emin-eref
  ENDIF

  IF (emax_input /= 0.0_DP) THEN
     ymax=emax_input 
  ELSE
     ymax=emax-eref
  ENDIF

  xscale=tpiba / bohr_radius_si / 1.D10
  IF (icode==1) THEN
     ylabel='Energy (eV)'
  ELSEIF (icode==2.OR.icode==4) THEN
     ylabel='Frequency (cm^{-1})'
  ELSEIF (icode==3) THEN
     ylabel='Mode-Gr\374neisen parameters  {/Symbol g}_{/Symbol n}&
                 &({/Helvetica-Bold q})'
  ENDIF

  CALL initialize_plot_dispersion(kx, tot_points, xscale, ymin, ymax, &
             eref, print_eref, nlines, start_point_eff, last_point_eff, &
             label_disp_q_eff, point_group_path, projective, lprojk_eff, &
             ylabel, gnu_filename, filenameps)

  IF (lprojpbs) CALL proj_band_structure(kx, e_eff, tot_points, nbnd, &
                    ymin, ymax, eref, e_rap, nrap, nbnd_rapk, start_rapk, &
                    nlines, start_point_eff, last_point_eff, &
                    nrap_plot_eff, rap_plot_eff )

  IF (identify_sur) CALL plot_surface_states(nbnd, tot_points, nlines, kx, &
                           e_rap, ymin, ymax, eref, nrap, nbnd_rapk, &
                           start_rapk, start_point_eff, last_point_eff, &
                           nrap_plot_eff, rap_plot_eff)

  IF ( nkz == 1 .OR. ( nkz > 1 .AND. .NOT. lprojpbs) .OR. force_bands ) THEN
     CALL plot_dispersion(nlines, nrap, has_points, dorap, with_lines, fileout)
  ELSE
     CALL gnuplot_print_objects()
  ENDIF

  CALL gnuplot_end()
  IF (lgnuplot.AND.ionode) &
     ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

!  IF (lgnuplot.AND.ionode) &
!     CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '//TRIM(gnu_filename),&
!                                                             WAIT=.FALSE.)
!
!  First deallocate the variables that define the path read from file
!
  DEALLOCATE (e)
  DEALLOCATE (k)

  DEALLOCATE (high_symmetry)
  DEALLOCATE (rap)
  DEALLOCATE (gcodek)
  DEALLOCATE (aux_ind)
  DEALLOCATE (gcodek_ext)
  DEALLOCATE (lprojk)
  DEALLOCATE (gaugek)
  DEALLOCATE (ptypek)
  DEALLOCATE (same_next)
!
!  then the variables that define the lines
!
  DEALLOCATE (start_point)
  DEALLOCATE (last_point)
  DEALLOCATE (point_group_path)
  DEALLOCATE (dorap)
  DEALLOCATE (is_in_range)
  DEALLOCATE (has_points)
  DEALLOCATE (with_lines)
  DEALLOCATE (projective)
  DEALLOCATE (nrap)
  DEALLOCATE (nbnd_rapk_min)
  DEALLOCATE (is_in_range_rap)
!
!  and the variable that define the effective path
!
  DEALLOCATE (k_eff)
  DEALLOCATE (e_eff)

  DEALLOCATE (rap_eff)
  DEALLOCATE (gcodek_eff)
  DEALLOCATE (aux_ind_eff)
  DEALLOCATE (gcodek_ext_eff)
  DEALLOCATE (ptypek_eff)
  DEALLOCATE (gaugek_eff)
  DEALLOCATE (lprojk_eff)
  DEALLOCATE (same_next_eff)
  DEALLOCATE (nrap_plot_eff)
  DEALLOCATE (rap_plot_eff)
  DEALLOCATE (label_disp_q_eff)

  DEALLOCATE (start_point_eff)
  DEALLOCATE (last_point_eff)
!
! The variables that rearrange the bands
!
  DEALLOCATE (nbnd_rapk)
  DEALLOCATE (start_rapk)
  DEALLOCATE (e_rap)
!
!  Additional working variables
!
  DEALLOCATE (nbnd_count)
  DEALLOCATE (rapin)
!
!  the x coordinate of the path and the k_rap
!
  DEALLOCATE (kx) 
  DEALLOCATE (k_rap)

  IF (identify_sur) THEN
     DEALLOCATE (lsurface_state_eff)
     DEALLOCATE (lsurface_state_rap)
  ENDIF

  RETURN
END SUBROUTINE plotband_sub

