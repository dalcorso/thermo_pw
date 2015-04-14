!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE thermo_readin()
  !-----------------------------------------------------------------------
  !
  !  This routine reads the input of pw.x and an additional namelist in which 
  !  we specify the calculations to do, the files where the data are read
  !  and written and the variables of the calculations.
  !
  USE kinds,                ONLY : DP
  USE thermo_mod,           ONLY : what, ngeo, step_ngeo, ntry
  USE control_thermo,       ONLY : outdir_thermo, flevdat,        &
                                   flfrc, flfrq, fldos, fltherm,  &
                                   flanhar, filband, flkeconv,    &
                                   flnkconv, flgrun
  USE temperature,          ONLY : tmin, tmax, deltat, ntemp
  USE ifc,                  ONLY : nq1_d, nq2_d, nq3_d, ndos_input, deltafreq, &
                                   zasr, freqmin_input, freqmax_input
  USE input_parameters,     ONLY : outdir, ibrav, forc_conv_thr
  USE read_input,           ONLY : read_input_file
  USE command_line_options, ONLY : input_file_ 
  USE control_paths,        ONLY : xqaux, wqaux, wqauxr, npk_label, letter, &
                                   label_list, nqaux, q_in_band_form, &
                                   q_in_cryst_coord, q2d, point_label_type, &
                                   disp_q, disp_wq, disp_nqs, npx, &
                                   label_disp_q, letter_path, nrap_plot_in, &
                                   rap_plot_in
  USE control_gnuplot,      ONLY : flgnuplot, flpsband, flpsdisp, &
                                   flpsmur, flpsdos, flpstherm, flpsanhar, &
                                   flpskeconv, flpsnkconv, flpsgrun,       &
                                   lgnuplot,           &
                                   flpbs, flprojlayer, gnuplot_command
  USE control_2d_bands,     ONLY : lprojpbs, nkz, gap_thr, sym_divide, &
                                   identify_sur, sur_layers, sur_thr, &
                                   force_bands, only_bands_plot, dump_states, &
                                   subtract_vacuum
  USE control_asy,          ONLY : flasy, lasymptote, asymptote_command
  USE control_bands,        ONLY : flpband, emin_input, emax_input, nbnd_bands,&
                                   lsym 
  USE control_grun,         ONLY : flpgrun, grunmin_input, grunmax_input
  USE control_conv,         ONLY : nke, deltake, nkeden, deltakeden, &
                                   nnk, deltank, nsigma, deltasigma
  USE control_mur,          ONLY : vmin_input, vmax_input, deltav, nvol
  USE control_elastic_constants, ONLY : at_save, tau_save, delta_epsilon, &
                                        ibrav_save, ngeo_strain, frozen_ions, &
                                        fl_el_cons, elastic_algorithm
  USE control_piezoelectric_tensor, ONLY : nosym_save
  USE piezoelectric_tensor, ONLY : nppl
  USE control_pwrun,        ONLY : celldm_save
  USE cell_base,            ONLY : at, bg, celldm
  USE ions_base,            ONLY : nat, tau
  USE symm_base,            ONLY : nosym
  USE mp_world,             ONLY : world_comm
  USE mp_images,            ONLY : nimage, my_image_id, root_image
  USE parser,               ONLY : read_line, parse_unit
  USE environment,          ONLY : environment_end
  USE mp_global,            ONLY : mp_global_end
  USE io_global,            ONLY : ionode, meta_ionode, meta_ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  INTEGER :: ios
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=512) :: dummy
  CHARACTER(LEN=256), ALLOCATABLE :: input(:)
  CHARACTER(LEN=256) :: asy_filename
  INTEGER, ALLOCATABLE :: iun_image(:)
  INTEGER :: image
  INTEGER :: iq, ipol, i, j, k, igeo, icont
  INTEGER :: iun_thermo, parse_unit_save

  INTEGER :: nch, nrp, ierr
  LOGICAL :: tend, terr, read_paths, set_internal_path, set_2d_path, exst
  CHARACTER(LEN=256) :: input_line, buffer
  REAL(DP) :: wq0
  !
  NAMELIST / input_thermo / what, ngeo, zasr,               &
                            flfrc, flfrq, fldos, fltherm,   &
                            flanhar, filband, flkeconv,     &
                            flnkconv,                       &
                            step_ngeo,                      &
                            nq1_d, nq2_d, nq3_d,            &
                            tmin, tmax, deltat, ntemp,      &
                            freqmin_input, freqmax_input,   &
                            ndos_input, deltafreq,          &
                            q2d, q_in_band_form,            &
                            q_in_cryst_coord,               &
                            point_label_type,               &
                            nbnd_bands,                     &
                            nppl,                           &
                            lsym,                           &
                            ntry,                           &
                            npx,                            &
                            lprojpbs, nkz, gap_thr,         &
                            only_bands_plot,                &
                            sur_thr, sur_layers,            &
                            sym_divide, identify_sur,       &
                            subtract_vacuum,                &
                            force_bands,                    &
                            dump_states,                    &
                            flevdat,                        &
                            flpband, flpgrun,               &
                            flgnuplot, flpsband,            &
                            flpsdisp, flpsdos, flpstherm,   &
                            flpsanhar, flpsmur, flpskeconv, &
                            flpsnkconv, flgrun,             &
                            fl_el_cons,                     &
                            flasy, asymptote_command,       &
                            lasymptote,                     &
                            emin_input, emax_input,         &
                            vmin_input, vmax_input, deltav, &
                            elastic_algorithm,              &
                            grunmin_input, grunmax_input,   &
                            delta_epsilon, ngeo_strain,     &
                            frozen_ions,                    &
                            nvol, nke, deltake,             &
                            nkeden, deltakeden,             &
                            nnk, deltank, nsigma, deltasigma, &
                            lgnuplot, gnuplot_command
  !
  !  First read the input of thermo. This input should be in a file
  !  called thermo_control
  !
  iun_thermo=2
  parse_unit_save=parse_unit
  parse_unit=iun_thermo
  IF (meta_ionode) &
     OPEN(UNIT=iun_thermo,FILE='thermo_control',STATUS='OLD', &
                               FORM='FORMATTED', ERR=10, IOSTAT=ios )
10  CALL mp_bcast(ios, meta_ionode_id, world_comm )
    CALL errore( 'thermo_readin', 'opening thermo_control file', ABS( ios ) )

  what=' '
  ngeo=0
  step_ngeo = 0.05_DP
  ntry=2

  nq1_d=128
  nq2_d=128
  nq3_d=128
  zasr='simple'
  freqmin_input=0.0_DP
  freqmax_input=0.0_DP
  deltafreq=1.0_DP
  ndos_input=1

  tmin=1.0_DP
  tmax=800.0_DP
  deltat=3.0_DP
  ntemp=1

  nbnd_bands=0
  emin_input=0.0_DP
  emax_input=0.0_DP
  lsym=.TRUE.

  vmin_input=0.0_DP
  vmax_input=0.0_DP
  deltav=0.0_DP
  nvol=1

  nke=5
  deltake=10.0_DP
  nkeden=1
  deltakeden=100.0_DP

  nnk=5
  deltank=2 
  nsigma=1  
  deltasigma=0.005_DP

  delta_epsilon=0.005_DP
  ngeo_strain=4
  frozen_ions=.FALSE.
  elastic_algorithm='standard'

  nppl=51

  npx=8

  lprojpbs = .TRUE.
  nkz = 1 
  gap_thr = 0.1_DP
  sym_divide = .FALSE.
  identify_sur =.FALSE.
  sur_layers=0
  sur_thr=0.0_DP
  subtract_vacuum=.TRUE.
  force_bands=.FALSE.
  only_bands_plot=.FALSE.
  dump_states=.FALSE.

  grunmin_input=0.0_DP
  grunmax_input=0.0_DP

  filband='output_band.dat'
  flpband='output_pband.dat'
  flpgrun='output_pgrun.dat'
  flfrc='output_frc.dat'
  flfrq='output_frq.dat'
  fldos='output_dos.dat'
  flkeconv='output_keconv.dat'
  flnkconv='output_nkconv.dat'
  fltherm='output_therm.dat'
  flanhar='output_anhar.dat'
  flgrun='output_grun.dat'
  flevdat='output_ev.dat'
  fl_el_cons='output_el_cons.dat'
  flpbs='output_pbs'
  flprojlayer='output_projlayer'

  flgnuplot='gnuplot.tmp'
  flpsmur='output_mur.ps'
  flpsband='output_band.ps'
  flpsdisp='output_disp.ps'
  flpsdos='output_dos.ps'
  flpstherm='output_therm.ps'
  flpsanhar='output_anhar.ps'
  flpskeconv='output_keconv.ps'
  flpsnkconv='output_nkconv.ps'
  flpsgrun='output_grun.ps'

  flasy='asy_tmp'
  lasymptote=.FALSE.
  asymptote_command='asy -f pdf -noprc'

  q2d=.FALSE.
  q_in_band_form=.TRUE.
  q_in_cryst_coord=.FALSE.
  point_label_type='SC'
  lgnuplot=.TRUE.
  gnuplot_command='gnuplot'

  IF (meta_ionode) READ( iun_thermo, input_thermo, IOSTAT = ios )
  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore( 'thermo_readin', 'reading input_thermo namelist', ABS( ios ) )
  !
  CALL bcast_thermo_input()
  !
  IF (what /= 'scf_2d_bands') THEN
     nkz=1
     lprojpbs=.FALSE.
  ENDIF

  read_paths=( what=='scf_bands'.OR.what=='scf_disp'.OR.what=='plot_bz'.OR. &
               what=='mur_lc_bands' .OR. what=='mur_lc_disp' .OR. &
               what=='mur_lc_t' .OR. what=='scf_2d_bands')

  IF ( ngeo==0 ) THEN
     IF (what(1:4) == 'scf_') ngeo=1
     IF (what(1:6) == 'mur_lc') ngeo=9
     IF (what(1:4) == 'elas') ngeo=1
     IF (what(1:4) == 'piez') ngeo=1
     IF (what(1:4) == 'pola') ngeo=1
  END IF

  nqaux=0
  set_internal_path=.FALSE.
  set_2d_path=.FALSE.
  IF ( read_paths ) THEN

     IF (meta_ionode) READ (iun_thermo, *, err=200, iostat = ios) nqaux

200  CALL mp_bcast(ios, meta_ionode_id, world_comm )
     IF (ios /= 0) THEN 
        IF (what=='scf_2d_bands') THEN
           set_2d_path=.TRUE.
        ELSE
           set_internal_path=.TRUE.
        ENDIF
        goto 70
     ENDIF
     CALL mp_bcast(nqaux, meta_ionode_id, world_comm )
!
!    Reads on input the k points
!
     ALLOCATE(xqaux(3,nqaux))
     ALLOCATE(wqaux(nqaux))
     ALLOCATE(letter(nqaux))
     ALLOCATE(letter_path(nqaux))
     ALLOCATE(label_list(nqaux))
     ALLOCATE(label_disp_q(nqaux))
     IF (.NOT.q_in_band_form) ALLOCATE(wqauxr(nqaux))
     ALLOCATE(nrap_plot_in(nqaux))
     ALLOCATE(rap_plot_in(12,nqaux))
     nrap_plot_in=0
     rap_plot_in=0

     npk_label=0
     letter_path='   '
     DO iq=1, nqaux
        IF (my_image_id==root_image) &
           CALL read_line( input_line, end_of_file = tend, error = terr )
        CALL mp_bcast(input_line, meta_ionode_id, world_comm)
        CALL mp_bcast(tend, meta_ionode_id, world_comm)
        CALL mp_bcast(terr,meta_ionode_id, world_comm)
        IF (tend) CALL errore('thermo_readin','Missing lines',1)
        IF (terr) CALL errore('thermo_readin','Error reading q points',1)
        DO j=1,256   ! loop over all characters of input_line
           IF ( (ICHAR(input_line(j:j)) < 58 .AND. &   ! a digit
                 ICHAR(input_line(j:j)) > 47)      &
             .OR.ICHAR(input_line(j:j)) == 43 .OR. &   ! the + sign
                 ICHAR(input_line(j:j)) == 45 .OR. &   ! the - sign
                 ICHAR(input_line(j:j)) == 46 ) THEN   ! a dot .
!
!   This is a digit, therefore this line contains the coordinates of the
!   k point. We read it and exit from the loop on characters
!
              IF (sym_divide) THEN
                 READ(input_line,*) xqaux(1,iq), xqaux(2,iq), &
                                    xqaux(3,iq), wq0, nrap_plot_in(iq)
                 IF (nrap_plot_in(iq)>0) THEN
                    nrp=nrap_plot_in(iq)
                    READ(input_line,*) xqaux(1,iq), xqaux(2,iq), &
                          xqaux(3,iq), wq0, nrap_plot_in(iq), &
                          (rap_plot_in(i,iq), i=1,nrp)
                 END IF   
              ELSE              
                 READ(input_line,*) xqaux(1,iq), xqaux(2,iq), &
                                    xqaux(3,iq), wq0
              ENDIF
              IF (q_in_band_form) THEN
                 wqaux(iq)=NINT(wq0)
              ELSE
                 wqauxr(iq)=wq0
                 wqaux(iq)=1
              ENDIF
!
!   search for a possible optional letter
!
              DO k=j,256
                 IF (ICHAR(input_line(k:k)) == 39) THEN
                    letter_path(iq)(1:1) = input_line(k+1:k+1)
                    IF (ICHAR(input_line(k+2:k+2)) /= 39) &
                       letter_path(iq)(2:2) = input_line(k+2:k+2)
                    IF (ICHAR(input_line(k+3:k+3)) /= 39) &
                       letter_path(iq)(3:3) = input_line(k+3:k+3)
                    EXIT
                 ENDIF
              ENDDO
              EXIT
           ELSEIF ((ICHAR(input_line(j:j)) < 123 .AND. &
                    ICHAR(input_line(j:j)) > 64))  THEN
!
!   This is a letter, not a space character. We read the next three 
!   characters and save them in the letter array, save also which q point
!   it is
!
              npk_label=npk_label+1
              READ(input_line(j:),'(a3)') letter(npk_label)
              label_list(npk_label)=iq
              letter_path(iq)=letter(npk_label)
!
!  now we remove the letters from input_line and read the number of points
!  of the line. The next two line should account for the case in which
!  there is only one space between the letter and the number of points.
!
              nch=3
              IF ( ICHAR(input_line(j+1:j+1))==32 .OR. &
                   ICHAR(input_line(j+2:j+2))==32 ) nch=2
              buffer=input_line(j+nch:)
              IF (sym_divide) THEN
                 READ(buffer,*,err=50,iostat=ios) wqaux(iq), &
                                      nrap_plot_in(iq)

                 IF (nrap_plot_in(iq)>0) THEN
                    nrp=nrap_plot_in(iq)
                    READ(buffer,*,err=50,iostat=ios) wqaux(iq), &
                          nrap_plot_in(iq), (rap_plot_in(i,iq), i=1,nrp)
                 END IF
              ELSE
                 READ(buffer,*,err=50,iostat=ios) wqaux(iq)
              ENDIF
50            IF (ios /=0) CALL errore('thermo_readin',&
                                     'problem reading number of points',1)
              EXIT
           ENDIF
        ENDDO
     ENDDO
  ENDIF
70  CONTINUE
  IF (meta_ionode) CLOSE( UNIT = iun_thermo, STATUS = 'KEEP' )
  parse_unit=parse_unit_save
  !
  !  Then open an input for each image and copy the input file
  !
  ALLOCATE(input(nimage))
  ALLOCATE(iun_image(nimage))

  DO image=1,nimage
     input(image)='_temporary_'//TRIM(int_to_char(image))
     iun_image(image)=100+image
  END DO

  IF (meta_ionode) THEN
     DO image=1,nimage
        OPEN(UNIT=iun_image(image),FILE=TRIM(input(image)),STATUS='UNKNOWN', &
             FORM='FORMATTED', ERR=30, IOSTAT=ios )
     ENDDO
     dummy=' '
     DO WHILE ( .TRUE. )
        READ (5,fmt='(A512)',END=20) dummy
        DO image=1,nimage
           WRITE (iun_image(image),'(A)') TRIM(dummy)
        ENDDO
     END DO
     !
20   DO image=1,nimage 
        CLOSE ( UNIT=iun_image(image), STATUS='KEEP' )
     END DO
  END IF 
30  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  !
  !  Read the input of pw.x and copy the result in the pw.x variables
  !  Note that each image reads its input file
  !
  CALL read_input_file('PW',TRIM(input(my_image_id+1)))
  outdir_thermo=outdir
  CALL iosys()
  at_save = at
  celldm_save=celldm
  ALLOCATE(tau_save(3,nat))
  tau_save=tau
!
!  bring tau_save in crystal coordinates. In strained geometries tau_save
!  is kept constant.
!
  CALL cryst_to_cart( nat, tau_save, bg, -1 )
  ibrav_save=ibrav
  nosym_save=nosym
  input_file_=input(my_image_id+1)
  IF (ionode) THEN
     INQUIRE( FILE=TRIM(input(my_image_id+1)), EXIST = exst )
     IF (exst) THEN
        OPEN(UNIT=iun_image(my_image_id+1),FILE=TRIM(input(my_image_id+1)), &
               STATUS='UNKNOWN', FORM='FORMATTED', ERR=40, IOSTAT=ios )
        CLOSE( UNIT = iun_image(my_image_id+1), STATUS = 'DELETE' )
40      CONTINUE
     ENDIF
  ENDIF

  IF (ibrav > 3 .AND. what=='mur_lc_t') CALL errore('thermo_pw','option not &
              & available for noncubic systems',1)

  IF (what=='piezoelectric_tensor' .OR. what=='mur_lc_piezoelectric_tensor') &
                                                                     THEN
     IF (.NOT.frozen_ions .AND. (forc_conv_thr > 5.d-5)) THEN
        WRITE(stdout,'(/,5x,"Force_conv_thr is too large for computing the &
                          &piezoelectric tensor ")')
        WRITE(stdout,'(5x,"5.d-5 or lower is required")')
        CALL errore('thermo_readin','Force_conv_thr too large',1) 
     ENDIF
  ENDIF

  IF (what=='scf_2d_bands'.AND.identify_sur) THEN
     IF (sur_layers==0) THEN
        sur_layers=MIN(2, nat/2)
     ENDIF
  ENDIF
     

     
  DEALLOCATE(input)
  DEALLOCATE(iun_image)

  IF (set_internal_path) CALL set_bz_path()
  IF (set_2d_path) CALL set_2d_bz_path()

  RETURN
END SUBROUTINE thermo_readin

!-----------------------------------------------------------------------
SUBROUTINE thermo_ph_readin()
  !-----------------------------------------------------------------------
  !
  !  This routine reads the input of the ph.x code when the calculation
  !  requires the use of the phonon routines.
  !
  USE thermo_mod, ONLY : what
  USE mp_world,   ONLY : world_comm
  USE io_global,  ONLY : meta_ionode, meta_ionode_id
  USE io_files,   ONLY : outdir_in_ph => outdir
  USE mp,         ONLY : mp_bcast
  USE output, ONLY : fildyn
  USE command_line_options, ONLY : input_file_ 
  USE input_parameters, ONLY : outdir
  !
  IMPLICIT NONE
  INTEGER :: ios
  CHARACTER(LEN=512) :: dummy
  !
  !  Only the meta_io_node reads the input and sends it to all images
  !  the input is read from file input_file_. This routine searches the
  !  string '---'. It is assumed that the input of the ph.x code is
  !  written after this string.
  !
  IF ( what == 'scf_ph' .OR. what== 'scf_disp' .OR. what == 'mur_lc_ph' &
     .OR. what== 'mur_lc_disp' .OR. what == 'mur_lc_t') THEN

     IF (meta_ionode) OPEN(unit=5, FILE='ph_control', STATUS='OLD', &
                 FORM='FORMATTED', ERR=20, IOSTAT=ios )
20   CALL mp_bcast(ios, meta_ionode_id, world_comm)
     CALL errore('thermo_ph_readin','error opening file '//'ph_control',&
                                                                     ABS(ios))
!
!    be sure that pw variables are completely deallocated
!
     CALL clean_pw(.TRUE.)
!
!    save the outdir. NB phq_readin can change the outdir, and it is
!    the responsability of the user to give the correct directory. Ideally
!    outdir should not be given in the input of thermo_pw.
!
     outdir_in_ph=TRIM(outdir)
     CALL phq_readin()
     IF (meta_ionode) CLOSE(5,STATUS='KEEP')
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE thermo_ph_readin

SUBROUTINE summarize_kpt(xqaux, wqaux, nqaux, letter_path )
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nqaux
  REAL(DP), INTENT(IN) :: xqaux(3,nqaux)
  INTEGER, INTENT(IN) :: wqaux(nqaux)
  CHARACTER(LEN=3), INTENT(IN) :: letter_path(nqaux)

  INTEGER :: ik
  
  WRITE(stdout,'(/,5x, "k points coordinates (2 pi / alat)")') 
  DO ik=1, nqaux
     WRITE(stdout, '(a3, 3f15.8,i5)') letter_path(ik), xqaux(:,ik), &
                                      wqaux(ik)
  ENDDO
  WRITE(stdout,*)
  WRITE(stdout,'(5x, "Input path: ")') 
  WRITE(stdout,'(i5)') nqaux
  DO ik=1, nqaux
     IF (letter_path(ik) == '   ') THEN
         WRITE(stdout, '(3x, 3f15.8,i5)')  xqaux(:,ik), wqaux(ik)
     ELSE
         WRITE(stdout, '(a3, i8)') letter_path(ik), wqaux(ik)
     ENDIF
  ENDDO
  WRITE(stdout,*)

RETURN
END SUBROUTINE summarize_kpt
