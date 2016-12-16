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
  USE constants,            ONLY : k_boltzmann_ry
  USE klist,                ONLY : degauss
  USE ktetra,               ONLY : ltetra
  USE thermo_mod,           ONLY : what, ngeo, step_ngeo, reduced_grid, &
                                   fact_ngeo, max_geometries, start_geo, &
                                   jump_geo
  USE control_thermo,       ONLY : outdir_thermo, after_disp, with_eigen,  &
                                   do_scf_relax, ltherm_dos, ltherm_freq,  &
                                   continue_zero_ibrav, find_ibrav
  USE data_files,           ONLY : flevdat, flfrc, flfrq, fldos, fltherm,  &
                                   flanhar, filband, flkeconv,             &
                                   flenergy, flpbs, flprojlayer,           &
                                   flnkconv, flgrun, flpgrun, fl_el_cons,  &
                                   flpband, flvec, flepsilon, fleldos,     &
                                   fleltherm, fldosfrq
  USE temperature,          ONLY : tmin, tmax, deltat, ntemp
  USE control_pressure,     ONLY : pressure
  USE ifc,                  ONLY : zasr
  USE control_dosq,         ONLY : nq1_d, nq2_d, nq3_d, ndos_input, deltafreq, &
                                   freqmin_input, freqmax_input, &
                                   phdos_sigma
  USE input_parameters,     ONLY : outdir,ibrav, forc_conv_thr, max_seconds, &
                                   calculation
  USE input_parameters, ONLY : a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units
  USE read_input,           ONLY : read_input_file
  USE command_line_options, ONLY : input_file_ 
  USE control_paths,        ONLY : xqaux, wqaux, wqauxr, npk_label, letter, &
                                   label_list, nqaux, q_in_band_form, &
                                   q_in_cryst_coord, q2d, point_label_type, &
                                   disp_q, disp_nqs, npx, &
                                   letter_path, nrap_plot_in, &
                                   label_disp_q, rap_plot_in, long_path
  USE control_gnuplot,      ONLY : flgnuplot, gnuplot_command, lgnuplot
  USE postscript_files,     ONLY : flpsband, flpsdisp, flpsmur, flpsdos, &
                                   flpstherm, flpsanhar, flpskeconv, &
                                   flpsnkconv, flpsgrun,  flpsenergy, &
                                   flpsepsilon, flpseldos, flpseltherm
  USE control_2d_bands,     ONLY : lprojpbs, nkz, gap_thr, sym_divide, &
                                   identify_sur, sur_layers, sur_thr, &
                                   force_bands, only_bands_plot, dump_states, &
                                   subtract_vacuum
  USE control_asy,          ONLY : flasy, lasymptote, asymptote_command
  USE control_bands,        ONLY : emin_input, emax_input, nbnd_bands, lsym, &
                                   enhance_plot
  USE control_eldos,        ONLY : deltae, ndose, nk1_d, nk2_d, nk3_d, &
                                   k1_d, k2_d, k3_d, sigmae, legauss
  USE control_grun,         ONLY : grunmin_input, grunmax_input, &
                                   temp_ph, volume_ph, celldm_ph, lv0_t, &
                                   lb0_t
  USE control_conv,         ONLY : nke, deltake, nkeden, deltakeden, &
                                   nnk, deltank, nsigma, deltasigma
  USE control_mur,          ONLY : vmin_input, vmax_input, deltav, nvol, &
                                   lmurn, celldm0 
  USE control_elastic_constants, ONLY : at_save, tau_save, delta_epsilon, &
                                        ngeo_strain, frozen_ions, &
                                        elastic_algorithm, poly_degree, &
                                        elcpvar, epsilon_0
  USE control_piezoelectric_tensor, ONLY : nosym_save
  USE control_xrdp,         ONLY : lambda, flxrdp, flpsxrdp, lformf, smin, &
                                   smax, nspoint, flformf, flpsformf, lcm, &
                                   lxrdp, lambda_elem
  USE xrdp_module,          ONLY : select_lambda
  USE control_energy_plot,  ONLY : ncontours, color_levels, ene_levels 
  USE control_quadratic_energy, ONLY : show_fit
  USE control_quartic_energy, ONLY : lquartic, lquartic_ph, lsolve
  USE piezoelectric_tensor, ONLY : nppl
  USE control_pwrun,        ONLY : celldm_save, ibrav_save, ityp_save, &
                                   amass_save
  USE lattices,             ONLY : find_ibrav_code
  USE control_ph,           ONLY : xmldyn
  USE output,               ONLY : fildyn
  USE cell_base,            ONLY : at, bg, celldm, cell_base_init
  USE ions_base,            ONLY : nat, tau, ntyp => nsp, ityp, amass, atm, &
                                   if_pos
  USE symm_base,            ONLY : nosym
  USE mp_world,             ONLY : world_comm
  USE mp_images,            ONLY : nimage, my_image_id, root_image
  USE parser,               ONLY : read_line, parse_unit
  USE environment,          ONLY : environment_end
  USE wrappers,             ONLY : f_mkdir_safe
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
  INTEGER :: iq, ipol, i, j, k, igeo, icont, ia, jpol
  INTEGER :: iun_thermo, parse_unit_save

  INTEGER :: nch, nrp, ierr
  LOGICAL :: tend, terr, read_paths, set_internal_path, set_2d_path, exst
  LOGICAL :: has_xml
  CHARACTER(LEN=256) :: input_line, buffer
  REAL(DP) :: wq0, max_seconds_
  REAL(DP) :: ur(3,3), global_s(3,3), omega, alat_save
  REAL(DP), PARAMETER :: eps1=1D-8
  !
  NAMELIST / input_thermo / what, ngeo, zasr,               &
                            flfrc, flfrq, fldos, fltherm,   &
                            flanhar, filband, flkeconv,     &
                            flnkconv,                       &
                            fleldos, fleltherm,             &
                            fldosfrq,                       &
                            fact_ngeo,                      &
                            step_ngeo,                      &
                            start_geo,                      &
                            jump_geo,                       &
                            max_geometries,                 &
                            lquartic,                       &
                            lquartic_ph,                    &
                            lsolve,                         &
                            show_fit,                       &
                            max_seconds,                    &
                            continue_zero_ibrav,            &
                            find_ibrav,                     &
                            nq1_d, nq2_d, nq3_d,            &
                            nk1_d, nk2_d, nk3_d,            &
                            k1_d, k2_d, k3_d,               &
                            deltae, ndose,                  &
                            sigmae, legauss,                &
                            tmin, tmax, deltat, ntemp,      &
                            pressure,                       &
                            freqmin_input, freqmax_input,   &
                            ndos_input, deltafreq,          &
                            phdos_sigma,                    &
                            q2d, q_in_band_form,            &
                            q_in_cryst_coord,               &
                            point_label_type,               &
                            nbnd_bands,                     &
                            nppl,                           &
                            lsym,                           &
                            enhance_plot,                   &
                            long_path,                      &
                            npx,                            &
                            lprojpbs, nkz, gap_thr,         &
                            reduced_grid,                   &
                            only_bands_plot,                &
                            sur_thr, sur_layers,            &
                            sym_divide, identify_sur,       &
                            subtract_vacuum,                &
                            force_bands,                    &
                            dump_states,                    &
                            ncontours,                      &
                            temp_ph, volume_ph, celldm_ph,  &
                            after_disp,                     &
                            with_eigen,                     &
                            do_scf_relax,                   &
                            ltherm_dos, ltherm_freq,        &
                            fildyn,                         &
                            flevdat,                        &
                            flpband, flpgrun,               &
                            flgnuplot, flpsband,            &
                            flpsdisp, flpsdos, flpstherm,   &
                            flpsanhar, flpsmur, flpskeconv, &
                            flpsnkconv, flgrun,             &
                            flpsenergy,                     &
                            flpsdos, flpseltherm,           &
                            flenergy,                       &
                            fl_el_cons,                     &
                            flvec,                          &
                            flxrdp,                         &
                            flasy, asymptote_command,       &
                            lasymptote,                     &
                            flpsxrdp, flpsxrdp,             &
                            flpsformf, flformf,             &
                            smin, smax, nspoint,            &
                            lformf, lcm, lxrdp,             &
                            emin_input, emax_input,         &
                            vmin_input, vmax_input, deltav, &
                            lmurn,                          &
                            elastic_algorithm,              &
                            poly_degree,                    &
                            grunmin_input, grunmax_input,   &
                            lv0_t, lb0_t,                   &
                            delta_epsilon, ngeo_strain,     &
                            epsilon_0,                      &
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
  max_seconds_=1.D8
  IF (meta_ionode) &
     OPEN(UNIT=iun_thermo,FILE='thermo_control',STATUS='OLD', &
                               FORM='FORMATTED', ERR=10, IOSTAT=ios )
10  CALL mp_bcast(ios, meta_ionode_id, world_comm )
    CALL errore( 'thermo_readin', 'opening thermo_control file', ABS( ios ) )

  what=' '
  ngeo=0
  fact_ngeo=1
  step_ngeo(1) = 0.05_DP
  step_ngeo(2) = 0.02_DP
  step_ngeo(3) = 0.02_DP
  step_ngeo(4) = 0.5_DP
  step_ngeo(5) = 0.5_DP
  step_ngeo(6) = 0.5_DP
  reduced_grid =.FALSE.
  start_geo=1
  jump_geo=1
  max_geometries=1000000
  lsolve=2
  lquartic=.TRUE.
  lquartic_ph=.FALSE.
  show_fit=.FALSE.
  max_seconds=1.D8

  nq1_d=192
  nq2_d=192
  nq3_d=192
  zasr='simple'
  freqmin_input=0.0_DP
  freqmax_input=0.0_DP
  deltafreq=1.0_DP
  phdos_sigma=2.0_DP
  ndos_input=1

  tmin=1.0_DP
  tmax=800.0_DP
  deltat=3.0_DP
  ntemp=1
  
  pressure=0.0_DP

  continue_zero_ibrav=.FALSE.
  find_ibrav=.FALSE.

  nbnd_bands=0
  emin_input=0.0_DP
  emax_input=0.0_DP
  long_path=.TRUE.
  lsym=.TRUE.
  enhance_plot=.FALSE.

  vmin_input=0.0_DP
  vmax_input=0.0_DP
  deltav=0.0_DP
  nvol=1
  lmurn=.TRUE.

  after_disp=.FALSE.
  with_eigen=.FALSE.
  fildyn=' '

  nke=5
  deltake=10.0_DP
  nkeden=1
  deltakeden=100.0_DP

  nnk=5
  deltank=2 
  nsigma=1  
  deltasigma=0.005_DP

  nk1_d=16
  nk2_d=16
  nk3_d=16
  k1_d=1
  k2_d=1
  k3_d=1
  deltae=0.0_DP
  ndose=0
  sigmae=0.0_DP
  legauss=.FALSE.

  delta_epsilon=0.005_DP
  epsilon_0=0.0_DP
  ngeo_strain=0
  frozen_ions=.FALSE.
  elastic_algorithm='standard'
  poly_degree=0

  nppl=51
  ncontours=0

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
  do_scf_relax=.FALSE.
  ltherm_dos=.TRUE.
  ltherm_freq=.TRUE.
  dump_states=.FALSE.

  grunmin_input=0.0_DP
  grunmax_input=0.0_DP
  temp_ph=0.0_DP
  volume_ph=0.0_DP
  celldm_ph=0.0_DP
  lv0_t=.TRUE.
  lb0_t=.TRUE.

  lambda=0.0_DP
  lambda_elem=' '
  smin=0.0_DP
  smax=1.0_DP
  nspoint=200
  lformf=.FALSE.
  lcm=.FALSE.
  lxrdp=.FALSE.

  filband='output_band.dat'
  flpband='output_pband.dat'
  flpgrun='output_pgrun.dat'
  flfrc='output_frc.dat'
  flfrq='output_frq.dat'
  fldosfrq='save_frequencies.dat'
  fldos='output_dos.dat'
  fleldos='output_eldos.dat'
  flvec='matdyn.modes'
  flkeconv='output_keconv.dat'
  flnkconv='output_nkconv.dat'
  fltherm='output_therm.dat'
  fleltherm='output_eltherm.dat'
  flanhar='output_anhar.dat'
  flgrun='output_grun.dat'
  flevdat='output_ev.dat'
  fl_el_cons='output_el_cons.dat'
  flxrdp='output_xrdp.dat'
  flformf='output_formf.dat'
  flpbs='output_pbs'
  flprojlayer='output_projlayer'
  flenergy='output_energy'
  flepsilon='epsilon'

  flgnuplot='gnuplot.tmp'
  flpsmur='output_mur.ps'
  flpsband='output_band.ps'
  flpsdisp='output_disp.ps'
  flpsdos='output_dos.ps'
  flpseldos='output_eldos.ps'
  flpstherm='output_therm.ps'
  flpseltherm='output_eltherm.ps'
  flpsanhar='output_anhar.ps'
  flpskeconv='output_keconv.ps'
  flpsnkconv='output_nkconv.ps'
  flpsgrun='output_grun.ps'
  flpsenergy='output_energy.ps'
  flpsepsilon='output_epsilon.ps'
  flpsxrdp='output_xrdp.ps'
  flpsformf='output_formf.ps'

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
  IF (nimage==1) max_seconds_=max_seconds
  !
  IF (what /= 'scf_2d_bands') THEN
     nkz=1
     lprojpbs=.FALSE.
  ENDIF

  IF (lmurn.AND.reduced_grid) CALL errore('thermo_readin',&
                             'lmurn and reduced_grid cannot be both .TRUE.',1)

  IF (after_disp) xmldyn=has_xml(fildyn)

  read_paths=( what=='scf_bands'.OR.what=='scf_disp'.OR.what=='plot_bz'.OR. &
               what=='mur_lc_bands' .OR. what=='mur_lc_disp' .OR. &
               what=='mur_lc_t' .OR. what=='scf_2d_bands')
!
!   ngeo_strain cannot be too small and the energy algorithm requires a
!   few more points
!
  IF (ngeo_strain<4) THEN
     ngeo_strain=4
     IF (elastic_algorithm=='energy') ngeo_strain=6
  ENDIF
!
!   The the default of the interpolation polynomial for elastic constants, if
!   not set in input, or if the input value is unreasonable
!
  IF (poly_degree < 2 ) THEN
     poly_degree = 3
     IF (elastic_algorithm=='energy') poly_degree=4
     IF (ngeo_strain < 6) THEN
       poly_degree = 2
       IF (elastic_algorithm=='energy') poly_degree=3
     ENDIF
  ENDIF
  elcpvar=poly_degree+1
  IF (ngeo_strain < elcpvar) CALL errore('thermo_readin','ngeo_strain is too small',1)


  IF (what(1:6)=='mur_lc') THEN
     IF (ncontours==0) THEN
        ncontours=9
        ALLOCATE(ene_levels(ncontours))
        ALLOCATE(color_levels(ncontours))
        ene_levels=-1000.0_DP
        color_levels='color_black'
     ELSE
        ALLOCATE(ene_levels(ncontours))
        ALLOCATE(color_levels(ncontours))
        ene_levels=-1000.0_DP
        color_levels='color_black'
        IF (meta_ionode) THEN
           DO icont=1,ncontours
              READ (iun_thermo, *, err=400, iostat = ios) ene_levels(icont), &
                                                          color_levels(icont)
           ENDDO
400        CONTINUE
        ENDIF
        CALL mp_bcast(ene_levels, meta_ionode_id, world_comm)
        CALL mp_bcast(color_levels, meta_ionode_id, world_comm)
     END IF
  END IF

  IF (volume_ph==0.0_DP.AND.celldm_ph(1)==0.0_DP.AND.temp_ph==0.0_DP) &
                                         temp_ph=tmin

  IF (lambda==0.0_DP) CALL select_lambda(lambda_elem,lambda)

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
  max_seconds=max_seconds_
  ALLOCATE(tau_save(3,nat))
  IF (ibrav==0.AND..NOT.continue_zero_ibrav) THEN
     alat_save=celldm(1)
     at=at*alat_save
     CALL find_ibrav_code(at(1,1),at(1,2),at(1,3),ibrav,celldm, &
                          ur,global_s,.FALSE.)
     a=0.0_DP
     b=0.0_DP
     c=0.0_DP
     cosab=0.0_DP
     cosac=0.0_DP
     cosbc=0.0_DP
     trd_ht=.FALSE.
     rd_ht=0.0_DP
     cell_units='alat'

     CALL cell_base_init(ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                      trd_ht, rd_ht, cell_units )

     tau_save=tau*alat_save
     tau=0.0_DP
     DO ia=1,nat
        DO ipol=1,3
           DO jpol=1,3
              tau(ipol,ia)=tau(ipol,ia) + global_s(jpol,ipol)*tau_save(jpol,ia)
           ENDDO
        ENDDO
     ENDDO
     tau=tau/celldm(1)

     WRITE(stdout,'(/,5x,"ibrav=0, please use:")')  
     WRITE(stdout,'(/,5x,"ibrav=",i3,",")') ibrav
     WRITE(stdout,'(5x,"celldm(1)= ", f15.10,",")') celldm(1)
     IF (ABS(celldm(2))>eps1) &
        WRITE(stdout,'(5x,"celldm(2)= ", f15.10,",")') celldm(2)
     IF (ABS(celldm(3))>eps1) &
        WRITE(stdout,'(5x,"celldm(3)= ", f15.10,",")') celldm(3)
     IF (ABS(celldm(4))>eps1) &
        WRITE(stdout,'(5x,"celldm(4)= ", f15.10,",")') celldm(4)
     IF (ABS(celldm(5))>eps1) &
        WRITE(stdout,'(5x,"celldm(5)= ", f15.10,",")') celldm(5)
     IF (ABS(celldm(6))>eps1) &
        WRITE(stdout,'(5x,"celldm(6)= ", f15.10,",")') celldm(6)

     WRITE(stdout,'(/,"ATOMIC COORDINATES (alat)")')
     DO ia=1,nat
        IF (if_pos(1,ia) /= 1 .OR. if_pos(2,ia) /= 1 .OR. if_pos(3,ia) /=1) &
           THEN
           WRITE(stdout,'(a3,3f17.10,3i5)') atm(ityp(ia)), &
                         (tau(ipol,ia),ipol=1,3), if_pos(:,ia)
        ELSE
           WRITE(stdout,'(a3,3f17.10)') atm(ityp(ia)), (tau(ipol,ia), ipol=1,3)
        END IF
     ENDDO
 

     IF (.NOT. find_ibrav) THEN
        WRITE(stdout,'(/,5x,"The code will now stop, modify the pw.x input")')
        WRITE(stdout,'(5x,"or set find_ibrav=.TRUE. in thermo_control to &
                                                   &continue ")')
        WRITE(stdout,'(5x,"with these modified coordinates.")')
        WRITE(stdout,'(/,5x,"Set continue_zero_ibrav=.TRUE. to continue &
                                     &with ibrav=0 (not recommended).")')

        CALL environment_end( 'THERMO_PW' )
        !
        CALL mp_global_end ()
        CALL do_stop( 0 )
     ENDIF
  ENDIF
  at_save = at
  celldm_save=celldm
  celldm0=celldm
  ALLOCATE(ityp_save(nat))
  ALLOCATE(amass_save(ntyp))
  ityp_save(:)=ityp(:)
  tau_save=tau
!
!  bring tau_save in crystal coordinates. In strained geometries tau_save
!  is kept constant.
!
  CALL cryst_to_cart( nat, tau_save, bg, -1 )
  ibrav_save=ibrav
  nosym_save=nosym
  input_file_=input(my_image_id+1)
  with_eigen=with_eigen.AND.(ibrav==1.OR.ibrav==2.OR.ibrav==3)


  IF (ionode) THEN
     INQUIRE( FILE=TRIM(input(my_image_id+1)), EXIST = exst )
     IF (exst) THEN
        OPEN(UNIT=iun_image(my_image_id+1),FILE=TRIM(input(my_image_id+1)), &
               STATUS='UNKNOWN', FORM='FORMATTED', ERR=40, IOSTAT=ios )
        CLOSE( UNIT = iun_image(my_image_id+1), STATUS = 'DELETE' )
40      CONTINUE
     ENDIF
  ENDIF

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

  IF (what(1:6)=='mur_lc'.AND..NOT.lmurn) THEN
     IF (calculation/='scf'.AND.calculation/='relax') &
        CALL errore('thermo_readin','thermo_pw requires scf or relax in &
                                               &pw input',1)
  ENDIF

  IF (what(1:7)=='scf_dos'.OR.what(1:10)=='mur_lc_dos') THEN
     IF (deltae==0.0_DP.AND.ndose==0.AND.(degauss>0.0_DP.OR.ltetra)) THEN
        deltae=1.5_DP * k_boltzmann_ry*MIN(4.0_DP, tmin) 
        WRITE(stdout,'(/,5x,"Deltae set to",f20.9," Ry")') deltae
     ENDIF
  END IF

  IF (what=='elastic_constants_t' .AND. elastic_algorithm/='standard') &
     CALL errore('thermo_readin','Only the standard algorithm is working &
                                          &in this case',1)

  IF ( ngeo(1)==0 ) THEN
     IF (what(1:4) == 'scf_') ngeo=1
     IF (what(1:6) == 'mur_lc'.OR.what=='elastic_constants_t') THEN
        IF (lmurn) THEN
           ngeo(1)=9
           DO igeo=2,6
              IF (ngeo(igeo)==0) ngeo(igeo)=1
           ENDDO
        ELSE
!
!   The default mesh is 5 in each crystallographic relevant direction,
!   except for cubic systems for which we take 9 point to make the
!   calculation compatible with lmurn=.TRUE.
!
           ngeo=5
           IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) ngeo(1)=9
        ENDIF
     ENDIF
  END IF
     
  CALL clean_ngeo(ngeo,fact_ngeo,ibrav)
     
  DEALLOCATE(input)
  DEALLOCATE(iun_image)

  IF (set_internal_path) CALL set_bz_path()
  IF (set_2d_path) CALL set_2d_bz_path()
!
!  create the restart directory
!
  ios=0
  IF (meta_ionode) ios = f_mkdir_safe( 'restart' )


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
  USE io_files,   ONLY : outdir_in_ph => tmp_dir
  USE mp,         ONLY : mp_bcast
  USE output, ONLY : fildyn
  USE control_ph, ONLY : ldisp
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
     CALL phq_readin_tpw()
     IF (meta_ionode) CLOSE(5,STATUS='KEEP')
     IF (.NOT.ldisp.AND. what /= 'scf_ph' .AND. what /= 'mur_lc_ph' ) &
        CALL errore('thermo_ph_reading','ldisp should be .TRUE.',1)
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
