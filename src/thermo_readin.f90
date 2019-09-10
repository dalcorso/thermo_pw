!
! Copyright (C) 2013-2016 Andrea Dal Corso
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
  !  If the calculation requires a path (see the variable read_path),
  !  and set_internal_path or set_2d_path are .FALSE. this routine allocates
  !  the variables of the path
  !  xqaux,  wqaux, wqauxr, letter, letter_path
  !  label_list, label_disp_q, nrap_plot_in, rap_plot_in
  !
  !  Moreover it create the restart directory
  !
  USE kinds,                ONLY : DP
!
!  variable read by this routine
!
  USE thermo_mod,           ONLY : what, ngeo, step_ngeo, reduced_grid,    &
                                   fact_ngeo, max_geometries,   &
                                   start_geometry, last_geometry, &
                                   ngeo_ph
  USE control_thermo,       ONLY : outdir_thermo, after_disp, with_eigen,  &
                                   do_scf_relax, ltherm_dos, ltherm_freq,  &
                                   continue_zero_ibrav, find_ibrav,        &
                                   set_internal_path, set_2d_path,         &
                                   all_geometries_together, max_seconds_tpw
  USE data_files,           ONLY : flevdat, flfrc, flfrq, fldos, fltherm,  &
                                   flanhar, filband, flkeconv, flenergy,   &
                                   flpbs, flprojlayer, flnkconv, flgrun,   &
                                   flpgrun, fl_el_cons, flpband, flvec,    &
                                   flepsilon, fleldos, fleltherm, fldosfrq
  USE temperature,          ONLY : tmin, tmax, deltat, ntemp
  USE control_pressure,     ONLY : pressure
  USE control_dosq,         ONLY : nq1_d, nq2_d, nq3_d, ndos_input, deltafreq, &
                                   freqmin_input, freqmax_input, phdos_sigma
  USE control_paths,        ONLY : xqaux, wqaux, wqauxr, npk_label, letter, &
                                   label_list, nqaux, q_in_band_form, &
                                   q2d, q_in_cryst_coord, point_label_type, &
                                   disp_q, disp_nqs, npx, &
                                   letter_path, nrap_plot_in, &
                                   label_disp_q, rap_plot_in, long_path, &
                                   old_path, path_fact, is_a_path
  USE control_gnuplot,      ONLY : flgnuplot, gnuplot_command, lgnuplot, flext
  USE postscript_files,     ONLY : flpsband, flpsdisp, flpsmur, flpsdos, &
                                   flpstherm, flpsanhar, flpskeconv, &
                                   flpsnkconv, flpsgrun,  flpsenergy, &
                                   flpsepsilon, flpseldos, flpseltherm
  USE control_2d_bands,     ONLY : lprojpbs, nkz, gap_thr, sym_divide, &
                                   identify_sur, sur_layers, sur_thr, sp_min, &
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
                                   lmurn
  USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain, &
                                   frozen_ions, elastic_algorithm, &
                                   poly_degree, epsilon_0, use_free_energy, &
                                   start_geometry_qha, last_geometry_qha
  USE control_xrdp,         ONLY : lambda, flxrdp, flpsxrdp, lformf, smin, &
                                   smax, nspoint, flformf, flpsformf, lcm, &
                                   lxrdp, lambda_elem
  USE control_energy_plot,  ONLY : ncontours, color_levels, ene_levels 
  USE control_quadratic_energy, ONLY : show_fit
  USE control_quartic_energy, ONLY : lquartic, poly_degree_ph,     &
                                   poly_degree_elc, poly_degree_cv, &
                                   poly_degree_bfact, lsolve
  USE piezoelectric_tensor, ONLY : nppl
  USE grun_anharmonic,      ONLY : poly_degree_grun
  USE images_omega, ONLY : omega_group
  USE control_qe,   ONLY : force_band_calculation, use_ph_images
  USE band_computation, ONLY : sym_for_diago
!
!  the QE variables needed here. max_seconds, zasr, xmldyn and fildyn
!  could be set by this routine 
!
  USE input_parameters,     ONLY : outdir, forc_conv_thr, max_seconds, &
                                   calculation
  USE control_ph,           ONLY : xmldyn
  USE ifc,                  ONLY : zasr
  USE cell_base,            ONLY : alat
  USE output,               ONLY : fildyn
  USE mp_world,             ONLY : world_comm
  USE mp_images,            ONLY : nimage, my_image_id, root_image
  USE parser,               ONLY : read_line, parse_unit
  USE read_input,           ONLY : read_input_file
  USE command_line_options, ONLY : input_file_ 
  USE check_stop,           ONLY : max_seconds_ => max_seconds
  USE io_global,            ONLY : ionode, meta_ionode, meta_ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  REAL(DP) :: wq0, save_max_seconds
  INTEGER, ALLOCATABLE :: iun_image(:)
  INTEGER :: image, iq, ipol, jpol, icont, iun_thermo, parse_unit_save, &
             nch, nrp, i, j, k, ios
  INTEGER :: iun_input
  INTEGER :: find_free_unit
  LOGICAL :: tend, terr, read_paths, exst, has_xml
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=512) :: dummy
  CHARACTER(LEN=256), ALLOCATABLE :: input(:)
  CHARACTER(LEN=256) :: input_line, buffer
  !
  NAMELIST / input_thermo / what,                           &
!
!  gnuplot
!
                            lgnuplot, gnuplot_command,      &
                            flgnuplot, flext,                     &
!
!  temperature and pressure
!
                            tmin, tmax, deltat, ntemp,      &
                            pressure,                       &
!
!   structure
!
                            continue_zero_ibrav,            &
                            find_ibrav,                     &
!
!   ke convergence
!
                            nke, deltake,                   &
                            nkeden, deltakeden,             &
                            flkeconv, flpskeconv,           &
!
!  nk convergence
!
                            nnk, deltank,                   &
                            nsigma, deltasigma,             &
                            flnkconv, flpsnkconv,           &
!
!  scf_bands
!
                            emin_input, emax_input,         &
                            nbnd_bands, lsym,               &
                            only_bands_plot,                &
                            enhance_plot, long_path,        &
                            old_path,                       &
                            path_fact,                      &
                            q2d, is_a_path,                 &
                            q_in_band_form,                 &
                            q_in_cryst_coord,               &
                            point_label_type,               &
                            filband, flpband,               &
                            flpsband,                       &
!
!   scf_2d_bands
!
                            lprojpbs, nkz, gap_thr,         &
                            sym_divide, identify_sur,       &
                            force_bands, dump_states,       &
                            sur_thr, sur_layers,            &
                            sp_min,                         &
                            subtract_vacuum,                &
                            flpbs, flprojlayer,             &
!
!   scf_dos
!
                            deltae, ndose,                  &
                            nk1_d, nk2_d, nk3_d,            &
                            k1_d, k2_d, k3_d,               &
                            sigmae, legauss,                &
                            fleldos, fleltherm,             &
                            flpsdos, flpseltherm,           &
!
!   plot_bz
!
                            lasymptote,                     &
                            flasy, asymptote_command,       &
                            npx,                            &
                            lambda, lambda_elem,            &
                            flxrdp, flpsxrdp,               &          
                            lxrdp, lformf,                  &
                            smin, smax, nspoint,            &
                            lcm, flformf, flpsformf,        &
!
!   scf_ph
!
                            flepsilon,                      &
                            flpsepsilon,                    &
                            force_band_calculation,         &
                            use_ph_images,                  &
                            sym_for_diago,                  &
!
!   scf_disp
!
                            freqmin_input, freqmax_input,   &
                            deltafreq, ndos_input,          &
                            nq1_d, nq2_d, nq3_d,            &
                            phdos_sigma,                    &
                            after_disp, fildyn,             &
                            zasr,                           &
                            ltherm_dos, ltherm_freq,        &
                            flfrc, flfrq, fldos, flvec,     &
                            fldosfrq, fltherm,              &
                            flpsdisp, flpsdos, flpstherm,   &
!
!   scf_elastic_constants
!
                            frozen_ions,                    &
                            ngeo_strain,                    &
                            elastic_algorithm,              &
                            delta_epsilon, epsilon_0,       &
                            poly_degree,                    &
                            fl_el_cons,                     &
!
!   scf_polarization
!
                            nppl,                           &
!
!   mur_lc
!
                            ngeo, step_ngeo,                &
                            lmurn,                          &
                            reduced_grid,                   &
                            show_fit,                       &
                            vmin_input, vmax_input, deltav, &
                            nvol,                           &
                            lquartic, lsolve,               &
                            flevdat,                        &
                            flpsmur,                        &
                            ncontours,                      &
                            flenergy, flpsenergy,           &
!
!   mur_lc_elastic_constants
!
                            do_scf_relax,                   &
!
!   mur_lc_t
!
                            grunmin_input, grunmax_input,   &
                            volume_ph, celldm_ph, temp_ph,  &
                            with_eigen,                     &
                            poly_degree_ph,                 &
                            poly_degree_cv,                 &
                            poly_degree_bfact,              &
                            poly_degree_elc,                &
                            lv0_t, lb0_t,                   &
                            poly_degree_grun,               &
                            flpgrun, flgrun, flpsgrun,      &
                            flanhar, flpsanhar,             &
                            fact_ngeo, ngeo_ph,             &
                            all_geometries_together,        &
!
!   elastic_constants_t
!
                            use_free_energy,                &
                            start_geometry_qha,             &
                            last_geometry_qha,              &
!
!   optical
!
                            omega_group,                    &
!
!   recover features
!
                            max_geometries,                 &
                            start_geometry,                 &
                            last_geometry,                  &
                            max_seconds                    

  !
  !  First read the input of thermo. This input should be in a file
  !  called thermo_control
  !
  parse_unit_save=parse_unit
  IF (ionode) iun_thermo=find_free_unit()
  parse_unit=iun_thermo
  save_max_seconds=1.D8
  IF (meta_ionode) THEN
     OPEN(UNIT=iun_thermo,FILE='thermo_control',STATUS='OLD', &
                               FORM='FORMATTED', ERR=10, IOSTAT=ios )
  ENDIF
10  CALL mp_bcast(ios, meta_ionode_id, world_comm )
    CALL errore( 'thermo_readin', 'opening thermo_control file', ABS( ios ) )

!
!  Default values of the input variables
!
  what=' '

  flgnuplot='gnuplot.tmp'
  lgnuplot=.TRUE.
  gnuplot_command='gnuplot'
  flext='.ps'

  tmin=1.0_DP
  tmax=800.0_DP
  deltat=3.0_DP
  ntemp=1
  pressure=0.0_DP

  continue_zero_ibrav=.FALSE.
  find_ibrav=.FALSE.

  nke=5
  deltake=10.0_DP
  nkeden=1
  deltakeden=100.0_DP
  flkeconv='output_keconv.dat'
  flpskeconv='output_keconv'

  nnk=5
  deltank=2 
  nsigma=1  
  deltasigma=0.005_DP
  flnkconv='output_nkconv.dat'
  flpsnkconv='output_nkconv'

  emin_input=0.0_DP
  emax_input=0.0_DP
  nbnd_bands=0
  only_bands_plot=.FALSE.
  lsym=.TRUE.
  enhance_plot=.FALSE.
  long_path=.TRUE.
  old_path=.FALSE.
  path_fact=1.0_DP
  filband='output_band.dat'
  flpband='output_pband.dat'
  flpsband='output_band'
  q2d=.FALSE.
  is_a_path=.TRUE.
  q_in_band_form=.TRUE.
  q_in_cryst_coord=.FALSE.
  point_label_type='SC'

  lprojpbs = .TRUE.
  nkz = 1 
  gap_thr = 0.1_DP
  sym_divide = .FALSE.
  identify_sur =.FALSE.
  dump_states=.FALSE.
  sur_layers=0
  sur_thr=0.0_DP
  sp_min=0.0_DP
  subtract_vacuum=.TRUE.
  force_bands=.FALSE.
  flpbs='output_pbs'
  flprojlayer='output_projlayer'

  deltae=0.0_DP
  ndose=0
  nk1_d=16
  nk2_d=16
  nk3_d=16
  k1_d=1
  k2_d=1
  k3_d=1
  sigmae=0.0_DP
  legauss=.FALSE.
  fleldos='output_eldos.dat'
  flpseldos='output_eldos'
  fleltherm='output_eltherm.dat'
  flpseltherm='output_eltherm'

  flasy='asy_tmp'
  lasymptote=.FALSE.
  asymptote_command='asy -f pdf -noprc'
  npx=8

  lambda=0.0_DP
  lambda_elem=' '
  flxrdp='output_xrdp.dat'
  flpsxrdp='output_xrdp'
  lxrdp=.FALSE.
  lformf=.FALSE.
  smin=0.0_DP
  smax=1.0_DP
  nspoint=200
  lcm=.FALSE.
  flformf='output_formf.dat'
  flpsformf='output_formf'

  flepsilon='epsilon'
  flpsepsilon='output_epsilon'
  force_band_calculation=.FALSE.
  IF (nimage>1) THEN
     use_ph_images=.FALSE.
  ELSE
     use_ph_images=.TRUE.
  ENDIF
  sym_for_diago=.FALSE.

  freqmin_input=0.0_DP
  freqmax_input=0.0_DP
  deltafreq=1.0_DP
  ndos_input=1
  nq1_d=192
  nq2_d=192
  nq3_d=192
  phdos_sigma=2.0_DP
  after_disp=.FALSE.
  fildyn=' '
  zasr='simple'
  ltherm_dos=.TRUE.
  ltherm_freq=.TRUE.
  flfrc='output_frc.dat'
  flfrq='output_frq.dat'
  flvec='matdyn.modes'
  fldos='output_dos.dat'
  fldosfrq='save_frequencies.dat'
  fltherm='output_therm.dat'
  flpsdisp='output_disp'
  flpsdos='output_dos'
  flpstherm='output_therm'

  frozen_ions=.FALSE.
  ngeo_strain=0
  elastic_algorithm='standard'
  delta_epsilon=0.005_DP
  epsilon_0=0.0_DP
  poly_degree=0
  fl_el_cons='output_el_cons.dat'

  poly_degree_grun=4

  nppl=51

  ngeo=0
  step_ngeo(1) = 0.05_DP
  step_ngeo(2) = 0.02_DP
  step_ngeo(3) = 0.02_DP
  step_ngeo(4) = 0.5_DP
  step_ngeo(5) = 0.5_DP
  step_ngeo(6) = 0.5_DP
  lmurn=.TRUE.
  reduced_grid =.FALSE.
  show_fit=.FALSE.
  vmin_input=0.0_DP
  vmax_input=0.0_DP
  deltav=0.0_DP
  nvol=1
  lquartic=.TRUE.
  lsolve=2
  flevdat='output_ev.dat'
  flpsmur='output_mur'
  ncontours=0
  flenergy='output_energy'
  flpsenergy='output_energy'

  do_scf_relax=.FALSE.

  grunmin_input=0.0_DP
  grunmax_input=0.0_DP
  volume_ph=0.0_DP
  celldm_ph=0.0_DP
  temp_ph=0.0_DP
  with_eigen=.FALSE.
  poly_degree_ph=4
  poly_degree_cv=4
  poly_degree_bfact=4
  poly_degree_elc=4
  lv0_t=.TRUE.
  lb0_t=.TRUE.
  flpgrun='output_pgrun.dat'
  flgrun='output_grun.dat'
  flpsgrun='output_grun'
  flanhar='output_anhar.dat'
  flpsanhar='output_anhar'
  fact_ngeo=1
  ngeo_ph=0
  omega_group=1
  all_geometries_together=.FALSE.

  use_free_energy=.FALSE.
  start_geometry_qha=1
  last_geometry_qha=1000000

  max_geometries=1000000
  start_geometry=1
  last_geometry=1000000
  max_seconds=1.D8


  IF (meta_ionode) READ( iun_thermo, input_thermo, ERR=100, END=100,  &
                                                   IOSTAT = ios )
100 CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore( 'thermo_readin', 'reading input_thermo namelist', ABS( ios ) )
!
  CALL bcast_thermo_input()
!
!   Here a few consistency check on the input variables
!
  IF (what==' ') CALL errore('thermo_readin','''what'' must be initialized',1)

  IF (what/='mur_lc_t'.AND.what/='elastic_constants_t'&
                             .AND.all_geometries_together) &
          CALL errore('thermo_readin','all_geometries_together requires &
                          &mur_lc_t or elastic_constants_t',1)

  IF (flext/='.pdf') flext='.ps'

  IF (max_geometries /= 1000000 .AND.all_geometries_together) &
          CALL errore('thermo_readin','all_geometries_together not compatible &
                                          &with max_geometries',1)
  !
  IF (lmurn.AND.reduced_grid) CALL errore('thermo_readin',&
                             'lmurn and reduced_grid cannot be both .TRUE.',1)
  IF (reduced_grid) THEN
     ltherm_dos=.FALSE.
     ltherm_freq=.FALSE.
     lv0_t=.FALSE.
     lb0_t=.FALSE.
  ENDIF

  IF (poly_degree_ph<1 .OR. poly_degree_ph>4) &
            CALL errore('thermo_readin','poly_degree_ph must be between &
                                                              & 1 and 4',1)
  IF (poly_degree_cv<1 .OR. poly_degree_cv>4) &
            CALL errore('thermo_readin','poly_degree_cv must be between & 
                                                               &1 and 4',1)
  IF (poly_degree_bfact<1 .OR. poly_degree_bfact>4) &
            CALL errore('thermo_readin','poly_degree_bcast must be between &
                                                               &1 and 4',1)
  IF (poly_degree_elc<1 .OR. poly_degree_elc>4) &
            CALL errore('thermo_readin','poly_degree_elc must be between &
                                                               &1 and 4',1)

  IF (what=='elastic_constants_t'.AND.elastic_algorithm/='energy_std' &
      .AND.elastic_algorithm/='energy'.AND.use_free_energy) &
     CALL errore('thermo_readin','Only the energy algorithms are available &
                                          &in this case',1)

  read_paths=( what=='scf_bands'.OR.what=='scf_disp'.OR.what=='plot_bz'.OR. &
               what=='mur_lc_bands' .OR. what=='mur_lc_disp' .OR. &
               what=='mur_lc_t' .OR. what=='scf_2d_bands'.OR. &
               what=='elastic_constants_t')

  IF (nimage==1) save_max_seconds=max_seconds
  max_seconds_tpw=max_seconds
  xmldyn=has_xml(fildyn)
  IF (q2d) is_a_path=.FALSE.
!
!   here read the contour levels
!
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
     ENDIF
  ENDIF
!
!   read the path if given in the input of thermo_pw
!
  nqaux=0
  set_internal_path=.FALSE.
  set_2d_path=.FALSE.
  IF ( read_paths ) THEN
     IF (meta_ionode) READ (iun_thermo, *, END=200, ERR=200, IOSTAT=ios) nqaux
200  CALL mp_bcast(ios, meta_ionode_id, world_comm )
     IF (ios /= 0) THEN 
        IF (what=='scf_2d_bands') THEN
           set_2d_path=.TRUE.
        ELSE
           set_internal_path=.TRUE.
        ENDIF
        GOTO 70
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
                    letter_path(iq)=''
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
                 READ(buffer,*,ERR=50,IOSTAT=ios) wqaux(iq), nrap_plot_in(iq)
                 IF (nrap_plot_in(iq)>0) THEN
                    nrp=nrap_plot_in(iq)
                    READ(buffer,*,ERR=50,IOSTAT=ios) wqaux(iq), &
                          nrap_plot_in(iq), (rap_plot_in(i,iq), i=1,nrp)
                 END IF
              ELSE
                 READ(buffer,*,ERR=50,IOSTAT=ios) wqaux(iq)
              ENDIF
50            IF (ios /=0) CALL errore('thermo_readin',&
                                     'problem reading number of points',1)
!
!   search for a possible optional letter
!
              DO k=j+nch,256
                 IF (ICHAR(input_line(k:k)) == 39) THEN
                    letter_path(iq)=''
                    letter_path(iq)(1:1) = input_line(k+1:k+1)
                    IF (ICHAR(input_line(k+2:k+2)) /= 39) &
                       letter_path(iq)(2:2) = input_line(k+2:k+2)
                    IF (ICHAR(input_line(k+3:k+3)) /= 39) &
                       letter_path(iq)(3:3) = input_line(k+3:k+3)
                    EXIT
                 ENDIF
              ENDDO
              EXIT
           ENDIF
        ENDDO
     ENDDO
  ENDIF
70  CONTINUE
  IF (meta_ionode) CLOSE( UNIT = iun_thermo, STATUS = 'KEEP' )
  parse_unit=parse_unit_save
!
!  Then open an input file for each image and copy there the input file
!
  ALLOCATE(input(nimage))
  ALLOCATE(iun_image(nimage))

  DO image=1,nimage
     input(image)='_temporary_'//TRIM(int_to_char(image))
     iun_image(image)=100+image
  END DO

  IF (meta_ionode) THEN
     IF (input_file_==' ') THEN
        iun_input=5
     ELSE
        iun_input=101+nimage
        OPEN(UNIT=iun_input,FILE=TRIM(input_file_),STATUS='OLD', &
          FORM='FORMATTED', ERR=30, IOSTAT=ios )
     END IF

     DO image=1,nimage
        OPEN(UNIT=iun_image(image),FILE=TRIM(input(image)),STATUS='UNKNOWN', &
             FORM='FORMATTED', ERR=30, IOSTAT=ios )
     ENDDO
     dummy=' '
     DO WHILE ( .TRUE. )
        READ (iun_input,fmt='(A512)',END=20) dummy
        DO image=1,nimage
           WRITE (iun_image(image),'(A)') TRIM(dummy)
        ENDDO
     ENDDO
     !
20   DO image=1,nimage 
        CLOSE ( UNIT=iun_image(image), STATUS='KEEP' )
     ENDDO
     IF (input_file_/=' ') THEN
        CLOSE(UNIT=iun_input, STATUS='KEEP')
     ENDIF
  ENDIF 
30  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore('thermo_readin','Creating input files',ios)
  !
  !  Read the input of pw.x and copy the result in the pw.x variables
  !  Note that each image reads its input file
  !
  CALL read_input_file('PW',TRIM(input(my_image_id+1)))
  outdir_thermo=outdir
  CALL iosys()
  max_seconds_=save_max_seconds

!
!   Now delete the temporary input files
!
  IF (ionode) THEN
     INQUIRE( FILE=TRIM(input(my_image_id+1)), EXIST = exst )
     IF (exst) THEN
        OPEN(UNIT=iun_image(my_image_id+1),FILE=TRIM(input(my_image_id+1)), &
               STATUS='UNKNOWN', FORM='FORMATTED', ERR=40, IOSTAT=ios )
        CLOSE( UNIT = iun_image(my_image_id+1), STATUS = 'DELETE' )
40      CONTINUE
     ENDIF
  ENDIF
!
!  here check the consistency of the input variables of thermo_pw with
!  those of pw.x
!
  IF (what=='piezoelectric_tensor' .OR. what=='mur_lc_piezoelectric_tensor') &
                                                                     THEN
     IF (.NOT.frozen_ions .AND. (forc_conv_thr > 5.d-5)) THEN
        WRITE(stdout,'(/,5x,"Force_conv_thr is too large for computing the &
                          &piezoelectric tensor ")')
        WRITE(stdout,'(5x,"5.d-5 or lower is required")')
        CALL errore('thermo_readin','Force_conv_thr too large',1) 
     ENDIF
  ENDIF

  IF (what(1:6)=='mur_lc'.AND..NOT.lmurn) THEN
     IF (calculation/='scf'.AND.calculation/='relax') &
        CALL errore('thermo_readin','thermo_pw requires scf or relax in &
                                               &pw input',1)
  ENDIF
!
!  default value of sp_min if not given in input
!
  IF (sp_min==0.0_DP) THEN
     sp_min=2.0_DP/alat
  ELSE
     sp_min=sp_min/alat
  ENDIF

  DEALLOCATE(input)
  DEALLOCATE(iun_image)

  RETURN
END SUBROUTINE thermo_readin

!-----------------------------------------------------------------------
SUBROUTINE thermo_ph_readin()
  !-----------------------------------------------------------------------
  !
  !  This routine reads the input of the ph.x code when the calculation
  !  requires the use of the phonon routines. It must be called by all the
  !  CPUs.
  !
  USE thermo_mod, ONLY : what
  USE mp_world,   ONLY : world_comm
  USE io_global,  ONLY : meta_ionode, meta_ionode_id
  USE io_files,   ONLY : outdir_in_ph => tmp_dir
  USE mp,         ONLY : mp_bcast
  USE control_ph, ONLY : ldisp
  USE internal_files_names, ONLY : fildyn_thermo
  USE output,     ONLY : fildyn
  USE input_parameters, ONLY : outdir
  !
  IMPLICIT NONE
  INTEGER :: ios
  CHARACTER(LEN=512) :: dummy
  !
  !  Only the meta_io_node reads the input and sends it to all images
  !  the input is read from file ph_control. This routine searches the
  !  string '---'. It is assumed that the input of the ph.x code is
  !  written after this string.
  !
  IF (meta_ionode) OPEN(UNIT=5, FILE='ph_control', STATUS='OLD', &
                 FORM='FORMATTED', ERR=20, IOSTAT=ios )
20   CALL mp_bcast(ios, meta_ionode_id, world_comm)
     CALL errore('thermo_ph_readin','error opening file '//'ph_control',&
                                                                     ABS(ios))
!
!    be sure that pw variables are completely deallocated
!
  CALL clean_all_pw()
!
!    save the outdir. NB phq_readin can change the outdir, and it is
!    the responsability of the user to give the correct directory. Ideally
!    outdir should not be given in the input of thermo_pw.
!
  outdir_in_ph=TRIM(outdir)
  CALL phq_readin_tpw()
  IF (meta_ionode) CLOSE(5,STATUS='KEEP')
  IF (.NOT.ldisp.AND. what /= 'scf_ph' .AND. what /= 'mur_lc_ph' ) &
        CALL errore('thermo_ph_readin','ldisp should be .TRUE.',1)
  fildyn_thermo="dynamical_matrices/"//TRIM(fildyn)
  !
  RETURN
  !
END SUBROUTINE thermo_ph_readin

