
! Copyright (C) 2013-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
PROGRAM thermo_pw
  !-----------------------------------------------------------------------
  !
  ! ... This is a driver for the calculation of thermodynamic quantities,
  ! ... using the harmonic and/or quasiharmonic approximation and the
  ! ... plane waves pseudopotential method.
  ! ... It reads the input of pwscf and an input that specifies
  ! ... which calculations to do and the parameters for these calculations.
  ! ... It checks the scratch directories to see what has been already
  ! ... calculated. The info for the quantities that have been already
  ! ... calculated is read inside the code. The others tasks are scheduled,
  ! ... their priorities determined, and distributed to the image driver.
  ! ... If there are several available images the different tasks are
  ! ... carried out in parallel. This driver can carry out a scf 
  ! ... calculation, a non scf calculation to determine the band structure,
  ! ... or a linear response calculation at a given q and for a given
  ! ... representation. Finally the root image can carry out several
  ! ... post processing tasks. The task currently implemented are:
  ! ... 
  ! ...   scf       : a single scf calculation to determine the total energy.
  ! ...   scf_ke    : many scf calculations at different cut-offs
  ! ...   scf_nk    : many scf calculations at different numbers of k points
  ! ...   scf_bands : a band structure calculation after a scf calcul.
  ! ...   scf_2d_bands : this is as scf_bands, but the cell is assumed to be 
  ! ...                  a slab and the default path is chosen in the 2d 
  ! ...                  Brillouin zone. This option can be used also to 
  ! ...                  calculate the projected bulk band structure.
  ! ...   scf_ph    : a phonon calculation after an scf run
  ! ...   scf_disp  : a phonon dispersion calculation after a scf run
  ! ...   scf_elastic_constants : elastic constants at zero temperature 
  ! ...   scf_piezoelectric_tensor : piezoelectric tensor at zero temperature
  ! ...
  ! ...   mur_lc    : lattice constant via Murnaghan equation or 
  ! ...               quadratic interpolation of the equation of state
  ! ...   mur_lc_bands  : a band structure calculation at the minimum or the
  ! ...               Murnaghan 
  ! ...   mur_lc_ph : a phonon calculation at the minimum of the Murnaghan
  ! ...   mur_lc_disp : a dispersion calculation at the minimum of the
  ! ...               Murnaghan with the possibility to compute the harmonic
  ! ...               thermodynamical quantities
  ! ...   mur_lc_elastic_constants : elastic constants at zero temperature 
  ! ...               at the minimum of the Murnaghan equation 
  ! ...   mur_lc_piezoelectric_tensor : piezoelectric tensor at zero temperature
  ! ...               at the minimum of the Murnaghan equation 
  ! ...
  ! ...   mur_lc_t  : lattice constant and bulk modulus as a function 
  ! ...               of temperature within the quasiharmonic approximation
  ! ...               for cubic systems or crystal parameters as a function
  ! ...               of temperature for tetragonal, hexagonal, trigonal,
  ! ...               and orthorombic systems. 
  ! ...
  USE kinds,            ONLY : DP
  USE check_stop,       ONLY : check_stop_init
  USE mp_global,        ONLY : mp_startup, mp_global_end
  USE mp_images,        ONLY : nimage, nproc_image, my_image_id, root_image
  USE environment,      ONLY : environment_start, environment_end
  USE mp_world,         ONLY : world_comm
  USE mp_asyn,          ONLY : with_asyn_images
  USE control_ph,       ONLY : with_ext_images, always_run, ldisp
  USE io_global,        ONLY : ionode, stdout, meta_ionode_id
  USE mp,               ONLY : mp_sum, mp_bcast
  USE control_thermo,   ONLY : lev_syn_1, lev_syn_2, lpwscf_syn_1,         &
                               lbands_syn_1, lph, outdir_thermo, lq2r,     &
                               lmatdyn, ldos, ltherm,                      &
                               lconv_ke_test, lconv_nk_test,               &
                               spin_component, after_disp, lelastic_const, &
                               lpiezoelectric_tensor, lpolarization, lpart2_pw
  USE postscript_files, ONLY : flpstherm, flpsdisp, flpsdos
  USE data_files,        ONLY : flevdat, fl_el_cons
  USE elastic_constants, ONLY : print_elastic_constants, &
                                compute_elastic_constants, epsilon_geo, &
                                sigma_geo, el_con, el_compliances, &
                                compute_elastic_compliances, &
                                print_elastic_compliances, read_elastic, &
                                write_elastic, print_macro_elasticity, &
                                compute_elastic_constants_adv, &
                                compute_elastic_constants_ene
  USE piezoelectric_tensor, ONLY : compute_piezo_tensor, &
                                compute_d_piezo_tensor, &
                                polar_geo, g_piezo_tensor, d_piezo_tensor, &
                                print_d_piezo_tensor, print_g_piezo_tensor
  USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions, &
                                elastic_algorithm, rot_mat, omega0
  USE internal_files_names,  ONLY : flfrq_thermo, flvec_thermo
  USE control_paths,    ONLY : nqaux
  USE control_gnuplot,  ONLY : flgnuplot
  USE control_bands,    ONLY : nbnd_bands
  USE control_pwrun,    ONLY : ibrav_save, do_punch
  USE thermo_sym,       ONLY : laue, code_group_save
  USE cell_base,        ONLY : omega
  USE wvfct,            ONLY : nbnd
  USE lsda_mod,         ONLY : nspin
  USE thermodynamics,   ONLY : phdos_save
  USE ph_freq_thermodynamics, ONLY: ph_freq_save
  USE temperature,      ONLY : ntemp
  USE phdos_module,     ONLY : destroy_phdos
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units, outdir
  USE control_mur,      ONLY : vmin, b0, b01, emin, celldm0, lmurn
  USE thermo_mod,       ONLY : what, ngeo, omega_geo, energy_geo, &
                               tot_ngeo, reduced_grid, ibrav_geo, celldm_geo, &
                               central_geo
  USE cell_base,        ONLY : ibrav_ => ibrav, celldm_ => celldm
  USE control_2d_bands, ONLY : only_bands_plot
  USE ph_restart,       ONLY : destroy_status_run
  USE save_ph,          ONLY : clean_input_variables
  USE output,           ONLY : fildyn
  USE io_files,         ONLY : tmp_dir, wfc_dir
  USE cell_base,        ONLY : cell_base_init
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, irr, ierr
  CHARACTER (LEN=9)   :: code = 'THERMO_PW'
  CHARACTER (LEN=256) :: auxdyn=' '
  CHARACTER (LEN=256) :: diraux=' '
  CHARACTER(LEN=6) :: int_to_char
  INTEGER :: part, nwork, igeom, itemp, nspin0, exit_status
  LOGICAL  :: exst, parallelfs
  LOGICAL :: check_file_exists, check_dyn_file_exists
  CHARACTER(LEN=256) :: file_dat
  ! Initialize MPI, clocks, print initial messages
  !
  CALL mp_startup ( start_images=.true. )
  CALL environment_start ( code )
  CALL start_clock( 'PWSCF' )
  with_asyn_images=(nimage > 1)
  !
  ! ... and begin with the initialization part
  !
  CALL thermo_readin()
  !
  CALL set_temperature()
  !
  CALL thermo_summary()
  !
  IF (my_image_id /= root_image) CALL check_stop_init()
  !
  part = 1
  !
  CALL initialize_thermo_work(nwork, part)
  file_dat='save_energy' 
  IF ( .NOT. check_file_exists(file_dat) ) THEN
     !
     !  In this part the images work asyncronously. No communication is
     !  allowed except though the master-workers mechanism
     !
     CALL run_thermo_asyncronously(nwork, part, 1, auxdyn)
     !
     !  In this part all images are syncronized and can communicate 
     !  their results thought the world_comm communicator
     !
     CALL mp_sum(energy_geo, world_comm)
     energy_geo=energy_geo / nproc_image
     CALL write_energy(nwork, file_dat)
  ELSE
     CALL read_energy(nwork, file_dat)
  ENDIF
!
!  In the kinetic energy test write the results
!
  IF (lconv_ke_test) THEN
     CALL write_e_ke()
     CALL plot_e_ke()
  ENDIF
!
! In the k-point test write the results
!
  IF (lconv_nk_test) THEN
     CALL write_e_nk()
     CALL plot_e_nk()
  ENDIF
!
!  In a murnaghan equation calculation determine the lattice constant,
!  bulk modulus and its derivative and write the results
!
  IF (lev_syn_1) THEN
     IF (lmurn) THEN
        CALL do_ev()
        CALL write_mur(vmin,b0,b01,emin)
        CALL plot_mur()
        CALL compute_celldm_geo(vmin, celldm0, &
                   celldm_geo(1,central_geo), omega_geo(central_geo))
     ELSE
        CALL write_gnuplot_energy(nwork)
        CALL quadratic_fit()
        IF (.NOT. reduced_grid) CALL plot_multi_energy()
        WRITE(stdout,'(5x,"The minimum energy is obtained for celldm")')
        WRITE(stdout,'(5x,6f12.5)') celldm0
     ENDIF
     CALL mp_bcast(celldm0, meta_ionode_id, world_comm)
     celldm=celldm0
     CALL cell_base_init ( ibrav, celldm0, a, b, c, cosab, cosac, cosbc, &
                      trd_ht, rd_ht, cell_units )
     CALL set_fft_mesh()
     omega0=omega
  END IF

  CALL deallocate_asyn()

  IF (lpwscf_syn_1) THEN
     outdir=TRIM(outdir_thermo)//'g1/'
     tmp_dir = TRIM ( outdir )
     wfc_dir = tmp_dir
     CALL check_tempdir ( tmp_dir, exst, parallelfs )

     IF (my_image_id==root_image) THEN
!
!   do the self consistent calculation at the new lattice constant
!
        do_punch=.TRUE.
        IF (.NOT.only_bands_plot) THEN
           WRITE(stdout,'(/,2x,76("+"))')
           WRITE(6,'(5x,"Doing a self-consistent calculation", i5)') 
           WRITE(stdout,'(2x,76("+"),/)')
           CALL do_pwscf(exit_status, .TRUE.)
        ENDIF
        IF (lbands_syn_1) THEN
!
!   do the band calculation after setting the path
!
           IF (.NOT.only_bands_plot) THEN
              CALL set_paths_disp()
              CALL set_k_points()
              IF (nbnd_bands > nbnd) nbnd = nbnd_bands
              WRITE(stdout,'(/,2x,76("+"))')
              WRITE(6,'(5x,"Doing a non self-consistent calculation", i5)') 
              WRITE(stdout,'(2x,76("+"),/)')
              CALL do_pwscf(exit_status, .FALSE.)
              nspin0=nspin
              IF (nspin==4) nspin0=1
              DO spin_component = 1, nspin0
                 CALL bands_sub()
                 CALL plotband_sub(1,1,' ')
              ENDDO
           ELSE
              CALL read_minimal_info(.TRUE.)
              CALL plotband_sub(1,1,' ')
           ENDIF
        ENDIF
     ENDIF
  END IF
     !
  IF (lpart2_pw) THEN
!
!   here the second part does not use the phonon code. This is for the
!   calculation of elastic constants
!
     part=2
     CALL initialize_thermo_work(nwork, part)
     !
     !  Asyncronous work starts again. No communication is
     !  allowed except though the master workers mechanism
     !
     CALL run_thermo_asyncronously(nwork, part, igeom, auxdyn)
     !
     ! here we return syncronized and calculate the elastic constants 
     ! from energy or stress 
     !

     IF (lelastic_const) THEN
       IF (elastic_algorithm == 'energy') THEN
!
!   save the energy calculated by all images
!
          CALL mp_sum(energy_geo, world_comm)
          energy_geo=energy_geo / nproc_image
       ELSE
!
!   save the stress tensors calculated by all images
!
          CALL mp_sum(sigma_geo, world_comm)
          sigma_geo=sigma_geo / nproc_image
!         CALL write_stress(nwork, file_dat)
      ENDIF
!
!  the elastic constants are calculated here
!
        IF (elastic_algorithm=='standard') THEN
           CALL compute_elastic_constants(sigma_geo, epsilon_geo, nwork, &
                               ngeo_strain, ibrav_save, laue)
        ELSE IF (elastic_algorithm=='advanced') THEN
           CALL compute_elastic_constants_adv(sigma_geo, epsilon_geo, &
                               nwork, ngeo_strain, ibrav_save, laue, rot_mat)
        ELSE IF (elastic_algorithm=='energy') THEN
           CALL compute_elastic_constants_ene(energy_geo, epsilon_geo, &
                               nwork, ngeo_strain, ibrav_save, laue, omega0)
        END IF
        CALL print_elastic_constants(el_con, frozen_ions)
!
!  now compute the elastic compliances and prints them
!
        CALL compute_elastic_compliances(el_con,el_compliances)
        CALL print_elastic_compliances(el_compliances, frozen_ions)
        CALL print_macro_elasticity( ibrav_save, code_group_save, el_con, &
                                     el_compliances)
!
!  save elastic constants and compliances on file
!
        IF (my_image_id==root_image) CALL write_elastic(fl_el_cons)
     ENDIF

     IF (lpiezoelectric_tensor) THEN
        CALL mp_sum(polar_geo, world_comm)
        polar_geo=polar_geo / nproc_image
!
!  the elastic constants are calculated here
!
        CALL compute_piezo_tensor(polar_geo, epsilon_geo, nwork, &
                               ngeo_strain, ibrav_save, code_group_save)
        CALL print_g_piezo_tensor(frozen_ions)

        IF (my_image_id==root_image) CALL read_elastic(fl_el_cons, exst)
        CALL mp_bcast(exst, meta_ionode_id, world_comm)
        IF (exst) THEN
           CALL mp_bcast(el_con, meta_ionode_id, world_comm)
           CALL mp_bcast(el_compliances, meta_ionode_id, world_comm)
           CALL compute_d_piezo_tensor(el_compliances)
           CALL print_d_piezo_tensor(frozen_ions)
        ENDIF
     END IF
     IF (lpolarization) THEN
        CALL mp_sum(polar_geo, world_comm)
        polar_geo=polar_geo / nproc_image
        CALL print_polarization(polar_geo(:,1), .TRUE. )
     ENDIF

     CALL deallocate_asyn()
  ENDIF

  IF (what(1:8) /= 'mur_lc_t') ngeo=1
!
!   This part makes now one or several phonon calculations, using the
!   image feature of this code and running asyncronously the images
!   different geometries are made in sequence. This should be improved,
!   there should be no need to resyncronize after each geometry
!
  IF (lph) THEN
     !
     ! ... reads the phonon input
     !
     with_ext_images=with_asyn_images
     always_run=.TRUE.
     CALL start_clock( 'PHONON' )
     DO igeom=1,tot_ngeo
        write(6,'(/,5x,40("%"))') 
        write(6,'(5x,"Computing geometry ", i5)') igeom
        write(6,'(5x,40("%"),/)') 
        outdir=TRIM(outdir_thermo)//'g'//TRIM(int_to_char(igeom))//'/'
        !
        IF (.NOT. after_disp) CALL thermo_ph_readin()
        IF (after_disp) ldisp=.TRUE.
        CALL set_files_names(igeom)

        IF (after_disp.AND.what=='mur_lc_t') THEN
!
!  The geometry is read by thermo_ph_readin from the output files of pw.x,
!  except in the case where after_disp=.TRUE.. In this case we have to
!  set it here.
!
           ibrav_=ibrav_geo(igeom)
           celldm_(:)=celldm_geo(:,igeom)
           CALL set_bz_path()
        ENDIF
!
!  Set the BZ path for the present geometry
!
        IF (nqaux > 0) CALL set_paths_disp()
        !
        ! ... Checking the status of the calculation and if necessary initialize
        ! ... the q mesh and all the representations
        !
        auxdyn=fildyn

        IF ( .NOT. check_dyn_file_exists(auxdyn)) THEN

           CALL check_initial_status(auxdyn)
           !
           part=2
           CALL initialize_thermo_work(nwork, part)
           !
           !  Asyncronous work starts again. No communication is
           !  allowed except though the master workers mechanism
           !
           CALL run_thermo_asyncronously(nwork, part, igeom, auxdyn)
           !  
           !   return to syncronous work. Collect the work of all images and
           !   writes the dynamical matrix
           !
           CALL collect_everything(auxdyn)
           !
        END IF
        IF (lq2r) THEN
!
!   Compute the interatomic force constants from the dynamical matrices
!   written on file
!
           CALL q2r_sub(auxdyn) 
!
!    compute interpolated dispersions
!
           IF (lmatdyn) THEN
              CALL matdyn_sub(0, igeom)
              CALL plotband_sub(2,igeom,' ')
           ENDIF
!
!    compute phonon dos
!
           IF (lmatdyn.AND.ldos) THEN
!
!    the phonon dos is calculated and the frequencies saved on the ph_freq_save 
!    structure
!
              IF (.NOT.ALLOCATED(phdos_save)) ALLOCATE(phdos_save(tot_ngeo))
              IF (.NOT.ALLOCATED(ph_freq_save)) ALLOCATE(ph_freq_save(tot_ngeo))
              CALL matdyn_sub(1,igeom)
              CALL plot_phdos()
           ENDIF
!
!    computes the thermodynamical properties
!
           IF (ldos.AND.ltherm) THEN
!
!    first from phonon dos
!
              CALL write_thermo(igeom)
              CALL write_thermo_ph(igeom)
              CALL plot_thermo(igeom)
           ENDIF
        ENDIF
        CALL deallocate_asyn()
        IF (.NOT. after_disp) THEN
           CALL clean_pw(.TRUE.)
           CALL close_phq(.FALSE.)
           CALL clean_input_variables()
           CALL destroy_status_run()
           CALL deallocate_part()
        ENDIF
     ENDDO

     CALL restore_files_names()
!
!     Here the Helmholtz free energy at each lattice constant is available.
!     We can write on file the free energy as a function of the volume at
!     any temperature. For each temperature we call the murnaghan equation
!     to fit the data. We save the minimum volume, the bulk modulus and its
!     pressure derivative for each temperature.
!
     IF (lev_syn_2) THEN
        IF (lmurn) THEN
           DO itemp = 1, ntemp
              CALL do_ev_t(itemp)
              CALL do_ev_t_ph(itemp)
           ENDDO
!
!    here we calculate several anharmonic quantities 
!
           CALL write_anharmonic()
           CALL write_ph_freq_anharmonic()
!
!    here we calculate and plot the gruneisen parameters along the given path.
!
           CALL write_gruneisen_band(flfrq_thermo,flvec_thermo)
           CALL plotband_sub(3,1,flfrq_thermo)
           CALL plotband_sub(4,1,flfrq_thermo)
!
!    here we compute the gruneisen parameters on the uniform mesh
!
           CALL compute_gruneisen()
!
!    here we calculate several anharmonic quantities and plot them.
!
           CALL write_grun_anharmonic()
           CALL plot_anhar() 
        ELSE
!
!    Anisotropic solid. Compute only the crystal parameters as a function
!    of temperature and the thermal expansion tensor
!
           DO itemp = 1, ntemp
              CALL quadratic_fit_t(itemp)
           ENDDO
           CALL write_anhar_anis()
!           CALL write_ph_freq_anhar_anis()
           CALL plot_anhar_anis()
!
!    here we calculate and plot the gruneisen parameters along the given path.
!
           CALL write_gruneisen_band_anis(flfrq_thermo,flvec_thermo)
           CALL plotband_sub(4,1,flfrq_thermo)
           CALL plot_gruneisen_band_anis(flfrq_thermo)
        ENDIF
     ENDIF

     IF (lmatdyn.AND.ldos) THEN
        DO igeom=1,tot_ngeo
           CALL destroy_phdos(phdos_save(igeom))
        ENDDO
        DEALLOCATE(phdos_save)
     ENDIF
  ENDIF
  !
  CALL deallocate_thermo()
  !
  CALL environment_end( 'THERMO_PW' )
  !
  CALL mp_global_end ()
  !
  STOP
  !
END PROGRAM thermo_pw

