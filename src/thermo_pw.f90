!
! Copyright (C) 2013 Andrea Dal Corso
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
  ! ...   scf_nk    : many scf calculations at different number of k points
  ! ...   scf_bands : a band structure calculation after a scf calcul.
  ! ...   scf_ph    : a phonon at a single q after a scf run
  ! ...   scf_disp  : a phonon dispersion calculation after a scf run
  ! ...   mur_lc    : lattice constant via murnaghan equation
  ! ...   mur_lc_bands  : a band structure calculation at the minimum or the
  ! ...               murnaghan
  ! ...   mur_lc_ph : a phonon calculation at the minimum of the murmaghan
  ! ...   mur_lc_disp : a dispersion calculation at the minimum of the
  ! ...               murnaghan with possibility to compute harmonic
  ! ...               thermodynamical quantities
  ! ...   mur_lc_t  : lattice constant and bulk modulus as a function 
  ! ...               of temperature within the quasiharmonic approximation
  ! ...   elastic_constants : elastic constants at zero temperature 
  ! ...   mur_lc_elastic_constants : elastic constants at zero temperature 
  ! ...               at the minimum of the Murnaghan equation 
  ! ...
  ! ...
  USE kinds,            ONLY : DP
  USE check_stop,       ONLY : check_stop_init
  USE mp_global,        ONLY : mp_startup, mp_global_end
  USE mp_images,        ONLY : nimage, nproc_image, my_image_id, root_image
  USE environment,      ONLY : environment_start, environment_end
  USE mp_world,         ONLY : world_comm
  USE mp_asyn,          ONLY : with_asyn_images
  USE control_ph,       ONLY : with_ext_images, always_run
  USE io_global,        ONLY : ionode, stdout, meta_ionode_id
  USE mp,               ONLY : mp_sum, mp_bcast
  USE control_thermo,   ONLY : lev_syn_1, lev_syn_2, lpwscf_syn_1,         &
                               lbands_syn_1, lph, outdir_thermo, lq2r,     &
                               lmatdyn, ldos, ltherm, flfrc, flfrq, fldos, &
                               fltherm, spin_component, flevdat,           &
                               lconv_ke_test, lconv_nk_test, compute_lc,   &
                               lelastic_const 
  USE ifc,              ONLY : freqmin, freqmax
  USE elastic_constants, ONLY : print_elastic_constants, &
                                compute_elastic_constants, epsilon_geo, &
                                sigma_geo, el_con, el_compliances, &
                                compute_elastic_compliances, &
                                print_elastic_compliances, read_elastic, &
                                write_elastic
  USE control_elastic_constants, ONLY : ibrav_save, ngeo_strain, frozen_ions, &
                                 fl_el_cons
  USE control_paths,    ONLY : nqaux
  USE control_gnuplot,  ONLY : flpsdos, flgnuplot, flpstherm, flpsdisp
  USE control_bands,    ONLY : flpband, nbnd_bands
  USE thermo_sym,       ONLY : laue, code_group_save
  USE wvfct,            ONLY : nbnd
  USE lsda_mod,         ONLY : nspin
  USE thermodynamics,   ONLY : phdos_save
  USE ph_freq_thermodynamics, ONLY: ph_freq_save
  USE temperature,      ONLY : ntemp
  USE phdos_module,     ONLY : destroy_phdos
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units, outdir
  USE control_mur,      ONLY : vmin, b0, b01, emin
  USE thermo_mod,       ONLY : what, ngeo, alat_geo, omega_geo, energy_geo, ntry
  USE ph_restart,       ONLY : destroy_status_run
  USE save_ph,          ONLY : clean_input_variables
  USE output,           ONLY : fildyn
  USE io_files,         ONLY : tmp_dir, wfc_dir
  USE cell_base,        ONLY : cell_base_init
  USE fft_base,         ONLY : dfftp, dffts
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, irr, ierr
  CHARACTER (LEN=9)   :: code = 'THERMO_PW'
  CHARACTER (LEN=256) :: auxdyn=' '
  CHARACTER (LEN=256) :: diraux=' '
  CHARACTER(LEN=6) :: int_to_char
  INTEGER :: part, nwork, igeom, itemp, nspin0, itry, exit_status
  REAL(DP) :: compute_alat_geo, fact
  LOGICAL :: all_done_asyn
  LOGICAL  :: exst, parallelfs
  LOGICAL :: check_file_exists, check_dyn_file_exists
  CHARACTER(LEN=256) :: fildyn_thermo, flfrc_thermo, flfrq_thermo, &
                        fldos_thermo, fltherm_thermo, flpband_thermo, &
                        flpsdos_thermo, flpstherm_thermo, flgnuplot_thermo, &
                        flpsdisp_thermo, file_dat
  !
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
  CALL check_stop_init()
  !
  DO itry = 1, ntry
     part = 1
     !
     CALL initialize_thermo_work(nwork, part)
     IF (.NOT. compute_lc) CYCLE
     !
     !  In this part the images work asyncronously. No communication is
     !  allowed except though the master-workers mechanism
     !
     file_dat='save_energy' // TRIM(int_to_char(itry))
     IF ( .NOT. check_file_exists(file_dat) ) THEN
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
        CALL do_ev(itry)
        CALL mur(vmin,b0,b01,emin)
        CALL plot_mur()
     ENDIF
     !
     CALL deallocate_asyn()
  END DO
  ! 
  !  This part is syncronized, only the first image makes a scf calculation
  !  and a band structure calculation, at the equilibrium lattice constant
  !  found previously.
  !
  !  set the lattice constant at the computed minimum of the Murnaghan
  !
  IF (lev_syn_1) THEN
     celldm(1)=compute_alat_geo(vmin, alat_geo(ngeo/2+1), &
                                              omega_geo(ngeo/2+1))
     CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                      trd_ht, rd_ht, cell_units )
     CALL clean_dfft()
  END IF

  IF (lpwscf_syn_1) THEN
     outdir=TRIM(outdir_thermo)//'g1/'
     tmp_dir = TRIM ( outdir )
     wfc_dir = tmp_dir
     CALL check_tempdir ( tmp_dir, exst, parallelfs )

     IF (my_image_id==root_image) THEN
!
!   do the self consistent calculation at the new lattice constant
!
        CALL do_pwscf(exit_status, .TRUE.)
        IF (lbands_syn_1) THEN
!
!   do the band calculation after setting the path
!
           CALL set_paths_disp()
           CALL set_k_points()
           IF (nbnd_bands > nbnd) nbnd = nbnd_bands
           CALL do_pwscf(exit_status, .FALSE.)
           nspin0=nspin
           IF (nspin==4) nspin0=1
           DO spin_component = 1, nspin0
              CALL bands_sub()
              CALL plotband_sub(1,1,' ')
           ENDDO
        ENDIF
     ENDIF
  END IF
  
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
     DO igeom=1,ngeo
        write(6,'(/,5x,40("%"))') 
        write(6,'(5x,"Computing geometry ", i5)') igeom
        write(6,'(5x,40("%"),/)') 
        outdir=TRIM(outdir_thermo)//'g'//TRIM(int_to_char(igeom))//'/'
        !
        CALL thermo_ph_readin()
        IF (igeom==1) fildyn_thermo=TRIM(fildyn)
        IF (igeom==1) flfrc_thermo=TRIM(flfrc)
        IF (igeom==1) flfrq_thermo=TRIM(flfrq)
        IF (igeom==1) fldos_thermo=TRIM(fldos)
        IF (igeom==1) flpsdos_thermo=TRIM(flpsdos)
        IF (igeom==1) fltherm_thermo=TRIM(fltherm)
        IF (igeom==1) flpstherm_thermo=TRIM(flpstherm)
        IF (igeom==1) flpband_thermo=TRIM(flpband)
        IF (igeom==1) flgnuplot_thermo=TRIM(flgnuplot)
        IF (igeom==1) flpsdisp_thermo=TRIM(flpsdisp)

        fildyn=TRIM(fildyn_thermo)//'.g'//TRIM(int_to_char(igeom))//'.'
        flfrc=TRIM(flfrc_thermo)//'.g'//TRIM(int_to_char(igeom))
        flfrq=TRIM(flfrq_thermo)//'.g'//TRIM(int_to_char(igeom))
        fldos=TRIM(fldos_thermo)//'.g'//TRIM(int_to_char(igeom))
        flpsdos=TRIM(flpsdos_thermo)//'.g'//TRIM(int_to_char(igeom))
        fltherm=TRIM(fltherm_thermo)//'.g'//TRIM(int_to_char(igeom))
        flpstherm=TRIM(flpstherm_thermo)//'.g'//TRIM(int_to_char(igeom))
        flpband=TRIM(flpband_thermo)//'.g'//TRIM(int_to_char(igeom))
        flgnuplot=TRIM(flgnuplot_thermo)//'.g'//TRIM(int_to_char(igeom))
        flpsdisp=TRIM(flpsdisp_thermo)//'.g'//TRIM(int_to_char(igeom))

        IF (nqaux > 0) CALL set_paths_disp()
        !
        ! ... Checking the status of the calculation and if necessary initialize
        ! ... the q mesh and all the representations
        !
        !
        CALL check_initial_status(auxdyn)
        !
        IF ( .NOT. check_dyn_file_exists(auxdyn)) THEN
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
           CALL q2r_sub(auxdyn) 
!
!    compute interpolated dispersions
!
           IF (lmatdyn) THEN
              CALL matdyn_sub(0, igeom)
              CALL plotband_sub(2,igeom,' ')
           ENDIF
!
!    computes phonon dos
!
           IF (lmatdyn.AND.ldos) THEN
!
!    the phonon dos is calculated and the frequencies saved on the ph_freq_save 
!    structure
!
              IF (.NOT.ALLOCATED(phdos_save)) ALLOCATE(phdos_save(ngeo))
              IF (.NOT.ALLOCATED(ph_freq_save)) ALLOCATE(ph_freq_save(ngeo))
              CALL matdyn_sub(1,igeom)
              CALL simple_plot('_dos', fldos, flpsdos, 'frequency (cm^{-1})', &
                       'DOS (states / cm^{-1} / cell)', 'red', freqmin, freqmax, &
                            0.0_DP, 0.0_DP)
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
        CALL clean_pw(.TRUE.)
        CALL close_phq(.FALSE.)
        CALL clean_input_variables()
        CALL destroy_status_run()
        CALL deallocate_part()
     ENDDO
     flgnuplot=TRIM(flgnuplot_thermo)
!
!     Here the Helmholtz free energy at each lattice constant is available.
!     We can write on file the free energy as a function of the volume at
!     any temperature. For each temperature we call the murnaghan equation
!     to fit the data. We save the minimum volume, the bulk modulus and its
!     pressure derivative for each temperature.
!
     IF (lev_syn_2) THEN
        diraux='evdir'
        CALL check_tempdir ( diraux, exst, parallelfs )
        flevdat=TRIM(diraux)//'/'//TRIM(flevdat)
        DO itemp = 1, ntemp
           CALL do_ev_t(itemp)
           CALL do_ev_t_ph(itemp)
        ENDDO
!
!    here we calculate and plot the gruneisen parameters along the given path.
!
        CALL write_gruneisen_band(flfrq_thermo)
        CALL plotband_sub(3,1,flfrq_thermo)
!
!    here we compute the gruneisen parameters on the uniform mesh
!
        CALL compute_gruneisen()
!
!    here we calculate several anharmonic quantities and plot them.
!
        CALL write_anharmonic()
        CALL write_ph_freq_anharmonic()
        CALL write_grun_anharmonic()
        CALL plot_anhar() 
     ENDIF

     IF (lmatdyn.AND.ldos) THEN
        DO igeom=1,ngeo
           CALL destroy_phdos(phdos_save(igeom))
        ENDDO
        DEALLOCATE(phdos_save)
     ENDIF
  ELSE
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

     IF (lelastic_const) THEN
!
!   save the stress tensors calculated by all images
!
        CALL mp_sum(sigma_geo, world_comm)
        sigma_geo=sigma_geo / nproc_image
!        CALL write_stress(nwork, file_dat)
!
!  the elastic constants are calculated here
!
        CALL compute_elastic_constants(sigma_geo, epsilon_geo, nwork, &
                               ngeo_strain, ibrav_save, laue)
        CALL print_elastic_constants(frozen_ions)
!
!  now compute the elastic compliances and prints them
!
        CALL compute_elastic_compliances(el_con,el_compliances)
        CALL print_elastic_compliances(frozen_ions)
!
!  save elastic constants and compliances on file
!
        IF (my_image_id==root_image) CALL write_elastic(fl_el_cons)
     ENDIF
     !
     CALL deallocate_asyn()
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
