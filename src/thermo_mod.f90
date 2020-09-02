!
! Copyright (C) 2013 - 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for the thermo_pw program
!
MODULE thermo_mod
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the type of calculations to do,
  ! ... the number of geometries and the parameters of each geometry
  ! ... Some space to save the results of pw.x at each geometry
  !
  SAVE

  CHARACTER(LEN=30) :: what             ! the type of calculation. See README.

  INTEGER  :: ngeo(6)                   ! number of different geometries 
                                        ! per celldm parameter

  INTEGER :: fact_ngeo(6)               ! factors used to reduce the
                                        ! size of the grid where phonons
                                        ! are interpolated

  INTEGER :: ngeo_ph(6)                 ! mesh where the phonon are calculated
                                        ! can be smaller than ngeo
                                        !
  LOGICAL, ALLOCATABLE :: no_ph(:)      ! decide in which geometries
                                        ! the phonons are calculated
  LOGICAL, ALLOCATABLE :: dynmat_on_file(:) ! set to .TRUE. if the dynmat
                                        ! of this geometry are on file

  INTEGER :: tot_ngeo                   ! total number of geometries for
                                        ! which we calculate the phonon
                                        ! dispersions

  REAL(DP), ALLOCATABLE :: omega_geo(:),    &   ! volume (geometry)
                           celldm_geo(:,:), &   ! The celldm for each geometry
                           energy_geo(:),   &   ! The total energy at each
                                                ! geometry
                           forces_geo(:,:,:), & ! the forces on atoms at each
                                                ! geometry
                           stress_geo(:,:,:)    ! the stress at each 
                                                ! geometry
  REAL(DP) ::              step_ngeo(6)         ! the difference of 
                                                ! parameters among
                                                ! different geometries.
  INTEGER, ALLOCATABLE :: ibrav_geo(:)          ! the Bravais lattice at
                                                ! each geometry
  INTEGER :: central_geo                        ! a reference geometry
  LOGICAL :: reduced_grid                       ! if .TRUE. use a reduced
                                                ! grid to interpolate the
                                                ! geometry
  INTEGER, ALLOCATABLE :: in_degree(:)          ! for each geometry says to
                                                ! which degree it belongs in 
                                                ! the reduced grid
  INTEGER :: red_central_geo                    ! With the reduced_grid option
                                                ! this is the central geometry
                                                ! used by all degrees
  REAL(DP) :: density                           ! the density of the solid

  INTEGER :: max_geometries                     ! This value controls the
                                                ! number of geometries for
                                                ! which phonons are calculated
                                                ! in each run
  INTEGER :: start_geometry, last_geometry      ! These variables control
                                                ! which geometries are done
                                                ! in each run

END MODULE thermo_mod

MODULE temperature
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the temperature
  !
  SAVE

  REAL(DP), ALLOCATABLE :: temp(:)      ! temperature (K), n

  REAL(DP) :: tmin, tmax                ! maximum and minimum temperature (K)
  REAL(DP) :: deltat                    ! delta T
  INTEGER  :: ntemp                     ! number of temperatures

END MODULE temperature

MODULE control_mur
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the calculation of the murnaghan
  !     equation
  !
  SAVE
  REAL(DP) :: vmin, b0, b01, emin   ! the minimum of the murnaghan at T=0
  REAL(DP) :: vmin_input, vmax_input, deltav   ! plot the fitted total energy
                                               ! and pressure from vmin_input
                                               ! to vnax_input in steps of i
                                               ! deltav
  INTEGER :: nvol                              ! the number of volumes for
                                               ! the plot
  LOGICAL :: lmurn             ! if .TRUE. makes a Murnaghan
                               ! fit of the energy as a function of the volume
                               ! otherwise makes a fit of the energy as a 
                               ! function of the celldm parameters, with a
                               ! quadratic or quartic polynomial of dimension 
                               ! up to 6

END MODULE control_mur

MODULE thermodynamics
  !
  USE kinds, ONLY: DP
  USE phdos_module, ONLY : phdos_type, gen_phdos_type
  !
  ! ... The variables needed to save the thermodynmac quantities 
  !     calculated from the phonon dos
  !
  SAVE

  TYPE(phdos_type), ALLOCATABLE :: phdos_save(:) ! phdos for each geometry
                                             ! geometry on a uniform mesh
  TYPE(gen_phdos_type) :: gen_phdos_save     ! generalized phdos
                                             ! geometry on a uniform mesh

  REAL(DP), ALLOCATABLE :: ph_ener(:,:)      ! phonon total energy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_free_ener(:,:) ! phonon free_energy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_entropy(:,:)   ! phonon entropy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_e0(:)          ! zero point energy, geometry
  REAL(DP), ALLOCATABLE :: ph_ce(:,:)        ! phonon heat capacity, T, geometry
  REAL(DP), ALLOCATABLE :: ph_b_fact(:,:,:,:,:)! atomic B factor

END MODULE thermodynamics

MODULE ph_freq_thermodynamics
  !
  USE kinds, ONLY: DP
  USE ph_freq_module, ONLY : ph_freq_type
  !
  ! ... The variables needed to save the thermodynmac quantities 
  !     calculated directly from the phonon frequencies and eigenvectors
  !     the f indicates that the quantities has been calculated 
  !     directly from phonon frequencies
  !
  SAVE

  TYPE(ph_freq_type), ALLOCATABLE :: ph_freq_save(:) ! frequencies for each 

  REAL(DP), ALLOCATABLE :: phf_ener(:,:)      ! phonon total energy, T, geometry
  REAL(DP), ALLOCATABLE :: phf_free_ener(:,:) ! phonon free_energy, T, geometry
  REAL(DP), ALLOCATABLE :: phf_entropy(:,:)   ! phonon entropy, T, geometry
  REAL(DP), ALLOCATABLE :: phf_e0(:)          ! zero point energy, geometry
  REAL(DP), ALLOCATABLE :: phf_ce(:,:)        ! phonon specific heat, T, geometry
  REAL(DP), ALLOCATABLE :: phf_b_fact(:,:,:,:,:) ! atomic B factor

END MODULE ph_freq_thermodynamics

MODULE anharmonic
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the free energy fitted by the Murnaghan
!
  USE kinds, ONLY: DP
  SAVE
!
!   here all the quantities obtained from the free energy calculated from
!   the phonon dos
!
  REAL(DP), ALLOCATABLE :: cp_t(:)     ! isobaric heat capacity (T)
  REAL(DP), ALLOCATABLE :: cv_t(:)     ! isocoric heat capacity (T)
  REAL(DP), ALLOCATABLE :: ener_t(:)   ! vibrational energy (T)
  REAL(DP), ALLOCATABLE :: free_ener_t(:)   ! vibrational free energy (T)
  REAL(DP), ALLOCATABLE :: entropy_t(:)! entropy (T)
  REAL(DP), ALLOCATABLE :: ce_t(:)     ! constant strain heat capacity (T)
  REAL(DP), ALLOCATABLE :: b0_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: alpha_t(:)  ! linear thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: beta_t(:)   ! volume thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: gamma_t(:)  ! average gruneisen parameter
  REAL(DP), ALLOCATABLE :: bths_t(:,:,:)   ! thermal stress
  REAL(DP), ALLOCATABLE :: ggamma_t(:,:,:) ! generalized average gruneisen 
                                       ! parameter

  REAL(DP), ALLOCATABLE :: vmin_t(:), b0_t(:), b01_t(:), free_e_min_t(:) 
                           ! the parameters of the minimum of the  
                           ! free energy calculated from dos at each temperature
  REAL(DP), ALLOCATABLE :: celldm_t(:,:) ! the lattice parameters as a 
                           ! function of temperature
  REAL(DP), ALLOCATABLE :: alpha_anis_t(:,:)  ! thermal expansion tensor

  REAL(DP), ALLOCATABLE :: cpmce_anis(:) ! difference cp-ce computed from
                                         ! elastic constants
  REAL(DP), ALLOCATABLE :: bfact_t(:,:,:)! b factor as a function of 
                                         ! temperature

  LOGICAL :: lelastic=.FALSE.            ! elastic constants available in
                                         ! some approximation
  REAL(DP), ALLOCATABLE :: el_cons_t(:,:,:) ! elastic constants as a function
                                         ! of temperature (constant T)
  REAL(DP), ALLOCATABLE :: el_comp_t(:,:,:) ! elastic compliances as a function
                                         ! of temperature (constant T)
  REAL(DP), ALLOCATABLE :: el_cons_s(:,:,:) ! elastic constants as a function
                                         ! of temperature (constant entropy)
  REAL(DP), ALLOCATABLE :: el_comp_s(:,:,:) ! elastic compliances as a function
                                         ! of temperature (constant entropy)
  REAL(DP), ALLOCATABLE :: macro_el_s(:,:) ! macroscopic elasticity as a 
                                           ! function of t (from el_cons_s)
  REAL(DP), ALLOCATABLE :: macro_el_t(:,:) ! macroscopic elasticity as a 
                                           ! function of t (from el_cons_t)
  REAL(DP), ALLOCATABLE :: v_s(:,:)        ! the sound velocities from
                                           ! macro_el_s: vp, vb, vg
  REAL(DP), ALLOCATABLE :: v_t(:,:)        ! the sound velocities from
                                           ! macro_el_t: vp, vb, vg
  REAL(DP), ALLOCATABLE :: el_con_geo_t(:,:,:,:) ! the temperature dependent
                                           ! elastic constants at all 
                                           ! geometries

END MODULE anharmonic

MODULE ph_freq_anharmonic
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the free energy fitted by the Murnaghan
!
  USE kinds, ONLY: DP
  SAVE
!
!   here all the quantities obtained from the free energy calculated from
!   the frequencies 
!
  REAL(DP), ALLOCATABLE :: cpf_t(:)     ! isobaric heat capacity (T)
  REAL(DP), ALLOCATABLE :: cvf_t(:)     ! isocoric heat capacity (T)
  REAL(DP), ALLOCATABLE :: enerf_t(:)   ! vibrational energy (T)
  REAL(DP), ALLOCATABLE :: free_enerf_t(:) ! vibrational free energy (T)
  REAL(DP), ALLOCATABLE :: entropyf_t(:)! entropy (T)
  REAL(DP), ALLOCATABLE :: cef_t(:)     ! constant strain heat capacity (T)
  REAL(DP), ALLOCATABLE :: b0f_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: alphaf_t(:)  ! linear thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: betaf_t(:)   ! volume thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: gammaf_t(:)  ! average gruneisen parameter
  REAL(DP), ALLOCATABLE :: bthsf_t(:,:,:)   ! thermal stress
  REAL(DP), ALLOCATABLE :: ggammaf_t(:,:,:) ! generalized average gruneisen 
                                       ! parameter

  REAL(DP), ALLOCATABLE :: vminf_t(:), b0f_t(:), b01f_t(:), free_e_minf_t(:) 
                           ! the parameters of the minimum of the  
                           ! free energy  
  REAL(DP), ALLOCATABLE :: celldmf_t(:,:) ! the lattice parameters as a 
                           ! function of temperature
  REAL(DP), ALLOCATABLE :: alphaf_anis_t(:,:)  ! thermal expansion tensor

  REAL(DP), ALLOCATABLE :: cpmcef_anis(:) ! difference cp-ce computed from
                                          ! elastic constants
  REAL(DP), ALLOCATABLE :: bfactf_t(:,:,:)! b factor as a function of 
                                          ! temperature

  LOGICAL :: lelasticf=.FALSE.           ! elastic constants available in
                                         ! some approximation
  REAL(DP), ALLOCATABLE :: el_consf_t(:,:,:) ! elastic constants as a function
                                          ! of temperature
  REAL(DP), ALLOCATABLE :: el_compf_t(:,:,:) ! elastic compliances as a function
                                         ! of temperature
  REAL(DP), ALLOCATABLE :: el_consf_s(:,:,:) ! elastic constants as a function
                                         ! of temperature (constant entropy)
  REAL(DP), ALLOCATABLE :: el_compf_s(:,:,:) ! elastic compliances as a function
                                         ! of temperature (constant entropy)
  REAL(DP), ALLOCATABLE :: macro_elf_s(:,:) ! macroscopic elasticity as a 
                                            ! function of t (from el_consf_s)
  REAL(DP), ALLOCATABLE :: macro_elf_t(:,:) ! macroscopic elasticity as a 
                                            ! function of t (from el_consf_t)
  REAL(DP), ALLOCATABLE :: vf_s(:,:)        ! the sound velocities from
                                            ! macro_elf_s: vp, vb, vg
  REAL(DP), ALLOCATABLE :: vf_t(:,:)        ! the sound velocities from
                                            ! macro_elf_t: vp, vb, vg
  REAL(DP), ALLOCATABLE :: el_conf_geo_t(:,:,:,:) ! the temperature dependent
                                            ! elastic constants at all 
                                            ! geometries
END MODULE ph_freq_anharmonic

MODULE grun_anharmonic
  USE kinds, ONLY: DP
  USE polynomial, ONLY : poly2
  USE ph_freq_module, ONLY : ph_freq_type
  !
  !  This module contains all the quantities calculated from the gruneisen
  !  parameters
  !
  SAVE

  REAL(DP), ALLOCATABLE :: betab(:)    ! volume thermal expansion multiplied
                                       ! by the bulk modulus
  REAL(DP), ALLOCATABLE :: alpha_an_g(:,:) ! thermal expansion tensor at each
                                        !    temperature
  REAL(DP), ALLOCATABLE :: cp_grun_t(:) ! isobaric heat capacity as a 
                                        ! function of temperature
  REAL(DP), ALLOCATABLE :: ce_grun_t(:) ! constant strain heat capacity 
                                        ! used for thermodynamic quantities 
                                        ! with gruneisen parameters
  REAL(DP), ALLOCATABLE :: cv_grun_t(:) ! isochoric heat capacity as a
                                        ! function of temperature
  REAL(DP), ALLOCATABLE :: b0_grun_s(:) ! isoentropic bulk modulus as a 
                                        ! function of temperature
  REAL(DP), ALLOCATABLE :: grun_gamma_t(:)  ! average gruneisen parameter
  REAL(DP), ALLOCATABLE :: grun_cpmce_anis(:) ! difference cp-c_epsilon 
                                          ! computed from elastic constants

  TYPE(poly2), ALLOCATABLE :: p_grun_p2(:,:) ! For each band and each q point
                                             ! these are the coefficients of 
                                             ! polynomial which fit the 
                                             ! frequencies as a function of
                                             ! the crystal parameters
  REAL(DP), ALLOCATABLE :: poly_grun(:,:,:)
  REAL(DP), ALLOCATABLE :: poly_grun_red(:,:,:,:) ! For each band, each q point
                                        ! and each crystal parameter
                                        ! these are the coefficients of 
                                        ! polynomial which fit the frequency
                                        ! as a function of crystal parameter
                                        ! this is for the reduced_grid case
  INTEGER :: poly_degree_grun           ! degree of the polynomial used to
                                        ! intepolate the frequencies
  LOGICAL :: done_grun=.FALSE.          ! the anharmonic quantities with
                                        ! Gruneisen parameters have been 
                                        ! calculated
  LOGICAL :: lelastic_grun=.FALSE.      ! elastic constants available in
                                        ! some approximation
  REAL(DP), ALLOCATABLE :: el_cons_grun_t(:,:,:) ! temperature dependent
                                ! elastic constant for thermodynamic quantities
                                ! with gruneisen parameters
  REAL(DP), ALLOCATABLE :: el_comp_grun_t(:,:,:) ! temperature dependent
                                ! elastic compliances for thermodynamic 
                                ! quantities with gruneisen parameters

END MODULE grun_anharmonic

MODULE ifc
  !
  USE kinds, ONLY: DP
  SAVE
  REAL(DP), ALLOCATABLE :: frc(:,:,:,:,:,:,:), zeu(:,:,:), &
               m_loc(:,:), wscache(:,:,:,:,:)
  REAL(DP) :: epsil_ifc(3,3)
  LOGICAL :: has_zstar
  ! frc : interatomic force constants in real space
  ! zeu : effective charges for the original cell
  ! m_loc: the magnetic moments of each atom
  ! wscache: the weight of each q point 
  !
  CHARACTER(LEN=3), ALLOCATABLE :: atm(:)

  CHARACTER(LEN=10) :: zasr      ! the type of asr
  ! 
END MODULE ifc

MODULE control_dosq
  !
  !   This module contains the quantities to make an integration
  !   over the Brilluin zone. The coordinates of the q points, their weights,
  !   if tetraedra are used also ntetra and tetra.
  !
  USE kinds, ONLY: DP
  SAVE

  INTEGER :: nq1_d, nq2_d, nq3_d    ! grid for phonon dos

  REAL(DP), ALLOCATABLE :: dos_q(:,:)  ! the cartesian coordinate of the q
                                       ! points for dos
  INTEGER :: ntetra
  REAL(DP), ALLOCATABLE ::  tetra(:,:,:)
  !
  INTEGER :: ndos_input          ! number of points in the dos plot
  REAL(DP) :: freqmin, freqmax   ! dos minimum and maximun frequency 
                                 ! at this geometry
  REAL(DP) :: freqmin_input, freqmax_input   ! dos minimum and maximun frequency 
                                 ! be given as input (in cm^{-1})
  REAL(DP) :: deltafreq          ! delta frequency in the dos plot. (in cm^{-1})
  REAL(DP) :: phdos_sigma        ! smearing for phdos calculation.

END MODULE control_dosq

MODULE control_thermo
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the type of calculation to do
  !
  SAVE
  !
  LOGICAL, ALLOCATABLE :: lpwscf(:),  & ! if .true. this work requires a scf calc.
                          lberry(:),  & ! if .true. this work requires 
                                        ! a berry_phase calculation
                          lstress(:), & ! if .true. this work computes stress
                          lphonon(:)    ! if .true. this work requires a phonon

  INTEGER, ALLOCATABLE :: geometry(:), & ! the geometry of a given phonon
                                         ! calculation for each work
                          iqw(:),      & ! the q vector of a given phonon
                                         ! calculation for each work
                          irrw(:)        ! the irrep of a given phonon 
                                         ! calculation for each work
  LOGICAL, ALLOCATABLE :: comp_f_iw(:,:), & ! for each work the list of
                                         ! frequencies to compute
                          comp_irr_iq_iw(:,:,:), & ! for each work the 
                                         ! comp_irr_iq of that work
                          comp_iq_iw(:,:) ! for each work the list
                                         ! of comp_iq of that work.
  LOGICAL :: read_paths      ! the paths for dispersion are read from input
  LOGICAL :: lev_syn_1=.FALSE. ! if .true. must calculate the murnaghan
                               ! at T=0
  LOGICAL :: lev_syn_2=.FALSE. ! if .true. must calculate the murnaghan
                               ! at all T
  LOGICAL :: lpwscf_syn_1=.FALSE. ! if .true. must calculate the synchronous pw
                                  ! for scf calculation
  LOGICAL :: lbands_syn_1=.FALSE. ! if .true. must calculate the synchronous pw
                                  ! for nscf calculation
  LOGICAL :: ldos_syn_1=.FALSE.   ! if .true. must calculate the synchronous pw
                                  ! for nscf calculation for dos
  LOGICAL :: lpart2_pw=.FALSE.    ! if .true. in the second part makes
                                  ! also pw calculations           
  LOGICAL :: with_eigen=.FALSE.   ! save phonon frequencies and eigenvectors

  LOGICAL :: ltherm_dos=.TRUE.    ! computes the phonon dos and the thermal
                                  ! properties from phonon dos
  LOGICAL :: ltherm_freq=.TRUE.   ! computes the phonon frequencies and 
                                  ! the thermal properties from phonon dos
                                  !
  LOGICAL :: lconv_ke_test=.FALSE.! if .true. this writes the ke test on file
  LOGICAL :: lconv_nk_test=.FALSE.! if .true. this writes the k-point on file
  LOGICAL :: lelastic_const=.FALSE. ! if .true. compute elastic constants
  LOGICAL :: lpiezoelectric_tensor=.FALSE. ! if .true. compute the piezoelectric 
                                  ! tensor
  LOGICAL :: lpolarization=.FALSE. ! if .true. compute the piezoelectric 
  LOGICAL :: lph=.FALSE.    ! if .true. must calculate phonon
  LOGICAL :: ltherm=.FALSE. ! if .true. the thermodynamic properties are
                            ! calculated
  LOGICAL :: after_disp=.FALSE. ! if .true. dynamical matrix files are supposed
                            ! to be already on disk. fildyn must be read
                            ! from thermo_control.
  LOGICAL :: lq2r=.FALSE.   ! if .true. the interatomic force constants 
                            ! are calculated and the phonon are interpolated

  LOGICAL :: do_scf_relax   ! to be used with cells with internal 
                            ! degrees of freedom. If .true. reoptimize
                            ! the atomic coordinates after finding 
                            ! the equilibrium crystal parameters. 
                            ! Not necessary for phonons, but used with
                            ! elastic constants and piezoelectric tensor
  LOGICAL :: continue_zero_ibrav ! if .TRUE. run with ibrav=0

  LOGICAL :: find_ibrav     ! if .TRUE. continue the calculation with the
                            ! converted input
  LOGICAL :: lectqha=.FALSE. ! if .true. compute the elastic constants 
                            ! within the qha at many geometries
  !
  CHARACTER(LEN=256) :: outdir_thermo ! the outdir read from the input
  !
  LOGICAL :: set_internal_path ! the path provided by thermo_pw is used
  !
  LOGICAL :: set_2d_path    ! the path provided by thermo_pw is used
  !
  LOGICAL :: all_geometries_together ! if true phonon dispersions for all
                            !  geometries are calculated together
  INTEGER :: spin_component ! the spin

  REAL(DP) :: max_seconds_tpw ! used to stop the run after max_seconds with
                              ! images
  !
END MODULE control_thermo

MODULE distribute_collection

SAVE

LOGICAL, ALLOCATABLE :: me_igeom(:) ! if .TRUE. this image must collect
                                    ! some work in this geometry

LOGICAL, ALLOCATABLE :: me_igeom_iq(:,:) ! if .TRUE. this image must collect
                                    ! some q point in this geometry

END MODULE distribute_collection

MODULE control_elastic_constants
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the calculation of the elastic 
  !     constants
  !
  SAVE
  REAL(DP) :: delta_epsilon               ! the distance between two strains (D)
  REAL(DP) :: epsilon_0                   ! a minimum strain (e_0). 
                                          ! For ngeo_strain even, the
                                          ! strains will be:
                                          ! -e_0-D/2, -e0-3D/2, ...
                                          ! e_0+D/2, e0+3D/2, ...
                                          ! For ngeo_strain odd:
                                          !  0, -e0-D, -e0-2D, ... 
                                          !      e0+D,  e0+2D
  REAL(DP), ALLOCATABLE :: rot_mat(:,:,:) ! rotation matrix between the
                                          ! cartesian coordinates of the
                                          ! strained and unstrained cell  
  INTEGER :: ngeo_strain        ! number of strain configurations

  LOGICAL :: frozen_ions        ! if .true. compute the elastic constant 
                                ! keeping the ions frozen at the strained
                                ! positions

  CHARACTER(LEN=12) :: elastic_algorithm ! can be standard, advanced, 
                                ! energy_std or energy
                                ! it chooses the routines to use to 
                                ! calculate elastic constants.
  LOGICAL :: elalgen=.FALSE.    ! elastic_algorithm requires energy
  INTEGER :: elcpvar            ! number of variables of the polynomial 
                                ! function that interpolates stress or energy
  INTEGER :: poly_degree        ! the order of the polynomial to interpolate
                                ! the stress-strain or the energy-strain
                                ! functions. Could differ from poly_degree_elc
                                ! which is the order of the polynomial
                                ! that interpolates the elastic constants
                                ! among different geometries

  LOGICAL :: el_cons_available=.FALSE.  ! when this flag becomes true it
                                ! means that the elastic constant have been
                                ! read from file and are available

  LOGICAL :: el_cons_t_available=.FALSE.  ! when this flag becomes true it
                                ! means that the elastic constants for each 
                                ! geometry are available and can be 
                                ! interpolated at each temperature
  LOGICAL :: el_cons_qha_available=.FALSE. ! when this flag becomes true 
                                ! the elastic constants computed 
                                ! within quasi harmonic approximation 
                                ! and BZ integrations with phdos are available 
                                ! at least for one geometry
  LOGICAL :: el_consf_qha_available=.FALSE. ! when this flag becomes true
                                ! the elastic constants computed 
                                ! within quasi harmonic approximation 
                                ! and BZ integrations are available at least
                                ! for one geometry
  LOGICAL :: el_cons_qha_geo_available=.FALSE. ! when this flag becomes true
                                ! the QHA elastic constants computed with
                                ! phdos are available for all geometries
  LOGICAL :: el_consf_qha_geo_available=.FALSE. ! when this flag becomes
                                ! true the QHA with BZ integration 
                                ! elastic constants are available for 
                                ! all geometries

  REAL(DP), ALLOCATABLE :: el_con_geo(:,:,:)  ! the elastic constants at
                                ! each geometry
  INTEGER, ALLOCATABLE :: el_con_ibrav_geo(:) ! the ibrav of each unperturbed
                                ! lattice.
  REAL(DP), ALLOCATABLE :: el_con_celldm_geo(:,:) ! the celldm of each 
                                ! unperturbed lattice.
  REAL(DP), ALLOCATABLE :: el_con_tau_crys_geo(:,:,:) ! the atomic positions of each
  REAL(DP), ALLOCATABLE :: el_con_omega_geo(:) ! the volume of each 
                                ! unperturbed cell.
  REAL(DP), ALLOCATABLE :: epsil_geo(:) ! strain amplitude for each geometry
                                !
  INTEGER :: ngeom=1            ! the number of geometries

  INTEGER :: work_base          ! number of works for one set of
                                ! elastic constants

  LOGICAL :: use_free_energy    ! if true makes the qha approximation 
                                ! otherwise the quasi-static one.
  INTEGER :: start_geometry_qha ! initial geometry calculated in this run

  INTEGER :: last_geometry_qha  ! last_geometry calculated in this run 

END MODULE control_elastic_constants
  !
MODULE control_conv
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the control of convergence
  !
  SAVE

  REAL(DP), ALLOCATABLE :: ke(:)    ! the kinetic energy that are tested
  REAL(DP), ALLOCATABLE :: keden(:) ! the kinetic energy of the charge
  REAL(DP), ALLOCATABLE :: sigma_test(:) ! the smearing value
  INTEGER,  ALLOCATABLE :: nk_test(:,:) ! the nk to calculate
  REAL(DP) :: deltake               ! the interval between kinetic energies
  INTEGER  :: nke                   ! the number of kinetic energies
  REAL(DP) :: deltakeden            ! the interval for kinetic energy for density
  INTEGER  :: nkeden                ! the number of kinetic energies for densiy
  INTEGER  :: ncutoffene            ! combined number of cut-off energies
  INTEGER  :: nnk                   ! the number of nk values to test k points
  INTEGER  :: deltank(3)            ! the step between nk values to test k points
  INTEGER  :: nsigma                ! the number of smearing values to test
  REAL(DP) :: deltasigma            ! the step between smaring values

END MODULE control_conv

MODULE control_paths
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the paths
  !
  SAVE
  !
  REAL(DP), ALLOCATABLE :: xqaux(:,:)     ! the initial and final points 
  REAL(DP), ALLOCATABLE :: wqauxr(:)      ! the weight for normal paths
  INTEGER, ALLOCATABLE  :: wqaux(:)       ! the number of points per line
  INTEGER, ALLOCATABLE  :: nrap_plot_in(:)! the number of representations 
                                          ! per line 
  INTEGER, ALLOCATABLE :: rap_plot_in(:,:)! which representations per line
  INTEGER, ALLOCATABLE :: nrap_plot(:)    ! number of representations per k
                                          ! point
  INTEGER, ALLOCATABLE :: rap_plot(:,:)   ! which representation per k point
  INTEGER               :: disp_nqs       ! total number of points to compute
  REAL(DP), ALLOCATABLE :: disp_q(:,:), disp_wq(:)  ! q path for interpolated
  INTEGER            :: npk_label         ! number of label points
  INTEGER            :: nqaux             ! number of auxiliary points
  INTEGER            :: npx               ! supercell size for BZ search
  CHARACTER(LEN=3), ALLOCATABLE :: letter(:) ! the labels
  CHARACTER(LEN=3), ALLOCATABLE :: letter_path(:) ! the optional labels
  INTEGER, ALLOCATABLE :: label_list(:)   ! correspondence label xqaux list
  INTEGER, ALLOCATABLE :: label_disp_q(:)   ! correspondence label disp_q list
  REAL(DP)          :: path_fact          ! a factor used to multiply the
                                          ! number of points of the default path
  REAL(DP)          :: dkmod_save         ! the typical distance between two
                                          ! points.
  CHARACTER(LEN=10) :: point_label_type   ! type of point labels
  LOGICAL, ALLOCATABLE :: high_sym_path(:)! high_symmetry points along the path
  LOGICAL  :: long_path      ! if .TRUE. use the complete path in the BZ
  LOGICAL  :: old_path      ! if .TRUE. use the alternative path in the BZ
  LOGICAL :: q_in_band_form
  LOGICAL :: q_in_cryst_coord
  LOGICAL :: q2d                 ! if .true. the k are in a 2d grid
  LOGICAL :: is_a_path=.TRUE.    ! if .true. the k point are in a path
  LOGICAL :: lbar_label=.FALSE.  ! if .TRUE. put a bar on each label 
                                 ! (needed for surface band plot)
END MODULE control_paths

MODULE control_bands
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the band plot
  !
  SAVE
  !
  INTEGER :: nbnd_bands                   ! number of bands for band calculation
  !
  REAL(DP) :: emin_input     ! minimum energy of the plot (eV)
  REAL(DP) :: emax_input     ! maximum energy of the plot (eV)
  LOGICAL  :: lsym           ! if .TRUE. does the symmetry analysis of the bands
  LOGICAL  :: enhance_plot   ! if .TRUE. an enhaced band plot is done
  LOGICAL  :: only_bands     ! if .TRUE. skip the scf calculation

END MODULE control_bands

MODULE control_2d_bands
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the two dimensional band structure
  !
  SAVE

  LOGICAL :: lprojpbs, &  ! if .TRUE. plot the projected band structure (PBS)
             sym_divide, & ! if .TRUE. the bands belonging to different 
                         ! irreps are ploted in different panels
             identify_sur, & ! if .TRUE. identify the surface states
             only_bands_plot ! if .TRUE. does not recalculate the bands

  INTEGER :: nkz         ! the number of values of k_z used for the PBS
  INTEGER :: sur_layers  ! number of surface layers
  INTEGER :: nlayers     ! number of layers identified by the code
  INTEGER :: surface1, surface2 ! tha surface layers
  REAL(DP) :: gap_thr, & ! energy gap in the PBS
              sur_thr, & ! minimum percentage of charge to be a surface state
              sp_min     ! minimum layer spacing. Given in input. 

  INTEGER, ALLOCATABLE  :: aux_ind_sur(:,:)
  REAL(DP), ALLOCATABLE :: averag(:,:,:,:)      ! charge density on each layer 
  REAL(DP), ALLOCATABLE :: vacuum(:,:,:)        ! charge on vacuum
  LOGICAL, ALLOCATABLE  :: lsurface_state(:,:)  ! if .TRUE. a given state is 
                                                ! a surface state
  LOGICAL, ALLOCATABLE  :: lsurface_state_rap(:,:)  ! the same info but
                                                ! divided for the different
                                                ! representations 
  LOGICAL :: force_bands ! if .TRUE. the bands are plotted in all cases
  LOGICAL :: dump_states ! if .TRUE. dumps the planar average on each state
                         !  on file
  LOGICAL :: subtract_vacuum ! if .TRUE. the charge density of each state
                         ! on vacuum is subtracted
END MODULE control_2d_bands

MODULE proj_rap_point_group

  USE kinds, ONLY: DP
  SAVE

  INTEGER ::  which_elem(48)   ! for each symmetry element says which 
                               ! element it is in the list of symmetry
                               ! operations of that group.
  INTEGER ::  group_desc(48)   ! for each symmetry of the list of
                               ! operation of a given group says which
                               ! element it is in the global list of
                               ! symmetry operations
  INTEGER :: nrap_proj         ! number of projective representations
  INTEGER :: nsym_proj         ! number of symmetry elements of the group
  INTEGER :: code_groupq_ext   ! the extended code of the point group
  INTEGER :: qptype(3)         ! the ptype of the present q
  REAL(DP) :: qgauge(48)        ! the gauge of the current q
  INTEGER :: lqproj            ! if .true. the present q is projective
  COMPLEX(DP) :: char_mat_proj(48,48)
  CHARACTER(LEN=45)  :: name_rap_proj(48)  ! the name of each projective
                                 ! irreducible projective representation
END MODULE proj_rap_point_group

MODULE control_eldos

  USE kinds, ONLY: DP
  SAVE

  REAL(DP) :: deltae        ! the distance between energy points in the plot
  INTEGER :: ndose               ! number of points in the dos plot
  INTEGER :: nk1_d, nk2_d, nk3_d ! the k-point mesh for dos
  INTEGER :: k1_d, k2_d, k3_d    ! the k-point mesh shift
  REAL(DP) :: sigmae        ! the optional smearing different from the one
                            ! used in the pw.x calculation
  LOGICAL :: legauss        ! if .true. use gaussian smearing for dos plot
  INTEGER :: save_ndos      ! number of points in the file 
  REAL(DP), ALLOCATABLE :: dos_k(:,:)  
  REAL(DP), ALLOCATABLE :: dos_wk(:)
  
END MODULE control_eldos

MODULE control_grun

  USE kinds, ONLY: DP
  SAVE
                                 
  REAL(DP) :: grunmin_input, &  ! minimum value for gruneisen plot
              grunmax_input     ! maximum value for gruneisen plot
                                !
  REAL(DP) :: volume_ph         ! volume at which the phonon frequencies
                                ! and the Gruneisen parameters dispersions
                                ! are interpolated.
                                ! 
  REAL(DP) :: celldm_ph(6)      ! cell parameters at which the phonon
                                ! frequencies and the Gruneisen parameters
                                ! dispersions are interpolated.
                                !
  REAL(DP) :: temp_ph           ! when volume_ph=0.0_DP or celldm_ph=0.0_DP
                                ! we use the volume or celldm at this 
                                ! temperature
                                !
  LOGICAL :: lv0_t,    &        ! if .TRUE. use the temperature dependent
                                ! volume for the calculation of anharmonic
                                ! quantities from gruneisen parameters
             lb0_t              ! if .TRUE. use the temperature dependent
                                ! bulk modulus for the calculation of 
                                ! anharmonic quantities from gruneisen 
                                ! parameters
                                !
  REAL(DP), ALLOCATABLE :: vgrun_t(:), &   ! temperature dependent volume for
                                ! thermodynamic quantities with gruneisen 
                                ! parameters 
                           b0_grun_t(:), & ! temperature dependent bulk modulus
                                ! for thermodynamic quantities with gruneisen
                                ! parameters
                           celldm_grun_t(:,:)  ! temperature_dependent celldm
                                ! for thermodynamic quantities with gruneisen
                                ! parameters
END MODULE control_grun

MODULE initial_conf
  USE kinds, ONLY: DP
  USE collect_info, ONLY : collect_info_type
  SAVE

  INTEGER  :: ibrav_save       ! save the Bravais lattice
  REAL(DP) :: celldm_save(6)   ! save the crystal parameters
  REAL(DP) :: omega_save       ! the volume of the initial configuration
  INTEGER, ALLOCATABLE :: ityp_save(:)  ! save the type of atoms. To be
                                        ! used after completely cleaning pw
  REAL(DP) :: at_save(3,3)     ! save the at of the configuration read by pw.x 
  REAL(DP) :: bg_save(3,3)     ! save the bg of the configuration read by pw.x 
  REAL(DP), ALLOCATABLE :: tau_save(:,:) ! save the atomic coordinates read
                               ! from pw.x input
  REAL(DP), ALLOCATABLE :: tau_save_crys(:,:) ! save the atomic coordinates read
                               ! from pw.x input in crystal coordinates
  INTEGER  :: nr1_save, nr2_save, nr3_save  ! save the fft dimensions
  LOGICAL :: nosym_save        ! save the input nosym

  TYPE(collect_info_type), ALLOCATABLE :: collect_info_save(:)

  LOGICAL :: epsil_save   ! save input epsil, changed in a dispersion run
  LOGICAL :: zeu_save     ! save input_zeu, changed in dispersion run
  LOGICAL :: zue_save     ! save input_zue, changed in dispersion run
  INTEGER :: start_q_save, last_q_save ! save start_q and last_q

END MODULE initial_conf

MODULE initial_param
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to save the initial parameters
  !
  SAVE

  REAL(DP) :: ecutwfc0  ! initial cutoff on wavefunction
  REAL(DP) :: ecutrho0  ! initial cutoff on the charge density

  REAL(DP) :: ethr0     ! the initial accuracy of the eigenvalues
  

END MODULE initial_param

MODULE equilibrium_conf
  USE kinds, ONLY: DP
  SAVE

  REAL(DP) :: celldm0(6)       ! equilibrium crystal parameters
  REAL(DP) :: omega0           ! the volume at the equilibrium configuration
  REAL(DP) :: at0(3,3)         ! save the equilibrium lattice vectors
  REAL(DP) :: bg0(3,3)         ! save the equilibrium reciprocal lattice vectors
  REAL(DP), ALLOCATABLE :: tau0(:,:) ! save the atomic coordinates read
                               ! from pw.x input
  REAL(DP), ALLOCATABLE :: tau0_crys(:,:) ! save the atomic coordinates 
                               ! in crystal coordinates
  INTEGER  :: nr1_0, nr2_0, nr3_0  ! the fft dimensions of the equilibrium
                                   !  configuration
END MODULE equilibrium_conf

MODULE control_pwrun

  USE kinds, ONLY: DP
  SAVE

  LOGICAL  :: do_punch=.TRUE.  ! set this variable to .FALSE. if pw has
                               ! not to save the punch files.
END MODULE control_pwrun

MODULE control_phrun

  SAVE

  CHARACTER(LEN=256)  :: auxdyn=""  ! the name of the dynamical matrix file

END MODULE control_phrun

MODULE control_energy_plot

  USE kinds, ONLY: DP
  SAVE

  INTEGER :: ncontours
  REAL(DP), ALLOCATABLE :: ene_levels(:)
  CHARACTER(LEN=12), ALLOCATABLE :: color_levels(:)

END MODULE control_energy_plot

MODULE control_debye

  USE kinds, ONLY: DP
  SAVE

  REAL(DP) :: debye_t, deb_e0

  REAL(DP), ALLOCATABLE :: deb_energy(:),       &  ! Debye vibrational energy
                           deb_free_energy(:),  &  ! Debye free energy
                           deb_entropy(:),      &  ! Debye entropy
                           deb_cv(:),           &  ! Debye cv
                           deb_b_fact(:,:,:,:), &  ! axiliary B factor
                           deb_bfact(:)            ! Debye atomic B factors

END MODULE control_debye

MODULE control_macro_elasticity

  USE kinds, ONLY: DP
  SAVE

  REAL(DP) :: macro_el(8)     ! the Voigt and Reuss approximations
  REAL(DP) :: vp, vb, vg      ! the sound velocities
  REAL(DP) :: approx_debye_t  ! approximate Debye temperature

END MODULE control_macro_elasticity

MODULE control_quadratic_energy

  USE kinds, ONLY: DP
  USE polynomial, ONLY : poly2
  SAVE

  INTEGER :: nvar                ! number of independent variables
  INTEGER :: ncoeff              ! number of coefficients of the polynomial fit
  INTEGER :: ncoeff2             ! number of coefficient of the quadratic
                                 ! part of the polynomial
  REAL(DP), ALLOCATABLE :: hessian_v(:,:), &   ! hessian eigenvectors
                           hessian_e(:),   &   ! hessian eigenvalues
                           x_pos_min(:)        ! coordinates of the minimum
  TYPE(poly2)           :: p2
  LOGICAL :: show_fit            ! if .TRUE. show the countour plots of the
                                 ! fitted polynomial
END MODULE control_quadratic_energy

MODULE control_quartic_energy

  USE kinds, ONLY: DP
  USE polynomial, ONLY : poly4
  SAVE

  LOGICAL :: lquartic                        ! if .TRUE. fit the energy/
                                             ! enthalpy
                                             ! with a quartic polynomial
  INTEGER :: poly_degree_ph                  ! degree of the polynomial that 
                                             ! fits the free energy.
  INTEGER :: poly_degree_thermo              ! degree of the polynomial that 
                                             ! fits the heat capacity, 
                                             ! vibrational energy and entropy
                                             ! and zero point energy
  INTEGER :: poly_degree_bfact               ! degree of the polynomial that 
                                             ! fits the b factor
  INTEGER :: poly_degree_elc                 ! degree of the polynomial that 
                                             ! fits the elastic constants
  INTEGER :: ncoeff4                         ! number of coefficients of 
                                             ! the polynomial fit
  REAL(DP), ALLOCATABLE :: x_min_4(:),   &   ! coordinates of the minimum
                           coeff4(:)         ! coefficients of quartic fit
  INTEGER :: lsolve                          ! 1, 2, 3 controls the method
                                             ! used to find the polynomial
                                             ! coefficients (Default 2)
  TYPE(poly4) :: p4                          ! coefficients of the polynomial

END MODULE control_quartic_energy

MODULE control_pressure
!
!   Some quantities can be calculated at finite pressure. This
!   module has the variables to control the external pressure
!
  USE kinds, ONLY: DP
  SAVE
!
  REAL(DP) :: pressure, & ! Some properties can be calculated
                          ! at fixed pressure. Read in kbar and converted
                          ! to Ry/(a.u.)^2
              pressure_kb ! the pressure in kbar.


END MODULE control_pressure

MODULE control_xrdp

  USE kinds, ONLY: DP
  SAVE

  REAL(DP) :: lambda              ! wavelenght in Angstrom of the X-rays
  CHARACTER(LEN=2) :: lambda_elem ! anode element

  CHARACTER(LEN=256) :: flxrdp    ! file where the data are written 
  CHARACTER(LEN=256) :: flpsxrdp  ! postscript file with the X-ray diffraction
                                  ! spectrum
  LOGICAL :: lformf               ! If true makes plot of the form factor
                                  ! of each atom type
  CHARACTER(LEN=256) :: flformf   ! file where the data are written 
  CHARACTER(LEN=256) :: flpsformf ! postscript file with the X-ray diffraction

  REAL(DP) :: smin, smax          ! the minimum and maximum s for the 
                                  ! atomic form factor plot
  INTEGER :: nspoint
  LOGICAL :: lcm                  ! if .true. use the Cromer-Mann coefficients
                                  ! to calculate the form factors, otherwise
                                  ! use Doyle-Turner or Smith-Burge parameters
  LOGICAL :: lxrdp                ! if .true. compute xrdp after the crystal
                                  ! parameters optimization

END MODULE control_xrdp

MODULE control_gnuplot
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the gnuplot. The files that
  !     contains the gnuplot scripts
  !
  SAVE

  CHARACTER(LEN=256) :: flgnuplot ! the name of file with the gnuplot script
  CHARACTER(LEN=256) :: gnuplot_command ! the gnuplot command
  CHARACTER(LEN=8) :: flext ! presently .ps or .pdf for output postscript of
                            ! pdf files

  LOGICAL :: lgnuplot    ! set to false not to use gnuplot

END MODULE control_gnuplot

MODULE data_files
  USE kinds,  ONLY : DP
  !
  ! ... The names of the files that contains the data. These files are
  !     usually readable and can be plotted directly by gnuplot, or
  !     can be read by the code plotband.f90 of QE and transformed in
  !     a format readable by gnuplot, or contains data saved by thermo_pw
  !     and must be read by this code.
  !
  SAVE

  CHARACTER(LEN=256) :: flfrc   ! file with the force constants
  CHARACTER(LEN=256) :: flfrq   ! file with the frequencies 
  CHARACTER(LEN=256) :: fldos   ! file with the phonon dos
  CHARACTER(LEN=256) :: fldosfrq   ! file with the frequencies for phonon dos
  CHARACTER(LEN=256) :: fleldos ! file with the electron dos
  CHARACTER(LEN=256) :: fltherm ! file with the harmonic thermodynamic 
  CHARACTER(LEN=256) :: fleltherm ! file with the electronic thermodynamic 
  CHARACTER(LEN=256) :: flanhar ! file with the anharmonic quantities
  CHARACTER(LEN=256) :: flevdat ! file with data for ev.x 
  CHARACTER(LEN=256) :: flenergy  ! the name of the file with the energy
                                  ! suited for gnuplot contour plots
  CHARACTER(LEN=256) :: flepsilon ! the name of the file with the dielectric
                                  ! constant
  CHARACTER(LEN=256) :: filband ! file with the bands, readable by plotband 
  CHARACTER(LEN=256) :: flpband ! the name of the file with the bands in 
                                ! a format readable by gnuplot
  CHARACTER(LEN=256) :: flkeconv ! file with the ke cut-off convergence
  CHARACTER(LEN=256) :: flnkconv ! file the k point convergence
  CHARACTER(LEN=256) :: flgrun   ! file with Gruneisen parameters
  CHARACTER(LEN=256) :: flpgrun  ! the name of the output file with the 
                                 ! gruneisen parameters in a format 
                                 ! readable by gnuplot
  CHARACTER(LEN=256) :: flpbs    ! the name of the file with the pbs
  CHARACTER(LEN=256) :: flvec    ! the name of the file with the eigenvectors
  CHARACTER(LEN=256) :: flprojlayer ! the name of the file with the projections
                                  ! of the wavefunctions on each layer
  CHARACTER(LEN=256) :: fl_el_cons ! the file where the elastic constants are
                                   ! written

END MODULE data_files

MODULE postscript_files
  USE kinds,  ONLY : DP
  !
  ! ... The names of the postscript files
  !
  SAVE

  CHARACTER(LEN=256) :: flpsmur   ! the name of the postscript file with E(V)
  CHARACTER(LEN=256) :: flpsband  ! the name of the output postscript file
  CHARACTER(LEN=256) :: flpsdisp  ! the name of the output postscript file
  CHARACTER(LEN=256) :: flpsdos   ! the name of the postscript file with dos
  CHARACTER(LEN=256) :: flpseldos ! the name of the postscript file with el dos
  CHARACTER(LEN=256) :: flpstherm ! the name of the postscript file with 
                                  ! thermodynamic quantities
  CHARACTER(LEN=256) :: flpseltherm ! the name of the postscript file with 
                                  ! electronic thermodynamic quantities
  CHARACTER(LEN=256) :: flpsanhar ! the name of the postscript file with 
                                  ! anharmonic quantities
  CHARACTER(LEN=256) :: flpskeconv! the name of the postscript file with 
                                  ! total energy at different cut-offs
  CHARACTER(LEN=256) :: flpsnkconv! the name of the postscript file with 
                                  ! total energy at different k points
  CHARACTER(LEN=256) :: flpsgrun  ! the name of the postscript file with 
                                  ! the gruneisen parameters
  CHARACTER(LEN=256) :: flpsenergy  ! the name of the postscript file with 
                                  ! the energy contours
  CHARACTER(LEN=256) :: flpsepsilon ! the name of the file with the dielectric
                                  ! constant
END MODULE postscript_files

MODULE internal_files_names
  USE kinds,  ONLY : DP
  !
  ! ... The names of the files before adding the geometry information
  !
  SAVE

  CHARACTER(LEN=256) :: fildyn_thermo    ! dynamical matrix file
  CHARACTER(LEN=256) :: flfrc_thermo     ! interatomic force constants file
  CHARACTER(LEN=256) :: flfrq_thermo     ! frequencies files for dispersion
  CHARACTER(LEN=256) :: fldos_thermo     ! dos file
  CHARACTER(LEN=256) :: fldosfrq_thermo  ! frequencies for dos 
  CHARACTER(LEN=256) :: fltherm_thermo   ! thermodynamic quantities files
  CHARACTER(LEN=256) :: flpband_thermo   ! dispersions in gnuplot format

  CHARACTER(LEN=256) :: flgnuplot_thermo  ! file with gnuplot script

  CHARACTER(LEN=256) :: flpsdisp_thermo   ! postscript file with dispersions
  CHARACTER(LEN=256) :: flpsdos_thermo    ! postscript file with dos
  CHARACTER(LEN=256) :: flpstherm_thermo  ! postscript file with thermodynamic
                                          ! quantities
  CHARACTER(LEN=256) :: flvec_thermo      ! the name of the file with the eigenvectors

END MODULE internal_files_names

MODULE control_asy
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the brillouin zone plot
  !
  SAVE
  !
  CHARACTER(LEN=256) :: flasy ! the name of file with the asymptote script

  CHARACTER(LEN=256) :: asymptote_command ! the asymptote command

  LOGICAL :: lasymptote ! if .true. asymptote is called within the thermo_pw
                        ! code

END MODULE control_asy

MODULE thermo_sym
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the symmetry plot
  !
  SAVE
  !
  INTEGER :: laue    ! The laue class of the solid
  INTEGER :: code_group_save ! the code of the point group of the unperturbed
                             ! solid
  LOGICAL :: ibrav_group_consistent ! if .true., point group and bravais lattice are
                             ! consistent 
  INTEGER :: fft_fact(3)     ! the factors that must be inside the fft
  INTEGER :: sg_number       ! the number of the space_group
  INTEGER :: aux_sg          ! the orientation of the space group

END MODULE thermo_sym

MODULE efermi_plot
  USE kinds,  ONLY : DP

  SAVE

  INTEGER :: n1, n2     ! the number of k points of the 2d plot

  REAL(DP) :: kxmin, kxmax, kymin, kymax ! the size of the plotted region

  LOGICAL, ALLOCATABLE :: has_ef(:) ! say which bands cross ef
END MODULE efermi_plot
