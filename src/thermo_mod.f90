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
                               ! quadratic function of dimension up to 6
  REAL(DP) :: celldm0(6)  ! the minimum celldm

END MODULE control_mur

MODULE thermodynamics
  !
  USE kinds, ONLY: DP
  USE phdos_module, ONLY : phdos_type
  !
  ! ... The variables needed to save the thermodynmac quantities 
  !     calculated from the phonon dos
  !
  SAVE

  TYPE(phdos_type), ALLOCATABLE :: phdos_save(:) ! phdos for each geometry
                                             ! geometry on a uniform mesh

  REAL(DP), ALLOCATABLE :: ph_ener(:,:)      ! phonon total energy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_free_ener(:,:) ! phonon free_energy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_entropy(:,:)   ! phonon entropy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_cv(:,:)        ! phonon specific heat, T, geometry

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
  REAL(DP), ALLOCATABLE :: phf_cv(:,:)        ! phonon specific heat, T, geometry

  REAL(DP), ALLOCATABLE :: vminf_t(:), b0f_t(:), b01f_t(:), free_e_minf_t(:) 
                           ! the parameters of the minimum of the  
                           ! free energy at each temperature

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
  REAL(DP), ALLOCATABLE :: cp_t(:)     ! isobaric specific heat (T)
  REAL(DP), ALLOCATABLE :: cv_t(:)     ! isocoric specific heat (T)
  REAL(DP), ALLOCATABLE :: b0_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: alpha_t(:)  ! linear thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: beta_t(:)   ! volume thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: gamma_t(:)  ! average gruneisen parameter

  REAL(DP), ALLOCATABLE :: vmin_t(:), b0_t(:), b01_t(:), free_e_min_t(:) 
                           ! the parameters of the minimum of the  
                           ! free energy calculated from dos at each temperature
  REAL(DP), ALLOCATABLE :: celldm_t(:,:) ! the lattice parameters as a 
                           ! function of temperature
  REAL(DP), ALLOCATABLE :: alpha_anis_t(:,:)  ! thermal expansion tensor

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
  REAL(DP), ALLOCATABLE :: cpf_t(:)     ! isobaric specific heat (T)
  REAL(DP), ALLOCATABLE :: cvf_t(:)     ! isocoric specific heat (T)
  REAL(DP), ALLOCATABLE :: b0f_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: alphaf_t(:)  ! linear thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: betaf_t(:)   ! volume thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: gammaf_t(:)  ! average gruneisen parameter

  REAL(DP), ALLOCATABLE :: vminf_t(:), b0f_t(:), b01f_t(:), free_e_minf_t(:) 
                           ! the parameters of the minimum of the  
                           ! free energy  
  REAL(DP), ALLOCATABLE :: celldmf_t(:,:) ! the lattice parameters as a 
                           ! function of temperature
  REAL(DP), ALLOCATABLE :: alphaf_anis_t(:,:)  ! thermal expansion tensor

END MODULE ph_freq_anharmonic

MODULE grun_anharmonic
  USE kinds, ONLY: DP
  USE ph_freq_module, ONLY : ph_freq_type
  !
  !  This module contains all the quantities calculated from the gruneisen
  !  parameters
  !
  SAVE
  REAL(DP), ALLOCATABLE :: betab(:)    ! volume thermal expansion multiplied
                                       ! by the bulk modulus
  REAL(DP), ALLOCATABLE :: grun_gamma_t(:)  ! average gruneisen parameter

  REAL(DP), ALLOCATABLE :: poly_grun(:,:,:) ! For each band and each q point
                                        ! these are the coefficients of 
                                        ! polynomial which fit the frequency
                                        ! as a function of volume
  INTEGER :: poly_order                 ! order of the polynomial + 1 

END MODULE grun_anharmonic

MODULE ifc
  !
  USE kinds, ONLY: DP
  SAVE
  REAL(DP), ALLOCATABLE :: frc(:,:,:,:,:,:,:), zeu(:,:,:), &
               m_loc(:,:), wscache(:,:,:,:,:)
  ! frc : interatomic force constants in real space
  ! tau : atomic positions for the original cell
  ! zeu : effective charges for the original cell
  ! m_loc: the magnetic moments of each atom
  ! wscache: the weight of each q point 
  !
  CHARACTER(LEN=3), ALLOCATABLE :: atm(:)
  !
  INTEGER :: nq1_d, nq2_d, nq3_d ! grid for phonon dos
  INTEGER :: ndos_input          ! number of points in the dos plot
  REAL(DP) :: freqmin, freqmax   ! dos minimum and maximun frequency 
                                 ! at this geometry
  REAL(DP) :: freqmin_input, freqmax_input   ! dos minimum and maximun frequency 
                                 ! be given as input (in cm^{-1})
  REAL(DP) :: deltafreq          ! delta frequency in the dos plot. (in cm^{-1})
  CHARACTER(LEN=10) :: zasr      ! the type of asr
  ! 
END MODULE ifc

MODULE thermo_priority
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the job priority
  !
  SAVE
  INTEGER, ALLOCATABLE :: npriority(:)  ! the number of jobs that must be
                                        ! calculated before this one
  INTEGER, ALLOCATABLE :: priority(:,:) ! the priority of the calculation
  !
  INTEGER :: max_priority
  !
END MODULE thermo_priority

MODULE control_thermo
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the type of calculation to do
  !
  SAVE
  !
  LOGICAL, ALLOCATABLE :: lpwscf(:),  & ! if .true. this work requires a scf calc.
                          lbands(:),  & ! if .true. this work requires a band
                          lberry(:),  & ! if .true. this work requires 
                                        ! a berry_phase calculation
                          lstress(:), & ! if .true. this work computes stress
                          lphonon(:)    ! if .true. this work requires a phonon

  LOGICAL :: read_paths      ! the paths for dispersion are read from input
  LOGICAL :: lev_syn_1=.FALSE. ! if .true. must calculate the murnaghan
                               ! at T=0
  LOGICAL :: lev_syn_2=.FALSE. ! if .true. must calculate the murnaghan
                               ! at all T
  LOGICAL :: lpwscf_syn_1=.FALSE. ! if .true. must calculate the syncronous pw
                                  ! for scf calculation
  LOGICAL :: lbands_syn_1=.FALSE. ! if .true. must calculate the syncronous pw
                                  ! for nscf calculation
  LOGICAL :: lpart2_pw=.FALSE.    ! if .true. in the second part makes
                                  ! also pw calculations           
  LOGICAL :: with_eigen=.FALSE.   ! save phonon frequencies and eigenvectors
                               !
  LOGICAL :: lconv_ke_test=.FALSE.! if .true. this writes the ke test on file
  LOGICAL :: lconv_nk_test=.FALSE.! if .true. this writes the k-point on file
  LOGICAL :: lelastic_const=.FALSE. ! if .true. compute elastic constants
  LOGICAL :: lpiezoelectric_tensor=.FALSE. ! if .true. compute the piezoelectric 
                                  ! tensor
  LOGICAL :: lpolarization=.FALSE. ! if .true. compute the piezoelectric 
  LOGICAL :: lph=.FALSE.    ! if .true. must calculate phonon
  LOGICAL :: ldos=.FALSE.   ! if .true. the phonon dos is calculated
  LOGICAL :: ltherm=.FALSE. ! if .true. the thermodynamical properties are
                            ! calculated
  LOGICAL :: after_disp=.FALSE. ! if .true. dynamical matrix files are supposed
                            ! to be already on disk. fildyn must be read
                            ! from thermo_control.
  LOGICAL :: lq2r=.FALSE.   ! if .true. the interatomic force constants 
                            ! are calculated
  LOGICAL :: lmatdyn=.FALSE.  ! if .true. the phonon are interpolated

  INTEGER :: spin_component
  !
  CHARACTER(LEN=256) :: outdir_thermo ! the outdir read from the input
  !
END MODULE control_thermo

MODULE control_elastic_constants
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the calculation of the elastic 
  !     constants
  !
  SAVE
  REAL(DP) :: delta_epsilon
  REAL(DP) :: at_save(3,3)
  REAL(DP), ALLOCATABLE :: tau_save(:,:)
  REAL(DP), ALLOCATABLE :: rot_mat(:,:,:) ! rotation matrix between the
                                          ! cartesian coordinates of the
                                          ! strained and unstrained cell  
  REAL(DP), ALLOCATABLE :: aap_mat(:,:,:) ! possible change of the definition
                                          ! of the direct lattice vectors
                                          ! a (old) in terms of a' (new)
  REAL(DP), ALLOCATABLE :: apa_mat(:,:,:) ! possible change of the definition
                                          ! of the direct lattice vectors
                                          ! a' (new) in terms of a (old)
  REAL(DP) :: omega0            ! for elastic constants not dependent 
                                ! on temperature this is the unperturbed 
                                ! volume
  INTEGER :: ngeo_strain        ! number of strain configurations

  LOGICAL :: frozen_ions        ! if .true. compute the elastic constant 
                                ! keeping the ions frozen at the strained
                                ! positions

  CHARACTER(LEN=12) :: elastic_algorithm ! can be standard, advanced or energy
                                ! it chooses the routines to use to 
                                ! calculate elastic constants.

END MODULE control_elastic_constants

MODULE control_piezoelectric_tensor
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the calculation of the elastic 
  !     constants
  !
  SAVE

  LOGICAL :: nosym_save
  !
END MODULE control_piezoelectric_tensor
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
  INTEGER,  ALLOCATABLE :: nk_test(:) ! the nk to calculate
  REAL(DP) :: deltake               ! the interval between kinetic energies
  INTEGER  :: nke                   ! the number of kinetic energies
  REAL(DP) :: deltakeden            ! the interval for kinetic energy for density
  INTEGER  :: nkeden                ! the number of kinetic energies for densiy
  INTEGER  :: nnk                   ! the number of nk values to test k points
  INTEGER  :: deltank               ! the step between nk values to test k points
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
  CHARACTER(LEN=10) :: point_label_type   ! type of point labels
  LOGICAL :: q_in_band_form
  LOGICAL :: q_in_cryst_coord
  LOGICAL :: q2d
                                                    ! phonon
  !
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
              sur_thr    ! minimum percentage of charge to be a surface state

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


MODULE control_grun

  USE kinds, ONLY: DP
  SAVE

                                 
  REAL(DP) :: grunmin_input, &  ! minimum value for gruneisen plot
              grunmax_input     ! maximum value for gruneisen plot

  REAL(DP) :: volume_ph         ! volume at which the phonon frequencies
                                ! and the Gruneisen parameters are 
                                ! interpolated

  REAL(DP) :: temp_ph           ! when volume_ph=0.0_DP we use the volume
                                ! corresponding to this temperature
END MODULE control_grun

MODULE control_pwrun

  USE kinds, ONLY: DP
  SAVE

  INTEGER  :: nr1_save, nr2_save, nr3_save  ! save the fft dimensions

  REAL(DP) :: celldm_save(6)   ! save the crystal parameters
  INTEGER  :: ibrav_save       ! save the Bravais lattice
  LOGICAL  :: do_punch=.TRUE.  ! set this variable to .FALSE. if pw has
                               ! not to save the punch files.

END MODULE control_pwrun

MODULE control_energy_plot

  USE kinds, ONLY: DP
  SAVE

  INTEGER :: ncontours
  REAL(DP), ALLOCATABLE :: ene_levels(:)
  CHARACTER(LEN=12), ALLOCATABLE :: color_levels(:)

END MODULE control_energy_plot

MODULE control_quadratic_energy

  USE kinds, ONLY: DP
  SAVE

  INTEGER :: degree              ! number of degrees of freedom
  INTEGER :: nvar                ! number of variables of the polynomial fit
  REAL(DP), ALLOCATABLE :: hessian_v(:,:), &   ! hessian eigenvectors
                           hessian_e(:),   &   ! hessian eigenvalues
                           x_pos_min(:),   &   ! coordinates of the minimum
                           coeff(:),       &   ! coefficients of quadratic fit
                           coeff_t(:,:)        ! coefficients at each
                                               ! temperature

END MODULE control_quadratic_energy

MODULE control_gnuplot
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the gnuplot. The files that
  !     contains the gnuplot scripts
  !
  SAVE

  CHARACTER(LEN=256) :: flgnuplot ! the name of file with the gnuplot script
  CHARACTER(LEN=256) :: gnuplot_command ! the gnuplot command

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
  CHARACTER(LEN=256) :: fltherm ! file with the harmonic thermodynamic 
  CHARACTER(LEN=256) :: flanhar ! file with the anharmonic quantities
  CHARACTER(LEN=256) :: flevdat ! file with data for ev.x 
  CHARACTER(LEN=256) :: flenergy  ! the name of the file with the energy
                                  ! suited for gnuplot contour plots
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
  CHARACTER(LEN=256) :: flpstherm ! the name of the postscript file with 
                                  ! thermodynamic quantities
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
END MODULE postscript_files

MODULE internal_files_names
  USE kinds,  ONLY : DP
  !
  ! ... The names of the files before adding the geometry information
  !
  SAVE

  CHARACTER(LEN=256) :: fildyn_thermo    ! dynamical matrix file
  CHARACTER(LEN=256) :: flfrc_thermo     ! interatomic force constants file
  CHARACTER(LEN=256) :: flfrq_thermo     ! frequencies files
  CHARACTER(LEN=256) :: fldos_thermo     ! dos file
  CHARACTER(LEN=256) :: fltherm_thermo   ! thermodynamical quantities files
  CHARACTER(LEN=256) :: flpband_thermo   ! dispersions in gnuplot format

  CHARACTER(LEN=256) :: flgnuplot_thermo  ! file with gnuplot script

  CHARACTER(LEN=256) :: flpsdisp_thermo   ! postscript file with dispersions
  CHARACTER(LEN=256) :: flpsdos_thermo    ! postscript file with dos
  CHARACTER(LEN=256) :: flpstherm_thermo  ! postscript file with thermodynamical
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

END MODULE thermo_sym
