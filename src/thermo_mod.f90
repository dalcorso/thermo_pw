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

  INTEGER  :: ngeo                      ! number of different geometries 
                                        ! calculated in this run

  REAL(DP), ALLOCATABLE :: alat_geo(:),     &   ! the lattice constant at
                                                ! each geometry
                           omega_geo(:),    &   ! volume (geometry)
                           celldm_geo(:,:), &   ! The celldm for each geometry
                           energy_geo(:),   &   ! The total energy at each
                                                ! geometry
                           forces_geo(:,:,:), & ! the forces on atoms at each
                                                ! geometry
                           stress_geo(:,:,:)    ! the stress at each 
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
  REAL(DP)     :: vmin, b0, b01, emin   ! the minimum of the murnaghan at T=0
  REAL(DP) :: vmin_input, vmax_input, deltav   ! plot the fitted total energy
                                               ! and pressure from vmin_input
                                               ! to vnax_input in steps of i
                                               ! deltav
  INTEGER :: nvol                              ! the number of volumes for
                                               ! the plot
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

END MODULE ph_freq_anharmonic

MODULE grun_anharmonic
  USE kinds, ONLY: DP
  USE ph_freq_module, ONLY : ph_freq_type
  !
  !  This module contains all the quantities calculated from the gruneisen
  !  parameters
  !
  SAVE
  TYPE(ph_freq_type) :: ph_grun        ! the gruneisen parameters
                                       ! on a mesh
  REAL(DP), ALLOCATABLE :: betab(:)    ! volume thermal expansion multiplied
                                       ! by the bulk modulus
  REAL(DP), ALLOCATABLE :: grun_gamma_t(:)  ! average gruneisen parameter

END MODULE grun_anharmonic

MODULE ifc
  !
  USE kinds, ONLY: DP
  SAVE
  REAL(DP), ALLOCATABLE :: frc(:,:,:,:,:,:,:), zeu(:,:,:), &
               m_loc(:,:)
  ! frc : interatomic force constants in real space
  ! tau : atomic positions for the original cell
  ! zeu : effective charges for the original cell
  ! m_loc: the magnetic moments of each atom
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
  LOGICAL :: lconv_ke_test=.FALSE.! if .true. this writes the ke test on file
  LOGICAL :: lconv_nk_test=.FALSE.! if .true. this writes the k-point on file
  LOGICAL :: lph=.FALSE.    ! if .true. must calculate phonon
  LOGICAL :: ldos=.FALSE.   ! if .true. the phonon dos is calculated
  LOGICAL :: ltherm=.FALSE. ! if .true. the thermodynamical properties are
                            ! calculated
  LOGICAL :: lq2r=.FALSE.   ! if .true. the interatomic force constants are calculated
  LOGICAL :: lmatdyn=.FALSE.     ! if .true. the phonon are interpolated
  
  CHARACTER(LEN=256) :: outdir_thermo, fildyn_thermo, flfrc, &
                        flfrq, fldos, fltherm, flanhar, flevdat, &
                        filband, flkeconv, flnkconv, flgrun
  INTEGER :: spin_component
  !
END MODULE control_thermo

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
  INTEGER, ALLOCATABLE  :: wqaux(:)       ! the number of points per line
  INTEGER               :: disp_nqs       ! total number of points to compute
  REAL(DP), ALLOCATABLE :: disp_q(:,:), disp_wq(:)  ! q path for interpolated
  INTEGER            :: npk_label         ! number of label points
  INTEGER            :: nqaux             ! number of auxiliary points
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
  CHARACTER(LEN=256) :: flpband ! the name of the output file
  REAL(DP) :: emin_input     ! minimum energy of the plot (eV)
  REAL(DP) :: emax_input     ! maximum energy of the plot (eV)

END MODULE control_bands

MODULE control_grun

  USE kinds, ONLY: DP
  SAVE

  CHARACTER(LEN=256) :: flpgrun ! the name of the output file with the 
                                ! gruneisen parameters in a format 
                                ! readable by gnuplot
                                !
  REAL(DP) :: grunmin_input, &  ! minimum value for gruneisen plot
              grunmax_input     ! maximum value for gruneisen plot

END MODULE control_grun

MODULE control_gnuplot
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the band plot
  !
  SAVE

  CHARACTER(LEN=256) :: flgnuplot ! the name of file with the gnuplot script
  CHARACTER(LEN=256) :: flpsmur   ! the name of the postscript file with E(V)
  CHARACTER(LEN=256) :: flpsband  ! the name of the output postscript file
  CHARACTER(LEN=256) :: flpsdisp  ! the name of the output postscript file
  CHARACTER(LEN=256) :: flpsdos ! the name of the postscript file with dos
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
  CHARACTER(LEN=256) :: gnuplot_command ! the gnuplot command

  LOGICAL :: lgnuplot    ! set to false not to use gnuplot

END MODULE control_gnuplot
