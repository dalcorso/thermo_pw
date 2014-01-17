!
! Copyright (C) 2013 - 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for the phonon program
!
MODULE thermo_mod
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the type of calculations to do
  !
  SAVE

  CHARACTER(LEN=30) :: what             ! the type of calculation. See README.
  REAL(DP)          :: vmin, b0, emin   ! the minimum of the murnaghan
  REAL(DP), ALLOCATABLE :: vmin_t(:), b0_t(:), free_e_min_t(:) ! the minimum 
                                      ! of the murnaghan at each temperature
                                      !
  REAL(DP), ALLOCATABLE :: alat_geo(:), energy_geo(:)
  !
END MODULE thermo_mod

MODULE thermodynamics
  !
  USE kinds, ONLY: DP
  USE phdos_module, ONLY : phdos_type
  SAVE

  TYPE(phdos_type), ALLOCATABLE :: phdos_save(:) ! phdos for each geometry

  REAL(DP), ALLOCATABLE :: temp(:,:)         ! temperature (K), n, geometry
  REAL(DP), ALLOCATABLE :: ph_ener(:,:)      ! phonon total energy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_free_ener(:,:) ! phonon free_energy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_entropy(:,:)   ! phonon entropy, T, geometry
  REAL(DP), ALLOCATABLE :: ph_cv(:,:)        ! phonon specific heat, T, geometry
  REAL(DP), ALLOCATABLE :: omegav(:)         ! volume (geometry)
  REAL(DP) :: tmin, tmax                ! maximum and minimum temperature (K)
  REAL(DP) :: deltat                    ! delta T
  INTEGER  :: ntemp                     ! number of temperatures
  INTEGER  :: ngeo                      ! number of geometries used for
                                        ! computing anharmonic properties
END MODULE thermodynamics

MODULE anharmonic
  USE kinds, ONLY: DP
  REAL(DP), ALLOCATABLE :: cp_t(:)     ! isobaric specific heat (T)
  REAL(DP), ALLOCATABLE :: cv_t(:)     ! isocoric specific heat (T)
  REAL(DP), ALLOCATABLE :: b0_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: alpha_t(:)  ! linear thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: beta_t(:)   ! volume thermal expansion coefficient

END MODULE anharmonic

MODULE ifc
  !
  USE kinds, ONLY: DP
  SAVE
  REAL(DP), ALLOCATABLE :: frc(:,:,:,:,:,:,:), tau_blk(:,:),  zeu(:,:,:), &
               m_loc(:,:)
  ! frc : interatomic force constants in real space
  ! tau_blk : atomic positions for the original cell
  ! zeu : effective charges for the original cell
  ! m_loc: the magnetic moments of each atom
  INTEGER, ALLOCATABLE  :: ityp_blk(:)
  ! ityp_blk : atomic types for each atom of the original cell
  !
  CHARACTER(LEN=3), ALLOCATABLE :: atm(:)
  !
  INTEGER :: disp_nqs
  REAL(DP), ALLOCATABLE :: disp_q(:,:), disp_wq(:)  ! q path for interpolated
                                                    ! phonon
  INTEGER :: nq1_d, nq2_d, nq3_d ! grid for phonon dos
  INTEGER :: ndos                ! number of points in the dos plot
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
  LOGICAL :: lev_syn_2=.FALSE. ! if .true. must calculate the murnaghan
  LOGICAL :: lpwscf_syn_1=.FALSE. ! if .true. must calculate the syncronous pw
                                  ! for scf calculation
  LOGICAL :: lbands_syn_1=.FALSE. ! if .true. must calculate the syncronous pw
                                  ! for nscf calculation
  LOGICAL :: lph=.FALSE. ! if .true. must calculate phonon
  LOGICAL :: ldos        ! if .true. the phonon dos is calculated
  LOGICAL :: ltherm      ! if .true. the thermodynamical properties are
                         ! calculated
  LOGICAL :: lq2r        ! if .true. the interatomic force constants are calculated
  LOGICAL :: lmatdyn     ! if .true. the phonon are interpolated
  
  CHARACTER(LEN=256) :: outdir_thermo, fildyn_thermo, flfrc, &
                        flfrq, fldos, fltherm, flanhar, flevdat, &
                        filband
  INTEGER :: spin_component
  !
END MODULE control_thermo

MODULE control_paths
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the paths
  !
  SAVE
  !
  REAL(DP), ALLOCATABLE :: xqaux(:,:)     ! the initial and final points 
  INTEGER, ALLOCATABLE :: wqaux(:)        ! the number of points per line
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
  INTEGER :: nbnd_bands                   ! number of bands for band calculation
  !
END MODULE control_paths

MODULE control_bands
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the band plot
  !
  SAVE
  !
  CHARACTER(LEN=256) :: flpband ! the name of the output file
  REAL(DP) :: emin_input     ! minimum energy of the plot (eV)
  REAL(DP) :: emax_input     ! maximum energy of the plot (eV)

END MODULE control_bands

MODULE control_gnuplot
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the band plot
  !
  SAVE

  CHARACTER(LEN=256) :: flgnuplot ! the name of file with the gnuplot script
  CHARACTER(LEN=256) :: flpsband  ! the name of the output postscript file
  CHARACTER(LEN=256) :: flpsdisp  ! the name of the output postscript file
  CHARACTER(LEN=256) :: flpsdos ! the name of the postscript file with dos
  CHARACTER(LEN=256) :: flpstherm ! the name of the postscript file with 
                                  ! thermodynamic quantities
  CHARACTER(LEN=256) :: flpsanhar ! the name of the postscript file with 
                                  ! anharmonic quantities

END MODULE control_gnuplot
