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
!----------------------------------------------------------------------------
MODULE thermo_mod
!----------------------------------------------------------------------------
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
                           stress_geo(:,:,:), & ! the stress at each 
                                                ! geometry
                           ef_geo(:)            ! save the fermi energy in
                                                ! metals

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
  REAL(DP) :: celldm_p1(6), celldm_m1(6)

  LOGICAL :: lcubic=.FALSE.                     ! .TRUE. when the bravais
                                                ! lattice is cubic

END MODULE thermo_mod

!----------------------------------------------------------------------------
MODULE temperature
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the temperature
  !
  SAVE

  REAL(DP), ALLOCATABLE :: temp(:)      ! temperature (K), n

  REAL(DP) :: tmin, tmax                ! maximum and minimum temperature (K)
  REAL(DP) :: deltat                    ! delta T
  INTEGER  :: ntemp                     ! number of temperatures
!
!  Some plots can be done as a function of pressure or of volume 
!  at several temperatures
!  These variables control which temperatures are plotted
!
  INTEGER  :: ntemp_plot                ! how many pressures are plot
  REAL(DP), ALLOCATABLE :: temp_plot(:) ! the values of the temperature 
                                        !                       (in kbar)
  INTEGER, ALLOCATABLE :: itemp_plot(:) ! the index in the list of temperatures
                                        ! generated by tmin and tmax
  INTEGER :: itemp300                   ! the index of the temperature closest
                                        ! to 300 K.
  REAL(DP), ALLOCATABLE :: temp_sigma(:), & ! Used when the total energy is
                           sigma_ry(:)      ! computed for serveral sigma 
                                            ! values
 
END MODULE temperature

!----------------------------------------------------------------------------
MODULE control_pressure
!----------------------------------------------------------------------------
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

!
!   Some plots can be done as a function of pressure. These variables
!   control these plots
!
  REAL(DP) :: pmin        ! minimal pressure in the pressure-volume plot
                          ! (lmurn=.FALSE.)
  REAL(DP) :: pmax        ! maximum pressure in the pressure-volume plot
                          ! (lmurn=.FALSE.)
  REAL(DP) :: deltap      ! step interval in the plot
  INTEGER  :: npress      ! number of computed pressures

  REAL(DP), ALLOCATABLE :: press(:)
!
!  Some plots can be done as a function of temperature at several pressures
!  These variables control which pressures are plotted
!
  INTEGER  :: npress_plot         ! how many pressures are plot
  REAL(DP), ALLOCATABLE :: press_plot(:) ! the values of the pressure (in kbar)
  INTEGER, ALLOCATABLE :: ipress_plot(:) ! the index in the list of pressures
                                         ! generated by pmin and pmax
END MODULE control_pressure

!----------------------------------------------------------------------------
MODULE control_mur
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the calculation of the murnaghan
  !     equation
  !
  SAVE
  REAL(DP) :: vmin, b0, b01, b02, emin   ! the minimum of the murnaghan at T=0
  LOGICAL :: lmurn             ! if .TRUE. makes a Murnaghan
                               ! fit of the energy as a function of the volume
                               ! otherwise makes a fit of the energy as a 
                               ! function of the celldm parameters, with a
                               ! quadratic or quartic polynomial of dimension 
                               ! up to 6

  REAL(DP) :: omegap0          ! omega at p=0. Estimated in all cases also
                               ! when enthalpy is minimized at finite pressure.

  REAL(DP), ALLOCATABLE :: p0(:) ! pressure versus volume at T=0 K

END MODULE control_mur

!----------------------------------------------------------------------------
MODULE control_mur_p
!----------------------------------------------------------------------------
   USE kinds,  ONLY : DP
   !
   ! ... The variables needed to control the calculation of the Murnaghan
   !     equation at any pressure
   !
   SAVE
   REAL(DP), ALLOCATABLE :: vmin_p(:), &   ! The minimum of the murnaghan
                            b0_p(:),   &   ! the bulk modulus
                            b01_p(:),  &   ! the derivative of the bulk modulus
                            b02_p(:),  &   ! the second derivative of the 
                                           ! bulk modulus
                            emin_p(:)      ! the minimum enthalpy 
                                           ! as a function of pressure
END MODULE control_mur_p

!----------------------------------------------------------------------------
MODULE control_vol
!----------------------------------------------------------------------------
   USE kinds,  ONLY : DP
   !
   ! ... The variables needed to control the calculation of 
   !     equation at any pressure
   !
   SAVE

   REAL(DP) :: vmin_input, vmax_input, deltav   ! plot the fitted total energy
                                                ! and pressure from vmin_input
                                                ! to vnax_input in steps of i
                                                ! deltav
   INTEGER :: nvol                              ! the number of volumes for
                                                ! the plot
   INTEGER  :: nvol_plot                ! how many pressures are plot
                                        !
   INTEGER, ALLOCATABLE :: ivol_plot(:) ! the index in the list of volume
                                        ! This must be in the list of ngeo
                                        ! computed geometries
END MODULE control_vol

!----------------------------------------------------------------------------
MODULE uniform_pressure
!----------------------------------------------------------------------------
   USE kinds,  ONLY : DP
   USE polynomial, ONLY : poly2, poly4
   !
   ! ... The variables needed to save the crystal parameters as a 
   !     function of pressure in the anisotropic case when lmurn=.FALSE. 
   !
   SAVE
   REAL(DP), ALLOCATABLE :: omega_p(:)    ! npress. The volume.
   REAL(DP), ALLOCATABLE :: bm_p(:)       ! npress. The bulk modulus
                                          !         from the derivative of
                                          !         omega_p with respect to p  
   REAL(DP), ALLOCATABLE :: celldm_p(:,:) ! (6 npress) the crystal parameters
   REAL(DP), ALLOCATABLE :: el_cons_p(:,:,:), & ! elastic constants
                            el_comp_p(:,:,:), & ! elastic compliances
                            b0ec_p(:),        & ! bulk modulus
                            macro_el_p(:,:),  & ! bulk modulus, shear modulus
                                                ! etc.
                            v_p(:,:),         & ! sound velocities
                            density_p(:)        ! density
                                                !
   TYPE(poly2), ALLOCATABLE :: p2_p(:)          ! the polynomial interpolating
                                                ! enthalpy (quadratic)
   TYPE(poly4), ALLOCATABLE :: p4_p(:)          ! (quartic)

   TYPE(poly2) :: p2_p_p1       ! the polynomial interpolating
                                ! enthalpy (quadratic)
   TYPE(poly4) :: p4_p_p1       ! (quartic) at pressure+deltap

   TYPE(poly2) :: p2_p_m1       ! the polynomial interpolating
                                ! enthalpy (quadratic)
   TYPE(poly4) :: p4_p_m1       ! (quartic) at pressure-deltap

END MODULE uniform_pressure
!
!----------------------------------------------------------------------------
MODULE thermodynamics
!----------------------------------------------------------------------------
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
  REAL(DP), ALLOCATABLE :: ph_t_debye(:,:)   ! debye temperature, T, geometry
  REAL(DP), ALLOCATABLE :: ph_e0(:)          ! zero point energy, geometry
  REAL(DP), ALLOCATABLE :: ph_ce(:,:)        ! phonon heat capacity, T, geometry
  REAL(DP), ALLOCATABLE :: ph_b_fact(:,:,:,:,:)! atomic B factor

END MODULE thermodynamics

!----------------------------------------------------------------------------
MODULE ph_freq_thermodynamics
!----------------------------------------------------------------------------
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
  REAL(DP), ALLOCATABLE :: phf_t_debye(:,:)   ! debye temperature, T, geometry
  REAL(DP), ALLOCATABLE :: phf_e0(:)          ! zero point energy, geometry
  REAL(DP), ALLOCATABLE :: phf_ce(:,:)        ! phonon specific heat, T, geometry
  REAL(DP), ALLOCATABLE :: phf_b_fact(:,:,:,:,:) ! atomic B factor

END MODULE ph_freq_thermodynamics

!----------------------------------------------------------------------------
MODULE el_thermodynamics
!----------------------------------------------------------------------------
  !
  USE kinds, ONLY: DP
  !
  ! ... The variables needed to save the electronic thermodynamic quantities 
  !     calculated from the electron dos at several geometries
  !
  SAVE

  REAL(DP), ALLOCATABLE :: el_ener(:,:)      ! electronic total energy, &
                                             ! T, geometry
  REAL(DP), ALLOCATABLE :: el_free_ener(:,:) ! electronic free_energy, T, &
                                             ! geometry
  REAL(DP), ALLOCATABLE :: el_entr(:,:)      ! electronic entropy, T, geometry
  REAL(DP), ALLOCATABLE :: el_mu(:,:)        ! electronic chemical potential, 
                                             ! T, geometry
  REAL(DP), ALLOCATABLE :: el_ce(:,:)        ! electronic heat capacity, &
                                             ! T, geometry
END MODULE el_thermodynamics
!
!----------------------------------------------------------------------------
MODULE emp_anharmonic
!----------------------------------------------------------------------------
!
!   here all the quantities obtained adding the free energy of the
!   electrons
!
  USE kinds, ONLY: DP
  SAVE
!
!  Interpolated empirical quantities (using phdos)
!
  REAL(DP), ALLOCATABLE :: emp_energy_t(:) ! empirical energy 
  REAL(DP), ALLOCATABLE :: emp_free_energy_t(:) ! empirical free energy
  REAL(DP), ALLOCATABLE :: emp_entropy_t(:) !  empirical entropy
  REAL(DP), ALLOCATABLE :: emp_ce_t(:) ! empirical heat capacity
!
!  Interpolated empirical quantities (using phonon frequencies)
!
  REAL(DP), ALLOCATABLE :: emp_energyf_t(:) ! empirical energy 
  REAL(DP), ALLOCATABLE :: emp_free_energyf_t(:) ! empirical free energy
  REAL(DP), ALLOCATABLE :: emp_entropyf_t(:) !  empirical entropy
  REAL(DP), ALLOCATABLE :: emp_cef_t(:) ! empirical heat capacity
!
!  Empirical quantities interpolated at the temperature dependent geometry 
!  at selected pressures 
!
  REAL(DP), ALLOCATABLE :: emp_free_ener_pt(:,:) ! empirical free energy
  REAL(DP), ALLOCATABLE :: emp_ener_pt(:,:) ! empirical energy
  REAL(DP), ALLOCATABLE :: emp_entr_pt(:,:) ! empirical entropy
  REAL(DP), ALLOCATABLE :: emp_ce_pt(:,:) ! empirical heat capacity at
                                          ! several pressures
!
!  Empirical quantities interpolated at the temperature dependent geometry 
!  at selected pressures 
!
  REAL(DP), ALLOCATABLE :: emp_free_enerf_pt(:,:) ! empirical free energy
  REAL(DP), ALLOCATABLE :: emp_enerf_pt(:,:) ! empirical energy
  REAL(DP), ALLOCATABLE :: emp_entrf_pt(:,:) ! empirical entropy
  REAL(DP), ALLOCATABLE :: emp_cef_pt(:,:)   ! empirical heat capacity at
                                             ! several pressures
!
!  Empirical quantities interpolated at the pressure dependent geometries
!  at selected temperatures
!
  REAL(DP), ALLOCATABLE :: emp_ce_ptt(:,:) ! empirical heat capacity at
                                           ! several temperatures
  REAL(DP), ALLOCATABLE :: emp_cef_ptt(:,:) ! empirical heat capacity at
                                           ! several temperatures
!
END MODULE emp_anharmonic

!----------------------------------------------------------------------------
MODULE anharmonic
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the free energy fitted by an equation of state.
!   Quantities obtained from the free energy calculated from
!   the phonon dos.
!
  USE kinds, ONLY: DP
  USE polynomial,  ONLY : poly1, poly2, poly3, poly4
  SAVE
!
!  The parameters of the interpolation at each temperature
!
  REAL(DP), ALLOCATABLE :: vmin_t(:)   ! minimum volume at each T
  REAL(DP), ALLOCATABLE :: b0_t(:)     ! bulk modulus at each T from EOS
  REAL(DP), ALLOCATABLE :: b0_ec_t(:)  ! bulk modulus at each T from EC
  REAL(DP), ALLOCATABLE :: b01_t(:)    ! pressure derivative of b0 at each T
  REAL(DP), ALLOCATABLE :: b02_t(:)    ! second pressure derivative of b0
                                       ! at each T
  REAL(DP), ALLOCATABLE :: free_e_min_t(:) ! total free energy at the minimum
                                       ! at each T
  REAL(DP), ALLOCATABLE :: a_t(:,:) ! the coefficients of the polynomial
                           ! that interpolates the free energy at each
                           ! temperature
  TYPE(poly1), ALLOCATABLE :: p1t_t(:)    ! coefficients of the polynomial
  TYPE(poly2), ALLOCATABLE :: p2t_t(:)    ! that interpolate the free energy
  TYPE(poly3), ALLOCATABLE :: p3t_t(:)    ! with multidimensional integral.
  TYPE(poly4), ALLOCATABLE :: p4t_t(:)
!
!  The interpolated harmonic quantities at each temperature
!
  REAL(DP), ALLOCATABLE :: ener_t(:)   ! vibrational energy (T)
  REAL(DP), ALLOCATABLE :: free_ener_t(:) ! vibrational free energy (T)
  REAL(DP), ALLOCATABLE :: entropy_t(:)! entropy (T)
  REAL(DP), ALLOCATABLE :: ce_t(:)     ! constant strain heat capacity (T)
  REAL(DP), ALLOCATABLE :: cv_t(:)     ! isocoric heat capacity (T)
  REAL(DP) :: e0_t                     ! zero point energy 
!
!  The calculated anharmonic quantities at each temperature
!
  REAL(DP), ALLOCATABLE :: cp_t(:)     ! isobaric heat capacity 
  REAL(DP), ALLOCATABLE :: b0_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: b0_ec_s(:)  ! constant entropy bulk modulus from ec
  REAL(DP), ALLOCATABLE :: alpha_t(:)  ! linear thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: beta_t(:)   ! volume thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: gamma_t(:)  ! average gruneisen parameter
!
!  Anharmonic quantities for anisotropic thermodynamic
!
  REAL(DP), ALLOCATABLE :: celldm_t(:,:) ! the lattice parameters as a 
                           ! function of temperature
  REAL(DP), ALLOCATABLE :: density_t(:)  ! the density
  REAL(DP), ALLOCATABLE :: alpha_anis_t(:,:)  ! thermal expansion tensor 
                                          ! (Voigt index)
  REAL(DP), ALLOCATABLE :: cpmce_anis(:) ! difference cp-ce computed from
                                         ! elastic constants
  REAL(DP), ALLOCATABLE :: bfact_t(:,:,:)! b factor as a function of 
                                         ! temperature
  REAL(DP), ALLOCATABLE :: csmct_t(:,:,:) ! difference of elastic constants
  REAL(DP), ALLOCATABLE :: bths_t(:,:,:)  ! thermal stress
  REAL(DP), ALLOCATABLE :: ggamma_t(:,:,:)! generalized average gruneisen 
                                       ! parameter
!
!  elastic constants and related quantities
!
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
  REAL(DP), ALLOCATABLE :: debye_macro_el_t(:) ! debye temperature from
                                           ! macro elasticity (isothermal)
  REAL(DP), ALLOCATABLE :: debye_macro_el_s(:) ! debye temperature from
                                           ! macro elasticity (isoentropic)
  REAL(DP), ALLOCATABLE :: el_con_geo_t(:,:,:,:) ! the temperature dependent
                                         ! elastic constants at all geometries
  REAL(DP), ALLOCATABLE :: dyde_t(:,:) !(nstep,ntemp) the derivative of the
                                       ! internal parameter with respect to 
                                       ! strain
!
!  The parameters of the interpolation neglecting the electronic exitation
!  contribution
!
  REAL(DP), ALLOCATABLE :: vmin_noe_t(:)  ! minimum volume at each T 
                                          ! without el.
  REAL(DP), ALLOCATABLE :: celldm_noe_t(:,:)  ! the lattice parameters as a 
                                          ! function of temperature
  REAL(DP), ALLOCATABLE :: density_noe_t(:)  ! the density
  REAL(DP), ALLOCATABLE :: b0_noe_t(:)    ! bulk modulus at each T 
                                          ! without elec.
  REAL(DP), ALLOCATABLE :: b01_noe_t(:)   ! pressure derivative of b0 at each T
                                          ! without electrons
  REAL(DP), ALLOCATABLE :: b02_noe_t(:)   ! second pressure derivative of b0
                                          ! at each T without electrons
  REAL(DP), ALLOCATABLE :: free_e_min_noe_t(:) ! total free energy at 
                                       ! the minimum at each T without el.
  REAL(DP), ALLOCATABLE :: ener_noe_t(:)   ! vibrational energy (T)
  REAL(DP), ALLOCATABLE :: free_ener_noe_t(:) ! vibrational free energy (T)
  REAL(DP), ALLOCATABLE :: entropy_noe_t(:)! entropy (T)
  REAL(DP) :: e0_noe_t                     ! zero point energy 

  REAL(DP), ALLOCATABLE :: a_noe_t(:,:) ! the coefficients of the polynomial
                           ! that interpolates the free energy at each
                           ! temperature
  TYPE(poly1), ALLOCATABLE :: p1t_noe_t(:) ! coefficients of the polynomial
  TYPE(poly2), ALLOCATABLE :: p2t_noe_t(:) ! that interpolate the free energy
  TYPE(poly3), ALLOCATABLE :: p3t_noe_t(:) ! with multidimensional integral.
  TYPE(poly4), ALLOCATABLE :: p4t_noe_t(:)
!
!  The interpolated harmonic quantities neglecting the electronic effect
!
  REAL(DP), ALLOCATABLE :: ce_noe_t(:)     ! constant strain heat capacity (T)
  REAL(DP), ALLOCATABLE :: cv_noe_t(:)     ! isocoric heat capacity (T)
!
!  The anharmonic quantities neglecting the electronic effect
!
  REAL(DP), ALLOCATABLE :: cp_noe_t(:)     ! isobaric heat capacity (T)
  REAL(DP), ALLOCATABLE :: b0_noe_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: beta_noe_t(:)   ! volume thermal expansion
  REAL(DP), ALLOCATABLE :: gamma_noe_t(:)  ! average gruneisen parameter

  LOGICAL :: noelcvg                       ! if .TRUE. the electronic
                                           ! contribution is not added to 
                                           ! cv when computing the average
                                           ! Gruneisen parameter
  REAL(DP), ALLOCATABLE :: celldm_t_p1(:,:)! the celldm at the pressure+dp
                                           ! as a function of T
  REAL(DP), ALLOCATABLE :: celldm_t_m1(:,:)! the celldm at the pressure-dp
                                           ! as a function of T
  REAL(DP), ALLOCATABLE :: celldm_noe_t_p1(:,:)! the celldm at the pressure+dp
                                           ! as a function of T
  REAL(DP), ALLOCATABLE :: celldm_noe_t_m1(:,:)! the celldm at the pressure-dp
                                           ! as a function of T
END MODULE anharmonic

!----------------------------------------------------------------------------
MODULE anharmonic_pt
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the Gibbs free energy fitted by an equation of state
!   plus a polynomial. These variables are computed at all temperatures 
!   for selected pressures.
!
  USE kinds, ONLY: DP
  SAVE
!
!  Parameters of the interpolation
!
  REAL(DP), ALLOCATABLE :: vmin_pt(:,:)  ! the volume at the minimum
  REAL(DP), ALLOCATABLE :: b0_pt(:,:)    ! bulk modulus for all T from EOS
  REAL(DP), ALLOCATABLE :: b0_ec_pt(:,:) ! bulk modulus for all T from EC
  REAL(DP), ALLOCATABLE :: b01_pt(:,:)   ! pressure derivative of b0
  REAL(DP), ALLOCATABLE :: b02_pt(:,:)   ! second pressure derivative of b0
  REAL(DP), ALLOCATABLE :: emin_pt(:,:)  ! Gibbs energy at the minimum
!
!  Interpolated harmonic quantities
!
  REAL(DP), ALLOCATABLE :: free_ener_pt(:,:)  ! the free energy for selected
                                         ! pressures
  REAL(DP), ALLOCATABLE :: ener_pt(:,:)  ! the energy for selected pressures
  REAL(DP), ALLOCATABLE :: entr_pt(:,:)  ! the entropy for selected pressures
  REAL(DP), ALLOCATABLE :: ce_pt(:,:)    ! the constant strain specific
                                         ! heat as a function of temperature
                                         ! at selected pressures
  REAL(DP), ALLOCATABLE :: cv_pt(:,:)    ! the constant volume specific
                                         ! heat as a function of temperature
                                         ! at selected pressures
!
!  The calculated anharmonic quantities
!
  REAL(DP), ALLOCATABLE :: cp_pt(:,:)    ! the isobaric specific heat 
                                         ! at several pressures
  REAL(DP), ALLOCATABLE :: beta_pt(:,:)  ! the thermal expansion as a
                                         ! function of temperature for 
                                         ! selected pressures
  REAL(DP), ALLOCATABLE :: b0_s_pt(:,:)  ! the bulk modulus at constant
                                         ! entropy 
  REAL(DP), ALLOCATABLE :: gamma_pt(:,:) ! the average Gruneisen parameter
!
! Anharmonic quantities for anisotropic thermodynamic
!
  REAL(DP), ALLOCATABLE :: celldm_pt(:,:,:)      ! crystal parameters
  REAL(DP), ALLOCATABLE :: density_pt(:,:)       ! the density
  REAL(DP), ALLOCATABLE :: alpha_anis_pt(:,:,:)  ! the thermal expansion
  REAL(DP), ALLOCATABLE :: cpmce_anis_pt(:,:)    ! the difference cp-ce
  REAL(DP), ALLOCATABLE :: bths_pt(:,:,:,:)      ! the thermal stress
  REAL(DP), ALLOCATABLE :: ggamma_pt(:,:,:,:)    ! the generalized Gruneisen
                                                 ! parameters
  REAL(DP), ALLOCATABLE :: csmct_pt(:,:,:,:) ! isoentropic elastic constants
                                             ! - isothermal ones
  REAL(DP), ALLOCATABLE :: el_cons_pt(:,:,:,:)   ! isothermal elast. cons.
  REAL(DP), ALLOCATABLE :: el_comp_pt(:,:,:,:)   ! isothermal elast. comp.
  REAL(DP), ALLOCATABLE :: el_cons_s_pt(:,:,:,:) ! isoentropic elast. cons.
  REAL(DP), ALLOCATABLE :: el_comp_s_pt(:,:,:,:) ! isoentropic elast. comp.
  REAL(DP), ALLOCATABLE :: macro_el_pt(:,:,:)    ! macroscopic elasticity
  REAL(DP), ALLOCATABLE :: macro_el_s_pt(:,:,:)  ! isoentropic macroscopic ela.
  REAL(DP), ALLOCATABLE :: debye_macro_el_pt(:,:)  ! isothermal debye temper.
  REAL(DP), ALLOCATABLE :: debye_macro_el_s_pt(:,:)  ! adiabatic debye temper.
  REAL(DP), ALLOCATABLE :: v_pt(:,:,:)           ! isothermal sound speed
  REAL(DP), ALLOCATABLE :: v_s_pt(:,:,:)         ! isoentropic sound speed

  REAL(DP), ALLOCATABLE :: celldm_pt_p1(:,:,:)   ! crystal parameters at p+dp
  REAL(DP), ALLOCATABLE :: celldm_pt_m1(:,:,:)   ! crystal parameters at p-dp

END MODULE anharmonic_pt
!----------------------------------------------------------------------------
MODULE ph_freq_anharmonic_pt
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the Gibbs free energy fitted by an equation of state
!   plus a polynomial. These variables are computed at all temperatures 
!   for selected pressures.
!
  USE kinds, ONLY: DP
  SAVE
!
!  Parameters of the interpolation
!
  REAL(DP), ALLOCATABLE :: vminf_pt(:,:)  ! the volume at the minimum
  REAL(DP), ALLOCATABLE :: b0f_pt(:,:)    ! bulk modulus for all T from EOS
  REAL(DP), ALLOCATABLE :: b0f_ec_pt(:,:) ! bulk modulus for all T from EC
  REAL(DP), ALLOCATABLE :: b01f_pt(:,:)   ! pressure derivative of b0
  REAL(DP), ALLOCATABLE :: b02f_pt(:,:)   ! second pressure derivative of b0
  REAL(DP), ALLOCATABLE :: eminf_pt(:,:)  ! Gibbs energy at the minimum
!
!  Interpolated harmonic quantities
!
  REAL(DP), ALLOCATABLE :: free_enerf_pt(:,:)  ! the free energy for selected
                                         ! pressures
  REAL(DP), ALLOCATABLE :: enerf_pt(:,:)  ! the energy for selected pressures
  REAL(DP), ALLOCATABLE :: entrf_pt(:,:)  ! the entropy for selected pressures
  REAL(DP), ALLOCATABLE :: cef_pt(:,:)    ! the constant strain specific
                                         ! heat as a function of temperature
                                         ! at selected pressures
  REAL(DP), ALLOCATABLE :: cvf_pt(:,:)    ! the constant volume specific
                                         ! heat as a function of temperature
                                         ! at selected pressures
!
!  The calculated anharmonic quantities
!
  REAL(DP), ALLOCATABLE :: cpf_pt(:,:)    ! the isobaric specific heat 
                                         ! at several pressures
  REAL(DP), ALLOCATABLE :: betaf_pt(:,:)  ! the thermal expansion as a
                                         ! function of temperature for 
                                         ! selected pressures
  REAL(DP), ALLOCATABLE :: b0f_s_pt(:,:)  ! the bulk modulus at constant
                                         ! entropy 
  REAL(DP), ALLOCATABLE :: gammaf_pt(:,:) ! the average Gruneisen parameter
!
! Anharmonic quantities for anisotropic thermodynamic
!
  REAL(DP), ALLOCATABLE :: celldmf_pt(:,:,:)      ! crystal parameters
  REAL(DP), ALLOCATABLE :: densityf_pt(:,:)       ! the density
  REAL(DP), ALLOCATABLE :: alphaf_anis_pt(:,:,:)  ! the thermal expansion
  REAL(DP), ALLOCATABLE :: cpmcef_anis_pt(:,:)    ! the difference cp-ce
  REAL(DP), ALLOCATABLE :: bthsf_pt(:,:,:,:)      ! the thermal stress
  REAL(DP), ALLOCATABLE :: ggammaf_pt(:,:,:,:)    ! the generalized Gruneisen
                                                 ! parameters
  REAL(DP), ALLOCATABLE :: csmctf_pt(:,:,:,:) ! isoentropic elastic constants
                                             ! - isothermal ones
  REAL(DP), ALLOCATABLE :: el_consf_pt(:,:,:,:)   ! isothermal elast. cons.
  REAL(DP), ALLOCATABLE :: el_compf_pt(:,:,:,:)   ! isothermal elast. comp.
  REAL(DP), ALLOCATABLE :: el_consf_s_pt(:,:,:,:) ! isoentropic elast. cons.
  REAL(DP), ALLOCATABLE :: el_compf_s_pt(:,:,:,:) ! isoentropic elast. comp.
  REAL(DP), ALLOCATABLE :: macro_elf_pt(:,:,:)    ! macroscopic elasticity
  REAL(DP), ALLOCATABLE :: macro_elf_s_pt(:,:,:)  ! isoentropic macroscopic ela.
  REAL(DP), ALLOCATABLE :: vf_pt(:,:,:)           ! isothermal sound speed
  REAL(DP), ALLOCATABLE :: vf_s_pt(:,:,:)         ! isoentropic sound speed
  REAL(DP), ALLOCATABLE :: debye_macro_elf_pt(:,:)  ! isothermal debye temper.
  REAL(DP), ALLOCATABLE :: debye_macro_elf_s_pt(:,:)  ! adiabatic debye temper.

  REAL(DP), ALLOCATABLE :: celldmf_pt_p1(:,:,:)   ! crystal parameters at p+dp
  REAL(DP), ALLOCATABLE :: celldmf_pt_m1(:,:,:)   ! crystal parameters at p-dp

END MODULE ph_freq_anharmonic_pt

!----------------------------------------------------------------------------
MODULE anharmonic_ptt
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the Gibbs free energy fitted by an equation of state
!   plus a polynomial. These variables are computed at all pressures 
!   for selected temperatures.
!
  USE kinds, ONLY: DP
  SAVE
!
!  The parameters of the interpolation
!
  REAL(DP), ALLOCATABLE :: vmin_ptt(:,:)  ! the volume at the minimum for all P
  REAL(DP), ALLOCATABLE :: vmin_ptt_p1(:,:) ! the volume at the minimum 
                                          ! for all P at T+deltat
  REAL(DP), ALLOCATABLE :: vmin_ptt_m1(:,:) ! the volume at the minimum 
                                          ! for all P at T-deltat
  REAL(DP), ALLOCATABLE :: b0_ptt(:,:)    ! the bulk modulus for all P from EOS
  REAL(DP), ALLOCATABLE :: b0_ec_ptt(:,:) ! the bulk modulus for all P from EC
  REAL(DP), ALLOCATABLE :: b01_ptt(:,:)   ! the pressure derivative of b0
  REAL(DP), ALLOCATABLE :: b02_ptt(:,:)   ! the second pressure derivative 
                                          ! of b0
  REAL(DP), ALLOCATABLE :: emin_ptt(:,:)  ! the Gibbs energy at the minimum
  REAL(DP), ALLOCATABLE :: emin_ptt_p1(:,:) ! the Gibbs energy at the minimum
                                          ! for all P at T+deltat
  REAL(DP), ALLOCATABLE :: emin_ptt_m1(:,:) ! the Gibbs energy at the minimum
                                          ! for all P at T-deltat
!
!  The interpolated harmonic quantities
!
  REAL(DP), ALLOCATABLE :: ener_ptt(:,:)  ! the energy
  REAL(DP), ALLOCATABLE :: entr_ptt(:,:)  ! the entropy
  REAL(DP), ALLOCATABLE :: ce_ptt(:,:)    ! the constant strain heat capacity
  REAL(DP), ALLOCATABLE :: cv_ptt(:,:)    ! the constant volume heat capacity
!
!  The anharmonic quantities
!
  REAL(DP), ALLOCATABLE :: beta_ptt(:,:)  ! the volume thermal expansion 
  REAL(DP), ALLOCATABLE :: cp_ptt(:,:)    ! the constant pressure heat 
                                          ! capacity
  REAL(DP), ALLOCATABLE :: b0_s_ptt(:,:)  ! the isoentropic bulk modulus
  REAL(DP), ALLOCATABLE :: gamma_ptt(:,:) ! the average Gruneisen parameter
!
! Anharmonic quantities for anisotropic thermodynamics
!
  REAL(DP), ALLOCATABLE :: celldm_ptt(:,:,:)      ! the crystal parameters 
  REAL(DP), ALLOCATABLE :: celldm_ptt_p1(:,:,:)   ! the crystal parameters at
                                                  ! T+deltat
  REAL(DP), ALLOCATABLE :: celldm_ptt_m1(:,:,:)   ! the crystal parameters at
                                                  ! T-deltat
  REAL(DP), ALLOCATABLE :: alpha_anis_ptt(:,:,:)  ! the thermal expansion 
                                                  ! tensor
  REAL(DP), ALLOCATABLE :: density_ptt(:,:)       ! the density
  REAL(DP), ALLOCATABLE :: cpmce_anis_ptt(:,:)    ! the difference cp-ce
  REAL(DP), ALLOCATABLE :: bths_ptt(:,:,:,:)      ! the thermal stress
  REAL(DP), ALLOCATABLE :: ggamma_ptt(:,:,:,:)    ! the generalized Gruneisen
                                                  ! parameters
  REAL(DP), ALLOCATABLE :: csmct_ptt(:,:,:,:) ! isoentropic elastic constants
                                             ! - isothermal ones
  REAL(DP), ALLOCATABLE :: el_cons_ptt(:,:,:,:)    ! isothermal elast. cons.
  REAL(DP), ALLOCATABLE :: el_comp_ptt(:,:,:,:)    ! isothermal elast. comp.
  REAL(DP), ALLOCATABLE :: el_cons_s_ptt(:,:,:,:)  ! isoentropic elast. cons.
  REAL(DP), ALLOCATABLE :: el_comp_s_ptt(:,:,:,:)  ! isoentropic elast. comp.
  REAL(DP), ALLOCATABLE :: macro_el_ptt(:,:,:)     ! macroscopic elasticity
  REAL(DP), ALLOCATABLE :: macro_el_s_ptt(:,:,:)   ! isoentropic macros. ela.
  REAL(DP), ALLOCATABLE :: debye_macro_el_ptt(:,:) ! isothermal debye temp.
  REAL(DP), ALLOCATABLE :: debye_macro_el_s_ptt(:,:) ! isoentropic debye temp.
  REAL(DP), ALLOCATABLE :: v_ptt(:,:,:)            ! isothermal sound speed
  REAL(DP), ALLOCATABLE :: v_s_ptt(:,:,:)          ! isoentropic sound speed

END MODULE anharmonic_ptt
!----------------------------------------------------------------------------
MODULE ph_freq_anharmonic_ptt
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the Gibbs free energy fitted by an equation of state
!   plus a polynomial. These variables are computed at all pressures 
!   for selected temperatures.
!
  USE kinds, ONLY: DP
  SAVE
!
!  The parameters of the interpolation
!
  REAL(DP), ALLOCATABLE :: vminf_ptt(:,:)! the volume at the minimum for all P
  REAL(DP), ALLOCATABLE :: vminf_ptt_p1(:,:) ! the volume at the minimum 
                                          ! for all P at T+deltat
  REAL(DP), ALLOCATABLE :: vminf_ptt_m1(:,:) ! the volume at the minimum 
                                          ! for all P at T-deltat
  REAL(DP), ALLOCATABLE :: b0f_ptt(:,:)    ! the bulk modulus for all P from EOS
  REAL(DP), ALLOCATABLE :: b0f_ec_ptt(:,:) ! the bulk modulus for all P from EC
  REAL(DP), ALLOCATABLE :: b01f_ptt(:,:)   ! the pressure derivative of b0
  REAL(DP), ALLOCATABLE :: b02f_ptt(:,:)   ! the second pressure derivative 
                                          ! of b0
  REAL(DP), ALLOCATABLE :: eminf_ptt(:,:)  ! the Gibbs energy at the minimum
  REAL(DP), ALLOCATABLE :: eminf_ptt_p1(:,:) ! the Gibbs energy at the minimum
                                          ! for all P at T+deltat
  REAL(DP), ALLOCATABLE :: eminf_ptt_m1(:,:) ! the Gibbs energy at the minimum
                                          ! for all P at T-deltat
!
!  The interpolated harmonic quantities
!
  REAL(DP), ALLOCATABLE :: enerf_ptt(:,:)  ! the energy
  REAL(DP), ALLOCATABLE :: entrf_ptt(:,:)  ! the entropy
  REAL(DP), ALLOCATABLE :: cef_ptt(:,:)    ! the constant strain heat capacity
  REAL(DP), ALLOCATABLE :: cvf_ptt(:,:)    ! the constant volume heat capacity
!
!  The anharmonic quantities
!
  REAL(DP), ALLOCATABLE :: betaf_ptt(:,:)  ! the volume thermal expansion 
  REAL(DP), ALLOCATABLE :: cpf_ptt(:,:)    ! the constant pressure heat 
                                          ! capacity
  REAL(DP), ALLOCATABLE :: b0f_s_ptt(:,:)  ! the isoentropic bulk modulus
  REAL(DP), ALLOCATABLE :: gammaf_ptt(:,:) ! the average Gruneisen parameter
!
! Anharmonic quantities for anisotropic thermodynamics
!
  REAL(DP), ALLOCATABLE :: celldmf_ptt(:,:,:)      ! the crystal parameters 
  REAL(DP), ALLOCATABLE :: celldmf_ptt_p1(:,:,:)   ! the crystal parameters at
                                                  ! T+deltat
  REAL(DP), ALLOCATABLE :: celldmf_ptt_m1(:,:,:)   ! the crystal parameters at
                                                  ! T-deltat
  REAL(DP), ALLOCATABLE :: alphaf_anis_ptt(:,:,:) ! the thermal expansion 
                                                  ! tensor
  REAL(DP), ALLOCATABLE :: densityf_ptt(:,:)       ! the density
  REAL(DP), ALLOCATABLE :: cpmcef_anis_ptt(:,:)    ! the difference cp-ce
  REAL(DP), ALLOCATABLE :: bthsf_ptt(:,:,:,:)      ! the thermal stress
  REAL(DP), ALLOCATABLE :: ggammaf_ptt(:,:,:,:)    ! the generalized Gruneisen
                                                  ! parameters
  REAL(DP), ALLOCATABLE :: csmctf_ptt(:,:,:,:) ! isoentropic elastic constants
                                             ! - isothermal ones
  REAL(DP), ALLOCATABLE :: el_consf_ptt(:,:,:,:)    ! isothermal elast. cons.
  REAL(DP), ALLOCATABLE :: el_compf_ptt(:,:,:,:)    ! isothermal elast. comp.
  REAL(DP), ALLOCATABLE :: el_consf_s_ptt(:,:,:,:)  ! isoentropic elast. cons.
  REAL(DP), ALLOCATABLE :: el_compf_s_ptt(:,:,:,:)  ! isoentropic elast. comp.
  REAL(DP), ALLOCATABLE :: macro_elf_ptt(:,:,:)     ! macroscopic elasticity
  REAL(DP), ALLOCATABLE :: macro_elf_s_ptt(:,:,:)   ! isoentropic macros. ela.
  REAL(DP), ALLOCATABLE :: vf_ptt(:,:,:)            ! isothermal sound speed
  REAL(DP), ALLOCATABLE :: vf_s_ptt(:,:,:)          ! isoentropic sound speed
  REAL(DP), ALLOCATABLE :: debye_macro_elf_ptt(:,:)  ! isothermal debye temper.
  REAL(DP), ALLOCATABLE :: debye_macro_elf_s_ptt(:,:)  ! adiabatic debye temper.

END MODULE ph_freq_anharmonic_ptt
!
!----------------------------------------------------------------------------
MODULE anharmonic_vt
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the Helmholtz free energy fitted by an equation 
!   of state plus a polynomial. These variables are computed at all 
!   temperatures for selected volumes.
!
  USE kinds, ONLY: DP
  SAVE
  REAL(DP), ALLOCATABLE :: press_vt(:,:) ! the thermal pressure as a function
                           ! of volume for selected temperatures
  REAL(DP), ALLOCATABLE :: press_vtt(:,:) ! the thermal pressure as a function
                           ! of temperature for selected volumes
END MODULE anharmonic_vt
!----------------------------------------------------------------------------
MODULE ph_freq_anharmonic_vt
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the Helmholtz free energy fitted by an equation 
!   of state plus a polynomial. These variables are computed at all 
!   temperatures for selected volumes.
!
  USE kinds, ONLY: DP
  SAVE
  REAL(DP), ALLOCATABLE :: pressf_vt(:,:) ! the thermal pressure as a function
                           ! of volume for selected temperatures
  REAL(DP), ALLOCATABLE :: pressf_vtt(:,:) ! the thermal pressure as a function
                           ! of temperature for selected volumes
END MODULE ph_freq_anharmonic_vt

!----------------------------------------------------------------------------
MODULE ph_freq_anharmonic
!----------------------------------------------------------------------------
!
!   The variables needed to describe the anharmonic quantities calculated
!   from the minimum of the free energy fitted by an equation of state.
!   Quantities obtained from the free energy calculated from all 
!   the frequencies.
!
  USE kinds, ONLY: DP
  USE polynomial,  ONLY : poly1, poly2, poly3, poly4
  SAVE
!
!  The parameters of the interpolation at each temperature
!
  REAL(DP), ALLOCATABLE :: vminf_t(:)  ! the minimum volume 
  REAL(DP), ALLOCATABLE :: b0f_t(:)    ! the bulk modulus B from EOS
  REAL(DP), ALLOCATABLE :: b0f_ec_t(:) ! the bulk modulus B from EC
  REAL(DP), ALLOCATABLE :: b01f_t(:)   ! pressure derivative of B
  REAL(DP), ALLOCATABLE :: b02f_t(:)   ! second pressure derivative of B
  REAL(DP), ALLOCATABLE :: free_e_minf_t(:) ! total free energy at the minimum
  REAL(DP), ALLOCATABLE :: af_t(:,:) ! the coefficients of the polynomial
                           ! that interpolates the free energy 
  TYPE(poly1), ALLOCATABLE :: p1tf_t(:)    ! coefficients of the polynomial
  TYPE(poly2), ALLOCATABLE :: p2tf_t(:)    ! that interpolate the free energy
  TYPE(poly3), ALLOCATABLE :: p3tf_t(:)    ! with multidimensional integral.
  TYPE(poly4), ALLOCATABLE :: p4tf_t(:)
!
!  The interpolated harmonic quantities
!
  REAL(DP), ALLOCATABLE :: enerf_t(:)   ! vibrational energy (T)
  REAL(DP), ALLOCATABLE :: free_enerf_t(:) ! vibrational free energy (T)
  REAL(DP), ALLOCATABLE :: entropyf_t(:)! entropy (T)
  REAL(DP), ALLOCATABLE :: cef_t(:)     ! constant strain heat capacity (T)
  REAL(DP), ALLOCATABLE :: cvf_t(:)     ! isocoric heat capacity (T)
  REAL(DP) :: e0f_t                     ! zero point energy 
!
!  The calculated anharmonic quantities
!
  REAL(DP), ALLOCATABLE :: cpf_t(:)     ! isobaric heat capacity (T)
  REAL(DP), ALLOCATABLE :: b0f_s(:)     ! constant entropy bulk modulus from eos
  REAL(DP), ALLOCATABLE :: b0f_ec_s(:)  ! constant entropy bulk modulus from ec
  REAL(DP), ALLOCATABLE :: alphaf_t(:)  ! linear thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: betaf_t(:)   ! volume thermal expansion coefficient
  REAL(DP), ALLOCATABLE :: gammaf_t(:)  ! average gruneisen parameter

!
!  Anharmonic quantities for anisotropic thermodynamic
!
  REAL(DP), ALLOCATABLE :: celldmf_t(:,:) ! the lattice parameters as a 
                           ! function of temperature
  REAL(DP), ALLOCATABLE :: densityf_t(:) ! the density
  REAL(DP), ALLOCATABLE :: alphaf_anis_t(:,:)  ! thermal expansion tensor

  REAL(DP), ALLOCATABLE :: cpmcef_anis(:) ! difference cp-ce computed from
                                          ! elastic constants
  REAL(DP), ALLOCATABLE :: bfactf_t(:,:,:)! b factor as a function of 
                                          ! temperature
  REAL(DP), ALLOCATABLE :: csmctf_t(:,:,:) ! difference of elastic constants
  REAL(DP), ALLOCATABLE :: bthsf_t(:,:,:)   ! thermal stress
  REAL(DP), ALLOCATABLE :: ggammaf_t(:,:,:) ! generalized average gruneisen 
                                          ! parameter

!
!  The parameters of the interpolation neglecting the electronic exitation
!  contribution
!
  REAL(DP), ALLOCATABLE :: vminf_noe_t(:)  ! minimum volume at each T 
                                          ! without el.
  REAL(DP), ALLOCATABLE :: celldmf_noe_t(:,:)  ! the lattice parameters as a 
                                          ! function of temperature
  REAL(DP), ALLOCATABLE :: densityf_noe_t(:)  ! the density
  REAL(DP), ALLOCATABLE :: b0f_noe_t(:)    ! bulk modulus at each T 
                                          ! without elec.
  REAL(DP), ALLOCATABLE :: b01f_noe_t(:)   ! pressure derivative of b0 at each T
                                          ! without electrons
  REAL(DP), ALLOCATABLE :: b02f_noe_t(:)   ! second pressure derivative of b0
                                          ! at each T without electrons
  REAL(DP), ALLOCATABLE :: free_e_minf_noe_t(:) ! total free energy at 
                                       ! the minimum at each T without el.
  REAL(DP), ALLOCATABLE :: enerf_noe_t(:)   ! vibrational energy (T)
  REAL(DP), ALLOCATABLE :: free_enerf_noe_t(:) ! vibrational free energy (T)
  REAL(DP), ALLOCATABLE :: entropyf_noe_t(:)! entropy (T)
  REAL(DP) :: e0f_noe_t                     ! zero point energy 

  REAL(DP), ALLOCATABLE :: af_noe_t(:,:) ! the coefficients of the polynomial
                           ! that interpolates the free energy at each
                           ! temperature
  TYPE(poly1), ALLOCATABLE :: p1tf_noe_t(:) ! coefficients of the polynomial
  TYPE(poly2), ALLOCATABLE :: p2tf_noe_t(:) ! that interpolate the free energy
  TYPE(poly3), ALLOCATABLE :: p3tf_noe_t(:) ! with multidimensional integral.
  TYPE(poly4), ALLOCATABLE :: p4tf_noe_t(:)
!
!  The interpolated harmonic quantities neglecting the electronic effect
!
  REAL(DP), ALLOCATABLE :: cef_noe_t(:)     ! constant strain heat capacity (T)
  REAL(DP), ALLOCATABLE :: cvf_noe_t(:)     ! isocoric heat capacity (T)
!
!  The anharmonic quantities neglecting the electronic effect
!
  REAL(DP), ALLOCATABLE :: cpf_noe_t(:)     ! isobaric heat capacity (T)
  REAL(DP), ALLOCATABLE :: b0f_noe_s(:)     ! constant entropy bulk modulus
  REAL(DP), ALLOCATABLE :: betaf_noe_t(:)   ! volume thermal expansion
  REAL(DP), ALLOCATABLE :: gammaf_noe_t(:)  ! average gruneisen parameter
!
!  elastic constants and related quantities
!
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
  REAL(DP), ALLOCATABLE :: dydef_t(:,:) !(nstep,ntemp) the derivative of the
                                       ! internal parameter with respect to 
                                       ! strain
  REAL(DP), ALLOCATABLE :: debye_macro_elf_t(:) ! debye temperature from
                                           ! macro elasticity (isothermal)
  REAL(DP), ALLOCATABLE :: debye_macro_elf_s(:) ! debye temperature from
                                           ! macro elasticity (isoentropic)
  REAL(DP), ALLOCATABLE :: el_conf_geo_t(:,:,:,:) ! the temperature dependent
                                            ! elastic constants at all 
                                            ! geometries
  REAL(DP), ALLOCATABLE :: celldmf_t_p1(:,:)! the celldm at the pressure+dp
                                           ! as a function of T
  REAL(DP), ALLOCATABLE :: celldmf_t_m1(:,:)! the celldm at the pressure-dp
                                           ! as a function of T
  REAL(DP), ALLOCATABLE :: celldmf_noe_t_p1(:,:)! the celldm at the pressure+dp
                                           ! as a function of T
  REAL(DP), ALLOCATABLE :: celldmf_noe_t_m1(:,:)! the celldm at the pressure-dp
                                           ! as a function of T

END MODULE ph_freq_anharmonic

!----------------------------------------------------------------------------
MODULE grun_anharmonic
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE el_anharmonic
!----------------------------------------------------------------------------
!
!   here all the quantities obtained adding the free energy of the
!   electrons
!
  USE kinds, ONLY: DP
  SAVE
!
!   The parameters of the interpolation. Used when the vibrational free
!   energy is not computed
!
  REAL(DP), ALLOCATABLE :: vmine_t(:),  &   ! volume at the minimum
                           b0e_t(:),    &   ! bulk modulus (B) at the minimum
                           b01e_t(:),   &   ! pressure derivative of B
                           b02e_t(:),   &   ! second pressure derivative of B
                           free_e_mine_t(:) ! free energy at the minimum
!
!   Interpolated electronic quantities
!
  REAL(DP), ALLOCATABLE :: el_energy_t(:) ! electronic energy 
  REAL(DP), ALLOCATABLE :: el_free_energy_t(:) ! electronic 
  REAL(DP), ALLOCATABLE :: el_entropy_t(:) !  electronic entropy
  REAL(DP), ALLOCATABLE :: el_ce_t(:) ! electronic heat capacity

  REAL(DP), ALLOCATABLE :: el_energyf_t(:) ! electronic energy 
  REAL(DP), ALLOCATABLE :: el_free_energyf_t(:) ! electronic free energy 
  REAL(DP), ALLOCATABLE :: el_entropyf_t(:) ! electronic entropy
  REAL(DP), ALLOCATABLE :: el_cef_t(:) ! electronic heat capacity

  REAL(DP), ALLOCATABLE :: el_beta_t(:) ! the electronic contribution to the
                                        ! thermal expansion (with phdos)
  REAL(DP), ALLOCATABLE :: el_betaf_t(:) ! the electronic contribution to the
                                        ! thermal expansion
  REAL(DP), ALLOCATABLE :: el_cp_t(:) ! the electronic contribution to the
                                        ! isobaric heat capacity
  REAL(DP), ALLOCATABLE :: el_cpf_t(:) ! the electronic contribution to the
                                        ! isobaric heat capacity
  REAL(DP), ALLOCATABLE :: el_gamma_t(:) ! the electronic contribution to the
                                        ! average gruneisen parameter
  REAL(DP), ALLOCATABLE :: el_gammaf_t(:) ! the electronic contribution to the
                                        ! average gruneisen parameter
  REAL(DP), ALLOCATABLE :: el_b0_t(:)   ! the electronic contribution to the
                                        ! bulk modulus
  REAL(DP), ALLOCATABLE :: el_b0f_t(:)   ! the electronic contribution to the
                                        ! bulk modulus
!
!  Harmonic quantities interpolated at the temperature dependent geometry 
!  at selected pressures 
!
  REAL(DP), ALLOCATABLE :: el_free_ener_pt(:,:) ! energy
  REAL(DP), ALLOCATABLE :: el_ener_pt(:,:) ! energy
  REAL(DP), ALLOCATABLE :: el_entr_pt(:,:) ! entropy
  REAL(DP), ALLOCATABLE :: el_ce_pt(:,:) ! heat capacity
!
!  Harmonic quantities interpolated at the temperature dependent geometry 
!  at selected pressures 
!
  REAL(DP), ALLOCATABLE :: el_free_enerf_pt(:,:) ! energy
  REAL(DP), ALLOCATABLE :: el_enerf_pt(:,:) ! energy
  REAL(DP), ALLOCATABLE :: el_entrf_pt(:,:) ! entropy
  REAL(DP), ALLOCATABLE :: el_cef_pt(:,:) ! heat capacity
!
!  Harmonic quantities interpolated at the pressure dependent geometries
!  at selected temperatures
!
  REAL(DP), ALLOCATABLE :: el_ce_ptt(:,:) ! heat capacity

  REAL(DP), ALLOCATABLE :: el_cef_ptt(:,:) ! heat capacity

END MODULE el_anharmonic

!----------------------------------------------------------------------------
MODULE ifc
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE control_dosq
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE control_thermo
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the type of calculation to do
  !
  SAVE
  !
  LOGICAL, ALLOCATABLE :: lpwscf(:),  & ! if .true. this work requires a scf calc.
                          lpwband(:), & ! if .true. this work makes a band
                                        ! calculation
                          lef(:),     & ! if .true. save the fermi level 
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
                          comp_iq_iw(:,:), & ! for each work the list
                                         ! of comp_iq of that work.
                          done_irr_iq_iw(:,:,:), & ! for each work the 
                                         ! done_irr_iq of that work
                          done_iq_iw(:,:) ! for each work the list
                                         ! of done_iq of that work.
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
  LOGICAL :: ltherm_glob=.FALSE.  ! if true computes the anharmonic properties
                                  ! by fitting the free energy versus volume
                                  ! at each temperature with an eos, 
                                  ! otherwise fit the temperature dependent 
                                  ! part of the free energy with a polynomial
  LOGICAL :: lhugoniot=.FALSE.    ! if .true. the code plots T(p) and V(p)
                                  ! along the Hugoniot
  LOGICAL :: lgeo_from_file=.FALSE. ! geometries for mur_lc read from file
  LOGICAL :: lgeo_to_file=.FALSE. ! geometries for mur_lc written to file
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

!----------------------------------------------------------------------------
MODULE distribute_collection
!----------------------------------------------------------------------------

SAVE

LOGICAL, ALLOCATABLE :: me_igeom(:) ! if .TRUE. this image must collect
                                    ! some work in this geometry

LOGICAL, ALLOCATABLE :: me_igeom_iq(:,:) ! if .TRUE. this image must collect
                                    ! some q point in this geometry

END MODULE distribute_collection

!----------------------------------------------------------------------------
MODULE control_elastic_constants
!----------------------------------------------------------------------------
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
                                !
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

  LOGICAL :: lelastic=.FALSE.   ! elastic constants as a function of &
                                ! temperature available in some approximation
  LOGICAL :: lelasticf=.FALSE.  ! elastic constants as a function of pressure
                                ! available in some approximation
  LOGICAL :: lelastic_p=.FALSE. ! elastic constants as a function of pressure
                                ! available
  LOGICAL :: lelastic_pt=.FALSE. ! elastic constants as a function of 
                                ! temperature for a few pressures available
  LOGICAL :: lelastic_ptt=.FALSE. ! elastic constants as a function of 
                                ! pressure for a few temperatures available
                                !
  LOGICAL :: el_cons_available=.FALSE.  ! when this flag becomes true it
                                ! means that the elastic constant have been
                                ! read from file and are available

  LOGICAL :: el_cons_geo_available=.FALSE.  ! when this flag becomes true it
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
  LOGICAL :: lelasticf_pt=.FALSE. ! elastic constants as a function of 
                                ! temperature for a few pressures available
  LOGICAL :: lelasticf_ptt=.FALSE. ! elastic constants as a function of 
                                ! pressure for a few temperatures available
                                !

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

  LOGICAL, ALLOCATABLE :: all_geometry_done_geo(:) ! for each unperturbed
                                ! geometry it is true if all dynamical
                                ! matrices of that geometry are computed

  INTEGER :: fndos_ec           ! number of qha elastic constants found on file
                                ! free energy computed from phonon dos

  INTEGER :: fnph_ec            ! number of qha elastic constants found on file
                                ! free energy computed from direct integral

  INTEGER :: f_geodos_ec        ! the first geometry that has qha ec

  INTEGER :: f_geoph_ec         ! the first geometry that has qha ec

  LOGICAL, ALLOCATABLE :: found_dos_ec(:) ! for each geometry it is true if
                                ! qha elastic constant have been found
  LOGICAL, ALLOCATABLE :: found_ph_ec(:) ! for each geometry it is true if
                                ! qha elastic constant have been found

  INTEGER :: nstep_ec           ! number of strain_types
                                !
  INTEGER :: nmove              ! number of atomic positions to study
                                !
  LOGICAL :: stype(21)          ! if .TRUE. this strain type needs also
                                ! atomic relaxations
  INTEGER :: move_at(21)        ! The atoms that move in any strain type

  REAL(DP) :: atom_dir(3,21)    ! The versor of the atomic displacement
  REAL(DP) :: atom_step(21)     ! The amount of displacement for each 
                                ! strain type

  REAL(DP), ALLOCATABLE :: tau_acc(:,:,:) ! a possible displacement of
                                ! the atoms with respect to the uniformely
                                ! strained configuration
                                !
  REAL(DP), ALLOCATABLE :: min_y(:,:,:) ! the minimum value of the internal
                                ! coordinate (ngeo_strain,21,ngeom)
  REAL(DP), ALLOCATABLE :: epsil_y(:,:,:) ! the strain amplitude
                                ! that has min_y internal coordinate 
                                ! as a minimum
  REAL(DP), ALLOCATABLE :: min_y_t(:,:,:,:) ! the minimum value of the internal
                                ! coordinate (ngeo_strain,21,ngeom,ntemp) at 
                                ! each temperature
  LOGICAL :: lcm_ec             ! if .true. the code moves the other atoms
                                ! so as to keep the center of mass of the
                                ! cell fixed.
  LOGICAL :: lzsisa             ! when .TRUE. and the previous calculations
                                ! have been executed makes the zsiza 
                                ! approximation. Only min_y is used at
                                ! each temperature
                                ! Default : .FALSE.
  LOGICAL :: lfp                ! when .TRUE. and the previous calculations
                                ! have been executed makes the frozen phonon
                                ! approximation, taking the min_y=0.0

  REAL(DP), ALLOCATABLE :: dyde(:,:,:) ! (21,ngeom,ntemp)

END MODULE control_elastic_constants
  !
!----------------------------------------------------------------------------
MODULE control_conv
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the control of convergence
  !
  SAVE

  REAL(DP), ALLOCATABLE :: ke(:)    ! the kinetic energies that are tested
  REAL(DP), ALLOCATABLE :: keden(:) ! the kinetic energies of the charge
  REAL(DP), ALLOCATABLE :: sigma_test(:) ! the smearing values
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

!----------------------------------------------------------------------------
MODULE control_paths
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE control_bands
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE control_2d_bands
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE proj_rap_point_group
!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------
MODULE control_eldos
!----------------------------------------------------------------------------

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

  LOGICAL :: lel_free_energy ! if .TRUE. add the electron free energy to
                            ! the total energy
  LOGICAL :: hot_electrons   ! if .TRUE. the code expects to find the 
                             ! energy for several values of sigma in 
                             ! directories called restart2, restart3 ...
                             ! and the electronic contribution is not
                             ! evaluated from eldos.
END MODULE control_eldos

!----------------------------------------------------------------------------
MODULE control_grun
!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------
MODULE control_ev
!----------------------------------------------------------------------------
  USE kinds, ONLY: DP
  SAVE

INTEGER :: ieos     ! the equation of state to use
INTEGER :: npt      ! the number of points for the fit

REAL(DP), ALLOCATABLE :: v0(:)   ! the volume of each point
REAL(DP), ALLOCATABLE :: e0(:)   ! the energy of each point

END MODULE control_ev

!----------------------------------------------------------------------------
MODULE initial_conf
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE initial_param
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to save the initial parameters
  !
  SAVE

  REAL(DP) :: ecutwfc0  ! initial cutoff on wavefunction
  REAL(DP) :: ecutrho0  ! initial cutoff on the charge density

  REAL(DP) :: ethr0     ! the initial accuracy of the eigenvalues
  

END MODULE initial_param

!----------------------------------------------------------------------------
MODULE equilibrium_conf
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE control_pwrun
!----------------------------------------------------------------------------

  USE kinds, ONLY: DP
  SAVE

  LOGICAL  :: do_punch=.TRUE.  ! set this variable to .FALSE. if pw has
                               ! not to save the punch files.
END MODULE control_pwrun

!----------------------------------------------------------------------------
MODULE control_phrun
!----------------------------------------------------------------------------

  SAVE

  CHARACTER(LEN=256)  :: auxdyn=""  ! the name of the dynamical matrix file

END MODULE control_phrun

!----------------------------------------------------------------------------
MODULE control_energy_plot
!----------------------------------------------------------------------------

  USE kinds, ONLY: DP
  SAVE

  INTEGER :: ncontours
  REAL(DP), ALLOCATABLE :: ene_levels(:)
  CHARACTER(LEN=12), ALLOCATABLE :: color_levels(:)

END MODULE control_energy_plot

!----------------------------------------------------------------------------
MODULE control_debye
!----------------------------------------------------------------------------

  USE kinds, ONLY: DP
  SAVE

  REAL(DP) :: debye_t, deb_e0

  REAL(DP), ALLOCATABLE :: deb_energy(:),       &  ! Debye vibrational energy
                           deb_free_energy(:),  &  ! Debye free energy
                           deb_entropy(:),      &  ! Debye entropy
                           deb_cv(:),           &  ! Debye cv
                           deb_b_fact(:,:,:,:), &  ! axiliary B factor
                           deb_bfact(:)            ! Debye atomic B factors

  INTEGER :: idebye                            ! find the debye temperature
                                               ! 1 free energy
                                               ! 2 energy
                                               ! 3 isochoric heat capacity
                                               ! 0 none

END MODULE control_debye

!----------------------------------------------------------------------------
MODULE control_macro_elasticity
!----------------------------------------------------------------------------

  USE kinds, ONLY: DP
  SAVE

  REAL(DP) :: macro_el(8)     ! the Voigt and Reuss approximations
  REAL(DP) :: vp, vb, vg      ! the sound velocities
  REAL(DP) :: approx_debye_t  ! approximate Debye temperature

END MODULE control_macro_elasticity

!----------------------------------------------------------------------------
MODULE control_quadratic_energy
!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------
MODULE control_quartic_energy
!----------------------------------------------------------------------------

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
  REAL(DP), ALLOCATABLE ::    &
             hessian4_v(:,:), &              ! hessian eigenvectors
             hessian4_e(:),   &              ! hessian eigenvalues
             x_min_4(:),      &              ! coordinates of the minimum
             coeff4(:)                       ! coefficients of quartic fit
  INTEGER :: lsolve                          ! 1, 2, 3 controls the method
                                             ! used to find the polynomial
                                             ! coefficients (Default 2)
  TYPE(poly4) :: p4                          ! coefficients of the polynomial

END MODULE control_quartic_energy


!----------------------------------------------------------------------------
MODULE control_xrdp
!----------------------------------------------------------------------------

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

!----------------------------------------------------------------------------
MODULE control_gnuplot
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE data_files
!----------------------------------------------------------------------------
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
  CHARACTER(LEN=256) :: flelanhar ! file with the anharmonic quantities
                                ! obtained adding only the electronic 
                                ! free energy
  CHARACTER(LEN=256) :: flevdat ! file with data for ev.x 
  CHARACTER(LEN=256) :: flenergy  ! the name of the file with the energy
                                  ! suited for gnuplot contour plots
  CHARACTER(LEN=256) :: flepsilon ! the name of the file with the dielectric
                                  ! constant
  CHARACTER(LEN=256) :: floptical ! the name of the file with the optical
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
  CHARACTER(LEN=256) :: flgeom   ! the file with the saved geometries 

END MODULE data_files

!----------------------------------------------------------------------------
MODULE postscript_files
!----------------------------------------------------------------------------
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
  CHARACTER(LEN=256) :: flpsoptical ! the name of the file with the optical
                                  ! constant
  CHARACTER(LEN=256) :: flpsepsilon ! the name of the file with the dielectric
                                  ! constant
  CHARACTER(LEN=256) :: flps_el_cons ! the name of the file with the 
                                  ! elastic constants
END MODULE postscript_files

!----------------------------------------------------------------------------
MODULE internal_files_names
!----------------------------------------------------------------------------
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

  CHARACTER(LEN=256) :: fleldos_thermo    ! the name of the file with electron
                                          ! dos
  CHARACTER(LEN=256) :: flpseldos_thermo  ! the name of the postscript file 
                                          ! with electron dos
  CHARACTER(LEN=256) :: fleltherm_thermo  ! the name of the file with the 
                                          ! electronic thermodynamic properties
  CHARACTER(LEN=256) :: flpseltherm_thermo  ! the name of the postscript
                                          ! file with the electronic 
                                          !thermodynamic properties

END MODULE internal_files_names

!----------------------------------------------------------------------------
MODULE control_asy
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to control the brillouin zone plot with asymptote
  !
  SAVE
  !
  CHARACTER(LEN=256) :: flasy ! the name of file with the asymptote script

  CHARACTER(LEN=256) :: asymptote_command ! the asymptote command

  LOGICAL :: lasymptote ! if .true. asymptote is called within the thermo_pw
                        ! code

END MODULE control_asy

!----------------------------------------------------------------------------
MODULE thermo_sym
!----------------------------------------------------------------------------
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

!----------------------------------------------------------------------------
MODULE efermi_plot
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP

  SAVE

  INTEGER :: n1, n2     ! the number of k points of the 2d plot

  REAL(DP) :: kxmin, kxmax, kymin, kymax ! the size of the plotted region

  LOGICAL, ALLOCATABLE :: has_ef(:) ! say which bands cross ef
END MODULE efermi_plot

!----------------------------------------------------------------------------
MODULE control_emp_free_ener
!----------------------------------------------------------------------------
  USE kinds,  ONLY : DP

  SAVE

  LOGICAL :: add_empirical ! if .TRUE. adds the empirical free energy term

  INTEGER :: efe        ! type of empirical free energy
                        ! 1) (alpha1+alpha2 * V) T^2
                        ! 2) -2/3 k_B nat alpha1 (v/v0p)^alpha2 T^2

  REAL(DP) :: alpha1, &  ! parameter of the empirical free energy
                         ! eV/K (efe=1) 1/K (efe=2)
              alpha2, &  ! parameter in eV/A^3/K (efe=1), adimensional (efe=2) 
              v0p        ! equilibrium volume in (a.u.)^3

  REAL(DP), ALLOCATABLE :: emp_ener(:,:)     ! empirical total energy, 
                                             ! T, geometry
  REAL(DP), ALLOCATABLE :: emp_free_ener(:,:) ! empirical  free_energy, T, 
                                             ! geometry
  REAL(DP), ALLOCATABLE :: emp_entr(:,:)     ! empirical entropy, T, geometry
  REAL(DP), ALLOCATABLE :: emp_ce(:,:)       ! empirical heat capacity, 
                                             ! T, geometry

END MODULE control_emp_free_ener

