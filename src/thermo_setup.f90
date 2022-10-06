!
! Copyright (C) 2016-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE thermo_setup()
  !-----------------------------------------------------------------------
  !
  !   This routine is called only once, after thermo_pw has read its
  !   input and the input of pw.x (and ph.x). It is used to save the
  !   initial configuration of the run and to set a few variables
  !   that need to be adjusted at the beginning of the calculation.
  !   It sets also the temperature and the pressure.
  !   Variables adjusted or calculated:
  !
  !   deltae      : for electronic dos calculation if not set
  !   ngeo        : if not set, or cleared to be compatible with the structure
  !   ngeo_strain : if too small or not set
  !   poly_degree : if too small or not set
  !   temp_ph     : if not set
  !   lambda      : if not set
  !   sur_layers  : if not set
  !   the path    : if not given in input is set here
  !
  !   celldm0     : set equal to celldm
  !   It changes the variables of QE only if ibrav=0 and allowed by the 
  !   user. In this case it changes ibrav, celldm, tau, at, and bg.
  !
  !   It sets the variable with_asyn_images if nimage>1
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : ry_kbar
  USE constants,            ONLY : k_boltzmann_ry
  USE thermo_mod,           ONLY : what, ngeo, fact_ngeo, ngeo_ph    
  USE temperature,          ONLY : tmin, ntemp_plot
  USE control_thermo,       ONLY : continue_zero_ibrav, find_ibrav, &
                                   set_internal_path, set_2d_path
  USE control_elastic_constants, ONLY : ngeo_strain, elastic_algorithm, &
                                  elcpvar, poly_degree, elalgen
  USE control_eldos,        ONLY : deltae, ndose, lel_free_energy
  USE control_mur,          ONLY : lmurn
  USE control_pressure,     ONLY : pmin, pmax, deltap, npress_plot
  USE equilibrium_conf,     ONLY : tau0, tau0_crys
  USE control_paths,        ONLY : npk_label, lbar_label
  USE control_grun,         ONLY : temp_ph, volume_ph, celldm_ph
  USE control_xrdp,         ONLY : lambda, lambda_elem
  USE control_2d_bands,     ONLY : lprojpbs, nkz, sur_layers, identify_sur, &
                                   sp_min
!
!  variables modified by this routine
!
  USE initial_conf,         ONLY : celldm_save, ibrav_save, ityp_save,       &
                                   nr1_save, nr2_save, nr3_save, &
                                   nosym_save, tau_save, tau_save_crys, &
                                   omega_save, at_save
  USE initial_param,        ONLY : ecutwfc0, ecutrho0, ethr0
  USE equilibrium_conf,     ONLY : nr1_0, nr2_0, nr3_0
  USE thermo_sym,           ONLY : code_group_save
!
!   helper codes from the library
!
  USE lattices,             ONLY : find_ibrav_code
  USE xrdp_module,          ONLY : select_lambda
  USE color_mod,            ONLY : set_colors
!
!  these are of QE that can be changed here:
!  ibrav, celldm, tau, dfftp parameters.
!  This routine sets also the point group variables of QE.
!
  USE cell_base,            ONLY : at, bg, ibrav, celldm, omega, &
                                   cell_base_init, alat
  USE ions_base,            ONLY : nat, tau, ntyp => nsp, ityp, amass, atm, &
                                   if_pos

  USE gvecw,                ONLY : ecutwfc
  USE gvect,                ONLY : ecutrho
  USE symm_base,            ONLY : nosym
  USE rap_point_group,      ONLY : code_group
  USE fft_base,             ONLY : dfftp
  USE klist,                ONLY : degauss, ltetra
  USE control_flags,        ONLY : ethr
  USE io_global,            ONLY : stdout
  USE mp_images,            ONLY : nimage
  USE mp_asyn,              ONLY : with_asyn_images

  !
  IMPLICIT NONE
  !
  INTEGER :: ia, ipol, jpol, igeo, ibrav_, code_group_ext
  REAL(DP), PARAMETER :: eps1=1D-8
  REAL(DP) :: ur(3,3), global_s(3,3), rd_ht(3,3), celldm_(6), alat_save, zero
  REAL(DP), ALLOCATABLE :: tau_aux(:,:)
  LOGICAL :: cubic
  !
  !
  with_asyn_images=(nimage>1)
  !
  !  
  cubic=(ibrav==1.OR.ibrav==2.OR.ibrav==3) 
  IF ((npress_plot>0.OR.ntemp_plot>0).AND.(.NOT.lmurn)) &
    CALL errore('thermo_setup','npress_plot and ntemp_plot need lmurn=.TRUE.',1) 
  !
  CALL set_temperature()
  !
  CALL set_pressure()
  !
  CALL set_pressure_kb()
  !
  !   ngeo_strain cannot be too small and the energy algorithm requires a
  !   few more points
  !
  elalgen=elastic_algorithm=='energy_std'.OR.elastic_algorithm=='energy'
  IF (ngeo_strain<4) THEN
     ngeo_strain=4
     IF (elalgen) ngeo_strain=6
  ENDIF
  !
  !   The default of the interpolation polynomial for elastic constants, if
  !   not set in input, or if the input value is unreasonable
  !
  IF (poly_degree < 2 ) THEN
     poly_degree = 3
     IF (elalgen) poly_degree=4
     IF (ngeo_strain < 6) THEN
       poly_degree = 2
       IF (elalgen) poly_degree=3
     ENDIF
  ENDIF
  elcpvar=poly_degree+1
  IF (ngeo_strain < elcpvar) CALL errore('thermo_setup','ngeo_strain is &
                                                                 &too small',1)
!
!  set a few default values
!
  IF (volume_ph==0.0_DP.AND.celldm_ph(1)==0.0_DP.AND.temp_ph==0.0_DP) &
                                                              temp_ph=tmin

  IF (lambda==0.0_DP) CALL select_lambda(lambda_elem,lambda)
!
!  the pressure must be in Ry/(a.u.)^3 but in input it was given in kbar
!
    pmax=pmax/ry_kbar
    pmin=pmin/ry_kbar
    deltap=deltap/ry_kbar
!
!   here deal with ibrav=0. The code finds ibrav and celldm and
!   writes them on output or change them and continue the calculation
!
  IF (ibrav==0.AND..NOT.continue_zero_ibrav) THEN
     ALLOCATE(tau_aux(3,nat))
     alat_save=celldm(1)
     at=at*alat_save
     code_group_ext=0
     CALL find_ibrav_code(at(1,1),at(1,2),at(1,3),ibrav_,celldm_, &
                                 code_group_ext, ur, global_s, .FALSE.)
     zero=0.0_DP
     rd_ht=0.0_DP

     CALL cell_base_init(ibrav_, celldm_, zero, zero, zero, zero, zero, zero, &
                      .FALSE., rd_ht, ' ' )

     tau_aux=tau*alat_save
     tau=0.0_DP
     DO ia=1,nat
        DO ipol=1,3
           DO jpol=1,3
              tau(ipol,ia)=tau(ipol,ia) + global_s(jpol,ipol)*tau_aux(jpol,ia)
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
!
!   by setting what to an empty string all the rest of thermo_pw is not
!   run and the code exits
!
        what=' '
     ENDIF
     DEALLOCATE(tau_aux)
  ENDIF

  IF (what=='mur_lc_t'.AND.(ibrav==5.OR.ibrav==91.OR.ABS(ibrav)==12&
                           .OR.ABS(ibrav)==13.OR.ibrav==14.OR.ibrav==0)) &
     CALL errore('thermo_setup','Thermal expansion not available',1)

  IF (set_internal_path.AND.(ABS(ibrav)==13.OR.ibrav==14.OR.ibrav==0)) &
     CALL errore('thermo_setup','internal path not available',1)

  IF (npk_label>0.AND.ibrav==0) &
     CALL errore('thermo_setup','ibrav=0 and path labels not available',1)
 
  IF (set_internal_path) CALL set_bz_path()
  IF (set_2d_path) CALL set_2d_bz_path()
!
! Save the initial configurantion in the initial variables
!
  ALLOCATE(ityp_save(nat))
  ALLOCATE(tau_save(3,nat))
  ALLOCATE(tau_save_crys(3,nat))

  ibrav_save=ibrav
  at_save = at
  omega_save=omega
  celldm_save=celldm
  nosym_save=nosym
  ityp_save(:)=ityp(:)
  tau_save=tau
!
!  bring tau_save in crystal coordinates. In strained geometries tau_save
!  is kept constant.
!
  tau_save_crys=tau_save
  CALL cryst_to_cart( nat, tau_save_crys, bg, -1 )
!
!  We now check the point group and the space group so that we can
!  find the optimal fft mesh for the run. Save also the code group
!  of the unpertubed configuration
!
  CALL find_fft_fact()

  nr1_save=dfftp%nr1
  nr2_save=dfftp%nr2
  nr3_save=dfftp%nr3

  code_group_save=code_group
  ecutwfc0=ecutwfc
  ecutrho0=ecutrho
  ethr0=ethr
!
! The equilibrium configuration is set here for the case in which we
! do not minimize the energy. In this case we keep the input geometry
!
  ALLOCATE(tau0(3,nat))
  ALLOCATE(tau0_crys(3,nat))

  CALL set_equilibrium_conf(celldm, tau, at, omega)
  nr1_0=nr1_save
  nr2_0=nr2_save
  nr3_0=nr3_save
!
!   Some initialization on ngeo
!
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
           IF (cubic) ngeo(1)=9
        ENDIF
     ENDIF
  END IF
  CALL clean_ngeo(ngeo,fact_ngeo,ngeo_ph,ibrav)
!
!  Initialize colors
!
  CALL set_colors()
!
!  setup for the electronic dos
!
  IF (what(1:7)=='scf_dos'.OR.what(1:10)=='mur_lc_dos'.OR.lel_free_energy) &
                                                                         THEN
     IF (deltae==0.0_DP.AND.ndose==0.AND.(degauss>0.0_DP.OR.ltetra)) THEN
        deltae=1.5_DP * k_boltzmann_ry*MIN(4.0_DP, tmin)
        WRITE(stdout,'(/,5x,"Deltae set to",f20.9," Ry")') deltae
     ENDIF
  END IF
!
!  surface specific initializations
!
  IF (what /= 'scf_2d_bands') THEN
     nkz=1
     lprojpbs=.FALSE.
  ELSE
     lbar_label=.TRUE.
  ENDIF
!
!  default value of sp_min if not given in input
!
  IF (sp_min==0.0_DP) THEN
     sp_min=2.0_DP/alat
  ELSE
     sp_min=sp_min/alat
  ENDIF

  IF (what=='scf_2d_bands'.AND.identify_sur) THEN
     IF (sur_layers==0) THEN
        sur_layers=MIN(2, nat/2)
     ENDIF
  ENDIF
!
!  Initialize here the internal names of the files for dealing with many
!  geometries
!
   CALL initialize_file_names()
   CALL initialize_el_file_names()

  RETURN
  !
  END SUBROUTINE thermo_setup
