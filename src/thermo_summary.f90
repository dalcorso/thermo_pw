!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE thermo_summary()
  !-----------------------------------------------------------------------
  !
  !  This routine writes a summary of the input and calculates a few
  !  quantities that can be deduced from the pw input and that can be
  !  useful for the thermo_pw calculation
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi
  USE thermo_mod,           ONLY : what, ngeo, density
  USE thermo_sym,           ONLY : laue, code_group_save, fft_fact, &
                                   ibrav_group_consistent
  USE cell_base,            ONLY : ibrav
  USE control_paths,        ONLY : xqaux, wqaux, npk_label, letter, &
                                   label_list, nqaux, point_label_type, &
                                   letter_path
  USE control_asy,          ONLY : flasy, lasymptote, asymptote_command
  USE control_elastic_constants, ONLY : frozen_ions, ngeo_strain, &
                                  elastic_algorithm
  USE control_xrdp,         ONLY : lambda, flxrdp, flpsxrdp, lformf, &
                                   flformf, smin, smax, nspoint, lcm
  USE control_2d_bands,     ONLY : lprojpbs
  USE space_groups,         ONLY : sg_name, find_space_group, set_fft_fact, &
                                   project_frac_tran
  USE point_group,          ONLY : sym_label, find_group_info_ext
  USE control_pwrun,        ONLY : nr1_save, nr2_save, nr3_save
  USE ktetra,               ONLY : tetra, ltetra
  USE control_flags,        ONLY : iverbosity
  USE noncollin_module,     ONLY : noncolin, m_loc
  USE fft_base,             ONLY : dfftp, dffts
  USE spin_orb,             ONLY : domag
  USE rap_point_group,      ONLY : code_group
  USE cell_base,            ONLY : at, bg, celldm, omega
  USE ions_base,            ONLY : tau, nat, ityp, nsp, amass, atm
  USE symm_base,            ONLY : irt, nsym, s, sr, ftau
  USE constants,            ONLY : amu_si, bohr_radius_si
  USE bz_2d_form,           ONLY : bz_2d, allocate_2d_bz, init_2d_bz
  USE xrdp_module,          ONLY : write_form_factor, compute_xrdp
  USE mp_world,             ONLY : world_comm
  USE mp_images,            ONLY : nimage, my_image_id, root_image
  USE environment,          ONLY : environment_end
  USE mp_global,            ONLY : mp_global_end
  USE io_global,            ONLY : ionode, stdout
  USE io_files,             ONLY : prefix
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  INTEGER :: ios
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: asy_filename, filename, xsf_filename
  CHARACTER(LEN=11) :: group_name
  REAL(DP) :: total_mass, total_expected_mass, current_mass, expected_mass, fact
  REAL(DP) :: atom_weight, celldm_2d(3), ft(3,48), ftpar(3,48), ftperp(3,48)
  REAL(DP) :: s01(3), s02(3), s01c(3), s02c(3), s01mod, s02mod
  REAL(DP), ALLOCATABLE :: xau(:,:)
  INTEGER :: atomic_number
  INTEGER :: laue_class
  INTEGER :: it, ia, na, ipol, jpol, unique, trig, aux_sg
  LOGICAL :: read_path, lelc, lpiezo, ltherm_expansion
  INTEGER :: ierr, iuout, system, sg_number
  LOGICAL :: check_group_ibrav
  CHARACTER(LEN=12) :: spaceg_name
  CHARACTER(LEN=11) :: gname
  TYPE(bz_2d) :: bz_2d_struc
  INTEGER :: ibz, i, group_desc(48), which_elem(48), isym, jsym, &
             code_group_ext, code_group1

  read_path=.FALSE.
  lelc = .FALSE.
  ltherm_expansion = .FALSE.
  lpiezo=.FALSE.
  WRITE(stdout,'(/)')
  SELECT CASE (TRIM(what))
     CASE ('plot_bz') 
          WRITE(stdout,'(5x,"Plotting the Brillouin Zone and k points path")')
          read_path=.TRUE.
     CASE ('scf') 
          WRITE(stdout,'(5x,"Doing a single scf calculation")')
     CASE ('scf_bands') 
          WRITE(stdout,'(5x,"Doing a band calculation")')
          WRITE(stdout,'(5x,"Use what=plot_bz to visualize the BZ path")')
          read_path=.TRUE.
     CASE ('scf_dos')
          WRITE(stdout,'(5x,"Doing a electronic bands dos calculation")')
     CASE ('scf_ph') 
          WRITE(stdout,'(5x,"Doing a phonon calculation")')
     CASE ('scf_disp')
          WRITE(stdout,'(5x,"Doing a phonon dispersion calculation")')
          WRITE(stdout,'(5x,"Use what=plot_bz to visualize the BZ path")')
          WRITE(stdout,'(5x,"Computing the harmonic thermodynamic quantities")')
          read_path=.TRUE.
     CASE ('mur_lc') 
          WRITE(stdout,'(5x,"Calculating the volume that minimizes the energy")')
     CASE ('mur_lc_bands') 
          WRITE(stdout,'(5x,"Calculating the bands at the Murnaghan minimum &
                                                 &volume")')
          WRITE(stdout,'(5x,"Use what=plot_bz to visualize the BZ path")')
          read_path=.TRUE.
     CASE ('mur_lc_dos')
          WRITE(stdout,'(5x,"Calculating the electronic bands dos at minimum &
                                                 &volume")')
     CASE ('mur_lc_ph') 
          WRITE(stdout,'(5x,"Doing a phonon calculation at the Murnaghan &
                                         &minimum volume")')
     CASE ('mur_lc_disp') 
          WRITE(stdout,'(5x,"Doing a phonon dispersion calculation at the &
                                           & minimum volume")')
          WRITE(stdout,'(5x,"Use what=plot_bz to visualize the BZ path")')
          WRITE(stdout,'(5x,"Computing the harmonic thermodynamic quantities")')
          read_path=.TRUE.
     CASE ('mur_lc_t') 
          WRITE(stdout,'(5x,"Computing the lattice constant and the bulk" )')
          WRITE(stdout,'(5x,"modulus as a function of temperature ")')
          read_path=.TRUE.
          ltherm_expansion = .TRUE.
     CASE ('scf_elastic_constants') 
          IF (frozen_ions) THEN
             WRITE(stdout,'(5x,"Computing the frozen ions elastic constants ")')
          ELSE
             WRITE(stdout,'(5x,"Computing the elastic constants ")')
          ENDIF
          lelc = .TRUE.
     CASE ('mur_lc_elastic_constants') 
          IF (frozen_ions) THEN
             WRITE(stdout,'(5x,"Computing the frozen ions elastic constants &
                         &at the Murnaghan minimum volume")')
          ELSE
             WRITE(stdout,'(5x,"Computing the elastic constants at the &
                                  &Murnaghan minimum volume ")')
          ENDIF
          lelc = .TRUE.
     CASE ('elastic_constants_t') 
          IF (frozen_ions) THEN
             WRITE(stdout,'(5x,"Computing the frozen ions elastic constants &
                         &at all geometries")')
          ELSE
             WRITE(stdout,'(5x,"Computing the elastic constants at &
                                  &all geometries ")')
          ENDIF
          lelc = .TRUE.
     CASE ('scf_piezoelectric_tensor') 
          IF (frozen_ions) THEN
             WRITE(stdout,'(5x,"Computing the frozen ions piezoelectric tensor")')
          ELSE
             WRITE(stdout,'(5x,"Computing the piezoelectric tensor")')
          ENDIF
          lpiezo = .TRUE.
     CASE ('mur_lc_piezoelectric_tensor') 
          IF (frozen_ions) THEN
             WRITE(stdout,'(5x,"Computing the frozen ions piezoelectric tensor &
                         &at the Murnaghan minimum volume")')
          ELSE
             WRITE(stdout,'(5x,"Computing the piezoelectric tensor at the &
                                  &Murnaghan minimum volume")')
          ENDIF
          lpiezo=.TRUE.
     CASE ('scf_polarization') 
          WRITE(stdout,'(5x,"Computing the spontaneous polarization")')
     CASE ('mur_lc_polarization') 
          WRITE(stdout,'(5x,"Computing the spontaneous polarization at the &
                                  &Murnaghan minimum volume")')
     CASE ('scf_nk')
          WRITE(stdout,'(5x,"Testing the total energy convergence with k &
                         &points sampling")')
     CASE ('scf_ke')
          WRITE(stdout,'(5x,"Testing the total energy convergence with kinetic &
                         &energy cutoff ")')
     CASE ('scf_2d_bands') 
          WRITE(stdout,'(5x,"Plotting the surface band structure")')
          IF (lprojpbs) WRITE(stdout,'(5x,"Only projected band structure")')
     CASE DEFAULT
        CALL errore('thermo_summary','what not programmed',1)
  END SELECT
!
!  We now check the point group and find the Laue class, so we write
!  on output the form of the tensor that is calculated
!
  dfftp%nr1=1536
  dfftp%nr2=1536
  dfftp%nr3=1536
  fft_fact=1
  CALL find_symmetry(fft_fact)
  CALL find_group(nsym,sr,gname,code_group)
  ibrav_group_consistent=check_group_ibrav(code_group, ibrav)

  IF ( ibrav_group_consistent ) THEN
     DO isym=1,nsym
        ft(1,isym)= -DBLE(ftau(1,isym)) / dfftp%nr1
        ft(2,isym)= -DBLE(ftau(2,isym)) / dfftp%nr2
        ft(3,isym)= -DBLE(ftau(3,isym)) / dfftp%nr3
     END DO

     CALL find_space_group(ibrav, nsym, sr, ft, at, bg, sg_number, aux_sg, &
                                                        s01, s02, .FALSE.)
     IF (sg_number > 0) THEN
        unique=0
        trig=0
        IF (ibrav==-12.OR.ibrav==-13) unique=1
        IF (ibrav==5) trig=1
        CALL set_fft_fact(sg_number, unique, trig, fft_fact)
        CALL clean_dfft()
        CALL find_symmetry(fft_fact)
     ENDIF
  ELSE
     CALL clean_dfft()
     CALL find_symmetry(fft_fact)
  ENDIF

  WRITE(stdout,'(/,5x,"FFT mesh: (",i5,",",i5,",",i5," )")') dfftp%nr1, &
                                              dfftp%nr2, dfftp%nr3
  nr1_save=dfftp%nr1
  nr2_save=dfftp%nr2
  nr3_save=dfftp%nr3

  code_group_save=code_group
!
!  Description of Bravais lattice
!
  WRITE(stdout,'(/,5x,"Bravais lattice:")')
  SELECT CASE (ibrav)
     CASE(1)  
         WRITE(stdout,'(/,5x, "ibrav=1: Simple cubic lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.",/)') celldm(1)
     CASE(2)  
         WRITE(stdout,'(/,5x, "ibrav=2 Face centered cubic lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.",/)') celldm(1)
     CASE(3)  
         WRITE(stdout,'(/,5x, "ibrav=3 Body centered cubic lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.",/)') celldm(1)
     CASE(4)  
         WRITE(stdout,'(/,5x, "ibrav=4 Hexagonal lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., c/a=",f11.6,/)') celldm(1), &
                                                            celldm(3)
     CASE(5)  
         WRITE(stdout,'(/,5x, "ibrav=5 Trigonal lattice")')
         IF (ltherm_expansion) &
            CALL errore('thermo_summary','Thermal expansion not available',1)
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., cos(alpha)=",f11.6,/)') &
                                        celldm(1), celldm(4)
     CASE(6)  
         WRITE(stdout,'(/,5x, "ibrav=6 Simple tetragonal lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., c/a=",f11.6,/)') &
                                  celldm(1), celldm(3)
     CASE(7)  
         WRITE(stdout,'(/,5x, "ibrav=7 Centered tetragonal lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., c/a=",f11.6,/)') celldm(1), &
                                                            celldm(3)
     CASE(8)  
         WRITE(stdout,'(/,5x, "ibrav=8 Simple orthorombic lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                          &", c/a=",f11.6,/)') celldm(1), celldm(2), celldm(3)
     CASE(9, -9)  
         WRITE(stdout,'(/,5x, "ibrav=9 One face (C) centered orthorombic lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                            &", c/a=",f11.6,/)') celldm(1), celldm(2), celldm(3)
     CASE(91)  
         WRITE(stdout,'(/,5x, "ibrav=91 One face (A) centered orthorombic lattice")')
         IF (ltherm_expansion) &
            CALL errore('thermo_summary','Thermal expansion not available',1)
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                      &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
     CASE(10)  
         WRITE(stdout,'(/,5x, "ibrav=10 Face centered orthorombic lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                        &", c/a=",f11.6,/)') celldm(1), celldm(2), celldm(3)
     CASE(11)  
         WRITE(stdout,'(/,5x, "ibrav=11 Body centered orthorombic lattice")')
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                  &", c/a=",f11.6,/)') celldm(1), celldm(2), celldm(3)
     CASE(12,-12)  
         IF (ibrav==12) THEN
            WRITE(stdout,'(/,5x, "ibrav=12 Monoclinic lattice (c unique)")')
            WRITE(stdout,'(5x,"Starting cell parameters:")')
            WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                   &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
            WRITE(stdout,'(/,5x,"cos(gamma)=",f11.6,/)') celldm(4)
            WRITE(stdout,'(5x,"gamma=",f11.6,/)') ACOS(celldm(4))*180.0_DP/pi
         ELSE
            WRITE(stdout,'(/,5x, "ibrav=12 Monoclinic lattice (b unique)")')
            WRITE(stdout,'(5x,"Starting cell parameters:")')
            WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                      &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
            WRITE(stdout,'(/,5x,"cos(beta)=",f11.6,/)') celldm(5)
            WRITE(stdout,'(5x,"beta=",f11.6,/)') ACOS(celldm(5))*180.0_DP/pi
         ENDIF
         IF (ltherm_expansion) &
            CALL errore('thermo_summary','Thermal expansion not available',1)
     CASE(13, -13)  
         IF (ibrav==13) THEN
            WRITE(stdout,'(/,5x, "ibrav=13 Centered monoclinic lattice &
                          &(c unique)")')
            WRITE(stdout,'(5x,"Starting cell parameters:")')
            WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                      &", c/a=",f11.6,/)') celldm(1), celldm(2), celldm(3)
            WRITE(stdout,'(/,5x,"cos(gamma)=",f11.6)') celldm(4)
            WRITE(stdout,'(5x,"gamma=",f11.6,/)') ACOS(celldm(4))*180.0_DP/pi
         ELSE
            WRITE(stdout,'(/,5x, "ibrav=-13 Centered monoclinic lattice &
                          &(b unique)")')
            WRITE(stdout,'(5x,"Starting cell parameters:")')
            WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                       &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
            WRITE(stdout,'(/,5x,"cos(beta)=",f11.6)') celldm(5)
            WRITE(stdout,'(5x,"beta=",f11.6,/)') ACOS(celldm(5))*180.0_DP/pi
         ENDIF
         IF (read_path) &
         WRITE(stdout,'(/,5x, "No Brillouin Zone support. You must provide the path ")')
         IF (ltherm_expansion) &
            CALL errore('thermo_summary','Thermal expansion not available',1)
     CASE(14)  
         WRITE(stdout,'(/,5x, "ibrav=14 Triclinic lattice")')
         IF (read_path) &
         WRITE(stdout,'(/,5x, "No Brillouin Zone support. You must provide &
                               &the path ")')
         IF (ltherm_expansion) &
            CALL errore('thermo_summary','Thermal expansion not available',1)
            WRITE(stdout,'(5x,"Starting cell parameters:")')
            WRITE(stdout,'(/,5x,"alat=",f11.6," a.u., b/a=",f11.6,&
                      &", c/a=",f11.6)') celldm(1), celldm(2), celldm(3)
            WRITE(stdout,'(/,5x,"cos(alpha)=",f11.6, &
                          &", cos(beta)=",f11.6,", cos(gamma)=",f11.6)') &
                          celldm(4), celldm(5), celldm(6)
            WRITE(stdout,'(5x,"alpha=",f11.6," deg,",2x,"beta=",f11.6,&
                                           &" deg,",2x,"gamma=",f11.6," deg",/)') &
                          ACOS(celldm(4))*180.0_DP/pi,  &
                          ACOS(celldm(5))*180.0_DP/pi,  &
                          ACOS(celldm(6))*180.0_DP/pi
     CASE(0)  
         WRITE(stdout,'(/,5x, "ibrav=0 user provided cell")')
         WRITE(stdout,'(/,5x, "Be careful many options do not work with ibrav=0")')
         IF (read_path) &
            WRITE(stdout,'(/,5x, "No Brillouin Zone support. You must provide&
                          & the path ")')
         IF (ltherm_expansion) &
            CALL errore('thermo_summary','Thermal expansion not available',1)
         WRITE(stdout,'(5x,"Starting cell parameters:")')
         WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.")') celldm(1)
  CASE DEFAULT
     CALL errore('thermo_summary','ibrav not programmed',1)
  END SELECT

  WRITE(stdout,'(5x,"Starting primitive lattice vectors:")')
  WRITE( stdout, '(5x, &
      &     "crystal axes: (cart. coord. in units of alat)",/,/, &
      &       3(15x,"a(",i1,") = (",3f11.6," )  ",/ ) )')  (jpol,  &
      (at (ipol, jpol) , ipol = 1, 3) , jpol = 1, 3)
  !
  WRITE(stdout,'(5x,"Starting reciprocal lattice vectors:")')
  WRITE( stdout, '(5x, &
      &   "reciprocal axes: (cart. coord. in units 2 pi/alat)",/,/, &
      &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (jpol,&
      &  (bg (ipol, jpol) , ipol = 1, 3) , jpol = 1, 3)

  WRITE(stdout,'(5x,"Starting atomic positions in Cartesian axes:")')
  WRITE( stdout, '(/,5x,"site n.     atom                  positions (alat units)")')

  WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
             (na, atm(ityp(na)), na, (tau(ipol,na), ipol=1,3), na=1,nat)
  !
  !   allocate work space
  !
  ALLOCATE (xau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of the
  !     direct lattice vectors
  !
  DO na = 1, nat
     DO ipol = 1, 3
        xau(ipol,na) = bg(1,ipol)*tau(1,na) + bg(2,ipol)*tau(2,na) + &
                       bg(3,ipol)*tau(3,na)
     ENDDO
  ENDDO

  WRITE( stdout, '(/,5x,"Starting atomic positions in crystallographic axes:")')
  WRITE( stdout, '(/,5x,"site n.     atom        ", &
       &             "          positions (cryst. coord.)")')

  WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f11.7,"  )")') &
        (na, atm(ityp(na)), na,  (xau(ipol,na), ipol=1,3), na=1,nat)
!
!   deallocate work space
!
  DEALLOCATE(xau)
!
! ----------------------------------------------------------------------
!  Information on the symmetry and the form of the physical quantities
!
! check the compatibility of point group and Bravais lattice. If we
! have compatibility, we can use symmetry to reduce the number 
! of independent components of tensors
!
  ibrav_group_consistent=check_group_ibrav(code_group, ibrav)

  IF (what=='scf_2d_band') GOTO 1000

  IF ( ibrav_group_consistent ) THEN

     CALL find_group_info_ext(nsym, sr, code_group1, code_group_ext, &
                                                    which_elem, group_desc) 

     WRITE(stdout,'(/,5x,"The point group",i4,1x,a," is compatible &
                           &with the Bravais lattice.")') code_group_ext, &
                                    TRIM(group_name(code_group))
                                   
     WRITE(stdout,'(/,5x,"The rotation matrices with the order used inside &
                          &thermo_pw are:")') 

     CALL print_symmetries_tpw ( 1, noncolin, domag )

     DO isym=1,nsym
        ft(1,isym)= -DBLE(ftau(1,isym)) / dfftp%nr1
        ft(2,isym)= -DBLE(ftau(2,isym)) / dfftp%nr2
        ft(3,isym)= -DBLE(ftau(3,isym)) / dfftp%nr3
     END DO

     CALL find_space_group(ibrav, nsym, sr, ft, at, bg, sg_number, aux_sg, &
                                                        s01,  s02, .TRUE.)
     CALL sg_name(sg_number, 1, spaceg_name)
     IF (sg_number > 0) THEN
        WRITE(stdout,'(/,5x,"Space group ",a,"   (group number",i4, ").")') &
                            TRIM(spaceg_name), sg_number
        s01mod=s01(1)**2 + s01(2)**2 + s01(3)**2
        s02mod=s02(1)**2 + s02(2)**2 + s02(3)**2
        s01c=s01
        s02c=s02
        CALL cryst_to_cart(1,s01c,at,1)
        CALL cryst_to_cart(1,s02c,at,1)
        IF (s01mod < 1.D-4.AND.s02mod>1.D3) THEN
           WRITE(stdout,'(5x,"The origin coincides with the ITA tables.")')
        ELSEIF (s01mod<1.D-4.AND.s02mod<1.D3) THEN
           WRITE(stdout,'(5x,"The origin coincides with &
                       &the ITA tables (origin choice 1).")')
           WRITE(stdout,'(/,5x,"To shift to origin choice 2 &
                              &subtract to all coordinates:")')
           WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
           WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s02(:), s02c(:)
        ELSEIF (s01mod > 1.D-4 .AND. s02mod<1.D-4 ) THEN
           WRITE(stdout,'(5x,"The origin coincides with &
                       & the ITA tables (origin choice 2).")')
           WRITE(stdout,'(/,5x,"To shift to origin choice 1 &
                              &subtract to all coordinates:")')
           WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
           WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s01(:), s01c(:)
        ELSEIF (s01mod > 1.D-4 .AND. s01mod<1.d3 .AND. s02mod > 1.d3 ) THEN
           WRITE(stdout,'(5x,"The origin does not coincide with &
                       &the ITA tables.")')
           WRITE(stdout,'(/,5x,"To have the origin as in the ITA &
                          &tables subtract to all coordinates:")')
           WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
           WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s01(:), s01c(:)
       ELSEIF (s01mod > 1.D-4 .AND. s02mod < 1.d3 ) THEN
           WRITE(stdout,'(5x,"The origin does not coincide with the ITA tables.")')
           WRITE(stdout,'(/,5x,"To shift to origin choice 1 &
                              &subtract to all coordinates:")')
           WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
           WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s01(:), s01c(:)
           WRITE(stdout,'(/,5x,"to shift to origin choice 2 &
                             &subtract to all coordinates: ")')
           WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
           WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s02(:), s02c(:)
        ENDIF
     ELSE
        WRITE(stdout,'(/,5x,"Unknown space group.")') 
     ENDIF
!
!  first rank tensors
!
     IF ( what=='polarization'.OR. what=='mur_lc_polarization'&
                              .OR. what=='plot_bz') THEN
        SELECT CASE (ibrav)
           CASE(1,2,3)  
!
!   cubic
!
              IF (code_group==29.OR. code_group==32) THEN
                 WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
                 WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
              ELSEIF (code_group==28.OR.code_group==30.OR.code_group==31) THEN
                 WRITE(stdout,'(/,5x, "This solid has not inversion but, &
                                       &first-rank tensors, such as ")')
                 WRITE(stdout,'(5x, "the spontaneous polarization, vanish.")')
              ENDIF
           CASE(4,5,6,7)  
!
!  hexagonal, trigonal, tetragonal
!
              IF (code_group==18.OR.code_group==22.OR.code_group==23&
                                .OR.code_group==25) THEN
                 WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
                 WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
              ELSEIF (code_group==9.OR.code_group==10.OR.code_group==17.OR. &
                      code_group==21.OR.code_group==24.OR.code_group==26 ) THEN
                 WRITE(stdout,'(/,5x, "This solid has not inversion but, &
                                       &first-rank tensors, such as ")')
                 WRITE(stdout,'(5x, "the spontaneous polarization, vanish.")')
              ELSEIF (code_group==5.OR.code_group==6.OR.code_group==7.OR. &
                      code_group==13.OR.code_group==14.OR.code_group==15) THEN
                 WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
                 WRITE(stdout,'(/,5x, "(  .   .   p3 )")')
              ENDIF
           CASE(8,9,10,11)  
!
!  orthorombic
!
              IF (code_group==20) THEN
                 WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
                 WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
              ELSEIF (code_group==8) THEN
                 WRITE(stdout,'(/,5x, "This solid has not inversion but, &
                                       &first-rank tensors, such as ")')
                 WRITE(stdout,'(5x, "the spontaneous polarization, vanish.")')
              ELSEIF (code_group==12) THEN
                 WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
                 WRITE(stdout,'(/,5x, "(  .   .   p3 )")')
              ENDIF
           CASE(12,13)  
!
!   monoclinic c orientation
!
              IF (code_group==12) THEN
                 WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
                 WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
              ELSEIF (code_group==4) THEN
                 WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
                 WRITE(stdout,'(/,5x, "(  .   .   p3 )")')
              ELSEIF (code_group==3) THEN
                 WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
                 WRITE(stdout,'(/,5x, "( p1   p2   . )")')
              ENDIF
           CASE(-12)  
!
!   monoclinic b orientation
!
              IF (code_group==12) THEN
                 WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
                 WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
              ELSEIF (code_group==4) THEN
                 WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
                 WRITE(stdout,'(/,5x, "(  .   p2   .  )")')
              ELSEIF (code_group==3) THEN
                 WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
                 WRITE(stdout,'(/,5x, "(  p1   .   p3 )")')
              ENDIF
           CASE(14)  
!
!  triclinc 
!
              IF (code_group==2) THEN
                 WRITE(stdout,'(/,5x, "This solid has inversion symmetry.")')
                 WRITE(stdout,'(5x, "First-rank tensors, such as the &
                                 &spontaneous polarization, vanish.")')
              ELSEIF (code_group==1) THEN
                 WRITE(stdout,'(/,5x, "First-rank tensors, such as the &
                               &spontaneous polarization, have the form:")')
                 WRITE(stdout,'(/,5x, "( p1   p2   p3 )")')
              ENDIF
           CASE DEFAULT 
        END SELECT
     ENDIF
!
!   second rank tensors
!
     IF (what=='mur_lc_t'.OR. what=='plot_bz') THEN
        SELECT CASE (ibrav)
           CASE(1,2,3)  
!
!   cubic
!
              WRITE(stdout,'(/,5x, "Second-rank tensors, such as the dielectric")')
              WRITE(stdout,'(5x, "tensor or the thermal expansion, have the form")')
              WRITE(stdout,'(/,5x, "( e11   .    .  )")')
              WRITE(stdout,'(5x, "(  .   e11   .  )")')
              WRITE(stdout,'(5x, "(  .    .   e11 )")')
           CASE(4,5,6,7)  
!
!  hexagonal, trigonal, tetragonal
!
              WRITE(stdout,'(/,5x, "Second-rank tensors, such as the dielectric")')
              WRITE(stdout,'(5x, "tensor or the thermal expansion, have the form")')
              WRITE(stdout,'(/,5x, "( e11   .    .  )")')
              WRITE(stdout,'(5x, "(  .   e11   .  )")')
              WRITE(stdout,'(5x, "(  .    .   e33 )")')
           CASE(8,9,10,11)  
!
!  orthorombic
!
              WRITE(stdout,'(/,5x, "Second-rank tensors, such as the dielectric")')
              WRITE(stdout,'(5x, "tensor or the thermal expansion, have the form")')
              WRITE(stdout,'(/,5x, "( e11   .    .  )")')
              WRITE(stdout,'(5x, "(  .   e22   .  )")')
              WRITE(stdout,'(5x, "(  .    .   e33 )")')
           CASE(12,13)  
!
!   monoclinic c orientation
!
              WRITE(stdout,'(/,5x, "Second-rank tensors, such as the dielectric")')
              WRITE(stdout,'(5x, "tensor or the thermal expansion, have the form")')
              WRITE(stdout,'(/,5x, "( e11  e12   .  )")')
              WRITE(stdout,'(5x, "( e12  e22   .  )")')
              WRITE(stdout,'(5x, "(  .    .   e33 )")')
           CASE(-12)  
!
!   monoclinic b orientation
!
              WRITE(stdout,'(/,5x, "Second-rank tensors, such as the dielectric")')
              WRITE(stdout,'(5x, "tensor or the thermal expansion, have the form")')
              WRITE(stdout,'(/,5x, "( e11   .   e13 )")')
              WRITE(stdout,'(5x, "(  .   e22   .  )")')
              WRITE(stdout,'(5x, "( e13   .   e33 )")')
           CASE(14)  
!
!  triclinc 
!
              WRITE(stdout,'(/,5x, "Second-rank tensors, such as the dielectric")')
              WRITE(stdout,'(5x, "tensor or the thermal expansion, have the form")')
              WRITE(stdout,'(/,5x, "( e11  e12  e13 )")')
              WRITE(stdout,'(5x, "( e12  e22  e23 )")')
              WRITE(stdout,'(5x, "( e13  e23  e33 )")')
           CASE DEFAULT 
        END SELECT
     ENDIF
!
!  third rank tensor, such as the piezoelectric tensor
!
     IF ( lpiezo .OR. what=='plot_bz') THEN
        WRITE(stdout,'(/,5x,"The piezoelectric tensor is defined as the derivative &
                          &of the polarization ")')
        WRITE(stdout,'(5x,"with respect to strain (in zero electric field).")')
        SELECT CASE (code_group) 
          CASE(2,16,18,19,20,22,23,25,27,29,32) 
             WRITE(stdout,'(/,5x,"Solid with inversion symmetry.")') 
             WRITE(stdout,'(5x,"Third-rank tensors, such as the piezoelectic&
                         & tensor, vanish.")')
          CASE(3)
!
!  C_s   Monoclinic
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             IF (ibrav==-12) THEN
                WRITE(stdout,'(/,5x,"( g11  g12  g13   .   g15   .  )")') 
                WRITE(stdout,'(5x,"(  .    .    .   g24   .   g26 )")') 
                WRITE(stdout,'(5x,"( g31  g32  g33   .   g35   .  )")') 
             ELSE
                WRITE(stdout,'(/,5x,"( g11  g12  g13   .    .   g16 )")') 
                WRITE(stdout,'(5x,"( g21  g22  g23   .    .   g26 )")') 
                WRITE(stdout,'(5x,"(  .    .    .   g16  g26   .  )")') 
             ENDIF
             WRITE(stdout,'(/,5x,"It requires five strains: e1, e2, e3, e4, &
                                                            &and e5")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          5*ngeo_strain

          CASE(4)
!
!  C_2   Monoclinic
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is:")')
             IF (ibrav==-12) THEN
                WRITE(stdout,'(/,5x,"(  .    .    .   g14   .   g16 )")') 
                WRITE(stdout,'(5x,"( g21  g22  g23   .   g25   .  )")') 
                WRITE(stdout,'(5x,"(  .    .    .   g34   .   g36 )")') 
             ELSE
                WRITE(stdout,'(/,5x,"(  .    .    .   g14  g15   .  )")') 
                WRITE(stdout,'(5x,"(  .    .    .   g24  g25   .  )")') 
                WRITE(stdout,'(5x,"( g31  g32  g33   .    .   g66 )")') 
             ENDIF
             WRITE(stdout,'(/,5x,"It requires all six strains")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          4*ngeo_strain

          CASE(6,7)
!
!  C_4, tetragonal, C_6 hexagonal
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"In this class the piezoelectric tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .   g14  g15   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .   g24 -g14   .  )")') 
             WRITE(stdout,'(5x,"( g31  g31  g33   .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires four strains: e1, e3, e4, and e5")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          4*ngeo_strain

          CASE(8)
!
!  D_2 (222) Orthorombic
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   g25   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .   g36 )")') 
             WRITE(stdout,'(/,5x,"It requires two strains: e4, e5, and e6")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          3*ngeo_strain

         CASE(9)
!
! D_3  Trigonal 
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"( g11 -g11   .   g14   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .  -g14 2g11 )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires two strains: e1 and e4")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          2*ngeo_strain

         CASE(10,11)
!
! D_4  tetragonal, D_6 hexagonal
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .  -g14   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires one strain: e4")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          ngeo_strain

         CASE(12)
!
! C_2v  Orthorombic
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .    .   g15   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .   g24   .    .  )")') 
             WRITE(stdout,'(5x,"( g31  g32  g33   .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires five strains: e1, e2, e3, e4,&
                               & and e5 ")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          5*ngeo_strain

         CASE(13)
!
! C_3v  Trigonal. Assuming m perpendicular to x1
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .    .   g15 -g21 )")') 
             WRITE(stdout,'(5x,"( g21 -g21   .   g15   .    .  )")') 
             WRITE(stdout,'(5x,"( g31  g31  g33   .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires three strain: e1, e3, and e4 ")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          3*ngeo_strain

         CASE(14,15)
!
! C_4v tetragonal, C_6v hexagonal
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .    .   g15   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .   g15   .    .  )")') 
             WRITE(stdout,'(5x,"( g31  g31  g33   .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires three strain: e1, e3, and e4 ")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          3*ngeo_strain

         CASE(17)
!
! C_3h hexagonal
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"( g11 -g11   .    .    .  -g12 )")') 
             WRITE(stdout,'(5x,"( g12 -g12   .    .    .   g11 )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires one strain: e1 ")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          ngeo_strain

         CASE(21)
!
! D_3h hexagonal
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .    .    .  -g12 )")') 
             WRITE(stdout,'(5x,"( g12 -g12   .    .    .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .    .  )")') 
             WRITE(stdout,'(/,5x,"It requires one strain: e1 ")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          ngeo_strain

         CASE(24)
!
! D_2d tetragonal: axis 2 || x1
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
!
             WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   g14   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .   g34 )")') 
             WRITE(stdout,'(/,5x,"It requires two strains: e4 and e6")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          2*ngeo_strain

         CASE(26)
!
! S_4 tetragonal
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .   g14  g15   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .  -g15  g14   .  )")') 
             WRITE(stdout,'(5x,"( g31 -g31   .    .    .   g36 )")') 
             WRITE(stdout,'(/,5x,"It requires three strains: e1, e4, and e6")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          3*ngeo_strain

         CASE(28,30)
!
! T, T_d cubic
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"(  .    .    .   g14   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   g14   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .   g14 )")') 
             WRITE(stdout,'(/,5x,"It requires one strain: e4")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          ngeo_strain

         CASE(31)
             WRITE(stdout,'(/,5x,"Solid with O symmetry. The &
                             &piezoelectic tensor vanishes")')

         CASE DEFAULT
!
!  C_1 
!
             WRITE(stdout,'(/,5x,"With this point group the piezoelectric &
                                                          &tensor is")')
             WRITE(stdout,'(/,5x,"( g11  g12  g13  g14  g15  g16 )")') 
             WRITE(stdout,'(5x,"( g21  g22  g23  g24  g25  g26 )")') 
             WRITE(stdout,'(5x,"( g31  g32  g33  g34  g35  g36 )")') 
             WRITE(stdout,'(/,5x,"It requires all six strains")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          6*ngeo_strain
       END SELECT
    ENDIF
!
!  Fourth rank tensors: these are the elastic constant
! 
    laue = laue_class(code_group)
    WRITE(stdout,'(/,5x,"The Laue class is ", a)') group_name(laue)

    IF (lelc.OR.what=='plot_bz') THEN
       WRITE(stdout,'(/,5x,"In this class the elastic tensor is")') 
       SELECT CASE (laue) 
          CASE (16)
!
!    monoclinic case, class C_2h (2/m) 
!
             IF (ibrav==-12) THEN
!
!    unique axis b
!
                WRITE(stdout,'(/,5x,"( c11  c12  c13   .   c15   .  )")') 
                WRITE(stdout,'(5x,"( c12  c22  c23   .   c25   .  )")') 
                WRITE(stdout,'(5x,"( c13  c23  c33   .   c35   .  )")') 
                WRITE(stdout,'(5x,"(  .    .    .   c44   .   c46 )")') 
                WRITE(stdout,'(5x,"( c15  c25  c35   .   c55   .  )")') 
                WRITE(stdout,'(5x,"(  .    .    .   c46   .   c66 )")') 
             ELSE
!
!   unique axis c
!
                WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .   c16 )")') 
                WRITE(stdout,'(5x,"( c12  c22  c23   .    .   c26 )")') 
                WRITE(stdout,'(5x,"( c13  c23  c33   .    .   c36 )")') 
                WRITE(stdout,'(5x,"(  .    .    .   c44   .   c46 )")') 
                WRITE(stdout,'(5x,"(  .    .    .    .   c55   .  )")') 
                WRITE(stdout,'(5x,"( c15  c25  c36  c46   .   c66 )")') 
             ENDIF 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires all six strains")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          6*ngeo_strain 
             ENDIF
          CASE (20)
!
!  orthorombic D_2h (mmm)
!
             WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c12  c22  c23   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c13  c23  c33   .    .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .   c44   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   c55   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .   c66 )")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires all six strains")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          6*ngeo_strain 
             ENDIF
          CASE (18)
!
!  tetragonal C_4h (4/m)
!
             WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .   c16 )")') 
             WRITE(stdout,'(5x,"( c12  c11  c13   .    .  -c16 )")') 
             WRITE(stdout,'(5x,"( c13  c13  c33   .    .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .   c44   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   c44   .  )")') 
             WRITE(stdout,'(5x,"( c16 -c16   .    .    .   c66 )")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
             WRITE(stdout,'(/,5x,"It requires four strains: e1, e3, e4, e6")') 
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          4*ngeo_strain 
             ENDIF

          CASE (22)
!
!  tetragonal D_4h (4/mmm)
!
             WRITE(stdout,'(/,5x,"( c11  c12  c13   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c12  c11  c13   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c13  c13  c33   .    .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .   c44   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   c44   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .   c66 )")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires four strains: e1, e3, &
                                                &e4, e6")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          4*ngeo_strain 
             ENDIF
          CASE (27)
!
!  trigonal S_6 (-3)
!
             WRITE(stdout,'(5x,"( c11  c12  c13  c14  c15   .  )")') 
             WRITE(stdout,'(5x,"( c12  c11  c13 -c14 -c15   .  )")') 
             WRITE(stdout,'(5x,"( c13  c13  c33   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c14 -c14   .   c44   .  -c15 )")') 
             WRITE(stdout,'(5x,"( c15 -c15   .    .   c44  c14 )")') 
             WRITE(stdout,'(5x,"(  .    .    .  -c1   c14   X  )")') 
             WRITE(stdout,'(5x,"X=(c11-c12)/2")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires three strains: e1, e3, &
                                                             &and e4")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          3*ngeo_strain 
             ENDIF
          CASE (25)
!
!  trigonal D_3d (-3m)
!

             WRITE(stdout,'(5x,"( c11  c12  c13  c14   .    .  )")') 
             WRITE(stdout,'(5x,"( c12  c11  c13 -c14   .    .  )")') 
             WRITE(stdout,'(5x,"( c13  c13  c33   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c14 -c14   .   c44   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   c44  c14 )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   c14   X  )")') 
             WRITE(stdout,'(5x,"X=(c11-c12)/2")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires three strains: e1, e3, &
                                                          &and e4")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          3*ngeo_strain 
             ENDIF

          CASE (19,23)
!
!  hexagonal C_6h (6/m), D_6h (6/mmm)
!

             WRITE(stdout,'(5x,"( c11  c12  c13   .    .    . )")') 
             WRITE(stdout,'(5x,"( c12  c11  c13   .    .    . )")') 
             WRITE(stdout,'(5x,"( c13  c13  c33   .    .    . )")') 
             WRITE(stdout,'(5x,"(  .    .    .   c44   .    . )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   c44   . )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .    X )")') 
             WRITE(stdout,'(5x,"X=(c11-c12)/2")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires three strains: e1, e3, &
                                                                   and e4")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          3*ngeo_strain 
             END IF
          CASE (29,32)
!
!  cubic T_h (m-3), O_h (m-3m)
!
             WRITE(stdout,'(/,5x,"( c11  c12  c12   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c12  c11  c12   .    .    .  )")') 
             WRITE(stdout,'(5x,"( c12  c12  c11   .    .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .   c44   .    .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .   c44   .  )")') 
             WRITE(stdout,'(5x,"(  .    .    .    .    .   c44 )")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires two strains: e1 and e4")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          2*ngeo_strain 
             ENDIF
          CASE DEFAULT
             IF (laue /=2) &
                WRITE(stdout,'(5x,"Laue class not programmed using C_i")')   
             WRITE(stdout,'(/,5x,"( c11  c12  c13  c14  c15  c16 )")') 
             WRITE(stdout,'(5x,"( c12  c22  c23  c24  c25  c26 )")') 
             WRITE(stdout,'(5x,"( c13  c23  c33  c34  c35  c36 )")') 
             WRITE(stdout,'(5x,"( c14  c24  c34  c44  c45  c46 )")') 
             WRITE(stdout,'(5x,"( c15  c25  c35  c45  c55  c56 )")') 
             WRITE(stdout,'(5x,"( c16  c26  c36  c46  c56  c66 )")') 
             IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
                WRITE(stdout,'(/,5x,"It requires all six strains")') 
                WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                          6*ngeo_strain 
             ENDIF
       END SELECT
    ENDIF
 ELSE
!
!  In this case the Bravais lattice and the point group are not consistent.
!  Usually this means that the user has used ibrav=0, or the solid has
!  too low symmetry. Since one cannot assume the direction of the z axis 
!  as the high symmetry axis, I skip completely the use of symmetry for
!  this case
!
    CALL print_symmetries_tpw ( 1, noncolin, domag )
    WRITE(stdout,'(/,5x,"ibrav=0 or Bravais lattice not compatible with &
                                    &the point group.")')
    WRITE(stdout,'(/,5x,"I will not use symmetry.")')
    WRITE(stdout,'(/,5x,"Cannot use the Laue class with ibrav=0, &
                                                &using laue=0")')
    laue=0

    IF ( what=='polarization'.OR. what=='mur_lc_polarization'&
                             .OR. what=='plot_bz') THEN
!
!  for first rank tensor one can still check for the existence
!  of inversion symmetry
!
       SELECT CASE (code_group) 
          CASE(2,16,18,19,20,22,23,25,27,29,32) 
             WRITE(stdout,'(/,5x,"Solid with inversion symmetry.")')
             WRITE(stdout,'(5x,"First-rank tensors, such as the spontaneous &
                           &polarization, vanish.")')
          CASE DEFAULT
             WRITE(stdout,'(/,5x,"Solid without inversion symmetry.")')
       END SELECT
    ENDIF
!
!  second rank tensor. All components are calculated
!
    IF (what=='mur_lc_t'.OR. what=='plot_bz') THEN

       WRITE(stdout,'(/,5x, "All components of second order tensors such as")')
       WRITE(stdout,'(5x, "the dielectric tensor or the thermal expansion are &
                                          &calculated ")')

       WRITE(stdout,'(/,5x, "( e11  e12  e13 )")')
       WRITE(stdout,'(5x, "( e12  e22  e23 )")')
       WRITE(stdout,'(5x, "( e13  e23  e33 )")')
    ENDIF
!
!  third rank tensor
!
    IF ( lpiezo .OR. what=='plot_bz') THEN
       SELECT CASE (code_group) 
          CASE(2,16,18,19,20,22,23,25,27,29,32) 
             WRITE(stdout,'(/,5x,"Solid with inversion symmetry.")')
             WRITE(stdout,'(5x,"Third-rank tensors, such as the piezoelectic &
                                &tensor, vanish.")')
          CASE DEFAULT

             WRITE(stdout,'(/,5x,"I will take a piezoelectric tensor of the &
                                                                      form")')
             WRITE(stdout,'(/,5x,"( g11  g12  g13  g14  g15  g16 )")') 
             WRITE(stdout,'(5x,"( g21  g22  g23  g24  g25  g26 )")') 
             WRITE(stdout,'(5x,"( g31  g32  g33  g34  g35  g36 )")') 
             WRITE(stdout,'(/,5x,"It requires all six strains")')
             WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
                       6*ngeo_strain
       END SELECT
    ENDIF
!
!  Fourth rank tensors
!
    IF (lelc.OR.what=='plot_bz') THEN
       WRITE(stdout,'(/,5x,"I will use elastic constants with the form")')
       WRITE(stdout,'(/,5x,"( c11  c12  c13  c14  c15  c16 )")') 
       WRITE(stdout,'(5x,"( c12  c22  c23  c24  c25  c26 )")') 
       WRITE(stdout,'(5x,"( c13  c23  c33  c34  c35  c36 )")') 
       WRITE(stdout,'(5x,"( c14  c24  c34  c44  c45  c46 )")') 
       WRITE(stdout,'(5x,"( c15  c25  c35  c45  c55  c56 )")') 
       WRITE(stdout,'(5x,"( c16  c26  c36  c46  c56  c66 )")') 
       IF (elastic_algorithm=='standard'.OR. &
                                    elastic_algorithm=='advanced') THEN
          WRITE(stdout,'(/,5x,"It requires all six strains")') 
          WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
             6*ngeo_strain 
       END IF
    ENDIF
END IF

IF (what(1:6)=='mur_lc') &
   WRITE(stdout,'(5x,"The Murnaghan relaxation will require", &
                     &i3, " scf calculations")') ngeo(1)*ngeo(2)*ngeo(3)* &
                            ngeo(4)*ngeo(5)*ngeo(6)

WRITE(stdout,'(/,5x,70("-"))')
IF (frozen_ions) THEN
   WRITE(stdout,'(5x,"Ions are not relaxed")')
ELSE
   WRITE(stdout,'(5x,"Ions are relaxed in each calculation")')
ENDIF
WRITE(stdout,'(5x,70("-"))')
!
!  ----------------------------------------------------------------------
!  Information on the density
!
  CALL compute_density(omega,density)
!
!  ----------------------------------------------------------------------
!  Brillouin zone plot
!
1000  CONTINUE
  IF (what=='plot_bz') THEN
     asy_filename=TRIM(flasy)//'.asy'
     IF ( my_image_id==root_image ) THEN
        IF (ibrav /= 13 .AND. ibrav /= -13 .AND. ibrav /=91 .AND. &
                        ibrav /=14 .AND. ibrav /=0) THEN 
           CALL plot_bz(ibrav, celldm, at, bg, point_label_type, &
                xqaux, wqaux, nqaux, letter, letter_path, npk_label, &
                label_list, asy_filename)

           IF (lasymptote.AND.ionode) &
               ierr=system(TRIM(asymptote_command)//' '//TRIM(asy_filename))
        ENDIF
!
!  Form factors of the atoms

        WRITE(stdout,'(/,5x,"Computing the X-ray powder diffraction intensities")')
        IF (lformf) THEN
           DO it=1,nsp
              IF (lcm) THEN
                WRITE(stdout,'(/,5x,"Form factors from Cromer-Mann &
                             &parameters")')
              ELSE
                 WRITE(stdout,'(/,5x,"Form factors from Doyle-Turner &
                             &or Smith-Burge parameters")')
              END IF
              filename=TRIM(flformf)//'.'//TRIM(int_to_char(it))
              CALL write_form_factor(atm(it),smin,smax,nspoint,lcm,filename)
           END DO
           CALL plot_form_factors()
        END IF

        CALL compute_xrdp(at,bg,celldm(1),nat,tau,nsp,ityp,atm, &
                               lambda,ibrav,lcm,flxrdp)
        CALL plot_xrdp('')
     ENDIF
     CALL summarize_kpt(xqaux, wqaux, nqaux, letter_path)

     iuout=35
     xsf_filename=TRIM(prefix)//'.xsf'
     OPEN(unit=iuout, file=xsf_filename, status='unknown', &
                                            form='formatted')

     CALL xsf_struct (celldm(1), at, nat, tau, atm, ityp, iuout)

     CLOSE(iuout)


     CALL environment_end( 'THERMO_PW' )
     !
     CALL mp_global_end ()
     CALL do_stop( 0 )
  ELSEIF (what=='scf_2d_bands') THEN
     !
     asy_filename=TRIM(flasy)//'.asy'
     IF ( my_image_id==root_image ) THEN

        CALL plot_2d_bz(ibrav, celldm, at, bg, xqaux, wqaux, nqaux, &
             letter, letter_path, npk_label, label_list, asy_filename)

        IF (lasymptote.AND.ionode) &
           ierr=system(TRIM(asymptote_command)//' '//TRIM(asy_filename))

     ENDIF
     !
  ENDIF

  RETURN
END SUBROUTINE thermo_summary

FUNCTION check_group_ibrav(code_group, ibrav)
!
!  This routine checks if the ibrav is compatible with the point group,
!  and if it is not if write a brief message saying which lattices are
!  compatible
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
USE lattices, ONLY : lattice_name
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, ibrav

INTEGER :: is_compatible(32,6)
LOGICAL :: check_group_ibrav
INTEGER :: i
CHARACTER(LEN=40) :: latt_name
CHARACTER(LEN=11) :: group_name

check_group_ibrav=.FALSE.
is_compatible=0
!
!   C_1, C_i triclinic
!
is_compatible(1,1)=14
is_compatible(2,1)=14
!
!   C_s, C_2, monoclinic
!
is_compatible(3,1)=12
is_compatible(3,2)=13
is_compatible(3,3)=-12
is_compatible(3,4)=-13
is_compatible(4,1)=12
is_compatible(4,2)=13
is_compatible(4,3)=-12
is_compatible(4,4)=-13
!
!   C_3, trigonal or hexagonal
!
is_compatible(5,1)=4
is_compatible(5,2)=5
is_compatible(5,3)=-5
!
!   C_4, tetragonal
!
is_compatible(6,1)=6
is_compatible(6,2)=7
!
!   C_6, hexagonal
!
is_compatible(7,1)=4
!
!   D_2, orthorombic
!
is_compatible(8,1)=8
is_compatible(8,2)=9
is_compatible(8,3)=10
is_compatible(8,4)=11
is_compatible(8,5)=-9
is_compatible(8,6)=91
!
!   D_3 trigonal or hexagonal
!
is_compatible(9,1)=4
is_compatible(9,2)=5
is_compatible(9,3)=-5
!
!   D_4 tetragonal
!
is_compatible(10,1)=6
is_compatible(10,2)=7
!
!   D_6 hexagonal
!
is_compatible(11,1)=4
!
!   C_2v orthorombic
!
is_compatible(12,1)=8
is_compatible(12,2)=9
is_compatible(12,3)=10
is_compatible(12,4)=11
is_compatible(12,5)=-9
is_compatible(12,6)=91
!
!   C_3v hexagonal or trigonal
!
is_compatible(13,1)=4
is_compatible(13,2)=5
is_compatible(13,3)=-5
!
!   C_4v tetragonal
!
is_compatible(14,1)=6
is_compatible(14,2)=7
!
!   C_6v hexagonal
!
is_compatible(15,1)=4
!
!   C_2h monoclinic
!
is_compatible(16,1)=12
is_compatible(16,2)=13
is_compatible(16,3)=-12
is_compatible(16,4)=-13
!
!  C_3h hexagonal
!
is_compatible(17,1)=4
!
!  C_4h tetragonal
!
is_compatible(18,1)=6
is_compatible(18,2)=7
!
!  C_6h hexagonal
!
is_compatible(19,1)=4
!
!  D_2h orthorombic
!
is_compatible(20,1)=8
is_compatible(20,2)=9
is_compatible(20,3)=10
is_compatible(20,4)=11
is_compatible(20,5)=-9
is_compatible(20,6)=91
!
!  D_3h hexagonal
!
is_compatible(21,1)=4
!
!  D_4h tetragonal
!
is_compatible(22,1)=6
is_compatible(22,2)=7
!
!  D_6h hexagonal
!
is_compatible(23,1)=4
!
!  D_2d tetragonal
!
is_compatible(24,1)=6
is_compatible(24,2)=7
!
!   D_3d hexagonal or trigonal
!
is_compatible(25,1)=4
is_compatible(25,2)=5
is_compatible(25,3)=-5
!
!   S_4 tetragonal
!
is_compatible(26,1)=6
is_compatible(26,2)=7
!
!   S_6 hexagonal or trigonal
!
is_compatible(27,1)=4
is_compatible(27,2)=5
is_compatible(27,3)=-5
!
!   T cubic
!
is_compatible(28,1)=1
is_compatible(28,2)=2
is_compatible(28,3)=3
!
!   T_h cubic
!
is_compatible(29,1)=1
is_compatible(29,2)=2
is_compatible(29,3)=3
!
!   T_d cubic
!
is_compatible(30,1)=1
is_compatible(30,2)=2
is_compatible(30,3)=3
!
!   O cubic
!
is_compatible(31,1)=1
is_compatible(31,2)=2
is_compatible(31,3)=3
!
!   O_h cubic
!
is_compatible(32,1)=1
is_compatible(32,2)=2
is_compatible(32,3)=3

IF (ibrav==0) THEN
   WRITE(stdout,'(5x,"ibrav=0, many features are not implemented")')
   WRITE(stdout,'(5x,"The point group ",a11," is compatible with:")') &
                                                      group_name(code_group)
   DO i=1,6
      IF (is_compatible(code_group,i) /=0) THEN
         CALL lattice_name(is_compatible(code_group,i),latt_name)
         WRITE(stdout,'(5x,a)') TRIM(latt_name)
      ENDIF
   ENDDO
ELSE
   DO i=1,6
      IF (is_compatible(code_group,i) == ibrav) THEN
         check_group_ibrav=.TRUE.
         GOTO 100
      ENDIF
   ENDDO
   CALL lattice_name(ibrav,latt_name)
   WRITE(stdout,'(/,5x,a," is incompatible with the ",&
                      & a," Bravais lattice")') &
                      TRIM(group_name(code_group)), TRIM(latt_name)
   WRITE(stdout,'(5x,"It is compatible with the ")') 
   DO i=1,6
      IF (is_compatible(code_group,i) /=0) THEN
         CALL lattice_name(is_compatible(code_group,i),latt_name)
         WRITE(stdout,'(5x,a," Bravais lattice; ibrav=",i5)') TRIM(latt_name), &
                   is_compatible(code_group,i)
      ENDIF
   ENDDO
   WRITE(stdout,'(/,5x,"You might want to change the Bravais lattice or to")') 
   WRITE(stdout,'(/,5x,"understand why the symmetries are wrong before&
                      &continuing")') 
   WRITE(stdout,'(5x,"The point group or the Laue class are not used to &
                      &reduce the number of ")')
   WRITE(stdout,'(5x,"computed tensor components")') 
100 CONTINUE
END IF

RETURN
END FUNCTION check_group_ibrav

SUBROUTINE find_fft_fact()
!
!  this routine finds the fft_fact that correspond to the space group
!  of the solid. It assumes that the tau and celldm are already known
!  It sets also the dimension of the fft that correspond to the cut-off.
!
USE kinds,            ONLY : DP
USE fft_base,         ONLY : dfftp
USE cell_base,        ONLY : ibrav, at, bg
USE thermo_sym,       ONLY : fft_fact, ibrav_group_consistent
USE rap_point_group,  ONLY : code_group
USE space_groups,     ONLY : find_space_group, set_fft_fact
USE symm_base,        ONLY : nsym, s, sr, ftau

IMPLICIT NONE
INTEGER :: sg_number
INTEGER :: unique, trig, isym, aux_sg
LOGICAL :: check_group_ibrav
REAL(DP) :: s01(3), s02(3), ft(3,48)
CHARACTER(LEN=12) :: spaceg_name
CHARACTER(LEN=11) :: gname


  dfftp%nr1=1536
  dfftp%nr2=1536
  dfftp%nr3=1536
  fft_fact=1
  CALL find_symmetry(fft_fact)
  CALL find_group(nsym,sr,gname,code_group)
  ibrav_group_consistent=check_group_ibrav(code_group, ibrav)

  IF ( ibrav_group_consistent ) THEN
     DO isym=1,nsym
        ft(1,isym)= -DBLE(ftau(1,isym)) / dfftp%nr1
        ft(2,isym)= -DBLE(ftau(2,isym)) / dfftp%nr2
        ft(3,isym)= -DBLE(ftau(3,isym)) / dfftp%nr3
     END DO
     CALL find_space_group(ibrav, nsym, sr, ft, at, bg, sg_number, aux_sg, &
                                  s01, s02, .FALSE.)

     IF (sg_number > 0) THEN
        unique=0
        trig=0
        IF (ibrav==-12.OR.ibrav==-13) unique=1
        IF (ibrav==5) trig=1
        CALL set_fft_fact(sg_number, unique, trig, fft_fact)
        CALL clean_dfft()
        CALL find_symmetry(fft_fact)
     ENDIF
  ELSE
     CALL clean_dfft()
     CALL find_symmetry(fft_fact)
  ENDIF

  RETURN
END SUBROUTINE find_fft_fact
