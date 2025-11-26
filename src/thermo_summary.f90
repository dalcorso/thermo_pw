!
! Copyright (C) 2014-2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE thermo_summary()
  !-----------------------------------------------------------------------
  !
  !  This routine writes a summary of the input. It writes what will
  !  be calculated, a few information on the reference configuration
  !  read in the pw.x input, it identifies the space group and writes the
  !  form of the physical tensor to calculate. When what='plot_bz', it 
  !  plots the Brillouin zone and gives information about the BZ path, 
  !  it writes an xsf file of the structure, and the X-ray power diffraction 
  !  spectrum.
  !
  !  This routine sets the following variables:
  !  laue      ! the laue class
  !  ibrav_group_consistent   ! consistency between bravais lattice and point
  !                             group
  !  sg_number ! the number of the space group
  !  aux_sg    ! the space group orientation
  !  
  !
  USE kinds,                ONLY : DP
!
!  summarized control variables
!
  USE thermo_mod,           ONLY : what, ngeo, density
  USE control_paths,        ONLY : nqaux, xqaux, wqaux, npk_label, letter, &
                                   letter_path, label_list, point_label_type
  USE control_asy,          ONLY : flasy, lasymptote, asymptote_command
  USE control_elastic_constants, ONLY : frozen_ions, ngeo_strain, &
                                   elastic_algorithm, elalgen, use_free_energy
  USE control_xrdp,         ONLY : lambda, flxrdp, lformf, &
                                   flformf, smin, smax, nspoint, lcm
  USE xrdp_module,          ONLY : write_form_factor, compute_xrdp
  USE control_2d_bands,     ONLY : lprojpbs
  USE control_thermo,       ONLY : with_eigen, ltherm_dos, ltherm_freq, &
                                   ltherm_glob
!
!  variables set by this routine
!
  USE thermo_sym,           ONLY : laue, ibrav_group_consistent, sg_number, &
                                   aux_sg
!
!  library modules
!
  USE space_groups,         ONLY : sg_name, find_space_group, sg_origin
  USE lattices,             ONLY : print_bravais_description
  USE point_group,          ONLY : find_group_info_ext
  USE nye,                  ONLY : print_vectors_shape, print_tensor2_shape, &
                                   print_el_cons_shape, print_piezo_shape, &
                                   print_b_fact_shape
  USE elastic_constants,    ONLY : print_el_cons_info
  USE piezoelectric_tensor, ONLY : print_piezo_info
!
!  pw modules and variables, set from input or in thermo_setup
!
  USE rap_point_group,      ONLY : code_group
  USE fft_base,             ONLY : dfftp
  USE noncollin_module,     ONLY : noncolin, domag
  USE cell_base,            ONLY : ibrav, at, bg, celldm, omega
  USE ions_base,            ONLY : tau, nat, ityp, nsp, atm
  USE symm_base,            ONLY : nsym, sr, ft

  USE mp_images,            ONLY : my_image_id, root_image
  USE io_global,            ONLY : ionode, stdout
  USE io_files,             ONLY : prefix
  !
  IMPLICIT NONE
  REAL(DP) :: ft_(3,48)
  REAL(DP) :: s01(3), s02(3)
  REAL(DP), ALLOCATABLE :: xau(:,:)
  INTEGER :: it, ia, na, ipol, jpol, iuout, &
             group_desc(48), which_elem(48), isym, code_group_ext,     &
             code_group1
  INTEGER :: laue_class, b_rest(nat)
  INTEGER :: find_free_unit
  INTEGER :: ierr, system
  LOGICAL :: lpolar, lelc, lpiezo, check_group_ibrav
  CHARACTER(LEN=12)  :: spaceg_name
  CHARACTER(LEN=11)  :: gname, group_name
  CHARACTER(LEN=6)   :: int_to_char
  CHARACTER(LEN=256) :: asy_filename, filename, xsf_filename

  IF (what==' ') RETURN
!
!  Summarize what will be calculated
!
  lelc = .FALSE.
  lpiezo=.FALSE.
  lpolar=.FALSE.
  WRITE(stdout,'(/)')
  SELECT CASE (TRIM(what))
     CASE ('plot_bz') 
          WRITE(stdout,'(5x,"Plotting the Brillouin Zone and k point path")')
     CASE ('scf') 
          WRITE(stdout,'(5x,"Doing a single scf calculation")')
     CASE ('scf_bands') 
          WRITE(stdout,'(5x,"Doing a band calculation")')
          WRITE(stdout,'(5x,"Use what=''plot_bz'' to visualize the BZ path")')
     CASE ('scf_dos')
          WRITE(stdout,'(5x,"Doing a electronic dos calculation")')
     CASE ('scf_ph') 
          WRITE(stdout,'(5x,"Doing a phonon calculation")')
     CASE ('scf_disp')
          WRITE(stdout,'(5x,"Doing a phonon dispersion calculation")')
          WRITE(stdout,'(5x,"Use what=''plot_bz'' to visualize the BZ path")')
          WRITE(stdout,'(5x,"Computing the harmonic thermodynamic properties")')
     CASE ('mur_lc') 
          WRITE(stdout,'(5x,"Calculating the volume that minimizes the &
                                                                  &energy")')
     CASE ('mur_lc_bands') 
          WRITE(stdout,'(5x,"Calculating the bands at the minimum volume")')
          WRITE(stdout,'(5x,"Use what=''plot_bz'' to visualize the BZ path")')
     CASE ('mur_lc_dos')
          WRITE(stdout,'(5x,"Calculating the electronic dos at the minimum &
                                                                   &volume")')
     CASE ('mur_lc_ph') 
          WRITE(stdout,'(5x,"Doing a phonon calculation at the minimum &
                                                                   &volume")')
     CASE ('mur_lc_disp') 
          WRITE(stdout,'(5x,"Doing a phonon dispersion calculation at the &
                                           & minimum volume")')
          WRITE(stdout,'(5x,"Use what=''plot_bz'' to visualize the BZ path")')
          WRITE(stdout,'(5x,"Computing the harmonic thermodynamic &
                                                               &properties")')
     CASE ('mur_lc_t') 
          WRITE(stdout,'(5x,"Computing the lattice constant and the bulk" )')
          WRITE(stdout,'(5x,"modulus as a function of temperature ")')
     CASE ('scf_elastic_constants') 
          WRITE(stdout,'(5x,"Computing the elastic constants ")')
          IF (frozen_ions) WRITE(stdout,'(5x,"The ions are frozen" )')
          lelc = .TRUE.
     CASE ('mur_lc_elastic_constants') 
          WRITE(stdout,'(5x,"Computing the elastic constants at the &
                                  &minimum volume ")')
          IF (frozen_ions) WRITE(stdout,'(5x,"The ions are frozen" )')
          lelc = .TRUE.
     CASE ('elastic_constants_geo')
          WRITE(stdout,'(5x,"Computing the elastic constants" )')
          IF (use_free_energy) THEN
             WRITE(stdout,'(5x,"within quasi-harmonic approximation" )')
          ELSE
             WRITE(stdout,'(5x,"at zero temperature" )')
          ENDIF
          WRITE(stdout,'(5x,"for several geometries" )')
          IF (frozen_ions) WRITE(stdout,'(5x,"The ions are frozen" )')
          lelc = .TRUE.
     CASE ('scf_piezoelectric_tensor') 
          WRITE(stdout,'(5x,"Computing the piezoelectric tensor")')
          IF (frozen_ions) WRITE(stdout,'(5x,"The ions are frozen" )')
          lpiezo = .TRUE.
     CASE ('mur_lc_piezoelectric_tensor') 
          WRITE(stdout,'(5x,"Computing the piezoelectric tensor at the &
                                  &minimum volume")')
          IF (frozen_ions) WRITE(stdout,'(5x,"The ions are frozen" )')
          lpiezo=.TRUE.
     CASE ('piezoelectric_tensor_geo')
          WRITE(stdout,'(5x,"Computing the piezoelectric tensor" )')
          WRITE(stdout,'(5x,"for several geometries" )')
          IF (frozen_ions) WRITE(stdout,'(5x,"The ions are frozen" )')
          lpiezo = .TRUE.
     CASE ('scf_polarization') 
          WRITE(stdout,'(5x,"Computing the spontaneous polarization")')
          lpolar=.TRUE.
     CASE ('mur_lc_polarization') 
          WRITE(stdout,'(5x,"Computing the spontaneous polarization at the &
                                  &minimum volume")')
          lpolar=.TRUE.
     CASE ('polarization_geo')
          WRITE(stdout,'(5x,"Computing the Berry phase polarization" )')
          WRITE(stdout,'(5x,"for several geometries" )')
          IF (frozen_ions) WRITE(stdout,'(5x,"The ions are frozen" )')
          lpolar = .TRUE.
     CASE ('scf_nk')
          WRITE(stdout,'(5x,"Testing the total energy convergence with the k &
                         &points mesh")')
     CASE ('scf_ke')
          WRITE(stdout,'(5x,"Testing the total energy convergence with the &
                         &kinetic energy cutoff ")')
     CASE ('scf_2d_bands') 
          WRITE(stdout,'(5x,"Plotting the surface band structure")')
          IF (lprojpbs) WRITE(stdout,'(5x,"Plotting the projected band &
                                                              &structure")')
     CASE DEFAULT
        CALL errore('thermo_summary','what not programmed',1)
  END SELECT

  WRITE(stdout,'(/,5x,"FFT mesh: (",i5,",",i5,",",i5," )")') dfftp%nr1, &
                                                  dfftp%nr2, dfftp%nr3
!
!  Description of the Bravais lattice and of the reciprocal lattice
!
  WRITE(stdout,'(/,5x,"Bravais lattice:")')
  IF (ibrav/=0) THEN
     CALL print_bravais_description(ibrav,celldm)
  ELSE  
     WRITE(stdout,'(/,5x, "ibrav=0 user provided cell")')
     WRITE(stdout,'(5x, "Be careful many options will not work &
                                                           &with ibrav=0")')
     WRITE(stdout,'(/,5x,"Starting cell parameters:")')
     WRITE(stdout,'(/,5x,"alat=",f11.6," a.u.")') celldm(1)
  ENDIF

  WRITE(stdout,'(/,5x,"Starting primitive lattice vectors:")')
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
!
!  Description of the atomic positions
!
  WRITE(stdout,'(5x,"Starting atomic positions in Cartesian axes:")')
  WRITE( stdout, '(/,5x,"site n.     atom                  positions (alat units)")')

  WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
             (na, atm(ityp(na)), na, (tau(ipol,na), ipol=1,3), na=1,nat)
!
!     Compute the coordinates of the atoms in crystal basis and write them
!
  ALLOCATE (xau(3,nat))
  xau=tau
  CALL cryst_to_cart(nat,xau,bg,-1)

  WRITE( stdout, '(/,5x,"Starting atomic positions in crystallographic axes:")')
  WRITE( stdout, '(/,5x,"site n.     atom        ", &
       &             "          positions (cryst. coord.)")')

  WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f11.7,"  )")') &
        (na, atm(ityp(na)), na,  (xau(ipol,na), ipol=1,3), na=1,nat)
  DEALLOCATE(xau)
!
!   Write a few information on the energy minimization
!
  IF (what(1:6)=='mur_lc') &
     WRITE(stdout,'(/,5x,"The energy minimization will require", &
                     &i3, " scf calculations")') ngeo(1)*ngeo(2)*ngeo(3)* &
                            ngeo(4)*ngeo(5)*ngeo(6)

!
! ----------------------------------------------------------------------
!  Information on the symmetry and the form of the physical tensors
!
!
! for surface band structure calculation the information on physical
! quantities is not written on output
!
  IF (what=='scf_2d_band') GOTO 1000

!
! check the compatibility of point group and Bravais lattice. If 
! they are compatibile, we use symmetry to reduce the number 
! of independent components of tensors
!
  ibrav_group_consistent=check_group_ibrav(code_group, ibrav)

  IF ( ibrav_group_consistent ) THEN

     CALL find_group_info_ext(nsym, sr, code_group1, code_group_ext, &
                                                    which_elem, group_desc) 

     WRITE(stdout,'(/,5x,"The point group",i4,1x,a," is compatible &
                           &with the Bravais lattice.")') code_group_ext, &
                                    TRIM(group_name(code_group))
                                   
     WRITE(stdout,'(/,5x,"The rotation matrices with the order used inside &
                          &thermo_pw are:")') 

     CALL print_symmetries_tpw ( 1, noncolin, domag )
!
!  The space group is identified another time here, because these
!  are the fft factors that pw.x will use. It should coincide with the
!  previously identified one.
!
     ft_(:,1:nsym)=-ft(:,1:nsym)

     CALL find_space_group(ibrav, nsym, sr, ft_, at, bg, sg_number, aux_sg, &
                                                        s01,  s02, .TRUE.)

     CALL sg_name(sg_number, 1, spaceg_name)
     CALL sg_origin(sg_number, spaceg_name, at, s01, s02 ) 
!
!---------------------------------------------------------------------------
!  Write information on the tensors to compute
!
!  first rank tensors
!
     IF ( lpolar .OR. what=='plot_bz') CALL print_vectors_shape(code_group,&
                                                                        ibrav)
!
!   second rank tensors
!
     IF (what=='mur_lc_t'.OR. what=='plot_bz') THEN
        WRITE(stdout,'(/,5x, "Second-rank tensors, such as the dielectric")')
        WRITE(stdout,'(5x, "tensor or the thermal expansion, have the form:")')
        CALL print_tensor2_shape(ibrav)
     ENDIF
!
!  third rank tensor, such as the piezoelectric tensor
!
     IF ( lpiezo .OR. what=='plot_bz') THEN
        WRITE(stdout,'(/,5x,"The piezoelectric tensor, defined as the &
                             &derivative of the polarization ")')
        WRITE(stdout,'(5x,"with respect to strain (in zero electric field) &
                                                                 &is:")')
        CALL print_piezo_shape(code_group,ibrav)
        CALL print_piezo_info(code_group,ibrav,ngeo_strain)
     ENDIF
!
!  Fourth rank tensors: elastic constants
! 
    laue = laue_class(code_group)
    WRITE(stdout,'(/,5x,"The Laue class is ", a)') group_name(laue)

    IF (lelc.OR.what=='plot_bz') THEN
       WRITE(stdout,'(/,5x,"In this class the elastic tensor is")') 
       CALL print_el_cons_shape(laue,ibrav)
       CALL print_el_cons_info(elastic_algorithm, laue, ibrav, ngeo_strain)
    ENDIF
 ELSE
!
!  In this case the Bravais lattice and the point group are not consistent.
!  Usually this means that the user has used ibrav=0, or the solid has
!  too low symmetry for the declared Bravais lattice and a different
!  lattice should be used. Since one cannot assume the direction of the z axis 
!  as the high symmetry axis, the code skips completely the use of symmetry for
!  this case
!
    CALL print_symmetries_tpw ( 1, noncolin, domag )
    WRITE(stdout,'(/,5x,"ibrav=0 or Bravais lattice not compatible with &
                                    &the point group.")')
    WRITE(stdout,'(/,5x,"The code will not use symmetry.")')
    WRITE(stdout,'(/,5x,"ibrav=0 use Laue class 2")')

    laue=2

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
                                          &calculated:")')

       CALL print_tensor2_shape(0)
    ENDIF
!
!  third rank tensor
!
    IF ( lpiezo .OR. what=='plot_bz') THEN
       SELECT CASE (code_group) 
          CASE(2,16,18,19,20,22,23,25,27,29,32) 
             CALL print_piezo_shape(code_group,ibrav)
          CASE DEFAULT
             CALL print_piezo_shape(code_group,ibrav)
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
       CALL print_el_cons_shape(laue,ibrav)
       IF (.NOT.elalgen) THEN
          WRITE(stdout,'(/,5x,"It requires all six strains")') 
          WRITE(stdout,'(5x,"for a total of",i3," scf calculations")') &
             6*ngeo_strain 
       ENDIF
    ENDIF
 ENDIF

 IF (what=='mur_lc_t'.AND.(ltherm_dos.OR.ltherm_freq)) THEN
    IF (ltherm_glob) THEN
       WRITE(stdout,'(/,5x,"Anharmonic properties from global interpolation &
                            &of the")') 
       WRITE(stdout,'(5x,"total free energy")') 
    ELSE
       WRITE(stdout,'(/,5x,"Anharmonic properties from separate interpolation &
                            &of the")')
       WRITE(stdout,'(5x,"total energy and of the free energy")') 
    ENDIF
 ENDIF
!
!  B factor matrix
! 
IF (( what=='scf_disp'.OR.what=='scf_ph'.OR.what=='mur_lc_ph'.OR. &
      what=='mur_lc_disp'.OR.what=='mur_lc_t'.OR.what=='plot_bz') &
      .AND.with_eigen ) THEN

   WRITE(stdout,'(/,5x,"B factor matrix:")')
   
   CALL b_factor_reduction(b_rest)
   
   DO na=1,nat
    
      WRITE(stdout,'(/,5x,"Site n.",i5,", atom",a6,":")') na, &
                                                          TRIM(atm(ityp(na)))
      !WRITE(stdout,'(/,5x,"case n.",i3)') b_rest(na)

      CALL print_b_fact_shape(b_rest(na))

   ENDDO
ENDIF

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
  CALL compute_density(omega,density,.TRUE.)
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

!           IF (lasymptote.AND.ionode) &
!              CALL EXECUTE_COMMAND_LINE(TRIM(asymptote_command)//' '&
!                                       //TRIM(asy_filename), WAIT=.FALSE.)
        ENDIF
!
!  Form factors of the atoms

        WRITE(stdout,'(/,5x,"Computing the X-ray powder diffraction &
                                                              &intensities")')
        IF (lformf) THEN
           DO it=1,nsp
              IF (lcm) THEN
                WRITE(stdout,'(/,5x,"Form factors from Cromer-Mann &
                             &parameters")')
              ELSE
                 WRITE(stdout,'(/,5x,"Form factors from Doyle-Turner &
                             &or Smith-Burge parameters")')
              ENDIF
              filename=TRIM(flformf)//'.'//TRIM(int_to_char(it))
              CALL write_form_factor(atm(it),smin,smax,nspoint,lcm,filename)
           ENDDO
           CALL plot_form_factors()
        ENDIF

        CALL compute_xrdp(at,bg,celldm(1),nat,tau,nsp,ityp,atm, &
                               lambda,ibrav,lcm,flxrdp)
        CALL plot_xrdp('')
!
!   write the xsf file for plotting the structure using xcrysden
!
        IF (ionode) THEN
           iuout=find_free_unit()
           xsf_filename=TRIM(prefix)//'.xsf'
           OPEN(UNIT=iuout, FILE=TRIM(xsf_filename), STATUS='unknown', &
                                                             FORM='formatted')
           CALL xsf_struct (celldm(1), at, nat, tau, atm, ityp, iuout)
           CLOSE(iuout)
        ENDIF
     ENDIF
     CALL summarize_kpt(xqaux, wqaux, nqaux, letter_path)

  ELSEIF (what=='scf_2d_bands') THEN
     !
     asy_filename=TRIM(flasy)//'.asy'
     IF ( my_image_id==root_image ) THEN

        CALL plot_2d_bz(ibrav, celldm, at, bg, xqaux, wqaux, nqaux, &
             letter, letter_path, npk_label, label_list, asy_filename)
                                                         

        IF (lasymptote.AND.ionode) &
           ierr=system(TRIM(asymptote_command)//' '//TRIM(asy_filename))

!        IF (lasymptote.AND.ionode) &
!           CALL EXECUTE_COMMAND_LINE(TRIM(asymptote_command)//' '&
!                                       //TRIM(asy_filename), WAIT=.FALSE.)
     ENDIF
     !
  ENDIF

  RETURN
END SUBROUTINE thermo_summary

!-----------------------------------------------------------------------
FUNCTION check_group_ibrav(code_group, ibrav)
!-----------------------------------------------------------------------
!
!  This routine checks if the ibrav is compatible with the point group,
!  and if it is not it writes a brief message saying which lattices are
!  compatible
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
USE lattices, ONLY : lattice_name, is_compatible_group_ibrav
USE noncollin_module, ONLY : noncolin, domag
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, ibrav

INTEGER :: is_compatible(32,6), ncomp(32)
LOGICAL :: check_group_ibrav
INTEGER :: i
CHARACTER(LEN=40) :: latt_name
CHARACTER(LEN=11) :: group_name

check_group_ibrav=is_compatible_group_ibrav(code_group, ibrav, &
                                                    is_compatible, ncomp )

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
   WRITE(stdout,'(/,5x,"understand why the symmetries are wrong before &
                      &continuing")') 
   WRITE(stdout,'(5x,"The point group or the Laue class are not used to &
                      &reduce the number of ")')
   WRITE(stdout,'(5x,"computed tensor components")') 
100 CONTINUE
ENDIF
IF (noncolin.AND.domag) check_group_ibrav=.FALSE.

RETURN
END FUNCTION check_group_ibrav

!-----------------------------------------------------------------------
SUBROUTINE find_fft_fact()
!-----------------------------------------------------------------------
!
!  This routine finds the fft_fact that corresponds to the space group
!  of the solid. It assumes that the tau and celldm are already known
!  It sets also the dimensions of the fft that correspond to the cut-off.
!
USE kinds,            ONLY : DP
USE fft_base,         ONLY : dfftp
USE cell_base,        ONLY : ibrav, at, bg
USE thermo_sym,       ONLY : fft_fact, ibrav_group_consistent
USE rap_point_group,  ONLY : code_group
USE space_groups,     ONLY : find_space_group, set_fft_fact
USE symm_base,        ONLY : nsym, sr, ft

IMPLICIT NONE
INTEGER :: sg_number
INTEGER :: unique, trig, isym, aux_sg
LOGICAL :: check_group_ibrav
REAL(DP) :: s01(3), s02(3), ft_(3,48)
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
   ft_(:,1:nsym)=-ft(:,1:nsym)
   CALL find_space_group(ibrav, nsym, sr, ft_, at, bg, sg_number, aux_sg,&
                                  s01, s02, .FALSE.)

   IF (sg_number > 0) THEN
      unique=0
      trig=0
      IF (ibrav==-12.OR.ibrav==-13) unique=1
      IF (ibrav==5) trig=1
      CALL set_fft_fact(sg_number, unique, trig, fft_fact)
      CALL add_origin_fact(s01,fft_fact)
      CALL clean_dfft()
      CALL find_symmetry(fft_fact)
   ENDIF
ELSE
   CALL clean_dfft()
   CALL find_symmetry(fft_fact)
ENDIF

RETURN
END SUBROUTINE find_fft_fact

!-----------------------------------------------------------------------
SUBROUTINE summarize_kpt(xqaux, wqaux, nqaux, letter_path )
!-----------------------------------------------------------------------
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
  
INTEGER, INTENT(IN)  :: nqaux
REAL(DP), INTENT(IN) :: xqaux(3,nqaux)
INTEGER, INTENT(IN)  :: wqaux(nqaux)
CHARACTER(LEN=3), INTENT(IN) :: letter_path(nqaux)

INTEGER :: ik
  
WRITE(stdout,'(/,5x, "k points coordinates (2 pi / alat)")') 
DO ik=1, nqaux
   WRITE(stdout, '(a3, 3f15.8,i5)') letter_path(ik), xqaux(:,ik), &
                                    wqaux(ik)
ENDDO
WRITE(stdout,*)
WRITE(stdout,'(5x, "Input path: ")') 
WRITE(stdout,'(i5)') nqaux
DO ik=1, nqaux
   IF (letter_path(ik) == '   ') THEN
      WRITE(stdout, '(3x, 3f15.8,i5)')  xqaux(:,ik), wqaux(ik)
   ELSE
      WRITE(stdout, '(a3, i8)') letter_path(ik), wqaux(ik)
   ENDIF
ENDDO
WRITE(stdout,*)

RETURN
END SUBROUTINE summarize_kpt

!-----------------------------------------------------------------------
SUBROUTINE add_origin_fact(s0,fft_fact)
!-----------------------------------------------------------------------
!
!   This routine receives the origin shift in crystal coordinates
!   and checks if the fft_factors contains the factors needed to represent
!   the origin shift. If not it adds them to the factors.
!
USE kinds, ONLY : DP
IMPLICIT NONE

REAL(DP), INTENT(IN)   :: s0(3)
INTEGER, INTENT(INOUT) :: fft_fact(3)

INTEGER :: ipol, i
!
!  try only the factors 2, 3, 4, 5, 6, 7, 8
!
DO ipol=1,3
   IF (ABS(s0(ipol)) > 1.D-8) THEN
factors:   DO i=2,8
         IF (fft_fact(ipol)==0 .OR. MOD(fft_fact(ipol),i)/=0) THEN
            IF (ABS(s0(ipol)*i-NINT(s0(ipol)*i))< 1.D-6) THEN
                fft_fact(ipol)= MAX(fft_fact(ipol),1)*i
                EXIT factors
            ENDIF
         ENDIF
      ENDDO  factors
   ENDIF
ENDDO

RETURN
END SUBROUTINE add_origin_fact
