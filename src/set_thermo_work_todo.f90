!
! Copyright (C) 2013-2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_thermo_work_todo(iwork, part, iq_point, irr_value)
!-----------------------------------------------------------------------
!
!  This routine receives from the asynchronous driver the number of the
!  work to do, iwork, and sets the input variables of pwscf or of the phonon 
!  according to iwork.
!  When this iwork and part correspond to a phonon calculation it gives 
!  the q point iq_point and the irreducible representation irr_value
!  on output if not using ph.x images. Otherwise the iq_point
!  gives the number of the image that would have done that work in ph.x 
!  and irr_value is not used.
!
  USE kinds,            ONLY : DP
  USE thermo_mod,       ONLY : what, ibrav_geo, celldm_geo, ef_geo
  USE control_thermo,   ONLY : outdir_thermo
  USE control_elastic_constants, ONLY : frozen_ions, use_free_energy
  USE control_bands,    ONLY : nbnd_bands
  USE control_conv,     ONLY : ke, keden, nk_test, sigma_test
  USE control_eldos,    ONLY : lel_free_energy
  USE initial_conf,     ONLY : ibrav_save, tau_save_crys
  USE equilibrium_conf, ONLY : at0, tau0
!
!  the library modules
!
  USE elastic_constants, ONLY : epsilon_geo
  USE strain_mod,        ONLY : apply_strain, print_strain
!
!  the pw variables that are set here or used to set the input
!
  USE input_parameters, ONLY : electron_maxstep, k_points, xk, wk, k1, k2, &
                               k3, nkstot
  USE control_flags,    ONLY : niter, lbands, tstress
  USE cell_base,        ONLY : cell_base_init, at
  USE ions_base,        ONLY : tau, nat
  USE gvecw,            ONLY : ecutwfc
  USE gvect,            ONLY : ecutrho
  USE gvecs,            ONLY : dual
  USE wvfct,            ONLY : nbnd
  USE ener,             ONLY : ef
  USE start_k,          ONLY : init_start_k
  USE klist,            ONLY : degauss

  USE io_global,        ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iwork, part
  INTEGER, INTENT(OUT) :: iq_point, irr_value

  INTEGER :: i, ia, nk1, nk2, nk3, ibrav, igeom, image
  REAL(DP) :: rd_ht(3,3), zero, celldm(6)
  CHARACTER(LEN=10) :: cell_units
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: outdir
  LOGICAL :: trd_ht
  !
  iq_point=0
  irr_value=0
  zero=0.0_DP

  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
!
!  Here we set the kinetic energy 
!
        CASE ('scf_ke')
           ecutwfc = ke(iwork)
           ecutrho = keden(iwork)
           dual = ecutrho / ecutwfc
!
!  recompute the FFT dimension
!
           CALL set_fft_mesh()
!
!  set the tmp_dir for this work
!
           outdir=TRIM(outdir_thermo)//'/ke'//TRIM(int_to_char(iwork))//'/'
           CALL set_tmp_dir( outdir)

        CASE ('scf_nk')
!
!  Here we set the number of k points and the degauss in metals
!
           IF (TRIM(k_points) /='automatic') &
              CALL errore('set_thermo_work_todo', &
                               'kpoint test requires automatic k point',1)
           degauss = sigma_test(iwork)
           nk1=nk_test(1,iwork)
           nk2=nk_test(2,iwork)
           nk3=nk_test(3,iwork)
!
!   for the shift we use the same parameters read in input
!
           CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, &
                               nkstot, xk, wk )
!
!  recompute the FFT dimension
!
           CALL set_fft_mesh()
!
!  set the tmp_dir for this work
!
           outdir=TRIM(outdir_thermo)//'/nk'//TRIM(int_to_char(iwork))//'/'
           CALL set_tmp_dir( outdir)
!
!   then all the cases that require many energies calculations at
!   different geometries. The geometry is set here in the pwscf variables
!
        CASE ('mur_lc',                      &
              'mur_lc_bands',                &
              'mur_lc_dos',                  &
              'mur_lc_ph',                   &
              'mur_lc_disp',                 &
              'mur_lc_t',                    &
              'mur_lc_elastic_constants',    &
              'mur_lc_piezoelectric_tensor', &
              'mur_lc_polarization')
!
!  Initialize the QE variables for the ionic relaxation. 
!
           CALL set_work_for_relaxation(iwork)
!
!    In the relaxed ion case we compute the stress but do use it
!    so lstress(iwork) is .FALSE. but here we set lstres=.TRUE.
!    To save time or if stress is not available comment this command.
!
           tstress=.NOT.frozen_ions
!
!   now set the celldm
!
           celldm(:)=celldm_geo(:,iwork)
           rd_ht=0.0_DP
           CALL cell_base_init ( ibrav_save, celldm, zero, zero, zero, zero, &
                                     zero, zero, .FALSE., rd_ht, ' ' )
!
!   recompute the fft mesh 
!
           CALL set_fft_mesh()
!
! strain uniformly the coordinates to the new celldm
!
           tau=tau_save_crys
           CALL cryst_to_cart( nat, tau, at, 1 )
!
! set the tmp_dir for this geometry
!
           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           CALL set_tmp_dir( outdir)
           IF (.NOT.frozen_ions) CALL clean_bfgs_history()
!
!    the case of quasi-harmonic elastic constants
!
        CASE ('elastic_constants_t')
           niter=electron_maxstep
           CALL set_work_for_relaxation(iwork)
           CALL set_work_for_elastic_const(iwork)
        CASE DEFAULT
           CALL errore('set_thermo_work_todo','unknown what',1)
     END SELECT
  ELSE IF (part==2) THEN
     SELECT CASE (TRIM(what))
!
!   here we set the representation and the q point to calculate in the
!   phonon calculation for the present geometry
!
        CASE ('scf_ph',         &
              'scf_disp',       &
              'mur_lc_ph',      &
              'mur_lc_disp',    &
              'mur_lc_t')       
           CALL set_work_for_ph(iwork, igeom, iq_point, irr_value)
        CASE('elastic_constants_t')
           IF (use_free_energy) THEN
              CALL set_work_for_ph(iwork, igeom, iq_point, irr_value)
           ELSEIF (lel_free_energy) THEN
              CALL set_work_for_bands(iwork)
           ENDIF
        CASE ('mur_lc')

           CALL set_work_for_bands(iwork)
!
!    Here the elastic constant calculation
!
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants')
           niter=electron_maxstep
           CALL set_work_for_relaxation(iwork)
           CALL set_work_for_elastic_const(iwork)
        CASE ('scf_piezoelectric_tensor', 'mur_lc_piezoelectric_tensor')
           niter=electron_maxstep
           ibrav=0
           DO i=1, 3
              CALL apply_strain(at0(1,i), at(1,i), epsilon_geo(1,1,iwork))
           ENDDO
           DO ia=1,nat
              CALL apply_strain(tau0(1,ia), tau(1,ia), epsilon_geo(1,1,iwork))
           ENDDO
           CALL print_strain(epsilon_geo(:,:,iwork))

           CALL set_work_for_relaxation(iwork)

           rd_ht = TRANSPOSE( at ) 
           trd_ht=.TRUE.
           cell_units='alat'
           CALL cell_base_init ( ibrav, celldm, zero, zero, zero, zero, &
                         zero, zero, trd_ht, rd_ht, cell_units )
           CALL set_fft_mesh()
           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           CALL set_tmp_dir( outdir )
           IF (.NOT.frozen_ions) CALL clean_bfgs_history()

        CASE ('scf_polarization','mur_lc_polarization')
     END SELECT
  ELSE
     CALL errore('set_thermo_work_todo','unknown part',1)
  END IF
  !
  RETURN
  !
END SUBROUTINE set_thermo_work_todo

!-----------------------------------------------------------------------
SUBROUTINE set_tmp_dir (outdir)
!-----------------------------------------------------------------------
!
!  This subroutine sets the tmp_dir and wfc_dir to the input outdir
!  and creates the directory if it does not exist. It stops
!  if the directory cannot be created or if some processor within
!  the same image cannot read or write on it.
!
USE io_files,    ONLY : tmp_dir, wfc_dir, check_tempdir

IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: outdir
LOGICAL :: exst, parallelfs

tmp_dir = TRIM ( outdir )
wfc_dir = tmp_dir
CALL check_tempdir ( tmp_dir, exst, parallelfs )
IF (.NOT. parallelfs) &
         CALL errore('set_tmp_dir','some node cannot read or write',1)
RETURN
END SUBROUTINE set_tmp_dir

!-----------------------------------------------------------------------
SUBROUTINE set_work_for_relaxation(iwork)
!-----------------------------------------------------------------------
!
!  This routine initializes the variables that usually are initialized
!  in iosys with calculation='relax'. This is needed because there is 
!  also the possibility that we relax the ions when in the pw.x input
!  calculation='scf'. Here we do a minimal initialization of the
!  pw.x variables. Only lbfgs relaxation is possible.
!  With frozen_ion=.TRUE. relaxation is disabled even if it was
!  requested in the input of pw.x.
!
USE control_thermo,   ONLY : lstress
USE control_elastic_constants, ONLY : frozen_ions 
USE input_parameters, ONLY : etot_conv_thr, forc_conv_thr
USE control_flags,    ONLY : lbfgs, nstep, lforce=>tprnfor, tstress
USE relax,            ONLY : epse, epsf

IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork

tstress=lstress(iwork)
IF (frozen_ions) THEN
   lbfgs=.FALSE.
ELSE
   lforce=.TRUE.
   lbfgs = .TRUE.
   IF (nstep==1) nstep = 20
   epse = etot_conv_thr
   epsf = forc_conv_thr
ENDIF
RETURN
END SUBROUTINE set_work_for_relaxation

!-----------------------------------------------------------------------
SUBROUTINE set_work_for_ph(iwork, igeom, iq_point, irr_value)
!-----------------------------------------------------------------------
!
!    This subroutine sets the comp_iq, comp_irr_iq, comp_f, that
!    corresponds to iwork. For doing this it uses the arrays
!    comp_iq_iw, comp_irr_iq_iw, comp_f_iw, prepared in
!    initialize_thermo_work. If we are doing all the geometries
!    together this routine initialize also the current geometry
!    reading the output of pw.x for the current geometry. 
!
USE control_thermo,   ONLY : all_geometries_together, geometry, irrw, &
                             iqw, comp_irr_iq_iw, comp_iq_iw, comp_f_iw, &
                             done_irr_iq_iw, done_iq_iw
USE ions_base,        ONLY : nat
USE grid_irr_iq,      ONLY : irr_iq, comp_irr_iq, done_irr_iq
USE disp,             ONLY : nqs, comp_iq, done_iq
USE control_ph,       ONLY : epsil, trans
USE control_thermo,   ONLY : lphonon
USE freq_ph,          ONLY : fpol
USE images_omega,     ONLY : comp_f
USE mp_images,        ONLY : nimage
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork
INTEGER, INTENT(OUT) :: igeom, iq_point, irr_value

INTEGER :: image
!
!  With recover this iwork could be already done and we return.
!
IF (.NOT.lphonon(iwork)) RETURN
!
igeom=1
iq_point=1
irr_value=0
!
!  If we are doing all geometries together at each task we must
!  read the output of pw.x. Get here the information.
!
IF (all_geometries_together) THEN
   igeom=geometry(iwork)
   CALL prepare_do_phonon(igeom)
ENDIF
!
!  set the default. Nothing is calculated if not explicitely set below
!
comp_irr_iq=.FALSE.
comp_iq=.FALSE.
done_irr_iq=.FALSE.
done_iq=.FALSE.
comp_f=.FALSE.
!
!  here sets the variables that controls the phonon run
!
IF (trans) THEN
   comp_irr_iq(:,:)=comp_irr_iq_iw(:,1:nqs,iwork)
   comp_iq(:)=comp_iq_iw(1:nqs,iwork)
   done_irr_iq(:,:)=done_irr_iq_iw(:,1:nqs,iwork)
   done_iq(:)=done_iq_iw(1:nqs,iwork)
   iq_point=iqw(iwork)
   irr_value=irrw(iwork)
ELSEIF (fpol) THEN
!
!   This is the frequency dependent case. There is only one q point 
!   and the 0 irrep to be done. We do all the frequencies set
!   into comp_f_iw for iwork
!
   comp_iq(1)=.TRUE.
   comp_irr_iq(0,1)=.TRUE.
   comp_f(:)=comp_f_iw(:,iwork)
ELSEIF (epsil) THEN
!
!  In this case we make only the electric field perturbation. There
!  is only one q point and the 0 irrep.
!
   comp_iq(1)=.TRUE.
   comp_irr_iq(0,1)=.TRUE.
ENDIF

RETURN
END SUBROUTINE set_work_for_ph

!-----------------------------------------------------------------------
SUBROUTINE set_work_for_bands(iwork)
!-----------------------------------------------------------------------
USE kinds,          ONLY : DP
USE ions_base,      ONLY : tau, nat
USE cell_base,      ONLY : cell_base_init, at
USE control_flags,  ONLY : lbands
USE control_bands,  ONLY : nbnd_bands
USE wvfct,          ONLY : nbnd
USE thermo_mod,     ONLY : what, ibrav_geo, celldm_geo, ef_geo
USE control_thermo, ONLY : outdir_thermo
USE klist,          ONLY : lgauss, ltetra
USE ener,           ONLY : ef
USE initial_conf,   ONLY : tau_save_crys

IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork
REAL(DP) :: rd_ht(3,3), zero, celldm(6)
CHARACTER(LEN=256) :: outdir
CHARACTER(LEN=6) :: int_to_char
!
!   now set the celldm
!
celldm(:)=celldm_geo(:,iwork)
rd_ht=0.0_DP
zero=0.0_DP
CALL cell_base_init ( ibrav_geo(iwork), celldm, zero, zero, &
                      zero, zero, zero, zero, .FALSE., rd_ht, ' ' )
CALL set_dos_kpoints()
lbands=.FALSE.
IF (lgauss.OR.ltetra) ef=ef_geo(iwork)
!
!    use nbnd_bands to control how many bands to compute. If it is zero
!    we use the number of bands of the self-consistent calculation
!
IF (nbnd_bands > nbnd) nbnd = nbnd_bands
!
!   recompute the fft mesh 
!
CALL set_fft_mesh()
!
! strain uniformly the coordinates to the new celldm
!
tau=tau_save_crys
CALL cryst_to_cart( nat, tau, at, 1 )
!
! set the tmp_dir for this geometry
!
outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
CALL set_tmp_dir(outdir)
RETURN
END SUBROUTINE set_work_for_bands
