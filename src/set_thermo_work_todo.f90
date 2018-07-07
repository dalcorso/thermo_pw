!
! Copyright (C) 2013-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_thermo_work_todo(iwork, part, iq_point, irr_value, auxdyn_loc)
  !-----------------------------------------------------------------------
!
!  This routine receives from the asynchronous driver the work to do in
!  the integer iwork and sets the input variables of pwscf or of the
!  phonon according to iwork before performing the actual calculation.
!  On output it gives the q point and the irreducible representation, if
!  this correspond to a phonon calculation.
!
  USE kinds,            ONLY : DP
  USE thermo_mod,       ONLY : what, celldm_geo, tot_ngeo, start_geometry, &
                               last_geometry
  USE control_thermo,   ONLY : outdir_thermo, lstress, lphonon, &
                               all_geometries_together
  USE control_elastic_constants, ONLY : frozen_ions
  USE control_conv,     ONLY : ke, keden, nk_test, sigma_test
  USE initial_conf,     ONLY : ibrav_save, tau_save_crys, collect_info_save, &
                               geometry
  USE equilibrium_conf, ONLY : at0, tau0
  USE images_omega,     ONLY : omega_group
  USE control_qe,       ONLY : use_ph_images
!
!  the library modules
!
  USE elastic_constants, ONLY : epsilon_geo
  USE strain_mod,        ONLY : apply_strain, print_strain
  USE mp_asyn,           ONLY : with_asyn_images
!
!  the pw variables that are set here or used to set the input
!
  USE input_parameters, ONLY : electron_maxstep, k_points, xk, wk, k1, k2, &
                               k3, nkstot, etot_conv_thr, forc_conv_thr
  USE control_flags,    ONLY : lbfgs, nstep, niter 
  USE cell_base,   ONLY : cell_base_init, at
  USE ions_base,   ONLY : tau, nat
  USE gvecw,       ONLY : ecutwfc
  USE gvect,       ONLY : ecutrho
  USE gvecs,       ONLY : dual
  USE force_mod,   ONLY : lforce, lstres
  USE relax,       ONLY : epse, epsf
  USE start_k,     ONLY : init_start_k
  USE klist,       ONLY : degauss
  USE freq_ph,     ONLY : fpol, nfs
  USE optical,     ONLY : start_freq, last_freq
  USE io_files,    ONLY : tmp_dir, wfc_dir, check_tempdir
!
!   the phonon variables set here or used to set the input
!
  USE grid_irr_iq, ONLY : irr_iq, comp_irr_iq, done_irr_iq
  USE disp,        ONLY : nqs, comp_iq, done_iq
  USE control_ph,       ONLY : recover
  USE images_omega, ONLY : comp_f
  USE mp_images,    ONLY : nimage

  USE io_global,   ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iwork, part
  INTEGER, INTENT(OUT) :: iq_point, irr_value
  CHARACTER(LEN=256), INTENT(INOUT) :: auxdyn_loc

  INTEGER :: jwork, irr, iq, i, j, ia, nk1, nk2, nk3, ibrav, start_omega, &
             igeom, startge, lastge, target_nqs, target_irr_iq, image
  REAL(DP) :: rd_ht(3,3), zero, celldm(6)
  CHARACTER(LEN=10) :: cell_units
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: outdir
  LOGICAL :: exst, parallelfs, trd_ht, something_to_do_all, std
  !
  iq_point=0
  irr_value=0
  zero=0.0_DP

  IF (part == 1) THEN
     SELECT CASE (TRIM(what))
!
!  first the cases where there is nothing to do
!
        CASE ('scf',              &
              'scf_bands',        &
              'scf_dos',          &
              'scf_ph',           &
              'scf_disp',         &
              'elastic_constants_t') 
!
!  then the cases in which we set the kinetic energy and the k points
!
        CASE ('scf_ke')
           ecutwfc = ke(iwork)
           ecutrho = keden(iwork)
           dual = ecutrho / ecutwfc
           CALL set_fft_mesh()
           outdir=TRIM(outdir_thermo)//'/ke'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
        CASE ('scf_nk')
           degauss = sigma_test(iwork)
           nk1=nk_test(1,iwork)
           nk2=nk_test(2,iwork)
           nk3=nk_test(3,iwork)
           IF (TRIM(k_points) /='automatic') &
              CALL errore('set_thermo_work_todo', &
                               'kpoint test requires automatic k point',1)
!
!   for the shift we use the same parameters read in input
!
           CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, &
                               nkstot, xk, wk )
           CALL set_fft_mesh()
           outdir=TRIM(outdir_thermo)//'/nk'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
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

           IF (frozen_ions) THEN
              lbfgs=.FALSE.
           ELSE
              lforce=.TRUE.
              lstres=.TRUE.
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
!
!   now set the celldm. 
!
           celldm(:)=celldm_geo(:,iwork)
           rd_ht=0.0_DP
           CALL cell_base_init ( ibrav_save, celldm, zero, zero, zero, zero, &
                                     zero, zero, .FALSE., rd_ht, ' ' )
           CALL set_fft_mesh()
!
! strain uniformly the coordinates to the new celldm
!
           tau=tau_save_crys
           CALL cryst_to_cart( nat, tau, at, 1 )

           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )

           IF (.NOT.frozen_ions) CALL clean_bfgs_history()

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
           IF (all_geometries_together) THEN
              igeom=geometry(iwork)
              std=something_to_do_all(iwork, igeom, iq_point, irr_value)
              IF (std) THEN
                 CALL initialize_ph_geometry(igeom, auxdyn_loc)
              ELSE
                 lphonon(iwork)=.FALSE.
                 RETURN
              ENDIF
           ELSE
              igeom=1
           ENDIF
           comp_irr_iq=.FALSE.
           comp_iq=.FALSE.
           comp_f=.FALSE.
           jwork=0
           IF (fpol) THEN
              comp_iq(1)=.TRUE.
              comp_irr_iq(0,1)=.TRUE.
              IF (with_asyn_images) THEN
                 start_omega=(iwork-1)*omega_group
                 DO i=1, omega_group
                    j=MIN(start_omega+i, nfs)
                    comp_f(j)=.TRUE.
                 ENDDO
              ELSE
                 comp_f=.TRUE.
              ENDIF
              lphonon(iwork)=.FALSE.
              DO i=start_freq,last_freq
                 lphonon(iwork)=lphonon(iwork).OR.comp_f(i)
              ENDDO
           ELSE
              IF (use_ph_images) THEN
                 image=MOD(iwork-1,nimage)+1
                 DO iq=1,collect_info_save(igeom)%nqs
                    DO irr=0, collect_info_save(igeom)%irr_iq(iq)
                       IF (collect_info_save(igeom)%&
                                          comp_irr_iq(irr,iq,image)==1) &
                          comp_irr_iq(irr,iq)=.TRUE.
                       IF (collect_info_save(igeom)%&
                                           done_irr_iq(irr,iq,image)==1) &
                          comp_irr_iq(irr,iq)=.FALSE.
                    ENDDO
                    IF (collect_info_save(igeom)%comp_iq(iq,image)==1) &
                          comp_iq(iq)=.TRUE.
                    IF (collect_info_save(igeom)%done_iq(iq,image)==1) &
                          comp_iq(iq)=.FALSE.
                 ENDDO
                 iq_point=iwork
                 irr_value=0
              ELSE
                 startge=1
                 lastge=1
                 target_nqs=nqs
                 IF (all_geometries_together) THEN
                    startge=start_geometry
                    lastge=last_geometry
                    target_nqs=collect_info_save(igeom)%nqs
                 ENDIF
                 DO igeom=startge,lastge
                    DO iq=1, target_nqs
                       target_irr_iq=irr_iq(iq)
                       IF (all_geometries_together) &
                          target_irr_iq=collect_info_save(igeom)%irr_iq(iq)
                       DO irr=0, target_irr_iq
                          jwork=jwork+1
                          IF (jwork==iwork) THEN
                             comp_irr_iq(irr,iq)=.TRUE.
                             comp_iq(iq)=.TRUE.
                             IF (all_geometries_together) THEN
                                IF (collect_info_save(igeom)%&
                                            done_irr_iq(irr,iq,1)==1) &
                                    comp_irr_iq(irr,iq)=.FALSE.
                                IF (collect_info_save(igeom)%done_iq(iq,1)==1)&
                                   comp_iq(iq)=.FALSE.
                             ELSE
                                IF (done_irr_iq(irr,iq)) &
                                   comp_irr_iq(irr,iq)=.FALSE.
                                IF (done_iq(iq)) comp_iq(iq)=.FALSE.
                             ENDIF
                             iq_point=iq
                             irr_value=irr
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
!
!    Here the elastic constant calculation
!
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants', &
                                                      'elastic_constants_t')
           niter = electron_maxstep
           IF (frozen_ions) THEN
              lstres=lstress(iwork)
              lbfgs=.FALSE.
           ELSE
              lforce=.TRUE.
              lstres=lstress(iwork)
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
           CALL set_work_for_elastic_const(iwork)
        CASE ('scf_piezoelectric_tensor', 'mur_lc_piezoelectric_tensor')
           ibrav=0
           niter = electron_maxstep
           DO i=1, 3
              CALL apply_strain(at0(1,i), at(1,i), epsilon_geo(1,1,iwork))
           ENDDO
           DO ia=1,nat
              CALL apply_strain(tau0(1,ia), tau(1,ia), epsilon_geo(1,1,iwork))
           ENDDO
           WRITE(stdout,'(/,2x,76("-"))')
           CALL print_strain(epsilon_geo(:,:,iwork))
           IF (frozen_ions) THEN
              lstres=.TRUE.
              lbfgs=.FALSE.
           ELSE
              lforce=.TRUE.
              lstres=.TRUE.
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
           rd_ht = TRANSPOSE( at ) 
           trd_ht=.TRUE.
           cell_units='alat'
           CALL cell_base_init ( ibrav, celldm, zero, zero, zero, zero, &
                         zero, zero, trd_ht, rd_ht, cell_units )
           CALL set_fft_mesh()

           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )

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
!
SUBROUTINE set_work_for_elastic_const(iwork)

USE kinds,            ONLY : DP
!
!  variables of thermo_pw
!
USE control_thermo,   ONLY : outdir_thermo
USE equilibrium_conf, ONLY : celldm0, at0, tau0_crys
USE thermo_mod,       ONLY : celldm_geo, ibrav_geo
USE control_elastic_constants, ONLY : frozen_ions, elastic_algorithm, rot_mat
!
!  library routines
!
USE elastic_constants, ONLY : epsilon_geo
USE strain_mod,       ONLY : apply_strain, print_strain
USE rotate,           ONLY : rotate_vect
!
!  variables of qe
!
USE cell_base,        ONLY : cell_base_init, at
USE ions_base,        ONLY : tau, nat
USE io_files,         ONLY : tmp_dir, wfc_dir, check_tempdir
USE io_global,        ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork

INTEGER  :: ivec, na, ipol, jpol, ibrav
REAL(DP), ALLOCATABLE :: tau_ocoord(:,:)
REAL(DP) :: rd_ht(3,3), zero, celldm_(6)
LOGICAL  :: exst, parallelfs, trd_ht
CHARACTER(LEN=256) :: outdir
CHARACTER(LEN=10)  :: cell_units
CHARACTER(LEN=6)   :: int_to_char

WRITE(stdout,'(/,2x,76("-"))')
CALL print_strain(epsilon_geo(:,:,iwork))
!
!  entering here we have:
!  
!  at0 that contains the unstrained at vectors in units of celldm0(1)
!
!  tau0_crys that contains the crystal coordinates of the atoms. In 
!  a uniform strain these coordinates do not change.
!
!  first strain the at0. 
!
DO ivec=1, 3
   CALL apply_strain(at0(1,ivec), at(1,ivec), epsilon_geo(1,1,iwork))
ENDDO
!
!   Now find tau strained in cartesian coordinates using the strained at
!
tau=tau0_crys
CALL cryst_to_cart( nat, tau, at, 1 )
!
zero=0.0_DP
IF (elastic_algorithm=='standard') THEN
   ibrav=0
   rd_ht = TRANSPOSE( at )
   trd_ht=.TRUE.
   cell_units='alat'
   CALL cell_base_init ( ibrav, celldm0, zero, zero, zero, zero, &
                     zero, zero, trd_ht, rd_ht, cell_units )
!
!  the atomic coordinates are strained uniformely and are not modified here
!
ELSEIF (elastic_algorithm=='advanced' .OR. &
                                        elastic_algorithm=='energy') THEN
!
!  compute the at on the basis of ibrav_geo and celldm_geo
!
   ibrav = ibrav_geo(iwork)
   celldm_(:)=celldm_geo(:,iwork)
   trd_ht=.FALSE.
   rd_ht=0.0_DP
   CALL cell_base_init ( ibrav, celldm_, zero, zero, zero, zero, &
                         zero, zero, .FALSE., rd_ht, ' ' )
!
!   In this scheme sometimes the cartesian axes of the strained 
!   and unstrained cells are different. We rotate all the atomic positions
!   already strained to the new axis.
!
   ALLOCATE(tau_ocoord(3,nat))
   tau_ocoord=tau

   CALL rotate_vect(rot_mat(1,1,iwork), nat, tau_ocoord, tau, 1)
   DEALLOCATE(tau_ocoord)
!
!  bring the tau in the correct units of the new alat
!
   tau=tau * celldm0(1) / celldm_(1)
!
!  find the optimal fft mesh
!
   CALL find_fft_fact()
ENDIF
CALL set_fft_mesh()

outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
tmp_dir = TRIM ( outdir )
wfc_dir = tmp_dir
CALL check_tempdir ( tmp_dir, exst, parallelfs )
IF (.NOT.frozen_ions) CALL clean_bfgs_history()

RETURN
END SUBROUTINE set_work_for_elastic_const

SUBROUTINE clean_bfgs_history()

USE io_files,    ONLY : tmp_dir, prefix, seqopn
USE io_global,   ONLY : ionode

IMPLICIT NONE
INTEGER :: iunupdate
INTEGER :: find_free_unit
LOGICAL :: exst
CHARACTER(LEN=256) :: filename

IF (ionode) THEN
   !
   !  clean the bfgs history
   !
   iunupdate=find_free_unit()
   CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
   CLOSE(iunupdate, STATUS='DELETE')
   filename = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs'
   OPEN( iunupdate, FILE=TRIM(filename), FORM='FORMATTED')
   CLOSE(iunupdate, STATUS='DELETE')
ENDIF
RETURN
END SUBROUTINE clean_bfgs_history

LOGICAL FUNCTION something_to_do_all(iwork, igeom, iq_point, irr_point)

USE thermo_mod,   ONLY : tot_ngeo, start_geometry, last_geometry
USE initial_conf, ONLY : collect_info_save
USE control_qe,   ONLY : use_ph_images
USE mp_images,    ONLY : nimage

IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork, igeom
INTEGER, INTENT(OUT) :: iq_point, irr_point

LOGICAL :: std
INTEGER :: iq, irr, jgeom, jwork, image

std=.FALSE.
IF (use_ph_images) THEN
   image=MOD(iwork-1, nimage) + 1 
   DO iq=1, collect_info_save(igeom)%nqs
      DO irr=0, collect_info_save(igeom)%irr_iq(iq)
         IF ((collect_info_save(igeom)%comp_irr_iq(irr,iq,image)==1).AND.&
             (collect_info_save(igeom)%done_irr_iq(irr,iq,image)==0)) std=.TRUE.
      ENDDO
   ENDDO
   irr_point=0
   iq_point=0
ELSE
   jwork=0
   DO jgeom=start_geometry, last_geometry
      DO iq=1, collect_info_save(jgeom)%nqs
         DO irr=0, collect_info_save(jgeom)%irr_iq(iq)
            jwork=jwork+1
            IF (jwork==iwork) THEN
               IF (collect_info_save(jgeom)%done_irr_iq(irr,iq,1)==0) std=.TRUE.
               irr_point=irr
               iq_point=iq
            ENDIF    
         ENDDO
      ENDDO
   ENDDO
ENDIF
something_to_do_all=std

RETURN
END FUNCTION something_to_do_all
