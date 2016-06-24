!
! Copyright (C) 2013-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_thermo_work_todo(iwork, part, iq_point, irr_value, igeo)
  !-----------------------------------------------------------------------
  !
  !  This routine receives from the asynchronous driver the work to do in
  !  the integer iwork and sets the input variables of pwscf or of the
  !  phonon according to iwork before performing the actual calculation.
  !
  USE kinds,            ONLY : DP
  USE thermo_mod,       ONLY : what, celldm_geo
  USE control_thermo,   ONLY : outdir_thermo, lstress
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units, outdir, &
                               electron_maxstep, &
                               k_points, xk, wk, nk1, nk2, nk3,  &
                               k1, k2, k3, nkstot, &
                               calculation, etot_conv_thr, forc_conv_thr
  USE control_conv, ONLY : ke, keden, nk_test, sigma_test
  USE control_elastic_constants, ONLY : at_save, tau_save, frozen_ions
  USE elastic_constants, ONLY : epsilon_geo, apply_strain, print_strain
  USE control_ph, ONLY : recover
  USE control_flags, ONLY : gamma_only, tstress, tprnfor, lbfgs, nstep, niter, &
                            lscf, lmd
  USE cellmd,        ONLY : lmovecell
  USE force_mod, ONLY : lforce, lstres
  USE relax,       ONLY : epse, epsf
  USE io_files,    ONLY : tmp_dir, wfc_dir, prefix, seqopn
  USE io_global,   ONLY : ionode
  USE cell_base,   ONLY : cell_base_init, at
  USE ions_base,   ONLY : tau, nat
  USE fft_base,    ONLY : dfftp, dffts
  USE gvecw,       ONLY : ecutwfc
  USE start_k,     ONLY : init_start_k
  USE klist,       ONLY : degauss
  USE gvect,       ONLY : ecutrho
  USE gvecs,       ONLY : dual
  USE grid_irr_iq, ONLY : irr_iq, comp_irr_iq, done_irr_iq
  USE disp,        ONLY : nqs, comp_iq, done_iq
  USE io_global,   ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: iwork, iq_point, irr_value
  INTEGER, INTENT(IN) :: part, igeo
  INTEGER :: jwork, irr, iq, i, ia, iunupdate
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename
  LOGICAL :: exst, parallelfs
  !
  iq_point=0
  irr_value=0

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
           CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, &
                               nkstot, xk, wk )
           CALL set_fft_mesh()
           gamma_only = ( k_points == 'gamma' )
           outdir=TRIM(outdir_thermo)//'/ke'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
!
!   then all the cases that require a many energies calculations at
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
              calculation='scf'
              lbfgs=.FALSE.
           ELSE
              IF (calculation=='vc-relax') lmovecell=.TRUE.
              lforce=.TRUE.
              lstres=.TRUE.
              lbfgs = .TRUE.
              nstep = 20
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ENDIF
!
!   now set the celldm. Both the input parameters and the variables in 
!   cell_base module
!
           celldm(:)=celldm_geo(:,iwork)
           CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                         trd_ht, rd_ht, cell_units )
!
!
!
           CALL set_fft_mesh()
!
! strain uniformly the coordinates to the new celldm
!
           tau=tau_save
           CALL cryst_to_cart( nat, tau, at, 1 )

           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )

           IF (.NOT.frozen_ions .AND. ionode) THEN
              !
              !  clean the bfgs history
              !
              iunupdate=2
              CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
              CLOSE(iunupdate, STATUS='DELETE')
              filename = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs'
              OPEN( iunupdate, FILE=TRIM(filename), FORM='FORMATTED')
              CLOSE(iunupdate, STATUS='DELETE')
           END IF

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
           comp_irr_iq=.FALSE.
           comp_iq=.FALSE.
           jwork=0
           DO iq=1,nqs
              DO irr=0, irr_iq(iq)
                 jwork=jwork+1
                 IF (jwork==iwork) THEN
                    IF (recover) THEN
                       comp_irr_iq(irr,iq)=.NOT.done_irr_iq(irr,iq)
                       comp_iq(iq)=.NOT.done_iq(iq)
                    ELSE
                       comp_irr_iq(irr,iq)=.TRUE.
                       comp_iq(iq)=.TRUE.
                    ENDIF
                    iq_point=iq
                    irr_value=irr
                 ENDIF
              ENDDO
           ENDDO
!
!    Here the elastic constant calculation
!
        CASE ('scf_elastic_constants', 'mur_lc_elastic_constants', &
                                                      'elastic_constants_t')
           tstress=.TRUE.
           tprnfor=.TRUE.
           niter = electron_maxstep
           IF (frozen_ions) THEN
              calculation='scf'
              lstres=lstress(iwork)
              lbfgs=.FALSE.
           ELSE
              calculation='relax'
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
           tstress=.TRUE.
           tprnfor=.TRUE.
           niter = electron_maxstep
           DO i=1, 3
              CALL apply_strain(at_save(1,i), at(1,i), epsilon_geo(1,1,iwork))
           ENDDO
           DO ia=1,nat
              CALL apply_strain(tau_save(1,ia), tau(1,ia), epsilon_geo(1,1,iwork))
           ENDDO
           WRITE(stdout,'(/,2x,76("-"))')
           CALL print_strain(epsilon_geo(:,:,iwork))
           IF (frozen_ions) THEN
              calculation='scf'
              lstres=.TRUE.
              lbfgs=.FALSE.
           ELSE
              calculation='relax'
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
           CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                         trd_ht, rd_ht, cell_units )
           CALL set_fft_mesh()
           outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
           IF (.NOT.frozen_ions .AND. ionode) THEN
              !
              !  clean the bfgs history
              !
              iunupdate=2
              CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
              CLOSE(iunupdate, STATUS='DELETE')
              filename = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs'
              OPEN( iunupdate, FILE=TRIM(filename), FORM='FORMATTED')
              CLOSE(iunupdate, STATUS='DELETE')
            END IF
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
!
SUBROUTINE set_work_for_elastic_const(iwork)
   USE kinds,      ONLY : DP
   USE thermo_mod, ONLY : celldm_geo, ibrav_geo
   USE cell_base,  ONLY : cell_base_init, at
   USE ions_base,  ONLY : tau, nat
   USE control_thermo,   ONLY : outdir_thermo
   USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                                trd_ht, rd_ht, cell_units, outdir
   USE control_elastic_constants, ONLY : at_save, tau_save, frozen_ions, &
                                elastic_algorithm, rot_mat, aap_mat, apa_mat
   USE elastic_constants, ONLY : epsilon_geo, apply_strain, print_strain
   USE control_mur,   ONLY : celldm0
   USE rotate,        ONLY : rotate_vect
   USE io_files,      ONLY : tmp_dir, wfc_dir, prefix, seqopn
   USE io_global,     ONLY : stdout, ionode

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: iwork

   INTEGER :: i, na, ipol, jpol, kpol, iunupdate
   REAL(DP), ALLOCATABLE :: tau_ocoord(:,:)
   REAL(DP) :: atp(3,3)
   LOGICAL :: exst, parallelfs
   CHARACTER(LEN=256) :: filename
   CHARACTER(LEN=6) :: int_to_char

   WRITE(stdout,'(/,2x,76("-"))')
   CALL print_strain(epsilon_geo(:,:,iwork))
!
!  entering here we have:
!  
!  at_save that contains the unstrained at vectors in units of celldm0(1)
!
!  tau_save that contains the crystal coordinates of the atoms. In 
!  a uniform strain these coordinates do not change.
!
!
!    first bring tau in the strained lattice
!
   IF (elastic_algorithm=='advanced' .OR. &
                              elastic_algorithm=='energy') THEN
      atp=0.0_DP
      DO ipol=1,3
         DO jpol=1,3
            DO kpol=1,3
            atp(ipol,jpol) = atp(ipol,jpol) + apa_mat(jpol,kpol,iwork)*&
                                         at_save(ipol,kpol)
            END DO
         END DO
      ENDDO
   ELSE
      atp(:,:)=at_save(:,:)
   ENDIF

   DO i=1, 3
      CALL apply_strain(atp(1,i), at(1,i), epsilon_geo(1,1,iwork))
   ENDDO
!
!  tau save are in crystal coordinates. A uniform strain of these coordinates
!  means to keep them constant. We just rotate them in case the direct
!  lattice vectors have changed
!
   IF (elastic_algorithm=='advanced' .OR. &
                              elastic_algorithm=='energy') THEN
      tau=0.0_DP
      DO na=1,nat
         DO ipol=1,3
            DO jpol=1,3
               tau(ipol,na) = tau(ipol,na) + aap_mat(jpol,ipol,iwork)*&
                                         tau_save(jpol,na)
            END DO
         END DO
      END DO
   ELSE
      tau=tau_save
   END IF
!
!  bring tau to cartesian coordinates
!
   CALL cryst_to_cart( nat, tau, at, 1 )

   IF (elastic_algorithm=='standard') THEN
      ibrav=0
      rd_ht = TRANSPOSE( at )
      trd_ht=.TRUE.
      cell_units='alat'
      CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, &
                     cosac, cosbc, trd_ht, rd_ht, cell_units )
      CALL set_fft_mesh()
   ELSEIF (elastic_algorithm=='advanced' .OR. &
                                           elastic_algorithm=='energy') THEN
      ibrav = ibrav_geo(iwork)
      celldm(:)=celldm_geo(:,iwork)
      cell_units='alat'
      trd_ht=.FALSE.
      rd_ht=0.0_DP
      CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, &
                         cosac, cosbc, trd_ht, rd_ht, cell_units )
      ALLOCATE(tau_ocoord(3,nat))
      tau_ocoord=tau
!
!   In this scheme sometimes the cartesian axes of the strained 
!   and unstrained cell are different. We rotate all the atomic positions
!   already strained to the new axis.
!
      CALL rotate_vect(rot_mat(1,1,iwork), nat, tau_ocoord, tau, 1)
      DEALLOCATE(tau_ocoord)
!
!  bring the tau in the correct units of the new alat
!
      tau=tau * celldm0(1) / celldm(1)
!
!  find the optimal fft mesh
!
      CALL find_fft_fact()
  ENDIF

  outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
  tmp_dir = TRIM ( outdir )
  wfc_dir = tmp_dir
  CALL check_tempdir ( tmp_dir, exst, parallelfs )
  IF (.NOT.frozen_ions .AND. ionode) THEN
     !
     !  clean the bfgs history
     !
     iunupdate=2
     CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
     CLOSE(iunupdate, STATUS='DELETE')
     filename = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs'
     OPEN(iunupdate, FILE=TRIM(filename), FORM='FORMATTED')
     CLOSE(iunupdate, STATUS='DELETE')
  END IF
  RETURN
  END SUBROUTINE set_work_for_elastic_const
