!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_thermo_work_todo(iwork, part, iq_point, irr_value, igeo)
  !-----------------------------------------------------------------------
  !
  !  This routine receives from the asyncronous driver the work to do in
  !  the integer iwork and sets the variables dependent from iwork before 
  !  performing the actual calculation.
  !
  USE kinds,            ONLY : DP
  USE thermo_mod,       ONLY : what, alat_geo
  USE control_thermo,   ONLY : outdir_thermo
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units, outdir
  USE input_parameters, ONLY : k_points, xk, wk, nk1, nk2, nk3,  &
                               k1, k2, k3, nkstot
  USE control_conv, ONLY : ke, keden, nk_test, sigma_test
  USE control_elastic_constants, ONLY : at_save, tau_save
  USE elastic_constants, ONLY : epsilon_geo, apply_strain
  USE control_flags, ONLY : gamma_only, tstress, tprnfor, lbfgs, nstep
  USE force_mod, ONLY : lforce
  USE relax,            ONLY : epse, epsf
  USE input_parameters, ONLY : calculation, etot_conv_thr, forc_conv_thr
  USE io_files,    ONLY : tmp_dir, wfc_dir, prefix, seqopn
  USE io_global,   ONLY : ionode
  USE cell_base,   ONLY : cell_base_init, at
  USE ions_base,   ONLY : tau, nat
  USE fft_base,    ONLY : dfftp, dffts
  USE wvfct,       ONLY : ecutwfc
  USE start_k,     ONLY : init_start_k
  USE klist,       ONLY : degauss
  USE gvect,       ONLY : ecutrho
  USE gvecs,       ONLY : dual
  USE grid_irr_iq, ONLY : irr_iq, comp_irr_iq
  USE disp,        ONLY : nqs, comp_iq
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
        CASE ('scf', 'scf_bands', 'scf_ph', 'scf_disp')
        CASE ('scf_ke')
           ecutwfc = ke(iwork)
           ecutrho = keden(iwork)
           dual = ecutrho / ecutwfc
           CALL clean_dfft()
           outdir=TRIM(outdir_thermo)//'ke'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
        CASE ('scf_nk')
           degauss = sigma_test(iwork)
           nk1=nk_test(iwork)
           nk2=nk1
           nk3=nk1
           CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, &
                               nkstot, xk, wk )
           CALL clean_dfft()
           gamma_only = ( k_points == 'gamma' )
           outdir=TRIM(outdir_thermo)//'ke'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
        CASE ('mur_lc', 'mur_lc_bands', 'mur_lc_ph', 'mur_lc_disp', &
              'mur_lc_t')
           celldm(1)=alat_geo(iwork)
           CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                         trd_ht, rd_ht, cell_units )
           CALL clean_dfft()
           outdir=TRIM(outdir_thermo)//'g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
        CASE ('elastic_constants', 'fi_elastic_constants')
           ibrav=0
           CALL clean_dfft()
           tstress=.TRUE.
           tprnfor=.TRUE.
           DO i=1, 3
              CALL apply_strain(at_save(1,i), at(1,i), epsilon_geo(1,1,iwork))
           ENDDO
           DO ia=1,nat
              write(6,*) 'straing atom', ia
              CALL apply_strain(tau_save(1,ia), tau(1,ia), epsilon_geo(1,1,iwork))
           ENDDO
           IF (what=='elastic_constants') THEN
              calculation='relax'
              lforce=.TRUE.
              lbfgs = .TRUE.
              nstep = 10
              epse = etot_conv_thr
              epsf = forc_conv_thr
           ELSE
              calculation='scf'
              lbfgs=.FALSE.
           ENDIF
           rd_ht = TRANSPOSE( at ) 
           trd_ht=.TRUE.
           cell_units='alat'
           CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                         trd_ht, rd_ht, cell_units )
           outdir=TRIM(outdir_thermo)//'g'//TRIM(int_to_char(iwork))//'/'
           tmp_dir = TRIM ( outdir )
           wfc_dir = tmp_dir
           CALL check_tempdir ( tmp_dir, exst, parallelfs )
           IF ( what=='elastic_constants' .AND. ionode) THEN
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
           CALL errore('set_thermo_work','unknown what',1)
     END SELECT
  ELSE IF (part==2) THEN
     SELECT CASE (TRIM(what))
        CASE ('scf_ph', 'scf_disp','mur_lc_ph','mur_lc_disp','mur_lc_t')
           comp_irr_iq=.FALSE.
           comp_iq=.FALSE.
           jwork=0
           DO iq=1,nqs
              DO irr=0, irr_iq(iq)
                 jwork=jwork+1
                 IF (jwork==iwork) THEN
                    comp_irr_iq(irr,iq)=.TRUE.
                    comp_iq(iq)=.TRUE.
                    iq_point=iq
                    irr_value=irr
                 ENDIF
              ENDDO
           ENDDO
     END SELECT
  ELSE
     CALL errore('set_thermo_work','unknown part',1)
  END IF
  !
  RETURN
  !
END SUBROUTINE set_thermo_work_todo
