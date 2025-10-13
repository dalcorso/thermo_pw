!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE run_nscf_tpw(do_band, iq)
  !-----------------------------------------------------------------------
  !
  !! This is the main driver of the \(\texttt{pwscf}\) program called from
  !! the \(\texttt{PHonon}\) code.
  !
  !
  USE control_flags,   ONLY : conv_ions, lforce=>tprnfor, tstress
  USE starting_scf,    ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir, seqopn
  USE lsda_mod,        ONLY : nspin
  USE control_flags,   ONLY : restart, lscf, io_level
  USE check_stop,      ONLY : check_stop_now
  USE fft_base,        ONLY : dffts, dfftp
  !!!
  USE fft_types, ONLY: fft_type_allocate
  USE cell_base, ONLY: at, bg, tpiba
  USE gvect,     ONLY: gcutm
  USE gvecs,     ONLY: gcutms
  USE gvect,     ONLY: eigts1, eigts2, eigts3
  USE many_k_mod, ONLY: eigts1_d, eigts2_d, eigts3_d
  USE klist,     ONLY: qnorm, nelec, nks, nkstot
  USE wvfct,     ONLY: nbnd, nbndx, npwx
  !!!
  USE disp,            ONLY : lgamma_iq
  USE control_ph,      ONLY : recover, tmp_dir_phq, &
                              ext_restart, bands_computed, newgrid, qplot, &
                              only_wfc
  USE control_lr,      ONLY : reduce_io
  USE io_global,       ONLY : stdout
  USE uspp,            ONLY : okvan
  USE grid_irr_iq,     ONLY : done_bands
  USE acfdtest,        ONLY : acfdt_is_active, acfdt_num_der, ir_point, delta_vrs
  USE scf,             ONLY : vrs
  USE mp_bands,        ONLY : intra_bgrp_comm, nyfft
  USE mp_pools,        ONLY : kunit
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_barrier

  USE lr_symm_base,    ONLY : minus_q, nsymq, invsymq
  USE control_lr,      ONLY : ethr_nscf
  USE qpoint,          ONLY : xq
  USE noncollin_module,ONLY : noncolin, domag, npol
  USE el_phon,         ONLY : elph_mat
  USE ahc,             ONLY : elph_ahc
  USE rism_module,     ONLY : lrism, rism_set_restart
  USE control_qe,      ONLY : many_k
  USE many_k_mod,      ONLY : deallocate_many_k, allocate_many_k, init_k_blocks
  USE two_chem,        ONLY : twochem
  USE input_parameters, ONLY : occupations
  !
 !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  INTEGER, INTENT(IN) :: iq
  !
  LOGICAL :: exst
  INTEGER :: nsolv
  !
  CALL start_clock( 'PWSCF' )
  !
  io_level=1
  ! FIXME: following section does not belong to this subroutine
  IF (done_bands(iq)) THEN
     WRITE (stdout,'(/,5x,"Bands found: reading from ",a)') TRIM(tmp_dir_phq)
     IF (many_k) CALL deallocate_many_k()
     CALL clean_pw_tpw( .TRUE. )
     CALL close_files(.true.)
     wfc_dir=tmp_dir_phq
     tmp_dir=tmp_dir_phq
     ! FIXME: kunit is set here: in this case we do not go through setup_nscf
     ! FIXME: and read_file calls divide_et_impera that needs kunit
     ! FIXME: qnorm (also set in setup_nscf) is needed by allocate_nlpot
     kunit = 2
     IF ( lgamma_iq(iq) ) kunit = 1
     IF (noncolin.AND.domag) THEN 
        kunit = 4
        IF (lgamma_iq(iq)) kunit=2
     ENDIF
     qnorm = SQRT(xq(1)**2+xq(2)**2+xq(3)**2) * tpiba
     CALL read_file()
     CALL set_small_group_of_q(nsymq,invsymq,minus_q)
     IF (many_k) THEN
        CALL init_k_blocks(npwx,npol,nbndx,nks,nkstot,nbnd,dfftp%nnr,&
                                                      dffts%nnr,okvan)
        nsolv=1
        IF (noncolin.AND.domag) nsolv=2
        CALL allocate_many_k(nsolv)
#if defined(__CUDA)
!
!  this should be taken care by read_file, but for some unknown reason
!  it is presently not updating these variables on the device
!
        eigts1_d=eigts1
        eigts2_d=eigts2
        eigts3_d=eigts3
#endif
     ENDIF
     RETURN
  ENDIF
  !
  IF (many_k) CALL deallocate_many_k()
  CALL clean_pw_tpw( .FALSE. )
  !
  CALL close_files(.true.)
  !
  ! From now on, work only on the _ph virtual directory
  !
  wfc_dir=tmp_dir_phq
  tmp_dir=tmp_dir_phq
  ! ... Setting the values for the nscf run
  !
  startingconfig = 'input'
  starting_pot   = 'file'
  starting_wfc   = 'atomic'
  lscf           = .FALSE.
  lforce         = .FALSE.
  tstress        = .FALSE.
  restart        = ext_restart
  conv_ions      = .true.
  ethr_nscf      = 1.0D-9 / nelec 
  !
  ! threshold for diagonalization ethr_nscf - should be good for all cases
  IF (lrism) CALL rism_set_restart()
  !
  CALL fft_type_allocate ( dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft)
  !
  CALL setup_nscf_tpw ( newgrid, xq, elph_mat .OR. elph_ahc )
  !
  if (twochem) occupations ='smearing' !this is needed to avoid init_twochem error check.
  CALL init_run_tpw()
!°°°°°°°°°°°°°°°°°°°° ACFDT TEST °°°°°°°°°°°°°°°°°°°°°°°°°
  IF (acfdt_is_active) THEN
    ! ACFDT mumerical derivative test: modify the potential
    IF (acfdt_num_der) vrs(ir_point,1)=vrs(ir_point,1) + delta_vrs
  ENDIF
!°°°°°°°°°°°°°°°°°END OF ACFDT TEST °°°°°°°°°°°°°°°°°°°°°°
!
  IF (do_band.AND..NOT.elph_mat) CALL non_scf_tpw ( )


  IF ( check_stop_now() ) THEN
!
!  In this case the code stops inside the band calculation. Save the
!  files and stop the pwscf run
!
     CALL punch( 'config' )
     CALL stop_run( -1 )
     CALL do_stop( 1 )
  ENDIF

  !
  IF (.NOT.reduce_io.and.do_band) THEN
     IF ( only_wfc ) THEN
        ! write wavefunctions to file in portable format
        CALL punch( 'all' )
     ELSE
        ! do not write wavefunctions: not sure why, I think
        ! they are written anyway in internal format - PG
        CALL punch( 'config' )
     END IF
  ENDIF
  !
  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  ext_restart=.FALSE.
  !
  IF (io_level > 0) THEN
     CALL close_files(.true.)
  ELSE
     CALL mp_barrier( intra_image_comm )  
  ENDIF
  !
  bands_computed=.TRUE.
  !
  CALL stop_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE run_nscf_tpw
