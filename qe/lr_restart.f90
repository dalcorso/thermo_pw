!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE lr_restart_tpw(iter_restart,rflag)
  !---------------------------------------------------------------------
  !
  ! Restart the Lanczos recursion
  !
  USE kinds,                ONLY : DP 
  USE io_global,            ONLY : stdout, ionode_id
  USE control_flags,        ONLY : gamma_only
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE io_files,             ONLY : tmp_dir, prefix, diropn, wfc_dir
  USE lr_lanczos,           ONLY : evc1, evc1_new, evc1_old, sevc1,     &
                                   beta_store, gamma_store, zeta_store, &
                                   iulanczos, iunrestart, nwordrestart, &
                                   lanczos_steps
  USE lr_global,            ONLY : rpert
  USE wvfct,                ONLY : nbnd, npwx
  USE io_global,            ONLY : ionode
  USE mp,                   ONLY : mp_bcast
  USE mp_images,            ONLY : intra_image_comm
  USE noncollin_module,     ONLY : npol
  USE qpoint,               ONLY : nksq

  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: iter_restart
  LOGICAL, INTENT(OUT) :: rflag
  !
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  CHARACTER (len=24) :: bgz_suffix
  ! local variables
  !
  INTEGER :: i, iunres
  INTEGER :: ik, ig, ip
  LOGICAL :: exst
  CHARACTER(len=256) :: tempfile, filename, tmp_dir_saved
  INTEGER :: pol_index
  REAL(kind=dp) :: norm0(3)
  !
  rflag=.TRUE.
  pol_index = 1
  !
  ! Optical case: recompute the kinetic-energy g2kin and 
  ! beta functions vkb (needed only in the US case).
  ! Note, this is done only in the gamma_only case,
  ! because in the k-points version all is recomputed
  ! on-the-fly for every k point.
  !
  !
  ! Reading Lanczos coefficients
  bgz_suffix = TRIM ( ".beta_gamma_z." )
  !
  filename = TRIM(prefix) // TRIM(bgz_suffix) // TRIM("dat")
  tempfile = TRIM(tmp_dir) // TRIM(filename)
  !
  INQUIRE (FILE = TRIM(tempfile), EXIST = exst)
  !
  IF (.NOT.exst) THEN
     !
     WRITE(stdout,*) "WARNING: " // TRIM(filename) // " does not exist"
     rflag = .FALSE.
     RETURN
     !
  ENDIF

  iunres = find_free_unit()
  !
  ! Ionode only reads the tempfile
  !
  IF (ionode) THEN
  !
  ! Read and broadcast beta, gamma, and zeta.
  !
     OPEN(iunres, file = TRIM(tempfile), form = 'formatted', status = 'old')
     !
     READ(iunres, *, end=301, err=303) iter_restart
     !
     IF ( iter_restart >= lanczos_steps ) iter_restart=lanczos_steps
     !
     READ(iunres, *, end=301, err=303) norm0(pol_index)
     !
     ! X. Ge: New structure for the pseudo-Hermitian algorithm.
     !
     beta_store(1) = norm0(pol_index)
     !
     DO i=1,7
        READ(iunres, *, end=301, err=303)
     ENDDO
     ! 
     DO i=1,iter_restart-1
        READ(iunres,*,end=301,err=303) beta_store(i+1)
        READ(iunres,*,end=301,err=303) gamma_store(i+1)
        READ(iunres,*,end=301,err=303) zeta_store (:,:,i)
     ENDDO
     !
     READ(iunres,*,end=301,err=303) beta_store(iter_restart)
     READ(iunres,*,end=301,err=303) gamma_store(iter_restart)
     READ(iunres,*,end=301,err=303) zeta_store (:,:,iter_restart)
     !
     CLOSE(iunres)
     !
  ENDIF
  CALL mp_bcast (iter_restart, ionode_id, intra_image_comm)
  CALL mp_bcast (norm0(pol_index), ionode_id, intra_image_comm)
  CALL mp_bcast (beta_store(:), ionode_id, intra_image_comm)
  CALL mp_bcast (gamma_store(:), ionode_id, intra_image_comm)
  CALL mp_bcast (zeta_store(:,:,:), ionode_id, intra_image_comm)
  !
  ! Optical case: read projection
  ! 
  iter_restart = iter_restart + 1 
  !
  ! Parallel reading
  ! Note: Restart files are always in outdir
  ! Reading Lanczos vectors
  !
  nwordrestart = 2 * nbnd * npwx * npol * nksq * rpert
  !
  CALL diropn( iunres, 'restart_lanczos.'//trim(int_to_char(1)), &
                                                         nwordrestart, exst)
  !
  CALL davcio(evc1(:,:,:,:), nwordrestart, iunres, 1, -1)
  CALL davcio(evc1_old(:,:,:,:), nwordrestart, iunres, 2, -1)
  CLOSE( unit = iunres)
  iunrestart=iunres
  !
  RETURN
  !
301 CALL errore ('lr_restart_tpw', 'A File is corrupted, file ended &
                                                          &unexpectedly', 1)
303 CALL errore ('lr_restart_tpw', 'A File is corrupted, error in reading &
                                                                   &data', 1)
  rflag=.FALSE.
  RETURN
  !
END SUBROUTINE lr_restart_tpw
