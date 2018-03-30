!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE lanczos_write_restart_tpw(lr_iteration)
  !---------------------------------------------------------------------
  ! 
  ! This subroutine writes the vectors necessary to
  ! restart the Lanczos recursion.
  !
 USE kinds,                ONLY : DP
 USE io_files,             ONLY : tmp_dir, prefix, diropn
 USE lr_lanczos,           ONLY : evc1, evc1_new, evc1_old, sevc1, beta_store, &
                                  gamma_store, zeta_store, iulanczos, & 
                                  iunrestart, nwordrestart
 USE lr_global,            ONLY : rpert, size_evc1
 USE wvfct,                ONLY : nbnd, npwx
 USE fft_base,             ONLY : dfftp
 USE io_global,            ONLY : ionode, stdout
 USE klist,                ONLY : nks, nelec
 USE noncollin_module,     ONLY : nspin_mag, noncolin, npol
 use lsda_mod,             ONLY : nspin
 USE cell_base,            ONLY : alat, omega
 USE qpoint,               ONLY : nksq, xq
 !
 IMPLICIT NONE
 CHARACTER(LEN=6), EXTERNAL :: int_to_char
 !
 ! local variables
 ! 
 REAL(DP) :: norm0(3)
 !
 INTEGER, EXTERNAL :: find_free_unit
 INTEGER, INTENT(IN) :: lr_iteration
 INTEGER :: i, j, pol_index, iunres
 CHARACTER (len=24) :: bgz_suffix
 CHARACTER(len=256) :: tempfile, filename
 LOGICAL :: exst
 real(kind=dp) :: degspin
 !
 ! If there is only one polarization dir, storage is one rank less.
 !
 pol_index = 1
 iunres = find_free_unit()
 ! 
 IF (ionode) THEN
    !
    ! Writing beta, gamma and zeta coefficients.
    !
    bgz_suffix = TRIM ( ".beta_gamma_z." )
    filename = trim(prefix) // trim(bgz_suffix) // trim("dat")
    tempfile = trim(tmp_dir) // trim(filename)
    !
    OPEN (iunres, file = tempfile, form = 'formatted', status = 'unknown')
    WRITE(iunres,*) LR_iteration
    !
    norm0(pol_index) = beta_store(1)
    WRITE(iunres,*) norm0(pol_index)
    !
    IF (nspin==2) THEN
       degspin = 1.0d0
    ELSE
       degspin = 2.0d0
    ENDIF
    IF (noncolin) degspin = 1.0d0
    !
    ! Write the degenaracy wrt spin
    !
    WRITE(iunres,*) degspin
    !
    ! ------ Needed for EELS ----------
    !
    ! Write the lattice parameter
    !
    WRITE(iunres,*) alat
    !
    ! Write the unit-cell volume
    !
    WRITE(iunres,*) omega
    !
    ! Write the number of valence (and semicore electrons) in the unit
    ! cell
    !
    WRITE(iunres,*) nelec
    !
    ! Write the components of the transferred momentum
    !
    WRITE(iunres,*) xq(1)
    WRITE(iunres,*) xq(2)
    WRITE(iunres,*) xq(3)
    !
    !-----------------------------------
    !
    DO i=1, lr_iteration-1
       !
       WRITE(iunres,*) beta_store(i+1)
       WRITE(iunres,*) gamma_store(i+1)
       !
       ! This is absolutely necessary for cross platform compatibility
       !
       WRITE(iunres,*) zeta_store (:,:,i)
       !
    ENDDO
    !
    ! X. Ge: Compatible with the old version. The beta & gamma will not be
    ! used in the spectrum calculation.
    !
    WRITE(iunres,*) beta_store(lr_iteration)             
    WRITE(iunres,*) gamma_store(lr_iteration)             
    WRITE(iunres,*) zeta_store (:,:,lr_iteration)        
    !
    CLOSE(iunres)
    !
 ENDIF
 !
 ! Parallel writing operations
 !
 ! Note: Restart files are writen in outdir.
 ! If you do not want them to be written,
 ! just disable restart saving completely.
 !
 ! Writing wavefuncion files for restart
 !
 nwordrestart = 2 * nbnd * npwx * npol * nksq * rpert
 !
 CALL diropn ( iunres,'restart_lanczos.'//trim(int_to_char(1)),&
                                                          nwordrestart, exst)
 !
 CALL davcio(evc1(:,:,:,:),nwordrestart,iunres,1,1)
 CALL davcio(evc1_old(:,:,:,:),nwordrestart,iunres,2,1)
 CLOSE(UNIT = iunres)
 iunrestart=iunres
 !
 RETURN
END SUBROUTINE lanczos_write_restart_tpw
