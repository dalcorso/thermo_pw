!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solve_linear_system(dvpsi, dpsi, h_diag, thresh, ik, isolv, lter)
!---------------------------------------------------------------------------
!
!   This routine is a driver that solves the linear system of the
!   phonon or electric field perturbation and possibly
!   changes the sign of the magnetic field in the Hamiltonian
!
USE kinds,        ONLY : DP
USE scf,          ONLY : vrs
#if defined(__CUDA)
  USE many_k_mod,   ONLY : vrs_d
#endif
USE klist,        ONLY : ngk
USE control_lr,   ONLY : nbnd_occ
USE wvfct,        ONLY : nbnd, npwx, et
USE noncollin_module, ONLY : npol
USE uspp,         ONLY : okvan, deeq_nc
USE qpoint,       ONLY : ikks, ikqs
USE qpoint_aux,   ONLY : ikmks
USE lr_nc_mag,    ONLY : deeq_nc_save
USE io_global,    ONLY : stdout

IMPLICIT NONE

INTEGER :: ik, isolv, lter
COMPLEX(DP) :: dvpsi (npwx*npol, nbnd)
COMPLEX(DP) ::  dpsi (npwx*npol, nbnd)
REAL(DP) :: h_diag (npwx*npol, nbnd)
REAL(DP) :: thresh

INTEGER :: ikk, ikq, ikmk, npwq
REAL(DP) ::  anorm
LOGICAL :: conv_root
EXTERNAL ch_psi_all, cg_psi

ikk=ikks(ik)
ikq=ikqs(ik)
npwq=ngk(ikq)
ikmk=ikk
IF (isolv==2) THEN
   ikmk=ikmks(ik)
   vrs(:,2:4)=-vrs(:,2:4)
#if defined(__CUDA)
   vrs_d = vrs
#endif
   IF (okvan) THEN
      deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,2)
      !$acc update device(deeq_nc)
   ENDIF
ENDIF

conv_root=.TRUE.
CALL cgsolve_all (ch_psi_all, cg_psi, et(1,ikmk), dvpsi, dpsi, &
                  h_diag, npwx, npwq, thresh, ik, lter, conv_root, &
                  anorm, nbnd_occ(ikk), npol )

IF (isolv==2) THEN
   vrs(:,2:4)=-vrs(:,2:4)
#if defined(__CUDA)
   vrs_d = vrs
#endif
   IF (okvan) THEN
      deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,1)
      !$acc update device(deeq_nc)
   ENDIF
ENDIF

IF (.NOT.conv_root) WRITE( stdout, '(5x,"kpoint",i4, " solve_linter: &
      &root not converged ",es10.3)') ik, anorm

RETURN
END SUBROUTINE solve_linear_system

