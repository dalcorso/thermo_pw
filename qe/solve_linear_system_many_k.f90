!
! Copyright (C) 2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE solve_linear_system_many_k(dvpsi, dpsi, h_diag, thresh, &
                                             lter, nk, npe, nsolv)
!----------------------------------------------------------------------------
!
!   This routine is a driver that solves the linear system of the
!   phonon or electric field perturbation and possibly
!   changes the sign of the magnetic field in the Hamiltonian
!
USE kinds,        ONLY : DP
USE scf,          ONLY : vrs
USE klist,        ONLY : nks, ngk
USE control_lr,   ONLY : nbnd_occ
USE wvfct,        ONLY : nbnd, npwx, et
USE noncollin_module, ONLY : npol
USE uspp,         ONLY : okvan, deeq_nc
USE qpoint,       ONLY : ikks, ikqs, nksq
USE qpoint_aux,   ONLY : ikmks
USE lr_nc_mag,    ONLY : deeq_nc_save
USE io_global,    ONLY : stdout
USE many_k_ph_mod, ONLY : current_ikb_ph, startkb_ph

IMPLICIT NONE

INTEGER :: nk, npe, nsolv
INTEGER :: lter(nk*npe*nsolv)
COMPLEX(DP) :: dvpsi (npwx*npol, nbnd*nk*npe*nsolv)
COMPLEX(DP) ::  dpsi (npwx*npol, nbnd*nk*npe*nsolv)
REAL(DP) :: h_diag (npwx*npol, nbnd*nk*npe*nsolv)
REAL(DP) :: thresh

INTEGER :: ikk, ikq, ikmk, id, ik, ik1, ipert, isolv, npwq
REAL(DP), ALLOCATABLE ::  anorm(:)
LOGICAL, ALLOCATABLE :: conv_root(:)
EXTERNAL ch_psi_all_many_k, cg_psi_many_k

#if defined(__CUDA)
ATTRIBUTES(DEVICE) :: h_diag, dvpsi, dpsi
#endif

ALLOCATE(conv_root(nk*npe*nsolv))
ALLOCATE(anorm(nk*npe*nsolv))

conv_root=.TRUE.
CALL cgsolve_all_many_k (ch_psi_all_many_k, cg_psi_many_k, et, dvpsi, dpsi, &
     h_diag, npwx, ngk, thresh, lter, conv_root, anorm, nbnd, npol, nk, &
     npe, nsolv, ikks, ikqs, nbnd_occ, nks, nksq)

DO ik1=1,nk
   ik=ik1+startkb_ph(current_ikb_ph)
   DO isolv=1, nsolv
      DO ipert=1,npe
         id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
         IF (.NOT.conv_root(id)) WRITE( stdout, '(5x,"kpoint-pert.",2i4, &
        &" solve_linter: root not converged ",es10.3)') ik, ipert, anorm(id)
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(anorm)
DEALLOCATE(conv_root)

RETURN
END SUBROUTINE solve_linear_system_many_k
