!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_dos_kpoints()
!
!  This subroutine computes the mesh of k point for electronic dos
!  calculation. 
!
  USE control_paths, ONLY : disp_nqs, disp_q, disp_wq
  USE cell_base,     ONLY : bg
  USE symm_base,     ONLY : s, nrot, t_rev, time_reversal
  USE control_dos,   ONLY : nk1_d, nk2_d, nk3_d, k1_d, k2_d, k3_d

  IMPLICIT NONE
  INTEGER :: nqx
 
  IF ( ALLOCATED(disp_q) )    DEALLOCATE (disp_q)
  IF ( ALLOCATED(disp_wq) )   DEALLOCATE (disp_wq)

  nqx=nk1_d*nk2_d*nk3_d
  ALLOCATE ( disp_q(3,nqx) )
  ALLOCATE ( disp_wq(nqx) )

  CALL kpoint_grid ( nrot, time_reversal, .FALSE., s, t_rev, bg, nqx, &
                     k1_d, k2_d, k3_d, nk1_d, nk2_d, nk3_d, disp_nqs, disp_q, &
                     disp_wq)

   RETURN
END SUBROUTINE set_dos_kpoints
