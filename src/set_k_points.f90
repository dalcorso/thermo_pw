!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_k_points()
!
!  This subroutine copies the k points path in the input k points of
!  the pw.x code.
!
  USE control_paths,    ONLY : disp_nqs, disp_q, disp_wq
  USE start_k,          ONLY : init_start_k
  IMPLICIT NONE
  INTEGER :: nk1, nk2, nk3, k1, k2, k3
  CHARACTER(len=80) :: k_points 

  nk1=0
  nk2=0
  nk3=0
  k1=0
  k2=0
  k3=0
  k_points='tpiba'
 
  CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, disp_nqs, disp_q, &
                                     disp_wq )
  RETURN
END SUBROUTINE set_k_points
