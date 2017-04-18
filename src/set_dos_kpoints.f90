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
  USE control_dosq,  ONLY : dos_q, dos_wq
  USE control_eldos, ONLY : nk1_d, nk2_d, nk3_d, k1_d, k2_d, k3_d
  USE start_k,       ONLY : init_start_k

  IMPLICIT NONE
  INTEGER :: nqx
 
  IF ( ALLOCATED(dos_q) )    DEALLOCATE (dos_q)
  IF ( ALLOCATED(dos_wq) )   DEALLOCATE (dos_wq)

  nqx=1
  ALLOCATE ( dos_q(3,nqx) )
  ALLOCATE ( dos_wq(nqx) )

  CALL init_start_k ( nk1_d, nk2_d, nk3_d, k1_d, k2_d, k3_d, 'automatic', 0,&
                               dos_q, dos_wq )

   RETURN
END SUBROUTINE set_dos_kpoints
