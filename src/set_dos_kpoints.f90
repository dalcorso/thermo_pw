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
  USE control_eldos, ONLY : nk1_d, nk2_d, nk3_d, k1_d, k2_d, k3_d, dos_k, &
                            dos_wk
  USE start_k,       ONLY : init_start_k

  IMPLICIT NONE
  INTEGER :: nqx
 
  IF ( ALLOCATED(dos_k) )    DEALLOCATE (dos_k)
  IF ( ALLOCATED(dos_wk) )   DEALLOCATE (dos_wk)

  nqx=1
  ALLOCATE ( dos_k(3,nqx) )
  ALLOCATE ( dos_wk(nqx) )

  CALL init_start_k ( nk1_d, nk2_d, nk3_d, k1_d, k2_d, k3_d, 'automatic', 0,&
                               dos_k, dos_wk )

   RETURN
END SUBROUTINE set_dos_kpoints
