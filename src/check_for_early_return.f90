!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE check_for_early_return(energy_geo, nwork, ngeom, lreturn)
!----------------------------------------------------------------------
!
USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : &
                              start_geometry_qha, last_geometry_qha
IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom
REAL(DP), INTENT(IN) :: energy_geo(nwork)
LOGICAL, INTENT(OUT) :: lreturn

INTEGER :: work_base, base_ind, igeom, iwork

lreturn=.FALSE.
work_base = nwork / ngeom
DO igeom=start_geometry_qha,last_geometry_qha
   base_ind=(igeom-1)*work_base
   DO iwork=1,work_base
      lreturn=lreturn.OR.(ABS(energy_geo(base_ind+iwork))<1.D-10)
   ENDDO
ENDDO

RETURN
END SUBROUTINE check_for_early_return
