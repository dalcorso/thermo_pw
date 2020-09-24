!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
SUBROUTINE print_polarization(polar_, flag)
!--------------------------------------------------------------
USE kinds,     ONLY : DP
USE constants, ONLY : electron_si, bohr_radius_si
USE cell_base, ONLY : alat, omega
USE io_global, ONLY : stdout
USE piezoelectric_tensor,  ONLY : polar_geo

USE mp,        ONLY : mp_sum
USE mp_world,  ONLY : world_comm
USE mp_images, ONLY : nproc_image

IMPLICIT NONE
REAL(DP), INTENT(IN) :: polar_(3)
LOGICAL, INTENT(IN) :: flag
REAL(DP) :: fact, polar(3)

IF (flag) THEN
   CALL mp_sum(polar_geo, world_comm)
   polar_geo=polar_geo / nproc_image
   polar=polar_geo(:,1)
ELSE
   polar=polar_
ENDIF

WRITE(stdout,'(/,5x,"The Berry phase polarization of this solid in &
                     cartesian coordinates is:")')

WRITE(stdout,'(/,5x,"(", 2(f10.5,","), f10.5, "   ) phase ")')  polar(:)
fact= alat / omega
WRITE(stdout,'(/,5x,"(", 2(f10.5,","), f10.5, "   ) e/(a.u.)^2")') &
                                   polar(:) * fact
fact= alat * electron_si / (bohr_radius_si)**2 / omega
WRITE(stdout,'(/,5x,"(", 2(f10.5,","), f10.5, "   ) C/m^2")') &
                                                 polar(:) * fact
IF (flag) THEN
   WRITE(stdout,'(/,5x,"Please note that only differences of polarization")')
   WRITE(stdout,'(5x,"have physical meaning. If you know that this vector")')
   WRITE(stdout,'(5x,"should be zero by symmetry (see above), you can use")')
   WRITE(stdout,'(5x,"this vector as a zero reference.")')
   WRITE(stdout,'(5x,"If your solid can have a spontaneous polarization")')
   WRITE(stdout,'(5x,"you need to subtract an appropriate zero reference")')
   WRITE(stdout,'(5x,"to this vector to find the spontaneous polarization.")')
ENDIF
WRITE(stdout,'(/,20x,40("-"),/)')

RETURN
END SUBROUTINE print_polarization
