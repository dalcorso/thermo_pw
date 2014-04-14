!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_elastic_cons_work( nwork, ngeo )
USE kinds, ONLY : DP
USE cell_base, ONLY : ibrav
USE control_elastic_constants, ONLY : delta_epsilon
USE elastic_constants, ONLY : epsilon_voigt, epsilon_geo, sigma_geo, &
                              trans_epsilon
IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeo
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min
INTEGER :: igeo, iwork

IF (ibrav==1 .OR. ibrav==2 .OR. ibrav==3) THEN
   nwork = 2 * ngeo
   ALLOCATE( epsilon_geo(3, 3, nwork) ) 
   ALLOCATE( sigma_geo(3, 3, nwork) )
   ALLOCATE( epsilon_voigt(6, nwork) )
   epsilon_voigt=0.0_DP
   sigma_geo=0.0_DP
   epsilon_geo=0.0_DP
   epsilon_min= - delta_epsilon * (ngeo / 2)
   DO igeo=1,ngeo
      epsilon_voigt(3, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 ) 
      epsilon_voigt(4, ngeo + igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
      epsilon_voigt(5, ngeo + igeo) = epsilon_voigt(4, ngeo + igeo)
      epsilon_voigt(6, ngeo + igeo) = epsilon_voigt(4, ngeo + igeo)
   ENDDO
   DO iwork = 1, nwork
      CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
   ENDDO
ELSE
   CALL errore('set_elastic_cons_work', & 
               'elastic constants for this ibrav not implemented', 1)
ENDIF

RETURN
END SUBROUTINE set_elastic_cons_work
