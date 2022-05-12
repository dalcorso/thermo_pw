!
! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE empirical_free_energy()

USE kinds,  ONLY : DP
USE thermo_mod, ONLY : ngeo, omega_geo
USE ions_base, ONLY : nat
USE temperature, ONLY : temp, ntemp
USE control_emp_free_ener, ONLY: emp_ener, emp_free_ener, emp_entr, emp_ce
USE control_emp_free_ener, ONLY : efe, v0p, alpha1, alpha2
USE model_free_energy, ONLY : emp_quadratic_t, emp_anhar_am

IMPLICIT NONE
REAL(DP) :: t, vol, free_ener
INTEGER :: itemp, igeom

free_ener=0.0_DP
IF (efe==1) THEN
   DO itemp=1,ntemp
      t=temp(itemp)
      DO igeom=1, ngeo(1)
         vol=omega_geo(igeom)
         CALL emp_quadratic_t(t, vol, alpha1, alpha2, &
            emp_free_ener(itemp,igeom), emp_ener(itemp,igeom), &
            emp_entr(itemp,igeom), emp_ce(itemp,igeom))
      ENDDO
   ENDDO
ELSEIF (efe==2) THEN
   DO itemp=1,ntemp
      t=temp(itemp)
      DO igeom=1, ngeo(1)
         vol=omega_geo(igeom)
         CALL emp_anhar_am(t, vol, alpha1, alpha2, v0p, nat, &
            emp_free_ener(itemp,igeom), emp_ener(itemp,igeom), &
            emp_entr(itemp,igeom), emp_ce(itemp,igeom))
      ENDDO
   ENDDO
ELSE
   CALL errore('empirical_free_energy','free energy not available',1)
ENDIF

RETURN
END SUBROUTINE empirical_free_energy
