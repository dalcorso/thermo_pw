!
! Copyright (C) 2014-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_cp(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)
!
!  This subroutine receives the constant volume heat capacity at
!  several volumes, the equilibrium volume as a function of temperature,
!  the bulk modulus as a function of the temperature and the 
!  thermal expansion as a function of temperature and computes the
!  isobaric heat capacity as a function of temperature, the
!  isoentropic bulk modulus as a function of temperature, 
!  the average gruneisen parameter as a function of temperature
!  cv_t and cp_t in Ry / K / cell
!  b0_t and b0_s in kbar
!  gamma_t adimensional
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : temp, ntemp
USE isoentropic,    ONLY : isobaric_heat_capacity, isoentropic_bulk_modulus, &
                           average_gruneisen

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: beta_t(ntemp), vmin_t(ntemp), b0_t(ntemp)
REAL(DP), INTENT(OUT) :: cv_t(ntemp), cp_t(ntemp), b0_s(ntemp), gamma_t(ntemp)
!
CALL isobaric_heat_capacity(vmin_t,b0_t,beta_t,temp,cp_t,ntemp)
!
!  The routine gives only the difference
!
cp_t(1:ntemp) = cp_t(1:ntemp)+cv_t(1:ntemp)

CALL isoentropic_bulk_modulus(vmin_t,b0_t,beta_t,cp_t,temp,b0_s,ntemp)
!
!  The routine gives only the difference
!
b0_s(1:ntemp) = b0_s(1:ntemp)+b0_t(1:ntemp)

CALL average_gruneisen(vmin_t,b0_t,beta_t,cv_t,temp,gamma_t,ntemp)

RETURN
END SUBROUTINE compute_cp
