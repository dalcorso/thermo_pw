!
! Copyright (C) 2014-2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_cp_bs_g(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)
!
!  This subroutine receives the thermal expansion, the equilibrium volume,
!  the bulk modulus and the constant volume heat capacity all as a function
!  of temperature and computes the isobaric heat capacity, the
!  isoentropic bulk modulus, and the average gruneisen parameter 
!  as a function of temperature. 
!  Units:
!  cv_t and cp_t in Ry / K / cell
!  b0_t and b0_s in kbar
!  gamma_t adimensional
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : temp, ntemp
USE isoentropic,    ONLY : isobaric_heat_capacity, isoentropic_bulk_modulus, &
                           average_gruneisen

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: beta_t(ntemp), vmin_t(ntemp), b0_t(ntemp), cv_t(ntemp)
REAL(DP), INTENT(OUT) :: cp_t(ntemp), b0_s(ntemp), gamma_t(ntemp)
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
END SUBROUTINE compute_cp_bs_g

SUBROUTINE compute_cv_bs_g(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)
!
!  This subroutine receives the thermal expansion, the equilibrium volume,
!  the bulk modulus and the isobaric heat capacity all as a function
!  of temperature and computes the isochoric heat capacity, the
!  isoentropic bulk modulus, and the average gruneisen parameter 
!  as a function of temperature. 
!  Units:
!  cv_t and cp_t in Ry / K / cell
!  b0_t and b0_s in kbar
!  gamma_t adimensional
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : temp, ntemp
USE isoentropic,    ONLY : isobaric_heat_capacity, isoentropic_bulk_modulus, &
                           average_gruneisen

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: beta_t(ntemp), vmin_t(ntemp), b0_t(ntemp), cp_t(ntemp)
REAL(DP), INTENT(OUT) :: cv_t(ntemp), b0_s(ntemp), gamma_t(ntemp)
!
CALL isobaric_heat_capacity(vmin_t,b0_t,beta_t,temp,cv_t,ntemp)
!
!  The routine gives C_p-C_v. We subtract this quantity from C_p to get C_v
!
cv_t(1:ntemp) = cp_t(1:ntemp)-cv_t(1:ntemp)

CALL isoentropic_bulk_modulus(vmin_t,b0_t,beta_t,cp_t,temp,b0_s,ntemp)
!
!  The routine gives only the difference
!
b0_s(1:ntemp) = b0_s(1:ntemp)+b0_t(1:ntemp)

CALL average_gruneisen(vmin_t,b0_t,beta_t,cv_t,temp,gamma_t,ntemp)

RETURN
END SUBROUTINE compute_cv_bs_g
