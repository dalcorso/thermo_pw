!
! Copyright (C) 2014-2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE compute_cp_bs_g(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t,&
                           subtract_el)
!---------------------------------------------------------------------------
!
!  This subroutine receives the thermal expansion, the equilibrium volume,
!  the bulk modulus and the constant volume heat capacity as a function
!  of temperature and computes the isobaric heat capacity, the
!  isoentropic bulk modulus, and the average Gruneisen parameter 
!  as a function of temperature. 
!
!  For testing purposes it is possible to compute the Gruneisen parameter
!  neglecting the electronic contribution on the heat capacity, even
!  if this contribution is included everywhere else (subtract_el=.TRUE.)
!
!  Units:
!  cv_t and cp_t in Ry / K / cell
!  b0_t and b0_s in kbar
!  gamma_t adimensional
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : temp, ntemp
USE isoentropic,    ONLY : isobaric_heat_capacity, isoentropic_bulk_modulus, &
                           average_gruneisen
USE el_anharmonic,  ONLY : el_ce_t

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: beta_t(ntemp), vmin_t(ntemp), b0_t(ntemp), cv_t(ntemp)
REAL(DP), INTENT(OUT) :: cp_t(ntemp), b0_s(ntemp), gamma_t(ntemp)
LOGICAL, INTENT(IN)   :: subtract_el
REAL(DP) :: aux(ntemp)
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

IF (subtract_el) THEN
   aux=cv_t-el_ce_t
ELSE
   aux=cv_t
ENDIF
CALL average_gruneisen(vmin_t,b0_t,beta_t,aux,temp,gamma_t,ntemp)

RETURN
END SUBROUTINE compute_cp_bs_g
!
!---------------------------------------------------------------------------
SUBROUTINE compute_cp_bs_gp(beta_ptt, vmin_ptt, b0_ptt, cv_ptt, cp_ptt, &
                   b0_s_ptt, gamma_ptt, el_ce_ptt, itemp, subtract_el)
!---------------------------------------------------------------------------
!
!  This subroutine receives the thermal expansion, the equilibrium volume,
!  the bulk modulus and the constant volume heat capacity all as a function
!  of pressure and computes the isobaric heat capacity, the
!  isoentropic bulk modulus, and the average gruneisen parameter 
!  as a function of pressure. 
!
!  For testing purposes it is possible to compute the Gruneisen parameter
!  neglecting the electronic contribution on the heat capacity, even
!  if this contribution is included everywhere else (subtract_el=.TRUE.)
!
!  Units:
!  cv_t and cp_t in Ry / K / cell
!  b0_t and b0_s in kbar
!  gamma_t adimensional
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE temperature,      ONLY : temp, ntemp
USE control_pressure, ONLY : npress
USE mp,               ONLY : mp_sum
USE mp_world,         ONLY : world_comm

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
LOGICAL, INTENT(IN) :: subtract_el
REAL(DP), INTENT(IN)  :: beta_ptt(npress), vmin_ptt(npress), b0_ptt(npress), &
                         cv_ptt(npress), el_ce_ptt(npress)
REAL(DP), INTENT(OUT) :: cp_ptt(npress), b0_s_ptt(npress), gamma_ptt(npress)

REAL(DP), ALLOCATABLE :: cpmcv(:), bsmbt(:)
REAL(DP) :: aux, auxc(npress)
INTEGER :: ipress, startp, lastp 

ALLOCATE(cpmcv(npress))
ALLOCATE(bsmbt(npress))

CALL divide(world_comm, npress, startp, lastp)

cpmcv=0.0_DP
DO ipress=startp,lastp
   cpmcv(ipress) = temp(itemp) * vmin_ptt(ipress) * &
                           beta_ptt(ipress)**2 * b0_ptt(ipress) / ry_kbar
ENDDO
CALL mp_sum(cpmcv,world_comm)
!
!  The routine gives only the difference
!
cp_ptt(1:npress) = cpmcv(1:npress)+cv_ptt(1:npress)

bsmbt=0.0_DP
DO ipress=startp,lastp
   aux = temp(itemp) * vmin_ptt(ipress) * beta_ptt(ipress)**2 * &
                             b0_ptt(ipress) / cp_ptt(ipress) / ry_kbar
   bsmbt(ipress) =  b0_ptt(ipress) * aux / ( 1.0_DP - aux )
END DO
CALL mp_sum(bsmbt,world_comm)
!
!  The routine gives only the difference
!
b0_s_ptt(1:npress) = bsmbt(1:npress)+b0_ptt(1:npress)
!
!  The average Gruneisen parameter
!
gamma_ptt=0.0_DP
IF (subtract_el) THEN
   auxc=cv_ptt-el_ce_ptt
ELSE
   auxc=cv_ptt
ENDIF

DO ipress=startp,lastp
   gamma_ptt(ipress)=beta_ptt(ipress) * b0_ptt(ipress) * vmin_ptt(ipress) &
                                      / auxc(ipress) / ry_kbar
ENDDO
CALL mp_sum(gamma_ptt,world_comm)

DEALLOCATE(cpmcv)
DEALLOCATE(bsmbt)

RETURN
END SUBROUTINE compute_cp_bs_gp

!---------------------------------------------------------------------------
SUBROUTINE compute_cv_bs_g(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)
!---------------------------------------------------------------------------
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
