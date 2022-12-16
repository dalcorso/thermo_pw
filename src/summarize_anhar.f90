!
! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE summarize_anhar()
!---------------------------------------------------------------------
!
!  This routine writes on output a few thermodynamic quantities
!  at T=0 K, at the lower and at the higher available temperatures 
!  and at the temperature closer to T=300 K.
!
USE kinds,                 ONLY : DP
USE constants,             ONLY : rydberg_si, avogadro, bohr_radius_si
USE temperature,           ONLY : temp, ntemp, itemp300, deltat
USE control_mur,           ONLY : vmin, b0, b01, b02, emin
USE control_ev,            ONLY : ieos
USE equilibrium_conf,      ONLY : celldm0
USE control_mur,           ONLY : lmurn
USE anharmonic,            ONLY : vmin_t, b0_t, b01_t, b02_t, &
                                  beta_t, cv_t, cp_t, gamma_t, b0_s, &
                                  celldm_t
USE io_global,             ONLY : stdout

IMPLICIT NONE
REAL(DP) :: cv_fact, dbetadt, dbsdt, dbtdt

cv_fact=rydberg_si*avogadro
!
!  T=0 K
!
WRITE(stdout,'(/,5x,"Thermodynamic properties at T = 0 K")') 
WRITE(stdout,'(5x,"V_0 =",f20.8," (a.u.)^3", f20.8, " A^3")') vmin, vmin * &
                                                    bohr_radius_si**3 * 1.D30
WRITE(stdout,'(5x,"a_0 =",f20.8," a.u.", f24.8," A")') celldm0(1), &
                                             celldm0(1)*bohr_radius_si * 1.D10
IF (lmurn) THEN
   WRITE(stdout,'(5x,"B_0 =",f17.5," kbar")') b0
   WRITE(stdout,'(5x,"dB_0/dp =",f13.5)') b01
   IF (ieos==2) WRITE(stdout,'(5x,"d^2B_0/dp^2 =",f9.5," 1/kbar")') b02
ENDIF
!
!  Effect of zero point motion
!
WRITE(stdout,'(/,5x,"Thermodynamic properties at T = ",f8.2," K")') temp(2)

WRITE(stdout,'(5x,"V_0 =",f20.8," (a.u.)^3",f20.8," A^3")') vmin_t(2), &
                         vmin_t(2) * bohr_radius_si**3 * 1.D30
WRITE(stdout,'(5x,"a_0 =",f20.8," a.u.", f24.8, " A")') celldm_t(1,2), &
                            celldm_t(1,2)* bohr_radius_si * 1.D10
WRITE(stdout,'(5x,"B_T =",f17.5," kbar")') b0_t(2)
IF (lmurn) THEN
   WRITE(stdout,'(5x,"dB_T/dp =",f13.5)') b01_t(2)
   IF (ieos==2) WRITE(stdout,'(5x,"d^2B_0/dp^2 =",f9.5," 1/kbar")') b02_t(2)
ENDIF
!
!  Room temperature
!
dbetadt = (beta_t(itemp300+1) - beta_t(itemp300-1) ) / 2.0_DP / deltat
dbtdt = (b0_t(itemp300+1) - b0_t(itemp300-1) ) / 2.0_DP / deltat
dbsdt = (b0_s(itemp300+1) - b0_s(itemp300-1) ) / 2.0_DP / deltat
WRITE(stdout,'(/,5x,"Thermodynamic properties at T = ",f8.2," K")') &
                                                           temp(itemp300)
WRITE(stdout,'(5x,"V_0 =",f20.8," (a.u.)^3", f20.8, " A^3")') &
              vmin_t(itemp300), vmin_t(itemp300)* bohr_radius_si**3 * 1.D30
WRITE(stdout,'(5x,"a_0 =",f20.8," a.u.",f24.8," A")') celldm_t(1,itemp300), &
                          celldm_t(1,itemp300)*bohr_radius_si * 1.D10
WRITE(stdout,'(5x,"B_T =",f17.5," kbar")') b0_t(itemp300)

IF (lmurn) THEN
   WRITE(stdout,'(5x,"dB_T/dp =",f13.5)') b01_t(itemp300)
   WRITE(stdout,'(5x,"dB_T/dT =",f13.5," kbar/K")') dbtdt
   IF (ieos==2) WRITE(stdout,'(5x,"d^2B_0/dp^2 =",f9.5," 1/kbar")') &
                                                           b02_t(itemp300)
ENDIF
WRITE(stdout,'(5x,"beta =",f19.8," x10^-6 1/K")') beta_t(itemp300)*1.D6
WRITE(stdout,'(5x,"d beta/dT=",f15.8," x10^-8 1/K^2")') dbetadt*1.D8
WRITE(stdout,'(5x,"d beta/dp=",f15.8," x10^-7 1/kbar/K")') &
                                           dbtdt*1.D7/b0_t(itemp300)**2
WRITE(stdout,'(5x,"B_s =",f17.5," kbar")') b0_s(itemp300)
WRITE(stdout,'(5x,"d B_s/dT =",f12.5," kbar / K")') dbsdt
WRITE(stdout,'(5x,"Cv =",f21.8," J / K / mol")') cv_t(itemp300) * cv_fact
WRITE(stdout,'(5x,"Cp =",f21.8," J / K / mol")') cp_t(itemp300) * cv_fact
WRITE(stdout,'(5x,"gamma =",f16.6)') gamma_t(itemp300)
WRITE(stdout,'(5x,"delta_T =",f14.6)') -dbtdt/beta_t(itemp300)/b0_t(itemp300)
WRITE(stdout,'(5x,"delta_S =",f14.6)') -dbsdt/beta_t(itemp300)/b0_s(itemp300)
!
!  The maximum calculated temperature
!
dbetadt = (beta_t(ntemp-1) - beta_t(ntemp-3) ) / 2.0_DP / deltat
dbtdt = (b0_t(ntemp-1) - b0_t(ntemp-3) ) / 2.0_DP / deltat
dbsdt = (b0_s(ntemp-1) - b0_s(ntemp-3) ) / 2.0_DP / deltat
WRITE(stdout,'(/,5x,"Thermodynamic properties at T = ",f8.2," K")') &
                                                           temp(ntemp-2)
WRITE(stdout,'(5x,"V_0 =",f20.8," (a.u.)^3",f20.8," A^3")') vmin_t(ntemp-2), &
                         vmin_t(ntemp-2) * bohr_radius_si**3 * 1.D30
WRITE(stdout,'(5x,"a_0 =",f20.8," a.u.", f24.8, " A")') celldm_t(1,ntemp-2), &
                            celldm_t(1,ntemp-2)* bohr_radius_si * 1.D10
WRITE(stdout,'(5x,"B_T =",f17.5," kbar")') b0_t(ntemp-2)
IF (lmurn) THEN
   WRITE(stdout,'(5x,"dB_T/dp =",f13.5)') b01_t(ntemp-2)
   WRITE(stdout,'(5x,"dB_T/dT =",f13.5," kbar/K")') dbtdt
   IF (ieos==2) WRITE(stdout,'(5x,"d^2B_0/dp^2 =",f9.5," 1/kbar")') &
                                                        b02_t(ntemp-2)
ENDIF

WRITE(stdout,'(5x,"beta =",f19.8," x10^-6 1/K")') beta_t(ntemp-2)*1.D6
WRITE(stdout,'(5x,"d beta/dT=",f15.8," x10^-8 1/K^2")') dbetadt*1.D8
WRITE(stdout,'(5x,"d beta/dp=",f15.8," x10^-7 1/kbar/K")') &
                                           dbtdt*1.D7/b0_t(ntemp-2)**2
WRITE(stdout,'(5x,"B_s =",f17.5," kbar")') b0_s(ntemp-2)
WRITE(stdout,'(5x,"d B_s/dT =",f12.5," kbar / K")') dbsdt
WRITE(stdout,'(5x,"Cv =",f21.8," J / K / mol")') cv_t(ntemp-2) * cv_fact
WRITE(stdout,'(5x,"Cp =",f21.8," J / K / mol")') cp_t(ntemp-2) * cv_fact
WRITE(stdout,'(5x,"gamma =",f16.6)') gamma_t(ntemp-2)
WRITE(stdout,'(5x,"delta_T =",f14.6)') -dbtdt/beta_t(ntemp-2)/b0_t(ntemp-2)
WRITE(stdout,'(5x,"delta_S =",f14.6)') -dbsdt/beta_t(ntemp-2)/b0_s(ntemp-2)

RETURN
END SUBROUTINE summarize_anhar
