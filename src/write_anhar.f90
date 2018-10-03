!
! Copyright (C) 2013-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_anharmonic()
!
!   This routine writes on output the anharmonic quantities calculated
!   at the volume that minimize the free energy (computed from phonon dos)
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp
USE thermodynamics, ONLY : ph_cv, ph_b_fact
USE anharmonic,     ONLY : alpha_t, beta_t, gamma_t, cp_t, cv_t, b0_s, &
                           vmin_t, free_e_min_t, b0_t, b01_t, bfact_t
USE data_files,     ONLY : flanhar
USE control_thermo, ONLY : with_eigen

IMPLICIT NONE
CHARACTER(LEN=256) :: filename

CALL compute_beta(vmin_t, beta_t, temp, ntemp)

alpha_t = beta_t / 3.0_DP

CALL interpolate_cv(vmin_t, ph_cv, cv_t) 
CALL compute_cp(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)
!
!   here we plot the quantities calculated from the phonon dos
!
filename="anhar_files/"//TRIM(flanhar)
CALL add_pressure(filename)

CALL write_ener_beta(temp, vmin_t, free_e_min_t, beta_t, ntemp, filename)
!
!   here the bulk modulus and the gruneisen parameter
!
filename="anhar_files/"//TRIM(flanhar)//'.bulk_mod'
CALL add_pressure(filename)

CALL write_bulk_anharm(temp, gamma_t, b0_t, b0_s, ntemp, filename)
!
!   here the derivative of the bulk modulus with respect to pressure
!
filename="anhar_files/"//TRIM(flanhar)//'.dbulk_mod'
CALL add_pressure(filename)

CALL write_dbulk_anharm(temp, b01_t, ntemp, filename)
!
!   here the heat capacities (constant strain and constant stress)
!
filename="anhar_files/"//TRIM(flanhar)//'.heat'
CALL add_pressure(filename)

CALL write_heat_anharm(temp, cv_t, cp_t, ntemp, filename)

!
!   here the b factors
!
IF (with_eigen) THEN
   CALL interpolate_b_fact(vmin_t, ph_b_fact, bfact_t)
   filename="anhar_files/"//TRIM(flanhar)
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfact_t, ntemp, filename)
END IF


RETURN
END SUBROUTINE write_anharmonic

SUBROUTINE write_ph_freq_anharmonic()
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp
USE control_thermo, ONLY : with_eigen
USE ph_freq_thermodynamics, ONLY : phf_cv, phf_b_fact
USE ph_freq_anharmonic, ONLY : alphaf_t, betaf_t, gammaf_t, cpf_t, cvf_t, &
                        b0f_s, free_e_minf_t, vminf_t, b0f_t, b01f_t, &
                        bfactf_t
USE data_files,     ONLY : flanhar

IMPLICIT NONE
CHARACTER(LEN=256) :: filename

CALL compute_beta(vminf_t, betaf_t, temp, ntemp)

alphaf_t = betaf_t / 3.0_DP

CALL interpolate_cv(vminf_t, phf_cv, cvf_t) 
CALL compute_cp(betaf_t, vminf_t, b0f_t, cvf_t, cpf_t, b0f_s, gammaf_t)

!
!   here we plot the quantities calculated from the phonon dos
!
filename="anhar_files/"//TRIM(flanhar)//'_ph'
CALL add_pressure(filename)

CALL write_ener_beta(temp, vminf_t, free_e_minf_t, betaf_t, ntemp, filename)
!
!   here the bulk modulus and the gruneisen parameter
!
filename="anhar_files/"//TRIM(flanhar)//'.bulk_mod_ph'
CALL add_pressure(filename)

CALL write_bulk_anharm(temp, gammaf_t, b0f_t, b0f_s, ntemp, filename)
!
!   here the derivative of the bulk modulus with respect to pressure
!
filename="anhar_files/"//TRIM(flanhar)//'.dbulk_mod_ph'
CALL add_pressure(filename)

CALL write_dbulk_anharm(temp, b01f_t, ntemp, filename)
!
!   here the heat capacities (constant strain and constant stress)
!
filename="anhar_files/"//TRIM(flanhar)//'.heat_ph'
CALL add_pressure(filename)

CALL write_heat_anharm(temp, cvf_t, cpf_t, ntemp, filename)
!
!   here the b factors
!
IF (with_eigen) THEN
   CALL interpolate_b_fact(vminf_t, phf_b_fact, bfactf_t)
   filename="anhar_files/"//TRIM(flanhar)//'_ph'
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfactf_t, ntemp, filename)
END IF


RETURN
END SUBROUTINE write_ph_freq_anharmonic

SUBROUTINE write_grun_anharmonic()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE ions_base,      ONLY : nat
USE temperature,    ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : ph_freq_save, phf_cv
USE anharmonic,     ONLY :  vmin_t, b0_t, cv_t
USE ph_freq_anharmonic,     ONLY :  vminf_t, cvf_t, b0f_t, cpf_t, b0f_s
USE grun_anharmonic, ONLY : betab, cp_grun_t, b0_grun_s, &
                            grun_gamma_t, poly_grun, poly_order
USE ph_freq_module, ONLY : thermal_expansion_ph, ph_freq_type,  &
                           destroy_ph_freq, init_ph_freq
USE control_grun,     ONLY : lv0_t, lb0_t
USE control_mur,    ONLY : vmin, b0
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_dosq,   ONLY : nq1_d, nq2_d, nq3_d
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm, i, nq, imode, iq, startq, lastq, iq_eff, nq_eff
TYPE(ph_freq_type) :: ph_freq    ! the frequencies at the volumes at
                                 ! which the gruneisen parameters are 
                                 ! calculated
TYPE(ph_freq_type) :: ph_grun    ! the gruneisen parameters recomputed
                                 ! at each temperature at the volume
                                 ! corresponding to that temperature
REAL(DP) :: vm, f, g
INTEGER :: find_free_unit

!
!  compute thermal expansion from gruneisen parameters. 
!  NB: betab is multiplied by the bulk modulus
!
nq=ph_freq_save(1)%nq
startq=ph_freq_save(1)%startq
lastq=ph_freq_save(1)%lastq
nq_eff = ph_freq_save(1)%nq_eff
CALL init_ph_freq(ph_grun, nat, nq1_d, nq2_d, nq3_d, nq_eff, startq, &
                                                             lastq, nq, .FALSE.)
CALL init_ph_freq(ph_freq, nat, nq1_d, nq2_d, nq3_d, nq_eff, startq, &
                                                             lastq, nq, .FALSE.)
DO iq=1, nq_eff
   ph_freq%wg(iq)=ph_freq_save(1)%wg(iq)
ENDDO

betab=0.0_DP
DO itemp = 1, ntemp
   IF (lv0_t) THEN
      IF (ltherm_freq) THEN
         vm=vminf_t(itemp)
      ELSEIF (ltherm_dos) THEN
         vm=vmin_t(itemp)
      ELSE
         vm=vmin
      ENDIF
   ELSE
      vm=vmin
   ENDIF
   ph_freq%nu= 0.0_DP
   ph_grun%nu= 0.0_DP
   iq_eff=0
   DO iq=startq, lastq
      iq_eff=iq_eff+1
      DO imode=1,3*nat
         f=poly_grun(poly_order,imode,iq)
         g=poly_grun(poly_order,imode,iq)*(poly_order-1.0_DP)
         DO i=poly_order-1,1,-1
            f = poly_grun(i,imode,iq) + f*vm
            g = poly_grun(i,imode,iq)*(i-1.0_DP) + g*vm
         END DO
!
!     g here is V d w / d V 
!
         ph_freq%nu(imode,iq_eff)=f
         IF ( f > 0.0_DP) THEN 
            ph_grun%nu(imode,iq_eff)=-g/vm/f
         ELSE
            ph_grun%nu(imode,iq_eff) = 0.0_DP
         ENDIF
      ENDDO
   ENDDO
   CALL thermal_expansion_ph(ph_freq, ph_grun, temp(itemp), betab(itemp))
   IF (lb0_t) THEN
      IF (ltherm_freq) THEN
         betab(itemp)=betab(itemp) * ry_kbar / b0f_t(itemp)
      ELSEIF(ltherm_dos) THEN 
         betab(itemp)=betab(itemp) * ry_kbar / b0_t(itemp)
      ELSE
         betab(itemp)=betab(itemp) * ry_kbar / b0
      ENDIF
   ELSE
      betab(itemp)=betab(itemp) * ry_kbar / b0
   ENDIF
END DO
CALL mp_sum(betab, world_comm)
!
IF (ltherm_freq) THEN
   CALL compute_cp(betab, vminf_t, b0f_t, cvf_t, cp_grun_t, &
                                                   b0_grun_s, grun_gamma_t)
ELSEIF (ltherm_dos) THEN
   CALL compute_cp(betab, vmin_t, b0_t, cv_t, cp_grun_t, &
                                                   b0_grun_s, grun_gamma_t)
ELSE
   vmin_t(1:ntemp)=vmin
   b0_t(1:ntemp)=b0
   CALL interpolate_cv(vmin_t, phf_cv, cv_t)
   CALL compute_cp(betab, vmin_t, b0_t, cv_t, cp_grun_t, &
                                              b0_grun_s, grun_gamma_t)
ENDIF

IF (meta_ionode) THEN
!
!   here quantities calculated from the gruneisen parameters
!
   filename="anhar_files/"//TRIM(flanhar)//'.aux_grun'
   CALL add_pressure(filename)

   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
   WRITE(iu_therm,'("#   T (K)   beta(T)x10^6    gamma(T)      &
                 &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )

   IF (ltherm_freq) THEN
      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(5e16.8)') temp(itemp), betab(itemp)*1.D6,    &
                        grun_gamma_t(itemp), cp_grun_t(itemp)-cvf_t(itemp),   &
                        b0_grun_s(itemp) - b0f_t(itemp)
      END DO
   ELSE
      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(5e16.8)') temp(itemp), betab(itemp)*1.D6,    &
                        grun_gamma_t(itemp), cp_grun_t(itemp) - cv_t(itemp), &
                        b0_grun_s(itemp) - b0_t(itemp)
      END DO
   ENDIF
   CLOSE(iu_therm)
END IF

CALL destroy_ph_freq(ph_freq)
CALL destroy_ph_freq(ph_grun)

RETURN
END SUBROUTINE write_grun_anharmonic

SUBROUTINE compute_beta(vmin_t, beta_t, temp, ntemp)
!
!  This routine receives as input the volume for ntemp temperatures
!  and computes the volume thermal expansion for ntemp-2 temperatures.
!  In the first and last point the thermal expansion is not computed.
!
USE kinds, ONLY : DP
USE mp,    ONLY : mp_sum
USE mp_world, ONLY : world_comm

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ntemp
  REAL(DP), INTENT(IN) :: vmin_t(ntemp), temp(ntemp)
  REAL(DP), INTENT(OUT) :: beta_t(ntemp) 

  INTEGER :: itemp, startt, lastt

  beta_t=0.0_DP
  CALL divide(world_comm, ntemp, startt, lastt)
!
!  just interpolate linearly
!
  DO itemp = max(2,startt), min(ntemp-1, lastt)
     beta_t(itemp) = (vmin_t(itemp+1)-vmin_t(itemp-1)) / &
                     (temp(itemp+1)-temp(itemp-1)) / vmin_t(itemp)
  END DO
  CALL mp_sum(beta_t, world_comm)

  RETURN
END SUBROUTINE compute_beta

SUBROUTINE write_ener_beta(temp, vmin, emin, beta, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), emin(ntemp), vmin(ntemp), beta(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# beta is the volume thermal expansion ")')
   WRITE(iu_therm,'("#   T (K)         V(T) (a.u.)^3          F (T) (Ry) &
                   &      beta (10^(-6) K^(-1))")' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,2e23.13,e18.8)') temp(itemp), &
                vmin(itemp), emin(itemp), beta(itemp)*1.D6
   END DO
  
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_ener_beta

SUBROUTINE write_bulk_anharm(temp, gammat, b0t, b0s, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), gammat(ntemp), b0t(ntemp), b0s(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("#  ")')
   WRITE(iu_therm,'("#   T (K)          gamma(T)            B_T(T) (kbar)     &
                                 & B_S(T)-B_T(T) (kbar) ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e22.13)') temp(itemp), &
                 gammat(itemp), b0t(itemp), b0s(itemp)-b0t(itemp)
   END DO
  
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_bulk_anharm

SUBROUTINE write_dbulk_anharm(temp, b01t, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), b01t(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("#  ")')
   WRITE(iu_therm,'("#   T (K)          dB/dp (T) ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e23.13)') temp(itemp), b01t(itemp)
   END DO
  
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_dbulk_anharm

SUBROUTINE write_heat_anharm(temp, cvt, cpt, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), cvt(ntemp), cpt(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# ")')
   WRITE(iu_therm,'("# T (K)        C_V(T) (Ry/cell/K)    C_P(T) (Ry/cell/K) &
                                 &  C_P-C_V(T) (Ry/cell/K)")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e22.13)') temp(itemp), &
                 cvt(itemp), cpt(itemp), cpt(itemp)-cvt(itemp)
   END DO
  
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_heat_anharm
!
! Copyright (C) 2018 Cristiano Malica
!
SUBROUTINE write_anharm_bfact(temp, bfact_t, ntemp, filename)

USE ions_base, ONLY : nat
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp 
REAL(DP), INTENT(IN) :: temp(ntemp), bfact_t(6,nat,ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: na, ijpol
INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

CHARACTER(LEN=6) :: int_to_char

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   
   DO na=1, nat

      OPEN(UNIT=iu_therm, FILE=TRIM(filename)//'.'//&
                                 TRIM(int_to_char(na))//'.dw', STATUS='UNKNOWN', &
                                                      FORM='FORMATTED')
      WRITE(iu_therm,'("# ")')

      WRITE(iu_therm,'(6x,"T",14x,"B_11",14x,"B_12",14x,"B_13",14x,"B_22",&
                           &14x,"B_23",14x,"B_33")')

      DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,6e18.8)') temp(itemp), (bfact_t(ijpol,na,itemp), ijpol=1,6)
      END DO

      CLOSE(UNIT=iu_therm, STATUS='KEEP')
   END DO

ENDIF
RETURN
END SUBROUTINE write_anharm_bfact


