!
! Copyright (C) 2013-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_anharmonic()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE temperature,    ONLY : ntemp, temp
USE thermodynamics, ONLY : ph_cv
USE anharmonic,     ONLY : alpha_t, beta_t, gamma_t, cp_t, cv_t, b0_s, &
                           vmin_t, b0_t, b01_t
USE thermo_mod,     ONLY : omega_geo
USE control_pressure, ONLY : pressure, pressure_kb
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char 

IF (my_image_id /= root_image) RETURN

CALL compute_beta(vmin_t, beta_t, temp, ntemp)

alpha_t = beta_t / 3.0_DP

CALL compute_cp(beta_t, vmin_t, b0_t, ph_cv, cv_t, cp_t, b0_s, gamma_t)

IF (ionode) THEN
!
!   here we plot the quantities calculated from the phonon dos
!
   filename="anhar_files/"//TRIM(flanhar)
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_murn_beta(temp, vmin_t, b0_t, b01_t, beta_t, ntemp, filename)
!
!   here auxiliary quantities calculated from the phonon dos
!
   filename="anhar_files/"//TRIM(flanhar)//'.aux'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_aux_anharm(temp, gamma_t, cv_t, cp_t, b0_t, b0_s, ntemp, filename)

END IF

RETURN
END SUBROUTINE write_anharmonic

SUBROUTINE write_ph_freq_anharmonic()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE temperature,    ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : phf_cv
USE ph_freq_anharmonic, ONLY : alphaf_t, betaf_t, gammaf_t, cpf_t, cvf_t, &
                        b0f_s, vminf_t, b0f_t, b01f_t
USE control_pressure, ONLY : pressure, pressure_kb
USE thermo_mod,     ONLY : omega_geo
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char 

IF (my_image_id /= root_image) RETURN

CALL compute_beta(vminf_t, betaf_t, temp, ntemp)

alphaf_t = betaf_t / 3.0_DP

CALL compute_cp(betaf_t, vminf_t, b0f_t, phf_cv, cvf_t, cpf_t, b0f_s, gammaf_t)

IF (ionode) THEN
!
!   here we plot the quantities calculated from the phonon dos
!
   filename="anhar_files/"//TRIM(flanhar)//'_ph'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_murn_beta(temp, vminf_t, b0f_t, b01f_t, betaf_t, ntemp, filename)

!
!   here auxiliary quantities calculated from the phonon dos
!
   filename="anhar_files/"//TRIM(flanhar)//'.aux_ph'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_aux_anharm(temp, gammaf_t, cvf_t, cpf_t, b0f_t, b0f_s, &
                                                          ntemp, filename)
END IF

RETURN
END SUBROUTINE write_ph_freq_anharmonic

SUBROUTINE write_grun_anharmonic()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : ngeo, omega_geo
USE temperature,    ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : ph_freq_save, phf_cv
USE ph_freq_anharmonic,     ONLY :  vminf_t, cvf_t, b0f_t, cpf_t, b0f_s
USE grun_anharmonic, ONLY : betab, grun_gamma_t, poly_grun, poly_order
USE ph_freq_module, ONLY : thermal_expansion_ph, ph_freq_type,  &
                           destroy_ph_freq, init_ph_freq
USE control_pressure, ONLY : pressure, pressure_kb
USE control_grun,     ONLY : lv0_t, lb0_t
USE control_mur,    ONLY : vmin, b0
USE ifc,            ONLY : nq1_d, nq2_d, nq3_d
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char 
INTEGER :: itemp, iu_therm, i, nq, imode, iq
TYPE(ph_freq_type) :: ph_freq    ! the frequencies at the volumes at
                                 ! which the gruneisen parameters are 
                                 ! calculated
TYPE(ph_freq_type) :: ph_grun    ! the gruneisen parameters recomputed
                                 ! at each temperature at the volume
                                 ! corresponding to that temperature
REAL(DP) :: vm

IF (my_image_id /= root_image) RETURN
!
!  compute thermal expansion from gruneisen parameters. 
!  NB: betab is multiplied by the bulk modulus
!
nq=ph_freq_save(1)%nq
CALL init_ph_freq(ph_grun, nat, nq1_d, nq2_d, nq3_d, nq, .FALSE.)
CALL init_ph_freq(ph_freq, nat, nq1_d, nq2_d, nq3_d, nq, .FALSE.)
DO itemp = 1, ntemp
   IF (lv0_t) THEN
      vm=vminf_t(itemp)
   ELSE
      vm=vmin
   ENDIF
   ph_freq%nu= 0.0_DP
   ph_freq%wg=ph_freq_save(1)%wg
   ph_grun%nu= 0.0_DP
   DO iq=1,nq
      DO imode=1,3*nat
         DO i=1,poly_order
            ph_freq%nu(imode,iq) = ph_freq%nu(imode,iq) + &
                  poly_grun(i,imode,iq) * vm**(i-1)
            ph_grun%nu(imode,iq) = ph_grun%nu(imode,iq) - &
                  poly_grun(i,imode,iq) * vm**(i-2) * (i-1.0_DP)
         END DO
         IF (ph_freq%nu(imode,iq) > 0.0_DP ) THEN
             ph_grun%nu(imode,iq) = ph_grun%nu(imode,iq) / &
                                    ph_freq%nu(imode,iq)
         ELSE
            ph_grun%nu(imode,iq) = 0.0_DP
         ENDIF
      END DO
   END DO
      
   CALL thermal_expansion_ph(ph_freq, ph_grun, temp(itemp), betab(itemp))
   IF (lb0_t) THEN
      betab(itemp)=betab(itemp) * ry_kbar / b0f_t(itemp)
   ELSE
      betab(itemp)=betab(itemp) * ry_kbar / b0
   ENDIF
END DO

CALL compute_cp(betab, vminf_t, b0f_t, phf_cv, cvf_t, cpf_t, b0f_s, &
                                                             grun_gamma_t)
IF (ionode) THEN
!
!   here quantities calculated from the gruneisen parameters
!
   filename="anhar_files/"//TRIM(flanhar)//'.aux_grun'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
   WRITE(iu_therm,'("#   T (K)       beta(T)      gamma(T)    C_v (T) ( Ry / N/ K)  " )' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e13.6,3e16.8)') temp(itemp), &
              betab(itemp)*1.D6, grun_gamma_t(itemp), cvf_t(itemp)
   END DO
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
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ntemp
  REAL(DP), INTENT(IN) :: vmin_t(ntemp), temp(ntemp)
  REAL(DP), INTENT(OUT) :: beta_t(ntemp) 

  INTEGER :: itemp

  beta_t=0.0_DP
!
!  just interpolate linearly
!
  DO itemp = 2, ntemp-1
     beta_t(itemp) = (vmin_t(itemp+1)-vmin_t(itemp-1)) / &
                     (temp(itemp+1)-temp(itemp-1)) / vmin_t(itemp)
  END DO

  RETURN
END SUBROUTINE compute_beta

SUBROUTINE write_murn_beta(temp, vmin, b0, b01, beta, ntemp, filename)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), vmin(ntemp), b0(ntemp), b01(ntemp), &
                        beta(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm

iu_therm=2
OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

WRITE(iu_therm,'("# beta is the volume thermal expansion ")')
WRITE(iu_therm,'("#   T (K)     V(T) (a.u.)^3   B (T) (kbar) &
                   & d B (T) / dP  beta (10^(-6) K^(-1))")' )

DO itemp = 2, ntemp-1
   WRITE(iu_therm, '(e12.5,e20.13,2e14.6,e18.8)') temp(itemp), &
                vmin(itemp), b0(itemp), b01(itemp), beta(itemp)*1.D6
END DO
  
CLOSE(iu_therm)
RETURN
END SUBROUTINE write_murn_beta

SUBROUTINE write_aux_anharm(temp, gamma_t, cv, cp, b0_t, b0_s, ntemp, filename)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), gamma_t(ntemp), cv(ntemp), cp(ntemp), &
                  b0_s(ntemp), b0_t(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm

iu_therm=2
OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
WRITE(iu_therm,'("#   T (K)       gamma(T)     C_v ( Ry / cell / K ) &
                 &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )

DO itemp = 2, ntemp-1
   WRITE(iu_therm, '(5e16.8)') temp(itemp),               &
                               gamma_t(itemp), cv(itemp), &
                               cp(itemp) - cv(itemp),   &
                               b0_s(itemp) - b0_t(itemp)
END DO
CLOSE(iu_therm)

RETURN
END SUBROUTINE write_aux_anharm
