!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_anharmonic()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo
USE temperature,    ONLY : ntemp, temp
USE thermodynamics, ONLY : ph_cv
USE anharmonic,     ONLY : alpha_t, beta_t, gamma_t, cp_t, cv_t, b0_s, &
                           vmin_t, b0_t, b01_t
USE control_thermo, ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm

IF (my_image_id /= root_image) RETURN

DO itemp = 2, ntemp-1
   beta_t(itemp) = (vmin_t(itemp+1)-vmin_t(itemp-1)) / &
                   (temp(itemp+1)-temp(itemp-1)) / vmin_t(itemp)
END DO

alpha_t = beta_t / 3.0_DP

CALL compute_cp(beta_t, vmin_t, b0_t, ph_cv, cv_t, cp_t, b0_s, gamma_t)

IF (ionode) THEN
!
!   here we plot the quantities calculated from the phonon dos
!
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(flanhar), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# beta is the volume thermal expansion ")')
   WRITE(iu_therm,'("#   T (K)     V(T) (a.u.)^3   B (T) (kbar) &
                      & d B (T) / dP     beta (10^(-6))")' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e20.13,2e14.6,e18.8)') temp(itemp), &
                   vmin_t(itemp), b0_t(itemp), b01_t(itemp), beta_t(itemp)*1.D6
   END DO
   CLOSE(iu_therm)
!
!   here auxiliary quantities calculated from the phonon dos
!
   filename=TRIM(flanhar)//'.aux'
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
   WRITE(iu_therm,'("#   T (K)       gamma(T)       C_v ( Ry / cell ) &
                    &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(5e16.8)') temp(itemp),                  &
                                  gamma_t(itemp), cv_t(itemp),  &
                                  cp_t(itemp) - cv_t(itemp),    &
                                  b0_s(itemp) - b0_t(itemp)
   END DO
   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_anharmonic

SUBROUTINE write_ph_freq_anharmonic()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo
USE temperature,    ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : phf_cv
USE ph_freq_anharmonic, ONLY : alphaf_t, betaf_t, gammaf_t, cpf_t, cvf_t, &
                        b0f_s, vminf_t, b0f_t, b01f_t
USE control_thermo, ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm

IF (my_image_id /= root_image) RETURN

DO itemp = 2, ntemp-1
   betaf_t(itemp) = (vminf_t(itemp+1)-vminf_t(itemp-1)) / &
                    (temp(itemp+1)-temp(itemp-1)) / vminf_t(itemp)
END DO
alphaf_t = betaf_t / 3.0_DP

CALL compute_cp(betaf_t, vminf_t, b0f_t, phf_cv, cvf_t, cpf_t, b0f_s, gammaf_t)

IF (ionode) THEN
!
!   here we plot the quantities calculated from the phonon dos
!
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(flanhar)//'_ph', STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# beta is the volume thermal expansion ")')
   WRITE(iu_therm,'("#   T (K)     V(T) (a.u.)^3   B (T) (kbar) &
                      & d B (T) / dP     beta (10^(-6))")' )


   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e20.13,2e14.6,e18.8)') temp(itemp), &
              vminf_t(itemp), b0f_t(itemp), b01f_t(itemp), betaf_t(itemp)*1.D6
   END DO
   CLOSE(iu_therm)
!
!   here auxiliary quantities calculated from the phonon dos
!
   filename=TRIM(flanhar)//'.aux_ph'
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
   WRITE(iu_therm,'("#   T (K)       gamma(T)       C_p ( Ry / cell ) &
                    &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(5e16.8)') temp(itemp),               &
                                  gammaf_t(itemp), cvf_t(itemp), &
                                  cpf_t(itemp) - cvf_t(itemp),   &
                                  b0f_s(itemp) - b0f_t(itemp)
   END DO
   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_ph_freq_anharmonic

SUBROUTINE write_grun_anharmonic()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo
USE temperature,    ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : ph_freq_save, phf_cv
USE ph_freq_anharmonic,     ONLY :  vminf_t, cvf_t, b0f_t
USE grun_anharmonic, ONLY : ph_grun, betab
USE ph_freq_module, ONLY : thermal_expansion_ph
USE control_thermo, ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm

IF (my_image_id /= root_image) RETURN

!
!  compute thermal expansion from gruneisen parameters. 
!  NB: betab is multiplied by the bulk modulus
!
DO itemp = 1, ntemp
   CALL thermal_expansion_ph(ph_freq_save(ngeo/2+1), ph_grun, &
                             temp(itemp), betab(itemp))
END DO

IF (ionode) THEN
!
!   here quantities calculated from the gruneisen parameters
!
   filename=TRIM(flanhar)//'.aux_grun'
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
   WRITE(iu_therm,'("#   T (K)       beta(T)      gamma(T)    beta(T) (B_0)  C_v (T) ( Ry / N/ K)  " )' )

   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e13.6,4e16.8)') temp(itemp), &
              betab(itemp)*1.D6*ry_kbar/b0f_t(itemp), &
              betab(itemp)*omega_geo(ngeo/2+1)/phf_cv(itemp,ngeo/2+1), &
              betab(itemp)*1.D6*ry_kbar*vminf_t(1)/vminf_t(itemp)/b0f_t(1), &
              cvf_t(itemp)
   END DO
   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_grun_anharmonic

