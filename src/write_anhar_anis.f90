!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_anhar_anis()
!
!   This routine computes the thermal expansion tensor for anisotropic
!   solids.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE temperature,    ONLY : ntemp, temp
USE thermodynamics, ONLY : ph_cv
USE anharmonic,     ONLY : alpha_anis_t, vmin_t, b0_t, celldm_t, beta_t
USE control_quadratic_energy, ONLY : nvar, degree, coeff_t
USE control_elastic_constants, ONLY : ibrav_save
USE cell_base,      ONLY : ibrav
USE control_thermo, ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm
REAL(DP) :: compute_omega_geo
REAL(DP) :: fact1, fact2, deriv1, deriv2

IF (my_image_id /= root_image) RETURN

alpha_anis_t=0.0_DP
SELECT CASE (ibrav) 
   CASE(1,2,3) 
      DO itemp = 2, ntemp-1
         alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
      END DO
   CASE(4,6,7)
      DO itemp = 2, ntemp-1
         alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
         alpha_anis_t(2,itemp) = ( celldm_t(3,itemp+1)*celldm_t(1,itemp+1)-   &
                                   celldm_t(3,itemp-1)*celldm_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldm_t(3,itemp)*  &
                          celldm_t(1,itemp)) 
      END DO
   CASE(5)
      DO itemp = 2, ntemp-1
         fact1 = 2.0_DP * ( 1.0_DP - celldm_t(4,itemp) )
         fact2 = 1.0_DP + 2.0_DP * celldm_t(4,itemp)
         deriv1 = ( celldm_t(1,itemp+1) - celldm_t(1,itemp-1) )/      &
                                     ( temp(itemp+1) - temp(itemp-1) )
         deriv2 = ( celldm_t(4,itemp+1) - celldm_t(4,itemp-1)) /      &
                                     ( temp(itemp+1) - temp(itemp-1) )
         alpha_anis_t(1,itemp) = deriv1 / celldm_t(1, itemp) - deriv2 / fact1 
         alpha_anis_t(2,itemp) = deriv1 / celldm_t(2, itemp) + deriv2 / fact2
      END DO
   CASE (8,9,10,11)
      DO itemp = 2, ntemp-1
         alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
         alpha_anis_t(2,itemp) = ( celldm_t(2,itemp+1)*celldm_t(1,itemp+1)-   &
                                   celldm_t(2,itemp-1)*celldm_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldm_t(2,itemp)*  &
                          celldm_t(1,itemp)) 
         alpha_anis_t(3,itemp) = ( celldm_t(3,itemp+1)*celldm_t(1,itemp+1)-   &
                                   celldm_t(3,itemp-1)*celldm_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldm_t(3,itemp)*  &
                          celldm_t(1,itemp)) 
      END DO
   CASE DEFAULT
!
!   In this case do nothing. The thermal expansion tensor is not written later
!
END SELECT
!
!  compute the volume as a function of temperature
!
DO itemp=1,ntemp
   vmin_t(itemp)=compute_omega_geo(ibrav_save, celldm_t(1,itemp))
ENDDO
!
!  compute the volume thermal expansion as the derivative of the volume
!  with respect to temperature
!
CALL compute_beta(vmin_t, beta_t, temp, ntemp)

!CALL compute_cp(beta_t, vmin_t, b0_t, ph_cv, cv_t, cp_t, b0_s, gamma_t)

IF (ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(flanhar), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# alpha is the linear thermal expansion ")')
   WRITE(iu_therm,'("#   T (K)        V(T) (a.u.)^3           beta (x10^6)  ")' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e20.13,e15.6)') temp(itemp), vmin_t(itemp), &
                                                           beta_t(itemp)*1.D6
   END DO
   CLOSE(iu_therm)
!
!   here auxiliary quantities calculated from the phonon dos
!
!   filename=TRIM(flanhar)//'.aux'
!   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
!   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
!   WRITE(iu_therm,'("#   T (K)       gamma(T)       C_v ( Ry / cell ) &
!                    &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )
!
!   DO itemp = 2, ntemp-1
!      WRITE(iu_therm, '(5e16.8)') temp(itemp),                  &
!                                  gamma_t(itemp), cv_t(itemp),  &
!                                  cp_t(itemp) - cv_t(itemp),    &
!                                  b0_s(itemp) - b0_t(itemp)
!   END DO
!   CLOSE(iu_therm)
!END IF
!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename=TRIM(flanhar)//'.celldm'
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (ibrav==1 .OR. ibrav==2 .OR. ibrav==3 ) THEN
      WRITE(iu_therm,'("#   T (K)      celldm(1)      alpha_xx(x10^6)")' )
      DO itemp = 1, ntemp-1
         WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldm_t(1,itemp), &
                                               alpha_anis_t(1,itemp)*1.D6
      END DO
   ELSEIF (ibrav==4 .OR. ibrav==6 .OR. ibrav==7 ) THEN
      WRITE(iu_therm,'("#   T (K)   celldm(1)   celldm(3)    alpha_xx(x10^6)   alpha_zz (x10^6")' )
      DO itemp = 1, ntemp-1
         WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldm_t(1,itemp), &
                                                        celldm_t(3,itemp), &
                                               alpha_anis_t(1,itemp)*1.D6, &
                                               alpha_anis_t(2,itemp)*1.D6
      END DO
   ELSEIF ( ibrav==5 ) THEN
      WRITE(iu_therm,'("#   T (K)   celldm(1)   celldm(4)    alpha_xx(x10^6)   alpha_zz (x10^6")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldm_t(1,itemp), &
                                                        celldm_t(4,itemp), &
                                               alpha_anis_t(1,itemp)*1.D6, &
                                               alpha_anis_t(2,itemp)*1.D6
      END DO
   ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)        celldm(2)        celldm(3)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,3e20.9)') temp(itemp), celldm_t(1,itemp), &
                                                        celldm_t(2,itemp), &
                                                        celldm_t(3,itemp)
      END DO
   ELSEIF (ibrav==12 .OR. ibrav==13) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)         celldm(2)        celldm(3)        celldm(4)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,4e17.9)') temp(itemp), celldm_t(1,itemp), &
                                                        celldm_t(2,itemp), &
                                                        celldm_t(3,itemp), &
                                                        celldm_t(4,itemp)
      END DO
   ELSEIF (ibrav==-12 .OR. ibrav==-13) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)         celldm(2)        celldm(3)        celldm(5)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,4e17.9)') temp(itemp), celldm_t(1,itemp), &
                                                        celldm_t(2,itemp), &
                                                        celldm_t(3,itemp), &
                                                        celldm_t(5,itemp)
      END DO
   ELSEIF (ibrav==14) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)         celldm(2)        &
                  &celldm(3)        celldm(4)        celldm(5)        celldm(6)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,6e15.7)') temp(itemp), celldm_t(1,itemp), &
                                                        celldm_t(2,itemp), &
                                                        celldm_t(3,itemp), &
                                                        celldm_t(4,itemp), &
                                                        celldm_t(5,itemp), &
                                                        celldm_t(6,itemp)
      END DO
   ELSE IF (ibrav==0) THEN
!
!  In this case we write nothing but do not stop
!
   ELSE
      CALL errore('write_anhar_anis','ibrav not programmed',1)
   END IF
   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_anhar_anis
!
SUBROUTINE write_ph_freq_anhar_anis()
!
!   This routine computes the thermal expansion tensor for anisotropic
!   solids using the free energy calculated from direct integration of
!   the phonon frequencies.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE temperature,    ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : phf_cv
USE ph_freq_anharmonic, ONLY : alphaf_anis_t, vminf_t, b0f_t, celldmf_t, &
                               betaf_t
USE control_quadratic_energy, ONLY : nvar, degree, coeff_t
USE control_elastic_constants, ONLY : ibrav_save
USE cell_base,      ONLY : ibrav
USE control_thermo, ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm
REAL(DP) :: compute_omega_geo
REAL(DP) :: fact1, fact2, deriv1, deriv2

IF (my_image_id /= root_image) RETURN

alphaf_anis_t=0.0_DP
SELECT CASE (ibrav) 
   CASE(1,2,3) 
      DO itemp = 2, ntemp-1
         alphaf_anis_t(1,itemp)=(celldmf_t(1,itemp+1)-celldmf_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldmf_t(1,itemp)
      END DO
   CASE(4,6,7)
      DO itemp = 2, ntemp-1
         alphaf_anis_t(1,itemp)=(celldmf_t(1,itemp+1)-celldmf_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldmf_t(1,itemp)
         alphaf_anis_t(2,itemp)=( celldmf_t(3,itemp+1)*celldmf_t(1,itemp+1)-   &
                                celldmf_t(3,itemp-1)*celldmf_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldmf_t(3,itemp)*  &
                          celldmf_t(1,itemp)) 
      END DO
   CASE(5)
      DO itemp = 2, ntemp-1
         fact1 = 2.0_DP * ( 1.0_DP - celldmf_t(4,itemp) )
         fact2 = 1.0_DP + 2.0_DP * celldmf_t(4,itemp)
         deriv1 = ( celldmf_t(1,itemp+1) - celldmf_t(1,itemp-1) )/      &
                                     ( temp(itemp+1) - temp(itemp-1) )
         deriv2 = ( celldmf_t(4,itemp+1) - celldmf_t(4,itemp-1)) /      &
                                     ( temp(itemp+1) - temp(itemp-1) )
         alphaf_anis_t(1,itemp) = deriv1 / celldmf_t(1, itemp) - deriv2 / fact1 
         alphaf_anis_t(2,itemp) = deriv1 / celldmf_t(2, itemp) + deriv2 / fact2
      END DO
   CASE (8,9,10,11)
      DO itemp = 2, ntemp-1
         alphaf_anis_t(1,itemp)=(celldmf_t(1,itemp+1)-celldmf_t(1,itemp-1)) / &
                         (temp(itemp+1)-temp(itemp-1)) / celldmf_t(1,itemp)
         alphaf_anis_t(2,itemp)=(celldmf_t(2,itemp+1)*celldmf_t(1,itemp+1)-   &
                                celldmf_t(2,itemp-1)*celldmf_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldmf_t(2,itemp)*  &
                          celldmf_t(1,itemp)) 
         alphaf_anis_t(3,itemp)=(celldmf_t(3,itemp+1)*celldmf_t(1,itemp+1)-   &
                                 celldmf_t(3,itemp-1)*celldmf_t(1,itemp-1) )/ &
                         (temp(itemp+1)-temp(itemp-1))/(celldmf_t(3,itemp)*  &
                          celldmf_t(1,itemp)) 
      END DO
   CASE DEFAULT
!
!   In this case do nothing. The thermal expansion tensor is not written later
!
END SELECT
!
!  compute the volume as a function of temperature
!
DO itemp=1,ntemp
   vminf_t(itemp)=compute_omega_geo(ibrav_save, celldmf_t(1,itemp))
ENDDO
!
!  compute the volume thermal expansion as the derivative of the volume
!  with respect to temperature
!
CALL compute_beta(vminf_t, betaf_t, temp, ntemp)

!DO itemp = 1, ntemp 
!   CALL compute_bulk_modulus(ibrav, celldm_t(1,itemp), coeff_t(1,itemp), &
!                                 nvar, degree, b0_t(itemp))
!   CALL compute_elcons_t(ibrav, celldm_t(1,itemp), coeff_t(1,itemp), &
!                                 nvar, degree, elcons_t(1,itemp))
!ENDDO

!CALL compute_cp(beta_t, vmin_t, b0_t, ph_cv, cv_t, cp_t, b0_s, gamma_t)

IF (ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(flanhar)//'_ph', STATUS='UNKNOWN', &
                                                  FORM='FORMATTED')
   WRITE(iu_therm,'("# alpha is the linear thermal expansion ")')
   WRITE(iu_therm,'("#   T (K)        V(T) (a.u.)^3           beta (x10^6)  ")' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e20.13,e15.6)') temp(itemp), vminf_t(itemp), &
                                              betaf_t(itemp)*1.D6
   END DO
   CLOSE(iu_therm)
!
!   here auxiliary quantities calculated from the phonon dos
!
!   filename=TRIM(flanhar)//'.aux'
!   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
!   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
!   WRITE(iu_therm,'("#   T (K)       gamma(T)       C_v ( Ry / cell ) &
!                    &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )
!
!   DO itemp = 2, ntemp-1
!      WRITE(iu_therm, '(5e16.8)') temp(itemp),                  &
!                                  gamma_t(itemp), cv_t(itemp),  &
!                                  cp_t(itemp) - cv_t(itemp),    &
!                                  b0_s(itemp) - b0_t(itemp)
!   END DO
!   CLOSE(iu_therm)
!END IF
!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename=TRIM(flanhar)//'.celldm_ph'
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (ibrav==1 .OR. ibrav==2 .OR. ibrav==3 ) THEN
      WRITE(iu_therm,'("#   T (K)      celldm(1)      alpha_xx(x10^6)")' )
      DO itemp = 1, ntemp-1
         WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                               alphaf_anis_t(1,itemp)*1.D6
      END DO
   ELSEIF (ibrav==4 .OR. ibrav==6 .OR. ibrav==7 ) THEN
      WRITE(iu_therm,'("#   T (K)   celldm(1)   celldm(3)    alpha_xx(x10^6)   alpha_zz (x10^6")' )
      DO itemp = 1, ntemp-1
         WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                        celldmf_t(3,itemp), &
                                            alphaf_anis_t(1,itemp)*1.D6, &
                                            alphaf_anis_t(2,itemp)*1.D6
      END DO
   ELSEIF ( ibrav==5 ) THEN
      WRITE(iu_therm,'("#   T (K)   celldm(1)   celldm(4)    alpha_xx(x10^6)   alpha_zz (x10^6")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                        celldmf_t(4,itemp), &
                                             alphaf_anis_t(1,itemp)*1.D6, &
                                             alphaf_anis_t(2,itemp)*1.D6
      END DO
   ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)        celldm(2)        celldm(3)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,3e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                        celldmf_t(2,itemp), &
                                                        celldmf_t(3,itemp)
      END DO
   ELSEIF (ibrav==12 .OR. ibrav==13) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)         celldm(2)        celldm(3)        celldm(4)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,4e17.9)') temp(itemp), celldmf_t(1,itemp), &
                                                        celldmf_t(2,itemp), &
                                                        celldmf_t(3,itemp), &
                                                        celldmf_t(4,itemp)
      END DO
   ELSEIF (ibrav==-12 .OR. ibrav==-13) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)         celldm(2)        celldm(3)        celldm(5)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,4e17.9)') temp(itemp), celldmf_t(1,itemp), &
                                                        celldmf_t(2,itemp), &
                                                        celldmf_t(3,itemp), &
                                                        celldmf_t(5,itemp)
      END DO
   ELSEIF (ibrav==14) THEN
      WRITE(iu_therm,'("#   T (K)       celldm(1)         celldm(2)        &
                  &celldm(3)        celldm(4)        celldm(5)        celldm(6)")' )
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,6e15.7)') temp(itemp), celldmf_t(1,itemp), &
                                                        celldmf_t(2,itemp), &
                                                        celldmf_t(3,itemp), &
                                                        celldmf_t(4,itemp), &
                                                        celldmf_t(5,itemp), &
                                                        celldmf_t(6,itemp)
      END DO
   ELSE IF (ibrav==0) THEN
!
!  In this case we write nothing but do not stop
!
   ELSE
      CALL errore('write_anhar_anis','ibrav not programmed',1)
   END IF
   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_ph_freq_anhar_anis
