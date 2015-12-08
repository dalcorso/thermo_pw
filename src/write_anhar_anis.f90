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
USE anharmonic,     ONLY : alpha_anis_t, vmin_t, b0_t, celldm_t, beta_t
USE control_pwrun,  ONLY : ibrav_save
USE control_pressure, ONLY : pressure, pressure_kb
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm
REAL(DP) :: compute_omega_geo
CHARACTER(LEN=8) :: float_to_char

IF (my_image_id /= root_image) RETURN

CALL compute_alpha_anis(celldm_t, alpha_anis_t, temp, ntemp, ibrav_save)
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

IF (ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   filename=flanhar
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
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
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_alpha_anis(ibrav_save, celldm_t, alpha_anis_t, temp, ntemp, filename )

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
USE control_pressure, ONLY : pressure, pressure_kb
USE ph_freq_anharmonic, ONLY : alphaf_anis_t, vminf_t, b0f_t, celldmf_t, &
                               betaf_t
USE control_pwrun,  ONLY : ibrav_save
USE control_pressure, ONLY : pressure, pressure_kb
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm
REAL(DP) :: compute_omega_geo
CHARACTER(LEN=8) :: float_to_char

IF (my_image_id /= root_image) RETURN

CALL compute_alpha_anis(celldmf_t, alphaf_anis_t, temp, ntemp, ibrav_save)
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

IF (ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   iu_therm=2
   filename=TRIM(flanhar)//'_ph'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', &
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
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_alpha_anis(ibrav_save, celldmf_t, alphaf_anis_t, temp, ntemp, filename )

ENDIF

RETURN
END SUBROUTINE write_ph_freq_anhar_anis

SUBROUTINE write_grun_anharmonic_anis()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE ions_base,      ONLY : nat
USE cell_base,      ONLY : ibrav
USE thermo_mod,     ONLY : ngeo
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure, pressure_kb
USE ph_freq_thermodynamics, ONLY : ph_freq_save, phf_cv
USE ph_freq_anharmonic,     ONLY : celldmf_t, vminf_t, cvf_t, b0f_t, cpf_t, &
                                   b0f_s
USE grun_anharmonic, ONLY : alpha_an_g, grun_gamma_t, poly_grun, done_grun
USE ph_freq_module, ONLY : thermal_expansion_ph, ph_freq_type,  &
                           destroy_ph_freq, init_ph_freq
USE elastic_constants, ONLY : read_elastic, el_compliances
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic,      &
                               evaluate_fit_grad_quadratic
USE ifc,            ONLY : nq1_d, nq2_d, nq3_d
USE data_files,     ONLY : flanhar, fl_el_cons
USE io_global,      ONLY : ionode, stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char
INTEGER :: itemp, iu_therm, i, nq, imode, iq, degree, nvar, nwork
INTEGER :: itens, jtens
TYPE(ph_freq_type) :: ph_freq    ! the frequencies at the volumes at
                                 ! which the gruneisen parameters are 
                                 ! calculated
TYPE(ph_freq_type), ALLOCATABLE :: ph_grun(:)  ! the gruneisen parameters 
                                 ! recomputed at each temperature at the 
                                 ! geometry corresponding to that temperature
REAL(DP) :: cm(6), aux(6), alpha_aux(6), alpha(6), f, vm
REAL(DP), ALLOCATABLE :: grad(:), x(:)
INTEGER :: compute_nwork
LOGICAL :: exst

done_grun=.FALSE.
IF (my_image_id /= root_image) RETURN
!
!  Not implemented cases
!
IF ( ibrav<1 .OR. ibrav>11 ) THEN
   WRITE(stdout,'(5x,"Thermal expansions from Gruneisen parameters &
                     & not available")' )
   RETURN
END IF
!
!  First check if the elastic constants are on file
!
CALL read_elastic(fl_el_cons, exst)
!
!  If the elastic constants are not available, this calculation cannot be
!  don_ge
!
IF (.NOT. exst) THEN
   WRITE(stdout,'(5x,"The elastic constants are needed to compute &
                    &thermal expansions from Gruneisen parameters")' )
   RETURN
ENDIF
!
!  compute thermal expansion from gruneisen parameters. 
!  NB: alpha_an calculated by thermal_expansion_ph is multiplied by the 
!      elastic compliances
!
nq=ph_freq_save(1)%nq
CALL compute_degree(ibrav,degree,nvar)
nwork=compute_nwork()

CALL init_ph_freq(ph_freq, nat, nq1_d, nq2_d, nq3_d, nq, .FALSE.)
ALLOCATE(ph_grun(degree))
ALLOCATE(grad(degree))
ALLOCATE(x(degree))
DO i=1, degree
   CALL init_ph_freq(ph_grun(i), nat, nq1_d, nq2_d, nq3_d, nq, .FALSE.)
END DO

DO itemp = 1, ntemp
   IF (MOD(itemp,30)==0) WRITE(stdout,'(5x,"Computing temperature T=",f10.4,&
                                                   &" K")') temp(itemp)
   cm(:)=celldmf_t(:,itemp)
   vm = vminf_t(itemp)
   CALL compute_x(cm,x,degree,ibrav)
   ph_freq%nu= 0.0_DP
   DO i=1,degree
      ph_grun(i)%nu= 0.0_DP
   END DO
   ph_freq%wg=ph_freq_save(1)%wg
!
!  compute the frequencies once
!
   DO iq=1,nq
      DO imode=1,3*nat
         CALL evaluate_fit_quadratic(degree,nvar,x,f,poly_grun(1,imode,iq))
         CALL evaluate_fit_grad_quadratic(degree,nvar,x,grad,&
                                                     poly_grun(1,imode,iq))
         ph_freq%nu(imode,iq) = f 
         IF (f > 0.0_DP ) THEN
            DO i=1,degree
               ph_grun(i)%nu(imode,iq)=grad(i) / f
            END DO
         ELSE
            DO i=1,degree
               ph_grun(i)%nu(imode,iq)=0.0_DP
            END DO
         END IF
      END DO
   END DO
!
!  Loop over the number of independent crystal parameters for this
!  Bravais lattice
!
   DO i=1,degree
      CALL thermal_expansion_ph(ph_freq, ph_grun(i), temp(itemp), &
                                           alpha_aux(i))
   END DO
!
!  Here convert from derivatives with respect to crystal parameters to
!  derivative with respect to strain. 
!  NB: This is just the derivative of stress with respect to temperature.
!
   alpha(:)=0.0_DP
   SELECT CASE(ibrav)
       CASE(1,2,3)
!
!  cubic 
!
           alpha(1)=alpha_aux(1) * cm(1) / 3.0_DP              
           alpha(2)=alpha_aux(1) * cm(1) / 3.0_DP            
           alpha(3)=alpha_aux(1) * cm(1) / 3.0_DP            
       CASE(4,6,7)
!
!  hexagonal or tetragonal
!
           alpha(1)=(alpha_aux(1) * cm(1) - alpha_aux(2) * cm(3) ) / 2.0_DP   
           alpha(2)=alpha(2)             
           alpha(3)=alpha_aux(2) * cm(3)             
       CASE(5)
!
!  trigonal 
!
           alpha(3)=( 1.0_DP + 2.0_DP * cm(4) ) * ( alpha_aux(1) * cm(1) + &
                      2.0_DP * ( 1.0_DP - cm(4) ) * alpha_aux(2) ) / 3.0_DP

           alpha(1)=(alpha_aux(1) * cm(1) - alpha(3) ) / 2.0_DP   
           alpha(2)=alpha(2)             
       CASE(8,9,10,11)
!
!   orthorombic case
!
           alpha(1)=alpha_aux(1) * cm(1) - alpha_aux(2) * cm(2) &
                                         - alpha_aux(3) * cm(3)
           alpha(2)=alpha_aux(2) * cm(2)             
           alpha(3)=alpha_aux(3) * cm(3)             
   CASE DEFAULT
       CALL errore('write_grun_anharmonic_anis','ibrav not programmed',1)
   END SELECT
!
!  The thermal expansion needs to be multiplied by the elastic compliances
!
   aux=0.0_DP
   DO itens=1,6
      DO jtens=1,6
         aux(itens)=aux(itens) + el_compliances(itens,jtens) &
                                              *alpha(jtens)
      END DO
      alpha_an_g(:,itemp) = -aux(:) * ry_kbar / vm
   END DO
END DO

!CALL compute_cp(betab, vminf_t, b0f_t, phf_cv, cvf_t, cpf_t, b0f_s, &
!                                                             grun_gamma_t)
IF (ionode) THEN
!
!   here quantities calculated from the gruneisen parameters
!
   filename=TRIM(flanhar)//'.aux_grun'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
   CALL write_alpha_anis(ibrav, celldmf_t, alpha_an_g, temp, ntemp, filename )
END IF
done_grun=.TRUE.

CALL destroy_ph_freq(ph_freq)
DO i=1,degree
   CALL destroy_ph_freq(ph_grun(i))
END DO
DEALLOCATE(ph_grun)
DEALLOCATE(grad)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_grun_anharmonic_anis

SUBROUTINE write_alpha_anis(ibrav, celldmf_t, alpha_t, temp, ntemp, filename)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, ntemp
REAL(DP), INTENT(IN) :: celldmf_t(6,ntemp), alpha_t(6,ntemp), &
                        temp(ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER :: itemp, iu_therm

iu_therm=2

OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

IF (ibrav==1 .OR. ibrav==2 .OR. ibrav==3 ) THEN
   WRITE(iu_therm,'("#   T (K)      celldm(1)      alpha_xx(x10^6)")' )
   DO itemp = 1, ntemp-1
      WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                            alpha_t(1,itemp)*1.D6
   END DO
ELSEIF (ibrav==4 .OR. ibrav==6 .OR. ibrav==7 ) THEN
   WRITE(iu_therm,'("#   T (K)         celldm(1)          celldm(3)          alpha_xx(x10^6)         alpha_zz (x10^6)")' )
   DO itemp = 1, ntemp-1
      WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(3,itemp), &
                                                     alpha_t(1,itemp)*1.D6, &
                                                     alpha_t(3,itemp)*1.D6
   END DO
ELSEIF ( ibrav==5 ) THEN
   WRITE(iu_therm,'("#   T (K)   celldm(1)   celldm(4)    alpha_xx(x10^6)   alpha_zz (x10^6)")' )
   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(4,itemp), &
                                                     alpha_t(1,itemp)*1.D6, &
                                                     alpha_t(3,itemp)*1.D6
   END DO
ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
   WRITE(iu_therm,'("#   T (K)   celldm(1)     celldm(2)    celldm(3) &
             &alpha_xx(x10^6)  alpha_yy(x10^6)   alpha_zz(x10^6)")' )
   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e12.5,6e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(2,itemp), &
                                                     celldmf_t(3,itemp), &
                                                     alpha_t(1,itemp)*1.D6, &
                                                     alpha_t(2,itemp)*1.D6, &
                                                     alpha_t(3,itemp)*1.D6  
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
   CALL errore('write_alpha_anis','ibrav not programmed',1)
END IF

CLOSE(iu_therm)

RETURN
END SUBROUTINE write_alpha_anis

SUBROUTINE compute_alpha_anis(celldm_t, alpha_anis_t, temp, ntemp, ibrav)

USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, ntemp
REAL(DP), INTENT(IN) :: celldm_t(6,ntemp), temp(ntemp)
REAL(DP), INTENT(INOUT) :: alpha_anis_t(6,ntemp)

INTEGER :: itemp
REAL(DP) :: fact1, fact2, deriv1, deriv2

   alpha_anis_t=0.0_DP
   SELECT CASE (ibrav) 
      CASE(1,2,3) 
         DO itemp = 2, ntemp-1
            alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                            (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
            alpha_anis_t(2,itemp) = alpha_anis_t(1,itemp)
            alpha_anis_t(3,itemp) = alpha_anis_t(1,itemp)
         END DO
      CASE(4,6,7)
          DO itemp = 2, ntemp-1
             alpha_anis_t(1,itemp) = (celldm_t(1,itemp+1)-celldm_t(1,itemp-1)) / &
                             (temp(itemp+1)-temp(itemp-1)) / celldm_t(1,itemp)
             alpha_anis_t(2,itemp) = alpha_anis_t(1,itemp)
             alpha_anis_t(3,itemp) = ( celldm_t(3,itemp+1)*celldm_t(1,itemp+1)-   &
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
            alpha_anis_t(2,itemp) = alpha_anis_t(1,itemp)
            alpha_anis_t(3,itemp) = deriv1 / celldm_t(2, itemp) + deriv2 / fact2
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

RETURN
END SUBROUTINE compute_alpha_anis
