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
USE anharmonic,     ONLY : alpha_anis_t, vmin_t, b0_t, celldm_t, beta_t, &
                           gamma_t, cv_t, cp_t, b0_s, cpmcv_anis, el_cons_t, &
                           lelastic
USE control_grun,   ONLY : lb0_t
USE control_pwrun,  ONLY : ibrav_save
USE control_pressure, ONLY : pressure, pressure_kb
USE control_macro_elasticity, ONLY : macro_el
USE control_elastic_constants, ONLY : el_cons_available, el_cons_t_available
USE elastic_constants, ONLY : el_con
USE isoentropic,    ONLY : isostress_heat_capacity
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

CALL interpolate_cv(vmin_t, ph_cv, cv_t)
IF (lelastic) THEN
   CALL compute_cp(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)
   CALL isostress_heat_capacity(vmin_t,el_cons_t,alpha_anis_t,temp,&
                                                         cpmcv_anis,ntemp)
ENDIF

IF (ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   filename='anhar_files/'//flanhar
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# beta is the volume thermal expansion ")')
   IF (el_cons_t_available) THEN
      WRITE(iu_therm,'("#   T (K)        V(T) (a.u.)^3       B_0(T) (kbar)    &
                                                           &beta (x10^6)  ")' )

      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(e12.5,e20.13,2e15.6)') temp(itemp), vmin_t(itemp), &
                                           b0_t(itemp), beta_t(itemp)*1.D6
      END DO
   ELSE
      WRITE(iu_therm,'("#   T (K)        V(T) (a.u.)^3         &
                                                        &  beta (x10^6)  ")' )

      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(e12.5,e20.13,e15.6)') temp(itemp), vmin_t(itemp), &
                                                           beta_t(itemp)*1.D6
      END DO
   END IF

   CLOSE(iu_therm)
!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_alpha_anis(ibrav_save, celldm_t, alpha_anis_t, temp, ntemp, &
                                                                   filename )
!
!   here auxiliary quantities calculated from the phonon dos
!
   IF (lelastic) THEN
      filename='anhar_files/'//TRIM(flanhar)//'.aux'
      OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                                             FORM='FORMATTED')
      WRITE(iu_therm,'("# gamma is the average Gruneisen parameter ")')
      WRITE(iu_therm,'("#   T (K)       gamma(T)       C_v ( Ry / cell ) &
                    &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )

      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(5e16.8)') temp(itemp),                  &
                                     gamma_t(itemp), cv_t(itemp),  &
                                     cp_t(itemp) - cv_t(itemp),    &
                                     b0_s(itemp) - b0_t(itemp)
      END DO
      CLOSE(iu_therm)
!
!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file aux
!
      filename='anhar_files/'//TRIM(flanhar)//'.anis'
      IF (pressure /= 0.0_DP) &
         filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
      OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                                          FORM='FORMATTED')
      WRITE(iu_therm,'("#   T (K)       (C_p - C_v)(T)  " )' )

      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(2e16.8)') temp(itemp), cpmcv_anis(itemp)
      END DO
      CLOSE(iu_therm)
   END IF
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
USE ph_freq_thermodynamics, ONLY : phf_cv
USE ph_freq_anharmonic, ONLY : alphaf_anis_t, vminf_t, b0f_t, celldmf_t, &
                               betaf_t, gammaf_t, cvf_t, cpf_t, b0f_s, &
                               cpmcvf_anis, el_consf_t, lelasticf
USE elastic_constants, ONLY : el_con
USE control_grun,   ONLY : lb0_t
USE control_pwrun,  ONLY : ibrav_save
USE control_pressure, ONLY : pressure, pressure_kb
USE control_macro_elasticity, ONLY : macro_el
USE control_elastic_constants, ONLY : el_cons_available
USE isoentropic,    ONLY : isostress_heat_capacity
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm
REAL(DP) :: compute_omega_geo, el_con_t(6,6,ntemp)
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

IF (lelasticf) THEN
   CALL interpolate_cv(vminf_t, phf_cv, cvf_t)
   CALL compute_cp(betaf_t, vminf_t, b0f_t, cvf_t, cpf_t, b0f_s, gammaf_t)
   CALL isostress_heat_capacity(vminf_t,el_consf_t,alphaf_anis_t,temp,&
                                                         cpmcvf_anis,ntemp)
ENDIF

IF (ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   iu_therm=2
   filename='anhar_files/'//TRIM(flanhar)//'_ph'
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
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm_ph'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))

   CALL write_alpha_anis(ibrav_save, celldmf_t, alphaf_anis_t, temp, ntemp, &
                                                               filename )
!
!   here auxiliary quantities calculated from the phonon dos
!
   IF (lelasticf) THEN
      filename='anhar_files/'//TRIM(flanhar)//'.aux_ph'
      OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
      WRITE(iu_therm,'("# gamma is the average Gruneisen parameter ")')
      WRITE(iu_therm,'("#   T (K)       gamma(T)       C_v ( Ry / cell ) &
                    &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )

      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(5e16.8)') temp(itemp),                  &
                                     gammaf_t(itemp), cvf_t(itemp),  &
                                     cpf_t(itemp) - cvf_t(itemp),    &
                                     b0f_s(itemp) - b0f_t(itemp)
      END DO
      CLOSE(iu_therm)
!
!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file aux
!
      filename='anhar_files/'//TRIM(flanhar)//'.anis_ph'
      IF (pressure /= 0.0_DP) &
         filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
      OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                                          FORM='FORMATTED')
      WRITE(iu_therm,'("#   T (K)       (C_p - C_v)(T)  " )' )

      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(2e16.8)') temp(itemp), cpmcvf_anis(itemp)
      END DO
      CLOSE(iu_therm)
   END IF

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
USE control_grun,   ONLY : lv0_t
USE control_mur,    ONLY : celldm0, vmin
USE thermodynamics,         ONLY : ph_cv
USE ph_freq_thermodynamics, ONLY : ph_freq_save, phf_cv
USE anharmonic,             ONLY : celldm_t, vmin_t, b0_t, cv_t, lelastic, &
                                   el_comp_t
USE ph_freq_anharmonic,     ONLY : celldmf_t, vminf_t, b0f_t, cvf_t, lelasticf
USE grun_anharmonic, ONLY : alpha_an_g, grun_gamma_t, poly_grun, done_grun, &
                            cp_grun_t, b0_grun_s, betab
USE ph_freq_module, ONLY : thermal_expansion_ph, ph_freq_type,  &
                           destroy_ph_freq, init_ph_freq
USE control_macro_elasticity, ONLY : macro_el
USE control_grun,     ONLY : lb0_t
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE elastic_constants, ONLY :  el_compliances
USE control_elastic_constants, ONLY : el_cons_available, el_cons_t_available, &
                                      omega0
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic,      &
                               evaluate_fit_grad_quadratic
USE control_dosq,   ONLY : nq1_d, nq2_d, nq3_d
USE data_files,     ONLY : flanhar, fl_el_cons
USE io_global,      ONLY : ionode, stdout
USE mp_images,      ONLY : my_image_id, root_image, intra_image_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char
INTEGER :: itemp, iu_therm, i, nq, imode, iq, degree, nvar, nwork
INTEGER :: itens, jtens, nstart, nlast
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
!  If the elastic constants are not available, this calculation cannot be
!  don_ge
!
IF ( .NOT.(lelastic.OR.lelasticf) ) THEN
   WRITE(stdout,'(5x,"The elastic constants are needed to compute ")')
   WRITE(stdout,'(5x,"thermal expansions from Gruneisen parameters")')
   RETURN
ENDIF

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the anharmonic properties from &
                                                   &Gruneisen parameters")')
WRITE(stdout,'(5x,"Writing on file anhar_files/",a)') TRIM(flanhar)// &
                                                                  '.aux_grun'
WRITE(stdout,'(2x,76("+"),/)')
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

CALL divide(intra_image_comm, ntemp, nstart, nlast)
alpha_an_g=0.0_DP
DO itemp = nstart, nlast
   IF (MOD(itemp-nstart+1,30)==0) &
             WRITE(6,'(5x,"Computing temperature ", i5 " / ",&
       & i5, 4x," T=",f12.2," K")') itemp-nstart+1, nlast-nstart+1, temp(itemp)
   IF (lv0_t.AND.ltherm_freq) THEN
      cm(:)=celldmf_t(:,itemp)
      vm = vminf_t(itemp)
   ELSEIF (lv0_t.AND.ltherm_dos) THEN
      cm(:)=celldm_t(:,itemp)
      vm = vmin_t(itemp)
   ELSE
      cm(:)=celldm0(:)
      vm = omega0
   END IF
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
           alpha(2)=alpha(1)             
           alpha(3)=alpha_aux(2) * cm(3)             
       CASE(5)
!
!  trigonal 
!
           alpha(3)=( 1.0_DP + 2.0_DP * cm(4) ) * ( alpha_aux(1) * cm(1) + &
                      2.0_DP * ( 1.0_DP - cm(4) ) * alpha_aux(2) ) / 3.0_DP

           alpha(1)=(alpha_aux(1) * cm(1) - alpha(3) ) / 2.0_DP   
           alpha(2)=alpha(1)             
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
   IF (el_cons_t_available) THEN
      DO itens=1,6
         DO jtens=1,6
            aux(itens)=aux(itens) + el_comp_t(itens,jtens,itemp) &
                                                 *alpha(jtens)
         END DO
      END DO
   ELSEIF (el_cons_available) THEN
      DO itens=1,6
         DO jtens=1,6
            aux(itens)=aux(itens) + el_compliances(itens,jtens)*alpha(jtens)
         END DO
      END DO
   END IF
   alpha_an_g(:,itemp) = -aux(:) * ry_kbar / vm
END DO
CALL mp_sum(alpha_an_g, intra_image_comm)

betab(:)=alpha_an_g(1,:)+alpha_an_g(2,:)+alpha_an_g(3,:)

IF (ltherm_freq) THEN
   CALL compute_cp(betab, vminf_t, b0f_t, cvf_t, cp_grun_t, b0_grun_s, &
                                                               grun_gamma_t)
ELSEIF (ltherm_dos) THEN
   CALL compute_cp(betab, vmin_t, b0_t, cv_t, cp_grun_t, b0_grun_s, &
                                                               grun_gamma_t)
ELSE
   CALL interpolate_cv(vmin_t, phf_cv, cv_t)
   CALL compute_cp(betab, vmin_t, b0_t, cv_t, cp_grun_t, &
                                                 b0_grun_s, grun_gamma_t)
END IF

IF (ionode) THEN
!
!   here quantities calculated from the gruneisen parameters
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm_grun'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
   IF (ltherm_freq) THEN
      CALL write_alpha_anis(ibrav, celldmf_t, alpha_an_g, temp, ntemp, &
                                                                filename )
   ELSEIF (ltherm_dos) THEN
      CALL write_alpha_anis(ibrav, celldm_t, alpha_an_g, temp, ntemp, &
                                                                filename )
   ELSE
!
!     here we plot the celldm that has been used to compute the thermal
!     expansion, not the celldm that results from the thermal expansion.
!
      DO itemp=1,ntemp
         celldm_t(:,itemp)=celldm0(:)
      ENDDO
      CALL write_alpha_anis(ibrav, celldm_t, alpha_an_g, temp, ntemp, &
                                                                filename )
      celldm_t=0.0_DP
   ENDIF


   filename="anhar_files/"//TRIM(flanhar)//'.aux_grun'
   IF (pressure /= 0.0_DP) &
      filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                                         FORM='FORMATTED')
   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
   WRITE(iu_therm,'("#   T (K)     beta(T)    gamma(T)      &
              &   (C_p - C_v)(T)      (B_S - B_T) (T) (kbar) " )' )

   IF (ltherm_freq) THEN
      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(5e16.8)') temp(itemp), betab(itemp)*1.D6,   &
                     grun_gamma_t(itemp), cp_grun_t(itemp) - cvf_t(itemp), &
                     b0_grun_s(itemp) - b0f_t(itemp)
      END DO
   ELSE
      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(5e16.8)') temp(itemp), betab(itemp)*1.D6,    &
                     grun_gamma_t(itemp), cp_grun_t(itemp) - cv_t(itemp), &
                     b0_grun_s(itemp) - b0_t(itemp)
      END DO
   END IF
   CLOSE(iu_therm)
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
   WRITE(iu_therm,'("#   T (K)        celldm(1)         alpha_xx(x10^6)")' )
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
