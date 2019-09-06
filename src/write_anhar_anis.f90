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
USE temperature,    ONLY : ntemp, temp
USE thermo_mod,     ONLY : ibrav_geo
USE thermo_sym,     ONLY : laue
USE thermodynamics, ONLY : ph_ce, ph_b_fact
USE anharmonic,     ONLY : alpha_anis_t, vmin_t, b0_t, celldm_t, beta_t, &
                           gamma_t, cv_t, ce_t, cp_t, b0_s, cpmce_anis,  &
                           el_cons_t, free_e_min_t, bths_t, ggamma_t,    &
                           bfact_t, lelastic, el_cons_s, el_comp_s
USE initial_conf,   ONLY : ibrav_save
USE control_thermo, ONLY : with_eigen
USE elastic_constants, ONLY : compute_elastic_compliances, &
                              write_el_cons_on_file
USE isoentropic,    ONLY : isostress_heat_capacity, thermal_stress,      &
                           gen_average_gruneisen, isoentropic_elastic_constants
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm, ibrav
INTEGER :: find_free_unit
REAL(DP) :: compute_omega_geo, csmct(6,6,ntemp)

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

CALL interpolate_cv(vmin_t, celldm_t, ph_ce, ce_t)
IF (lelastic) THEN
   CALL isostress_heat_capacity(vmin_t,el_cons_t,alpha_anis_t,temp, &
                                                         cpmce_anis,ntemp)
   cp_t=ce_t+cpmce_anis
   CALL compute_cv_bs_g(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t)

   CALL thermal_stress(el_cons_t,alpha_anis_t,bths_t,ntemp)
   CALL gen_average_gruneisen(vmin_t,bths_t,cv_t,temp,ggamma_t,ntemp)

   CALL isoentropic_elastic_constants(vmin_t,bths_t,cv_t,temp,csmct,ntemp)
   el_cons_s=el_cons_t + csmct
   DO itemp=2, ntemp-1
      CALL compute_elastic_compliances(el_cons_s(:,:,itemp), &
                                                   el_comp_s(:,:,itemp))
   ENDDO
ENDIF

ibrav=ibrav_geo(1)
IF (meta_ionode) THEN
!
!   here we write the anharmonic quantities calculated from the phonon dos
!   on disk
!
   filename='anhar_files/'//flanhar
   CALL add_pressure(filename)

   CALL write_ener_beta(temp, vmin_t, free_e_min_t, beta_t, ntemp, filename)
!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm'
   CALL add_pressure(filename)

   CALL write_alpha_anis(ibrav_save, celldm_t, alpha_anis_t, temp, ntemp, &
                                                                   filename )
!
!   here auxiliary quantities calculated from the phonon dos
!
   IF (lelastic) THEN
      !
      !   here the bulk modulus and the gruneisen parameter
      !
      filename="anhar_files/"//TRIM(flanhar)//'.bulk_mod'
      CALL add_pressure(filename)

      CALL write_bulk_anharm(temp, b0_t, b0_s, ntemp, filename)
!
!   the heat capacities
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat'
      CALL add_pressure(filename)

      CALL write_heat_anharm(temp, ce_t, cv_t, cp_t, ntemp, filename)

!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file heat
!
      filename='anhar_files/'//TRIM(flanhar)//'.heat_anis'
      CALL add_pressure(filename)
      CALL write_heat_anharm_anis(temp, ce_t, cv_t, cp_t, ntemp, filename)
!
!   Here the thermal stresses
!
      filename='anhar_files/'//TRIM(flanhar)//'.tstress'
      CALL add_pressure(filename)
      CALL write_thermal_stress(temp, bths_t, ntemp, filename)
!
!   Here the average Gruneisen paramater and the quantities that contribute
!   to it
!
      filename='anhar_files/'//TRIM(flanhar)//'.gamma'
      CALL add_pressure(filename)
      CALL write_gamma_anharm(temp, gamma_t, cv_t, beta_t, b0_t, ntemp, &
                                                                 filename)
!
!   Here the generalized Gruneisen parameters
!
      filename='anhar_files/'//TRIM(flanhar)//'.ggamma'
      CALL add_pressure(filename)
      CALL write_generalized_gamma(temp, ggamma_t, ntemp, filename)
   ELSE
!
!   only the interpolated heat capacity is available
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat'
      CALL add_pressure(filename)
      CALL write_heat_anharm_small(temp, ce_t, ntemp, filename)
   ENDIF
ENDIF
IF (lelastic) THEN
!
!   Here the elastic constants at constant entropy
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_cons_s'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_cons_s, &
                                                  b0_s, filename, 0)
!
!  and here the elastic compliances
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_comp_s'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_comp_s, &
                                                     b0_s, filename, 1)
ENDIF

IF (with_eigen) THEN
   CALL interpolate_b_fact_anis(celldm_t, ph_b_fact, bfact_t)
   filename="anhar_files/"//TRIM(flanhar)//'.anis'
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfact_t, ntemp, filename)
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
USE temperature,    ONLY : ntemp, temp
USE thermo_mod,     ONLY : ibrav_geo
USE thermo_sym,     ONLY : laue
USE ph_freq_thermodynamics, ONLY : phf_ce, phf_b_fact
USE ph_freq_anharmonic, ONLY : alphaf_anis_t, vminf_t, b0f_t, celldmf_t, &
                               betaf_t, gammaf_t, cvf_t, cef_t, cpf_t, b0f_s, &
                               cpmcef_anis, el_consf_t, lelasticf,       &
                               free_e_minf_t, bthsf_t, ggammaf_t, bfactf_t, &
                               el_consf_s, el_compf_s
USE elastic_constants, ONLY : compute_elastic_compliances, &
                              write_el_cons_on_file
USE initial_conf,   ONLY : ibrav_save
USE control_thermo, ONLY : with_eigen
USE isoentropic,    ONLY : isostress_heat_capacity, thermal_stress,      &
                           gen_average_gruneisen, isoentropic_elastic_constants
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit, ibrav
REAL(DP) :: compute_omega_geo, el_con_t(6,6,ntemp), csmct(6,6,ntemp)

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

CALL interpolate_cv(vminf_t, celldmf_t, phf_ce, cef_t)
IF (lelasticf) THEN
   CALL isostress_heat_capacity(vminf_t,el_consf_t,alphaf_anis_t,temp,&
                                                         cpmcef_anis,ntemp)
   cpf_t=cef_t+cpmcef_anis
   CALL compute_cv_bs_g(betaf_t, vminf_t, b0f_t, cvf_t, cpf_t, b0f_s, gammaf_t)
   CALL thermal_stress(el_consf_t,alphaf_anis_t,bthsf_t,ntemp)
   CALL gen_average_gruneisen(vminf_t,bthsf_t,cvf_t,temp,ggammaf_t,ntemp)

   CALL isoentropic_elastic_constants(vminf_t,bthsf_t,cvf_t,temp,csmct,ntemp)
   el_consf_s=el_consf_t + csmct
   DO itemp=2, ntemp-1
      CALL compute_elastic_compliances(el_consf_s(:,:,itemp), &
                                                        el_compf_s(:,:,itemp))
   ENDDO
ENDIF

ibrav=ibrav_geo(1)
IF (meta_ionode) THEN
!
!   here we plot the anharmonic quantities calculated from the phonon dos
!
   iu_therm=find_free_unit()
   filename='anhar_files/'//TRIM(flanhar)//'_ph'
   CALL add_pressure(filename)

   CALL write_ener_beta(temp, vminf_t, free_e_minf_t, betaf_t, ntemp, filename)

!
!  Here we write on output the celldm parameters and their derivative
!  with respect to temperature. 
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm_ph'
   CALL add_pressure(filename)

   CALL write_alpha_anis(ibrav_save, celldmf_t, alphaf_anis_t, temp, ntemp, &
                                                               filename )
!
!   here auxiliary quantities calculated from the phonon dos
!
   IF (lelasticf) THEN
      !
      !   here the bulk modulus and the gruneisen parameter
      !
      filename="anhar_files/"//TRIM(flanhar)//'.bulk_mod_ph'
      CALL add_pressure(filename)

      CALL write_bulk_anharm(temp, b0f_t, b0f_s, ntemp, filename)
!
!   the heat capacities
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat_ph'
      CALL add_pressure(filename)

      CALL write_heat_anharm(temp, cef_t, cvf_t, cpf_t, ntemp, filename)
!
!   Here the thermal stresses
!
      filename='anhar_files/'//TRIM(flanhar)//'.tstress_ph'
      CALL add_pressure(filename)
      CALL write_thermal_stress(temp, bthsf_t, ntemp, filename)
!
!   Here the average Gruneisen paramater and the quantities that contribute
!   to it
!
      filename='anhar_files/'//TRIM(flanhar)//'.gamma_ph'
      CALL add_pressure(filename)
      CALL write_gamma_anharm(temp, gammaf_t, cvf_t, betaf_t, b0f_t, ntemp, &
                                                                 filename)
!
!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file aux
!
      filename='anhar_files/'//TRIM(flanhar)//'.heat_anis_ph'
      CALL add_pressure(filename)
      CALL write_heat_anharm_anis(temp, cef_t, cvf_t, cpf_t, ntemp, filename)
!
!   Here the generalized Gruneisen parameters
!
      filename='anhar_files/'//TRIM(flanhar)//'.ggamma_ph'
      CALL add_pressure(filename)
      CALL write_generalized_gamma(temp, ggammaf_t, ntemp, filename)
   ELSE
!
!   only the interpolated heat capacity is available
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat_ph'
      CALL add_pressure(filename)
      CALL write_heat_anharm_small(temp, cef_t, ntemp, filename)
   ENDIF
ENDIF

IF (lelasticf) THEN
!
!   Here the elastic constants at constant entropy
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_cons_s_ph'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_consf_s, b0f_s, &
                                                              filename, 0)
!
!   and here the elastic compliances at constant entropy
!
   filename='anhar_files/'//TRIM(flanhar)//'.el_comp_s_ph'
   CALL add_pressure(filename)
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_compf_s, b0f_s, &
                                                              filename, 1)
ENDIF

IF (with_eigen) THEN
   CALL interpolate_b_fact_anis(celldmf_t, phf_b_fact, bfactf_t)
   filename="anhar_files/"//TRIM(flanhar)//'.anis_ph'
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfactf_t, ntemp, filename)
END IF

RETURN
END SUBROUTINE write_ph_freq_anhar_anis

SUBROUTINE write_grun_anhar_anis()
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE ions_base,      ONLY : nat
USE cell_base,      ONLY : ibrav
USE thermo_mod,     ONLY : ngeo, no_ph, reduced_grid
USE temperature,    ONLY : ntemp, temp
USE control_grun,   ONLY : vgrun_t, celldm_grun_t, b0_grun_t, lb0_t
USE ph_freq_thermodynamics, ONLY : ph_freq_save
USE grun_anharmonic, ONLY : alpha_an_g, grun_gamma_t, done_grun,          &
                           cp_grun_t, b0_grun_s, betab, grun_cpmce_anis,  &
                           cv_grun_t, ce_grun_t, lelastic_grun,           &
                           el_cons_grun_t, el_comp_grun_t, lelastic_grun, &
                           poly_order, p_grun_p2, poly_grun_red
USE polyfit_mod,    ONLY : compute_poly, compute_poly_deriv
USE ph_freq_module, ONLY : thermal_expansion_ph, ph_freq_type,  &
                           destroy_ph_freq, init_ph_freq
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic,      &
                               evaluate_quadratic_grad
USE isoentropic,    ONLY : isostress_heat_capacity
USE control_dosq,   ONLY : nq1_d, nq2_d, nq3_d
USE polynomial,     ONLY : clean_poly
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode, stdout
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm, i, nq, imode, iq, nvar, nwork
INTEGER :: itens, jtens, startq, lastq, iq_eff
TYPE(ph_freq_type) :: ph_freq    ! the frequencies at the volumes at
                                 ! which the gruneisen parameters are 
                                 ! calculated
TYPE(ph_freq_type), ALLOCATABLE :: ph_grun(:)  ! the gruneisen parameters 
                                 ! recomputed at each temperature at the 
                                 ! geometry corresponding to that temperature
REAL(DP) :: cm(6), aux(6), alpha_aux(6), alpha(6), f, g, vm
REAL(DP), ALLOCATABLE :: grad(:), x(:)
INTEGER :: compute_nwork, central_geo, find_free_unit

done_grun=.FALSE.
!
!  Not implemented cases
!
IF ( ibrav<1 .OR. ibrav>11 ) THEN
   WRITE(stdout,'(5x,"Thermal expansions from Gruneisen parameters &
                     & not available")' )
   RETURN
END IF
!
!  If the elastic constants are not available, this calculation cannot be done
!
IF (.NOT.lelastic_grun ) THEN
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
! divide the q points among the processors. Each processor has only a part
! of the q points and computes the contribution of these points to the
! anharmonic properties
!
CALL find_central_geo(ngeo, no_ph, central_geo)
nq=ph_freq_save(central_geo)%nq
startq=ph_freq_save(central_geo)%startq
lastq=ph_freq_save(central_geo)%lastq
CALL init_ph_freq(ph_freq, nat, nq1_d, nq2_d, nq3_d, startq, lastq,  &
                                                               nq, .FALSE.)
ph_freq%wg=ph_freq_save(central_geo)%wg
!
! now allocate space for each set of gruneisen parameters
!
nvar=crystal_parameters(ibrav)
nwork=compute_nwork()
ALLOCATE(ph_grun(nvar))
ALLOCATE(grad(nvar))
ALLOCATE(x(nvar))
DO i=1, nvar
   CALL init_ph_freq(ph_grun(i), nat, nq1_d, nq2_d, nq3_d, startq, &
                                                    lastq, nq, .FALSE.)
END DO
!
!  computes the anharmonic quantities at each temperature
!
alpha_an_g=0.0_DP
DO itemp = 1, ntemp
   IF (MOD(itemp,30)==0) &
             WRITE(6,'(5x,"Computing temperature ", i5 " / ",&
       & i5, 4x," T=",f12.2," K")') itemp, ntemp, temp(itemp)
!
!  Computes the volume at which the Gruneisen parameters and the frequencies
!  are interpolated
!
   cm(:)=celldm_grun_t(:,itemp)
   vm = vgrun_t(itemp)

   CALL compress_celldm(cm,x,nvar,ibrav)
!
!  compute the frequencies and the gruneisen parameters from the interpolating 
!  polynomial
!
   ph_freq%nu= 0.0_DP
   DO i=1,nvar
      ph_grun(i)%nu= 0.0_DP
   END DO
   iq_eff=0
   DO iq=startq, lastq
      iq_eff=iq_eff+1
      IF (reduced_grid) THEN
         DO i=1,nvar
            DO imode=1,3*nat
               CALL compute_poly(x(i), poly_order, &
                                            poly_grun_red(1,imode,i,iq),f)
!
!  this function gives the derivative with respect to x(i) multiplied by x(i)
!
               CALL compute_poly_deriv(x(i), poly_order, &
                                             poly_grun_red(1,imode,i,iq),g)

               ph_freq%nu(imode,iq_eff) = f
               IF (f > 0.0_DP ) THEN
                  ph_grun(i)%nu(imode,iq_eff)= - g / f 
               ELSE
                  ph_grun(i)%nu(imode,iq_eff)=0.0_DP
               END IF
            ENDDO
         ENDDO
      ELSE
         DO imode=1,3*nat
            CALL evaluate_fit_quadratic(nvar,x,f,p_grun_p2(imode,iq))
            CALL evaluate_quadratic_grad(nvar,x,grad,p_grun_p2(imode,iq))
            ph_freq%nu(imode,iq_eff) = f 
            IF (f > 0.0_DP ) THEN
               DO i=1,nvar
                  ph_grun(i)%nu(imode,iq_eff)=-grad(i) / f
               END DO
            ELSE
               DO i=1,nvar
                  ph_grun(i)%nu(imode,iq_eff)=0.0_DP
               END DO
            END IF
         END DO
      ENDIF
   END DO
!
!  Compute thermal expansion from Gruneisen parameters. Loop over the number 
!  of independent crystal parameters for this Bravais lattice
!  alpha calculated by thermal_expansion_ph is not multiplied by the elastic 
!  compliances
!
   DO i=1,nvar
      CALL thermal_expansion_ph(ph_freq, ph_grun(i), temp(itemp), &
                                alpha_aux(i), ce_grun_t(itemp))
   END DO
!
!  Here convert from derivatives with respect to crystal parameters to
!  derivatives with respect to strain. 
!
   CALL convert_ac_alpha(alpha_aux, alpha, cm, ibrav)
!
!  To get the thermal expansion we need to multiply by the elastic compliances
!
   aux=0.0_DP
   DO itens=1,6
      DO jtens=1,6
         aux(itens)=aux(itens) + el_comp_grun_t(itens,jtens,itemp)*alpha(jtens)
      END DO
   END DO
   alpha_an_g(:,itemp) = aux(:) * ry_kbar / vm
END DO
CALL mp_sum(alpha_an_g, world_comm)
CALL mp_sum(ce_grun_t, world_comm)
!
!  compute the volume thermal expansion as the trace of the thermal expansion 
!  tensor
!
betab(:)=alpha_an_g(1,:)+alpha_an_g(2,:)+alpha_an_g(3,:)
!
!  computes the other anharmonic quantities

CALL isostress_heat_capacity(vgrun_t,el_cons_grun_t,alpha_an_g,temp,&
                                                    grun_cpmce_anis,ntemp)
cp_grun_t = ce_grun_t + grun_cpmce_anis
CALL compute_cv_bs_g(betab, vgrun_t, b0_grun_t, cv_grun_t, &
                                     cp_grun_t, b0_grun_s, grun_gamma_t)

IF (meta_ionode) THEN
!
!   here quantities calculated from the gruneisen parameters
!
   filename='anhar_files/'//TRIM(flanhar)//'.celldm_grun'
   CALL add_pressure(filename)

   CALL write_alpha_anis(ibrav, celldm_grun_t, alpha_an_g, temp, ntemp, &
                                                                filename )
!
!  Here the average Gruneisen parameter and the quantities that form it
!
   filename="anhar_files/"//TRIM(flanhar)//'.aux_grun'
   CALL add_pressure(filename)
   CALL write_aux_grun(temp, betab, cp_grun_t, cv_grun_t, b0_grun_s, &
                                               b0_grun_t, ntemp, filename)
!
!  Here the average Gruneisen paramater and the quantities that contribute
!  to it
!
   filename='anhar_files/'//TRIM(flanhar)//'.gamma_grun'
   CALL add_pressure(filename)
   CALL write_gamma_anharm(temp, grun_gamma_t, cv_grun_t, betab, &
                                               b0_grun_t, ntemp, filename)
!
!  Here we write on output the anharmonic properties computed for
!  anisotropic solids, using the thermal expansion tensor, as opposed
!  to the volume thermal expansion used in the file aux
!
   filename='anhar_files/'//TRIM(flanhar)//'.heat_anis_grun'
   CALL add_pressure(filename)
   CALL write_heat_anharm_anis(temp, ce_grun_t, cv_grun_t, cp_grun_t, &
                                                           ntemp, filename)
ENDIF

done_grun=.TRUE.

CALL destroy_ph_freq(ph_freq)
DO i=1,nvar
   CALL destroy_ph_freq(ph_grun(i))
END DO
DEALLOCATE(ph_grun)

IF (.NOT.reduced_grid) THEN
   DO imode=1,3*nat
      DO iq=startq, lastq
         CALL clean_poly(p_grun_p2(imode,iq))
      ENDDO
   ENDDO
   DEALLOCATE(p_grun_p2)
ENDIF

DEALLOCATE(grad)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_grun_anhar_anis

SUBROUTINE write_alpha_anis(ibrav, celldmf_t, alpha_t, temp, ntemp, filename)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, ntemp
REAL(DP), INTENT(IN) :: celldmf_t(6,ntemp), alpha_t(6,ntemp), temp(ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

iu_therm=find_free_unit()

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
   DO itemp = 1, ntemp-1
      WRITE(iu_therm, '(e12.5,4e20.9)') temp(itemp), celldmf_t(1,itemp), &
                                                     celldmf_t(4,itemp), &
                                                     alpha_t(1,itemp)*1.D6, &
                                                     alpha_t(3,itemp)*1.D6
   END DO
ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
   WRITE(iu_therm,'("#   T (K)   celldm(1)     celldm(2)    celldm(3) &
             &alpha_xx(x10^6)  alpha_yy(x10^6)   alpha_zz(x10^6)")' )
   DO itemp = 1, ntemp-1
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
         alpha_anis_t(3,itemp) = deriv1 / celldm_t(1, itemp) + deriv2 / fact2
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

SUBROUTINE convert_ac_alpha(alpha_aux, alpha, cm, ibrav)
!
!  this subroutine receives the thermal expansion calculated with the
!  gruneisen parameters which are derivatives of the frequencies with 
!  respect to the crystal parameters and transforms it into a thermal
!  expansion tensor.
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER :: ibrav
REAL(DP) :: alpha(6), alpha_aux(6), cm(6)

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

      alpha(1)=(alpha_aux(1)*cm(1) - alpha(3) ) / 2.0_DP   
      alpha(2)=alpha(1)             
   CASE(8,9,10,11)
!
!   orthorhombic case
!
      alpha(1)=alpha_aux(1)*cm(1) - alpha_aux(2)*cm(2) - alpha_aux(3)*cm(3)
      alpha(2)=alpha_aux(2)*cm(2)             
      alpha(3)=alpha_aux(3)*cm(3)             
CASE DEFAULT
   CALL errore('convert_ac_alpha','ibrav not programmed',1)
END SELECT

RETURN
END SUBROUTINE convert_ac_alpha

SUBROUTINE write_heat_anharm_anis(temp, cet, cvt, cpt, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), cet(ntemp), cvt(ntemp), cpt(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# All heat capacities in (Ry/cell/K)")')
   WRITE(iu_therm,'("# T (K)      (C_s-C_e)(T)              C_P(T) &
                              &                C_V-C_e(T) ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e22.13)') temp(itemp), &
                 cpt(itemp)-cet(itemp), cpt(itemp), cvt(itemp)-cet(itemp)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_heat_anharm_anis

SUBROUTINE write_heat_anharm_small(temp, cet, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), cet(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# ")')
   WRITE(iu_therm,'("# T (K)        C_e(T) (Ry/cell/K) ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,e22.13)') temp(itemp), cet(itemp)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_heat_anharm_small

SUBROUTINE write_thermal_stress(temp, bths_t, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), bths_t(3,3,ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# Thermal stresses in (kbar)      ")')
   WRITE(iu_therm,'("# T (K)",5x," b_11", 9x," b_12", 9x," b_13", 9x,&
                       " b_22", 9x, " b_23", 9x, " b_33")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,6e14.6)') temp(itemp), bths_t(1,1,itemp),   &
                  bths_t(1,2,itemp), bths_t(1,3,itemp), bths_t(2,2,itemp), &
                  bths_t(2,3,itemp), bths_t(3,3,itemp) 
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_thermal_stress

SUBROUTINE write_generalized_gamma(temp, ggamma_t, ntemp, filename)
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode

IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), ggamma_t(3,3,ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# Generalized Gruneisen parameter      ")')
   WRITE(iu_therm,'("# T (K)",5x," g_11", 9x," g_12", 9x," g_13", 9x,&
                       " g_22", 9x, " g_23", 9x, " g_33")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,6e14.6)') temp(itemp),   ggamma_t(1,1,itemp), &
             ggamma_t(1,2,itemp), ggamma_t(1,3,itemp), ggamma_t(2,2,itemp), &
             ggamma_t(2,3,itemp), ggamma_t(3,3,itemp) 
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_generalized_gamma

SUBROUTINE set_elastic_grun()

USE anharmonic,         ONLY : lelastic
USE ph_freq_anharmonic, ONLY : lelasticf
USE grun_anharmonic,    ONLY : el_cons_grun_t, el_comp_grun_t, lelastic_grun
USE temperature,        ONLY : ntemp
USE anharmonic,         ONLY : el_cons_t, el_comp_t
USE ph_freq_anharmonic, ONLY : el_consf_t, el_compf_t

IMPLICIT NONE

INTEGER :: itemp

lelastic_grun=.FALSE.
IF (lelasticf) THEN
   DO itemp=1, ntemp
      el_cons_grun_t(:,:,itemp)=el_consf_t(:,:,itemp)
      el_comp_grun_t(:,:,itemp)=el_compf_t(:,:,itemp)
   ENDDO
   lelastic_grun=.TRUE.
ELSEIF(lelastic) THEN
   DO itemp=1, ntemp
      el_cons_grun_t(:,:,itemp)=el_cons_t(:,:,itemp)
      el_comp_grun_t(:,:,itemp)=el_comp_t(:,:,itemp)
   ENDDO
   lelastic_grun=.TRUE.
ENDIF

RETURN
END SUBROUTINE set_elastic_grun
