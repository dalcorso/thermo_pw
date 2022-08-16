!
! Copyright (C) 2013-2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_anhar()
!-----------------------------------------------------------------------
!
!   This routine writes on output the anharmonic quantities calculated
!   at the volume that minimizes the free energy (computed from phonon dos)
!   It writes beta, cp, cp-cv, bs, bs-bt, gamma
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp
USE anharmonic,     ONLY : alpha_t, beta_t, gamma_t, cp_t, cv_t,         &
                           ener_t, free_ener_t, entropy_t, b0_s, vmin_t, &
                           free_e_min_t, b0_t, b01_t, b02_t, bfact_t,    &
                           e0_t, beta_noe_t, vmin_noe_t, b0_noe_s,       &
                           b0_noe_t, cp_noe_t, cv_noe_t, gamma_noe_t,    &
                           noelcvg
USE el_anharmonic,  ONLY : el_b0_t, el_beta_t, el_gamma_t, el_cp_t,      &
                           el_b0_t
USE control_eldos,  ONLY : lel_free_energy
USE control_thermo, ONLY : with_eigen
USE data_files,     ONLY : flanhar

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
LOGICAL :: subtract_el
!
!  start by computing the thermal exspansion
!
CALL compute_beta(vmin_t, beta_t, temp, ntemp)
alpha_t = beta_t / 3.0_DP
!
!  here compute the isobaric heat capacity, the isoentropic bulk modulus
!  and the Gruneisen parameter. Remove the electronic contribution
!  from cv_t before computing the gruneisen parameter if required on input.
!
subtract_el=(lel_free_energy.AND.noelcvg)
CALL compute_cp_bs_g(beta_t, vmin_t, b0_t, cv_t, cp_t, b0_s, gamma_t, &
                      subtract_el)
!
!  Determine the electronic contribution to beta, cp, and gamma, computing
!  the same quantities without the electronic contribution
!
IF (lel_free_energy) THEN
   CALL compute_beta(vmin_noe_t, beta_noe_t, temp, ntemp)
   CALL compute_cp_bs_g(beta_noe_t, vmin_noe_t, b0_noe_t, cv_noe_t, cp_noe_t, &
   b0_noe_s, gamma_noe_t, .FALSE.)
   el_beta_t=beta_t - beta_noe_t
   el_b0_t = b0_t - b0_noe_t
   el_cp_t = cp_t - cp_noe_t 
   el_gamma_t=gamma_t-gamma_noe_t
ENDIF
!
!   here we plot the quantities calculated from the phonon dos
!
filename="anhar_files/"//TRIM(flanhar)
CALL add_pressure(filename)

CALL write_ener_beta(temp, vmin_t, free_e_min_t, beta_t, ntemp, filename)
!
!   here the bulk modulus 
!
filename="anhar_files/"//TRIM(flanhar)//'.bulk'
CALL add_pressure(filename)

CALL write_bulk_anharm(temp, b0_t, b0_s, ntemp, filename)
!
!   here the derivative of the bulk modulus with respect to pressure
!
filename="anhar_files/"//TRIM(flanhar)//'.dbulk'
CALL add_pressure(filename)

CALL write_dbulk_anharm(temp, b01_t, b02_t, ntemp, filename)
!
!   here the heat capacities (constant strain and constant stress)
!
filename="anhar_files/"//TRIM(flanhar)//'.heat'
CALL add_pressure(filename)

CALL write_heat_anhar(temp, cv_t, cv_t, cp_t, ntemp, filename)
!
!  here the average Gruneisen parameter and the quantities that form it
!
filename="anhar_files/"//TRIM(flanhar)//'.gamma'
CALL add_pressure(filename)

CALL write_gamma_anharm(temp, gamma_t, cv_t, beta_t, b0_t, ntemp, filename)
!
!   here all the interpolated harmonic quantities
!
filename="anhar_files/"//TRIM(flanhar)//'.therm'
CALL add_pressure(filename)
CALL write_thermo_anharm(temp, ntemp, e0_t, ener_t, free_ener_t, &
                                            entropy_t, cv_t, filename)
!
!   here the b factors
!
IF (with_eigen) THEN
   filename="anhar_files/"//TRIM(flanhar)
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfact_t, ntemp, filename)
END IF
!
!   Write also the electronic contribution to the anharmonic quantities
!
IF (lel_free_energy) THEN
   filename="anhar_files/"//TRIM(flanhar)//'.el_anhar'
   CALL add_pressure(filename)
   CALL write_el_anhar(temp, el_gamma_t, el_cp_t, el_beta_t, &
                                        el_b0_t, ntemp, filename)
ENDIF

RETURN
END SUBROUTINE write_anhar
!
!-------------------------------------------------------------------------
SUBROUTINE write_ph_freq_anhar()
!-------------------------------------------------------------------------
!
!   This routine writes on output the anharmonic quantities calculated
!   at the volume that minimizes the free energy (computed from Brillouin
!   zone integration). It writes beta, cp, cp-cv, bs, bs-bt, gamma.
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp
USE ph_freq_anharmonic, ONLY : alphaf_t, betaf_t, gammaf_t, cpf_t, cvf_t,   &
                           enerf_t, free_enerf_t, entropyf_t, b0f_s,        &
                           free_e_minf_t, vminf_t, b0f_t, b01f_t, b02f_t,   &
                           bfactf_t, e0f_t
USE anharmonic,     ONLY : vmin_noe_t, beta_noe_t, b0_noe_t, cv_noe_t,      &
                           cp_noe_t, b0_noe_s, gamma_noe_t, noelcvg
USE el_anharmonic,  ONLY : el_b0f_t, el_betaf_t, el_cpf_t, el_gammaf_t
USE control_eldos,  ONLY : lel_free_energy
USE control_thermo, ONLY : with_eigen
USE data_files,     ONLY : flanhar

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
LOGICAL :: subtract_el
!
!  start by computing the thermal exspansion
!
CALL compute_beta(vminf_t, betaf_t, temp, ntemp)

alphaf_t = betaf_t / 3.0_DP

subtract_el=(lel_free_energy.AND.noelcvg)
CALL compute_cp_bs_g(betaf_t, vminf_t, b0f_t, cvf_t, cpf_t, b0f_s, gammaf_t,&
                                             subtract_el)
!
IF (lel_free_energy) THEN
   CALL compute_beta(vmin_noe_t, beta_noe_t, temp, ntemp)
   CALL compute_cp_bs_g(beta_noe_t, vmin_noe_t, b0_noe_t, cv_noe_t, cp_noe_t, &
   b0_noe_s, gamma_noe_t,.FALSE.)

   el_betaf_t=betaf_t - beta_noe_t
   el_b0f_t = b0f_t - b0_noe_t
   el_cpf_t = cpf_t - cp_noe_t 
   el_gammaf_t=gammaf_t-gamma_noe_t
ENDIF
!
!   here we plot the quantities calculated from the explicit sum over 
!   frequencies
!
filename="anhar_files/"//TRIM(flanhar)//'_ph'
CALL add_pressure(filename)

CALL write_ener_beta(temp, vminf_t, free_e_minf_t, betaf_t, ntemp, filename)
!
!   here the bulk modulus 
!
filename="anhar_files/"//TRIM(flanhar)//'.bulk_ph'
CALL add_pressure(filename)

CALL write_bulk_anharm(temp, b0f_t, b0f_s, ntemp, filename)
!
!   here the derivative of the bulk modulus with respect to pressure
!
filename="anhar_files/"//TRIM(flanhar)//'.dbulk_ph'
CALL add_pressure(filename)

CALL write_dbulk_anharm(temp, b01f_t, b02f_t, ntemp, filename)
!
!   here the heat capacities (constant strain and constant stress)
!
filename="anhar_files/"//TRIM(flanhar)//'.heat_ph'
CALL add_pressure(filename)

CALL write_heat_anhar(temp, cvf_t, cvf_t, cpf_t, ntemp, filename)
!
!  here the average Gruneisen parameter and the quantities that form it
!
filename="anhar_files/"//TRIM(flanhar)//'.gamma_ph'
CALL add_pressure(filename)

CALL write_gamma_anharm(temp, gammaf_t, cvf_t, betaf_t, b0f_t, ntemp, filename)
!
!  here the interpolated harmonic quantities
!
filename="anhar_files/"//TRIM(flanhar)//'.therm_ph'
CALL add_pressure(filename)
CALL write_thermo_anharm(temp, ntemp, e0f_t, enerf_t, free_enerf_t, & 
                                            entropyf_t, cvf_t, filename)
!
!   here the b factors
!
IF (with_eigen) THEN
   filename="anhar_files/"//TRIM(flanhar)//'_ph'
   CALL add_pressure(filename)
   CALL write_anharm_bfact(temp, bfactf_t, ntemp, filename)
END IF
!
!   Write also the electronic contribution to the anharmonic quantities
!
IF (lel_free_energy) THEN
   filename="anhar_files/"//TRIM(flanhar)//'.elf_anhar'
   CALL add_pressure(filename)
   CALL write_el_anhar(temp, el_gammaf_t, el_cpf_t, el_betaf_t, &
                                        el_b0f_t, ntemp, filename)
ENDIF

RETURN
END SUBROUTINE write_ph_freq_anhar
!
!-------------------------------------------------------------------------
SUBROUTINE write_grun_anharmonic()
!-------------------------------------------------------------------------
!
!  This routine computes the anharmonic properties from the Gruneisen
!  parameters. Used when lmurn=.TRUE.. (It sssumes to have the volume and
!  the bulk modulus as a function of temperature).
!  
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE ions_base,      ONLY : nat
USE temperature,    ONLY : ntemp, temp
USE thermo_mod,     ONLY : ngeo, no_ph
USE ph_freq_thermodynamics, ONLY : ph_freq_save
USE grun_anharmonic, ONLY : betab, cp_grun_t, b0_grun_s, &
                           grun_gamma_t, poly_grun, poly_degree_grun, &
                           ce_grun_t, cv_grun_t
USE ph_freq_module, ONLY : thermal_expansion_ph, ph_freq_type,  &
                           destroy_ph_freq, init_ph_freq
USE el_anharmonic,  ONLY : el_beta_t, el_cef_t
USE anharmonic,     ONLY : noelcvg
USE polyfit_mod,    ONLY : compute_poly, compute_poly_deriv
USE control_grun,   ONLY : vgrun_t, b0_grun_t
USE control_dosq,   ONLY : nq1_d, nq2_d, nq3_d
USE control_eldos,  ONLY : lel_free_energy
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
LOGICAL :: subtract_el
INTEGER :: itemp, nq, imode, iq, startq, lastq, iq_eff
TYPE(ph_freq_type) :: ph_freq    ! the frequencies at the volumes at
                                 ! which the gruneisen parameters are 
                                 ! calculated
TYPE(ph_freq_type) :: ph_grun    ! the gruneisen parameters recomputed
                                 ! at each temperature at the volume
                                 ! corresponding to that temperature
REAL(DP) :: vm, f, g
INTEGER :: central_geo
!
!  divide the q vectors among the processors. Allocate space to save the
!  frequencies and the gruneisen parameter on a part of the mesh of q points.
!
CALL find_central_geo(ngeo,no_ph,central_geo)
nq=ph_freq_save(central_geo)%nq
startq=ph_freq_save(central_geo)%startq
lastq=ph_freq_save(central_geo)%lastq

CALL init_ph_freq(ph_grun, nat, nq1_d, nq2_d, nq3_d, startq, &
                                                         lastq, nq, .FALSE.)
CALL init_ph_freq(ph_freq, nat, nq1_d, nq2_d, nq3_d, startq, &
                                                         lastq, nq, .FALSE.)
ph_freq%wg(:)=ph_freq_save(central_geo)%wg(:)

betab=0.0_DP
DO itemp = 1, ntemp
!
!  Set the volume at this temperature.
!
   vm=vgrun_t(itemp)
!
!  Use the fitting polynomials to find the frequencies and the gruneisen
!  parameters at this volume.
!
   ph_freq%nu= 0.0_DP
   ph_grun%nu= 0.0_DP
   iq_eff=0
   DO iq=startq, lastq
      iq_eff=iq_eff+1
      DO imode=1,3*nat
         CALL compute_poly(vm, poly_degree_grun, poly_grun(:,imode,iq),f)
         CALL compute_poly_deriv(vm, poly_degree_grun, poly_grun(:,imode,iq),g)
!
!     g here is d w / d V 
!     ph_grun%nu will contain the gruneisen parameter divided by the volume
!     as requested by the thermal_expansion_ph routine
!
         ph_freq%nu(imode,iq_eff)=f
         IF ( f > 0.0_DP) THEN 
            ph_grun%nu(imode,iq_eff)=-g/f
         ELSE
            ph_grun%nu(imode,iq_eff) = 0.0_DP
         ENDIF
      ENDDO
   ENDDO
!
!  computes thermal expansion from gruneisen parameters. 
!  this routine gives betab multiplied by the bulk modulus
!
   CALL thermal_expansion_ph(ph_freq, ph_grun, temp(itemp), betab(itemp), &
                                                            ce_grun_t(itemp))
!
!  divide by the bulk modulus
!
   betab(itemp)=betab(itemp) * ry_kbar / b0_grun_t(itemp)

END DO
CALL mp_sum(betab, world_comm)
CALL mp_sum(ce_grun_t, world_comm)
cv_grun_t=ce_grun_t
!
!  Add the electronic quantities 
!
IF (lel_free_energy) THEN
   betab=betab+el_beta_t
   ce_grun_t=ce_grun_t+el_cef_t
   cv_grun_t=ce_grun_t
ENDIF
!
!  computes the other anharmonic quantities
!
subtract_el=(lel_free_energy.AND.noelcvg)
CALL compute_cp_bs_g(betab, vgrun_t, b0_grun_t, cv_grun_t, cp_grun_t, &
                                  b0_grun_s, grun_gamma_t, subtract_el)
IF (meta_ionode) THEN
!
!  Here the average Gruneisen parameter and the quantities that form it
!
   filename="anhar_files/"//TRIM(flanhar)//'.gamma_grun'
   CALL add_pressure(filename)
   CALL write_gamma_anharm(temp, grun_gamma_t, ce_grun_t, betab, &
                                               b0_grun_t, ntemp, filename)
!
!   here quantities calculated from the mode Gruneisen parameters
!   We add the electronic contribution before writing if needed
!
   filename="anhar_files/"//TRIM(flanhar)//'.aux_grun'
   CALL add_pressure(filename)
   CALL write_aux_grun(temp, betab, cp_grun_t, ce_grun_t, b0_grun_s, &
                                               b0_grun_t, ntemp, filename) 
ENDIF

CALL destroy_ph_freq(ph_freq)
CALL destroy_ph_freq(ph_grun)

RETURN
END SUBROUTINE write_grun_anharmonic
!
!-------------------------------------------------------------------------
SUBROUTINE write_el_fit_harmonic()
!-------------------------------------------------------------------------
!
USE control_ev,     ONLY : ieos
USE data_files,     ONLY : flelanhar
USE el_anharmonic,  ONLY : vmine_t, b0e_t, b01e_t, b02e_t, free_e_mine_t
USE io_global,      ONLY : meta_ionode
USE temperature,    ONLY : temp, ntemp

IMPLICIT NONE
INTEGER :: iu_therm
INTEGER :: itemp
INTEGER :: find_free_unit
CHARACTER(LEN=256) :: filename

filename='anhar_files/'//TRIM(flelanhar)
IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("#  ")')
   IF (ieos==2) THEN
      WRITE(iu_therm,'("#   T (K)        V_T(T) (a.u.)^3        B_0(T) (kbar) &
                &    dB/dp      d^2B/dp^2 (1/kbar)    free energy (Ry) ")')
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,5e22.13)') temp(itemp), &
                 vmine_t(itemp), b0e_t(itemp), b01e_t(itemp), &
                 b02e_t(itemp), free_e_mine_t(itemp)
      ENDDO
   ELSE
      WRITE(iu_therm,'("#   T (K)        V_T(T) (a.u.)^3        B_0(T) (kbar) &
                               &    dB0/dp       free energy (Ry) ")')
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e12.5,4e22.13)') temp(itemp), &
                 vmine_t(itemp), b0e_t(itemp), b01e_t(itemp), &
                 free_e_mine_t(itemp)
      ENDDO
   ENDIF

   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_el_fit_harmonic
!
!-------------------------------------------------------------------------
SUBROUTINE compute_beta(vmin_t, beta_t, temp, ntemp)
!-------------------------------------------------------------------------
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
!
!------------------------------------------------------------------------
SUBROUTINE write_ener_beta(temp, vmin, emin, beta, ntemp, filename)
!------------------------------------------------------------------------
!
!  This routine writes on file the minimum volume, the minimum energy,
!  and the thermal expansion as a function of temperature. 
!
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
                   &      beta x 10^6 (1/K)")' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,2e23.13,e18.8)') temp(itemp), &
                vmin(itemp), emin(itemp), beta(itemp)*1.D6
   END DO
  
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_ener_beta
!
!------------------------------------------------------------------------
SUBROUTINE write_ener_beta_t(press, vmin, vminv0, emin, beta, npress, &
                             itemp, filename)
!------------------------------------------------------------------------
!
!  This routine writes on output the minimum volume, the minimum energy,
!  and the thermal expansion as a function of pressure at several 
!  temperatures
!
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
USE temperature, ONLY : temp

IMPLICIT NONE
INTEGER, INTENT(IN) :: npress, itemp
REAL(DP), INTENT(IN) :: press(npress), emin(npress), vmin(npress), &
                        vminv0(npress), beta(npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_press
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   OPEN(UNIT=iu_press, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                                         FORM='FORMATTED')
   WRITE(iu_press,'("#   p (kbar)      Volume(p) ((a.u.)^3)      V/V(300K) &
        &    E_min (p) (Ry)         beta(p) x 10^6 (1/K)  T=",f16.4," K")' ) &
                                                               temp(itemp)
   DO ipress=1,npress
      WRITE(iu_press,'(5e20.10)') press(ipress), vmin(ipress), vminv0(ipress),&
                                  emin(ipress), beta(ipress)*1.D6
   ENDDO
   CLOSE(UNIT=iu_press, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_ener_beta_t
!
!------------------------------------------------------------------------
SUBROUTINE write_ener_beta_vol(temp, vmin, vminv0, emin, beta, ntemp, &
                                                               filename)
!------------------------------------------------------------------------
!
!  This routine writes on output the minimum volume, the minimum volume
!  divided by the volume at 300 K, the minimum energy,
!  and the thermal expansion as a function of temperature.
!
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), emin(ntemp), vmin(ntemp), &
                        vminv0(ntemp), beta(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# beta is the volume thermal expansion ")')
   WRITE(iu_therm,'("#   T (K)         V(T) (a.u.)^3   V(T)/V0   F (T) (Ry) &
                   &      beta x 10^6 (1/K)")' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e23.13,e18.8)') temp(itemp), &
            vmin(itemp), vminv0(itemp), emin(itemp), beta(itemp)*1.D6
   END DO
  
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_ener_beta_vol
!
!------------------------------------------------------------------------
SUBROUTINE write_bulk_anharm(temp, b0t, b0s, ntemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the isothermal bulk modulus, the
!   isoentropic bulk modulus and their difference as a function of 
!   temperature.
!
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), b0t(ntemp), b0s(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("#  ")')
   WRITE(iu_therm,'("#   T (K)        B_T(T) (kbar)        B_S(T) (kbar) &
                                 &       B_S(T)-B_T(T) (kbar) ")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e22.13)') temp(itemp), &
                 b0t(itemp), b0s(itemp), b0s(itemp)-b0t(itemp)
   END DO
  
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_bulk_anharm
!
!------------------------------------------------------------------------
SUBROUTINE write_bulk_anharm_t(press, b0t, b0s, npress, itemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the isothermal bulk modulus, the
!   isoentropic bulk modulus, and their difference as a function of 
!   pressure for several temperatures.
!
USE kinds,       ONLY : DP
USE io_global,   ONLY : meta_ionode
USE temperature, ONLY : temp

IMPLICIT NONE
INTEGER, INTENT(IN) :: npress, itemp
REAL(DP), INTENT(IN) :: press(npress), b0t(npress), b0s(npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_press
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_press=find_free_unit()
   OPEN(UNIT=iu_press, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                              FORM='FORMATTED')

   WRITE(iu_press,'("#  ")')
   WRITE(iu_press,'("# p(kbar)        B_T(p) (kbar)        B_S(p) (kbar) &
                   &       B_S(p)-B_T(p) (kbar)  T=",f16.4," K")') temp(itemp)

   DO ipress = 1, npress
      WRITE(iu_press, '(e12.5,3e22.13)') press(ipress), &
                 b0t(ipress), b0s(ipress), b0s(ipress)-b0t(ipress)
   END DO
  
   CLOSE(iu_press)
ENDIF
RETURN
END SUBROUTINE write_bulk_anharm_t
!
!------------------------------------------------------------------------
SUBROUTINE write_dbulk_anharm(temp, b01t, b02t, ntemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the first (and second)
!   pressure derivative of the isothermal bulk modulus as a function of 
!   temperature for several pressures.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode
USE control_ev, ONLY : ieos
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), b01t(ntemp), b02t(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("#  ")')
   IF (ieos==2) THEN
      WRITE(iu_therm,'("#   T (K)",10x,"dB/dp (T)",10x,"d^2B/dp^2 (1/kbar)")')
      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(e12.5,2e23.13)') temp(itemp), b01t(itemp), &
                                                         b02t(itemp)
      ENDDO
   ELSE
      WRITE(iu_therm,'("#   T (K)          dB/dp (T) ")')
      DO itemp = 2, ntemp-1
         WRITE(iu_therm, '(e12.5,e23.13)') temp(itemp), b01t(itemp)
      END DO
   ENDIF 
   CLOSE(iu_therm)
ENDIF
RETURN
END SUBROUTINE write_dbulk_anharm
!
!------------------------------------------------------------------------
SUBROUTINE write_dbulk_anharm_t(press, b01, b02, npress, itemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the first (and second)
!   pressure derivative of the isothermal bulk modulus as a function of 
!   pressures for several temperatures.
!
!
USE kinds,       ONLY : DP
USE io_global,   ONLY : meta_ionode
USE control_ev,  ONLY : ieos
USE temperature, ONLY : temp
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress, itemp
REAL(DP), INTENT(IN) :: press(npress), b01(npress), b02(npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_press
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_press=find_free_unit()
   OPEN(UNIT=iu_press, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_press,'("#  ")')
   IF (ieos==2) THEN
      WRITE(iu_press,'("#  p (kbar)",10x,"dB/dp (T)",10x,&
                       &"d^2B/dp^2 (1/kbar), T=",f16.4," K")') temp(itemp)
      DO ipress = 1, npress
         WRITE(iu_press, '(e12.5,2e23.13)') press(ipress), b01(ipress), &
                                                           b02(ipress)
      ENDDO
   ELSE
      WRITE(iu_press,'("#   p (kbar)         dB/dp (T)   T=",f16.4," K")') &
                                                                 temp(itemp)
      DO ipress = 1, npress
         WRITE(iu_press, '(e12.5,e23.13)') press(ipress), b01(ipress)
      END DO
   ENDIF 
   CLOSE(iu_press)
ENDIF
RETURN
END SUBROUTINE write_dbulk_anharm_t
!
!------------------------------------------------------------------------
SUBROUTINE write_gamma_anharm(temp, gammat, cvt, beta, b0t, ntemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the average Gruneisen parameter,
!   the isochoric heat capacity and the product of the thermal expansion
!   and the bulk modulus as a function of temperature.
!
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), gammat(ntemp), cvt(ntemp), beta(ntemp), &
                                     b0t(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# ")')
   WRITE(iu_therm,'("# T (K)          gamma(T)             C_V(T) &
                              &(Ry/cell/K)    beta B_T (kbar/K)")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e22.13)') temp(itemp), gammat(itemp), &
                                         cvt(itemp), beta(itemp)*b0t(itemp)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_gamma_anharm
!
!------------------------------------------------------------------------
SUBROUTINE write_gamma_anharm_t(press, gammat, cvt, beta, b0t, npress, &
                                                        itemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the average Gruneisen parameter,
!   the isochoric heat capacity and the product of the thermal expansion
!   and the bulk modulus as a function of pressure for several 
!   temperatures.
!
USE kinds,       ONLY : DP
USE io_global,   ONLY : meta_ionode
USE temperature, ONLY : temp

IMPLICIT NONE
INTEGER,  INTENT(IN) :: npress, itemp
REAL(DP), INTENT(IN) :: press(npress), gammat(npress), cvt(npress), &
                        beta(npress), b0t(npress)
CHARACTER(LEN=*) :: filename

INTEGER :: ipress, iu_press
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_press=find_free_unit()
   OPEN(UNIT=iu_press, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_press,'("# ")')
   WRITE(iu_press,'("# p (kbar)       gamma(p)             C_V(p) &
                              &(Ry/cell/K)    beta B_T (kbar/K)")')

   DO ipress = 1, npress
      WRITE(iu_press, '(e12.5,3e22.13)') press(ipress), gammat(ipress), &
                                 cvt(ipress), beta(ipress)*b0t(ipress)
   ENDDO
   CLOSE(iu_press)
ENDIF

RETURN
END SUBROUTINE write_gamma_anharm_t
!
!------------------------------------------------------------------------
SUBROUTINE write_heat_anhar(temp, cet, cvt, cpt, ntemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the isochoric heat capacity,
!   the isobaric heat capacity and difference between the two
!   as a function of temperature.
!
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

   WRITE(iu_therm,'("# ")')
   WRITE(iu_therm,'("# T (K)        C_e(T) (Ry/cell/K)    (C_P-C_V)(T) &
                              &(Ry/cell/K)     C_e+C_P-C_V(T) (Ry/cell/K)")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,3e22.13)') temp(itemp), &
                   cet(itemp), cpt(itemp)-cvt(itemp), cet(itemp)+cpt(itemp)&
                                                    -cvt(itemp)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_heat_anhar
!
!------------------------------------------------------------------------
SUBROUTINE write_heat_anhar_t(press, cet, cvt, cpt, npress, &
                                                         itemp, filename)
!------------------------------------------------------------------------
!
!   This routine writes on output the isochoric heat capacity,
!   the isobaric heat capacity and difference between the two
!   as a function of pressure for several temperatures.
!
USE kinds,       ONLY : DP
USE io_global,   ONLY : meta_ionode
USE temperature, ONLY : temp

IMPLICIT NONE
INTEGER,  INTENT(IN) :: npress, itemp
REAL(DP), INTENT(IN) :: press(npress), cet(npress), cvt(npress), cpt(npress)
CHARACTER(LEN=*)     :: filename

INTEGER :: ipress, iu_press
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_press=find_free_unit()
   OPEN(UNIT=iu_press, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_press,'("# ")')
   WRITE(iu_press,'("# press(kbar)  C_e(p) (Ry/cell/K)    (C_P-C_V)(p) &
          &(Ry/cell/K)     C_e+C_P-C_V(p) (Ry/cell/K), T=",f16.4," K")') &
                                     temp(itemp)
   DO ipress = 1, npress
      WRITE(iu_press, '(e12.5,3e22.13)') press(ipress), &
           cet(ipress), cpt(ipress)-cvt(ipress), &
                       cet(ipress)+cpt(ipress)-cvt(ipress)
   ENDDO
   CLOSE(iu_press)
ENDIF

RETURN
END SUBROUTINE write_heat_anhar_t
!
!------------------------------------------------------------------------
SUBROUTINE write_el_anhar(temp, gammat, cpt, beta, b0t, ntemp, filename)
!------------------------------------------------------------------------
!
!  This routine writes on output the electronic contribution to the
!  average Gruneisen parameter, to the isobaric heat capacity, to the
!  thermal expansion and to the isothermal bulk modulus as a function
!  of temperature.
!
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER,  INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), gammat(ntemp), cpt(ntemp), beta(ntemp), &
                                     b0t(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# ")')
   WRITE(iu_therm,'("# T (K)          gamma(T)             C_P(T) &
                    &(Ry/cell/K)    betax10^6 (1/K)    B_T (kbar/K)")')

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,4e22.13)') temp(itemp), gammat(itemp), &
                           cpt(itemp), beta(itemp)*1.D6, b0t(itemp)
   ENDDO
   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_el_anhar
!
!------------------------------------------------------------------------
SUBROUTINE write_aux_grun(temp, betab, cp_grun_t, cv_grun_t, b0_grun_s, b0_t, &
                                                           ntemp, filename) 
!------------------------------------------------------------------------
!
!  This routine writes on output the anharmonic quantities computed
!  using the Gruneisen parameters: thermal expansion, difference
!  between isobaric and isochoric heat capacity, and difference
!  between isothermal and isobaric bulk moculus.
!
USE kinds,     ONLY : DP
USE io_global, ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), betab(ntemp), cp_grun_t(ntemp), &
                        cv_grun_t(ntemp), b0_grun_s(ntemp), b0_t(ntemp)

CHARACTER(LEN=*) :: filename

INTEGER :: itemp, iu_therm
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# gamma is the average gruneisen parameter ")')
   WRITE(iu_therm,'("#   T (K)         beta(T)x10^6 (1/K)&
            &   (C_p - C_v)(T)   (B_S - B_T) (T) (kbar) " )' )

   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(4e16.8)') temp(itemp), betab(itemp)*1.D6,    &
                        cp_grun_t(itemp)-cv_grun_t(itemp),   &
                        b0_grun_s(itemp) - b0_t(itemp)
   ENDDO
   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_aux_grun
!
!------------------------------------------------------------------------
SUBROUTINE set_volume_b0_grun()
!------------------------------------------------------------------------
!
!  This routines sets the volume, the celldm and the bulk modulus to
!  be used in the calculations with Gruneisen parameters. 
!
USE control_grun,       ONLY : lv0_t, lb0_t, vgrun_t, celldm_grun_t, &
                               b0_grun_t
USE control_thermo,     ONLY : ltherm_dos, ltherm_freq
USE temperature,        ONLY : ntemp
USE anharmonic,         ONLY : vmin_t, celldm_t, b0_t
USE ph_freq_anharmonic, ONLY : vminf_t, celldmf_t, b0f_t
USE equilibrium_conf,   ONLY : celldm0, omega0
USE control_mur,        ONLY : b0

IMPLICIT NONE

INTEGER :: itemp

DO itemp=1, ntemp
   IF (lv0_t.AND.ltherm_freq) THEN
      vgrun_t(itemp)=vminf_t(itemp)
      celldm_grun_t(:,itemp)=celldmf_t(:,itemp)
   ELSEIF (lv0_t.AND.ltherm_dos) THEN
      vgrun_t(itemp)=vmin_t(itemp)
      celldm_grun_t(:,itemp)=celldm_t(:,itemp)
   ELSE
      vgrun_t(itemp)=omega0
      celldm_grun_t(:,itemp)=celldm0(:)
   ENDIF

   IF (lb0_t.AND.ltherm_freq) THEN
      b0_grun_t(itemp)=b0f_t(itemp)
   ELSEIF(lb0_t.AND.ltherm_dos) THEN
      b0_grun_t(itemp)=b0_t(itemp)
   ELSE
      b0_grun_t(itemp)=b0
   ENDIF
ENDDO

RETURN
END SUBROUTINE set_volume_b0_grun
!
!  Interpolate the thermodynamic quantities at the temperature dependent
!  volume
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_harmonic()
!-----------------------------------------------------------------------
!
!   This subroutine interpolates the computed harmonic quantities at the
!   temperature dependent volume. This version is for quantities computed
!   from phonon dos.
!
USE kinds,          ONLY : DP
USE thermodynamics, ONLY : ph_e0, ph_ce, ph_b_fact, ph_ener, ph_free_ener, &
                                                             ph_entropy
USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_ce
USE control_emp_free_ener,  ONLY : add_empirical, emp_ener, emp_free_ener, &
                           emp_entr, emp_ce
USE el_anharmonic,  ONLY : el_energy_t, el_free_energy_t, el_entropy_t,    &
                           el_ce_t
USE emp_anharmonic,  ONLY : emp_energy_t, emp_free_energy_t, emp_entropy_t, &
                           emp_ce_t
USE anharmonic,     ONLY : vmin_t, celldm_t, cv_t, ce_t, ener_t,           &
                           free_ener_t, entropy_t, e0_t, bfact_t
USE control_thermo, ONLY : with_eigen
USE control_eldos,  ONLY : lel_free_energy

IMPLICIT NONE
!
!   here the vibrational energy, entropy, zero point energy 
!
CALL interpolate_thermo(vmin_t, celldm_t, ph_ener, ener_t)

CALL interpolate_thermo(vmin_t, celldm_t, ph_free_ener, free_ener_t)

CALL interpolate_thermo(vmin_t, celldm_t, ph_entropy, entropy_t)

CALL interpolate_e0(vmin_t, celldm_t, ph_e0, e0_t) 

CALL interpolate_thermo(vmin_t, celldm_t, ph_ce, ce_t) 
cv_t=ce_t

IF (with_eigen) CALL interpolate_b_fact(vmin_t, ph_b_fact, bfact_t)
!
!  If available interpolate the electronic part and add it
!  to the vibrational one
!
IF (lel_free_energy) THEN
   CALL interpolate_thermo(vmin_t, celldm_t, el_ener, el_energy_t)
   ener_t = ener_t + el_energy_t 
   CALL interpolate_thermo(vmin_t, celldm_t, el_free_ener, el_free_energy_t)
   free_ener_t = free_ener_t + el_free_energy_t 
   CALL interpolate_thermo(vmin_t, celldm_t, el_entr, el_entropy_t)
   entropy_t = entropy_t + el_entropy_t
   CALL interpolate_thermo(vmin_t, celldm_t, el_ce, el_ce_t) 
   ce_t=ce_t+el_ce_t
   cv_t=ce_t
ELSE
   el_energy_t=0.0_DP
   el_free_energy_t=0.0_DP
   el_entropy_t=0.0_DP
   el_ce_t=0.0_DP
ENDIF

IF (add_empirical) THEN
   CALL interpolate_thermo(vmin_t, celldm_t, emp_ener, emp_energy_t)
   ener_t = ener_t + emp_energy_t
   CALL interpolate_thermo(vmin_t, celldm_t, emp_free_ener, emp_free_energy_t)
   free_ener_t = free_ener_t + emp_free_energy_t
   CALL interpolate_thermo(vmin_t, celldm_t, emp_entr, emp_entropy_t)
   entropy_t = entropy_t + emp_entropy_t
   CALL interpolate_thermo(vmin_t, celldm_t, emp_ce, emp_ce_t)
   ce_t=ce_t+emp_ce_t
   cv_t=ce_t
ELSE
   emp_energy_t=0.0_DP
   emp_free_energy_t=0.0_DP
   emp_entropy_t=0.0_DP
   emp_ce_t=0.0_DP
ENDIF

RETURN
END SUBROUTINE interpolate_harmonic

!-----------------------------------------------------------------------
SUBROUTINE interpolate_harmonic_ph()
!-----------------------------------------------------------------------
!
!   This subroutine interpolates the computed harmonic quantities at the
!   temperature dependent volume. This version is for quantities computed
!   from Brillouin zone integration.
!
USE kinds,    ONLY : DP
USE ph_freq_thermodynamics, ONLY : phf_e0, phf_ce, phf_b_fact, phf_ener, &
                                   phf_free_ener, phf_entropy
USE el_thermodynamics,  ONLY : el_ener, el_free_ener, el_entr, el_ce
USE control_emp_free_ener,  ONLY : add_empirical, emp_ener, emp_free_ener, &
                               emp_entr, emp_ce
USE ph_freq_anharmonic, ONLY : vminf_t, celldmf_t, cvf_t, cef_t, enerf_t,&
                               free_enerf_t, entropyf_t, e0f_t, bfactf_t
USE el_anharmonic,  ONLY : el_energyf_t, el_free_energyf_t, el_entropyf_t, &
                           el_cef_t
USE emp_anharmonic,  ONLY : emp_energyf_t, emp_free_energyf_t, &
                           emp_entropyf_t, emp_cef_t
USE control_thermo, ONLY : with_eigen
USE control_eldos,  ONLY : lel_free_energy

IMPLICIT NONE
!
!   interpolate the vibrational energy, entropy, zero point energy,
!   heat capacity and b-factor
!
CALL interpolate_thermo(vminf_t, celldmf_t, phf_ener, enerf_t) 

CALL interpolate_thermo(vminf_t, celldmf_t, phf_free_ener, free_enerf_t)

CALL interpolate_thermo(vminf_t, celldmf_t, phf_entropy, entropyf_t)

CALL interpolate_e0(vminf_t, celldmf_t, phf_e0, e0f_t)

CALL interpolate_thermo(vminf_t, celldmf_t, phf_ce, cef_t) 
cvf_t=cef_t

IF (with_eigen) CALL interpolate_b_fact(vminf_t, phf_b_fact, bfactf_t)
!
!  If available interpolate the electronic part and add it
!  to the vibrational one
!
IF (lel_free_energy) THEN
   CALL interpolate_thermo(vminf_t, celldmf_t, el_ener, el_energyf_t)
   enerf_t = enerf_t + el_energyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, el_free_ener, el_free_energyf_t)
   free_enerf_t = free_enerf_t + el_free_energyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, el_entr, el_entropyf_t)
   entropyf_t = entropyf_t + el_entropyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, el_ce, el_cef_t)
   cef_t = cef_t + el_cef_t
   cvf_t=cef_t
ELSE
   el_energyf_t=0.0_DP
   el_free_energyf_t=0.0_DP
   el_entropyf_t=0.0_DP
   el_cef_t=0.0_DP
ENDIF

IF (add_empirical) THEN
   CALL interpolate_thermo(vminf_t, celldmf_t, emp_ener, emp_energyf_t)
   enerf_t = enerf_t + emp_energyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, emp_free_ener, emp_free_energyf_t)
   free_enerf_t = free_enerf_t + emp_free_energyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, emp_entr, emp_entropyf_t)
   entropyf_t = entropyf_t + emp_entropyf_t
   CALL interpolate_thermo(vminf_t, celldmf_t, emp_ce, emp_cef_t)
   cef_t=cef_t+emp_cef_t
   cvf_t=cef_t
ELSE
   emp_energyf_t=0.0_DP
   emp_free_energyf_t=0.0_DP
   emp_entropyf_t=0.0_DP
   emp_cef_t=0.0_DP
ENDIF

RETURN
END SUBROUTINE interpolate_harmonic_ph
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_harmonic_pt()
!-----------------------------------------------------------------------
!
!   This subroutine interpolates the computed harmonic quantities at the
!   temperature and pressure dependent volume. This version is for 
!   quantities computed from phonon dos. Presently only the heat capacity
!   is interpolated.
!
USE kinds,          ONLY : DP
USE thermodynamics, ONLY : ph_ener, ph_free_ener, ph_entropy, ph_ce
USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_ce
USE control_emp_free_ener, ONLY : add_empirical, emp_ener, emp_ce, emp_entr
USE el_anharmonic,  ONLY : el_ener_pt, el_free_ener_pt, el_entr_pt, el_ce_pt
USE emp_anharmonic, ONLY : emp_free_ener_pt, emp_ener_pt, emp_ce_pt, &
                           emp_entr_pt
USE anharmonic,     ONLY : celldm_t ! this is not used yet, must become pt
USE anharmonic_pt,  ONLY : vmin_pt, ener_pt, free_ener_pt, entr_pt, cv_pt, &
                           ce_pt
USE control_eldos,  ONLY : lel_free_energy
USE control_pressure, ONLY : npress_plot

IMPLICIT NONE
INTEGER :: ipressp

IF (npress_plot==0) RETURN

DO ipressp=1,npress_plot
   CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, ph_ener, &
                                                         ener_pt(:,ipressp))
   CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, ph_free_ener, &
                                                    free_ener_pt(:,ipressp))
   CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, ph_entropy, &
                                                         entr_pt(:,ipressp))
   CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, ph_ce, &
                                                         ce_pt(:,ipressp))
ENDDO
cv_pt=ce_pt
!
!  If available interpolate the electronic part and add it
!  to the vibrational one
!
IF (lel_free_energy) THEN
   DO ipressp=1,npress_plot
      CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, el_ener, &
                                                    el_ener_pt(:,ipressp))
      CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, el_free_ener, &
                                                 el_free_ener_pt(:,ipressp))
      CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, el_entr, &
                                                    el_entr_pt(:,ipressp))
      CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, el_ce, &
                                                    el_ce_pt(:,ipressp))
   ENDDO
   free_ener_pt=free_ener_pt+el_free_ener_pt
   ener_pt=ener_pt+el_ener_pt
   entr_pt=entr_pt+el_entr_pt
   ce_pt=ce_pt+el_ce_pt
   cv_pt=ce_pt
ELSE
   el_free_ener_pt=0.0_DP
   el_ener_pt=0.0_DP
   el_entr_pt=0.0_DP
   el_ce_pt=0.0_DP
ENDIF
!
!  If requested interpolate the empirical part and add it
!  to the vibrational one
!
IF (add_empirical) THEN
   DO ipressp=1,npress_plot
      CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, emp_ener, &
                                                    emp_ener_pt(:,ipressp))
      CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, emp_entr, &
                                                    emp_entr_pt(:,ipressp))
      CALL interpolate_thermo(vmin_pt(:,ipressp), celldm_t, emp_ce, &
                                                    emp_ce_pt(:,ipressp))
   ENDDO
   ener_pt=ener_pt+emp_ener_pt
   entr_pt=entr_pt+emp_entr_pt
   ce_pt=ce_pt+emp_ce_pt
   cv_pt=ce_pt
ELSE
   emp_ener_pt=0.0_DP
   emp_entr_pt=0.0_DP
   emp_ce_pt=0.0_DP
ENDIF

RETURN
END SUBROUTINE interpolate_harmonic_pt
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_harmonic_noe_t()
!-----------------------------------------------------------------------
!
!   This subroutine interpolates the computed harmonic quantities at the
!   temperature and pressure dependent volume neglecting the effect of
!   electron free energy. This version is for quantities computed from 
!   phonon dos. Presently only the heat capacity is interpolated.
!
USE kinds,          ONLY : DP
USE thermodynamics, ONLY : ph_ce
USE anharmonic,     ONLY : celldm_t ! this is not used yet, must become pt
USE anharmonic,     ONLY : vmin_noe_t, cv_noe_t, ce_noe_t
USE control_eldos,  ONLY : lel_free_energy

IMPLICIT NONE

IF (.NOT.lel_free_energy) RETURN

CALL interpolate_thermo(vmin_noe_t, celldm_t, ph_ce, ce_noe_t)
cv_noe_t=ce_noe_t

RETURN
END SUBROUTINE interpolate_harmonic_noe_t
!
!-----------------------------------------------------------------------
SUBROUTINE interpolate_harmonic_ptt()
!-----------------------------------------------------------------------
!
!   This subroutine interpolates the computed harmonic quantities at the
!   temperature and pressure dependent volume. This version is for 
!   quantities computed from phonon dos. Presently only the heat capacity
!   is interpolated.
!
USE kinds,             ONLY : DP
USE thermodynamics,    ONLY : ph_ce
USE el_thermodynamics, ONLY : el_ce
USE el_anharmonic,     ONLY : el_ce_ptt
USE control_emp_free_ener, ONLY : add_empirical, emp_ce
USE anharmonic,        ONLY : celldm_t ! this is not used yet, must become ptt
USE anharmonic_ptt,    ONLY : vmin_ptt, cv_ptt, ce_ptt
USE emp_anharmonic,    ONLY : emp_ce_ptt
USE temperature,       ONLY : ntemp_plot, itemp_plot
USE control_eldos,     ONLY : lel_free_energy

IMPLICIT NONE
INTEGER :: itempp, itemp

IF (ntemp_plot==0) RETURN

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   CALL interpolate_thermo_p(vmin_ptt(:,itempp), ph_ce, ce_ptt(:,itempp), &
                                                                   itemp)
ENDDO
cv_ptt=ce_ptt
!
!  If available interpolate the electronic part and add it
!  to the vibrational one
!
IF (lel_free_energy) THEN
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_thermo_p(vmin_ptt(:,itempp), el_ce, &
                                              el_ce_ptt(:,itempp), itemp)
   ENDDO
   ce_ptt=ce_ptt+el_ce_ptt
   cv_ptt=ce_ptt
ELSE
   el_ce_ptt=0.0_DP
ENDIF
!
!  If requested in input interpolate the empirical part and add it
!  to the vibrational one
!
IF (add_empirical) THEN
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_thermo_p(vmin_ptt(:,itempp), emp_ce, &
                                              emp_ce_ptt(:,itempp), itemp)
   ENDDO
   ce_ptt=ce_ptt+emp_ce_ptt
   cv_ptt=ce_ptt
ELSE
   emp_ce_ptt=0.0_DP
ENDIF

RETURN
END SUBROUTINE interpolate_harmonic_ptt

!-----------------------------------------------------------------------
SUBROUTINE write_anhar_p()
!-----------------------------------------------------------------------
!
!   This routine writes on output the anharmonic quantities calculated
!   at the volume that minimize the gibbs free energy (computed from 
!   phonon dos).
!   In this routine the pressure is a parameter and the quantities
!   are calculated as a function of temperature at selected pressures
!
USE kinds,          ONLY : DP
USE temperature,    ONLY : ntemp, temp, itemp300
USE anharmonic_pt,  ONLY : vmin_pt, emin_pt, beta_pt, b0_pt, b01_pt, &
                           b02_pt, cv_pt, ce_pt, cp_pt, gamma_pt, b0_s_pt,  &
                           free_ener_pt, ener_pt, entr_pt
USE anharmonic,     ONLY : noelcvg, e0_t  
USE control_eldos,  ONLY : lel_free_energy
USE data_files,     ONLY : flanhar
USE control_pressure, ONLY : press, ipress_plot, npress_plot

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, ipress, ipressp
LOGICAL :: subtract_el

REAL(DP) :: aux(ntemp)
!

DO ipressp=1,npress_plot
!
!  Compute beta at several pressures if required in input
!  vmin_pt must have been computed before
!
   CALL compute_beta(vmin_pt(:,ipressp), beta_pt(:,ipressp), temp, ntemp)

   ipress=ipress_plot(ipressp)

   IF (itemp300 > 0) THEN
      aux(:)=vmin_pt(:,ipressp)/vmin_pt(itemp300,ipressp)
   ELSE
      aux(:)=vmin_pt(:,ipressp)
   ENDIF
!
!  Compute the isobaric heat capacity, the isoentropic bulk modulus and
!  the Gruneisen parameter. ce_pt must have been interpolated outside
!  the routine
!
   subtract_el=(lel_free_energy.AND.noelcvg)
   CALL compute_cp_bs_g(beta_pt(:,ipressp), vmin_pt(:,ipressp),        &
                 b0_pt(:,ipressp), ce_pt(:,ipressp), cp_pt(:,ipressp), &
                 b0_s_pt(:,ipressp), gamma_pt(:,ipressp), subtract_el)
!
!  write on output the temperature, the volume, the free energy and beta.
!
   filename="anhar_files/"//TRIM(flanhar)//'.press'
   CALL add_value(filename, press(ipress))
   CALL write_ener_beta_vol(temp, vmin_pt(:,ipressp), aux(:), &
                  emin_pt(:,ipressp), beta_pt(:,ipressp), ntemp, filename)

!
!  The bulk modulus and its derivative as a function of temperature at
!  several pressures
!
   filename="anhar_files/"//TRIM(flanhar)//'.bulk_press'
   CALL add_value(filename, press(ipress))

   CALL write_bulk_anharm(temp, b0_pt(:,ipressp), b0_s_pt(:,ipressp), &
                                                       ntemp, filename)
   filename="anhar_files/"//TRIM(flanhar)//'.dbulk_press'
   CALL add_value(filename, press(ipress))

   CALL write_dbulk_anharm(temp, b01_pt(:,ipressp), b02_pt(:,ipressp), &
                                                    ntemp, filename)
!
!   here the heat capacities (constant strain and constant stress)
!
   filename="anhar_files/"//TRIM(flanhar)//'.heat_press'
   CALL add_value(filename, press(ipress))

   CALL write_heat_anhar(temp, ce_pt(:,ipressp), ce_pt(:,ipressp), &
                cp_pt(:,ipressp), ntemp, filename)
!
!   here the Gruneisen parameter and the terms that form it
!
   filename="anhar_files/"//TRIM(flanhar)//'.gamma_press'
   CALL add_value(filename, press(ipress))

   CALL write_gamma_anharm(temp, gamma_pt(:,ipressp), ce_pt(:,ipressp), &
              beta_pt(:,ipressp), b0_pt(:,ipressp), ntemp, filename)
!
!   here all the interpolated harmonic quantities
!
   filename="anhar_files/"//TRIM(flanhar)//'.therm_press'
   CALL add_value(filename, press(ipress))
   CALL write_thermo_anharm(temp, ntemp, e0_t, ener_pt(:,ipressp), &
             free_ener_pt(:,ipressp), entr_pt(:,ipressp), cv_pt(:,ipressp), &
             filename)

ENDDO

RETURN
END SUBROUTINE write_anhar_p
!
!-------------------------------------------------------------------------
SUBROUTINE write_anhar_t()
!-------------------------------------------------------------------------
!
!  This routine computes the volume thermal expansion, the bulk modulus,
!  the heat capacity and the average Gruneisen parameter as a function 
!  of pressure for selected temperatures.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE data_files,     ONLY : flanhar
USE temperature,    ONLY : temp, deltat, itemp_plot, ntemp_plot, itemp300
USE control_eldos,  ONLY : lel_free_energy
USE anharmonic,     ONLY : vmin_t, a_t, noelcvg
USE anharmonic_ptt, ONLY : beta_ptt, vmin_ptt, b0_ptt, b01_ptt, b02_ptt, &
                           gamma_ptt, ce_ptt, cp_ptt, b0_s_ptt, emin_ptt
USE el_anharmonic,  ONLY : el_ce_ptt
USE anharmonic_vt,  ONLY : press_vt
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_pressure,  ONLY : press, npress
USE control_mur,    ONLY : vmin, b0, b01, b02, emin
USE control_mur_p,  ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
USE io_global,      ONLY : meta_ionode
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE

CHARACTER(LEN=256) :: filename
LOGICAL :: subtract_el
REAL(DP) :: vmm1, vmp1
REAL(DP), ALLOCATABLE :: aux(:,:)
INTEGER :: m1, ipress, itemp, itempp, startp, lastp

IF (ntemp_plot==0) RETURN

m1=poly_degree_ph+1

ALLOCATE(aux(npress,ntemp_plot))
press_vt=0.0_DP
beta_ptt=0.0_DP
cp_ptt=0.0_DP
b0_s_ptt=0.0_DP
gamma_ptt=0.0_DP

CALL divide(world_comm, npress, startp, lastp)

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   CALL write_mur_pol(vmin, b0, b01, b02, emin, a_t(:,itemp), m1, itempp)
   CALL write_thermal_press(a_t(:,itemp),m1, itempp)

   DO ipress=startp,lastp
      CALL find_min_mur_pol(vmin_p(ipress), b0_p(ipress) / ry_kbar,    &
                         b01_p(ipress), b02_p(ipress)* ry_kbar,        &
                         a_t(:,itemp-1), m1, vmm1)
      CALL find_min_mur_pol(vmin_p(ipress), b0_p(ipress) / ry_kbar,    &
                         b01_p(ipress), b02_p(ipress) * ry_kbar,       &
                         a_t(:,itemp+1), m1, vmp1)
      beta_ptt(ipress, itempp) = (vmp1 - vmm1) / 2.0_DP / deltat /     &
                                        vmin_ptt(ipress,itempp)
   ENDDO 
ENDDO
CALL mp_sum(press_vt, world_comm)
CALL mp_sum(beta_ptt, world_comm)

subtract_el=(lel_free_energy.AND.noelcvg)
DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   CALL compute_cp_bs_gp(beta_ptt(:,itempp), vmin_ptt(:,itempp),      &
              b0_ptt(:,itempp), ce_ptt(:,itempp), cp_ptt(:,itempp),   &
              b0_s_ptt(:,itempp), gamma_ptt(:,itempp), el_ce_ptt(:,itempp), &
                      itemp, subtract_el)
   IF (itemp300 > 0) THEN
      aux(:,itempp)=vmin_ptt(:,itempp)/vmin_t(itemp300)
   ELSE
      aux(:,itempp)=vmin_ptt(:,itempp)
   ENDIF
ENDDO

IF (meta_ionode) THEN
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
!
!   volume, energy, beta
!
      filename="anhar_files/"//TRIM(flanhar)//'.temp'
      CALL add_value(filename, temp(itemp))
      CALL write_ener_beta_t(press, vmin_ptt(:,itempp), aux(:,itempp), & 
             emin_ptt(:,itempp), beta_ptt(:,itempp), npress, itemp, filename)
!
!   bulk modulus
!
      filename="anhar_files/"//TRIM(flanhar)//'.bulk_temp'
      CALL add_value(filename, temp(itemp))
      CALL write_bulk_anharm_t(press, b0_ptt(:,itempp), b0_s_ptt(:,itempp),  &
                                                     npress, itemp, filename)
!
!   pressure derivative of the bulk modulus
!

      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp'
      CALL add_value(filename, temp(itemp))
      CALL write_dbulk_anharm_t(press, b01_ptt(:,itempp), b02_ptt(:,itempp),  &
                                                     npress, itemp, filename)
!
!   heat capacities 
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat_temp'
      CALL add_value(filename, temp(itemp))

      CALL write_heat_anhar_t(press, ce_ptt(:,itempp), ce_ptt(:,itempp), &
              cp_ptt(:,itempp), npress, itemp, filename)
!
!   Average Gruneisen parameter
!
      filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp'
      CALL add_value(filename, temp(itemp))

      CALL write_gamma_anharm_t(press, gamma_ptt(:,itempp), &
            ce_ptt(:,itempp), beta_ptt(:,itempp), b0_ptt(:,itempp), &
                                           npress, itemp, filename)
   ENDDO
ENDIF

DEALLOCATE(aux)

RETURN
END SUBROUTINE write_anhar_t
!
!-----------------------------------------------------------------------
SUBROUTINE write_anhar_v()
!-----------------------------------------------------------------------
!
!  This routine writes on output the thermal pressure as a function of
!  the temperature for several volumes
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : omega_geo
USE control_vol,    ONLY : nvol_plot, ivol_plot
USE temperature,    ONLY : temp, ntemp
USE anharmonic_vt,  ONLY : press_vtt
USE data_files,     ONLY : flanhar
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8) :: float_to_char
INTEGER :: find_free_unit, iu_therm
INTEGER :: itemp, ivol, ivolp
REAL(DP) :: omega

DO ivolp=1,nvol_plot
   ivol=ivol_plot(ivolp)

   omega=omega_geo(ivol)

   IF (meta_ionode) THEN
      iu_therm=find_free_unit()
      filename="anhar_files/"//TRIM(flanhar)//'.vol.'//& 
                            TRIM(float_to_char(omega,2))
      OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                                        FORM='FORMATTED')
      WRITE(iu_therm,'("#",5x,"  T (K) ", 16x, " pressure (kbar)",10x, &
               &"pressure(T=0 K) (kbar)",5x,"at V=", f15.6, " (a.u.)^3")') &
                                                                     omega
      DO itemp=1,ntemp
         WRITE(iu_therm,"(f20.8,2f25.12)") temp(itemp), &
                                  press_vtt(itemp,ivolp), press_vtt(1,ivolp)         
      ENDDO

      CLOSE(unit=iu_therm, STATUS='KEEP')
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_anhar_v

!-----------------------------------------------------------------------
SUBROUTINE write_anhar_el()
!-----------------------------------------------------------------------
!
!  This routine writes on output the electronic harmonic quantities 
!  interpolated at the temperature dependent volume.
!
USE kinds,       ONLY : DP
USE data_files,  ONLY : flanhar
USE temperature, ONLY : temp, ntemp
USE io_global,   ONLY : meta_ionode
USE el_anharmonic,  ONLY : el_energy_t, el_free_energy_t, el_entropy_t,    &
                           el_ce_t
!
IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm

INTEGER :: find_free_unit

filename="anhar_files/"//TRIM(flanhar)//".el_press"
CALL add_pressure(filename)

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("#",5x,"   T (K)", 10x, " Cv (Ry/K)", 9x,    &
                  &" energy (Ry)", 5x, "  free energy (Ry)", 5x, &
                  &" entropy (Ry/K)"  )') 
   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,4e20.12)') temp(itemp), &
             el_ce_t(itemp), el_energy_t(itemp), el_free_energy_t(itemp), &
             el_entropy_t(itemp)
   END DO

   CLOSE(iu_therm)
ENDIF

RETURN
END SUBROUTINE write_anhar_el
!
!-----------------------------------------------------------------------
SUBROUTINE write_anhar_el_pt()
!-----------------------------------------------------------------------
!
!  This routine writes on output the electronic harmonic quantities 
!  interpolated at the temperature and pressure dependent volume.
!  The quantities are written as a function of temperature for 
!  several pressures. Presently only the heat capacity is interpolated
!
USE kinds,       ONLY : DP
USE data_files,  ONLY : flanhar
USE temperature, ONLY : temp, ntemp
USE io_global,   ONLY : meta_ionode
USE el_anharmonic,  ONLY : el_ce_pt
USE control_pressure, ONLY : press, npress_plot, ipress_plot, pressure_kb
!
IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: itemp, iu_therm, ipress, ipressp

INTEGER :: find_free_unit

DO ipressp=1,npress_plot
   ipress=ipress_plot(ipressp)
   IF (press(ipress)==pressure_kb) CYCLE
   filename="anhar_files/"//TRIM(flanhar)//".el_press"
   CALL add_value(filename,press(ipress))
   CALL add_pressure(filename)

   IF (meta_ionode) THEN
      iu_therm=find_free_unit()
      OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                                  FORM='FORMATTED')
      WRITE(iu_therm,'("#",5x,"   T (K)", 10x, " C_v (Ry/K)")')
      DO itemp = 1, ntemp
         WRITE(iu_therm, '(e16.8,e20.12)') temp(itemp), el_ce_pt(itemp,ipressp)
      ENDDO

      CLOSE(iu_therm)
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_anhar_el_pt
!
!-----------------------------------------------------------------------
SUBROUTINE fit_free_energy()
!-----------------------------------------------------------------------
!
!  This subroutine fits the vibrational (and possibily electronic) free
!  energy as a function of the volume with a polynomial of degree 
!  poly_degree_ph. This is done for all the mesh of temperatures.
!  This version use the free energy computed from phonon dos.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_eldos,  ONLY : lel_free_energy
USE anharmonic,     ONLY : a_t
USE thermodynamics, ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE control_emp_free_ener, ONLY : add_empirical, emp_free_ener
USE temperature,    ONLY : ntemp
USE polyfit_mod,    ONLY : polyfit
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER :: idata, ndata, itemp, startt, lastt
REAL(DP) :: x(ngeo(1)), y(ngeo(1))

a_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp = startt, lastt
   ndata=0
   DO idata=1,ngeo(1)
      IF (no_ph(idata)) CYCLE
      ndata=ndata+1
      x(ndata)=omega_geo(idata)
      y(ndata)=ph_free_ener(itemp,idata)
      IF (lel_free_energy) y(ndata)=y(ndata)+el_free_ener(itemp,idata)
      IF (add_empirical) y(ndata)=y(ndata)+emp_free_ener(itemp,idata)
!      IF (MOD(itemp,100)==0) WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
   ENDDO
   CALL polyfit(x, y, ndata, a_t(:,itemp), poly_degree_ph)
ENDDO
CALL mp_sum(a_t, world_comm)

RETURN
END SUBROUTINE fit_free_energy
!
!-----------------------------------------------------------------------
SUBROUTINE fit_free_energy_noe()
!-----------------------------------------------------------------------
!
!  This subroutine fits the vibrational free energy as a function of 
!  the volume with a polynomial of degree poly_degree_ph. 
!  This is done for all the mesh of temperatures/
!  This routine is used only when lel_free_energy=.TRUE., because 
!  when lel_free_energy=.FALSE. a_noe_t would be equal to a_t.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_eldos,  ONLY : lel_free_energy
USE anharmonic,     ONLY : a_noe_t
USE thermodynamics, ONLY : ph_free_ener
USE temperature,    ONLY : ntemp
USE polyfit_mod,    ONLY : polyfit
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER :: idata, ndata, itemp, startt, lastt
REAL(DP) :: x(ngeo(1)), y(ngeo(1))

IF (.NOT.lel_free_energy) RETURN

a_noe_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp = startt, lastt
   ndata=0
   DO idata=1,ngeo(1)
      IF (no_ph(idata)) CYCLE
      ndata=ndata+1
      x(ndata)=omega_geo(idata)
      y(ndata)=ph_free_ener(itemp,idata)
!      IF (MOD(itemp,100)==0) WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
   ENDDO
   CALL polyfit(x, y, ndata, a_noe_t(:,itemp), poly_degree_ph)
ENDDO
CALL mp_sum(a_noe_t, world_comm)

RETURN
END SUBROUTINE fit_free_energy_noe
!
!-----------------------------------------------------------------------
SUBROUTINE fit_free_energy_ph() 
!-----------------------------------------------------------------------
!
!  This subroutine fits the vibrational (and possibily electronic) free
!  energy as a function of the volume with a polynomial of degree 
!  poly_degree_ph. This is done for all the mesh of temperatures.
!  This version use the free energy computed from Brillouin zone
!  integration.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph
USE control_quartic_energy, ONLY : poly_degree_ph
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics,      ONLY : el_free_ener
USE control_emp_free_ener, ONLY : add_empirical, emp_free_ener
USE ph_freq_anharmonic,     ONLY : af_t
USE control_eldos,  ONLY : lel_free_energy
USE temperature,    ONLY : ntemp
USE polyfit_mod,    ONLY : polyfit
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER :: idata, ndata, itemp, startt, lastt
REAL(DP) :: x(ngeo(1)), y(ngeo(1))

af_t=0.0_DP
CALL divide(world_comm, ntemp, startt, lastt)

DO itemp=startt, lastt
   ndata=0
   DO idata=1,ngeo(1)
      IF (no_ph(idata)) CYCLE
      ndata=ndata+1
      x(ndata)=omega_geo(idata)
      y(ndata)=phf_free_ener(itemp,idata) 
      IF (lel_free_energy) y(ndata)=y(ndata)+el_free_ener(itemp,idata)
      IF (add_empirical) y(ndata)=y(ndata)+emp_free_ener(itemp,idata)
!     WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
   ENDDO
   CALL polyfit(x, y, ndata, af_t(:,itemp), poly_degree_ph)
ENDDO
CALL mp_sum(af_t, world_comm)

RETURN
END SUBROUTINE fit_free_energy_ph
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_t()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus and
!  its derivative at a given temperature. It uses the free energy 
!  computed by phonon dos.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : celldm_geo, central_geo, omega_geo
USE control_mur,    ONLY : vmin, b0, b01, b02, emin
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, b02_t, free_e_min_t, a_t, &
                           celldm_t
USE temperature,    ONLY : ntemp
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER ::  m1, itemp, startt, lastt
REAL(DP) :: vm, aux, aux1, aux2, aux3

m1=poly_degree_ph+1

vmin_t=0.0_DP
b0_t=0.0_DP
b01_t=0.0_DP
b02_t=0.0_DP
free_e_min_t=0.0_DP
celldm_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp=startt,lastt
   CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, b02 * ry_kbar, &
                                                   a_t(:,itemp), m1, vm)
   CALL eos_energy_pol(ieos, vm, aux3, vmin, b0/ry_kbar, b01, b02*ry_kbar, &
                                                   a_t(:,itemp), m1)
   CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin, b0/ry_kbar, b01, &
                                    b02 * ry_kbar, a_t(:,itemp), m1)
   vmin_t(itemp)=vm
   b0_t(itemp)=aux * ry_kbar
   b01_t(itemp)=aux1 
   b02_t(itemp)=aux2 / ry_kbar
   free_e_min_t(itemp)=emin+ aux3
   CALL compute_celldm_geo(vm, celldm_t(1,itemp), &
                   celldm_geo(1,central_geo), omega_geo(central_geo))
ENDDO

CALL mp_sum(vmin_t, world_comm)
CALL mp_sum(b0_t,   world_comm)
CALL mp_sum(b01_t,  world_comm)
CALL mp_sum(b02_t,  world_comm)
CALL mp_sum(free_e_min_t, world_comm)
CALL mp_sum(celldm_t, world_comm)

RETURN
END SUBROUTINE anhar_ev_t
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_noe_t()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus, and
!  its derivative at a given temperature. It uses the free energy 
!  computed by phonon dos. The electronic part is not added.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE control_mur,    ONLY : vmin, b0, b01, b02, emin
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE anharmonic,     ONLY : vmin_noe_t, b0_noe_t, b01_noe_t, b02_noe_t, &
                           free_e_min_noe_t, a_noe_t
USE temperature,    ONLY : ntemp
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER ::  m1, itemp, startt, lastt
REAL(DP) :: vm, aux, aux1, aux2, aux3

m1=poly_degree_ph+1

vmin_noe_t=0.0_DP
b0_noe_t=0.0_DP
b01_noe_t=0.0_DP
b02_noe_t=0.0_DP
free_e_min_noe_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp=startt,lastt
   CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, b02 * ry_kbar, &
                                                   a_noe_t(:,itemp), m1, vm)
   CALL eos_energy_pol(ieos, vm, aux3, vmin, b0/ry_kbar, b01, b02*ry_kbar, &
                                                   a_noe_t(:,itemp), m1)
   CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin, b0/ry_kbar, b01, &
                                    b02 * ry_kbar, a_noe_t(:,itemp), m1)
   vmin_noe_t(itemp)=vm
   b0_noe_t(itemp)=aux * ry_kbar
   b01_noe_t(itemp)=aux1 
   b02_noe_t(itemp)=aux2 / ry_kbar
   free_e_min_noe_t(itemp)=emin+ aux3
ENDDO

CALL mp_sum(vmin_noe_t, world_comm)
CALL mp_sum(b0_noe_t,   world_comm)
CALL mp_sum(b01_noe_t,  world_comm)
CALL mp_sum(b02_noe_t,  world_comm)
CALL mp_sum(free_e_min_noe_t, world_comm)

RETURN
END SUBROUTINE anhar_ev_noe_t
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_t_ph() 
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!  at a given temperature. It uses the free energy computed by the
!  Brillouin zone integration.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE ph_freq_anharmonic,     ONLY : vminf_t, b0f_t, b01f_t, b02f_t, &
                                   free_e_minf_t, af_t
USE control_mur,    ONLY : emin, vmin, b0, b01, b02
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE temperature,    ONLY : ntemp
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER :: m1, itemp, startt, lastt
REAL(DP) :: vm, aux, aux1, aux2, aux3

m1=poly_degree_ph+1

vminf_t=0.0_DP
b0f_t=0.0_DP
b01f_t=0.0_DP
b02f_t=0.0_DP
free_e_minf_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp=startt,lastt
   CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, b02*ry_kbar, &
                                        af_t(:,itemp), m1, vm)
   CALL eos_energy_pol(ieos, vm, aux3, vmin, b0/ry_kbar, b01, b02*ry_kbar, &
                                        af_t(:,itemp), m1)
   CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin, b0/ry_kbar, &
                              b01, b02 * ry_kbar, af_t(:,itemp), m1)

   vminf_t(itemp)=vm
   b0f_t(itemp)=aux * ry_kbar
   b01f_t(itemp)=aux1
   b02f_t(itemp)=aux2 / ry_kbar
   free_e_minf_t(itemp)=emin + aux3
   !
ENDDO
CALL mp_sum(vminf_t, world_comm)
CALL mp_sum(b0f_t, world_comm)
CALL mp_sum(b01f_t, world_comm)
CALL mp_sum(b02f_t, world_comm)
CALL mp_sum(free_e_minf_t, world_comm)

RETURN
END SUBROUTINE anhar_ev_t_ph
!
!-----------------------------------------------------------------------
SUBROUTINE summarize_anhar_param()
!-----------------------------------------------------------------------
!
!  This routine writes on output the parameters of the equation of state
!  computed for all temperatures (phonon dos version).
!
USE kinds,        ONLY : DP
USE thermo_mod,   ONLY : central_geo, omega_geo, celldm_geo
USE anharmonic,   ONLY : vmin_t, b0_t, b01_t, b02_t, free_e_min_t
USE temperature,  ONLY : ntemp, temp
USE io_global,    ONLY : stdout

IMPLICIT NONE
INTEGER  :: itemp
REAL(DP) :: celldm_(6)

!
!   This routine is for debugging purposes only. Comment the RETURN if you
!   want a long output
!
RETURN
DO itemp=1,ntemp
   CALL compute_celldm_geo(vmin_t(itemp), celldm_, &
                        celldm_geo(1,central_geo), omega_geo(central_geo))
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x, "free energy from phonon dos, at T= ", f12.6)') &
                                                                temp(itemp)
   CALL summarize_mur(celldm_(1), b0_t(itemp), b01_t(itemp), b02_t(itemp), &
                                                    free_e_min_t(itemp))

   WRITE(stdout,'(2x,76("-"),/)')
ENDDO

RETURN
END SUBROUTINE summarize_anhar_param
!
!-----------------------------------------------------------------------
SUBROUTINE summarize_anhar_param_ph()
!-----------------------------------------------------------------------
!
!  This routine writes on output the parameters of the equation of state
!  computed for all temperatures (Brillouin zone integration version).
!
USE kinds,        ONLY : DP
USE thermo_mod,   ONLY : central_geo, omega_geo, celldm_geo
USE ph_freq_anharmonic, ONLY : vminf_t, b0f_t, b01f_t, b02f_t, free_e_minf_t
USE temperature,  ONLY : ntemp, temp
USE io_global,    ONLY : stdout

IMPLICIT NONE
INTEGER  :: itemp
REAL(DP) :: celldm_(6)

!
!   This routine is for debugging purposes only. Comment the RETURN if you
!   want a long output
!
RETURN
DO itemp=1,ntemp
   CALL compute_celldm_geo(vminf_t(itemp), celldm_, &
                        celldm_geo(1,central_geo), omega_geo(central_geo))
   WRITE(stdout,'(/,2x,76("+"))')
   WRITE(stdout,'(5x, "free energy from phonon frequencies, &
                                       &at T= ", f12.6)') temp(itemp)
   CALL summarize_mur(celldm_(1), b0f_t(itemp), b01f_t(itemp), b02f_t(itemp), &
                                                    free_e_minf_t(itemp))

   WRITE(stdout,'(2x,76("+"),/)')
ENDDO

RETURN
END SUBROUTINE summarize_anhar_param_ph

!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_pt()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus, and its
!  derivatives as a function of temperature for the npress_plot pressures 
!  that we want to plot in output.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE control_pressure, ONLY : npress_plot, ipress_plot
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_mur_p,  ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
USE temperature,    ONLY : ntemp
USE anharmonic,     ONLY : a_t
USE anharmonic_pt,  ONLY : vmin_pt, b0_pt, b01_pt, b02_pt, emin_pt
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER  :: itemp, m1, ipress, ipressp, startt, lastt
REAL(DP) :: vm, aux, aux1, aux2, aux3

IF (npress_plot==0) RETURN

m1=poly_degree_ph+1

vmin_pt=0.0_DP
b0_pt=0.0_DP
b01_pt=0.0_DP
b02_pt=0.0_DP
emin_pt=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO ipressp=1, npress_plot
   ipress=ipress_plot(ipressp)
   DO itemp=startt,lastt
      CALL find_min_mur_pol(vmin_p(ipress), b0_p(ipress) / ry_kbar, &
                 b01_p(ipress), b02_p(ipress)*ry_kbar, a_t(:,itemp), m1, vm)

      CALL eos_energy_pol(ieos, vm, aux3, vmin_p(ipress), &
         b0_p(ipress)/ry_kbar, b01_p(ipress), b02_p(ipress)*ry_kbar, &
                                                  a_t(:,itemp), m1)

      CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin_p(ipress), &
         b0_p(ipress)/ry_kbar, b01_p(ipress), b02_p(ipress) * ry_kbar, &
                                                       a_t(:,itemp), m1)
      vmin_pt(itemp,ipressp)=vm
      b0_pt(itemp,ipressp)=aux * ry_kbar
      b01_pt(itemp,ipressp)=aux1
      b02_pt(itemp,ipressp)=aux2 / ry_kbar
      emin_pt(itemp,ipressp)=emin_p(ipress)+aux3
   ENDDO
   !
ENDDO

CALL mp_sum(vmin_pt, world_comm)
CALL mp_sum(b0_pt,   world_comm)
CALL mp_sum(b01_pt,  world_comm)
CALL mp_sum(b02_pt,  world_comm)
CALL mp_sum(emin_pt, world_comm)

RETURN
END SUBROUTINE anhar_ev_pt
!
!-------------------------------------------------------------------------
SUBROUTINE anhar_ev_ptt()
!-------------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus, and its
!  derivatives as a function of pressure for the ntemp_plot temperatures 
!  that we want to plot in output.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_ev,     ONLY : ieos
USE anharmonic,     ONLY : a_t
USE anharmonic_ptt, ONLY : vmin_ptt, emin_ptt, b0_ptt, b01_ptt, b02_ptt
USE control_quartic_energy, ONLY : poly_degree_ph
USE temperature,    ONLY : ntemp_plot, itemp_plot
USE control_pressure,  ONLY : npress
USE control_mur_p,  ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm


IMPLICIT NONE

REAL(DP) :: vm, aux, aux1, aux2, aux3
INTEGER :: m1, ipress, itempp, itemp, startp, lastp

IF (ntemp_plot==0) RETURN

m1=poly_degree_ph+1
vmin_ptt=0.0_DP
b0_ptt=0.0_DP
b01_ptt=0.0_DP
b02_ptt=0.0_DP
emin_ptt=0.0_DP
CALL divide(world_comm, npress, startp, lastp)

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   DO ipress=startp,lastp
      CALL find_min_mur_pol(vmin_p(ipress), b0_p(ipress) / ry_kbar, &
                         b01_p(ipress), b02_p(ipress) * ry_kbar,    &
                         a_t(:,itemp), m1, vm)
      CALL eos_energy_pol(ieos, vm, aux, vmin_p(ipress),            &
                        b0_p(ipress)/ry_kbar, b01_p(ipress),        &
                        b02_p(ipress) * ry_kbar, a_t(:,itemp), m1)
      CALL eos_bulk_pol(ieos, vm, aux1, aux2, aux3, vmin_p(ipress), &
                        b0_p(ipress)/ry_kbar, b01_p(ipress),        &
                        b02_p(ipress) * ry_kbar, a_t(:,itemp), m1)
      vmin_ptt(ipress,itempp) = vm
      b0_ptt(ipress,itempp) = aux1 * ry_kbar
      b01_ptt(ipress,itempp) = aux2
      b02_ptt(ipress,itempp) = aux3 / ry_kbar
      emin_ptt(ipress,itempp)= aux
   ENDDO
ENDDO

CALL mp_sum(vmin_ptt, world_comm)
CALL mp_sum(b0_ptt, world_comm)
CALL mp_sum(b01_ptt, world_comm)
CALL mp_sum(b02_ptt, world_comm)
CALL mp_sum(emin_ptt, world_comm)

RETURN
END SUBROUTINE anhar_ev_ptt
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_vt()
!-----------------------------------------------------------------------
!
!  This subroutine computes the thermal pressure as a function of 
!  temperature for the nvol_plot volumes required in input.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : omega_geo
USE control_quartic_energy, ONLY : poly_degree_ph
USE anharmonic,       ONLY : a_t
USE anharmonic_vt,    ONLY : press_vtt
USE temperature,      ONLY : ntemp
USE polyfit_mod,      ONLY : compute_poly_deriv
USE control_vol,      ONLY : nvol_plot, ivol_plot
USE mp,               ONLY : mp_sum
USE mp_world,         ONLY : world_comm

IMPLICIT NONE
REAL(DP) :: press, omega
INTEGER  :: itemp, ivolp, igeo, startt, lastt

IF (nvol_plot==0) RETURN

press_vtt=0.0_DP
CALL divide(world_comm, ntemp, startt, lastt)
DO ivolp=1,nvol_plot
   igeo = ivol_plot(ivolp)
   omega = omega_geo(igeo)
   DO itemp=startt,lastt
      CALL compute_poly_deriv(omega, poly_degree_ph, a_t(:,itemp), press)
      press = -press
      press_vtt(itemp,ivolp) = press * ry_kbar
   ENDDO
ENDDO
CALL mp_sum(press_vtt, world_comm)

RETURN
END SUBROUTINE anhar_ev_vt
!
!-----------------------------------------------------------------------
SUBROUTINE write_hugoniot()
!-----------------------------------------------------------------------
!
!  This subroutine computes for each volume in the mesh of volumes
!  the Hugoniot pressure and temperature if the temperature is within
!  the range of the temperature mesh.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE control_mur,      ONLY : vmin, b0, b01, b02, emin
USE control_quartic_energy, ONLY : poly_degree_ph
USE anharmonic,       ONLY : vmin_t, ener_t, a_t, celldm_t
USE control_ev,       ONLY : ieos
USE control_vol,      ONLY : vmin_input, deltav, nvol
USE thermodynamics,   ONLY : ph_ener
USE el_thermodynamics, ONLY : el_ener
USE el_anharmonic,    ONLY : el_energy_t
USE temperature,      ONLY : temp, ntemp, itemp300, deltat
USE control_eldos,    ONLY : lel_free_energy
USE control_quartic_energy, ONLY : poly_degree_ph
USE eos,              ONLY : eos_press_pol, eos_energy
USE data_files,       ONLY : flanhar
USE io_global,        ONLY : meta_ionode
USE mp,               ONLY : mp_sum
USE mp_world,         ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filename
REAL(DP) :: press, omega, ener0, ele300, t_hugo, p_hugo
REAL(DP), ALLOCATABLE :: omegat(:), enert(:), el_enert(:), presst(:), hugo(:)
INTEGER  :: m1, itemp, itemp_hugo, iu_hugo, ivol, startt, lastt
INTEGER  :: find_free_unit 

m1=poly_degree_ph+1

ALLOCATE(omegat(ntemp))
ALLOCATE(enert(ntemp))
ALLOCATE(el_enert(ntemp))
ALLOCATE(hugo(ntemp))
ALLOCATE(presst(ntemp))
!
!   Open the file where the Hugoniot pressure and temperature will be written
!
filename="anhar_files/"//TRIM(flanhar)//'.hugoniot'
CALL add_pressure(filename)
IF (meta_ionode) THEN
   iu_hugo=find_free_unit()
   OPEN(UNIT=iu_hugo, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_hugo,'("# Volume (a.u.)^3      Pressure (kbar)          T  (K)")')
ENDIF
!
!  For each volume compute the Hugoniot pressure and temperature
!
ele300=0.0_DP
el_enert=0.0_DP
CALL divide(world_comm, ntemp, startt, lastt)
DO ivol=1,nvol
   omega = vmin_input + deltav * (ivol-1)
   omegat=omega
!
!  Interpolate at the volume omega the cold equation of state
!
   CALL eos_energy(ieos, omega, ener0, vmin, b0/ry_kbar, b01, b02*ry_kbar)
!
!  For all temperatures interpolate the harmonic vibrational energy
!  and (if available) the electronic energy
!
   CALL interpolate_thermo(omegat, celldm_t, ph_ener, enert)
   IF (lel_free_energy) THEN
      CALL interpolate_thermo(omegat, celldm_t, el_ener, el_enert)
      ele300=el_energy_t(itemp300)
   ENDIF
!
!  Now compute the difference between 1/2 P_H (V_0-V) and E_H-E_0
!
   presst=0.0_DP
   hugo=0.0_DP
   DO itemp=startt,lastt
      CALL eos_press_pol(ieos, omega, presst(itemp), vmin, b0/ry_kbar, &
                                     b01, b02*ry_kbar, a_t(:,itemp), m1)

      hugo(itemp)= presst(itemp)*0.5_DP*(vmin_t(itemp300)-omega)  - &
               (ener0 + enert(itemp) + el_enert(itemp) - ener_t(itemp300) - &
                                         ele300)
   ENDDO
   CALL mp_sum(presst, world_comm)
   CALL mp_sum(hugo, world_comm)
!
!  Find the point in the temperature mesh in which difference changes sign. 
!   
   itemp_hugo=0
   DO itemp=1,ntemp
      IF (hugo(itemp)>0.0_DP) itemp_hugo=itemp
   ENDDO
!
!   Write the Hugoniot pressure and temperature only if the temperature
!   is in the range of the temperature mesh
!
   IF (meta_ionode) THEN
      IF (itemp_hugo > 0 .AND. itemp_hugo < ntemp) THEN
         t_hugo=temp(itemp_hugo) - deltat * hugo(itemp_hugo) /            &
                             (hugo(itemp_hugo+1)-hugo(itemp_hugo))
         p_hugo=presst(itemp_hugo) + (presst(itemp_hugo+1)-               &
                                               presst(itemp_hugo)) *      &
                             (t_hugo-temp(itemp_hugo)) / deltat
         WRITE(iu_hugo,'(3e20.10)') omega, p_hugo*ry_kbar, t_hugo
      ENDIF
   ENDIF
ENDDO
IF (meta_ionode) CLOSE (UNIT=iu_hugo, STATUS='KEEP')

DEALLOCATE(enert)
DEALLOCATE(el_enert)
DEALLOCATE(omegat)
DEALLOCATE(hugo)
DEALLOCATE(presst)

RETURN
END SUBROUTINE write_hugoniot
!
!-----------------------------------------------------------------------
SUBROUTINE find_min_mur_pol(v0, b0, b01, b02, a, m1, vm)
!-----------------------------------------------------------------------
!
!  This routine minimizes a function equal to an equation os state
!  (Birch-Murnaghan third or fourth order, or Munaghan)
!  with parameters v0, b0, b01, and b02 and a polynomial of degree m1-1 of
!  the form a(1) + a(2) * v + a(3) * v**2 +... a(m1) * v**(m1-1)
!  NB: b0 must be in atomic units not kbar
!  NB: b02 must be in atomic units not 1/kbar
!  The polynomial gives the energy in Ry.
!  On output vm is in a.u.
!
USE kinds, ONLY : DP 
USE io_global, ONLY : stdout
USE eos, ONLY : eos_press_pol, eos_dpress_pol
USE control_ev, ONLY : ieos
IMPLICIT NONE
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, a(m1)
REAL(DP), INTENT(OUT) :: vm
REAL(DP), PARAMETER :: tol=1.D-14
INTEGER, PARAMETER :: maxstep=100
REAL(DP) :: fx, fxp, v1, v1old
INTEGER :: istep

v1=v0
v1old=v1
DO istep=1,maxstep
   CALL eos_press_pol(ieos, v1, fx, v0, b0, b01, b02, a, m1)
   CALL eos_dpress_pol(ieos, v1, fxp, v0, b0, b01, b02, a, m1)
!
!  Newton method
!
   v1 = v1 + fx/fxp
!   WRITE(stdout,'(5x,"Step", i4, " V1=", f20.12, " f= ", f20.12)') istep, v1, -fx
!   FLUSH(stdout)
   IF (ABS(v1-v1old) < tol .OR. ABS(fx) < tol ) GOTO 100
   v1old=v1
ENDDO
CALL errore('find_min_mur_pol','minimum not found',1)
100 CONTINUE
vm=v1
!WRITE(stdout,'("Vmin", 3f20.12)') vm

RETURN
END SUBROUTINE find_min_mur_pol
!
!----------------------------------------------------------------------
SUBROUTINE write_mur_pol(omega0, b0, b01, b02, emin, a_t, m1, itempp)
!----------------------------------------------------------------------
!
!  This routine writes on file the energy versus volume
!  curve, together with the pressure versus volume curve. Depending
!  on which equation of state has been used to interpolate the
!  T=0 K energy it calls the appropriate routine.
!  It receives also the coefficients of the polynomial which interpolates
!  the phonon (+ electron if available) free energy.
!  It receives the parameters of the equation of state:
!
! in input emin in Ry, omega0 in (a.u.)**3, b0 in kbar, b01 adimensional
! b02 in 1/kbar
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,            ONLY : DP
USE data_files,       ONLY : flanhar
USE thermo_mod,       ONLY : ngeo, omega_geo
USE control_mur,      ONLY : p0
USE control_vol,      ONLY : nvol, vmin_input, vmax_input, deltav
USE control_ev,       ONLY : ieos
USE anharmonic_vt,    ONLY : press_vt
USE polyfit_mod,      ONLY : compute_poly, compute_poly_deriv
USE eos,              ONLY : eos_energy_pol, eos_press_pol
USE control_pressure, ONLY : pressure_kb
USE temperature,      ONLY : temp, itemp_plot
USE constants,        ONLY : ry_kbar
USE io_global,        ONLY : meta_ionode

IMPLICIT NONE

REAL(DP), INTENT(IN) :: emin, omega0, b0, b01, b02, a_t(m1)
INTEGER, INTENT(IN) :: itempp
CHARACTER(LEN=256)   :: filename, filename1
CHARACTER(LEN=8) :: float_to_char
REAL(DP) :: omega, free, pres, e, p
INTEGER  :: i, j, m1, iu_mur, itemp
INTEGER  :: find_free_unit

itemp=itemp_plot(itempp)
IF (itemp<1) RETURN

filename="anhar_files/"//TRIM(flanhar)//'.mur_temp'
CALL add_value(filename,temp(itemp))
CALL add_pressure(filename)

IF (meta_ionode) THEN
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   CALL write_mur_start_line(itemp, iu_mur)
   DO i=1,nvol
      omega= vmin_input + deltav * (i-1)
      CALL eos_energy_pol(ieos, omega, e, omega0, b0/ry_kbar, b01, &
                                                  b02*ry_kbar, a_t, m1)
      e=e+emin
      CALL eos_press_pol(ieos, omega, p, omega0, b0/ry_kbar, b01, &
                                                 b02*ry_kbar, a_t, m1)
      press_vt(i,itempp)=p*ry_kbar + pressure_kb
      WRITE(iu_mur,'(f18.10,4f20.10)') omega, e, e+p*omega, &
                                        press_vt(i,itempp), p0(i)
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_mur_pol
!
!----------------------------------------------------------------------
SUBROUTINE write_thermal_press(a_t, m1, itempp)
!----------------------------------------------------------------------
!
!  This routine receives the coefficients of the polynomial which interpolates
!  the phonon (+ electron if available) free energy and
!  writes on file the interpolated free energy and
!  the interpolated thermal pressure versus volume.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE data_files,       ONLY : flanhar
USE control_vol,      ONLY : nvol, vmin_input, deltav
USE polyfit_mod,      ONLY : compute_poly, compute_poly_deriv
USE eos,              ONLY : eos_energy_pol, eos_press_pol
USE temperature,      ONLY : temp, itemp_plot
USE io_global,        ONLY : meta_ionode

IMPLICIT NONE
INTEGER, INTENT(IN) :: itempp, m1
REAL(DP), INTENT(IN) :: a_t(m1)
CHARACTER(LEN=256)   :: filename
REAL(DP) :: omega, free, pres
INTEGER  :: ivol, iu_mur, itemp
INTEGER  :: find_free_unit

itemp=itemp_plot(itempp)
IF (itemp<1) RETURN

filename="anhar_files/"//TRIM(flanhar)//'.poly_free_temp'
CALL add_value(filename,temp(itemp))

IF (meta_ionode) THEN
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_mur,'("#",2x,"Volume (a.u.)^3",6x,"free energy",6x,&
       & "thermal pressure (kbar) ")')
   DO ivol=1,nvol
      omega= vmin_input + deltav * (ivol-1)
      CALL compute_poly(omega, m1-1, a_t, free)
      CALL compute_poly_deriv(omega, m1-1, a_t, pres)
      WRITE(iu_mur,'(f18.10,2e20.10)') omega, free, -pres * ry_kbar
   ENDDO
   CLOSE (UNIT=iu_mur, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_thermal_press
