!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE write_pyro_t_qha()
!---------------------------------------------------------------------------
!
!  This routine computes the temperature dependent pyroelectric coefficient
!  at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at temperature T.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE ions_base,  ONLY : nat
USE initial_conf, ONLY : ibrav_save
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_atomic_pos, ONLY : dtau_duint, iconstr_internal, &
                               linternal_thermo, max_nint_var, nint_var
USE control_piezoelectric_tensor, ONLY : lpiezo, lpiezof
USE control_pyroelectric_tensor, ONLY : lpyro, lpyrof
USE polarization_vector, ONLY : write_pyro_on_file, compute_pyro, &
                                piezo_pyro
USE equilibrium_conf, ONLY : celldm0, omega0
USE anharmonic, ONLY : e_piezo_tensor_t, uint_t, pyro_t, piezo_pyro_t, &
                       alpha_anis_t, zeu_t, alpha_int_t, celldm_t, &
                       uint_zsisa_t, alpha_int_zsisa_t
USE ph_freq_anharmonic, ONLY : e_piezo_tensorf_t, uintf_t, pyrof_t,    &
                               piezo_pyrof_t, alphaf_anis_t,           &
                               zeuf_t, alphaf_int_t, celldmf_t,        &
                               uintf_zsisa_t, alphaf_int_zsisa_t
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE voigt,       ONLY : to_voigt2
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filepyro, filename
CHARACTER(LEN=80) :: astring
INTEGER :: igeo, nvar, ndata, jdata
REAL(DP) :: alpha_t(3,3), omega, alpha_aux(6)
REAL(DP) :: compute_omega_geo

INTEGER :: i, j, idata, itemp, ipol, jpol, startt, lastt
INTEGER :: compute_nwork 

ALLOCATE(dtau_duint(3,nat,max_nint_var))

pyro_t=0.0_DP
pyrof_t=0.0_DP
piezo_pyro_t=0.0_DP
piezo_pyrof_t=0.0_DP

IF (iconstr_internal==0) RETURN

IF (ltherm_dos) THEN
  filename='anhar_files/'//TRIM(flanhar)//'.uint'
  CALL add_pressure(filename)
  CALL read_uint_anis(uint_t, alpha_int_t, nint_var, temp, ntemp, filename)
ENDIF

IF (ltherm_freq) THEN
  filename='anhar_files/'//TRIM(flanhar)//'.uint_ph'
  CALL add_pressure(filename)
  CALL read_uint_anis(uintf_t, alphaf_int_t, nint_var, temp, ntemp, filename)
  filename='anhar_files/'//TRIM(flanhar)//'.uint_zsisa_ph'
  CALL add_pressure(filename)
  CALL read_uint_anis(uintf_zsisa_t, alphaf_int_zsisa_t, nint_var, temp, &
                                                          ntemp, filename)
ENDIF

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp=startt,lastt
   IF (itemp==1.OR.itemp==ntemp) CYCLE
   
   IF (lpiezo.AND.ltherm_dos) THEN
      CALL dtau_dinternal(celldm_t(1,itemp), dtau_duint, nat, nint_var, &
                                                         iconstr_internal) 
      alpha_aux(:)=alpha_int_t(:,itemp)*uint_t(:,itemp) - &
                   alpha_int_zsisa_t(:,itemp)*uint_zsisa_t(:,itemp)
      CALL compute_pyro(zeu_t(1,1,1,itemp), dtau_duint, &
          alpha_aux, pyro_t(1,itemp), nat, max_nint_var, nint_var)     
      omega=compute_omega_geo(ibrav_save, celldm_t(1,itemp))
      pyro_t(:,itemp) = pyro_t(:,itemp) * celldmf_t(1,itemp) / omega


      CALL to_voigt2(alpha_anis_t(1,itemp),alpha_t,.FALSE.)


      CALL piezo_pyro(e_piezo_tensor_t(1,1,itemp), alpha_t, &
                     piezo_pyro_t(1,itemp))
   END IF

   IF (lpiezof.AND.ltherm_freq) THEN
      CALL dtau_dinternal(celldmf_t(1,itemp), dtau_duint, nat, nint_var, &
                                                         iconstr_internal) 
      alpha_aux(:)=alphaf_int_t(:,itemp)*uintf_t(:,itemp) -  &
                   alphaf_int_zsisa_t(:,itemp)*uintf_zsisa_t(:,itemp)
      CALL compute_pyro(zeuf_t(1,1,1,itemp), dtau_duint,  &
         alpha_aux, pyrof_t(1,itemp), nat, max_nint_var, nint_var)
      omega=compute_omega_geo(ibrav_save, celldmf_t(1,itemp))
      pyrof_t(:,itemp) = pyrof_t(:,itemp) * celldmf_t(1,itemp) / omega
      CALL to_voigt2(alphaf_anis_t(1,itemp),alpha_t,.FALSE.)
      CALL piezo_pyro(e_piezo_tensorf_t(1,1,itemp), alpha_t, &
                     piezo_pyrof_t(1,itemp))
   END IF
ENDDO

IF (ltherm_dos) THEN
   CALL mp_sum(pyro_t, world_comm)
   CALL mp_sum(piezo_pyro_t, world_comm)
   lpyro=.TRUE.
   filepyro='anhar_files/'//TRIM(flanhar)//'.pyro'
   astring="#    Pyroelectric tensor (e/(a.u.)^2/K)"
   CALL write_pyro_on_file(temp, ntemp, pyro_t, piezo_pyro_t, &
                                       astring, filepyro, 0)
ENDIF
!
IF (ltherm_freq) THEN
   CALL mp_sum(pyrof_t, world_comm)
   CALL mp_sum(piezo_pyrof_t, world_comm)
   lpyrof=.TRUE.
   filepyro='anhar_files/'//TRIM(flanhar)//'.pyro_ph'
   astring="#    Pyroelectric tensor (e/(a.u.)^2/K)"
   CALL write_pyro_on_file(temp, ntemp, pyrof_t, piezo_pyrof_t, &
                                       astring, filepyro, 0)
ENDIF

DEALLOCATE(dtau_duint)
!
RETURN
END SUBROUTINE write_pyro_t_qha
!
!---------------------------------------------------------------------------
SUBROUTINE write_pyro_pt_qha()
!---------------------------------------------------------------------------
!
!  This routine computes the temperature dependent pyroelectric tensor
!  at the crystal parameters that minimize the Helmholtz (or Gibbs) 
!  free energy at temperature T for several pressures.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE initial_conf, ONLY : ibrav_save
USE ions_base,  ONLY : nat
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_pressure, ONLY : ipress_plot, npress_plot, press
USE control_piezoelectric_tensor, ONLY : lpiezo_pt, lpiezof_pt
USE control_pyroelectric_tensor, ONLY : lpyro_pt, lpyrof_pt
USE control_atomic_pos, ONLY : dtau_duint, iconstr_internal, &
                               linternal_thermo, max_nint_var, nint_var
USE polarization_vector, ONLY : write_pyro_on_file, compute_pyro, &
                                piezo_pyro
USE anharmonic_pt, ONLY : e_piezo_tensor_pt, uint_pt, pyro_pt, piezo_pyro_pt, &
                       alpha_anis_pt, zeu_pt, alpha_int_pt, celldm_pt, &
                       uint_zsisa_pt, alpha_int_zsisa_pt
USE ph_freq_anharmonic_pt, ONLY : e_piezo_tensorf_pt, uintf_pt, pyrof_pt, &
                               piezo_pyrof_pt, alphaf_anis_pt,            & 
                               zeuf_pt, alphaf_int_pt, celldmf_pt,        &
                               uintf_zsisa_pt, alphaf_int_zsisa_pt
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE voigt,       ONLY : to_voigt2
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filepyro, filename
CHARACTER(LEN=80) :: astring
INTEGER :: igeo, nvar, ndata, jdata
REAL(DP) :: alpha_t(3,3), omega, alpha_aux(6)
REAL(DP) :: compute_omega_geo
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)

INTEGER :: i, j, idata, itemp, startt, lastt, ipressp, ipress
INTEGER :: compute_nwork 

ALLOCATE(dtau_duint(3,nat,max_nint_var))
pyro_pt=0.0_DP
pyrof_pt=0.0_DP
piezo_pyro_pt=0.0_DP
piezo_pyrof_pt=0.0_DP

IF (iconstr_internal==0) RETURN

CALL divide(world_comm, ntemp, startt, lastt)

DO ipressp=1,npress_plot
   ipress=ipress_plot(ipressp)

   IF (ltherm_dos) THEN
      filename='anhar_files/'//TRIM(flanhar)//'.uint_press'
      CALL add_value(filename, press(ipress))
      CALL read_uint_anis(uint_pt, alpha_int_pt, nint_var, temp, &
                                                ntemp, filename)
      filename='anhar_files/'//TRIM(flanhar)//'.uint_zsisa_ph_press'
      CALL add_value(filename, press(ipress))
      CALL read_uint_anis(uint_zsisa_pt, alpha_int_zsisa_pt, nint_var, &
                                                     temp, ntemp, filename)
   ENDIF

   IF (ltherm_freq) THEN
      filename='anhar_files/'//TRIM(flanhar)//'.uint_ph_press'
      CALL add_value(filename, press(ipress))
      CALL read_uint_anis(uintf_pt(1,1,ipressp), alphaf_int_pt(1,1,ipressp), &
                       nint_var, temp, ntemp, filename)
      filename='anhar_files/'//TRIM(flanhar)//'.uint_zsisa_ph_press'
      CALL add_value(filename, press(ipress))
      CALL read_uint_anis(uintf_zsisa_pt(1,1,ipressp), &
           alphaf_int_zsisa_pt(1,1,ipressp), nint_var, temp, ntemp, filename)
   ENDIF

   DO itemp=startt,lastt
      IF (itemp==1.OR.itemp==ntemp) CYCLE
      IF (lpiezo_pt.AND.ltherm_dos) THEN
         CALL dtau_dinternal(celldm_pt(1,itemp,ipressp), dtau_duint, nat, &
                                nint_var, iconstr_internal)
         alpha_aux(:)=alpha_int_pt(:,itemp,ipressp)*uint_pt(:,itemp,ipressp)&
                 -alpha_int_zsisa_pt(:,itemp,ipressp)* &
                  uint_zsisa_pt(:,itemp,ipressp)
         CALL compute_pyro(zeu_pt(1,1,1,itemp,ipressp), dtau_duint, &
             alpha_aux, pyro_pt(1,itemp,ipressp), nat, max_nint_var, nint_var)
         omega=compute_omega_geo(ibrav_save, celldm_pt(1,itemp,ipressp))
         pyro_pt(:,itemp,ipressp) = pyro_pt(:,itemp,ipressp) * &
                                   celldmf_pt(1,itemp,ipressp) / omega
         CALL to_voigt2(alpha_anis_pt(1,itemp,ipressp),alpha_t,.FALSE.)

         CALL piezo_pyro(e_piezo_tensor_pt(1,1,itemp,ipressp), alpha_t, &
                     piezo_pyro_pt(1,itemp,ipressp))
      END IF

      IF (lpiezof_pt.AND.ltherm_freq) THEN
         CALL dtau_dinternal(celldmf_pt(1,itemp,ipressp), dtau_duint, nat, &
                                                nint_var, iconstr_internal)
         alpha_aux(:)=alphaf_int_pt(:,itemp,ipressp)*uintf_pt(:,itemp,ipressp)&
                 -alphaf_int_zsisa_pt(:,itemp,ipressp)* &
                  uintf_zsisa_pt(:,itemp,ipressp)
         CALL compute_pyro(zeuf_pt(1,1,1,itemp,ipressp), dtau_duint, &
             alpha_aux, pyrof_pt(1,itemp,ipressp), nat, max_nint_var, nint_var)
         omega=compute_omega_geo(ibrav_save, celldmf_pt(1,itemp,ipressp))
         pyrof_pt(:,itemp,ipressp) = pyrof_pt(:,itemp,ipressp) * &
                                   celldmf_pt(1,itemp,ipressp) / omega
         CALL to_voigt2(alphaf_anis_pt(1,itemp,ipressp),alpha_t,.FALSE.)

         CALL piezo_pyro(e_piezo_tensorf_pt(1,1,itemp,ipressp), alpha_t, &
                     piezo_pyrof_pt(1,itemp,ipressp))
      END IF
   ENDDO

   IF (ltherm_dos) THEN
      CALL mp_sum(pyro_pt(:,:,ipressp), world_comm)
      CALL mp_sum(piezo_pyro_pt(:,:,ipressp), world_comm)
      lpyro_pt=.TRUE.
      filepyro='anhar_files/'//TRIM(flanhar)//'.pyro_press'
      CALL add_value(filepyro, press(ipress))
      astring="#    Pyroelectric tensor (e/(a.u.)^2/K)"
      CALL write_pyro_on_file(temp, ntemp, pyro_pt(1,1,ipressp), &
                           piezo_pyro_pt(1,1,ipressp), astring, filepyro, 0)
   ENDIF
!
   IF (ltherm_freq) THEN
      CALL mp_sum(pyrof_pt(:,:,ipressp), world_comm)
      CALL mp_sum(piezo_pyrof_pt(:,:,ipressp), world_comm)
      lpyrof_pt=.TRUE.
      filepyro='anhar_files/'//TRIM(flanhar)//'.pyro_ph_press'
      CALL add_value(filepyro, press(ipress))
      astring="#    Pyroelectric tensor (e/(a.u.)^2/K)"
      CALL write_pyro_on_file(temp, ntemp, pyrof_pt(1,1,ipressp), &
                    piezo_pyrof_pt(1,1,ipressp), astring, filepyro, 0)
   ENDIF
ENDDO
DEALLOCATE(dtau_duint)
!
RETURN
END SUBROUTINE write_pyro_pt_qha
!
!---------------------------------------------------------------------------
SUBROUTINE write_pyro_ptt_qha()
!---------------------------------------------------------------------------
!
!  This routine computes the pressure dependent pyroelectric tensor
!  at the crystal parameters that minimize the Gibbs free energy at 
!  pressure p for several temperatures. 
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE initial_conf, ONLY : ibrav_save
USE ions_base,  ONLY : nat
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_pressure, ONLY : ipress_plot, npress_plot, press
USE control_piezoelectric_tensor, ONLY : lpiezo_ptt, lpiezof_ptt
USE control_pyroelectric_tensor, ONLY : lpyro_ptt, lpyrof_ptt
USE control_atomic_pos, ONLY : dtau_duint, iconstr_internal, &
                               linternal_thermo, max_nint_var, nint_var
USE polarization_vector, ONLY : write_pyro_on_file, compute_pyro, &
                                piezo_pyro
USE anharmonic_ptt, ONLY : e_piezo_tensor_ptt, uint_ptt, pyro_ptt, &
                       piezo_pyro_ptt, &
                       alpha_anis_ptt, zeu_ptt, alpha_int_ptt, celldm_ptt, &
                       uint_zsisa_ptt, alpha_int_zsisa_ptt
USE ph_freq_anharmonic_ptt, ONLY : e_piezo_tensorf_ptt, uintf_ptt, pyrof_ptt, &
                               piezo_pyrof_ptt, alphaf_anis_ptt,            & 
                               zeuf_ptt, alphaf_int_ptt, celldmf_ptt,       &
                               uintf_zsisa_ptt, alphaf_int_zsisa_ptt
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp, ntemp_plot, itemp_plot
USE control_pressure, ONLY : npress
USE voigt,       ONLY : to_voigt2
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filepyro, filename
CHARACTER(LEN=80) :: astring
INTEGER :: igeo, nvar, ndata, jdata
REAL(DP) :: alpha_t(3,3), omega, alpha_aux(6)
REAL(DP) :: compute_omega_geo

INTEGER :: i, j, idata, itemp, startp, lastp, itempp, ipress
INTEGER :: compute_nwork 

ALLOCATE(dtau_duint(3,nat,max_nint_var))
pyro_ptt=0.0_DP
pyrof_ptt=0.0_DP
piezo_pyro_ptt=0.0_DP
piezo_pyrof_ptt=0.0_DP

IF (iconstr_internal==0) RETURN

CALL divide(world_comm, npress, startp, lastp)

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)

   IF (ltherm_dos) THEN
      filename='anhar_files/'//TRIM(flanhar)//'.uint_temp'
      CALL add_value(filename, temp(itemp))
      CALL read_uint_anis(uint_ptt(1,1,itempp), alpha_int_ptt(1,1,itempp), &
                                     nint_var, press, npress, filename)
      filename='anhar_files/'//TRIM(flanhar)//'.uint_zsisa_ph_temp'
      CALL add_value(filename, temp(itemp))
      CALL read_uint_anis(uint_zsisa_ptt(1,1,itempp),                      &
                     alpha_int_zsisa_ptt(1,1,itempp), nint_var,            &
                                                  press, npress, filename)
   ENDIF

   IF (ltherm_freq) THEN
      filename='anhar_files/'//TRIM(flanhar)//'.uint_ph_temp'
      CALL add_value(filename, temp(itemp))
      CALL read_uint_anis(uintf_ptt(1,1,itempp), alphaf_int_ptt(1,1,itempp), &
                       nint_var, press, npress, filename)
      filename='anhar_files/'//TRIM(flanhar)//'.uint_zsisa_ph_temp'
      CALL add_value(filename, temp(itemp))
      CALL read_uint_anis(uintf_zsisa_ptt(1,1,itempp), &
           alphaf_int_zsisa_ptt(1,1,itempp), nint_var, press, npress, filename)
   ENDIF

   DO ipress=startp,lastp
      IF (ipress==npress) CYCLE
      IF (lpiezo_ptt.AND.ltherm_dos) THEN
         CALL dtau_dinternal(celldm_ptt(1,ipress,itempp), dtau_duint, nat, &
                                nint_var, iconstr_internal)
         alpha_aux(:)=alpha_int_ptt(:,ipress,itempp)*uint_ptt(:,ipress,itempp)&
                 -alpha_int_zsisa_ptt(:,ipress,itempp)* &
                  uint_zsisa_ptt(:,ipress,itempp)
         CALL compute_pyro(zeu_ptt(1,1,1,ipress,itempp), dtau_duint, &
             alpha_aux, pyro_ptt(1,ipress,itempp), nat, max_nint_var, nint_var)
         omega=compute_omega_geo(ibrav_save, celldm_ptt(1,ipress,itempp))
         pyro_ptt(:,ipress,itempp) = pyro_ptt(:,ipress,itempp) * &
                                   celldmf_ptt(1,ipress,itempp) / omega
         CALL to_voigt2(alpha_anis_ptt(1,ipress,itempp),alpha_t,.FALSE.)

         CALL piezo_pyro(e_piezo_tensor_ptt(1,1,ipress,itempp), alpha_t, &
                     piezo_pyro_ptt(1,ipress,itempp))
      END IF

      IF (lpiezof_ptt.AND.ltherm_freq) THEN
         CALL dtau_dinternal(celldmf_ptt(1,ipress,itempp), dtau_duint, nat, &
                                                nint_var, iconstr_internal)
         alpha_aux(:)=alphaf_int_ptt(:,ipress,itempp)*   &
                             uintf_ptt(:,ipress,itempp)  &
                 -alphaf_int_zsisa_ptt(:,ipress,itempp)* &
                  uintf_zsisa_ptt(:,ipress,itempp)
         CALL compute_pyro(zeuf_ptt(1,1,1,ipress,itempp), dtau_duint, &
             alpha_aux, pyrof_ptt(1,ipress,itempp), nat, max_nint_var, nint_var)
         omega=compute_omega_geo(ibrav_save, celldmf_ptt(1,ipress,itempp))
         pyrof_ptt(:,ipress,itempp) = pyrof_ptt(:,ipress,itempp) * &
                                   celldmf_ptt(1,ipress,itempp) / omega
         CALL to_voigt2(alphaf_anis_ptt(1,ipress,itempp),alpha_t,.FALSE.)

         CALL piezo_pyro(e_piezo_tensorf_ptt(1,1,ipress,itempp), alpha_t, &
                     piezo_pyrof_ptt(1,ipress,itempp))
      END IF
   ENDDO

   IF (ltherm_dos) THEN
      CALL mp_sum(pyro_ptt(:,:,itempp), world_comm)
      CALL mp_sum(piezo_pyro_ptt(:,:,itempp), world_comm)
      lpyro_ptt=.TRUE.
      filepyro='anhar_files/'//TRIM(flanhar)//'.pyro_temp'
      CALL add_value(filepyro, temp(itemp))
      astring="#    Pyroelectric tensor (e/(a.u.)^2/K)"
      CALL write_pyro_on_file(press, npress, pyro_ptt(1,1,itempp), &
                           piezo_pyro_ptt(1,1,itempp), astring, filepyro, 2)
   ENDIF
!
   IF (ltherm_freq) THEN
      CALL mp_sum(pyrof_ptt(:,:,itempp), world_comm)
      CALL mp_sum(piezo_pyrof_ptt(:,:,itempp), world_comm)
      lpyrof_ptt=.TRUE.
      filepyro='anhar_files/'//TRIM(flanhar)//'.pyro_ph_temp'
      CALL add_value(filepyro, temp(itemp))
      astring="#    Pyroelectric tensor (e/(a.u.)^2/K)"
      CALL write_pyro_on_file(press, npress, pyrof_ptt(1,1,itempp), &
                    piezo_pyrof_ptt(1,1,itempp), astring, filepyro, 2)
   ENDIF
ENDDO
DEALLOCATE(dtau_duint)
!
RETURN
END SUBROUTINE write_pyro_ptt_qha
