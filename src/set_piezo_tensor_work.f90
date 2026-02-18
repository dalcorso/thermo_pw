!
! Copyright (C) 2014-2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE set_piezo_tensor_work( ngeom, nwork )
!------------------------------------------------------------------------
USE kinds,               ONLY : DP
USE cell_base,           ONLY : ibrav
USE ions_base,           ONLY : nat, ityp, amass
USE thermo_mod,          ONLY : ibrav_geo, celldm_geo, tau_geo, at_geo, &
                                uint_geo, iwho, epsilon_infty_geo,      &
                                zeu_geo, freq_geo, z_geo,               &
                                epsilon_zero_geo, epsilon_zerom1_geo,   &
                                omega_geo
USE control_epsilon_infty, ONLY : lzeu_geo, lepsilon_infty_geo
USE control_thermo,      ONLY : ltau_from_file
USE initial_conf,        ONLY : ibrav_save
USE equilibrium_conf,    ONLY : celldm0, tau0_crys
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain, epsil_geo,&
                                      work_base, el_con_celldm_geo,        &
                                      el_con_ibrav_geo, rot_mat,           &
                                      iconstr_internal_ec, nint_var_ec,    &
                                      max_nint_var, nmove, atom_step,      &
                                      atom_dir, move_at, stype,            &
                                      nstep_ec, min_y, lcm_ec, epsil_y, old_ec,&
                                      stypec, int_ngeo_ec, int_step_ngeo_ec,&
                                      ninternal_ec, elastic_algorithm,     &
                                      tau_acc, el_con_tau_crys_geo
USE control_piezoelectric_tensor, ONLY : e_piezo_tensor_relax_geo,     &
                                         e_piezo_tensor_fi_geo,        &
                                         dtau_dint_pt, dint_depsilon_geo, &
                                         decompose_piezo, piezo_zeu_geo
USE elastic_constants,   ONLY : epsilon_voigt, epsilon_geo, sigma_geo
USE strain_mod,          ONLY : trans_epsilon, set_strain_adv
USE piezoelectric_tensor, ONLY : allocate_piezo
USE rap_point_group,     ONLY : code_group
IMPLICIT NONE
INTEGER, INTENT(IN)  :: ngeom
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsil, tot_mass, sumdisp
INTEGER :: igeo, iwork, istep, nstep, nstep_tot, igeom, istart, na, ivar, &
           iwork_base, ipol, imove
CHARACTER(LEN=2) :: strain_list(21)
LOGICAL :: check_group_ibrav, flag
REAL(DP) :: compute_omega_geo

epsilon_min= - delta_epsilon * (ngeo_strain - 1) / 2.0_DP
IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group)
     CASE(2,16,18,19,20,22,23,25,27,29,32)
!
!  Point groups with inversion, the piezoelectric tensor vanishes
!  C_i, C_2h, C_4h, C_6h, D_2h, D_4h, D_6h, D_3d, S_6, T_h, O_h
!
        nstep=0
     CASE(6,7)
!
!   C_4, C_6
!
        nstep = 4
        strain_list(1) = 'C '
        strain_list(2) = 'E '
        strain_list(3) = 'I '
        strain_list(4) = 'H '

     CASE(8)
!
!   D_2
!
        nstep = 3
        strain_list(1) = 'I '
        strain_list(2) = 'H '
        strain_list(3) = 'G '

     CASE(9)
!
!   D_3
!
        nstep = 2
        strain_list(1) = 'C '
        strain_list(2) = 'I '

     CASE(10,11,28,30)
!
!  D_4, D_6, T, T_d
!
        nstep = 1
        strain_list(1) = 'I '

     CASE(12)
!
!   C_2v
!
        nstep = 5
        strain_list(1) = 'C '
        strain_list(2) = 'D '
        strain_list(3) = 'E '
        strain_list(4) = 'I '
        strain_list(5) = 'H '

     CASE(13,14,15)
!
!   C_3v, C_4v, C_6v
!
        nstep = 3
        strain_list(1) = 'B '
        strain_list(2) = 'E '
        strain_list(3) = 'H '

     CASE(17,21)
!
!  C_3h, D_3h
!
        nstep=1
        strain_list(1) = 'C '
     CASE(24)
!
!  D_2d
!
        nstep = 2
        strain_list(1) = 'I '
        strain_list(2) = 'H '
     CASE(26)
!
!  S_4
!
        nstep = 4 
        strain_list(1) = 'C '
        strain_list(2) = 'I '
        strain_list(3) = 'H '
        strain_list(4) = 'G '

     CASE DEFAULT
        nstep = 6 
        strain_list(1) = 'C '
        strain_list(2) = 'D '
        strain_list(3) = 'E '
        strain_list(4) = 'I '
        strain_list(5) = 'H '
        strain_list(6) = 'G '
   END SELECT
ELSE
   nstep = 6
   strain_list(1) = 'C '
   strain_list(2) = 'D '
   strain_list(3) = 'E '
   strain_list(4) = 'I '
   strain_list(5) = 'H '
   strain_list(6) = 'G '
ENDIF
!
!  Compute the total number of perturbed geometries per unperturbed
!  one not accounting for ngeo_strain, but accounting for the possibility
!  that in some strain types we do a nmove or ninternal_ec(istep) 
!  internal parameter displacements.
!
nstep_tot=0
DO istep=1,nstep
   IF (stype(istep)) THEN
      nstep_tot=nstep_tot+nmove
      int_ngeo_ec(:,istep)=1
   ELSEIF (stypec(istep)) THEN
      ninternal_ec(istep)=1
      CALL clean_int_ngeo(int_ngeo_ec(1,istep), nint_var_ec(istep))
      DO ivar=1,nint_var_ec(istep)
         ninternal_ec(istep)=ninternal_ec(istep)*int_ngeo_ec(ivar,istep)
      ENDDO
      nstep_tot=nstep_tot+ninternal_ec(istep)
   ELSE
      nstep_tot=nstep_tot+1
      int_ngeo_ec(:,istep)=1
   ENDIF
ENDDO
!
!   Total number of distorted geometries
!
nwork = nstep_tot * ngeo_strain * ngeom
!
work_base= nstep_tot * ngeo_strain

nstep_ec= nstep
CALL allocate_piezo(nwork)
ALLOCATE( epsilon_voigt(6, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
ALLOCATE( epsil_geo(nwork) )
ALLOCATE( e_piezo_tensor_relax_geo(3,6,ngeom) )
ALLOCATE( e_piezo_tensor_fi_geo(3,6,ngeom) )
ALLOCATE( sigma_geo(3, 3, nwork) )
ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )
ALLOCATE( at_geo(3,3,nwork) )
ALLOCATE( tau_geo(3,nat,nwork) )
ALLOCATE( rot_mat(3,3,nwork) )
ALLOCATE( uint_geo(max_nint_var,nwork) )
ALLOCATE( tau_acc(3,nat,nwork) )
ALLOCATE( min_y(max_nint_var,ngeo_strain,21,ngeom) )
ALLOCATE( epsil_y(ngeo_strain,21,ngeom) )
ALLOCATE( epsilon_infty_geo(3,3,nwork) )
ALLOCATE( zeu_geo(3,3,nat,nwork) )
ALLOCATE( lepsilon_infty_geo(nwork) )
ALLOCATE( lzeu_geo(nwork) )
lepsilon_infty_geo=.FALSE.
lzeu_geo=.FALSE.
ALLOCATE( epsilon_zero_geo(3,3,nwork) )
ALLOCATE( epsilon_zerom1_geo(3,3,nwork) ) 
ALLOCATE( freq_geo(3*nat,nwork) )
ALLOCATE( z_geo(3*nat,3*nat,nwork) )
ALLOCATE( omega_geo(nwork) )

IF (ngeom==1) THEN
   IF (.NOT. ALLOCATED(el_con_ibrav_geo)) ALLOCATE(el_con_ibrav_geo(1))
   IF (.NOT. ALLOCATED(el_con_celldm_geo)) ALLOCATE(el_con_celldm_geo(6,1))
   IF (.NOT. ALLOCATED(el_con_tau_crys_geo)) &
                            ALLOCATE(el_con_tau_crys_geo(3,nat,1))
   el_con_ibrav_geo(1)=ibrav_save
   el_con_celldm_geo(:,1)=celldm0(:)
   el_con_tau_crys_geo(:,:,1)=tau0_crys(:,:)
ENDIF

flag=(elastic_algorithm=='standard')
iwork=0
tau_acc=0.0_DP
celldm_geo=0.0_DP
DO igeom=1, ngeom
   DO istep=1,nstep
      DO igeo=1,ngeo_strain
         iwork=iwork+1
         epsil = epsilon_min + delta_epsilon * ( igeo - 1 )
         epsil_geo(iwork) = epsil 
         CALL set_strain_adv(strain_list(istep), el_con_ibrav_geo(igeom),  &
              el_con_celldm_geo(1,igeom), epsil, &
              epsilon_voigt(1,iwork), ibrav_geo(iwork), &
              celldm_geo(1,iwork), rot_mat(1,1,iwork), flag )
!
!   set celldm_geo when we use ibrav_geo=0 needed for tau_acc if used
!
         IF (ibrav_geo(iwork)==0) celldm_geo(:,iwork)=celldm0(:)
         IF (stype(istep)) THEN
            iwork_base=iwork
            tau_acc(:,move_at(istep),iwork)=(1.0_DP-(nmove+1.0_DP)/2.0_DP)* &
                                         atom_step(istep)*atom_dir(:,istep)
            DO imove=2, nmove
               iwork=iwork+1
               epsil_geo(iwork)=epsil_geo(iwork_base)
               epsilon_voigt(:,iwork)=epsilon_voigt(:,iwork_base)
               ibrav_geo(iwork)=ibrav_geo(iwork_base)
               celldm_geo(:,iwork)=celldm_geo(:,iwork_base)
               rot_mat(:,:,iwork)=rot_mat(:,:,iwork_base)
               tau_acc(:,move_at(istep),iwork)=(imove-(nmove+1.0_DP)/2.0_DP)* &
                                         atom_step(istep)*atom_dir(:,istep)
            ENDDO
         ENDIF
         IF (stypec(istep)) THEN
            iwork_base=iwork
            DO imove=2, ninternal_ec(istep)
               iwork=iwork+1
               epsil_geo(iwork)=epsil_geo(iwork_base)
               epsilon_voigt(:,iwork)=epsilon_voigt(:,iwork_base)
               ibrav_geo(iwork)=ibrav_geo(iwork_base)
               celldm_geo(:,iwork)=celldm_geo(:,iwork_base)
               rot_mat(:,:,iwork)=rot_mat(:,:,iwork_base)
            ENDDO
         ENDIF
      ENDDO
!
!   set separately the atomic positions for all ngeo_strain*ninternal_ec(istep)
!   geometries
!
      IF (stypec(istep)) THEN
         istart=iwork-ngeo_strain*ninternal_ec(istep)+1
         CALL set_tau_acc(celldm_geo(1,istart), tau_acc(1,1,istart), &
                  uint_geo(1,istart), ngeo_strain*ninternal_ec(istep), nat, &
                  nint_var_ec(istep), istep)
      ENDIF
   ENDDO
ENDDO
!
!  Correct the displacement so that it does not change the center of mass
!  of the system
!
tot_mass=0.0_DP
DO na=1,nat
   tot_mass=tot_mass+amass(ityp(na))
ENDDO
IF ((tot_mass>0.0_DP).AND.lcm_ec) THEN
   DO iwork=1,nwork
      DO ipol=1,3
         sumdisp=0.0_DP
         DO na=1,nat
            sumdisp=sumdisp+tau_acc(ipol,na,iwork) * amass(ityp(na))
         ENDDO
         DO na=1,nat
            tau_acc(ipol,na,iwork)=tau_acc(ipol,na,iwork)-sumdisp/tot_mass
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
!  tau_acc in units of celldm0(1)
!
tau_acc=tau_acc/celldm0(1)

epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

DO iwork=1,nwork
   IF (ibrav_geo(iwork)==0) THEN
      omega_geo(iwork)=compute_omega_geo(ibrav_save, celldm_geo(:,iwork))
   ELSE
      omega_geo(iwork)=compute_omega_geo(ibrav_geo(iwork), celldm_geo(:,iwork))
   ENDIF
ENDDO


IF (decompose_piezo) THEN
   nstep=nwork/ngeo_strain/ngeom
   ALLOCATE(piezo_zeu_geo(3,3,nat,ngeom))
   ALLOCATE(dtau_dint_pt(3,nat,max_nint_var,nstep,ngeom))
   ALLOCATE(dint_depsilon_geo(max_nint_var,nstep,ngeom))
   IF (ltau_from_file) THEN
      DO iwork=1,nwork
         CALL check_geometry_exist(iwork,2,iwho)
      ENDDO
   ENDIF
   DO igeom=1,ngeom
      DO istep=1, nstep
         CALL dtau_dinternal(el_con_celldm_geo(1,igeom),        &
                dtau_dint_pt(1,1,1,istep,igeom), nat, nint_var_ec(istep),&
                        iconstr_internal_ec(istep))
      ENDDO
   ENDDO
ENDIF

RETURN
END SUBROUTINE set_piezo_tensor_work
