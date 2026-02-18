!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE manage_piezo_tensor(nwork,ngeom)
!----------------------------------------------------------------------------
!
USE kinds,                ONLY : DP
USE initial_conf,         ONLY : ibrav_save
USE constants,            ONLY : electron_si, bohr_radius_si
USE thermo_mod,           ONLY : energy_geo
USE thermo_sym,           ONLY : code_group_save
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions, epsil_geo, &
                                 el_con_ibrav_geo, el_con_celldm_geo,      &
                                 el_con_at_geo, el_con_celldm_geo,         &
                                 el_con_omega_geo, epsil_geo,              &
                                 start_geometry_qha, last_geometry_qha

USE control_piezoelectric_tensor, ONLY : g_piezo_tensor_geo, &
                                 eg_piezo_tensor_geo, e_piezo_tensor_geo, &
                                 d_piezo_tensor_geo, polar0_geo,          &
                                 decompose_piezo
USE piezoelectric_tensor, ONLY : compute_improper_piezo_tensor, &
                                 compute_d_piezo_tensor,        &
                                 polar_strain, print_piezo_tensor,&
                                 e_piezo_tensor, tot_b_phase,     &
                                 eg_piezo_tensor, g_piezo_tensor, &
                                 d_piezo_tensor,                  &
                                 compute_proper_piezo_tensor,     &
                                 compute_polarization_equil,      &
                                 proper_improper_piezo,           &
                                 clean_piezo_tensor,              &  
                                 write_piezo_tensor
                                    
USE elastic_constants,    ONLY : epsilon_geo, el_con, el_compliances, &
                                 read_elastic
USE data_files,           ONLY : fl_el_cons, fl_piezo

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : my_image_id, root_image, nproc_image
USE mp,                   ONLY : mp_bcast, mp_sum
USE io_global,            ONLY : stdout, meta_ionode_id 

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom
LOGICAL :: exst

INTEGER :: iwork, igeom, base_ind, work_base, work_base_eff, nwork_eff
REAL(DP), ALLOCATABLE :: epsilon_geo_eff(:,:,:), energy_geo_eff(:), &
                         epsil_geo_eff(:), polar_strain_eff(:,:),   &
                         tot_b_phase_eff(:,:) 
CHARACTER(LEN=256) :: filepiezo, filelastic
CHARACTER(LEN=80)  :: label
CHARACTER(LEN=6)   :: int_to_char
LOGICAL :: lreturn
REAL(DP) :: fact
!
!  First collect the total energies
!
CALL mp_sum(energy_geo, world_comm)
energy_geo=energy_geo / nproc_image
!
!  the elastic constants are calculated here if we have the energies
!  of all geometries
!
lreturn=.FALSE.
work_base = nwork / ngeom
DO igeom=start_geometry_qha,last_geometry_qha
   base_ind=(igeom-1)*work_base
   DO iwork=1,work_base
      lreturn=lreturn.OR.(ABS(energy_geo(base_ind+iwork))<1.D-10)
   ENDDO
ENDDO
IF (lreturn) RETURN
!
!  First collect the polarization among all images
!
CALL mp_sum(polar_strain, world_comm)
polar_strain=polar_strain / nproc_image
CALL mp_sum(tot_b_phase, world_comm)
tot_b_phase=tot_b_phase / nproc_image

ALLOCATE(energy_geo_eff(nwork))
ALLOCATE(epsilon_geo_eff(3,3,nwork))
ALLOCATE(epsil_geo_eff(nwork))
ALLOCATE(polar_strain_eff(3,nwork))
ALLOCATE(tot_b_phase_eff(3,nwork))

CALL redefine_energies(energy_geo, epsilon_geo, epsil_geo, nwork,  &
                       energy_geo_eff, epsilon_geo_eff, nwork_eff)

work_base_eff = nwork_eff / ngeom

CALL redefine_polar(polar_strain, tot_b_phase, nwork, polar_strain_eff, &
                    tot_b_phase_eff, epsil_geo, epsil_geo_eff)

DO igeom=start_geometry_qha, last_geometry_qha

   WRITE(stdout,'(2x,76("#"),/)')
   WRITE(stdout,'(5x,"Computing the piezoelectric tensor for the equilibrium &
                                         &geometry=",i4,/)') igeom
   WRITE(stdout,'(5x,i3,6f10.5,/)') el_con_ibrav_geo(igeom), &
                                   el_con_celldm_geo(:,igeom)
   WRITE(stdout,'(2x,76("#"),/)')
   base_ind= (igeom-1)*work_base_eff
   !
   !  First compute the polarization of the unstrained state
   !
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Polarization of the equilibrium geometry")')
   WRITE(stdout,'(2x,76("-"),/)')
   CALL compute_polarization_equil(polar_strain_eff(:,base_ind+1), &
          epsil_geo_eff(base_ind+1), polar0_geo(:,igeom), work_base_eff, &
          ngeo_strain)
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the improper piezoelectric tensor")')
   WRITE(stdout,'(2x,76("-"),/)')
!
!  the piezoelectric tensor is calculated here. First the improper one
!
   CALL compute_improper_piezo_tensor(polar_strain_eff(:,base_ind+1), &
                epsilon_geo_eff(:,:,base_ind+1), work_base_eff, ngeo_strain, &
                ibrav_save, code_group_save)
   g_piezo_tensor_geo(:,:,igeom)=g_piezo_tensor(:,:)

   label="Improper total piezoelectric tensor gamma_ij [ 10^{-2} e/(a.u.)^2 ]"
   fact=100.0_DP
   CALL print_piezo_tensor(g_piezo_tensor, fact, label, frozen_ions)
   label="Improper total piezoelectric tensor gamma_ij [ C/m^2 ]"
   fact= electron_si / (bohr_radius_si)**2
   CALL print_piezo_tensor(g_piezo_tensor, fact, label, frozen_ions)

   CALL proper_improper_piezo(polar0_geo(:,igeom), g_piezo_tensor, &
                                                      eg_piezo_tensor, -1)
   CALL clean_piezo_tensor(eg_piezo_tensor, ibrav_save, code_group_save)

   eg_piezo_tensor_geo(:,:,igeom)=eg_piezo_tensor(:,:)

   label="Proper total transformed piezoelectric tensor gamma_ij [C/m^2] "
   fact= electron_si / (bohr_radius_si)**2
   CALL print_piezo_tensor( eg_piezo_tensor, fact, label, frozen_ions)

   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the proper piezoelectric tensor")')
   WRITE(stdout,'(2x,76("-"),/)')
!
!  and then the proper one (this can be compared with experiments)
!
   CALL compute_proper_piezo_tensor(tot_b_phase_eff(:,base_ind+1),  &
           epsilon_geo_eff(:,:,base_ind+1), work_base_eff, ngeo_strain, &
           ibrav_save, code_group_save, el_con_at_geo(:,:,igeom) )

   e_piezo_tensor=e_piezo_tensor * el_con_celldm_geo(1,igeom) / &
                                           el_con_omega_geo(igeom)
   e_piezo_tensor_geo(:,:,igeom)=e_piezo_tensor(:,:)
   label="Proper total piezoelectric tensor gamma_ij [ 10^{-2} e/(a.u.)^2 ]"
   fact=100.0_DP
   CALL print_piezo_tensor(e_piezo_tensor, fact, label, frozen_ions)

   label="Proper total piezoelectric tensor gamma_ij [ C/m^2 ]"
   fact= electron_si / (bohr_radius_si)**2
   CALL print_piezo_tensor(e_piezo_tensor, fact, label, frozen_ions)
!
!  If a file with the elastic constants exists read them and compute the
!  strain piezoelectric tensor.
!
   filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.g'//&
                                                     TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) CALL read_elastic(filelastic, exst)
   CALL mp_bcast(exst, meta_ionode_id, world_comm)
   IF (exst) THEN
      CALL mp_bcast(el_con, meta_ionode_id, world_comm)
      CALL mp_bcast(el_compliances, meta_ionode_id, world_comm)
      CALL compute_d_piezo_tensor(el_compliances)
      d_piezo_tensor_geo(:,:,igeom)=d_piezo_tensor(:,:)
      label=" Strain total piezoelectric tensor d_ij [pC/N]"
!
!  the factor 10000.0 comes from the fact that the compliances were in 1/kbar
!  that is in 1/10^8 Pa, to have pC/N we need to multiply by 10^12 and
!  10^12/10^8=10^4=10000.0
!
      fact= electron_si / (bohr_radius_si)**2 * 1.D4
      CALL print_piezo_tensor(d_piezo_tensor, fact, label, frozen_ions)
   ENDIF
!
!   Now write all piezo tensors on file
!
   filepiezo='elastic_constants/'//TRIM(fl_piezo)
   IF (frozen_ions) filepiezo=TRIM(filepiezo)//'.fi'
   filepiezo=TRIM(filepiezo)//'.g'//TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) CALL write_piezo_tensor(filepiezo,&
                                              polar0_geo(:,igeom))

   IF (decompose_piezo) CALL write_piezo_decomposition(nwork_eff,igeom)
ENDDO

DEALLOCATE(energy_geo_eff)
DEALLOCATE(epsilon_geo_eff)
DEALLOCATE(epsil_geo_eff)
DEALLOCATE(polar_strain_eff)
DEALLOCATE(tot_b_phase_eff)

RETURN
END SUBROUTINE manage_piezo_tensor

!----------------------------------------------------------------------------
SUBROUTINE write_piezo_decomposition(nwork, igeom)
!----------------------------------------------------------------------------
!
USE kinds,                ONLY : DP
USE constants,            ONLY : electron_si, bohr_radius_si
USE thermo_mod,           ONLY : celldm_geo, tau_geo, uint_geo, ibrav_geo, &
                                 at_geo
USE ions_base,            ONLY : nat
USE initial_conf,         ONLY : ibrav_save
USE thermo_sym,           ONLY : code_group_save
USE control_atomic_pos,   ONLY : max_nint_var
USE control_elastic_constants, ONLY : ngeom, ngeo_strain, nint_var_ec, &
                                      iconstr_internal_ec, epsil_geo,  &
                                      el_con_celldm_geo, el_con_omega_geo
USE control_piezoelectric_tensor, ONLY : e_piezo_tensor_relax_geo, &
                                         e_piezo_tensor_fi_geo,    &
                                         piezo_zeu_geo,            &
                                         dtau_dint_pt, dint_depsilon_geo
USE piezoelectric_tensor, ONLY: e_piezo_tensor, compute_relax_piezo,   &
                                print_piezo_tensor, e_piezo_tensor_fi, &
                                read_piezo_tensor_fi
USE dielectric_constant, ONLY : read_dielectric_properties_from_file
USE data_files,           ONLY : fl_piezo, fl_dielectric
USE mp_images,            ONLY : my_image_id, root_image
USE io_global,            ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, igeom

INTEGER :: work_base, nstep, base_ind, igeo_strain, istep, iwork, na, ipol, &
           ivar, times, i
REAL(DP) :: polar0_geo(3), epsilon_infty(3,3)
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filepiezo
CHARACTER(LEN=80) :: label
REAL(DP) :: fact
LOGICAL :: exist

WRITE(stdout,'(2x,76("*"))')
WRITE(stdout,'(5x,"Decomposition of the piezoelectric tensor")')
WRITE(stdout,'(2x,76("*"),/)')
work_base= nwork / ngeom
nstep= work_base / ngeo_strain
base_ind= (igeom-1)*work_base
!
! first find the internal variables corresponding to each strain
!
DO istep=1,nstep
   DO igeo_strain=1,ngeo_strain
      iwork=base_ind + (istep-1) * ngeo_strain + igeo_strain
      CALL internal_to_tau(ibrav_geo(iwork), celldm_geo(1,iwork),   &
          tau_geo(1,1,iwork), at_geo(1,1,iwork), uint_geo(1,iwork), &
          nat, nint_var_ec(istep), iconstr_internal_ec(istep), 2)
   ENDDO
ENDDO
!
! compute the derivative of the internal variables with respect to strain
!
CALL compute_duint_depsilon(epsil_geo(base_ind+1), uint_geo(1,base_ind+1), &
                            dint_depsilon_geo(1,1,igeom), nint_var_ec,     &
                            nat, max_nint_var, nstep, work_base, ngeo_strain) 
!
! and print them on output
!
WRITE(stdout,'(5x,"Derivatives of the internal variables with &
                   &respect to strain")')
WRITE(label,'(5x,"strain type ")') 
times=MAXVAL(nint_var_ec(1:nstep))
DO i=1, times
   WRITE(label,'(a,"    du(",i1,")/de")') TRIM(label), i
ENDDO
WRITE(stdout,'(a)') TRIM(label)
DO istep=1,nstep
   WRITE(stdout,'(5x,i5,3x,4f15.8)') istep, (dint_depsilon_geo(ivar,istep,&
                                        igeom), ivar=1,nint_var_ec(istep) )
ENDDO
!
! Now read the Born effective charges of this equilibrium geometry from file
!
filepiezo='elastic_constants/'//TRIM(fl_dielectric)//'.g'//&
                                              TRIM(int_to_char(igeom))
IF (my_image_id==root_image) &
   CALL read_dielectric_properties_from_file(filepiezo, nat, epsilon_infty, &
                                         piezo_zeu_geo(1,1,1,igeom))
!
!  and computes the ionic contribution to the piezoelectric tensor
!
CALL compute_relax_piezo(ibrav_save, code_group_save, nat, max_nint_var,   &
                         nint_var_ec, nstep,                               &
                         e_piezo_tensor_relax_geo(1,1,igeom),              &
                         piezo_zeu_geo(1,1,1,igeom),                       &
                         dtau_dint_pt(1,1,1,1,igeom),                         &
                         dint_depsilon_geo(1,1,igeom))
!
!  put the correct units on the ionic contribution
!
e_piezo_tensor_relax_geo(:,:,igeom)= e_piezo_tensor_relax_geo(:,:,igeom) * &
                 el_con_celldm_geo(1,igeom) / el_con_omega_geo(igeom)
!
!  read from file the clamped ion piezoelectric tensor
!
filepiezo='elastic_constants/'//TRIM(fl_piezo)//'.fi'//'.g'//&
                                                TRIM(int_to_char(igeom))
IF (my_image_id==root_image) CALL read_piezo_tensor_fi(filepiezo,&
                                         polar0_geo,exist)

e_piezo_tensor_fi_geo(:,:,igeom)=e_piezo_tensor_fi(:,:)
!
!  Print the two contributions: Clamped ion
!
label= "Proper clamped-ions piezoelectric tensor gamma_ij [ C/m^2 ]"
fact= electron_si / (bohr_radius_si)**2 
CALL print_piezo_tensor(e_piezo_tensor_fi_geo(1,1,igeom), fact, label, .TRUE.)
!
!  Print the two contributions: Ionic 
!
label= "Ionic contribution to the piezoelectric tensor gamma_ij [ C/m^2 ]"
fact= electron_si / (bohr_radius_si)**2 
CALL print_piezo_tensor(e_piezo_tensor_relax_geo(1,1,igeom), fact, &
                                                          label, .FALSE.)
!
!  Print the total proper piezoelectric tensor 
!
e_piezo_tensor(:,:)=e_piezo_tensor_fi_geo(:,:,igeom)+   &
                    e_piezo_tensor_relax_geo(:,:,igeom)

label= "Proper total piezoelectric tensor gamma_ij [ 10^{-2} e/(a.u.)^2 ]"
fact= 100.0_DP
CALL print_piezo_tensor(e_piezo_tensor, fact, label, .FALSE.)

label= "Proper total piezoelectric tensor gamma_ij [ C/m^2 ]"
fact= electron_si / (bohr_radius_si)**2
CALL print_piezo_tensor(e_piezo_tensor, fact, label, .FALSE.)
!
RETURN
END SUBROUTINE write_piezo_decomposition
!
!----------------------------------------------------------------------------
SUBROUTINE compute_duint_depsilon(epsil_geo, uint_geo, duint_depsilon,    &
                                  nint_var_ec, nat, max_nint_var,         &
                                  nstep, work_base, ngeo_strain) 
!----------------------------------------------------------------------------

!
USE kinds, ONLY : DP
USE polyfit_mod, ONLY : polyfit
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat, nstep, work_base, ngeo_strain, max_nint_var 
INTEGER, INTENT(IN) :: nint_var_ec(nstep)
REAL(DP), INTENT(IN) :: epsil_geo(work_base), uint_geo(max_nint_var,work_base) 
REAL(DP), INTENT(OUT) :: duint_depsilon(max_nint_var,nstep)

INTEGER :: ivar, igeo, istep, iwork, base_ind
REAL(DP), ALLOCATABLE :: x(:), f(:)
REAL(DP) :: a(2)

ALLOCATE(x(ngeo_strain))
ALLOCATE(f(ngeo_strain))

iwork=0
base_ind=0
DO istep=1, nstep
   DO ivar=1, nint_var_ec(istep)
      iwork=base_ind
      DO igeo=1,ngeo_strain
         iwork=iwork+1
         x(igeo)=epsil_geo(iwork)
         f(igeo)=uint_geo(ivar,iwork)
      ENDDO
      CALL polyfit(x,f,ngeo_strain,a,1)
      duint_depsilon(ivar,istep)=a(2)
   ENDDO
   base_ind=base_ind+ngeo_strain
ENDDO

DEALLOCATE(x)
DEALLOCATE(f)

RETURN
END SUBROUTINE compute_duint_depsilon
