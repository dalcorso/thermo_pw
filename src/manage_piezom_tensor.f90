!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE manage_piezom_tensor(nwork,ngeom)
!----------------------------------------------------------------------------
!
USE kinds,                ONLY : DP
USE initial_conf,         ONLY : ibrav_save
USE thermo_mod,           ONLY : energy_geo
USE thermo_sym,           ONLY : b_birss_code, b_birss_ext_code
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions, epsil_geo, &
                                 el_con_ibrav_geo, el_con_celldm_geo,      &
                                 el_con_at_geo, el_con_celldm_geo,         &
                                 el_con_omega_geo, epsil_geo,              &
                                 start_geometry_qha, last_geometry_qha

USE control_piezomagnetic_tensor, ONLY : piezom_tensor_geo, mag_strain,    &
                                 mag0_geo
USE piezomagnetic_tensor, ONLY : compute_piezom_tensor, print_piezom_tensor, &
                                 piezom_tensor,                     &
                                 compute_magnetization_equil
                                    
USE elastic_constants,    ONLY : epsilon_geo
USE data_files,           ONLY : fl_piezom

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : my_image_id, root_image, nproc_image
USE mp,                   ONLY : mp_bcast, mp_sum
USE io_global,            ONLY : stdout, meta_ionode_id 

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom
LOGICAL :: exst

INTEGER :: iwork, igeom, base_ind, work_base, work_base_eff, nwork_eff
REAL(DP), ALLOCATABLE :: epsilon_geo_eff(:,:,:), energy_geo_eff(:),    &
                         epsil_geo_eff(:), mag_strain_eff(:,:)
CHARACTER(LEN=256) :: filepiezom, filelastic
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
!  the piezomagnetic tensor is calculated here if we have the energies
!  of all geometries
!
CALL check_for_early_return(energy_geo, nwork, ngeom, lreturn)
IF (lreturn) RETURN
!
!  First collect the magnetization among all images
!
CALL mp_sum(mag_strain, world_comm)
mag_strain = mag_strain / nproc_image
!
!  Redefine energies and magnetization if we have computed them for 
!  several atomic positions
!
ALLOCATE(energy_geo_eff(nwork))
ALLOCATE(epsilon_geo_eff(3,3,nwork))
ALLOCATE(epsil_geo_eff(nwork))
ALLOCATE(mag_strain_eff(3,nwork))
!
CALL redefine_energies(energy_geo, epsilon_geo, epsil_geo, nwork,  &
                       energy_geo_eff, epsilon_geo_eff, nwork_eff)

work_base_eff = nwork_eff / ngeom

CALL redefine_mag(mag_strain, nwork, mag_strain_eff, epsil_geo, &
                                                        epsil_geo_eff)

DO igeom=start_geometry_qha, last_geometry_qha

   WRITE(stdout,'(2x,76("#"),/)')
   WRITE(stdout,'(5x,"Computing the piezomagnetic tensor for the equilibrium &
                                         &geometry=",i4,/)') igeom
   WRITE(stdout,'(5x,i3,6f10.5,/)') el_con_ibrav_geo(igeom), &
                                   el_con_celldm_geo(:,igeom)
   WRITE(stdout,'(2x,76("#"),/)')
   base_ind= (igeom-1)*work_base_eff
   !
   !  First compute the magnetization of the unstrained state
   !
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Magnetization of the equilibrium geometry")')
   WRITE(stdout,'(2x,76("-"),/)')
   CALL compute_magnetization_equil(mag_strain_eff(:,base_ind+1), &
          epsil_geo_eff(base_ind+1), mag0_geo(:,igeom), work_base_eff, &
          ngeo_strain)
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Computing the piezomagnetic tensor")')
   WRITE(stdout,'(2x,76("-"),/)')
!
!  the piezomagnetic tensor is calculated here.
!
   CALL compute_piezom_tensor(mag_strain_eff(:,base_ind+1), &
                epsilon_geo_eff(:,:,base_ind+1), work_base_eff, &
                ngeo_strain, ibrav_save, b_birss_code, b_birss_ext_code)

   piezom_tensor_geo(:,:,igeom)=piezom_tensor(:,:)

   label="Piezomagnetic tensor h_ij [ Bohr magneton ]"
   fact=1.0_DP
   CALL print_piezom_tensor(piezom_tensor, fact, label)

ENDDO

DEALLOCATE(energy_geo_eff)
DEALLOCATE(epsilon_geo_eff)
DEALLOCATE(epsil_geo_eff)
DEALLOCATE(mag_strain_eff)

RETURN
END SUBROUTINE manage_piezom_tensor

SUBROUTINE check_for_early_return(energy_geo, nwork, ngeom, lreturn)

USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : &
                              start_geometry_qha, last_geometry_qha

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom
REAL(DP), INTENT(IN) :: energy_geo(nwork)
LOGICAL, INTENT(OUT) :: lreturn

INTEGER :: work_base, base_ind, igeom, iwork

lreturn=.FALSE.
work_base = nwork / ngeom
DO igeom=start_geometry_qha,last_geometry_qha
   base_ind=(igeom-1)*work_base
   DO iwork=1,work_base
      lreturn=lreturn.OR.(ABS(energy_geo(base_ind+iwork))<1.D-10)
   ENDDO
ENDDO
RETURN
END SUBROUTINE check_for_early_return
