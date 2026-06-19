!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE manage_magnetoelectric_tensor(nwork,ngeom)
!----------------------------------------------------------------------------
!
USE kinds,                ONLY : DP
USE initial_conf,         ONLY : ibrav_save
USE thermo_mod,           ONLY : energy_geo, bfield_geo
USE thermo_sym,           ONLY : a_birss_code
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions, &
                                 start_geometry_qha, last_geometry_qha

USE piezoelectric_tensor, ONLY : polar_strain, tot_b_phase
USE magnetoelectric_tensor, ONLY : compute_magnetoelectric_tensor, &
                                   write_magnetoelectric_tensor
                                    
USE data_files,           ONLY : fl_magnetoelec

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : my_image_id, root_image, nproc_image
USE mp,                   ONLY : mp_bcast, mp_sum
USE io_global,            ONLY : stdout, meta_ionode_id 

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom
LOGICAL :: exst

INTEGER :: iwork, igeom, base_ind, work_base, work_base_eff, nwork_eff, &
           ipol, jpol, istep
REAL(DP) :: magnetoelectric(3,3, nwork)
REAL(DP), ALLOCATABLE :: polar_strain_aux(:,:)
CHARACTER(LEN=256) :: filemagnetoelec
CHARACTER(LEN=6) :: labp(3)
CHARACTER(LEN=1) :: labp1
CHARACTER(LEN=6)   :: int_to_char
LOGICAL :: lreturn
REAL(DP) :: fact
!
!  First collect the total energies
!
CALL mp_sum(energy_geo, world_comm)
energy_geo=energy_geo / nproc_image
!
!  the magnetoelectric tensor is calculated here if we have the energies
!  of all geometries
!
CALL check_for_early_return(energy_geo, nwork, ngeom, lreturn)
IF (lreturn) RETURN
!
!  First collect the polarization and the Berry phase among all images
!
CALL mp_sum(polar_strain, world_comm)
polar_strain=polar_strain / nproc_image
CALL mp_sum(tot_b_phase, world_comm)
tot_b_phase=tot_b_phase / nproc_image

work_base = nwork / ngeom

DO igeom=start_geometry_qha, last_geometry_qha

   WRITE(stdout,'(/,2x,76("#"),/)')
   WRITE(stdout,'(5x,"Computing the magnetoelectric tensor for &
                        &geometry=",i4,/)') igeom
   WRITE(stdout,'(2x,76("#"),/)')
   base_ind= (igeom-1)*work_base
!
!  the magnetoelectric tensor is calculated here.
!
   CALL compute_magnetoelectric_tensor(a_birss_code,            &
           bfield_geo(:,base_ind+1), polar_strain(:,base_ind+1),   &
           work_base, magnetoelectric(1,1,igeom))
!
!
!
   WRITE(stdout,'(/,5x,"Magnetoelectric tensor (ps/m)")')
   WRITE(stdout,'(/,20x,"B_x",17x,"B_y",17x,"B_z")')
   labp(1)=' P_x ('
   labp(2)=' P_y ('
   labp(3)=' P_z ('
   labp1=')'
   DO ipol=1,3
      WRITE(stdout,'(2x,a,3f20.10,1x,a1)') labp(ipol), &
             (magnetoelectric(ipol, jpol, igeom), jpol=1,3), labp1
   ENDDO
!
!   Now write the magnetoelectric tensor on file in the elastic_constants
!   directory
!
   filemagnetoelec='elastic_constants/'//TRIM(fl_magnetoelec)
   IF (frozen_ions) filemagnetoelec=TRIM(filemagnetoelec)//'.fi'
   filemagnetoelec=TRIM(filemagnetoelec)//'.g'//TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) &
        CALL write_magnetoelectric_tensor(filemagnetoelec,magnetoelectric)

ENDDO

RETURN
END SUBROUTINE manage_magnetoelectric_tensor

