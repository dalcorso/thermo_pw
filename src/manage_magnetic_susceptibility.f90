!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE manage_magnetic_susceptibility(nwork,ngeom)
!----------------------------------------------------------------------------
!
USE kinds,                ONLY : DP
USE thermo_mod,           ONLY : energy_geo, bfield_geo
USE control_thermo,       ONLY : lmagnetoelec

USE control_piezomagnetic_tensor, ONLY : mag_strain
USE control_elastic_constants, ONLY : frozen_ions
USE control_bfield,       ONLY : ngeo_b                                    
USE thermo_sym,           ONLY : a_birss_code
USE initial_conf,         ONLY : ibrav_save
USE equilibrium_conf,     ONLY : omega0
USE magnetization_vector, ONLY : compute_magnetic_susceptibility, &
                       compute_magnetic_susceptibility_from_magnetoelectric, &
                                 write_magnetic_susceptibility
USE data_files,           ONLY : fl_magsusc

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : nproc_image, my_image_id, root_image
USE mp,                   ONLY : mp_sum
USE io_global,            ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, ngeom

CHARACTER(LEN=256) :: filemagsusc
CHARACTER(LEN=6) :: int_to_char, labm(3)
CHARACTER(LEN=1) :: labm1
INTEGER :: iwork, ipol, jpol, istart, igeom, iwork_tot
LOGICAL :: lreturn
REAL(DP) :: magnetic_susceptibility(3,3,ngeom)

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(/,5x,"Computing magnetic susceptibility: chi_{ij} &
                                               &= mu0 d M_i / d B_j")')
WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(/,16x,"i j",7x,"B_j (Ry)",8x,"M_i (Bohr mag.)")') 
!
!  First collect the total energies (when not already done)
!
IF (.NOT.lmagnetoelec) THEN
   CALL mp_sum(energy_geo, world_comm)
   energy_geo=energy_geo / nproc_image
!
!  the piezomagnetic tensor is calculated here if we have the energies
!  of all geometries
!
   CALL check_for_early_return(energy_geo, nwork, ngeom, lreturn)
   IF (lreturn) RETURN
!
!  Collect the magnetization among all images
!
   CALL mp_sum(mag_strain, world_comm)
   mag_strain = mag_strain / nproc_image
ENDIF
!
iwork_tot=0
DO igeom=1, ngeom
   istart= 3 * ngeo_b * (igeom-1)+1
   IF (lmagnetoelec) THEN 
      CALL compute_magnetic_susceptibility_from_magnetoelectric( &
                   mag_strain(1,istart), bfield_geo(1,istart),   &
                   nwork/ngeom, ibrav_save, a_birss_code, omega0, &
                   magnetic_susceptibility(1,1,igeom))
   ELSE
      CALL compute_magnetic_susceptibility(mag_strain(1,istart),  &
                   bfield_geo(1,istart), nwork/ngeom, ibrav_save, &
                   omega0, magnetic_susceptibility(1,1,igeom))
   ENDIF
   WRITE(stdout,'(/,5x,"Magnetic susceptibility of geometry ",i5,&
                       &" (SI adimensional)")') igeom
   WRITE(stdout,'(/,20x,"B_x",17x,"B_y",17x,"B_z")') 
   labm(1)=' M_x ('
   labm(2)=' M_y ('
   labm(3)=' M_z ('
   labm1=')'
   DO ipol=1,3
      WRITE(stdout,'(2x,a,3f20.10,1x,a1)') labm(ipol), &
             (magnetic_susceptibility(ipol, jpol, igeom), jpol=1,3), labm1
   ENDDO
   filemagsusc='elastic_constants/'//TRIM(fl_magsusc)
   IF (frozen_ions) filemagsusc=TRIM(filemagsusc)//'.fi'
   filemagsusc=TRIM(filemagsusc)//'.g'//TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) &
      CALL write_magnetic_susceptibility(filemagsusc,&
                             magnetic_susceptibility(1,1,igeom))
ENDDO

!
RETURN
END SUBROUTINE manage_magnetic_susceptibility

