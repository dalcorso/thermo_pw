!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE manage_polarization(nwork)
!----------------------------------------------------------------------------
!
USE thermo_mod,           ONLY : celldm_geo, energy_geo, start_geometry, &
                                 ibrav_geo, last_geometry
USE piezoelectric_tensor, ONLY : polar_strain,       &
                                 tot_b_phase
USE polarization_vector,  ONLY : write_polarization
USE data_files,           ONLY : fl_polar

USE mp_world,             ONLY : world_comm
USE mp_images,            ONLY : my_image_id, root_image, nproc_image
USE mp,                   ONLY : mp_bcast, mp_sum
USE io_global,            ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork
LOGICAL :: exst

INTEGER :: iwork, igeom
CHARACTER(LEN=256) :: filepolar
CHARACTER(LEN=6)   :: int_to_char
LOGICAL :: lreturn
!
!  First collect the total energies
!
CALL mp_sum(energy_geo, world_comm)
energy_geo=energy_geo / nproc_image
!
!  the polarization is printed here if we have the energies
!  of all geometries required in the run
!
lreturn=.FALSE.
DO iwork=start_geometry,last_geometry
   lreturn=lreturn.OR.(ABS(energy_geo(iwork))<1.D-10)
ENDDO
IF (lreturn) RETURN
!
!  First collect the polarization among all images
!
CALL mp_sum(polar_strain, world_comm)
polar_strain=polar_strain / nproc_image
CALL mp_sum(tot_b_phase, world_comm)
tot_b_phase=tot_b_phase / nproc_image

DO igeom=start_geometry, last_geometry
   WRITE(stdout,'(/,2x,76("-"))')
   WRITE(stdout,'(5x,"Polarization of geometry",i5)') igeom
   WRITE(stdout,'(2x,76("-"),/)')
   WRITE(stdout,'(5x,i5,": ", i3,6f10.5)') igeom, ibrav_geo(igeom), &
                                   celldm_geo(:,igeom)
   WRITE(stdout,'(/,5x,"Total Berry phases computed along the three &
                            &reciprocal lattice vectors")') 
   WRITE(stdout,'(3f18.10)') tot_b_phase(1:3,igeom)
   CALL print_polarization(polar_strain(1,igeom), .FALSE. )
   filepolar='elastic_constants/'//TRIM(fl_polar)//'.g'//&
                                                     TRIM(int_to_char(igeom))
   IF (my_image_id==root_image) CALL write_polarization(filepolar,&
                                    polar_strain(1,igeom),tot_b_phase(1,igeom))
ENDDO

RETURN
END SUBROUTINE manage_polarization
