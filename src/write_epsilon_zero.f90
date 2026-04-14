!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE write_epsilon_zero(igeom)
!-------------------------------------------------------------------------
USE kinds,            ONLY : DP
USE constants,        ONLY : amu_ry, ry_to_cmm1
USE thermo_mod,       ONLY : z_geo, freq_geo, omega_geo, epsilon_zero_geo, &
                             epsilon_infty_geo, zeu_geo, epsilon_zerom1_geo
USE ions_base,        ONLY : nat, ityp, amass
USE ifc,              ONLY : has_zstar, epsil_ifc, zeu
USE control_epsilon_infty, ONLY : lepsilon_infty_geo, lzeu_geo, &
                                  lepsilon_zero_geo
USE dielectric_constant, ONLY : write_dielectric_properties_to_file, &
                                polar_mode_permittivity_tpw
USE matrix_inversion, ONLY : invmat
USE data_files,       ONLY : fl_dielectric
USE mp_images,        ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6) :: int_to_char
COMPLEX(DP), ALLOCATABLE :: z(:,:)
REAL(DP), ALLOCATABLE :: w2(:)
INTEGER :: imode, na, nta, ipol

IF (.NOT. has_zstar) RETURN

epsilon_infty_geo(:,:,igeom)=epsil_ifc(:,:)
zeu_geo(:,:,:,igeom)=zeu(:,:,:)
lepsilon_infty_geo(igeom)=.TRUE.
lzeu_geo(igeom)=.TRUE.
!
!   Compute epsilon_zero. The routine requires the displacement
!
ALLOCATE(w2(3*nat))
ALLOCATE(z(3*nat,3*nat))
DO imode = 1, 3*nat
   DO na = 1,nat
      nta = ityp(na)
      DO ipol = 1,3
         z((na-1)*3+ipol,imode) = z_geo((na-1)*3+ipol,imode,igeom)/ &
                               SQRT(amu_ry*amass(nta))
      ENDDO
   ENDDO
   w2(imode)=SIGN((freq_geo(imode,igeom) /  ry_to_cmm1)**2, &
                                               freq_geo(imode,igeom) )
ENDDO
!
!  actual calculation
!
CALL polar_mode_permittivity_tpw( nat, epsilon_infty_geo(1,1,igeom), &
            z, zeu_geo(1,1,1,igeom), w2, omega_geo(igeom), &
            epsilon_zero_geo(1,1,igeom))
!
!  compute the inverse
!
CALL invmat(3, epsilon_zero_geo(:,:,igeom), epsilon_zerom1_geo(:,:,igeom))

DEALLOCATE(z)
DEALLOCATE(w2)

filename='elastic_constants/'//TRIM(fl_dielectric)//'.g'//&
                                              TRIM(int_to_char(igeom))
IF (my_image_id==root_image) &
   CALL write_dielectric_properties_to_file(filename, nat,               &
             epsilon_infty_geo(1,1,igeom), epsilon_zero_geo(1,1,igeom),  &
             zeu_geo(1,1,1,igeom) )


lepsilon_zero_geo(igeom)=.TRUE.
RETURN
END SUBROUTINE write_epsilon_zero
