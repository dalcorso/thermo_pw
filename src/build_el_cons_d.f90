! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE build_el_cons_d()
  !-----------------------------------------------------------------------
  !
  !  This routine computes the inverse of the dielectric constant and
  !  if the piezoelectric tensor for all geometries is available
  !  and the elastic constants calculated at constan E
  !  are available, it computes the elastic constant at
  !  constant D.
  !
  !
  USE kinds,             ONLY : DP
  USE thermo_mod,        ONLY : epsilon_zerom1_geo
  USE io_global,         ONLY : stdout

  USE control_piezoelectric_tensor, ONLY : piezo_available,  &
                                        piezo_geo_available, &
                                        e_piezo_tensor_geo

  USE control_elastic_constants, ONLY : el_con_geo, el_con_d_geo, &
                      el_cons_available, el_cons_geo_available, el_con_geo, &
                      frozen_ions
  USE elastic_constants, ONLY : ece_to_ecd, print_elastic_constants

  USE thermo_mod,        ONLY : tot_ngeo
  !
  IMPLICIT NONE
  INTEGER            :: igeom
  !
  IF (.NOT.(piezo_available).AND..NOT.(piezo_geo_available)) RETURN
  IF (.NOT.(el_cons_available).AND..NOT.(el_cons_geo_available)) RETURN
  
  DO igeom=1,tot_ngeo
     CALL ece_to_ecd(e_piezo_tensor_geo(1,1,igeom), &
                     epsilon_zerom1_geo(1,1,igeom), &
                     el_con_geo(1,1,igeom), el_con_d_geo(1,1,igeom)) 
     CALL print_elastic_constants(el_con_d_geo(1,1,igeom), frozen_ions)
  ENDDO

  RETURN
END SUBROUTINE build_el_cons_d
!
