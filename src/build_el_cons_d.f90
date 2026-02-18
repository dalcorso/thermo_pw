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
!---------------------------------------------------------------------
SUBROUTINE build_epsilon_zerom1()
  !-----------------------------------------------------------------------
  !
  !  This routine computes the inverse of the dielectric constant 
  !  if epsilon_zero is available
  !
  !
  USE kinds,             ONLY : DP
  USE thermo_mod,        ONLY : epsilon_zero_geo, epsilon_zerom1_geo
  USE io_global,         ONLY : stdout

  USE control_epsilon_infty,  ONLY : epsilon_infty_geo_available, &
                                     lepsilon_infty_geo
  USE control_dielectric_constant, ONLY : epsilonm1_geo_available
  USE matrix_inversion, ONLY : invmat

  USE thermo_mod,       ONLY : tot_ngeo
  !
  IMPLICIT NONE
  INTEGER            :: igeom, ipol
  !
  epsilon_infty_geo_available=.TRUE.
  DO igeom=1,tot_ngeo
     epsilon_infty_geo_available=epsilon_infty_geo_available.AND. &
                                        lepsilon_infty_geo(igeom)
  ENDDO

  IF (.NOT.(epsilon_infty_geo_available)) RETURN
  
  DO igeom=1,tot_ngeo
     CALL invmat(3, epsilon_zero_geo(:,:,igeom), epsilon_zerom1_geo(:,:,igeom))
     WRITE(6,*) 'epsilon0m1 geometry', igeom
     DO ipol=1,3
        WRITE(6,*) epsilon_zerom1_geo(ipol,:,igeom)
     ENDDO
     FLUSH(6)
  ENDDO
  epsilonm1_geo_available=.TRUE.

  RETURN
END SUBROUTINE build_epsilon_zerom1
