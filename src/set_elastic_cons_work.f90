!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_elastic_cons_work( nwork )
USE kinds, ONLY : DP
USE cell_base, ONLY : ibrav
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain
USE elastic_constants, ONLY : epsilon_voigt, epsilon_geo, sigma_geo, &
                              trans_epsilon
USE thermo_sym, ONLY : laue
IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsilon_min_off
INTEGER :: igeo, iwork, i, j

epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP
epsilon_min_off= - delta_epsilon * (ngeo_strain - 1 ) 
SELECT CASE (laue) 
   CASE(29,32)
!
!  cubic case
!
      nwork = 2 * ngeo_strain
      ALLOCATE( epsilon_voigt(6, nwork) )
      epsilon_voigt=0.0_DP
      DO igeo=1,ngeo_strain
         epsilon_voigt(3, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 ) 
         epsilon_voigt(4, ngeo_strain + igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
!        epsilon_voigt(5, ngeo_strain + igeo) = &
!                       epsilon_voigt(4, ngeo_strain + igeo)
!        epsilon_voigt(6, ngeo_strain + igeo) = &
!                       epsilon_voigt(4, ngeo_strain + igeo)
      ENDDO
   CASE(19,23,25,27)
!
!  hexagonal or trigonal cases
!
      nwork = 3 * ngeo_strain
      ALLOCATE( epsilon_voigt(6, nwork) )
      epsilon_voigt=0.0_DP
      DO igeo=1,ngeo_strain
         epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 ) 
         epsilon_voigt(3, ngeo_strain+igeo) = epsilon_min + &
                                              delta_epsilon * ( igeo - 1 ) 
         epsilon_voigt(4, 2*ngeo_strain + igeo) = epsilon_min_off + &
                                      2.0_DP * delta_epsilon * ( igeo - 1 )
      ENDDO
   CASE(18,22)
!
!  tetragonal case
!
      nwork = 4 * ngeo_strain
      ALLOCATE( epsilon_voigt(6, nwork) )
      epsilon_voigt=0.0_DP
      DO igeo=1,ngeo_strain
         epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
         epsilon_voigt(3, ngeo_strain+igeo) = epsilon_min + &
                                             delta_epsilon * ( igeo - 1 )
         epsilon_voigt(4, 2*ngeo_strain + igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
         epsilon_voigt(6, 3*ngeo_strain + igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
      ENDDO
CASE DEFAULT
!
!   generic case, used when ibrav=0 or for triclinic, monoclinic, and 
!   orthorombic systems, no information is deduced from symmetry 
!   in the first three cases, all 21 elements of the elastic constants 
!   matrix are computed. 
!   Requires 6 * ngeo_strain self consistent calculations. In the monoclinic
!   and orthorombic case, some elements vanish by symmetry, but we need to
!   make six independent strains in any case.
!
   nwork = 6 * ngeo_strain
   ALLOCATE( epsilon_voigt(6, nwork) )
   epsilon_voigt=0.0_DP
   DO igeo=1,ngeo_strain
      epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
      epsilon_voigt(2, ngeo_strain+igeo) = epsilon_min + &
                                             delta_epsilon * ( igeo - 1 )
      epsilon_voigt(3, 2*ngeo_strain+igeo) = epsilon_min + &
                                             delta_epsilon * ( igeo - 1 )
      epsilon_voigt(4, 3*ngeo_strain + igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
      epsilon_voigt(5, 4*ngeo_strain + igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
      epsilon_voigt(6, 5*ngeo_strain + igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
   ENDDO
END SELECT

ALLOCATE( sigma_geo(3, 3, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
sigma_geo=0.0_DP
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

RETURN
END SUBROUTINE set_elastic_cons_work
