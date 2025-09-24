!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE set_piezo_tensor_work( nwork )
!------------------------------------------------------------------------
USE kinds,               ONLY : DP
USE cell_base,           ONLY : ibrav
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain
USE elastic_constants,   ONLY : epsilon_voigt, epsilon_geo
USE strain_mod,          ONLY : trans_epsilon
USE piezoelectric_tensor, ONLY : polar_geo
USE rap_point_group,     ONLY : code_group
IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsilon_min_off
INTEGER :: igeo, iwork
LOGICAL :: check_group_ibrav

epsilon_min= - delta_epsilon * (ngeo_strain - 1) / 2.0_DP
epsilon_min_off= - delta_epsilon * (ngeo_strain - 1) 
IF (check_group_ibrav(code_group, ibrav)) THEN
   SELECT CASE (code_group)
     CASE(2,16,18,19,20,22,23,25,27,29,32)
!
!  Point groups with inversion, the piezoelectric tensor vanishes
!  C_i, C_2h, C_4h, C_6h, D_2h, D_4h, D_6h, D_3d, S_6, T_h, O_h
!
     CASE(6,7)
!
!   C_4, C_6
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
           epsilon_voigt(5, 3*ngeo_strain + igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(8)
!
!   D_2
!
        nwork = 3 * ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(4, igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
           epsilon_voigt(5, ngeo_strain + igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
           epsilon_voigt(6, 2*ngeo_strain + igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(9)
!
!   D_3
!
        nwork = 2 * ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(1, igeo) = epsilon_min +  delta_epsilon * ( igeo - 1 )
           epsilon_voigt(4, ngeo_strain + igeo) = epsilon_min_off + &
                                          2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(10,11,28,30)
!
!  D_4, D_6, T, T_d
!
        nwork = ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(4, igeo) = epsilon_min_off + &
                                          2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(12)
!
!   C_2v
!
        nwork = 5 * ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
           epsilon_voigt(2, ngeo_strain + igeo) = epsilon_min + &
                                                  delta_epsilon * ( igeo - 1 )
           epsilon_voigt(3, 2*ngeo_strain + igeo) = epsilon_min + &
                                                  delta_epsilon * ( igeo - 1 )
           epsilon_voigt(4, 3*ngeo_strain + igeo) = epsilon_min_off + &
                                         2.0_DP * delta_epsilon * ( igeo - 1 )
           epsilon_voigt(5, 4*ngeo_strain + igeo) = epsilon_min_off + &
                                         2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(13,14,15)
!
!   C_3v, C_4v, C_6v
!
        nwork = 3 * ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
           epsilon_voigt(3, ngeo_strain + igeo) = epsilon_min + &
                                                  delta_epsilon * ( igeo - 1 )
           epsilon_voigt(5, 2*ngeo_strain + igeo) = epsilon_min_off + &
                                         2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(17,21)
!
!  C_3h, D_3h
!
        nwork = ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(24)
!
!  D_2d
!
        nwork = 2 * ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(4, igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
           epsilon_voigt(6, ngeo_strain + igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE(26)
!
!  S_4
!
        nwork = 4 * ngeo_strain
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        DO igeo=1,ngeo_strain
           epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
           epsilon_voigt(4, ngeo_strain + igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
           epsilon_voigt(5, 2 * ngeo_strain + igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
           epsilon_voigt(6, 3 * ngeo_strain + igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
        ENDDO
     CASE DEFAULT
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
ELSE
   nwork = 6 * ngeo_strain
   ALLOCATE( epsilon_voigt(6, nwork) )
   epsilon_voigt=0.0_DP
   DO igeo=1,ngeo_strain
      epsilon_voigt(1, igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
      epsilon_voigt(2, ngeo_strain+igeo) = epsilon_min + &
                                             delta_epsilon * ( igeo - 1 )
      epsilon_voigt(3, 2*ngeo_strain + igeo) = epsilon_min + &
                                             delta_epsilon * ( igeo - 1 )
      epsilon_voigt(4, 3*ngeo_strain + igeo) = epsilon_min_off + &
                                    2.0_DP * delta_epsilon * ( igeo - 1 )
      epsilon_voigt(5, 4*ngeo_strain + igeo) = epsilon_min_off + &
                                    2.0_DP * delta_epsilon * ( igeo - 1 )
      epsilon_voigt(6, 5*ngeo_strain + igeo) = epsilon_min_off + &
                                    2.0_DP * delta_epsilon * ( igeo - 1 )
   ENDDO
ENDIF

ALLOCATE( polar_geo(3, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
polar_geo=0.0_DP
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

RETURN
END SUBROUTINE set_piezo_tensor_work
