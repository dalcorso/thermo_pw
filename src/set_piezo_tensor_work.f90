!
! Copyright (C) 2014-2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE set_piezo_tensor_work( ngeom, nwork )
!------------------------------------------------------------------------
USE kinds,               ONLY : DP
USE cell_base,           ONLY : ibrav
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain, epsil_geo,&
                                      work_base
USE elastic_constants,   ONLY : epsilon_voigt, epsilon_geo
USE strain_mod,          ONLY : trans_epsilon
USE piezoelectric_tensor, ONLY : allocate_piezo
USE rap_point_group,     ONLY : code_group
IMPLICIT NONE
INTEGER, INTENT(IN)  :: ngeom
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsilon_min_off, epsil
INTEGER :: igeo, iwork, istep, nstep, igeom, istart
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
        nstep = 4
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1, ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(1, istart+igeo) = epsilon_min + delta_epsilon &
                                                             * ( igeo - 1 )
              epsilon_voigt(3, ngeo_strain+istart+igeo) = epsilon_min + &
                                                 delta_epsilon * ( igeo - 1 )
              epsilon_voigt(4, 2*ngeo_strain+istart+igeo) = epsilon_min_off &
                                  + 2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(5, 3*ngeo_strain+istart+igeo) = epsilon_min_off& 
                             + 2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(8)
!
!   D_2
!
        nstep = 3
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1, ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(4, istart+igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(5, ngeo_strain+istart+igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(6, 2*ngeo_strain+istart+igeo) = epsilon_min_off &
                                      + 2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(9)
!
!   D_3
!
        nstep = 2
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1, ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(1, istart+igeo) = epsilon_min + delta_epsilon*&
                                                               ( igeo - 1 )
              epsilon_voigt(4, ngeo_strain+istart+igeo) = epsilon_min_off + &
                                          2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(10,11,28,30)
!
!  D_4, D_6, T, T_d
!
        nstep = 1
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1,ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(4, istart+igeo) = epsilon_min_off + &
                                          2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(12)
!
!   C_2v
!
        nstep = 5
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1,ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(1, istart+igeo) = epsilon_min + delta_epsilon * &
                                                               ( igeo - 1 )
              epsilon_voigt(2, ngeo_strain+istart+igeo) = epsilon_min + &
                                              delta_epsilon * ( igeo - 1 )
              epsilon_voigt(3, 2*ngeo_strain+istart+igeo) = epsilon_min + &
                                              delta_epsilon * ( igeo - 1 )
              epsilon_voigt(4, 3*ngeo_strain+istart+igeo) = epsilon_min_off + &
                                     2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(5, 4*ngeo_strain+istart+igeo) = epsilon_min_off + &
                                     2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(13,14,15)
!
!   C_3v, C_4v, C_6v
!
        nstep = 3
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1,ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(1, istart+igeo) = epsilon_min + delta_epsilon * &
                                                            ( igeo - 1 )
              epsilon_voigt(2, istart+igeo) = epsilon_voigt(1, igeo)
              epsilon_voigt(3, ngeo_strain+istart+igeo) = epsilon_min + &
                                                  delta_epsilon * ( igeo - 1 )
              epsilon_voigt(5, 2*ngeo_strain+istart+igeo) = epsilon_min_off + &
                                         2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(17,21)
!
!  C_3h, D_3h
!
        nstep=1
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1,ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(1, istart+igeo) = epsilon_min + delta_epsilon*&
                                                                   (igeo-1)
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(24)
!
!  D_2d
!
        nstep = 2
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1,ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(4, istart+igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(6, ngeo_strain+istart+igeo) = epsilon_min_off + &
                                        2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE(26)
!
!  S_4
!
        nstep = 4 
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1,ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(1, istart+igeo) = epsilon_min + delta_epsilon * &
                                                                ( igeo - 1 )
              epsilon_voigt(4, ngeo_strain+istart+igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(5, 2*ngeo_strain+istart+igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(6, 3*ngeo_strain+istart+igeo) = epsilon_min_off + &
                                       2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
     CASE DEFAULT
        nstep = 6 
        nwork = nstep * ngeo_strain * ngeom
        ALLOCATE( epsilon_voigt(6, nwork) )
        epsilon_voigt=0.0_DP
        istart=0
        DO igeom=1, ngeom
           DO igeo=1,ngeo_strain
              epsilon_voigt(1,istart+igeo) = epsilon_min + delta_epsilon * &
                                                                  ( igeo - 1 )
              epsilon_voigt(2,ngeo_strain+istart+igeo) = epsilon_min + &
                                                  delta_epsilon * ( igeo - 1 )
              epsilon_voigt(3,2*ngeo_strain+istart+igeo) = epsilon_min + &
                                                  delta_epsilon * ( igeo - 1 )
              epsilon_voigt(4,3*ngeo_strain+istart+igeo) = epsilon_min_off + &
                                         2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(5,4*ngeo_strain +istart+igeo) = epsilon_min_off + &
                                         2.0_DP * delta_epsilon * ( igeo - 1 )
              epsilon_voigt(6,5*ngeo_strain+istart+igeo) = epsilon_min_off + &
                                         2.0_DP * delta_epsilon * ( igeo - 1 )
           ENDDO
           istart=istart + nstep * ngeo_strain
        ENDDO
   END SELECT
ELSE
   nstep = 6
   nwork = nstep * ngeo_strain * ngeom
   ALLOCATE( epsilon_voigt(6, nwork) )
   epsilon_voigt=0.0_DP
   istart=0
   DO igeom=1,ngeom
      DO igeo=1,ngeo_strain
         epsilon_voigt(1, istart+igeo) = epsilon_min + delta_epsilon * ( igeo - 1 )
         epsilon_voigt(2, ngeo_strain+istart+igeo) = epsilon_min + &
                                             delta_epsilon * ( igeo - 1 )
         epsilon_voigt(3, 2*ngeo_strain+istart+igeo) = epsilon_min + &
                                             delta_epsilon * ( igeo - 1 )
         epsilon_voigt(4, 3*ngeo_strain +istart+igeo) = epsilon_min_off + &
                                    2.0_DP * delta_epsilon * ( igeo - 1 )
         epsilon_voigt(5, 4*ngeo_strain +istart+igeo) = epsilon_min_off + &
                                    2.0_DP * delta_epsilon * ( igeo - 1 )
         epsilon_voigt(6, 5*ngeo_strain +istart+igeo) = epsilon_min_off + &
                                    2.0_DP * delta_epsilon * ( igeo - 1 )
      ENDDO
      istart=istart + nstep * ngeo_strain
   ENDDO
ENDIF

work_base= nstep * ngeo_strain

CALL allocate_piezo(nwork)
ALLOCATE( epsilon_geo(3, 3, nwork) )
ALLOCATE( epsil_geo(nwork) )
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

iwork=0
DO igeom=1,ngeom
   DO istep=1,nstep
      DO igeo=1,ngeo_strain
         iwork=iwork+1
         epsil = epsilon_min + delta_epsilon * ( igeo - 1 )
         epsil_geo(iwork) = epsil
      ENDDO
   ENDDO
ENDDO
RETURN
END SUBROUTINE set_piezo_tensor_work
