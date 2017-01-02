!
! Copyright (C) 2014-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_elastic_cons_work( nwork )

USE kinds,      ONLY : DP
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain, epsilon_0

USE thermo_sym, ONLY : laue
!
!   library helper routine
!
USE elastic_constants, ONLY : epsilon_voigt, epsilon_geo, sigma_geo, &
                       trans_epsilon

IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsilon_min_off
INTEGER  :: igeo, iwork

epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP - epsilon_0
epsilon_min_off= - delta_epsilon * (ngeo_strain - 1 ) - 2.0_DP * epsilon_0

IF (ALLOCATED(epsilon_voigt)) DEALLOCATE(epsilon_voigt)
IF (ALLOCATED(sigma_geo))     DEALLOCATE(sigma_geo)
IF (ALLOCATED(epsilon_geo))   DEALLOCATE(epsilon_geo)
IF (ALLOCATED(ibrav_geo))     DEALLOCATE(ibrav_geo)
IF (ALLOCATED(celldm_geo))    DEALLOCATE(celldm_geo)

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
         IF (igeo> ngeo_strain/2) THEN
            epsilon_voigt(3, igeo)=epsilon_voigt(3, igeo)+2.0_DP*epsilon_0
            epsilon_voigt(4, ngeo_strain+igeo)= &
                         epsilon_voigt(4, ngeo_strain+igeo)+4.0_DP*epsilon_0
         ENDIF
         IF (MOD(ngeo_strain,2)==1 .AND. igeo==ngeo_strain/2 + 1) THEN 
            epsilon_voigt(3, igeo)=epsilon_voigt(3, igeo)-epsilon_0
            epsilon_voigt(4, ngeo_strain+igeo)= &
                         epsilon_voigt(4, ngeo_strain+igeo)-2.0_DP*epsilon_0
         ENDIF
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
         IF (igeo> ngeo_strain/2) THEN
            epsilon_voigt(1, igeo)=epsilon_voigt(1, igeo) + 2.0_DP*epsilon_0
            epsilon_voigt(3, ngeo_strain+igeo)=&
                     epsilon_voigt(3, ngeo_strain+igeo) + 2.0_DP*epsilon_0
            epsilon_voigt(4, 2*ngeo_strain + igeo)= &
                     epsilon_voigt(4, 2*ngeo_strain + igeo) + 4.0_DP*epsilon_0
         ENDIF
         IF (MOD(ngeo_strain,2)==1 .AND. igeo==ngeo_strain/2 + 1) THEN 
            epsilon_voigt(1, igeo)=epsilon_voigt(1, igeo) - epsilon_0
            epsilon_voigt(3, ngeo_strain+igeo)=&
                     epsilon_voigt(3, ngeo_strain+igeo) - epsilon_0
            epsilon_voigt(4, 2*ngeo_strain + igeo)= &
                     epsilon_voigt(4, 2*ngeo_strain + igeo) - 2.0_DP*epsilon_0
         ENDIF
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
         IF (igeo> ngeo_strain/2) THEN
            epsilon_voigt(1, igeo)=epsilon_voigt(1, igeo)+2.0_DP*epsilon_0
            epsilon_voigt(3, ngeo_strain + igeo)=&
                          epsilon_voigt(3, ngeo_strain+igeo)+2.0_DP*epsilon_0
            epsilon_voigt(4, 2*ngeo_strain + igeo)=&
                          epsilon_voigt(4, 2*ngeo_strain+igeo)+4.0_DP*epsilon_0
            epsilon_voigt(6, 3*ngeo_strain+igeo)=&
                          epsilon_voigt(6, 3*ngeo_strain+igeo)+4.0_DP*epsilon_0
         ENDIF
         IF (MOD(ngeo_strain,2)==1 .AND. igeo==ngeo_strain/2 + 1) THEN
            epsilon_voigt(1, igeo)=epsilon_voigt(1, igeo) - epsilon_0
            epsilon_voigt(3, ngeo_strain + igeo)=&
                          epsilon_voigt(3, ngeo_strain+igeo) - epsilon_0
            epsilon_voigt(4, 2*ngeo_strain + igeo)=&
                          epsilon_voigt(4, 2*ngeo_strain+igeo)-2.0_DP*epsilon_0
            epsilon_voigt(6, 3*ngeo_strain+igeo)=&
                          epsilon_voigt(6, 3*ngeo_strain+igeo)-2.0_DP*epsilon_0
         ENDIF
      ENDDO
CASE DEFAULT
!
!   generic case, used when ibrav=0 or for triclinic, monoclinic, and 
!   orthorhombic systems, no information is deduced from symmetry 
!   in the first three cases, all 21 elements of the elastic constants 
!   matrix are computed. 
!   Requires 6 * ngeo_strain self consistent calculations. In the monoclinic
!   and orthorhombic case, some elements vanish by symmetry, but we need to
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
      IF (igeo> ngeo_strain/2) THEN
         epsilon_voigt(1, igeo)=epsilon_voigt(1, igeo)+2.0_DP*epsilon_0
         epsilon_voigt(1, ngeo_strain+igeo)=&
                       epsilon_voigt(1, ngeo_strain+igeo)+2.0_DP*epsilon_0
         epsilon_voigt(3, 2*ngeo_strain + igeo)=&
                       epsilon_voigt(3, 2*ngeo_strain+igeo)+2.0_DP*epsilon_0
         epsilon_voigt(4, 3*ngeo_strain + igeo)=&
                       epsilon_voigt(4, 3*ngeo_strain+igeo)+4.0_DP*epsilon_0
         epsilon_voigt(5, 4*ngeo_strain+igeo)=&
                       epsilon_voigt(5, 4*ngeo_strain+igeo)+4.0_DP*epsilon_0
         epsilon_voigt(6, 5*ngeo_strain+igeo)=&
                       epsilon_voigt(6, 5*ngeo_strain+igeo)+4.0_DP*epsilon_0
      ENDIF
      IF (MOD(ngeo_strain,2)==1 .AND. igeo==ngeo_strain/2 + 1) THEN
         epsilon_voigt(1, igeo)=epsilon_voigt(1, igeo) - epsilon_0
         epsilon_voigt(1, ngeo_strain+igeo)=&
                       epsilon_voigt(1, ngeo_strain+igeo) - epsilon_0
         epsilon_voigt(3, 2*ngeo_strain + igeo)=&
                       epsilon_voigt(3, 2*ngeo_strain+igeo) - epsilon_0
         epsilon_voigt(4, 3*ngeo_strain + igeo)=&
                       epsilon_voigt(4, 3*ngeo_strain+igeo)-2.0_DP*epsilon_0
         epsilon_voigt(5, 4*ngeo_strain+igeo)=&
                       epsilon_voigt(5, 4*ngeo_strain+igeo)-2.0_DP*epsilon_0
         epsilon_voigt(6, 5*ngeo_strain+igeo)=&
                       epsilon_voigt(6, 5*ngeo_strain+igeo)-2.0_DP*epsilon_0
      ENDIF
   ENDDO
END SELECT

ALLOCATE( sigma_geo(3, 3, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
sigma_geo=0.0_DP
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )

RETURN
END SUBROUTINE set_elastic_cons_work

SUBROUTINE set_elastic_cons_work_adv( nwork )
!
!  With the option advanced this routine computes the work to do,
!  the Bravais lattice and the new celldm parameters for each strain.
!  Moreover it can calculate the rotation matrix between the unstrained and
!  strained cartesian axes.
!  Moreover it sets the strain matrix for each geometery, both in 
!  Voigt notation and in strain form.
!
USE kinds,      ONLY : DP
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain, rot_mat, &
                              aap_mat, apa_mat, elastic_algorithm, epsilon_0


USE initial_conf, ONLY : ibrav_save
USE equilibrium_conf, ONLY : celldm0
USE thermo_sym, ONLY : laue
!
!  library helper modules
!
USE elastic_constants, ONLY : epsilon_voigt, sigma_geo, epsilon_geo, &
                              trans_epsilon
USE strain_mod, ONLY : apply_strain_adv

IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsil
INTEGER :: igeo, iwork, base_ind, i, j, istep, nstep
CHARACTER(LEN=2) :: strain_list(9)

epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP - epsilon_0
SELECT CASE (laue) 
   CASE(29,32)
!
!  cubic system
!
      IF (ibrav_save==1.OR.ibrav_save==2.OR.ibrav_save==3) THEN    
         IF (elastic_algorithm=='advanced') THEN
            nstep = 2
            nwork = nstep * ngeo_strain
            strain_list(1) = 'C '
            strain_list(2) = 'F '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 3
            nwork = nstep * ngeo_strain
            strain_list(1) = 'A '
            strain_list(2) = 'C '
            strain_list(3) = 'F '
         ENDIF
      ELSE
         CALL errore('set_elastic_cons_work_adv',&
                                   'Bravais lattice not available',1)
      ENDIF
   CASE (23)
!
!  hexagonal system
!
      IF (ibrav_save==4) THEN    
         IF (elastic_algorithm=='advanced') THEN
            nstep = 3
            nwork = nstep * ngeo_strain
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'H '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 5
            nwork = nstep * ngeo_strain
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'B1'
            strain_list(4) = 'A '
            strain_list(5) = 'H '
         ENDIF
      ELSE
         CALL errore('set_elastic_cons_work_adv',&
                                   'Bravais lattice not available',ibrav_save)
      ENDIF
   CASE (18,22)
!
!   tetragonal system
!
      IF (ibrav_save==6.OR.ibrav_save==7) THEN
         IF (elastic_algorithm=='advanced') THEN
            nstep = 4
            nwork = nstep * ngeo_strain
            strain_list(1) = 'E '
            strain_list(2) = 'C '
            strain_list(3) = 'G '
            strain_list(4) = 'H '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 6
            nwork = nstep * ngeo_strain
            strain_list(1) = 'E '
            strain_list(2) = 'C '
            strain_list(3) = 'B '
            strain_list(4) = 'B1'
            strain_list(5) = 'G '
            strain_list(6) = 'H '
         ENDIF
      ELSE
         CALL errore('set_elastic_cons_work_adv',&
                                   'Bravais lattice not correct',ibrav_save)
      ENDIF
   CASE (20,0)
!
!   orthorhombic system
!
      IF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10.OR.ibrav_save==11) THEN
         IF (elastic_algorithm=='advanced') THEN
            nstep = 6
            nwork = nstep * ngeo_strain
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'G '
            strain_list(5) = 'H '
            strain_list(6) = 'I '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 9
            nwork = nstep * ngeo_strain
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'B '
            strain_list(5) = 'B1'
            strain_list(6) = 'B2'
            strain_list(7) = 'G '
            strain_list(8) = 'H '
            strain_list(9) = 'I '
         ENDIF
      ELSE
         CALL errore('set_elastic_cons_work_adv',&
                    'Bravais lattice and Laue class not compatible',ibrav_save)
      ENDIF
   CASE DEFAULT
      CALL errore('set_elastic_cons_work_adv',&
                              'Bravais lattice not available',1)
END SELECT

ALLOCATE( epsilon_voigt(6, nwork) )
ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )
ALLOCATE( rot_mat(3,3,nwork) )
ALLOCATE( aap_mat(3,3,nwork) )
ALLOCATE( apa_mat(3,3,nwork) )

DO istep=1,nstep
   base_ind = (istep-1) * ngeo_strain
   DO igeo=1,ngeo_strain
      epsil=epsilon_min + delta_epsilon * ( igeo - 1 )
      IF (igeo > ngeo_strain/2) epsil=epsil + 2.0_DP*epsilon_0
      IF (MOD(ngeo_strain,2)==1 .AND. igeo==(ngeo_strain/2 + 1)) epsil=epsil &
                                                          -epsilon_0
      CALL apply_strain_adv(strain_list(istep), ibrav_save, celldm0, &
           epsil, ibrav_geo(base_ind+igeo), celldm_geo(1,base_ind+igeo), &
           epsilon_voigt(1,base_ind+igeo), rot_mat(1,1,base_ind+igeo), &
           aap_mat(1,1,base_ind+igeo), apa_mat(1,1,base_ind+igeo) )
   END DO
END DO

ALLOCATE( sigma_geo(3, 3, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
sigma_geo=0.0_DP
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

RETURN
END SUBROUTINE set_elastic_cons_work_adv

