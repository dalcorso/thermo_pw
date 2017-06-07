!
! Copyright (C) 2014-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_elastic_cons_work( nwork )
!
!  This routine sets the work to do for computing the elastic constants for
!  each Laue class. 
!  With the standard algorithm it computes only the strain (both in Voigt 
!  and in matrix form), while with the advanced and energy algorithms 
!  it computes also the Bravais lattice and the new celldm parameters 
!  for each strain (only for the lattices for which this option is available).
!  In the latter case it gives also the rotation matrix that brings the
!  coordinates of a point in the unstrained cartesian axis into the
!  coordinates in the cartesian axis of the strained system.
!  The strain is given in the coordinates axis of the unstrained system.
!
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : ibrav_geo, celldm_geo
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain, rot_mat, &
                               elastic_algorithm, epsilon_0
USE initial_conf,      ONLY : ibrav_save
USE equilibrium_conf,  ONLY : celldm0
USE thermo_sym,        ONLY : laue
!
!  library helper modules
!
USE elastic_constants, ONLY : epsilon_voigt, sigma_geo, epsilon_geo
USE strain_mod,        ONLY : set_strain_adv, trans_epsilon

IMPLICIT NONE
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsil
INTEGER  :: igeo, iwork, base_ind, istep, nstep
CHARACTER(LEN=2) :: strain_list(9)

nstep=0
SELECT CASE (laue) 
   CASE(29,32)
!
!  cubic system (T_h or O_h Laue classes)
!
      IF (ibrav_save==1.OR.ibrav_save==2.OR.ibrav_save==3) THEN    
         IF (elastic_algorithm=='standard'.OR.&
                              elastic_algorithm=='advanced') THEN
            nstep = 2
            strain_list(1) = 'E '
            strain_list(2) = 'F3'
            IF (ibrav_save==1) strain_list(2) = 'F '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 3
            strain_list(1) = 'A '
            strain_list(2) = 'E '
            strain_list(3) = 'F3'
            IF (ibrav_save==1) strain_list(3) = 'F '
         ENDIF
      ELSE
         CALL errore('set_elastic_cons_work',&
                              'Uncorrect lattice for cubic system',1)
      ENDIF
   CASE (19,23)
!
!  hexagonal system (C_6h or D_6h Laue classes)
!
      IF (ibrav_save==4) THEN    
         IF (elastic_algorithm=='standard'.OR.&
                              elastic_algorithm=='advanced') THEN
            nstep = 3
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'H '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 5
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'B1'
            strain_list(4) = 'A '
            strain_list(5) = 'H '
         ENDIF
      ELSE
         CALL errore('set_elastic_cons_work',&
                       'Uncorrect lattice for hexagonal system',ibrav_save)
      ENDIF
   CASE (18,22)
!
!   tetragonal system, (C_4h or D_4h Laue classes)
!
      IF (ibrav_save==6.OR.ibrav_save==7) THEN
         IF (elastic_algorithm=='standard'.OR.&
                              elastic_algorithm=='advanced') THEN
            nstep = 4
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'H '
            strain_list(4) = 'G '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 6
            strain_list(1) = 'E '
            strain_list(2) = 'C '
            strain_list(3) = 'B '
            strain_list(4) = 'B1'
            strain_list(5) = 'G '
            strain_list(6) = 'H '
         ENDIF
      ELSE
         CALL errore('set_elastic_cons_work',&
                       'Uncorrect lattice for tetragonal system',ibrav_save)
      ENDIF
   CASE (20)
!
!   orthorhombic system, (D_2h Laue class)
!
      IF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10&
                                                    .OR.ibrav_save==11) THEN
         IF (elastic_algorithm=='standard'.OR.&
                              elastic_algorithm=='advanced') THEN
            nstep = 6
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'G '
            strain_list(5) = 'H '
            strain_list(6) = 'I '
         ELSEIF (elastic_algorithm=='energy') THEN
            nstep = 9
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
         CALL errore('set_elastic_cons_work',&
                    'Uncorrect lattice for orthorombic system',ibrav_save)
      ENDIF
   CASE (25,27)
!
!   D_3d and S_6 Trigonal system. Only the standard algorithm is available
!
      IF (ibrav_save==5) THEN
         IF (elastic_algorithm=='standard'&
                       .OR.elastic_algorithm=='advanced') THEN
            nstep = 3
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'I '
         END IF
      ELSE
         CALL errore('set_elastic_cons_work',&
                    'Uncorrect lattice for trigonal system',ibrav_save)
      END IF
   CASE (16)
!
!   C_2h. Monoclinic system. Only the standard or advanced algorithm 
!         are available
!
      IF (ibrav_save==12.OR.ibrav_save==13.OR.ibrav_save==-12 &
                                          .OR.ibrav_save==-13) THEN
         ! both b and c-unique
         IF (elastic_algorithm=='standard' &
                .OR.elastic_algorithm=='advanced') THEN
            nstep = 6
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'G '
            strain_list(5) = 'H '
            strain_list(6) = 'I '
         END IF
      ELSE
         CALL errore('set_elastic_cons_work',&
                    'Uncorrect lattice for monoclinic system',ibrav_save)
      END IF

   CASE (2)
!
!   C_i   Triclinic system. Only the standard algorithm is available
!
      IF (ibrav_save==14.OR.ibrav_save==0) THEN
         IF (elastic_algorithm=='standard') THEN
            nstep = 6
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'G '
            strain_list(5) = 'H '
            strain_list(6) = 'I '
         END IF
      ELSE
         CALL errore('set_elastic_cons_work',&
                    'Uncorrect lattice for triclinic system',ibrav_save)
      END IF
   CASE DEFAULT
      CALL errore('set_elastic_cons_work','Laue class not available',1)
END SELECT
IF (nstep<1) CALL errore('set_elastic_cons_work', 'Uncorrect nstep, &
                                             &use standard algorithm',1)
nwork = nstep * ngeo_strain

ALLOCATE( epsilon_voigt(6, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
ALLOCATE( sigma_geo(3, 3, nwork) )
ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )
ALLOCATE( rot_mat(3,3,nwork) )

epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP - epsilon_0
DO istep=1,nstep
   base_ind = (istep-1) * ngeo_strain
   DO igeo=1,ngeo_strain
      epsil=epsilon_min + delta_epsilon * ( igeo - 1 )
      IF (igeo > ngeo_strain/2) epsil=epsil + 2.0_DP*epsilon_0
      IF (MOD(ngeo_strain,2)==1 .AND. igeo==(ngeo_strain/2 + 1)) epsil=epsil &
                                                          -epsilon_0

      CALL set_strain_adv(strain_list(istep), ibrav_save, celldm0, &
           epsil, epsilon_voigt(1,base_ind+igeo), ibrav_geo(base_ind+igeo), &
           celldm_geo(1,base_ind+igeo), rot_mat(1,1,base_ind+igeo) )
   ENDDO
ENDDO

sigma_geo=0.0_DP
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

RETURN
END SUBROUTINE set_elastic_cons_work
