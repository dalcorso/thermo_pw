!
! Copyright (C) 2014-2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE initialize_elastic_cons( ngeom, nwork )
!----------------------------------------------------------------------
!
!  This routine sets the work to do for computing the elastic constants for
!  each Laue class. 
!  With the standard or energy_std algorithm it computes only the strain 
!  (both in Voigt and in matrix form), while with the advanced and energy 
!  algorithms it computes also the Bravais lattice and the new celldm 
!  parameters for each strain (only for the lattices for which this 
!  option is available).
!  In the latter case it gives also the rotation matrix that brings the
!  coordinates of a point in the unstrained cartesian axis into the
!  coordinates in the cartesian axis of the strained system.
!  The strain is given in the coordinates axis of the unstrained system.
!  It sets the following variables
!  nstep_ec    ! number of strain types
!  nwork       ! total number of scf calculations to do 
!  work_base   ! total number of scf calculation for one equilibrium geometry
!  epsilon_voigt(6,nwork) ! the strain of each run in voigt notation
!  epsilon_geo(3,3,nwork) ! the strain of each run
!  epsil_geo(nwork) ! the amplitude of the strain of each run
!  
!  advanced or energy algorithm
!  ibrav_geo(nwork)  ! the bravais lattice code of each run
!  celldm_geo(6,nwork) ! the crystal parameters of each run
!  rot_mat(3,3,nwork) ! the possible rotation with respect to the orientation
!                     ! of the equilibrium cell 
!
!  If also atoms are optimized in each run
!  tau_acc(3,nat,nwork)   ! the atomic coordinates of each run
!                         ! this is supposed to be only the displacement
!                         ! with respect to the uniformely strained atomic
!                         ! positions. These displacements are given in
!                         ! cartesian coordinates in units of the 
!                         ! unperturbed celldm(1).
!
!  Moreover it allocates space for
!  sigma_geo(3,3,nwork)   ! the stress of each run
!  min_y(max_nint_var,ngeo_strain,21,ngeom) ! the internal coordinates 
                               ! that minimize the energy for each strain 
!  epsil_y(ngeo_strain,21,ngeom) ! the strain amplitude of each run.
!  
!  
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : ibrav_geo, celldm_geo, at_geo, tau_geo, &
                              uint_geo, freq_geo, z_geo, epsilon_zero_geo, &
                              epsilon_zerom1_geo
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain, rot_mat, &
                               elastic_algorithm,                          &
                               el_con_ibrav_geo, el_con_celldm_geo,        &
                               work_base, elalgen, epsil_geo, tau_acc,     &
                               nmove, atom_step, atom_dir, move_at, stype, &
                               nstep_ec, min_y, lcm_ec, epsil_y, old_ec,   &
                               nint_var_ec,                                &
                               stypec, int_ngeo_ec, int_step_ngeo_ec,      &
                               ninternal_ec 
USE control_atomic_pos, ONLY : max_nint_var
USE initial_conf,      ONLY : ibrav_save
USE equilibrium_conf,  ONLY : celldm0
USE thermo_sym,        ONLY : laue
!
!  library helper modules
!
USE elastic_constants, ONLY : epsilon_voigt, sigma_geo, epsilon_geo
USE strain_mod,        ONLY : set_strain_adv, trans_epsilon

USE ions_base,         ONLY : nat, ityp, amass

IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeom
INTEGER, INTENT(OUT) :: nwork
REAL(DP) :: epsilon_min, epsil, tot_mass, sumdisp
INTEGER  :: igeo, iwork, igeom, istep, nstep, iwork_base, imove, &
            nstep_tot, na, ipol, istart, ivar
CHARACTER(LEN=2) :: strain_list(21)
LOGICAL :: flag, lcons

nstep=0
SELECT CASE (laue) 
   CASE(29,32)
!
!  cubic system (T_h or O_h Laue classes)
!
      IF (ibrav_save==1.OR.ibrav_save==2.OR.ibrav_save==3) THEN    
         IF (.NOT.elalgen) THEN
            nstep = 2
            strain_list(1) = 'E '
            strain_list(2) = 'F3'
            IF (ibrav_save==1) strain_list(2) = 'F '
         ELSE
            nstep = 3
            strain_list(1) = 'A '
            strain_list(2) = 'E '
            strain_list(3) = 'F3'
            IF (ibrav_save==1) strain_list(3) = 'F '
         ENDIF
      ELSE
         CALL errore('initialize_elastic_cons',&
                              'Incorrect lattice for cubic system',1)
      ENDIF
   CASE (19,23)
!
!  hexagonal system (C_6h or D_6h Laue classes)
!
      IF (ibrav_save==4) THEN    
         IF (.NOT.elalgen) THEN
            nstep = 3
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'H '
         ELSE
            nstep = 5
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'B1'
            IF (old_ec==0) THEN
               strain_list(4) = 'A '
            ELSEIF (old_ec==1) THEN
               strain_list(4) = 'B '
            ELSEIF (old_ec==2) THEN
               strain_list(3)='B '
               strain_list(4)='A '
            ENDIF
            strain_list(5) = 'H '
         ENDIF
      ELSE
         CALL errore('initialise_elastic_cons',&
                       'Incorrect lattice for hexagonal system',ibrav_save)
      ENDIF
   CASE (18,22)
!
!   tetragonal system, (C_4h or D_4h Laue classes)
!
      IF (ibrav_save==6.OR.ibrav_save==7) THEN
         IF (.NOT.elalgen) THEN
            nstep = 4
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'H '
            strain_list(4) = 'G '
         ELSE
            nstep = 6
            strain_list(1) = 'E '
            strain_list(2) = 'C '
            strain_list(3) = 'B '
            strain_list(4) = 'B1'
            strain_list(5) = 'G '
            strain_list(6) = 'H '
            IF (laue==18) THEN
               nstep = 7
               strain_list(7) = 'CG'
               IF (elastic_algorithm=='energy') nstep=0
            ENDIF
         ENDIF
      ELSE
         CALL errore('initialize_elastic_cons',&
                       'Incorrect lattice for tetragonal system',ibrav_save)
      ENDIF
   CASE (20)
!
!   orthorhombic system, (D_2h Laue class)
!
      IF (ibrav_save==8.OR.ibrav_save==9.OR.ibrav_save==10&
                                                    .OR.ibrav_save==11) THEN
         IF (.NOT.elalgen) THEN
            nstep = 6
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'G '
            strain_list(5) = 'H '
            strain_list(6) = 'I '
         ELSE
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
         CALL errore('initialize_elastic_cons',&
                    'Incorrect lattice for orthorombic system',ibrav_save)
      ENDIF
   CASE (25,27)
!
!   D_3d and S_6 Trigonal system (hexagonal or trigonal lattice). 
!
      IF (ibrav_save==4.OR.ibrav_save==5) THEN
         IF (.NOT.elalgen) THEN
            nstep = 3
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'H '
            IF (ibrav_save==5) strain_list(3) = 'I '
         ELSE
            nstep = 6
            strain_list(1) = 'C '
            strain_list(2) = 'E '
            strain_list(3) = 'B1'
            strain_list(4) = 'A '
            strain_list(5) = 'H '
            strain_list(6) = 'CI'
            IF (laue==27) THEN
               nstep = 7
               strain_list(7) = 'CG'
            ENDIF
            IF (elastic_algorithm=='energy') nstep=0
         END IF
      ELSE
         CALL errore('initialize_elastic_cons',&
                    'Incorrect lattice for trigonal system',ibrav_save)
      END IF
   CASE (16)
!
!   C_2h. Monoclinic system. Only the standard or advanced algorithm 
!         are available
!
      IF (ibrav_save==12.OR.ibrav_save==13.OR.ibrav_save==-12 &
                                          .OR.ibrav_save==-13) THEN
         ! both b and c-unique
         IF (.NOT.elalgen) THEN
            nstep = 6
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'G '
            strain_list(5) = 'H '
            strain_list(6) = 'I '
         ELSEIF (elastic_algorithm=='energy_std') THEN
            nstep=13
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'B '
            strain_list(5) = 'B1'
            strain_list(6) = 'B2'
            strain_list(7) = 'G '
            strain_list(8) = 'H '
            strain_list(9) = 'I '
            IF (ibrav_save>0) THEN
!
!   c unique
!
               strain_list(10) = 'CG'
               strain_list(11) = 'DG'
               strain_list(12) = 'EG'
               strain_list(13) = 'IH'
            ELSE
!
!   b unique
!
               strain_list(10) = 'CH'
               strain_list(11) = 'DH'
               strain_list(12) = 'EH'
               strain_list(13) = 'GI'
            ENDIF
         END IF
      ELSE
         CALL errore('initialize_elastic_cons',&
                    'Incorrect lattice for monoclinic system',ibrav_save)
      END IF

   CASE (2)
!
!   C_i   Triclinic system. 
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
         ELSEIF (elastic_algorithm=='energy_std') THEN
            nstep = 21
            strain_list(1) = 'C '
            strain_list(2) = 'D '
            strain_list(3) = 'E '
            strain_list(4) = 'B '
            strain_list(5) = 'B1'
            strain_list(6) = 'B2'
            strain_list(7) = 'G '
            strain_list(8) = 'H '
            strain_list(9) = 'I '
            strain_list(10) = 'CG'
            strain_list(11) = 'CH'
            strain_list(12) = 'CI'
            strain_list(13) = 'DG'
            strain_list(14) = 'DH'
            strain_list(15) = 'DI'
            strain_list(16) = 'EG'
            strain_list(17) = 'EH'
            strain_list(18) = 'EI'
            strain_list(19) = 'GH'
            strain_list(20) = 'IH'
            strain_list(21) = 'GI'
         END IF
      ELSE
         CALL errore('initialize_elastic_cons',&
                    'Incorrect lattice for triclinic system',ibrav_save)
      END IF
   CASE DEFAULT
      CALL errore('initialize_elastic_cons','Laue class not available',1)
END SELECT
IF (nstep<1.AND.elastic_algorithm=='advanced') &
   CALL errore('initialize_elastic_cons', 'Incorrect nstep, &
                                             &use standard algorithm',1)
IF (nstep<1.AND.elastic_algorithm=='energy') &
   CALL errore('initialize_elastic_cons', 'Incorrect nstep, &
                                             &use energy_std algorithm',1)
IF (nstep<1) &
   CALL errore('initialize_elastic_cons', 'Incorrect nstep, &
                                             &check elastic_algorithm',1)
!
!  Compute the total number of perturbed geometries per unperturbed
!  one not accounting for ngeo_strain, but accounting for the possibility
!  that in some strain types we do a nmove or ninternal_ec(istep) 
!  internal parameter displacements.
!
nstep_tot=0
DO istep=1,nstep
   IF (stype(istep)) THEN
      nstep_tot=nstep_tot+nmove
      int_ngeo_ec(:,istep)=1
   ELSEIF (stypec(istep)) THEN
      ninternal_ec(istep)=1
      CALL clean_int_ngeo(int_ngeo_ec(1,istep), nint_var_ec(istep))
      DO ivar=1,nint_var_ec(istep)
         ninternal_ec(istep)=ninternal_ec(istep)*int_ngeo_ec(ivar,istep)
      ENDDO
      nstep_tot=nstep_tot+ninternal_ec(istep)
   ELSE
      nstep_tot=nstep_tot+1
      int_ngeo_ec(:,istep)=1
   ENDIF
ENDDO
!
!   Total number of distorted geometries
!
nwork = nstep_tot * ngeo_strain * ngeom
!
!   Total number of distorted geometries per unperturbed geometry
!
work_base = nstep_tot * ngeo_strain
nstep_ec=nstep
ALLOCATE( epsilon_voigt(6, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
ALLOCATE( epsil_geo(nwork) )
ALLOCATE( sigma_geo(3, 3, nwork) )
ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )
ALLOCATE( at_geo(3,3,nwork) )
ALLOCATE( tau_geo(3,nat,nwork) )
ALLOCATE( rot_mat(3,3,nwork) )
ALLOCATE( uint_geo(max_nint_var,nwork) )
ALLOCATE( tau_acc(3,nat,nwork) )
ALLOCATE( min_y(max_nint_var,ngeo_strain,21,ngeom))
ALLOCATE( epsil_y(ngeo_strain,21,ngeom))
ALLOCATE( epsilon_zero_geo(3,3,nwork) )
ALLOCATE( epsilon_zerom1_geo(3,3,nwork) )
ALLOCATE( freq_geo(3*nat,nwork) )
ALLOCATE( z_geo(3*nat,3*nat,nwork) )

IF (ngeom==1) THEN
   IF (.NOT. ALLOCATED(el_con_ibrav_geo)) ALLOCATE(el_con_ibrav_geo(1))
   IF (.NOT. ALLOCATED(el_con_celldm_geo)) ALLOCATE(el_con_celldm_geo(6,1))
   el_con_ibrav_geo(1)=ibrav_save
   el_con_celldm_geo(:,1)=celldm0(:)
ENDIF
flag=(elastic_algorithm=='standard'.OR.elastic_algorithm=='energy_std')
epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP 
iwork=0
tau_acc=0.0_DP
DO igeom=1, ngeom
   DO istep=1,nstep
      DO igeo=1,ngeo_strain
         iwork=iwork+1
         epsil = epsilon_min + delta_epsilon * ( igeo - 1 )
         epsil_geo(iwork) = epsil 
         CALL set_strain_adv(strain_list(istep), el_con_ibrav_geo(igeom),  &
              el_con_celldm_geo(1,igeom), epsil, &
              epsilon_voigt(1,iwork), ibrav_geo(iwork), &
              celldm_geo(1,iwork), rot_mat(1,1,iwork), flag )
         IF (stype(istep)) THEN
            iwork_base=iwork
            tau_acc(:,move_at(istep),iwork)=(1.0_DP-(nmove+1.0_DP)/2.0_DP)* &
                                         atom_step(istep)*atom_dir(:,istep)
            DO imove=2, nmove
               iwork=iwork+1
               epsil_geo(iwork)=epsil_geo(iwork_base)
               epsilon_voigt(:,iwork)=epsilon_voigt(:,iwork_base)
               ibrav_geo(iwork)=ibrav_geo(iwork_base)
               celldm_geo(:,iwork)=celldm_geo(:,iwork_base)
               rot_mat(:,:,iwork)=rot_mat(:,:,iwork_base)
               tau_acc(:,move_at(istep),iwork)=(imove-(nmove+1.0_DP)/2.0_DP)* &
                                         atom_step(istep)*atom_dir(:,istep)
            ENDDO
         ENDIF
         IF (stypec(istep)) THEN
            iwork_base=iwork
            DO imove=2, ninternal_ec(istep)
               iwork=iwork+1
               epsil_geo(iwork)=epsil_geo(iwork_base)
               epsilon_voigt(:,iwork)=epsilon_voigt(:,iwork_base)
               ibrav_geo(iwork)=ibrav_geo(iwork_base)
               celldm_geo(:,iwork)=celldm_geo(:,iwork_base)
               rot_mat(:,:,iwork)=rot_mat(:,:,iwork_base)
            ENDDO
         ENDIF
      ENDDO
!
!   set separately the atomic positions for all ngeo_strain*ninternal_ec(istep)
!   geometries
!
      IF (stypec(istep)) THEN
         istart=iwork-ngeo_strain*ninternal_ec(istep)+1
         CALL set_tau_acc(celldm_geo(1,istart), tau_acc(1,1,istart), &
                  uint_geo(1,istart), ngeo_strain*ninternal_ec(istep), nat, &
                  nint_var_ec(istep), istep)
      ENDIF
   ENDDO
ENDDO
!
!  Correct the displacement so that it does not change the center of mass
!  of the system
!
tot_mass=0.0_DP
DO na=1,nat
   tot_mass=tot_mass+amass(ityp(na))
ENDDO
IF ((tot_mass>0.0_DP).AND.lcm_ec) THEN
   DO iwork=1,nwork
      DO ipol=1,3
         sumdisp=0.0_DP
         DO na=1,nat
            sumdisp=sumdisp+tau_acc(ipol,na,iwork) * amass(ityp(na))
         ENDDO
         DO na=1,nat
            tau_acc(ipol,na,iwork)=tau_acc(ipol,na,iwork)-sumdisp/tot_mass
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
!  tau_acc in units of celldm0(1)
!
tau_acc=tau_acc/celldm0(1)
!
!  Compute the strain tensor from epsilon_voigt
!
sigma_geo=0.0_DP
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

IF (elastic_algorithm=='standard'.OR.elastic_algorithm=='energy_std') THEN
   ibrav_geo=0
   celldm_geo=0.0_DP
ENDIF

RETURN
END SUBROUTINE initialize_elastic_cons
