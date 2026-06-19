!
! Copyright (C) 2026 Andrea Dal Corso and A. Ahmed
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE energy_for_ec(ngeom_)
!----------------------------------------------------------------------
!
!  This routine sets the distorted geometries to calculate the elastic 
!  constants for each Laue class. If the distorted geometry has the
!  same Bravais lattice of the equilibrium configuration the energy
!  of the distorted lattice can be obtained from the polynomial that
!  interpolates the energies on the grid. This routine uses this polynomial
!  and writes on file the distorted energies. These files are in
!  the restart directories and can be used by a calculation of 
!  elastic constant instead of doing another scf calculation.
!  This routines works only for elastic_algorithm='energy' or
!  elastic_algorithm='energy_std'.
!  It sets the following variables (that remain however local to this
!  routine and get lost as soon as we exit from here).
!  with_same_lattice ! strain types for which this routine can be used
!  nwork       ! total number of energies to calculate  
!              ! note: this is number account for all distorted geometries
!              !       and can be used as an index to the name of the file
!              !       where the interpolated energies are saved.
!  epsilon_voigt(6,nwork) ! the strain of each run in voigt notation
!  epsilon_geo(3,3,nwork) ! the strain of each run
!  epsil_geo(nwork) ! the amplitude of the strain of each run
!  
!  for energy algorithm
!  ibrav_geo(nwork)  ! the bravais lattice code of each run
!  celldm_geo(6,nwork) ! the crystal parameters of each run
!  rot_mat(3,3,nwork) ! the possible rotation with respect to the orientation
!                     ! of the equilibrium cell 
!  
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : ibrav_geo_=>ibrav_geo, celldm_geo_=>celldm_geo
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain,          &
                               elastic_algorithm, elalgen, old_ec
USE geometry_file,     ONLY : ngeo_file, celldm_geo_file
USE control_quadratic_energy, ONLY : p2
USE control_quartic_energy,   ONLY : p4
USE control_thermo,    ONLY : lgeo_to_file
USE initial_conf,      ONLY : ibrav_save
USE equilibrium_conf,  ONLY : celldm0
USE thermo_sym,        ONLY : laue
!
!  library helper modules
!
USE elastic_constants, ONLY : set_strain_for_ec
USE strain_mod,        ONLY : set_strain_adv, trans_epsilon


USE io_global,         ONLY : ionode, ionode_id
USE mp,                ONLY : mp_bcast
USE mp_images,         ONLY : intra_image_comm

IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeom_
REAL(DP) :: epsilon_min, epsil, omega
REAL(DP), ALLOCATABLE :: epsilon_voigt(:,:), epsilon_geo(:,:,:), &
                         epsil_geo(:), celldm_geo(:,:),          &
                         el_con_celldm_geo(:,:), rot_mat(:,:,:), &
                         energy_geo(:), el_con_at_geo(:,:,:), at_geo(:,:,:)
INTEGER, ALLOCATABLE :: ibrav_geo(:), el_con_ibrav_geo(:)
INTEGER  :: igeo, iwork, igeom, istep, nstep, nstep_tot, part, &
            iu_ene, nwork, ios, ngeom
INTEGER  :: find_free_unit
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=2) :: strain_list(21)
INTEGER :: with_same_lattice(21)
LOGICAL :: flag

IF (elastic_algorithm/='energy'.AND.elastic_algorithm/='energy_std') &
   CALL errore('energy_for_ec', 'Incorrect algorithm', 1)

CALL set_strain_for_ec(strain_list, nstep, with_same_lattice, laue, &
                        ibrav_save, old_ec, elalgen, elastic_algorithm)

ngeom=ngeom_
IF (lgeo_to_file) ngeom=ngeo_file

IF (nstep<1.AND.elastic_algorithm=='energy') &
   CALL errore('energy_for_ec', 'Incorrect nstep, use energy_std algorithm',1)

IF (nstep<1) &
   CALL errore('energy_for_ec', 'Incorrect nstep, &
                                             &check elastic_algorithm',1)
!
nstep_tot=nstep
!
!   Total number of distorted geometries
!
nwork = nstep_tot * ngeo_strain * ngeom
!
!   Total number of distorted geometries per unperturbed geometry
!
ALLOCATE( energy_geo(nwork) )
ALLOCATE( epsilon_voigt(6, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
ALLOCATE( epsil_geo(nwork) )
ALLOCATE( at_geo(3, 3, nwork) )
ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )
ALLOCATE( rot_mat(3,3,nwork) )

IF (ngeom==1) THEN
   IF (.NOT. ALLOCATED(el_con_ibrav_geo)) ALLOCATE(el_con_ibrav_geo(1))
   IF (.NOT. ALLOCATED(el_con_celldm_geo)) ALLOCATE(el_con_celldm_geo(6,1))
   IF (.NOT. ALLOCATED(el_con_at_geo)) ALLOCATE(el_con_at_geo(3,3,1))
   el_con_ibrav_geo(1)=ibrav_save
   el_con_celldm_geo(:,1)=celldm0(:)
   CALL latgen(ibrav_save,celldm0,at_geo(1,1,1),at_geo(1,2,1),at_geo(1,3,1),&
                                                                     omega)
   el_con_at_geo(:,:,1)= at_geo(:,:,1)
ELSE
   ALLOCATE(el_con_ibrav_geo(ngeom))
   ALLOCATE(el_con_celldm_geo(6,ngeom))
   ALLOCATE(el_con_at_geo(3,3,ngeom))
   DO igeom=1, ngeom
      IF (lgeo_to_file) THEN
         el_con_ibrav_geo(igeom)=ibrav_save
         el_con_celldm_geo(:,igeom)=celldm_geo_file(:,igeom)
      ELSE
         el_con_ibrav_geo(igeom)=ibrav_geo_(igeom)
         el_con_celldm_geo(:,igeom)=celldm_geo_(:,igeom)
      ENDIF
      CALL latgen(el_con_ibrav_geo(igeom), el_con_celldm_geo(1,igeom), &
           at_geo(1,1,igeom),at_geo(1,2,igeom), at_geo(1,3,igeom), omega)
      el_con_at_geo(:,:,igeom)=at_geo(:,:,igeom)
   ENDDO
ENDIF
flag=(elastic_algorithm=='energy_std')
epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP 
iwork=0
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
      ENDDO
   ENDDO
ENDDO
!
!  Compute the strain tensor from epsilon_voigt
!
epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO

iwork=0
part=1
DO igeom=1, ngeom
   DO istep=1,nstep
      DO igeo=1,ngeo_strain
         iwork=iwork+1
         IF (with_same_lattice(istep)==1) THEN
            IF (elastic_algorithm=='energy_std') THEN
!
!    here it should be easier to compute the celldm_geo because we are
!    considering only strain that do not change the lattice
!
               CALL compute_celldm_geo_from_strain(el_con_at_geo(1,1,igeom), &
                     epsilon_geo(1,1,iwork), celldm_geo(1,iwork),ibrav_save)
               ibrav_geo(iwork)=ibrav_save
            ENDIF
            WRITE(6,'(2i5,6f15.6)') iwork, ibrav_geo(iwork), &
                                          celldm_geo(1:6,iwork)
            CALL compute_energy_from_poly(ibrav_geo(iwork),&
                       celldm_geo(1,iwork), energy_geo(iwork),p2,p4)
            IF (iwork > 0.AND.ionode) THEN
               iu_ene=find_free_unit()
               filename='energy_dist/e_work_part.'//TRIM(int_to_char(iwork))//&
                                            '.'//TRIM(int_to_char(part))
               OPEN(UNIT=iu_ene, FILE=TRIM(filename), STATUS='UNKNOWN', &
                                 FORM='FORMATTED', ERR=20, IOSTAT=ios)

               WRITE(iu_ene,*) energy_geo(iwork)

               CLOSE(iu_ene,status='KEEP')
            ENDIF
20          CALL mp_bcast(ios, ionode_id, intra_image_comm)
            CALL errore('energy_for_ec','Opening energy file',ios)
         ENDIF
      ENDDO
   ENDDO
ENDDO 

DEALLOCATE( energy_geo )
DEALLOCATE( epsilon_voigt )
DEALLOCATE( epsilon_geo )
DEALLOCATE( epsil_geo )
DEALLOCATE( at_geo )
DEALLOCATE( ibrav_geo )
DEALLOCATE( celldm_geo )
DEALLOCATE( rot_mat )

RETURN
END SUBROUTINE energy_for_ec

!----------------------------------------------------------------------------
SUBROUTINE compute_energy_from_poly(ibrav, celldm, energy, p2, p4)
!----------------------------------------------------------------------------
!
!   This routine receives an ibrav and a celldm and computes the
!   values of the polynomial at celldm saving it into energy
!
USE kinds, ONLY : DP
USE control_quartic_energy, ONLY : lquartic
USE polynomial, ONLY : poly2, poly4 
USE control_mur, ONLY : lmurn
USE lattices,    ONLY : compress_celldm, crystal_parameters
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic
USE quartic_surfaces, ONLY : evaluate_fit_quartic
IMPLICIT NONE

INTEGER, INTENT(IN) :: ibrav
REAL(DP), INTENT(IN) :: celldm(6)
REAL(DP), INTENT(OUT) :: energy

TYPE(poly2), INTENT(IN) :: p2
TYPE(poly4), INTENT(IN) :: p4

INTEGER :: nvar
REAL(DP), ALLOCATABLE :: x(:)

nvar=crystal_parameters(ibrav)
ALLOCATE(x(nvar))

IF (lmurn) THEN
   x(1)=celldm(1)
ELSE
   CALL compress_celldm(celldm, x, nvar, ibrav)
ENDIF

IF (lquartic) THEN
   CALL evaluate_fit_quartic(nvar, x, energy, p4)
ELSE
   CALL evaluate_fit_quadratic(nvar, x, energy, p2)
ENDIF

DEALLOCATE(x)

RETURN
END SUBROUTINE compute_energy_from_poly

!-----------------------------------------------------------------------------
SUBROUTINE compute_celldm_geo_from_strain(el_cons_at, &
                     epsilon_geo, celldm_geo, ibrav_save)
!-----------------------------------------------------------------------------
USE kinds,            ONLY : DP
USE strain_mod,       ONLY : apply_strain
USE lattices,         ONLY : compute_conventional, lattice_parameters, &
                             set_celldm

IMPLICIT NONE
REAL(DP), INTENT(IN) :: el_cons_at(3,3), epsilon_geo(3,3)
REAL(DP), INTENT(OUT) :: celldm_geo(6)
INTEGER, INTENT(IN) :: ibrav_save

REAL(DP) :: ats(3,3), atp(3,3)
REAL(DP) :: a, b, c, alpha, beta, ggamma
INTEGER :: ivec

DO ivec=1,3
   CALL apply_strain(el_cons_at(1,ivec), ats(1,ivec), epsilon_geo)
ENDDO

CALL compute_conventional(ats, atp, ibrav_save )

CALL lattice_parameters(atp, a, b, c, alpha, beta, ggamma)

CALL set_celldm(ibrav_save, a, b, c, alpha, beta, ggamma, celldm_geo)

RETURN
END SUBROUTINE compute_celldm_geo_from_strain
