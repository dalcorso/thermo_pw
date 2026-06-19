!
! Copyright (C) 2026 Andrea Dal Corso and A. Ahmed
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE free_energy_for_ec(ngeom_)
!----------------------------------------------------------------------
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : ibrav_geo_=>ibrav_geo, celldm_geo_=>celldm_geo, &
                              omega_geo_=>omega_geo, no_ph, tot_ngeo
USE ions_base,         ONLY : nat
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain,          &
                               elastic_algorithm, elalgen, old_ec
USE thermodynamics,    ONLY : ph_ener_eos, ph_entropy_eos, ph_ce_eos, ph_e0
USE control_mur,       ONLY : lmurn
USE geometry_file,     ONLY : ngeo_file, celldm_geo_file
USE control_thermo,    ONLY : lgeo_to_file
USE initial_conf,      ONLY : ibrav_save
USE equilibrium_conf,  ONLY : celldm0
USE thermo_sym,        ONLY : laue
USE temperature,       ONLY : ntemp, temp
USE data_files,        ONLY : fltherm


USE elastic_constants, ONLY : set_strain_for_ec
USE strain_mod,        ONLY : set_strain_adv, trans_epsilon

! Polynomial evaluators
USE polynomial,        ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE linear_surfaces,   ONLY : evaluate_fit_linear, fit_multi_linear      
USE quadratic_surfaces,ONLY : evaluate_fit_quadratic, fit_multi_quadratic   
USE cubic_surfaces,    ONLY : evaluate_fit_cubic, fit_multi_cubic        
USE quartic_surfaces,  ONLY : evaluate_fit_quartic, fit_multi_quartic    
USE control_quartic_energy, ONLY : lsolve, poly_degree_ph

USE anharmonic, ONLY : p1t_t, p2t_t, p3t_t, p4t_t

USE control_mur,       ONLY : lmurn
USE lattices,          ONLY : compress_celldm, crystal_parameters

USE io_global,         ONLY : meta_ionode

IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeom_

REAL(DP) :: epsilon_min, epsil, omega, e0
REAL(DP), ALLOCATABLE :: epsilon_voigt(:,:), epsilon_geo(:,:,:), &
                         epsil_geo(:), celldm_geo(:,:),          &
                         el_con_celldm_geo(:,:), rot_mat(:,:,:), &
                         el_con_at_geo(:,:,:), at_geo(:,:,:)
INTEGER, ALLOCATABLE :: ibrav_geo(:), el_con_ibrav_geo(:)

INTEGER  :: igeo, iwork, igeom, istep, nstep, nstep_tot, part, &
            iu_out, nwork, ios, ngeom, itemp, nvar, ndata
INTEGER  :: find_free_unit
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=2) :: strain_list(21)
INTEGER :: with_same_lattice(21)
LOGICAL :: flag

REAL(DP), ALLOCATABLE :: ener(:), freee(:), entr(:), cv(:)
!
!   variables for interpolation
!
REAL(DP), ALLOCATABLE :: x(:,:), x1(:), ene(:), ent(:), cvh(:), e0_fit(:)
REAL(DP) :: tot_states
TYPE(poly1) :: p1_e0
TYPE(poly2) :: p2_e0
TYPE(poly3) :: p3_e0
TYPE(poly4) :: p4_e0
TYPE(poly1), ALLOCATABLE :: p1(:,:)
TYPE(poly2), ALLOCATABLE :: p2(:,:)
TYPE(poly3), ALLOCATABLE :: p3(:,:)
TYPE(poly4), ALLOCATABLE :: p4(:,:)

WRITE(6,*) 'entered free_energy_for_ec'
FLUSH(6)
e0=1.0_DP
IF (.NOT.elalgen) &
   CALL errore('free_energy_for_ec', 'Incorrect algorithm', 1)
!
!  We start intepolating with a polynomial the energy, entropy and heat 
!  capacity. The free energy is already interpolated.
!
nvar=crystal_parameters(ibrav_save)
IF (lmurn) nvar=1
ngeom=ngeom_
IF (lgeo_to_file) ngeom=ngeo_file
WRITE(6,*) 'free energy ngeom', ngeom
!
! Allocate the polynomials.
!
IF (poly_degree_ph==1) THEN
   ALLOCATE(p1(ntemp,3))
ELSEIF (poly_degree_ph==2) THEN
   ALLOCATE(p2(ntemp,3))
ELSEIF (poly_degree_ph==3) THEN
   ALLOCATE(p3(ntemp,3))
ELSEIF (poly_degree_ph==4) THEN
   ALLOCATE(p4(ntemp,3))
ENDIF
!
!  Allocate space to write the data to interpolate
!
WRITE(6,*) 'tot_ngeo nvar', tot_ngeo, nvar
FLUSH(6)
ALLOCATE(x(nvar,tot_ngeo))
ALLOCATE(ene(tot_ngeo))
ALLOCATE(ent(tot_ngeo))
ALLOCATE(cvh(tot_ngeo))
ALLOCATE(e0_fit(tot_ngeo))
!
!  here generate the x variable of the function, a vector of dimension nvar
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   IF (lmurn) THEN
      x(1,igeo)=omega_geo_(igeo)
   ELSE
      CALL compress_celldm(celldm_geo_(1,igeo), x(1,igeo), nvar, ibrav_save)
   ENDIF
ENDDO
!
!   The interpolation of the zero point energy.
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   ndata=ndata+1
   e0_fit(ndata)= ph_e0(igeo)
ENDDO
IF (poly_degree_ph == 4) THEN
   CALL init_poly(nvar,p4_e0)
   CALL fit_multi_quartic(ndata, nvar, lsolve, x, e0_fit, p4_e0)
ELSEIF (poly_degree_ph == 3) THEN
   CALL init_poly(nvar,p3_e0)
   CALL fit_multi_cubic(ndata, nvar, lsolve, x, e0_fit, p3_e0)
ELSEIF (poly_degree_ph == 2) THEN
   CALL init_poly(nvar,p2_e0)
   CALL fit_multi_quadratic(ndata, nvar, lsolve, x, e0_fit, p2_e0)
ELSEIF (poly_degree_ph == 1) THEN
   CALL init_poly(nvar,p1_e0)
   CALL fit_multi_linear(ndata, nvar, lsolve, x, e0_fit, p1_e0)
ENDIF
!
!   for each temperature interpolate energy, entropy and C_V.
!   The interpolation of the free energy is already available.
! 
DO itemp=1,ntemp
!
!  collect the computed data at this temperature
!
   ndata=0
   DO igeo=1, tot_ngeo
      IF (no_ph(igeo)) CYCLE
      ndata=ndata+1
      ene(ndata)= ph_ener_eos(itemp,igeo)
      ent(ndata)= ph_entropy_eos(itemp,igeo)
      cvh(ndata)= ph_ce_eos(itemp,igeo)
   ENDDO
!
!  Interpolate the computed data with a polynomial.
!
   IF (poly_degree_ph == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
      CALL init_poly(nvar,p4(itemp,1))
      CALL init_poly(nvar,p4(itemp,2))
      CALL init_poly(nvar,p4(itemp,3))
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, ene, p4(itemp,1))
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, ent, p4(itemp,2))
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, cvh, p4(itemp,3))
   ELSEIF (poly_degree_ph == 3) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p3(itemp,1))
      CALL init_poly(nvar,p3(itemp,2))
      CALL init_poly(nvar,p3(itemp,3))
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, ene, p3(itemp,1))
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, ent, p3(itemp,2))
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, cvh, p3(itemp,3))
   ELSEIF (poly_degree_ph == 2) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p2(itemp,1))
      CALL init_poly(nvar,p2(itemp,2))
      CALL init_poly(nvar,p2(itemp,3))
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, ene, p2(itemp,1))
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, ent, p2(itemp,2))
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, cvh, p2(itemp,3))
   ELSEIF (poly_degree_ph == 1) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p1(itemp,1))
      CALL init_poly(nvar,p1(itemp,2))
      CALL init_poly(nvar,p1(itemp,3))
      CALL fit_multi_linear(ndata, nvar, lsolve, x, ene, p1(itemp,1))
      CALL fit_multi_linear(ndata, nvar, lsolve, x, ent, p1(itemp,2))
      CALL fit_multi_linear(ndata, nvar, lsolve, x, cvh, p1(itemp,3))
   ENDIF
ENDDO  ! on temperatures
!
!  free the interpolate functions
!
DEALLOCATE(x)
DEALLOCATE(ene)
DEALLOCATE(ent)
DEALLOCATE(cvh)
DEALLOCATE(e0_fit)
!
!   Now generate all the strained structures for which the Bravais
!   lattice type does not change and evaluate in these strained 
!   structures the energy, free energy, entropy and heat capacity.
!
CALL set_strain_for_ec(strain_list, nstep, with_same_lattice, laue, &
                        ibrav_save, old_ec, elalgen, elastic_algorithm)
IF (nstep<1.AND.elastic_algorithm=='energy') &
   CALL errore('free_energy_for_ec', 'Incorrect nstep, use energy_std algorithm',1)

IF (nstep<1) &
   CALL errore('free_energy_for_ec', 'Incorrect nstep, check elastic_algorithm',1)

nstep_tot=nstep

! Total number of distorted geometries
nwork = nstep_tot * ngeo_strain * ngeom
!
! Allocate space for describing strain and the strained structures
!
ALLOCATE( epsilon_voigt(6, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
ALLOCATE( epsil_geo(nwork) )
ALLOCATE( at_geo(3, 3, nwork) )
ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )
ALLOCATE( rot_mat(3,3,nwork) )
!
!   Now compute the strained geometries on which we have to interpolate
!   energy, free energy, entropy and heat capacity.
!   Start with the unperturbed geometries
!
IF (ngeom==1) THEN
   IF (.NOT. ALLOCATED(el_con_ibrav_geo)) ALLOCATE(el_con_ibrav_geo(1))
   IF (.NOT. ALLOCATED(el_con_celldm_geo)) ALLOCATE(el_con_celldm_geo(6,1))
   IF (.NOT. ALLOCATED(el_con_at_geo)) ALLOCATE(el_con_at_geo(3,3,1))
   el_con_ibrav_geo(1)=ibrav_save
   el_con_celldm_geo(:,1)=celldm0(:)
   CALL latgen(ibrav_save,celldm0,at_geo(1,1,1),at_geo(1,2,1),at_geo(1,3,1), &
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
      CALL latgen(ibrav_geo_(igeom),celldm_geo_(1,igeom), &
            at_geo(1,1,igeom),at_geo(1,2,igeom), at_geo(1,3,igeom), omega)
      el_con_at_geo(:,:,igeom)=at_geo(:,:,igeom)
   ENDDO
ENDIF

flag=(elastic_algorithm=='energy_std')
epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP

! Build strained geometries 
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

epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO
!
!  And interpolate there the four quantities. At each temperature write 
!  the file.
!
iwork=0
tot_states=3*nat
ALLOCATE(x1(nvar))
ALLOCATE(freee(ntemp))
ALLOCATE(ener(ntemp))
ALLOCATE(entr(ntemp))
ALLOCATE(cv(ntemp))
DO igeom=1, ngeom
   DO istep=1,nstep
      DO igeo=1,ngeo_strain
         iwork=iwork+1
         IF (with_same_lattice(istep)==1) THEN
            IF (elastic_algorithm=='energy_std') THEN
               CALL compute_celldm_geo_from_strain(el_con_at_geo(1,1,igeom), &
                     epsilon_geo(1,1,iwork), celldm_geo(1,iwork), ibrav_save)
               ibrav_geo(iwork)=ibrav_save
            ENDIF
            IF (lmurn) THEN
                x1(1)=celldm_geo(1,iwork)
            ELSE
               CALL compress_celldm(celldm_geo(1,iwork), x1, nvar, &
                                                    ibrav_geo(iwork))
            ENDIF

            DO itemp=1, ntemp
               ! Evaluate temperature-dependent polynomials at this geometry:
               IF (poly_degree_ph==1) THEN
                  ! poly1
                  CALL evaluate_fit_linear  (nvar, x1, freee(itemp), &
                                                                p1t_t(itemp))
                  CALL evaluate_fit_linear  (nvar, x1, ener(itemp), &
                                                                p1(itemp,1))
                  CALL evaluate_fit_linear  (nvar, x1, entr(itemp), &
                                                                p1(itemp,2))
                  CALL evaluate_fit_linear  (nvar, x1, cv(itemp), &
                                                                p1(itemp,3))
               ELSEIF (poly_degree_ph==2) THEN
                  ! poly2
                  CALL evaluate_fit_quadratic (nvar, x1, freee(itemp), &
                                                                p2t_t(itemp))
                  CALL evaluate_fit_quadratic (nvar, x1, ener(itemp), &
                                                                p2(itemp,1))
                  CALL evaluate_fit_quadratic (nvar, x1, entr(itemp), &
                                                                p2(itemp,2))
                  CALL evaluate_fit_quadratic (nvar, x1, cv(itemp), &
                                                                p2(itemp,3))
               ELSEIF (poly_degree_ph==3) THEN
                  !  poly3
                  CALL evaluate_fit_cubic (nvar, x1, freee(itemp), &
                                                               p3t_t(itemp))
                  CALL evaluate_fit_cubic (nvar, x1, ener(itemp), &
                                                               p3(itemp,1))
                  CALL evaluate_fit_cubic (nvar, x1, entr(itemp), &
                                                               p3(itemp,2))
                  CALL evaluate_fit_cubic (nvar, x1, cv(itemp), p3(itemp,3))
               ELSEIF (poly_degree_ph==4) THEN
                  !   poly4
                  CALL evaluate_fit_quartic (nvar, x1, freee(itemp), &
                                                                 p4t_t(itemp))
                  CALL evaluate_fit_quartic (nvar, x1, ener(itemp), &
                                                                 p4(itemp,1))
                  CALL evaluate_fit_quartic (nvar, x1, entr(itemp), &
                                                                 p4(itemp,2))
                  CALL evaluate_fit_quartic (nvar, x1, cv(itemp), &
                                                                 p4(itemp,3))
               ENDIF
            ENDDO
!
!   Interpolation of zero points energy
!
            IF (poly_degree_ph==1) THEN
               CALL evaluate_fit_linear (nvar, x1, e0, p1_e0)
            ELSEIF (poly_degree_ph==2) THEN
               CALL evaluate_fit_quadratic (nvar, x1, e0, p2_e0)
            ELSEIF (poly_degree_ph==3) THEN
               CALL evaluate_fit_cubic (nvar, x1, e0, p3_e0)
            ELSEIF (poly_degree_ph==4) THEN
               CALL evaluate_fit_quartic (nvar, x1, e0, p4_e0)
            ENDIF
            IF (meta_ionode) THEN
               iu_out = find_free_unit()
               filename='energy_dist/'//TRIM(fltherm)//'.g'//&
                                          TRIM(int_to_char(iwork))
               CALL write_thermo_info(e0, tot_states, ntemp, temp,  &
                   ener, freee, entr, cv, 1,filename)

            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO
DEALLOCATE(x1)
DEALLOCATE(freee)
DEALLOCATE(ener)
DEALLOCATE(entr)
DEALLOCATE(cv)

DEALLOCATE( epsilon_voigt )
DEALLOCATE( epsilon_geo )
DEALLOCATE( epsil_geo )
DEALLOCATE( at_geo )
DEALLOCATE( ibrav_geo )
DEALLOCATE( celldm_geo )
DEALLOCATE( rot_mat )

IF (poly_degree_ph==1) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p1(itemp,1))
      CALL clean_poly(p1(itemp,2))
      CALL clean_poly(p1(itemp,3))
   ENDDO
   DEALLOCATE(p1)
   CALL clean_poly(p1_e0)
ELSEIF (poly_degree_ph==2) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p2(itemp,1))
      CALL clean_poly(p2(itemp,2))
      CALL clean_poly(p2(itemp,3))
   ENDDO
   DEALLOCATE(p2)
   CALL clean_poly(p2_e0)
ELSEIF (poly_degree_ph==3) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p3(itemp,1))
      CALL clean_poly(p3(itemp,2))
      CALL clean_poly(p3(itemp,3))
   ENDDO
   DEALLOCATE(p3)
   CALL clean_poly(p3_e0)
ELSEIF (poly_degree_ph==4) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p4(itemp,1))
      CALL clean_poly(p4(itemp,2))
      CALL clean_poly(p4(itemp,3))
   ENDDO
   DEALLOCATE(p4)
   CALL clean_poly(p4_e0)
ENDIF
!


RETURN
END SUBROUTINE free_energy_for_ec
!
!----------------------------------------------------------------------
SUBROUTINE free_energy_for_ec_ph(ngeom_)
!----------------------------------------------------------------------
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : ibrav_geo_=>ibrav_geo, celldm_geo_=>celldm_geo, &
                              omega_geo_=>omega_geo, no_ph, tot_ngeo
USE ions_base,         ONLY : nat
USE control_elastic_constants, ONLY : delta_epsilon, ngeo_strain,          &
                               elastic_algorithm, elalgen, old_ec
USE ph_freq_thermodynamics,    ONLY : phf_ener_eos, phf_entropy_eos,       &
                               phf_ce_eos, phf_e0
USE control_mur,       ONLY : lmurn
USE geometry_file,     ONLY : ngeo_file, celldm_geo_file
USE control_thermo,    ONLY : lgeo_to_file
USE initial_conf,      ONLY : ibrav_save
USE equilibrium_conf,  ONLY : celldm0
USE thermo_sym,        ONLY : laue
USE temperature,       ONLY : ntemp, temp
USE data_files,        ONLY : fltherm


USE elastic_constants, ONLY : set_strain_for_ec
USE strain_mod,        ONLY : set_strain_adv, trans_epsilon

! Polynomial evaluators
USE polynomial,        ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE linear_surfaces,   ONLY : evaluate_fit_linear, fit_multi_linear      
USE quadratic_surfaces,ONLY : evaluate_fit_quadratic, fit_multi_quadratic   
USE cubic_surfaces,    ONLY : evaluate_fit_cubic, fit_multi_cubic        
USE quartic_surfaces,  ONLY : evaluate_fit_quartic, fit_multi_quartic    
USE control_quartic_energy, ONLY : lsolve, poly_degree_ph

USE ph_freq_anharmonic, ONLY : p1tf_t, p2tf_t, p3tf_t, p4tf_t

USE control_mur,       ONLY : lmurn
USE lattices,          ONLY : compress_celldm, crystal_parameters

USE io_global,         ONLY : meta_ionode

IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeom_

REAL(DP) :: epsilon_min, epsil, omega, e0
REAL(DP), ALLOCATABLE :: epsilon_voigt(:,:), epsilon_geo(:,:,:), &
                         epsil_geo(:), celldm_geo(:,:),          &
                         el_con_celldm_geo(:,:), rot_mat(:,:,:), &
                         el_con_at_geo(:,:,:), at_geo(:,:,:)
INTEGER, ALLOCATABLE :: ibrav_geo(:), el_con_ibrav_geo(:)

INTEGER  :: igeo, iwork, igeom, istep, nstep, nstep_tot, part, &
            iu_out, nwork, ios, ngeom, itemp, nvar, ndata
INTEGER  :: find_free_unit
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=2) :: strain_list(21)
INTEGER :: with_same_lattice(21)
LOGICAL :: flag

REAL(DP), ALLOCATABLE :: ener(:), freee(:), entr(:), cv(:)
!
!   variables for interpolation
!
REAL(DP), ALLOCATABLE :: x(:,:), x1(:), ene(:), ent(:), cvh(:), e0_fit(:)
REAL(DP) :: tot_states
TYPE(poly1) :: p1_e0
TYPE(poly2) :: p2_e0
TYPE(poly3) :: p3_e0
TYPE(poly4) :: p4_e0
TYPE(poly1), ALLOCATABLE :: p1(:,:)
TYPE(poly2), ALLOCATABLE :: p2(:,:)
TYPE(poly3), ALLOCATABLE :: p3(:,:)
TYPE(poly4), ALLOCATABLE :: p4(:,:)

IF (.NOT.elalgen) &
   CALL errore('free_energy_for_ec', 'Incorrect algorithm', 1)
!
!  We start intepolating with a polynomial the energy, entropy and heat 
!  capacity. The free energy is already interpolated.
!
nvar=crystal_parameters(ibrav_save)
IF (lmurn) nvar=1
ngeom=ngeom_
IF (lgeo_to_file) ngeom=ngeo_file
!
! Allocate the polynomials.
!
IF (poly_degree_ph==1) THEN
   ALLOCATE(p1(ntemp,3))
ELSEIF (poly_degree_ph==2) THEN
   ALLOCATE(p2(ntemp,3))
ELSEIF (poly_degree_ph==3) THEN
   ALLOCATE(p3(ntemp,3))
ELSEIF (poly_degree_ph==4) THEN
   ALLOCATE(p4(ntemp,3))
ENDIF
!
!  Allocate space to write the data to interpolate
!
ALLOCATE(x(nvar,tot_ngeo))
ALLOCATE(ene(tot_ngeo))
ALLOCATE(ent(tot_ngeo))
ALLOCATE(cvh(tot_ngeo))
ALLOCATE(e0_fit(tot_ngeo))
!
!  here generate the x variable of the function, a vector of dimension nvar
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   IF (lmurn) THEN
      x(1,igeo)=omega_geo_(igeo)
   ELSE
      CALL compress_celldm(celldm_geo_(1,igeo), x(1,igeo), nvar, ibrav_save)
   ENDIF
ENDDO
!
!   The interpolation of the zero point energy.
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   ndata=ndata+1
   e0_fit(ndata)= phf_e0(igeo)
ENDDO
IF (poly_degree_ph == 4) THEN
   CALL init_poly(nvar,p4_e0)
   CALL fit_multi_quartic(ndata, nvar, lsolve, x, e0_fit, p4_e0)
ELSEIF (poly_degree_ph == 3) THEN
   CALL init_poly(nvar,p3_e0)
   CALL fit_multi_cubic(ndata, nvar, lsolve, x, e0_fit, p3_e0)
ELSEIF (poly_degree_ph == 2) THEN
   CALL init_poly(nvar,p2_e0)
   CALL fit_multi_quadratic(ndata, nvar, lsolve, x, e0_fit, p2_e0)
ELSEIF (poly_degree_ph == 1) THEN
   CALL init_poly(nvar,p1_e0)
   CALL fit_multi_linear(ndata, nvar, lsolve, x, e0_fit, p1_e0)
ENDIF
!
!   for each temperature interpolate energy, entropy and C_V.
!   The interpolation of the free energy is already available.
! 
DO itemp=1,ntemp
!
!  collect the computed data at this temperature
!
   WRITE(6,*) 'working on temperature', itemp, ntemp
   FLUSH(6)
   ndata=0
   DO igeo=1, tot_ngeo
      IF (no_ph(igeo)) CYCLE
      ndata=ndata+1
      ene(ndata)= phf_ener_eos(itemp,igeo)
      ent(ndata)= phf_entropy_eos(itemp,igeo)
      cvh(ndata)= phf_ce_eos(itemp,igeo)
   ENDDO
!
!  Interpolate the computed data with a polynomial.
!
   IF (poly_degree_ph == 4) THEN
!
!   compute the coefficients of the quartic polynomial
!
      CALL init_poly(nvar,p4(itemp,1))
      CALL init_poly(nvar,p4(itemp,2))
      CALL init_poly(nvar,p4(itemp,3))
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, ene, p4(itemp,1))
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, ent, p4(itemp,2))
      CALL fit_multi_quartic(ndata, nvar, lsolve, x, cvh, p4(itemp,3))
   ELSEIF (poly_degree_ph == 3) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p3(itemp,1))
      CALL init_poly(nvar,p3(itemp,2))
      CALL init_poly(nvar,p3(itemp,3))
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, ene, p3(itemp,1))
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, ent, p3(itemp,2))
      CALL fit_multi_cubic(ndata, nvar, lsolve, x, cvh, p3(itemp,3))
   ELSEIF (poly_degree_ph == 2) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p2(itemp,1))
      CALL init_poly(nvar,p2(itemp,2))
      CALL init_poly(nvar,p2(itemp,3))
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, ene, p2(itemp,1))
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, ent, p2(itemp,2))
      CALL fit_multi_quadratic(ndata, nvar, lsolve, x, cvh, p2(itemp,3))
   ELSEIF (poly_degree_ph == 1) THEN
!
!   compute the coefficients of the cubic polynomial
!
      CALL init_poly(nvar,p1(itemp,1))
      CALL init_poly(nvar,p1(itemp,2))
      CALL init_poly(nvar,p1(itemp,3))
      CALL fit_multi_linear(ndata, nvar, lsolve, x, ene, p1(itemp,1))
      CALL fit_multi_linear(ndata, nvar, lsolve, x, ent, p1(itemp,2))
      CALL fit_multi_linear(ndata, nvar, lsolve, x, cvh, p1(itemp,3))
   ENDIF
ENDDO  ! on temperatures
!
!  free the interpolate functions
!
DEALLOCATE(x)
DEALLOCATE(ene)
DEALLOCATE(ent)
DEALLOCATE(cvh)
DEALLOCATE(e0_fit)
!
!   Now generate all the strained structures for which the Bravais
!   lattice type does not change and evaluate in these strained 
!   structures the energy, free energy, entropy and heat capacity.
!
CALL set_strain_for_ec(strain_list, nstep, with_same_lattice, laue, &
                        ibrav_save, old_ec, elalgen, elastic_algorithm)
IF (nstep<1.AND.elastic_algorithm=='energy') &
   CALL errore('free_energy_for_ec', 'Incorrect nstep, use energy_std algorithm',1)

IF (nstep<1) &
   CALL errore('free_energy_for_ec', 'Incorrect nstep, check elastic_algorithm',1)

nstep_tot=nstep

! Total number of distorted geometries
nwork = nstep_tot * ngeo_strain * ngeom
!
! Allocate space for describing strain and the strained structures
!
ALLOCATE( epsilon_voigt(6, nwork) )
ALLOCATE( epsilon_geo(3, 3, nwork) )
ALLOCATE( epsil_geo(nwork) )
ALLOCATE( at_geo(3, 3, nwork) )
ALLOCATE( ibrav_geo(nwork) )
ALLOCATE( celldm_geo(6,nwork) )
ALLOCATE( rot_mat(3,3,nwork) )
!
!   Now compute the strained geometries on which we have to interpolate
!   energy, free energy, entropy and heat capacity.
!   Start with the unperturbed geometries
!
IF (ngeom==1) THEN
   IF (.NOT. ALLOCATED(el_con_ibrav_geo)) ALLOCATE(el_con_ibrav_geo(1))
   IF (.NOT. ALLOCATED(el_con_celldm_geo)) ALLOCATE(el_con_celldm_geo(6,1))
   IF (.NOT. ALLOCATED(el_con_at_geo)) ALLOCATE(el_con_at_geo(3,3,1))
   el_con_ibrav_geo(1)=ibrav_save
   el_con_celldm_geo(:,1)=celldm0(:)
   CALL latgen(ibrav_save,celldm0,at_geo(1,1,1),at_geo(1,2,1),at_geo(1,3,1), &
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
      CALL latgen(ibrav_geo_(igeom),celldm_geo_(1,igeom), &
            at_geo(1,1,igeom),at_geo(1,2,igeom), at_geo(1,3,igeom), omega)
      el_con_at_geo(:,:,igeom)=at_geo(:,:,igeom)
   ENDDO
ENDIF

flag=(elastic_algorithm=='energy_std')
epsilon_min= - delta_epsilon * (ngeo_strain - 1 ) / 2.0_DP

! Build strained geometries 
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

epsilon_geo=0.0_DP
DO iwork = 1, nwork
   CALL trans_epsilon(epsilon_voigt(1,iwork), epsilon_geo(1,1,iwork), 1)
ENDDO
!
!  And interpolate there the four quantities. At each temperature write 
!  the file.
!
iwork=0
tot_states=3*nat
ALLOCATE(x1(nvar))
ALLOCATE(freee(ntemp))
ALLOCATE(ener(ntemp))
ALLOCATE(entr(ntemp))
ALLOCATE(cv(ntemp))
DO igeom=1, ngeom
   DO istep=1,nstep
      DO igeo=1,ngeo_strain
         iwork=iwork+1
         IF (with_same_lattice(istep)==1) THEN
            IF (elastic_algorithm=='energy_std') THEN
               CALL compute_celldm_geo_from_strain(el_con_at_geo(1,1,igeom), &
                     epsilon_geo(1,1,iwork), celldm_geo(1,iwork), ibrav_save)
               ibrav_geo(iwork)=ibrav_save
            ENDIF
            IF (lmurn) THEN
                x1(1)=celldm_geo(1,iwork)
            ELSE
               CALL compress_celldm(celldm_geo(1,iwork), x1, nvar, &
                                                    ibrav_geo(iwork))
            ENDIF

            DO itemp=1, ntemp
               ! Evaluate temperature-dependent polynomials at this geometry:
               IF (poly_degree_ph==1) THEN
                  ! poly1
                  CALL evaluate_fit_linear  (nvar, x1, freee(itemp), &
                                                               p1tf_t(itemp))
                  CALL evaluate_fit_linear  (nvar, x1, ener(itemp), &
                                                               p1(itemp,1))
                  CALL evaluate_fit_linear  (nvar, x1, entr(itemp), &
                                                               p1(itemp,2))
                  CALL evaluate_fit_linear  (nvar, x1, cv(itemp), &
                                                               p1(itemp,3))
               ELSEIF (poly_degree_ph==2) THEN
                  ! poly2
                  CALL evaluate_fit_quadratic (nvar, x1, freee(itemp), &
                                                               p2tf_t(itemp))
                  CALL evaluate_fit_quadratic (nvar, x1, ener(itemp), &
                                                                 p2(itemp,1))
                  CALL evaluate_fit_quadratic (nvar, x1, entr(itemp), &
                                                                 p2(itemp,2))
                  CALL evaluate_fit_quadratic (nvar, x1, cv(itemp),  &
                                                                 p2(itemp,3))
               ELSEIF (poly_degree_ph==3) THEN
                  !  poly3
                  CALL evaluate_fit_cubic (nvar, x1, freee(itemp), &
                                                                p3tf_t(itemp))
                  CALL evaluate_fit_cubic (nvar, x1, ener(itemp), p3(itemp,1))
                  CALL evaluate_fit_cubic (nvar, x1, entr(itemp), p3(itemp,2))
                  CALL evaluate_fit_cubic (nvar, x1, cv(itemp), p3(itemp,3))
               ELSEIF (poly_degree_ph==4) THEN
                  !   poly4
                  CALL evaluate_fit_quartic (nvar, x1, freee(itemp), &
                                                                 p4tf_t(itemp))
                  CALL evaluate_fit_quartic (nvar, x1, ener(itemp), &
                                                                 p4(itemp,1))
                  CALL evaluate_fit_quartic (nvar, x1, entr(itemp), &
                                                                 p4(itemp,2))
                  CALL evaluate_fit_quartic (nvar, x1, cv(itemp), p4(itemp,3))
               ENDIF
            ENDDO
!
!   Interpolation of zero points energy
!
            IF (poly_degree_ph==1) THEN
               CALL evaluate_fit_linear (nvar, x1, e0, p1_e0)
            ELSEIF (poly_degree_ph==2) THEN
               CALL evaluate_fit_quadratic (nvar, x1, e0, p2_e0)
            ELSEIF (poly_degree_ph==3) THEN
               CALL evaluate_fit_cubic (nvar, x1, e0, p3_e0)
            ELSEIF (poly_degree_ph==4) THEN
               CALL evaluate_fit_quartic (nvar, x1, e0, p4_e0)
            ENDIF
            IF (meta_ionode) THEN
               iu_out = find_free_unit()
               filename='energy_dist/'//TRIM(fltherm)//'.g'//&
                                          TRIM(int_to_char(iwork))//'_ph'
               CALL write_thermo_info(e0, tot_states, ntemp, temp,  &
                   ener, freee, entr, cv, 1,filename)

            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(x1)
DEALLOCATE(freee)
DEALLOCATE(ener)
DEALLOCATE(entr)
DEALLOCATE(cv)

DEALLOCATE( epsilon_voigt )
DEALLOCATE( epsilon_geo )
DEALLOCATE( epsil_geo )
DEALLOCATE( at_geo )
DEALLOCATE( ibrav_geo )
DEALLOCATE( celldm_geo )
DEALLOCATE( rot_mat )

IF (poly_degree_ph==1) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p1(itemp,1))
      CALL clean_poly(p1(itemp,2))
      CALL clean_poly(p1(itemp,3))
   ENDDO
   DEALLOCATE(p1)
   CALL clean_poly(p1_e0)
ELSEIF (poly_degree_ph==2) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p2(itemp,1))
      CALL clean_poly(p2(itemp,2))
      CALL clean_poly(p2(itemp,3))
   ENDDO
   DEALLOCATE(p2)
   CALL clean_poly(p2_e0)
ELSEIF (poly_degree_ph==3) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p3(itemp,1))
      CALL clean_poly(p3(itemp,2))
      CALL clean_poly(p3(itemp,3))
   ENDDO
   DEALLOCATE(p3)
   CALL clean_poly(p3_e0)
ELSEIF (poly_degree_ph==4) THEN
   DO itemp=1,ntemp
      CALL clean_poly(p4(itemp,1))
      CALL clean_poly(p4(itemp,2))
      CALL clean_poly(p4(itemp,3))
   ENDDO
   DEALLOCATE(p4)
   CALL clean_poly(p4_e0)
ENDIF
!


RETURN
END SUBROUTINE free_energy_for_ec_ph
!
!-----------------------------------------------------------------------------
