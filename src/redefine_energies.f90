!
! Copyright (C) 2024 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE redefine_energies(energy_geo, epsilon_geo, epsil_geo, nwork,  &
                             energy_geo_eff, epsilon_geo_eff, nwork_eff)
!------------------------------------------------------------------------
!
! This routine is used when the elastic constants are calculated by
! moving atoms in some strain types. 
! It receives the list of energies and strain and remove all the energies
! that correspond to the same strain but different atomic positions and
! substitutes them with the energy at the minimum. To obtain it, it fits
! the energies as a function of the atomic position with a second
! order polynomial and finds the minimum.
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : uint_geo
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain, &
                                      stype, nmove, atom_step, min_y, &
                                      epsil_y, lfp, start_geometry_qha, &
                                      last_geometry_qha, ninternal_ec,  &
                                      nint_var_ec, stypec
USE quadratic_surfaces,       ONLY : fit_multi_quadratic, &
                                     find_quadratic_extremum
USE quartic_surfaces,         ONLY : fit_multi_quartic, &
                                     find_quartic_extremum, &
                                     evaluate_fit_quartic
USE control_quartic_energy,  ONLY : lsolve
USE polyfit_mod, ONLY : polyfit, compute_poly
USE polynomial, ONLY : poly2, poly4, clean_poly, init_poly
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
INTEGER, INTENT(OUT)  :: nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(OUT) :: energy_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo(nwork)
REAL(DP), INTENT(IN)  :: epsilon_geo(3,3,nwork)
REAL(DP), INTENT(OUT) :: epsilon_geo_eff(3,3,nwork)

INTEGER :: igeom, istep, igeo, iwork, jwork, imov, ivar
REAL(DP) :: a(3), emin, xmin(1)
REAL(DP), ALLOCATABLE :: ene(:), y(:), yv(:,:), xmin_v(:), ene_v(:)
TYPE(poly2) :: p2
TYPE(poly4) :: p4


ALLOCATE(ene(nmove))
ALLOCATE(y(nmove))

iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index
!
!  Increase the iwork and jwork for the geometries not computed
!
DO igeom=1, start_geometry_qha-1
   DO istep=1, nstep_ec
      DO igeo=1, ngeo_strain
         IF (stype(istep)) THEN
            DO imov=1, nmove
               jwork=jwork+1
            ENDDO
         ELSEIF (stypec(istep)) THEN
            DO imov=1, ninternal_ec(istep)
               jwork=jwork+1
            ENDDO
         ELSE
            jwork=jwork+1
         ENDIF
         iwork=iwork+1
      ENDDO
   ENDDO
ENDDO

DO igeom=start_geometry_qha, last_geometry_qha
   DO istep=1, nstep_ec
      DO igeo=1, ngeo_strain
         IF (stype(istep)) THEN
            DO imov=1, nmove
               jwork=jwork+1
               ene(imov)=energy_geo(jwork)
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
               y(imov)=(imov-(nmove+1.0_DP)/2.0_DP)*atom_step(istep)
            ENDDO
!          
!           Fit the data with a parabola find the minimum,
!           then refit with a quartic and find the minimum starting 
!           from the one of the parabola
!
            CALL fit_multi_quartic(nmove,1,lsolve,y,ene,p4)
            IF (lfp) THEN
               xmin(1)=0.0_DP
               CALL evaluate_fit_quartic(1,xmin,emin,p4)
            ELSE
               CALL fit_multi_quadratic(nmove,1,lsolve,y,ene,p2)
               CALL find_quadratic_extremum(1,xmin,emin,p2)
               CALL find_quartic_extremum(1,xmin,emin,p4)
            ENDIF
            iwork=iwork+1
            energy_geo_eff(iwork)= emin
            epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)   
            min_y(1,igeo,istep,igeom)=xmin(1)
            epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
         ELSEIF (stypec(istep)) THEN
            ALLOCATE(yv(nint_var_ec(istep),ninternal_ec(istep)))
            ALLOCATE(ene_v(ninternal_ec(istep)))
            ALLOCATE(xmin_v(nint_var_ec(istep)))
            CALL init_poly(nint_var_ec(istep),p2)
            CALL init_poly(nint_var_ec(istep),p4)

            DO imov=1, ninternal_ec(istep)
               jwork=jwork+1
               ene_v(imov)=energy_geo(jwork)
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
               DO ivar=1,nint_var_ec(istep)
                  yv(ivar, imov)= uint_geo(ivar, jwork)
               ENDDO
!              WRITE(6,*) yv(ivar, imov), ene(imov)
            ENDDO
            CALL fit_multi_quartic(ninternal_ec(istep),nint_var_ec(istep),&
                                   lsolve,yv,ene_v,p4)
            IF (lfp) THEN
               xmin_v(1:nint_var_ec(istep))=0.0_DP
               CALL evaluate_fit_quartic(nint_var_ec(istep),xmin_v,emin,p4)
            ELSE
               CALL fit_multi_quadratic(ninternal_ec(istep),nint_var_ec(istep),&
                           lsolve,yv,ene_v,p2)
               CALL find_quadratic_extremum(nint_var_ec(istep),xmin_v,emin,p2)
               CALL find_quartic_extremum(nint_var_ec(istep),xmin_v,emin,p4)
            ENDIF
            iwork=iwork+1
            energy_geo_eff(iwork)= emin
            epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)
            min_y(1:nint_var_ec(istep),igeo,istep,igeom)=xmin_v(:)
            epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
            DEALLOCATE(xmin_v)
            DEALLOCATE(ene_v)
            DEALLOCATE(yv)
            CALL clean_poly(p2)
            CALL clean_poly(p4)
         ELSE
            jwork=jwork+1
            iwork=iwork+1
            energy_geo_eff(iwork)=energy_geo(jwork)
            epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!  Increase the iwork and jwork for the geometries not computed
!
DO igeom=last_geometry_qha+1, ngeom
   DO istep=1, nstep_ec
      DO igeo=1, ngeo_strain
         IF (stype(istep)) THEN
            DO imov=1, nmove
               jwork=jwork+1
            ENDDO
         ELSEIF (stypec(istep)) THEN
            DO imov=1, ninternal_ec(istep)
               jwork=jwork+1
            ENDDO
         ELSE
            jwork=jwork+1
         ENDIF
         iwork=iwork+1
      ENDDO
   ENDDO
ENDDO
nwork_eff=iwork

DEALLOCATE(ene)
DEALLOCATE(y)

RETURN
END SUBROUTINE redefine_energies

!--------------------------------------------------------------------------
SUBROUTINE write_min_y()
!--------------------------------------------------------------------------
USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain, &
                                      stype, min_y_t, epsil_y,      &
                                      dyde, all_geometry_done_geo,  &
                                      start_geometry_qha, last_geometry_qha, &
                                      stypec, nint_var_ec
USE data_files,                ONLY : flanhar
USE temperature,               ONLY : temp, ntemp
USE polyfit_mod,               ONLY : polyfit
USE io_global,                 ONLY : meta_ionode

IMPLICIT NONE
INTEGER :: igeom, istep, itemp, igeo, iu_rel, ivar, nvar
INTEGER :: find_free_unit
REAL(DP), ALLOCATABLE ::  x(:), f(:)
REAL(DP) :: a(2)
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename

IF (.NOT.ANY(stype)) RETURN
iu_rel=find_free_unit()
ALLOCATE(x(ngeo_strain))
ALLOCATE(f(ngeo_strain))

DO igeom=start_geometry_qha, last_geometry_qha
   IF (.NOT.all_geometry_done_geo(igeom)) CYCLE
   CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.int_rel', &
                                filename, igeom)
   DO istep=1,nstep_ec
      IF (stype(istep).OR.stypec(istep)) THEN
!
!      compute the derivative of u with respect to epsilon
!      Interpolate with a straight line and find the slope
!
         IF (stype(istep)) THEN
            nvar=1
         ELSEIF(stypec(istep)) THEN
            nvar=nint_var_ec(istep)
         ENDIF
            
         DO itemp=1,ntemp
            DO ivar=1, nvar
               DO igeo=1,ngeo_strain
                  x(igeo)=epsil_y(igeo,istep,igeom)
                  f(igeo)=min_y_t(ivar,igeo,istep,igeom,itemp)
               ENDDO
               CALL polyfit(x,f,ngeo_strain,a,1)
               dyde(ivar,istep,igeom,itemp)=a(2)
            ENDDO
         ENDDO
         IF (meta_ionode) THEN
            DO ivar=1, nvar
               filename=TRIM(filename)//'.'//int_to_char(istep)//'.'&
                                                          //int_to_char(ivar)
               OPEN (UNIT=iu_rel, FILE=TRIM(filename), STATUS='unknown', &
                                                    FORM='formatted')
               WRITE(iu_rel,'("#Temp (K)   dy/de", 15f15.7)') (epsil_y(igeo,&
                                       istep,igeom),igeo=1,ngeo_strain)
               DO itemp=1, ntemp
                  WRITE(iu_rel,'(15f15.7)') temp(itemp),  dyde(ivar,istep,&
                             igeom,itemp), (min_y_t(ivar,igeo,istep, &
                             igeom,itemp),igeo=1,ngeo_strain)
               ENDDO
               CLOSE(UNIT=iu_rel, STATUS='KEEP')
            ENDDO
         ENDIF
      ENDIF
   ENDDO   
ENDDO

DEALLOCATE(x)
DEALLOCATE(f)

RETURN
END SUBROUTINE write_min_y

!--------------------------------------------------------------------------
SUBROUTINE read_dyde(igeom)
!--------------------------------------------------------------------------
!
!  This routine reads the file that contains the derivative of the
!  displacement with respect to strain as a function of temperature.
!  It is used when what='mur_lc_t' if stype is set to .TRUE.
!  in input as when the elastic constants are calculated.
!
USE kinds, ONLY : DP
USE thermo_mod,                ONLY : tot_ngeo
USE control_elastic_constants, ONLY : dyde, stype, all_geometry_done_geo
USE data_files,                ONLY : flanhar
USE temperature,               ONLY : temp, ntemp
USE mp_world,                  ONLY : world_comm
USE mp,                        ONLY : mp_sum
USE io_global,                 ONLY : meta_ionode

IMPLICIT NONE
INTEGER :: igeom, istep, itemp, iu_rel
INTEGER :: find_free_unit
LOGICAL :: exst
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename

IF (.NOT.ANY(stype)) RETURN
iu_rel=find_free_unit()

dyde(:,:,igeom,:)=0.0_DP
IF (meta_ionode) THEN
   CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.int_rel', &
                                filename, igeom)
   DO istep=1,21
      IF (stype(istep)) THEN
         filename=TRIM(filename)//'.'//int_to_char(istep)
         INQUIRE(FILE=TRIM(filename),EXIST=exst)
         IF (exst) THEN
            OPEN (UNIT=iu_rel, FILE=TRIM(filename), STATUS='old', &
                                                 FORM='formatted')
            READ(iu_rel,*)
            DO itemp=1, ntemp
               READ(iu_rel,'(15f15.7)') temp(itemp), dyde(1,istep,igeom,itemp)
            ENDDO
            CLOSE(UNIT=iu_rel, STATUS='KEEP')
         ENDIF    
      ENDIF
   ENDDO   
ENDIF
CALL mp_sum(dyde(:,:,igeom,:),world_comm)

RETURN
END SUBROUTINE read_dyde
!
!------------------------------------------------------------------------
SUBROUTINE redefine_energies_qua_t(energy_geo, epsilon_geo, epsil_geo,  &
                nwork, energy_geo_eff, epsilon_geo_eff, min_y_t,        &
                nwork_eff, igeom)
!------------------------------------------------------------------------
!
! This routine is used when the elastic constants are calculated by
! moving atoms in some strain types. 
! It receives the list of free energies and strain and remove all the 
! free energies that correspond to the same strain but different 
! atomic positions and substitutes them with the free energy at the minimum. 
! To obtain it, it fits the energies as a function of the atomic position 
! with a fourth order polynomial and finds the minimum.
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : uint_geo
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain,      &
                                      stype, nmove, atom_step,           &
                                      epsil_y, min_y, lzsisa, lfp,       &
                                      stypec, nint_var_ec, ninternal_ec, &
                                      max_nint_var
USE control_quartic_energy, ONLY : lsolve
USE polynomial,   ONLY : init_poly, clean_poly, poly2, poly4
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             evaluate_fit_quartic
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
INTEGER, INTENT(IN)   :: igeom
INTEGER, INTENT(OUT)  :: nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(OUT) :: energy_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo(nwork)
REAL(DP), INTENT(IN)  :: epsilon_geo(3,3,nwork)
REAL(DP), INTENT(OUT) :: epsilon_geo_eff(3,3,nwork)
REAL(DP), INTENT(OUT) :: min_y_t(max_nint_var, ngeo_strain, 21)

INTEGER :: istep, igeo, iwork, jwork, imov, nvar, ivar
REAL(DP) :: a(3), emin, xmin(1)
REAL(DP), ALLOCATABLE :: ene(:), y(:), yv(:,:), ene_v(:), xmin_v(:)
TYPE(poly2) :: p2, p2c
TYPE(poly4) :: p4, p4c

nvar=1
CALL init_poly(nvar,p2)
CALL init_poly(nvar,p4)

iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index
DO istep=1, nstep_ec
   DO igeo=1, ngeo_strain
      IF (stype(istep)) THEN
         ALLOCATE(ene(nmove))
         ALLOCATE(y(nmove))
         DO imov=1, nmove
            jwork=jwork+1
            ene(imov)=energy_geo(jwork)
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
            y(imov)=(imov-(nmove+1.0_DP)/2.0_DP)*atom_step(istep)
         ENDDO
!          
!           Fit the data with a quartic and put the minimum free energy
!           or the free energy at the minimum of the energy (lsize) or
!           the free energy at the uniformely strained positions
!
         CALL fit_multi_quartic(nmove,1,lsolve,y,ene,p4)
         IF (lzsisa) THEN
            xmin(1)=min_y(1,igeo,istep,igeom)
            CALL evaluate_fit_quartic(1,xmin,emin,p4)
         ELSEIF (lfp) THEN
            xmin(1)=0.0_DP
            CALL evaluate_fit_quartic(1,xmin(1),emin,p4)
         ELSE
            CALL fit_multi_quadratic(nmove,1,lsolve,y,ene,p2)
            CALL find_quadratic_extremum(1,xmin,emin,p2)
            CALL find_quartic_extremum(1,xmin,emin,p4)
         ENDIF
         iwork=iwork+1
         energy_geo_eff(iwork)= emin
         epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)   
         min_y_t(1,igeo,istep)=xmin(1)
         epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
         DEALLOCATE(y)
         DEALLOCATE(ene)
      ELSEIF (stypec(istep)) THEN
         ALLOCATE(yv(nint_var_ec(istep), ninternal_ec(istep)))
         ALLOCATE(ene_v(ninternal_ec(istep)))
         ALLOCATE(xmin_v(nint_var_ec(istep)))
         CALL init_poly(nint_var_ec(istep),p2c)
         CALL init_poly(nint_var_ec(istep),p4c)
         DO imov=1, ninternal_ec(istep)
            jwork=jwork+1
            ene_v(imov)=energy_geo(jwork)
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
            DO ivar=1, nint_var_ec(istep)
               yv(ivar, imov)= uint_geo(ivar, jwork)
            ENDDO
         ENDDO
!          
!           Fit the data with a quartic and put the minimum free energy
!           or the free energy at the minimum of the energy (lsize) or
!           the free energy at the uniformely strained positions
!
         CALL fit_multi_quartic(ninternal_ec(istep),nint_var_ec(istep),&
                           lsolve,yv,ene_v,p4c)
         IF (lzsisa) THEN
            xmin_v(1:nint_var_ec(istep))=min_y(1:nint_var_ec(istep),igeo,&
                                                               istep,igeom)
            CALL evaluate_fit_quartic(nint_var_ec(istep),xmin_v,emin,p4c)
         ELSEIF (lfp) THEN
            xmin_v=0.0_DP
            CALL evaluate_fit_quartic(nint_var_ec(istep),xmin_v,emin,p4c)
         ELSE
            CALL fit_multi_quadratic(ninternal_ec(istep),nint_var_ec(istep),&
                              lsolve, yv,ene_v,p2c)
            CALL find_quadratic_extremum(nint_var_ec(istep),xmin_v,emin,p2c)
            CALL find_quartic_extremum(nint_var_ec(istep),xmin_v,emin,p4c)
         ENDIF
         iwork=iwork+1
         energy_geo_eff(iwork)= emin
         epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)   
         min_y_t(1:nint_var_ec(istep),igeo,istep)=xmin_v(:)
         epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
         DEALLOCATE(yv)
         DEALLOCATE(ene_v)
         DEALLOCATE(xmin_v)
         CALL clean_poly(p2c)
         CALL clean_poly(p4c)
      ELSE
         jwork=jwork+1
         iwork=iwork+1
         energy_geo_eff(iwork)=energy_geo(jwork)
         epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)
      ENDIF
   ENDDO
ENDDO
nwork_eff=iwork

CALL clean_poly(p2)
CALL clean_poly(p4)

RETURN
END SUBROUTINE redefine_energies_qua_t
!
!------------------------------------------------------------------------
SUBROUTINE redefine_energies_min(energy_geo, celldm_geo, uint_geo, nwork, &
                       energy_geo_eff, celldm_geo_eff, min_u, nwork_eff)
!------------------------------------------------------------------------
!
! This routine is used when the energies are calculated both as a 
! function of internal and external parameters. 
! It receives the list of energies and internal parameters 
! that correspond to the same geometry but different atomic positions and
! substitutes them with the energy at the minimum. To obtain it, it fits
! the energies as a function of the atomic position with a fourth
! order polynomial and finds the minimum.
!
USE kinds, ONLY : DP
USE control_atomic_pos, ONLY : nint_var, ninternal, uint0
USE control_elastic_constants, ONLY : lfp
USE control_quartic_energy, ONLY : lsolve
USE polynomial,   ONLY : init_poly, clean_poly, poly2, poly4
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             evaluate_fit_quartic 
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
INTEGER, INTENT(OUT)  :: nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(IN)  :: celldm_geo(6,nwork)
REAL(DP), INTENT(IN)  :: uint_geo(nint_var,nwork)
REAL(DP), INTENT(OUT) :: energy_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: celldm_geo_eff(6,nwork)
REAL(DP), INTENT(OUT) :: min_u(nint_var,nwork)

INTEGER :: igeom, istep, igeo, iwork, jwork, imov, nvar, nmove
REAL(DP) :: emin, xmin(nint_var)
REAL(DP), ALLOCATABLE :: ene(:), y(:,:)
TYPE(poly2) :: p2
TYPE(poly4) :: p4

nwork_eff=nwork/ninternal
nvar=nint_var
nmove=ninternal
ALLOCATE(ene(nmove))
ALLOCATE(y(nint_var,nmove))

iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index

CALL init_poly(nvar,p2)
CALL init_poly(nvar,p4)
DO igeo=1, nwork_eff
   DO imov=1, nmove
      jwork=jwork+1
      ene(imov)=energy_geo(jwork)
!
!  y is the internal parameter that determine the atomic positions.
!
      y(1:nint_var,imov)=uint_geo(1:nint_var,jwork)
   ENDDO
!          
!           Fit the data with a quartic polynomial and find the minimum,
!           to find the minimum first fit with a parabola, find the minimum
!           and start from that guess.
!
   CALL fit_multi_quartic(nmove,nint_var,lsolve,y,ene,p4)
   IF (lfp) THEN
!
!  frozen phonon uint does not change
!
      xmin(:)=uint0(:)
      CALL evaluate_fit_quartic(nint_var,xmin,emin,p4)
   ELSE
      CALL fit_multi_quadratic(nmove,nint_var,lsolve,y,ene,p2)
      CALL find_quadratic_extremum(nint_var,xmin,emin,p2)
      CALL find_quartic_extremum(nint_var,xmin,emin,p4)
   ENDIF
   iwork=iwork+1
   energy_geo_eff(iwork)= emin
   celldm_geo_eff(:,iwork)=celldm_geo(:,jwork)   
   min_u(:,iwork)=xmin(:)
ENDDO

DEALLOCATE(ene)
DEALLOCATE(y)
CALL clean_poly(p2)
CALL clean_poly(p4)

RETURN
END SUBROUTINE redefine_energies_min
!
!------------------------------------------------------------------------
SUBROUTINE redefine_u_min_t(energy_geo, ph_free_ener, uint_geo, &
                                         nwork, min_u, nwork_eff)
!------------------------------------------------------------------------
!
! This routine is used when the energies are calculated both as a 
! function of internal and external parameters. 
! It receives the list of energies and free energies as a function
! of both internal and external parameters.
! By fixing the external parameters and minimizing the free energy
! with respect to the internal parameters it finds the minimum value
! of the internal parameters at each temperature.
! If lel_free_energy is true, it assumes to have also the electronic
! contribution and adds it to the free energy before finding the minimum.
!
! The output of this routine is, for each temperature, the min_u the
! internal parameters that for each external parameter minimizes the 
! free energy.
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : no_ph
USE control_atomic_pos, ONLY : nint_var, ninternal, uint0
USE control_elastic_constants, ONLY : lfp
USE control_eldos, ONLY : lel_free_energy
USE control_quartic_energy, ONLY : lsolve
USE el_thermodynamics, ONLY : el_free_ener
USE temperature, ONLY : ntemp
USE polynomial,   ONLY : init_poly, clean_poly, poly2, poly4
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             evaluate_fit_quartic 
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork, nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(IN)  :: ph_free_ener(ntemp, nwork)
REAL(DP), INTENT(IN)  :: uint_geo(nint_var,nwork)
REAL(DP), INTENT(OUT) :: min_u(nint_var,ntemp,nwork_eff)

INTEGER :: igeom, istep, igeo, iwork, jwork, imov, nvar, nmove, nwork_eff_, &
           itemp
REAL(DP) :: emin, xmin(nint_var)
REAL(DP), ALLOCATABLE :: ene(:), y(:,:)
TYPE(poly2) :: p2
TYPE(poly4) :: p4

nwork_eff_=nwork/ninternal
IF (nwork_eff_ /= nwork_eff) CALL errore('redefine_free_ener_min', &
                               'problem with nwork_eff', 1)
nvar=nint_var
nmove=ninternal
ALLOCATE(ene(nmove))
ALLOCATE(y(nint_var,nmove))

CALL init_poly(nvar,p2)
CALL init_poly(nvar,p4)
DO itemp=1, ntemp
   iwork=0     ! run on current energy index
   jwork=0     ! run on previous energy index
   DO igeo=1, nwork_eff
      DO imov=1, nmove
         jwork=jwork+1
         ene(imov)=energy_geo(jwork)+ph_free_ener(itemp,jwork)
         IF (lel_free_energy) ene(imov)=ene(imov)+el_free_ener(itemp,jwork)
!
!  y is the internal parameter that determine the atomic positions.
!
         y(1:nint_var,imov)=uint_geo(1:nint_var,jwork)
      ENDDO
!          
!     Fit the data with a quartic polynomial and find the minimum,
!     to find the minimum first fit with a parabola, find the minimum
!     and start from that guess.
!
      CALL fit_multi_quartic(nmove,nint_var,lsolve,y,ene,p4)
      IF (lfp) THEN
!
!  frozen phonon uint does not change
!
         xmin(:)=uint0(:)
         CALL evaluate_fit_quartic(nint_var,xmin,emin,p4)
      ELSE
         IF (itemp==1) THEN
            CALL fit_multi_quadratic(nmove,nint_var,lsolve,y,ene,p2)
            CALL find_quadratic_extremum(nint_var,xmin,emin,p2)
         ENDIF
         CALL find_quartic_extremum(nint_var,xmin,emin,p4)
      ENDIF
      iwork=iwork+1
      min_u(1:nint_var,itemp,iwork)=xmin(1:nint_var)
   ENDDO
ENDDO

DEALLOCATE(ene)
DEALLOCATE(y)
CALL clean_poly(p2)
CALL clean_poly(p4)

RETURN
END SUBROUTINE redefine_u_min_t
!
!------------------------------------------------------------------------
SUBROUTINE redefine_thermo_t(ph_free_ener, ph_e0, ph_ener, ph_entropy, &
                             ph_ce, uint_geo, nwork, ph_free_ener_eos, &
                             ph_e0_eos, ph_ener_eos, ph_entropy_eos, &
                             ph_ce_eos, min_u, nwork_eff)
!------------------------------------------------------------------------
!
!  This routine interpolates the vibrational free energy, energy, 
!  entropy and heat capacity with a polynomial of the internal parameters
!  and computes them at the minimum of the internal parameters.
!  The minimum of the internal parameters for each external geometry
!  is already contained in min_u.
!
USE kinds, ONLY : DP
USE control_atomic_pos, ONLY : nint_var, ninternal, uint0
USE control_elastic_constants, ONLY : lfp
USE control_quartic_energy, ONLY : lsolve
USE temperature, ONLY : ntemp
USE polynomial,   ONLY : init_poly, clean_poly, poly2, poly4
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             evaluate_fit_quartic 
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork, nwork_eff
REAL(DP), INTENT(IN)  :: ph_free_ener(ntemp, nwork)
REAL(DP), INTENT(IN)  :: ph_e0(nwork)
REAL(DP), INTENT(IN)  :: ph_ener(ntemp, nwork)
REAL(DP), INTENT(IN)  :: ph_entropy(ntemp, nwork)
REAL(DP), INTENT(IN)  :: ph_ce(ntemp, nwork)
REAL(DP), INTENT(IN)  :: uint_geo(nint_var,nwork)
REAL(DP), INTENT(OUT) :: ph_free_ener_eos(ntemp, nwork_eff)
REAL(DP), INTENT(OUT) :: ph_e0_eos(nwork_eff)
REAL(DP), INTENT(OUT) :: ph_ener_eos(ntemp, nwork_eff)
REAL(DP), INTENT(OUT) :: ph_entropy_eos(ntemp, nwork_eff)
REAL(DP), INTENT(OUT) :: ph_ce_eos(ntemp, nwork_eff)
REAL(DP), INTENT(OUT) :: min_u(nint_var,ntemp,nwork_eff)

INTEGER :: igeom, istep, igeo, iwork, jwork, imov, nvar, nmove, nwork_eff_, &
           itemp
REAL(DP) :: emin, xmin(nint_var)
REAL(DP), ALLOCATABLE :: ener(:), free_ener(:), entropy(:), ce(:), y(:,:), &
                         e0(:)
TYPE(poly4) :: p4

nwork_eff_=nwork/ninternal
IF (nwork_eff_ /= nwork_eff) CALL errore('redefine_free_ener_min', &
                               'problem with nwork_eff', 1)
nvar=nint_var
nmove=ninternal
ALLOCATE(free_ener(nmove))
ALLOCATE(ener(nmove))
ALLOCATE(entropy(nmove))
ALLOCATE(ce(nmove))
ALLOCATE(e0(nmove))
ALLOCATE(y(nint_var,nmove))

CALL init_poly(nvar,p4)
DO itemp=1, ntemp
   iwork=0     ! run on current energy index
   jwork=0     ! run on previous energy index
   DO igeo=1, nwork_eff
      DO imov=1, nmove
         jwork=jwork+1
         free_ener(imov)=ph_free_ener(itemp,jwork)
         ener(imov)=ph_ener(itemp,jwork)
         entropy(imov)=ph_entropy(itemp,jwork)
         ce(imov)=ph_ce(itemp,jwork)
         IF (itemp==1) e0(imov)=ph_e0(jwork)
!
!  y is the internal parameter that determine the atomic positions.
!
         y(1:nint_var,imov)=uint_geo(1:nint_var,jwork)
      ENDDO
!          
!     Fit the data with a quartic polynomial and find the value of
!     the polynomial at the minimum internal parameter.
!
      iwork=iwork+1
      IF (itemp==1) THEN
         CALL fit_multi_quartic(nmove,nint_var,lsolve,y,e0,p4)
         CALL evaluate_fit_quartic(nint_var,min_u(1,itemp,iwork),emin,p4)
         ph_e0_eos(iwork)= emin
      ENDIF
      CALL fit_multi_quartic(nmove,nint_var,lsolve,y,free_ener,p4)
      CALL evaluate_fit_quartic(nint_var,min_u(1,itemp,iwork),emin,p4)
      ph_free_ener_eos(itemp,iwork)= emin
      CALL fit_multi_quartic(nmove,nint_var,lsolve,y,ener,p4)
      CALL evaluate_fit_quartic(nint_var,min_u(1,itemp,iwork),emin,p4)
      ph_ener_eos(itemp,iwork)= emin
      CALL fit_multi_quartic(nmove,nint_var,lsolve,y,entropy,p4)
      CALL evaluate_fit_quartic(nint_var,min_u(1,itemp,iwork),emin,p4)
      ph_entropy_eos(itemp,iwork)= emin
      CALL fit_multi_quartic(nmove,nint_var,lsolve,y,ce,p4)
      CALL evaluate_fit_quartic(nint_var,min_u(1,itemp,iwork),emin,p4)
      ph_ce_eos(itemp,iwork)= emin
   ENDDO
ENDDO

DEALLOCATE(free_ener)
DEALLOCATE(ener)
DEALLOCATE(entropy)
DEALLOCATE(e0)
DEALLOCATE(ce)
DEALLOCATE(y)
CALL clean_poly(p4)

RETURN
END SUBROUTINE redefine_thermo_t

!------------------------------------------------------------------------
SUBROUTINE redefine_polar(polar_strain, tot_b_phase, nwork, &
               polar_strain_eff, tot_b_phase_eff, epsil_geo, epsil_geo_eff)
!------------------------------------------------------------------------
!
! This routine is used when the piezoelectric tensor is calculated by
! moving atoms in some strain types. 
! It receives the list of polarization and berry phase and remove 
! all the energies that correspond to the same strain but different 
! atomic positions and substitutes them with the energy at the minimum. 
! To obtain it, it fits the energies as a function of the atomic 
! position with a second order polynomial and finds the minimum.
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : uint_geo
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain,     &
                                      stype, nmove, atom_step, min_y,   &
                                      epsil_y, lfp, start_geometry_qha, &
                                      last_geometry_qha, ninternal_ec,  &
                                      nint_var_ec, stypec
USE quadratic_surfaces,       ONLY : fit_multi_quadratic, &
                                     find_quadratic_extremum
USE quartic_surfaces,         ONLY : fit_multi_quartic, &
                                     find_quartic_extremum, &
                                     evaluate_fit_quartic
USE control_quartic_energy,  ONLY : lsolve
USE polynomial, ONLY : poly4, clean_poly, init_poly
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
REAL(DP), INTENT(IN)  :: polar_strain(3,nwork)
REAL(DP), INTENT(OUT) :: polar_strain_eff(3,nwork)
REAL(DP), INTENT(IN)  :: tot_b_phase(3,nwork)
REAL(DP), INTENT(OUT) :: tot_b_phase_eff(3,nwork)
REAL(DP), INTENT(IN)  :: epsil_geo(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo_eff(nwork)

INTEGER :: igeom, istep, igeo, iwork, jwork, imov, ivar, ipol
REAL(DP) :: pol(3), bph(3)
REAL(DP), ALLOCATABLE :: y(:), yv(:,:), pstrain(:,:), bphase(:,:)
TYPE(poly4) :: p4(3), p4_b(3)



iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index
!
!  Increase the iwork and jwork for the geometries not computed
!
DO igeom=1, start_geometry_qha-1
   DO istep=1, nstep_ec
      DO igeo=1, ngeo_strain
         IF (stype(istep)) THEN
            DO imov=1, nmove
               jwork=jwork+1
            ENDDO
         ELSEIF (stypec(istep)) THEN
            DO imov=1, ninternal_ec(istep)
               jwork=jwork+1
            ENDDO
         ELSE
            jwork=jwork+1
         ENDIF
         iwork=iwork+1
      ENDDO
   ENDDO
ENDDO

DO igeom=start_geometry_qha, last_geometry_qha
   DO istep=1, nstep_ec
      DO igeo=1, ngeo_strain
         IF (stype(istep)) THEN
            ALLOCATE(pstrain(nmove,3))
            ALLOCATE(bphase(nmove,3))
            ALLOCATE(y(nmove))
            DO ipol=1,3
               CALL init_poly(1,p4(ipol))
               CALL init_poly(1,p4_b(ipol))
            ENDDO
            DO imov=1, nmove
               jwork=jwork+1
               DO ipol=1,3 
                  pstrain(imov,ipol)=polar_strain(ipol,jwork)
                  bphase(imov,ipol)=tot_b_phase(ipol,jwork)
               ENDDO
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
               y(imov)=(imov-(nmove+1.0_DP)/2.0_DP)*atom_step(istep)
            ENDDO
!          
!           Fit the data with a quartic and evaluate them at the minimum
!           of the internal parameter that should be already available.
!
            DO ipol=1,3
               CALL fit_multi_quartic(nmove,1,lsolve,y,pstrain(1,ipol),p4(ipol))
               CALL fit_multi_quartic(nmove,1,lsolve,y,bphase(1,ipol),&
                                                                     p4_b(ipol))
               CALL evaluate_fit_quartic(1,min_y(1,igeo,istep,igeom),pol(ipol),&
                                                       p4(ipol))
               CALL evaluate_fit_quartic(1,min_y(1,igeo,istep,igeom),bph(ipol),&
                                                       p4_b(ipol))
            ENDDO
            iwork=iwork+1
            polar_strain_eff(iwork,:)= pol(:)
            tot_b_phase_eff(iwork,:)= bph(:)
            epsil_geo_eff(iwork)=epsil_geo(jwork)
            DO ipol=1,3
               CALL clean_poly(p4(ipol))
               CALL clean_poly(p4_b(ipol))
            ENDDO
            DEALLOCATE(pstrain)
            DEALLOCATE(bphase)
            DEALLOCATE(y)
         ELSEIF (stypec(istep)) THEN
            ALLOCATE(yv(nint_var_ec(istep),ninternal_ec(istep)))
            ALLOCATE(pstrain(ninternal_ec(istep),3))
            ALLOCATE(bphase(ninternal_ec(istep),3))
            DO ipol=1,3
               CALL init_poly(nint_var_ec(istep),p4(ipol))
               CALL init_poly(nint_var_ec(istep),p4_b(ipol))
            ENDDO
            DO imov=1, ninternal_ec(istep)
               jwork=jwork+1
               DO ipol=1,3 
                  pstrain(imov,ipol)=polar_strain(ipol,jwork)
                  bphase(imov,ipol)=tot_b_phase(ipol,jwork)
               ENDDO
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
               DO ivar=1,nint_var_ec(istep)
                  yv(ivar, imov)= uint_geo(ivar, jwork)
               ENDDO
!              WRITE(6,*) yv(ivar, imov), ene(imov)
            ENDDO
!
!   fit the polarization and the berry phase with a quartic and interpolate
!   at the internal parameter at the minimum energy (already available)
!
            DO ipol=1,3
               CALL fit_multi_quartic(ninternal_ec(istep),nint_var_ec(istep),&
                                   lsolve,yv,pstrain(1,ipol),p4(ipol))
               CALL fit_multi_quartic(ninternal_ec(istep),nint_var_ec(istep),&
                                   lsolve,yv,bphase(1,ipol),p4_b(ipol))
               CALL evaluate_fit_quartic(nint_var_ec(istep),&
                    min_y(1:nint_var_ec(istep),igeo,istep,igeom),pol(ipol),&
                                                               p4(ipol))
               CALL evaluate_fit_quartic(nint_var_ec(istep),&
                    min_y(1:nint_var_ec(istep),igeo,istep,igeom),bph(ipol),&
                                                              p4_b(ipol))
            ENDDO
            iwork=iwork+1
            polar_strain_eff(:,iwork)=pol(:)
            tot_b_phase_eff(:,iwork)= bph(:)
            epsil_geo_eff(iwork)=epsil_geo(jwork)
            DEALLOCATE(pstrain)
            DEALLOCATE(bphase)
            DEALLOCATE(yv)
            DO ipol=1,3
               CALL clean_poly(p4(ipol))
               CALL clean_poly(p4_b(ipol))
            ENDDO
         ELSE
            jwork=jwork+1
            iwork=iwork+1
            polar_strain_eff(:,iwork)=polar_strain(:,jwork)
            tot_b_phase_eff(:,iwork)=tot_b_phase(:,jwork)
            epsil_geo_eff(iwork)=epsil_geo(jwork)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!  Increase the iwork and jwork for the geometries not computed
!
DO igeom=last_geometry_qha+1, ngeom
   DO istep=1, nstep_ec
      DO igeo=1, ngeo_strain
         IF (stype(istep)) THEN
            DO imov=1, nmove
               jwork=jwork+1
            ENDDO
         ELSEIF (stypec(istep)) THEN
            DO imov=1, ninternal_ec(istep)
               jwork=jwork+1
            ENDDO
         ELSE
            jwork=jwork+1
         ENDIF
         iwork=iwork+1
      ENDDO
   ENDDO
ENDDO


RETURN
END SUBROUTINE redefine_polar

!------------------------------------------------------------------------
SUBROUTINE redefine_polar_qha_t(polar_strain, tot_b_phase, nwork, &
               polar_strain_eff, tot_b_phase_eff, epsil_geo,      &
               epsil_geo_eff, min_y_t)
!------------------------------------------------------------------------
!
! This routine is used when the piezoelectric tensor is calculated by
! moving atoms in some strain types. 
! It receives the list of polarization and berry phase and remove 
! all the energies that correspond to the same strain but different 
! atomic positions and substitutes them with the energy at the minimum. 
! To obtain it, it fits the energies as a function of the atomic 
! position with a second order polynomial and finds the minimum.
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : uint_geo
USE control_elastic_constants, ONLY : nstep_ec, ngeo_strain,            &
                                      stype, nmove, atom_step,          &
                                      epsil_y, lfp, start_geometry_qha, &
                                      last_geometry_qha, ninternal_ec,  &
                                      nint_var_ec, stypec, max_nint_var
USE quadratic_surfaces,       ONLY : fit_multi_quadratic, &
                                     find_quadratic_extremum
USE quartic_surfaces,         ONLY : fit_multi_quartic, &
                                     find_quartic_extremum, &
                                     evaluate_fit_quartic
USE control_quartic_energy,  ONLY : lsolve
USE polynomial, ONLY : poly4, clean_poly, init_poly
IMPLICIT NONE

INTEGER,  INTENT(IN)  :: nwork
REAL(DP), INTENT(IN)  :: polar_strain(3,nwork)
REAL(DP), INTENT(OUT) :: polar_strain_eff(3,nwork)
REAL(DP), INTENT(IN)  :: tot_b_phase(3,nwork)
REAL(DP), INTENT(OUT) :: tot_b_phase_eff(3,nwork)
REAL(DP), INTENT(IN)  :: epsil_geo(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: min_y_t(max_nint_var, ngeo_strain, 21)

INTEGER :: istep, igeo, iwork, jwork, imov, ivar, ipol
REAL(DP) :: pol(3), bph(3)
REAL(DP), ALLOCATABLE :: y(:), yv(:,:), pstrain(:,:), bphase(:,:)
TYPE(poly4) :: p4(3), p4_b(3)


iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index

DO istep=1, nstep_ec
   DO igeo=1, ngeo_strain
      IF (stype(istep)) THEN
         ALLOCATE(pstrain(nmove,3))
         ALLOCATE(bphase(nmove,3))
         ALLOCATE(y(nmove))
         DO ipol=1,3
            CALL init_poly(1,p4(ipol))
            CALL init_poly(1,p4_b(ipol))
         ENDDO
         DO imov=1, nmove
            jwork=jwork+1
            DO ipol=1,3 
               pstrain(imov,ipol)=polar_strain(ipol,jwork)
               bphase(imov,ipol)=tot_b_phase(ipol,jwork)
            ENDDO
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
            y(imov)=(imov-(nmove+1.0_DP)/2.0_DP)*atom_step(istep)
         ENDDO
!          
!           Fit the data with a quartic and evaluate them at the minimum
!           of the internal parameter that should be already available.
!
         DO ipol=1,3
            CALL fit_multi_quartic(nmove,1,lsolve,y,pstrain(1,ipol),p4(ipol))
            CALL fit_multi_quartic(nmove,1,lsolve,y,bphase(1,ipol), p4_b(ipol))
            CALL evaluate_fit_quartic(1,min_y_t(1,igeo,istep), &
                          pol(ipol),p4(ipol))
            CALL evaluate_fit_quartic(1,min_y_t(1,igeo,istep), &
                                      bph(ipol), p4_b(ipol))
         ENDDO
         iwork=iwork+1
         polar_strain_eff(iwork,:)= pol(:)
         tot_b_phase_eff(iwork,:)= bph(:)
         epsil_geo_eff(iwork)=epsil_geo(jwork)
         DO ipol=1,3
            CALL clean_poly(p4(ipol))
            CALL clean_poly(p4_b(ipol))
         ENDDO
         DEALLOCATE(pstrain)
         DEALLOCATE(bphase)
         DEALLOCATE(y)
      ELSEIF (stypec(istep)) THEN
         ALLOCATE(yv(nint_var_ec(istep),ninternal_ec(istep)))
         ALLOCATE(pstrain(ninternal_ec(istep),3))
         ALLOCATE(bphase(ninternal_ec(istep),3))
         DO ipol=1,3
            CALL init_poly(nint_var_ec(istep),p4(ipol))
            CALL init_poly(nint_var_ec(istep),p4_b(ipol))
         ENDDO
         DO imov=1, ninternal_ec(istep)
            jwork=jwork+1
            DO ipol=1,3 
               pstrain(imov,ipol)=polar_strain(ipol,jwork)
               bphase(imov,ipol)=tot_b_phase(ipol,jwork)
            ENDDO
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
            DO ivar=1,nint_var_ec(istep)
               yv(ivar, imov)= uint_geo(ivar, jwork)
            ENDDO
!            WRITE(6,*) yv(ivar, imov), ene(imov)
         ENDDO
!
!   fit the polarization and the berry phase with a quartic and interpolate
!   at the internal parameter at the minimum energy (already available)
!
         DO ipol=1,3
            CALL fit_multi_quartic(ninternal_ec(istep),nint_var_ec(istep),&
                                lsolve,yv,pstrain(1,ipol),p4(ipol))
            CALL fit_multi_quartic(ninternal_ec(istep),nint_var_ec(istep),&
                                lsolve,yv,bphase(1,ipol),p4_b(ipol))
            CALL evaluate_fit_quartic(nint_var_ec(istep),                 &
                 min_y_t(1:nint_var_ec(istep),igeo,istep), pol(ipol), p4(ipol))
            CALL evaluate_fit_quartic(nint_var_ec(istep),                 &
                 min_y_t(1:nint_var_ec(istep),igeo,istep), bph(ipol), &
                                                                    p4_b(ipol))
         ENDDO
         iwork=iwork+1
         polar_strain_eff(:,iwork)=pol(:)
         tot_b_phase_eff(:,iwork)= bph(:)
         epsil_geo_eff(iwork)=epsil_geo(jwork)
         DEALLOCATE(pstrain)
         DEALLOCATE(bphase)
         DEALLOCATE(yv)
         DO ipol=1,3
            CALL clean_poly(p4(ipol))
            CALL clean_poly(p4_b(ipol))
         ENDDO
      ELSE
         jwork=jwork+1
         iwork=iwork+1
         polar_strain_eff(:,iwork)=polar_strain(:,jwork)
         tot_b_phase_eff(:,iwork)=tot_b_phase(:,jwork)
         epsil_geo_eff(iwork)=epsil_geo(jwork)
      ENDIF
   ENDDO
ENDDO

RETURN
END SUBROUTINE redefine_polar_qha_t
