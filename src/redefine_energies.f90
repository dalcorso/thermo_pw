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
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain, &
                                      stype, nmove, atom_step, min_y, &
                                      epsil_y, lfp
USE polyfit_mod, ONLY : polyfit, compute_poly
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
INTEGER, INTENT(OUT)  :: nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(OUT) :: energy_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo(nwork)
REAL(DP), INTENT(IN)  :: epsilon_geo(3,3,nwork)
REAL(DP), INTENT(OUT) :: epsilon_geo_eff(3,3,nwork)

INTEGER :: igeom, istep, igeo, iwork, jwork, imov
REAL(DP) :: a(3), emin, xmin
REAL(DP), ALLOCATABLE :: ene(:), y(:)

ALLOCATE(ene(nmove))
ALLOCATE(y(nmove))

iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index
DO igeom=1, ngeom  
   WRITE(6,*) 'compressing geometry', igeom
   DO istep=1, nstep_ec
      DO igeo=1, ngeo_strain
         IF (stype(istep)) THEN
            DO imov=1, nmove
               jwork=jwork+1
               WRITE(6,*) 'imov, jwork', imov, jwork
               ene(imov)=energy_geo(jwork)
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
               y(imov)=(imov-(nmove+1.0_DP)/2.0_DP)*atom_step(istep)
!               WRITE(6,*) y(imov), ene(imov)
            ENDDO
!          
!           Fit the data with a parabola and put the minimum energy
!
            CALL polyfit(y,ene,nmove,a,2)
!
!   In the frozen phonon case take xmin=0.0_DP
!
            IF (lfp) THEN
               xmin=0.0_DP
            ELSE
               xmin= - a(2) / 2.0_DP / a(3)
            ENDIF
            CALL compute_poly(xmin, 2, a, emin)
            iwork=iwork+1
            WRITE(6,'("igeom, istep, igeo, iwork, jwork",5i7)') igeom, istep, igeo, iwork, jwork
            energy_geo_eff(iwork)= emin
            epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)   
            min_y(igeo,istep,igeom)=xmin
            epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
         ELSE
            jwork=jwork+1
            iwork=iwork+1
            WRITE(6,'("igeom, istep, igeo, iwork, jwork",5i7)') igeom, istep, igeo, iwork, jwork
            energy_geo_eff(iwork)=energy_geo(jwork)
            epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)
         ENDIF
      ENDDO
   ENDDO
ENDDO
nwork_eff=iwork

DEALLOCATE(ene)
DEALLOCATE(y)

RETURN
END SUBROUTINE redefine_energies

!------------------------------------------------------------------------
SUBROUTINE redefine_energies_t(energy_geo, epsilon_geo, epsil_geo, nwork,  &
                energy_geo_eff, epsilon_geo_eff, nwork_eff, igeom, itemp)
!------------------------------------------------------------------------
!
! This routine is used when the elastic constants are calculated by
! moving atoms in some strain types. 
! It receives the list of free energies and strain and remove all the 
! free energies that correspond to the same strain but different 
! atomic positions and substitutes them with the free energy at the minimum. 
! To obtain it, it fits the energies as a function of the atomic position 
! with a second order polynomial and finds the minimum.
!
USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain, &
                                      stype, nmove, atom_step, min_y_t, &
                                      epsil_y, min_y, lzsisa, lfp
USE polyfit_mod, ONLY : polyfit, compute_poly
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
INTEGER, INTENT(IN)   :: igeom
INTEGER, INTENT(IN)   :: itemp
INTEGER, INTENT(OUT)  :: nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(OUT) :: energy_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo(nwork)
REAL(DP), INTENT(IN)  :: epsilon_geo(3,3,nwork)
REAL(DP), INTENT(OUT) :: epsilon_geo_eff(3,3,nwork)

INTEGER :: istep, igeo, iwork, jwork, imov
REAL(DP) :: a(3), emin, xmin
REAL(DP), ALLOCATABLE :: ene(:), y(:)

ALLOCATE(ene(nmove))
ALLOCATE(y(nmove))

iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index
WRITE(6,*) 'compressing geometry', igeom
DO istep=1, nstep_ec
   DO igeo=1, ngeo_strain
      IF (stype(istep)) THEN
         DO imov=1, nmove
            jwork=jwork+1
            WRITE(6,*) 'imov, jwork', imov, jwork
            ene(imov)=energy_geo(jwork)
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
            y(imov)=(imov-(nmove+1.0_DP)/2.0_DP)*atom_step(istep)
!               WRITE(6,*) y(imov), ene(imov)
         ENDDO
!          
!           Fit the data with a parabola and put the minimum energy
!
         
         CALL polyfit(y,ene,nmove,a,2)
         IF (lzsisa) THEN
            xmin=min_y(igeo,istep,igeom)
         ELSEIF (lfp) THEN
            xmin=0.0_DP
         ELSE
            xmin= - a(2) / 2.0_DP / a(3)
         ENDIF
         CALL compute_poly(xmin, 2, a, emin)
         iwork=iwork+1
         energy_geo_eff(iwork)= emin
         epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)   
         min_y_t(igeo,istep,igeom,itemp)=xmin
         epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
      ELSE
         jwork=jwork+1
         iwork=iwork+1
         energy_geo_eff(iwork)=energy_geo(jwork)
         epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)
      ENDIF
   ENDDO
ENDDO
nwork_eff=iwork

DEALLOCATE(ene)
DEALLOCATE(y)

RETURN
END SUBROUTINE redefine_energies_t

!--------------------------------------------------------------------------
SUBROUTINE write_min_y()
!--------------------------------------------------------------------------
USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain, &
                                      stype, min_y_t, epsil_y,      &
                                      dyde, all_geometry_done_geo
USE data_files,                ONLY : flanhar
USE temperature,               ONLY : temp, ntemp
USE polyfit_mod,               ONLY : polyfit

IMPLICIT NONE
INTEGER :: igeom, istep, itemp, igeo, iu_rel
INTEGER :: find_free_unit
REAL(DP), ALLOCATABLE ::  x(:), f(:)
REAL(DP) :: a(2)
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename

IF (.NOT.ANY(stype)) RETURN
iu_rel=find_free_unit()
ALLOCATE(x(ngeo_strain))
ALLOCATE(f(ngeo_strain))

DO igeom=1, ngeom
   IF (.NOT.all_geometry_done_geo(igeom)) CYCLE
   CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.int_rel', &
                                filename, igeom)
   DO istep=1,nstep_ec
      IF (stype(istep)) THEN
!
!      compute the derivative of u with respect to epsilon
!      Interpolate with a straight line and find the slope
!
         DO itemp=1,ntemp
            DO igeo=1,ngeo_strain
               x(igeo)=epsil_y(igeo,istep,igeom)
               f(igeo)=min_y_t(igeo,istep,igeom,itemp)
            ENDDO
            CALL polyfit(x,f,ngeo_strain,a,1)
            dyde(istep,igeom,itemp)=a(2)
         ENDDO

         filename=TRIM(filename)//'.'//int_to_char(istep)
         OPEN (UNIT=iu_rel, FILE=TRIM(filename), STATUS='unknown', &
                                                 FORM='formatted')
         WRITE(iu_rel,'("#Temp (K)   dy/de", 15f15.7)') (epsil_y(igeo,&
                                    istep,igeom),igeo=1,ngeo_strain)
         DO itemp=1, ntemp
            WRITE(iu_rel,'(15f15.7)') temp(itemp),  dyde(istep,igeom,itemp), &
                         (min_y_t(igeo,istep,igeom,itemp),igeo=1,ngeo_strain)
         ENDDO
         CLOSE(UNIT=iu_rel, STATUS='KEEP')
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

dyde(:,igeom,:)=0.0_DP
IF (meta_ionode) THEN
   CALL add_geometry_number('anhar_files/', TRIM(flanhar)//'.int_rel', &
                                filename, igeom)
   DO istep=1,21
      WRITE(6,*) 'check for istep', istep, igeom, TRIM(filename)
      IF (stype(istep)) THEN
         filename=TRIM(filename)//'.'//int_to_char(istep)
         INQUIRE(FILE=TRIM(filename),EXIST=exst)
         WRITE(6,*) 'check for istep ', istep, TRIM(filename), exst
         IF (exst) THEN
            OPEN (UNIT=iu_rel, FILE=TRIM(filename), STATUS='old', &
                                                 FORM='formatted')
            READ(iu_rel,*)
            DO itemp=1, ntemp
               READ(iu_rel,'(15f15.7)') temp(itemp), dyde(istep,igeom,itemp)
               IF (itemp==5) WRITE(6,'(15f15.7)') temp(itemp), &
                                               dyde(istep,igeom,itemp)
            ENDDO
            CLOSE(UNIT=iu_rel, STATUS='KEEP')
         ENDIF    
      ENDIF
   ENDDO   
ENDIF
CALL mp_sum(dyde(:,igeom,:),world_comm)

RETURN
END SUBROUTINE read_dyde

!------------------------------------------------------------------------
SUBROUTINE redefine_energies_qua(energy_geo, epsilon_geo, epsil_geo, nwork,  &
                             energy_geo_eff, epsilon_geo_eff, nwork_eff)
!------------------------------------------------------------------------
!
! This routine is used when the elastic constants are calculated by
! moving atoms in some strain types. 
! It receives the list of energies and strain and remove all the energies
! that correspond to the same strain but different atomic positions and
! substitutes them with the energy at the minimum. To obtain it, it fits
! the energies as a function of the atomic position with a fourth
! order polynomial and finds the minimum.
!
USE kinds, ONLY : DP
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain, &
                                      stype, nmove, atom_step, min_y, &
                                      epsil_y, lzsisa, lfp
USE control_quartic_energy, ONLY : lsolve
USE polynomial,   ONLY : init_poly, clean_poly, poly2, poly4
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             evaluate_fit_quartic

IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
INTEGER, INTENT(OUT)  :: nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(OUT) :: energy_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo(nwork)
REAL(DP), INTENT(IN)  :: epsilon_geo(3,3,nwork)
REAL(DP), INTENT(OUT) :: epsilon_geo_eff(3,3,nwork)

INTEGER :: igeom, istep, igeo, iwork, jwork, imov, nvar
REAL(DP) :: emin, xmin(1)
REAL(DP), ALLOCATABLE :: ene(:), y(:)
TYPE(poly2) :: p2
TYPE(poly4) :: p4


nvar=1
ALLOCATE(ene(nmove))
ALLOCATE(y(nmove))

iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index
CALL init_poly(nvar,p2)
CALL init_poly(nvar,p4)
DO igeom=1, ngeom  
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
               WRITE(6,*) y(imov), ene(imov)
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
            min_y(igeo,istep,igeom)=xmin(1)
            epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
         ELSE
            jwork=jwork+1
            iwork=iwork+1
            energy_geo_eff(iwork)=energy_geo(jwork)
            epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)
         ENDIF
      ENDDO
   ENDDO
ENDDO
nwork_eff=iwork

DEALLOCATE(ene)
DEALLOCATE(y)
CALL clean_poly(p2)
CALL clean_poly(p4)

RETURN
END SUBROUTINE redefine_energies_qua
!
!------------------------------------------------------------------------
SUBROUTINE redefine_energies_qua_t(energy_geo, epsilon_geo, epsil_geo,  &
                nwork, energy_geo_eff, epsilon_geo_eff, nwork_eff, igeom, & 
                itemp)
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
USE control_elastic_constants, ONLY : ngeom, nstep_ec, ngeo_strain, &
                                      stype, nmove, atom_step, min_y_t, &
                                      epsil_y, min_y, lzsisa, lfp
USE control_quartic_energy, ONLY : lsolve
USE polynomial,   ONLY : init_poly, clean_poly, poly2, poly4
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_quadratic_extremum
USE quartic_surfaces, ONLY : fit_multi_quartic, find_quartic_extremum, &
                             evaluate_fit_quartic
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nwork
INTEGER, INTENT(IN)   :: igeom
INTEGER, INTENT(IN)   :: itemp
INTEGER, INTENT(OUT)  :: nwork_eff
REAL(DP), INTENT(IN)  :: energy_geo(nwork)
REAL(DP), INTENT(OUT) :: energy_geo_eff(nwork)
REAL(DP), INTENT(OUT) :: epsil_geo(nwork)
REAL(DP), INTENT(IN)  :: epsilon_geo(3,3,nwork)
REAL(DP), INTENT(OUT) :: epsilon_geo_eff(3,3,nwork)

INTEGER :: istep, igeo, iwork, jwork, imov, nvar
REAL(DP) :: a(3), emin, xmin(1)
REAL(DP), ALLOCATABLE :: ene(:), y(:)
TYPE(poly2) :: p2
TYPE(poly4) :: p4

nvar=1
ALLOCATE(ene(nmove))
ALLOCATE(y(nmove))
CALL init_poly(nvar,p2)
CALL init_poly(nvar,p4)

iwork=0     ! run on current energy index
jwork=0     ! run on previous energy index
WRITE(6,*) 'compressing geometry', igeom
DO istep=1, nstep_ec
   DO igeo=1, ngeo_strain
      IF (stype(istep)) THEN
         DO imov=1, nmove
            jwork=jwork+1
            WRITE(6,*) 'imov, jwork', imov, jwork
            ene(imov)=energy_geo(jwork)
!
!  y is the displacement (in a.u.) of the atom with respect to the uniformely
!  strained position
!
            y(imov)=(imov-(nmove+1.0_DP)/2.0_DP)*atom_step(istep)
            WRITE(6,*) y(imov), ene(imov)
         ENDDO
!          
!           Fit the data with a quartic and put the minimum free energy
!           or the free energy at the minimum of the energy (lsize) or
!           the free energy at the uniformely strained positions
!
         CALL fit_multi_quartic(nmove,1,lsolve,y,ene,p4)
         IF (lzsisa) THEN
            xmin(1)=min_y(igeo,istep,igeom)
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
         WRITE(6,'("igeom, istep, igeo, iwork, jwork",5i7)') igeom, istep, &
                                                          igeo, iwork, jwork
         min_y_t(igeo,istep,igeom,itemp)=xmin(1)
         epsil_y(igeo,istep,igeom)=epsil_geo(jwork)
      ELSE
         jwork=jwork+1
         iwork=iwork+1
         WRITE(6,'("igeom, istep, igeo, iwork, jwork",5i7)') igeom, istep, &
                                                       igeo, iwork, jwork
         energy_geo_eff(iwork)=energy_geo(jwork)
         epsilon_geo_eff(:,:,iwork)=epsilon_geo(:,:,jwork)
      ENDIF
   ENDDO
ENDDO
nwork_eff=iwork

DEALLOCATE(ene)
DEALLOCATE(y)
CALL clean_poly(p2)
CALL clean_poly(p4)

RETURN
END SUBROUTINE redefine_energies_qua_t

