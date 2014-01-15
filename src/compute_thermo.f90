SUBROUTINE compute_thermo(igeom)

USE kinds,          ONLY : DP
USE phdos_module,   ONLY : phdos_type, read_phdos_data, zero_point_energy, &
                          free_energy, vib_energy, vib_entropy, &
                          specific_heat_cv, integrated_dos
USE thermodynamics, ONLY : tmin, tmax, deltat, ntemp, ngeo, temp, ph_ener, &
                          ph_free_ener, ph_entropy, ph_cv, ph_cp, phdos_save
USE mp_images,      ONLY : root_image, my_image_id
USE io_global,      ONLY : ionode, stdout
USE control_thermo, ONLY : fltherm

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

INTEGER  :: i, ios
REAL(DP) :: e0, tot_states
INTEGER  :: itemp
INTEGER  :: iu_therm
!
IF (my_image_id /= root_image) RETURN

IF ( igeom < 1 .OR. igeom > ngeo ) CALL errore('print_thermo', & 
                                               'Too many geometries',1)
WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamical properties")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(fltherm)
WRITE(stdout,'(2x,76("+"),/)')

!
!  Allocate thermodynamic quantities
!
IF (deltat <= 0.0_8) CALL errore('print_thermo','Negative deltat',1)
ntemp=1+NINT((tmax-tmin)/deltat)

IF (.NOT.ALLOCATED(temp)) ALLOCATE(temp(ntemp,ngeo))

IF (.NOT.ALLOCATED(ph_free_ener)) ALLOCATE(ph_free_ener(ntemp,ngeo))
IF (.NOT.ALLOCATED(ph_ener)) ALLOCATE(ph_ener(ntemp,ngeo))
IF (.NOT.ALLOCATED(ph_entropy)) ALLOCATE(ph_entropy(ntemp,ngeo))
IF (.NOT.ALLOCATED(ph_cv)) ALLOCATE(ph_cv(ntemp,ngeo))
IF (.NOT.ALLOCATED(ph_cp)) ALLOCATE(ph_cp(ntemp,ngeo))

CALL zero_point_energy(phdos_save(igeom), e0)
CALL integrated_dos(phdos_save(igeom), tot_states)

DO itemp = 1, ntemp
   temp(itemp, igeom) = tmin + (itemp-1) * deltat
   CALL free_energy(phdos_save(igeom), temp(itemp,igeom), &
                                       ph_free_ener(itemp,igeom))
   CALL vib_energy(phdos_save(igeom), temp(itemp,igeom), &
                                       ph_ener(itemp, igeom))
   CALL vib_entropy(phdos_save(igeom), temp(itemp,igeom), &
                                       ph_entropy(itemp, igeom))
   CALL specific_heat_cv(phdos_save(igeom), temp(itemp,igeom), &
                                       ph_cv(itemp, igeom))
   ph_free_ener(itemp,igeom)=ph_free_ener(itemp,igeom)+e0
   ph_ener(itemp,igeom)=ph_ener(itemp,igeom)+e0
END DO

IF (ionode) THEN
   iu_therm=2
   OPEN (UNIT=iu_therm, FILE=TRIM(fltherm), STATUS='unknown',&
                                                     FORM='formatted')
   WRITE(iu_therm,'("# Zero point energy is:", f15.5, " Ry")') e0
   WRITE(iu_therm,'("# Total number of states is:", f15.5, " Ry")') tot_states
   WRITE(iu_therm,'("# Temperature T in K, ")')
   WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
   WRITE(iu_therm,'("# Entropy in Ry / (cell K)")')
   WRITE(iu_therm,'("# Specific heat Cv in Ry / (cell K)")')
   WRITE(iu_therm,'("# Multiply by 13.6058 to have energies in &
                       & eV / cell etc.")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  & energies in cal / mol")')
   WRITE(iu_therm,'("# assuming that N_A cells contain 1 mole &
                  & (N_A is the Avogadro")')
   WRITE(iu_therm,'("# number). If N_A cells contain N moles &
                  & (for instance in ")')
   WRITE(iu_therm,'("# silicon N=2) divide by N to have energies in &
                  & cal / mol etc. ")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  & have energies in J / mol")')
   WRITE(iu_therm,'("#",5x,"   T  ", 7x, " energy ", 4x, "  free energy ",&
                  & 4x, " entropy ", 7x, " Cv ")') 

   DO itemp = 1, ntemp
      WRITE(iu_therm, '(5e15.6)') temp(itemp,igeom), &
                    ph_ener(itemp,igeom), ph_free_ener(itemp,igeom), &
                    ph_entropy(itemp,igeom), ph_cv(itemp,igeom)
   END DO

   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE compute_thermo
