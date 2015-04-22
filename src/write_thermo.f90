!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_thermo(igeom)

USE kinds,          ONLY : DP
USE phdos_module,   ONLY : phdos_type, read_phdos_data, zero_point_energy, &
                           free_energy, vib_energy, vib_entropy, &
                           specific_heat_cv, integrated_dos
USE thermo_mod,     ONLY : tot_ngeo
USE temperature,    ONLY : tmin, tmax, deltat, ntemp, temp
USE thermodynamics, ONLY : ph_ener, ph_free_ener, ph_entropy, ph_cv, phdos_save
USE mp_images,      ONLY : root_image, my_image_id, intra_image_comm
USE mp,             ONLY : mp_bcast
USE io_global,      ONLY : ionode, ionode_id, stdout
USE control_thermo, ONLY : fltherm

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

INTEGER  :: i, ios
REAL(DP) :: e0, tot_states
INTEGER  :: itemp
INTEGER  :: iu_therm, idum
LOGICAL  :: check_file_exists, do_read
!
do_read=.FALSE.
IF ( check_file_exists(fltherm) ) do_read=.TRUE.
IF (my_image_id /= root_image) RETURN

IF (do_read) THEN
   IF (ionode) THEN
      iu_therm=2
      OPEN (UNIT=iu_therm, FILE=TRIM(fltherm), STATUS='unknown',&
                                                     FORM='formatted')
      DO idum=1,12
         READ(iu_therm,*)
      ENDDO
      DO itemp = 1, ntemp
         READ(iu_therm, '(e16.8,4e20.12)') temp(itemp), &
                    ph_ener(itemp,igeom), ph_free_ener(itemp,igeom), &
                    ph_entropy(itemp,igeom), ph_cv(itemp,igeom)
      END DO
      CLOSE(iu_therm)
   END IF
   CALL mp_bcast(ph_ener(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(ph_free_ener(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(ph_entropy(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(ph_cv(:,igeom), ionode_id, intra_image_comm)
   RETURN
END IF


IF ( igeom < 1 .OR. igeom > tot_ngeo ) CALL errore('print_thermo', & 
                                               'Too many geometries',1)
WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamical properties from phonon dos")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(fltherm)
WRITE(stdout,'(2x,76("+"),/)')

!
!  Allocate thermodynamic quantities
!

CALL zero_point_energy(phdos_save(igeom), e0)
CALL integrated_dos(phdos_save(igeom), tot_states)

DO itemp = 1, ntemp
   CALL free_energy(phdos_save(igeom), temp(itemp), ph_free_ener(itemp,igeom))
   CALL vib_energy(phdos_save(igeom), temp(itemp), ph_ener(itemp,igeom))
   CALL vib_entropy(phdos_save(igeom), temp(itemp), ph_entropy(itemp, igeom))
   CALL specific_heat_cv(phdos_save(igeom), temp(itemp), ph_cv(itemp, igeom))
   ph_free_ener(itemp,igeom)=ph_free_ener(itemp,igeom)+e0
   ph_ener(itemp,igeom)=ph_ener(itemp,igeom)+e0
END DO

IF (ionode) THEN
   iu_therm=2
   OPEN (UNIT=iu_therm, FILE=TRIM(fltherm), STATUS='unknown',&
                                                     FORM='formatted')
   WRITE(iu_therm,'("# Zero point energy is:", f9.5, " Ry/cell,", f9.5, &
                    &" kJ*N/mol,", f9.5, " kcal*N/mol")') e0, &
                       e0 * 1313.313_DP, e0 * 313.7545_DP 
   WRITE(iu_therm,'("# Total number of states is:", f15.5,",")') tot_states
   WRITE(iu_therm,'("# Temperature T in K, ")')
   WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
   WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
   WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
   WRITE(iu_therm,'("# Multiply by 13.6058 to have energies in &
                       &eV/cell etc..")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  &energies in cal*N/mol.")')
   WRITE(iu_therm,'("# We assume that N_A cells contain N moles &
                  &(N_A is the Avogadro number).")')
   WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                   &energies in cal/mol etc. ")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  &have energies in J/mol.")')
   WRITE(iu_therm,'("#",5x,"   T  ", 7x, " energy ", 4x, "  free energy ",&
                  & 4x, " entropy ", 7x, " Cv ")') 

   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,4e20.12)') temp(itemp), &
                    ph_ener(itemp,igeom), ph_free_ener(itemp,igeom), &
                    ph_entropy(itemp,igeom), ph_cv(itemp,igeom)
   END DO

   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_thermo

SUBROUTINE write_thermo_ph(igeom)

USE kinds,            ONLY : DP
USE ph_freq_module,   ONLY : ph_freq_type, zero_point_energy_ph, &
                          free_energy_ph, vib_energy_ph, vib_entropy_ph, &
                          specific_heat_cv_ph
USE temperature,      ONLY : tmin, tmax, deltat, ntemp, temp
USE ph_freq_thermodynamics, ONLY : phf_ener, phf_free_ener, phf_entropy, &
                           phf_cv, ph_freq_save
USE thermo_mod,       ONLY : tot_ngeo
USE mp_images,        ONLY : root_image, my_image_id, intra_image_comm
USE mp,               ONLY : mp_bcast
USE io_global,        ONLY : ionode, ionode_id, stdout
USE control_thermo,   ONLY : fltherm

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: filename
LOGICAL :: check_file_exists, do_read

INTEGER  :: i, ios
REAL(DP) :: e0
INTEGER  :: itemp
INTEGER  :: iu_therm, idum
!
do_read=.FALSE.
filename=TRIM(fltherm)//'_ph'
IF ( check_file_exists(filename) ) do_read=.TRUE.

IF (my_image_id /= root_image) RETURN

IF (do_read) THEN
   IF (ionode) THEN
      iu_therm=2
      OPEN (UNIT=iu_therm, FILE=TRIM(filename), STATUS='old',&
                                                     FORM='formatted')
      DO idum=1,11
         READ(iu_therm,*)
      ENDDO
      DO itemp = 1, ntemp
         READ(iu_therm, '(e16.8,4e20.12)') temp(itemp), &
                    phf_ener(itemp,igeom), phf_free_ener(itemp,igeom), &
                    phf_entropy(itemp,igeom), phf_cv(itemp,igeom)
      END DO
      CLOSE(iu_therm)
   END IF
   CALL mp_bcast(phf_ener(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(phf_free_ener(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(phf_entropy(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(phf_cv(:,igeom), ionode_id, intra_image_comm)
   RETURN
END IF


IF ( igeom < 1 .OR. igeom > tot_ngeo ) CALL errore('print_thermo', & 
                                               'Too many geometries',1)
WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamical properties from frequencies")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filename)
WRITE(stdout,'(2x,76("+"),/)')

!
!  Allocate thermodynamic quantities
!
CALL zero_point_energy_ph(ph_freq_save(igeom), e0)

DO itemp = 1, ntemp
   IF (MOD(itemp,30)==0) WRITE(6,'(5x,"Computing temperature ", i5)') itemp
   CALL free_energy_ph(ph_freq_save(igeom), temp(itemp), &
                                       phf_free_ener(itemp,igeom))
   CALL vib_energy_ph(ph_freq_save(igeom), temp(itemp), &
                                       phf_ener(itemp, igeom))
   CALL specific_heat_cv_ph(ph_freq_save(igeom), temp(itemp), &
                                       phf_cv(itemp, igeom))
   phf_free_ener(itemp,igeom)=phf_free_ener(itemp,igeom)+e0
   phf_ener(itemp,igeom)=phf_ener(itemp,igeom)+e0
END DO
phf_entropy(:,igeom)=(phf_ener(:, igeom) - phf_free_ener(:, igeom))/  &
                        temp(:)
IF (ionode) THEN
   iu_therm=2
   OPEN (UNIT=iu_therm, FILE=TRIM(filename), STATUS='unknown',&
                                                     FORM='formatted')
   WRITE(iu_therm,'("# Zero point energy is:", f9.5, " Ry/cell,", f9.5, &
                    &" kJ*N/mol,", f9.5, " kcal*N/mol")') e0, &
                       e0 * 1313.313_DP, e0 * 313.7545_DP 
   WRITE(iu_therm,'("# Temperature T in K, ")')
   WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
   WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
   WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
   WRITE(iu_therm,'("# Multiply by 13.6058 to have energies in &
                       &eV/cell etc..")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  &energies in cal*N/mol.")')
   WRITE(iu_therm,'("# We assume that N_A cells contain N moles &
                  &(N_A is the Avogadro number).")')
   WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                   &energies in cal/mol etc. ")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  &have energies in J/mol.")')
   WRITE(iu_therm,'("#",5x,"   T  ", 7x, " energy ", 4x, "  free energy ",&
                  & 4x, " entropy ", 7x, " Cv ")') 

   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,4e20.12)') temp(itemp), &
                    phf_ener(itemp,igeom), phf_free_ener(itemp,igeom), &
                    phf_entropy(itemp,igeom), phf_cv(itemp,igeom)
   END DO

   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_thermo_ph
