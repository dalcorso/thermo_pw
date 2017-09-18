!
! Copyright (C) 2014-15 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_thermo(igeom)
!
!  This routine writes on file the harmonic thermodynamical quantities
!
USE kinds,          ONLY : DP
USE phdos_module,   ONLY : zero_point_energy, fecv, integrated_dos
USE thermo_mod,     ONLY : tot_ngeo
USE temperature,    ONLY : ntemp, temp
USE thermodynamics, ONLY : ph_ener, ph_free_ener, ph_entropy, ph_cv, phdos_save
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_bcast, mp_sum
USE io_global,      ONLY : meta_ionode, meta_ionode_id, stdout
USE data_files,     ONLY : fltherm

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

INTEGER  :: idum, itemp, startt, lastt, iu_therm
INTEGER  :: find_free_unit
REAL(DP) :: e0, tot_states
CHARACTER(LEN=256) :: filetherm
LOGICAL  :: check_file_exists, do_read
!
do_read=.FALSE.
filetherm="therm_files/"//TRIM(fltherm)
IF ( check_file_exists(filetherm) ) do_read=.TRUE.
IF ( igeom < 1 .OR. igeom > tot_ngeo ) CALL errore('write_thermo', & 
                                               'Too many geometries',1)
IF (do_read) THEN
   IF (meta_ionode) THEN
      iu_therm=find_free_unit()
      OPEN (UNIT=iu_therm, FILE=TRIM(filetherm), STATUS='old',&
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
   CALL mp_bcast(ph_ener(:,igeom), meta_ionode_id, world_comm)
   CALL mp_bcast(ph_free_ener(:,igeom), meta_ionode_id, world_comm)
   CALL mp_bcast(ph_entropy(:,igeom), meta_ionode_id, world_comm)
   CALL mp_bcast(ph_cv(:,igeom), meta_ionode_id, world_comm)
   RETURN
END IF

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamic properties from phonon dos")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filetherm)
WRITE(stdout,'(2x,76("+"),/)')

CALL zero_point_energy(phdos_save(igeom), e0)
CALL integrated_dos(phdos_save(igeom), tot_states)
!
!  Divide the temperatures among processors
!
CALL divide (world_comm, ntemp, startt, lastt)
ph_free_ener(:,igeom)=0.0_DP
ph_ener(:,igeom)=0.0_DP
ph_cv(:,igeom)=0.0_DP
ph_entropy(:,igeom)=0.0_DP
DO itemp = startt, lastt
   IF (MOD(itemp-startt+1,30)==0) &
                     WRITE(6,'(5x,"Computing temperature ", i5 " / ",&
       & i5, 4x," T=",f12.2," K")') itemp-startt+1, lastt-startt+1, temp(itemp)

   CALL fecv(phdos_save(igeom), temp(itemp), ph_free_ener(itemp, igeom), &
                                  ph_ener(itemp,igeom), ph_cv(itemp, igeom))
   ph_free_ener(itemp,igeom)=ph_free_ener(itemp,igeom)+e0
   ph_ener(itemp,igeom)=ph_ener(itemp,igeom)+e0
   ph_entropy(itemp,igeom)=(ph_ener(itemp, igeom)-ph_free_ener(itemp,igeom))/&
                            temp(itemp)
END DO
!
!  and collect the results
!
CALL mp_sum(ph_free_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(ph_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(ph_cv(1:ntemp,igeom),world_comm)
CALL mp_sum(ph_entropy(1:ntemp,igeom),world_comm)

IF (meta_ionode) &
   CALL write_thermo_info(e0, tot_states, ntemp, temp, ph_ener(1,igeom), &
              ph_free_ener(1,igeom), ph_entropy(1,igeom), ph_cv(1,igeom),&
                                                             1,filetherm)

RETURN
END SUBROUTINE write_thermo

SUBROUTINE write_thermo_ph(igeom)

USE kinds,            ONLY : DP
USE ph_freq_module,   ONLY : ph_freq_type, zero_point_energy_ph, fecv_ph, &
                             free_energy_ph, vib_energy_ph, specific_heat_cv_ph
USE temperature,      ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : phf_ener, phf_free_ener, phf_entropy, &
                             phf_cv, ph_freq_save
USE thermo_mod,       ONLY : tot_ngeo
USE mp,               ONLY : mp_bcast, mp_sum
USE mp_world,         ONLY : world_comm
USE io_global,        ONLY : meta_ionode, meta_ionode_id, stdout
USE data_files,       ONLY : fltherm

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: filename
LOGICAL :: check_file_exists, do_read

INTEGER  :: idum, itemp, iu_therm
INTEGER  :: find_free_unit
REAL(DP) :: e0
!
do_read=.FALSE.
filename="therm_files/"//TRIM(fltherm)//'_ph'
IF ( check_file_exists(filename) ) do_read=.TRUE.

IF ( igeom < 1 .OR. igeom > tot_ngeo ) CALL errore('write_thermo', & 
                                               'Too many geometries',1)

IF (do_read) THEN
   IF (meta_ionode) THEN
      iu_therm=find_free_unit()
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
   CALL mp_bcast(phf_ener(:,igeom), meta_ionode_id, world_comm)
   CALL mp_bcast(phf_free_ener(:,igeom), meta_ionode_id, world_comm)
   CALL mp_bcast(phf_entropy(:,igeom), meta_ionode_id, world_comm)
   CALL mp_bcast(phf_cv(:,igeom), meta_ionode_id, world_comm)
   RETURN
END IF

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamic properties from frequencies")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filename)
WRITE(stdout,'(2x,76("+"),/)')
!
!  Allocate thermodynamic quantities
!
CALL zero_point_energy_ph(ph_freq_save(igeom), e0)
!
!  Divide the temperatures among processors
!
phf_free_ener(:,igeom)=0.0_DP
phf_ener(:,igeom)=0.0_DP
phf_cv(:,igeom)=0.0_DP
phf_entropy(:,igeom)=0.0_DP
DO itemp = 1, ntemp
   IF (MOD(itemp,30)==0) &
        WRITE(stdout,'(5x,"Computing temperature ", i5 " / ",&
        & i5,4x," T=",f12.2," K")') itemp, ntemp, temp(itemp)
   CALL fecv_ph(ph_freq_save(igeom), temp(itemp), phf_free_ener(itemp,igeom), &
                      phf_ener(itemp, igeom), phf_cv(itemp, igeom))
   phf_free_ener(itemp,igeom)=phf_free_ener(itemp,igeom)+e0
   phf_ener(itemp,igeom)=phf_ener(itemp,igeom)+e0
   phf_entropy(itemp,igeom)=(phf_ener(itemp, igeom) - &
                             phf_free_ener(itemp, igeom))/temp(itemp)
END DO
!
!  In ph_freq_save the frequencies are distributed among all processors
!  so the resulting thermodynamical quantities must be collected
!
CALL mp_sum(phf_free_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(phf_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(phf_cv(1:ntemp,igeom),world_comm)
CALL mp_sum(phf_entropy(1:ntemp,igeom),world_comm)

IF (meta_ionode) &
   CALL write_thermo_info(e0, 0.0_DP, ntemp, temp, phf_ener(1,igeom), &
              phf_free_ener(1,igeom), phf_entropy(1,igeom), phf_cv(1,igeom),& 
                                                            2,filename)
RETURN
END SUBROUTINE write_thermo_ph

SUBROUTINE write_thermo_debye(igeom)

USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE debye_module,     ONLY : debye_e0, debye_vib_energy, debye_free_energy, &
                             debye_entropy, debye_cv
USE control_debye,    ONLY : deb_e0, deb_cv, deb_entropy, deb_energy, &
                             deb_free_energy, debye_t
USE temperature,      ONLY : ntemp, temp
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode, stdout
USE data_files,       ONLY : fltherm

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6) :: int_to_char

INTEGER  :: idum
INTEGER  :: itemp
!
filename='therm_files/'//TRIM(fltherm)//'_debye.g'//TRIM(int_to_char(igeom))

IF (my_image_id /= root_image) RETURN

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamic properties from elastic constants")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filename)
WRITE(stdout,'(2x,76("+"),/)')
!
!  Compute thermodynamic quantities
!
CALL debye_e0 (debye_t, nat, deb_e0)
CALL debye_cv (debye_t, temp, ntemp, nat, deb_cv)
CALL debye_vib_energy (debye_t, temp, ntemp, nat, deb_energy)
CALL debye_free_energy (debye_t, temp, ntemp, nat, deb_free_energy)
CALL debye_entropy (debye_t, temp, ntemp, nat, deb_entropy)
!
!  Add the zero point energy
!
DO itemp=1,ntemp
   deb_energy(itemp) = deb_energy(itemp) + deb_e0
   deb_free_energy(itemp) = deb_free_energy(itemp) + deb_e0
END DO
!
!  Write on file
!
 IF (ionode) &
    CALL write_thermo_info(deb_e0, debye_t, ntemp, temp, deb_energy, &
               deb_free_energy, deb_entropy, deb_cv, 3, filename)

RETURN
END SUBROUTINE write_thermo_debye

SUBROUTINE write_thermo_info(e0, tot_states, ntemp, temp, energy, &
                                 free_energy, entropy, cv, iflag, filename)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, iflag
REAL(DP), INTENT(IN) :: e0, tot_states, temp(ntemp), energy(ntemp),  &
                        free_energy(ntemp), entropy(ntemp), cv(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: iu_therm, itemp
INTEGER :: find_free_unit

iu_therm=find_free_unit()
OPEN (UNIT=iu_therm, FILE=TRIM(filename), STATUS='unknown',&
                                                     FORM='formatted')
WRITE(iu_therm,'("# Zero point energy:", f8.5, " Ry/cell,", f9.5, &
                 &" kJ/(N mol),", f9.5, " kcal/(N mol)")') e0, &
                    e0 * 1313.313_DP, e0 * 313.7545_DP 
IF (iflag==3) THEN
   WRITE(iu_therm,'("# Temperature T in K, Debye temperature=",f12.3, &
                                          &" K,")') tot_states
ELSE
   WRITE(iu_therm,'("# Temperature T in K, ")')
ENDIF
IF (iflag==1) &
WRITE(iu_therm,'("# Total number of states is:", f15.5,",")') tot_states
WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
WRITE(iu_therm,'("# Multiply by 13.6058 to have energies in &
                       &eV/cell etc..")')
WRITE(iu_therm,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  &energies in cal/(N mol).")')
WRITE(iu_therm,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  &have energies in J/(N mol).")')
WRITE(iu_therm,'("# N is the number of formula units per cell.")')
WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                &energies in cal/mol etc. ")')
WRITE(iu_therm,'("#",5x,"   T  ", 10x, " energy ", 9x, "  free energy ",&
               & 9x, " entropy ", 12x, " Cv ")')

DO itemp = 1, ntemp
   WRITE(iu_therm, '(e16.8,4e20.12)') temp(itemp), energy(itemp), &
                 free_energy(itemp), entropy(itemp), cv(itemp)
END DO
CLOSE(iu_therm)

RETURN
END SUBROUTINE write_thermo_info
