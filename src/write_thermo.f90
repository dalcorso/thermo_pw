!
! Copyright (C) 2014-15 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_thermo(igeom)
!-----------------------------------------------------------------------
!
!  This routine computes and writes on file the harmonic thermodynamic 
!  quantities calculated using the phonon dos
!
USE kinds,          ONLY : DP
USE ions_base,      ONLY : nat, nsp, ityp, amass
USE phdos_module,   ONLY : zero_point_energy, fecv, integrated_dos, &
                           phdos_debye_factor
USE thermo_mod,     ONLY : tot_ngeo
USE temperature,    ONLY : ntemp, temp
USE thermodynamics, ONLY : ph_ener, ph_free_ener, ph_entropy, ph_ce, &
                           ph_e0, ph_t_debye, ph_b_fact, phdos_save, &
                           gen_phdos_save
USE control_thermo, ONLY : with_eigen
USE control_debye,  ONLY : idebye
USE data_files,     ONLY : fltherm
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum
USE io_global,      ONLY : meta_ionode, stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

INTEGER  :: itemp, startt, lastt, iu_therm
INTEGER  :: find_free_unit
REAL(DP) :: e0, tot_states
CHARACTER(LEN=256) :: filetherm
LOGICAL  :: do_read
!
IF ( igeom < 1 .OR. igeom > tot_ngeo ) CALL errore('write_thermo', & 
                                               'Too many geometries',1)

filetherm='therm_files/'//TRIM(fltherm)
CALL read_thermo(ntemp, temp, ph_ener(1,igeom), ph_free_ener(1,igeom), &
                 ph_entropy(1,igeom), ph_ce(1,igeom), do_read, filetherm)
IF (with_eigen.AND.do_read) &
   CALL read_b_factor(ntemp, ph_b_fact(1,1,1,1,igeom),filetherm) 
IF (do_read) RETURN

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamic properties from phonon dos")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filetherm)
WRITE(stdout,'(2x,76("+"),/)')

CALL zero_point_energy(phdos_save(igeom), e0)
ph_e0(igeom)=e0
CALL integrated_dos(phdos_save(igeom), tot_states)
!
!  Divide the temperatures among processors
!
CALL divide (world_comm, ntemp, startt, lastt)
ph_free_ener(:,igeom)=0.0_DP
ph_ener(:,igeom)=0.0_DP
ph_ce(:,igeom)=0.0_DP
ph_entropy(:,igeom)=0.0_DP
ph_b_fact(:,:,:,:,igeom)=0.0_DP
DO itemp = startt, lastt
   IF (MOD(itemp-startt+1,30)==0) &
                     WRITE(6,'(5x,"Computing temperature ", i5, " / ",&
       & i5, 4x," T=",f12.2," K")') itemp-startt+1, lastt-startt+1, temp(itemp)

   CALL fecv(phdos_save(igeom), temp(itemp), ph_free_ener(itemp, igeom), &
                                  ph_ener(itemp,igeom), ph_ce(itemp, igeom))
   ph_free_ener(itemp,igeom)=ph_free_ener(itemp,igeom)!+e0
   ph_ener(itemp,igeom)=ph_ener(itemp,igeom)!+e0
   ph_entropy(itemp,igeom)=(ph_ener(itemp, igeom)-ph_free_ener(itemp,igeom))/&
                            temp(itemp)
   IF (with_eigen) THEN
      CALL phdos_debye_factor(gen_phdos_save, temp(itemp), &
                          ph_b_fact(1,1,1,itemp,igeom), nat, amass, nsp, ityp)
   ENDIF

END DO
!
!  and collect the results
!
CALL mp_sum(ph_free_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(ph_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(ph_ce(1:ntemp,igeom),world_comm)
CALL mp_sum(ph_entropy(1:ntemp,igeom),world_comm)
IF (with_eigen) CALL mp_sum(ph_b_fact(1:3,1:3,1:nat,1:ntemp,igeom),world_comm)
IF (idebye==1) THEN
   CALL find_t_debye(ph_free_ener(1:ntemp,igeom), temp, ntemp, &
                                             ph_t_debye(1:ntemp,igeom))
ELSEIF(idebye==2) THEN
   CALL find_t_debye_ene(ph_ener(1:ntemp,igeom), temp, ntemp, &
                                             ph_t_debye(1:ntemp,igeom))
ELSEIF(idebye==3) THEN
   CALL find_t_debye_cv(ph_ce(1:ntemp,igeom), temp, ntemp, &
                                             ph_t_debye(1:ntemp,igeom))
ENDIF
ph_free_ener(1:ntemp,igeom) = ph_free_ener(1:ntemp,igeom)+e0
ph_ener(1:ntemp,igeom) = ph_ener(1:ntemp,igeom)+e0

IF (meta_ionode) THEN
   CALL write_thermo_info(e0, tot_states, ntemp, temp, ph_ener(1,igeom), &
              ph_free_ener(1,igeom), ph_entropy(1,igeom), ph_ce(1,igeom),&
                                                             1,filetherm)
   IF (idebye/=0) CALL write_thermo_info_deb(temp, ph_t_debye(1,igeom),  &
                                                          ntemp, filetherm)
ENDIF

IF (meta_ionode.AND.with_eigen) &
   CALL write_dw_info(ntemp, temp, ph_b_fact(1,1,1,1,igeom), nat, filetherm)

RETURN
END SUBROUTINE write_thermo
!
!-----------------------------------------------------------------------
SUBROUTINE write_thermo_ph(igeom)
!-----------------------------------------------------------------------
!  This routine computes and writes on file the harmonic thermodynamic 
!  quantities calculated using the direct sum over the phonon frequencies
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat, nsp, ityp, amass
USE symme,            ONLY : symtensor
USE ph_freq_module,   ONLY : ph_freq_type, zero_point_energy_ph, fecv_ph, &
                             free_energy_ph, vib_energy_ph, &
                             specific_heat_cv_ph, debye_waller_factor
USE temperature,      ONLY : ntemp, temp
USE ph_freq_thermodynamics, ONLY : phf_ener, phf_free_ener, phf_entropy, &
                             phf_e0, phf_ce, phf_b_fact, ph_freq_save,   &
                             phf_t_debye
USE control_debye,    ONLY : idebye
USE thermo_mod,       ONLY : tot_ngeo
USE control_thermo,   ONLY : with_eigen
USE control_debye,    ONLY : idebye
USE data_files,       ONLY : fltherm
USE mp,               ONLY : mp_sum
USE mp_world,         ONLY : world_comm
USE io_global,        ONLY : meta_ionode, stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: filename
LOGICAL :: do_read

INTEGER  :: itemp, iu_therm, ipol, jpol, na
INTEGER  :: find_free_unit
REAL(DP) :: e0
!
IF ( igeom < 1 .OR. igeom > tot_ngeo ) CALL errore('write_thermo', & 
                                               'Too many geometries',1)
filename='therm_files/'//TRIM(fltherm)//'_ph'
CALL read_thermo(ntemp, temp, phf_ener(1,igeom), phf_free_ener(1,igeom), &
                 phf_entropy(1,igeom), phf_ce(1,igeom), do_read, filename)
IF (with_eigen.AND.do_read) &
   CALL read_b_factor(ntemp, phf_b_fact(1,1,1,1,igeom),filename)
IF (do_read) RETURN

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamic properties from frequencies")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filename)
WRITE(stdout,'(2x,76("+"),/)')
!
!  Compute the zero point energy
!
CALL zero_point_energy_ph(ph_freq_save(igeom), e0)
phf_e0(igeom)=e0
!
!  Now compute the other thermodynamic quantities. Note that the each
!  processor adds its own q points so the thermodynamic quantities need
!  to be collected among processors.
!
phf_free_ener(:,igeom)=0.0_DP
phf_ener(:,igeom)=0.0_DP
phf_ce(:,igeom)=0.0_DP
phf_entropy(:,igeom)=0.0_DP
DO itemp = 1, ntemp
   IF (MOD(itemp,30)==0) &
        WRITE(stdout,'(5x,"Computing temperature ", i5, " / ",&
        & i5,4x," T=",f12.2," K")') itemp, ntemp, temp(itemp)
   CALL fecv_ph(ph_freq_save(igeom), temp(itemp), phf_free_ener(itemp,igeom), &
                      phf_ener(itemp, igeom), phf_ce(itemp, igeom))
   phf_free_ener(itemp,igeom)=phf_free_ener(itemp,igeom)!+e0
   phf_ener(itemp,igeom)=phf_ener(itemp,igeom)!+e0
   phf_entropy(itemp,igeom)=(phf_ener(itemp, igeom) - &
                             phf_free_ener(itemp, igeom))/temp(itemp)
   IF (with_eigen) THEN
      CALL debye_waller_factor(ph_freq_save(igeom), temp(itemp), &
                          phf_b_fact(1,1,1,itemp,igeom), nat, amass, nsp, ityp)
      CALL symtensor(nat, phf_b_fact(1,1,1,itemp,igeom))
   ENDIF
END DO
!
!  In ph_freq_save the frequencies are distributed among all processors
!  so the resulting thermodynamical quantities must be collected
!
CALL mp_sum(phf_e0(igeom),world_comm)
CALL mp_sum(phf_free_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(phf_ener(1:ntemp,igeom),world_comm)
CALL mp_sum(phf_ce(1:ntemp,igeom),world_comm)
CALL mp_sum(phf_entropy(1:ntemp,igeom),world_comm)
CALL mp_sum(e0,world_comm)
IF (with_eigen) CALL mp_sum(phf_b_fact(1:3,1:3,1:nat,1:ntemp,igeom),world_comm)

IF (idebye==1) THEN
   CALL find_t_debye(phf_free_ener(1:ntemp,igeom), temp, ntemp, &
                                             phf_t_debye(1:ntemp,igeom))
ELSEIF(idebye==2) THEN
   CALL find_t_debye_ene(phf_ener(1:ntemp,igeom), temp, ntemp, &
                                             phf_t_debye(1:ntemp,igeom))
ELSEIF(idebye==3) THEN
   CALL find_t_debye_cv(phf_ce(1:ntemp,igeom), temp, ntemp, &
                                             phf_t_debye(1:ntemp,igeom))
ENDIF

phf_free_ener(1:ntemp,igeom) = phf_free_ener(1:ntemp,igeom)+e0
phf_ener(1:ntemp,igeom) = phf_ener(1:ntemp,igeom)+e0

IF (meta_ionode) THEN
   CALL write_thermo_info(phf_e0(igeom), 0.0_DP, ntemp, temp, &
              phf_ener(1,igeom), &
              phf_free_ener(1,igeom), phf_entropy(1,igeom), phf_ce(1,igeom),& 
                                                            2,filename)
   IF (idebye/=0) CALL write_thermo_info_deb(temp, phf_t_debye(1,igeom), &
                                                   ntemp, filename)
ENDIF
IF (meta_ionode.AND.with_eigen) &
   CALL write_dw_info(ntemp, temp, phf_b_fact(1,1,1,1,igeom), nat, filename)

RETURN
END SUBROUTINE write_thermo_ph

!-----------------------------------------------------------------------
SUBROUTINE write_thermo_debye(igeom)
!-----------------------------------------------------------------------
! 
!  This routine computes and writes on file the harmonic thermodynamic 
!  quantities calculated using the Debye model. It receives as input only
!  the debye temperature
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat, nsp, amass, ityp
USE debye_module,     ONLY : debye_e0, debye_vib_energy, debye_free_energy, &
                             debye_entropy, debye_cv, debye_b_factor
USE control_debye,    ONLY : deb_e0, deb_cv, deb_entropy, deb_energy, &
                             deb_free_energy, deb_b_fact, deb_bfact, debye_t
USE temperature,      ONLY : ntemp, temp
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode, stdout
USE data_files,       ONLY : fltherm

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6) :: int_to_char

INTEGER  :: itemp
LOGICAL  :: do_read
!
filename='therm_files/'//TRIM(fltherm)//'_debye.g'//TRIM(int_to_char(igeom))

CALL read_thermo(ntemp, temp, deb_energy, deb_free_energy, deb_entropy, &
                        deb_cv, do_read, filename)
IF (do_read) RETURN

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
IF (nsp==1) CALL debye_b_factor(debye_t, temp, ntemp, amass(1), deb_bfact)
!
!  Add the zero point energy
!
DO itemp=1,ntemp
   deb_energy(itemp) = deb_energy(itemp) + deb_e0
   deb_free_energy(itemp) = deb_free_energy(itemp) + deb_e0
END DO
!
!  Write on file the information
!
IF (ionode) THEN
   CALL write_thermo_info(deb_e0, debye_t, ntemp, temp, deb_energy, &
               deb_free_energy, deb_entropy, deb_cv, 3, filename)
   IF (nsp==1) CALL write_dw_debye(ntemp, temp, deb_bfact, filename)
END IF

RETURN
END SUBROUTINE write_thermo_debye
!
!-----------------------------------------------------------------------
SUBROUTINE write_thermo_info(e0, tot_states, ntemp, temp, energy, &
                                 free_energy, entropy, cv, iflag, filename)
!-----------------------------------------------------------------------
!
! This routine writes on file the text inside the files with the 
! thermodynamic quantities that explain the units and the quantities
! contained in the file
!
USE kinds,        ONLY : DP
USE constants,    ONLY : rytoev, electronvolt_si, rydberg_si, avogadro
IMPLICIT NONE
REAL(DP), PARAMETER :: caltoj=4.184_DP
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
             e0*rydberg_si*avogadro/1.D3, e0*rydberg_si*avogadro/caltoj/1.D3
IF (iflag==3) THEN
   WRITE(iu_therm,'("# Temperature T in K, Debye temperature=",f12.3, &
                                          &" K,")') tot_states
ELSE
   WRITE(iu_therm,'("# Temperature T in K, ")')
ENDIF
IF (iflag==1) THEN
   WRITE(iu_therm,'("# Total number of states is:", f15.5,",")') tot_states
ELSE
   WRITE(iu_therm,'("# ")') 
ENDIF
WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
WRITE(iu_therm,'("# Multiply by ",f7.4," to have energies in &
                       &eV/cell etc..")') rytoev
WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to have &
                  &energies in cal/(N mol).")') rytoev, electronvolt_si &
                         * avogadro / caltoj, rydberg_si*avogadro/caltoj
WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to &
                  &have energies in J/(N mol).")') rytoev, electronvolt_si&
                         * avogadro, rydberg_si*avogadro
WRITE(iu_therm,'("# N is the number of formula units per cell.")')
WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                &energies in cal/mol etc. ")')
WRITE(iu_therm,'("#",5x,"   T  ", 10x, " energy ", 9x, "  free energy ",&
               & 9x, " entropy ", 12x, " Cv ")')

DO itemp = 1, ntemp
   WRITE(iu_therm, '(e16.8,4e20.12)') temp(itemp), energy(itemp), &
                 free_energy(itemp), entropy(itemp), cv(itemp)
END DO
CLOSE(UNIT=iu_therm, STATUS='KEEP')

RETURN
END SUBROUTINE write_thermo_info
!
!-----------------------------------------------------------------------
SUBROUTINE write_thermo_info_deb(temp, debye_t, ntemp, filename)
!-----------------------------------------------------------------------
!
!
USE kinds,        ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), debye_t(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: iu_therm, itemp
INTEGER :: find_free_unit

iu_therm=find_free_unit()
OPEN (UNIT=iu_therm, FILE=TRIM(filename)//".tdeb", STATUS='unknown', FORM='formatted')
!
!  Debye temperature on the first point is not reliable.
!
DO itemp = 2, ntemp
   WRITE(iu_therm, '(e16.8,e20.12)') temp(itemp), debye_t(itemp)
ENDDO

CLOSE(UNIT=iu_therm, STATUS='KEEP')

RETURN
END SUBROUTINE write_thermo_info_deb

!-----------------------------------------------------------------------
SUBROUTINE write_dw_info(ntemp, temp, b_fact, nat, filename)
!-----------------------------------------------------------------------
!
!  This routine writes on files the debye-waller factor as a function of
!  temperature
!
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, nat
REAL(DP), INTENT(IN) :: temp(ntemp), b_fact(3,3,nat,ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: iu_therm, itemp, na, ipol, jpol
INTEGER :: find_free_unit
CHARACTER(LEN=6) :: int_to_char

iu_therm=find_free_unit()
DO na=1,nat
   OPEN (UNIT=iu_therm, FILE=TRIM(filename)//'.'//&
                             TRIM(int_to_char(na))//'.dw', STATUS='unknown',&
                                               FORM='formatted')
   WRITE(iu_therm,'("#",6x,"T",14x,"B_11",14x,"B_12",14x,"B_13",14x,"B_22",&
                        &14x,"B_23",14x,"B_33")')

   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,6e18.8)') temp(itemp), &
                        ((b_fact(ipol,jpol,na,itemp), jpol=ipol,3), ipol=1,3)
   END DO
   CLOSE(UNIT=iu_therm, STATUS='KEEP')
END DO

RETURN
END SUBROUTINE write_dw_info
!
!-----------------------------------------------------------------------
SUBROUTINE read_thermo(ntemp, temp, ph_ener, ph_free_ener, ph_entropy, &
                       ph_ce, do_read, filetherm)
!-----------------------------------------------------------------------
!
!  This routine reads a file that contains the harmonic thermodynamic
!  quantities. It must be called by all processors, only the meta_ionode
!  reads and broadcast the thermodynamic quantities.
!  When with_eigen is true it reads also the atomic b factor.
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE io_global,        ONLY : meta_ionode, meta_ionode_id
USE control_thermo,   ONLY : with_eigen
USE mp_world,         ONLY : world_comm
USE mp,               ONLY : mp_bcast

IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(INOUT) :: temp(ntemp), ph_ener(ntemp), ph_free_ener(ntemp),&
                           ph_entropy(ntemp), ph_ce(ntemp)
LOGICAL, INTENT(OUT) :: do_read

CHARACTER(LEN=*) :: filetherm
CHARACTER(LEN=256) :: filename_loc
INTEGER :: iu_therm, idum, itemp, na, find_free_unit
LOGICAL :: check_file_exists
CHARACTER(LEN=6) :: int_to_char

do_read=.FALSE.
IF ( check_file_exists(filetherm) ) do_read=.TRUE.
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
                    ph_ener(itemp), ph_free_ener(itemp), &
                    ph_entropy(itemp), ph_ce(itemp)
      END DO
      CLOSE(UNIT=iu_therm, STATUS='KEEP')
   END IF
   CALL mp_bcast(ph_ener, meta_ionode_id, world_comm)
   CALL mp_bcast(ph_free_ener, meta_ionode_id, world_comm)
   CALL mp_bcast(ph_entropy, meta_ionode_id, world_comm)
   CALL mp_bcast(ph_ce, meta_ionode_id, world_comm)

END IF

RETURN
END SUBROUTINE read_thermo

!---------------------------------------------------------------------------
SUBROUTINE read_b_factor(ntemp, b_fact, filetherm)
!---------------------------------------------------------------------------
!
!   Reads the b factor and broacast to all processors
!
USE kinds,      ONLY : DP
USE ions_base,  ONLY : nat
USE io_global,  ONLY : meta_ionode, meta_ionode_id
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast

IMPLICIT NONE
INTEGER, INTENT(IN)     :: ntemp
REAL(DP), INTENT(INOUT) :: b_fact(3,3,nat,ntemp)
CHARACTER(LEN=*)        :: filetherm

INTEGER :: na, itemp, ipol, jpol, iu_therm, find_free_unit
CHARACTER(LEN=256) :: filename_loc
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: rdum

DO na=1,nat
   filename_loc=TRIM(filetherm)//'.'//TRIM(int_to_char(na))//'.dw'
   IF (meta_ionode) THEN
      iu_therm=find_free_unit()
      OPEN (UNIT=iu_therm, FILE=TRIM(filename_loc), STATUS='unknown',&
                                             FORM='formatted')
      READ(iu_therm, *)
      DO itemp = 1, ntemp
         READ(iu_therm, '(e16.8,6e18.8)') rdum, &
                     ((b_fact(ipol,jpol,na,itemp), jpol=ipol,3), ipol=1,3)
      END DO
      CLOSE(UNIT=iu_therm, STATUS='KEEP')
   ENDIF
ENDDO
CALL mp_bcast(b_fact, meta_ionode_id, world_comm)

RETURN
END SUBROUTINE read_b_factor

!
!  Copyright (C) 2018 Cristiano Malica
!
!-----------------------------------------------------------------------
SUBROUTINE write_dw_debye(ntemp, temp, deb_bfact, filename)
!-----------------------------------------------------------------------
!
! This routine writes on file the atomic b factor computed using the
! debye model.
!
USE kinds, ONLY : DP

IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), deb_bfact(ntemp)
CHARACTER(LEN=*) :: filename

INTEGER :: iu_therm, itemp
INTEGER :: find_free_unit

iu_therm=find_free_unit()

OPEN (UNIT=iu_therm, FILE=TRIM(filename)//'.dw', STATUS='unknown', FORM='formatted')
WRITE(iu_therm,'("#",6x,"T",14x,"B_33")')

DO itemp = 1, ntemp
   WRITE(iu_therm, '(e16.8,6e18.8)') temp(itemp), deb_bfact(itemp) 
END DO
CLOSE(UNIT=iu_therm, STATUS='KEEP')

RETURN
END SUBROUTINE write_dw_debye

