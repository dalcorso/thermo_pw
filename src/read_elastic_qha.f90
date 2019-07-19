! Copyright (C) 2019 Cristiano Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE read_elastic_qha()

USE kinds,            ONLY : DP
USE cell_base,        ONLY : ibrav
USE thermo_sym,       ONLY : laue
USE temperature,      ONLY : ntemp, temp
USE control_elastic_constants,  ONLY : el_cons_qha_available
USE anharmonic,         ONLY : lelastic, el_cons_t, el_comp_t, b0_t
USE ph_freq_anharmonic, ONLY : lelasticf, el_consf_t, el_compf_t, b0f_t 
USE data_files,       ONLY : flanhar
USE control_thermo,   ONLY : ltherm_freq, ltherm_dos

IMPLICIT NONE

REAL(DP) :: k0_t(ntemp)
CHARACTER(LEN=256) :: filelastic, filelastic_comp
LOGICAL  :: check_file_exists

el_cons_qha_available=.FALSE.

IF (ltherm_freq) THEN
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_ph'
   filelastic_comp='anhar_files/'//TRIM(flanhar)//'.el_comp_ph'
   IF (check_file_exists(filelastic).AND. &
      check_file_exists(filelastic_comp)) THEN
      el_cons_qha_available=.TRUE.
      lelasticf=.TRUE.
      CALL read_el_cons_from_file(temp, ntemp, ibrav, laue, &
                                  el_consf_t, b0f_t, filelastic)
      CALL read_el_cons_from_file(temp, ntemp, ibrav, laue, & 
                                  el_compf_t, k0_t, filelastic_comp)
   ENDIF 

ELSEIF (ltherm_dos) THEN
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
   filelastic_comp='anhar_files/'//TRIM(flanhar)//'.el_comp'
   IF (check_file_exists(filelastic).AND. &
      check_file_exists(filelastic_comp)) THEN
      el_cons_qha_available=.TRUE.
      lelastic=.TRUE.
      CALL read_el_cons_from_file(temp, ntemp, ibrav, laue, &
                                  el_cons_t, b0_t, filelastic)
      CALL read_el_cons_from_file(temp, ntemp, ibrav, laue, & 
                                  el_comp_t, k0_t, filelastic_comp)
   ENDIF   
ENDIF
RETURN
END SUBROUTINE read_elastic_qha
