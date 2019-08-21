!
! Copyright (C) 2019 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_elastic_constants_t()
!
!  This routine sets the temperature dependent elastic constants
!  depending on what has been found in the files and the user requests.
!  If possible it will use the quasi-harmonic elastic constants found
!  in the anhar_files directory (quasi-harmonic approximation), otherwise 
!  it will use the elastic constants at several geometries found in the 
!  elastic_constants directory (quasi-static approximation). 
!
USE kinds, ONLY : DP
USE control_thermo,        ONLY : ltherm_dos, ltherm_freq
USE control_grun,          ONLY : lb0_t
USE control_mur,           ONLY : b0
USE anharmonic,            ONLY : el_cons_t, el_comp_t, b0_t, lelastic
USE ph_freq_anharmonic,    ONLY : el_consf_t, el_compf_t, b0f_t, lelasticf
USE control_macro_elasticity, ONLY: macro_el
USE temperature,           ONLY : ntemp
USE elastic_constants,     ONLY : el_con, el_compliances
USE control_elastic_constants, ONLY : el_cons_available, el_cons_t_available, &
                                  el_cons_qha_available,                 &
                                  el_cons_qha_geo_available,             &
                                  el_consf_qha_available,                &
                                  el_consf_qha_geo_available
IMPLICIT NONE

INTEGER :: itemp

IF (lb0_t) THEN
   IF (el_cons_qha_geo_available.OR.el_consf_qha_geo_available) THEN
      CALL write_elastic_t_qha()
   ELSEIF (el_cons_t_available) THEN
      CALL write_elastic_t()
   ELSEIF(el_cons_qha_available.OR.el_consf_qha_available) THEN
!
!  In this case do nothing because el_cons_t and el_comp_t or
!  el_consf_t and el_compf_t have been already set 
!
   ELSE
      CALL errore('manage_anhar_anis','Temperature dependent elastic &
                   &constants not available',-1)
   ENDIF
ELSEIF(el_cons_qha_available.OR.el_consf_qha_available) THEN
   IF (ltherm_dos.AND.el_cons_qha_available) THEN
      DO itemp=1,ntemp
         el_cons_t(:,:,itemp)=el_cons_t(:,:,2)
         el_comp_t(:,:,itemp)=el_comp_t(:,:,2)
      ENDDO
      lelastic=.TRUE.
   ENDIF
   IF (ltherm_freq.AND.el_consf_qha_available) THEN
      DO itemp=1,ntemp
         el_consf_t(:,:,itemp)=el_consf_t(:,:,2)
         el_compf_t(:,:,itemp)=el_compf_t(:,:,2)
      ENDDO
      lelasticf=.TRUE.
   ENDIF
ELSEIF(el_cons_available) THEN
   b0=macro_el(5)
   b0_t=macro_el(5)
   b0f_t=macro_el(5)
   DO itemp=1,ntemp
      el_cons_t(:,:,itemp)=el_con(:,:)
      el_comp_t(:,:,itemp)=el_compliances(:,:)
      el_consf_t(:,:,itemp)=el_con(:,:)
      el_compf_t(:,:,itemp)=el_compliances(:,:)
   ENDDO
   lelastic=.TRUE.
   lelasticf=.TRUE.
ENDIF

RETURN
END SUBROUTINE set_elastic_constants_t
