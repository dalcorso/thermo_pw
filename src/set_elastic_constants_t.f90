!
! Copyright (C) 2019 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE set_elastic_constants_t()
!----------------------------------------------------------------------
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
USE anharmonic,            ONLY : el_cons_t, el_comp_t, b0_t
USE ph_freq_anharmonic,    ONLY : el_consf_t, el_compf_t, b0f_t
USE control_macro_elasticity, ONLY: macro_el
USE temperature,           ONLY : ntemp
USE elastic_constants,     ONLY : el_con, el_compliances
USE control_elastic_constants, ONLY : el_cons_available,                 &
                                  el_cons_geo_available,                 &
                                  el_cons_qha_available,                 &
                                  el_cons_qha_geo_available,             &
                                  el_consf_qha_available,                &
                                  el_consf_qha_geo_available, lelastic,  &
                                  lelasticf, stype
IMPLICIT NONE

INTEGER :: itemp, istep

IF (lb0_t) THEN
   IF (el_cons_qha_geo_available.OR.el_consf_qha_geo_available) THEN
      CALL write_elastic_t_qha()
      CALL write_elastic_pt_qha()
      CALL write_elastic_ptt_qha()
      DO istep=1,21
         IF (stype(istep)) THEN
            CALL write_dyde_t_qha(istep)
         ENDIF
      ENDDO
   ELSEIF (el_cons_geo_available) THEN
      CALL write_elastic_t()
      CALL write_elastic_pt()
      CALL write_elastic_ptt()
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

!----------------------------------------------------------------------
SUBROUTINE set_piezo_tensor_t()
!----------------------------------------------------------------------
!
!  This routine sets the temperature dependent piezoelectric tensor
!  depending on what has been found in the files and the user requests.
!  It will use the piezoelectric tensor at several geometries found in the 
!  elastic_constants directory (quasi-static approximation). 
!
USE control_piezoelectric_tensor, ONLY : piezo_geo_available
IMPLICIT NONE

IF (piezo_geo_available) THEN
   CALL write_piezo_t()
   CALL write_piezo_pt()
   CALL write_piezo_ptt()
ENDIF

RETURN
END SUBROUTINE set_piezo_tensor_t

!----------------------------------------------------------------------
SUBROUTINE set_epsilon_infty_t()
!----------------------------------------------------------------------
!
!  This routine sets the temperature dependent dielectric constants
!  depending on what has been found in the files and the user requests.
!  It will use the dielectric constant at several geometries found in the 
!  dynamical matrix file (quasi-static approximation). 
!
USE thermo_mod, ONLY : tot_ngeo
USE control_epsilon_infty, ONLY : epsilon_infty_geo_available, &
                                  lepsilon_infty_geo
IMPLICIT NONE

INTEGER :: igeom

epsilon_infty_geo_available=.TRUE.
DO igeom=1,tot_ngeo
   epsilon_infty_geo_available=epsilon_infty_geo_available.AND. &
                                         lepsilon_infty_geo(igeom)
ENDDO


IF (epsilon_infty_geo_available) THEN
   CALL write_epsilon_infty_t()
   CALL write_epsilon_infty_pt()
   CALL write_epsilon_infty_ptt()
ENDIF

RETURN
END SUBROUTINE set_epsilon_infty_t

!----------------------------------------------------------------------
SUBROUTINE set_zeu_t()
!----------------------------------------------------------------------
!
!  This routine sets the temperature dependent born effective charge
!  depending on what has been found in the files and the user requests.
!  It will use the born effective charges at several geometries found in the 
!  dynamical matrix file (quasi-static approximation). 
!
USE thermo_mod, ONLY : tot_ngeo
USE control_epsilon_infty, ONLY : zeu_geo_available, lzeu_geo
IMPLICIT NONE

INTEGER :: igeom

zeu_geo_available=.TRUE.
DO igeom=1,tot_ngeo
   zeu_geo_available=zeu_geo_available.AND.lzeu_geo(igeom)
ENDDO

IF (zeu_geo_available) THEN
   CALL write_zeu_t()
   CALL write_zeu_pt()
   CALL write_zeu_ptt()
ENDIF

RETURN
END SUBROUTINE set_zeu_t
