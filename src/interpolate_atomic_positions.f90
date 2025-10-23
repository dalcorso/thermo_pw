!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE interpolate_atomic_positions()
!----------------------------------------------------------------------
!
USE thermo_mod,         ONLY : uint_geo_eos, tot_ngeo
USE thermodynamics,     ONLY : uint_geo_eos_t
USE anharmonic,         ONLY : celldm_t, uint_t, tau_t, uint_zsisa_t,     &
                               tau_zsisa_t
USE temperature,        ONLY : ntemp_plot, itemp_plot
USE anharmonic_pt,      ONLY : celldm_pt, uint_pt, tau_pt, uint_zsisa_pt, &
                               tau_zsisa_pt
USE anharmonic_ptt,     ONLY : celldm_ptt, uint_ptt, tau_ptt,             &
                               uint_zsisa_ptt, tau_zsisa_ptt,             &
                               uint_ptt_p1, uint_ptt_m1,                  &
                               tau_ptt_p1, tau_ptt_m1
USE control_pressure,   ONLY : npress_plot
USE control_atomic_pos, ONLY : linternal_thermo, linterpolate_tau

IMPLICIT NONE
INTEGER :: ipressp, itempp, itemp, igeo
 
IF (linternal_thermo) THEN
   CALL interpolate_uint_t(celldm_t, uint_geo_eos_t, uint_t, tau_t)
   DO ipressp=1,npress_plot
      CALL interpolate_uint_t(celldm_pt(1,1,ipressp), uint_geo_eos_t, &
                             uint_pt(1,1,ipressp), tau_pt(1,1,1,ipressp))
   ENDDO
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_uint_ptt(celldm_ptt(1,1,itempp), &
                      uint_geo_eos_t, uint_ptt(1,1,itempp), &
                      tau_ptt(1,1,1,itempp), uint_ptt_p1(1,1,itempp), &
                      tau_ptt_p1(1,1,1,itempp), uint_ptt_m1(1,1,itempp), &
                      tau_ptt_m1(1,1,1,itempp), itemp)
   ENDDO
ELSEIF (linterpolate_tau) THEN
   CALL interpolate_uint(celldm_t, uint_geo_eos, &
                                           uint_zsisa_t, tau_zsisa_t)
   DO ipressp=1,npress_plot
      CALL interpolate_uint(celldm_pt(1,1,ipressp), uint_geo_eos, &
                   uint_zsisa_pt(1,1,ipressp), tau_zsisa_pt(1,1,1,ipressp))
   ENDDO
   DO itempp=1,ntemp_plot
      CALL interpolate_uint_zptt(celldm_ptt(1,1,itempp), uint_geo_eos, &
               uint_zsisa_ptt(1,1,itempp), tau_zsisa_ptt(1,1,1,itempp))
   ENDDO
ENDIF

RETURN
END SUBROUTINE interpolate_atomic_positions

!----------------------------------------------------------------------
SUBROUTINE interpolate_atomic_positionsf()
!----------------------------------------------------------------------
!
USE thermo_mod,             ONLY : uint_geo_eos
USE ph_freq_thermodynamics, ONLY : uintf_geo_eos_t
USE ph_freq_anharmonic,     ONLY : celldmf_t, uintf_t, tauf_t,        &
                                   uintf_zsisa_t, tauf_zsisa_t
USE temperature,            ONLY : ntemp_plot, itemp_plot
USE ph_freq_anharmonic_pt,  ONLY : celldmf_pt, uintf_pt, tauf_pt,     &
                                   uintf_zsisa_pt, tauf_zsisa_pt
USE ph_freq_anharmonic_ptt, ONLY : celldmf_ptt, uintf_ptt, tauf_ptt,  &
                                   uintf_zsisa_ptt, tauf_zsisa_ptt,   &
                                   uintf_ptt_p1, uintf_ptt_m1,        &
                                   tauf_ptt_p1, tauf_ptt_m1
USE control_pressure,       ONLY : npress_plot
USE control_atomic_pos,     ONLY : linternal_thermo, linterpolate_tau

IMPLICIT NONE
INTEGER :: ipressp, itempp, itemp

IF (linternal_thermo) THEN 
   CALL interpolate_uint_t(celldmf_t, uintf_geo_eos_t, uintf_t, tauf_t)
   DO ipressp=1,npress_plot
      CALL interpolate_uint_t(celldmf_pt(1,1,ipressp), uintf_geo_eos_t, &
                              uintf_pt(1,1,ipressp), tauf_pt(1,1,1,ipressp))
   ENDDO
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_uint_ptt(celldmf_ptt(1,1,itempp), &
                      uintf_geo_eos_t, uintf_ptt(1,1,itempp), &
                      tauf_ptt(1,1,1,itempp), uintf_ptt_p1(1,1,itempp), &
                      tauf_ptt_p1(1,1,1,itempp), uintf_ptt_m1(1,1,itempp), &
                      tauf_ptt_m1(1,1,1,itempp), itemp)
   ENDDO
ELSEIF (linterpolate_tau) THEN
   CALL interpolate_uint(celldmf_t, uint_geo_eos, &
                                           uintf_zsisa_t, tauf_zsisa_t)
   DO ipressp=1,npress_plot
      CALL interpolate_uint(celldmf_pt(1,1,ipressp), uint_geo_eos, &
                   uintf_zsisa_pt(1,1,ipressp), tauf_zsisa_pt(1,1,1,ipressp))
   ENDDO
   DO itempp=1,ntemp_plot
      CALL interpolate_uint_zptt(celldmf_ptt(1,1,itempp), uint_geo_eos, &
            uintf_zsisa_ptt(1,1,itempp), tauf_zsisa_ptt(1,1,1,itempp))
   ENDDO
ENDIF

RETURN
END SUBROUTINE interpolate_atomic_positionsf
