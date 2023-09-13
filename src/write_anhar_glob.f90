!
! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_t()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus and
!  its derivative as a function of temperature. It uses the free energy 
!  computed by phonon dos.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, central_geo, celldm_geo, omega_geo, &
                           energy_geo
USE control_ev,     ONLY : npt, v0, e0
USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, b02_t, free_e_min_t, celldm_t
USE thermodynamics, ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE control_eldos,  ONLY : lel_free_energy
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER ::  itemp, startt, lastt, ipt

vmin_t=0.0_DP
b0_t=0.0_DP
b01_t=0.0_DP
b02_t=0.0_DP
free_e_min_t=0.0_DP
celldm_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itemp=startt,lastt
   DO ipt=1,npt
      e0(ipt)=energy_geo(ipt) + pressure * v0(ipt) + ph_free_ener(itemp,ipt)
      IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp,ipt)
   ENDDO
   CALL ev_sub_nodisk(vmin_t(itemp), b0_t(itemp), b01_t(itemp), &
                                  b02_t(itemp), free_e_min_t(itemp) )
   CALL compute_celldm_geo(vmin_t(itemp), celldm_t(1,itemp), &
                   celldm_geo(1,central_geo), omega_geo(central_geo))
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vmin_t, world_comm)
CALL mp_sum(b0_t,   world_comm)
CALL mp_sum(b01_t,  world_comm)
CALL mp_sum(b02_t,  world_comm)
CALL mp_sum(free_e_min_t, world_comm)
CALL mp_sum(celldm_t, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_t
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_noe_t()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus and
!  its derivative as a function of temperature. It uses the free energy 
!  computed by phonon dos. The electronic part is not added.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, central_geo, omega_geo, energy_geo, celldm_geo
USE control_ev,     ONLY : npt, v0, e0
USE anharmonic,     ONLY : vmin_noe_t, b0_noe_t, b01_noe_t, b02_noe_t, &
                           free_e_min_noe_t, a_noe_t, celldm_noe_t
USE thermodynamics, ONLY : ph_free_ener
USE control_eldos,  ONLY : lel_free_energy
USE control_pressure, ONLY : pressure
USE temperature,    ONLY : ntemp
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER ::  itemp, ipt, startt, lastt

vmin_noe_t=0.0_DP
b0_noe_t=0.0_DP
b01_noe_t=0.0_DP
b02_noe_t=0.0_DP
free_e_min_noe_t=0.0_DP
celldm_noe_t=0.0_DP

IF (.NOT.lel_free_energy) RETURN

CALL divide(world_comm, ntemp, startt, lastt)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itemp=startt,lastt
   DO ipt=1,npt
      e0(ipt)=energy_geo(ipt) + pressure * v0(ipt) + ph_free_ener(itemp,ipt)
   ENDDO
   CALL ev_sub_nodisk(vmin_noe_t(itemp), b0_noe_t(itemp), b01_noe_t(itemp), &
                               b02_noe_t(itemp), free_e_min_noe_t(itemp) )
   CALL compute_celldm_geo(vmin_noe_t(itemp), celldm_noe_t(1,itemp), &
                   celldm_geo(1,central_geo), omega_geo(central_geo))
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vmin_noe_t, world_comm)
CALL mp_sum(b0_noe_t,   world_comm)
CALL mp_sum(b01_noe_t,  world_comm)
CALL mp_sum(b02_noe_t,  world_comm)
CALL mp_sum(free_e_min_noe_t, world_comm)
CALL mp_sum(celldm_noe_t, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_noe_t

!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_t_ph()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus and
!  its derivative as a function of temperature. It uses the free energy 
!  computed by the Brillouin zone integration.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, central_geo, omega_geo, energy_geo, celldm_geo
USE control_ev,     ONLY : npt, v0, e0
USE ph_freq_anharmonic,  ONLY : vminf_t, b0f_t, b01f_t, b02f_t, &
                           free_e_minf_t, celldmf_t
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE control_eldos,  ONLY : lel_free_energy
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER ::  itemp, startt, lastt, ipt

vminf_t=0.0_DP
b0f_t=0.0_DP
b01f_t=0.0_DP
b02f_t=0.0_DP
free_e_minf_t=0.0_DP
celldmf_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itemp=startt,lastt
   DO ipt=1,npt
      e0(ipt)=energy_geo(ipt) + pressure * v0(ipt) + phf_free_ener(itemp,ipt)
      IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp,ipt)
   ENDDO
   CALL ev_sub_nodisk(vminf_t(itemp), b0f_t(itemp), b01f_t(itemp), &
                                  b02f_t(itemp), free_e_minf_t(itemp) )
   CALL compute_celldm_geo(vminf_t(itemp), celldmf_t(1,itemp), &
                   celldm_geo(1,central_geo), omega_geo(central_geo))
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vminf_t, world_comm)
CALL mp_sum(b0f_t,   world_comm)
CALL mp_sum(b01f_t,  world_comm)
CALL mp_sum(b02f_t,  world_comm)
CALL mp_sum(free_e_minf_t, world_comm)
CALL mp_sum(celldmf_t, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_t_ph
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_noe_t_ph()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus and
!  its derivative as a function of temperature. It uses the free energy 
!  computed by phonon dos. The electronic part is not added.
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, central_geo, omega_geo, energy_geo, celldm_geo
USE control_ev,     ONLY : npt, v0, e0
USE ph_freq_anharmonic,  ONLY : vminf_noe_t, b0f_noe_t, b01f_noe_t, &
                           b02f_noe_t, free_e_minf_noe_t, af_noe_t, &
                           celldmf_noe_t
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE control_eldos,  ONLY : lel_free_energy
USE control_pressure, ONLY : pressure
USE temperature,    ONLY : ntemp
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER ::  itemp, ipt, startt, lastt

vminf_noe_t=0.0_DP
b0f_noe_t=0.0_DP
b01f_noe_t=0.0_DP
b02f_noe_t=0.0_DP
free_e_minf_noe_t=0.0_DP
celldmf_noe_t=0.0_DP

IF (.NOT.lel_free_energy) RETURN

CALL divide(world_comm, ntemp, startt, lastt)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itemp=startt,lastt
   DO ipt=1,npt
      e0(ipt)=energy_geo(ipt) + pressure * v0(ipt) + phf_free_ener(itemp,ipt)
   ENDDO
   CALL ev_sub_nodisk(vminf_noe_t(itemp), b0f_noe_t(itemp), b01f_noe_t(itemp),&
                               b02f_noe_t(itemp), free_e_minf_noe_t(itemp) )
   CALL compute_celldm_geo(vminf_noe_t(itemp), celldmf_noe_t(1,itemp), &
                   celldm_geo(1,central_geo), omega_geo(central_geo))
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vminf_noe_t, world_comm)
CALL mp_sum(b0f_noe_t,   world_comm)
CALL mp_sum(b01f_noe_t,  world_comm)
CALL mp_sum(b02f_noe_t,  world_comm)
CALL mp_sum(free_e_minf_noe_t, world_comm)
CALL mp_sum(celldmf_noe_t, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_noe_t_ph
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_pt()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus, and its
!  derivatives as a function of temperature for the npress_plot pressures 
!  that we want to plot in output.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, central_geo, omega_geo, energy_geo, &
                           celldm_geo
USE control_pressure, ONLY : npress, press, npress_plot, ipress_plot
USE control_ev,     ONLY : npt, v0, e0
USE temperature,    ONLY : ntemp
USE thermodynamics, ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE anharmonic_pt,  ONLY : vmin_pt, b0_pt, b01_pt, b02_pt, emin_pt, &
                           celldm_pt
USE control_eldos,  ONLY : lel_free_energy
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER  :: itemp, ipt, ipress, ipressp, startt, lastt
REAL(DP) :: vm

IF (npress_plot==0) RETURN

vmin_pt=0.0_DP
b0_pt=0.0_DP
b01_pt=0.0_DP
b02_pt=0.0_DP
emin_pt=0.0_DP
celldm_pt=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO ipressp=1, npress_plot
   ipress=ipress_plot(ipressp)
   DO itemp=startt,lastt
      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                                   ph_free_ener(itemp,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp,ipt)
      ENDDO
      CALL ev_sub_nodisk(vmin_pt(itemp,ipressp), b0_pt(itemp,ipressp), &
       b01_pt(itemp,ipressp), b02_pt(itemp,ipressp), emin_pt(itemp,ipressp) )
      CALL compute_celldm_geo(vmin_pt(itemp,ipressp), &
            celldm_pt(1,itemp,ipressp), celldm_geo(1,central_geo), &
                                        omega_geo(central_geo))
   ENDDO
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vmin_pt, world_comm)
CALL mp_sum(b0_pt,   world_comm)
CALL mp_sum(b01_pt,  world_comm)
CALL mp_sum(b02_pt,  world_comm)
CALL mp_sum(emin_pt, world_comm)
CALL mp_sum(celldm_pt, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_pt
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_ph_pt()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus, and its
!  derivatives as a function of temperature for the npress_plot pressures 
!  that we want to plot in output.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, central_geo, omega_geo, energy_geo, celldm_geo
USE control_pressure, ONLY : npress, press, npress_plot, ipress_plot
USE control_ev,     ONLY : npt, v0, e0
USE temperature,    ONLY : ntemp
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE ph_freq_anharmonic_pt,  ONLY : vminf_pt, b0f_pt, b01f_pt, b02f_pt, &
                           eminf_pt, celldmf_pt
USE control_eldos,  ONLY : lel_free_energy
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE
INTEGER  :: itemp, ipt, ipress, ipressp, startt, lastt
REAL(DP) :: vm

IF (npress_plot==0) RETURN

vminf_pt=0.0_DP
b0f_pt=0.0_DP
b01f_pt=0.0_DP
b02f_pt=0.0_DP
eminf_pt=0.0_DP
celldmf_pt=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO ipressp=1, npress_plot
   ipress=ipress_plot(ipressp)
   DO itemp=startt,lastt
      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                                   phf_free_ener(itemp,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp,ipt)
      ENDDO
      CALL ev_sub_nodisk(vminf_pt(itemp,ipressp), b0f_pt(itemp,ipressp), &
       b01f_pt(itemp,ipressp), b02f_pt(itemp,ipressp), eminf_pt(itemp,ipressp))
      CALL compute_celldm_geo(vminf_pt(itemp,ipressp), &
            celldmf_pt(1,itemp,ipressp), celldm_geo(1,central_geo), &
                                        omega_geo(central_geo))
   ENDDO
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vminf_pt, world_comm)
CALL mp_sum(b0f_pt,   world_comm)
CALL mp_sum(b01f_pt,  world_comm)
CALL mp_sum(b02f_pt,  world_comm)
CALL mp_sum(eminf_pt, world_comm)
CALL mp_sum(celldmf_pt, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_ph_pt
!
!-------------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_ptt()
!-------------------------------------------------------------------------
!
!  This routine computes the volume thermal expansion, the bulk modulus,
!  the heat capacity and the average Gruneisen parameter as a function 
!  of pressure for selected temperatures.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, central_geo, omega_geo, energy_geo, &
                           celldm_geo
USE control_ev,     ONLY : npt, e0, v0
USE thermodynamics, ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE control_eldos,  ONLY : lel_free_energy
USE anharmonic_ptt, ONLY : vmin_ptt, emin_ptt, b0_ptt, b01_ptt, b02_ptt, &
                           celldm_ptt
USE temperature,    ONLY : ntemp_plot, itemp_plot
USE control_pressure,  ONLY : press, npress
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE

REAL(DP) :: vm
INTEGER ::  ipt, ipress, itempp, itemp, startp, lastp

IF (ntemp_plot==0) RETURN

vmin_ptt=0.0_DP
b0_ptt=0.0_DP
b01_ptt=0.0_DP
b02_ptt=0.0_DP
emin_ptt=0.0_DP
celldm_ptt=0.0_DP
CALL divide(world_comm, npress, startp, lastp)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   DO ipress=startp,lastp
      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                                   ph_free_ener(itemp,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp,ipt)
      ENDDO
      CALL ev_sub_nodisk(vmin_ptt(ipress,itempp), b0_ptt(ipress,itempp), &
     b01_ptt(ipress,itempp), b02_ptt(ipress,itempp), emin_ptt(ipress,itempp) )
      CALL compute_celldm_geo(vmin_ptt(ipress,itempp), &
            celldm_ptt(1,ipress,itempp), celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
   ENDDO
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vmin_ptt, world_comm)
CALL mp_sum(b0_ptt, world_comm)
CALL mp_sum(b01_ptt, world_comm)
CALL mp_sum(b02_ptt, world_comm)
CALL mp_sum(emin_ptt, world_comm)
CALL mp_sum(celldm_ptt, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_ptt
!
!-------------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_ph_ptt()
!-------------------------------------------------------------------------
!
!  This routine computes the volume thermal expansion, the bulk modulus,
!  the heat capacity and the average Gruneisen parameter as a function 
!  of pressure for selected temperatures.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, central_geo, omega_geo, energy_geo, &
                           celldm_geo
USE control_ev,     ONLY : npt, e0, v0
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE control_eldos,  ONLY : lel_free_energy
USE ph_freq_anharmonic_ptt, ONLY : vminf_ptt, eminf_ptt, b0f_ptt, &
                           b01f_ptt, b02f_ptt, celldmf_ptt
USE temperature,    ONLY : ntemp_plot, itemp_plot
USE control_pressure,  ONLY : press, npress
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE

REAL(DP) :: vm
INTEGER ::  ipt, ipress, itempp, itemp, startp, lastp

IF (ntemp_plot==0) RETURN

vminf_ptt=0.0_DP
b0f_ptt=0.0_DP
b01f_ptt=0.0_DP
b02f_ptt=0.0_DP
eminf_ptt=0.0_DP
celldmf_ptt=0.0_DP
CALL divide(world_comm, npress, startp, lastp)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   DO ipress=startp,lastp
      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                                   phf_free_ener(itemp,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp,ipt)
      ENDDO
      CALL ev_sub_nodisk(vminf_ptt(ipress,itempp), b0f_ptt(ipress,itempp),&
     b01f_ptt(ipress,itempp), b02f_ptt(ipress,itempp), &
                              eminf_ptt(ipress,itempp) )
      CALL compute_celldm_geo(vminf_ptt(ipress,itempp), &
            celldmf_ptt(1,ipress,itempp), celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
   ENDDO
ENDDO
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_sum(vminf_ptt, world_comm)
CALL mp_sum(b0f_ptt, world_comm)
CALL mp_sum(b01f_ptt, world_comm)
CALL mp_sum(b02f_ptt, world_comm)
CALL mp_sum(eminf_ptt, world_comm)
CALL mp_sum(celldmf_ptt, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_ph_ptt
!
!-----------------------------------------------------------------------
SUBROUTINE anhar_ev_glob_vt()
!-----------------------------------------------------------------------
!
!  This subroutine computes the thermal pressure as a function of 
!  temperature for the nvol_plot volumes required in input.
!  It requires p0, the pressure at T=0 K as a function of the volume
!  and the parameters of the equation of state that interpolates the
!  free energy at each temperature.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : omega_geo
USE anharmonic_vt,    ONLY : press_vtt
USE control_mur,      ONLY : p0
USE control_ev,       ONLY : ieos
USE temperature,      ONLY : ntemp
USE control_vol,      ONLY : nvol_plot, ivol_plot
USE anharmonic,       ONLY : vmin_t, b0_t, b01_t, b02_t
USE eos,              ONLY : eos_press
USE io_global,        ONLY : meta_ionode_id
USE mp,               ONLY : mp_sum, mp_bcast
USE mp_world,         ONLY : world_comm

IMPLICIT NONE
REAL(DP) :: press, omega, p00
INTEGER  :: itemp, ivolp, igeo, startt, lastt

IF (nvol_plot==0) RETURN
!
!  p0 is known only to the meta_ionode that has written the file,
!  we now need it in all nodes
!
CALL mp_bcast(p0, meta_ionode_id, world_comm)
press_vtt=0.0_DP
CALL divide(world_comm, ntemp, startt, lastt)
DO ivolp=1,nvol_plot
   igeo = ivol_plot(ivolp)
   omega = omega_geo(igeo)
   CALL interpolate_p0(p00,omega)
   DO itemp=startt,lastt
      CALL eos_press(ieos, omega, press, vmin_t(itemp), b0_t(itemp)/ry_kbar, &
                     b01_t(itemp), b02_t(itemp)*ry_kbar)
      press_vtt(itemp,ivolp) = press * ry_kbar - p00
   ENDDO
ENDDO
CALL mp_sum(press_vtt, world_comm)

RETURN
END SUBROUTINE anhar_ev_glob_vt
!
!-----------------------------------------------------------------------
SUBROUTINE ph_freq_anhar_ev_glob_vt()
!-----------------------------------------------------------------------
!
!  This subroutine computes the thermal pressure as a function of 
!  temperature for the nvol_plot volumes required in input.
!  It requires p0, the pressure at T=0 K as a function of the volume
!  and the parameters of the equation of state that interpolates the
!  free energy at each temperature.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE thermo_mod,       ONLY : omega_geo
USE ph_freq_anharmonic_vt, ONLY : pressf_vtt
USE control_mur,      ONLY : p0
USE control_ev,       ONLY : ieos
USE temperature,      ONLY : ntemp
USE control_vol,      ONLY : nvol_plot, ivol_plot
USE ph_freq_anharmonic, ONLY : vminf_t, b0f_t, b01f_t, b02f_t
USE eos,              ONLY : eos_press
USE io_global,        ONLY : meta_ionode_id
USE mp,               ONLY : mp_sum, mp_bcast
USE mp_world,         ONLY : world_comm

IMPLICIT NONE
REAL(DP) :: press, omega, p00
INTEGER  :: itemp, ivolp, igeo, startt, lastt

IF (nvol_plot==0) RETURN
!
!  p0 is known only to the meta_ionode that has written the file,
!  we now need it in all nodes
!
CALL mp_bcast(p0, meta_ionode_id, world_comm)
pressf_vtt=0.0_DP
CALL divide(world_comm, ntemp, startt, lastt)
DO ivolp=1,nvol_plot
   igeo = ivol_plot(ivolp)
   omega = omega_geo(igeo)
   CALL interpolate_p0(p00,omega)
   DO itemp=startt,lastt
      CALL eos_press(ieos, omega, press, vminf_t(itemp), b0f_t(itemp)/ry_kbar,&
                     b01f_t(itemp), b02f_t(itemp)*ry_kbar)
      pressf_vtt(itemp,ivolp) = press * ry_kbar - p00
   ENDDO
ENDDO
CALL mp_sum(pressf_vtt, world_comm)

RETURN
END SUBROUTINE ph_freq_anhar_ev_glob_vt
!
!----------------------------------------------------------------------
SUBROUTINE write_mur_pol_glob(omega0, b0, b01, b02, emin, itempp, filename)
!----------------------------------------------------------------------
!
!  This routine writes on file the energy versus volume
!  curve, together with the pressure versus volume curve. Depending
!  on which equation of state has been used to interpolate the
!  T=0 K energy, it calls the appropriate routine.
!  It receives the parameters of the equation of state:
!
! in input emin in Ry, omega0 in (a.u.)**3, b0 in kbar, b01 adimensional
! b02 in 1/kbar
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE data_files,       ONLY : flanhar
USE thermo_mod,       ONLY : ngeo, omega_geo
USE control_mur,      ONLY : p0
USE control_vol,      ONLY : nvol, vmin_input, vmax_input, deltav
USE control_ev,       ONLY : ieos
USE anharmonic_vt,    ONLY : press_vt
USE polyfit_mod,      ONLY : compute_poly, compute_poly_deriv
USE eos,              ONLY : eos_energy, eos_press
USE control_pressure, ONLY : pressure_kb
USE temperature,      ONLY : temp, itemp_plot
USE io_global,        ONLY : meta_ionode

IMPLICIT NONE

REAL(DP), INTENT(IN) :: emin, omega0, b0, b01, b02
INTEGER, INTENT(IN) :: itempp
CHARACTER(LEN=256)   :: filename
CHARACTER(LEN=8) :: float_to_char
REAL(DP) :: omega, free, pres, e, p
INTEGER  :: i, j, iu_mur, itemp
INTEGER  :: find_free_unit

itemp=itemp_plot(itempp)
IF (itemp<1) RETURN

CALL add_value(filename,temp(itemp))
CALL add_pressure(filename)

IF (meta_ionode) THEN
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   CALL write_mur_start_line(itemp, iu_mur)
   DO i=1,nvol
      omega= vmin_input + deltav * (i-1)
      CALL eos_energy(ieos, omega, e, omega0, b0/ry_kbar, b01, &
                                                  b02*ry_kbar)
      e=e+emin
      CALL eos_press(ieos, omega, p, omega0, b0/ry_kbar, b01, &
                                                 b02*ry_kbar)
      press_vt(i,itempp)=p*ry_kbar + pressure_kb
      WRITE(iu_mur,'(f18.10,4f20.10)') omega, e, e+p*omega, &
                                        press_vt(i,itempp), p0(i)
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_mur_pol_glob
!
!
!-------------------------------------------------------------------------
SUBROUTINE write_anhar_glob_ptt()
!-------------------------------------------------------------------------
!
!  This routine computes the volume thermal expansion, the bulk modulus,
!  the heat capacity and the average Gruneisen parameter as a function 
!  of pressure for selected temperatures.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo
USE control_ev,     ONLY : npt, v0, e0
USE data_files,     ONLY : flanhar
USE temperature,    ONLY : temp, deltat, itemp_plot, ntemp_plot, itemp300
USE thermodynamics, ONLY : ph_ce
USE anharmonic,     ONLY : vmin_t, a_t, b0_t, b01_t, b02_t, free_e_min_t, &
                           noelcvg
USE anharmonic_ptt, ONLY : beta_ptt, vmin_ptt, b0_ptt, b01_ptt, b02_ptt, &
                           gamma_ptt, ce_ptt, cp_ptt, b0_s_ptt, emin_ptt
USE el_anharmonic,  ONLY : el_ce_ptt
USE anharmonic_vt,  ONLY : press_vt
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_pressure,  ONLY : press, npress
USE control_eldos,  ONLY : lel_free_energy
USE el_thermodynamics, ONLY : el_free_ener
USE thermodynamics, ONLY : ph_free_ener
USE control_vol,    ONLY : vmin_input, deltav, nvol
USE io_global,      ONLY : meta_ionode
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE

CHARACTER(LEN=256) :: filename
LOGICAL :: subtract_el
REAL(DP) :: vmm1, vmp1, omega, bulk
REAL(DP) :: b0m1, b01m1, b02m1, b0p1, b01p1, b02p1, enem1, enep1
REAL(DP), ALLOCATABLE :: aux(:,:)
INTEGER :: ivol, i1, m1, ipt, iu_mur, ipress, itemp, itempp, startp, lastp
INTEGER :: find_free_unit

IF (ntemp_plot==0) RETURN

m1=poly_degree_ph+1

ALLOCATE(aux(npress,ntemp_plot))
press_vt=0.0_DP
beta_ptt=0.0_DP
cp_ptt=0.0_DP
b0_s_ptt=0.0_DP
gamma_ptt=0.0_DP

CALL divide(world_comm, npress, startp, lastp)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   filename="anhar_files/"//TRIM(flanhar)//'.mur_temp'
   CALL write_mur_pol_glob(vmin_t(itemp), b0_t(itemp), b01_t(itemp), &
           b02_t(itemp), free_e_min_t(itemp), itempp, filename)
!   CALL write_thermal_press(a_t(:,itemp), m1, itempp)

   DO ipress=startp,lastp
      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                   ph_free_ener(itemp-1,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp-1,ipt)
      ENDDO
      CALL ev_sub_nodisk(vmm1, b0m1, b01m1, b02m1, enem1 )

      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                   ph_free_ener(itemp+1,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp+1,ipt)
      ENDDO
      CALL ev_sub_nodisk(vmp1, b0p1, b01p1, b02p1, enep1 )
      beta_ptt(ipress, itempp) = (vmp1 - vmm1) / 2.0_DP / deltat /     &
                                        vmin_ptt(ipress,itempp)
   ENDDO 
ENDDO
!CALL mp_sum(press_vt, world_comm)
CALL mp_sum(beta_ptt, world_comm)
DEALLOCATE(v0)
DEALLOCATE(e0)

subtract_el=(lel_free_energy.AND.noelcvg)
DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   CALL compute_cp_bs_gp(beta_ptt(:,itempp), vmin_ptt(:,itempp),      &
              b0_ptt(:,itempp), ce_ptt(:,itempp), cp_ptt(:,itempp),   &
              b0_s_ptt(:,itempp), gamma_ptt(:,itempp),                &
              el_ce_ptt(:,itempp),itemp,subtract_el)

   IF (itemp300 > 0) THEN
      aux(:,itempp)=vmin_ptt(:,itempp)/vmin_t(itemp300)
   ELSE
      aux(:,itempp)=vmin_ptt(:,itempp)
   ENDIF
ENDDO

IF (meta_ionode) THEN
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      iu_mur=find_free_unit()
      filename="anhar_files/"//TRIM(flanhar)//'.temp'
      CALL add_value(filename, temp(itemp))
      CALL write_ener_beta_ptt(press, vmin_ptt(:,itempp), aux(:,itempp), & 
             emin_ptt(:,itempp), beta_ptt(:,itempp), npress, itemp, filename)

      filename="anhar_files/"//TRIM(flanhar)//'.bulk_temp'
      CALL add_value(filename, temp(itemp))
      CALL write_bulk_anharm_ptt(press, b0_ptt(:,itempp), b0_s_ptt(:,itempp), &
                                                     npress, itemp, filename)

      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_temp'
      CALL add_value(filename, temp(itemp))
      CALL write_dbulk_anharm_t(press, b01_ptt(:,itempp), b02_ptt(:,itempp),  &
                                                     npress, itemp, filename)
!
!   here the heat capacities (constant strain and constant stress)
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat_temp'
      CALL add_value(filename, temp(itemp))

      CALL write_heat_anhar_t(press, ce_ptt(:,itempp), ce_ptt(:,itempp), &
              cp_ptt(:,itempp), npress, itemp, filename)

      filename="anhar_files/"//TRIM(flanhar)//'.gamma_temp'
      CALL add_value(filename, temp(itemp))

      CALL write_gamma_anharm_t(press, gamma_ptt(:,itempp), &
            ce_ptt(:,itempp), beta_ptt(:,itempp), b0_ptt(:,itempp), &
                                           npress, itemp, filename)

   ENDDO
ENDIF

DEALLOCATE(aux)

RETURN
END SUBROUTINE write_anhar_glob_ptt

!-------------------------------------------------------------------------
SUBROUTINE write_ph_freq_anhar_glob_ptt()
!-------------------------------------------------------------------------
!
!  This routine computes the volume thermal expansion, the bulk modulus,
!  the heat capacity and the average Gruneisen parameter as a function 
!  of pressure for selected temperatures.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo
USE control_ev,     ONLY : npt, v0, e0
USE data_files,     ONLY : flanhar
USE temperature,    ONLY : temp, deltat, itemp_plot, ntemp_plot, itemp300
USE thermodynamics, ONLY : ph_ce
USE ph_freq_anharmonic,     ONLY : vminf_t, af_t, b0f_t, b01f_t, &
                            b02f_t, free_e_minf_t
USE anharmonic,             ONLY : noelcvg
USE ph_freq_anharmonic_ptt, ONLY : betaf_ptt, vminf_ptt, b0f_ptt, b01f_ptt, &
                           b02f_ptt, gammaf_ptt, cef_ptt, cpf_ptt, b0f_s_ptt, &
                           eminf_ptt
USE el_anharmonic,  ONLY : el_cef_ptt
USE ph_freq_anharmonic_vt,  ONLY : pressf_vt
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_pressure,  ONLY : press, npress
USE control_eldos,  ONLY : lel_free_energy
USE el_thermodynamics, ONLY : el_free_ener
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE control_vol,    ONLY : vmin_input, deltav, nvol
USE io_global,      ONLY : meta_ionode
USE mp,             ONLY : mp_sum
USE mp_world,       ONLY : world_comm

IMPLICIT NONE

CHARACTER(LEN=256) :: filename
LOGICAL :: subtract_el
REAL(DP) :: vmm1, vmp1, omega, bulk
REAL(DP) :: b0m1, b01m1, b02m1, b0p1, b01p1, b02p1, enem1, enep1
REAL(DP), ALLOCATABLE :: aux(:,:)
INTEGER :: ivol, i1, m1, ipt, iu_mur, ipress, itemp, itempp, startp, lastp
INTEGER :: find_free_unit

IF (ntemp_plot==0) RETURN

m1=poly_degree_ph+1

ALLOCATE(aux(npress,ntemp_plot))
pressf_vt=0.0_DP
betaf_ptt=0.0_DP
cpf_ptt=0.0_DP
b0f_s_ptt=0.0_DP
gammaf_ptt=0.0_DP

CALL divide(world_comm, npress, startp, lastp)

npt=ngeo(1)
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_geo(ipt)
ENDDO

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   filename="anhar_files/"//TRIM(flanhar)//'.mur_ph_temp'
   CALL write_mur_pol_glob(vminf_t(itemp), b0f_t(itemp), b01f_t(itemp), &
           b02f_t(itemp), free_e_minf_t(itemp), itempp, filename)
!   CALL write_thermal_press(a_t(:,itemp), m1, itempp)

   DO ipress=startp,lastp
      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                   phf_free_ener(itemp-1,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp-1,ipt)
      ENDDO
      CALL ev_sub_nodisk(vmm1, b0m1, b01m1, b02m1, enem1 )

      DO ipt=1,npt
         e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar + &
                                   phf_free_ener(itemp+1,ipt)
         IF (lel_free_energy) e0(ipt)=e0(ipt)+el_free_ener(itemp+1,ipt)
      ENDDO
      CALL ev_sub_nodisk(vmp1, b0p1, b01p1, b02p1, enep1 )
      betaf_ptt(ipress, itempp) = (vmp1 - vmm1) / 2.0_DP / deltat /     &
                                        vminf_ptt(ipress,itempp)
   ENDDO 
ENDDO
!CALL mp_sum(press_vt, world_comm)
CALL mp_sum(betaf_ptt, world_comm)
DEALLOCATE(v0)
DEALLOCATE(e0)

subtract_el=(lel_free_energy.AND.noelcvg)
DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   CALL compute_cp_bs_gp(betaf_ptt(:,itempp), vminf_ptt(:,itempp),      &
              b0f_ptt(:,itempp), cef_ptt(:,itempp), cpf_ptt(:,itempp),   &
              b0f_s_ptt(:,itempp), gammaf_ptt(:,itempp),                &
              el_cef_ptt(:,itempp),itemp,subtract_el)

   IF (itemp300 > 0) THEN
      aux(:,itempp)=vminf_ptt(:,itempp)/vminf_t(itemp300)
   ELSE
      aux(:,itempp)=vminf_ptt(:,itempp)
   ENDIF
ENDDO

IF (meta_ionode) THEN
   DO itempp=1,ntemp_plot
      itemp=itemp_plot(itempp)
      iu_mur=find_free_unit()
      filename="anhar_files/"//TRIM(flanhar)//'.ph_temp'
      CALL add_value(filename, temp(itemp))
      CALL write_ener_beta_ptt(press, vminf_ptt(:,itempp), aux(:,itempp), & 
             eminf_ptt(:,itempp), betaf_ptt(:,itempp), npress, itemp, filename)

      filename="anhar_files/"//TRIM(flanhar)//'.bulk_ph_temp'
      CALL add_value(filename, temp(itemp))
      CALL write_bulk_anharm_ptt(press, b0f_ptt(:,itempp), &
                    b0f_s_ptt(:,itempp), npress, itemp, filename)

      filename="anhar_files/"//TRIM(flanhar)//'.dbulk_ph_temp'
      CALL add_value(filename, temp(itemp))
      CALL write_dbulk_anharm_t(press, b01f_ptt(:,itempp), &
                                 b02f_ptt(:,itempp), npress, itemp, filename)
!
!   here the heat capacities (constant strain and constant stress)
!
      filename="anhar_files/"//TRIM(flanhar)//'.heat_ph_temp'
      CALL add_value(filename, temp(itemp))

      CALL write_heat_anhar_t(press, cef_ptt(:,itempp), cef_ptt(:,itempp), &
              cpf_ptt(:,itempp), npress, itemp, filename)

      filename="anhar_files/"//TRIM(flanhar)//'.gamma_ph_temp'
      CALL add_value(filename, temp(itemp))

      CALL write_gamma_anharm_t(press, gammaf_ptt(:,itempp), &
            cef_ptt(:,itempp), betaf_ptt(:,itempp), b0f_ptt(:,itempp), &
                                           npress, itemp, filename)

   ENDDO
ENDIF

DEALLOCATE(aux)

RETURN
END SUBROUTINE write_ph_freq_anhar_glob_ptt
