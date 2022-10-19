!
! Copyright (C) 2013-2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------
SUBROUTINE manage_energy_minimum(nwork)
!------------------------------------------------------
!
!  This routine coordinates the work to find the minimum
!  of the total energy as a function of the crystal
!  parameters.
!

USE kinds,            ONLY : DP
USE thermo_mod,       ONLY : celldm_geo, omega_geo, central_geo, density, &
                             energy_geo
USE control_mur,      ONLY : vmin, b0, b01, b02, emin, lmurn
USE control_thermo,   ONLY : lgeo_from_file
USE geometry_file,    ONLY : compute_celldm_geo_file
USE equilibrium_conf, ONLY : celldm0, omega0, tau0, at0, tau0_crys
USE initial_conf,     ONLY : ibrav_save, tau_save_crys
USE control_xrdp,     ONLY : lxrdp
USE control_pressure, ONLY : pressure_kb

USE ions_base,        ONLY : nat, tau
USE cell_base,        ONLY : at, bg, omega, cell_base_init

USE io_global,        ONLY : meta_ionode_id
USE mp_images,        ONLY : nproc_image
USE mp_world,         ONLY : world_comm
USE mp,               ONLY : mp_sum, mp_bcast

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork

REAL(DP)            :: rd_ht(3,3), zero
LOGICAL             :: trd_ht
CHARACTER(LEN=10)   :: cell_units
!
!  First collect the total energies
!
CALL mp_sum(energy_geo, world_comm)
energy_geo=energy_geo / nproc_image

IF (lmurn) THEN
   CALL do_ev()
   CALL write_mur(vmin,b0,b01,b02,emin)
   CALL write_mur_p()
   CALL plot_mur()
   CALL plot_mur_p()
   IF (lgeo_from_file) THEN
      CALL compute_celldm_geo_file(vmin, celldm0, pressure_kb)
   ELSE
      CALL compute_celldm_geo(vmin, celldm0, &
                   celldm_geo(1,central_geo), omega_geo(central_geo))
   ENDIF
ELSE
   CALL write_gnuplot_energy(nwork)
   CALL quadratic_fit()
   CALL write_quadratic()
   CALL plot_multi_energy()
   CALL write_e_omega()
   CALL plot_mur()
   CALL plot_geo_p()
ENDIF
CALL mp_bcast(celldm0, meta_ionode_id, world_comm)
CALL mp_bcast(emin, meta_ionode_id, world_comm)
CALL write_minimum_energy_data()

IF (.NOT.lmurn) THEN
!
!   If npress>0 the crystal parameters as a function of pressure have
!   been calculated. Interpolate the elastic constants at those
!   geometries.
!
   CALL write_elastic_p()
!
ELSEIF (lgeo_from_file) THEN
   CALL write_elastic_mur_p()
ENDIF
CALL plot_elastic_t1(2, .FALSE.)
CALL plot_elastic_t1(3, .FALSE.)
CALL plot_macro_elastic_t1()
CALL plot_sound_speed_t1()

trd_ht=.FALSE.
rd_ht=0.0_DP
zero=0.0_DP
cell_units='alat'
CALL cell_base_init ( ibrav_save, celldm0, zero, zero, zero, zero, &
                                      zero, zero, trd_ht, rd_ht, cell_units )
CALL set_fft_mesh()
omega0=omega
at0(:,:)=at(:,:)
!
!   strain the initial tau uniformly to the new at and set the equilibrium tau
!   to the strained tau
!
CALL adjust_tau(tau_save_crys, tau, at)
tau0=tau
!
!  bring the strained tau in crystal coordinates
!
tau0_crys=tau
CALL cryst_to_cart(nat, tau0_crys, bg, -1)
!
!  recompute the density at the minimum volume
!
     CALL compute_density(omega,density,.TRUE.)
!
!  compute the xrdp at the minimum volume if required by the user
!
     IF (lxrdp) CALL manage_xrdp(' ')

RETURN
END SUBROUTINE manage_energy_minimum
