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
                             energy_geo, uint_geo, tau_geo, celldm_geo_eos, &
                             omega_geo_eos, energy_geo_eos, uint_geo_eos,   &
                             tau_geo_eos, tot_ngeo, tot_ngeo_eos
USE control_mur,      ONLY : vmin, b0, b01, b02, emin, lmurn
USE control_thermo,   ONLY : lgeo_from_file
USE control_atomic_pos, ONLY : linternal_thermo, iconstr_internal, nint_var, &
                               tau_eq, linterpolate_tau
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
INTEGER, INTENT(INOUT) :: nwork

REAL(DP)            :: rd_ht(3,3), zero
REAL(DP), ALLOCATABLE :: energy_geo_eff(:), celldm_geo_eff(:,:), min_u(:,:)
REAL(DP) :: compute_omega_geo
LOGICAL             :: trd_ht, lreturn
CHARACTER(LEN=10)   :: cell_units
INTEGER :: nwork_eff, igeom, iwork
!
!  First collect the total energies
!
CALL mp_sum(energy_geo, world_comm)
energy_geo=energy_geo / nproc_image
!
!   First check that all geometries have been computed
!
lreturn=.FALSE.
DO iwork=1,nwork
   lreturn=lreturn.OR.(ABS(energy_geo(iwork))<1.D-10)
ENDDO
IF (lreturn) RETURN

ALLOCATE(energy_geo_eos(nwork))
ALLOCATE(omega_geo_eos(nwork))
ALLOCATE(celldm_geo_eos(6,nwork))

IF (linternal_thermo) THEN
   ALLOCATE(uint_geo_eos(nint_var,tot_ngeo_eos))
   ALLOCATE(tau_geo_eos(3,nat,tot_ngeo_eos))
   CALL redefine_energies_min(energy_geo, celldm_geo, uint_geo, nwork, &
        energy_geo_eos, celldm_geo_eos, uint_geo_eos, nwork_eff)
   DO igeom = 1, tot_ngeo_eos
      omega_geo_eos(igeom)=compute_omega_geo(ibrav_save,&
                                                   celldm_geo_eos(:,igeom))
      CALL internal_to_tau(celldm_geo_eos(1,igeom), tau_geo_eos(1,1,igeom), &
                          uint_geo_eos(1,igeom), nat, nint_var,      &
                                             iconstr_internal, 1)
   ENDDO
   CALL summarize_eos_geometries(tot_ngeo_eos)
ELSE
   IF (linterpolate_tau) THEN
      ALLOCATE(uint_geo_eos(nint_var,nwork))
      ALLOCATE(tau_geo_eos(3,nat,nwork))
   ENDIF
   energy_geo_eos(1:nwork) = energy_geo(1:nwork)
   omega_geo_eos(1:nwork) = omega_geo(1:nwork)
   celldm_geo_eos(1:6,1:nwork)=celldm_geo(1:6,1:nwork)
ENDIF

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
                   celldm_geo_eos(1,central_geo), omega_geo_eos(central_geo))
   ENDIF
ELSE
   CALL write_gnuplot_energy(tot_ngeo_eos)
   CALL quadratic_fit()
   CALL write_quadratic()
   CALL write_e_omega()
   CALL compute_celldm_pm()
   CALL compute_bulk_modulus()
   CALL plot_multi_energy()
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
CALL plot_macro_elastic_p()
CALL plot_sound_speed_p()

trd_ht=.FALSE.
rd_ht=0.0_DP
zero=0.0_DP
cell_units='alat'
CALL cell_base_init ( ibrav_save, celldm0, zero, zero, zero, zero, &
                                      zero, zero, trd_ht, rd_ht, cell_units )
CALL set_fft_mesh()
omega0=omega
at0(:,:)=at(:,:)
IF (linternal_thermo) THEN
   tau=tau_eq
ELSE
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
ENDIF
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

!-----------------------------------------------------------------------
SUBROUTINE summarize_eos_geometries(nwork)
!-----------------------------------------------------------------------
USE constants,     ONLY : bohr_radius_si
USE thermo_mod,    ONLY : ibrav_geo, omega_geo_eos, celldm_geo_eos, &
                          tau_geo_eos, uint_geo_eos
USE initial_conf,  ONLY : celldm_save, ibrav_save, tau_save, atm_save, &
                          ityp_save
USE control_atomic_pos, ONLY : nint_var, uint0 
USE io_global,     ONLY : stdout
USE ions_base,     ONLY : nat

IMPLICIT NONE
INTEGER :: nwork
INTEGER :: igeo, ipol, na, i

WRITE(stdout,'(/,5x,70("-"))')
WRITE(stdout,'(/,5x,"Number of geometries for EOS:", i5)') nwork
WRITE(stdout,'(/,5x,"Input ibrav and celldm:")')
WRITE(stdout,'(12x, i3, 6f10.5,/)') ibrav_save, celldm_save(:)

WRITE(stdout,'(/,5x,"Input atomic positions: uint=",2f15.8)') (uint0(i), &
                                                   i=1, nint_var)
WRITE( stdout, '(/,3x,"Cartesian axes")')
WRITE( stdout, '(/,5x,"site n.     atom                  positions (alat units)")')
WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
     (na, atm_save(ityp_save(na)), na, (tau_save(ipol,na), ipol=1,3), &
                                                               na=1,nat)
WRITE( stdout, *) 

WRITE(stdout,'(/,5x,"List of the geometries for EOS :")')

DO igeo=1,nwork
   WRITE(stdout,'(5x,i5,": ", i3,6f10.5)') igeo, ibrav_geo(igeo), &
                                   celldm_geo_eos(:,igeo)
   WRITE( stdout,'(5x,"Atomic positions found: uint=",2f15.8)') &
                                 (uint_geo_eos(i,igeo), i=1, nint_var)
   WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")')       &
        (na, atm_save(ityp_save(na)), na, (tau_geo_eos(ipol,na,igeo),        &
                                                   ipol=1,3),  na=1,nat)
   WRITE( stdout, *) 
ENDDO

WRITE(stdout,'(/,5x,"Volumes: ",10x,"(a.u.)^3",10x,"(A)^3")')
DO igeo=1,nwork
   WRITE(stdout,'(5x,i5,2f20.10)') igeo, omega_geo_eos(igeo), &
                            omega_geo_eos(igeo)*(bohr_radius_si)**3/1.D-30
ENDDO
WRITE(stdout,'(/,5x,70("-"))')

RETURN
END SUBROUTINE summarize_eos_geometries
