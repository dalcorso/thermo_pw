!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_ev_input(file_dat, pressure)
  !-----------------------------------------------------------------------
  !
  !  This routine receives the summary of the cell volumes and of
  !  the energy at each volume and writes the input file 
  !  for the ev_sub subroutine.
  !
  !
  USE kinds,            ONLY : DP
  USE mp_images,        ONLY : my_image_id, root_image
  USE thermo_mod,       ONLY : omega_geo, energy_geo, ngeo
  USE io_global,        ONLY : ionode

  IMPLICIT NONE
  CHARACTER(LEN=256) :: file_dat
  CHARACTER(LEN=256) :: filedata
  REAL(DP) ::  pressure
  INTEGER :: iu_ev, igeom
  INTEGER :: find_free_unit
  !
  IF (my_image_id /= root_image) RETURN
  !
  filedata=TRIM(file_dat)
  CALL write_ev_driver(filedata) 
  !
  IF (ionode) THEN
     iu_ev=find_free_unit()
     OPEN(UNIT=iu_ev, FILE=TRIM(filedata), STATUS='UNKNOWN', FORM='FORMATTED')
     DO igeom=1,ngeo(1)
        WRITE(iu_ev,'(2e30.15)') omega_geo(igeom), energy_geo(igeom) + &
                                        pressure * omega_geo(igeom)
     ENDDO
     !
     CLOSE(UNIT=iu_ev,STATUS='keep')
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE write_ev_input

!-----------------------------------------------------------------------
SUBROUTINE write_ev_driver(file_dat)
!-----------------------------------------------------------------------

USE io_global,        ONLY : ionode
USE control_ev,       ONLY : ieos

IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: file_dat
CHARACTER(LEN=256) :: filename
INTEGER :: iu_ev
INTEGER :: find_free_unit
!
IF (ionode) THEN
   iu_ev=find_free_unit()
   filename='energy_files/input_ev'
   CALL add_pressure(filename)
     
   OPEN(UNIT=iu_ev, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_ev,'("au")')
   WRITE(iu_ev,'("hex")')
   WRITE(iu_ev,'(i3)') ieos
   WRITE(iu_ev,'(a)') TRIM(file_dat)
   WRITE(iu_ev,'(a)') TRIM(file_dat)//'.ev.out'
   CLOSE(iu_ev)
ENDIF
!
RETURN
END SUBROUTINE write_ev_driver

!-----------------------------------------------------------------------
SUBROUTINE do_ev()
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!
USE constants,   ONLY : ry_kbar
USE thermo_mod,  ONLY : ngeo, omega_geo, energy_geo
USE control_mur, ONLY : vmin, b0, b01, b02, emin
USE control_mur_p,  ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
USE control_ev,  ONLY : ieos, npt, v0, e0
USE temperature, ONLY : ntemp_plot
USE data_files,  ONLY : flevdat
USE control_pressure, ONLY : pressure, press, npress, npress_plot
USE io_global,   ONLY : meta_ionode_id, stdout
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_bcast

IMPLICIT NONE
CHARACTER(LEN=256) :: file_dat, filename
INTEGER :: ipress, ipt

  file_dat="energy_files/"//TRIM(flevdat) 
  CALL add_pressure(file_dat)
  filename='energy_files/input_ev'
  CALL add_pressure(filename)

  CALL write_ev_input(file_dat, pressure)
  CALL ev_sub(vmin, b0, b01, b02, emin, filename )

  CALL mp_bcast(vmin, meta_ionode_id, world_comm)
  CALL mp_bcast(b0, meta_ionode_id, world_comm)
  CALL mp_bcast(b01, meta_ionode_id, world_comm)
  CALL mp_bcast(b02, meta_ionode_id, world_comm)
  CALL mp_bcast(emin, meta_ionode_id, world_comm)
!
!    Compute also the Murnaghan parameters as a function of pressure 
!
  IF (ntemp_plot>0.OR.npress_plot>0) THEN
     npt=ngeo(1)
     ALLOCATE(v0(npt))
     ALLOCATE(e0(npt))
     DO ipt=1, npt
        v0(ipt)=omega_geo(ipt)
     ENDDO
     DO ipress=1,npress
        DO ipt=1,npt
           e0(ipt)=energy_geo(ipt) + press(ipress) * v0(ipt) / ry_kbar
        ENDDO
        CALL ev_sub_nodisk(vmin_p(ipress), b0_p(ipress), b01_p(ipress), &
                                  b02_p(ipress), emin_p(ipress) )
     ENDDO
     DEALLOCATE(e0)
     DEALLOCATE(v0)
  ENDIF
  !
  RETURN
END SUBROUTINE do_ev

!-----------------------------------------------------------------------
SUBROUTINE do_ev_t(itemp)
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus and
!  its derivative at a given temperature. It uses the free energy 
!  computed by phonon dos.
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo, no_ph
USE control_mur,    ONLY : vmin, b0, b01, b02, emin
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_eldos,  ONLY : lel_free_energy
USE thermodynamics, ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE anharmonic,     ONLY : vmin_t, b0_t, b01_t, b02_t, free_e_min_t, a_t
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE data_files,     ONLY : flevdat
USE polyfit_mod,    ONLY : polyfit
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER  :: m1
INTEGER :: idata, ndata, i1
REAL(DP) :: a(poly_degree_ph+1), x(ngeo(1)), y(ngeo(1)), aux, aux1, aux2, aux3


  IF (my_image_id /= root_image) RETURN

  m1=poly_degree_ph+1
!  WRITE(stdout,*) 
  ndata=0
  DO idata=1,ngeo(1)
     IF (no_ph(idata)) CYCLE
     ndata=ndata+1
     x(ndata)=omega_geo(idata)
     y(ndata)=ph_free_ener(itemp,idata)
     IF (lel_free_energy) y(ndata)=y(ndata)+el_free_ener(itemp,idata)
!    WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
  ENDDO
  CALL polyfit(x, y, ndata, a, poly_degree_ph)

  CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, b02 * ry_kbar, a, m1, vm)
  CALL eos_energy_pol(ieos, vm, aux3, vmin, b0/ry_kbar, b01, b02*ry_kbar, a, m1)
  CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin, b0/ry_kbar, b01, &
                                                           b02 * ry_kbar, a, m1)
  vmin_t(itemp)=vm
  b0_t(itemp)=aux * ry_kbar
  b01_t(itemp)=aux1 
  b02_t(itemp)=aux2 / ry_kbar
  free_e_min_t(itemp)=emin+ aux3
  a_t(1:m1,itemp)=a(1:m1)
  !
  CALL compute_celldm_geo(vmin_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("-"))')
  WRITE(stdout,'(5x, "free energy from phonon dos, at T= ", f12.6)') &
                                                                temp(itemp)
  CALL summarize_mur(celldm_(1), b0_t(itemp), b01_t(itemp), b02_t(itemp), &
                                                     free_e_min_t(itemp))

  WRITE(stdout,'(2x,76("-"),/)')

  RETURN
END SUBROUTINE do_ev_t
!

!-----------------------------------------------------------------------
SUBROUTINE do_ev_t_ph(itemp) 
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!  at a given temperature. It uses the free energy computed by the
!  sum over the phonon frequencies
!
USE kinds,          ONLY : DP
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo, no_ph
USE constants,      ONLY : ry_kbar
USE ph_freq_thermodynamics, ONLY : phf_free_ener
USE el_thermodynamics,      ONLY : el_free_ener
USE ph_freq_anharmonic,     ONLY : vminf_t, b0f_t, b01f_t, b02f_t, free_e_minf_t
USE control_mur,    ONLY : emin, vmin, b0, b01, b02
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_eldos,  ONLY : lel_free_energy
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE polyfit_mod,    ONLY : polyfit
USE data_files,     ONLY : flevdat
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER :: m1
INTEGER :: idata, ndata, i1
REAL(DP) :: a(poly_degree_ph+1), x(ngeo(1)), y(ngeo(1)), aux, aux1, aux2, aux3
!REAL(DP) :: compute_mur_fun


  IF (my_image_id /= root_image) RETURN

  m1=poly_degree_ph+1
!  WRITE(stdout,*) 
  ndata=0
  DO idata=1,ngeo(1)
     IF (no_ph(idata)) CYCLE
     ndata=ndata+1
     x(ndata)=omega_geo(idata)
     y(ndata)=phf_free_ener(itemp,idata) 
     IF (lel_free_energy) y(ndata)=y(ndata)+el_free_ener(itemp,idata)
!    WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
  ENDDO
  CALL polyfit(x, y, ndata, a, poly_degree_ph)

  CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, b02*ry_kbar, a, m1, vm)
  CALL eos_energy_pol(ieos, vm, aux3, vmin, b0/ry_kbar, b01, b02*ry_kbar, a, m1)
  CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin, b0/ry_kbar, &
                              b01, b02 * ry_kbar, a, m1)

  vminf_t(itemp)=vm
  b0f_t(itemp)=aux * ry_kbar
  b01f_t(itemp)=aux1
  b02f_t(itemp)=aux2 / ry_kbar
  free_e_minf_t(itemp)=emin + aux3
  !
  CALL compute_celldm_geo(vminf_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x, "free energy from phonon frequencies at T=",f12.4)') &
                                                               temp(itemp)
  CALL summarize_mur(celldm_(1), b0f_t(itemp), b01f_t(itemp), b02f_t(itemp),&
                                                      free_e_minf_t(itemp))

  WRITE(stdout,'(2x,76("+"),/)')

  RETURN
END SUBROUTINE do_ev_t_ph

!-----------------------------------------------------------------------
SUBROUTINE do_ev_t_el(itemp)
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume and bulk modulus
!  at a given temperature. It uses the free energy computed by phonon dos
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE thermo_mod,     ONLY : ngeo, omega_geo, energy_geo, celldm_geo, &
                           central_geo
USE control_mur,    ONLY : vmin, b0, b01, b02, emin
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE el_thermodynamics, ONLY : el_free_ener
USE el_anharmonic,  ONLY : vmine_t, b0e_t, b01e_t, b02e_t, free_e_mine_t
USE temperature,    ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE data_files,     ONLY : flevdat
USE polyfit_mod,    ONLY : polyfit
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
INTEGER :: igeom, iu_ev
REAL(DP) :: free_e, vm, celldm_(6)
INTEGER  :: m1
INTEGER :: idata, ndata, i1
REAL(DP) :: a(poly_degree_ph+1), x(ngeo(1)), y(ngeo(1)), aux, aux1, aux2, aux3

  IF (my_image_id /= root_image) RETURN

  m1=poly_degree_ph+1
  WRITE(stdout,*) ngeo(1), central_geo
  ndata=0
  DO idata=1,ngeo(1)
     ndata=ndata+1
     x(ndata)=omega_geo(idata)
     y(ndata)=el_free_ener(itemp,idata)
     WRITE(stdout,'(2f25.14)') x(ndata), y(ndata)
  ENDDO
  WRITE(stdout,*)
  CALL polyfit(x, y, ndata, a, poly_degree_ph)

  CALL find_min_mur_pol(vmin, b0 / ry_kbar, b01, b02 * ry_kbar, a, m1, vm)

  CALL eos_energy_pol(ieos, vm, aux3, vmin, b0/ry_kbar, b01, b02*ry_kbar, a, m1)
  CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin, b0/ry_kbar, b01, &
                                        b02 * ry_kbar, a, m1)

  vmine_t(itemp)=vm
  b0e_t(itemp)=aux
  b01e_t(itemp)=aux1
  b02e_t(itemp)=aux2
  free_e_mine_t(itemp) = emin + aux3
  !
  CALL compute_celldm_geo(vmine_t(itemp), celldm_, celldm_geo(1,central_geo), &
                                         omega_geo(central_geo))
  WRITE(stdout,'(/,2x,76("-"))')
  WRITE(stdout,'(5x, "free energy from electron dos, at T= ", f12.6)') &
                                                                temp(itemp)
  IF (pressure_kb /= 0.0_DP) &
     WRITE(stdout, '(5x,"pressure = ",f15.6," kbar")') pressure_kb
  WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
                                 &" a.u.")') celldm_(1)
  WRITE(stdout,'(5x, "The bulk modulus is ",24x,f12.3,"  kbar")')  b0e_t(itemp)
  WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",&
                               &f9.3)')  b01e_t(itemp)
  IF (ieos==2) &
  WRITE(stdout,'(5x, "The second pressure derivative of the bulk modulus is ",&
                               &f9.3)')  b02e_t(itemp)
  IF (pressure_kb /= 0.0_DP) THEN
     WRITE(stdout,'(5x,"The Gibbs energy at the minimum is    ",&
                                    6x,f20.9," Ry")') free_e_mine_t(itemp)
  ELSE
     WRITE(stdout,'(5x,"The free energy at the minimum is",6x,f20.9," Ry")') &
                                                 free_e_mine_t(itemp)
  END IF

  WRITE(stdout,'(2x,76("-"),/)')

  RETURN
END SUBROUTINE do_ev_t_el

!-----------------------------------------------------------------------
SUBROUTINE do_ev_pt(itemp)
!-----------------------------------------------------------------------
!
!  This subroutine computes the equilibrium volume, bulk modulus, and its
!  derivatives as a function of temperature for the npress_plot pressures 
!  that we want to plot in output
!
USE kinds,          ONLY : DP
USE constants,      ONLY : ry_kbar
USE control_pressure, ONLY : npress, npress_plot, ipress_plot
USE control_ev,     ONLY : ieos
USE eos,            ONLY : eos_bulk_pol, eos_energy_pol
USE control_quartic_energy, ONLY : poly_degree_ph
USE control_mur_p,  ONLY : vmin_p, b0_p, b01_p, b02_p, emin_p
USE anharmonic,     ONLY : a_t
USE anharmonic_pt,  ONLY : vmin_pt, b0_pt, b01_pt, b02_pt, emin_pt
USE io_global,      ONLY : stdout
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
REAL(DP) :: free_e, vm
INTEGER  :: m1, i1, ipress, ipressp
REAL(DP) :: aux, aux1, aux2, aux3

IF (my_image_id /= root_image) RETURN

IF (npress_plot==0) RETURN

m1=poly_degree_ph+1

DO ipressp=1, npress_plot
   ipress=ipress_plot(ipressp)
   CALL find_min_mur_pol(vmin_p(ipress), b0_p(ipress) / ry_kbar, &
                 b01_p(ipress), b02_p(ipress)*ry_kbar, a_t(:,itemp), m1, vm)

   CALL eos_energy_pol(ieos, vm, aux3, vmin_p(ipress), b0_p(ipress)/ry_kbar,  &
                 b01_p(ipress), b02_p(ipress)*ry_kbar, a_t(:,itemp), m1)

   CALL eos_bulk_pol(ieos, vm, aux, aux1, aux2, vmin_p(ipress), &
                 b0_p(ipress)/ry_kbar, b01_p(ipress), b02_p(ipress) * ry_kbar, &
                                                       a_t(:,itemp), m1)
   vmin_pt(itemp,ipressp)=vm
   b0_pt(itemp,ipressp)=aux * ry_kbar
   b01_pt(itemp,ipressp)=aux1
   b02_pt(itemp,ipressp)=aux2 / ry_kbar
   emin_pt(itemp,ipressp)=emin_p(ipress)+aux3
   !
ENDDO

RETURN
END SUBROUTINE do_ev_pt
!
!-----------------------------------------------------------------------
SUBROUTINE do_ev_vt(itemp)
!-----------------------------------------------------------------------
!
!  This subroutine computes the thermal pressure as a function of 
!  temperature for the volumes specified in input.
!
USE kinds,            ONLY : DP
USE thermo_mod,       ONLY : omega_geo
USE constants,        ONLY : ry_kbar
USE control_pressure, ONLY : npress, npress_plot, ipress_plot
USE control_quartic_energy, ONLY : poly_degree_ph
USE anharmonic,       ONLY : a_t
USE anharmonic_vt,    ONLY : press_vt
USE polyfit_mod,      ONLY : compute_poly_deriv
USE control_vol,      ONLY : nvol_plot, ivol_plot
USE io_global,        ONLY : stdout
USE mp_images,        ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: itemp
REAL(DP) :: press, omega
INTEGER  :: m1, ivolp, igeo

IF (my_image_id /= root_image) RETURN

IF (nvol_plot==0) RETURN

m1=poly_degree_ph+1

DO ivolp=1,nvol_plot
   igeo = ivol_plot(ivolp)
   omega = omega_geo(igeo)
   CALL compute_poly_deriv(omega, m1-1, a_t(:,itemp), press)
   press = -press
   press_vt(itemp,ivolp) = press * ry_kbar
ENDDO

RETURN
END SUBROUTINE do_ev_vt

!-----------------------------------------------------------------------
SUBROUTINE find_min_mur_pol(v0, b0, b01, b02, a, m1, vm)
!-----------------------------------------------------------------------
!
!  This routine minimizes a function equal to an equation os state
!  (Birch-Murnaghan third or fourth order, or Munaghan)
!  with parameters v0, b0, b01, and b02 and a polynomial of degree m1-1 of
!  the form a(1) + a(2) * v + a(3) * v**2 +... a(m1) * v**(m1-1)
!  NB: b0 must be in atomic units not kbar
!  NB: b02 must be in atomic units not 1/kbar
!  The polynomial gives the energy in Ry.
!  On output vm is in a.u.
!
USE kinds, ONLY : DP 
USE io_global, ONLY : stdout
USE eos, ONLY : eos_press_pol, eos_dpress_pol
USE control_ev, ONLY : ieos
IMPLICIT NONE
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, a(m1)
REAL(DP), INTENT(OUT) :: vm
REAL(DP), PARAMETER :: tol=1.D-12
INTEGER, PARAMETER :: maxstep=100
REAL(DP) :: fx, fxp, v1, v1old
INTEGER :: istep

v1=v0
v1old=v1
DO istep=1,maxstep
   CALL eos_press_pol(ieos, v1, fx, v0, b0, b01, b02, a, m1)
   CALL eos_dpress_pol(ieos, v1, fxp, v0, b0, b01, b02, a, m1)
!
!  Newton method
!
   v1 = v1 + fx/fxp
!   WRITE(stdout,'(5x,"Step", i4, " V1=", f20.12, " f= ", f20.12)') istep, v1, -fx
!   FLUSH(stdout)
   IF (ABS(v1-v1old) < tol .OR. ABS(fx) < tol ) GOTO 100
   v1old=v1
ENDDO
CALL errore('find_min_mur_pol','minimum not found',1)
100 CONTINUE
vm=v1
!WRITE(stdout,'("Vmin", 3f20.12)') vm

RETURN
END SUBROUTINE find_min_mur_pol
!
!-------------------------------------------------------------------------
SUBROUTINE summarize_mur(celldm, b0, b01, b02, free_e_min)
!-------------------------------------------------------------------------
!
USE kinds,     ONLY : DP
USE control_ev, ONLY : ieos
USE control_pressure, ONLY : pressure_kb
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) ::  celldm, b0, b01, b02, free_e_min

IF (pressure_kb /= 0.0_DP) &
   WRITE(stdout, '(5x,"pressure = ",f15.6," kbar")') pressure_kb
IF (ieos==1) THEN
   WRITE(stdout,'(5x,"Using Birch-Murnaghan equation of order 3")')
ELSEIF(ieos==2) THEN
   WRITE(stdout,'(5x,"Using Birch-Murnaghan equation of order 4")')
ELSEIF(ieos==4) THEN
   WRITE(stdout,'(5x,"Using the Murnaghan equation")')
ENDIF
WRITE(stdout,'(5x, "The equilibrium lattice constant is ",9x,f12.4,&
                                 &" a.u.")') celldm
WRITE(stdout,'(5x, "The bulk modulus is ",24x,f12.3,"  kbar")')  b0
WRITE(stdout,'(5x, "The pressure derivative of the bulk modulus is ",&
                               &f9.3)')  b01
IF (ieos==2) &
   WRITE(stdout,'(5x, "The second derivative of the bulk &
                     &modulus is ", f13.5," 1/kbar")')  b02
IF (pressure_kb /= 0.0_DP) THEN
   WRITE(stdout,'(5x,"The Gibbs energy at the minimum is    ",&
                                    6x,f20.9," Ry")') free_e_min
ELSE
   WRITE(stdout,'(5x,"The free energy at the minimum is",6x,f20.9," Ry")') &
                                                 free_e_min
END IF

RETURN
END SUBROUTINE summarize_mur
!
