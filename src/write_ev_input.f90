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
USE control_ev,  ONLY : npt, v0, e0
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
     CALL mp_bcast(vmin_p, meta_ionode_id, world_comm)
     CALL mp_bcast(b0_p, meta_ionode_id, world_comm)
     CALL mp_bcast(b01_p, meta_ionode_id, world_comm)
     CALL mp_bcast(b02_p, meta_ionode_id, world_comm)
     CALL mp_bcast(emin_p, meta_ionode_id, world_comm)
  ENDIF
  !
  RETURN
END SUBROUTINE do_ev

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
  b0e_t(itemp)=aux * ry_kbar
  b01e_t(itemp)=aux1 
  b02e_t(itemp)=aux2 / ry_kbar
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
