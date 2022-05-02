! Copyright (C) 2021 Andrea Dal Corso and Xuejun Gong
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE eos
!
!   this module contains the support routines for the equations of
!   state and their derivatives. Some of these formulas were contained
!   in the routine ev.f90 of the QE distribution, but cannot
!   be used because they are internal function to the ev.f90 routine.
!
!   The modulus provides for each eos some routines. All these routines
!   receive as input the parameters: 
!   ieos: the required equation of state
!         1 - Birch-Murnaghan 3 order
!         2 - Birch-Murnaghan 4 order
!         4 - Murnaghan
!   
!   v0: the equilibrium volume
!   b0: the bulk modulus at equilibrium
!   b01: the derivative of the bulk modulus with respect to pressure
!   b02: the second derivative of the bulk modulus with
!        respect to pressure (used only for ieos=2)
!   
!   The routines provided by this module are:
!   eos_energy gives as output the energy at the given input volume
!   eos_press  gives as output the pressure at the given input volume
!   eos_dpress gives as output minus the derivative of the pressure with
!              respect to volume (second derivative of energy with respect
!              to the volume)
!   eos_bulk   gives as output the bulk modulus and its first derivative
!              with respect to pressure at the given input volume
!              when ieos=2 gives also the second derivative of the
!              bulk modulus with respect to pressure
!
!   All routines have also a version that in addition to the parameters
!   of the eos receives the coefficients of a polynomial which is
!   supposed to interpolate the free energy and gives the given quantity 
!   assuming the eos is the given function plus the polynomial.
!
!   In input the volume is in (a.u.)^3, the bulk modulus is in kbar
!   and the output energy is in Ry.
!
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC  eos_energy, eos_press, eos_dpress, eos_bulk,      &
          eos_energy_pol, eos_press_pol, eos_dpress_pol, eos_bulk_pol

CONTAINS
!-------------------------------------------------------------------
SUBROUTINE eos_energy(ieos, vm, eout, v0, b0, b01, b02)
!-------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN)   :: ieos
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02
REAL(DP), INTENT(OUT) :: eout

IF (ieos==1) THEN
   CALL birch3_energy(v0, b0, b01, vm, eout)
ELSEIF (ieos==2) THEN
   CALL birch4_energy(v0, b0, b01, b02, vm, eout)
ELSEIF (ieos==4) THEN
   CALL murnaghan_energy(v0, b0, b01, vm, eout)
ELSE
   CALL errore('eos_energy','ieos not programmed',1)
ENDIF

RETURN
END SUBROUTINE eos_energy
!
!-------------------------------------------------------------------
SUBROUTINE eos_press(ieos, vm, press, v0, b0, b01, b02)
!-------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN)   :: ieos
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02
REAL(DP), INTENT(OUT) :: press

IF (ieos==1) THEN
   CALL birch3_press(v0, b0, b01, vm, press)
ELSEIF (ieos==2) THEN
   CALL birch4_press(v0, b0, b01, b02, vm, press)
ELSEIF (ieos==4) THEN
   CALL murnaghan_press(v0, b0, b01, vm, press)
ELSE
   CALL errore('eos_press','ieos not programmed',1)
ENDIF

RETURN
END SUBROUTINE eos_press
!
!-------------------------------------------------------------------
SUBROUTINE eos_dpress(ieos, vm, dpress, v0, b0, b01, b02)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)   :: ieos
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02
REAL(DP), INTENT(OUT) :: dpress

IF (ieos==1) THEN
   CALL birch3_dpress(v0, b0, b01, vm, dpress)
ELSEIF (ieos==2) THEN
   CALL birch4_dpress(v0, b0, b01, b02, vm, dpress)
ELSEIF (ieos==4) THEN
   CALL murnaghan_dpress(v0, b0, b01, vm, dpress)
ELSE
   CALL errore('eos_dpress','ieos not programmed',1)
ENDIF

RETURN
END SUBROUTINE eos_dpress
!
!-------------------------------------------------------------------
SUBROUTINE eos_bulk(ieos, vm, bulk, dbulk, d2bulk, v0, b0, b01, b02)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)   :: ieos
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02
REAL(DP), INTENT(OUT) :: bulk, dbulk, d2bulk

d2bulk=0.0_DP
IF (ieos==1) THEN
   CALL birch3_bulk(v0, b0, b01, vm, bulk, dbulk)
ELSEIF (ieos==2) THEN
   CALL birch4_bulk(v0, b0, b01, b02, vm, bulk, dbulk, d2bulk)
ELSEIF (ieos==4) THEN
   CALL murnaghan_bulk(v0, b0, b01, vm, bulk, dbulk)
ELSE
   CALL errore('eos_bulk','ieos not programmed',1)
ENDIF

RETURN
END SUBROUTINE eos_bulk
!
!-------------------------------------------------------------------
SUBROUTINE eos_energy_pol(ieos, vm, eout, v0, b0, b01, b02, a, m1)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)   :: ieos, m1
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02, a(m1)
REAL(DP), INTENT(OUT) :: eout

IF (ieos==1) THEN
   CALL birch3_energy_pol(v0, b0, b01, vm, a, m1, eout)
ELSEIF (ieos==2) THEN
   CALL birch4_energy_pol(v0, b0, b01, b02, vm, a, m1, eout)
ELSEIF (ieos==4) THEN
   CALL murnaghan_energy_pol(v0, b0, b01, vm, a, m1, eout)
ELSE
   CALL errore('eos_energy_pol','ieos not programmed',1)
ENDIF

RETURN
END SUBROUTINE eos_energy_pol
!
!-------------------------------------------------------------------
SUBROUTINE eos_press_pol(ieos, vm, press, v0, b0, b01, b02, a, m1)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)   :: ieos, m1
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02, a(m1)
REAL(DP), INTENT(OUT) :: press

IF (ieos==1) THEN
   CALL birch3_press_pol(v0, b0, b01, vm, a, m1, press)
ELSEIF (ieos==2) THEN
   CALL birch4_press_pol(v0, b0, b01, b02, vm, a, m1, press)
ELSEIF (ieos==4) THEN
   CALL murnaghan_press_pol(v0, b0, b01, vm, a, m1, press)
ELSE
   CALL errore('eos_press_pol','ieos not programmed',1)
ENDIF

RETURN
END SUBROUTINE eos_press_pol
!
!-------------------------------------------------------------------
SUBROUTINE eos_dpress_pol(ieos, vm, dpress, v0, b0, b01, b02, a, m1)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)   :: ieos, m1
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02, a(m1)
REAL(DP), INTENT(OUT) :: dpress

IF (ieos==1) THEN
   CALL birch3_dpress_pol(v0, b0, b01, vm, a, m1, dpress)
ELSEIF (ieos==2) THEN
   CALL birch4_dpress_pol(v0, b0, b01, b02, vm, a, m1, dpress)
ELSEIF (ieos==4) THEN
   CALL murnaghan_dpress_pol(v0, b0, b01, vm, a, m1, dpress)
ELSE
   CALL errore('eos_dpress_pol','ieos not programmed',1)
ENDIF

RETURN
END SUBROUTINE eos_dpress_pol
!
!-------------------------------------------------------------------
SUBROUTINE eos_bulk_pol(ieos, vm, bulk, dbulk, d2bulk, v0, b0, &
                              b01, b02, a, m1)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)   :: ieos, m1
REAL(DP), INTENT(IN)  :: vm, v0, b0, b01, b02, a(m1)
REAL(DP), INTENT(OUT) :: bulk, dbulk, d2bulk

d2bulk=0.0_DP
IF (ieos==1) THEN
   CALL birch3_bulk_pol(v0, b0, b01, vm, a, m1, bulk, dbulk)
ELSEIF (ieos==2) THEN
   CALL birch4_bulk_pol(v0, b0, b01, b02, vm, a, m1, bulk, dbulk, d2bulk)
ELSEIF (ieos==4) THEN
   CALL murnaghan_bulk_pol(v0, b0, b01, vm, a, m1, bulk, dbulk)
ELSE
   CALL errore('eos_bulk_pol','ieos not programmed',1)
ENDIF
RETURN
END SUBROUTINE eos_bulk_pol
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_energy(v0, b0, b01, vm, eout)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: eout

eout = b0 * v0 * (( v0 / vm )**(b01-1.0_DP) / (b01-1.0_DP) + (vm/v0) )/b01 &
           - v0 * b0 / (b01-1.0_DP) 

RETURN
END SUBROUTINE murnaghan_energy
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_energy_pol(v0, b0, b01, vm, a, m1, eout)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: eout

REAL(DP) :: aux

CALL murnaghan_energy(v0, b0, b01, vm, eout)
CALL compute_poly(vm, m1-1, a, aux)
eout = eout + aux

RETURN
END SUBROUTINE murnaghan_energy_pol
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_press(v0, b0, b01, vm, press)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: press

press = b0 * ((v0 / vm)**b01 - 1.0_DP) / b01

RETURN
END SUBROUTINE murnaghan_press
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_press_pol(v0, b0, b01, vm, a, m1, press)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: press

REAL(DP) :: aux

CALL murnaghan_press(v0, b0, b01, vm, press)
CALL compute_poly_deriv(vm, m1-1, a, aux)

press= press - aux

RETURN
END SUBROUTINE murnaghan_press_pol
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_dpress(v0, b0, b01, vm, dpress)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: dpress

dpress = b0 * (v0 / vm)**b01 / vm

RETURN
END SUBROUTINE murnaghan_dpress
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_dpress_pol(v0, b0, b01, vm, a, m1, dpress)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv2
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: dpress

REAL(DP) :: aux

CALL murnaghan_dpress(v0, b0, b01, vm, dpress)
CALL compute_poly_deriv2(vm, m1-1, a, aux)

dpress = dpress + aux

RETURN
END SUBROUTINE murnaghan_dpress_pol
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_bulk(v0, b0, b01, vm, bulk, dbulk)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: bulk, dbulk

bulk = b0 * (v0 / vm)**b01 

dbulk= b01

RETURN
END SUBROUTINE murnaghan_bulk
!
!-------------------------------------------------------------------
SUBROUTINE murnaghan_bulk_pol(v0, b0, b01, vm, a, m1, bulk, dbulk)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv2, compute_poly_deriv3
!
IMPLICIT NONE
INTEGER,  INTENT(IN)  :: m1
REAL(DP), INTENT(IN)  :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: bulk, dbulk

REAL(DP) :: aux, aux1

CALL murnaghan_bulk(v0, b0, b01, vm, bulk, dbulk)

dbulk=dbulk*bulk
CALL compute_poly_deriv2(vm, m1-1, a, aux)
bulk = bulk + aux * vm
CALL compute_poly_deriv3(vm, m1-1, a, aux1)
dbulk = dbulk - vm * ( aux + vm * aux1)
dbulk = dbulk / bulk

RETURN
END SUBROUTINE murnaghan_bulk_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch3_energy(v0, b0, b01, vm, eout)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: eout
REAL(DP) :: x, y, x2, c1

x=v0/vm
x2=x*x
y=x2**(1.0_DP/3.0_DP)
c1 = 3.d0*(b01-4.d0)/8.d0

eout = 4.5d0*b0*v0*( (-0.5d0+c1)*y     &
      +(0.25_DP-c1)*y**2               &
      +c1*x2/3.d0+(1.d0/4.d0-c1/3.d0) )

RETURN
END SUBROUTINE birch3_energy
!
!-------------------------------------------------------------------
SUBROUTINE birch3_energy_pol(v0, b0, b01, vm, a, m1, eout)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: eout

REAL(DP) :: aux

CALL birch3_energy(v0, b0, b01, vm, eout)

CALL compute_poly(vm, m1-1, a, aux)
eout = eout + aux

RETURN
END SUBROUTINE birch3_energy_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch3_press(v0, b0, b01, vm, press)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: press
REAL(DP) :: x, x2, x3, y, c1

x=v0/vm
x2=x*x
x3=x2*x
y=x2**(1.0_DP/3.0_DP)
c1 = 3.d0*(b01-4.d0)/8.d0

press=3.d0*b0*((-0.5d0+c1)*y*x   &
           +(0.5d0-2.d0*c1)*y*y*x  &
           + c1*x3 )
RETURN
END SUBROUTINE birch3_press
!
!-------------------------------------------------------------------
SUBROUTINE birch3_press_pol(v0, b0, b01, vm, a, m1, press)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: press

REAL(DP) :: aux

CALL birch3_press(v0, b0, b01, vm, press)
CALL compute_poly_deriv(vm, m1-1, a, aux)
press= press - aux

RETURN
END SUBROUTINE birch3_press_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch3_dpress(v0, b0, b01, vm, dpress)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: dpress
REAL(DP) :: x, x2, x3, y, c1, dbdv

x=v0/vm
x2=x*x
x3=x2*x
y=x2**(1.0_DP/3.0_DP)
c1 = 3.d0*(b01-4.d0)/8.d0
 
dpress = b0 * ( 9.0_DP * c1 * x3 &
          + (- 14.0_DP * c1 + 3.5_DP)*y*y*x &
          + (5.0_DP * c1 - 2.5_DP)*y*x ) / vm

RETURN
END SUBROUTINE birch3_dpress
!
!-------------------------------------------------------------------
SUBROUTINE birch3_dpress_pol(v0, b0, b01, vm, a, m1, dpress)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv2
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: dpress

REAL(DP) :: aux

CALL birch3_dpress(v0, b0, b01, vm, dpress)
CALL compute_poly_deriv2(vm, m1-1, a, aux)

dpress = dpress + aux

RETURN
END SUBROUTINE birch3_dpress_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch3_bulk(v0, b0, b01, vm, bulk, dbulk)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, vm
REAL(DP), INTENT(OUT) :: bulk, dbulk

REAL(DP) :: x, x2, x3, y, c1, dbdv

x=v0/vm
x2=x*x
x3=x2*x
y=x2**(1.0_DP/3.0_DP)
c1 = 3.d0*(b01-4.d0)/8.d0

bulk = b0 * ( 9.0_DP * c1 * x3 &
          + (- 14.0_DP * c1 + 3.5_DP)*y*y*x &
          + (5.0_DP * c1 - 2.5_DP)*y*x )

dbdv = b0 * ( -27.0_DP * c1 * x2 * x2 &
        + ( 98.0_DP * c1 / 3.0_DP - 49.0_DP / 6.0_DP )*x2*y*y &
        + ( 25.0_DP / 6.0_DP - 25.0_DP * c1 / 3.0_DP )*x2*y ) /v0

dbulk = -vm * dbdv / bulk

RETURN
END SUBROUTINE birch3_bulk
!
!-------------------------------------------------------------------
SUBROUTINE birch3_bulk_pol(v0, b0, b01, vm, a, m1, bulk, dbulk)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv2, compute_poly_deriv3
!
IMPLICIT NONE
!
INTEGER,  INTENT(IN)  :: m1
REAL(DP), INTENT(IN)  :: v0, b0, b01, vm, a(m1)
REAL(DP), INTENT(OUT) :: bulk, dbulk

REAL(DP) :: aux, aux1

CALL birch3_bulk(v0, b0, b01, vm, bulk, dbulk)

dbulk=dbulk*bulk
CALL compute_poly_deriv2(vm, m1-1, a, aux)
bulk = bulk + aux * vm
CALL compute_poly_deriv3(vm, m1-1, a, aux1)
dbulk = dbulk - vm * ( aux + vm * aux1)
dbulk = dbulk / bulk

RETURN
END SUBROUTINE birch3_bulk_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch4_energy(v0, b0, b01, b02, vm, eout)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, vm
REAL(DP), INTENT(OUT) :: eout
REAL(DP) :: x, x2, c0, c1, y

x=v0/vm
x2=x*x
y=x2**(1.0_DP/3.0_DP)
c0 = (9.d0*b0*b02 + 9.d0*b01**2 - 63.d0*b01 + 143.d0)/48.d0
c1 = 3.d0*(b01-4.d0)/8.d0

eout = 4.5d0*b0*v0*( (-0.5d0 + c1 - c0)*y        &
      +(0.25_DP-c1 + 3.0_DP * c0 * 0.5_DP)*y*y   &
      +(c1/3.d0-c0)*x2 + c0 * 0.25_DP * x2 * y   &
      +(1.d0/4.d0 - c1/3.d0 + c0*0.25_DP) )

RETURN
END SUBROUTINE birch4_energy
!
!-------------------------------------------------------------------
SUBROUTINE birch4_energy_pol(v0, b0, b01, b02, vm, a, m1, eout)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, vm, a(m1)
REAL(DP), INTENT(OUT) :: eout

REAL(DP) :: aux

CALL birch4_energy(v0, b0, b01, b02, vm, eout)
CALL compute_poly(vm, m1-1, a, aux)
eout = eout + aux

RETURN
END SUBROUTINE birch4_energy_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch4_press(v0, b0, b01, b02, vm, press)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, vm
REAL(DP), INTENT(OUT) :: press
REAL(DP) :: x, x2, x3, c0, c1, y

x=v0/vm
x2=x*x
x3=x2*x
y=x2**(1.0_DP/3.0_DP)
c0 = (9.d0*b0*b02 + 9.d0*b01**2 - 63.d0*b01 + 143.d0)/48.d0
c1 = 3.d0*(b01-4.d0)/8.d0

press=3.d0*b0*( (-0.5d0+ c1 - c0 )*y*x           &
           +( 0.5d0-2.d0*c1 + 3.0_DP*c0 )*y*y*x  &
           +( c1 - 3.0_DP * c0 )*x3              &
           + c0*x3*y )

RETURN
END SUBROUTINE birch4_press
!
!-------------------------------------------------------------------
SUBROUTINE birch4_press_pol(v0, b0, b01, b02, vm, a, m1, press)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, vm, a(m1)
REAL(DP), INTENT(OUT) :: press

REAL(DP) :: aux

CALL birch4_press(v0, b0, b01, b02, vm, press)
CALL compute_poly_deriv(vm, m1-1, a, aux)
press= press - aux

RETURN
END SUBROUTINE birch4_press_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch4_dpress(v0, b0, b01, b02, vm, dpress)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, vm
REAL(DP), INTENT(OUT) :: dpress
REAL(DP) :: x, x2, x3, x4, c0, c1, y

x=v0/vm
x2=x*x
x3=x2*x
x4=x2*x2
y=x2**(1.0_DP/3.0_DP)
c0 = (9.d0*b0*b02 + 9.d0*b01**2 - 63.d0*b01 + 143.d0)/48.d0
c1 = 3.d0*(b01-4.d0)/8.d0

dpress = b0 * ( 11.0_DP * c0 * x3 * y + (9.0_DP * c1 - 27.0_DP *c0) * x3 &
          + (- 14.0_DP * c1 +21.0_DP * c0 +3.5_DP) *y*y*x               &
          + (5.0_DP * c1 - 5.0_DP * c0 - 2.5_DP) *y*x ) / vm


RETURN
END SUBROUTINE birch4_dpress
!
!-------------------------------------------------------------------
SUBROUTINE birch4_dpress_pol(v0, b0, b01, b02, vm, a, m1, dpress)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv2
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: m1
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, vm, a(m1)
REAL(DP), INTENT(OUT) :: dpress

REAL(DP) :: aux

CALL birch4_dpress(v0, b0, b01, b02, vm, dpress)
CALL compute_poly_deriv2(vm, m1-1, a, aux)
dpress = dpress + aux

RETURN
END SUBROUTINE birch4_dpress_pol
!
!-------------------------------------------------------------------
SUBROUTINE birch4_bulk(v0, b0, b01, b02, vm, bulk, dbulk, d2bulk)
!-------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: v0, b0, b01, b02, vm
REAL(DP), INTENT(OUT) :: bulk, dbulk, d2bulk
REAL(DP) :: x, x2, x3, x4, c0, c1, y, dbdv, dbd2v

x=v0/vm
x2=x*x
x3=x2*x
x4=x2*x2
y=x2**(1.0_DP/3.0_DP)
c0 = (9.d0*b0*b02 + 9.d0*b01**2 - 63.d0*b01 + 143.d0)/48.d0
c1 = 3.d0*(b01-4.d0)/8.d0

bulk = b0 * ( 11.0_DP * c0 * x3 * y + (9.0_DP * c1 - 27.0_DP *c0) * x3 &
          + (- 14.0_DP * c1 +21.0_DP * c0 +3.5_DP) *y*y*x               &
          + (5.0_DP * c1 - 5.0_DP * c0 - 2.5_DP) *y*x )

dbdv = b0 * ( -121.0_DP * c0 / 3.0_DP * x4 * y + &
          (-27.0_DP * c1 + 81.0_DP * c0) * x4    &
        + ( 98.0_DP * c1 / 3.0_DP -49.0_DP * c0 - 49.0_DP / 6.0_DP)*x2*y*y &
        + ( 25.0_DP / 6.0_DP + 25.0_DP * (c0 - c1) / 3.0_DP ) * x2 * y ) /v0


dbulk=- vm * dbdv / bulk

dbd2v = b0 * ( 1694._DP / 9.0_DP * c0 * x2 * x3 * y             &
           + 324.0_DP * ( c1 / 3.0_DP - c0) * x2 * x3           &
           + (-980.0_DP / 9.0_DP * c1 + 490.0_DP * c0 / 3.0_DP  &
           + 245.0_DP / 9.0_DP ) * x3 * y * y                   &
           + (200.0_DP / 9.0_DP * ( c1 - c0 ) - 100.0_DP/9.0_DP)*x3*y)/v0**2

!
!d2bulk= vm * dbdv / bulk **2 - vm**2 * dbdv**2 / bulk**3 + &
!                               (vm/bulk)**2 * dbd2v

d2bulk=  - (1.0_DP + dbulk) * dbulk / bulk + (vm/bulk)**2 * dbd2v

RETURN
END SUBROUTINE birch4_bulk
!
!-------------------------------------------------------------------
SUBROUTINE birch4_bulk_pol(v0, b0, b01, b02, vm, a, m1, bulk, &
                                                 dbulk, d2bulk)
!-------------------------------------------------------------------
!
USE polyfit_mod, ONLY : compute_poly_deriv2, compute_poly_deriv3, &
                        compute_poly_deriv4
!
IMPLICIT NONE
!
INTEGER,  INTENT(IN)  :: m1
REAL(DP), INTENT(IN)  :: v0, b0, b01, b02, vm, a(m1)
REAL(DP), INTENT(OUT) :: bulk, dbulk, d2bulk

REAL(DP) :: x, x2, x3, x4, c0, c1, y, dbdv, dbd2v, aux, aux1, aux2

x=v0/vm
x2=x*x
x3=x2*x
x4=x2*x2
y=x2**(1.0_DP/3.0_DP)
c0 = (9.d0*b0*b02 + 9.d0*b01**2 - 63.d0*b01 + 143.d0)/48.d0
c1 = 3.d0*(b01-4.d0)/8.d0

bulk = b0 * ( 11.0_DP * c0 * x3 * y + (9.0_DP * c1 - 27.0_DP *c0) * x3 &
          + (- 14.0_DP * c1 +21.0_DP * c0 +3.5_DP) *y*y*x               &
          + (5.0_DP * c1 - 5.0_DP * c0 - 2.5_DP) *y*x )

CALL compute_poly_deriv2(vm, m1-1, a, aux)
bulk = bulk + aux * vm

dbdv = b0 * ( -121.0_DP * c0 / 3.0_DP * x4 * y + &
          (-27.0_DP * c1 + 81.0_DP * c0) * x4  + &
          ( 98.0_DP * c1 / 3.0_DP -49.0_DP * c0 - 49.0_DP / 6.0_DP)*x2*y*y + &
          ( 25.0_DP / 6.0_DP + 25.0_DP * (c0 - c1) / 3.0_DP ) * x2 * y ) /v0

CALL compute_poly_deriv3(vm, m1-1, a, aux1)
dbdv = dbdv + aux + aux1 * vm 

dbulk = - vm * dbdv / bulk

dbd2v = b0 * ( 1694._DP / 9.0_DP * c0 * x2 * x3 * y             &
           + 324.0_DP * ( c1 / 3.0_DP - c0) * x2 * x3           &
           + (-980.0_DP / 9.0_DP * c1 + 490.0_DP * c0 / 3.0_DP  &
           + 245.0_DP / 9.0_DP ) * x3 * y * y                   &
           + (200.0_DP / 9.0_DP * (c1 - c0) - 100.0_DP/9.0_DP)*x3*y)/v0**2

CALL compute_poly_deriv4(vm, m1-1, a, aux2)

dbd2v= dbd2v + 2.0_DP * aux1 + vm * aux2

!
!d2bulk= vm * dbdv / bulk **2 - vm**2 * dbdv**2 / bulk**3 + &
!                               (vm/bulk)**2 * dbd2v
!
d2bulk=  - (1.0_DP + dbulk) * dbulk / bulk + (vm/bulk)**2 * dbd2v

RETURN
END SUBROUTINE birch4_bulk_pol

END MODULE eos
