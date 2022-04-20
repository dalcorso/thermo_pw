!
! Copyright (C) 2022 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE thermodynamics_mod
!
!  This module provide methods to calculate various thermodynamic quantities. 
!  Presently it provides:
!
!  entropy_from_f   ! receives as input the Helmholtz free energy as a 
!                   ! function of temperature and gives as output the entropy
!  energy_from_f    ! receives as input the Helmholtz free energy and
!                   ! entropy as a function of tempereture and gives as 
!                   ! output the internal energy (or receives Gibbs energy
!                   ! and entropy and gives as output enthalpy
!  f_from_entropy   ! receives as input the entropy and gives as output
!                   ! the Helmholtz free energy 
!  f_from_pressure  ! receives as input the pressure and gives as output
!                   ! the Helmholtz free energy 
!  hc_from_entropy  ! receives as input the entropy as a function of 
!                   ! temperature and gives a output the heat capacity
!  beta_from_v      ! receives as input the volume as a function of temperature
!                   ! and gives as output the thermal expansion
!  b_from_v         ! receives as input the volume as a function of pressure
!                   ! and gives as output the bulk modulus
!  b_from_p         ! receives as input the pressure as a function of volume
!                   ! and gives as output the bulk modulus
!  comp_from_v      ! receives as input the volume as a function of pressure
!                   ! and gives as output the compressibility
!  v_from_g         ! receives the Gibbs energy as a function of pressure
!                   ! and computes the volume as a function of pressure
!  p_from_f         ! receives the Helmholtz energy as a function of volume
!                   ! and computes the pressure as a function of volume
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
SAVE
PRIVATE

PUBLIC entropy_from_f, energy_from_f, f_from_entropy, f_from_pressure, &
       hc_from_entropy, beta_from_v, p_from_f, b_from_v, b_from_p, comp_from_v

CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE entropy_from_f(hfe,temp,ntemp,entr)
!-----------------------------------------------------------------------
!
!  Note: the entropy in the first and last point is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: hfe(ntemp), temp(ntemp)
REAL(DP), INTENT(OUT) :: entr(ntemp)

INTEGER :: i

entr=0.0_DP
DO i=2,ntemp-1
   entr(i) = - ( hfe(i+1)-hfe(i-1) ) / (temp(i+1) - temp(i-1))
ENDDO

RETURN
END SUBROUTINE entropy_from_f

!-----------------------------------------------------------------------
SUBROUTINE energy_from_f(hfe,temp,ntemp,entr,ener)
!-----------------------------------------------------------------------
!
!  Note: the energy in the first and last point is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: hfe(ntemp), temp(ntemp), entr(ntemp)
REAL(DP), INTENT(OUT) :: ener(ntemp)

INTEGER :: i

ener=0.0_DP
DO i=2,ntemp-1
   ener(i) = hfe(i) + temp(i) * entr(i)
ENDDO

RETURN
END SUBROUTINE energy_from_f

!-----------------------------------------------------------------------
SUBROUTINE f_from_entropy(hfe,temp,ntemp,entr)
!-----------------------------------------------------------------------
!
!  Note: the entropy in the first and last point is not computed
!  This routine provides f-f_0, where f is the Helmholtz free energy
!  if the entropy is a fixed volume, the Gibbs free energy if entropy
!  is at fixed pressure.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), entr(ntemp)
REAL(DP), INTENT(OUT) :: hfe(ntemp)

INTEGER :: i

hfe(1)=0.0_DP
hfe(2)=0.0_DP
DO i=2,ntemp-1
   hfe(i+1) = hfe(i-1) - entr(i) * (temp(i+1)-temp(i-1))
ENDDO

RETURN
END SUBROUTINE f_from_entropy
!
!-----------------------------------------------------------------------
SUBROUTINE f_from_pressure(hfe,vol,nvol,press)
!-----------------------------------------------------------------------
!
!  Note: the pressure in the first and last point is not computed
!  This routine provides f-f_0, where f is the Helmholtz free energy
!  if the pressure is at fixed temperature.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvol
REAL(DP), INTENT(IN) :: vol(nvol), press(nvol)
REAL(DP), INTENT(OUT) :: hfe(nvol)

INTEGER :: i

hfe(1)=0.0_DP
DO i=2,nvol-1
   hfe(i) = hfe(i-1) - press(i) * (vol(i+1)-vol(i-1))
ENDDO

RETURN
END SUBROUTINE f_from_pressure

!-----------------------------------------------------------------------
SUBROUTINE hc_from_entropy(hc,temp,ntemp,entr)
!-----------------------------------------------------------------------
!
!  Note: the entropy in the first and last point is not computed
!        and the heat capacity is not computed in the first and
!        last two points
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), entr(ntemp)
REAL(DP), INTENT(OUT) :: hc(ntemp)

INTEGER :: i

hc(1:2)=0.0_DP
DO i=3,ntemp-2
   hc(i) = (entr(i+1) - entr(i-1)) / (temp(i+1)-temp(i-1)) * temp(i)
ENDDO

RETURN
END SUBROUTINE hc_from_entropy
!
!-----------------------------------------------------------------------
SUBROUTINE beta_from_v(v,temp,ntemp,beta)
!-----------------------------------------------------------------------
!
!  Note: the thermal expansion in the first and last point is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), v(ntemp)
REAL(DP), INTENT(OUT) :: beta(ntemp)

INTEGER :: i

beta=0.0_DP
DO i=2,ntemp-1
  beta(i) = (v(i+1) - v(i-1)) / (temp(i+1)-temp(i-1)) / v(i)
ENDDO

RETURN
END SUBROUTINE beta_from_v
!
!-----------------------------------------------------------------------
SUBROUTINE b_from_v(v,press,npress,b)
!-----------------------------------------------------------------------
!
!  Note: the bulk modulus in the first and last point is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), v(npress)
REAL(DP), INTENT(OUT) :: b(npress)

INTEGER :: i

b=0.0_DP
DO i=2,npress-1
  b(i) = -(v(i+1) - v(i-1)) / (press(i+1)-press(i-1)) / v(i)
  b(i) = 1.0_DP / b(i)
ENDDO

RETURN
END SUBROUTINE b_from_v
!
!-----------------------------------------------------------------------
SUBROUTINE b_from_p(press,v,nvol,b)
!-----------------------------------------------------------------------
!
!  Note: the bulk modulus in the first and last point is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvol
REAL(DP), INTENT(IN) :: press(nvol), v(nvol)
REAL(DP), INTENT(OUT) :: b(nvol)

INTEGER :: i

b=0.0_DP
DO i=2,nvol-1
  b(i) = -(press(i+1) - press(i-1)) / (v(i+1)-v(i-1)) * v(i)
ENDDO

RETURN
END SUBROUTINE b_from_p

!-----------------------------------------------------------------------
SUBROUTINE comp_from_v(v,press,npress,comp)
!-----------------------------------------------------------------------
!
!  Note: the compressibility in the first and last point is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), v(npress)
REAL(DP), INTENT(OUT) :: comp(npress)

INTEGER :: i

comp=0.0_DP
DO i=2,npress-1
  comp(i) = -(v(i+1) - v(i-1)) / (press(i+1)-press(i-1)) / v(i)
ENDDO

RETURN
END SUBROUTINE comp_from_v
!
!-----------------------------------------------------------------------
SUBROUTINE v_from_g(gibbs,press,npress,v)
!-----------------------------------------------------------------------
!
!  Note: the volume in the first and last point is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), gibbs(npress)
REAL(DP), INTENT(OUT) :: v(npress)

INTEGER :: i

v=0.0_DP
DO i=2,npress-1
  v(i) = -(gibbs(i+1) - gibbs(i-1)) / (press(i+1)-press(i-1)) 
ENDDO

RETURN
END SUBROUTINE v_from_g
!
!-----------------------------------------------------------------------
SUBROUTINE p_from_f(hf,v,nvol,press)
!-----------------------------------------------------------------------
!
!  Note: the pressure in the first and last volume is not computed
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvol
REAL(DP), INTENT(IN) :: v(nvol), hf(nvol)
REAL(DP), INTENT(OUT) :: press(nvol)

INTEGER :: i

press=0.0_DP
DO i=2,nvol-1
  press(i) = -(hf(i+1) - hf(i-1)) / (v(i+1)-v(i-1)) 
ENDDO

RETURN
END SUBROUTINE p_from_f

END MODULE thermodynamics_mod

