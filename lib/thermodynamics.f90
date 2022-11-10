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
!  thermo_from_hf   ! receives the free energy on a mesh of volumes and
!                   ! temperatures and gives on the same mesh
!                   ! the thermodynamic quantities
!  fit_at_p_target  ! transforms a function of F(V,T) into F(\bar P,T)
!                   ! at a target pressure \bar P. Needs p(V,T)
!  fit_at_p         ! transforms a function of F(V,T) into F(P,T).
!                   ! Needs P(V,T).
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
SAVE
PRIVATE

PUBLIC entropy_from_f, energy_from_f, f_from_entropy, f_from_pressure, &
       hc_from_entropy, beta_from_v, p_from_f, b_from_v, b_from_p, &
       comp_from_v, thermo_from_hf, fit_at_p_target, fit_at_p

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
!  Note: the bulk modulus in the first and last point is linearly 
!  extrapolated
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
b(1)=b(2) + (b(3) - b(2)) * (press(1)-press(2))/ (press(3)-press(2)) 
b(npress)=b(npress-1) + (b(npress-2) - b(npress-1)) *     &
                        (press(npress)-press(npress-1))   &
                      / (press(npress-2)-press(npress-1))
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
!
!-----------------------------------------------------------------------
SUBROUTINE thermo_from_hf(hf, v, nvol, temp, ntemp, p, bp, bs, & 
                          cv, cp, betat, gammat, entr, ener)
!-----------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER :: ntemp, nvol
REAL(DP) :: hf(nvol,ntemp), v(nvol), temp(ntemp) 
REAL(DP) :: p(nvol, ntemp), bp(nvol,ntemp), bs(nvol,ntemp), &
            cv(nvol, ntemp), cp(nvol,ntemp), betat(nvol, ntemp), &
            gammat(nvol,ntemp), entr(nvol,ntemp), ener(nvol,ntemp)

INTEGER :: itemp, ivol
REAL(DP) :: aux(ntemp), aux1(ntemp), aux2(ntemp), abt(nvol,ntemp)
REAL(DP) :: deltat

deltat=temp(3)-temp(2)

DO itemp=1, ntemp
   CALL p_from_f(hf(1,itemp),v,nvol,p(1,itemp))
   CALL b_from_p(p(2,itemp),v(2),nvol-2,bp(2,itemp))
ENDDO

DO ivol=1, nvol
   aux(1:ntemp)=hf(ivol,1:ntemp)
   CALL entropy_from_f(aux,temp,ntemp,aux1)
   CALL energy_from_f(aux,temp,ntemp,aux1,aux2)
   entr(ivol,1:ntemp)=aux1(1:ntemp)
   ener(ivol,1:ntemp)=aux2(1:ntemp)
   DO itemp=3,ntemp-2
      cv(ivol,itemp) = (ener(ivol,itemp+1) - ener (ivol,itemp-1))/2.0_DP/deltat
      abt(ivol,itemp) = ( p(ivol,itemp+1) - p(ivol,itemp-1) ) / 2.0_DP / deltat
      betat(ivol,itemp) = abt(ivol,itemp) / bp(ivol,itemp)
      gammat(ivol,itemp) = abt(ivol,itemp) * v(ivol) / cv(ivol,itemp)
      cp(ivol,itemp) = cv(ivol,itemp) + abt(ivol,itemp) * betat(ivol,itemp) &
                                      * temp(itemp) * v(ivol)
      bs(ivol,itemp) = bp(ivol,itemp) + abt(ivol,itemp) * temp(itemp) &
                                        * gammat(ivol, itemp)
   ENDDO
   cv(ivol,2)=2.0_DP * cv(ivol,3) - cv(ivol,4)
   cv(ivol,ntemp-1)=2.0_DP * cv(ivol,ntemp-2) - cv(ivol,ntemp-3)
   betat(ivol,2)=2.0_DP * betat(ivol,3) - betat(ivol,4)
   betat(ivol,ntemp-1)=2.0_DP * betat(ivol,ntemp-2) - betat(ivol,ntemp-3)
   gammat(ivol,2)=2.0_DP * gammat(ivol,3) - gammat(ivol,4)
   gammat(ivol,ntemp-1)=2.0_DP * gammat(ivol,ntemp-2) - gammat(ivol,ntemp-3)
   cp(ivol,2)=2.0_DP * cp(ivol,3) - cp(ivol,4)
   cp(ivol,ntemp-1)=2.0_DP * cp(ivol,ntemp-2) - cp(ivol,ntemp-3)
   bs(ivol,2)=2.0_DP * bs(ivol,3) - bs(ivol,4)
   bs(ivol,ntemp-1)=2.0_DP * bs(ivol,ntemp-2) - bs(ivol,ntemp-3)
   bp(ivol,2)=2.0_DP * bp(ivol,3) - bp(ivol,4)
   bp(ivol,ntemp-1)=2.0_DP * bp(ivol,ntemp-2) - bp(ivol,ntemp-3)
ENDDO

RETURN
END SUBROUTINE thermo_from_hf
!
!----------------------------------------------------------------------
SUBROUTINE fit_at_p_target(f, ntemp, v, nvol, p, p_target, vp, fp)
!----------------------------------------------------------------------
!
!  This function receives a function f of volume and temperature,
!  the pressure as a function of volume and temperature
!  and gives the function f at p_target as a function of 
!  temperature e il volume v(T) alla pression p_target
!
!  It is supposed that a linear interpolation of the values at the
!  two closest volumes is sufficient.
!
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, nvol

REAL(DP), INTENT(IN) :: f(nvol, ntemp), p(nvol,ntemp), v(nvol)
REAL(DP), INTENT(INOUT) :: fp(ntemp), vp(ntemp)
REAL(DP), INTENT(IN) :: p_target

INTEGER  :: itemp, i
REAL(DP) :: aux1, aux2, p1, p2, v1, v2

DO itemp=2,ntemp-1
   DO i=3,nvol-2
      IF (p(i,itemp)>p_target) THEN
         p1=p(i,itemp)
         aux1=f(i,itemp)
         v1=v(i)
      ELSEIF (p(i-1,itemp)>p_target) THEN
         p2=p(i,itemp)
         aux2=f(i,itemp)
         v2=v(i)
      ENDIF
   ENDDO
   vp(itemp)=v1+(p_target-p1) * (v2 - v1) / (p2 - p1)
   fp(itemp)=aux1+( aux2 - aux1 )* (vp(itemp) - v1) / (v2 - v1)
ENDDO

RETURN
END SUBROUTINE fit_at_p_target
!
!----------------------------------------------------------------------
SUBROUTINE fit_at_p(f, ntemp, v, nvol, p, press, npress, vp, fp, itemp)
!----------------------------------------------------------------------
!
!  This function receives a function f of volume and temperature,
!  the pressure as a function of volume and temperature
!  and a mesh of pressures press(npress) and gives  
!  the function f(npress,ntemp) at all pressures press(npress)
!
!  It is supposed that a linear interpolation of the values at the
!  two closest volumes is sufficient.
!
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, npress, nvol, itemp

REAL(DP), INTENT(IN) :: f(nvol, ntemp), p(nvol,ntemp), v(nvol), press(npress)
REAL(DP), INTENT(INOUT) :: fp(npress), vp(npress)

INTEGER  :: i, ipress
REAL(DP) :: p_target
REAL(DP) :: aux1, aux2, v1, v2, p1, p2

DO ipress=1,npress
   p_target=press(ipress)
   DO i=3,nvol-2
      IF (p(i,itemp)>p_target) THEN
         p1=p(i,itemp)
         aux1=f(i,itemp)
         v1=v(i)
      ELSEIF (p(i-1,itemp)>p_target) THEN
         p2=p(i,itemp)
         aux2=f(i,itemp)
         v2=v(i)
      ENDIF
   ENDDO
   vp(ipress)=v1+(p_target-p1) * (v2 - v1) / (p2 - p1)
   fp(ipress)=aux1+( aux2 - aux1 )* (vp(ipress) - v1) / (v2 - v1)
ENDDO

RETURN
END SUBROUTINE fit_at_p

END MODULE thermodynamics_mod

