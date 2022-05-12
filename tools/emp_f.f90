!
! Copyright (C) 2022 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM emp_f
!
!   This program thermodynamic properties of W starting from the 
!   empirical expression of the free energy given in the paper
!   P.I. Dorogokupets and A.R. Oganov, Phys. Rev. B 75, 024115 (2007).
!   or those given in the paper
!   K.D. Litasov et al. Journ. of Appl. Phys. 113, 093507 (2013).
!   The parameters for the different elements are given in input.
!   There is also the possibility to use the Mie-Gruneisen-Debye model
!   (See the same paper of Litasov et. al.)
!
USE kinds,       ONLY : DP
USE constants,   ONLY : avogadro, bohr_radius_angs, rydberg_si, ry_kbar
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end
USE thermodynamics_mod, ONLY : p_from_f, b_from_p, entropy_from_f, &
                        energy_from_f, f_from_pressure
USE debye_module,ONLY : debye_vib_energy
USE eos,         ONLY : eos_energy
USE io_global,   ONLY : stdout
IMPLICIT NONE
CHARACTER(LEN=9) :: code="EMP_F"
CHARACTER(LEN=256) :: filename, string
CHARACTER(LEN=8) :: float_to_char
INTEGER :: nvol, ntemp, npress, ntemp_plot, nvol_plot, npress_plot

REAL(DP) :: v0, b0, b01, b02, vmin, vmax, deltav, x, y, eta,            &
            t, r, es, hs, as, ms, expb, v1, v2, p1, p2, b1, b2,         &
            bs1, bs2, cv1, cv2, cp1, cp2, alp1, alp2, gam1, gam2,       &
            q, t0, alpha0, alpha1, aux0, temp0
REAL(DP) :: mb(2), thetab(2), db(2), b(2), g, me(2), thetae(2),         &
            g0, gi, beta, vs, theta, db0dt
REAL(DP), ALLOCATABLE :: v_t(:), alpha_t(:), b0_t(:), b0_s(:), cv_t(:), &
            cp_t(:), gamma_t(:), gam(:), k0(:), auxv(:)
REAL(DP), ALLOCATABLE :: v(:), e(:)
REAL(DP), ALLOCATABLE :: temp(:), vol(:), vol_plot(:), temp_plot(:),    &
            f(:,:), p(:,:), bp(:,:), cv(:,:), abt(:,:), betat(:,:),     &
            gammat(:,:), entr(:,:), ener(:,:), cp(:,:), bs(:,:),        &
            press_plot(:), aux(:), aux1(:), aux2(:), p_aux(:,:)
REAL(DP), ALLOCATABLE :: v_p(:), b0_p(:), b0_sp(:), cv_p(:), cp_p(:),   &
            alpha_p(:), gamma_p(:), press(:)
REAL(DP) :: tmin, tmax, deltat, dt, dtmin
REAL(DP) :: pmin, pmax, deltap
REAL(DP) :: ep, ge, ne, dt_min, dv_min, vau, vtau, k0kbar
INTEGER  :: i, i1, itemp, itempp, ivol, ivolp, which_f, ipressp, iu_p, it, iv,&
            itemp300, ipress, leng
INTEGER, ALLOCATABLE  :: itemp_plot(:), ivol_plot(:)
!

CALL mp_startup ( start_images=.TRUE. )
CALL environment_start ( code )

WRITE(stdout,'(5x,"Dorogokupets (1), Litasov (2), MGD (3) or  &
                                             & HTBM (4) free energy?")') 
READ(5,*) which_f
IF (which_f < 1 .OR. which_f > 4 ) CALL errore('emp_f','unknown free energy',1) 
WRITE(6,'(5x,"Volume (m^3), bulk modulus b0 (kbar), &
                                &pressure derivative of b0?")') 
READ(5,*) v0, b0, b01

b0=b0*1.D8
eta=1.5_DP * ( b01 - 1.0_DP )

IF (which_f==1) THEN
   WRITE(6,'(5x,"mb(1), d(2), thetab(1)?")') 
   READ(5,*) mb(1), db(1), thetab(1)

   WRITE(6,'(5x,"mb(2), d(2), thetab(2)?")') 
   READ(5,*) mb(2), db(2), thetab(2)
ENDIF

IF (which_f==1.OR.which_f==2) THEN

   WRITE(6,'(5x,"me(1), thetae(1)?")') 
   READ(5,*) me(1), thetae(1)

   WRITE(6,'(5x,"me(2), thetae(2)?")') 
   READ(5,*) me(2), thetae(2)

   WRITE(6,'(5x,"gamma_0, gamma_infty, beta?")') 
   READ(5,*) g0, gi, beta

   WRITE(6,'(5x,"a, m?")') 
   READ(5,*) as, ms

   WRITE(6,'(5x,"e, g?")') 
   READ(5,*) ep, ge

   IF (which_f==1) THEN
      WRITE(6,'(5x,"H, S?")') 
      READ(5,*) hs, es
   ENDIF
ELSEIF (which_f==3) THEN
   WRITE(6,'(5x,"gamma_0, gamma_infinity and q (or beta)?")') 
   READ(5,*)  g0, gi, q
   WRITE(6,'(5x,"theta_0?")') 
   READ(5,*)  t0
ELSEIF (which_f==4) THEN
   WRITE(6,'(5x,"thermal expansion (1/K)?")') 
   READ(5,*)  alpha0
   WRITE(6,'(5x,"temperature derivative of the thermal expansion(1/K^2)?")') 
   READ(5,*)  alpha1
   WRITE(6,'(5x,"temperature derivative of b0 (Pa/K)?")') 
   READ(5,*)  db0dt
ENDIF
!
!  The gas constant
!
r=8.3144626181532_DP
!
! Temperature mesh 
!
WRITE(6,'(5x,"Tmin, Tmax, ntemp")') 
READ(5,*) tmin, tmax, ntemp

ALLOCATE(temp(ntemp))
ALLOCATE(v_t(ntemp))
ALLOCATE(k0(ntemp))
ALLOCATE(alpha_t(ntemp)) 
ALLOCATE(b0_t(ntemp)) 
ALLOCATE(b0_s(ntemp)) 
ALLOCATE(cv_t(ntemp))
ALLOCATE(cp_t(ntemp))
ALLOCATE(gamma_t(ntemp))

deltat= (tmax - tmin) / (ntemp - 1.0_DP)

DO itemp=1,ntemp
   temp(itemp)= tmin + deltat * (itemp-1.0_DP)
ENDDO

itemp300=ntemp
dtmin=1.D10
DO itemp=1, ntemp
   dt=ABS(temp(itemp) - 298.15D0)
   IF (dt<dtmin) THEN
      itemp300=itemp
      dtmin=dt
   ENDIF
ENDDO
IF (itemp300==ntemp) itemp300=0

WRITE(6,'(5x,"Number of temperatures to plot")') 
READ(5,*) ntemp_plot

ALLOCATE(temp_plot(ntemp_plot))
ALLOCATE(itemp_plot(ntemp_plot))

DO itemp=1,ntemp_plot
   WRITE(6,'(5x,"Temperature ",i5)') itemp 
   READ(5,*) temp_plot(itemp)
ENDDO

DO itempp=1,ntemp_plot
   dt_min=1.D20
   DO itemp=1,ntemp
      IF (ABS(temp(itemp)-temp_plot(itempp)) < dt_min ) THEN
         dt_min=ABS(temp(itemp)-temp_plot(itempp))
         itemp_plot(itempp)=itemp
      ENDIF
   ENDDO
ENDDO

!DO itempp=1,ntemp_plot
!   WRITE(6,*) itempp, temp(itemp_plot(itempp))
!ENDDO
!
! Volume mesh 
!
WRITE(6,'(5x,"Vmin, Vmax, nvol")') 
READ(5,*) vmin, vmax, nvol

ALLOCATE(vol(nvol))
ALLOCATE(v(nvol))
ALLOCATE(e(nvol))

vmin=vmin*v0
vmax=vmax*v0
deltav= (vmax - vmin ) / (nvol - 1.0_DP)

DO ivol=1,nvol
   v(ivol)= vmin + deltav * (ivol-1.0_DP)
ENDDO

WRITE(6,'(5x,"Pmin, Pmax, npress (kbar)")') 
READ(5,*) pmin, pmax, npress

pmin=pmin*1.D8
pmax=pmax*1.D8

ALLOCATE(press(npress))

deltap= (pmax - pmin ) / (npress - 1.0_DP)

DO ipress=1,npress
   press(ipress)= pmin + deltap * (ipress-1.0_DP)
ENDDO

ALLOCATE(f(nvol,ntemp))
ALLOCATE(p(nvol,ntemp))
ALLOCATE(bp(nvol,ntemp))
ALLOCATE(cv(nvol,ntemp))
ALLOCATE(abt(nvol,ntemp))
ALLOCATE(betat(nvol,ntemp))
ALLOCATE(gammat(nvol,ntemp)) 
ALLOCATE(entr(nvol,ntemp)) 
ALLOCATE(ener(nvol,ntemp)) 
ALLOCATE(cp(nvol,ntemp)) 
ALLOCATE(bs(nvol,ntemp))
ALLOCATE(p_aux(nvol,ntemp))
ALLOCATE(gam(nvol))
ALLOCATE(aux(ntemp))
ALLOCATE(auxv(nvol))
ALLOCATE(aux1(ntemp))
ALLOCATE(aux2(ntemp))
ALLOCATE(v_p(npress))
ALLOCATE(b0_p(npress))
ALLOCATE(b0_sp(npress))
ALLOCATE(cv_p(npress))
ALLOCATE(cp_p(npress))
ALLOCATE(alpha_p(npress))
ALLOCATE(gamma_p(npress))

WRITE(6,'(5x,"Number of volumes to plot")') 
READ(5,*) nvol_plot

ALLOCATE(vol_plot(nvol_plot))
ALLOCATE(ivol_plot(nvol_plot))

DO ivol=1,nvol_plot
   WRITE(6,'(5x,"Volume ",i5)') ivol
   READ(5,*) vol_plot(ivol)
   vol_plot(ivol)=vol_plot(ivol)*v0
ENDDO

DO ivolp=1,nvol_plot
   dv_min=1.D20
   DO ivol=1,nvol
      IF (ABS(v(ivol)-vol_plot(ivolp)) < dv_min ) THEN
         dv_min = ABS(v(ivol)-vol_plot(ivolp))
         ivol_plot(ivolp)=ivol
      ENDIF
   ENDDO
ENDDO
!
!
WRITE(6,'(5x,"Number of pressures to plot")') 
READ(5,*) npress_plot

ALLOCATE(press_plot(npress_plot))

DO ipressp=1,npress_plot
   WRITE(6,'(5x,"Pressure ",i5," in kbar")') ipressp
   READ(5,*) press_plot(ipressp)
   press_plot(ipressp)=press_plot(ipressp)*1.D8
ENDDO

DO i=1, nvol
   v(i)= vmin + (i-1) * deltav
   x= v(i) / v0
   y= x**(1.0_DP / 3.0_DP)
   e(i)= 9.0_DP * b0 * v0 *(1.0_DP - (1.0_DP - eta * (1.0_DP - y)) &
                                         * exp(eta * (1.0_DP-y))) / eta**2
!   WRITE(6,'(3e20.9)') v(i), v(i)/v0, e(i)
ENDDO

IF (which_f==1.OR.which_f==2) THEN
   f=0.0_DP
   DO itemp=1, ntemp
      t=temp(itemp)
      DO i=1, nvol
         x= v(i) / v0
         y= x**(1.0_DP / 3.0_DP)
         vs=x**(-gi)*exp((g0-gi)/beta*(1.0_DP - x**beta))
         DO i1=1,2
            g = db(i1) * log ( 1.0_DP + thetab(i1) * vs / t / db(i1) )
            b(i1) = 1.0_DP / ( exp(g)-1.0_DP )
         ENDDO
         DO i1=1,2
!
!  quasi harmonic term
!
!
            IF (which_f==1) THEN
               f(i,itemp) = f(i,itemp) + mb(i1) * r * ( (db(i1)-1.0_DP) * &
                             thetab(i1)*vs / & 
                             2.0_DP / db(i1) - t * log(1.0_DP + b(i1)) )
               f(i,itemp) = f(i,itemp) + me(i1) * r * ( thetae(i1)*vs * 0.5_DP &
                       + t * log ( 1.0_DP - exp( -thetae(i1)*vs/t ) ))
            ELSEIF (which_f==2) THEN
               f(i,itemp) = f(i,itemp) + me(i1) * r * (  &
                         thetae(i1)*vs * 0.5_DP + &
                         t * log ( 1.0_DP - exp( -thetae(i1)*vs/t ) ))

            ENDIF
!
!  anharmonic term
!
            IF (which_f==1) THEN
               expb = exp( thetab(i1)* vs / t )
               f(i,itemp) = f(i,itemp) + mb(i1) * r * as * x**ms *(  &
                     (thetab(i1)* vs * (0.5_DP + 1.0_DP /      &
                     (expb-1.0_DP) ) )**2 +                    &
                     2.0_DP * (thetab(i1)*vs / t )**2 * expb * t**2 / &
                     (expb - 1.0_DP)**2 ) / 6.0_DP
                expb = exp( thetae(i1) * vs / t )
                f(i,itemp) = f(i,itemp) + me(i1) * r * as * x**ms *(  &
                    (thetae(i1)* vs * (0.5_DP + 1.0_DP /       &
                    (expb-1.0_DP) ) )**2 +                     &
                    2.0_DP * (thetae(i1)*vs / t )**2 * expb * t**2 / &
                    (expb - 1.0_DP)**2 ) / 6.0_DP
            ENDIF
         ENDDO
!
!   Anharmonic contribution for Litasov functional 
!
         IF (which_f==2) THEN
            f(i,itemp) = f(i,itemp) - 1.5_DP * r * as * x**ms * t**2
         ENDIF
!
!  Electronic contribution
!
!         f(i,itemp) = - 1.5_DP * r * ep * x**ge * t**2 
         f(i,itemp) = f(i,itemp) - 1.5_DP * r * ep * x**ge * t**2 
!
!  Defects contribution
!
         IF (which_f==1) THEN
            f(i,itemp) = f(i,itemp)     &
                    - 1.5_DP * r * t * exp (es / x - hs / t / x / x)
         ENDIF
      ENDDO
   ENDDO
ELSEIF (which_f==3) THEN
   f=0.0_DP
   DO i=1, nvol
      gam(i)=g0 * ( v(i) / v0 )**q
      theta = ( g0 - gam(i) ) / q
      theta = t0 * EXP (theta)
      CALL debye_vib_energy(theta, temp, ntemp, 1, aux)
      DO itemp=1,ntemp
         f(i, itemp) =avogadro * rydberg_si * aux(itemp)/3.0_DP &
                  + 3.0_DP * r * ( 3.0_DP * theta / 8.0_DP +   &
                  temp(itemp) *log(1.0_DP - exp(-theta/temp(itemp)))) 
      ENDDO
   ENDDO
ELSEIF (which_f==4) THEN
   DO itemp=1,ntemp
      temp0=300.0_DP
      k0(itemp)= b0 + db0dt * (temp(itemp)-temp0)
      v_t(itemp)= v0 * exp(alpha0 *(temp(itemp)-temp0) + alpha1 * 0.5_DP * &
                                             (temp(itemp)**2 - temp0 **2))
      DO i=1, nvol
!
!  ieos=1 means Birch-Murnaghan of third order
!
         vau = v(i) *1.D30/avogadro/bohr_radius_angs**3
         vtau = v_t(itemp) *1.D30/avogadro/bohr_radius_angs**3
         k0kbar=k0(itemp)/1.D8 /ry_kbar
         CALL eos_energy(1, vau, aux0, vtau, k0kbar, b01, b02)
         f(i,itemp)= aux0 * avogadro * rydberg_si
         WRITE(6,*) 'i, itemp, e(i), aux0', i, e(i), aux0 * avogadro * &
                                              rydberg_si
      ENDDO
   ENDDO
ENDIF
DO i=1,nvol
   auxv(i)=f(i,itemp300)
ENDDO
DO itemp=1,ntemp
   DO i=1,nvol
      f(i,itemp)=e(i) + f(i,itemp)-auxv(i)
   ENDDO
ENDDO
!
!   Compute pressure and bulk modulus as a function of volume
!
DO itemp=1, ntemp
   CALL p_from_f(f(1,itemp),v,nvol,p(1,itemp))
   CALL b_from_p(p(2,itemp),v(2),nvol-2,bp(2,itemp))
ENDDO
!
!   Compute entropy and internal energy as a function of volume and temperature
!
DO ivol=1, nvol
   aux(1:ntemp)=f(ivol,1:ntemp)
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
ENDDO
!
!  Write the computed quatities on output
!  Free energy as a function of V for a few T:
!
iu_p=58
OPEN(UNIT=iu_p, FILE='free_energy_v.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"V/V0",11x," T=",8(f8.1," K",6x,"  T="))') &
                                  (temp_plot(it),it=1,ntemp_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))
DO ivol=3,nvol-2
   WRITE(iu_p,'(8e20.10)') v(ivol)/v0, (f(ivol,itemp_plot(it)),&
                                                      it=1,ntemp_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Free energy as a function of T for a few volumes 
!
iu_p=58
OPEN(UNIT=iu_p, FILE='free_energy_t.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (f(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Energy
!
iu_p=58
OPEN(UNIT=iu_p, FILE='energy.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (ener(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Entropy
!
iu_p=58
OPEN(UNIT=iu_p, FILE='entropy.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (entr(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  C_V
!
iu_p=58
OPEN(UNIT=iu_p, FILE='cv.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (cv(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  beta
!
iu_p=58
OPEN(UNIT=iu_p, FILE='beta.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))


DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (betat(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  gamma
!
iu_p=58
OPEN(UNIT=iu_p, FILE='gamma.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (gammat(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  beta B_T
!
iu_p=58
OPEN(UNIT=iu_p, FILE='beta_bt.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (abt(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  C_P
!
iu_p=58
OPEN(UNIT=iu_p, FILE='cp.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (cp(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  B_T
!
iu_p=58
OPEN(UNIT=iu_p, FILE='bt.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (bp(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  B_S
!
iu_p=58
OPEN(UNIT=iu_p, FILE='bs.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(iv)*1.D30/avogadro/bohr_radius_angs**3,iv=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (bs(ivol_plot(iv),itemp), &
                                                        iv=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Pressure:
!
iu_p=58
OPEN(UNIT=iu_p, FILE='press.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"V/V0",11x," T=",8(f8.1," K",6x,"  T="))') &
                                  (temp_plot(it),it=1,ntemp_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO i=3,nvol-2
   WRITE(iu_p,'(8e20.10)') v(i)/v0, (p(i,itemp_plot(it)),it=1,ntemp_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Bulk modulus:
!
OPEN(UNIT=iu_p, FILE='bulk.dat', STATUS='unknown', FORM='formatted')

DO i=3,nvol-2
   WRITE(iu_p,'(8e20.10)') v(i)/v0, (bp(i,itemp_plot(it)),it=1,ntemp_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!   Compute v_t(itemp) at zero pressure and interpolates the other quantities
!   at this volume.
!
DO ipressp=1,npress_plot
   filename='const_p.dat.'//TRIM(float_to_char(press_plot(ipressp)/1.D8,1))
   OPEN(UNIT=iu_p, FILE=filename, STATUS='unknown', FORM='formatted')
   WRITE(iu_p,'("#t, v_t, b0_t, b0_s, cv_t, cp_t, alpha_t, gamma_t")')
   DO itemp=3,ntemp-2
      DO i=3,nvol-2
         IF (p(i,itemp)>press_plot(ipressp)) THEN
            p1=p(i,itemp)
            b1=bp(i,itemp)
            bs1=bs(i,itemp)
            cv1=cv(i,itemp)
            cp1=cp(i,itemp)
            alp1=betat(i,itemp)
            gam1=gammat(i,itemp)
            v1=v(i)
         ELSEIF (p(i-1,itemp)>press_plot(ipressp)) THEN
            p2=p(i,itemp)
            b2=bp(i,itemp)
            bs2=bs(i,itemp)
            cv2=cv(i,itemp)
            cp2=cp(i,itemp)
            alp2=betat(i,itemp)
            gam2=gammat(i,itemp)
            v2=v(i)
         ENDIF
      ENDDO
      v_t(itemp)=v1+(press_plot(ipressp)-p1) * (v2 - v1) / (p2 - p1)
      b0_t(itemp)=b1+( b2 - b1 )* (v_t(itemp) - v1) / (v2 - v1)
      b0_s(itemp)=bs1+( bs2 - bs1 )* (v_t(itemp) - v1) / (v2 - v1)
      cv_t(itemp)=cv1+( cv2 - cv1 )* (v_t(itemp) - v1) / (v2 - v1)
      cp_t(itemp)=cp1+( cp2 - cp1 )* (v_t(itemp) - v1) / (v2 - v1)
      alpha_t(itemp)=alp1+( alp2 - alp1 )* (v_t(itemp) - v1) / (v2 - v1)
      gamma_t(itemp)=gam1+( gam2 - gam1 )* (v_t(itemp) - v1) / (v2 - v1)
      WRITE(iu_p,'(8e20.10)') temp(itemp), v_t(itemp), b0_t(itemp), &
                        b0_s(itemp), cv_t(itemp), cp_t(itemp), &
                        alpha_t(itemp), gamma_t(itemp)
   ENDDO
   CLOSE(UNIT=iu_p, STATUS='keep')
ENDDO

DO itempp=1,ntemp_plot
   filename='const_t.dat.'//TRIM(float_to_char(temp_plot(itempp),1))
   OPEN(UNIT=iu_p, FILE=filename, STATUS='unknown', FORM='formatted')
   WRITE(iu_p,'("#p, v_p, b0_p, b0_sp, cv_p, cp_p, alpha_p, gamma_p")')
   itemp=itemp_plot(itempp)
   DO ipress=1,npress
      DO i=3,nvol-2
         IF (p(i,itemp)>press(ipress)) THEN
            p1=p(i,itemp)
            b1=bp(i,itemp)
            bs1=bs(i,itemp)
            cv1=cv(i,itemp)
            cp1=cp(i,itemp)
            alp1=betat(i,itemp)
            gam1=gammat(i,itemp)
            v1=v(i)
         ELSEIF (p(i-1,itemp)>press(ipress)) THEN
            p2=p(i,itemp)
            b2=bp(i,itemp)
            bs2=bs(i,itemp)
            cv2=cv(i,itemp)
            cp2=cp(i,itemp)
            alp2=betat(i,itemp)
            gam2=gammat(i,itemp)
            v2=v(i)
         ENDIF
      ENDDO
      v_p(ipress)=v1+(press(ipress)-p1) * (v2 - v1) / (p2 - p1)
      b0_p(ipress)=b1+( b2 - b1 )* (v_p(ipress) - v1) / (v2 - v1)
      b0_sp(ipress)=bs1+( bs2 - bs1 )* (v_p(ipress) - v1) / (v2 - v1)
      cv_p(ipress)=cv1+( cv2 - cv1 )* (v_p(ipress) - v1) / (v2 - v1)
      cp_p(ipress)=cp1+( cp2 - cp1 )* (v_p(ipress) - v1) / (v2 - v1)
      alpha_p(ipress)=alp1+( alp2 - alp1 )* (v_p(ipress) - v1) / (v2 - v1)
      gamma_p(ipress)=gam1+( gam2 - gam1 )* (v_p(ipress) - v1) / (v2 - v1)
      WRITE(iu_p,'(8e20.10)') press(ipress)/1.D8, v_p(ipress), b0_p(ipress), &
                        b0_sp(ipress), cv_p(ipress), cp_p(ipress), &
                        alpha_p(ipress), gamma_p(ipress)
   ENDDO
   CLOSE(UNIT=iu_p, STATUS='keep')
ENDDO

DEALLOCATE(p_aux)
DEALLOCATE(v)
DEALLOCATE(e)
DEALLOCATE(ivol_plot)
DEALLOCATE(itemp_plot)
DEALLOCATE(temp_plot)
DEALLOCATE(press_plot)
DEALLOCATE(temp)
DEALLOCATE(v_t)
DEALLOCATE(alpha_t) 
DEALLOCATE(b0_t) 
DEALLOCATE(b0_s) 
DEALLOCATE(cv_t)
DEALLOCATE(cp_t)
DEALLOCATE(gamma_t)
DEALLOCATE(vol)

DEALLOCATE(auxv)
DEALLOCATE(aux)
DEALLOCATE(aux1)
DEALLOCATE(aux2)
DEALLOCATE(gam)
DEALLOCATE(f)
DEALLOCATE(p)
DEALLOCATE(bp)
DEALLOCATE(cv)
DEALLOCATE(abt)
DEALLOCATE(betat)
DEALLOCATE(gammat) 
DEALLOCATE(entr) 
DEALLOCATE(ener) 
DEALLOCATE(cp) 
DEALLOCATE(bs)
DEALLOCATE(k0)
DEALLOCATE(v_p)
DEALLOCATE(b0_p)
DEALLOCATE(b0_sp)
DEALLOCATE(cv_p)
DEALLOCATE(cp_p)
DEALLOCATE(alpha_p)
DEALLOCATE(gamma_p)
DEALLOCATE(press)

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM emp_f
