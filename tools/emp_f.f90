!
! Copyright (C) 2022 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM emp_f
!
!   This program computes the thermodynamic properties of a few
!   materials from the empirical expression of the free energy 
!   given in the papers:
!   1) P.I. Dorogokupets and A.R. Oganov, Phys. Rev. B 75, 024115 (2007).
!   2) K.D. Litasov et al. Journ. of Appl. Phys. 113, 093507 (2013).
!   3) T.S. Sokolova, P.I. Dorogokupets, K.D. Litasov, Russian Geology and
!   geophysics 54, 181 (2013).
!   4) P.I. Dorogokupets, T.S. Sokolova, B.S. Danilov, K.D. Litasov
!      Geodynamics and tectonophysics 3, 129 (2012).
!
!   Other set of parameters can be found in 
!   5) K.D. Litasov et al. Journ. of Appl. Phys. 113, 133505 (2013).
!
!   The T=0 K equation of state is given by ieos with the following 
!   implemented options
!   ieos
!    5           Vinet equation.
!    6           Kunc equation. kappa must be given in input.
!    7           Holzapfel equation. n and Z must be given in input.
!    The thermal free energy is given by
!    which_f
!      1        Ref. 1)
!      2        Ref. 2)
!      3        Ref. 3)
!      4        High temperature Birch-Murnaghan equation. See Ref. 2)
!      5        Debye, Mie, Gruneisen model. See Ref. 2)         
!
USE kinds,       ONLY : DP
USE constants,   ONLY : avogadro, bohr_radius_angs, rydberg_si, ry_kbar, eps12
USE thermodynamics_mod, ONLY : p_from_f, b_from_p, entropy_from_f, &
                        energy_from_f, f_from_pressure
USE debye_module,ONLY : debye_vib_energy
USE eos,         ONLY : eos_energy
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end
USE io_global,   ONLY : stdout

IMPLICIT NONE
CHARACTER(LEN=9)   :: code="EMP_F"
CHARACTER(LEN=256) :: filename, string
CHARACTER(LEN=8)   :: float_to_char
!
!  Parameters of the functionals
!
REAL(DP) :: v0, b0, b01, b02,                  &   ! volume, bulk modulus
                                                   ! and pressure derivatives
            mb(2), thetab(2), db(2),           &   ! frequencies and weights 
            me(2), thetae(2),                  &   
            g0, gi, beta,                      &   ! Gruneisen parameter
            t0_p, delta0_p,                    &   ! t and delta for Gruneisen
            es, hs,                            &   ! defects contribution
            as, ms,                            &   ! anharmonic term
            ep, ge,                            &   ! electronic contribution   
            alpha0, alpha1,                    &   ! thermal expansion
            t0, q,                             &   ! Debye temperature and q
            theta, db0dt,                      &   ! theta and db0dt
            kappa, zeta, enne                      ! parameters of eos
!
!  Volume dependent variables, T=0 K
!
REAL(DP), ALLOCATABLE ::    &
            v(:),           &  ! volume
            e(:),           &  ! energy
            p0(:),          &  ! pressure
            k0(:),          &  ! bulk modulus
            k01(:),         &  ! pressure derivative of bulk modulus
            gamma0(:),      &  ! Gruneisen parameter
            auxv(:)            ! auxiliary variable 
!
!   Volume depedent variables, T given in input
!
REAL(DP), ALLOCATABLE ::    &
            p_v(:)            ! Pressure at the target temperature
!
!   Temperature dependent variables
!
REAL(DP), ALLOCATABLE ::    &
            temp(:),        &  ! temperature
            v_t(:),         &  ! volume 
            b0_t(:),        &  ! bulk modulus
            aux(:),         &  ! auxiliary arrays
            aux1(:),        &  ! 
            aux2(:)
!
!  Pressure dependent variables 
!
REAL(DP), ALLOCATABLE ::    &
            press(:)           ! pressure
!
!   Variables on a mesh of volumes and temperatures
!
REAL(DP), ALLOCATABLE ::    &
            f(:,:),         &  ! Helmholtz free energy
            entr(:,:),      &  ! entropy
            ener(:,:),      &  ! energy
            p(:,:),         &  ! pressure
            bp(:,:),        &  ! bulk modulus
            bs(:,:),        &  ! isoentropic bulk modulus
            cv(:,:),        &  ! isochoric heat capacity
            cp(:,:),        &  ! isobaric  heat capacity
            betat(:,:),     &  ! thermal expansion
            gammat(:,:),    &  ! Gruneisen parameter
            abt(:,:)           ! product of bulk modulus and thermal expansion
!
!   Control of volume, temperature, and pressure
!
REAL(DP) ::                          &
            vmin, vmax, deltav,      &
            tmin, tmax, deltat,      &
            pmin, pmax, deltap

INTEGER ::                           &
            nvol,                    &   ! number of volumes
            ntemp,                   &   ! number of temperatures
            npress                       ! number of pressures

INTEGER ::  itemp300                     ! temperature on the mesh closest to
                                         ! room temperature
!
!   Auxiliary for interpolation at a given volume
!
REAL(DP) ::                 &
            v1, v2,         &  !  Interpolation of the volume
            p1, p2,         &  !  pressure
            b1, b2,         &  !  isothermal bulk modulus
            bs1, bs2,       &  !  adiabatic bulk modulus
            cv1, cv2,       &  !  isochoric heat capacity
            cp1, cp2,       &  !  isobaric heat capacity
            bet1, bet2,     &  !  thermal expansion
            gam1, gam2,     &  !  Gruneisen parameter
            entr1, entr2       !  entropy
!
!  Values at the given target temperature for the two volumes closest
!  to vtarget, or to the volume corresponding to ptarget
!
REAL(DP) ::                      &
            p_tar(2),            &  ! pressure at the target volume
            bp_tar(2),           &  ! isothermal bulk modulus
            bs_tar(2),           &  ! isoentropic bulk modulus
            cv_tar(2),           &  ! isochoric heat capacity
            cp_tar(2),           &  ! isobaric heat capacity
            beta_tar(2),         &  ! thermal expansion
            gamma_tar(2),        &  ! Gruneisen parameter
            entr_tar(2)             ! entropy
!
!  Values at the given target temperature and volume, or pressure
!
REAL(DP) ::                 &
            ttarget,        &  ! temperature
            ptarget,        &  ! pressure
            vtarget,        &  ! volume 
            deltatar,       &  ! distance of the target temperature from mesh
            v_o,            &  ! volume at the target pressure
            p_o,            &  ! pressure at the target volume
            b0_o,           &  ! isothermal bulk modulus
            bs_o,           &  ! isoentropic bulk modulus
            cv_o,           &  ! isochoric heat capacity
            cp_o,           &  ! isobaric heat capacity
            beta_o,         &  ! thermal expansion
            gamma_o,        &  ! Gruneisen parameter
            entr_o             ! entropy
!
!   Control of the output
!
INTEGER::                            &
            nvol_plot,               &   ! number of volumes to plot
            ntemp_plot,              &   ! number of temperatures to plot
            npress_plot                  ! number of pressures to plot

INTEGER, ALLOCATABLE  ::             &   !
            ivol_plot(:),            &   ! volume to plot on the mesh
            itemp_plot(:)                ! temperature to plot on the mesh

REAL(DP), ALLOCATABLE ::             &
            vol_plot(:),             &   ! volumes to write on output
            temp_plot(:),            &   ! temperatures to write on output
            press_plot(:)                ! pressure to write on output
!
!  Control variables
!
INTEGER ::                           &
            which_f,                 &   ! choice of the functional 
            ieos,                    &   ! choice of the equation of state
            iwrite                       ! choose what to write on output
!
!   Auxiliary variables
!
REAL(DP) ::         &
            dt,     &  ! delta t
            dtmin,  &  ! minimum value of delta t
            dvmin,  &  ! minumum distance of volume to plot from volume mesh
            x, y, z, & ! auxiliary for the equation of state 
            eta,    &  !  
            c0, c2, &  ! 
            auxf, dauxdy, & !
            b(2), g,     &  !
            expat, vs,   &  !
            intg0, intg, &  ! integrals of the average Gruneisen parameter 
                            ! divided by V.
            t,        &     ! used to indicate the current temperature
            expb,     &     ! used to compute the anharmonic term
            r,        &     ! the gas constant
            aux0, temp0, &
            vau, vtau,   &  ! used to pass volumes in (a.u.)^3
            k0kbar,      &  ! used to pass bulk modulus in Ry/(a.u.)^3
            deltav_o        ! distance of the target volume from the mesh
!
!  Counters
!
INTEGER ::                    &
            i, j, i1,         &   !
            itemp, itempp,    &   !  counters on temperature
            ivol, ivolp,      &   !  counters on volume
            ipress, ipressp,  &   !  counters on pressure
            itemp_tar, i_tar      !

INTEGER ::  iu_p,             &   ! output unit
            leng                  ! string length
!
CALL mp_startup ( start_images=.TRUE. )
CALL environment_start ( code )
!
!    Read input variables
!    
!    Functional
!
WRITE(stdout,'(/,5x,"Which free energy do you want to use?")')
WRITE(stdout,'(7x,"1) Dorogokupets et al.")')
WRITE(stdout,'(7x,"2) Litasov et al.")')
WRITE(stdout,'(7x,"3) Sokolova et al.")') 
WRITE(stdout,'(7x,"4) High Tempertare Birch-Murnaghan")') 
WRITE(stdout,'(7x,"5) Mie-Gruneisen-Debye")') 
READ(5,*) which_f
IF (which_f < 1 .OR. which_f > 5) CALL errore('emp_f','unknown free energy',1) 
!
!    Cold equation of state
!
WRITE(stdout,'(5x,"Which cold equation of state?")') 
WRITE(stdout,'(7x,"5) Vinet ")') 
WRITE(stdout,'(7x,"6) Kunc")') 
WRITE(stdout,'(7x,"7) Holzapfel")') 
READ(5,*) ieos

IF (ieos==6) THEN
   WRITE(stdout,'(5x,"Value of k")') 
   READ(5,*) kappa
ELSEIF (ieos==7) THEN
   WRITE(stdout,'(5x,"Enne, zeta?")') 
   READ(5,*) enne, zeta
ENDIF
IF (ieos < 5 .OR. which_f > 7) CALL errore('emp_f',&
                                           'unknown equation of state',1) 
!
!   parameters common to all equations of state
!
WRITE(stdout,'(5x,"Volume (m^3), bulk modulus b0 (kbar), &
                                &pressure derivative of b0?")') 
READ(5,*) v0, b0, b01
!
!  In input b0 in kbar, here we use the SI units (Pa = N/m^2)
!
b0=b0*1.D8
!
!  parameters for the acoustic modes (only Ref. 1)
!
IF (which_f==1) THEN
   WRITE(stdout,'(5x,"mb(1), d(2), thetab(1)?")') 
   READ(5,*) mb(1), db(1), thetab(1)

   WRITE(stdout,'(5x,"mb(2), d(2), thetab(2)?")') 
   READ(5,*) mb(2), db(2), thetab(2)
ENDIF
!
!  parameters for optical modes (Refs. 1,2,3)
!
IF (which_f==1.OR.which_f==2.OR.which_f==3) THEN

   WRITE(stdout,'(5x,"me(1), thetae(1)?")') 
   READ(5,*) me(1), thetae(1)

   WRITE(stdout,'(5x,"me(2), thetae(2)?")') 
   READ(5,*) me(2), thetae(2)
!
!  Parameters for the volume dependence of the Gruneisen parameter 
!
   IF (which_f==3) THEN
      WRITE(stdout,'(5x,"t and delta?")') 
      READ(5,*) t0_p, delta0_p
   ELSE
      WRITE(stdout,'(5x,"gamma_0, gamma_infty, beta?")') 
      READ(5,*) g0, gi, beta
   ENDIF
!
!  Parameters for the anharmonic part (all references)
!
   WRITE(stdout,'(5x,"a, m?")') 
   READ(5,*) as, ms
!
!  Parameters for the electronic part (all references)
!
   WRITE(stdout,'(5x,"e, g?")') 
   READ(5,*) ep, ge
!
!  Parameters for the defect term (only Ref.1)
!
   IF (which_f==1) THEN
      WRITE(stdout,'(5x,"H, S?")') 
      READ(5,*) hs, es
   ENDIF
ELSEIF (which_f==4) THEN
!
!  The High Temperature Birch Murnaghan model requires in addition
!  the parameters to interpolate the thermal expansion and 
!  the derivative of the bulk modulus with respect to temperature.
!
   WRITE(stdout,'(5x,"thermal expansion (1/K)?")') 
   READ(5,*)  alpha0
   WRITE(stdout,'(5x,"temperature derivative of the thermal expansion(1/K^2)?")') 
   READ(5,*)  alpha1
   WRITE(stdout,'(5x,"temperature derivative of b0 (Pa/K)?")') 
   READ(5,*)  db0dt
ELSEIF (which_f==5) THEN
!
!   This is the Mie-Gruneisen-Debye model. It requires the parameters
!   for the Gruneisen parameter and the Debye temperature.
!
   WRITE(stdout,'(5x,"gamma_0, gamma_infinity and q (or beta)?")') 
   READ(5,*)  g0, gi, q
   WRITE(stdout,'(5x,"theta_0?")') 
   READ(5,*)  t0
ENDIF
!
!  The gas constant
!
r=8.3144626181532_DP
!
! Temperature mesh 
!
WRITE(stdout,'(5x,"Tmin, Tmax, ntemp")') 
READ(5,*) tmin, tmax, ntemp
!
!  Allocation of variables depending only on temperature
!
ALLOCATE(temp(ntemp))
ALLOCATE(v_t(ntemp))
ALLOCATE(b0_t(ntemp)) 
ALLOCATE(aux(ntemp))
ALLOCATE(aux1(ntemp))
ALLOCATE(aux2(ntemp))
!
!  mesh of temperatures
!
deltat= (tmax - tmin) / (ntemp - 1.0_DP)
DO itemp=1,ntemp
   temp(itemp)= tmin + deltat * (itemp-1.0_DP)
ENDDO
!
!  Find on the mesh the value of temperature closest to room temperature
!
itemp300=ntemp
dtmin=1.D10
DO itemp=1, ntemp
   dt=ABS(temp(itemp) - 298.15D0)
   IF (dt<dtmin) THEN
      itemp300=itemp
      dtmin=dt
   ENDIF
ENDDO
IF (itemp300==ntemp) CALL errore('emp_f','Temperature mesh must contain &
                                             room temperature',1)
!
!  Read at which temperatures to make the plot
!
WRITE(stdout,'(5x,"Number of temperatures to plot")') 
READ(5,*) ntemp_plot

ALLOCATE(temp_plot(ntemp_plot))
ALLOCATE(itemp_plot(ntemp_plot))

DO itemp=1,ntemp_plot
   WRITE(stdout,'(5x,"Temperature ",i5)') itemp 
   READ(5,*) temp_plot(itemp)
ENDDO
!
!  find the temperature closest to temp_plot and set the index itemp_plot
!
DO itempp=1,ntemp_plot
   dtmin=1.D20
   DO itemp=1,ntemp
      IF (ABS(temp(itemp)-temp_plot(itempp)) < dtmin ) THEN
         dtmin=ABS(temp(itemp)-temp_plot(itempp))
         itemp_plot(itempp)=itemp
      ENDIF
   ENDDO
ENDDO
!
DO itempp=1,ntemp_plot
   WRITE(stdout,'(5x,"Temp. #",i5," Requested T=",f15.5,&
                  &" K, Found T=",f15.5," K")') itempp, temp_plot(itempp), &
                                                temp(itemp_plot(itempp))
ENDDO
!
! Volume mesh 
!
WRITE(stdout,'(5x,"Vmin, Vmax, nvol")') 
READ(5,*) vmin, vmax, nvol
!
!  Allocate the variable depending only on the volume
!
ALLOCATE(v(nvol))
ALLOCATE(e(nvol))
ALLOCATE(p0(nvol))
ALLOCATE(k0(nvol))
ALLOCATE(k01(nvol))
ALLOCATE(gamma0(nvol))
ALLOCATE(auxv(nvol))
!
!  Compute the mesh of volumes. Note that in input vmin and vmax are in
!  units of v0.
!
vmin=vmin*v0
vmax=vmax*v0
deltav= (vmax - vmin ) / (nvol - 1.0_DP)

DO ivol=1,nvol
   v(ivol)= vmin + deltav * (ivol-1.0_DP)
ENDDO
!
!  Read the number of volumes to plot, which volumes and find 
!  the closest on the mesh of volumes.
!
WRITE(stdout,'(5x,"Number of volumes to plot")') 
READ(5,*) nvol_plot

ALLOCATE(vol_plot(nvol_plot))
ALLOCATE(ivol_plot(nvol_plot))

DO ivol=1,nvol_plot
   WRITE(stdout,'(5x,"Volume (in v0 units)",i5)') ivol
   READ(5,*) vol_plot(ivol)
   vol_plot(ivol)=vol_plot(ivol)*v0
ENDDO
!
!  find the volume closest to vol_plot and set the index ivol_plot
!
DO ivolp=1,nvol_plot
   dvmin=1.D20
   DO ivol=1,nvol
      IF (ABS(v(ivol)-vol_plot(ivolp)) < dvmin ) THEN
         dvmin = ABS(v(ivol)-vol_plot(ivolp))
         ivol_plot(ivolp)=ivol
      ENDIF
   ENDDO
ENDDO
!
! Plot the information on the volumes
!
DO ivolp=1,nvol_plot
   WRITE(stdout,'(5x,"Volume. #",i5," Requested V/V0=",f13.5,        &
                  &", Found V/V0=",f13.7)') itempp, vol_plot(ivolp)/v0, &
                                                v(ivol_plot(ivolp))/v0
ENDDO
!
!   Allocate the variables that depend on volume and temperature.
!
ALLOCATE(f(nvol,ntemp))
ALLOCATE(entr(nvol,ntemp)) 
ALLOCATE(ener(nvol,ntemp)) 
ALLOCATE(p(nvol,ntemp))
ALLOCATE(bp(nvol,ntemp))
ALLOCATE(bs(nvol,ntemp))
ALLOCATE(cv(nvol,ntemp))
ALLOCATE(cp(nvol,ntemp)) 
ALLOCATE(betat(nvol,ntemp))
ALLOCATE(gammat(nvol,ntemp)) 
ALLOCATE(abt(nvol,ntemp))
!
!   read the pressure mesh
!
WRITE(stdout,'(5x,"Pmin, Pmax, npress (kbar)")') 
READ(5,*) pmin, pmax, npress
!
!   1.D8 transforms kbar in Pa
!
pmin=pmin*1.D8
pmax=pmax*1.D8
!
!  Allocate pressure dependent variables
!
ALLOCATE(press(npress))
!
!  Compute the values of the pressure on the mesh
!
deltap= (pmax - pmin ) / (npress - 1.0_DP)
DO ipress=1,npress
   press(ipress)= pmin + deltap * (ipress-1.0_DP)
ENDDO
!
!  Read the number of pressure to plot, which pressures and find 
!  the closest on the mesh of pressures.
!
!
WRITE(stdout,'(5x,"Number of pressures to plot")') 
READ(5,*) npress_plot

ALLOCATE(press_plot(npress_plot))

DO ipressp=1,npress_plot
   WRITE(stdout,'(5x,"Pressure ",i5," in kbar")') ipressp
   READ(5,*) press_plot(ipressp)
   press_plot(ipressp)=press_plot(ipressp)*1.D8
ENDDO
!
!-----------------------------------------------------------------------
!  Now compute the T=0K energy as function of the volume
!
DO i=1, nvol
   v(i)= vmin + (i-1) * deltav
   x= v(i) / v0
   y= x**(1.0_DP / 3.0_DP)
!
!    Chose here the equation of state for fitting the energy
!
   IF (ieos==6) THEN
!
!    Kunc equation, kappa is given input. (kappa=2) is equivalent to 
!    Vinet equation
!
      eta=1.5_DP * b01 - kappa + 0.5_DP
      z=exp(eta * (1.0_DP-y))
      p0(i) = 3.0_DP * b0 * (1.0_DP - y) * z / y**kappa
      k0(i)=b0 *( kappa +(1.0_DP - kappa + eta)*y - eta * y**2 ) * z/y**kappa
      k01(i)= (kappa**2 + (2.0_DP*kappa *(1.0_DP+eta)-kappa**2-eta-1.0_DP)*y &
               +(eta**2-2.0_DP*kappa*eta+3.0_DP*eta)*y**2-eta**2 *y**3)/ &
                    (kappa +(1.0_DP - kappa + eta) * y - eta * y**2 ) / 3.0_DP
   ELSEIF (ieos==7) THEN
!
!   Holzapfel equation, enne and zeta required in input. enne is the number
!   of atoms per cell, Z is the atomic number. v0 in m^3/mol
!   In this formula b0 must be in GPa, since in the code it is in Pa
!
      c0=-log(3.0_DP * b0 / 1003.6_DP / (enne*zeta)**(5.0_DP/3.0_DP) / 1.D9 * &
                                              (v0*1.d6)**(5.0_DP/3.0_DP))
      c2=(b01-3.0_DP)*1.5_DP - c0
      z=exp(c0 * (1.0_DP-y))
      p0(i) = 3.0_DP * b0 * (1.0_DP - y) * (1.0_DP+c2*y*(1-y)) * z / y**5
      auxf=5.0_DP + (-4.0_DP + 4.0_DP*c2+c0)*y &
              + (-6.0_DP*c2-c0+c0*c2)*y**2 + (2.0_DP*c2-2.0_DP*c0*c2)*y**3 &
              +c0*c2*y**4
      dauxdy= -4.0_DP + 4.0_DP*c2+c0 + 2.0_DP *(-6.0_DP*c2-c0+c0*c2)*y + &
              3.0_DP * (2.0_DP*c2-2.0_DP*c0*c2)*y**2 + &
              4.0_DP * c0*c2*y**3
      k0(i) = b0 * auxf * z / y**5
      k01(i)= (5.0_DP + (c0 - dauxdy / auxf) * y ) /3.0_DP
   ELSE
!
!  Default: Vinet equation.
!
      eta=1.5_DP * (b01 - 1.0_DP)
      z=exp(eta * (1.0_DP-y))
      e(i)= 9.0_DP * b0 * v0 *(1.0_DP - (1.0_DP - eta * (1.0_DP - y)) &
                                         * z) / eta**2
      p0(i) = 3.0_DP * b0 * (1.0_DP - y) * z / y**2
      k0(i) = b0 *( 1.0_DP +(1.0_DP + eta * y) * (1.0_DP-y)) * z / y**2
      k01(i)= (2.0_DP + eta * y + (y *(1.0_DP - eta) + 2.0_DP * y**2 * eta) / &
                    (1.0_DP + (1.0_DP - y) * (1.0_DP+ y * eta))) / 3.0_DP
   ENDIF
!   WRITE(stdout,'(5e20.9)') v(i), v(i)/v0, k0(i), k01(i), e(i)
ENDDO
!
!  Compute the energy from the integral of the pressure if not
!  given as an analytic expression
!
IF (ieos == 6 .OR. ieos==7) THEN
   CALL f_from_pressure(e,v,nvol,p0)
   e=e/2.0_DP
ENDIF
!
!  the functional of Ref.3 requires the average Gruneisen parameter as
!  a function of the volume, computed from the T=0K EOS and two 
!  additional parameters.
!
IF (which_f==3) THEN
!
!   The original reference Burakovsky and Preston, J. Phys. Chem. Solid.
!   65, 1581 (2004) has a 3.0 in this formula, while
!   the paper of Sokolova et. al. has a 2.0, so comment and
!   uncomment the lines with the value that you prefer. 
!
   DO i=1, nvol
      gamma0(i)= (k01(i)*0.5_DP - 1.0_DP / 6.0_DP - t0_p *  &
               (1.0_DP-p0(i)/k0(i)/2.0_DP)/3.0_DP) /        &
!               (1.0_DP-p0(i)/k0(i)/3.0_DP)/3.0_DP) /        &
               (1.0_DP - 2.0_DP * t0_p * p0(i)              &
               /3.0_DP / k0(i)) + delta0_p
   ENDDO
!
!  Compute the integral of gamma0 from V_1 (the first volume of the mesh)
!  to V_0
!
   intg0=0.0_DP
   DO i=1,nvol
      IF (v(i) < (v0-eps12)) intg0=intg0+gamma0(i)*MIN(deltav,v0-v(i))/v(i)
   ENDDO
ENDIF
!
!   Now compute the free energy on the mesh of volumes and temperatures
!
IF (which_f==1.OR.which_f==2.OR.which_f==3) THEN
   f=0.0_DP
   expat=1.0_DP
   DO itemp=1, ntemp
      t=temp(itemp)
      DO i=1, nvol
         x= v(i) / v0
         y= x**(1.0_DP / 3.0_DP)
         IF (which_f==1.OR.which_f==2) THEN
            vs=x**(-gi)*exp((g0-gi)/beta*(1.0_DP - x**beta))
            DO i1=1,2
               g = db(i1) * log ( 1.0_DP + thetab(i1) * vs / t / db(i1) )
               b(i1) = 1.0_DP / ( exp(g)-1.0_DP )
            ENDDO
         ELSEIF (which_f==3) THEN
            IF (as /= 0.0_DP) expat = exp (0.5_DP * as * x**ms * t)
            intg=0.0_DP
            DO j=1,i-1
               intg=intg + gamma0(j) * deltav/ v(j)
            ENDDO
            vs = exp(-intg+intg0) * expat
         ENDIF
!
!  quasi harmonic term
!
         DO i1=1,2
            IF (which_f==1) THEN
               f(i,itemp) = f(i,itemp) + mb(i1) * r * ( (db(i1)-1.0_DP) * &
                             thetab(i1)*vs / & 
                             2.0_DP / db(i1) - t * log(1.0_DP + b(i1)) )
               f(i,itemp) = f(i,itemp) + me(i1)*r*( thetae(i1)*vs * 0.5_DP &
                       + t * log ( 1.0_DP - exp( -thetae(i1)*vs/t ) ))
            ELSEIF (which_f==2.OR.which_f==3) THEN
               f(i,itemp) = f(i,itemp) + me(i1) * r * (  &
                         thetae(i1)*vs * 0.5_DP + &
                         t * log ( 1.0_DP - exp( -thetae(i1)*vs/t ) ))
            ENDIF
!
!  anharmonic term for the functional of Ref.1).
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
!   Anharmonic contribution for the functionals of Ref.2) and 3) 
!
         IF (which_f==2.OR.which_f==3) THEN
            f(i,itemp) = f(i,itemp) - 1.5_DP * r * as * x**ms * t**2
         ENDIF
!
!  Electronic contribution (all functionals)
!
         f(i,itemp) = f(i,itemp) - 1.5_DP * r * ep * x**ge * t**2 
!
!  Defects contribution, only in Ref.1)
!
         IF (which_f==1) THEN
            f(i,itemp) = f(i,itemp)     &
                    - 1.5_DP * r * t * exp (es / x - hs / t / x / x)
         ENDIF
      ENDDO
   ENDDO
ELSEIF (which_f==4) THEN
!
!  This is high temperature Birch Murnaghan model
!
   DO itemp=1,ntemp
      temp0=300.0_DP
      b0_t(itemp)= b0 + db0dt * (temp(itemp)-temp0)
      v_t(itemp)= v0 * exp(alpha0 *(temp(itemp)-temp0) + alpha1 * 0.5_DP * &
                                             (temp(itemp)**2 - temp0 **2))
      DO i=1, nvol
!
!  The first argument of the routine eos_energy 1 means Birch-Murnaghan 
!  of third order. The routine requires volumes in (a.u.)^3 and
!  bulk modulus in Ry/(a.u.)^3. b01 is assumed indipendent from
!  temperature. b02 is not used.
!
         vau = v(i)*1.D30/avogadro/bohr_radius_angs**3
         vtau = v_t(itemp)*1.D30/avogadro/bohr_radius_angs**3
         k0kbar=b0_t(itemp)/1.D8 /ry_kbar
         CALL eos_energy(1, vau, aux0, vtau, k0kbar, b01, b02)

         f(i,itemp)= aux0 * avogadro * rydberg_si
!
!        WRITE(6,*) 'i, itemp, e(i), aux0', i, e(i), f(i,itemp) 
!
      ENDDO
   ENDDO
ELSEIF (which_f==5) THEN
!
!  This is the Mie-Gruneisen-Debye model
!
   f=0.0_DP
   DO i=1, nvol
      gamma0(i)=g0 * ( v(i) / v0 )**q
      theta = ( g0 - gamma0(i) ) / q
      theta = t0 * EXP (theta)
      CALL debye_vib_energy(theta, temp, ntemp, 1, aux)
      DO itemp=1,ntemp
         f(i, itemp) =avogadro * rydberg_si * aux(itemp)/3.0_DP &
                  + 3.0_DP * r * ( 3.0_DP * theta / 8.0_DP +    &
                  temp(itemp) *log(1.0_DP - exp(-theta/temp(itemp)))) 
      ENDDO
   ENDDO
ENDIF
!
!  Since the input parameters are at room temperature we subtract here
!  the free energy at room temperature
!
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

DO ivol=1, nvol
!
!   At each volume, compute entropy and internal energy for all temperatures
!   Then with a loop on temperature computes cv, the product of bulk
!   modulus and thermal expansion, then thermal expansion, average
!   gruneisen parameter, cp, and bs.
!
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
!-------------------------------------------------------------------------
!  Write the computed quatities on output
!  Free energy as a function of V for a few T:
!
iu_p=58
OPEN(UNIT=iu_p, FILE='free_energy_v.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"V/V0",11x," T=",8(f8.1," K",6x,"  T="))') &
                                  (temp_plot(itempp),itempp=1,ntemp_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))
DO ivol=3,nvol-2
   WRITE(iu_p,'(8e20.10)') v(ivol)/v0, (f(ivol,itemp_plot(itempp)),&
                                                      itempp=1,ntemp_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Free energy as a function of T for a few volumes 
!
iu_p=58
OPEN(UNIT=iu_p, FILE='free_energy_t.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3, &
              ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (f(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Energy as a function of temperature for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='energy.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3,         &
                     ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (ener(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Entropy as a function of temperature for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='entropy.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3,         &
              ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (entr(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  C_V as a function of T for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='cv.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3, &
                                                     ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (cv(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  thermal expansion as a function of T for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='beta.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3,&
                                                         ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (betat(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Average gruneisen parameter as a function of temperature for a few
!  volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='gamma.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3, &
                                                           ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (gammat(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  The product of thermal expansion and bulk modulus as a function of 
!  temperature for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='beta_bt.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3,&
                                                       ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (abt(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Isobaric heat capacity as a function of T for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='cp.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3, &
                                                          ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (cp(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Isothermal bulk modulus as a function of T for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='bt.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
             (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3, &
                                                       ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (bp(ivol_plot(ivolp),itemp), &
                                                     ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Isoentropic bulk modulus as a function of T for a few volumes
!
iu_p=58
OPEN(UNIT=iu_p, FILE='bs.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"T (K)",11x," V=",8(f8.1," (a.u.)^3 ",2x,"  V="))') &
        (vol_plot(ivolp)*1.D30/avogadro/bohr_radius_angs**3,ivolp=1,nvol_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO itemp=3,ntemp-2
   WRITE(iu_p,'(8e20.10)') temp(itemp), (bs(ivol_plot(ivolp),itemp), &
                                                        ivolp=1,nvol_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Pressure as a function of V/V0 for a few temperatures
!
iu_p=58
OPEN(UNIT=iu_p, FILE='press.dat', STATUS='unknown', FORM='formatted')

WRITE(string,'("#",8x,"V/V0",11x," T=",8(f8.1," K",6x,"  T="))') &
                                  (temp_plot(itempp),itempp=1,ntemp_plot)
leng=LEN(TRIM(string))-5
WRITE(iu_p,'(a)') TRIM(string(1:leng))

DO i=3,nvol-2
   WRITE(iu_p,'(8e20.10)') v(i)/v0, (p(i,itemp_plot(itempp)),&
                                                     itempp=1,ntemp_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!  Isothermal bulk modulus as a function of V/V0 for a few temperatures
!
OPEN(UNIT=iu_p, FILE='bulk.dat', STATUS='unknown', FORM='formatted')

DO i=3,nvol-2
   WRITE(iu_p,'(8e20.10)') v(i)/v0, (bp(i,itemp_plot(itempp)),&
                                                     itempp=1,ntemp_plot) 
ENDDO

CLOSE(UNIT=iu_p, STATUS='keep')
!
!   Compute v_t(itemp) for a set of pressures and interpolates the 
!   other quantities at this volume.
!
DO ipressp=1,npress_plot
   filename='const_p.dat.'//TRIM(float_to_char(press_plot(ipressp)/1.D8,1))
   OPEN(UNIT=iu_p, FILE=filename, STATUS='unknown', FORM='formatted')
   WRITE(iu_p,'("#t, v_t, b0_t, b0_s, cv_t, cp_t, alpha_t, gamma_t")')
   DO itemp=3,ntemp-2
      i_tar=0
      DO i=nvol-2, 3, -1
         IF (p(i,itemp)<press_plot(ipressp)) i_tar=i
      ENDDO
      i_tar=i_tar-1
      IF (i_tar==-1) THEN
         WRITE(stdout,'(5x,"At temperature ",i5,f15.6)') itemp, temp(itemp) 
         WRITE(stdout,'(5x,"press_plot ",f15.6)') press_plot(ipressp)/1.D8
         CALL errore('emp_f','press_plot is out of range',1)
      ENDIF
!
!   prepare the linear interpolation
!
      p1=p(i_tar,itemp)
      p2=p(i_tar+1,itemp)
      b1=bp(i_tar,itemp)
      b2=bp(i_tar+1,itemp)
      bs1=bs(i_tar,itemp)
      bs2=bs(i_tar+1,itemp)
      cv1=cv(i_tar,itemp)
      cv2=cv(i_tar+1,itemp)
      cp1=cp(i_tar,itemp)
      cp2=cp(i_tar+1,itemp)
      bet1=betat(i_tar,itemp)
      bet2=betat(i_tar+1,itemp)
      gam1=gammat(i_tar,itemp)
      gam2=gammat(i_tar+1,itemp)
      v1=v(i_tar)
      v2=v(i_tar+1)
!
!   interpolate
!
      v_o=v1+(press_plot(ipressp)-p1) * (v2 - v1) / (p2 - p1)
      deltav_o=(v_o-v1) / deltav
      b0_o=b1+( b2 - b1 )* deltav_o
      bs_o=bs1+( bs2 - bs1 )* deltav_o
      cv_o=cv1+( cv2 - cv1 )* deltav_o
      cp_o=cp1+( cp2 - cp1 )* deltav_o
      beta_o=bet1+( bet2 - bet1 )* deltav_o
      gamma_o=gam1+( gam2 - gam1 )* deltav_o
      WRITE(iu_p,'(8e20.10)') temp(itemp), v_o, b0_o, bs_o, cv_o, &
                              cp_o, beta_o, gamma_o
   ENDDO
   CLOSE(UNIT=iu_p, STATUS='keep')
ENDDO
!
!  For a set of temperatures, required in input,
!  interpolate at the volume that correspond to a given pressure 
!
DO itempp=1,ntemp_plot
   filename='const_t.dat.'//TRIM(float_to_char(temp_plot(itempp),1))
   OPEN(UNIT=iu_p, FILE=filename, STATUS='unknown', FORM='formatted')
   WRITE(iu_p,'("#p, v_p, b0_p, b0_sp, cv_p, cp_p, beta_p, gamma_p")')
   itemp=itemp_plot(itempp)
   DO ipress=1,npress
      i_tar=0
      DO i=nvol-2, 3, -1
         IF (p(i,itemp)<press(ipress)) i_tar=i
      ENDDO
      i_tar=i_tar-1
      IF (i_tar==-1) THEN
         WRITE(stdout,'(5x,"At temperature ",i5,f15.6)') itemp, temp(itemp)
         WRITE(stdout,'(5x,"press ",f15.6)') press(ipress)/1.D8
         CALL errore('emp_f','press is out of range',1)
      ENDIF
!
!   prepare the interpolation
!
      p1=p(i_tar,itemp)
      p2=p(i_tar+1,itemp)
      b1=bp(i_tar,itemp)
      b2=bp(i_tar+1,itemp)
      bs1=bs(i_tar,itemp)
      bs2=bs(i_tar+1,itemp)
      cv1=cv(i_tar,itemp)
      cv2=cv(i_tar+1,itemp)
      cp1=cp(i_tar,itemp)
      cp2=cp(i_tar+1,itemp)
      bet1=betat(i_tar,itemp)
      bet2=betat(i_tar+1,itemp)
      gam1=gammat(i_tar,itemp)
      gam2=gammat(i_tar+1,itemp)
      v1=v(i_tar)
      v2=v(i_tar+1)

      v_o=v1+(press(ipress)-p1) * (v2 - v1) / (p2 - p1)
      deltav_o=(v_o-v1) / deltav
      b0_o=b1+( b2 - b1 )* deltav_o
      bs_o=bs1+( bs2 - bs1 )* deltav_o
      cv_o=cv1+( cv2 - cv1 )* deltav_o
      cp_o=cp1+( cp2 - cp1 )* deltav_o
      beta_o=bet1+( bet2 - bet1 )* deltav_o
      gamma_o=gam1+( gam2 - gam1 )* deltav_o
      WRITE(iu_p,'(8e20.10)') press(ipress)/1.D8, v_o, b0_o, bs_o, cv_o, &
                              cp_o, beta_o, gamma_o
   ENDDO
   CLOSE(UNIT=iu_p, STATUS='keep')
ENDDO

WRITE(stdout,'(/,5x,"Termodynamic quantities have been calculated.")')
WRITE(stdout,'(5x,"What do you want to write?")')
WRITE(stdout,'(7x,"1) Write at T and P.")')
WRITE(stdout,'(7x,"2) Write at T and V.")')
WRITE(stdout,'(7x,"3) Nothing.")')
READ(5,*) iwrite
!
!   Read the target temperature, pressure or volume
!
IF (iwrite==1) THEN
   WRITE(stdout,'(/,5x,"Enter T (K) and P (kbar):")')
   READ(5,*) ttarget, ptarget

   WRITE(stdout,'(/,5x,"Thermodynamic quantities at T=",f10.2," K,",   &
                             " p=",f10.2," kbar:")') ttarget, ptarget
!
!  Pressure from kbar -> Pa
!
   ptarget=ptarget*1.D8

ELSEIF (iwrite==2) THEN
   WRITE(stdout,'(/,5x,"Enter T (K) and Volume (V0 units):")')
   READ(5,*) ttarget, vtarget

   WRITE(stdout,'(/,5x,"Thermodynamic quantities at T=",f10.2," K,",   &
                                   " v/v0=",f10.2)') ttarget, vtarget
   vtarget=vtarget*v0
!
!  Find the volume closest to vtarget, but smaller
!
   i_tar=0
   DO i=3,nvol-2
      IF (v(i)<vtarget) i_tar=i
   ENDDO
   IF (i_tar==0.OR. i_tar==(nvol-2)) &
      CALL errore('emp_f','Target volume is out of range',1)
ENDIF

IF (iwrite==1.OR.iwrite==2) THEN
!
!  find the temperature closest to ttarget, but lower
!
   itemp_tar=ntemp
   DO itemp=1, ntemp
      IF (temp(itemp)<ttarget) itemp_tar=itemp
   ENDDO
   IF (itemp_tar==ntemp) CALL errore('emp_f','Target temperature is        &
                                                            &out of range',1)
!
!   Interpolate the pressure at the temperature itemp_tar, for all volumes
!
   deltatar=(ttarget-temp(itemp_tar))/deltat

   IF (iwrite==1) THEN
      ALLOCATE(p_v(nvol))
      DO ivol=3,nvol-2
         p_v(ivol)= p(ivol, itemp_tar) + deltatar *          &
                     (p(ivol,itemp_tar+1)- p(ivol,itemp_tar))  
      ENDDO
!
!  Find volume for which the pressure is slightly larger 
!  that ptarget. Note that pressure decreases with volume
!
      i_tar=0
      DO i=3,nvol-2
         IF (p_v(i)>ptarget) i_tar=i
      ENDDO
      IF (i_tar==0.OR. i_tar==(nvol-2)) &
         CALL errore('emp_f','Target pressure is out of range',1)
!
!  The volume at ptarget (iwrite=1, T and P are given)
!
      vtarget=v(i_tar)+(ptarget-p_v(i_tar)) * deltav / &
                                            (p_v(i_tar+1) - p_v(i_tar))
      DEALLOCATE(p_v)
   ENDIF
!
!   Interpolate linearly all quantities at ttarget, for the two volumes
!   of the mesh that contain vtarget
!
   DO ivol=i_tar,i_tar+1
      p_tar(ivol-i_tar+1)= p(ivol, itemp_tar) + deltatar*          &
                     (p(ivol,itemp_tar+1)- p(ivol,itemp_tar))  
      bp_tar(ivol-i_tar+1)= bp(ivol, itemp_tar) + deltatar*        &
                     (bp(ivol,itemp_tar+1)- bp(ivol,itemp_tar)) 
      bs_tar(ivol-i_tar+1)= bs(ivol, itemp_tar) + deltatar*        &
                     (bs(ivol,itemp_tar+1)- bs(ivol,itemp_tar))  
      cv_tar(ivol-i_tar+1)= cv(ivol, itemp_tar) + deltatar*        & 
                     (cv(ivol,itemp_tar+1)- cv(ivol,itemp_tar))  
      cp_tar(ivol-i_tar+1)= cp(ivol, itemp_tar) + deltatar*        &
                     (cp(ivol,itemp_tar+1)-cp(ivol,itemp_tar))  
      beta_tar(ivol-i_tar+1)= betat(ivol, itemp_tar) + deltatar*  &
                  (betat(ivol,itemp_tar+1)- betat(ivol,itemp_tar))  
      gamma_tar(ivol-i_tar+1)= gammat(ivol, itemp_tar) + deltatar* &
                  (gammat(ivol,itemp_tar+1)- gammat(ivol,itemp_tar))  
      entr_tar(ivol-i_tar+1)= entr(ivol, itemp_tar) + deltatar* &
                  (entr(ivol,itemp_tar+1)- entr(ivol,itemp_tar))  
   ENDDO

   deltatar=ttarget-temp(itemp_tar)
!
!  The pressure at vtarget (iwrite=2, T and V are given)
!
   IF (iwrite==2) &
      ptarget=p_tar(1)+(vtarget-v(i_tar))*(p_tar(2)-p_tar(1))/deltav

   deltav_o=(vtarget-v(i_tar)) / deltav
!
!  the other quantities at the target volume
!
   b0_o=bp_tar(1)+(bp_tar(2) - bp_tar(1))* deltav_o
   bs_o=bs_tar(1)+(bs_tar(2) - bs_tar(1))* deltav_o
   cv_o=cv_tar(1)+(cv_tar(2) - cv_tar(1))* deltav_o
   cp_o=cp_tar(1)+(cp_tar(2) - cp_tar(1))* deltav_o
   beta_o=beta_tar(1)+(beta_tar(2) - beta_tar(1))*  deltav_o
   gamma_o=gamma_tar(1)+(gamma_tar(2) - gamma_tar(1))* deltav_o
   entr_o=entr_tar(1)+(entr_tar(2) - entr_tar(1))* deltav_o
!
!   Write on output the interpolated quantities
!
   WRITE(stdout,'(/,5x,"Volume=",f14.3," cm^3/mol",f16.3," A^3",&
      &f14.3," (a.u.)^3")') vtarget*1.D6, vtarget/avogadro/1.D-30, &
                           vtarget/avogadro/1.D-30/bohr_radius_angs**3
   IF (iwrite==1) &
      WRITE(stdout,'(/,5x,"V/V0=",27x,f14.3)') vtarget/v0
   IF (iwrite==2) &
      WRITE(stdout,'(/,5x,"Pressure= ",15x,f21.3," kbar")') ptarget/1.D8

   WRITE(stdout,'(5x,"Isothermal bulk modulus= ",f21.3," kbar")') b0_o/1.D8
   WRITE(stdout,'(5x,"Isoentropic bulk modulus= ",f20.3," kbar")') bs_o/1.D8
   WRITE(stdout,'(5x,"Isochoric heat capacity= ",f21.3," J/K/mol")') cv_o
   WRITE(stdout,'(5x,"Isobaric heat capacity= ",f22.3," J/K/mol")') cp_o
   WRITE(stdout,'(5x,"Thermal expansion x 10^6= ",f20.3, " 1/K")') beta_o*1.D6
   WRITE(stdout,'(5x,"Average Gruneisen parameter= ",f17.3)') gamma_o
   WRITE(stdout,'(5x,"Entropy = ",16x,f20.3," J/K/mol")') entr_o

ENDIF
!
!  Deallocate all variables
!
DEALLOCATE(press_plot)
DEALLOCATE(press)

DEALLOCATE(abt)
DEALLOCATE(gammat) 
DEALLOCATE(betat)
DEALLOCATE(cp) 
DEALLOCATE(cv)
DEALLOCATE(bs)
DEALLOCATE(bp)
DEALLOCATE(p)
DEALLOCATE(ener) 
DEALLOCATE(entr) 
DEALLOCATE(f)

DEALLOCATE(ivol_plot)
DEALLOCATE(vol_plot)

DEALLOCATE(auxv)
DEALLOCATE(gamma0)
DEALLOCATE(k01)
DEALLOCATE(k0)
DEALLOCATE(p0)
DEALLOCATE(e)
DEALLOCATE(v)

DEALLOCATE(itemp_plot)
DEALLOCATE(temp_plot)

DEALLOCATE(aux2)
DEALLOCATE(aux1)
DEALLOCATE(aux)
DEALLOCATE(b0_t) 
DEALLOCATE(v_t)
DEALLOCATE(temp)

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM emp_f
