!
! Copyright (C) 2022 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM emp_g
!
!   This small program prints the fit of the experimental 
!   properties of W from 300 K to 3700 K.
!   P. Gustavson, Int. Jour. of Thermophysics 6, 395 (1985).
!
USE kinds,          ONLY : DP
USE mp_global,      ONLY : mp_startup, mp_global_end
USE environment,    ONLY : environment_start, environment_end
USE thermodynamics_mod, ONLY : entropy_from_f, hc_from_entropy, beta_from_v

IMPLICIT NONE

CHARACTER(LEN=9) :: code='EMP_G'
CHARACTER(LEN=8) :: float_to_char
CHARACTER(LEN=256) :: filename
REAL(DP) :: tmin, deltat, t
REAL(DP) :: pmin, deltap, p
REAL(DP) :: am, bm, cm, dm, em, fm, gm, a, alpha0, alpha1, k0, k1, k2, n
REAL(DP) :: aux, aux1
REAL(DP), ALLOCATABLE :: g(:,:), press(:), v(:,:), b(:,:), s(:,:), cp(:,:), &
                         beta(:,:), temp(:)
INTEGER  :: i, ip, npt, npp, iu_out

CALL mp_startup ( start_images=.TRUE. )
CALL environment_start ( code )

deltat=5.0_DP
tmin=280._DP
npt=680

pmin=0._DP * 1.D8      ! 1.D8 conversione kbar -> Pa. Pressure in Pa.
deltap=500.0_DP * 1.D8
npp=7

am=-7647.26_DP
bm=130.4_DP
cm=-24.1_DP
dm=-1.936D-3
em=2.07D-7
fm=-5.33D-11
gm=4.45D4
a=9.5168D-6
alpha0=9.386D-6
alpha1=5.51D-9
k0=3.1575D-12
k1=1.6D-16
k2=3.1D-20
n=4.0_DP

ALLOCATE(press(npp))
ALLOCATE(temp(npt))
ALLOCATE(g(npt,npp))
ALLOCATE(v(npt,npp))
ALLOCATE(b(npt,npp))
ALLOCATE(beta(npt,npp))
ALLOCATE(cp(npt,npp))
ALLOCATE(s(npt,npp))

DO ip=1,npp
   press(ip)= pmin + deltap * (ip-1.0_DP)
ENDDO

DO i=1,npt
   temp(i)= tmin + deltat * (i-1.0_DP)
ENDDO

DO i=1,npt
   t=temp(i)
   DO ip=1, npp
      p=press(ip)
      g(i,ip)= am+ bm * t + cm * t * log(t) + dm * t**2 + &
               em * t**3 + fm * t**4 + gm / t
      aux=exp(alpha0 * t + 0.5_DP * alpha1 * t**2)
      aux1=k0 + k1*t + k2 * t**2
      g(i,ip)= g(i,ip) + a * aux/ aux1 / (n-1.0_DP) * &
                            ((1.0_DP+n*p*aux1)**(1.0_DP-1.0_DP/n)-1.0_DP)

      v(i,ip)= a * aux/ (1.0_DP+n*p*aux1)**(1.0_DP/n)
      b(i,ip)= (1.0_DP + n*p*aux1) / aux1
   ENDDO
ENDDO

!WRITE(6,'("#           p (kbar)      G (J/mol)         V/V0          B (kbar)")') 
!WRITE(6,'("#           T (K)         G (J/mol)         H-H0(J/mol)     S &
!                     &(J/mol/K)       C_p (J/mol/K) ")') 

DO ip=1,npp
   CALL entropy_from_f(g(1,ip),temp,npt,s(1,ip))
   CALL hc_from_entropy(cp(1,ip),temp,npt,s(1,ip))
   CALL beta_from_v(v(1,ip),temp,npt,beta(1,ip))
ENDDO

iu_out=25
DO ip=1,npp
   filename='data_p'//float_to_char(press(ip)/1.D8,0)
   OPEN(UNIT=iu_out, FILE=TRIM(filename), STATUS='unknown', FORM='formatted')
   WRITE(iu_out,'("# T (K)",8x,"G (J/mol)",8x,"H (J/mol)",8x,&
                "S (J/mol/K)",8x,"V (m^3)",8x,"B (kbar)",8x,&
                "beta*10^6",8x,"Cp (J/mol/K)" )') 
   DO i=3,npt-2
      WRITE(iu_out,'(8e19.10)') temp(i), g(i,ip), g(i,ip)+temp(i)*s(i,ip),  &
                           s(i,ip), v(i,ip), b(i,ip)/1.D8, beta(i,ip)*1.D6, &
                           cp(i,ip)
   ENDDO
   CLOSE(UNIT=iu_out, STATUS='KEEP')
ENDDO

DEALLOCATE(beta)
DEALLOCATE(cp)
DEALLOCATE(s)
DEALLOCATE(b)
DEALLOCATE(v)
DEALLOCATE(g)
DEALLOCATE(press)
DEALLOCATE(temp)

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM emp_g
