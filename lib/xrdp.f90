!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE xrdp_module
!
!  This module provides routines to compute the X-ray powder diffraction 
!  pattern. It has routines to:
!  
!  compute_form_factor : receive s and the element name and gives the form
!                        factor
!  compute_param_ff : receive the element name and sets the parameters
!  
!  compute_xrdp : recieves the direct and reciprocal lattice vectors of
!                 the crystal, the atomic positions, and the X-rays 
!                 wavelength and computes the scattering angles of the
!                 Bragg peaks and their relative intensities.
!
!  Useful discussions with F. Zadra during the development of this module
!  are gratefully acknowledged.
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
USE mp_images, ONLY : intra_image_comm
USE mp, ONLY : mp_sum
IMPLICIT NONE
SAVE
PRIVATE

PUBLIC compute_form_factor, write_form_factor, compute_xrdp, select_lambda

CONTAINS

SUBROUTINE write_form_factor(element,smin,smax,npoint,lcm,filename)

USE kinds, ONLY : DP
USE io_global, ONLY : ionode
IMPLICIT NONE
CHARACTER(LEN=2), INTENT(IN) :: element
REAL(DP), INTENT(INOUT) :: smin, smax
LOGICAL :: lcm
INTEGER, INTENT(IN) :: npoint
CHARACTER(LEN=*) :: filename

REAL(DP) :: deltas
REAL(DP), ALLOCATABLE :: ff(:), s(:)
INTEGER :: ipoint

IF (smin < 0.0_DP) smin=0.0_DP
IF (smax == 0.0_DP) smax=1.0_DP
deltas=(smax-smin) / (npoint -1)

ALLOCATE(s(npoint))
ALLOCATE(ff(npoint))

DO ipoint=1,npoint
   s(ipoint) = deltas * (ipoint-1) 
   IF (lcm) THEN
      CALL cromermann_form_factor(element, s(ipoint), ff(ipoint))
   ELSE
      CALL compute_form_factor(element, s(ipoint), ff(ipoint))
   ENDIF
END DO

IF (ionode) THEN
   OPEN (UNIT=47, FILE=TRIM(filename), STATUS='unknown',FORM='formatted')
   DO ipoint=1,npoint
      WRITE(47,'(2f16.8)') s(ipoint), ff(ipoint)
   ENDDO
   CLOSE(UNIT=47)
ENDIF

DEALLOCATE(s)
DEALLOCATE(ff)

RETURN
END SUBROUTINE write_form_factor

SUBROUTINE compute_xrdp(at,bg,alat,nat,tau,ntyp,ityp,atm,&
                                           lambda,ibrav,lcm,in_filename)
!
! This routine generates all the reciprocal lattice vectors contained in
! a sphere of radius 
! |G| = 1/ (4 pi \lambda)
! The value of lambda can be given also as input to this routine. 
! If \lambda is zero, the routine uses the value
! \lambda=1.541838 A which is the Cu K_alpha X ray wavelenght. 
! alat is supposed to be in a.u.
!
! This is faster than collecting all the G vectors from the different
! processors and selecting only the smaller ones.
! The routine opens a file called filename and writes there its output.
! The format of the file is the following
!
! miller_indeces (ijk) multiplicity  2\theta  d_ijk  intensity
!
! The intensity is calculated multiplying by the polarization factor, 
! the Lorentz factor, and the factor that account for the
! Debye-Scherrer geometry
! (1+cos(2 theta)^2)/(sin(theta)^2 cos(theta))
!
! If lcm is .TRUE. use the Cromer-Mann coefficients for the atomic
! form factors, otherwise uses the ... coefficients.
!
USE kinds, ONLY : DP
USE constants, ONLY : pi, bohr_radius_si
USE io_global, ONLY : ionode
USE lattices,  ONLY : compute_conventional

IMPLICIT NONE
INTEGER, INTENT(IN) :: nat, ntyp, ibrav
INTEGER, INTENT(IN) :: ityp(nat)
REAL(DP), INTENT(IN) :: tau(3,nat)
REAL(DP), INTENT(IN)   :: at(3,3), bg(3,3), alat, lambda
LOGICAL, INTENT(IN) :: lcm
CHARACTER(LEN=3), INTENT(IN) :: atm(ntyp)
CHARACTER(LEN=*), INTENT(IN) :: in_filename

INTEGER                :: nx1, nx2, nx3

INTEGER :: i, j, k, ngm, ntheta, na, nt, ipol
INTEGER, ALLOCATABLE :: igm(:), multiplicity(:)
REAL(DP) :: gvec(3), gmod, gcutm, gcutm2, lambda0, gtau, ff, intmax, atc(3,3)
REAL(DP), ALLOCATABLE :: g(:,:), gg(:), theta(:), s(:), gi(:,:), intensity(:), &
                         intensity_g(:), sgg(:), thetag(:)
INTEGER, ALLOCATABLE :: miller_g(:,:), miller(:,:)
COMPLEX(DP) :: factor
COMPLEX(DP) :: phase
CHARACTER(LEN=6) :: int_to_char
!
!  convert lambda in atomic units
!
lambda0=lambda * 1.D-10/ bohr_radius_si
!
!  set the maximum value of G
!
!
!   gcutm is in units of 2 \pi / a
! 
gcutm=alat * 2.0_DP / lambda0
gcutm2=gcutm**2
!
!  Determine the size of a mesh that contains the sphere
!
nx1=int(gcutm*sqrt(at(1,1)**2+at(2,1)**2+at(3,1)**2))+1
nx2=int(gcutm*sqrt(at(1,2)**2+at(2,2)**2+at(3,2)**2))+1
nx3=int(gcutm*sqrt(at(1,3)**2+at(2,3)**2+at(3,3)**2))+1
!
! count the total number of g-vectors such that |g|^2 <= gcutm2
! call this number ngm
!
ngm=0
DO i=-nx1, nx1
   DO j=-nx2, nx2
      DO k=-nx3, nx3
         gvec(:) = i * bg(:,1) + j * bg(:,2) + k * bg(:,3)
         gmod=gvec(1)**2 + gvec(2)**2 + gvec(3)**2 
         IF (gmod <= gcutm2) ngm=ngm+1
      ENDDO
   ENDDO
ENDDO 
!
! allocate memory for g-vectors
!
ALLOCATE(g(3,ngm))
ALLOCATE(gg(ngm))
ALLOCATE(igm(ngm))
ALLOCATE(miller_g(3,ngm))
ALLOCATE(intensity_g(ngm))
ALLOCATE(sgg(ngm))
ALLOCATE(thetag(ngm))
!
! re-calculate the g-vectors and their squared norms, store and order them
!
ngm=0
DO i=-nx1, nx1
   DO j=-nx2, nx2
      DO k=-nx3, nx3
         gvec(:) = i * bg(:,1) + j * bg(:,2) + k * bg(:,3)
         gmod=gvec(1)**2 + gvec(2)**2 + gvec(3)**2 
         IF (gmod <= gcutm2) THEN
            ngm=ngm+1
            g(:,ngm)=gvec(:)
            gg(ngm)=gmod
            igm(ngm) = ngm
            miller_g(1,ngm)=i
            miller_g(2,ngm)=j
            miller_g(3,ngm)=k
         ENDIF
      ENDDO
   ENDDO
ENDDO 
!
CALL hpsort_eps( ngm, gg, igm, 1.D-8 )

DO i=1,ngm
   sgg(i)=SQRT(gg(i)) / ( 2.0_DP * alat * bohr_radius_si * 1.D10 )
   thetag(i) = 2.0_DP*ASIN(SQRT(gg(i))*lambda0/2.0_DP/alat) * 180.0_DP / pi
   factor=(0.0_DP,0.0_DP)
   DO nt=1,ntyp
      IF (lcm) THEN
         CALL cromermann_form_factor(atm(nt), sgg(i), ff)
      ELSE
         CALL compute_form_factor(atm(nt), sgg(i), ff)
      ENDIF
      DO na=1,nat
         IF (ityp(na)==nt) THEN
            gtau = g(1,igm(i))*tau(1,na)+g(2,igm(i))*tau(2,na)+g(3,igm(i))*tau(3,na) 
            gtau = gtau * 2.0_DP * pi 
            phase = CMPLX(COS(gtau), SIN(gtau))
            factor= factor + ff * phase
         END IF 
      END DO
   END DO
   intensity_g(i)=ABS(factor)**2 *(1.0_DP + COS(thetag(i)*pi/180_DP)**2)* &
         0.5_DP / SIN(thetag(i)*pi/360_DP)**2 / COS(thetag(i)*pi/360_DP)
END DO
!
!  Now count how many different distances there are
!
ntheta=0
gmod=0.0_DP
DO i=2,ngm
   IF (ABS(gg(i)-gmod)> 1.D-8) THEN
      ntheta=ntheta+1
      gmod=gg(i)
   ENDIF
ENDDO

ALLOCATE(theta(ntheta))
ALLOCATE(s(ntheta))
ALLOCATE(gi(3,ntheta))
ALLOCATE(miller(3,ntheta))
ALLOCATE(intensity(ntheta))
ALLOCATE(multiplicity(ntheta))
!
!  And save the angle and s and add all intensities of the peaks with the same
!  modulus of g.
!
ntheta=0
gmod=0.0_DP
multiplicity=1
DO i=2,ngm
   IF (ABS(gg(i)-gmod)> 1.D-8) THEN
      ntheta=ntheta+1
      gmod = gg(i)
      theta(ntheta) = thetag(i)
      s(ntheta) = sgg(i)
      gi(:,ntheta) = g(:,igm(i))
      miller(:,ntheta)=miller_g(:,igm(i))
      intensity(ntheta)=intensity_g(i)
   ELSE
      intensity(ntheta)=intensity(ntheta)+intensity_g(i)
      multiplicity(ntheta)=multiplicity(ntheta)+1
   ENDIF
ENDDO
!
!  Find the peak with largest intensity
!
intmax = 0.0_DP
DO i=1,ntheta
   IF (intensity(i) > intmax) intmax=intensity(i)
ENDDO
!
!  Find the lattice vectors of the primitive lattice given the centered ones
!
CALL compute_conventional(at, atc, ibrav)
!
!  Change the Miller indeces for the centered lattices, so that they
!  are the Miller indeces of the primitive lattices
!
DO i=1, ntheta
   DO ipol=1,3
      miller(ipol,i) = gi(1,i) * atc(1,ipol)  &
                     + gi(2,i) * atc(2,ipol)  &
                     + gi(3,i) * atc(3,ipol)
   END DO
END DO
!
!  Open a file and write inside all the computed quantities
!
IF (ionode) THEN
   OPEN (UNIT=47, FILE=TRIM(in_filename), STATUS='unknown',FORM='formatted')

   WRITE(47,'("# theta calculated for lambda=",f15.7," A")') lambda0 &
                                                   * 1.D10 * bohr_radius_si
   WRITE(47,'("# Miller indeces multiplicity  2 x theta        d (A)     intensity" )')

   DO i=1,ntheta
      WRITE(47,'(3i5,5x,i4,f16.4,f18.8,f12.4)') miller(1:3,i), &
         multiplicity(i), theta(i), 1.0_DP/2.0_DP/s(i), intensity(i)*&
                                    100.0_DP / intmax
   ENDDO
   CLOSE(UNIT=47)
END IF

DEALLOCATE( g )
DEALLOCATE( gg )
DEALLOCATE( igm )
DEALLOCATE( intensity_g )
DEALLOCATE( thetag )
DEALLOCATE( sgg )

DEALLOCATE( theta )
DEALLOCATE( s )
DEALLOCATE( gi )
DEALLOCATE( miller )
DEALLOCATE( miller_g )
DEALLOCATE( intensity )
DEALLOCATE( multiplicity )

RETURN
END SUBROUTINE compute_xrdp

SUBROUTINE select_lambda(elem_name,lambda)
!
!  This routine sets the X-ray wavelengh from the anode element. The output values
!  are taken from J.A. Bearden "X-Ray Wavelengths". Rev. of Mod. Phys. 39, 78–124, (1967) and are in Angstrom.
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=2), INTENT(IN) :: elem_name
REAL(DP), INTENT(OUT) :: lambda

SELECT CASE (elem_name)
   CASE ('Cr')
      lambda=2.291_DP
   CASE ('Fe')
      lambda=1.93736_DP
   CASE ('Co')
      lambda=1.79026_DP
   CASE ('Cu')
      lambda=1.541838_DP
   CASE ('Mo')
      lambda=0.71073_DP
   CASE DEFAULT
!
!   If the name is not understood put the Cu line wavelength
!
      lambda=1.541838_DP
END SELECT

RETURN
END SUBROUTINE select_lambda



SUBROUTINE cromermann_form_factor(element, s, ff)
USE kinds, ONLY : DP

IMPLICIT NONE
CHARACTER(LEN=2), INTENT(IN) :: element
REAL(DP), INTENT(IN) :: s
REAL(DP), INTENT(OUT) :: ff

INTEGER :: ipar

REAL(DP) :: par_a_cm(4), par_b_cm(4), par_c_cm

CALL set_cromermann_coefficients(element, par_a_cm, par_b_cm, par_c_cm)

ff=0.0_DP
DO ipar=1,4
   ff= ff + par_a_cm(ipar) * EXP( -par_b_cm(ipar) * s**2 )
END DO

ff=ff+par_c_cm

RETURN
END SUBROUTINE cromermann_form_factor 

SUBROUTINE set_cromermann_coefficients(element,par_a_cm,par_b_cm,par_c_cm)
!
!  This routine set the parameters of Cromer and Mann that can be
!  found in the Intenational Tables of Crystallography Vol. C.
!
!  The parameters have been taken from the file f0_CromerMann.dat available
!  on the web. The file belongs to the DABAX library. More information on
!  DABAX can be found at:
!  http://www.esrf.fr/computing/expg/subgroups/theory/DABAX/dabax.html
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=2), INTENT(IN) :: element
REAL(DP),INTENT(INOUT) :: par_a_cm(4), par_b_cm(4), par_c_cm

par_a_cm=0.0_DP
par_b_cm=0.0_DP
par_c_cm=0.0_DP

SELECT CASE (element)
   CASE ('Ac')
      par_a_cm(1)=35.6597_DP   
      par_a_cm(2)=23.1032_DP   
      par_a_cm(3)=12.5977_DP   
      par_a_cm(4)=4.08655_DP   
      par_c_cm=13.5266_DP  
      par_b_cm(1)=0.589092_DP   
      par_b_cm(2)=3.65155_DP   
      par_b_cm(3)=18.599_DP  
      par_b_cm(4)=117.02_DP
   CASE ('Ag')
      par_a_cm(1)=19.2808_DP   
      par_a_cm(2)=16.6885_DP   
      par_a_cm(3)=4.8045_DP   
      par_a_cm(4)=1.0463_DP   
      par_c_cm=5.1790_DP  
      par_b_cm(1)=0.64460_DP   
      par_b_cm(2)=7.4726_DP
      par_b_cm(3)=24.6605_DP 
      par_b_cm(4)=99.8156_DP
   CASE ('Al')
      par_a_cm(1)=6.4202_DP   
      par_a_cm(2)=1.9002_DP   
      par_a_cm(3)=1.5936_DP   
      par_a_cm(4)=1.9646_DP   
      par_c_cm=1.1151_DP  
      par_b_cm(1)=3.0387_DP   
      par_b_cm(2)=0.7426_DP
      par_b_cm(3)=31.5472_DP 
      par_b_cm(4)=85.0886_DP
   CASE ('Am')
      par_a_cm(1)=36.6706_DP   
      par_a_cm(2)=24.0992_DP   
      par_a_cm(3)=17.3415_DP   
      par_a_cm(4)=3.49331_DP   
      par_c_cm=13.3592_DP  
      par_b_cm(1)=0.483629_DP   
      par_b_cm(2)=3.20647_DP   
      par_b_cm(3)=14.3136_DP  
      par_b_cm(4)=102.273_DP
   CASE ('Ar')
      par_a_cm(1)=7.4845_DP   
      par_a_cm(2)=6.7723_DP   
      par_a_cm(3)=0.6539_DP   
      par_a_cm(4)=1.6442_DP   
      par_c_cm=1.4445_DP  
      par_b_cm(1)=0.9072_DP   
      par_b_cm(2)=14.8407_DP
      par_b_cm(3)=43.8983_DP 
      par_b_cm(4)=33.3929_DP
   CASE ('As')
      par_a_cm(1)=16.6723_DP   
      par_a_cm(2)=6.0701_DP   
      par_a_cm(3)=3.4313_DP   
      par_a_cm(4)=4.2779_DP   
      par_c_cm=2.531_DP  
      par_b_cm(1)=2.6345_DP   
      par_b_cm(2)=0.2647_DP
      par_b_cm(3)=12.9479_DP  
      par_b_cm(4)=47.7972_DP
   CASE ('At')
      par_a_cm(1)=35.3163_DP   
      par_a_cm(2)=19.0211_DP   
      par_a_cm(3)=9.49887_DP   
      par_a_cm(4)=7.42518_DP   
      par_c_cm=13.7108_DP  
      par_b_cm(1)=0.68587_DP   
      par_b_cm(2)=3.97458_DP   
      par_b_cm(3)=11.3824_DP  
      par_b_cm(4)=45.47155_DP
   CASE ('Au')
      par_a_cm(1)=16.8819_DP   
      par_a_cm(2)=18.5913_DP   
      par_a_cm(3)=25.5582_DP   
      par_a_cm(4)=5.860_DP   
      par_c_cm=12.0658_DP  
      par_b_cm(1)=0.4611_DP   
      par_b_cm(2)=8.6216_DP   
      par_b_cm(3)=1.4826_DP  
      par_b_cm(4)=36.3956_DP
   CASE ('B')
      par_a_cm(1)=2.0545_DP   
      par_a_cm(2)=1.3326_DP   
      par_a_cm(3)=1.0979_DP   
      par_a_cm(4)=0.7068_DP   
      par_c_cm=-0.1932_DP  
      par_b_cm(1)=23.2185_DP   
      par_b_cm(2)=1.021_DP
      par_b_cm(3)=60.3498_DP  
      par_b_cm(4)=0.1403_DP
   CASE ('Ba')
      par_a_cm(1)=20.3361_DP   
      par_a_cm(2)=19.297_DP   
      par_a_cm(3)=10.888_DP   
      par_a_cm(4)=2.6959_DP   
      par_c_cm=2.7731_DP  
      par_b_cm(1)=3.2160_DP   
      par_b_cm(2)=0.2756_DP
      par_b_cm(3)=20.2073_DP 
      par_b_cm(4)=167.202_DP
   CASE ('Be')
      par_a_cm(1)=1.5919_DP  
      par_a_cm(2)=1.1278_DP  
      par_a_cm(3)=0.5391_DP  
      par_a_cm(4)=0.7029_DP  
      par_c_cm=3.8500D-02  
      par_b_cm(1)=43.6427_DP  
      par_b_cm(2)=1.8623_DP
      par_b_cm(3)=103.483_DP 
      par_b_cm(4)=0.542_DP
   CASE ('Bi')
      par_a_cm(1)=33.3689_DP   
      par_a_cm(2)=12.9510_DP   
      par_a_cm(3)=16.5877_DP   
      par_a_cm(4)=6.46920_DP   
      par_c_cm=13.5782_DP  
      par_b_cm(1)=0.7040_DP   
      par_b_cm(2)=2.9238_DP   
      par_b_cm(3)=8.7937_DP  
      par_b_cm(4)=48.0093_DP
   CASE ('Bk')
      par_a_cm(1)=36.7881_DP   
      par_a_cm(2)=24.7736_DP   
      par_a_cm(3)=17.8919_DP   
      par_a_cm(4)=4.23284_DP   
      par_c_cm=13.2754_DP  
      par_b_cm(1)=0.451018_DP   
      par_b_cm(2)=3.04619_DP   
      par_b_cm(3)=12.8946_DP  
      par_b_cm(4)=86.0030_DP
   CASE ('Br')
      par_a_cm(1)=17.1789_DP   
      par_a_cm(2)=5.23580_DP   
      par_a_cm(3)=5.6377_DP   
      par_a_cm(4)=3.9851_DP   
      par_c_cm=2.955700_DP  
      par_b_cm(1)=2.1723_DP   
      par_b_cm(2)=16.5796_DP
      par_b_cm(3)=0.2609_DP  
      par_b_cm(4)=41.4328_DP
   CASE ('C')
      par_a_cm(1)=2.31_DP   
      par_a_cm(2)=1.02_DP   
      par_a_cm(3)=1.5886_DP   
      par_a_cm(4)=0.865_DP   
      par_c_cm= 0.2156_DP  
      par_b_cm(1)=20.8439_DP   
      par_b_cm(2)=10.2075_DP
      par_b_cm(3)=0.5687_DP  
      par_b_cm(4)=51.6512_DP
   CASE ('Ca')
      par_a_cm(1)=8.6266_DP   
      par_a_cm(2)=7.3873_DP   
      par_a_cm(3)=1.5899_DP   
      par_a_cm(4)=1.0211_DP   
      par_c_cm=1.3751_DP  
      par_b_cm(1)=10.4421_DP   
      par_b_cm(2)=0.6599_DP
      par_b_cm(3)=85.7484_DP  
      par_b_cm(4)=178.437_DP
   CASE ('Cd')
      par_a_cm(1)=19.2214_DP   
      par_a_cm(2)=17.6444_DP   
      par_a_cm(3)=4.461_DP   
      par_a_cm(4)=1.6029_DP   
      par_c_cm=5.0694_DP  
      par_b_cm(1)=0.5946_DP   
      par_b_cm(2)=6.9089_DP
      par_b_cm(3)=24.7008_DP 
      par_b_cm(4)=87.4825_DP
   CASE ('Ce')
      par_a_cm(1)=21.1671_DP   
      par_a_cm(2)=19.7695_DP   
      par_a_cm(3)=11.8513_DP   
      par_a_cm(4)=3.33049_DP   
      par_c_cm=1.86264_DP  
      par_b_cm(1)=2.81219_DP   
      par_b_cm(2)=0.226836_DP
      par_b_cm(3)=17.6083_DP 
      par_b_cm(4)=127.113_DP
   CASE ('Cf')
      par_a_cm(1)=36.9185_DP   
      par_a_cm(2)=25.1995_DP  
      par_a_cm(3)=18.3317_DP   
      par_a_cm(4)=4.24391_DP   
      par_c_cm=13.2674_DP  
      par_b_cm(1)=0.437533_DP   
      par_b_cm(2)=3.00775_DP   
      par_b_cm(3)=12.4044_DP  
      par_b_cm(4)=83.7881_DP
   CASE ('Cl')
      par_a_cm(1)=11.4604_DP   
      par_a_cm(2)=7.1964_DP   
      par_a_cm(3)=6.2556_DP   
      par_a_cm(4)=1.6455_DP   
      par_c_cm=-9.5574_DP  
      par_b_cm(1)=1.040D-02   
      par_b_cm(2)=1.1662_DP
      par_b_cm(3)=18.5194_DP 
      par_b_cm(4)=47.7784_DP
   CASE ('Cm')
      par_a_cm(1)=36.6488_DP   
      par_a_cm(2)=24.4096_DP   
      par_a_cm(3)=17.399_DP   
      par_a_cm(4)=4.21665_DP   
      par_c_cm=13.2887_DP  
      par_b_cm(1)=0.465154_DP   
      par_b_cm(2)=3.08997_DP   
      par_b_cm(3)=13.4346_DP  
      par_b_cm(4)=88.4834_DP
   CASE ('Co')
      par_a_cm(1)=12.2841_DP   
      par_a_cm(2)=7.3409_DP   
      par_a_cm(3)=4.0034_DP   
      par_a_cm(4)=2.3488_DP   
      par_c_cm=1.0118_DP  
      par_b_cm(1)=4.2791_DP   
      par_b_cm(2)=0.2784_DP
      par_b_cm(3)=13.5359_DP  
      par_b_cm(4)=71.1692_DP
   CASE ('Cr')
      par_a_cm(1)=10.6406_DP   
      par_a_cm(2)=7.35370_DP   
      par_a_cm(3)=3.324000_DP   
      par_a_cm(4)=1.4922_DP   
      par_c_cm=1.1832_DP  
      par_b_cm(1)=6.1038_DP   
      par_b_cm(2)=0.392_DP
      par_b_cm(3)=20.2626_DP 
      par_b_cm(4)=98.7399_DP
   CASE ('Cs')
      par_a_cm(1)=20.3892_DP   
      par_a_cm(2)=19.1062_DP   
      par_a_cm(3)=10.6620_DP   
      par_a_cm(4)=1.4953_DP   
      par_c_cm=3.3352_DP  
      par_b_cm(1)=3.569_DP   
      par_b_cm(2)=0.3107_DP
      par_b_cm(3)=24.3879_DP 
      par_b_cm(4)=213.904_DP
   CASE ('Cu')
      par_a_cm(1)=13.338_DP   
      par_a_cm(2)=7.1676_DP   
      par_a_cm(3)=5.6158_DP   
      par_a_cm(4)=1.6735_DP   
      par_c_cm=1.191_DP  
      par_b_cm(1)=3.5828_DP   
      par_b_cm(2)=0.247_DP
      par_b_cm(3)=11.3966_DP  
      par_b_cm(4)=64.8126_DP
   CASE ('Dy')
      par_a_cm(1)=26.507_DP   
      par_a_cm(2)=17.6383_DP   
      par_a_cm(3)=14.5596_DP   
      par_a_cm(4)=2.96577_DP   
      par_c_cm=4.29728_DP  
      par_b_cm(1)=2.1802_DP   
      par_b_cm(2)=0.202172_DP
      par_b_cm(3)=12.1899_DP 
      par_b_cm(4)=111.874_DP
   CASE ('Eu')
      par_a_cm(1)=24.6274_DP   
      par_a_cm(2)=19.0886_DP   
      par_a_cm(3)=13.7603_DP   
      par_a_cm(4)=2.9227_DP   
      par_c_cm=2.5745_DP  
      par_b_cm(1)=2.3879_DP   
      par_b_cm(2)=0.1942_DP
      par_b_cm(3)=13.7546_DP 
      par_b_cm(4)=123.174_DP
   CASE ('Er')
      par_a_cm(1)=27.6563_DP   
      par_a_cm(2)=16.4285_DP   
      par_a_cm(3)=14.9779_DP   
      par_a_cm(4)=2.98233_DP   
      par_c_cm=5.92046_DP  
      par_b_cm(1)=2.07356_DP   
      par_b_cm(2)=0.223545_DP
      par_b_cm(3)=11.3604_DP 
      par_b_cm(4)=105.703_DP
   CASE ('F')
      par_a_cm(1)=3.5392_DP   
      par_a_cm(2)=2.6412_DP   
      par_a_cm(3)=1.517_DP   
      par_a_cm(4)=1.0243_DP  
      par_c_cm=0.2776_DP  
      par_b_cm(1)=10.2825_DP   
      par_b_cm(2)=4.2944_DP
      par_b_cm(3)=0.2615_DP  
      par_b_cm(4)=26.1476_DP
   CASE ('Fe')
      par_a_cm(1)=11.7695_DP   
      par_a_cm(2)=7.3573_DP   
      par_a_cm(3)=3.5222_DP   
      par_a_cm(4)=2.3045_DP   
      par_c_cm=1.0369_DP  
      par_b_cm(1)=4.7611_DP   
      par_b_cm(2)=0.3072_DP
      par_b_cm(3)=15.3535_DP 
      par_b_cm(4)=76.8805_DP
   CASE ('Fr')
      par_a_cm(1)=35.9299_DP  
      par_a_cm(2)=23.0547_DP  
      par_a_cm(3)=12.1439_DP  
      par_a_cm(4)=2.11253_DP  
      par_c_cm=13.7247_DP 
      par_b_cm(1)=0.646453_DP  
      par_b_cm(2)=4.17619_DP  
      par_b_cm(3)=23.1052_DP 
      par_b_cm(4)=150.645_DP
   CASE ('Ga')
      par_a_cm(1)=15.2354_DP   
      par_a_cm(2)=6.7006_DP   
      par_a_cm(3)=4.3591_DP   
      par_a_cm(4)=2.9623_DP   
      par_c_cm=1.7189_DP  
      par_b_cm(1)=3.0669_DP   
      par_b_cm(2)=0.2412_DP
      par_b_cm(3)=10.7805_DP  
      par_b_cm(4)=61.4135_DP
   CASE ('Gd')
      par_a_cm(1)=25.0709_DP   
      par_a_cm(2)=19.0798_DP   
      par_a_cm(3)=13.8518_DP   
      par_a_cm(4)=3.54545_DP   
      par_c_cm=2.4196_DP  
      par_b_cm(1)=2.25341_DP   
      par_b_cm(2)=0.181951_DP
      par_b_cm(3)=12.9331_DP 
      par_b_cm(4)=101.398_DP
   CASE ('Ge')
      par_a_cm(1)=16.0816_DP   
      par_a_cm(2)=6.3747_DP   
      par_a_cm(3)=3.7068_DP   
      par_a_cm(4)=3.683_DP   
      par_c_cm=2.1313_DP  
      par_b_cm(1)=2.8509_DP   
      par_b_cm(2)=0.2516_DP
      par_b_cm(3)=11.4468_DP  
      par_b_cm(4)=54.7625_DP
   CASE ('H')
      par_a_cm(1)=0.493002_DP   
      par_a_cm(2)=0.322912_DP   
      par_a_cm(3)=0.140191_DP   
      par_a_cm(4)=4.081D-02   
      par_c_cm=3.0380001D-03  
      par_b_cm(1)=10.5109_DP
      par_b_cm(2)=26.1257_DP  
      par_b_cm(3)=3.14236_DP 
      par_b_cm(4)=57.7997_DP
   CASE ('He')
      par_a_cm(1)=0.8734_DP   
      par_a_cm(2)=0.6309_DP   
      par_a_cm(3)=0.3112_DP   
      par_a_cm(4)=0.178_DP   
      par_c_cm=6.3999998D-03  
      par_b_cm(1)=9.1037_DP
      par_b_cm(2)=3.3568_DP   
      par_b_cm(3)=22.9276_DP  
      par_b_cm(4)=0.9821_DP
   CASE ('Hf')
      par_a_cm(1)=29.144_DP   
      par_a_cm(2)=15.1726_DP   
      par_a_cm(3)=14.7586_DP   
      par_a_cm(4)=4.30013_DP   
      par_c_cm=8.58154_DP  
      par_b_cm(1)=1.83262_DP   
      par_b_cm(2)=9.5999_DP
      par_b_cm(3)=0.275116_DP  
      par_b_cm(4)=72.029_DP
   CASE ('Hg')
      par_a_cm(1)=20.6809_DP   
      par_a_cm(2)=19.0417_DP   
      par_a_cm(3)=21.6575_DP   
      par_a_cm(4)=5.9676_DP   
      par_c_cm=12.6089_DP  
      par_b_cm(1)=0.545_DP   
      par_b_cm(2)=8.4484_DP   
      par_b_cm(3)=1.5729_DP  
      par_b_cm(4)=38.3246_DP
   CASE ('Ho')
      par_a_cm(1)=26.9049_DP   
      par_a_cm(2)=17.294_DP   
      par_a_cm(3)=14.5583_DP   
      par_a_cm(4)=3.63837_DP   
      par_c_cm=4.56796_DP  
      par_b_cm(1)=2.07051_DP   
      par_b_cm(2)=0.19794_DP
      par_b_cm(3)=11.4407_DP 
      par_b_cm(4)=92.6566_DP
   CASE ('In')
      par_a_cm(1)=19.1624_DP   
      par_a_cm(2)=18.5596_DP   
      par_a_cm(3)=4.2948_DP   
      par_a_cm(4)=2.0396_DP   
      par_c_cm=4.9391_DP  
      par_b_cm(1)=0.5476_DP   
      par_b_cm(2)=6.3776_DP
      par_b_cm(3)=25.8499_DP  
      par_b_cm(4)=92.8029_DP
   CASE ('Ir')
      par_a_cm(1)=27.3049_DP   
      par_a_cm(2)=16.7296_DP   
      par_a_cm(3)=15.6115_DP   
      par_a_cm(4)=5.83377_DP   
      par_c_cm=11.4722_DP  
      par_b_cm(1)=1.59279_DP   
      par_b_cm(2)=8.86553_DP
      par_b_cm(3)=0.417916_DP  
      par_b_cm(4)=45.0011_DP
   CASE ('La')
      par_a_cm(1)=20.578_DP   
      par_a_cm(2)=19.599_DP   
      par_a_cm(3)=11.3727_DP   
      par_a_cm(4)=3.28719_DP   
      par_c_cm=2.14678_DP  
      par_b_cm(1)=2.94817_DP   
      par_b_cm(2)=0.244475_DP
      par_b_cm(3)=18.7726_DP 
      par_b_cm(4)=133.124_DP
   CASE ('Lu')
      par_a_cm(1)=28.9476_DP   
      par_a_cm(2)=15.2208_DP   
      par_a_cm(3)=15.1_DP   
      par_a_cm(4)=3.71601_DP   
      par_c_cm=7.97628_DP  
      par_b_cm(1)=1.90182_DP   
      par_b_cm(2)=9.98519_DP
      par_b_cm(3)=0.261033_DP  
      par_b_cm(4)=84.3298_DP
   CASE ('K')
      par_a_cm(1)=8.2186_DP   
      par_a_cm(2)=7.4398_DP   
      par_a_cm(3)=1.0519_DP   
      par_a_cm(4)=0.8659_DP   
      par_c_cm=1.4228_DP  
      par_b_cm(1)=12.7949_DP   
      par_b_cm(2)=0.7748_DP
      par_b_cm(3)=213.187_DP 
      par_b_cm(4)=41.6841_DP
   CASE ('Kr')
      par_a_cm(1)=17.3555_DP   
      par_a_cm(2)=6.7286_DP   
      par_a_cm(3)=5.5493_DP   
      par_a_cm(4)=3.5375_DP   
      par_c_cm=2.825_DP  
      par_b_cm(1)=1.9384_DP   
      par_b_cm(2)=16.5623_DP
      par_b_cm(3)=0.2261_DP  
      par_b_cm(4)=39.3972_DP
   CASE ('I')
      par_a_cm(1)=20.1472_DP  
      par_a_cm(2)=18.9949_DP  
      par_a_cm(3)=7.5138_DP  
      par_a_cm(4)=2.2735_DP  
      par_c_cm=4.0712_DP 
      par_b_cm(1)=4.347_DP  
      par_b_cm(2)=0.3814_DP
      par_b_cm(3)=27.766_DP 
      par_b_cm(4)=66.8776_DP
   CASE ('Li')
      par_a_cm(1)=1.1282_DP   
      par_a_cm(2)=0.7508_DP   
      par_a_cm(3)=0.6175_DP   
      par_a_cm(4)=0.4653_DP   
      par_c_cm=3.7700001D-02  
      par_b_cm(1)=3.9546_DP
      par_b_cm(2)=1.0524_DP  
      par_b_cm(3)=85.3905_DP 
      par_b_cm(4)=168.261_DP
   CASE ('Mg')
      par_a_cm(1)=5.4204_DP   
      par_a_cm(2)=2.1735_DP   
      par_a_cm(3)=1.2269_DP   
      par_a_cm(4)=2.3073_DP   
      par_c_cm=0.8584_DP  
      par_b_cm(1)=2.8275_DP   
      par_b_cm(2)=79.2611_DP
      par_b_cm(3)=0.3808_DP  
      par_b_cm(4)=7.1937_DP
   CASE ('Mn')
      par_a_cm(1)=11.2819_DP   
      par_a_cm(2)=7.3573_DP   
      par_a_cm(3)=3.0193_DP   
      par_a_cm(4)=2.2441_DP   
      par_c_cm=1.0896_DP  
      par_b_cm(1)=5.3409_DP   
      par_b_cm(2)=0.3432_DP
      par_b_cm(3)=17.8674_DP 
      par_b_cm(4)=83.7543_DP
   CASE ('Mo')
      par_a_cm(1)=3.7025_DP   
      par_a_cm(2)=17.2356_DP   
      par_a_cm(3)=12.8876_DP   
      par_a_cm(4)=3.7429_DP   
      par_c_cm=4.3875_DP  
      par_b_cm(1)=0.2772_DP   
      par_b_cm(2)=1.0958_DP
      par_b_cm(3)=11.004_DP 
      par_b_cm(4)=61.6584_DP
   CASE ('N')
      par_a_cm(1)=12.2126_DP   
      par_a_cm(2)=3.1322_DP   
      par_a_cm(3)=2.0125_DP   
      par_a_cm(4)=1.1663_DP   
      par_c_cm=  -11.529_DP  
      par_b_cm(1)=5.7000001D-03   
      par_b_cm(2)=9.8933_DP
      par_b_cm(3)=28.9975_DP  
      par_b_cm(4)=0.5826_DP
   CASE ('Na')
      par_a_cm(1)=4.7626_DP   
      par_a_cm(2)=3.1736_DP   
      par_a_cm(3)=1.2674_DP   
      par_a_cm(4)=1.1128_DP   
      par_c_cm=0.676_DP  
      par_b_cm(1)=3.285_DP   
      par_b_cm(2)=8.8422_DP
      par_b_cm(3)=0.3136_DP  
      par_b_cm(4)=129.424_DP
   CASE ('Nb')
      par_a_cm(1)=17.6142_DP   
      par_a_cm(2)=12.0144_DP   
      par_a_cm(3)=4.04183_DP   
      par_a_cm(4)=3.53346_DP   
      par_c_cm=3.75591_DP  
      par_b_cm(1)=1.18865_DP   
      par_b_cm(2)=11.766_DP
      par_b_cm(3)=0.204785_DP  
      par_b_cm(4)=69.7957_DP
   CASE ('Nd')
      par_a_cm(1)=22.6845_DP   
      par_a_cm(2)=19.6847_DP   
      par_a_cm(3)=12.774_DP   
      par_a_cm(4)=2.85137_DP   
      par_c_cm=1.98486_DP  
      par_b_cm(1)=2.66248_DP   
      par_b_cm(2)=0.210628_DP
      par_b_cm(3)=15.885_DP 
      par_b_cm(4)=137.903_DP
   CASE ('Ne')
      par_a_cm(1)=3.9553_DP   
      par_a_cm(2)=3.1125_DP   
      par_a_cm(3)=1.4546_DP   
      par_a_cm(4)=1.1251_DP   
      par_c_cm=0.3515_DP  
      par_b_cm(1)=8.4042_DP   
      par_b_cm(2)=3.4262_DP
      par_b_cm(3)=0.2306_DP  
      par_b_cm(4)=21.7184_DP
   CASE ('Ni')
      par_a_cm(1)=12.8376_DP   
      par_a_cm(2)=7.292_DP   
      par_a_cm(3)=4.4438_DP   
      par_a_cm(4)=2.38_DP   
      par_c_cm=1.0341_DP  
      par_b_cm(1)=3.8785_DP   
      par_b_cm(2)=0.2565_DP
      par_b_cm(3)=12.1763_DP  
      par_b_cm(4)=66.3421_DP
   CASE ('Np')
      par_a_cm(1)=36.1874_DP   
      par_a_cm(2)=23.5964_DP   
      par_a_cm(3)=15.6402_DP   
      par_a_cm(4)=4.18550_DP   
      par_c_cm=13.3573_DP  
      par_b_cm(1)=0.511929_DP   
      par_b_cm(2)=3.25396_DP   
      par_b_cm(3)=15.3622_DP  
      par_b_cm(4)=97.4908_DP
   CASE ('O')
      par_a_cm(1)=3.0485_DP   
      par_a_cm(2)=2.2868_DP   
      par_a_cm(3)=1.5463_DP   
      par_a_cm(4)=0.867_DP   
      par_c_cm=0.2508_DP  
      par_b_cm(1)=13.2771_DP   
      par_b_cm(2)=5.7011_DP
      par_b_cm(3)=0.3239_DP  
      par_b_cm(4)=32.9089_DP
   CASE ('Os')
      par_a_cm(1)=28.1894_DP   
      par_a_cm(2)=16.155_DP   
      par_a_cm(3)=14.9305_DP   
      par_a_cm(4)=5.67589_DP   
      par_c_cm=11.0005_DP  
      par_b_cm(1)=1.62903_DP   
      par_b_cm(2)=8.97948_DP
      par_b_cm(3)=0.382661_DP  
      par_b_cm(4)=48.1647_DP
   CASE ('P')
      par_a_cm(1)=6.4345_DP   
      par_a_cm(2)=4.1791_DP   
      par_a_cm(3)=1.78_DP   
      par_a_cm(4)=1.4908_DP   
      par_c_cm=1.1149_DP  
      par_b_cm(1)=1.9067_DP   
      par_b_cm(2)=27.157_DP
      par_b_cm(3)=0.526_DP  
      par_b_cm(4)=68.1645_DP
   CASE ('Pa')
      par_a_cm(1)=35.8847_DP   
      par_a_cm(2)=23.2948_DP   
      par_a_cm(3)=14.1891_DP   
      par_a_cm(4)=4.17287_DP   
      par_c_cm=13.4287_DP  
      par_b_cm(1)=0.547751_DP   
      par_b_cm(2)=3.41519_DP   
      par_b_cm(3)=16.9235_DP  
      par_b_cm(4)=105.251_DP
   CASE ('Pb')
      par_a_cm(1)=31.0617_DP  
      par_a_cm(2)=13.0637_DP  
      par_a_cm(3)=18.442_DP  
      par_a_cm(4)=5.9696_DP  
      par_c_cm=13.4118_DP 
      par_b_cm(1)=0.6902_DP  
      par_b_cm(2)=2.3576_DP  
      par_b_cm(3)=8.618_DP 
      par_b_cm(4)=47.257_DP
   CASE ('Pd')
      par_a_cm(1)=19.3319_DP   
      par_a_cm(2)=15.5017_DP   
      par_a_cm(3)=5.29537_DP   
      par_a_cm(4)=0.605844_DP   
      par_c_cm=5.26593_DP  
      par_b_cm(1)=0.698655_DP   
      par_b_cm(2)=7.98929_DP
      par_b_cm(3)=25.2052_DP 
      par_b_cm(4)=76.8986_DP
   CASE ('Po')
      par_a_cm(1)=34.6726_DP   
      par_a_cm(2)=15.4733_DP   
      par_a_cm(3)=13.1138_DP   
      par_a_cm(4)=7.02588_DP   
      par_c_cm=13.677_DP  
      par_b_cm(1)=0.700999_DP   
      par_b_cm(2)=3.55078_DP   
      par_b_cm(3)=9.55642_DP  
      par_b_cm(4)=47.0045_DP
   CASE ('Pm')
      par_a_cm(1)=23.3405_DP   
      par_a_cm(2)=19.6095_DP   
      par_a_cm(3)=13.1235_DP   
      par_a_cm(4)=2.87516_DP   
      par_c_cm=2.02876_DP  
      par_b_cm(1)=2.5627_DP   
      par_b_cm(2)=0.202088_DP
      par_b_cm(3)=15.1009_DP 
      par_b_cm(4)=132.721_DP
   CASE ('Pr')
      par_a_cm(1)=22.044_DP   
      par_a_cm(2)=19.6697_DP   
      par_a_cm(3)=12.3856_DP   
      par_a_cm(4)=2.82428_DP   
      par_c_cm=2.0583_DP  
      par_b_cm(1)=2.77393_DP   
      par_b_cm(2)=0.222087_DP
      par_b_cm(3)=16.7669_DP 
      par_b_cm(4)=143.644_DP
   CASE ('Pt')
      par_a_cm(1)=27.0059_DP   
      par_a_cm(2)=17.7639_DP   
      par_a_cm(3)=15.7131_DP   
      par_a_cm(4)=5.7837_DP   
      par_c_cm=11.6883_DP  
      par_b_cm(1)=1.51293_DP   
      par_b_cm(2)=8.81174_DP   
      par_b_cm(3)=0.424593_DP  
      par_b_cm(4)=38.6103_DP
   CASE ('Pu')
      par_a_cm(1)=35.5103_DP   
      par_a_cm(2)=22.5787_DP   
      par_a_cm(3)=12.7766_DP   
      par_a_cm(4)=4.92159_DP   
      par_c_cm=13.2116_DP  
      par_b_cm(1)=0.498626_DP   
      par_b_cm(2)=2.96627_DP   
      par_b_cm(3)=11.9484_DP  
      par_b_cm(4)=22.7502_DP
   CASE ('Ra')
      par_a_cm(1)=35.763_DP   
      par_a_cm(2)=22.9064_DP   
      par_a_cm(3)=12.4739_DP   
      par_a_cm(4)=3.21097_DP   
      par_c_cm=13.6211_DP  
      par_b_cm(1)=0.616341_DP   
      par_b_cm(2)=3.87135_DP   
      par_b_cm(3)=19.9887_DP  
      par_b_cm(4)=142.325_DP
   CASE ('Rb')
      par_a_cm(1)=17.1784_DP   
      par_a_cm(2)=9.64350_DP   
      par_a_cm(3)=5.1399_DP   
      par_a_cm(4)=1.5292_DP   
      par_c_cm=3.4873_DP  
      par_b_cm(1)=1.7888_DP   
      par_b_cm(2)=17.3151_DP
      par_b_cm(3)=0.2748_DP  
      par_b_cm(4)=164.934_DP
   CASE ('Re')
      par_a_cm(1)=28.7621_DP   
      par_a_cm(2)=15.7189_DP   
      par_a_cm(3)=14.5564_DP   
      par_a_cm(4)=5.44174_DP   
      par_c_cm=10.472_DP  
      par_b_cm(1)=1.67191_DP   
      par_b_cm(2)=9.09227_DP
      par_b_cm(3)=0.3505_DP  
      par_b_cm(4)=52.0861_DP
   CASE ('Rh')
      par_a_cm(1)=19.2957_DP   
      par_a_cm(2)=14.3501_DP   
      par_a_cm(3)=4.73425_DP   
      par_a_cm(4)=1.28918_DP   
      par_c_cm=5.328_DP  
      par_b_cm(1)=0.751536_DP   
      par_b_cm(2)=8.21758_DP
      par_b_cm(3)=25.8749_DP 
      par_b_cm(4)=98.6062_DP
   CASE ('Ru')
      par_a_cm(1)=19.2674_DP   
      par_a_cm(2)=12.9182_DP   
      par_a_cm(3)=4.86337_DP   
      par_a_cm(4)=1.56756_DP   
      par_c_cm=5.37874_DP  
      par_b_cm(1)=0.80852_DP   
      par_b_cm(2)=8.43467_DP
      par_b_cm(3)=24.7997_DP 
      par_b_cm(4)=94.2928_DP
   CASE ('S')
      par_a_cm(1)=6.9053_DP   
      par_a_cm(2)=5.2034_DP   
      par_a_cm(3)=1.4379_DP   
      par_a_cm(4)=1.5863_DP   
      par_c_cm=0.8669_DP  
      par_b_cm(1)=1.4679_DP   
      par_b_cm(2)=22.2151_DP
      par_b_cm(3)=0.2536_DP  
      par_b_cm(4)=56.172_DP
   CASE ('Sb')
      par_a_cm(1)=19.6418_DP   
      par_a_cm(2)=19.0455_DP   
      par_a_cm(3)=5.0371_DP   
      par_a_cm(4)=2.6827_DP   
      par_c_cm=4.5909_DP  
      par_b_cm(1)=5.3034_DP   
      par_b_cm(2)=0.4607_DP
      par_b_cm(3)=27.9074_DP 
      par_b_cm(4)=75.2825_DP
   CASE ('Sc')
      par_a_cm(1)=9.189_DP   
      par_a_cm(2)=7.3679_DP   
      par_a_cm(3)=1.6409_DP   
      par_a_cm(4)=1.468_DP   
      par_c_cm=1.3329_DP  
      par_b_cm(1)=9.0213_DP   
      par_b_cm(2)=0.5729_DP
      par_b_cm(3)=136.108_DP  
      par_b_cm(4)=51.3531_DP
   CASE ('Se')
      par_a_cm(1)=17.0006_DP   
      par_a_cm(2)=5.8196_DP   
      par_a_cm(3)=3.9731_DP   
      par_a_cm(4)=4.3543_DP   
      par_c_cm=2.8409_DP  
      par_b_cm(1)=2.4098_DP   
      par_b_cm(2)=0.2726_DP
      par_b_cm(3)=15.2372_DP  
      par_b_cm(4)=43.8163_DP
   CASE ('Si')
      par_a_cm(1)=5.66269_DP   
      par_a_cm(2)=3.07164_DP   
      par_a_cm(3)=2.62446_DP   
      par_a_cm(4)=1.3932_DP   
      par_c_cm=1.24707_DP  
      par_b_cm(1)=2.6652_DP   
      par_b_cm(2)=38.6634_DP
      par_b_cm(3)=0.916946_DP
      par_b_cm(4)=93.5458_DP
   CASE ('Sm')
      par_a_cm(1)=24.0042_DP   
      par_a_cm(2)=19.4258_DP   
      par_a_cm(3)=13.4396_DP   
      par_a_cm(4)=2.89604_DP   
      par_c_cm=2.20963_DP  
      par_b_cm(1)=2.47274_DP   
      par_b_cm(2)=0.196451_DP
      par_b_cm(3)=14.3996_DP 
      par_b_cm(4)=128.007_DP
   CASE ('Sn')
      par_a_cm(1)=19.1889_DP   
      par_a_cm(2)=19.1005_DP   
      par_a_cm(3)=4.4585_DP   
      par_a_cm(4)=2.4663_DP   
      par_c_cm=4.7821_DP  
      par_b_cm(1)=5.8303_DP   
      par_b_cm(2)=0.5031_DP
      par_b_cm(3)=26.8909_DP 
      par_b_cm(4)=83.9571_DP
   CASE ('Sr')
      par_a_cm(1)=17.5663_DP   
      par_a_cm(2)=9.8184_DP   
      par_a_cm(3)=5.422_DP   
      par_a_cm(4)=2.6694_DP   
      par_c_cm=2.5064_DP  
      par_b_cm(1)=1.5564_DP   
      par_b_cm(2)=14.0988_DP
      par_b_cm(3)=0.1664_DP  
      par_b_cm(4)=132.376_DP
   CASE ('Ta')
      par_a_cm(1)=29.2024_DP  
      par_a_cm(2)=15.2293_DP  
      par_a_cm(3)=14.5135_DP  
      par_a_cm(4)=4.76492_DP  
      par_c_cm=9.24354_DP 
      par_b_cm(1)=1.77333_DP  
      par_b_cm(2)=9.37046_DP
      par_b_cm(3)=0.295977_DP 
      par_b_cm(4)=63.3644_DP
   CASE ('Tb')
      par_a_cm(1)=25.8976_DP   
      par_a_cm(2)=18.2185_DP   
      par_a_cm(3)=14.3167_DP   
      par_a_cm(4)=2.95354_DP   
      par_c_cm=3.58324_DP  
      par_b_cm(1)=2.24256_DP   
      par_b_cm(2)=0.196143_DP
      par_b_cm(3)=12.6648_DP 
      par_b_cm(4)=115.362_DP
   CASE ('Tc')
      par_a_cm(1)=19.1301_DP   
      par_a_cm(2)=11.0948_DP   
      par_a_cm(3)=4.64901_DP   
      par_a_cm(4)=2.71263_DP   
      par_c_cm=5.40428_DP  
      par_b_cm(1)=0.864132_DP   
      par_b_cm(2)=8.14487_DP
      par_b_cm(3)=21.5707_DP 
      par_b_cm(4)=86.8472_DP
   CASE ('Te')
      par_a_cm(1)=19.9644_DP   
      par_a_cm(2)=19.0138_DP   
      par_a_cm(3)=6.14487_DP   
      par_a_cm(4)=2.5239_DP   
      par_c_cm=4.352_DP  
      par_b_cm(1)=4.81742_DP   
      par_b_cm(2)=0.420885_DP
      par_b_cm(3)=28.5284_DP 
      par_b_cm(4)=70.8403_DP
   CASE ('Ti')
      par_a_cm(1)=9.7595_DP   
      par_a_cm(2)=7.3558_DP   
      par_a_cm(3)=1.6991_DP   
      par_a_cm(4)=1.9021_DP   
      par_c_cm=1.2807_DP  
      par_b_cm(1)=7.8508_DP   
      par_b_cm(2)=0.5_DP
      par_b_cm(3)=35.6338_DP 
      par_b_cm(4)=116.105_DP
   CASE ('Th')
      par_a_cm(1)=35.5645_DP   
      par_a_cm(2)=23.4219_DP   
      par_a_cm(3)=12.7473_DP   
      par_a_cm(4)=4.80703_DP   
      par_c_cm=13.4314_DP  
      par_b_cm(1)=0.563359_DP   
      par_b_cm(2)=3.46204_DP   
      par_b_cm(3)=17.8309_DP  
      par_b_cm(4)=99.1722_DP
   CASE ('Tl')
      par_a_cm(1)=27.5446_DP   
      par_a_cm(2)=19.1584_DP   
      par_a_cm(3)=15.538_DP   
      par_a_cm(4)=5.52593_DP   
      par_c_cm=13.1746_DP  
      par_b_cm(1)=0.65515_DP   
      par_b_cm(2)=8.70751_DP   
      par_b_cm(3)=1.96347_DP  
      par_b_cm(4)=45.8149_DP
   CASE ('Tm')
      par_a_cm(1)=28.1819_DP   
      par_a_cm(2)=15.8851_DP   
      par_a_cm(3)=15.1542_DP   
      par_a_cm(4)=2.98706_DP   
      par_c_cm=6.75621_DP  
      par_b_cm(1)=2.02859_DP   
      par_b_cm(2)=0.238849_DP
      par_b_cm(3)=10.9975_DP 
      par_b_cm(4)=102.961_DP
   CASE ('U')
      par_a_cm(1)=36.0228_DP   
      par_a_cm(2)=23.4128_DP   
      par_a_cm(3)=14.9491_DP   
      par_a_cm(4)=4.188_DP   
      par_c_cm=13.3966_DP  
      par_b_cm(1)=0.5293_DP   
      par_b_cm(2)=3.3253_DP   
      par_b_cm(3)=16.0927_DP  
      par_b_cm(4)=100.613_DP
   CASE ('V')
      par_a_cm(1)=10.2971_DP   
      par_a_cm(2)=7.3511_DP   
      par_a_cm(3)=2.0703_DP   
      par_a_cm(4)=2.0571_DP   
      par_c_cm=1.2199_DP  
      par_b_cm(1)=6.8657_DP   
      par_b_cm(2)=0.4385_DP
      par_b_cm(3)=26.8938_DP 
      par_b_cm(4)=102.478_DP
   CASE ('W')
      par_a_cm(1)=29.0818_DP   
      par_a_cm(2)=15.4300_DP   
      par_a_cm(3)=14.4327_DP   
      par_a_cm(4)=5.11982_DP   
      par_c_cm=9.8875_DP  
      par_b_cm(1)=1.72029_DP   
      par_b_cm(2)=9.2259_DP
      par_b_cm(3)=0.321703_DP  
      par_b_cm(4)=57.056_DP
   CASE ('Xe')
      par_a_cm(1)=20.2933_DP   
      par_a_cm(2)=19.0298_DP   
      par_a_cm(3)=8.9767_DP   
      par_a_cm(4)=1.99_DP   
      par_c_cm=3.7118_DP  
      par_b_cm(1)=3.9282_DP   
      par_b_cm(2)=0.344_DP
      par_b_cm(3)=26.4659_DP  
      par_b_cm(4)=64.2658_DP
   CASE ('Y')
      par_a_cm(1)=17.776_DP   
      par_a_cm(2)=10.2946_DP  
      par_a_cm(3)=5.72629_DP  
      par_a_cm(4)=3.26588_DP  
      par_c_cm=1.91213_DP 
      par_b_cm(1)=1.4029_DP  
      par_b_cm(2)=12.8006_DP
      par_b_cm(3)=0.125599_DP
      par_b_cm(4)=104.354_DP
   CASE ('Yb')
      par_a_cm(1)=28.6641_DP   
      par_a_cm(2)=15.4345_DP  
      par_a_cm(3)=15.3087_DP  
      par_a_cm(4)=2.98963_DP  
      par_c_cm=7.56672_DP 
      par_b_cm(1)=1.9889_DP  
      par_b_cm(2)=0.257119_DP
      par_b_cm(3)=10.6647_DP
      par_b_cm(4)=100.417_DP
   CASE ('Zn')
      par_a_cm(1)=14.0743_DP   
      par_a_cm(2)=7.0318_DP   
      par_a_cm(3)=5.1652_DP   
      par_a_cm(4)=2.41_DP   
      par_c_cm=1.3041_DP  
      par_b_cm(1)=3.2655_DP   
      par_b_cm(2)=0.2333_DP
      par_b_cm(3)=10.3163_DP  
      par_b_cm(4)=58.7097_DP
   CASE ('Zr')
      par_a_cm(1)=17.8765_DP   
      par_a_cm(2)=10.9480_DP   
      par_a_cm(3)=5.41732_DP   
      par_a_cm(4)=3.65721_DP   
      par_c_cm=2.06929_DP  
      par_b_cm(1)=1.27618_DP   
      par_b_cm(2)=11.916_DP
      par_b_cm(3)=0.117622_DP  
      par_b_cm(4)=87.6627_DP
CASE DEFAULT
      WRITE(stdout, '(" Warning: Atomic name not recognized:",a,& 
                          &".This atom will not scatter X-ray")') element
END SELECT

RETURN
END SUBROUTINE set_cromermann_coefficients

SUBROUTINE compute_form_factor(element, s, ff)
!
!  This routine receives the name of the elements and the parameter 
!  s = |G| / (4 pi) in A^-1 and gives as output the form factor
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=2), INTENT(IN) :: element
REAL(DP), INTENT(IN) :: s
REAL(DP), INTENT(OUT) :: ff

REAL(DP) :: par_ff_z, par_ff_a(4), par_ff_b(4)
INTEGER :: ipar

CALL compute_param_ff(element, par_ff_z, par_ff_a, par_ff_b)

ff=0.0_DP
DO ipar=1,4 
   ff= ff - 41.78214_DP * s**2 * par_ff_a(ipar) * EXP( -par_ff_b(ipar) * s**2 )
END DO

ff=ff+par_ff_z

RETURN
END SUBROUTINE compute_form_factor

SUBROUTINE compute_param_ff(element, par_ff_z, par_ff_a, par_ff_b)
!
!  This routine sets the form factor parameters. The parameters are
!  taken from Table 12.1 of
!  M. de Graef and M.E. McHenry Structure of materials, Cambridge (2007)
!
!  The parameters are from 
!  G.H. Smith and R.E. Burge, Acta Cryst. 15, 182 (1962).
!  P.A. Doyle and P.S. Turner, Acta Cryst. A24, 390 (1968).
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=2) :: element
REAL(DP), INTENT(INOUT) :: par_ff_z, par_ff_a(4), par_ff_b(4)

par_ff_z=0.0_DP
par_ff_a=0.0_DP
par_ff_b=0.0_DP

SELECT CASE (element) 

   CASE ('Ac') 
      par_ff_z=89.0_DP
      par_ff_a(1)=6.278_DP 
      par_ff_a(2)=5.195_DP 
      par_ff_a(3)=2.321_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.323_DP 
      par_ff_b(2)=4.949_DP 
      par_ff_b(3)=0.557_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Ag')
      par_ff_z=47.0_DP
      par_ff_a(1)=2.036_DP 
      par_ff_a(2)=3.272_DP 
      par_ff_a(3)=2.511_DP 
      par_ff_a(4)=0.837_DP
      par_ff_b(1)=61.497_DP 
      par_ff_b(2)=11.824_DP 
      par_ff_b(3)=2.846_DP 
      par_ff_b(4)=0.327_DP
   CASE ('Al')
      par_ff_z=13.0_DP
      par_ff_a(1)=2.276_DP 
      par_ff_a(2)=2.428_DP 
      par_ff_a(3)=0.858_DP 
      par_ff_a(4)=0.317_DP
      par_ff_b(1)=72.322_DP 
      par_ff_b(2)=19.773_DP 
      par_ff_b(3)=3.080_DP 
      par_ff_b(4)=0.408_DP
   CASE ('Am')
      par_ff_z=95.0_DP
      par_ff_a(1)=6.378_DP 
      par_ff_a(2)=5.495_DP 
      par_ff_a(3)=2.495_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.156_DP 
      par_ff_b(2)=5.102_DP 
      par_ff_b(3)=0.565_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Ar')
      par_ff_z=18.0_DP
      par_ff_a(1)=1.274_DP 
      par_ff_a(2)=2.190_DP 
      par_ff_a(3)=0.793_DP 
      par_ff_a(4)=0.326_DP
      par_ff_b(1)=26.682_DP 
      par_ff_b(2)=8.813_DP 
      par_ff_b(3)=2.219_DP 
      par_ff_b(4)=0.307_DP
   CASE ('As')
      par_ff_z=33.0_DP
      par_ff_a(1)=2.399_DP 
      par_ff_a(2)=2.790_DP 
      par_ff_a(3)=1.529_DP 
      par_ff_a(4)=0.594_DP
      par_ff_b(1)=45.718_DP 
      par_ff_b(2)=12.817_DP 
      par_ff_b(3)=2.280_DP 
      par_ff_b(4)=0.328_DP
   CASE ('At')
      par_ff_z=85.0_DP
      par_ff_a(1)=6.133_DP 
      par_ff_a(2)=5.031_DP 
      par_ff_a(3)=2.239_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.047_DP 
      par_ff_b(2)=4.957_DP 
      par_ff_b(3)=0.558_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Au')
      par_ff_z=79.0_DP
      par_ff_a(1)=2.388_DP 
      par_ff_a(2)=4.226_DP 
      par_ff_a(3)=2.689_DP 
      par_ff_a(4)=1.255_DP
      par_ff_b(1)=42.866_DP 
      par_ff_b(2)=9.743_DP 
      par_ff_b(3)=2.264_DP 
      par_ff_b(4)=0.307_DP
   CASE ('B')
      par_ff_z=5.0_DP
      par_ff_a(1)=0.945_DP 
      par_ff_a(2)=1.312_DP 
      par_ff_a(3)=0.419_DP 
      par_ff_a(4)=0.116_DP
      par_ff_b(1)=46.444_DP 
      par_ff_b(2)=14.178_DP 
      par_ff_b(3)=3.223_DP 
      par_ff_b(4)=0.377_DP
   CASE ('Ba')
      par_ff_z=56.0_DP
      par_ff_a(1)=7.821_DP 
      par_ff_a(2)=6.004_DP 
      par_ff_a(3)=3.280_DP 
      par_ff_a(4)=1.103_DP
      par_ff_b(1)=117.657_DP 
      par_ff_b(2)=18.778_DP 
      par_ff_b(3)=3.263_DP 
      par_ff_b(4)=0.376_DP
   CASE ('Be')
      par_ff_z=4.0_DP
      par_ff_a(1)=1.250_DP 
      par_ff_a(2)=1.334_DP 
      par_ff_a(3)=0.360_DP 
      par_ff_a(4)=0.106_DP
      par_ff_b(1)=60.804_DP 
      par_ff_b(2)=18.591_DP 
      par_ff_b(3)=3.653_DP 
      par_ff_b(4)=0.416_DP
   CASE ('Bi')
      par_ff_z=83.0_DP
      par_ff_a(1)=3.841_DP 
      par_ff_a(2)=4.679_DP 
      par_ff_a(3)=3.192_DP 
      par_ff_a(4)=1.363_DP
      par_ff_b(1)=50.261_DP 
      par_ff_b(2)=11.999_DP 
      par_ff_b(3)=2.560_DP 
      par_ff_b(4)=0.318_DP
   CASE ('Bk')
      par_ff_z=97.0_DP
      par_ff_a(1)=6.502_DP 
      par_ff_a(2)=5.478_DP 
      par_ff_a(3)=2.510_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.375_DP 
      par_ff_b(2)=4.975_DP 
      par_ff_b(3)=0.561_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Br')
      par_ff_z=35.0_DP
      par_ff_a(1)=2.166_DP 
      par_ff_a(2)=2.904_DP 
      par_ff_a(3)=1.395_DP 
      par_ff_a(4)=0.589_DP
      par_ff_b(1)=33.899_DP 
      par_ff_b(2)=10.497_DP 
      par_ff_b(3)=2.041_DP 
      par_ff_b(4)=0.307_DP
   CASE ('C')
      par_ff_z=6.0_DP
      par_ff_a(1)=0.731_DP 
      par_ff_a(2)=1.195_DP 
      par_ff_a(3)=0.456_DP 
      par_ff_a(4)=0.125_DP
      par_ff_b(1)=36.995_DP 
      par_ff_b(2)=11.297_DP 
      par_ff_b(3)=2.814_DP 
      par_ff_b(4)=0.346_DP
   CASE ('Ca')
      par_ff_z=20.0_DP
      par_ff_a(1)=4.470_DP 
      par_ff_a(2)=2.971_DP 
      par_ff_a(3)=1.970_DP 
      par_ff_a(4)=0.482_DP
      par_ff_b(1)=99.523_DP 
      par_ff_b(2)=22.696_DP 
      par_ff_b(3)=4.195_DP 
      par_ff_b(4)=0.417_DP
   CASE ('Cd')
      par_ff_z=48.0_DP
      par_ff_a(1)=2.574_DP 
      par_ff_a(2)=3.259_DP 
      par_ff_a(3)=2.547_DP 
      par_ff_a(4)=0.838_DP
      par_ff_b(1)=55.675_DP 
      par_ff_b(2)=11.838_DP 
      par_ff_b(3)=2.784_DP 
      par_ff_b(4)=0.322_DP
   CASE ('Ce')
      par_ff_z=58.0_DP
      par_ff_a(1)=5.007_DP 
      par_ff_a(2)=3.980_DP 
      par_ff_a(3)=1.678_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.283_DP 
      par_ff_b(2)=5.183_DP 
      par_ff_b(3)=0.589_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Cf')
      par_ff_z=98.0_DP
      par_ff_a(1)=6.548_DP 
      par_ff_a(2)=5.526_DP 
      par_ff_a(3)=2.520_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.461_DP 
      par_ff_b(2)=4.965_DP 
      par_ff_b(3)=0.557_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Cl')
      par_ff_z=17.0_DP
      par_ff_a(1)=1.452_DP 
      par_ff_a(2)=2.292_DP 
      par_ff_a(3)=0.787_DP 
      par_ff_a(4)=0.322_DP
      par_ff_b(1)=30.935_DP 
      par_ff_b(2)=9.980_DP 
      par_ff_b(3)=2.234_DP 
      par_ff_b(4)=0.323_DP
   CASE ('Cm')
      par_ff_z=96.0_DP
      par_ff_a(1)=6.460_DP 
      par_ff_a(2)=5.469_DP 
      par_ff_a(3)=2.471_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.396_DP 
      par_ff_b(2)=4.970_DP 
      par_ff_b(3)=0.554_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Co')
      par_ff_z=27.0_DP
      par_ff_a(1)=2.367_DP 
      par_ff_a(2)=2.236_DP 
      par_ff_a(3)=1.724_DP 
      par_ff_a(4)=0.515_DP
      par_ff_b(1)=61.431_DP 
      par_ff_b(2)=14.180_DP 
      par_ff_b(3)=2.725_DP 
      par_ff_b(4)=0.344_DP
   CASE ('Cr')
      par_ff_z=24.0_DP
      par_ff_a(1)=2.307_DP 
      par_ff_a(2)=2.334_DP 
      par_ff_a(3)=1.823_DP 
      par_ff_a(4)=0.490_DP
      par_ff_b(1)=78.405_DP 
      par_ff_b(2)=15.785_DP 
      par_ff_b(3)=3.157_DP 
      par_ff_b(4)=0.364_DP
   CASE ('Cs')
      par_ff_z=55.0_DP
      par_ff_a(1)=6.062_DP 
      par_ff_a(2)=5.986_DP 
      par_ff_a(3)=3.303_DP 
      par_ff_a(4)=1.096_DP
      par_ff_b(1)=155.837_DP 
      par_ff_b(2)=19.695_DP 
      par_ff_b(3)=3.335_DP 
      par_ff_b(4)=0.379_DP
   CASE ('Cu')
      par_ff_z=29.0_DP
      par_ff_a(1)=1.579_DP 
      par_ff_a(2)=1.820_DP 
      par_ff_a(3)=1.658_DP 
      par_ff_a(4)=0.532_DP
      par_ff_b(1)=62.940_DP 
      par_ff_b(2)=12.453_DP 
      par_ff_b(3)=2.504_DP 
      par_ff_b(4)=0.333_DP
   CASE ('Dy')
      par_ff_z=66.0_DP
      par_ff_a(1)=5.332_DP 
      par_ff_a(2)=4.370_DP 
      par_ff_a(3)=1.863_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.888_DP 
      par_ff_b(2)=5.198_DP 
      par_ff_b(3)=0.581_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Er')
      par_ff_z=68.0_DP
      par_ff_a(1)=5.436_DP 
      par_ff_a(2)=4.437_DP 
      par_ff_a(3)=1.891_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.655_DP 
      par_ff_b(2)=5.117_DP 
      par_ff_b(3)=0.577_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Eu')
      par_ff_z=63.0_DP
      par_ff_a(1)=6.267_DP 
      par_ff_a(2)=4.844_DP 
      par_ff_a(3)=3.202_DP 
      par_ff_a(4)=1.200_DP
      par_ff_b(1)=100.298_DP 
      par_ff_b(2)=16.066_DP 
      par_ff_b(3)=2.980_DP 
      par_ff_b(4)=0.367_DP
   CASE ('F')
      par_ff_z=9.0_DP
      par_ff_a(1)=0.387_DP 
      par_ff_a(2)=0.811_DP 
      par_ff_a(3)=0.475_DP 
      par_ff_a(4)=0.146_DP
      par_ff_b(1)=20.239_DP 
      par_ff_b(2)=6.609_DP 
      par_ff_b(3)=1.931_DP 
      par_ff_b(4)=0.279_DP
   CASE ('Fe')
      par_ff_z=26.0_DP
      par_ff_a(1)=2.544_DP 
      par_ff_a(2)=2.343_DP 
      par_ff_a(3)=1.759_DP 
      par_ff_a(4)=0.506_DP
      par_ff_b(1)=64.424_DP 
      par_ff_b(2)=14.880_DP 
      par_ff_b(3)=2.854_DP 
      par_ff_b(4)=0.350_DP
   CASE ('Fr')
      par_ff_z=87.0_DP
      par_ff_a(1)=6.201_DP 
      par_ff_a(2)=5.121_DP 
      par_ff_a(3)=2.275_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.200_DP 
      par_ff_b(2)=4.954_DP 
      par_ff_b(3)=0.556_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Ga')
      par_ff_z=31.0_DP
      par_ff_a(1)=2.321_DP 
      par_ff_a(2)=2.486_DP 
      par_ff_a(3)=1.688_DP 
      par_ff_a(4)=0.599_DP
      par_ff_b(1)=65.602_DP 
      par_ff_b(2)=15.458_DP 
      par_ff_b(3)=2.581_DP 
      par_ff_b(4)=0.351_DP
   CASE ('Gd')
      par_ff_z=64.0_DP
      par_ff_a(1)=5.225_DP 
      par_ff_a(2)=4.314_DP 
      par_ff_a(3)=1.827_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.158_DP 
      par_ff_b(2)=5.259_DP 
      par_ff_b(3)=0.586_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Ge')
      par_ff_z=32.0_DP
      par_ff_a(1)=2.447_DP 
      par_ff_a(2)=2.702_DP 
      par_ff_a(3)=1.616_DP 
      par_ff_a(4)=0.601_DP
      par_ff_b(1)=55.893_DP 
      par_ff_b(2)=14.393_DP 
      par_ff_b(3)=2.446_DP 
      par_ff_b(4)=0.342_DP
   CASE ('H')
      par_ff_z=1.0_DP
      par_ff_a(1)=0.202_DP 
      par_ff_a(2)=0.244_DP 
      par_ff_a(3)=0.082_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=30.868_DP 
      par_ff_b(2)=8.544_DP 
      par_ff_b(3)=1.273_DP 
      par_ff_b(4)=0._DP
   CASE ('He')
      par_ff_z=2.0_DP
      par_ff_a(1)=0.091_DP 
      par_ff_a(2)=0.181_DP 
      par_ff_a(3)=0.110_DP 
      par_ff_a(4)=0.036_DP
      par_ff_b(1)=18.183_DP 
      par_ff_b(2)=6.212_DP 
      par_ff_b(3)=1.803_DP 
      par_ff_b(4)=0.284_DP
   CASE ('Hf')
      par_ff_z=72.0_DP
      par_ff_a(1)=5.588_DP 
      par_ff_a(2)=4.619_DP 
      par_ff_a(3)=1.997_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.001_DP 
      par_ff_b(2)=5.164_DP 
      par_ff_b(3)=0.579_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Hg')
      par_ff_z=80.0_DP
      par_ff_a(1)=2.682_DP 
      par_ff_a(2)=4.241_DP 
      par_ff_a(3)=2.755_DP 
      par_ff_a(4)=1.270_DP
      par_ff_b(1)=42.822_DP 
      par_ff_b(2)=9.856_DP 
      par_ff_b(3)=2.295_DP 
      par_ff_b(4)=0.307_DP
   CASE ('Ho')
      par_ff_z=67.0_DP
      par_ff_a(1)=5.376_DP 
      par_ff_a(2)=4.403_DP 
      par_ff_a(3)=1.884_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.773_DP 
      par_ff_b(2)=5.174_DP 
      par_ff_b(3)=0.582_DP 
      par_ff_b(4)=0.0_DP
   CASE ('I')
      par_ff_z=53.0_DP
      par_ff_a(1)=3.473_DP 
      par_ff_a(2)=4.060_DP 
      par_ff_a(3)=2.522_DP 
      par_ff_a(4)=0.840_DP
      par_ff_b(1)=39.441_DP 
      par_ff_b(2)=11.816_DP 
      par_ff_b(3)=2.415_DP 
      par_ff_b(4)=0.298_DP
   CASE ('In')
      par_ff_z=49.0_DP
      par_ff_a(1)=3.153_DP 
      par_ff_a(2)=3.557_DP 
      par_ff_a(3)=2.818_DP 
      par_ff_a(4)=0.884_DP
      par_ff_b(1)=66.649_DP 
      par_ff_b(2)=14.449_DP 
      par_ff_b(3)=2.976_DP 
      par_ff_b(4)=0.335_DP
   CASE ('Ir')
      par_ff_z=77.0_DP
      par_ff_a(1)=5.754_DP 
      par_ff_a(2)=4.851_DP 
      par_ff_a(3)=2.096_DP 
      par_ff_a(4)=0._DP
      par_ff_b(1)=29.159_DP 
      par_ff_b(2)=5.152_DP 
      par_ff_b(3)=0.570_DP 
      par_ff_b(4)=0._DP
   CASE ('K')
      par_ff_z=19.0_DP
      par_ff_a(1)=3.951_DP 
      par_ff_a(2)=2.545_DP 
      par_ff_a(3)=1.980_DP 
      par_ff_a(4)=0.482_DP
      par_ff_b(1)=137.075_DP 
      par_ff_b(2)=22.402_DP 
      par_ff_b(3)=4.532_DP 
      par_ff_b(4)=0.434_DP
   CASE ('Kr')
      par_ff_z=36.0_DP
      par_ff_a(1)=2.034_DP 
      par_ff_a(2)=2.927_DP 
      par_ff_a(3)=1.342_DP 
      par_ff_a(4)=0.589_DP
      par_ff_b(1)=29.999_DP 
      par_ff_b(2)=9.598_DP 
      par_ff_b(3)=1.952_DP 
      par_ff_b(4)=0.299_DP
   CASE ('La')
      par_ff_z=57.0_DP
      par_ff_a(1)=4.940_DP 
      par_ff_a(2)=3.968_DP 
      par_ff_a(3)=1.663_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.716_DP 
      par_ff_b(2)=5.245_DP 
      par_ff_b(3)=0.594_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Li')
      par_ff_z=3.0_DP
      par_ff_a(1)=1.611_DP 
      par_ff_a(2)=1.246_DP 
      par_ff_a(3)=0.326_DP 
      par_ff_a(4)=0.099_DP
      par_ff_b(1)=107.638_DP 
      par_ff_b(2)=30.480_DP 
      par_ff_b(3)=4.533_DP 
      par_ff_b(4)=0.495_DP
   CASE ('Lu')
      par_ff_z=71.0_DP
      par_ff_a(1)=5.553_DP 
      par_ff_a(2)=4.580_DP 
      par_ff_a(3)=1.969_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.907_DP 
      par_ff_b(2)=5.160_DP 
      par_ff_b(3)=0.577_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Mg')
      par_ff_z=12.0_DP
      par_ff_a(1)=2.268_DP 
      par_ff_a(2)=1.803_DP 
      par_ff_a(3)=0.839_DP 
      par_ff_a(4)=0.289_DP
      par_ff_b(1)=73.670_DP 
      par_ff_b(2)=20.175_DP 
      par_ff_b(3)=3.013_DP 
      par_ff_b(4)=0.405_DP
   CASE ('Mn')
      par_ff_z=25.0_DP
      par_ff_a(1)=2.747_DP 
      par_ff_a(2)=2.456_DP 
      par_ff_a(3)=1.792_DP 
      par_ff_a(4)=0.498_DP
      par_ff_b(1)=67.786_DP 
      par_ff_b(2)=15.674_DP 
      par_ff_b(3)=3.000_DP 
      par_ff_b(4)=0.357_DP
   CASE ('Mo')
      par_ff_z=42.0_DP
      par_ff_a(1)=3.120_DP 
      par_ff_a(2)=3.906_DP 
      par_ff_a(3)=2.361_DP 
      par_ff_a(4)=0.850_DP
      par_ff_b(1)=72.464_DP 
      par_ff_b(2)=14.642_DP 
      par_ff_b(3)=3.237_DP 
      par_ff_b(4)=0.366_DP
   CASE ('N')
      par_ff_z=7.0_DP
      par_ff_a(1)=0.572_DP 
      par_ff_a(2)=1.043_DP 
      par_ff_a(3)=0.465_DP 
      par_ff_a(4)=0.131_DP
      par_ff_b(1)=28.847_DP 
      par_ff_b(2)=9.054_DP 
      par_ff_b(3)=2.421_DP 
      par_ff_b(4)=0.317_DP
   CASE ('Na')
      par_ff_z=11.0_DP
      par_ff_a(1)=2.241_DP 
      par_ff_a(2)=1.333_DP 
      par_ff_a(3)=0.907_DP 
      par_ff_a(4)=0.286_DP
      par_ff_b(1)=108.004_DP 
      par_ff_b(2)=24.505_DP 
      par_ff_b(3)=3.391_DP 
      par_ff_b(4)=0.435_DP
   CASE ('Nb')
      par_ff_z=41.0_DP
      par_ff_a(1)=4.237_DP 
      par_ff_a(2)=3.105_DP 
      par_ff_a(3)=1.234_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=27.415_DP 
      par_ff_b(2)=5.074_DP 
      par_ff_b(3)=0.593_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Nd')
      par_ff_z=60.0_DP
      par_ff_a(1)=5.151_DP 
      par_ff_a(2)=4.075_DP 
      par_ff_a(3)=1.683_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.304_DP 
      par_ff_b(2)=5.073_DP 
      par_ff_b(3)=0.571_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Ne')
      par_ff_z=10.0_DP
      par_ff_a(1)=0.303_DP 
      par_ff_a(2)=0.720_DP 
      par_ff_a(3)=0.475_DP 
      par_ff_a(4)=0.153_DP
      par_ff_b(1)=17.640_DP 
      par_ff_b(2)=5.860_DP 
      par_ff_b(3)=1.762_DP 
      par_ff_b(4)=0.266_DP
   CASE ('Ni')
      par_ff_z=28.0_DP
      par_ff_a(1)=2.210_DP 
      par_ff_a(2)=2.134_DP 
      par_ff_a(3)=1.689_DP 
      par_ff_a(4)=0.524_DP
      par_ff_b(1)=58.727_DP 
      par_ff_b(2)=13.553_DP 
      par_ff_b(3)=2.609_DP 
      par_ff_b(4)=0.339_DP
   CASE ('Np')
      par_ff_z=93.0_DP
      par_ff_a(1)=6.323_DP 
      par_ff_a(2)=5.414_DP 
      par_ff_a(3)=2.453_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.142_DP 
      par_ff_b(2)=5.096_DP 
      par_ff_b(3)=0.568_DP 
      par_ff_b(4)=0.0_DP
   CASE ('O')
      par_ff_z=8.0_DP
      par_ff_a(1)=0.455_DP 
      par_ff_a(2)=0.917_DP 
      par_ff_a(3)=0.472_DP 
      par_ff_a(4)=0.138_DP
      par_ff_b(1)=23.780_DP 
      par_ff_b(2)=7.622_DP 
      par_ff_b(3)=2.144_DP 
      par_ff_b(4)=0.296_DP
   CASE ('Os')
      par_ff_z=76.0_DP
      par_ff_a(1)=5.750_DP 
      par_ff_a(2)=4.773_DP 
      par_ff_a(3)=2.079_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.933_DP 
      par_ff_b(2)=5.139_DP 
      par_ff_b(3)=0.573_DP 
      par_ff_b(4)=0.0_DP
   CASE ('P')
      par_ff_z=15.0_DP
      par_ff_a(1)=1.888_DP 
      par_ff_a(2)=2.469_DP 
      par_ff_a(3)=0.805_DP 
      par_ff_a(4)=0.320_DP
      par_ff_b(1)=44.876_DP 
      par_ff_b(2)=13.538_DP 
      par_ff_b(3)=2.642_DP 
      par_ff_b(4)=0.361_DP
   CASE ('Pa')
      par_ff_z=91.0_DP
      par_ff_a(1)=6.306_DP 
      par_ff_a(2)=5.303_DP 
      par_ff_a(3)=2.386_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.688_DP 
      par_ff_b(2)=5.026_DP 
      par_ff_b(3)=0.561_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Pb')
      par_ff_z=82.0_DP
      par_ff_a(1)=3.510_DP 
      par_ff_a(2)=4.552_DP 
      par_ff_a(3)=3.154_DP 
      par_ff_a(4)=1.359_DP
      par_ff_b(1)=52.914_DP 
      par_ff_b(2)=11.884_DP 
      par_ff_b(3)=2.571_DP 
      par_ff_b(4)=0.321_DP
   CASE ('Pd')
      par_ff_z=46.0_DP
      par_ff_a(1)=4.436_DP 
      par_ff_a(2)=3.454_DP 
      par_ff_a(3)=1.383_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.670_DP 
      par_ff_b(2)=5.269_DP 
      par_ff_b(3)=0.595_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Pm')
      par_ff_z=61.0_DP
      par_ff_a(1)=5.201_DP 
      par_ff_a(2)=4.094_DP 
      par_ff_a(3)=1.719_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.079_DP 
      par_ff_b(2)=5.081_DP 
      par_ff_b(3)=0.576_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Po')
      par_ff_z=84.0_DP
      par_ff_a(1)=6.070_DP 
      par_ff_a(2)=4.997_DP 
      par_ff_a(3)=2.232_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.075_DP 
      par_ff_b(2)=4.999_DP 
      par_ff_b(3)=0.563_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Pr')
      par_ff_z=59.0_DP
      par_ff_a(1)=5.085_DP 
      par_ff_a(2)=4.043_DP 
      par_ff_a(3)=1.684_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.588_DP 
      par_ff_b(2)=5.143_DP 
      par_ff_b(3)=0.581_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Pt')
      par_ff_z=78.0_DP
      par_ff_a(1)=5.803_DP 
      par_ff_a(2)=4.870_DP 
      par_ff_a(3)=2.127_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.016_DP 
      par_ff_b(2)=5.150_DP 
      par_ff_b(3)=0.572_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Pu')
      par_ff_z=94.0_DP
      par_ff_a(1)=6.415_DP 
      par_ff_a(2)=5.419_DP 
      par_ff_a(3)=2.449_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.836_DP 
      par_ff_b(2)=5.022_DP 
      par_ff_b(3)=0.561_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Ra')
      par_ff_z=88.0_DP
      par_ff_a(1)=6.215_DP 
      par_ff_a(2)=5.170_DP 
      par_ff_a(3)=2.316_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.382_DP 
      par_ff_b(2)=5.002_DP 
      par_ff_b(3)=0.562_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Rb')
      par_ff_z=37.0_DP
      par_ff_a(1)=4.776_DP 
      par_ff_a(2)=3.859_DP 
      par_ff_a(3)=2.234_DP 
      par_ff_a(4)=0.868_DP
      par_ff_b(1)=140.782_DP 
      par_ff_b(2)=18.991_DP 
      par_ff_b(3)=3.701_DP 
      par_ff_b(4)=0.419_DP
   CASE ('Re')
      par_ff_z=75.0_DP
      par_ff_a(1)=5.695_DP 
      par_ff_a(2)=4.740_DP 
      par_ff_a(3)=2.064_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.968_DP 
      par_ff_b(2)=5.156_DP 
      par_ff_b(3)=0.575_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Rh')
      par_ff_z=45.0_DP
      par_ff_a(1)=4.431_DP 
      par_ff_a(2)=3.343_DP 
      par_ff_a(3)=1.345_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=27.911_DP 
      par_ff_b(2)=5.153_DP 
      par_ff_b(3)=0.592_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Rn')
      par_ff_z=86.0_DP
      par_ff_a(1)=4.078_DP 
      par_ff_a(2)=4.978_DP 
      par_ff_a(3)=3.096_DP 
      par_ff_a(4)=1.326_DP
      par_ff_b(1)=38.406_DP 
      par_ff_b(2)=11.020_DP 
      par_ff_b(3)=2.355_DP 
      par_ff_b(4)=0.299_DP
   CASE ('Ru')
      par_ff_z=44.0_DP
      par_ff_a(1)=4.358_DP 
      par_ff_a(2)=3.298_DP 
      par_ff_a(3)=1.323_DP 
      par_ff_a(4)=0._DP
      par_ff_b(1)=27.881_DP 
      par_ff_b(2)=5.179_DP 
      par_ff_b(3)=0.594_DP 
      par_ff_b(4)=0.0_DP
   CASE ('S')
      par_ff_z=16.0_DP
      par_ff_a(1)=1.659_DP 
      par_ff_a(2)=2.386_DP 
      par_ff_a(3)=0.790_DP 
      par_ff_a(4)=0.321_DP
      par_ff_b(1)=36.650_DP 
      par_ff_b(2)=11.488_DP 
      par_ff_b(3)=2.469_DP 
      par_ff_b(4)=0.340_DP
   CASE ('Sb')
      par_ff_z=51.0_DP
      par_ff_a(1)=3.564_DP 
      par_ff_a(2)=3.844_DP 
      par_ff_a(3)=2.687_DP 
      par_ff_a(4)=0.864_DP
      par_ff_b(1)=50.487_DP 
      par_ff_b(2)=13.316_DP 
      par_ff_b(3)=2.691_DP 
      par_ff_b(4)=0.316_DP
   CASE ('Sc')
      par_ff_z=21.0_DP
      par_ff_a(1)=3.966_DP 
      par_ff_a(2)=2.917_DP 
      par_ff_a(3)=1.925_DP 
      par_ff_a(4)=0.480_DP
      par_ff_b(1)=88.960_DP 
      par_ff_b(2)=20.606_DP 
      par_ff_b(3)=3.856_DP 
      par_ff_b(4)=0.399_DP
   CASE ('Se')
      par_ff_z=34.0_DP
      par_ff_a(1)=2.298_DP 
      par_ff_a(2)=2.854_DP 
      par_ff_a(3)=1.456_DP 
      par_ff_a(4)=0.590_DP
      par_ff_b(1)=38.830_DP 
      par_ff_b(2)=11.536_DP 
      par_ff_b(3)=2.146_DP 
      par_ff_b(4)=0.316_DP
   CASE ('Si')
      par_ff_z=14.0_DP
      par_ff_a(1)=2.129_DP 
      par_ff_a(2)=2.533_DP 
      par_ff_a(3)=0.835_DP 
      par_ff_a(4)=0.322_DP
      par_ff_b(1)=57.775_DP 
      par_ff_b(2)=16.476_DP 
      par_ff_b(3)=2.880_DP 
      par_ff_b(4)=0.386_DP
   CASE ('Sm')
      par_ff_z=62.0_DP
      par_ff_a(1)=5.255_DP 
      par_ff_a(2)=4.113_DP 
      par_ff_a(3)=1.743_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.016_DP 
      par_ff_b(2)=5.037_DP 
      par_ff_b(3)=0.577_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Sn')
      par_ff_z=50.0_DP
      par_ff_a(1)=3.450_DP 
      par_ff_a(2)=3.735_DP 
      par_ff_a(3)=2.118_DP 
      par_ff_a(4)=0.877_DP
      par_ff_b(1)=59.104_DP 
      par_ff_b(2)=14.179_DP 
      par_ff_b(3)=2.855_DP 
      par_ff_b(4)=0.327_DP
   CASE ('Sr')
      par_ff_z=38.0_DP
      par_ff_a(1)=5.848_DP 
      par_ff_a(2)=4.003_DP 
      par_ff_a(3)=2.342_DP 
      par_ff_a(4)=0.880_DP
      par_ff_b(1)=104.972_DP 
      par_ff_b(2)=19.367_DP 
      par_ff_b(3)=3.737_DP 
      par_ff_b(4)=0.414_DP
   CASE ('Ta')
      par_ff_z=73.0_DP
      par_ff_a(1)=5.659_DP 
      par_ff_a(2)=4.630_DP 
      par_ff_a(3)=2.014_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.807_DP 
      par_ff_b(2)=5.114_DP 
      par_ff_b(3)=0.578_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Tb')
      par_ff_z=65.0_DP
      par_ff_a(1)=5.272_DP 
      par_ff_a(2)=4.347_DP 
      par_ff_a(3)=1.844_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.046_DP 
      par_ff_b(2)=5.226_DP 
      par_ff_b(3)=0.585_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Tc')
      par_ff_z=43.0_DP
      par_ff_a(1)=4.318_DP 
      par_ff_a(2)=3.270_DP 
      par_ff_a(3)=1.287_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.246_DP 
      par_ff_b(2)=5.148_DP 
      par_ff_b(3)=0.590_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Te')
      par_ff_z=52.0_DP
      par_ff_a(1)=4.785_DP 
      par_ff_a(2)=3.688_DP 
      par_ff_a(3)=1.500_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=27.999_DP 
      par_ff_b(2)=5.083_DP 
      par_ff_b(3)=0.581_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Th')
      par_ff_z=90.0_DP
      par_ff_a(1)=6.264_DP 
      par_ff_a(2)=5.263_DP 
      par_ff_a(3)=2.367_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.651_DP 
      par_ff_b(2)=5.030_DP 
      par_ff_b(3)=0.563_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Ti')
      par_ff_z=22.0_DP
      par_ff_a(1)=3.565_DP 
      par_ff_a(2)=2.818_DP 
      par_ff_a(3)=1.893_DP 
      par_ff_a(4)=0.483_DP
      par_ff_b(1)=81.982_DP 
      par_ff_b(2)=19.049_DP 
      par_ff_b(3)=3.590_DP 
      par_ff_b(4)=0.386_DP
   CASE ('Tl')
      par_ff_z=81.0_DP
      par_ff_a(1)=5.932_DP 
      par_ff_a(2)=4.972_DP 
      par_ff_a(3)=2.195_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.086_DP 
      par_ff_b(2)=5.126_DP 
      par_ff_b(3)=0.572_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Tm')
      par_ff_z=69.0_DP
      par_ff_a(1)=5.441_DP 
      par_ff_a(2)=4.510_DP 
      par_ff_a(3)=1.956_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=29.149_DP 
      par_ff_b(2)=5.264_DP 
      par_ff_b(3)=0.590_DP 
      par_ff_b(4)=0.0_DP
   CASE ('U')
      par_ff_z=92.0_DP
      par_ff_a(1)=6.767_DP 
      par_ff_a(2)=6.729_DP 
      par_ff_a(3)=4.014_DP 
      par_ff_a(4)=1.561_DP
      par_ff_b(1)=85.951_DP 
      par_ff_b(2)=15.642_DP 
      par_ff_b(3)=2.936_DP 
      par_ff_b(4)=0.335_DP
   CASE ('V')
      par_ff_z=23.0_DP
      par_ff_a(1)=3.245_DP 
      par_ff_a(2)=2.698_DP 
      par_ff_a(3)=1.860_DP 
      par_ff_a(4)=0.486_DP
      par_ff_b(1)=76.379_DP 
      par_ff_b(2)=17.726_DP 
      par_ff_b(3)=3.363_DP 
      par_ff_b(4)=0.374_DP
   CASE ('W')
      par_ff_z=74.0_DP
      par_ff_a(1)=5.709_DP 
      par_ff_a(2)=4.677_DP 
      par_ff_a(3)=2.019_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.782_DP 
      par_ff_b(2)=5.084_DP 
      par_ff_b(3)=0.572_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Xe')
      par_ff_z=54.0_DP
      par_ff_a(1)=3.366_DP 
      par_ff_a(2)=4.147_DP 
      par_ff_a(3)=2.443_DP 
      par_ff_a(4)=0.829_DP
      par_ff_b(1)=35.509_DP 
      par_ff_b(2)=11.117_DP 
      par_ff_b(3)=2.294_DP 
      par_ff_b(4)=0.289_DP
   CASE ('Y')
      par_ff_z=39.0_DP
      par_ff_a(1)=4.129_DP 
      par_ff_a(2)=3.012_DP 
      par_ff_a(3)=1.179_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=27.548_DP 
      par_ff_b(2)=5.088_DP 
      par_ff_b(3)=0.591_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Yb')
      par_ff_z=70.0_DP
      par_ff_a(1)=5.529_DP 
      par_ff_a(2)=4.533_DP 
      par_ff_a(3)=1.945_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.927_DP 
      par_ff_b(2)=5.144_DP 
      par_ff_b(3)=0.578_DP 
      par_ff_b(4)=0.0_DP
   CASE ('Zn')
      par_ff_z=30.0_DP
      par_ff_a(1)=1.942_DP 
      par_ff_a(2)=1.950_DP 
      par_ff_a(3)=1.619_DP 
      par_ff_a(4)=0.543_DP
      par_ff_b(1)=54.162_DP 
      par_ff_b(2)=12.518_DP 
      par_ff_b(3)=2.416_DP 
      par_ff_b(4)=0.330_DP
   CASE ('Zr')
      par_ff_z=40.0_DP
      par_ff_a(1)=4.105_DP 
      par_ff_a(2)=3.144_DP 
      par_ff_a(3)=1.229_DP 
      par_ff_a(4)=0.0_DP
      par_ff_b(1)=28.492_DP 
      par_ff_b(2)=5.277_DP 
      par_ff_b(3)=0.601_DP 
      par_ff_b(4)=0.0_DP
   CASE DEFAULT
        WRITE(stdout, '(" Warning: Atomic name not recognized:",a,& 
                          &".This atom will not scatter X-ray")') element
END SELECT

RETURN
END SUBROUTINE compute_param_ff
!

END MODULE xrdp_module
