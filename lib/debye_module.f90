!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE debye_module
!
!  This module provides routines to compute the thermodynamical quantities
!  within the Debye model.
!  It has routines to:
!  
!  compute the mean sound velocity from the elastic constants
!  compute the Debye temperature using the elastic constants
!  compute the specific heat given the Debye temperature
!  compute the vibrational energy given the Debye temperature
!  compute the vibrational free energy given the Debye temperature
!  compute the vibrational entropy given the Debye temperature
!
!
USE kinds, ONLY : DP
USE constants, ONLY : k_boltzmann_si, h_planck_si, bohr_radius_si, pi, &
                      k_boltzmann_ry
USE io_global, ONLY : stdout
USE mp_images, ONLY : intra_image_comm
USE mp, ONLY : mp_sum
IMPLICIT NONE
SAVE
PRIVATE
INTEGER, PARAMETER :: npt=1000  ! the number of points used to make
                                ! the Debye integral

PUBLIC  compute_debye_temperature, compute_average_sound, debye_vib_energy, &
        compute_debye_temperature_macro_el,                   &
        debye_free_energy, debye_entropy, debye_cv, debye_e0, &
        compute_debye_temperature_poisson, debye_b_factor, &
        debye_free_energy_0d, debye_cv_0d, debye_energy_0d, &
        write_debye_on_file

CONTAINS

!--------------------------------------------------------------
SUBROUTINE compute_average_sound(elconv, density, average_sound_speed)
!--------------------------------------------------------------
!
!  This routine receives as input the elastic constants in Voigt notation 
!  C_{ij} and the density of the solid. It gives as output the average 
!  sound speed.
!
USE kinds, ONLY : DP
USE constants, ONLY : pi
USE elastic_constants, ONLY : compute_sound
USE voigt, ONLY : to_voigt4
IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: elconv(6,6)
REAL(DP), INTENT(IN) :: density

REAL(DP), INTENT(OUT) :: average_sound_speed

REAL(DP):: sound_speed(3)
REAL(DP):: sound_disp(3,3)
REAL(DP) :: elcon(3,3,3,3)
REAL(DP) :: qvec(3), theta, phi, delta_theta, delta_phi, integral
INTEGER :: npt_theta, npt_phi
INTEGER :: i, j, k, l
!
!  Transform the elastic constants in the tensor notation
!
CALL to_voigt4(elconv, elcon, .FALSE.)
!
!  This is a two dimensional integral
!
npt_theta=40
npt_phi=40
delta_theta= pi / npt_theta
delta_phi= 2.0_DP * pi / npt_phi
integral=0.0_DP
DO i=1,npt_theta
   theta = delta_theta * (i-1) + delta_theta / 2.0_DP
   DO j=1,npt_phi
      phi = delta_phi * (j-1)
      qvec(1) = SIN(theta) * COS(phi)
      qvec(2) = SIN(theta) * SIN(phi)
      qvec(3) = COS(theta) 
!      WRITE(stdout,*) 'qvec', qvec(1), qvec(2), qvec(3)
      CALL compute_sound(elcon, qvec, density, sound_speed, sound_disp)
!      WRITE(stdout,*) 'sound speed', sound_speed(1), sound_speed(2), 
!                       sound_speed(3)
      IF (sound_speed(1)>0.0_DP.AND.sound_speed(2)>0.0_DP.AND. &
          sound_speed(3)>0.0_DP) &
      integral = integral + SIN(theta) * ( 1.0_DP / sound_speed(1)**3 + &
                  1.0_DP / sound_speed(2)**3 + 1.0_DP / sound_speed(3)**3 )
   END DO
END DO
integral= integral * delta_theta * delta_phi / 12.0_DP / pi
IF (integral > 0.0_DP) THEN
   average_sound_speed=1.0_DP / integral**(1.0_DP/3.0_DP)
   WRITE(stdout,'(/,20x,40("-"),/)')
   WRITE(stdout,'(5x,"Average Debye sound velocity = ", f12.3, " m/s")') &
                                                    average_sound_speed 
ELSE
   WRITE(stdout,'(5x,"Average Debye sound velocity not available")') 
ENDIF
RETURN
END SUBROUTINE compute_average_sound

!-------------------------------------------------------------------------
SUBROUTINE compute_debye_temperature(el_con, density, nat, omega, debye_t)
!-------------------------------------------------------------------------
!
!  This routine receives as input: the elastic constants in Voigt form
!  and in kbar units, the density in kg/m^3, omega in (a.u.)^3, the
!  number of atoms per cell and gives as output the Debye temperature
!  in K. The Debye temperature is also written in output.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: el_con(6,6)
REAL(DP), INTENT(IN) :: density
REAL(DP), INTENT(IN) :: omega
REAL(DP), INTENT(OUT) :: debye_t

REAL(DP) :: fact
REAL(DP) :: average_sound_speed
!
!  first compute the average sound speed
!
CALL compute_average_sound(el_con, density, average_sound_speed)

fact = h_planck_si / k_boltzmann_si / bohr_radius_si
debye_t = average_sound_speed * fact *    &
                       (3.0_DP / pi /4.0_DP * nat / omega)**(1.0_DP/3.0_DP)

WRITE(stdout,'(/,5x,"Debye temperature = ", f12.3, " K")') debye_t

RETURN
END SUBROUTINE compute_debye_temperature

!-------------------------------------------------------------------------
SUBROUTINE compute_debye_temperature_macro_el(vp, vs, density, nat, &
                                              omega, debye_t)
!-------------------------------------------------------------------------
!
!  This routine receives as input: the compressional and the shear
!  sound velocities in m/s units, the density in kg/m^3, omega in (a.u.)^3, 
!  the number of atoms per cell and gives as output the Debye temperature
!  in K. The Debye temperature is also written in output.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(IN) :: vp, vs
REAL(DP), INTENT(IN) :: density
REAL(DP), INTENT(IN) :: omega
REAL(DP), INTENT(OUT) :: debye_t

REAL(DP) :: fact
REAL(DP) :: average_sound_speed
!
!  first compute the average sound speed
!
average_sound_speed = (( 1.0_DP / vp **3 + 2.0_DP / vs**3)&
                                                  /3.0_DP)**(-1.0_DP/3.0_DP)
!
!  then the debye temperature
!
fact = h_planck_si / k_boltzmann_si / bohr_radius_si
debye_t = average_sound_speed * fact *    &
                       (3.0_DP / pi /4.0_DP * nat / omega)**(1.0_DP/3.0_DP)

!WRITE(stdout,'(/,5x,"Debye temperature = ", f12.3, " K")') debye_t

RETURN
END SUBROUTINE compute_debye_temperature_macro_el

!-------------------------------------------------------------------------
SUBROUTINE debye_vib_energy(debye_t, temp, ntemp, nat, deb_energy)
!-------------------------------------------------------------------------
!
! This routine receives in input the Debye temperature and a set of
! temperatures and computes the Debye vibrational energy in this set
! of temperatures.
! NB: the output energy has not the zero point contribution which
!     is calculated separately by debye_e0
!
IMPLICIT NONE
INTEGER :: ntemp, nat
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: deb_energy(ntemp)

REAL(DP) :: theta_over_t(ntemp)
INTEGER :: i

DO i=1,ntemp
   theta_over_t(i) = debye_t / temp(i)
ENDDO

DO i=1,ntemp
   deb_energy(i) = 3.0_DP * nat * deb_int_ene(theta_over_t(i)) * &
                   k_boltzmann_ry * temp(i)
END DO

RETURN
END SUBROUTINE debye_vib_energy
!
!-------------------------------------------------------------------------
SUBROUTINE debye_energy_0d(debye_t, tt, nat, deb_energy)
!-------------------------------------------------------------------------
!
! This routine receives in input the Debye temperature and a set of
! temperatures and computes the Debye vibrational energy in this set
! of temperatures.
! NB: the output energy has not the zero point contribution which
!     is calculated separately by debye_e0
!
IMPLICIT NONE
INTEGER :: nat
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: tt
REAL(DP), INTENT(INOUT) :: deb_energy

REAL(DP) :: theta_over_t

theta_over_t = debye_t / tt
deb_energy = 3.0_DP * nat * deb_int_ene(theta_over_t) * k_boltzmann_ry * tt

RETURN
END SUBROUTINE debye_energy_0d

!-------------------------------------------------------------------------
SUBROUTINE debye_free_energy(debye_t, temp, ntemp, nat, deb_free_energy)
!-------------------------------------------------------------------------
!
! This routine receives in input the Debye temperature and a set of
! temperatures and computes the Debye vibrational free energy in this set
! of temperatures.
! NB: the output free energy has not the zero point contribution which
!     is calculated separately by debye_e0
!
IMPLICIT NONE
INTEGER :: ntemp, nat
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: deb_free_energy(ntemp)

REAL(DP) :: theta_over_t(ntemp)
INTEGER :: i

DO i=1,ntemp
   theta_over_t(i) = debye_t / temp(i)
ENDDO

DO i=1,ntemp
!   deb_free_energy(i) = 9.0_DP * nat * deb_int_free_ene(theta_over_t(i)) * &
!                   k_boltzmann_ry * temp(i)
!
!  An alternative expression in terms of the energy integral 
!  See M.A. Blanco, E. Francisco, V. Luana, Comp. Phys. Comm. 158, 57 (2004)
!  can be used for debugging purposes.
!
   deb_free_energy(i) = nat * k_boltzmann_ry * temp(i) * &
                        ( - deb_int_ene(theta_over_t(i)) &
           + 3.0_DP * LOG ( 1.0_DP - EXP (- theta_over_t(i)) ) )

END DO

RETURN
END SUBROUTINE debye_free_energy

!-------------------------------------------------------------------------
SUBROUTINE debye_free_energy_0d(debye_t, tt, nat, deb_free_energy)
!-------------------------------------------------------------------------
!
! This routine receives the Debye temperature and one temperature T
! and computes the Debye vibrational free energy at this temperature.
!
! NB: the output free energy has not the zero point contribution which
!     is calculated separately by debye_e0
!
IMPLICIT NONE
INTEGER :: nat
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: tt
REAL(DP), INTENT(INOUT) :: deb_free_energy

REAL(DP) :: theta_over_t

theta_over_t = debye_t / tt
deb_free_energy = nat * k_boltzmann_ry * tt * (- deb_int_ene(theta_over_t)  &
                + 3.0_DP * LOG (1.0_DP - EXP (-theta_over_t)))

RETURN
END SUBROUTINE debye_free_energy_0d

!-------------------------------------------------------------------------
SUBROUTINE debye_entropy(debye_t, temp, ntemp, nat, deb_entropy)
!-------------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER :: ntemp, nat
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: deb_entropy(ntemp)

!REAL(DP) :: deb_energy(ntemp), deb_free_energy(ntemp)
REAL(DP) :: theta_over_t(ntemp)
INTEGER :: i

DO i=1,ntemp
   theta_over_t(i) = debye_t / temp(i)
ENDDO

!  From Comp. Phys. Comm. 158, 57 (2004)

DO i=1,ntemp
   deb_entropy(i) = nat * k_boltzmann_ry * &
           (4.0_DP * deb_int_ene(theta_over_t(i)) &
           -3.0_DP * LOG ( 1.0_DP - EXP (- theta_over_t(i))) )
END DO
!
! For debug purposes the entropy can be calculated also from the
! difference between energy and free energy. Comment the above lines
! and decomment these ones if you want to check an alternative expression.
!
!CALL debye_vib_energy(debye_t, temp, ntemp, nat, deb_energy)
!CALL debye_free_energy(debye_t, temp, ntemp, nat, deb_free_energy)

!DO i=1,ntemp
!   IF (temp(i) > 0.0_DP) THEN
!      deb_entropy(i) = ( deb_energy(i) - deb_free_energy(i) ) / temp(i)
!   ELSE
!      deb_entropy(i) = 0.0_DP
!   ENDIF
!END DO

RETURN
END SUBROUTINE debye_entropy

!-------------------------------------------------------------------------
SUBROUTINE debye_cv (debye_t, temp, ntemp, nat, deb_cv)
!-------------------------------------------------------------------------
!
!  The Debye temperature is in K, the temperature is in K, deb_cv is
!  the heat capacity in Ry / K / cell 
!
IMPLICIT NONE
INTEGER :: ntemp, nat
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(INOUT) :: deb_cv(ntemp)

INTEGER :: i
REAL(DP) :: theta_over_t(ntemp)

DO i=1,ntemp
   theta_over_t(i) = debye_t / temp(i)
END DO

DO i=1,ntemp
!   deb_cv(i) = 9.0_DP * nat * deb_int_cv(theta_over_t(i)) * k_boltzmann_ry
!
!  Alternative expression for the specific heat that use 
!  deb_int_ene without introducing deb_int_cv. Only for debug purposes.
!  From Comp. Phys. Comm. 158, 57 (2004).
!
   deb_cv(i) = nat * k_boltzmann_ry * (                          &
               12.0_DP * deb_int_ene(theta_over_t(i)) - 9.0_DP * &
               theta_over_t(i) / (EXP(theta_over_t(i))-1.0_DP) )
END DO

RETURN
END SUBROUTINE debye_cv
!
!-------------------------------------------------------------------------
SUBROUTINE debye_cv_0d (debye_t, tt, nat, deb_cv)
!-------------------------------------------------------------------------
!
!  The Debye temperature is in K, the temperature is in K, deb_cv is
!  the heat capacity in Ry / K / cell 
!
IMPLICIT NONE
INTEGER :: ntemp, nat
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: tt
REAL(DP), INTENT(INOUT) :: deb_cv

INTEGER :: i
REAL(DP) :: theta_over_t

theta_over_t = debye_t / tt

!deb_cv = 9.0_DP * nat * deb_int_cv(theta_over_t) * k_boltzmann_ry
!
!  Alternative expression for the specific heat that use 
!  deb_int_ene without introducing deb_int_cv. Only for debug purposes.
!  From Comp. Phys. Comm. 158, 57 (2004)
!
   deb_cv = nat * k_boltzmann_ry * (                          &
               12.0_DP * deb_int_ene(theta_over_t) - 9.0_DP * &
               theta_over_t / (EXP(theta_over_t)-1.0_DP) )

RETURN
END SUBROUTINE debye_cv_0d

!-------------------------------------------------------------------------
SUBROUTINE debye_e0(debye_t, nat, deb_e0)
!-------------------------------------------------------------------------
!

IMPLICIT NONE
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(OUT) :: deb_e0
INTEGER :: nat

deb_e0 = 9.0_DP * nat * k_boltzmann_ry * debye_t / 8.0_DP

RETURN
END SUBROUTINE debye_e0

!-------------------------------------------------------------------------
FUNCTION deb_int_cv(y)
!-------------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP) :: y, deb_int_cv

REAL(DP) :: deltax, integral, din, x
INTEGER :: i
!
deltax=y/npt

integral=0.0_DP
DO i=1,npt
   x = deltax * i
   din= x**4 / EXP(x) / ( 1.0_DP - EXP(-x ) )**2
   integral = integral + din
END DO
integral=integral - 0.5_DP * din

deb_int_cv = integral * deltax / y**3

RETURN 
END FUNCTION deb_int_cv


!-------------------------------------------------------------------------
FUNCTION deb_int_ene(y)
!-------------------------------------------------------------------------
!
IMPLICIT NONE
REAL(DP) :: y, deb_int_ene

REAL(DP) :: deltax, integral, din, x, expo
INTEGER :: i
!
IF (y>40.0_DP) THEN
!
!  In this case we compute analytically the integral
!
   deb_int_ene= pi**4/5.0_DP/y**3 - exp(-y) * 3.0_DP *     &
              (6.0_DP/y**3 + 6.0_DP/y**2 +3.0_DP/y+1.0_DP)
   RETURN
ENDIF
deltax=y/npt

integral=0.0_DP
DO i=1,npt
   x = deltax * i
   expo=EXP( -x)
   din= x**3 * expo / ( 1.0_DP - expo )
   integral = integral + din
END DO
integral= integral-0.5_DP*din

deb_int_ene = 3.0_DP * integral * deltax / y**3 

RETURN 
END FUNCTION deb_int_ene

!-------------------------------------------------------------------------
FUNCTION deb_int_free_ene(y)
!-------------------------------------------------------------------------
IMPLICIT NONE
REAL(DP) :: y, deb_int_free_ene

REAL(DP) :: deltax, integral, din, x
INTEGER :: i
!
deltax=y/npt

integral=0.0_DP
DO i=1,npt
   x = deltax * i
   din=x**2 *(LOG( 1.0_DP - EXP( -x ) ) )
   integral = integral + din
END DO
integral=integral-0.5_DP * din

deb_int_free_ene = integral * deltax / y**3 

RETURN 
END FUNCTION deb_int_free_ene

!-------------------------------------------------------------------------
SUBROUTINE compute_debye_temperature_poisson(in_poisson, bulk_modulus, &
                                          density, nat, omega, debye_t)
!-------------------------------------------------------------------------
!
!  This routine receives as input the poisson ratio, the bulk modulus
!  in kbar, the density in kg/m^3, the number of atoms per cell 
!  and the cell volume and gives as output the debye temperature in K.
!  It is simpler to use than compute_debye_temperature because it does not 
!  require the complete elastic constants tensor, but assumes isotropic solid
!  so it can be used for polycrystalline averages or cubic systems.
!  If you do not know the poisson ratio set it to 0.0_DP. In this case
!  the routine uses 0.25.
!  This expression is taken from Comp. Phys. Comm. 158, 57 (2004)
!  
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(IN) :: in_poisson, bulk_modulus, density, omega
REAL(DP), INTENT(OUT) :: debye_t

REAL(DP) :: fact, fun, fun1, poisson

fact = h_planck_si / k_boltzmann_si / bohr_radius_si
fact = fact * (3.0_DP / pi /4.0_DP * nat / omega)**(1.0_DP/3.0_DP) &
       * SQRT(bulk_modulus * 1.D8 /density)
!
!  Default value
!
poisson=in_poisson
IF (in_poisson==0.0_DP) poisson=0.25_DP

fun = 2.0_DP*(1.0_DP + poisson) / 3.0_DP / ( 1.0_DP - 2.0_DP * poisson)
IF (fun>0.0_DP) THEN
   fun = SQRT (fun**3.0_DP)
ELSE
   WRITE(stdout, '(/,5x,"The approximate Debye temperature not available" )')
   RETURN
ENDIF
fun1= (1.0_DP + poisson) / 3.0_DP / ( 1.0_DP - poisson)
IF (fun1>0.0_DP) THEN
   fun1= SQRT (fun1**3.0_DP)
ELSE
   WRITE(stdout, '(/,5x,"The approximate Debye temperature not available" )')
   RETURN
ENDIF
fun = 3.0_DP / ( 2.0_DP * fun + fun1)
fun = fun**(1.0_DP / 3.0_DP)

debye_t = fun * fact

WRITE(stdout, '(/,5x,"The approximate Debye temperature is ",f12.3," K" )') &
     debye_t

RETURN
END SUBROUTINE compute_debye_temperature_poisson

!-------------------------------------------------------------------------
SUBROUTINE write_debye_on_file(temp, ntemp, debye_t, debye_s, filename, &
                               iflag)
!-------------------------------------------------------------------------
!
! This routine creates a file with sound velocities as a function 
! of temperature or of pressure.
! flag=0 temp contains the temperature
! flag=1 temp contains the pressure
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, iflag
REAL(DP), INTENT(IN) :: temp(ntemp), debye_t(ntemp), debye_s(ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_debye, ios
INTEGER :: find_free_unit
CHARACTER(LEN=7) :: label

iu_debye=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_debye, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_debye_on_file','opening debye file',ABS(ios))

IF (iflag==0) THEN
   label='T (K)  '
ELSE
   label='p(kbar)'
ENDIF

IF (meta_ionode) THEN
   WRITE(iu_debye,'("#",2x,"debye_t: computed with isothermal speeds,&
                      & debye_s: computed with adiabatic speeds")')
   WRITE(iu_debye,'("#",2x,a7, 11x, "debye_t (K) ", 6x, "debye_s (K)")') label
   DO itemp=2,ntemp-1
      WRITE(iu_debye,'(e16.8, 8e18.10)') temp(itemp), debye_t(itemp), &
                                                debye_s(itemp)
   ENDDO

   CLOSE(iu_debye)
ENDIF

RETURN
END SUBROUTINE write_debye_on_file

!
! Copyright (C) 2018 Cristiano Malica
!
!-------------------------------------------------------------------------
SUBROUTINE debye_b_factor(debye_t, temp, ntemp, amass, deb_bfact)
!-------------------------------------------------------------------------
!
!  The Debye temperature is in K, the temperature is in K, deb_bfact is
!  the B-factor in Angstrom^2 
!

USE constants, ONLY : amu_si

IMPLICIT NONE
INTEGER :: ntemp, na, i
REAL(DP), INTENT(IN) :: debye_t
REAL(DP), INTENT(IN) :: temp(ntemp), amass(1)
REAL(DP), INTENT(INOUT) :: deb_bfact(ntemp)

REAL(DP) :: theta_over_t(ntemp), fact

fact = 1.D20 * 1.5_DP * h_planck_si**2 / k_boltzmann_si / debye_t / amass(1) &
                                                        / amu_si
DO i=1,ntemp
   theta_over_t(i) = debye_t / temp(i)
END DO

DO i=1,ntemp
   deb_bfact(i) = fact * (1.0_DP + 4.0_DP * deb_int_bfact(theta_over_t(i))) 
END DO

RETURN
END SUBROUTINE debye_b_factor

!-------------------------------------------------------------------------
FUNCTION deb_int_bfact(y)
!-------------------------------------------------------------------------

IMPLICIT NONE
REAL(DP) :: y, deb_int_bfact

REAL(DP) :: deltax, integral, x
INTEGER :: i
!
deltax=y/npt

integral=0.0_DP
DO i=1,npt
   x = deltax * i
   integral = integral + x / ( EXP(x) - 1.0_DP) 
END DO

deb_int_bfact = integral * deltax / y**2

RETURN
END FUNCTION deb_int_bfact



END MODULE debye_module
