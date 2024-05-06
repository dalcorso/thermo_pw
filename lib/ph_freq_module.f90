!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ph_freq_module
!
!  This module provide methods to calculate their vibrational contribution 
!  to the free energy, energy, entropy and constant volume specific heat.
!  It defines a type ph_freq_type that contains the phonon frequencies on 
!  a uniform mesh of k points.
!
!
USE kinds, ONLY : DP
USE constants, ONLY : k_boltzmann_ry, ry_to_cmm1
USE io_global, ONLY : stdout
IMPLICIT NONE
SAVE
PRIVATE

REAL(DP), PARAMETER :: kb=k_boltzmann_ry ! Boltzmann constant in Ry/K
REAL(DP), PARAMETER :: kb1=1.0_DP/kb/ry_to_cmm1 ! inverse Boltzmann 
                                                ! constant in cm^{-1}/K

REAL(DP), PARAMETER :: thr_ph=1.D-3   ! a phonon with frequency smaller than
                                      ! this is considered of zero frequency.
!
! NB: we need this threshold to decide which are the three acoustic 
! frequencies. For numerical reason they are never exactly 0 and if 
! introduced in the formulas of this module they give large errors 
! in the average Gruneisen parameters at small temperature.
! All phonon frequencies smaller than thr_ph are removed from the calculation
! A warning is written on output if they are more or less than 3. 
!
REAL(DP), PARAMETER :: thr_taylor=1.D-3   ! When the argument of the functions
                                          ! is smaller than this threshold use
                                          ! a Taylor expansion expression 
                                          ! to avoid the numerical error due
                                          ! to the calculation of 
                                          ! 1.0-(1.0-epsilon) 

TYPE ph_freq_type
   INTEGER :: number_of_points        ! the total number of frequencies
   INTEGER :: nqtot                   ! the total number of q points
   INTEGER :: nq                      ! points reduced by symmetry
   INTEGER :: nq_eff                  ! the points reduced by symmetry and
                                      ! belonging to this processor
   INTEGER :: startq, lastq           ! position of the q in the global q list
   INTEGER :: nq1, nq2, nq3           ! the number of points in the three
                                      ! reciprocal lattice directions
   INTEGER :: nat                     ! the number of atoms
   LOGICAL :: with_eigen              ! if true allocate and save the
                                      ! eigenvectors
   REAL(DP), ALLOCATABLE :: nu(:,:)   ! the frequencies (cm-1)
   COMPLEX(DP), ALLOCATABLE :: displa(:,:,:)  ! the displacements
   REAL(DP), ALLOCATABLE :: wg(:)     ! the weight of each point
END TYPE ph_freq_type


PUBLIC :: ph_freq_type, zero_point_energy_ph, free_energy_ph, vib_energy_ph, &
          vib_entropy_ph, specific_heat_cv_ph, init_ph_freq, destroy_ph_freq, &
          thermal_expansion_ph, read_ph_freq_data, write_ph_freq_data, &
          fecv_ph, debye_waller_factor

CONTAINS

!--------------------------------------------------------------------
SUBROUTINE init_ph_freq(ph_freq, nat, nq1, nq2, nq3, startq, lastq, nq, flag)
!--------------------------------------------------------------------
!
!  If flag is true save also the eigenvectors
!
IMPLICIT NONE

TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq
INTEGER, INTENT(IN) :: nq1, nq2, nq3, nat, startq, lastq, nq
LOGICAL, INTENT(IN) :: flag
INTEGER :: ndiv, nqtot, nq_eff

nqtot = nq1 * nq2 * nq3
ph_freq%nqtot=nqtot
ph_freq%nq_eff=lastq-startq+1
ph_freq%startq=startq
ph_freq%lastq=lastq
ph_freq%nq=nq
ph_freq%number_of_points = 3 * nat * nqtot
ph_freq%nq1 = nq1
ph_freq%nq2 = nq2
ph_freq%nq3 = nq3
ph_freq%nat = nat
ph_freq%with_eigen=flag
ndiv = ph_freq%number_of_points 
nq_eff=ph_freq%nq_eff
ALLOCATE(ph_freq%nu(3*nat, nq_eff))
IF (flag) ALLOCATE(ph_freq%displa(3*nat, 3*nat, nq_eff))
ALLOCATE(ph_freq%wg(nq_eff))

RETURN
END SUBROUTINE init_ph_freq

!--------------------------------------------------------------------
SUBROUTINE destroy_ph_freq(ph_freq)
!--------------------------------------------------------------------
IMPLICIT NONE
TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq

IF (ALLOCATED(ph_freq%nu)) DEALLOCATE(ph_freq%nu)
IF (ALLOCATED(ph_freq%wg)) DEALLOCATE(ph_freq%wg)
IF (ALLOCATED(ph_freq%displa)) DEALLOCATE(ph_freq%displa)

RETURN
END SUBROUTINE destroy_ph_freq

!--------------------------------------------------------------------
SUBROUTINE read_ph_freq_data(ph_freq, filename)
!--------------------------------------------------------------------
!
!  This subroutine reads the ph_freq data from a file. It allocates space
!  for the frequencies and their representation
!
IMPLICIT NONE
TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq
CHARACTER(LEN=256), INTENT(IN) :: filename
INTEGER :: iunit, ios
REAL(DP), ALLOCATABLE :: nu(:), dos(:)
INTEGER :: nat, nq1, nq2, nq3, nq_eff, startq, lastq, nq
INTEGER :: iq, imode, jmode, ndiv, nqtot
INTEGER :: find_free_unit
LOGICAL :: with_eigen

iunit=find_free_unit()
OPEN(file=TRIM(filename), unit=iunit, status='old', &
     form='unformatted', err=100, iostat=ios)
100 CONTINUE
IF (ios /= 0) CALL errore('read_ph_freq_data','opening file',ios)

READ(iunit) nat, nq1, nq2, nq3, nq_eff, startq, lastq, nq
READ(iunit) with_eigen

nqtot = nq1 * nq2 * nq3
ph_freq%nqtot=nqtot
ph_freq%nq_eff=nq_eff
ph_freq%startq=startq
ph_freq%lastq=lastq
ph_freq%nq=nq
ph_freq%number_of_points = 3 * nat * nqtot
ph_freq%nq1 = nq1
ph_freq%nq2 = nq2
ph_freq%nq3 = nq3
ph_freq%nat = nat
ph_freq%with_eigen=with_eigen
ndiv = ph_freq%number_of_points
ALLOCATE(ph_freq%nu(3*nat, nq_eff))
IF (with_eigen) ALLOCATE(ph_freq%displa(3*nat, 3*nat, nq_eff))
ALLOCATE(ph_freq%wg(nq_eff))

READ(iunit, END=20, ERR=10, IOSTAT=ios) ph_freq%wg
    ! nu(i) = frequencies (cm^{-1}) 
READ(iunit, END=20, ERR=10, IOSTAT=ios) ph_freq%nu
IF (with_eigen) &
   READ(iunit, END=20, ERR=10, IOSTAT=ios) ph_freq%displa
20 CONTINUE
   ios=0
10 IF (ios /= 0 ) CALL errore('read_ph_freq_data', 'problem reading phonon &
                                                             &frequencies', 1)
CLOSE(iunit)

RETURN
END SUBROUTINE read_ph_freq_data

!--------------------------------------------------------------------
SUBROUTINE write_ph_freq_data(ph_freq, filename)
!--------------------------------------------------------------------
!
!  This subroutine writes the ph_freq data from a file. It allocates space
!  for the frequencies
!
IMPLICIT NONE
TYPE(ph_freq_type), INTENT(INOUT) :: ph_freq
CHARACTER(LEN=256), INTENT(IN) :: filename
INTEGER :: iunit, ios
INTEGER :: find_free_unit

iunit=find_free_unit()
OPEN(FILE=TRIM(filename), UNIT=iunit, STATUS='unknown', &
     FORM='unformatted', ERR=100, IOSTAT=ios)
100 CONTINUE
IF (ios /= 0) CALL errore('write_ph_freq_data','opening file',ios)

WRITE(iunit) ph_freq%nat, ph_freq%nq1, ph_freq%nq2, &
             ph_freq%nq3, ph_freq%nq_eff, ph_freq%startq, ph_freq%lastq, &
             ph_freq%nq

WRITE(iunit) ph_freq%with_eigen
WRITE(iunit, ERR=10, IOSTAT=ios) ph_freq%wg
   ! nu(i) = frequencies (cm^{-1}) 
WRITE(iunit, ERR=10, IOSTAT=ios) ph_freq%nu 
IF (ph_freq%with_eigen) &
   WRITE(iunit, ERR=10, IOSTAT=ios) ph_freq%displa
10 CONTINUE
IF (ios /= 0 ) CALL errore('write_ph_freq_data', 'problem reading phdos', 1)

CLOSE(iunit)

RETURN
END SUBROUTINE write_ph_freq_data

!--------------------------------------------------------------------
SUBROUTINE zero_point_energy_ph(ph_freq, ener)
!--------------------------------------------------------------------
!
!  This subroutine receives as input a set of phonon frequencies and computes 
!  the zero point energy that corresponds to that frequencies. The output 
!  energy is in Ry.
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(OUT) :: ener
REAL(DP) :: wg, nu
INTEGER :: nq_eff, iq, imode, nat

nq_eff=ph_freq%nq_eff
nat=ph_freq%nat

ener=0.0_DP
DO iq=1, nq_eff
   wg= ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode, iq)
      IF (nu > thr_ph) ener = ener + 0.5_DP * nu * wg
   ENDDO
ENDDO
! result is in cm^{-1}, bring it to Ry
ener = ener / ry_to_cmm1 

RETURN
END SUBROUTINE zero_point_energy_ph

!--------------------------------------------------------------------
SUBROUTINE free_energy_ph(ph_freq, temp, free_ener)
!--------------------------------------------------------------------
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational free energy at that temperature. ener contains
!  only the vibrational contribution WITHOUT the zero point energy. 
!  
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: free_ener

INTEGER :: nq_eff, iq, nat, imode, counter
REAL(DP) :: nu, arg, wg, temp1, onesixth, one24, one120

free_ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
nq_eff=ph_freq%nq_eff
nat=ph_freq%nat
onesixth=1.0_DP / 6.0_DP
one24=1.0_DP /24.0_DP
counter=0
DO iq=1,nq_eff
   wg=ph_freq%wg(iq)
   DO imode=1,3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      IF (arg > thr_taylor) THEN
         free_ener = free_ener +kb*temp*wg*(LOG(1.0_DP-EXP(-arg)))
      ELSEIF (nu > thr_ph) THEN
         free_ener = free_ener+kb*temp*wg*(LOG( arg - arg**2*0.5_DP + &
                                          arg**3 * onesixth - arg**4 * one24 )) 
      ELSE 
         counter=counter+1
      ENDIF
   ENDDO
ENDDO
IF (counter > 3) WRITE(stdout,'(5x,"WARNING: Too many acoustic modes ",i8)') &
                 counter-3

RETURN
END SUBROUTINE free_energy_ph

!--------------------------------------------------------------------
SUBROUTINE vib_energy_ph(ph_freq, temp, ener)
!--------------------------------------------------------------------
!
!  This routine receives as input a set of phonon frequencies and a 
!  temperature and gives as 
!  output the vibrational energy at that temperature. ener contains
!  the energy WITHOUT the zero point energy. 
!  
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: ener

INTEGER :: nq_eff, iq, imode, nat
REAL(DP) :: nu, temp1, arg, wg, onesixth, one24, earg

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
onesixth=1.0_DP / 6.0_DP
one24=1.0_DP /24.0_DP
nq_eff=ph_freq%nq_eff
nat=ph_freq%nat
DO iq=1,nq_eff
   wg=ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      earg = EXP(-arg)
      IF (arg > thr_taylor) THEN
        ener = ener +  wg*nu*earg / ( 1.0_DP - earg ) 
      ELSEIF ( nu > thr_ph ) THEN
         ener = ener +  wg*nu*earg/(arg-arg**2*0.5_DP+arg**3*onesixth &
                                                    - arg**4 * one24)
      ENDIF
   ENDDO
ENDDO
ener = ener / ry_to_cmm1

RETURN
END SUBROUTINE vib_energy_ph

!--------------------------------------------------------------------
SUBROUTINE vib_entropy_ph(ph_freq, temp, entr)
!--------------------------------------------------------------------
!
!  This routine receives as input a phdos and a temperature and gives as 
!  output the vibrational entropy at that temperature. 
!  
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: entr
REAL(DP) :: ener, free_ener

CALL free_energy_ph(ph_freq, temp, free_ener)
CALL vib_energy_ph(ph_freq, temp, ener)

IF (temp > 0.0_DP) THEN
   entr = ( ener - free_ener ) / temp
ELSE
   entr = 0.0_DP
ENDIF

RETURN
END SUBROUTINE vib_entropy_ph

!--------------------------------------------------------------------
SUBROUTINE specific_heat_cv_ph(ph_freq, temp, cv)
!--------------------------------------------------------------------
!
!  This routine receives as input a set of phonon frequencies and a 
!  temperature and gives as output the constant volume specific heat at 
!  that temperature. 
!  The output cv is in Ry / K.
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: cv

INTEGER :: nq_eff, iq, imode, nat
REAL(DP) :: nu, temp1, arg, wg, onesixth, one12, one24, earg, one120, one360

cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
onesixth=1.0_DP / 6.0_DP
one12 = 1.0_DP / 12.0_DP
one24 = 1.0_DP / 24.0_DP
one120= 1.0_DP/ 120.0_DP
one360= 1.0_DP/ 360.0_DP
temp1 = 1.0_DP / temp
nq_eff=ph_freq%nq_eff
nat=ph_freq%nat
DO iq=1,nq_eff
   wg=ph_freq%wg(iq)
   DO imode=1,3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      earg=EXP(-arg)
      IF (arg > thr_taylor ) THEN
          cv = cv + wg * earg * ( arg / ( 1.0_DP - earg  ) ) ** 2 
      ELSEIF (nu > thr_ph) THEN
          cv = cv + wg / (1.0_DP + arg**2 *(one12 + arg**2 * one360))
!          cv = cv + wg * earg *  &
!                  (arg /( arg - arg**2*0.5_DP + arg**3 * onesixth &
!                             - arg**4 * one24 + arg**5 * one120))**2
      ENDIF
   ENDDO
ENDDO
cv = cv * kb

RETURN
END SUBROUTINE specific_heat_cv_ph

!--------------------------------------------------------------------
SUBROUTINE thermal_expansion_ph(ph_freq, ph_grun, temp, betab, cv)
!--------------------------------------------------------------------
!
!  This routine receives as input a set of phonon frequencies
!  on a uniform q mesh, the gruneisen parameters divided by the cell volume
!  on the same mesh and a temperature and gives as output the volume 
!  thermal expansion multiplied by the bulk modulus at that temperature 
!  and the isostrain heat capacity.
!  The output is in Ry / (a.u.)^3 / K. To obtain the volume thermal expansion 
!  divide by the bulk modulus in Ry / (a.u.)^3.
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq, ph_grun
REAL(DP), INTENT(IN)  :: temp
REAL(DP), INTENT(OUT) :: betab, cv

INTEGER :: nq_eff, iq, imode, nat
REAL(DP) :: nu, temp1, arg, gamman, wg, onesixth, one24, earg, aux

betab=0.0_DP
cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
onesixth=1.0_DP/6.0_DP
one24=1.0_DP/24.0_DP
temp1 = 1.0_DP / temp
nq_eff=ph_freq%nq_eff
nat=ph_freq%nat
DO iq=1, nq_eff
   wg=ph_freq%wg(iq)
   DO imode=1, 3*nat
      nu=ph_freq%nu(imode,iq)
      gamman=ph_grun%nu(imode,iq) 
      arg= kb1 * nu * temp1
      earg=EXP(-arg)
      IF (arg > thr_taylor ) THEN
         aux= wg * earg * ( arg / ( 1.0_DP - earg  ) ) ** 2
         cv = cv + aux
         betab = betab +  gamman * aux
      ELSEIF (nu > thr_ph) THEN
         aux= wg * earg * (1.0_DP/( 1.0_DP + arg*(-0.5_DP + arg*(onesixth - &
                        arg*one24))))**2
         cv = cv +  aux 
         betab = betab + gamman * aux 
      ENDIF
   ENDDO
ENDDO
cv = cv * kb
betab = betab * kb

RETURN
END SUBROUTINE thermal_expansion_ph

!--------------------------------------------------------------------
SUBROUTINE fecv_ph(ph_freq, temp, free_ener, ener, cv)
!--------------------------------------------------------------------
!
!  This routine receives as input a set of frequencies 
!  and a temperature and gives as output the vibrational free energy,
!  energy and heat capacity at that temperature. free_ener and ener contains
!  only the vibrational contribution WITHOUT the zero point energy. 
!  
!
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(OUT) :: free_ener, ener, cv

INTEGER :: nq_eff, iq, nat, imode, counter
REAL(DP) :: nu, arg, earg, wg, temp1, onesixth, one24, one12, one360

free_ener=0.0_DP
ener=0.0_DP
cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
nq_eff=ph_freq%nq_eff
nat=ph_freq%nat
onesixth=1.0_DP / 6.0_DP
one12=1.0_DP / 12.0_DP
one360=1.0_DP / 360.0_DP
one24=1.0_DP /24.0_DP
counter=0
DO iq=1,nq_eff
   wg=ph_freq%wg(iq)
   DO imode=1,3*nat
      nu=ph_freq%nu(imode,iq)
      arg= kb1 * nu * temp1
      earg = EXP(-arg)
      IF (arg > thr_taylor) THEN
         free_ener = free_ener +temp*wg*LOG(1.0_DP-earg)
         ener = ener +  wg*nu*earg / ( 1.0_DP - earg ) 
         cv = cv + wg * earg * ( arg / ( 1.0_DP - earg  ) ) ** 2 
      ELSEIF (nu > thr_ph) THEN
         free_ener = free_ener+temp*wg*(LOG( arg* (1.0_DP + arg*(-0.5_DP + &
                                  arg* (onesixth - arg * one24))) )) 
         cv = cv + wg / (1.0_DP + arg**2 *(one12 + arg**2 * one360))
!         cv = cv + wg * earg *  &
!                  (1.0_DP /( 1.0_DP + arg*(-0.5_DP + arg *( onesixth &
!                             - arg* one24 ))))**2 / arg
         ener = ener +  wg*nu*earg/arg/(1.0_DP+arg*(-0.5_DP+arg*(onesixth &
                                                    - arg * one24)))
      ELSE 
         counter=counter+1
      ENDIF
   ENDDO
ENDDO
IF (counter > 3) WRITE(stdout,'(5x,"WARNING: Too many acoustic modes",i8)') &
                 counter-3

free_ener=free_ener*kb
ener = ener / ry_to_cmm1
cv = cv * kb

RETURN
END SUBROUTINE fecv_ph

!--------------------------------------------------------------------
SUBROUTINE debye_waller_factor(ph_freq, temp, b_fact, nat, amass, ntyp, ityp)
!--------------------------------------------------------------------
!
! This subroutine receives as input a set of phonon frequencies and
! respective displacements (eigenvalues and eigenvectors of dynamical
! matrix) and computes the Debye Waller matrix for each atom.
!      
USE constants, ONLY : h_planck_si, c_si, bohr_radius_si, pi, amu_si

IMPLICIT NONE
TYPE(ph_freq_type), INTENT(IN) :: ph_freq
INTEGER, INTENT(IN) :: nat, ntyp, ityp(nat)
REAL(DP), INTENT(IN) :: temp, amass(ntyp)
REAL(DP), INTENT(INOUT) :: b_fact(3,3,nat) 
REAL(DP) :: wg, nu, temp1, arg, expt, tfact, ufact, fact
COMPLEX(DP) :: u1, u2
INTEGER :: nq_eff, iq, ipol, jpol, indi, indj, imode, na

nq_eff=ph_freq%nq_eff
IF (nat/=ph_freq%nat) CALL errore('debye_waller_factor','incompatible nat',1)
temp1 = 1.0_DP / temp

b_fact=0.0_DP
!
!   The input frequencies are supposed to be in cm^-1.
!   c_si * 100 is the speed of light in cm / sec, and its multiplication
!   for the frequency in cm^-1 trasforms it in Hz.
!   amu_si converts the mass containd in amass (suppose to be in amu) in Kg 
!   8 pi^2 comes from a factor 2 in the denominator of the formula,
!   a 2 pi converts h_planck_si in hbar_si, and a 2pi 
!   transforms the frequency into angular frequency.
!
!   NB: the B-factor is defined as 8 pi**2 the mean-square displacement,
!       so we multiply by this factor and write on output the B-factor
!       Switch the comments in the two following lines if you want
!       to output the mean-square displacement.
!       
!
!fact= h_planck_si * 1.D20 / c_si / 800.0_DP / amu_si / pi**2
fact= h_planck_si * 1.D20 / c_si / 100.0_DP / amu_si 
!
!  The output B-factor is in Angstrom^2. Uncomment the following line
!  if you want it in (a.u)^2
!
! fact= fact / 1.D20 / (bohr_radius_si)**2

DO na=1, nat
   DO ipol=1, 3
      indi=3*(na-1)+ipol
      DO jpol=ipol, 3
         indj=3*(na-1)+jpol
         DO iq=1, nq_eff
            wg=ph_freq%wg(iq)
            DO imode=1, 3*nat
               nu=ph_freq%nu(imode, iq)
               u1=ph_freq%displa(indi,imode,iq)
               u2=ph_freq%displa(indj,imode,iq)
               ufact = DREAL(u1*CONJG(u2))
               arg = kb1 * nu * temp1 
               IF (arg > thr_taylor) THEN
                  expt = EXP(-arg)
                  tfact = (1.0_DP+expt)/(1.0_DP-expt)
                  b_fact(ipol,jpol,na)=b_fact(ipol,jpol,na)+ &
                                     wg*ufact*tfact/nu
               ELSEIF (nu > thr_ph) THEN
!                  b_fact(ipol,jpol,na)=b_fact(ipol,jpol,na)+ &
!                      wg*ufact * (2.0_DP - arg + arg**2*0.5_DP-arg**3/6.0_DP) &
!                      /arg/(1.0_DP-arg*0.5_DP+arg**2/6.0_DP) / nu
                  b_fact(ipol,jpol,na)=b_fact(ipol,jpol,na)+ &
                      wg*ufact *(2.0_DP/arg+arg/6.0_DP-arg**3/360.0_DP)/nu
               ENDIF
            ENDDO
         ENDDO
         IF (ipol/=jpol) b_fact(jpol,ipol,na)=b_fact(ipol,jpol,na)
      ENDDO
   ENDDO
!
!   Mass is in amu and the conversion factor is inside fact
!
   b_fact(:,:,na)=b_fact(:,:,na)*fact/amass(ityp(na))
ENDDO   

RETURN
END SUBROUTINE debye_waller_factor

END MODULE ph_freq_module
