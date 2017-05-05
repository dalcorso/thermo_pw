!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE strain_mod
!
!   This module contains the support routines for the application of strain
!   to a crystal lattice. This module defines the following strains
!
!   A   e  0  0     B   e  0  0    B1  e  0  0    B2  0  0  0   
!       0  e  0         0  e  0        0  0  0        0  e  0
!       0  0  e         0  0  0        0  0  e        0  0  e
!
!   C   e  0  0     D   0  0  0    E   0  0  0    
!       0  0  0         0  e  0        0  0  0
!       0  0  0         0  0  0        0  0  e
!
!   F   0  e  e     F1  0 -e  e    F2  0  e -e    F3  0  e  e
!       e  0  e        -e  0  e        e  0  e        e  0 -e
!       e  e  0         e  e  0       -e  e  0        e -e  0
!
!   G   0  e  0     H   0  0  e    I   0  0  0
!       e  0  0         0  0  0        0  0  e
!       0  0  0         e  0  0        0  e  0
!
!   N   e  0  0     O   e  0  0    P   0  0  0
!       0 -e  0         0  0  0        0  e  0
!       0  0  0         0  0 -e        0  0 -e
!
! Not all strains are available for all Bravais lattices. 
!
! The main routine is apply_strain_adv that receives as input:
! strain_code Two characters code with the strain to apply. (the letter above)
! ibrav       The Bravais lattice index of the unstrained lattice 
! celldm      The lattice parameters of the unstrained lattice 
! epsil       The value of the strain
! 
! As output it gives:
! epsil_voigt    ! the strain in Voigt notation
!
! In general the strained lattice has a different ibrav and a different
! celldm with respect to the unstrained Bravais lattice. In general the
! strained Bravais lattice is described with respect to a different set of
! cartesian axis and in some cases the primitive vectors of the strained
! lattice are a linear combination of those obtained by applying epsilon
! to the unstrained vectors and rotating them. Note that the epsil_voigt
! describe the strain tensor in the cartesian axis of the unstrained 
! Bravais lattice.
!
! If the adv flag is set to .TRUE. it gives also
! 
! ibrav_strain   ! The Bravais lattice index of the strained lattice
! celldm_strain  ! The lattice parameters of the strained lattice
! rot            ! the 3x3 rotation matrix that transforms the cartesian
!                ! coordinates of the unstrained solid in the cartesian 
!                ! coordinates of the strained solid 
!
! The ibrav code is the one defined by the routine latgen of QE.
!
! Main limitations: trigonal, monoclinic, and triclinic systems are 
!                   defined only for the standard algorithm
!
!
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC  apply_strain_adv, trans_epsilon, apply_strain, print_strain

CONTAINS
!

SUBROUTINE apply_strain_adv(strain_code, ibrav, celldm, epsil, epsilon_voigt, &
                     ibrav_strain, celldm_strain, rot )

USE matrix_inversion, ONLY : invmat
IMPLICIT NONE
INTEGER, INTENT(IN)  :: ibrav
REAL(DP), INTENT(IN) :: celldm(6), epsil
CHARACTER(LEN=2),INTENT(IN) :: strain_code

INTEGER, INTENT(OUT) :: ibrav_strain
REAL(DP), INTENT(INOUT) :: celldm_strain(6), epsilon_voigt(6), rot(3,3)

REAL(DP), PARAMETER :: sqrt2=SQRT(2.0_DP), sqrt3=SQRT(3.0_DP), &
                                           sqrt6=SQRT(6.0_DP)
REAL(DP) :: phi, den, aepsilon
!
! some defaults
!
rot(:,:)=0.0_DP
rot(1,1)=1.0_DP
rot(2,2)=1.0_DP
rot(3,3)=1.0_DP
epsilon_voigt=0.0_DP
celldm_strain=0.0_DP

IF (ibrav==1) THEN
   IF (strain_code=='A ') THEN
      ibrav_strain=1
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * ( 1.0_DP + epsil )
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain=6
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      celldm_strain(1) = celldm(1) * ( 1.0_DP + epsil )
      celldm_strain(3) = 1.0_DP / ( 1.0_DP + epsil )
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=6
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = 1.0_DP + epsil
   ELSEIF (strain_code=='F ') THEN 
      ibrav_strain=5
      epsilon_voigt(4) = 2.0_DP*epsil
      epsilon_voigt(5) = 2.0_DP*epsil
      epsilon_voigt(6) = 2.0_DP*epsil
      aepsilon = 1.0_DP + 2.0_DP * epsil**2
      celldm_strain(1) = celldm(1) * SQRT(aepsilon)
      celldm_strain(4) = (epsil**2 + 2.0_DP*epsil) / aepsilon
      rot=0.0_DP
      rot(1,1)=1.0_DP/sqrt2
      rot(1,3)=-1.0_DP/sqrt2
      rot(2,1)=-1.0_DP/sqrt6
      rot(2,2)=2.0_DP/sqrt6
      rot(2,3)=-1.0_DP/sqrt6
      rot(3,1)=1.0_DP/sqrt3
      rot(3,2)=1.0_DP/sqrt3
      rot(3,3)=1.0_DP/sqrt3
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF(ibrav==2) THEN
   IF (strain_code=='A ') THEN
      ibrav_strain=2
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain=7
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      celldm_strain(1) = celldm(1) * ( 1.0_DP + epsil ) / sqrt2
      celldm_strain(3) = sqrt2 / ( 1.0_DP + epsil )
      rot=0.0_DP
      rot(1,1)= 1.0_DP / sqrt2
      rot(1,2)=-1.0_DP / sqrt2
      rot(2,1)= 1.0_DP / sqrt2
      rot(2,2)= 1.0_DP / sqrt2
      rot(3,3)= 1.0_DP
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=7
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) / sqrt2
      celldm_strain(3) = sqrt2 * (1.0_DP + epsil)
      rot=0.0_DP
      rot(1,1)= 1.0_DP / sqrt2
      rot(1,2)=-1.0_DP / sqrt2
      rot(2,1)= 1.0_DP / sqrt2
      rot(2,2)= 1.0_DP / sqrt2
      rot(3,3)= 1.0_DP
   ELSEIF (strain_code=='F3') THEN
      ibrav_strain=5
!
!  These sign of epsilon_voigt are necessary for the choice of the
!  primitive vectors of the fcc lattice in latgen. These vectors must
!  becomes the three vector of the trigonal cell.
!
      epsilon_voigt(4) = 2.0_DP * epsil
      epsilon_voigt(5) = -2.0_DP * epsil
      epsilon_voigt(6) = -2.0_DP * epsil
      aepsilon = 6.0_DP*epsil**2+4.0_DP*epsil+2.0_DP
      celldm_strain(1)=celldm(1)*&
                                            SQRT(aepsilon)*0.5_DP
      celldm_strain(4) = (5.0_DP*epsil**2 &
                                           + 6.0_DP*epsil+1.0_DP)/aepsilon
      rot=0.0_DP
      rot(1,1)=1.0_DP/sqrt2
      rot(1,2)=1.0_DP/sqrt2
      rot(2,1)=-1.0_DP/sqrt6
      rot(2,2)=1.0_DP/sqrt6
      rot(2,3)=-2.0_DP/sqrt6
      rot(3,1)=-1.0_DP/sqrt3
      rot(3,2)=1.0_DP/sqrt3
      rot(3,3)=1.0_DP/sqrt3
   ELSEIF (strain_code=='G ') THEN
      ibrav_strain=13
!
      epsilon_voigt(6) = 2.0_DP * epsil

      celldm_strain(1)=celldm(1)*SQRT(1.0_DP+epsil**2)
      celldm_strain(2)=(1.0_DP+epsil)/SQRT(2.0_DP*(1.0_DP+epsil**2))
      celldm_strain(3)=1.0_DP/SQRT(1.0_DP+epsil**2)
      celldm_strain(4) = (1.0_DP+epsil)/SQRT(2.0_DP*(1.0_DP+epsil**2))
      phi=ATAN(epsil)
      rot=0.0_DP
      rot(3,3)=1.0_DP
      rot(1,1)= COS(phi)
      rot(2,1)= -SIN(phi)
      rot(1,2)= SIN(phi)
      rot(2,2)= COS(phi)
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==3) THEN
   IF (strain_code=='A ') THEN
      ibrav_strain=3
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain=7
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = 1.0_DP / (1.0_DP + epsil)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=7
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = 1.0_DP + epsil
   ELSEIF (strain_code=='F3') THEN
!
!  These sign of epsilon_voigt are necessary for the choice of the
!  primitive vectors of the fcc lattice in latgen. These vectors must
!  becomes the three vector of the trigonal cell.
!
      ibrav_strain=5
      epsilon_voigt(4) =-2.0_DP*epsil
      epsilon_voigt(5) = 2.0_DP*epsil
      epsilon_voigt(6) = 2.0_DP*epsil
      aepsilon = 3.0_DP + 4.0_DP * epsil**2 + 4.0_DP * epsil
      celldm_strain(1)=celldm(1)*SQRT(aepsilon)*0.5_DP
      celldm_strain(4)=-(4.0_DP*epsil+1.0_DP)/aepsilon
      rot=0.0_DP
      rot(1,1)=1.0_DP/sqrt2
      rot(1,2)=1.0_DP/sqrt2
      rot(2,1)=-1.0_DP/sqrt6
      rot(2,2)=1.0_DP/sqrt6
      rot(2,3)=-2.0_DP/sqrt6
      rot(3,1)=-1.0_DP/sqrt3
      rot(3,2)=1.0_DP/sqrt3
      rot(3,3)=1.0_DP/sqrt3
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==4) THEN
   IF (strain_code=='C ') THEN
      ibrav_strain=9
      epsilon_voigt(1) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = sqrt3 / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=4
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='H ') THEN
      ibrav_strain=-13
      epsilon_voigt(5) = 2.0_DP * epsil
      celldm_strain(1) = celldm(1) *  SQRT(1.0_DP + epsil**2)
      celldm_strain(2) = sqrt3 / SQRT(1.0_DP + epsil**2)
      celldm_strain(3) = celldm(3)
      celldm_strain(5) = 2.0_DP * epsil / (1.0_DP + epsil**2)
      rot(:,:)=0.0_DP
      phi=-ATAN(epsil)
      rot(2,2)=1.0_DP
      rot(1,1)= COS(phi)
      rot(1,3)= -SIN(phi)
      rot(3,1)= SIN(phi)
      rot(3,3)= COS(phi)
   ELSEIF (strain_code=='B1') THEN
      ibrav_strain=9
      epsilon_voigt(1) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = sqrt3 / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='A ') THEN
      ibrav_strain=4
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==5) THEN
   IF (strain_code=='A ') THEN
      ibrav_strain=5
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      epsilon_voigt(3) = epsil
   ELSEIF (strain_code=='C ') THEN
      epsilon_voigt(1) = epsil
   ELSEIF (strain_code=='E ') THEN
      epsilon_voigt(3) = epsil
   ELSEIF (strain_code=='I ') THEN
      epsilon_voigt(4) = 2.0_DP * epsil
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==6.OR.ibrav==7) THEN
   IF (strain_code=='E ') THEN
      ibrav_strain = ibrav
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='C ') THEN
      IF (ibrav==6) THEN
         ibrav_strain = 8
      ELSE
         ibrav_strain = 11
      ENDIF
      epsilon_voigt(1) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = 1.0_DP / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain = ibrav
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='B1') THEN
      IF (ibrav==6) THEN
         ibrav_strain = 8
      ELSE
         ibrav_strain=11
      ENDIF
      epsilon_voigt(1) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = 1.0_DP / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='G ') THEN
      IF (ibrav==6) THEN
         ibrav_strain=9
      ELSE
         ibrav_strain=10
      ENDIF
      epsilon_voigt(6) = 2.0_DP * epsil
      celldm_strain(1) = celldm(1) * sqrt2 * (1.0_DP - epsil)
      celldm_strain(2) = (1.0_DP + epsil) / (1.0_DP - epsil)
      celldm_strain(3) = celldm(3) / sqrt2 / (1.0_DP - epsil)
      rot(1,1) = 1.0_DP / sqrt2
      rot(1,2) =-1.0_DP / sqrt2
      rot(2,1) = 1.0_DP / sqrt2
      rot(2,2) = 1.0_DP / sqrt2
   ELSEIF (strain_code=='H ') THEN
      IF (ibrav==6) THEN
         ibrav_strain = -12
         epsilon_voigt(5) = 2.0_DP * epsil
         celldm_strain(1) = celldm(1) * (1.0_DP + epsil**2)
         celldm_strain(2) = 1.0_DP / (1.0_DP + epsil**2)
         celldm_strain(3) = celldm(3)
         celldm_strain(5) = 2.0_DP * epsil / (1.0_DP + epsil**2)
         phi=-ATAN(epsil)
         rot(:,:)=0.0_DP
         rot(2,2)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(3,1)= SIN(phi)
         rot(3,3)= COS(phi)
      ELSE
         ibrav_strain=-13
         epsilon_voigt(5) = 2.0_DP * epsil
         den= SQRT((1.0_DP+celldm(3)**2)   &
                    *(1.0_DP + epsil**2) - 4.0_DP * epsil * celldm(3))
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = 1.0_DP / den
         celldm_strain(3) = SQRT(1.0_DP + epsil**2)/den
         celldm_strain(5) = ( (1.0_DP + epsil**2)-      &
                              2.0_DP*epsil*celldm(3) ) /      &
                                 SQRT(1.0_DP + epsil**2) / den
         phi=-ATAN((epsil - celldm(3)) / &
                     (1.0_DP - epsil * celldm(3)))
         rot(:,:)=0.0_DP
         rot(2,2)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(3,1)= SIN(phi)
         rot(3,3)= COS(phi)
      ENDIF
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==8.OR.ibrav==9.OR.ibrav==10.OR.ibrav==11) THEN
   IF (strain_code=='C ') THEN
      ibrav_strain = ibrav
      epsilon_voigt(1) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='D ') THEN
      ibrav_strain = ibrav
      epsilon_voigt(2) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain = ibrav
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain = ibrav
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='B1') THEN
      ibrav_strain = ibrav
      epsilon_voigt(1) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='B2') THEN
      ibrav_strain = ibrav
      epsilon_voigt(2) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='G ') THEN
      IF (ibrav==8) THEN
         ibrav_strain = 12
         epsilon_voigt(6) = 2.0_DP*epsil
         celldm_strain(1) = celldm(1)*SQRT(1.0_DP + epsil**2)
         celldm_strain(2) = celldm(2)
         celldm_strain(3) = celldm(3)/SQRT(1.0_DP + epsil**2)
         celldm_strain(4) = 2.0_DP * epsil/ (1.0_DP + epsil**2)
         phi=-ATAN(epsil)
         rot(:,:)=0.0_DP
         rot(3,3)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,2)= -SIN(phi)
         rot(2,1)= SIN(phi)
         rot(2,2)= COS(phi)
      ELSEIF( ibrav==9) THEN
         ibrav_strain = 12
         epsilon_voigt(6) = 2.0_DP * epsil
         den= SQRT((1.0_DP+celldm(2)**2)   &
                    *(1.0_DP + epsil**2) + 4.0_DP * epsil * celldm(2))
         celldm_strain(1) = celldm(1) * SQRT(1.0_DP + epsil**2)
         celldm_strain(2) = 0.5_DP * den/SQRT(1.0_DP + epsil**2)
         celldm_strain(3) = celldm(3) / &
                                         SQRT(1.0_DP + epsil**2)
         celldm_strain(4) = ( (1.0_DP + epsil**2)+&
                                 2.0_DP*epsil*celldm(2) ) / &
                                 SQRT(1.0_DP + epsil**2) / den
         phi=-ATAN(epsil)
         rot(:,:)=0.0_DP
         rot(3,3)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,2)= -SIN(phi)
         rot(2,1)= SIN(phi)
         rot(2,2)= COS(phi)
      ELSEIF( ibrav==10) THEN
         ibrav_strain = 13
         epsilon_voigt(6) = 2.0_DP * epsil
         den= SQRT((1.0_DP+celldm(2)**2)   &
                    *(1.0_DP + epsil**2) + 4.0_DP * epsil * celldm(2))
         celldm_strain(1) = celldm(1) * SQRT(1.0_DP + epsil**2)
         celldm_strain(2) = 0.5_DP * den/SQRT(1.0_DP + epsil**2)
         celldm_strain(3) = celldm(3) / SQRT(1.0_DP + epsil**2)
         celldm_strain(4) = ( (1.0_DP + epsil**2)+&
                                 2.0_DP*epsil*celldm(2) ) / &
                                 SQRT(1.0_DP + epsil**2) / den

         phi=-ATAN(epsil)
         rot(:,:)=0.0_DP
         rot(3,3)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,2)= -SIN(phi)
         rot(2,1)= SIN(phi)
         rot(2,2)= COS(phi)
      ELSEIF( ibrav==11) THEN
         ibrav_strain=13
         epsilon_voigt(6) = 2.0_DP * epsil
         den= SQRT((1.0_DP+celldm(2)**2)   &
                  *(1.0_DP + epsil**2) - 4.0_DP * epsil * celldm(2))
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = SQRT(1.0_DP + epsil**2)/den
         celldm_strain(3) = celldm(3) / den
         celldm_strain(4) = ( (1.0_DP + epsil**2)-&
                              2.0_DP*epsil*celldm(2) ) / &
                              SQRT(1.0_DP + epsil**2) / den
         phi=-ATAN((epsil - celldm(2)) / &
                  (1.0_DP - epsil * celldm(2)))
         rot(:,:)=0.0_DP
         rot(3,3)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,2)= -SIN(phi)
         rot(2,1)= SIN(phi)
         rot(2,2)= COS(phi)
      ENDIF
   ELSEIF (strain_code=='H ') THEN
      IF (ibrav==8.OR.ibrav==9) THEN
         IF (ibrav==8) THEN
            ibrav_strain=-12
         ELSE
            ibrav_strain=-13
         ENDIF
         epsilon_voigt(5) = epsil * 2.0_DP
         celldm_strain(1)=celldm(1)*SQRT(1.0_DP + epsil**2)
         celldm_strain(2)=celldm(2)/SQRT(1.0_DP + epsil**2)
         celldm_strain(3)=celldm(3)
         celldm_strain(5) = 2.0_DP * epsil/ (1.0_DP + epsil**2)
         phi=-ATAN(epsil)
         rot(:,:)=0.0_DP
         rot(2,2)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(3,1)= SIN(phi)
         rot(3,3)= COS(phi)
      ELSEIF (ibrav==10) THEN
         ibrav_strain=-13
         epsilon_voigt(5) = epsil * 2.0_DP
         den= SQRT((1.0_DP+celldm(3)**2)   &
                 *(1.0_DP + epsil**2) + 4.0_DP * epsil * celldm(3))
         celldm_strain(1) = celldm(1) * SQRT(1.0_DP + epsil**2)
         celldm_strain(2) = celldm(2) / SQRT(1.0_DP + epsil**2)
         celldm_strain(3) = 0.5_DP*den/SQRT(1.0_DP + epsil**2)
         celldm_strain(5) = ( (1.0_DP + epsil**2)+&
                               2.0_DP*epsil*celldm(3) ) / &
                                 SQRT(1.0_DP + epsil**2) / den
         phi=-ATAN(epsil)
         rot(:,:)=0.0_DP
         rot(2,2)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(3,1)= SIN(phi)
         rot(3,3)= COS(phi)
      ELSEIF (ibrav==11) THEN
         ibrav_strain=-13
         epsilon_voigt(5) = 2.0_DP * epsil
         den= SQRT((1.0_DP+celldm(3)**2)   &
                    *(1.0_DP + epsil**2) - 4.0_DP * epsil * celldm(3))
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = celldm(2) / den
         celldm_strain(3) = SQRT(1.0_DP + epsil**2)/den
         celldm_strain(5) = ( (1.0_DP + epsil**2)-&
                                 2.0_DP*epsil*celldm(3) ) / &
                                 SQRT(1.0_DP + epsil**2) / den
         phi=-ATAN((epsil - celldm(3)) / &
                  (1.0_DP - epsil * celldm(3)))
         rot(:,:)=0.0_DP
         rot(2,2)=1.0_DP
         rot(1,1)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(3,1)= SIN(phi)
         rot(3,3)= COS(phi)
      ENDIF
   ELSEIF (strain_code=='I ') THEN
      IF (ibrav==8.OR.ibrav==9) THEN
         IF (ibrav==8) THEN
            ibrav_strain=12
         ELSE
            ibrav_strain=13
         ENDIF
         epsilon_voigt(4) = epsil * 2.0_DP
         celldm_strain(1) = celldm(1) * celldm(2) * SQRT(1.0_DP + epsil**2)
         celldm_strain(2) = celldm(3)/celldm(2)
         celldm_strain(3) = 1.0_DP / celldm(2) / SQRT(1.0_DP + epsil**2)
         celldm_strain(4) = 2.0_DP * epsil/ (1.0_DP + epsil**2)
         phi=-ATAN(epsil)
         rot(:,:)=0.0_DP
         rot(3,1)=1.0_DP
         rot(1,2)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(2,2)= SIN(phi)
         rot(2,3)= COS(phi)
      ELSEIF (ibrav==10) THEN
         ibrav_strain=13
         epsilon_voigt(4) = epsil * 2.0_DP
         den=SQRT((celldm(2)**2+celldm(3)**2)*(1.0_DP + epsil**2)&
                      + 4.0_DP * epsil * celldm(2) * celldm(3))
         celldm_strain(1) = celldm(1) * celldm(2) *  SQRT(1.0_DP + epsil**2)
         celldm_strain(2) =  0.5_DP * den/celldm(2) / SQRT(1.0_DP + epsil**2)
         celldm_strain(3) = 1.0_DP / celldm(2) / SQRT(1.0_DP + epsil**2)
         celldm_strain(4) = ((1.0_DP + epsil**2)*celldm(2) &
                         + 2.0_DP * epsil * celldm(3) ) &
                             / den / SQRT(1.0_DP + epsil**2)
         phi=-ATAN( epsil )
         rot(:,:)=0.0_DP
         rot(3,1)=1.0_DP
         rot(1,2)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(2,2)= SIN(phi)
         rot(2,3)= COS(phi)
      ELSEIF (ibrav==11) THEN
         ibrav_strain=13
         epsilon_voigt(4) = epsil * 2.0_DP
         den=SQRT((celldm(2)**2+celldm(3)**2)*(1.0_DP + epsil**2)&
                      - 4.0_DP * epsil * celldm(2) * celldm(3))
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = celldm(2)/den &
                                            * SQRT(1.0_DP + epsil**2)
         celldm_strain(3) = 1.0_DP / den
         celldm_strain(4) = ((1.0_DP + epsil**2) - 2.0_DP * &
                             epsil * celldm(3)/celldm(2) )/den/&
                                       SQRT(1.0_DP + epsil**2)
         phi=-ATAN( (epsil - celldm(3)/celldm(2)) / &
                    (1.0_DP - epsil * celldm(3)/celldm(2)))
         rot(:,:)=0.0_DP
         rot(3,1)=1.0_DP
         rot(1,2)= COS(phi)
         rot(1,3)= -SIN(phi)
         rot(2,2)= SIN(phi)
         rot(2,3)= COS(phi)
      ENDIF
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==12.OR.ibrav==13.OR.ibrav==-12.OR.ibrav==-13) THEN
   IF (strain_code=='C ') THEN
      epsilon_voigt(1) = epsil 
   ELSEIF (strain_code=='D ') THEN
      epsilon_voigt(2) = epsil 
   ELSEIF (strain_code=='E ') THEN
      epsilon_voigt(3) = epsil 
   ELSEIF (strain_code=='G ') THEN
      epsilon_voigt(6) = 2.0_DP * epsil 
   ELSEIF (strain_code=='H ') THEN
      epsilon_voigt(5) = 2.0_DP * epsil 
   ELSEIF (strain_code=='I ') THEN
      epsilon_voigt(4) = 2.0_DP * epsil 
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ABS(ibrav))
   ENDIF
ELSEIF (ibrav==14) THEN
   IF (strain_code=='C ') THEN
      epsilon_voigt(1) = epsil 
   ELSEIF (strain_code=='D ') THEN
      epsilon_voigt(2) = epsil 
   ELSEIF (strain_code=='E ') THEN
      epsilon_voigt(3) = epsil 
   ELSEIF (strain_code=='G ') THEN
      epsilon_voigt(6) = 2.0_DP * epsil 
   ELSEIF (strain_code=='H ') THEN
      epsilon_voigt(5) = 2.0_DP * epsil 
   ELSEIF (strain_code=='I ') THEN
      epsilon_voigt(4) = 2.0_DP * epsil 
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ABS(ibrav))
   ENDIF
ELSE
   CALL errore('apply_strain_adv','bravais lattice not programmed',1)
ENDIF

RETURN
END SUBROUTINE apply_strain_adv
!
SUBROUTINE trans_epsilon(eps_voigt, eps_tensor, flag)
!
!  This routine transforms the strain tensor from a one dimensional
!  vector with 6 components, to a 3 x 3 symmetric tensor and viceversa
!
!  flag
!   1      6   --> 3 x 3
!  -1    3 x 3 --> 6 
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP) :: eps_voigt(6), eps_tensor(3,3)
INTEGER :: flag
INTEGER :: i, j

IF ( flag == 1 ) THEN
   eps_tensor(1,1) = eps_voigt(1)
   eps_tensor(1,2) = eps_voigt(6) * 0.5_DP
   eps_tensor(1,3) = eps_voigt(5) * 0.5_DP
   eps_tensor(2,2) = eps_voigt(2) 
   eps_tensor(2,3) = eps_voigt(4) * 0.5_DP
   eps_tensor(3,3) = eps_voigt(3) 
   DO i = 1, 3
      DO j = 1, i-1
         eps_tensor(i,j) = eps_tensor(j,i)
      END DO
   END DO
ELSEIF ( flag == -1 ) THEN
  eps_voigt(1) = eps_tensor(1,1)
  eps_voigt(2) = eps_tensor(2,2)
  eps_voigt(3) = eps_tensor(3,3)
  eps_voigt(4) = eps_tensor(2,3) * 2.0_DP
  eps_voigt(5) = eps_tensor(1,3) * 2.0_DP
  eps_voigt(6) = eps_tensor(1,2) * 2.0_DP
ELSE
   CALL errore('trans_epsilon', 'unknown flag', 1)
ENDIF

RETURN
END SUBROUTINE trans_epsilon

SUBROUTINE apply_strain(a, b, epsiln)
!
!  This routine receives as input a vector a and a strain tensor \epsil and
!  gives as output a vector b = a + \epsil a 
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: a(3), epsiln(3,3)
REAL(DP), INTENT(OUT) :: b(3)
INTEGER :: i

b = a 
DO i=1,3
   b(:)=b(:) + epsiln(:,i) * a(i)
ENDDO

RETURN
END SUBROUTINE apply_strain

SUBROUTINE print_strain(strain)
IMPLICIT NONE
REAL(DP), INTENT(IN) :: strain(3,3)
INTEGER :: i, j

WRITE(stdout,'(/,5x,"Applying the following strain")')  
DO i=1,3
   WRITE(stdout,'(5x, "(",2(f12.5,","),f12.5,"  )")') (strain(i,j), j=1,3)
ENDDO
WRITE(stdout,'(/)')

RETURN
END SUBROUTINE print_strain

END MODULE strain_mod
