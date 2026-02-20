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

!  and some strain that are the composition of two of the previous 
!
!   CG  e  e  0     CH  e  0  e    CI  e  0  0 
!       e  0  0         0  0  0        0  0  e
!       0  0  0         e  0  0        0  e  0
!
!   DG  0  e  0     DH  0  0  e    DI  0  0  0
!       e  e  0         0  e  0        0  e  e 
!       0  0  0         e  0  0        0  e  0
!
!   EG  0  e  0     EH  0  0  e    EI  0  0  0
!       e  0  0         0  0  0        0  0  e
!       0  0  e         e  0  e        0  e  e
!
!   GH  0  e  e     GI  0  e  0    IH  0  0  e
!       e  0  0         e  0  e        0  0  e
!       e  0  0         0  e  0        e  e  0
!
!
! Not all strains are available for all Bravais lattices. 
!
! The main routine is set_strain_adv that receives as input:
! strain_code Two characters code with the strain to apply. (the letter above)
! ibrav       The Bravais lattice index of the unstrained lattice 
! celldm      The lattice parameters of the unstrained lattice 
! epsil       The value of the strain
! 
! As output it gives:
! epsil_voigt    ! the strain in Voigt notation
!
! In general the strained lattice has a different ibrav and a different
! celldm with respect to the unstrained Bravais lattice. The
! strained Bravais lattice is described with respect to a different set of
! cartesian axis and in some cases the primitive vectors of the strained
! lattice are a linear combination of those obtained by applying epsilon
! to the unstrained vectors and rotating them. Note that the epsil_voigt
! describe the strain tensor in the cartesian axis of the unstrained 
! lattice.
!
! In many cases it gives also
! 
! ibrav_strain   ! The Bravais lattice index of the strained lattice
! celldm_strain  ! The lattice parameters of the strained lattice
! rot            ! the 3x3 rotation matrix that transforms the cartesian
!                ! coordinates of the unstrained solid in the cartesian 
!                ! coordinates of the strained solid 
!
! The ibrav code is the one defined by the routine latgen of QE.
!
! The module provides the following routines
!
! set_strain : is the main routine that receives a code of the strain,
!              and sets epsilon_voigt
! set_strain_adv : calls set_strain, but in addition receives the Bravais
!                  lattice index and the celldm of the unperturbed lattice
!                  and sets the ibrav and the celldm of the perturbed lattice.
! trans_epsilon  : converts from Voigt to standard notation and viceversa
! apply_strain   : reveives a 3x3 strain matrix and apply it to a vector
! print_strain   : prints the strain tensor
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

  PUBLIC  set_strain, set_strain_adv, trans_epsilon, apply_strain, print_strain

CONTAINS
!

!-----------------------------------------------------------------------
SUBROUTINE set_strain(strain_code, epsil, epsilon_voigt)
!-----------------------------------------------------------------------
IMPLICIT NONE
REAL(DP), INTENT(IN) :: epsil
REAL(DP), INTENT(OUT) :: epsilon_voigt(6)
CHARACTER(LEN=2),INTENT(IN) :: strain_code

epsilon_voigt=0.0_DP
IF (strain_code=='A ') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(2) = epsil
   epsilon_voigt(3) = epsil
ELSEIF (strain_code=='B ') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(2) = epsil
ELSEIF (strain_code=='B1') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(3) = epsil
ELSEIF (strain_code=='B2') THEN
   epsilon_voigt(2) = epsil
   epsilon_voigt(3) = epsil
ELSEIF (strain_code=='C ') THEN
   epsilon_voigt(1) = epsil
ELSEIF (strain_code=='D ') THEN
   epsilon_voigt(2) = epsil
ELSEIF (strain_code=='E ') THEN
   epsilon_voigt(3) = epsil
ELSEIF (strain_code=='F ') THEN
   epsilon_voigt(4) = 2.0_DP * epsil
   epsilon_voigt(5) = 2.0_DP * epsil
   epsilon_voigt(6) = 2.0_DP * epsil
ELSEIF (strain_code=='F1') THEN
   epsilon_voigt(4) = 2.0_DP * epsil
   epsilon_voigt(5) = 2.0_DP * epsil
   epsilon_voigt(6) =-2.0_DP * epsil
ELSEIF (strain_code=='F2') THEN
   epsilon_voigt(4) = 2.0_DP * epsil
   epsilon_voigt(5) =-2.0_DP * epsil
   epsilon_voigt(6) = 2.0_DP * epsil
ELSEIF (strain_code=='F3') THEN
   epsilon_voigt(4) =-2.0_DP * epsil
   epsilon_voigt(5) = 2.0_DP * epsil
   epsilon_voigt(6) = 2.0_DP * epsil
ELSEIF (strain_code=='G ') THEN
   epsilon_voigt(6) = 2.0_DP * epsil
ELSEIF (strain_code=='H ') THEN
   epsilon_voigt(5) = 2.0_DP * epsil
ELSEIF (strain_code=='I ') THEN
   epsilon_voigt(4) = 2.0_DP * epsil
ELSEIF (strain_code=='N ') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(2) = -epsil
ELSEIF (strain_code=='O ') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(3) = -epsil
ELSEIF (strain_code=='P ') THEN
   epsilon_voigt(2) = epsil
   epsilon_voigt(3) = -epsil
ELSEIF (strain_code=='CG') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(6) = epsil
ELSEIF (strain_code=='CH') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(5) = epsil
ELSEIF (strain_code=='CI') THEN
   epsilon_voigt(1) = epsil
   epsilon_voigt(4) = epsil
ELSEIF (strain_code=='DG') THEN
   epsilon_voigt(2) = epsil
   epsilon_voigt(6) = epsil
ELSEIF (strain_code=='DH') THEN
   epsilon_voigt(2) = epsil
   epsilon_voigt(5) = epsil
ELSEIF (strain_code=='DI') THEN
   epsilon_voigt(2) = epsil
   epsilon_voigt(4) = epsil
ELSEIF (strain_code=='EG') THEN
   epsilon_voigt(3) = epsil
   epsilon_voigt(6) = epsil
ELSEIF (strain_code=='EH') THEN
   epsilon_voigt(3) = epsil
   epsilon_voigt(5) = epsil
ELSEIF (strain_code=='EI') THEN
   epsilon_voigt(3) = epsil
   epsilon_voigt(4) = epsil
ELSEIF (strain_code=='GH') THEN
   epsilon_voigt(4) = epsil
   epsilon_voigt(5) = epsil
ELSEIF (strain_code=='GI') THEN
   epsilon_voigt(4) = epsil
   epsilon_voigt(6) = epsil
ELSEIF (strain_code=='IH') THEN
   epsilon_voigt(5) = epsil
   epsilon_voigt(6) = epsil
ELSE
   WRITE(stdout,'(a2)') strain_code 
   CALL errore('set_strain','strain not programmed',1)
ENDIF

RETURN
END SUBROUTINE set_strain

!-----------------------------------------------------------------------
SUBROUTINE set_strain_adv(strain_code, ibrav, celldm, epsil, epsilon_voigt, &
                     ibrav_strain, celldm_strain, rot, flag )
!-----------------------------------------------------------------------
!
!  This routine sets for every Bravais lattice the strain and
!  the new Bravais lattice and celldm of the strained system.
!  If needed it sets also the rotation matrix from the cartesian
!  axis of the unstrained lattice to those of the strained lattice.
!  
!  TODO:      ibrav        strain
!               5           B1, CI, CG
!               6,7         CG
!               12,13       B B1 B2 CG DG EG HI
!              -12,-13      B B1 B2 CH DH EH GI
!               14          B B1 B2 C  D  E  G
!                           H I  CG CH CI DG DH DI
!                           EG EH EI GH IH IG
!
!  If flag is .TRUE. sets only epsilon_voigt

USE constants,        ONLY : pi
USE rotate,           ONLY : set_rot_xyz
IMPLICIT NONE
INTEGER, INTENT(IN)  :: ibrav
REAL(DP), INTENT(IN) :: celldm(6), epsil
CHARACTER(LEN=2),INTENT(IN) :: strain_code
LOGICAL, INTENT(IN) :: flag

INTEGER, INTENT(OUT) :: ibrav_strain
REAL(DP), INTENT(INOUT) :: celldm_strain(6), epsilon_voigt(6), rot(3,3)

REAL(DP), PARAMETER :: sqrt2=SQRT(2.0_DP), sqrt3=SQRT(3.0_DP), &
                                           sqrt6=SQRT(6.0_DP)
REAL(DP) :: phi, den, den1, den2, sing, sint, cost, aepsilon, rotcr(3,3)
!
! some defaults
!

rot(:,:)=0.0_DP
rot(1,1)=1.0_DP
rot(2,2)=1.0_DP
rot(3,3)=1.0_DP
epsilon_voigt=0.0_DP
celldm_strain=0.0_DP

rotcr(1,1)=1.0_DP/sqrt2
rotcr(1,2)=1.0_DP/sqrt2
rotcr(1,3)=0.0_DP
rotcr(2,1)=-1.0_DP/sqrt6
rotcr(2,2)=1.0_DP/sqrt6
rotcr(2,3)=-2.0_DP/sqrt6
rotcr(3,1)=-1.0_DP/sqrt3
rotcr(3,2)=1.0_DP/sqrt3
rotcr(3,3)=1.0_DP/sqrt3
!
!  set the epsilon_voigt
!
CALL set_strain(strain_code,epsil,epsilon_voigt)
ibrav_strain=0
celldm_strain=0.0_DP
IF (ibrav==14 .OR. flag) RETURN

IF (ibrav==1) THEN
!
!  simple cubic 
!
   IF (strain_code=='A ') THEN
      ibrav_strain=1
      celldm_strain(1) = celldm(1) * ( 1.0_DP + epsil )
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain=6
      celldm_strain(1) = celldm(1) * ( 1.0_DP + epsil )
      celldm_strain(3) = 1.0_DP / ( 1.0_DP + epsil )
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=6
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = 1.0_DP + epsil
   ELSEIF (strain_code=='F ') THEN 
      ibrav_strain=5
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
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF(ibrav==2) THEN
!
!   fcc
!
   IF (strain_code=='A ') THEN
      ibrav_strain=2
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain=7
      celldm_strain(1) = celldm(1) * ( 1.0_DP + epsil ) / sqrt2
      celldm_strain(3) = sqrt2 / ( 1.0_DP + epsil )
      phi = pi / 4.0_DP
      CALL set_rot_xyz(rot, 3, phi)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=7
      celldm_strain(1) = celldm(1) / sqrt2
      celldm_strain(3) = sqrt2 * (1.0_DP + epsil)
      phi = pi / 4.0_DP
      CALL set_rot_xyz(rot, 3, phi)
   ELSEIF (strain_code=='F3') THEN
      ibrav_strain=5
!
!  The choice of F3 instead of F is necessary for the choice of the
!  primitive vectors of the fcc lattice in latgen. These vectors must
!  becomes the three vector of the trigonal cell.
!
      aepsilon = 6.0_DP*epsil**2-4.0_DP*epsil+2.0_DP
      celldm_strain(1)=celldm(1)*SQRT(aepsilon)*0.5_DP
      celldm_strain(4) = (5.0_DP*epsil**2-6.0_DP*epsil+1.0_DP)/aepsilon
      rot=rotcr
!   ELSEIF (strain_code=='G ') THEN
!      ibrav_strain=13
!      celldm_strain(1)=celldm(1)*SQRT(1.0_DP+epsil**2)
!      celldm_strain(2)=(1.0_DP+epsil)/SQRT(2.0_DP*(1.0_DP+epsil**2))
!      celldm_strain(3)=1.0_DP/SQRT(1.0_DP+epsil**2)
!      celldm_strain(4) = (1.0_DP+epsil)/SQRT(2.0_DP*(1.0_DP+epsil**2))
!      phi=ATAN(epsil)
!      CALL set_rot_xyz(rot, 3, phi)
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==3) THEN
!
!  bcc
!
   IF (strain_code=='A ') THEN
      ibrav_strain=3
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain=7
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = 1.0_DP / (1.0_DP + epsil)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=7
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = 1.0_DP + epsil
   ELSEIF (strain_code=='F3') THEN
!
!  The choice of F3 instead of F is necessary for the choice of the
!  primitive vectors of the bcc lattice in latgen. These vectors must
!  become the three vectors of the rhombohedral cell.
!
      ibrav_strain=5
      aepsilon = 3.0_DP + 4.0_DP * epsil**2 + 4.0_DP * epsil
      celldm_strain(1)=celldm(1)*SQRT(aepsilon)*0.5_DP
      celldm_strain(4)=-(4.0_DP*epsil+1.0_DP)/aepsilon
      rot=rotcr
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==4) THEN
!
!   hexagonal
!
   IF (strain_code=='A ') THEN
      ibrav_strain=4
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain=4
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='C ') THEN
      ibrav_strain=9
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = sqrt3 / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=4
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='H ') THEN
      ibrav_strain=-13
      celldm_strain(1) = celldm(1) *  SQRT(1.0_DP + epsil**2)
      celldm_strain(2) = sqrt3 / SQRT(1.0_DP + epsil**2)
      celldm_strain(3) = celldm(3)
      celldm_strain(5) = 2.0_DP * epsil / (1.0_DP + epsil**2)
      phi=-ATAN(epsil)
      CALL set_rot_xyz(rot, 2, phi)
   ELSEIF (strain_code=='B1') THEN
      ibrav_strain=9
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = sqrt3 / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==5) THEN
!
!  rhombohedral lattice
!
   cost=SQRT(1.0_dp + 2.0_dp*celldm(4)) / sqrt3
   sint=SQRT(1.0_dp - celldm(4)) * sqrt2 / sqrt3

   IF (strain_code=='A ') THEN
      ibrav_strain=ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(4) = celldm(4)
   ELSEIF (strain_code=='C ') THEN
      ibrav_strain=-13
      den = SQRT(1.0_DP + 3.0_DP * cost**2)
      celldm_strain(1) = celldm(1) * den
      celldm_strain(2) = sqrt3 * (1.0_DP + epsil) * sint / den
      celldm_strain(3) = 1.0_DP / den
      celldm_strain(5) = (3.0_DP * cost**2 - 1.0_DP) / den
      CALL set_rot_h(rot, sint, cost)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=ibrav
      den= SQRT(sint**2 + (1.0_DP + epsil)**2 * cost**2)
      celldm_strain(1) = celldm(1) * den
      celldm_strain(4) = ((1.0_DP + epsil)**2*cost**2 - 0.5_DP *sint**2) &
                                                   / den / den 
   ELSEIF (strain_code=='I ') THEN
      ibrav_strain=-13
      den= SQRT((1.0_DP+epsil**2)*(1.0_DP + 3.0_DP * cost**2)- &
                              8.0_DP * epsil * sint * cost )
      den1 = SQRT(1.0_DP + epsil**2 + 4.0_DP * epsil * sint * cost )

      celldm_strain(1) = celldm(1) * den
      celldm_strain(2) = sqrt3 * sint / den
      celldm_strain(3) = den1 / den
      celldm_strain(5) = ((1.0_DP + epsil**2)*(3.0_DP* cost**2 - 1.0_DP) &
                        + 2.0_DP * epsil * sint * cost) / den / den1
      CALL set_rot_i(rot, epsil, sint, cost)
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==6.OR.ibrav==7) THEN
!
!   tetragonal or centered tetragonal
!
   IF (strain_code=='A ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) 
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='B1') THEN
      IF (ibrav==6) THEN
         ibrav_strain = 8
      ELSE
         ibrav_strain=11
      ENDIF
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = 1.0_DP / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='C ') THEN
      IF (ibrav==6) THEN
         ibrav_strain = 8
      ELSE
         ibrav_strain = 11
      ENDIF
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = 1.0_DP / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='G ') THEN
      IF (ibrav==6) THEN
         ibrav_strain=9
      ELSE
         ibrav_strain=10
      ENDIF
      celldm_strain(1) = celldm(1) * sqrt2 * (1.0_DP - epsil)
      celldm_strain(2) = (1.0_DP + epsil) / (1.0_DP - epsil)
      celldm_strain(3) = celldm(3) / sqrt2 / (1.0_DP - epsil)
      phi= pi / 4.0_DP
      CALL set_rot_xyz(rot, 3, phi)
   ELSEIF (strain_code=='H ') THEN
      IF (ibrav==6) THEN
         ibrav_strain = -12
         den = SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = 1.0_DP / den
         celldm_strain(3) = celldm(3)
         celldm_strain(5) = 2.0_DP * epsil / den / den
         phi=-ATAN(epsil)
         CALL set_rot_xyz(rot, 2, phi)
      ELSE
!
!    This is option 2, written as before
!
         ibrav_strain=-13
         den=SQRT((1.0_DP+celldm(3)**2)   &
                    *(1.0_DP + epsil**2) - 4.0_DP * epsil * celldm(3))
         den1=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = 1.0_DP / den
         celldm_strain(3) = den1 / den
         celldm_strain(5) = ( (1.0_DP + epsil**2) &
                                 -2.0_DP*epsil*celldm(3) ) / den1 / den
         phi=-ATAN((epsil - celldm(3)) / &
                     (1.0_DP - epsil * celldm(3)))
!
!    This is option 1. The two options differs for the choice of the
!    monoclinic vectors
!

!         ibrav_strain=-13
!         den= SQRT((1.0_DP+celldm(3)**2)   &
!                    *(1.0_DP + epsil**2) + 4.0_DP * epsil * celldm(3))
!         den1=SQRT(1.0_DP + epsil**2)
!         celldm_strain(1) = celldm(1) * den
!         celldm_strain(2) = 1.0_DP / den
!         celldm_strain(3) = celldm(3) * den1 / den
!         celldm_strain(5) = ( (1.0_DP + epsil**2)*celldm(3) &
!                                 +2.0_DP*epsil ) / den1 / den
!         phi=-ATAN((epsil + celldm(3)) / &
!                     (1.0_DP + epsil * celldm(3)))
         CALL set_rot_xyz(rot, 2, phi)
      ENDIF
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==8.OR.ibrav==9.OR.ibrav==10.OR.ibrav==11) THEN
!
!  orthorombic
!
   IF (strain_code=='A ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2)
      celldm_strain(3) = celldm(3) 
   ELSEIF (strain_code=='B ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='B1') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='B2') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='C ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
   ELSEIF (strain_code=='D ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain = ibrav
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2)
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
   ELSEIF (strain_code=='G ') THEN
      IF (ibrav==8) THEN
         ibrav_strain = 12
         den = SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1)*den
         celldm_strain(2) = celldm(2)
         celldm_strain(3) = celldm(3)/den
         celldm_strain(4) = 2.0_DP * epsil/ den / den
         phi=-ATAN(epsil)
         CALL set_rot_xyz(rot, 3, phi)
      ELSEIF( ibrav==9) THEN
         ibrav_strain = 12
         den= SQRT((1.0_DP+celldm(2)**2)   &
                    *(1.0_DP + epsil**2) + 4.0_DP * epsil * celldm(2))
         den1=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den1
         celldm_strain(2) = 0.5_DP * den / den1
         celldm_strain(3) = celldm(3) / den1
         celldm_strain(4) = ( (1.0_DP + epsil**2)+&
                        2.0_DP*epsil*celldm(2) ) / den1 / den
         phi=-ATAN(epsil)
         CALL set_rot_xyz(rot, 3, phi)
      ELSEIF( ibrav==10) THEN
         ibrav_strain = 13
         den= SQRT((1.0_DP+celldm(2)**2)   &
                    *(1.0_DP + epsil**2) + 4.0_DP * epsil * celldm(2))
         den1=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den1
         celldm_strain(2) = 0.5_DP * den/den1
         celldm_strain(3) = celldm(3) / den1
         celldm_strain(4) = ( (1.0_DP + epsil**2)+&
                                 2.0_DP*epsil*celldm(2) ) / den1 / den
         phi=-ATAN(epsil)
         CALL set_rot_xyz(rot, 3, phi)
      ELSEIF( ibrav==11) THEN
         ibrav_strain=13
         den= SQRT((1.0_DP+celldm(2)**2)   &
                  *(1.0_DP + epsil**2) - 4.0_DP * epsil * celldm(2))
         den1= SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = den1 / den
         celldm_strain(3) = celldm(3) / den
         celldm_strain(4) = ( (1.0_DP + epsil**2)-&
                              2.0_DP*epsil*celldm(2) ) / den1 / den
         phi=-ATAN((epsil - celldm(2)) / &
                  (1.0_DP - epsil * celldm(2)))
         CALL set_rot_xyz(rot, 3, phi)
      ENDIF
   ELSEIF (strain_code=='H ') THEN
      IF (ibrav==8.OR.ibrav==9) THEN
         IF (ibrav==8) THEN
            ibrav_strain=-12
         ELSE
            ibrav_strain=-13
         ENDIF
         den=SQRT(1.0_DP + epsil**2)
         celldm_strain(1)=celldm(1) * den
         celldm_strain(2)=celldm(2)/ den
         celldm_strain(3)=celldm(3)
         celldm_strain(5) = 2.0_DP * epsil/ den / den
         phi=-ATAN(epsil)
         CALL set_rot_xyz(rot, 2, phi)
      ELSEIF (ibrav==10) THEN
         ibrav_strain=-13
         den= SQRT((1.0_DP+celldm(3)**2)   &
                 *(1.0_DP + epsil**2) + 4.0_DP * epsil * celldm(3))
         den1=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den1
         celldm_strain(2) = celldm(2) / den1
         celldm_strain(3) = 0.5_DP*den/ den1
         celldm_strain(5) = ( (1.0_DP + epsil**2)+&
                               2.0_DP*epsil*celldm(3) ) / den1 / den
         phi=-ATAN(epsil)
         CALL set_rot_xyz(rot, 2, phi)
      ELSEIF (ibrav==11) THEN
         ibrav_strain=-13
         den= SQRT((1.0_DP+celldm(3)**2)   &
                    *(1.0_DP + epsil**2) - 4.0_DP * epsil * celldm(3))
         den1=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = celldm(2) / den
         celldm_strain(3) = den1/den
         celldm_strain(5) = ( (1.0_DP + epsil**2)-&
                                 2.0_DP*epsil*celldm(3) ) / den1 / den
         phi=-ATAN((epsil - celldm(3)) / &
                  (1.0_DP - epsil * celldm(3)))

         CALL set_rot_xyz(rot, 2, phi)
      ENDIF
   ELSEIF (strain_code=='I ') THEN
      IF (ibrav==8.OR.ibrav==9) THEN
         IF (ibrav==8) THEN
            ibrav_strain=12
         ELSE
            ibrav_strain=13
         ENDIF
         den=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * celldm(2) * den
         celldm_strain(2) = celldm(3)/celldm(2)
         celldm_strain(3) = 1.0_DP / celldm(2) / den
         celldm_strain(4) = 2.0_DP * epsil/ den / den 
         phi=-ATAN(epsil)
         CALL set_rot_f(rot, phi)
      ELSEIF (ibrav==10) THEN
         ibrav_strain=13
         den=SQRT((celldm(2)**2+celldm(3)**2)*(1.0_DP + epsil**2)&
                      + 4.0_DP * epsil * celldm(2) * celldm(3))
         den1=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * celldm(2) * den1
         celldm_strain(2) =  0.5_DP * den/celldm(2) / den1
         celldm_strain(3) = 1.0_DP / celldm(2) / den1
         celldm_strain(4) = ((1.0_DP + epsil**2)*celldm(2) &
                         + 2.0_DP * epsil * celldm(3) ) / den / den1
         phi=-ATAN( epsil )
         CALL set_rot_f(rot, phi)
      ELSEIF (ibrav==11) THEN
         ibrav_strain=13
         den=SQRT((celldm(2)**2+celldm(3)**2)*(1.0_DP + epsil**2)&
                      - 4.0_DP * epsil * celldm(2) * celldm(3))
         den1=SQRT(1.0_DP + epsil**2)
         celldm_strain(1) = celldm(1) * den
         celldm_strain(2) = celldm(2) * den1/den 
         celldm_strain(3) = 1.0_DP / den
         celldm_strain(4) = ((1.0_DP + epsil**2) - 2.0_DP * &
                             epsil * celldm(3)/celldm(2) )/den/den1
         phi=-ATAN( (epsil - celldm(3)/celldm(2)) / &
                    (1.0_DP - epsil * celldm(3)/celldm(2)))
         CALL set_rot_f(rot, phi)
      ENDIF
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==12.OR.ibrav==13) THEN
!
!   monoclinic c-unique
!
   IF (strain_code=='A ') THEN
      ibrav_strain=ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) 
      celldm_strain(3) = celldm(3)
      celldm_strain(4) = celldm(4)
   ELSEIF (strain_code=='C ') THEN
      ibrav_strain=ibrav
      den = SQRT(1.0_DP + (2.0_DP * epsil + epsil**2) * celldm(4)**2)
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) * den / (1.0_DP + epsil) 
      celldm_strain(3) = celldm(3) / (1.0_DP + epsil)
      celldm_strain(4) = celldm(4) * (1.0_DP + epsil) / den 
   ELSEIF (strain_code=='D ') THEN
      ibrav_strain=ibrav
      sing = SQRT ( 1.0_DP - celldm(4)**2 ) 
      den = SQRT(1.0_DP + (2.0_DP * epsil + epsil**2) * sing**2 )
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) * den 
      celldm_strain(3) = celldm(3) 
      celldm_strain(4) = celldm(4) / den 
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=ibrav
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) 
      celldm_strain(3) = celldm(3) * (1.0_DP + epsil)
      celldm_strain(4) = celldm(4)  
   ELSEIF (strain_code=='G ') THEN
      ibrav_strain=ibrav
      sing = SQRT (1.0_DP - celldm(4)**2) 
      den  = SQRT (1.0_DP + epsil**2 + 4.0_DP * epsil * celldm(4) * sing)
      den1 = SQRT (1.0_DP + epsil**2)
      celldm_strain(1) = celldm(1) * den1
      celldm_strain(2) = celldm(2) * den / den1
      celldm_strain(3) = celldm(3) / den1
      celldm_strain(4) = (celldm(4)*(1.0_DP+epsil**2)+2.0_DP*epsil*sing) &
                                     / den / den1
      phi=-ATAN(epsil)
      CALL set_rot_xyz(rot, 3, phi)
   ELSEIF (strain_code=='H ') THEN
      ibrav_strain=14
      sing = SQRT (1.0_DP - celldm(4)**2) 
      den  = SQRT (1.0_DP + epsil**2*celldm(4)**2)
      den1 = SQRT (1.0_DP + epsil**2)
      celldm_strain(1) = celldm(1) * den1
      celldm_strain(2) = celldm(2) * den / den1
      celldm_strain(6) = celldm(4) * den1 / den
      IF (ibrav==12) THEN
         celldm_strain(3) = celldm(3) 
         celldm_strain(4) = 2.0_DP * celldm(4)* epsil / den / den1
         celldm_strain(5) = 2.0_DP * epsil**2 / den1 / den1
      ELSE
         den2 = SQRT((1.0_DP + celldm(3)**2)*(1.0_DP + epsil**2) + 4.0_DP * &
                                       epsil*celldm(3))
         celldm_strain(3) = 0.5_DP * den2 / den1
         celldm_strain(4) = celldm(4)* &
                         (1.0_DP+epsil**2+2.0_DP*epsil*celldm(3))/den/den2
         celldm_strain(5) =(1.0_DP+epsil**2+2.0_DP*epsil*celldm(3))/den1/den2
      ENDIF
      phi=-ATAN(epsil)
      CALL set_rot_xyz(rot, 2, phi)
   ELSEIF (strain_code=='I ') THEN
      ibrav_strain=14
      sing = SQRT (1.0_DP - celldm(4)**2) 
      den  = SQRT (1.0_DP + epsil**2*sing**2)
      den1 = SQRT (1.0_DP + epsil**2)
      celldm_strain(1) = celldm(1) 
      celldm_strain(2) = celldm(2) * den
      celldm_strain(6) = celldm(4) / den
      IF (ibrav==12) THEN
         celldm_strain(3) = celldm(3) * den1 
         celldm_strain(4) = 2.0_DP * sing * epsil / den / den1
         celldm_strain(5) = 0.0_DP 
      ELSEIF (ibrav==13) THEN
         den2 = SQRT(1.0_DP + celldm(3)**2*(1.0_DP + epsil**2))
         celldm_strain(3) = celldm(3) * den2 * 0.5_DP 
         celldm_strain(4) = (celldm(4) + 2*epsil*celldm(3)*sing) / den / den2 
         celldm_strain(5) = 1.0_DP / den2 
      ENDIF
      phi=-ATAN(epsil)
      CALL set_rot_xyz(rot, 1, phi)
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSEIF (ibrav==-12.OR.ibrav==-13) THEN
!
!  Monoclinic b-unique
!
   IF (strain_code=='A ') THEN
      ibrav_strain=ibrav
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) 
      celldm_strain(3) = celldm(3)
      celldm_strain(5) = celldm(5)
   ELSEIF (strain_code=='C ') THEN
      ibrav_strain=ibrav
      den = SQRT(1.0_DP + (2.0_DP * epsil + epsil**2) * celldm(5)**2)
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = celldm(2) / (1.0_DP + epsil) 
      celldm_strain(3) = celldm(3) * den / (1.0_DP + epsil)
      celldm_strain(5) = celldm(5) * (1.0_DP + epsil) / den 
   ELSEIF (strain_code=='D ') THEN
      ibrav_strain=ibrav
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) * (1.0_DP + epsil)
      celldm_strain(3) = celldm(3) 
      celldm_strain(5) = celldm(5)  
   ELSEIF (strain_code=='E ') THEN
      ibrav_strain=ibrav
      sing = SQRT ( 1.0_DP - celldm(5)**2 ) 
      den = SQRT(1.0_DP + (2.0_DP * epsil + epsil**2) * sing**2 )
      celldm_strain(1) = celldm(1)
      celldm_strain(2) = celldm(2) 
      celldm_strain(3) = celldm(3) * den
      celldm_strain(5) = celldm(5) / den 
   ELSEIF (strain_code=='G ') THEN
      ibrav_strain=14
      sing = SQRT (1.0_DP - celldm(5)**2) 
      den  = SQRT (1.0_DP + epsil**2*celldm(5)**2)
      den1 = SQRT (1.0_DP + epsil**2)
      celldm_strain(1) = celldm(1) * den1
      celldm_strain(3) = celldm(3)* den / den1
      celldm_strain(5) = celldm(5)* den1 / den 
      IF (ibrav==-12) THEN
         celldm_strain(2) = celldm(2) 
         celldm_strain(4) = 2.0_DP * celldm(5)* epsil / den / den1
         celldm_strain(6) = 2.0_DP * epsil**2 / den1 / den1
      ELSEIF (ibrav==-13) THEN
         den2 = SQRT((1.0_DP + celldm(2)**2)*(1.0_DP + epsil**2) + 4.0_DP * &
                                       epsil*celldm(2))
         celldm_strain(2) = 0.5_DP * den2 / den1
         celldm_strain(4) = celldm(5)* &
                         (1.0_DP+epsil**2+2.0_DP*epsil*celldm(2))/den/den2
         celldm_strain(6) =(1.0_DP+epsil**2+2.0_DP*epsil*celldm(2))/den1/den2
      ENDIF
      phi=-ATAN(epsil)
      CALL set_rot_xyz(rot, 3, phi)
   ELSEIF (strain_code=='H ') THEN
      ibrav_strain=ibrav
      sing = SQRT (1.0_DP - celldm(5)**2) 
      den  = SQRT (1.0_DP + epsil**2 + 4.0_DP * epsil * celldm(5) * sing)
      den1 = SQRT (1.0_DP + epsil**2)
      celldm_strain(1) = celldm(1) * den1
      celldm_strain(2) = celldm(2) / den1
      celldm_strain(3) = celldm(3) * den / den1
      celldm_strain(5) = (celldm(5)*(1.0_DP+epsil**2)+2.0_DP*epsil*sing) &
                                     / den / den1
      phi=-ATAN(epsil)
      CALL set_rot_xyz(rot, 2, phi)
   ELSEIF (strain_code=='I ') THEN
      ibrav_strain=14
      sing = SQRT (1.0_DP - celldm(5)**2) 
      den  = SQRT (1.0_DP + epsil**2*sing**2)
      den1 = SQRT (1.0_DP + epsil**2)
      celldm_strain(1) = celldm(1) 
      celldm_strain(3) = celldm(3) * den
      celldm_strain(5) = celldm(5) / den 
      IF (ibrav==-12) THEN
         celldm_strain(2) = celldm(2) * den1
         celldm_strain(4) = 2.0_DP * sing * epsil / den / den1
         celldm_strain(6) = 0.0_DP
      ELSEIF (ibrav==-13) THEN
         den2 = SQRT(1.0_DP + celldm(2)**2*(1.0_DP + epsil**2))
         celldm_strain(2) = celldm(2) * den2 * 0.5_DP 
         celldm_strain(4) = (celldm(5) + 2*epsil*celldm(2)*sing) / den / den2 
         celldm_strain(6) = 1.0_DP / den2 
      ENDIF
      phi=-ATAN(epsil)
      CALL set_rot_xyz(rot, 1, phi)
   ELSE
      CALL errore('set_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSE
   CALL errore('set_strain_adv','bravais lattice not programmed',1)
ENDIF

RETURN
END SUBROUTINE set_strain_adv
!
!-----------------------------------------------------------------------
SUBROUTINE trans_epsilon(eps_voigt, eps_tensor, flag)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
SUBROUTINE apply_strain(a, b, epsiln)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
SUBROUTINE print_strain(strain)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
SUBROUTINE set_rot_f(rot, phi)
!-----------------------------------------------------------------------

IMPLICIT NONE
REAL(DP), INTENT(IN) :: phi
REAL(DP), INTENT(OUT) :: rot(3,3)

rot(:,:)=0.0_DP
rot(3,1)=1.0_DP
rot(1,2)= COS(phi)
rot(1,3)= -SIN(phi)
rot(2,2)= SIN(phi)
rot(2,3)= COS(phi)

RETURN
END SUBROUTINE set_rot_f

!-----------------------------------------------------------------------
SUBROUTINE set_rot_h(rot, sint, cost)
!-----------------------------------------------------------------------

IMPLICIT NONE
REAL(DP), INTENT(IN) :: sint, cost
REAL(DP), INTENT(OUT) :: rot(3,3)

REAL(DP) :: den

rot(:,:)=0.0_DP
rot(2,1)=1.0_DP
den= SQRT(1.0_DP + 3.0_DP * cost**2)

rot(1,2)= -sint / den
rot(1,3)= 2.0_DP * cost / den
rot(3,2)= 2.0_DP * cost / den
rot(3,3)= sint / den

RETURN
END SUBROUTINE set_rot_h

!-----------------------------------------------------------------------
SUBROUTINE set_rot_i(rot, epsil, sint, cost)
!-----------------------------------------------------------------------

IMPLICIT NONE
REAL(DP), INTENT(IN) :: epsil, sint, cost
REAL(DP), INTENT(OUT) :: rot(3,3)

REAL(DP) :: den, a, b

rot(:,:)=0.0_DP
rot(2,1)=1.0_DP
den= SQRT((1.0_DP+epsil**2)*(1.0_DP + 3.0_DP * cost**2)- 8.0_DP * epsil * &
                        sint * cost )

a= 2.0_DP * epsil * cost - sint
b= 2.0_DP * cost - epsil * sint

rot(1,2) = a / den
rot(1,3) = b / den
rot(3,2) = b / den
rot(3,3) = -a / den

RETURN
END SUBROUTINE set_rot_i

END MODULE strain_mod

