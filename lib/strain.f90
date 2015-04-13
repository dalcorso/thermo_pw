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
!   to a crystal lattice. This module uses the following strains
!
!   A   e  0  0     B  e  0  0    B1  e  0  0    B2  0  0  0   
!       0  e  0        0  e  0        0  0  0        0  e  0
!       0  0  e        0  0  0        0  0  e        0  0  e
!
!   C   e  0  0     D  0  0  0    E   0  0  0    
!       0  0  0        0  e  0        0  0  0
!       0  0  0        0  0  0        0  0  e
!
!   F   0  e  e     G  0  e  0    H   0  0  e    I   0  0  0
!       e  0  e        e  0  0        0  0  0        0  0  e
!       e  e  0        0  0  0        e  0  0        0  e  0
!
!   N   e  0  0     O  e  0  0    P   0  0  0
!       0 -e  0        0  0  0        0  e  0
!       0  0  0        0  0 -e        0  0 -e
!
! Not all strains are available for all Bravais lattices. 
!
! The main routine is apply_strain_adv that receives as input:
! strain_code Two characters code with the strain to apply. (the letter above)
! ibrav       The Bravais lattice index of the lattice that is strained
! celldm      The lattice parameters of the lattice that is strained
! epsil       The value of the strain
! 
! As output it gives:
! ibrav_strain   ! The Bravais lattice index of the strained lattice
! celldm_strain  ! The lattice parameters of the strained lattice
! rot            ! the 3x3 rotation matrix that transform the original cartesian
!                ! coordintes in cartesian coordinates of the strained solid 
! epsil_voigt    ! the strain in Voigt notation
!
! The following two matrices express the relation between primitive vectors
! of the strained and unstrained lattice. Note that the strain is not introduced
! here. This is only a rotation between vectors that in different lattices
! are defined in different ways.
!
! aap            ! The matrix that express the original primitive vector in
!                ! terms of the vectors of the strained lattice
! apa            ! The matrix that express the strained primitive vectors in
!                ! terms of the vectors of the unstrained lattice
!
! The ibrav and the a and ap primitive vectors are those defined by the
! routine latgen of QE.
!
! Main limitations: trigonal, monoclinic, and triclinic systems are not
!                   supported
!
!
  USE kinds,     ONLY : DP
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC  apply_strain_adv

CONTAINS
!

SUBROUTINE apply_strain_adv(strain_code, ibrav, celldm, epsil, ibrav_strain, &
                            celldm_strain, epsilon_voigt, rot, aap, apa )

IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav
REAL(DP), INTENT(IN) :: celldm(6), epsil
CHARACTER(LEN=2),INTENT(IN) :: strain_code

INTEGER, INTENT(OUT) :: ibrav_strain
REAL(DP), INTENT(INOUT) :: celldm_strain(6), epsilon_voigt(6), rot(3,3), aap(3,3), &
                           apa(3,3)

REAL(DP), PARAMETER :: sqrt2=SQRT(2.0_DP), sqrt3=SQRT(3.0_DP), sqrt6=SQRT(6.0_DP)
REAL(DP) :: phi, den, aepsilon
!
! some defaults
!
rot(:,:)=0.0_DP
rot(1,1)=1.0_DP
rot(2,2)=1.0_DP
rot(3,3)=1.0_DP
aap(:,:)=rot(:,:)
apa(:,:)=rot(:,:)
epsilon_voigt=0.0_DP
celldm_strain=0.0_DP

IF (ibrav==1) THEN
   IF (strain_code=='A ') THEN
      ibrav_strain=1
      epsilon_voigt(1) = epsil
      epsilon_voigt(2) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * ( 1.0_DP + epsil )
   ELSEIF (strain_code=='C ') THEN
      ibrav_strain=6
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = 1.0_DP + epsil
   ELSEIF (strain_code=='F ') THEN 
      ibrav_strain=6
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
   ELSEIF (strain_code=='C ') THEN
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
   ELSEIF (strain_code=='F ') THEN
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
   ELSEIF (strain_code=='C ') THEN
      ibrav_strain=7
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1)
      celldm_strain(3) = 1.0_DP + epsil
   ELSEIF (strain_code=='F ') THEN
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
      aap(:,:)= 0.0_DP
      aap(1,1) = 1.0_DP
      aap(2,1) = 1.0_DP
      aap(2,2) = 1.0_DP
      aap(2,3) = 1.0_DP
      aap(3,3) = 1.0_DP
      apa(:,:) = 0.0_DP
      apa(1,1) = 1.0_DP
      apa(2,1) =-1.0_DP
      apa(2,2) = 1.0_DP
      apa(2,3) =-1.0_DP
      apa(3,3) = 1.0_DP
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
      aap(1,2)=-1.0_DP
      apa(1,2)=1.0_DP
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
      aap(:,:)=0.0_DP
      aap(1,1)=1.0_DP
      aap(1,2)=1.0_DP
      aap(2,1)=-1.0_DP
      aap(3,3)=1.0_DP
      apa(:,:)=0.0_DP
      apa(1,2)=-1.0_DP
      apa(2,1)=1.0_DP
      apa(2,2)=1.0_DP
      apa(3,3)=1.0_DP
   ELSEIF (strain_code=='B1') THEN
      ibrav_strain=9
      epsilon_voigt(1) = epsil
      epsilon_voigt(3) = epsil
      celldm_strain(1) = celldm(1) * (1.0_DP + epsil)
      celldm_strain(2) = sqrt3 / (1.0_DP + epsil)
      celldm_strain(3) = celldm(3)
      aap(1,2)=-1.0_DP
      apa(1,2)=1.0_DP
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
         aap(2,1)=1.0_DP
         aap(2,2)=-1.0_DP
         aap(2,3)=1.0_DP
         apa(:,:)=aap(:,:)
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
         aap(2,1)=1.0_DP
         aap(2,2)=-1.0_DP
         aap(2,3)=1.0_DP
         apa(:,:)=aap(:,:)
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
         aap=0.0_DP
         aap(1,1)=1.0_DP
         aap(2,3)=1.0_DP
         aap(3,1)=1.0_DP
         aap(3,2)=-1.0_DP
         apa=0.0_DP
         apa(1,1)=1.0_DP
         apa(2,1)=1.0_DP
         apa(2,3)=-1.0_DP
         apa(3,2)=1.0_DP
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
         aap=0.0_DP
         aap(1,2)=-1.0_DP
         aap(1,3)=1.0_DP
         aap(2,1)=-1.0_DP
         aap(2,3)=1.0_DP
         aap(3,2)=-1.0_DP
         apa=0.0_DP
         apa(1,1)=1.0_DP
         apa(1,2)=-1.0_DP
         apa(1,3)=-1.0_DP
         apa(2,3)=-1.0_DP
         apa(3,1)=1.0_DP
         apa(3,3)=-1.0_DP
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
         apa(:,:) = 0.0_DP
         apa(1,1) = 1.0_DP
         apa(1,2) =-1.0_DP
         apa(2,1) = 1.0_DP
         apa(3,3) = 1.0_DP
         aap(:,:) = 0.0_DP
         aap(1,2) = 1.0_DP
         aap(2,1) =-1.0_DP
         aap(2,2) = 1.0_DP
         aap(3,3) = 1.0_DP
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
         apa(:,:) = 0.0_DP
         apa(1,2) = 1.0_DP
         apa(1,3) =-1.0_DP
         apa(2,2) = 1.0_DP
         apa(3,1) = 1.0_DP
         aap(:,:) = 0.0_DP
         aap(1,3) = 1.0_DP
         aap(2,2) = 1.0_DP
         aap(3,1) =-1.0_DP
         aap(3,2) = 1.0_DP
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
         apa(:,:) = 0.0_DP
         apa(1,2) =-1.0_DP
         apa(2,1) = 1.0_DP
         apa(2,2) =-1.0_DP
         apa(3,1) = 1.0_DP
         apa(3,2) =-1.0_DP
         apa(3,3) = 1.0_DP
         aap(:,:) = 0.0_DP
         aap(1,1) = -1.0_DP
         aap(1,2) = 1.0_DP
         aap(2,1) =-1.0_DP
         aap(3,2) =-1.0_DP
         aap(3,3) = 1.0_DP
      ENDIF
   ELSEIF (strain_code=='H ') THEN
      IF (ibrav==8.OR.ibrav==9) THEN
         IF (ibrav==8) THEN
            ibrav_strain=-12
         ELSE
            ibrav_strain=-13
            apa(:,:) = 0.0_DP
            apa(1,2) =-1.0_DP
            apa(2,1) = 1.0_DP
            apa(3,3) = 1.0_DP
            aap(:,:) = 0.0_DP
            aap(1,2) = 1.0_DP
            aap(2,1) =-1.0_DP
            aap(3,3) = 1.0_DP
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
         apa(:,:) = 0.0_DP
         apa(1,2) =-1.0_DP
         apa(2,1) = 1.0_DP
         apa(3,3) = 1.0_DP
         aap(:,:) = 0.0_DP
         aap(1,2) = 1.0_DP
         aap(2,1) =-1.0_DP
         aap(3,3) = 1.0_DP
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
         aap=0.0_DP
         aap(1,1)=-1.0_DP
         aap(1,3)=1.0_DP
         aap(2,1)=-1.0_DP
         aap(3,2)=-1.0_DP
         apa=0.0_DP
         apa(1,2)=-1.0_DP
         apa(2,3)=-1.0_DP
         apa(3,1)=1.0_DP
         apa(3,2)=-1.0_DP
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
         aap=0.0_DP
         aap(1,3)=1.0_DP
         aap(2,1)=1.0_DP
         aap(3,2)=1.0_DP
         apa=0.0_DP
         apa(1,2)=1.0_DP
         apa(2,3)=1.0_DP
         apa(3,1)=1.0_DP
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
         aap=0.0_DP
         aap(1,1)=-1.0_DP
         aap(1,2)= 1.0_DP
         aap(2,3)= 1.0_DP
         aap(3,2)= 1.0_DP
         apa=0.0_DP
         apa(1,1)=-1.0_DP
         apa(1,3)= 1.0_DP
         apa(2,3)= 1.0_DP
         apa(3,2)= 1.0_DP
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

         aap=0.0_DP
         aap(1,1)=-1.0_DP
         aap(1,2)= 1.0_DP
         aap(2,2)=1.0_DP
         aap(2,3)=-1.0_DP
         aap(3,3)=-1.0_DP
         apa=0.0_DP
         apa(1,1)=-1.0_DP
         apa(1,2)=1.0_DP
         apa(1,3)=-1.0_DP
         apa(2,2)=1.0_DP
         apa(2,3)=-1.0_DP
         apa(3,3)=-1.0_DP
      ENDIF
   ELSE
      CALL errore('apply_strain_adv','strain not programmed',ibrav)
   ENDIF
ELSE
   CALL errore('apply_strain_adv','bravais lattice not programmed',1)
ENDIF


RETURN
END SUBROUTINE apply_strain_adv

END MODULE strain_mod
