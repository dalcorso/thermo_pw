!
! Copyright (C) 2014-2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM supercell_pos
!
!  This program reads a space group number and a set of inequivalent
!  positions and finds all the positions inside a unit cell.
!  Furthermore it can write the coordinates of n1 x n2 x n3 cells. 
!  For centered lattices it can write the coordinates of the
!  atoms in the primitive unit cell of the centered lattice or the 
!  coordinates inside the conventional unit cell. The latter has 2, 2, 3, 4 
!  times the number of atoms of the primitice unit cell, in the base-centered, 
!  body-centered, hexagonal, and face-centered cells, respectively.
!  Input variables:
!  which_input, integer. 1: requires the space group number and the 
!                           inequivalent atomic positions
!                        2: requires the Bravais lattice index ibrav as in QE
!                           and all atomic positions inside the unit cell.
!  If which_input=1
!     space_group, integer : space group number, as on ITA tables
!     origin, unique axis, trigonal or hexagonal, integers: 
!                           origin choice (1 or 2),
!                           unique axis b or c for monoclinic lattices (1 or 2),
!                           trigonal or hexagonal coordinates (1 or 2)     
!
!  NB: For the trigonal groups the coordinates can be given in the 
!      rhombohedral or hexagonal basis and the celldm must correspond to
!      the coordinates. So coordinates in rhombohedral axes require celldm(1)
!      and celldm(4), coordinates in hexagonal axes, celldm(1) and celldm(3).
!      On output the code will give the coordinates in the trigonal cell
!      if iconv=0 or in the hexagonal cell if iconv=1. 
!      The hexagonal cell contains 3 times more atoms. 
!
!  If which_input=2
!     ibrav, integer    : Bravais lattice index (ibrav as in QE)
!     units             : 1 alat units, 2 crystal coordinates
!  --- in both cases
!  celldm,  real(8), dimension(6) : dimensions of the unit cell. Same
!                       conventions as QE, see INPUT_PW for more info.
!  n1, n2, n3, integer : number of cells along a_1, a_2, a_3
!  iconv, integer : 1   for centered cells convert atomic positions to 
!                       the conventional unit cell
!  nat,  integer : number of atoms
!  'atm', tau(1), tau(2), tau(3) ! nat lines with atom name and coordinates.
!
USE kinds, ONLY : DP
USE wyckoff,     ONLY : nattot, tautot, ityptot, extfortot, &
                        clean_spacegroup, sup_spacegroup

IMPLICIT NONE
REAL(DP), PARAMETER :: onet = 1.0_DP / 3.0_DP, twot = 2.0_DP * onet
!
!  variables for dealing with space group
!
INTEGER :: space_group_code
INTEGER :: or, unb, trig
INTEGER :: ibrav_sg, origin_choice
LOGICAL :: uniqueb, rhombohedral
INTEGER :: ineq_nat           ! the number of inequivalent atoms
REAL(DP), ALLOCATABLE :: ineq_tau(:,:)  ! the coordinates of inequivalent atoms
INTEGER, ALLOCATABLE :: ineq_ityp(:)    ! the type of each inequivalent atom
INTEGER, ALLOCATABLE :: if_pos(:,:)     ! auxiliary compatibility with QE
REAL(DP), ALLOCATABLE :: rd_for(:,:)    ! auxiliary compatibility with QE
CHARACTER(LEN=3), ALLOCATABLE :: label(:)  ! the name of each inequivalent atom
!
!   variables to specify the geometry
!
INTEGER :: n1, n2, n3 
!
!   the type of atoms
!
INTEGER :: ntyp         ! number of types of atoms  
CHARACTER(LEN=3), ALLOCATABLE :: atm(:)  ! the name of each type
!
!  all atoms in the primitive unit cell
!
INTEGER :: nat                     ! number of atoms
REAL(DP), ALLOCATABLE :: tau(:,:)  ! coordinates
INTEGER, ALLOCATABLE :: ityp(:)    ! the type of each atom
!
!  all atoms in the conventional unit cell
!
INTEGER :: nat_new   ! number of atoms
REAL(DP), ALLOCATABLE :: tau_new(:,:)   ! coordinates of the atom
INTEGER, ALLOCATABLE :: ityp_new(:)     ! type of each atom
!
!  all atoms in the supercell 
!
INTEGER :: all_nat                      ! number of atoms
REAL(DP), ALLOCATABLE :: all_tau(:,:)   ! coordinates of the atom
INTEGER, ALLOCATABLE :: all_ityp(:)     ! type of each atom
!
!   Description of the cell size and shape
!
INTEGER  :: ibrav
REAL(DP) :: celldm(6), omega, at(3,3), bg(3,3)
!
! counters and auxiliary variables
!
INTEGER :: which_input, iconv, units
INTEGER :: na, iat, ia, natoms, nt, nb, i1, i2, i3, iuout
REAL(DP) :: a, cg
LOGICAL :: found
!
!  Part 1 read input variables
!
WRITE(6,'(5x," Space group (1) or standard coordinates (2)? ")')
READ(5,*) which_input

IF (which_input==1) THEN
   WRITE(6,'(5x," Space group number? ")')
   READ(5,*) space_group_code

   WRITE(6,'(5x," Origin, unique axis b, trigonal or hexagonal? (default 1 1 1) ")')
   READ(5,*) or, unb, trig
ELSE
   WRITE(6,'(5x," Bravais lattice code ibrav (as in QE)? ")')
   READ(5,*) ibrav
   WRITE(6,'(5x," Units (alat (1) or crystal coordinates (2)) ? ")')
   READ(5,*) units
ENDIF

WRITE(6,'(5x,"celldm ? (For instance 1.0 0.0 0.0 0.0 0.0 0.0)  ")')
READ(5,*) celldm(1), celldm(2), celldm(3), celldm(4), celldm(5), celldm(6)

WRITE(6,'(5x,"n1, n2, n3? (for instance 1 1 1) ")') 
READ(5,*) n1, n2, n3

WRITE(6,'(5x,"Transform to conventional cell? (1=Yes, 0=No) ")') 
READ(5,*) iconv
!
!  now read atomic positions and compute all atoms in the primitive unit cell 
!  count also the number of types of atoms and the label of each type
!  after this part the code must have
!  nat
!  tau
!  ntyp
!  ityp
!  atm
!
WRITE(6,'(5x," Number of atoms and atomic coordinates? ")')
IF (which_input==1) THEN
   READ(5,*) ineq_nat

   ALLOCATE(ineq_tau(3,ineq_nat))
   ALLOCATE(label(ineq_nat))
   ALLOCATE(ineq_ityp(ineq_nat))
   ALLOCATE(if_pos(3,ineq_nat))
   ALLOCATE(rd_for(3,ineq_nat))

   DO na=1, ineq_nat
      READ(5,*) label(na), ineq_tau(1,na), ineq_tau(2,na), ineq_tau(3,na)
   ENDDO

   ntyp=0
   DO na=1, ineq_nat
      found=.TRUE.
      DO nb=1, na-1
         IF (label(na) == label(nb).AND.found) THEN
            found=.FALSE.
            ineq_ityp(na)=ineq_ityp(nb)
         ENDIF
      ENDDO
      IF (found) THEN
         ntyp=ntyp+1
         ineq_ityp(na)=ntyp
      ENDIF
   ENDDO

   ALLOCATE(atm(ntyp))
   DO nt=1,ntyp
      DO na=1,ineq_nat
         IF (ineq_ityp(na)==nt) THEN
            atm(nt)=label(na)
         ENDIF
      ENDDO
   ENDDO
   
   uniqueb=.FALSE.
   rhombohedral=.FALSE.
   IF (space_group_code > 2 .AND. space_group_code < 16 ) THEN 
      uniqueb=(unb==1)
   ELSEIF ( space_group_code > 142 .AND. space_group_code < 168) THEN
      rhombohedral=(trig==1)
   ENDIF
   origin_choice=or
   rd_for=0.0_DP
   if_pos=0

   CALL sup_spacegroup(ineq_tau, ineq_ityp, rd_for, if_pos, space_group_code, &
        ineq_nat, uniqueb, rhombohedral, origin_choice, ibrav)

   nat=nattot
   ALLOCATE(tau(3,nat))
   ALLOCATE(ityp(nat))

   tau(:,:)=tautot(:,:)
   ityp(:) = ityptot(:)
   CALL clean_spacegroup()

   DEALLOCATE(rd_for)
   DEALLOCATE(if_pos)
   DEALLOCATE(ineq_ityp)
   DEALLOCATE(label)
   DEALLOCATE(ineq_tau)
ELSE
   uniqueb=.FALSE.
   rhombohedral = (ibrav==5)
   READ(5,*) nat
   ALLOCATE(tau(3,nat))
   ALLOCATE(ityp(nat))
   ALLOCATE(label(nat))
   DO na=1, nat
      READ(5,*) label(na), tau(1,na), tau(2,na), tau(3,na)
   ENDDO
!
!   count the number of types and set ityp
!
   ntyp=0
   DO na=1, nat
      found=.TRUE.
      DO nb=1, na-1
         IF (label(na) == label(nb).AND.found) THEN
            found=.FALSE.
            ityp(na)=ityp(nb)
         ENDIF
      ENDDO
      IF (found) THEN
         ntyp=ntyp+1
         ityp(na)=ntyp
      ENDIF
   ENDDO

   ALLOCATE(atm(ntyp))
   DO nt=1,ntyp
      DO na=1,nat
         IF (ityp(na)==nt) THEN
            atm(nt)=label(na)
         ENDIF
      ENDDO
   ENDDO
   !
   !  If coordinates are in alat units, transform to crystal coordinates
   !
   IF (units==1) THEN
      CALL latgen(ibrav, celldm, at(1,1), at(1,2), at(1,3), omega)
      at= at / celldm(1)
      CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
      CALL cryst_to_cart( nat, tau, bg, -1 )
   END IF 
   DEALLOCATE(label)
ENDIF

IF (ibrav==5 .AND. .NOT. rhombohedral) THEN
!
!  the hexagonal celldm have been given, transform to rhombohedral
!
   a = celldm(1) * SQRT(3.0_DP+celldm(3)**2) / 3.0_DP
   cg = ( celldm(3)**2 - 1.5_DP ) / ( celldm(3)**2 + 3.0_DP )
   celldm(1)=a
   celldm(3)=0.0_DP
   celldm(4)=cg
ENDIF
!
!  In this part we set up the conventional cell, or copy the data if
!  the conventional cell is not requested
!  nat_new
!  tau_new
!  ityp_new
!
IF (iconv==1) THEN
   SELECT CASE(ibrav)
      CASE(2)
        ibrav=1
        CALL transform_fcc_cubic(tau,nat) 
        nat_new = 4 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) 
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        tau_new(1,2*nat+1:3*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,2*nat+1:3*nat)=tau(2,1:nat) 
        tau_new(3,2*nat+1:3*nat)=tau(3,1:nat) + 0.5_DP
        tau_new(1,3*nat+1:4*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,3*nat+1:4*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,3*nat+1:4*nat)=tau(3,1:nat) 
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
        ityp_new(2*nat+1:3*nat)=ityp(1:nat)
        ityp_new(3*nat+1:4*nat)=ityp(1:nat)
      CASE(3)
        ibrav=1
        CALL transform_bcc_cubic(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
      CASE(5)
        ibrav=4
        CALL transform_trig_hex(tau,nat) 
        nat_new = 3 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + twot
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + onet
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + onet
        tau_new(1,2*nat+1:3*nat)=tau(1,1:nat) + onet
        tau_new(2,2*nat+1:3*nat)=tau(2,1:nat) + twot
        tau_new(3,2*nat+1:3*nat)=tau(3,1:nat) + twot
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
        ityp_new(2*nat+1:3*nat)=ityp(1:nat)

        a=celldm(1)
        cg=celldm(4)
        celldm(1)=a*SQRT(2.0_DP*(1.0_DP-cg))
        celldm(3)=a*SQRT(3.0_DP*(1.0_DP+2.0_DP*cg)) / celldm(1)
        celldm(4)=0.0_DP
      CASE(7)
        ibrav=6
        CALL transform_ct_tet(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
      CASE(9)
        ibrav=8
        CALL transform_ofco_orth_c(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) 
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
      CASE(91)
        ibrav=8
        CALL transform_ofco_orth_a(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) 
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
      CASE(10)
        ibrav=8
        CALL transform_fco_orth(tau,nat) 
        nat_new = 4 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) 
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        tau_new(1,2*nat+1:3*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,2*nat+1:3*nat)=tau(2,1:nat) 
        tau_new(3,2*nat+1:3*nat)=tau(3,1:nat) + 0.5_DP
        tau_new(1,3*nat+1:4*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,3*nat+1:4*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,3*nat+1:4*nat)=tau(3,1:nat) 
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
        ityp_new(2*nat+1:3*nat)=ityp(1:nat)
        ityp_new(3*nat+1:4*nat)=ityp(1:nat)
      CASE(11)
        ibrav=8
        CALL transform_bco_orth(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
      CASE(13)
        ibrav=12
        CALL transform_ofc_mon_uniq_c(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) 
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
      CASE(-13)
        ibrav=-12
        CALL transform_ofc_mon_uniq_b(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) 
        ityp_new(1:nat)=ityp(1:nat)
        ityp_new(nat+1:2*nat)=ityp(1:nat)
      CASE DEFAULT
        nat_new=nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(ityp_new(nat_new))
        ityp_new(:)=ityp(:)
        tau_new(:,:)=tau(:,:)
   END SELECT
ELSE
   nat_new=nat
   ALLOCATE(tau_new(3,nat_new))
   ALLOCATE(ityp_new(nat_new))
   tau_new(:,:)=tau(:,:)
   ityp_new(:)=ityp(:)
ENDIF
!
!  now setup the supercell. After this part we have
!  all_nat
!  all_tau
!  all_ityp
!  
!
CALL latgen(ibrav, celldm, at(1,1), at(1,2), at(1,3), omega)
all_nat=nat_new * n1 * n2 * n3

WRITE(6,'(5x, "celldm(1)= ",f15.9)') celldm(1) * n1
IF (celldm(2) > 0.0_DP ) &
   WRITE(6,'(5x, "celldm(2)= ",f15.9)') (celldm(2) * n2) / n1
IF (celldm(3) > 0.0_DP ) &
   WRITE(6,'(5x, "celldm(3)= ",f15.9)') (celldm(3) * n3) / n1
IF (celldm(4) /= 0.0_DP ) &
   WRITE(6,'(5x, "celldm(4)= ",f15.9)') celldm(4) 
IF (celldm(5) /= 0.0_DP ) &
   WRITE(6,'(5x, "celldm(5)= ",f15.9)') celldm(5) 
IF (celldm(6) /= 0.0_DP ) &
   WRITE(6,'(5x, "celldm(6)= ",f15.9)') celldm(6) 

WRITE(6,'(5x,"ibrav= ",i3)') ibrav
WRITE(6,'(5x,"nat= ",i5)') all_nat

ALLOCATE(all_tau(3,all_nat))
ALLOCATE(all_ityp(all_nat))

all_nat=0
DO i1=-(n1-1)/2, n1/2
   DO i2=-(n2-1)/2, n2/2
      DO i3=-(n3-1)/2, n3/2
         DO na=1,nat_new
            all_nat=all_nat+1
            all_tau(1,all_nat) = (tau_new(1,na) + i1) / n1 
            all_tau(2,all_nat) = (tau_new(2,na) + i2) / n2
            all_tau(3,all_nat) = (tau_new(3,na) + i3) / n3
            all_ityp(all_nat) = ityp_new(na)
         ENDDO
      ENDDO
   ENDDO
ENDDO

DO na=1,all_nat
   WRITE(6,'(a,3f18.10)') TRIM(atm(all_ityp(na))), all_tau(1,na), &
                                                   all_tau(2,na), &
                                                   all_tau(3,na) 
ENDDO

iuout=35
OPEN(unit=iuout, file='supercell.xsf', status='unknown', &
                                              form='formatted')

at=at/celldm(1)
at(:,1)=at(:,1)*n1
at(:,2)=at(:,2)*n2
at(:,3)=at(:,3)*n3
CALL cryst_to_cart( all_nat, all_tau, at, 1 )

CALL xsf_struct (celldm(1), at, all_nat, all_tau, atm, all_ityp, iuout)

CLOSE(iuout)

DEALLOCATE(atm)
DEALLOCATE(ityp)
DEALLOCATE(tau)
DEALLOCATE(ityp_new)
DEALLOCATE(tau_new)
DEALLOCATE(all_tau)
DEALLOCATE(all_ityp)

END PROGRAM supercell_pos

SUBROUTINE transform_fcc_cubic(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) = - 0.5_DP * ( tau_new(1,na) + tau_new(3,na) )
   tau(2,na) = 0.5_DP * ( tau_new(2,na) + tau_new(3,na) )
   tau(3,na) = 0.5_DP * ( tau_new(1,na) + tau_new(2,na) )
END DO

RETURN
END SUBROUTINE transform_fcc_cubic

SUBROUTINE transform_bcc_cubic(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) = - 0.5_DP * ( tau_new(1,na) - tau_new(2,na) - tau_new(3,na) )
   tau(2,na) =   0.5_DP * ( tau_new(1,na) + tau_new(2,na) - tau_new(3,na) )
   tau(3,na) =   0.5_DP * ( tau_new(1,na) + tau_new(2,na) + tau_new(3,na) )
END DO

RETURN
END SUBROUTINE transform_bcc_cubic

SUBROUTINE transform_trig_hex(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), PARAMETER :: onet = 1.0_DP / 3.0_DP
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) = onet*( tau_new(1,na) + tau_new(2,na) - 2.0_DP * tau_new(3,na) )
   tau(2,na) = onet*(-tau_new(1,na) + 2.0_DP * tau_new(2,na) - tau_new(3,na) )
   tau(3,na) =-onet*( tau_new(1,na) + tau_new(2,na) + tau_new(3,na) )
END DO

RETURN
END SUBROUTINE transform_trig_hex

SUBROUTINE transform_ct_tet(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) =  0.5_DP * (   tau_new(1,na) + tau_new(2,na) - tau_new(3,na) )
   tau(2,na) =  0.5_DP * ( - tau_new(1,na) + tau_new(2,na) - tau_new(3,na) )
   tau(3,na) =  0.5_DP * (   tau_new(1,na) + tau_new(2,na) + tau_new(3,na) )
END DO

RETURN
END SUBROUTINE transform_ct_tet

SUBROUTINE transform_ofco_orth_c(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) =  0.5_DP * ( tau_new(1,na) - tau_new(2,na) )
   tau(2,na) =  0.5_DP * ( tau_new(1,na) + tau_new(2,na) )
   tau(3,na) =  tau_new(3,na) 
END DO

RETURN
END SUBROUTINE transform_ofco_orth_c

SUBROUTINE transform_ofco_orth_a(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) =  tau_new(1,na) 
   tau(2,na) =  0.5_DP * ( tau_new(2,na) - tau_new(3,na) )
   tau(3,na) =  0.5_DP * ( tau_new(2,na) + tau_new(3,na) )
END DO

RETURN
END SUBROUTINE transform_ofco_orth_a

SUBROUTINE transform_fco_orth(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) =  0.5_DP * ( tau_new(1,na) + tau_new(2,na) )
   tau(2,na) =  0.5_DP * ( tau_new(2,na) + tau_new(3,na) )
   tau(3,na) =  0.5_DP * ( tau_new(1,na) + tau_new(3,na) )
END DO

RETURN
END SUBROUTINE transform_fco_orth

SUBROUTINE transform_bco_orth(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) = - 0.5_DP * ( tau_new(1,na) - tau_new(2,na) - tau_new(3,na) )
   tau(2,na) =   0.5_DP * ( tau_new(1,na) + tau_new(2,na) + tau_new(3,na) )
   tau(3,na) =   0.5_DP * ( tau_new(1,na) + tau_new(2,na) + tau_new(3,na) )
END DO

RETURN
END SUBROUTINE transform_bco_orth

SUBROUTINE transform_ofc_mon_uniq_c(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) =  0.5_DP * ( tau_new(1,na) + tau_new(3,na) )
   tau(2,na) =  tau_new(2,na)  
   tau(3,na) =  0.5_DP * ( - tau_new(1,na) + tau_new(3,na) )
END DO

RETURN
END SUBROUTINE transform_ofc_mon_uniq_c

SUBROUTINE transform_ofc_mon_uniq_b(tau,nat) 
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) =  0.5_DP * ( tau_new(1,na) - tau_new(2,na) )
   tau(2,na) =  0.5_DP * ( tau_new(1,na) + tau_new(2,na) ) 
   tau(3,na) =  tau_new(3,na)
END DO

RETURN
END SUBROUTINE transform_ofc_mon_uniq_b
