!
! Copyright (C) 2014-2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM supercell_pos
!
!  This program reads the code of a space group and a set of inequivalent
!  positions and finds all the positions inside a unit cell of the solid.
!  Furthermore it can produce n1 x n2 x n3 cells. It could be useful to
!  study defects. Note however that the defect must be created by
!  the user removing some atoms from the output.
!  For centered lattices it can produce both the coordinates of the
!  atoms in the unit cell of the centered lattice or the coordinates inside the
!  conventional unit cell. The latter cell has 2, 2, 3, 4 times the number of
!  atoms of the primitice unit cell, in the base-centered, body-centered,
!  rhombohedral, and faces centered cells, respectively.
!  Input variables:
!  which_input, integer. If 1 the space group and inequivalent atomic positions
!                        are given, if 2 the Bravais lattice and all atomic
!                        position inside the unit cell.
!  If which_input=1
!     space_group, integer : number of the space group, as on ITA tables
!     origin, unique axis, trigonal or hexagonal, integers: origin choice,
!                           unique axis b or c for monoclinic lattices,
!                           trigonal or or hexagonal coordinates      
!  If which_input=2
!     ibrav, integer    : bravais lattice index
!     units             : 1 alat units, 2 crystal coordinates
!  --- in both cases
!  celldm,  real(8), dimension(6) : dimensions of the unit cell
!  n1, n2, n3, integer : number of cells along a_1, a_2, a_3
!  iconv, integer : 1 for centered cell convert to conventional unit cell
!  nat,  integer : number of atoms
!  'atm', tau(1), tau(2), tau(3) ! atom name and coordinates.
!
USE kinds, ONLY : DP
USE wyckoff,     ONLY : nattot, tautot, ityptot, extfortot, &
                        clean_spacegroup, sup_spacegroup

IMPLICIT NONE
INTEGER :: space_group_code
INTEGER :: ineq_nat           ! the number of inequivalent atoms
REAL(DP), ALLOCATABLE :: ineq_tau(:,:), tau(:,:), tau_new(:,:), rd_for(:,:)
INTEGER :: na, iat, ia, natoms, nat, or, unb, trig
CHARACTER(LEN=3), ALLOCATABLE :: label(:), label_tau(:), label_tau_new(:), &
                                 ineq_label(:)
REAL(DP) :: celldm(6), omega, at(3,3), bg(3,3), a, cg
REAL(DP), PARAMETER :: onet = 1.0_DP / 3.0_DP, twot = 2.0_DP * onet
INTEGER, ALLOCATABLE :: ityp(:), ineq_ityp(:), if_pos(:,:)
LOGICAL :: uniqueb, rhombohedral, found
INTEGER :: n1, n2, n3, nb, ntyp, ibrav, which_input, units
INTEGER :: i1, i2, i3, iconv, nat_new, ibrav_sg, origin_choice

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

WRITE(6,'(5x,"celldm? (default 1.0 0.0 0.0 0.0 0.0 0.0) ")')
READ(5,*) celldm(1), celldm(2), celldm(3), celldm(4), celldm(5), celldm(6)

WRITE(6,'(5x,"n1, n2, n3? (default 1 1 1) ")') 
READ(5,*) n1, n2, n3

WRITE(6,'(5x,"Transform to conventional cell? (default 1=Yes, (0=No)) ")') 
READ(5,*) iconv

WRITE(6,'(5x," Number of atoms and atomic coordinates? ")')
IF (which_input==1) THEN
   READ(5,*) ineq_nat

   ALLOCATE(ineq_tau(3,ineq_nat))
   ALLOCATE(ineq_label(ineq_nat))
   ALLOCATE(label(ineq_nat))
   ALLOCATE(ineq_ityp(ineq_nat))
   ALLOCATE(rd_for(3,ineq_nat))
   ALLOCATE(if_pos(3,ineq_nat))

   DO na=1, ineq_nat
      READ(5,*) label(na), ineq_tau(1,na), ineq_tau(2,na), ineq_tau(3,na)
   ENDDO

   ntyp=0
   DO na=1, ineq_nat
      found=.TRUE.
      DO nb=1, na-1
         IF (label(na) == label(nb)) THEN
            found=.FALSE.
            ineq_ityp(na)=ineq_ityp(nb)
         ENDIF
      ENDDO
      IF (found) THEN
         ntyp=ntyp+1
         ineq_ityp(na)=ntyp
         ineq_label(ntyp)=label(na)
      ENDIF
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

   CALL sup_spacegroup(ineq_tau,ineq_ityp,rd_for,if_pos,&
                                     space_group_code,ineq_nat,uniqueb,&
     rhombohedral,origin_choice,ibrav)

   nat=nattot
   ALLOCATE(tau(3,nat))
   ALLOCATE(ityp(nat))
   ALLOCATE(label_tau(nat))

   tau(:,:)=tautot(:,:)
   ityp(:) = ityptot(:)
   CALL clean_spacegroup()

   DO na=1,nat
      label_tau(na)=ineq_label(ityp(na))
   ENDDO
   DEALLOCATE(ityp)
ELSE
   READ(5,*) nat
   ALLOCATE(tau(3,nat))
   ALLOCATE(label_tau(nat))
   DO na=1, nat
      READ(5,*) label_tau(na), tau(1,na), tau(2,na), tau(3,na)
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
ENDIF

IF (iconv==1) THEN
   SELECT CASE(ibrav)
      CASE(2)
        ibrav=1
        CALL transform_fcc_cubic(tau,nat) 
        nat_new = 4 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
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
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
        label_tau_new(2*nat+1:3*nat)=label_tau(1:nat)
        label_tau_new(3*nat+1:4*nat)=label_tau(1:nat)
      CASE(3)
        ibrav=1
        CALL transform_bcc_cubic(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
      CASE(5)
        ibrav=4
        CALL transform_trig_hex(tau,nat) 
        nat_new = 3 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + twot
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + onet
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + onet
        tau_new(1,2*nat+1:3*nat)=tau(1,1:nat) + onet
        tau_new(2,2*nat+1:3*nat)=tau(2,1:nat) + twot
        tau_new(3,2*nat+1:3*nat)=tau(3,1:nat) + twot
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
        label_tau_new(2*nat+1:3*nat)=label_tau(1:nat)
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
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
      CASE(9)
        ibrav=8
        CALL transform_ofco_orth_c(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) 
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
      CASE(91)
        ibrav=8
        CALL transform_ofco_orth_a(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) 
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
      CASE(10)
        ibrav=8
        CALL transform_fco_orth(tau,nat) 
        nat_new = 4 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
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
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
        label_tau_new(2*nat+1:3*nat)=label_tau(1:nat)
        label_tau_new(3*nat+1:4*nat)=label_tau(1:nat)
      CASE(11)
        ibrav=8
        CALL transform_bco_orth(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
      CASE(13)
        ibrav=12
        CALL transform_ofc_mon_uniq_c(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) 
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) + 0.5_DP
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
      CASE(-13)
        ibrav=-12
        CALL transform_ofc_mon_uniq_b(tau,nat) 
        nat_new = 2 * nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,1:nat)=tau(:,1:nat)
        tau_new(1,nat+1:2*nat)=tau(1,1:nat) + 0.5_DP
        tau_new(2,nat+1:2*nat)=tau(2,1:nat) + 0.5_DP
        tau_new(3,nat+1:2*nat)=tau(3,1:nat) 
        label_tau_new(1:nat)=label_tau(1:nat)
        label_tau_new(nat+1:2*nat)=label_tau(1:nat)
      CASE DEFAULT
        nat_new=nat
        ALLOCATE(tau_new(3,nat_new))
        ALLOCATE(label_tau_new(nat_new))
        tau_new(:,:)=tau(:,:)
        label_tau_new(1:nat)=label_tau(1:nat)
   END SELECT
ELSE
   nat_new=nat
   ALLOCATE(tau_new(3,nat_new))
   ALLOCATE(label_tau_new(nat_new))
   tau_new(:,:)=tau(:,:)
   label_tau_new(1:nat)=label_tau(1:nat)
ENDIF

CALL latgen(ibrav, celldm, at(1,1), at(1,2), at(1,3), omega)

WRITE(6,'(5x, "celldm(1)= ",f12.6)') celldm(1) * n1
IF (celldm(2) > 0.0_DP ) &
   WRITE(6,'(5x, "celldm(2)= ",f12.6)') (celldm(2) * n2) / n1
IF (celldm(3) > 0.0_DP ) &
   WRITE(6,'(5x, "celldm(3)= ",f12.6)') (celldm(3) * n3) / n1
IF (celldm(4) /= 0.0_DP ) &
   WRITE(6,'(5x, "celldm(4)= ",f12.6)') celldm(4) 
IF (celldm(5) /= 0.0_DP ) &
   WRITE(6,'(5x, "celldm(5)= ",f12.6)') celldm(5) 
IF (celldm(6) /= 0.0_DP ) &
   WRITE(6,'(5x, "celldm(6)= ",f12.6)') celldm(6) 

WRITE(6,'(5x,"ibrav= ",i3)') ibrav
WRITE(6,'(5x,"nat= ",i5)') nat_new * n1 * n2 * n3

DO i1=-(n1-1)/2, n1/2
   DO i2=-(n2-1)/2, n2/2
      DO i3=-(n3-1)/2, n3/2
         DO na=1,nat_new
            WRITE(6,'(a,3f18.10)') TRIM(label_tau_new(na)),          &
                 (tau_new(1,na)+i1)/n1, (tau_new(2,na) + i2)/n2, &
                 (tau_new(3,na) + i3)/n3
         ENDDO
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(ineq_tau)
DEALLOCATE(label)
DEALLOCATE(tau)
DEALLOCATE(tau_new)
DEALLOCATE(label_tau)
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
REAL(DP), PARAMETER :: onet = 1.0_DP / 3.0_DP, twot = 2.0_DP * onet
INTEGER, INTENT(IN) :: nat
REAL(DP), INTENT(INOUT) :: tau(3,nat)
REAL(DP) :: tau_new(3,nat)
INTEGER :: na

tau_new=tau
DO na=1,nat
   tau(1,na) = onet*( tau_new(1,na) + tau_new(2,na) - 2.0_DP * tau_new(3,na) )
   tau(2,na) = onet*(-tau_new(1,na) + 2.0_DP * tau_new(2,na) - tau_new(3,na) )
   tau(3,na) = onet*( tau_new(1,na) + tau_new(2,na) + tau_new(3,na) )
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
