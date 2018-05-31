!
! Copyright (C) 2014-2016 Andrea Dal Corso 
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
!  --- in both cases (if ibrav/=0)
!  celldm,  real(8), dimension(6) : dimensions of the unit cell. Same
!                       conventions as QE, see INPUT_PW for more info.
!  --- if ibrav=0 
!      at units        real        ! a real number. All at are multiplied by
!                                  ! this number
!      at(1,1), at(2,1), at(3,1)   ! first primitive vector
!      at(1,2), at(2,2), at(3,2)   ! second primitive vector
!      at(1,3), at(2,3), at(3,3)   ! third primitive vector
!  --- in all cases
!  iconv, integer : 1   for centered cells convert atomic positions to 
!                       the conventional unit cell
!  format                          ! 1 to give n1, n2, n3
!                                  ! 2 to give the matrix m_ij
!                                  ! 3 to give vectors A_1, A_2, A_3 in 
!                                  !   cartesian coordinates
!  case format=1
!  n1, n2, n3                      ! 
!                                  ! The supercell is defined by
!                                    A_1 = n1 a_1
!                                    A_2 = n2 a_2
!                                    A_3 = n3 a_3
!
!  case format=2
!  m11, m12, m13, integer : 
!  m21, m22, m23, integer : 
!  m31, m32, m33, integer : 
!                               The supercell is defined by the 
!                               three primitive vectors
!
!                               A_1 = m11 a_1 + m12 a_2 + m13 a_3
!                               A_2 = m21 a_1 + m22 a_2 + m23 a_3
!                               A_3 = m31 a_1 + m32 a_2 + m33 a_3
!
!  case format=3
!  A_1x, A_1y, A_1z
!  A_2x, A_2y, A_2z
!  A_3x, A_3y, A_3z
!
!                              Note that a_1, a_2, a_3 are those of the
!                              input lattice if iconv=0, or those of
!                              the conventional cell if iconv=1. 
!                              In case format=3 A_1, A_2, and A_3 given
!                              in input must have integer components
!                              on a_1, a_2 and a_3.
!
!  icenter, integer : 0 no change to atomic coordinates
!                   : 1 output atomic coordinates between -0.5 and 0.5 in
!                       crystal coordinates
!                   : 2 output atomic coordinates between 0 and 1 in
!                       crystal coordinates
!  nat,  integer : number of atoms
!  'atm', tau(1), tau(2), tau(3) ! nat lines with atom name and coordinates.
!
USE kinds,       ONLY : DP
USE wyckoff,     ONLY : nattot, tautot, ityptot, extfortot, &
                        clean_spacegroup, sup_spacegroup
USE wy_pos,      ONLY : wypos
USE parser,      ONLY : read_line, get_field, field_count
USE atomic_pos,  ONLY : find_ityp
USE lattices,    ONLY : find_ibrav_code
USE wrappers,    ONLY : feval_infix
USE io_global,   ONLY : ionode, stdout
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end

IMPLICIT NONE
INTEGER, PARAMETER :: ntypx = 10
REAL(DP), PARAMETER :: onet = 1.0_DP / 3.0_DP, twot = 2.0_DP * onet
REAL(DP), PARAMETER :: eps=1.D-9
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
INTEGER :: mm(3,3)   !  m_ij  see above for the definition of the supercell
!
!   the type of atoms
!
INTEGER :: ntyp         ! number of types of atoms  
CHARACTER(LEN=3)   :: atm(ntypx)  ! the name of each type
!
!  all atoms in the primitive unit cell
!
INTEGER :: nat                     ! number of atoms
REAL(DP), ALLOCATABLE :: tau(:,:)  ! coordinates
REAL(DP), ALLOCATABLE :: stau(:,:) ! save the coordinates
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
REAL(DP) :: celldm(6), omega, at(3,3), bg(3,3), ur(3,3), global_s(3,3)
!
!   Description of the supercell size and shape
!
INTEGER :: ibravs
REAL(DP) :: celldms(6), ats(3,3), bgs(3,3)
!
! counters and auxiliary variables
!
INTEGER :: which_input, iconv, icenter, units
INTEGER :: na, iat, ia, natoms, nt, nb, i1, i2, i3, iuout, nfield, idx, k, &
           ivec, jvec, ipol, jpol, ierr, code_group_ext, copies, idet, iv, &
           n1, n2, n3, iformat
INTEGER :: find_free_unit
REAL(DP) :: a, cg, inp(3), prod1, prod2, prod3, rmu(3), radius, rmod, v(3,7), &
            rmm(3,3), check
LOGICAL :: found
CHARACTER (LEN=256) :: input_line, field_str, wp
CHARACTER(LEN=9) :: code='SUPERCELL'
!
!  Part 1 read input variables
!
CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'(5x," Space group (1) or standard coordinates (2)? ")')
READ(5,*) which_input
WRITE(stdout,'(i5)') which_input

IF (which_input==1) THEN
   units=2
   WRITE(stdout,'(5x," Space group number? ")')
   READ(5,*) space_group_code
   WRITE(stdout,'(i5)') space_group_code

   WRITE(stdout,'(5x," Origin, unique axis b, trigonal or hexagonal? (default 1 1 1) ")')
   READ(5,*) or, unb, trig

   WRITE(stdout,'(3i5)') or, unb, trig
   uniqueb=.FALSE.
   rhombohedral=.FALSE.
   IF (space_group_code > 2 .AND. space_group_code < 16 ) THEN 
      uniqueb=(unb==1)
   ELSEIF ( space_group_code > 142 .AND. space_group_code < 168) THEN
      rhombohedral=(trig==1)
   ENDIF
   origin_choice=or
   WRITE(stdout,'(5x,"celldm ? (For instance 1.0 0.0 0.0 0.0 0.0 0.0)  ")')
   READ(5,*) celldm(1), celldm(2), celldm(3), celldm(4), celldm(5), celldm(6)
   WRITE(stdout,'(6f12.6)') celldm(1), celldm(2), celldm(3), celldm(4), &
                            celldm(5), celldm(6)
ELSE
   WRITE(stdout,'(5x," Bravais lattice code ibrav (as in QE)? ")')
   READ(5,*) ibrav
   WRITE(stdout, '(i5)') ibrav
   WRITE(stdout,'(5x," Units (alat (1) or crystal coordinates (2)) ? ")')
   READ(5,*) units
   WRITE(stdout, '(i5)') units
   IF (ibrav/=0) THEN
      WRITE(stdout,'(5x,"celldm ? (For instance 1.0 0.0 0.0 0.0 0.0 0.0)  ")')
      READ(5,*) celldm(1), celldm(2), celldm(3), celldm(4), celldm(5), celldm(6)
      WRITE(stdout,'(6f12.6)') celldm(1), celldm(2), celldm(3), celldm(4), &
                               celldm(5), celldm(6)
   ELSE
      celldm=0.0_DP
      WRITE(stdout,'(5x,"Units of at (at are multiplied by this number)? ")')
      READ(5,*) celldm(1)
      WRITE(stdout,'(f12.6)') celldm(1)
      WRITE(stdout,'(5x,"at ? ")')
      DO ivec=1,3
         READ(5,*) at(:,ivec)
         WRITE(stdout,'(3f12.6)') at(:,ivec)
      ENDDO 
      at=at*celldm(1)
   ENDIF
ENDIF

WRITE(stdout,'(5x,"Transform to conventional cell? (1=Yes, 0=No) ")') 
READ(5,*) iconv
WRITE(stdout,'(i5)') iconv
WRITE(stdout,'(5x,"Format of the input supercell? (1, 2, or 3) ")') 
READ(5,*) iformat
WRITE(stdout,'(i5)') iformat
IF (iformat==1) THEN
   WRITE(stdout,'(5x,"n1, n2, n3? (for instance 1 1 1)")')
   READ(5,*) n1, n2, n3
   WRITE(stdout,'(3i5)') n1, n2, n3
ELSEIF(iformat==2) THEN
   WRITE(stdout,'(5x,"m11, m12, m13? (for instance 1 0 0)")') 
   READ(5,*) (mm(1, ipol), ipol=1,3)
   WRITE(stdout,'(3i6)') (mm(1, ipol), ipol=1,3)
   WRITE(stdout,'(5x,"m21, m22, m23? (for instance 0 1 0)")') 
   READ(5,*) (mm(2, ipol), ipol=1,3)
   WRITE(stdout,'(3i6)') (mm(2, ipol), ipol=1,3)
   WRITE(stdout,'(5x,"m31, m32, m33? (for instance 0 0 1)")') 
   READ(5,*) (mm(3, ipol), ipol=1,3)
   WRITE(stdout,'(3i6)') (mm(3, ipol), ipol=1,3)
ELSEIF(iformat==3) THEN
   WRITE(stdout,'(5x,"R_1x, R_1y, R_1z ? (for instance 1.0 0.0 0.0)")') 
   READ(5,*) (ats(1, ipol), ipol=1,3)
   WRITE(stdout,'(3f15.7)') (ats(1, ipol), ipol=1,3)
   WRITE(stdout,'(5x,"R_2x, R_2y, R_2z ? (for instance 0.0 1.0 0.0)")') 
   READ(5,*) (ats(2, ipol), ipol=1,3)
   WRITE(stdout,'(3f15.7)') (ats(2, ipol), ipol=1,3)
   WRITE(stdout,'(5x,"R_3x, R_3y, R_3z ? (for instance 0.0 0.0 1.0)")') 
   READ(5,*) (ats(3, ipol), ipol=1,3)
   WRITE(stdout,'(3f15.7)') (ats(3, ipol), ipol=1,3)
ENDIF
WRITE(stdout,'(5x,"Centered output crystal coordinates? &
                 &(0=No, 1=Yes -0.5,0.5, 2=Yes 0,1) ")') 
READ(5,*) icenter
WRITE(stdout,'(i5)') icenter
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
WRITE(stdout,'(5x," Number of atoms and atomic coordinates? ")')
IF (which_input==1) THEN
   READ(5,*) ineq_nat

   ALLOCATE(ineq_tau(3,ineq_nat))
   ALLOCATE(label(ineq_nat))
   ALLOCATE(ineq_ityp(ineq_nat))
   ALLOCATE(if_pos(3,ineq_nat))
   ALLOCATE(rd_for(3,ineq_nat))

   DO na=1, ineq_nat
      CALL read_line( input_line )
      WRITE(stdout,*) TRIM(input_line)
!      READ(5,*) label(na), ineq_tau(1,na), ineq_tau(2,na), ineq_tau(3,na)
      CALL field_count( nfield, input_line )
      !
      ! read atom symbol (column 1)
      !
      CALL get_field(1, label(na), input_line)
      label(na) = TRIM(label(na))
      !
      !
      ! read field 2 (atom X coordinate or Wyckoff position symbol)
      !
      CALL get_field(2, field_str, input_line)
      !     
      ! Check if position na is expressed in wyckoff positions
      !
      idx = LEN_TRIM(field_str)
      IF ( (idx < 4) .AND. &
          ( IACHAR(field_str(idx:idx)) > 64 .AND. &
            IACHAR(field_str(idx:idx)) < 123 ) ) THEN

         IF ( nfield < 3 .and. nfield > 8 ) &
         CALL errore( 'supercell', 'wrong number of columns ' // &
                        & 'in ATOMIC_POSITIONS', na )
         wp=field_str
         inp(:)=1.d5
         !
         DO k = 3,MIN(nfield,5)
            ! read k-th field (coordinate k-2)
            CALL get_field(k, field_str, input_line)
            inp(k-2) = feval_infix(ierr, field_str )
            CALL errore('supercell', 'error reading field', ierr)
         ENDDO

         CALL wypos(ineq_tau(1,na),wp,inp,space_group_code, &
                 uniqueb,rhombohedral,origin_choice)

      ELSE
         ! 
         ! no wyckoff positions 
         !
         IF ( nfield /= 4 .and. nfield /= 7 ) &
         CALL errore( 'supercell', 'wrong number of columns ' // &
                        & 'in ATOMIC_POSITIONS', ia )
         !
         ! field just read is coordinate X
         !
         ineq_tau(1,na) = feval_infix(ierr, field_str )
         CALL errore('supercell','error reading field', ierr)
         DO k = 3,4
            ! read fields 3 and 4 (atom Y and Z coordinate)
            CALL get_field(k, field_str, input_line)
            ineq_tau(k-1,na) = feval_infix(ierr, field_str )
            CALL errore('supercell', 'error reading field', ierr)
         END DO
      ENDIF
   ENDDO

   CALL find_ityp(ineq_nat, label, ntyp, ineq_ityp, atm, ntypx)
   
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
   ALLOCATE(stau(3,nat))
   ALLOCATE(ityp(nat))
   ALLOCATE(label(nat))
   DO na=1, nat
      READ(5,*) label(na), tau(1,na), tau(2,na), tau(3,na)
   ENDDO


   CALL find_ityp(nat, label, ntyp, ityp, atm, ntypx)
   !
   !  If coordinates are in alat units, transform to crystal coordinates
   !
   IF (units==1) THEN
      IF (ibrav/=0) THEN
         CALL latgen(ibrav, celldm, at(1,1), at(1,2), at(1,3), omega)
         at= at / celldm(1)
      ENDIF
      CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
      CALL cryst_to_cart( nat, tau, bg, -1 )
   END IF 

   IF (ibrav==0) THEN
!
!   the at have been given using ibrav=0, we find the Bravais
!   lattice code and convert the coordinates to our conventions
!
      code_group_ext=0
      CALL find_ibrav_code(at(1,1), at(1,2), at(1,3), ibrav, celldm, &
                                          code_group_ext, ur, global_s,.TRUE.)
      WRITE(stdout,*) 'ur'
      DO ipol=1,3
         WRITE(stdout,*) (ur(ipol,jpol), jpol=1,3)
      END DO

      stau=0.0_DP
      DO na=1,nat
         DO ivec=1,3
            DO jvec=1,3
               stau(ivec,na)=stau(ivec,na) + ur(jvec, ivec) * tau(jvec,na)
            ENDDO
         ENDDO
         WRITE(stdout,*) (stau(ipol,na), ipol=1,3)
      ENDDO
      tau=stau
   ENDIF
   DEALLOCATE(stau)
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
CALL latgen(ibrav, celldm, at(1,1), at(1,2), at(1,3), omega)
at=at/celldm(1)
!
!   The primitive vectors of the bulk are
!
WRITE(stdout,'(/,5x,"Primitive lattice vectors")')
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') at(:,1)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') at(:,2)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') at(:,3)

CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )

IF (iformat==1) THEN
   ats(:,1) = n1 * at(:,1)
   ats(:,2) = n2 * at(:,2)
   ats(:,3) = n3 * at(:,3)
   mm=0
   mm(1,1) = n1
   mm(2,2) = n2
   mm(3,3) = n3
ELSEIF (iformat==2) THEN
   ats=0.0_DP
   DO ipol=1,3
      DO jpol=1,3
         ats(:,ipol) = ats(:,ipol) + mm(ipol,jpol) * at(:, jpol)
      ENDDO
   ENDDO
ELSEIF (iformat==3) THEN
!
!  bring the ats to crystal coordinates
!
   rmm=ats
   CALL cryst_to_cart(3, rmm, bg, -1)

   mm=NINT(rmm)
   check=0.0_DP
   DO ipol=1,3
      DO jpol=1,3
         check=check+ABS(rmm(ipol, jpol)-mm(ipol,jpol))
      ENDDO
   ENDDO

   IF (check > 1.D-8) CALL errore('supercell','input vectors not Bravais &
                                                    &lattice vectors',1)
   WRITE(stdout, '(/,5x,"Supercell primitive lattice vectors &
                                            &in crystal coordinates ")') 
   DO ipol=1,3
      WRITE(stdout,'(3i6)') (mm(ipol, jpol), jpol=1,3)
   ENDDO
ENDIF
!
!  compute the at of the supercell
!
WRITE(stdout,'(/,5x,"Supercell primitive lattice vectors")')
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') ats(:,1)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') ats(:,2)
WRITE(stdout,'(5x,"(",2(f15.6,","),f15.6,")")') ats(:,3)
!
!  idet is the number of bulk cells inside the supercell
!
idet = mm(1,1) * ( mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3) )-   &
       mm(1,2) * ( mm(2,1) * mm(3,3) - mm(3,1) * mm(2,3) )+   &
       mm(1,3) * ( mm(2,1) * mm(3,2) - mm(3,1) * mm(2,2) )

IF (idet < 0) CALL errore('supercell','Left handed supercell vectors &
                                                             &not allowed',1)
IF (idet == 0) CALL errore('supercell','Parallel or zero supercell vectors',1)

WRITE(stdout,'(/,5x,"The supercell contains",i5," bulk cells")') idet

IF (iformat==1) THEN
   ibravs=ibrav
   celldms(1)=n1 
   celldms(2)=( n2 * celldm(2) ) / n1
   celldms(3)=( n3 * celldm(3) ) / n1
   celldms(4)=celldm(4)
   celldms(5)=celldm(5)
   celldms(6)=celldm(6)
ELSE
!
!  Now find the Bravais lattice and the celldm of the new cell
!
   code_group_ext=0
   CALL find_ibrav_code(ats(1,1),ats(1,2),ats(1,3),ibravs,celldms, &
                                 code_group_ext, ur, global_s, .FALSE.)
!
ENDIF

CALL recips( ats(1,1), ats(1,2), ats(1,3), bgs(1,1), bgs(1,2), bgs(1,3) )


IF (ionode) THEN
   iuout=find_free_unit()
   OPEN(unit=iuout, file='supercell.pwin', status='unknown', &
                                              form='formatted')
   CALL write_cell_info ( iuout, ibravs, celldm(1), celldms)
ENDIF
CALL write_cell_info ( stdout, ibravs, celldm(1), celldms)
!
!   bring tau_new to cartesian components
!
CALL cryst_to_cart( nat_new, tau_new, at, 1 )
!
!   compute the total number of atoms
!
all_nat=nat_new * idet

WRITE(stdout,'(5x,"nat= ",i5)') all_nat

ALLOCATE(all_tau(3,all_nat))
ALLOCATE(all_ityp(all_nat))

IF (iformat==1) THEN
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
ELSE
   all_nat=nat_new
   iat=nat
   copies = 1
   all_ityp(1:nat_new)=ityp_new(1:nat_new)
   all_tau(:,1:nat_new)=tau_new(:,1:nat_new)
!
!  Build the vertices of the supercell
!
   v(:,1)=ats(:,1)
   v(:,2)=ats(:,2)
   v(:,3)=ats(:,3)
   v(:,4)=ats(:,1) + ats(:,2)
   v(:,5)=ats(:,1) + ats(:,3)
   v(:,6)=ats(:,2) + ats(:,3)
   v(:,7)=ats(:,1) + ats(:,2) + ats(:,3)
!
!  And find the most distant from the origin
!
   radius=0.0_DP
   DO iv=1,7
      rmod = v(1,iv) ** 2 + v(2,iv) ** 2 + v(3,iv) ** 2
      IF (rmod>radius) radius=rmod
   ENDDO
!
!  find the minimum n1, n2, n3 that contain the sphere
!
   n1=INT(sqrt(radius)*sqrt(bg(1,1)**2+bg(2,1)**2+bg(3,1)**2))+1
   n2=INT(sqrt(radius)*sqrt(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2))+1
   n3=INT(sqrt(radius)*sqrt(bg(1,3)**2+bg(2,3)**2+bg(3,3)**2))+1
!
!  and search all the points within this cell
!
   DO i1= -n1, n1
      DO i2= -n2, n2
         DO i3= -n3, n3
            IF ( i1==0 .AND. i2==0 .AND. i3==0) CYCLE
            rmu(:) = i1 * at(:,1) + i2 * at(:,2) + i3 * at(:,3)
            prod1 = rmu(1) * bgs(1,1) + rmu(2) * bgs(2,1) + rmu(3) * bgs(3,1)
            prod2 = rmu(1) * bgs(1,2) + rmu(2) * bgs(2,2) + rmu(3) * bgs(3,2)
            prod3 = rmu(1) * bgs(1,3) + rmu(2) * bgs(2,3) + rmu(3) * bgs(3,3)
            IF ( prod1 >=-eps .AND. prod1 < (1.0_DP-eps) .AND. &
                 prod2 >=-eps .AND. prod2 < (1.0_DP-eps) .AND. &
                 prod3 >=-eps .AND. prod3 < (1.0_DP-eps) ) THEN

               copies=copies+1
               DO na=1,nat_new
                  all_nat=all_nat+1
                  all_tau(:,all_nat) = tau_new(:,na) + rmu(:) 
                  all_ityp(all_nat) = ityp_new(na)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   IF (copies /= idet) CALL errore('supercell','Some problem with copies',1)
!
!   Bring to crystal coordinates
!
   CALL cryst_to_cart( all_nat, all_tau, bgs, -1 )
ENDIF

IF (icenter == 1.OR.icenter==2) THEN
!
!   Bring all crystal coordinates between -0.5 and 0.5 or between 0 and 1
!
   DO na=1,all_nat
      DO ipol=1,3
         all_tau(ipol,na)=all_tau(ipol,na)-NINT(all_tau(ipol,na))
         IF (all_tau(ipol,na)<0.0_DP.AND.icenter==2) &
                             all_tau(ipol,na)=all_tau(ipol,na)+1.0_DP
      ENDDO
   ENDDO
ENDIF
!
!  and write on output the coordinates
!
WRITE(stdout,'("ATOMIC_POSITIONS (crystal)")')
IF (ionode) WRITE(iuout,'("ATOMIC_POSITIONS (crystal)")')
DO na=1,all_nat
   WRITE(stdout,'(a,3f21.13)') TRIM(atm(all_ityp(na))), all_tau(1,na), &
                                                        all_tau(2,na), &
                                                        all_tau(3,na) 
   IF (ionode) WRITE(iuout,'(a,3f21.13)') TRIM(atm(all_ityp(na))), &
                                                all_tau(1,na), &
                                                all_tau(2,na), &
                                                all_tau(3,na) 
ENDDO
IF (ionode) CLOSE(iuout)

IF (ionode) THEN
   OPEN(unit=iuout, file='supercell.xsf', status='unknown', &
                                              form='formatted')

   CALL cryst_to_cart( all_nat, all_tau, ats, 1 )
   CALL xsf_struct (celldm(1), ats, all_nat, all_tau, atm, all_ityp, iuout)
   CLOSE(iuout)
ENDIF

DEALLOCATE(ityp)
DEALLOCATE(tau)
DEALLOCATE(ityp_new)
DEALLOCATE(tau_new)
DEALLOCATE(all_tau)
DEALLOCATE(all_ityp)

CALL environment_end( code )
CALL mp_global_end ()

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
   tau(2,na) = - 0.5_DP * ( tau_new(2,na) - tau_new(3,na) )
   tau(3,na) =   0.5_DP * ( tau_new(2,na) + tau_new(3,na) )
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
   tau(1,na) =   0.5_DP * ( tau_new(1,na) - tau_new(2,na) - tau_new(3,na) )
   tau(2,na) =   0.5_DP * ( tau_new(1,na) + tau_new(2,na) - tau_new(3,na) )
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

SUBROUTINE write_cell_info ( iunout, ibrav, units, celldm)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: iunout, ibrav
REAL(DP), INTENT(IN) :: units, celldm(6)

WRITE(iunout,'(5x,"ibrav= ",i3)') ibrav
WRITE(iunout,'(5x, "celldm(1)= ",f15.9)') celldm(1) * units
IF (celldm(2) > 0.0_DP ) &
   WRITE(iunout,'(5x, "celldm(2)= ",f15.9)') celldm(2) 
IF (celldm(3) > 0.0_DP ) &
   WRITE(iunout,'(5x, "celldm(3)= ",f15.9)') celldm(3)
IF (celldm(4) /= 0.0_DP ) &
   WRITE(iunout,'(5x, "celldm(4)= ",f15.9)') celldm(4) 
IF (celldm(5) /= 0.0_DP ) &
   WRITE(iunout,'(5x, "celldm(5)= ",f15.9)') celldm(5) 
IF (celldm(6) /= 0.0_DP ) &
   WRITE(iunout,'(5x, "celldm(6)= ",f15.9)') celldm(6) 

RETURN

END SUBROUTINE write_cell_info
