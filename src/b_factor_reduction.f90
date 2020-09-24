! Copyright (C) 2018 Cristiano Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE b_factor_reduction(b_rest)
!---------------------------------------------------------------------
!
!  This is a driver to find the simmetry restrictions
!  to B factor matrix
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat, tau, ityp
USE symme,            ONLY : symtensor
USE cell_base,        ONLY : ibrav
USE symm_base,        ONLY : s, t_rev, irt, nrot, nsym, invsym, nosym, &
                             set_sym_bl, find_sym
USE noncollin_module, ONLY : m_loc, noncolin
USE random_numbers,   ONLY : randy
USE spin_orb,         ONLY : domag
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode, stdout

IMPLICIT NONE
INTEGER  :: ierr, system, na
INTEGER  :: b_rest(nat)
CHARACTER(LEN=6) :: int_to_char
REAL(DP) :: random_bmatrix(3,3,nat), rnumb, toll, b11, b22, b33, b12, b13, b23 
REAL(DP) :: db12, db23, db1211, db1323, db1213, db2212, db2313 ! d means difference
REAL(DP) :: sb1323, sb1213 ! s means sum
REAL(DP) :: arg ! for random number generation
!INTEGER  :: time(8), seed(12), ipol, jpol
INTEGER  :: ipol, jpol
LOGICAL  :: magnetic_sym

toll=10E-8

magnetic_sym=noncolin.AND.domag
ALLOCATE( m_loc( 3, nat ) )
m_loc=0.0_DP

CALL set_sym_bl ( )
!
! ... eliminate rotations that are not symmetry operations
!
CALL find_sym ( nat, tau, ityp, magnetic_sym, m_loc )

!CALL date_and_time(VALUES=time)
!seed(1) = time(4)*(36000*time(5)+6000*time(6)+100*time(7)+time(8))
!CALL random_seed(PUT=seed)
arg=randy(0)

DO na=1, nat
   DO ipol=1, 3
      DO jpol=ipol, 3
!         CALL random_number(rnumb)
         rnumb=2.0_DP*randy()-1.0_DP
         random_bmatrix(ipol,jpol,na)=rnumb
         IF (ipol/=jpol) random_bmatrix(jpol,ipol,na)=random_bmatrix(ipol,jpol,na)
         !WRITE(stdout,'(/,5x,"Matrix element:",e12.5)') random_bmatrix(ipol,jpol,na)
      END DO
   END DO
END DO

CALL symtensor(nat, random_bmatrix)

DO na=1, nat
   
   b11=random_bmatrix(1,1,na)
   b22=random_bmatrix(2,2,na)
   b33=random_bmatrix(3,3,na)
   b12=random_bmatrix(1,2,na)
   b13=random_bmatrix(1,3,na)
   b23=random_bmatrix(2,3,na)
   db12=ABS(b11-b22)
   db23=ABS(b22-b33)
   db1211=ABS(2*b12-b11) ! Here, a factor 2 is present
   db1323=ABS(b13-b23)
   db1213=ABS(b12-b13)
   db2212=ABS(b22-2*b12) ! Here, a factor 2 is present
   db2313=ABS(b23-2*b13) ! Here, a factor 2 is present
   sb1323=ABS(b13+b23)
   sb1213=ABS(b12+b13)

   IF (db12.LT.toll) THEN

      IF (db23.LT.toll) THEN
        
         IF (ABS(b12).LT.toll) THEN
            b_rest(na) = 17

         ELSE
            b_rest(na) = 18

         END IF

      ELSE IF (db1211.LT.toll) THEN
         b_rest(na) = 16

      ELSE IF (ABS(b12).LT.toll) THEN
         b_rest(na) = 8
 
      ELSE IF (ABS(b13).LT.toll) THEN
         b_rest(na) = 5

      ELSE IF (db1323.LT.toll) THEN
         b_rest(na) = 6

      ELSE IF (sb1323.LT.toll) THEN
         b_rest(na) = 7
      
      END IF

   ELSE IF (db23.LT.toll) THEN
      
      IF (ABS(b12).LT.toll) THEN
         
         IF (ABS(b23).LT.toll) THEN
            b_rest(na) = 12

         ELSE
            b_rest(na) = 9

         END IF
      
      ELSE IF (db1213.LT.toll) THEN
         b_rest(na) = 10

      ELSE IF (sb1213.LT.toll) THEN
         b_rest(na) = 11

      END IF

   ELSE IF (db2212.LT.toll) THEN

      IF (ABS(b13).LT.toll) THEN
         b_rest(na) = 14

      ELSE IF (ABS(b23).LT.toll) THEN
         b_rest(na) = 13

      ELSE IF (db2313.LT.toll) THEN
         b_rest(na) = 15

      END IF

   ELSE IF (ABS(b12).LT.toll) THEN
      
      IF (ABS(b13).LT.toll) THEN
      
         IF (ABS(b23).LT.toll) THEN
            b_rest(na) = 4

         ELSE
            b_rest(na) = 3

         END IF 

      ELSE
         b_rest(na) = 1

      END IF

   ELSE IF (ABS(b13).LT.toll) THEN
      b_rest(na) = 2
    
   END IF

END DO

DEALLOCATE(m_loc)
DEALLOCATE(irt)

RETURN
END SUBROUTINE b_factor_reduction
