!
! Copyright (C) 2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM kovalev
!
!  This program contains a few useful tools to help the reading of the
!  Kovalev tables for space group representations
!
USE kinds, ONLY : DP
USE constants, ONLY : pi
USE rotate, ONLY : euler_to_su2
USE point_group, ONLY : sym_label, find_double_product_table_from_sym
USE io_global, ONLY : stdout
IMPLICIT NONE

REAL(DP), PARAMETER :: tmp = 3.0_DP * pi / 2.0_DP, pm = pi / 2.0_DP, &
                       pt = pi / 3.0_DP, dtp = 2.0_DP * pt, &
                       qpt = 4.0_DP * pt, cpt = 5.0_DP * pt

REAL(DP) :: kovalev_euler_cubic_psi(24)
DATA kovalev_euler_cubic_psi / 0.0_DP, 0.0_DP, 0.0_DP, &
                                   pi,     pi, 0.0_DP, &
                                   pi, 0.0_DP,     pm, &
                                  tmp,    tmp,    tmp, &
                                  tmp,     pm,    tmp, &
                                   pm,   0.0_DP,   pi, &
                               0.0_DP,       pi,  tmp, &
                                  tmp,       pm,   pm   /

REAL(DP) :: kovalev_euler_cubic_theta(24)
DATA kovalev_euler_cubic_theta / 0.0_DP,     pi,   pi, &
                                 0.0_DP,     pm,   pm, &
                                     pm,     pm,   pm, &
                                    tmp,     pm,   pm, &
                                     pi, 0.0_DP, 0.0_DP, &
                                     pi,     pm,   pm,   &
                                     pm,     pm,   pm,   &
                                     pm,     pm,   pm   /

REAL(DP) :: kovalev_euler_cubic_phi(24)
DATA kovalev_euler_cubic_phi / 0.0_DP,   0.0_DP,     pi,  &
                               0.0_DP,       pm,    tmp,  &
                                  tmp,       pm, 0.0_DP,  &
                               0.0_DP,       pi, 0.0_DP,  &
                               0.0_DP,   0.0_DP, 0.0_DP,  &
                               0.0_DP,       pi, 0.0_DP,  &
                               0.0_DP,       pi,    tmp,  &
                                   pm,       pm,    tmp /

REAL(DP) :: kovalev_euler_hex_psi(12)
DATA kovalev_euler_hex_psi /   0.0_DP,       pt,    dtp,  &
                                   pi,      qpt,    cpt, &
                               0.0_DP,       pt,    dtp,  &
                                   pi,       qpt,   cpt /

REAL(DP) :: kovalev_euler_hex_theta(12)
DATA kovalev_euler_hex_theta /   0.0_DP,    0.0_DP,    0.0_DP,  &
                                 0.0_DP,    0.0_DP,    0.0_DP,  &
                                   pi,        pi,        pi,  &
                                   pi,        pi,        pi  /
REAL(DP) :: kovalev_euler_hex_phi(12)

INTEGER :: kovalev_cubic(48)
DATA kovalev_cubic / 1,   4,  3,  2, 17, 19, 18, 20, 21, 22, 24, 23, &
                     6,   8,  7,  5, 14, 13, 16, 15, 10, 12,  9, 11, &
                     33, 36, 35, 34, 49, 51, 50, 52, 53, 54, 56, 55, &
                     38, 40, 39, 37, 46, 45, 47, 48, 42, 44, 41, 43  /

INTEGER :: kovalev_hexagonal(24)
DATA kovalev_hexagonal / 1,  25, 27,  2, 28, 26, 4,  30, 32,  3, 31, 29,  &
                        33,  57, 59, 34, 60, 58, 36, 62, 64, 35, 63, 61 /

!  DATA kovalev_hexagonal / 1,  25, 27,  2, 28, 26, 31, 29,  4, 30, 32, 3, &
!                          33,  57, 59, 34, 60, 58, 63, 61, 36, 62, 64, 35 /


COMPLEX(DP) :: a, b
REAL(DP) :: psi, theta, phi, arg
INTEGER :: work_choice, isym, jsym, start, last, nsym, ntables, itables, &
           group_index_ext
INTEGER :: stdin
INTEGER :: prd(48,48), epos(48,48), group_desc(48)


kovalev_euler_hex_phi=0.0_DP

WRITE(stdout,'(/,5x,"Choose what to write")')
WRITE(stdout,'(5x,"1) Write Kovalev symmetry operations list")')
!WRITE(stdout,'(5x,"2) Write Kovalev cubic product table")')
!WRITE(stdout,'(5x,"3) Write Kovalev hexagonal product table")')
!WRITE(stdout,'(5x,"4) Write Kovalev su2 matrices")')

stdin=5
READ(stdin,*) work_choice


IF (work_choice == 1 ) THEN
WRITE(stdout,'(/,5x, "Kovalev symmetries - cubic groups",/)')
DO isym=1,48
   WRITE(stdout,'(5x,i5," - ",i3,3x,a8)') isym, kovalev_cubic(isym), &
                                        sym_label(kovalev_cubic(isym))
ENDDO
WRITE(stdout,'(/,5x, "Kovalev symmetries - hexagonal groups",/)')
DO isym=1,24
   WRITE(stdout,'(5x,i5," - ",i3,3x,a8)') isym, kovalev_hexagonal(isym), &
                                        sym_label(kovalev_hexagonal(isym))
ENDDO

!ELSEIF (work_choice == 2 .OR. work_choice==3 ) THEN
!IF (work_choice==2) THEN
!   nsym=48
!   group_index_ext=136
!   group_desc=kovalev_cubic
!ELSE
!   nsym=24
!   group_index_ext=111
!   group_desc(1:24)=kovalev_hexagonal
!ENDIF
!
!CALL find_double_product_table_from_sym_kov(prd, epos, group_desc, nsym)
!WRITE(stdout,'(/,5x, "The double group product table",/)')
!ntables= nsym / 8
!IF (MOD(nsym,8) /= 0 ) ntables=ntables+1
!DO itables=1,ntables
!   start=(itables-1)*8+1
!   last=MIN(itables*8,nsym)
!   WRITE(stdout,'(8x,8(1x,a8))') (sym_label(group_desc(jsym)),jsym=start,last)
!   DO isym=1,nsym
!      WRITE(stdout,'(a8,i3,2x,7(i7,2x))') sym_label(group_desc(isym)), &
!     (prd(isym, jsym)*epos(isym,jsym),jsym=start,last)
!   ENDDO
!   WRITE(stdout,*)
!ENDDO
!WRITE(stdout,'(/,5x,"Row x column multiplication table")')
!WRITE(stdout,'(5x,"- means multiplication by -E:")')
!WRITE(stdout,'(5x, "Point group product table: neglect the minus signs")')
!
!ELSEIF (work_choice == 4 ) THEN
!
!   WRITE(stdout,*) 'Cayley-Klein a and b parameters -- cubic'  
!   DO isym=1,24
!      psi=kovalev_euler_cubic_psi(isym)
!      theta=kovalev_euler_cubic_theta(isym)
!      phi=kovalev_euler_cubic_phi(isym)
!
!      CALL euler_to_su2(psi,theta,phi,a,b)
!      WRITE(stdout,'(2i5,3x,a8,4f13.5)') isym, kovalev_cubic(isym), &
!                              sym_label(kovalev_cubic(isym)), a, b
!   ENDDO
!
!   WRITE(stdout,*) 'Cayley-Klein a and b parameters -- hexagonal'  
!
!   DO isym=1,12
!      psi=kovalev_euler_hex_psi(isym)
!      theta=kovalev_euler_hex_theta(isym)
!      phi=kovalev_euler_hex_phi(isym)
!
!      CALL euler_to_su2(psi,theta,phi,a,b)
!      WRITE(stdout,'(2i5,3x,a8,4f13.5)') isym, kovalev_hexagonal(isym), &
!                              sym_label(kovalev_hexagonal(isym)), a, b
!   ENDDO

ENDIF

CONTAINS

!-------------------------------------------------------------------------
SUBROUTINE set_sym_su2_kov(sym_num, smat, sinv)
!-------------------------------------------------------------------------
!
!  This routine uses the Cayley-Klein parameters to set the su2 rotation
!  matrices for the 32 proper rotations defined in the module. 
!  See sym_label for the name of each rotation. The integer sinv is set 
!  to 1 if the operation has not inversion, to -1 if the operation is to 
!  be multiplied by inversion. Note that inversion does not act on spin, 
!  so the su2 matrices remain the same. To multiply by E' just change
!  the sign of the matrix. The choices of the su2 matrix corresponding
!  to each rotation is done following Kovalev Irreducible representations
!  of the space groups.
!  The order of the operation is that of the point_group module. Use
!  the appropriate array to obtain the kovalev order.
!
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sym_num
  INTEGER, INTENT(OUT) :: sinv
  COMPLEX(DP), INTENT(OUT) :: smat(2,2)

  INTEGER :: snum
  COMPLEX(DP) :: a(32), b(32)
  REAL(DP) :: sqrt2=SQRT(2.0_DP), sqrt3=SQRT(3.0_DP)

  IF (sym_num < 1 .OR. sym_num > 64) CALL errore('set_sym_su2', &
              'problem with symmetry number',1)

!
!  the entries with a ! differ from those of the point group module
!
  a=(0.0_DP,0.0_DP)
  b=(0.0_DP,0.0_DP)

  a(1)=(1.0_DP,0.0_DP)
  a(2)=(0.0_DP,-1.0_DP)
!
  b(3)=(1.0_DP,0.0_DP)
  b(4)=(0.0_DP,-1.0_DP)
  b(5)=(-1.0_DP,-1.0_DP) / sqrt2
!
  b(6)=(-1.0_DP,1.0_DP) / sqrt2
!
  a(7)=(-1.0_DP,-1.0_DP) / sqrt2
  a(8)=(1.0_DP,-1.0_DP) / sqrt2
  a(9)=(0.0_DP,-1.0_DP) / sqrt2
  b(9)=(0.0_DP,-1.0_DP) / sqrt2
!
  a(10)=(0.0_DP,1.0_DP) / sqrt2
  b(10)=(0.0_DP,-1.0_DP) / sqrt2
!
  a(11)=(-1.0_DP,0.0_DP) / sqrt2
  b(11)=(1.0_DP,0.0_DP) / sqrt2
!
  a(12)=(-1.0_DP,0.0_DP) / sqrt2
  b(12)=(-1.0_DP,0.0_DP) / sqrt2
  a(13)=(0.0_DP,-1.0_DP) / sqrt2
  b(13)=(-1.0_DP,0.0_DP) / sqrt2
  a(14)=(0.0_DP,-1.0_DP) / sqrt2
  b(14)=(1.0_DP,0.0_DP) / sqrt2
!
  a(15)=(-1.0_DP,0.0_DP) / sqrt2
  b(15)=(0.0_DP,-1.0_DP) / sqrt2
  a(16)=(1.0_DP,0.0_DP) / sqrt2
  b(16)=(0.0_DP,-1.0_DP) / sqrt2
!
  a(17)=(-1.0_DP,-1.0_DP) / 2.0_DP
  b(17)=(-1.0_DP,-1.0_DP) / 2.0_DP
!
  a(18)=(-1.0_DP,1.0_DP) / 2.0_DP
  b(18)=(1.0_DP,-1.0_DP) / 2.0_DP
!
  a(19)=(-1.0_DP,-1.0_DP) / 2.0_DP
  b(19)=(1.0_DP,1.0_DP) / 2.0_DP
  a(20)=(1.0_DP,-1.0_DP) / 2.0_DP
  b(20)=(1.0_DP,-1.0_DP) / 2.0_DP
  a(21)=(1.0_DP,-1.0_DP) / 2.0_DP
  b(21)=(-1.0_DP,-1.0_DP) / 2.0_DP
  a(22)=(1.0_DP,1.0_DP) / 2.0_DP
  b(22)=(-1.0_DP,1.0_DP) / 2.0_DP
!
  a(23)=(-1.0_DP,-1.0_DP) / 2.0_DP
  b(23)=(-1.0_DP,1.0_DP) / 2.0_DP
!
  a(24)=(-1.0_DP,1.0_DP) / 2.0_DP
  b(24)=(-1.0_DP,-1.0_DP) / 2.0_DP
  a(25)=CMPLX(sqrt3,-1.0_DP) / 2.0_DP

  a(26)=-CMPLX(sqrt3,1.0_DP) / 2.0_DP
  a(27)=CMPLX(1.0_DP,-sqrt3) / 2.0_DP
!
  a(28)=-CMPLX(1.0_DP,sqrt3) / 2.0_DP
!
  b(29)=-CMPLX(1.0_DP,-sqrt3) / 2.0_DP
  b(30)=CMPLX(-1.0_DP,-sqrt3) / 2.0_DP
!
  b(31)=-CMPLX(sqrt3,-1.0_DP) / 2.0_DP
  b(32)=CMPLX(-sqrt3,-1.0_DP) / 2.0_DP

  sinv=1
  snum=sym_num
  IF (sym_num > 32) THEN
     sinv=-1
     snum=sym_num-32
  ENDIF
  smat(1,1) = a(snum)
  smat(2,2) = CONJG(a(snum))
  smat(1,2) = b(snum)
  smat(2,1) = -CONJG(b(snum))

  RETURN
  END SUBROUTINE set_sym_su2_kov

!-------------------------------------------------------------------------
 SUBROUTINE find_double_product_table_from_sym_kov(prd, epos, group_desc, nsym)
!-------------------------------------------------------------------------
!
!  This routine provides the product table prd of a given double group.
!  Each entry prd(i,j) is the index of the operation that results from
!  the product of S_i S_j (in this order).
!  The symmetry operations are only those without E'. If necessary
!  the table can be extended to the entire group using EE=E, EE'=E'E=E', 
!  E'E'=E.
!
!  The routine provides also the factor system epos of the point group
!  considered as a projective representation of the double group.
!
 USE kinds, ONLY : DP
 IMPLICIT NONE
 INTEGER, INTENT(OUT) :: prd(48,48), epos(48,48)

 INTEGER :: group_desc(48), sinv(48), isym, jsym, ksym, lsym, nsym
 COMPLEX(DP) :: group_mat(2,2,48)
  
 epos=1
 prd=0
 DO isym=1, nsym
    DO jsym=1, nsym
       CALL product_sym_su2_kov(group_desc(isym), group_desc(jsym), &
                            ksym, epos(isym,jsym))
       DO lsym=1,nsym
          IF (group_desc(lsym)==ksym) prd(isym,jsym)=lsym
       ENDDO
    END DO
 END DO

 RETURN
 END SUBROUTINE find_double_product_table_from_sym_kov

!-------------------------------------------------------------------------
 SUBROUTINE product_sym_su2_kov(isym, jsym, prd, epos)
!-------------------------------------------------------------------------
 !
 !  This routine recives the indeces of two symmetry operations, the
 !  list of symmetry operations in su2 form and gives the index of the
 !  product inside the group list and the possible -E operation
 !
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: isym, jsym
 INTEGER, INTENT(OUT) :: prd, epos

 COMPLEX(DP) :: a_mat(2,2), b_mat(2,2), c_mat(2,2), group_mat(2,2,64)
 INTEGER :: sinv(64)
 REAL(DP) :: diff, diff1
 INTEGER :: sinv1, sinv2, sinvp
 INTEGER :: pcount, ksym

 DO ksym=1,64
    CALL set_sym_su2_kov(ksym, group_mat(1,1,ksym), sinv(ksym))
!    WRITE(stdout,*) 'ksym', ksym, group_desc(ksym)
 ENDDO

 a_mat(:,:)=group_mat(:,:,isym)
 sinv1=sinv(isym)
 b_mat(:,:)=group_mat(:,:,jsym)
 sinv2=sinv(jsym)

 sinvp = sinv1*sinv2
 c_mat(1,1)=a_mat(1,1) * b_mat(1,1) + a_mat(1,2) * b_mat(2,1)
 c_mat(1,2)=a_mat(1,1) * b_mat(1,2) + a_mat(1,2) * b_mat(2,2)

 epos=1
 prd=0
 pcount=0
 DO ksym = 1, 64
    IF (sinv(ksym)==sinvp) THEN
       diff=ABS(c_mat(1,1)-group_mat(1,1,ksym))+ &
            ABS(c_mat(1,2)-group_mat(1,2,ksym))
       IF (diff < 1.D-6) THEN
          prd=ksym
          pcount=pcount+1
       ELSE
          diff1=ABS(c_mat(1,1)+group_mat(1,1,ksym))+ &
                ABS(c_mat(1,2)+group_mat(1,2,ksym))
          IF (diff1 < 1.D-6) THEN
             prd=ksym
             epos=-1
             pcount=pcount+1
          ENDIF
       ENDIF
    END IF
 END DO
 IF (pcount/=1) &
    CALL errore('product_sym_su2','The product of these matrices is &
                                   &not among the allowed symmetries',1)
 RETURN
 END SUBROUTINE product_sym_su2_kov
!

END PROGRAM kovalev


