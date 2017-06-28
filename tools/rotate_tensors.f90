!
! Copyright (C) 2017 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM rotate_tensors
!
!  This program reads a rotation matrix that transform the coordinates
!  of a point in reference system 1 into those of a point in reference system 2.
!  Then it reveives a tensor of a rank from 1 to 4 which is
!  supposed to be written in the reference system 1 or 2. It gives as
!  output the tensor components in the other reference system.
!  The input of the code is the following: 
!
!  R_11, R_12, R_13  : the first row of the rotation matrix
!  R_21, R_22, R_23  : the second row of the rotation matrix
!  R_31, R_32, R_33  : the third row of the rotation matrix
!  rank     : the rank of the tensor
!  The format of the tensor depends on its rank and is the following:
!  rank 1
!  a_1, a_2, a_3
!  rank 2 
!  a_11, a_12, a_13, 
!  a_21, a_22, a_23
!  a_31, a_32, a_33
!  rank 3
!  d_11  d_12  d_13   d_14   d_15   d_16
!  d_21  d_22  d_23   d_24   d_25   d_26
!  d_31  d_32  d_33   d_34   d_35   d_36
!  rank 4
!  c_11  c_12  c_13   c_14   c_15   c_16
!  c_21  c_22  c_23   c_24   c_25   c_26
!  c_31  c_32  c_33   c_34   c_35   c_36
!  c_41  c_42  c_43   c_44   c_45   c_46
!  c_51  c_52  c_53   c_54   c_55   c_56
!  c_61  c_62  c_63   c_64   c_65   c_66
! 
!  1 or 2  the tensor is in the coordinates of reference 1 (1) or 2 (2)
!        R is applied in case 1, R^-1 in case 2.
!
USE kinds,  ONLY : DP
USE rotate, ONLY : is_rotation, rotate_vect, rotate_tensors2, &
                   rotate_tensors3, rotate_tensors4
USE voigt,  ONLY : to_voigt3, to_voigt4

USE mp_global,     ONLY : mp_startup, mp_global_end
USE environment,   ONLY : environment_start, environment_end
USE io_global,     ONLY : stdout

IMPLICIT NONE
INTEGER :: nat
REAL(DP) :: a(3), rot(3,3), tensor1(3), tensor2(3,3), tensor3(3,3,3),  &
            tensor4(3,3,3,3), rtensor1(3), rtensor2(3,3), rtensor3(3,3,3), &
            rtensor4(3,3,3,3), tensor3v(3,6), tensor4v(6,6)
INTEGER :: na, ipol, jpol, rank

CHARACTER(LEN=9) :: code='r_tensor' 

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'(5x,"Rotation matrix?")') 
DO ipol=1,3
   READ(5,*) (rot(ipol,jpol), jpol=1,3)
   WRITE(stdout,'(3f15.8)') (rot(ipol,jpol), jpol=1,3)
ENDDO

IF (.NOT.is_rotation(rot)) THEN
   WRITE(stdout,'(/,5x,"WARNING: input matrix not a rotation")')
END IF

WRITE(stdout,'(5x,"Rank of the tensor")')
READ(5,*) rank
WRITE(stdout,'(i5)') rank

IF (rank==1) THEN
   WRITE(stdout,'(5x,"Enter the tensor components")')
   READ(5,*) (tensor1(ipol), ipol=1,3)
   WRITE(stdout,'(3f15.8)') (tensor1(ipol), ipol=1,3)
ELSEIF (rank==2) THEN
   WRITE(stdout,'(5x,"Enter the tensor components")')
   DO ipol=1,3
      READ(5,*) (tensor2(ipol, jpol), jpol=1,3)
      WRITE(stdout,'(3f15.8)') (tensor2(ipol,jpol), jpol=1,3)
   ENDDO
ELSEIF (rank==3) THEN
   WRITE(stdout,'(5x,"Enter the tensor components")')
   DO ipol=1,3
      READ(5,*) (tensor3v(ipol, jpol), jpol=1,6)
      WRITE(stdout,'(6f12.8)') (tensor3v(ipol,jpol), jpol=1,6)
   ENDDO
ELSEIF (rank==4) THEN
   WRITE(stdout,'(5x,"Enter the tensor components")')
   DO ipol=1,6
      READ(5,*) (tensor4v(ipol, jpol), jpol=1,6)
      WRITE(stdout,'(6f12.5)') (tensor4v(ipol,jpol), jpol=1,6)
   ENDDO
ENDIF

IF (rank==1) THEN
   CALL rotate_vect(rot, 1, tensor1, rtensor1, 1) 
ELSEIF(rank==2) THEN
   CALL rotate_tensors2(rot, 1, tensor2, rtensor2, 1) 
ELSEIF(rank==3) THEN
   CALL to_voigt3( tensor3v, tensor3, .FALSE. )
   CALL rotate_tensors3(rot, 1, tensor3, rtensor3, 1) 
   CALL to_voigt3( tensor3v, tensor3, .FALSE. )

   WRITE(stdout,'(/,5x,"Rotated tensor of rank 3")') 
   DO ipol=1,3
      WRITE(stdout,'(6f12.5)') (tensor3v(ipol,jpol), jpol=1,6)
   ENDDO
ELSEIF(rank==4) THEN
   CALL to_voigt4( tensor4v, tensor4, .FALSE. )
   CALL rotate_tensors4( rot, 1, tensor4, rtensor4, 1 )  
   CALL to_voigt4( tensor4v, rtensor4, .TRUE. )

   WRITE(stdout,'(/,5x,"Rotated tensor of rank 4")') 
   DO ipol=1,6
      WRITE(stdout,'(6f12.5)') (tensor4v(ipol,jpol), jpol=1,6)
   ENDDO
ENDIF

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM rotate_tensors
