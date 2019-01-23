!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_inverse_s(nsym, s, invs)
!
!  This routine is identical to the one contained in symm_base, but
!  has explicit arguments. Given the number of symmetries and the
!  symmetry matrices it fills invs with the index of the inverse of 
!  each matrix
!
IMPLICIT NONE
INTEGER :: s(3,3,48), nsym, invs(48)

INTEGER :: isym, jsym, ss(3,3)
LOGICAL :: found
!
!   Find the inverse of each matrix
!
DO isym = 1, nsym
   found = .FALSE.
   DO jsym = 1, nsym
      !
      ss = MATMUL (s(:,:,jsym),s(:,:,isym))
      ! s(:,:,1) is the identity
      IF ( ALL ( s(:,:,1) == ss(:,:) ) ) THEN
         invs (isym) = jsym
         found = .TRUE.
      ENDIF
   ENDDO
   IF ( .NOT.found) CALL errore ('find_inverse_s', ' Not a group', 1)
ENDDO

RETURN
END SUBROUTINE find_inverse_s

