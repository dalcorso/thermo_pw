!
! Copyright (C) 2018 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE atomic_pos
!
! This module contains routines to manipulate atomic positions
!
! Presently it contains a routine that receives a list of atoms and
! their label and finds the number of types, sets the indeces ityp
! that gives the type of each atom and atm the label of each type.
!
!
PRIVATE
SAVE

PUBLIC find_ityp

CONTAINS
!------------------------------------------------------------------
   SUBROUTINE find_ityp(nat, label, ntyp, ityp, atm, ntypx)
!------------------------------------------------------------------

   IMPLICIT NONE 

   INTEGER, INTENT(IN) :: nat, ntypx
   CHARACTER(LEN=3), INTENT(IN) :: label(nat)

   INTEGER, INTENT(INOUT) :: ityp(nat)
   INTEGER, INTENT(OUT) :: ntyp

   CHARACTER(LEN=3), INTENT(INOUT) :: atm(ntypx)

   INTEGER :: na, nb, nt
   LOGICAL :: found
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
         IF (ntyp > ntypx) CALL errore('find_ityp','too many types',1)
         ityp(na)=ntyp
         atm(ntyp)=label(na)
      ENDIF
   ENDDO

   RETURN
   END SUBROUTINE find_ityp

END MODULE atomic_pos
