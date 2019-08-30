!
! Copyright (C) 2013-2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE clean_ngeo(ngeo,fact_ngeo,ngeo_ph,ibrav)
!-----------------------------------------------------------------------
!
!  This routine cleans the ngeo variable, setting 1 on all the values
!  that correspond to crystallographic parameters that are fixed in
!  a given ibrav. Moreover if ngeo is zero for some parameter needed
!  by ibrav, the routine sets a default value of ndefault, presently set
!  to 5
!
USE control_mur, ONLY : lmurn
IMPLICIT NONE

INTEGER, INTENT(IN)    :: ibrav
INTEGER, INTENT(INOUT) :: ngeo(6), fact_ngeo(6), ngeo_ph(6)

INTEGER :: i, ndefault, ngeo_aux(6), fact_ngeo_aux(6)
LOGICAL :: clean_fact

ngeo_aux=1
ngeo_aux(1)=ngeo(1)
fact_ngeo_aux=1
fact_ngeo_aux(1)=fact_ngeo(1)
ndefault=5

IF (.NOT.lmurn) THEN
   SELECT CASE (ibrav)
      CASE(1,2,3)
      CASE (4,6,7)
         IF (ngeo(3) /= 0) THEN
             ngeo_aux(3)=ngeo(3)
         ELSE
             ngeo_aux(3)=ndefault
         ENDIF
         fact_ngeo_aux(3)=fact_ngeo(3)
      CASE (5)
         IF (ngeo(4) /= 0) THEN
            ngeo_aux(4)=ngeo(4)
         ELSE
            ngeo_aux(4)=ndefault
         ENDIF
         fact_ngeo_aux(4)=fact_ngeo(4)
      CASE(8,9,-9,91,10,11)
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=ndefault
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=ndefault
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
      CASE(12,13)
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=ndefault
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=ndefault
         ENDIF
         IF (ngeo(4) /= 0) THEN
            ngeo_aux(4)=ngeo(4)
         ELSE
            ngeo_aux(4)=ndefault
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
         fact_ngeo_aux(4)=fact_ngeo(4)
      CASE(-12,-13)   
         IF (ngeo(2) /= 0) THEN
            ngeo_aux(2)=ngeo(2)
         ELSE
            ngeo_aux(2)=ndefault
         ENDIF
         IF (ngeo(3) /= 0) THEN
            ngeo_aux(3)=ngeo(3)
         ELSE
            ngeo_aux(3)=ndefault
         ENDIF
         IF (ngeo(5) /= 0) THEN
            ngeo_aux(5)=ngeo(5)
         ELSE
            ngeo_aux(5)=ndefault
         ENDIF
         fact_ngeo_aux(2)=fact_ngeo(2)
         fact_ngeo_aux(3)=fact_ngeo(3)
         fact_ngeo_aux(5)=fact_ngeo(5)
   CASE DEFAULT

!  If the Bravais lattice is unkown, 14 or 0 we let the user choose
!
         ngeo_aux=ngeo
         DO i=2,6
            IF (ngeo_aux(i)==0) ngeo_aux(i)=ndefault
         ENDDO
         fact_ngeo_aux=fact_ngeo
   END SELECT
ENDIF

ngeo=ngeo_aux
fact_ngeo=fact_ngeo_aux
!
!  if ngeo_ph has been set, check if it is compatible with ngeo and
!  clean fact_ngeo which is not used.
!
clean_fact=.FALSE.
DO i=1,6
   IF (ngeo_ph(i)<=0) ngeo_ph(i)=ngeo(i)
   IF (ngeo_ph(i)>ngeo(i)) ngeo_ph(i)=ngeo(i)
   IF (MOD(ngeo_ph(i),2) /= MOD(ngeo(i),2)) &
               CALL errore('clean_ngeo','ngeo_ph incompatible with ngeo',1)
   IF (ngeo_ph(i)/=ngeo(i)) clean_fact=.TRUE.
ENDDO
IF (clean_fact) fact_ngeo=1

RETURN
END SUBROUTINE clean_ngeo

