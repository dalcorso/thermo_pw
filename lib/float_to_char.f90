!
! Copyright (C) 2015-2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
FUNCTION float_to_char(a,n)
!
!  This function transforms a float number into a character string.
!  n gives the number of digits of the decimal part.
!  The length of the string is of 8 characters, so the number cannot
!  be larger than 999999. (one place is reserved for the sign).
!  n is reduced if the string has not sufficient space
!
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=8) :: float_to_char
REAL(DP), INTENT(IN) :: a
INTEGER, INTENT(IN) :: n

CHARACTER(LEN=8) :: buffer
INTEGER :: i, ai, neff, leng
CHARACTER(LEN=6) :: as, fas, int_to_char
REAL(DP) :: aeff, aux
LOGICAL  :: lsign
!
!   If negative set lsign to .TRUE.
!
lsign=.FALSE.
aux=a
IF (a<0.0_DP) THEN
   lsign=.TRUE.
   aux=-a
ENDIF

IF (aux > 999999.) CALL errore('float_to_char','float too large',1)
!
!  find the effective number of decimal digits that can be written in the
!  string
!
DO i = 0,n
   IF (aux*10**i <= 9999999.) neff=i
ENDDO
!
!  round a to that number of digits
!
aeff=NINT(aux*10**neff)/10.0_DP**neff
!
!  and transform it to a string, first converting the integer part and
!  then the decimal part
!
ai=INT(aeff)
as=int_to_char(ai)
ai=NINT((aeff-ai)*10**neff)
fas=int_to_char(ai)
!
!  if the length of the decimal part string has less than neff digits,
!  there are leading 0 that must be inserted
!
leng=LEN(TRIM(fas))
DO i=1,neff-leng
   fas='0'//TRIM(fas)
ENDDO
buffer=TRIM(as)
!
!  finally join the two strings
!
IF (neff>0) buffer=TRIM(as)//'.'//TRIM(fas)
!
!  and put the negative sign if necessary
!
IF (lsign) buffer='-'//TRIM(buffer)
float_to_char=TRIM(buffer)

RETURN
END FUNCTION float_to_char
