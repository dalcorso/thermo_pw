!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE clean_dfft()

USE fft_base, ONLY : dfftp, dffts

IMPLICIT NONE

dfftp%nr1=0
dfftp%nr2=0
dfftp%nr3=0
dffts%nr1=0
dffts%nr2=0
dffts%nr3=0

RETURN
END SUBROUTINE clean_dfft

