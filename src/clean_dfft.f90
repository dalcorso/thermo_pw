!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE clean_dfft()
!----------------------------------------------------------------------
!
USE fft_base, ONLY : dfftp, dffts
USE input_parameters, ONLY : nr1, nr2, nr3, nr1s, nr2s, nr3s

IMPLICIT NONE

dfftp%nr1=nr1
dfftp%nr2=nr2
dfftp%nr3=nr3
dffts%nr1=nr1s
dffts%nr2=nr2s
dffts%nr3=nr3s

RETURN
END SUBROUTINE clean_dfft

