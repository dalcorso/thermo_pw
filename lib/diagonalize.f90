!
! Copyright (C) 2019 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE diagonalize
!
!   This module contains the routines to diagonalize a matrix
!   not available in QE.
!
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: diagonalize_r

CONTAINS

!-------------------------------------------------------------------------
SUBROUTINE diagonalize_r (npw,nbnd,h,e,v)
!-------------------------------------------------------------------------
!
! f90 interface to LAPACK routine ZHEEVX which calculates
! nbnd eigenvalues and eigenvectors of a complex hermitean matrix h. 
! 
IMPLICIT NONE

INTEGER, INTENT(IN) :: npw
INTEGER, INTENT(IN) :: nbnd

REAL(DP), INTENT(INOUT) :: h(npw,npw)   ! matrix to be diagonalized
REAL(DP), INTENT(OUT)   :: e(nbnd)      ! eigenvalues
REAL(DP), INTENT(OUT)   :: v(npw,nbnd)  ! eigenvectors (column-wise)

INTEGER  :: lwork,  &! auxiliary variable
            info,   &! flag saying if LAPACK execution was OK
            m        ! number of eigenvalues

CHARACTER(LEN=1) :: jobz, range, uplo  ! select the task in LAPACK

REAL(DP), ALLOCATABLE    :: work(:)      ! as above
INTEGER, ALLOCATABLE     :: iwork(:)     !    "
INTEGER, ALLOCATABLE     :: ifail(:)     !    "
REAL(DP)                 :: rdummy, zero ! dummy variable, zero
REAL(DP), ALLOCATABLE    :: ee(:)        ! axiliary space for eigenvalues
!
!   Initialize flags
!
jobz  = 'V' ! compute eigenvalues and eigenvectors
uplo  = 'U' ! LAPACK routines use the upper triangle of the input matrix
range = 'I' ! compute bands from 1 to nbnd

zero = 0.0_DP
v(:,:) = 0.0_DP
!
! allocate arrays of known size
!
ALLOCATE( ee(npw) )
ALLOCATE( iwork(5*npw) )
ALLOCATE( ifail(npw) )
ALLOCATE( work(16*npw) )
lwork=16*npw
!
! and diagonalize the matrix
!
CALL dsyevx(jobz, range, uplo, npw, h, npw, rdummy, rdummy, 1, nbnd, zero, &
            m, ee, v, npw, work, lwork, iwork, ifail, info)
!
IF (ABS(info) /= 0) THEN
   WRITE(stdout,'("Error in the diagonalization, info= ", i5)') info
   STOP 1
ENDIF
!
!
!  NB: h is overwritten by this routine. We save the eigenvalues on e
!      for plotting
!
e(1:nbnd)=ee(1:nbnd) 
       
DEALLOCATE(work)
DEALLOCATE(iwork)
DEALLOCATE(ifail)
DEALLOCATE(ee)

RETURN
END SUBROUTINE diagonalize_r

END MODULE diagonalize
