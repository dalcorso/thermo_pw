!
! Copyright (C) 2023 Andrea Dal Corso and Xuejun Gong
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE diago_dev_interf
   ATTRIBUTES(GLOBAL) SUBROUTINE diago_dev(ntime,nsizex,nsizek,nvec,outk,&
                       a, b, w, x, work, lwork, adiag, bdiag, rwork, lrwork, &
                       iwork, liwork, info, ifail, m)
   USE cudafor
   USE kinds,   ONLY : DP
   IMPLICIT NONE

   INTEGER, VALUE :: nsizex
   INTEGER, VALUE :: ntime
   INTEGER, VALUE :: lrwork
   INTEGER, VALUE :: liwork
   INTEGER, VALUE :: lwork
   INTEGER, VALUE :: nvec
   LOGICAL, DEVICE :: outk(ntime)
   INTEGER, DEVICE :: nsizek(ntime)
   INTEGER, DEVICE :: info(ntime)
   INTEGER, DEVICE :: ifail(nsizex,ntime)
   INTEGER, DEVICE :: iwork(liwork,ntime)
   INTEGER, DEVICE :: m(ntime)
   COMPLEX(DP), DEVICE :: a(nsizex,nsizex,ntime), b(nsizex,nsizex,ntime), &
                          x(nsizex,nsizex,ntime), work(lwork,ntime)
   REAL(DP), DEVICE :: rwork(lrwork,ntime)
   REAL(DP), DEVICE :: w(nsizex,ntime), adiag(nsizex,ntime), &
                                        bdiag(nsizex,ntime)

   END SUBROUTINE diago_dev
END INTERFACE diago_dev_interf
#endif
