!
! Copyright (C) 2024 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE print_gpu_memory()
!-------------------------------------------------------------------------
#if defined(__CUDA)
USE cudafor
IMPLICIT NONE
INTEGER :: ierr
INTEGER(KIND=cuda_count_kind) :: total, free

    ierr=cudaMemGetInfo(free, total)
    WRITE(6,'(5x,"used= ", f12.2, " Mb, free=", i20, &
         & " total=", f12.2, " GB")') (total-free)/1.D6, free, total/1.D9
#endif
RETURN
END SUBROUTINE print_gpu_memory
