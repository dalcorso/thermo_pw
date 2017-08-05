!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE psym_dmageq (dvtosym)
  !-----------------------------------------------------------------------
  !
  ! ...  p-symmetrize the magnetization change due to an electric field.
  !
  USE kinds,     ONLY : DP
  USE lsda_mod,   ONLY : nspin
  USE mp_bands,  ONLY : me_bgrp
  USE fft_base,  ONLY : dfftp
  USE scatter_mod, ONLY : cgather_sym
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym (dfftp%nnr, nspin)
    ! the potential to symmetrize
    !-local variable
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, ir3, ioff, ioff_tg, nxyp

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:)
  ! the potential to symm

  CALL start_clock ('psym_dmageq')

  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin))
  DO is = 1, nspin
     CALL cgather_sym (dfftp, dvtosym (:, is), ddvtosym (:, is) )
  ENDDO

  CALL sym_dmageq (ddvtosym)
  nxyp = dfftp%nr1x * dfftp%my_nr2p
  DO is = 1, nspin
     DO ir3 = 1, dfftp%my_nr3p
        ioff    = dfftp%nr1x * dfftp%my_nr2p * (ir3-1)
        ioff_tg = dfftp%nr1x * dfftp%nr2x* (dfftp%my_i0r3p+ir3-1) + &
                               dfftp%nr1x * dfftp%my_i0r2p
        CALL zcopy (nxyp, ddvtosym (ioff_tg+1, is), 1, dvtosym (ioff+1, is), 1)
     ENDDO
  ENDDO
  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psym_dmageq')

#else

  CALL sym_dmageq (dvtosym)

#endif

  RETURN

END SUBROUTINE psym_dmageq
