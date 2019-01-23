!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE apply_trev_r(psic)
!
!  This routine applies the time reversal operator to the wavefunctions
!  evc at the k point ikk_evc and puts the output in evc with the order
!  of G vectors of ikk_tevc
!
USE kinds,     ONLY : DP
USE wvfct,     ONLY : nbnd
USE fft_base,  ONLY:  dfftp
USE noncollin_module,  ONLY : npol

IMPLICIT NONE
COMPLEX(DP) :: psic(dfftp%nnr,npol,nbnd)
COMPLEX(DP), ALLOCATABLE :: aux2(:,:)

INTEGER :: ibnd, npt

ALLOCATE(aux2(dfftp%nnr,npol))
npt=dfftp%nnr

DO ibnd=1,nbnd
   aux2 (1:npt, 1) = - psic(1:npt, 2, ibnd)
   aux2 (1:npt, 2) = psic(1:npt, 1, ibnd)
   aux2=CONJG(aux2)
   psic(:,:,ibnd)=aux2(:,:)
ENDDO

DEALLOCATE(aux2)

RETURN
END SUBROUTINE apply_trev_r
