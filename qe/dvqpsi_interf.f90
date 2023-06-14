!
! Copyright (C) 2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
INTERFACE  dvqpsi_interf
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE setpsi_zero(nbndk, st, npol, psicr, nbnd, &
                                          nnrs, nk, npe, nsolv)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)

  END SUBROUTINE setpsi_zero

! -----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE setpsi_evc(nbndk, st, npol, psicr, nbnd, nnr, &
                                         nk, npe, nsolv, ikb)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,  ONLY : DP

  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnr, nbnd, ikb
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnr, npol, nbnd*nk*npe*nsolv)

END SUBROUTINE setpsi_evc
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvlocpsi_dev(nbndk, st, npol, psicr, dpsicr, &
                                           dvloc, nbnd, nnrs, nk, npe, nsolv)
!----------------------------------------------------------------------- 
  !  
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  COMPLEX(DP), DEVICE, INTENT(IN) :: dvloc(nnrs, npe)

  END SUBROUTINE dvlocpsi_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvscfinpsi_gpu(nbndk, st, npol, psicr, dpsicr, &
                   dvscfins, nspin_mag, nbnd, nnrs, nk, npe, nsolv, ikb)
!----------------------------------------------------------------------- 
  !  
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd, nspin_mag, ikb
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  COMPLEX(DP), DEVICE, INTENT(IN) :: dvscfins(nnrs, nspin_mag, npe)

  END SUBROUTINE dvscfinpsi_gpu

!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE set_dvpsik_dev(lda, npw, nbndk, st, ikt,  &
                npol, dvpsik, psicr, nbnd, nnrs, nk, npe, nsolv) 
!----------------------------------------------------------------------- 
!  
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd
  INTEGER, INTENT(IN), DEVICE :: ikt(nk)
  INTEGER, INTENT(IN), DEVICE :: npw(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: dvpsik(lda*npol, nbnd*nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(IN) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !
  END SUBROUTINE set_dvpsik_dev
!
!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE adddvscf_dev(lda, npw, nbndk, st, ikt,  &
                npol, dvpsik, psicr, nbnd, nnrs, nk, npe, nsolv) 
!----------------------------------------------------------------------- 
!  
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  INTEGER, INTENT(IN), VALUE :: nk, nnrs, nbnd, npe, nsolv
  INTEGER, INTENT(IN), DEVICE :: ikt(nk)
  INTEGER, INTENT(IN), DEVICE :: npw(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  INTEGER, INTENT(IN), VALUE :: npol
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: dvpsik(lda*npol, nbnd*nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(IN) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !
  END SUBROUTINE adddvscf_dev

!------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvqpsi_us_dev0(ikt, st, npol, current_ikb_ph, &
                                      nbnd, nk, npe, nsolv, imode0, nspin)
!------------------------------------------------------------------------
!
USE cudafor
USE kinds, ONLY : DP

IMPLICIT NONE

INTEGER, VALUE :: nk, npe, nsolv, nbnd, current_ikb_ph, imode0, nspin
INTEGER, DEVICE :: ikt(nk)
INTEGER, DEVICE :: st(nk*npe*nsolv)
INTEGER, VALUE  :: npol
END SUBROUTINE dvqpsi_us_dev0
!
!------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvqpsi_us_dev1(ndmx, st, ikt, dvpsi,      &
                                         npol, nkb, nbnd, nk, npe, nsolv)
!------------------------------------------------------------------------
!
USE cudafor
USE kinds, ONLY : DP

IMPLICIT NONE

INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, nkb

INTEGER, DEVICE :: st(nk*npe*nsolv)
INTEGER, DEVICE :: ikt(nk)

COMPLEX(DP), DEVICE :: dvpsi(ndmx*npol,nbnd*nk*npe*nsolv)

END SUBROUTINE dvqpsi_us_dev1
!
!------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvqpsi_us_dev2(ndmx, st, ikt, dvpsi, npol, &
                                             nkb, nbnd, nk, npe, nsolv)
!------------------------------------------------------------------------
!
USE cudafor
USE kinds, ONLY : DP

IMPLICIT NONE

INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, nkb

INTEGER, DEVICE :: st(nk*npe*nsolv)
INTEGER, DEVICE :: ikt(nk)

COMPLEX(DP), DEVICE :: dvpsi(ndmx*npol,nbnd*nk*npe*nsolv)
END SUBROUTINE dvqpsi_us_dev2
!
END INTERFACE dvqpsi_interf
#endif
