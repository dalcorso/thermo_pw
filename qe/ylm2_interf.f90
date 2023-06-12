INTERFACE ylm2_interf
 ATTRIBUTES(GLOBAL) SUBROUTINE ylm2_dev(ylm_d, lmaxkb, npwx, nk, ikt)
    !----------------------------------------------------------------------
    !
    USE cudafor
    USE kinds,        ONLY : DP

    IMPLICIT NONE

    INTEGER, VALUE, INTENT(IN) :: lmaxkb, npwx, nk
    INTEGER, INTENT(IN), DEVICE :: ikt(nk)
    REAL(DP), DEVICE, INTENT(OUT) :: ylm_d((lmaxkb + 1) **2, npwx, nk)

END SUBROUTINE ylm2_dev
END INTERFACE ylm2_interf
