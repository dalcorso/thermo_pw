INTERFACE ylmr2_interf
   ATTRIBUTES(DEVICE) SUBROUTINE ylmr2_dev (lmax2, g, gg, ylm)
   !-----------------------------------------------------------------------
   USE upf_kinds, ONLY : DP
   !
   IMPLICIT NONE
   !
   INTEGER, VALUE, INTENT(IN) :: lmax2
   REAL(DP), VALUE, INTENT(IN) :: gg
   REAL(DP), INTENT(IN) :: g(3)
   !
   REAL(DP), DEVICE, INTENT(OUT) :: ylm (lmax2)
   !
   END SUBROUTINE ylmr2_dev
END INTERFACE ylmr2_interf
