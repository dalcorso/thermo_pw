!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE ev_mod

   USE kinds,  ONLY : DP
   USE constants, ONLY : ry_kbar
   !
   !   This module contains some auxiliary routines to compute the
   !   equation of state. These routine were originally on the ev.f90
   !   routine of QE.
   !
   SAVE
   PRIVATE 

   INTEGER, PARAMETER :: nmaxpt=100 
   INTEGER :: npt             ! number of points 
   INTEGER :: ieos            ! index of the equation of state
   REAL(DP) :: emin           ! the minimum energy
   REAL(DP) :: v0(nmaxpt), efit(nmaxpt), etot(nmaxpt)

   PUBLIC v0, etot, emin, npt, eosdiff, eqstate, &
          find_minimum, initialize_data_ev

CONTAINS
!-----------------------------------------------------------------------
   SUBROUTINE initialize_data_ev(v0_, etot_, npt_, ieos_)
!-----------------------------------------------------------------------
!
!  This routine copy the external data into the variables of the
!  module. It is used to request the equation of state, to set
!  the number of data points npt and for each point to set the
!  volume and the energy.
!
   IMPLICIT NONE
   INTEGER :: npt_, ieos_
   REAL(DP), INTENT(IN) :: v0_(npt_), etot_(npt_)

   ieos=ieos_
   npt=npt_
   v0(1:npt)=v0_(1:npt_) 
   etot(1:npt)=etot_(1:npt_) 

   RETURN
   END SUBROUTINE initialize_data_ev
!
!-----------------------------------------------------------------------
   SUBROUTINE eqstate(npar,par,chisq,ediff)
!-----------------------------------------------------------------------
!
!  This routine receives the parameters par and computes the chisq
!  and the differences between the energy computed from the equation
!  of state and that set in etot by initialize_data_ev.
!
   USE eos,  ONLY : eos_energy
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: npar
   REAL(DP), INTENT(IN) :: par(npar)
   REAL(DP), INTENT(OUT):: chisq
   REAL(DP), OPTIONAL, INTENT(OUT):: ediff(npt)
   INTEGER :: i
   REAL(DP) :: k0, dk0, d2k0, vol0
!
   vol0 = par(1)
   k0   = par(2)/ry_kbar ! converts k0 to Ry atomic units...
   dk0  = par(3)
   d2k0 = par(4)*ry_kbar ! and d2k0/dp2 to (Ry a.u.)^(-1)
!
   DO i=1,npt
      CALL eos_energy(ieos, v0(i), efit(i), vol0, k0, dk0, d2k0)
   ENDDO
!
!      emin = equilibrium energy obtained by minimizing chi**2
!
   emin = 0.0d0
   DO i = 1,npt
      emin = emin + etot(i)-efit(i)
   ENDDO
   emin = emin/npt
!
!  Computation of chisq
!
   chisq = 0.0d0
   DO i = 1,npt
      efit(i) = efit(i)+emin
      chisq   = chisq + (etot(i)-efit(i))**2
      IF (PRESENT(ediff)) ediff(i) = efit(i)-etot(i)
   ENDDO
   chisq = chisq/npt
!
   RETURN
   END SUBROUTINE eqstate
!
!----------------------------------------------------------------------------
    SUBROUTINE eosdiff(m_, n_, par_, f_, i_)
!----------------------------------------------------------------------------
!
!  This routine is a driver that call eqstate. The input parameters
!  of this routine match those required by the minpack routine lmdif0.
!
    IMPLICIT NONE
    INTEGER, INTENT(in)  :: m_, n_
    INTEGER, INTENT(inout)   :: i_
    REAL(DP),INTENT(in)    :: par_(n_)
    REAL(DP),INTENT(out)   :: f_(m_)
    REAL(DP) :: chisq_
       !
    CALL eqstate(n_,par_,chisq_, f_)

    RETURN
    END SUBROUTINE eosdiff
!
!-----------------------------------------------------------------------
    SUBROUTINE find_minimum(npar,par,chisq)
!-----------------------------------------------------------------------
!
!   This is a driver that calls the routine lmdif0 and minimizes the
!   chisq varying the parameters par. The routine eosdiff must
!   provide the energy differences given the parameters par.
!   In input the variable par must contain an estimate of their values.
!

    USE lmdif_module, ONLY : lmdif0
    IMPLICIT NONE
    INTEGER ,INTENT(IN)  :: npar
    REAL(DP),INTENT(OUT) :: par(npar)
    REAL(DP),INTENT(OUT) :: chisq
    !
    REAL(DP) :: ediff(npt)
    INTEGER :: info
    !      
    CALL lmdif0(eosdiff, npt, npar, par, ediff, 1.d-12, info)
    !
    IF (info>0 .AND. info<5) THEN
!        PRINT*, "Minimization succeeded"
    ELSEIF(info>=5) THEN
       CALL errore("find_minimum", "Minimization stopped before convergence",1)
    ELSEIF(info<=0) THEN 
       CALL errore("find_minimum", "Minimization error", 1)
    ENDIF
      !
    CALL eqstate(npar,par,chisq)

    RETURN
    END SUBROUTINE find_minimum

END MODULE ev_mod
