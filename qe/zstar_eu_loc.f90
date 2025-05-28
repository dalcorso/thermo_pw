!
! Copyright (C) 2017 Dal Corso Andrea 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE zstar_eu_loc (drhoscf, zstareu0)
  !-----------------------------------------------------------------------
  !    This subroutine computes the contribution of the local
  !    potential to the electronic term <psi|dv-e ds|dpsi> of the effective 
  !    charges. It works in all cases, NC, US and PAW.
  !    On output the effective charges are distributed among processors
  !    and pools.
  !
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat
  USE fft_base,  ONLY : dffts
  USE cell_base, ONLY : omega
  USE modes,     ONLY : u
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN) :: drhoscf (dffts%nnr, nspin_mag, 3)
  COMPLEX(DP), INTENT(INOUT) :: zstareu0 (3, 3 * nat)
  ! the change of density due to perturbations
  ! auxiliary array where the result is added

  INTEGER :: ipert, is, nu_j
  ! counter on perturbations
  ! counter on spin polarizations
  ! counter on the j modes

  COMPLEX(DP), EXTERNAL :: zdotc
  COMPLEX(DP), ALLOCATABLE :: dvloc (:)
  ! d Vloc / dtau

  ALLOCATE (dvloc( dffts%nnr))    
  !
  DO nu_j = 1, 3 * nat
     CALL compute_dvloc (u(1,nu_j), .FALSE., dvloc)
     DO ipert = 1, 3
        DO is = 1, nspin_lsda
           zstareu0 (ipert, nu_j) = zstareu0 (ipert, nu_j) - &
                zdotc (dffts%nnr, drhoscf (1, is, ipert), 1, dvloc, 1) * &
                  omega / (dffts%nr1 * dffts%nr2 * dffts%nr3)
        ENDDO
     ENDDO
  ENDDO

  DEALLOCATE(dvloc)
  RETURN
END SUBROUTINE zstar_eu_loc
