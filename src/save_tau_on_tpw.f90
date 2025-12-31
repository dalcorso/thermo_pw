!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE save_tau_on_tpw(igeom)
!-------------------------------------------------------------------------
!
!  This routine saves the atomic positions read from the dynamical matrices
!  in the thermo_pw variables. It also uses the contraint to compute the
!  values of the internal parameters. It assumes that the variables to
!  be saved are in the Quantum ESPRESSO variables.
!
USE kinds, ONLY : DP
USE thermo_mod, ONLY : uint_geo_eos, tau_geo_eos
USE ions_base,  ONLY : tau, nat
USE cell_base,  ONLY : celldm, ibrav, at
USE control_atomic_pos, ONLY : nint_var, iconstr_internal
IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom
REAL(DP) :: uint(nint_var)
!
!  get the constraint from the atomic positions
!
CALL internal_to_tau(ibrav, celldm, tau, at, uint, nat, nint_var, &
                     iconstr_internal, 2)
!
!  and copy it in the variables of thermo_pw
!
uint_geo_eos(1:nint_var,igeom)=uint(1:nint_var)
!
!  copy also the atomic positions 
!
tau_geo_eos(1:3,1:nat,igeom)=tau(1:3,1:nat)
!
RETURN
END SUBROUTINE save_tau_on_tpw
