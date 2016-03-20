! Copyright (C) 2016-present Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
FUNCTION scal_prod(ndx, ndim, psi1, psi2)
!
!  This function computes the scalar product between two wavefunctions.
!  It knows if the calculation has spinors or scalar wavefunctions, 
!  if the functions are distributed among several processors or if 
!  the calculation is gamma_only and it acts differently on the 
!  wavefunctions depending on these settings. 
!  It must be called by all processors that belong to the same band group.
!
USE kinds, ONLY : DP
USE control_flags,      ONLY : gamma_only
USE noncollin_module,   ONLY : noncolin, npol
USE gvect,              ONLY : gstart

USE mp, ONLY : mp_sum
USE mp_bands, ONLY : intra_bgrp_comm

IMPLICIT NONE
COMPLEX(DP) :: scal_prod
INTEGER :: ndx, ndim

COMPLEX(DP) :: psi1(ndx), psi2(ndx)
COMPLEX(DP) :: aux
REAL(DP) :: ddot
COMPLEX(DP) :: zdotc

INTEGER :: ndmx

IF (noncolin) THEN
   ndmx = ndx / npol
ELSE
   ndmx = ndx
ENDIF

IF (gamma_only) THEN
   aux=2.0_DP*ddot(2*ndmx*npol,psi1,1,psi2,1)
   IF  (gstart==2) aux=aux-DBLE(psi1(1))*DBLE(psi2(1))
ELSE
  IF (noncolin) THEN
     aux= ZDOTC (ndim, psi1, 1, psi2, 1)
     aux= aux+ ZDOTC (ndim, psi1(ndmx+1), 1, psi2(ndmx+1), 1)
  ELSE
     aux= ZDOTC (ndim, psi1, 1, psi2, 1)
  ENDIF
ENDIF
CALL mp_sum(aux, intra_bgrp_comm )

scal_prod = aux
RETURN
END FUNCTION scal_prod
