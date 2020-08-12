!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE manage_mixing(dvscfout, dvscfin, dbecsum, mixin, npe, iter, kter, &
                                                            dr2, convt )
!
!   This routines is called by all the solve_xx routines to mix the
!   induced potential and possibly the derivatives of the partial
!   occupations in the PAW case.
!
USE kinds,             ONLY : DP
USE ions_base,         ONLY : nat
USE paw_variables,     ONLY : okpaw
USE fft_base,          ONLY : dfftp
USE noncollin_module,  ONLY : noncolin, npol, nspin_mag
USE uspp_param,        ONLY : nhm
USE control_ph,        ONLY : nmix_ph, tr2_ph, alpha_mix, flmixdpot

IMPLICIT NONE
INTEGER :: npe, iter, kter

COMPLEX(DP) :: dvscfout (dfftp%nnr, nspin_mag, npe)
COMPLEX(DP) :: dvscfin (dfftp%nnr, nspin_mag, npe)
COMPLEX(DP) :: dbecsum ((nhm * (nhm + 1))/2, nat, nspin_mag, npe)
COMPLEX(DP) :: mixin(dfftp%nnr*nspin_mag*npe+(nhm*(nhm+1)*nat*nspin_mag*npe)/2)
REAL(DP) :: dr2
LOGICAL :: convt

INTEGER :: ndim
COMPLEX(DP), ALLOCATABLE :: mixout(:)

IF (okpaw) THEN
   ALLOCATE (mixout(dfftp%nnr*nspin_mag*npe+(nhm*(nhm+1)*nat*nspin_mag*npe)/2))
   !
   !  In this case we mix also dbecsum
   !
   CALL setmixout(npe*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*npe)/2, &
                    mixout, dvscfout, dbecsum, ndim, -1 )
   CALL mix_potential (2*npe*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, & 
               alpha_mix(kter), dr2, npe*tr2_ph/npol, iter, nmix_ph,  &
                                                      flmixdpot, convt)
   CALL setmixout(npe*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*npe)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )
   !
   DEALLOCATE(mixout)
ELSE
   CALL mix_potential (2*npe*dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
                   alpha_mix(kter), dr2, npe*tr2_ph/npol, iter, &
                   nmix_ph, flmixdpot, convt)
ENDIF

RETURN
END SUBROUTINE manage_mixing

