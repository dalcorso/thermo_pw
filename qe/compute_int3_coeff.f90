!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_int3_coeff(dvscfin, dbecsum, npe)

USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE fft_base,         ONLY : dfftp
USE uspp_param,       ONLY : nhm
USE scf,              ONLY : rho
USE noncollin_module, ONLY : noncolin, nspin_mag, domag
USE uspp,             ONLY : okvan
USE paw_variables,    ONLY : okpaw
USE paw_onecenter,    ONLY : paw_dpotential
USE lrus,             ONLY : int3, int3_paw, int3_nc
USE nc_mag_aux,       ONLY : int3_save

IMPLICIT NONE

INTEGER :: npe
COMPLEX(DP) :: dvscfin ( dfftp%nnr, nspin_mag, npe)
COMPLEX(DP) :: dbecsum ( (nhm * (nhm + 1))/2, nat, nspin_mag, npe)

INTEGER :: ipert

IF (.NOT.okvan) RETURN

IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
!
!   with the new change of the potential we compute the integrals
!   of the change of potential and Q
!
CALL newdq (dvscfin, npe)
!
!  In the noncollinear magnetic case computes the int3 coefficients with
!  the opposite sign of the magnetic field. They are saved in int3_save,
!  that must have been allocated by the calling routine 
!
IF (noncolin.AND.domag) THEN
   int3_save(:,:,:,:,:,1)=int3_nc(:,:,:,:,:)
   IF (okpaw) rho%bec(:,:,2:4)=-rho%bec(:,:,2:4)
   DO ipert=1,npe
      dvscfin(:,2:4,ipert)=-dvscfin(:,2:4,ipert)
      IF (okpaw) dbecsum(:,:,2:4,ipert)=-dbecsum(:,:,2:4,ipert)
   ENDDO
!
!   if needed recompute the paw coeffients with the opposite sign of
!   the magnetic field
!
   IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
!
!   here compute the int3 integrals
!
   CALL newdq (dvscfin, npe)
   int3_save(:,:,:,:,:,2)=int3_nc(:,:,:,:,:)
!
!  restore the correct sign of the magnetic field.
!
   DO ipert=1,npe
      dvscfin(:,2:4,ipert)=-dvscfin(:,2:4,ipert)
      IF (okpaw) dbecsum(:,:,2:4,ipert)=-dbecsum(:,:,2:4,ipert)
   ENDDO
   IF (okpaw) rho%bec(:,:,2:4)=-rho%bec(:,:,2:4)
!
!  put into int3_nc the coefficient with +B
!
   int3_nc(:,:,:,:,:)=int3_save(:,:,:,:,:,1)
ENDIF

RETURN
END SUBROUTINE compute_int3_coeff
