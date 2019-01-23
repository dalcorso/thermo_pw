!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE add_dvscf_rhs( dvscfins, isolv, ipert, ik, npe )
!
! This routine calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
! dvscf_q is given in input on the smooth mesh
!
! In the US/PAW case it assumes also tha int3 integral are calculated
! and that vkb contains the beta functions. It applies also the term with 
! int3.
!
! The routine is task group aware. 
! isolv=1 means that the potential is the one given in input
! isolv=2 reverses the sign of the magnetic field 
!
! it works only for the perturbation ipert and the k point ik.
! npe is the number of perturbations (dimension of dvscfins)
!
! The result is added to dvpsi (note that this array is not set to zero)
!
USE kinds,          ONLY : DP
USE fft_base,       ONLY : dffts
USE wvfct,          ONLY : nbnd, npwx
USE wavefunctions_module, ONLY : evc
USE noncollin_module,  ONLY : noncolin, npol, nspin_mag
USE control_lr,     ONLY : nbnd_occ
USE spin_orb,       ONLY : domag
USE lsda_mod,       ONLY : lsda, current_spin, isk
USE uspp,           ONLY : okvan
USE lrus,           ONLY : int3_nc, becp1
USE qpoint_aux,     ONLY : becpt 
USE nc_mag_aux,     ONLY : int3_save
USE fft_helper_subroutines, ONLY : fftx_ntgrp
USE qpoint,         ONLY : ikks 
USE eqv,            ONLY : dvpsi

IMPLICIT NONE

INTEGER :: isolv, ipert, ik, npe
COMPLEX(DP) :: dvscfins(dffts%nnr, nspin_mag, npe)

INTEGER :: ibnd, incr, ipol, ikk, v_siz
COMPLEX(DP), ALLOCATABLE :: tg_dv(:,:), tg_psic(:,:), aux1(:,:), aux2(:,:)
 
CALL start_clock ('vpsifft')
!
!  allocate space for task groups if necessary
!
incr=1
IF ( dffts%has_task_groups ) THEN
   !
   v_siz =  dffts%nnr_tg
   ALLOCATE( tg_dv( v_siz, nspin_mag ) )
   ALLOCATE( tg_psic( v_siz, npol ) )
   incr = fftx_ntgrp(dffts)
   !
ENDIF
!
!  allocate axiliary space and compute the needed index
!
ALLOCATE (aux1(dffts%nnr, npol))
ALLOCATE (aux2(npwx*npol, nbnd))
ikk = ikks(ik)
IF (lsda) current_spin = isk (ikk)
!
!  change the sign of the magnetic field if required
!
IF (isolv==2) THEN
   dvscfins(:,2:4,ipert)=-dvscfins(:,2:4,ipert)
   IF (okvan) int3_nc(:,:,:,:,ipert)=int3_save(:,:,:,:,ipert,2)
ENDIF
!
!  Set the potential for task groups
!
IF (dffts%has_task_groups) THEN
   IF (noncolin) THEN
      CALL tg_cgather(dffts, dvscfins(:,1,ipert), tg_dv(:,1))
      IF (domag) THEN
         DO ipol=2,4
            CALL tg_cgather(dffts, dvscfins(:,ipol,ipert), tg_dv(:,ipol))
         ENDDO
      ENDIF
   ELSE
      CALL tg_cgather(dffts, dvscfins(:,current_spin,ipert), tg_dv(:,1))
   ENDIF
ENDIF
!
!  apply the potential
!
aux2=(0.0_DP,0.0_DP)
DO ibnd = 1, nbnd_occ (ikk), incr
   IF (dffts%has_task_groups) THEN
      CALL cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
      CALL apply_dpot(v_siz, tg_psic, tg_dv, 1)
      CALL cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, nbnd_occ(ikk))
   ELSE
      CALL cft_wave (ik, evc (1, ibnd), aux1, +1)
      CALL apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipert), current_spin)
      CALL cft_wave (ik, aux2 (1, ibnd), aux1, -1)
   ENDIF
ENDDO
!
!  add the computed term to dvpsi
!
CALL ZAXPY(npwx*npol*nbnd, (1.0_DP,0.0_DP), aux2, 1, dvpsi, 1)

CALL stop_clock ('vpsifft')
!
!  In the case of US pseudopotentials there is an additional
!  selfconsist term which comes from the dependence of D on
!  V_{eff} on the bare change of the potential
!
IF (isolv==1) THEN
   CALL adddvscf_tpw (ipert, ik, becp1)
ELSE
   CALL adddvscf_tpw (ipert, ik, becpt)
ENDIF
!
!  reset the original magnetic field if it was changed
!
IF (isolv==2) THEN
   dvscfins(:,2:4,ipert)=-dvscfins(:,2:4,ipert)
   IF (okvan) int3_nc(:,:,:,:,ipert)=int3_save(:,:,:,:,ipert,1)
ENDIF

DEALLOCATE(aux2)
DEALLOCATE(aux1)
!
! deallocate task group memory
!
IF ( dffts%has_task_groups ) THEN
   DEALLOCATE( tg_psic )
   DEALLOCATE( tg_dv )
ENDIF

RETURN
END SUBROUTINE add_dvscf_rhs
