!
! Copyright (C) 2018 Andrea Dal Corso 
! This routine has been obtained by modifyng the routine solve_e 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_int3q(irr, imode0, rpert, drhoscf, int3_paw0, dvscfin)
  !-----------------------------------------------------------------------
  !
  !  This routine computes the parts of the right hand side of the 
  !  linear system when the perturbed wavefunctions dpsi vanish. 
  !  In addition to -P^+_c dV_ext/d u computed elsewhere there is
  !  a nonvanishing term with NLCC or in the US or PAW case. 
  !  In output drhoscf contains the change of the charge, dvscfin the 
  !  change of the Hxc potential.
  !  The int3 integral due to this induced potential is computed here.
  ! 
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat
  USE qpoint,           ONLY : xq
  USE gvecs,            ONLY : doublegrid
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : fft_interpolate
  USE lr_symm_base,     ONLY : irotmq, minus_q, nsymq, rtau
  USE modes,            ONLY : npertx, u, t, tmq
  USE buffers,          ONLY : get_buffer
  USE units_ph,         ONLY : iudrhous, lrdrhous
  USE control_ph,       ONLY : lgamma_gamma
  USE noncollin_module, ONLY : noncolin, nspin_mag, domag
  USE scf,              ONLY : rho
  USE uspp,             ONLY : okvan, nlcc_any
  USE uspp_param,       ONLY : nhm
  USE paw_variables,    ONLY : okpaw
  USE phus,             ONLY : becsumort
  USE lrus,             ONLY : int3_paw
  USE paw_onecenter,    ONLY : paw_dpotential
  USE paw_symmetry,     ONLY : paw_dusymmetrize, paw_dumqsymmetrize
  USE dv_of_drho_lr,    ONLY : dv_of_drho
  USE mp,               ONLY : mp_sum
  USE mp_pools,         ONLY : inter_pool_comm

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: irr, imode0, rpert
  COMPLEX(DP), INTENT(INOUT) :: drhoscf(dfftp%nnr, nspin_mag, rpert)
  COMPLEX(DP), INTENT(INOUT) :: dvscfin(dfftp%nnr, nspin_mag, rpert)
  COMPLEX(DP), INTENT(INOUT) :: int3_paw0 (nhm, nhm, nat, nspin_mag, rpert)

  COMPLEX(DP), ALLOCATABLE ::          &
                   dbecsum(:,:,:,:),   &  ! the change of becsum
                   drhoc(:)            ! The change of the core charge

  INTEGER :: ipol, mu, is

  CALL start_clock ('set_int3q')

  ALLOCATE (drhoc(dfftp%nnr))
  IF (okpaw) ALLOCATE(dbecsum((nhm*(nhm+1))/2, nat, nspin_mag, rpert))
  !
  !  Put on drhoscf the change of the augmentation charge computed
  !  in drho. In the paw case we need becsumort.
  !
  do ipol = 1, rpert
     mu = imode0 + ipol
     IF (okvan) THEN
        CALL get_buffer (drhoscf(1,1,ipol), lrdrhous, iudrhous, mu)
     ELSE
        drhoscf(:,:,ipol)=(0.0_DP,0.0_DP)
     ENDIF
     IF (okpaw) dbecsum(:,:,:,ipol)=becsumort(:,:,:,mu)
  ENDDO

  CALL mp_sum ( drhoscf, inter_pool_comm )
  !
  !   Symmetrize the charge
  !
  IF (.NOT.lgamma_gamma) THEN
     CALL psymdvscf (rpert, irr, drhoscf)
     IF ( noncolin.and.domag ) CALL psym_dmag( rpert, irr, drhoscf)
     IF (okpaw) THEN
        IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,rpert,irr, &
                                             npertx,irotmq,rtau,xq,tmq)
        CALL PAW_dusymmetrize(dbecsum,rpert,irr,npertx,nsymq,rtau,xq,t)
     ENDIF
  ENDIF
  !
  !   calculate the corresponding linear potential response
  !
  dvscfin = drhoscf
  DO ipol=1,rpert
     IF (imode0+ipol > 0) THEN
        CALL addcore (imode0+ipol, drhoc)
     ELSE
        drhoc(:) = (0.0_DP,0.0_DP)
     ENDIF
     CALL dv_of_drho (dvscfin (1, 1, ipol), .TRUE., drhoc)
  ENDDO
!
!   In the PAW case computes the change of the D coefficients
!
   IF (okpaw) THEN
      CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,rpert)
      int3_paw0=int3_paw
   ENDIF
!
!   In the US and PAW case computes the integral of dV_Hxc and the 
!   augmentation function. This quantity is needed in adddvscf.
!
  CALL newdq(dvscfin,rpert)

  IF (okpaw) DEALLOCATE(dbecsum)
  DEALLOCATE(drhoc)

  CALL stop_clock ('set_int3q')

  RETURN
END SUBROUTINE set_int3q
