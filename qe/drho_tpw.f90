!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

SUBROUTINE drho_tpw
  !-----------------------------------------------------------------------
  !
  !! Here we compute, for each mode the change of the charge density
  !! due to the displacement, at fixed wavefunctions. These terms
  !! are saved on disk. The orthogonality part is included in the
  !! computed change.
  !
  !
  USE kinds,      ONLY : DP
  USE gvecs,      ONLY : doublegrid
  USE fft_base,   ONLY : dfftp, dffts
  USE lsda_mod,   ONLY : nspin
  USE cell_base,  ONLY : omega
  USE ions_base,  ONLY : nat
  USE buffers,    ONLY : save_buffer
  USE noncollin_module, ONLY : noncolin, npol, nspin_lsda, nspin_mag
  USE uspp_param, ONLY : upf, nhm
  USE uspp,       ONLY : okvan, nkb
  USE wvfct,      ONLY : nbnd
  USE paw_variables,    ONLY : okpaw
  USE control_ph, ONLY : ldisp, all_done, rec_code_read

  USE lrus,       ONLY : becp1
  USE klist,      ONLY : lgauss
  USE two_chem,   ONLY : twochem
  USE qpoint,     ONLY : nksq
  USE control_lr, ONLY : lgamma

  USE dynmat,     ONLY : dyn00
  USE modes,      ONLY : npertx, npert, nirr, u
  USE phus,       ONLY : becsumort, alphap
  USE units_ph,   ONLY : lrdrhous, iudrhous

  USE mp_pools,   ONLY : inter_pool_comm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE becmod,     ONLY : bec_type, allocate_bec_type, deallocate_bec_type
  USE fft_interfaces, ONLY : fft_interpolate

  IMPLICIT NONE

  INTEGER :: mode, is, ir, irr, iper, npe, nrstot, nu_i, nu_j, ik, &
             ipol
  ! counter on modes
  ! counter on atoms and polarizations
  ! counter on atoms
  ! counter on spin
  ! counter on perturbations
  ! the number of points
  ! counter on modes
  ! counter on k-point
  ! counter on coordinates

  REAL(DP), ALLOCATABLE :: wgg (:,:,:)
  ! the weight of each point

  COMPLEX(DP) :: wdyn (3 * nat, 3 * nat)
  TYPE (bec_type), POINTER :: becq(:), alpq(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvlocin (:), drhous (:,:,:),&
       drhoust (:,:,:), dbecsum(:,:,:,:), dbecsum_nc(:,:,:,:,:)
  ! auxiliary to store bec at k+q
  ! auxiliary to store alphap at
  ! the change of the local potential
  ! the change of the charge density
  ! the change of the charge density
  ! the derivative

!
!  The PAW case requires dbecsumort so we recalculate this starting part
!  This will be changed soon
!
  IF (all_done) RETURN
  IF ((rec_code_read >=-20 .AND..NOT.okpaw)) RETURN

  dyn00(:,:) = (0.d0,0.d0)
  IF (.NOT.okvan) RETURN
  CALL start_clock ('drho')
  !
  !    first compute the terms needed for the change of the charge density
  !    due to the displacement of the augmentation charge
  !
  CALL compute_becsum_ph()
  if(twochem.and.lgamma.and.lgauss) call compute_becsum_ph_cond()
  !
  CALL compute_alphasum()
  if(twochem.and.lgamma.and.lgauss) call compute_alphasum_cond()
  !
  !    then compute the weights
  !
  call start_clock('compute_wgg')
  ALLOCATE (wgg (nbnd, nbnd, nksq))
  IF (lgamma) THEN
     becq => becp1
     alpq => alphap
  ELSE
     ALLOCATE (becq ( nksq))
     ALLOCATE (alpq ( 3, nksq))
     DO ik =1,nksq
        CALL allocate_bec_type (  nkb, nbnd, becq(ik))
        DO ipol=1,3
           CALL allocate_bec_type (  nkb, nbnd, alpq(ipol,ik))
        ENDDO
     ENDDO
  ENDIF
  CALL compute_weight (wgg)
  call stop_clock('compute_wgg')
  !
  !    becq and alpq are sufficient to compute the part of C^3 (See Eq. 37
  !    which does not contain the local potential
  !
  IF (.NOT.lgamma) CALL compute_becalp (becq, alpq)
  call start_clock('nldyntot')
  CALL compute_nldyn (dyn00, wgg, becq, alpq)
  call stop_clock('nldyntot')
  !
  !   now we compute the change of the charge density due to the change of
  !   the orthogonality constraint
  !
  ALLOCATE (drhous ( dfftp%nnr, nspin_mag , 3 * nat))
  ALLOCATE (dbecsum( nhm * (nhm + 1) /2, nat, nspin_mag, 3 * nat))
  dbecsum=(0.d0,0.d0)
  IF (noncolin) THEN
     ALLOCATE (dbecsum_nc( nhm, nhm, nat, nspin, 3 * nat))
     dbecsum_nc=(0.d0,0.d0)
     CALL compute_drhous_nc_tpw (drhous, dbecsum_nc, wgg, becq, alpq)
  ELSE
     CALL compute_drhous_tpw (drhous, dbecsum, wgg, becq, alpq)
  ENDIF

  IF (.NOT.lgamma) THEN
     DO ik=1,nksq
        CALL deallocate_bec_type(becq(ik))
        DO ipol=1,3
           CALL deallocate_bec_type(alpq(ipol,ik))
        ENDDO
     END DO
     DEALLOCATE (becq)
     DEALLOCATE (alpq)
  ENDIF
  DEALLOCATE (wgg)
  !
  !  The part of C^3 (Eq. 37) which contain the local potential can be
  !  evaluated with an integral of this change of potential and drhous
  !
  ALLOCATE (dvlocin(dffts%nnr))

  wdyn (:,:) = (0.d0, 0.d0)
  nrstot = dffts%nr1 * dffts%nr2 * dffts%nr3
  DO nu_i = 1, 3 * nat
     CALL compute_dvloc (u(1,nu_i), .FALSE., dvlocin)
     DO nu_j = 1, 3 * nat
        DO is = 1, nspin_lsda
        ! FIXME: use zgemm instead of dot_product
           wdyn (nu_j, nu_i) = wdyn (nu_j, nu_i) + &
                dot_product (drhous(1:dffts%nnr,is,nu_j), dvlocin) * &
                omega / DBLE (nrstot)
        ENDDO
     ENDDO
  ENDDO
  !
  ! collect contributions from all pools (sum over k-points)
  !
  CALL mp_sum ( dyn00, inter_pool_comm )
  CALL mp_sum ( wdyn, inter_pool_comm )
  !
  ! collect contributions from nodes of a pool (sum over G & R space)
  !
  CALL mp_sum ( wdyn, intra_bgrp_comm )

  CALL zaxpy (3 * nat * 3 * nat, (1.d0, 0.d0), wdyn, 1, dyn00, 1)
  !
  !     force this term to be hermitean
  !
  DO nu_i = 1, 3 * nat
     DO nu_j = 1, nu_i
        dyn00(nu_i,nu_j) = 0.5d0*( dyn00(nu_i,nu_j) + CONJG(dyn00(nu_j,nu_i)))
        dyn00(nu_j,nu_i) = CONJG(dyn00(nu_i,nu_j))
     ENDDO
  ENDDO
  !      call tra_write_matrix('drho dyn00',dyn00,u,nat)
  !
  !    add the augmentation term to the charge density and save it
  !
  ALLOCATE (drhoust(dfftp%nnr, nspin_mag, npertx))
  drhoust=(0.d0,0.d0)
  !
  !  The calculation of dbecsum is distributed across processors (see addusdbec)
  !  Sum over processors the contributions coming from each slice of bands
  !
  IF (noncolin) THEN
     CALL mp_sum ( dbecsum_nc, intra_bgrp_comm )
  ELSE
     CALL mp_sum ( dbecsum, intra_bgrp_comm )
  ENDIF

  IF (noncolin.AND.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 3*nat)

  mode = 0
  IF (okpaw) becsumort=(0.0_DP,0.0_DP)
  DO irr = 1, nirr
     npe = npert (irr)
     IF (doublegrid) THEN
        DO is = 1, nspin_mag
           DO iper = 1, npe
              CALL fft_interpolate (dffts, drhous(:,is,mode+iper), dfftp, &
                                                  drhoust(:,is,iper))
           ENDDO
        ENDDO
     ELSE
        CALL zcopy (dfftp%nnr*nspin_mag*npe, drhous(1,1,mode+1), 1, drhoust, 1)
     ENDIF

     CALL dscal (2*dfftp%nnr*nspin_mag*npe, 0.5d0, drhoust, 1)

     CALL addusddens (drhoust, dbecsum(1,1,1,mode+1), mode, npe, 1)
     DO iper = 1, npe
        nu_i = mode+iper
        CALL save_buffer (drhoust (1, 1, iper), lrdrhous, iudrhous, nu_i)
     ENDDO
     mode = mode+npe
  ENDDO
  !
  !  Collect the sum over k points in different pools.
  !
  IF (okpaw) CALL mp_sum ( becsumort, inter_pool_comm )

  DEALLOCATE (drhoust)
  DEALLOCATE (dvlocin)
  DEALLOCATE (dbecsum)
  IF (noncolin) DEALLOCATE (dbecsum_nc)
  DEALLOCATE (drhous)

  CALL stop_clock ('drho')
  RETURN
END SUBROUTINE drho_tpw
