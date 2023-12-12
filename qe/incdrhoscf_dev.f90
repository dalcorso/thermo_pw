!
! Copyright (C) 2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
!----------------------------------------------------------------------
SUBROUTINE incdrhoscf_dev(outk_d, npwk_d, nbndk_d, nbndk, st_d, st, npol, &
           drhoscf, dbecsum, dpsik, psicrm, dpsicrm, nbnd, nksb_ph, npe,  &
           nnrs, nnr)
!----------------------------------------------------------------------
!
!  This is a host routine that calls several global routines to 
!  update the induced charge density and in the US case also dbecsum
!  with the contribution of the block of k points currently calculated.
!  It assumes that the wavefunctions are already in a real space smooth
!  grid on the device and are contained in psicrm.
!  dpsicrm has sufficient space on device to contain real space 
!  grids for all the induced wavefunction, that are in input here 
!  in dpsik and Fourier transformed into dpsicrm.
!
USE cudafor
USE kinds,         ONLY : DP
USE uspp,          ONLY : nkb
USE lrus,          ONLY : becp1
USE wvfct,         ONLY : npwx
USE cell_base,     ONLY : omega
USE ions_base,     ONLY : nat
USE fft_base,      ONLY : dffts, dfftp
USE uspp,          ONLY : okvan
USE uspp_param,    ONLY : nhm
USE many_k_mod,    ONLY : becpk_d
USE many_k_ph_mod, ONLY : current_ikb_ph, startkb_ph, nksbx_ph, dbecsum_d
USE noncollin_module, ONLY : nspin_mag
IMPLICIT NONE

#include<h_psi_interf.f90>
#include<incdrhoscf_interf.f90>

INTEGER :: nksb_ph, nbnd, npe, npol, nnrs, nnr
! input:: the number of k points in the current block
! input:: the number of bands
! input:: the number of perturbations
! input:: the number of components of each wavefunctions
! input:: the size of the smooth fft mesh
! inpit:: the size of the thick fft mesh
LOGICAL, DEVICE :: outk_d(nksb_ph * npe)
! input: Allow to skip some set of k points (in this routine 
! is should be .false. for all sets. It is required by the fft routine).
INTEGER, DEVICE :: npwk_d(nksb_ph * npe)
! input: the number of plane waves for each set
INTEGER :: nbndk(nksb_ph * npe)
INTEGER, DEVICE :: nbndk_d(nksb_ph * npe)
! input: the number of bands to sum in each set, on host and on device
INTEGER :: st(nksb_ph * npe)
INTEGER, DEVICE :: st_d(nksb_ph * npe)
! input: the starting point of the calculation on psicrm and dpsicrm 
! on host and on device
COMPLEX(DP) :: dbecsum((nhm*(nhm+1))/2, nat, nspin_mag, npe)
! inp/out: the dbecsum to accumulate (on host)
COMPLEX(DP) :: drhoscf (nnr, nspin_mag, npe)
! inp/out: the induced charge density to accumulate (on host)
COMPLEX(DP), DEVICE :: dpsik(npwx*npol,nbnd*nksbx_ph*npe)
! input: the change of wavefunctions
REAL(DP), DEVICE :: psicrm(2,nnrs,npol,nbnd*nksb_ph*npe)
! input: the Fourier transform of the wavefunctions
REAL(DP), DEVICE :: dpsicrm(2,nnrs,npol,nbnd*nksb_ph*npe)
! space for the Fourier transform of the induced wavefunctions

INTEGER :: nr1s, nr2s, nr3s, nr1sx, nr2sx, adim, ikb, ierr
INTEGER :: ik1, id, ipert, ijh, st_, ibnd, na
COMPLEX(DP) :: dbecsum_h((nhm*(nhm+1))/2, nat, nspin_mag, nbnd*nksbx_ph*npe)

nnr=dfftp%nnr
nr1s=dffts%nr1
nr2s=dffts%nr2
nr3s=dffts%nr3
nr1sx=dffts%nr1x
nr2sx=dffts%nr2x
adim=MAX(nr1s, nr2s, nr3s)
ikb=current_ikb_ph
!
!  set the work space for fft psicrm and dpsicrm to zero
!
CALL setdpsi_zero<<<dim3(nksb_ph*npe,nbnd,nnrs/64+1),dim3(1,1,64)>>> &
        (nbndk_d, st_d, npol, dpsicrm, nbnd, nnrs, nksb_ph, npe)
ierr=cudaDeviceSynchronize()
!
!  copy there all dpsi
!
CALL set_dpsi_grid<<<dim3(nksb_ph,npe,nbnd),dim3(1,1,1)>>>(npwx, nbndk_d, &
      st_d, npol, dpsik, dpsicrm, nbnd, nnrs, nksb_ph, npe, current_ikb_ph)
ierr=cudaDeviceSynchronize()
!
!  Fourier transform of dpsi for all k and perturbations
!
CALL fft1inv_dev<<<dim3(nksb_ph*npe,nbnd,nr1s/32+1),dim3(1,1,32)>>>    &
        (outk_d, nbndk_d, st_d, npol, dpsicrm, nbnd,                   &
         nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph*npe, adim)
ierr=cudaDeviceSynchronize()
CALL fft2inv_dev<<<dim3(nksb_ph*npe,nbnd,nr1s/32+1),dim3(1,1,32)>>>    &
        (outk_d, nbndk_d, st_d, npol, dpsicrm, nbnd,                   &
         nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph*npe, adim)
ierr=cudaDeviceSynchronize()
CALL fft3inv_dev<<<dim3(nksb_ph*npe,nbnd/32+1,1),dim3(1,32,1)>>>       &
        (outk_d, nbndk_d, st_d, npol, dpsicrm, nbnd,                   &
         nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph*npe, adim)
ierr=cudaDeviceSynchronize()
!
!  final addition of this block of wavefunctions to the drho
!
!$acc enter data copyin(drhoscf(1:nnr,1:nspin_mag,1:npe))
!$acc host_data use_device(drhoscf)
CALL incdrho_dev<<<dim3(1,npe,nnrs/64+1),dim3(1,1,64)>>>(nbndk_d, st_d,    &
         npol, drhoscf, psicrm, dpsicrm, nnrs, nnr, nksb_ph,               &
         npe, nbnd, current_ikb_ph, nspin_mag, omega)
ierr=cudaDeviceSynchronize()
!$acc end host_data 
!$acc exit data copyout(drhoscf)
IF (okvan) THEN
   CALL incdrho_calbec<<<dim3(nksb_ph*npe,nkb,nbnd),dim3(1,1,1)>>>(npwx,  &
        st_d, nbndk_d, npwk_d, dpsik, current_ikb_ph, npol, nksb_ph, npe, &
        nbnd)
   ierr=cudaDeviceSynchronize()
   CALL incdrho_addusdbec<<<dim3(nksb_ph*npe,1,nbnd),dim3(1,1,1)>>>(st_d, &
                nbndk_d, current_ikb_ph, npol, nksb_ph, npe)
   ierr=cudaDeviceSynchronize()
!
!  This last sum is made on the host
!
   dbecsum_h=dbecsum_d
   DO na=1,nat
      DO ijh=1,(nhm * (nhm + 1))/2
         DO ik1=1, nksb_ph
            DO ipert=1, npe
               id=ik1+(ipert-1)*nksb_ph
               st_=st(id)
               DO ibnd=1, nbndk(id)
                  dbecsum(ijh,na,1,ipert)=dbecsum(ijh,na,1,ipert) + &
                         dbecsum_h(ijh,na,1,st_+ibnd)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF

RETURN
END SUBROUTINE incdrhoscf_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE setdpsi_zero(nbndk, st, npol, dpsicr, nbnd, &
                                           nnrs, nk, npe)
!-----------------------------------------------------------------------
  !
  !  This routine sets to zero the psicr variable. Completely parallelized
  !  on the GPU.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nnrs, nbnd
  !! input: the number of k point
  !! input: the number of perturbations
  !! input: the number of points of the smooth fft mesh
  !! input: the maximum number of bands per set
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe)
  !! input: the number of bands per set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe)
  !! input: the starting point of each set in the list of dpsicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunctions
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd * nk * npe)
  !! input: fft space to fourier tranform dpsi
  !
  INTEGER :: k0, k2, id, i, mb, st_, ipol
  !
  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe) RETURN
  st_=st(id)
  mb=nbndk(id)

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_
  i=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (i>nnrs) RETURN
!
!  zero all the dpsic
!
  DO ipol=1,npol
     dpsicr(1,i,ipol,k2)=0.0_DP
     dpsicr(2,i,ipol,k2)=0.0_DP
  ENDDO

  RETURN

  END SUBROUTINE setdpsi_zero

! -----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE set_dpsi_grid( lda, nbndk, st, npol, dpsik, &
                                        dpsicr, nbnd, nnrs, nk, npe, ikb)
!-----------------------------------------------------------------------
!
  ! This routine copy the change of the wavefunctions on the fft mesh dpsicr,
  ! in parallel for all k points and bands. dpsicr must be zero in input.
  !
  USE cudafor
  USE util_param,  ONLY : DP
  USE many_k_mod,  ONLY : nl_d
  USE many_k_ph_mod, ONLY : ikqs => ikqs_d, startkb_ph => startkb_ph_d
  USE klist,       ONLY : igk_k_d, ngk => ngk_d

  IMPLICIT NONE

  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nk, npe, nnrs, nbnd, ikb
  !! input: the number of k point
  !! input: the number of perturbations
  !! input: the number of points of the smooth fft mesh
  !! input: the maximum number of bands per set
  !! input: the current block of k points
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe)
  !! input: the number of bands of each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe)
  !! input: the starting point of each set in dpsik and dpsicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of dpsik 
  COMPLEX(DP), DEVICE, INTENT(IN) :: dpsik(lda*npol, nbnd * nk * npe)
  !! input: the dpsi in reciprocal space
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd * nk * npe)
  !! inp/out: dpsi in real space grid (in input it must be zero)
  !
  !  ... local variables
  !
  INTEGER :: k0, k1, k2, i, ik, ik1, ikq, npwq, kstart, ipert, &
             id, mb, st_, iv

  ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik1>nk) RETURN
  ik=ik1 + startkb_ph(ikb)
  ikq=ikqs(ik)
  npwq = ngk(ikq)
  kstart=nbnd*(ik1-1)

  ipert=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (ipert>npe) RETURN

  id = ik1 + (ipert-1)*nk
  st_= st(id)
  mb = nbndk(id)

  k0=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (k0>mb) RETURN
  k1 = kstart + k0
  k2 = st_ + k0
!
!  Distribute the dpsi of this thread in its FFT mesh
!
  DO i=1,npwq
     iv=nl_d(igk_k_d(i,ikq))
     dpsicr(1,iv,1,k2)= REAL(dpsik(i,k2))
     dpsicr(2,iv,1,k2)= AIMAG(dpsik(i,k2))
  ENDDO

  IF (npol==2) THEN
     DO i=1,npwq
        iv=nl_d(igk_k_d(i,ikq))
        dpsicr(1,iv,2,k2)= REAL(dpsik(lda+i,k2))
        dpsicr(2,iv,2,k2)= AIMAG(dpsik(lda+i,k2))
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE set_dpsi_grid

! -----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE incdrho_dev( nbndk, st, npol, drhoscf, &
           psicr, dpsicr, nnrs, nnr, nk, npe, nbnd, ikb, nspin_mag,  &
           omega)
!-----------------------------------------------------------------------
!
  USE cudafor
  USE kinds, ONLY : DP
  USE many_k_mod, ONLY : wk => wk_d, isk => isk_d, lsda => lsda_d
  USE many_k_ph_mod, ONLY : ikks => ikks_d, &
                            startkb_ph => startkb_ph_d
  IMPLICIT NONE

  REAL(DP), INTENT(IN), VALUE :: omega
  !! input: the volume of the unit cell
  INTEGER, INTENT(IN), VALUE :: nk, npe, nnrs, nnr, nbnd, ikb, nspin_mag
  !! input: the number of k point
  !! input: the number of perturbations
  !! input: the number of points of the smooth fft mesh
  !! input: the number of points of the thick fft mesh
  !! input: the maximum number of bands per set
  !! input: the current block of k points
  !! input: the number of components of drhoscf
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe)
  !! input: the number of bands of each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe)
  !! input: the starting point of each set in psicr and dpsicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of psicr and dpsicr
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: drhoscf(nnr, nspin_mag, npe)
  !! inp/out: the induced charge updated by this routine
  REAL(DP), DEVICE, INTENT(IN) :: psicr(2, nnrs, npol, nbnd * nk * npe)
  !! input: the evc in the real space smooth grid
  REAL(DP), DEVICE, INTENT(IN) :: dpsicr(2, nnrs, npol, nbnd * nk * npe)
  !! input: the dpsi in the real space smooth grid
  !
  !  ... local variables
  !
  INTEGER  :: ik1, ik, ikk, id, ibnd, j, ir, st_, is_, ipert
  REAL(DP) :: wgt

  COMPLEX(DP) :: asum

  ipert=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (ipert> npe) RETURN

  ir=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ir> nnrs) RETURN
  DO ik1=1, nk
     ik=ik1 + startkb_ph(ikb)
     ikk=ikks(ik)
     wgt = 2.d0 * wk(ikk) / omega
     id = ik1 + (ipert-1) * nk
     st_=st(id)
     is_=1
     IF (lsda) is_=isk(ikk)
     asum=(0.0_DP,0.0_DP)
     DO ibnd=1,nbndk(id)
        j=st_+ibnd
        asum = asum + CMPLX(            &
          psicr(1,ir,1,j)*dpsicr(1,ir,1,j)+psicr(2,ir,1,j)*dpsicr(2,ir,1,j), &
          psicr(1,ir,1,j)*dpsicr(2,ir,1,j)-psicr(2,ir,1,j)*dpsicr(1,ir,1,j), &
                 KIND=DP)
     ENDDO
     drhoscf (ir,is_,ipert)=drhoscf (ir,is_,ipert) + wgt * asum 
  ENDDO

  RETURN
END SUBROUTINE incdrho_dev
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE incdrho_calbec(ndmx, st, nbndk, npwk, dpsi, &
                                      current_ikb_ph, npol, nk, npe, nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_mod,    ONLY : vkbk_d, nkb => nkb_d
 USE many_k_ph_mod, ONLY : dbecq_d

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nbnd, current_ikb_ph
 !! input: the number of k point
 !! input: the number of perturbations
 !! input: the number of points of the smooth fft mesh
 !! input: the number of points of the thick fft mesh
 !! input: the maximum number of bands per set
 !! input: the current block of k points
 INTEGER, DEVICE :: npwk(nk*npe)
 !! input: the number of plane waves for each set
 INTEGER, DEVICE :: nbndk(nk*npe)
 !! input: the number of bands of each set
 INTEGER, DEVICE :: st(nk*npe)
 !! input: the starting point of each set in psicr and dpsicr
 INTEGER, VALUE :: npol
 !! input: the number of components of dpsi
 COMPLEX(DP), DEVICE :: dpsi(ndmx*npol,nbnd*nk*npe)
 !! input: the change of the wavefunctions dpsi

 INTEGER :: ik1, ik, ikk, ipert, id, st_, ibnd, ig, kstart, nbd, k, j, n, i0
 COMPLEX(DP) :: asum

 id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (id>nk*npe) RETURN
 ipert=(id-1)/nk+1
 ik1=MOD(id-1,nk)+1

 kstart=nkb*(ik1-1)
 n=npwk(id)

 st_=st(id)

 i0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (i0>nkb) RETURN

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbndk(id)) RETURN

 asum=(0.0_DP,0.0_DP)
 DO j=1, n
    asum=asum+CONJG(vkbk_d(j,kstart+i0))*dpsi(j,st_+ibnd)
 ENDDO
 dbecq_d(i0,ibnd,id)=asum

 RETURN

END SUBROUTINE incdrho_calbec
!
!--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE incdrho_addusdbec(st, nbndk, current_ikb_ph, &
                               npol, nk, npe)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP
 USE many_k_mod,    ONLY : ntyp => ntyp_d, nat => nat_d, tvanp => tvanp_d, &
                           ityp => ityp_d, nh => nh_d, wk => wk_d,         &
                           nhm => nhm_d
 USE many_k_ph_mod, ONLY : dbecq_d, dbecsum_d, ikks => ikks_d,             &
                           startkb_ph => startkb_ph_d, becp1k_d
 USE uspp,          ONLY : ijtoh => ijtoh_d

 IMPLICIT NONE

 INTEGER, VALUE :: nk, npe, npol, current_ikb_ph

 INTEGER, DEVICE :: nbndk(nk*npe)
 !! input: the number of bands of each set
 INTEGER, DEVICE :: st(nk*npe)
 !! input: the starting point of each set in psicr and dpsicr

 INTEGER :: ik1, ik, ikk, ipert, id, st_, ibnd, ijkb0, nt, na, ih, jh, ijh, &
            ikb, jkb, mb 
 REAL(DP) :: wgt

 id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (id>nk*npe) RETURN
 ipert=(id-1)/nk+1
 ik1=MOD(id-1,nk)+1

 ik=startkb_ph(current_ikb_ph)+ik1 
 ikk = ikks(ik)
 mb = nbndk(id)

 wgt = wk(ikk) 

 st_=st(id)

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>mb) RETURN

 dbecsum_d(1:(nhm*(nhm+1))/2, 1:nat, 1, st_+ibnd)=CMPLX(0.0_DP,0.0_DP,KIND=DP)

 ijkb0 = 0
 DO nt = 1, ntyp
    IF (tvanp(nt) ) THEN
       DO na = 1, nat
          IF (ityp (na)==nt) THEN
             !
             DO ih = 1, nh(nt)
                ikb = ijkb0 + ih
                ijh=ijtoh(ih,ih,nt)
                dbecsum_d (ijh, na, 1, st_+ibnd) = &
                   dbecsum_d(ijh, na, 1, st_+ibnd) + wgt * &
                       (CONJG(becp1k_d(ikb,1,ibnd,ik))*dbecq_d(ikb,ibnd,id))
                DO jh = ih + 1, nh (nt)
                   ijh=ijtoh(ih,jh,nt)
                   jkb = ijkb0 + jh
                   dbecsum_d (ijh, na, 1, st_+ibnd) =                        &
                           dbecsum_d (ijh, na, 1, st_+ibnd) +                &
                    wgt*(CONJG(becp1k_d(ikb,1,ibnd,ik))*dbecq_d(jkb,ibnd,id)+&
                         CONJG(becp1k_d(jkb,1,ibnd,ik))*dbecq_d(ikb,ibnd,id) )
                ENDDO
                ijh = ijh + 1
             ENDDO
             ijkb0 = ijkb0 + nh (nt)
          ENDIF
       ENDDO
    ELSE
       DO na = 1, nat
          IF (ityp (na)==nt) ijkb0 = ijkb0 + nh (nt)
       ENDDO
    ENDIF
 ENDDO

 RETURN
END SUBROUTINE incdrho_addusdbec
#endif
