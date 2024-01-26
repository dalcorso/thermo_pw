!
! Copyright (C) 2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
!------------------------------------------------------------------------
SUBROUTINE dvqpsi_dev(iter, ikb, nk, npe, nsolv, nnrs, nspin_mag, npol,     &
                      imode0, outk_d, npwk_d, nbndk_d, st_d, ikt_d, psicrm, &
                      dpsicrm, dvloc_d, dvscfins)
!------------------------------------------------------------------------
!
USE cudafor
USE kinds,     ONLY : DP
USE wvfct,     ONLY : nbnd, npwx
USE cell_base, ONLY : tpiba
USE many_k_ph_mod, ONLY : nksb_ph, dvpsik_d
USE uspp,      ONLY : okvan, nkb
USE lsda_mod,  ONLY : nspin
USE fft_base,  ONLY : dffts
!
!   When iter=1 this subroutines applies dv_bare to the wavefunctions
!   and save it in dvpsik_d, which is overwritten.
!   When iter>1 it applies dvscfins to the wavefunctions and adds it
!   to dvpsik_d.
!   It does this for all current k points and all current perturbations
!   and for all bands. It uses global routines so the calculation is
!   done on GPU by different threads, each thread computing one band
!   one k point and one perturbation. In some cases also the G vectors,
!   or FFT grid points are distributed among threads.
!
IMPLICIT NONE
#include<dvqpsi_interf.f90>
#include<h_psi_interf.f90>

INTEGER :: iter, ikb, nk, npe, nsolv, nnrs, imode0, nspin_mag, npol
! input: the current iteration
! input: the number of the block of k points 
! input: the number of k points
! input: the number of perturbations
! input: nsolv=2 in the noncollinear magnetic case, nsolv=1 otherwise
! input: the number of points of the smooth fft mesh
! input: the starting mode of these perturbations
! input: the number of components of dvscfins
! input: the number of components of the wavefunctions
LOGICAL, DEVICE :: outk_d(nk * npe * nsolv)
!  input: auxiliary array needed in the fft routine. Must be .FALSE. here.
INTEGER, DEVICE :: npwk_d(nk * npe * nsolv)
!  input: number of plane waves for each set of functions
INTEGER, DEVICE :: nbndk_d(nk * npe * nsolv)
!  input: number of bands to compute in each set
INTEGER, DEVICE :: ikt_d(nk)
!  input: the global k point index for each k point
INTEGER, DEVICE :: st_d(nk * npe * nsolv)
!  input: the starting point of each set in psicrm and dpsicrm
REAL(DP), DEVICE :: psicrm(2, nnrs, npol, nbnd*nk*npe*nsolv)
!  inp/out: A smooth fft grid for saving the wavefunctions in real space.
REAL(DP), DEVICE :: dpsicrm(2, nnrs, npol, nbnd*nk*npe*nsolv)
!  inp/out: A smooth fft grid for saving the product of the wavefunctions and
!           the potential.
COMPLEX(DP), DEVICE :: dvloc_d (nnrs, npe)
!  input: the bare potential in the smooth fft mesh for each perturbation
COMPLEX(DP), DEVICE :: dvscfins ( nnrs, nspin_mag, npe)
!  input: the Hxc induced potential in the smooth fft mesh 
INTEGER :: nr1s, nr2s, nr3s, nr1sx, nr2sx, adim
INTEGER :: ierr, i

nr1s=dffts%nr1
nr2s=dffts%nr2
nr3s=dffts%nr3
nr1sx=dffts%nr1x
nr2sx=dffts%nr2x
adim=MAX(nr1s, nr2s, nr3s)

!
!  Clean the work space for fft psicrm 
!
   CALL setpsi_zero<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd,nnrs/32+1),&
      dim3(1,1,32)>>>(nbndk_d, st_d, npol, psicrm, nbnd, nnrs,    &
      nksb_ph(ikb), npe, nsolv)
   ierr=cudaDeviceSynchronize()
!
!  copy the evc for all k on the FFT mesh
!
   CALL setpsi_evc<<<dim3(nksb_ph(ikb),npe*nsolv,nbnd),dim3(1,1,1)>>>   &
           (nbndk_d, st_d, npol, psicrm, nbnd, nnrs, nksb_ph(ikb), npe, &
            nsolv, ikb)
   ierr=cudaDeviceSynchronize()
!
!  Fourier transform of evc
!
   CALL fft1inv_dev<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd,nr1s/32+1),          &
        dim3(1,1,32)>>>(outk_d, nbndk_d, st_d, npol, psicrm, nbnd,          &
        nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph(ikb)*npe*nsolv, adim)
   ierr=cudaDeviceSynchronize()
   CALL fft2inv_dev<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd,nr1s/32+1),          &
        dim3(1,1,32)>>>(outk_d, nbndk_d, st_d, npol, psicrm, nbnd,          &
        nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph(ikb)*npe*nsolv, adim)
   ierr=cudaDeviceSynchronize()
   CALL fft3inv_dev<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd/32+1,1),             &
        dim3(1,32,1)>>>(outk_d, nbndk_d, st_d, npol, psicrm, nbnd,          &
        nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph(ikb)*npe*nsolv, adim)
   ierr=cudaDeviceSynchronize()
 !
 !   ... apply potential to psi
 !
   IF (iter==1) THEN
      !
      ! at the first iteration apply dv_bare
      !
      CALL dvlocpsi_dev<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd,nnrs/32+1),     &
            dim3(1,1,32)>>> (nbndk_d, st_d, npol, psicrm, dpsicrm,         &
                             dvloc_d, nbnd, nnrs, nksb_ph(ikb), npe, nsolv)
      ierr=cudaDeviceSynchronize()
   ELSE
      !
      ! at the other iterations apply dvscfins
      !
      CALL dvscfinpsi_gpu<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd,nnrs/32+1),   &
           dim3(1,1,32)>>> (nbndk_d, st_d, npol, psicrm, dpsicrm,          &
           dvscfins, nspin_mag, nbnd, nnrs, nksb_ph(ikb), npe, nsolv, ikb)
      ierr=cudaDeviceSynchronize()
   ENDIF
   !
   !   Now return to reciprocal space
   !
   !
   !   ... direct fft in the direction x
   !
   CALL fft3fwd_dev<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd/32+1,1),         &
            dim3(1,32,1)>>>(outk_d, nbndk_d, st_d, npol, dpsicrm, nbnd, &
            nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph(ikb)*npe*nsolv, adim)
   ierr=cudaDeviceSynchronize()
   !
   !   ... direct fft in the direction y
   !
   CALL fft2fwd_dev<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd,nr1s/32+1),      &
          dim3(1,1,32)>>>(outk_d, nbndk_d, st_d, npol, dpsicrm, nbnd,   &
            nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph(ikb)*npe*nsolv, adim)
   ierr=cudaDeviceSynchronize()
   !
   !   ... direct fft in the direction z
   !
   CALL fft1fwd_dev<<<dim3(nksb_ph(ikb)*npe*nsolv,nbnd,nr1s/32+1),        &
             dim3(1,1,32)>>> (outk_d, nbndk_d, st_d, npol, dpsicrm, nbnd, &
             nr1s, nr2s, nr3s, nr1sx, nr2sx, nnrs, nksb_ph(ikb)*npe*nsolv,& 
             adim)
   ierr=cudaDeviceSynchronize()
   !
   IF (iter==1) THEN
      !
      !   at the first iteration dpsicmr contains dv_bare psi and
      !   this must be copied into dvpsi. 
      !
      CALL set_dvpsik_dev<<<dim3(nksb_ph(ikb),npe*nsolv,nbnd),dim3(1,1,1)>>>&
          ( npwx, npwk_d, nbndk_d, st_d, ikt_d, npol, dvpsik_d, dpsicrm,   &
            nbnd, nnrs, nksb_ph(ikb), npe, nsolv)
      ierr=cudaDeviceSynchronize()
      !
      !  add the terms that used to be in dvqpsi_us_only
      !
      CALL dvqpsi_us_dev0<<<dim3(nksb_ph(ikb)*npe*nsolv,1,nbnd),dim3(1,1,1)>>>&
               (ikt_d, st_d, npol, ikb, nbnd, nksb_ph(ikb), npe, &
                nsolv, imode0, nspin)
      ierr=cudaDeviceSynchronize()
      CALL dvqpsi_us_dev1<<<dim3(nksb_ph(ikb)*npe*nsolv,npwx,nbnd),         &
                dim3(1,1,1)>>>(npwx, st_d, ikt_d, dvpsik_d,                 &
                               npol, nkb, nbnd, nksb_ph(ikb), npe, nsolv)
      ierr=cudaDeviceSynchronize()
      CALL dvqpsi_us_dev2<<<dim3(nksb_ph(ikb)*npe*nsolv,npwx,nbnd),         &
                dim3(1,1,1)>>>(npwx, st_d, ikt_d, dvpsik_d, npol,           &
                nkb, nbnd, nksb_ph(ikb), npe, nsolv)
      ierr=cudaDeviceSynchronize()
   ELSE
      !
      ! at the other iterations psicrm contains dvscf psi that
      ! must be added to dvpsi 
      !
      CALL adddvscf_dev<<<dim3(nksb_ph(ikb), npe*nsolv, nbnd),dim3(1,1,1)>>> &
          ( npwx, npwk_d, nbndk_d, st_d, ikt_d, npol, dvpsik_d, dpsicrm,    &
            nbnd, nnrs, nksb_ph(ikb), npe, nsolv)
      ierr=cudaDeviceSynchronize()
   ENDIF   
RETURN
END SUBROUTINE dvqpsi_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE setpsi_zero(nbndk, st, npol, psicr, nbnd, &
                                          nnrs, nk, npe, nsolv)
!-----------------------------------------------------------------------
  !
  !  This routine sets to zero the psicr variable in parallel on the GPU.
  !  Each thread computes one k point, one perturbation, one band and
  !  one FFT grid point.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of linear systems to solve
  !! input: the number of points in the smooth FFT mesh
  !! input: the number of bands
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: the number of bands to compute for each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! inp/out: the real space grid to be set to zero

  INTEGER :: k0, k2, id, i, mb, ipol, st_
  !
  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe*nsolv) RETURN
  st_=st(id)
  mb=nbndk(id)

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_

  i=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (i>nnrs) RETURN
!
!  zero all the psicr and dpsic
!
  DO ipol=1, npol
     psicr(1,i,ipol,k2)=0.0_DP
     psicr(2,i,ipol,k2)=0.0_DP
  ENDDO

  RETURN
  END SUBROUTINE setpsi_zero

! -----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE setpsi_evc(nbndk, st, npol, psicr, nbnd, &
                                         nnrs, nk, npe, nsolv, ikb)
!-----------------------------------------------------------------------
!
!  This routine distributes the wavefunctions evc on the fft mesh.
!  It is parallel on the GPU, each thread computing one k point,
!  one perturbation, and one band.
!
  USE cudafor
  USE util_param,  ONLY : DP
  USE many_k_mod,  ONLY : evck => evck_d, nl_d, npwx => npwx_d
  USE many_k_ph_mod, ONLY : ikks => ikks_d, startkb_ph => startkb_ph_d
  USE klist,       ONLY : igk_k_d, ngk => ngk_d

  IMPLICIT NONE
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd, ikb
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of linear systems to solve
  !! input: the number of points in the smooth FFT mesh
  !! input: the number of the k point block
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: the number of bands of each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! inp/out: the real space grid set to zero
  !
  !  ... local variables
  !
  INTEGER :: k0, k1, k2, ik, ik1, ikk, npw, kstart, ipert, id, mb, st_, iv, &
             isp, isolv, i

  ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik1>nk) RETURN
  ik=ik1 + startkb_ph(ikb)
  ikk=ikks(ik)
  npw=ngk(ikk)
  kstart=nbnd*(ik1-1)

  isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (isp>npe*nsolv) RETURN
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1

  id = ik1 + (ipert-1)*nk + (isolv-1)*nk*npe
  st_= st(id)
  mb = nbndk(id)

  k0=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (k0>mb) RETURN
  k1 = kstart + k0 + (isolv-1)*nk*nbnd
  k2 = st_ + k0
!
!  Distribute the evc of this thread in its FFT mesh
!
  DO i=1,npw
     iv=nl_d(igk_k_d(i,ikk))
     psicr(1,iv,1,k2)= DBLE(evck(i,k1))
     psicr(2,iv,1,k2)= AIMAG(evck(i,k1))
  ENDDO

  IF (npol==2) THEN
     DO i=1,npw
        iv=nl_d(igk_k_d(i,ikk))
        psicr(1,iv,2,k2)= DBLE(evck(npwx+i,k1))
        psicr(2,iv,2,k2)= AIMAG(evck(npwx+i,k1))
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE setpsi_evc
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvlocpsi_dev(nbndk, st, npol, psicr, dpsicr, &
                              dvloc, nbnd, nnrs, nk, npe, nsolv)
!----------------------------------------------------------------------- 
  !  
  ! This routine apply the local potential to the wavefunctions 
  ! in parallel for all k points, perturbation, band, and smooth FFT
  ! grid points.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of points in the smooth FFT mesh
  !! input: the number of bands
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: the number of bands of each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! input: the wavefunctions to which the potential is applied
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! output: the result. dpsicr is completely overwritten.
  COMPLEX(DP), DEVICE, INTENT(IN) :: dvloc(nnrs, npe)
  !! input: the potential to apply
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, i, ik1, id, ipert, isolv, isp, mb, st_

  COMPLEX(DP) :: aux
  !
  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe*nsolv) RETURN
  isp=(id-1)/nk+1
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1

  st_=st(id)
  mb=nbndk(id) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  i=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (i>nnrs) RETURN

  aux=CMPLX(psicr(1,i,1,k2), psicr(2,i,1,k2), KIND=DP)
  aux = aux * dvloc(i,ipert)

  dpsicr(1,i,1,k2) = DBLE(aux)
  dpsicr(2,i,1,k2) = AIMAG(aux)

  IF (npol==2) THEN
     aux=CMPLX(psicr(1,i,2,k2), psicr(2,i,2,k2), KIND=DP)
     aux = aux * dvloc(i,ipert)

     dpsicr(1,i,2,k2) = DBLE(aux)
     dpsicr(2,i,2,k2) = AIMAG(aux)
  ENDIF


  RETURN
  !
END SUBROUTINE dvlocpsi_dev
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvscfinpsi_gpu(nbndk, st, npol, psicr, dpsicr, &
                   dvscfins, nspin_mag, nbnd, nnrs, nk, npe, nsolv, ikb)
!----------------------------------------------------------------------- 
  !  
  ! This routine apply the local potential to the wavefunctions 
  ! in parallel for all k points and bands.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : isk => isk_d, lsda => lsda_d, noncolin => &
                             noncolin_d, domag => domag_d
  USE many_k_ph_mod,  ONLY : startkb_ph => startkb_ph_d, ikks => ikks_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd, nspin_mag, ikb
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of linear systems to solve 
  !! input: the number of points in the smooth FFT mesh
  !! input: the number of bands
  !! input: the number of components of the Hxc induced potential
  !! input: the number of the k point block
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: the number of bands of each set
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! input: the wavefunctions to which the potential is applied
  REAL(DP), DEVICE, INTENT(INOUT) :: dpsicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! output: the result. dpsicr is completely overwritten.
  COMPLEX(DP), DEVICE, INTENT(IN) :: dvscfins(nnrs, nspin_mag, npe)
  !! input: the induced Hxc potential for each perturbation
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, i, ik, ik1, ikk, id, isp, isolv, ipert, mb, st_, is_

  COMPLEX(DP) :: aux, aux1, sup, sdwn
  !
  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe*nsolv) RETURN
  isp=(id-1)/nk+1
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1

  ik1=MOD(id-1,nk)+1
  ik=ik1+startkb_ph(ikb)
  ikk=ikks(ik)
  st_=st(id)
  mb=nbndk(id) 
  is_= 1
  IF (lsda) is_=isk(ikk)

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  i=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (i>nnrs) RETURN

  aux=CMPLX(psicr(1,i,1,k2), psicr(2,i,1,k2), KIND=DP)
  IF (noncolin) THEN
     aux1=CMPLX(psicr(1,i,2,k2), psicr(2,i,2,k2), KIND=DP)
     IF (domag) THEN
        IF (isolv==1) THEN
           sup=aux*(dvscfins(i,1,ipert)+dvscfins(i,4,ipert))+ &
               aux1*(dvscfins(i,2,ipert)-(0.d0,1.d0)*dvscfins(i,3,ipert))
           sdwn=aux1*(dvscfins(i,1,ipert)-dvscfins(i,4,ipert)) + &
                aux*(dvscfins(i,2,ipert)+(0.d0,1.d0)*dvscfins(i,3,ipert))
        ELSE
           sup=aux*(dvscfins(i,1,ipert)-dvscfins(i,4,ipert))+ &
              -aux1*(dvscfins(i,2,ipert)-(0.d0,1.d0)*dvscfins(i,3,ipert))
           sdwn=aux1*(dvscfins(i,1,ipert)+dvscfins(i,4,ipert)) + &
               -aux*(dvscfins(i,2,ipert)+(0.d0,1.d0)*dvscfins(i,3,ipert))
        ENDIF
        aux=sup
        aux1=sdwn
     ELSE
        aux1=CMPLX(psicr(1,i,2,k2), psicr(2,i,2,k2), KIND=DP)
        aux=aux*dvscfins(i,1,ipert)
        aux1=aux1*dvscfins(i,1,ipert)
     ENDIF
     dpsicr(1,i,2,k2) = DBLE(aux1)
     dpsicr(2,i,2,k2) = AIMAG(aux1)
  ELSE
     aux = aux * dvscfins(i,is_,ipert)
  ENDIF
  dpsicr(1,i,1,k2) = DBLE(aux)
  dpsicr(2,i,1,k2) = AIMAG(aux)

  RETURN
  !
END SUBROUTINE dvscfinpsi_gpu
!
!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE set_dvpsik_dev(lda, npw, nbndk, st, ikt,  &
                npol, dvpsik, psicr, nbnd, nnrs, nk, npe, nsolv) 
!----------------------------------------------------------------------- 
!  
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : nl_d
  USE klist,          ONLY : igk_k_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of linear systems to solve
  !! input: the number of points in the smooth FFT mesh
  !! input: the number of bands
  INTEGER, INTENT(IN), DEVICE :: ikt(nk)
  !! input: the global k point index for each k point
  INTEGER, INTENT(IN), DEVICE :: npw(nk*npe*nsolv)
  !! input: number of plane waves for each set of functions
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: number of bands for each set of functions
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: dvpsik(lda*npol, nbnd*nk*npe*nsolv)
  !! output: the result. dvpsik is completely overwritten.
  REAL(DP), DEVICE, INTENT(IN) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! input: the fft smooth mesh to save in dvpsik
  !
  !  ... local variables
  !
  INTEGER :: k0, k1, ik, ik1, ikq, isp, i, id, np, mb, st_, iv, isolv, ipert
  !
  ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik1>nk) RETURN
  ikq=ikt(ik1)

  isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (isp>npe*nsolv) RETURN
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1

  id=ik1+(ipert-1)*nk + (isolv-1)*npe*nk
  np = npw(id)
  st_=st(id)
  mb=nbndk(id) 

  k0=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (k0>mb) RETURN
  k1 = k0 + st_ 
!
!  Collect hpsi_d of this thread in its FFT mesh
!
  DO i=1,np
     iv=nl_d(igk_k_d(i,ikq))
     dvpsik(i,k1)=CMPLX(psicr(1,iv,1,k1), psicr(2,iv,1,k1), KIND=DP)
  ENDDO

  DO i=np+1,lda
     dvpsik(i,k1)=(0.0_DP,0.0_DP)
  ENDDO

  IF (npol==2) THEN
     DO i=1,np
        iv=nl_d(igk_k_d(i,ikq))
        dvpsik(lda+i,k1)=CMPLX(psicr(1,iv,2,k1), psicr(2,iv,2,k1), KIND=DP)
     ENDDO
     DO i=np+1,lda
        dvpsik(lda+i,k1)=(0.0_DP,0.0_DP)
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE set_dvpsik_dev

!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE adddvscf_dev(lda, npw, nbndk, st, ikt,  &
                npol, dvpsik, psicr, nbnd, nnrs, nk, npe, nsolv) 
!----------------------------------------------------------------------- 
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : nl_d
  USE klist,          ONLY : igk_k_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nk, npe, nsolv, nnrs, nbnd
  !! input: the number of k points
  !! input: the number of perturbations
  !! input: the number of linear systems to solve
  !! input: the number of points in the smooth FFT mesh
  !! input: the number of bands
  INTEGER, INTENT(IN), DEVICE :: ikt(nk)
  !! input: the global k point index for each k point
  INTEGER, INTENT(IN), DEVICE :: npw(nk*npe*nsolv)
  !! input: number of plane waves for each set of functions
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !! input: number of bands for each set of functions
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  !! input: the starting point of the wavefunctions for each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: dvpsik(lda*npol, nbnd*nk*npe*nsolv)
  !! output: the result. The result is added to dvpsik.
  REAL(DP), DEVICE, INTENT(IN) :: psicr(2, nnrs, npol, nbnd*nk*npe*nsolv)
  !! input: the fft smooth mesh to save in dvpsik
  !
  !  ... local variables
  !
  INTEGER :: k0, k1, ik, ik1, ikq, i, id, np, mb, st_, iv, ipert, isolv, isp
  !
  ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik1>nk) RETURN
  ikq=ikt(ik1)

  isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (isp > npe*nsolv) RETURN
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1

  id=ik1 + (ipert-1)*nk + (isolv-1) * nk * npe
  np = npw(id)
  st_=st(id)
  mb=nbndk(id) 

  k0=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (k0>mb) RETURN
  k1 = k0 + st_ 
!
!  Collect dvpsik of this thread
!
  DO i=1,np
     iv=nl_d(igk_k_d(i,ikq))
     dvpsik(i,k1)=dvpsik(i,k1)+ &
                     CMPLX(psicr(1,iv,1,k1), psicr(2,iv,1,k1), KIND=DP)
  ENDDO

  DO i=np+1,lda
     dvpsik(i,k1)=(0.0_DP,0.0_DP)
  ENDDO

  IF (npol==2) THEN
     DO i=1,np
        iv=nl_d(igk_k_d(i,ikq))
        dvpsik(lda+i,k1)=dvpsik(lda+i,k1)+ &
                     CMPLX(psicr(1,iv,2,k1), psicr(2,iv,2,k1), KIND=DP)
     ENDDO

     DO i=np+1,lda
        dvpsik(lda+i,k1)=(0.0_DP,0.0_DP)
     ENDDO

  ENDIF
  RETURN
  END SUBROUTINE adddvscf_dev
!
!------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvqpsi_us_dev0(ikt, st, npol, &
                    current_ikb_ph, nbnd, nk, npe, nsolv, imode0, nspin)
!------------------------------------------------------------------------
!
!  This routines computes the coefficients ps1k_d, and ps2k_d,
!  or ps1k_nc_d and ps2k_nc_d in the noncollinear case.
!  These coefficients where originally computed in add_dvqpsi_us_only.f90
!
USE cudafor
USE kinds, ONLY : DP
USE many_k_ph_mod, ONLY : ps1k_d, ps2k_d, int1_d, int2_d, ikks => ikks_d, &
                 startkb_ph => startkb_ph_d, u => u_d, deff_d, alphak_d,  &
                 becp1k_d, deff_nc_d, int1_nc_save_d, int2_so_d, ps1k_nc_d, &
                 ps2k_nc_d, alphatk_d, becptk_d, int1_nc => int1_nc_d
USE many_k_mod, ONLY : vkbk_d, nat => nat_d, ntyp => ntyp_d, &
                ityp => ityp_d, nh => nh_d, isk => isk_d, nhm => nhm_d,   &
                qq_at => qq_at_d, deeq =>deeq_d, noncolin => noncolin_d,  &
                lspinorb => lspinorb_d, lsda => lsda_d, nkb => nkb_d,     &
                okvan => okvan_d, tpiba => tpiba_d 
USE wvfct_gpum, ONLY : et => et_d
 
IMPLICIT NONE

#include<compute_deff_interf.f90>

INTEGER, VALUE :: nk, npe, nsolv, nbnd, current_ikb_ph, imode0, nspin
!! input: the number of k points
!! input: the number of perturbations
!! input: the number of linear systems to solve
!! input: the number of bands
!! input: the number of the k points block
!! input: the current starting mode
INTEGER, DEVICE :: ikt(nk)
!! input: the global k point index for each k point
INTEGER, DEVICE :: st(nk*npe*nsolv)
!! input: the starting point of the wavefunctions for each set in psicr
INTEGER, VALUE :: npol
!! input: the number of components of each wavefunction

REAL(DP), PARAMETER :: eps = 1.d-12
INTEGER  :: iv, ig, ik, ikk, ik1, id, ipert, ikq, ibnd
INTEGER  :: ijkb0, nt, na, nb, mu, ih, jh, ikb, jkb, ipol, nu, st_, &
            isp, isolv, current_spin
INTEGER :: is, js, ijs
COMPLEX(DP) :: asum
COMPLEX(DP) :: uact(3*nat)

id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
IF (id>nk*npe*nsolv) RETURN
isp=(id-1)/nk+1
isolv=(isp-1)/npe+1
ipert=MOD(isp-1,npe)+1

ik1=MOD(id-1,nk)+1
ik=ik1 + startkb_ph(current_ikb_ph)
ikk=ikks(ik)
ikq=ikt(ik1)
current_spin=1
IF (lsda) current_spin=isk(ikk)
st_=st(id)

ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ibnd>nbnd) RETURN

IF (noncolin) THEN
   CALL compute_deff_nc_dev(nhm, nat, nspin, isolv, npol, okvan, nsolv, &
                                deff_nc_d(1,1,1,1,st_+ibnd), et(ibnd,ikk))
ELSE
   CALL compute_deff_dev(nhm, nat, current_spin, okvan, &
                                  deff_d(1,1,1,st_+ibnd), et(ibnd,ikk))
ENDIF

IF (noncolin) THEN
   ps1k_nc_d(:,:,ibnd,id)=(0.0_DP, 0.0_DP)
   ps2k_nc_d(:,:,ibnd,:,id)=(0.0_DP, 0.0_DP)
ELSE
   ps1k_d(:,ibnd,id)=(0.0_DP, 0.0_DP)
   ps2k_d(:,ibnd,:,id)=(0.0_DP, 0.0_DP)
ENDIF

DO na=1,nat
   mu = 3 * (na - 1)
   DO ipol=1,3 
      uact(mu+ipol)=u(mu+ipol, imode0+ipert)
   ENDDO
ENDDO

ijkb0 = 0
DO nt = 1, ntyp
   DO na = 1, nat
      IF (ityp (na)==nt) THEN
         mu = 3 * (na - 1)
         DO ih = 1, nh (nt)
            ikb = ijkb0 + ih
            DO jh = 1, nh (nt)
               jkb = ijkb0 + jh
               DO ipol = 1, 3
                  IF ( ABS (uact (mu + 1) ) + &
                       ABS (uact (mu + 2) ) + &
                       ABS (uact (mu + 3) ) > eps) THEN
                     IF (noncolin) THEN
                        ijs=0
                        DO is=1,npol
                           DO js=1,npol
                              ijs=ijs+1
                              IF (isolv==1) THEN
                                 ps1k_nc_d(ikb,is,ibnd,id)=                 &
                                      ps1k_nc_d(ikb,is,ibnd,id) +           &
                                      deff_nc_d(ih,jh,na,ijs,st_+ibnd) *    &
!                                      alphak_d(jkb,js,ibnd,ipol,ik1)*       &
                                      alphak_d(jkb,js,ibnd,ipol,ik)*       &
                                       uact(mu + ipol)
                                 ps2k_nc_d(ikb,is,ibnd,ipol,id)=            &
                                      ps2k_nc_d(ikb,is,ibnd,ipol,id)+       &
                                         deff_nc_d(ih,jh,na,ijs,st_+ibnd) * &
!                                          becp1k_d(jkb,js,ibnd,ik1) *       &
                                          becp1k_d(jkb,js,ibnd,ik) *       &
                                       (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                              ELSE
                                 ps1k_nc_d(ikb,is,ibnd,id)=                 &
                                      ps1k_nc_d(ikb,is,ibnd,id) +           &
                                      deff_nc_d(ih,jh,na,ijs,st_+ibnd) *    &
!                                      alphatk_d(jkb,js,ibnd,ipol,ik1)*       &
                                      alphatk_d(jkb,js,ibnd,ipol,ik)*       &
                                       uact(mu + ipol)
                                 ps2k_nc_d(ikb,is,ibnd,ipol,id)=            &
                                      ps2k_nc_d(ikb,is,ibnd,ipol,id)+       &
                                         deff_nc_d(ih,jh,na,ijs,st_+ibnd) * &
!                                          becptk_d(jkb,js,ibnd,ik1) *       &
                                          becptk_d(jkb,js,ibnd,ik) *       &
                                       (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                              ENDIF
                           END DO
                        END DO
                     ELSE
                        ps1k_d (ikb, ibnd, id) = ps1k_d (ikb, ibnd, id) +    &
                                deff_d(ih, jh, na, st_+ibnd) *               &
!                           alphak_d(jkb, 1, ibnd, ipol, ik1) * uact (mu + ipol)
                           alphak_d(jkb, 1, ibnd, ipol, ik) * uact (mu + ipol)
                        ps2k_d (ikb, ibnd, ipol, id) =                       &
                             ps2k_d (ikb, ibnd, ipol, id) +                  &
!                         deff_d(ih,jh,na,st_+ibnd)*becp1k_d(jkb,1,ibnd,ik1)* &
                         deff_d(ih,jh,na,st_+ibnd)*becp1k_d(jkb,1,ibnd,ik)* &
                        (0.0_DP,-1.0_DP) * uact (mu + ipol) * tpiba
                     ENDIF
                     IF (okvan) THEN
                        IF (noncolin) THEN
                           ijs=0
                           DO is=1,npol
                              DO js=1,npol
                                 ijs=ijs+1
                                 IF (nsolv==2) THEN
                                    IF (isolv==1) THEN
                                       ps1k_nc_d(ikb,is,ibnd,id)=         &
                                         ps1k_nc_d(ikb,is,ibnd,id)+       &
                                      int1_nc_save_d(ih,jh,ipol,na,ijs,isolv)*&
!                                          becp1k_d(jkb,js,ibnd,ik1)*      &
                                          becp1k_d(jkb,js,ibnd,ik)*      &
                                          uact(mu+ipol)
                                    ELSE
                                       ps1k_nc_d(ikb,is,ibnd,id)=            &
                                       ps1k_nc_d(ikb,is,ibnd,id)+            &
                                  int1_nc_save_d(ih,jh,ipol,na,ijs,isolv) *  &
!                                      becptk_d(jkb,js,ibnd,ik1)*             &
                                      becptk_d(jkb,js,ibnd,ik)*             &
                                      uact(mu+ipol)
                                    ENDIF
                                 ELSE
                                    ps1k_nc_d(ikb,is,ibnd,id)=         &
                                      ps1k_nc_d(ikb,is,ibnd,id)+       &
                                   int1_nc(ih,jh,ipol,na,ijs) *  &
!                                       becp1k_d(jkb,js,ibnd,ik1)*      &
                                       becp1k_d(jkb,js,ibnd,ik)*      &
                                       uact(mu+ipol)
                                 ENDIF
                              END DO
                           END DO
                        ELSE
                           ps1k_d (ikb, ibnd, id) = ps1k_d(ikb, ibnd, id) +  &
                             (int1_d (ih, jh, ipol, na, current_spin) *      &
!                             becp1k_d(jkb, 1, ibnd, ik1) ) * uact (mu +ipol)
                             becp1k_d(jkb, 1, ibnd, ik) ) * uact (mu +ipol)
                        ENDIF
                     ENDIF
                  END IF  ! uact>0
                  IF (okvan) THEN
                     DO nb = 1, nat
                        nu = 3 * (nb - 1)
                        IF (noncolin) THEN
                           IF (lspinorb) THEN
                              ijs=0
                              DO is=1,npol
                                 DO js=1,npol
                                    ijs=ijs+1
                                    IF (isolv==1) THEN
                                       ps1k_nc_d(ikb,is,ibnd,id)=            &
                                       ps1k_nc_d(ikb,is,ibnd,id)+            &
                                       int2_so_d(ih,jh,ipol,nb,na,ijs)*      &
!                                       becp1k_d(jkb,js,ibnd,ik1)*uact(nu+ipol)
                                       becp1k_d(jkb,js,ibnd,ik)*uact(nu+ipol)
                                    ELSE
                                       ps1k_nc_d(ikb,is,ibnd,id)=            &
                                       ps1k_nc_d(ikb,is,ibnd,id)+            &
                                       int2_so_d(ih,jh,ipol,nb,na,ijs)*      &
!                                       becptk_d(jkb,js,ibnd,ik1)*uact(nu+ipol)
                                       becptk_d(jkb,js,ibnd,ik)*uact(nu+ipol)
                                    ENDIF
                                 END DO
                              END DO
                           ELSE
                              DO is=1,npol
                                 IF (isolv==1) THEN
                                    ps1k_nc_d(ikb,is,ibnd,id)=              &
                                    ps1k_nc_d(ikb,is,ibnd,id)+              &
                                    int2_d(ih,jh,ipol,nb,na) *              &
!                                    becp1k_d(jkb,is,ibnd,ik1)*uact(nu+ipol)
                                    becp1k_d(jkb,is,ibnd,ik)*uact(nu+ipol)
                                 ELSE
                                    ps1k_nc_d(ikb,is,ibnd,id)=              &
                                    ps1k_nc_d(ikb,is,ibnd,id)+              &
                                    int2_d(ih,jh,ipol,nb,na) *              &
!                                    becptk_d(jkb,is,ibnd,ik1)*uact(nu+ipol)
                                    becptk_d(jkb,is,ibnd,ik)*uact(nu+ipol)
                                 ENDIF
                              END DO
                           END IF
                        ELSE
                           ps1k_d (ikb, ibnd, id) = ps1k_d (ikb, ibnd, id) + &
                             (int2_d (ih, jh, ipol, nb, na) *                &
!                              becp1k_d(jkb, 1, ibnd, ik1) ) * uact (nu + ipol)
                              becp1k_d(jkb, 1, ibnd, ik) ) * uact (nu + ipol)
                        ENDIF
                     ENDDO
                  ENDIF  ! okvan
               ENDDO ! ipol
            ENDDO ! jh
         ENDDO ! ih
         ijkb0 = ijkb0 + nh (nt)
      ENDIF
   ENDDO  ! na
ENDDO ! nt

RETURN
END SUBROUTINE dvqpsi_us_dev0
!
!------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvqpsi_us_dev1(ndmx, st, ikt, dvpsi, &
                                   npol, nkb, nbnd, nk, npe, nsolv)
!------------------------------------------------------------------------
!
!  This routines add to dvpsi, the first term calculated in 
!  add_dvqpsi_us_only.f90
!
!
USE cudafor
USE kinds, ONLY : DP
USE many_k_ph_mod, ONLY : ps1k_d, ps1k_nc_d
USE many_k_mod, ONLY : vkbk_d, noncolin => noncolin_d
USE klist, ONLY : ngk => ngk_d

IMPLICIT NONE

INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, nkb

INTEGER, DEVICE :: st(nk*npe*nsolv)
INTEGER, DEVICE :: ikt(nk)

COMPLEX(DP), DEVICE :: dvpsi(ndmx*npol,nbnd*nk*npe*nsolv)

INTEGER  :: iv, ig, ik1, ikb, ipol, id, isolv, ipert, isp, kstart, st_, &
            npwq, ikq, ibnd
COMPLEX(DP) :: asum

id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
IF (id>nk*npe*nsolv) RETURN
isp=(id-1)/nk+1
isolv=(isp-1)/npe+1
ipert=MOD(isp-1,npe)+1

ik1=MOD(id-1,nk)+1
ikq=ikt(ik1)
npwq=ngk(ikq)

kstart=nkb*(ik1-1)

st_=st(id)

ig=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
IF (ig>npwq) RETURN

ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ibnd>nbnd) RETURN

IF (noncolin) THEN
   asum=(0.0_DP,0.0_DP)
   DO ikb=1,nkb
      asum=asum + ps1k_nc_d(ikb,1,ibnd,id) * vkbk_d(ig,kstart+ikb) 
   ENDDO
   dvpsi(ig, st_+ibnd)=dvpsi(ig, st_+ibnd) + asum
   asum=(0.0_DP,0.0_DP)
   DO ikb=1,nkb
      asum=asum + ps1k_nc_d(ikb,2,ibnd,id) * vkbk_d(ig,kstart+ikb) 
   ENDDO
   dvpsi(ndmx+ig, st_+ibnd)=dvpsi(ndmx+ig, st_+ibnd) + asum
ELSE
   asum=(0.0_DP,0.0_DP)
   DO ikb=1,nkb
      asum=asum + ps1k_d(ikb,ibnd,id) * vkbk_d(ig,kstart+ikb) 
   ENDDO
   dvpsi(ig, st_+ibnd)=dvpsi(ig, st_+ibnd) + asum
ENDIF

RETURN
END SUBROUTINE dvqpsi_us_dev1
!
!------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE dvqpsi_us_dev2(ndmx, st, ikt, dvpsi, npol, &
                                             nkb, nbnd, nk, npe, nsolv)
!------------------------------------------------------------------------
!
!  This routines add to dvpsi, the second term calculated in 
!  add_dvqpsi_us_only.f90
!
USE cudafor
USE kinds, ONLY : DP
USE many_k_ph_mod, ONLY : ps2k_d, ps2k_nc_d
USE many_k_mod, ONLY : vkbk_d, xk => xk_d, noncolin => noncolin_d, g => g_d
USE klist, ONLY : igk_k => igk_k_d, ngk => ngk_d

IMPLICIT NONE

INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, nkb

INTEGER, DEVICE :: st(nk*npe*nsolv)
INTEGER, DEVICE :: ikt(nk)

COMPLEX(DP), DEVICE :: dvpsi(ndmx*npol,nbnd*nk*npe*nsolv)

REAL(DP), PARAMETER :: eps = 1.D-12
LOGICAL :: ok
INTEGER  :: iv, ig, ik1, ikb, ipol, id, ipert, isolv, isp, kstart, st_, &
            npwq, ikq, ibnd
COMPLEX(DP) :: asum

id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
IF (id>nk*npe*nsolv) RETURN
isp=(id-1)/nk+1
isolv=(isp-1)/npe+1
ipert=MOD(isp-1,npe)+1
ik1=MOD(id-1,nk)+1
ikq=ikt(ik1)
npwq=ngk(ikq)

kstart=nkb*(ik1-1)

st_=st(id)

ig=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
IF (ig>npwq) RETURN

ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
IF (ibnd>nbnd) RETURN

ok=.FALSE.
IF (noncolin) THEN
   DO ikb=1,nkb
      DO ipol=1,3
         ok=ok.OR.(ABS(ps2k_nc_d(ikb,1,ibnd,ipol,id))>eps) &
              .OR.(ABS(ps2k_nc_d(ikb,2,ibnd,ipol,id))>eps)
      ENDDO
   ENDDO
ELSE
   DO ikb=1,nkb
      DO ipol=1,3
         ok=ok.OR.(ABS(ps2k_d(ikb,ibnd,ipol,id))>eps)
      ENDDO
   ENDDO
ENDIF
IF (.NOT.ok) RETURN

iv=igk_k(ig,ikq)
IF (noncolin) THEN
   asum=(0.0_DP,0.0_DP)
   DO ikb=1,nkb
      DO ipol=1,3
         asum=asum+ps2k_nc_d(ikb,1,ibnd,ipol,id) *  &
                      vkbk_d(ig,kstart+ikb) * (xk(ipol, ikq) + g(ipol,iv))
      ENDDO
   ENDDO
   dvpsi(ig, st_+ibnd)=dvpsi(ig, st_+ibnd) + asum
   asum=(0.0_DP,0.0_DP)
   DO ikb=1,nkb
      DO ipol=1,3
         asum=asum+ps2k_nc_d(ikb,2,ibnd,ipol,id) *  &
                      vkbk_d(ig,kstart+ikb) * (xk(ipol, ikq) + g(ipol,iv))
      ENDDO
   ENDDO
   dvpsi(ndmx+ig, st_+ibnd)=dvpsi(ndmx+ig, st_+ibnd) + asum
ELSE
   asum=(0.0_DP,0.0_DP)
   DO ikb=1,nkb
      DO ipol=1,3
         asum=asum+ps2k_d(ikb,ibnd,ipol,id) *  &
                      vkbk_d(ig,kstart+ikb) * (xk(ipol, ikq) + g(ipol,iv))
      ENDDO
   ENDDO
   dvpsi(ig, st_+ibnd)=dvpsi(ig, st_+ibnd) + asum
ENDIF
RETURN
END SUBROUTINE dvqpsi_us_dev2
#endif
