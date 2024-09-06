! Copyright (C) 2023 Andrea Dal Corso 
!                                 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routines derives from the coresponding parts contained
! in Quantum espresso, but computes entirely the product of
! the Hamiltonian with a wavefunction psi on the GPU device.
! It does this in parallel for many k points and bands. 
!
#if defined(__CUDA)
!-------------------------------------------------------------------------
  SUBROUTINE h_s_psik_dev(npwx, outk_d, kdimk_d, npw_d, nveck_d, nb1k_d, &
             st_d, stx_d, ikt_d, ikblk_d, npol, psi, hpsi, spsi, psicmr, &
             nvec, nvecx, nnr, nset, minus_b)
!-------------------------------------------------------------------------
!
!  This subroutine applies H and S to a set of wavefunctions contained in
!  psi. The array psi whose leading dimension is npwx*npol, contains
!  nvecx*nset wavefunctions. For each set the number of psi to calculate
!  is given by nveck_d, while the results is saved in hpsi
!  and spsi. The psi to calculate start at the nb1k_d position and
!  are saved in the same positions in hpsi and spsi.
!  For each set ikt_d is the index of the k point in the
!  complete list of k points, while ikblk_d is the index of the k point in 
!  the currently calculated block of k points. nset is the number of
!  k points to compute in pw.x, but it is the number of k points multiplied
!  the number of perturbations in the phonon.
!  The routine receives also the address psicmr of a large size
!  of device memory sufficient to contain up to nvec*nset smooth FFT meshes.
!  This array must be already allocated on device. 
!  The array outk allows to skip entirely all calculations for some set.
!  If the logical value minus_b is .TRUE. and we are doing a magnetic 
!  noncollinear calculation, nset is divided in two. The first half is 
!  calculated with B, the second set with -B.
!
  USE cudafor
  USE kinds,       ONLY : DP
  USE uspp,        ONLY : nkb
  USE fft_base,    ONLY : dffts
  IMPLICIT NONE

#include<h_psi_interf.f90>

  INTEGER, INTENT(IN) :: npwx
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN) :: nset, nvec, nvecx, nnr, npol
  !! input: the number of sets
  !! input: the number of vectors in each set for psicmr
  !! input: the number of vectors in each set for psi, hpsi, and spsi.
  !! input: the fft dimension.
  !! input: the number of components of each wavefunction.
  LOGICAL, INTENT(IN) :: minus_b
  !! input: if true the second half of sets is computed with -B_xc
  !!        (only noncollinear magnetic case).
  LOGICAL, INTENT(IN), DEVICE :: outk_d(nset)
  !! input: when .TRUE. the set is not calculated
  INTEGER, INTENT(IN), DEVICE :: kdimk_d(nset), npw_d(nset)
  !! input: kdim is npwx*npol in the noncollinear case, or equal to npw_d
  !!        in the collinear case
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: ikt_d(nset)
  !! input: index of k the global list of k points for this set
  INTEGER, INTENT(IN), DEVICE :: ikblk_d(nset)
  !! input: index of k in the current block of k points
  INTEGER, INTENT(IN), DEVICE :: nveck_d(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k_d(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx_d(nset) 
  !! input: start of each set in psi, hpsi, and spsi
  INTEGER, INTENT(IN), DEVICE :: st_d(nset)
  !! input: start of each set in psicmr
  COMPLEX(DP), DEVICE :: psi(npwx*npol, nvecx * nset)
  !! inp: the psi vector
  COMPLEX(DP), DEVICE :: hpsi(npwx*npol, nvecx * nset)
  !! out: the hpsi vector
  COMPLEX(DP), DEVICE :: spsi(npwx*npol, nvecx * nset)
  !! out: the spsi vector
  REAL(DP), DEVICE :: psicmr(2, nnr, npol, nvec * nset)
  !! used for internal purpose: FFT space

  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, adim, ierr

  CALL start_clock('h_psi_dev')
  nr1=dffts%nr1
  nr2=dffts%nr2
  nr3=dffts%nr3
  nr1x=dffts%nr1x
  nr2x=dffts%nr2x
  adim=MAX(nr1, nr2, nr3)
!
!  compute the kinetic energy
!
  CALL h_psi_ke<<<dim3(npwx/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>>&
       ( npwx, outk_d, npw_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol, psi, &
         hpsi, nvecx, nset )
  ierr=cudaDeviceSynchronize()
!
!  compute calbec
!
  CALL h_psi_calbec<<<dim3(nkb/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>>&
       ( npwx, outk_d, npw_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol, psi, &
         nvecx, nset )
  ierr=cudaDeviceSynchronize()
!
!  compute the products of beta and the coefficients of the nonlocal PP
!
  IF (npol==1) THEN
     CALL compute_ps_gpu<<<dim3(nset,nvec,1),dim3(1,1,1)>>>(outk_d, &
                                          nveck_d, ikt_d, nset)
     ierr=cudaDeviceSynchronize()
  ELSE
     CALL compute_ps_nc_gpu<<<dim3(nset,nvec,1),dim3(1,1,1)>>>(outk_d, &
                                nveck_d, nset, minus_b)
     ierr=cudaDeviceSynchronize()
  ENDIF
!
!  initialize s
!
  CALL copy_s_gpu<<<dim3((npwx*npol)/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>>  &
       ( npwx, outk_d, kdimk_d, nveck_d, nb1k_d, stx_d, npol, psi, spsi, &
                                                        nvecx, nset )
 ierr=cudaDeviceSynchronize()
!
! computation of veff * psi
!
!    ... set to zero the auxiliary space psic
!
 CALL vlocpsi_gpu_setpsi_zero<<<dim3(nset,nvec,nnr/64+1),&
                  dim3(1,1,64)>>>( outk_d, nveck_d, st_d, npol, &
                  psicmr, nvec, nnr, nset)
 ierr=cudaDeviceSynchronize()
 !
 !   ... copy psi into psicr
 !
 CALL put_psi_on_grid<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>> &
          ( npwx, outk_d, npw_d, nveck_d, &
            nb1k_d, st_d, stx_d, ikt_d, npol, psi, psicmr, nvec, &
            nvecx, nnr, nset)
 ierr=cudaDeviceSynchronize()
!
!   ... inverse fft in the direction z
!
 CALL fft1inv_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>&
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2, &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... inverse fft in the direction y
!
 CALL fft2inv_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>&
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2, &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... inverse fft in the direction x
!
 CALL fft3inv_dev<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>>&
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2, &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... apply potential to psi
!
 CALL vlocpsi_gpu_vp<<<dim3(nset,nvec,nnr/64+1),dim3(1,1,64)>>> &
      (outk_d, nveck_d, st_d, ikt_d, npol, psicmr, nvec, nnr, nset, minus_b)
 ierr=cudaDeviceSynchronize()
!
!   ... direct fft in the direction x
!
 CALL fft3fwd_dev<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>>&
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2, &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... direct fft in the direction y
!
 CALL fft2fwd_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>&
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2, &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... direct fft in the direction z
!
 CALL fft1fwd_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>&
             (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2, &
             nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... copy of psicmr into h_psi 
!
 CALL add_grid_to_hpsi<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>> &
          ( npwx, outk_d, npw_d, nveck_d, nb1k_d, st_d, stx_d, ikt_d, &
            npol, hpsi, psicmr, nvec, nvecx, nnr, nset) 
 ierr=cudaDeviceSynchronize()
!
!  Add the contribution of the nonlocal pseudopotential to hpsi and the
!  contribution of qq to spsi
!
 CALL add_vnlpsi_us_gpu<<<dim3(npwx/32+1,nvec,nset),dim3(32,1,1)>>> &
       ( npwx, outk_d, npw_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol, &
                 hpsi, spsi, nvecx, nset )
 ierr=cudaDeviceSynchronize()
 CALL stop_clock('h_psi_dev')
 RETURN
END SUBROUTINE h_s_psik_dev

!-------------------------------------------------------------------------
  SUBROUTINE h_psik_dev(npwx, outk_d, kdimk_d, npw_d, nveck_d, nb1k_d, &
             st_d, stx_d, ikt_d, ikblk_d, npol, psi, hpsi, psicmr, nvec, &
             nvecx, nnr, nset, minus_b)
!-------------------------------------------------------------------------
!
!  This subroutine applies H to a set of wavefunctions contained in
!  psi. The array psi whose leading dimension is npwx*npol, contains
!  nvecx*nset wavefunctions. For each set the number of psi to calculate
!  is given by nveck_d, while the results is saved in hpsi.
!  The psi to calculate start at the nb1k_d position and
!  are saved in the same positions in hpsi.
!  For each set ikt_d is the index of the k point in the
!  complete list of k points, while ikblk_d is the index of the k point in 
!  the currently calculated block of k points. nset is the number of
!  k points to compute in pw.x, but it is the number of k points multiplied
!  the number of perturbations in the phonon.
!  The routine receives also the address psicmr of a large size
!  of device memory sufficient to contain up to nvec*nset smooth FFT meshes.
!  This array must be already allocated on device. 
!  The array outk allows to skip entirely all calculations for some sets.
!  If the logical value minus_b is .TRUE. and we are doing a magnetic 
!  noncollinear calculation, nset is divided in two. The first half is 
!  calculated with B_xc, the second set with -B_xc.
!
  USE cudafor
  USE kinds,       ONLY : DP
  USE uspp,        ONLY : nkb
  USE fft_base,    ONLY : dffts
!
!  dffts gives information on the size of the mesh
!
  IMPLICIT NONE

#include<h_psi_interf.f90>

  INTEGER, INTENT(IN) :: npwx
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN) :: nset, nvec, nvecx, nnr, npol
  !! input: the number of sets
  !! input: the number of vectors in each set for psicmr
  !! input: the number of vectors in each set for psi, hpsi, and spsi.
  !! input: the fft dimension.
  !! input: the number of components of each wavefunction.
  LOGICAL, INTENT(IN) :: minus_b
  !! input: if true the second half of sets is computed with -B_xc
  !!        (only noncollinear magnetic case).
  LOGICAL,INTENT(IN), DEVICE :: outk_d(nset)
  !! input: when .TRUE. the k point is not calculated
  INTEGER, INTENT(IN), DEVICE :: kdimk_d(nset), npw_d(nset)
  !! input: kdim is npwx*npol in the noncollinear case, or equal to npw_d
  !!        in the collinear case
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: ikt_d(nset)
  !! input: index of k the global list of k points for this set
  INTEGER, INTENT(IN), DEVICE :: ikblk_d(nset)
  !! input: index of k in the current block of k points
  INTEGER, INTENT(IN), DEVICE :: nveck_d(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k_d(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx_d(nset)
  !! input: start of each set in psi and hpsi
  INTEGER, INTENT(IN), DEVICE :: st_d(nset)
  !! input: start of each set in psicmr
  COMPLEX(DP), DEVICE :: psi(npwx*npol, nvecx * nset)
  !! inp: the psi vector
  COMPLEX(DP), DEVICE :: hpsi(npwx*npol, nvecx * nset)
  !! out: the hpsi vector
  REAL(DP), DEVICE :: psicmr(2, nnr, npol, nvec * nset)
  !! used for internal purpose: FFT space

  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, adim, ierr

  CALL start_clock('h_psi_dev')
  nr1=dffts%nr1
  nr2=dffts%nr2
  nr3=dffts%nr3
  nr1x=dffts%nr1x
  nr2x=dffts%nr2x
  adim=MAX(nr1, nr2, nr3)
!
!  compute the kinetic energy
!
  CALL h_psi_ke<<<dim3(npwx/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>>&
       ( npwx, outk_d, kdimk_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol, psi, &
         hpsi, nvecx, nset )
  ierr=cudaDeviceSynchronize()
!
!  compute calbec
!
  CALL h_psi_calbec<<<dim3(nkb/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>>&
       ( npwx, outk_d, npw_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol, psi, &
         nvecx, nset )
  ierr=cudaDeviceSynchronize()
!
!  compute the products of beta and the coefficients of the nonlocal PP
!
  IF (npol==1) THEN
     CALL compute_ps_gpu<<<dim3(nset,nvec,1),dim3(1,1,1)>>>(outk_d, &
                                       nveck_d, ikt_d, nset)
     ierr=cudaDeviceSynchronize()
  ELSE
     CALL compute_ps_nc_gpu<<<dim3(nset,nvec,1),dim3(1,1,1)>>>(outk_d, &
                                nveck_d, nset, minus_b)
     ierr=cudaDeviceSynchronize()
  ENDIF
!
! computation of veff * psi
!
!    ... set to zero the auxiliary space psic
!
 CALL vlocpsi_gpu_setpsi_zero<<<dim3(nset/32+1,nvec,nnr/32+1),    &
                  dim3(32,1,32)>>>( outk_d, nveck_d, st_d, npol,  &
                  psicmr, nvec, nnr, nset)
 ierr=cudaDeviceSynchronize()
 !
 !   ... copy psi into psicr
 !
 CALL put_psi_on_grid<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>>  &
          ( npwx, outk_d, npw_d, nveck_d, nb1k_d, st_d, stx_d,     &
            ikt_d, npol, psi, psicmr, nvec, nvecx, nnr, nset)
 ierr=cudaDeviceSynchronize()
!
!   ... inverse fft in the direction z
!
 CALL fft1inv_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>        &
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2,   &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... inverse fft in the direction y
!
 CALL fft2inv_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>        &
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2,   &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... inverse fft in the direction x
!
 CALL fft3inv_dev<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>>          &
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2,   &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... apply potential to psi
!
 CALL vlocpsi_gpu_vp<<<dim3(nset/32+1,nvec,nnr/32+1),dim3(32,1,32)>>> &
      (outk_d, nveck_d, st_d, ikt_d, npol, psicmr, nvec, nnr, nset, minus_b)
 ierr=cudaDeviceSynchronize()
!
!   ... direct fft in the direction x
!
 CALL fft3fwd_dev<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>>           &
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2,    &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... direct fft in the direction y
!
 CALL fft2fwd_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>         &
            (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2,    &
            nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... direct fft in the direction z
!
 CALL fft1fwd_dev<<<dim3(nset,nvec,nr1/32+1),dim3(1,1,32)>>>         &
             (outk_d, nveck_d, st_d, npol, psicmr, nvec, nr1, nr2,   &
             nr3, nr1x, nr2x, nnr, nset, adim)
 ierr=cudaDeviceSynchronize()
!
!   ... copy of psicmr into h_psi 
!
 CALL add_grid_to_hpsi<<<dim3(nset,nvec/32+1,1),dim3(1,32,1)>>>      &
          (npwx, outk_d, npw_d, nveck_d, nb1k_d, st_d, stx_d, ikt_d, &
           npol, hpsi, psicmr, nvec, nvecx, nnr, nset)
 ierr=cudaDeviceSynchronize()
!
!  Add the contribution of the nonlocal pseudopotential to hpsi and the
!  contribution of qq to spsi
!
 CALL add_vnlpsi_gpu<<<dim3(npwx/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>> &
       ( npwx, outk_d, npw_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol,     &
         hpsi, nvecx, nset)
 ierr=cudaDeviceSynchronize()
 CALL stop_clock('h_psi_dev')
 RETURN
END SUBROUTINE h_psik_dev
!
!-------------------------------------------------------------------------
  SUBROUTINE s_psik_dev(npwx, outk_d, kdimk_d, npw_d, nveck_d, nb1k_d, &
             st_d, stx_d, ikblk_d, npol, psi, spsi, nvec, nvecx, nset)
!-------------------------------------------------------------------------
!
!  This subroutine applies S to a set of wavefunctions contained in
!  psi. The array psi whose leading dimension is npwx*npol, contains
!  nvecx*nset wavefunctions. For each set the number of psi to calculate
!  is given by nveck_d, while the results is saved in spsi. The 
!  psi to calculate start at the nb1k_d position and
!  are saved in the same positions in spsi.
!  while ikblk_d is the index of the k point in 
!  the currently calculated block of k points. nset is the number of
!  k points to compute pw.x, but it is the number of k points multiplied
!  the number of perturbations in the phonon.
!  The array outk allows to skip entirely all calculations for some set.
!
  USE cudafor
  USE kinds,       ONLY : DP
  USE uspp,        ONLY : nkb, okvan
  IMPLICIT NONE

#include<h_psi_interf.f90>

  INTEGER, INTENT(IN) :: npwx
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN) :: nset, nvec, nvecx, npol
  !! input: the number of sets
  !! input: the number of vectors in each set for psicmr
  !! input: the number of vectors in each set for psi, hpsi, and spsi.
  !! input: the number of components of each wavefunction.
  LOGICAL, INTENT(IN), DEVICE :: outk_d(nset)
  !! input: when .TRUE. the set is not calculated
  INTEGER, INTENT(IN), DEVICE :: kdimk_d(nset), npw_d(nset)
  !! input: kdim is npwx*npol in the noncollinear case, or equal to npw_d
  !!        in the collinear case
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: ikblk_d(nset)
  !! input: index of k in the current block of k points
  INTEGER, INTENT(IN), DEVICE :: nveck_d(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k_d(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx_d(nset) 
  !! input: start of each set in psi, hpsi, and spsi
  INTEGER, INTENT(IN), DEVICE :: st_d(nset)
  !! input: start of each set in psicmr
  COMPLEX(DP), DEVICE :: psi(npwx*npol, nvecx * nset)
  !! inp: the psi vector
  COMPLEX(DP), DEVICE :: spsi(npwx*npol, nvecx * nset)
  !! out: the spsi vector

  INTEGER :: ierr

  CALL start_clock('s_psi_dev')
!
!  initialize s
!
  CALL copy_s_gpu<<<dim3((npwx*npol)/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>>  &
       ( npwx, outk_d, kdimk_d, nveck_d, nb1k_d, stx_d, npol, psi, spsi, &
                                                        nvecx, nset )
 ierr=cudaDeviceSynchronize()
 IF (.NOT.okvan) THEN
    CALL stop_clock('s_psi_dev')
    RETURN
 ENDIF
!
!  compute calbec
!
  CALL h_psi_calbec<<<dim3(nkb/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>>&
       ( npwx, outk_d, npw_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol, psi, &
         nvecx, nset )
  ierr=cudaDeviceSynchronize()
!
!  compute the products of beta and the coefficients of the nonlocal PP
!
  IF (npol==1) THEN
     CALL compute_ps_s_gpu<<<dim3(nset,nvec,1),dim3(1,1,1)>>>(outk_d, &
                                          nveck_d, nset)
     ierr=cudaDeviceSynchronize()
  ELSE
     CALL compute_ps_s_nc_gpu<<<dim3(nset,nvec,1),dim3(1,1,1)>>>(outk_d, &
                                nveck_d, nset)
     ierr=cudaDeviceSynchronize()
  ENDIF
 !
 CALL add_spsi_gpu<<<dim3(npwx/4+1,nvec/4+1,nset/32+1),dim3(4,4,32)>>> &
       ( npwx, outk_d, npw_d, nveck_d, nb1k_d, stx_d, ikblk_d, npol, &
                 spsi, nvecx, nset )
 ierr=cudaDeviceSynchronize()
 CALL stop_clock('s_psi_dev')
 RETURN
END SUBROUTINE s_psik_dev
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE h_psi_ke( lda, outk, npw, nveck, nb1k, stx, &
                                   ikblk, npol, psi_d, hpsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  ! This routine applies the kinetic energy to the wavefunctions
  ! It runs in parallel on the GPU. Each thread computes one set,
  ! one band and one G vector.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY: g2kink_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  !! input: the number of sets
  !! input: the dimension of each sets in psi and hpsi
  LOGICAL,INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  !! input: the index of the k point of this set in the block
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in psi and hpsi
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components of each wavefunction.
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx * nset)
  !! inp: the psi vector
  COMPLEX(DP), DEVICE :: hpsi_d(lda*npol, nvecx * nset)
  !! out: the hpsi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, ik1, i, k, n, mb, nb1, stx_
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  ik1=ikblk(ik)
  n = npw(ik)
  nb1=nb1k(ik)
  stx_=stx(ik)
  mb=nveck(ik) 

  i=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i>lda) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k = k0 + stx_ + nb1 - 1

  IF (i <= n) THEN
     hpsi_d (i, k) = g2kink_d (i,ik1) * psi_d (i, k)
  ELSE
     hpsi_d (i, k) = (0.0_DP, 0.0_DP)
  END IF
  IF ( npol==2 ) THEN
     IF (i <= n) THEN
        hpsi_d (lda+i, k) = g2kink_d (i,ik1) * psi_d (lda+i, k)
     ELSE
        hpsi_d (lda+i, k) = (0.0_DP, 0.0_DP)
     END IF
  ENDIF
  !
  RETURN
  !
END SUBROUTINE h_psi_ke
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE h_psi_calbec( lda, outk, npw, nveck, &
                             nb1k, stx, ikblk, npol, psi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  ! This routine computes becp for all sets.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : becpk_d, vkbk_d, nkb => nkb_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  !! input: the number of sets
  !! input: the number of vectors in each set for psi
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  !! input: the index of the k point of this set in the block
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in psi 
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx * nset)
  !! input: the input wavefunctions.
  !
  !  ... local variables
  !
  INTEGER :: k0, i0, ik, ik1, i, j, k, n, m, nb1, stx_
  !
  COMPLEX(DP) :: asum
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  ik1=ikblk(ik)
  n = npw(ik)
  stx_=stx(ik)
  nb1=nb1k(ik)
  m=nveck(ik) 

  i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i0>nkb) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>m) RETURN
  i=i0 + nkb*(ik1-1)
  k = k0 + stx_ + nb1 - 1

  asum=(0.0_DP,0.0_DP)
  DO j=1, n
     asum=asum+CONJG(vkbk_d(j,i))*psi_d(j,k)
  ENDDO
  becpk_d(i0,1,k0,ik)=asum

  IF (npol==2) THEN
     asum=(0.0_DP,0.0_DP)
     DO j=1, n
        asum=asum+CONJG(vkbk_d(j,i))*psi_d(lda+j,k)
     ENDDO
     becpk_d(i0,2,k0,ik)=asum
  ENDIF
  !
  RETURN
  !
END SUBROUTINE h_psi_calbec
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE copy_s_gpu( lda, outk, kdimk, nveck, &
                        nb1k, stx, npol, psi_d, spsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  ! This routine copies the function psi in s_psi for all sets.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi_d and spsi_d
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  !! input: the number of sets
  !! input: the number of vectors in each set for psi
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: kdimk(nset)
  !! input: kdimk is npwx*npol in the noncollinear case, and
  !!        the number of plane waves in the other cases.
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in psi_d and spsi_d
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  COMPLEX(DP), DEVICE :: psi_d(lda*npol, nvecx*nset)
  !! inp: the psi vector
  COMPLEX(DP), DEVICE :: spsi_d(lda*npol, nvecx*nset)
  !! output: the spsi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, i0, ik, k, n, m, nb1, stx_
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  n = kdimk(ik)
  stx_=stx(ik)
  nb1=nb1k(ik)
  m=nveck(ik) 

  i0=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i0>n) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>m) RETURN
  k = k0 + stx_ + nb1 - 1

  spsi_d(i0,k)=psi_d(i0,k)
  !
  RETURN
  !
END SUBROUTINE copy_s_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE add_vnlpsi_us_gpu( lda, outk, npw, nveck, &
                     nb1k, stx, ikblk, npol, hpsi_d, spsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  ! This routine adds the contribution of the non local part to hpsi and 
  ! spsi. It assumes that psk_d and pssk_d are already filled with the
  ! sum_n becp(n,ibnd) * D_mn and sum_n becp(n,ibnd) * q_mn respectively.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : psk_d, pssk_d, vkbk_d, nkb => nkb_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of hpsi_d and spsi_d
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  !! input: the number of sets
  !! input: the number of vectors in each set for psi
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  !! input: the index of the k point of this set in the block
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in hpsi_d and spsi_d
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  COMPLEX(DP), DEVICE :: hpsi_d(lda*npol, nvecx*nset)
  !! in/out: the hpsi vector
  COMPLEX(DP), DEVICE :: spsi_d(lda*npol, nvecx*nset)
  !! in/out: the spsi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, ik1, i, j, j1, k, n, mb, nb1, stx_, sh
  !
  COMPLEX(DP) :: asum
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  ik1=ikblk(ik)
  n = npw(ik)
  stx_=stx(ik)
  nb1=nb1k(ik)
  mb=nveck(ik) 

  i=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i>n) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k = k0 + stx_ + nb1 - 1
  sh= nkb*(ik1-1)

  asum=(0.0_DP,0.0_DP)
  DO j=1, nkb
     j1=j+sh
     asum=asum+vkbk_d(i,j1)*psk_d(j,1,k0,ik)
  ENDDO
  hpsi_d(i,k)=hpsi_d(i,k) + asum

  IF (npol==2) THEN
     asum=(0.0_DP,0.0_DP)
     DO j=1, nkb
        j1=j+sh
        asum=asum+vkbk_d(i,j1)*psk_d(j,2,k0,ik)
     ENDDO
     hpsi_d(lda+i,k)=hpsi_d(lda+i,k) + asum
  ENDIF

  asum=(0.0_DP,0.0_DP)
  DO j=1, nkb
     j1=j+sh
     asum=asum+vkbk_d(i,j1)*pssk_d(j,1,k0,ik)
  ENDDO
  spsi_d(i,k)=spsi_d(i,k) + asum
  !
  IF (npol==2) THEN
     asum=(0.0_DP,0.0_DP)
     DO j=1, nkb
        j1=j+sh
        asum=asum+vkbk_d(i,j1)*pssk_d(j,2,k0,ik)
     ENDDO
     spsi_d(lda+i,k)=spsi_d(lda+i,k) + asum
  ENDIF
  !
  RETURN
  !
END SUBROUTINE add_vnlpsi_us_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE add_vnlpsi_gpu( lda, outk, npw, nveck, &
                     nb1k, stx, ikblk, npol, hpsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  ! This routine adds the contribution of the non local part to hpsi 
  ! It assumes that psk_d are already filled with sum_n becp(n,ibnd) * D_mn.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : psk_d, vkbk_d, nkb => nkb_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of psi
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  !! input: the number of sets
  !! input: the number of vectors in each set for hpsi_d
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  !! input: the index of the k point of this set in the block
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in hpsi_d and spsi_d
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  COMPLEX(DP), DEVICE :: hpsi_d(lda*npol, nvecx*nset)
  !! in/out: the hpsi_d vector
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, ik1, i, j, j1, k, n, mb, nb1, stx_, sh
  !
  COMPLEX(DP) :: asum
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  ik1=ikblk(ik)
  n =npw(ik)
  stx_=stx(ik)
  nb1=nb1k(ik)
  mb=nveck(ik) 

  i=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i>n) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k = k0 + stx_ + nb1 - 1
  sh=nkb * (ik1-1)

  asum=(0.0_DP,0.0_DP)
  DO j=1, nkb
     j1=j+sh
     asum=asum+vkbk_d(i,j1)*psk_d(j,1,k0,ik)
  ENDDO
  hpsi_d(i,k)=hpsi_d(i,k) + asum
  
  IF (npol==2) THEN
     asum=(0.0_DP,0.0_DP)
     DO j=1, nkb
        j1=j+sh
        asum=asum+vkbk_d(i,j1)*psk_d(j,2,k0,ik)
     ENDDO
     hpsi_d(lda+i,k)=hpsi_d(lda+i,k) + asum
  ENDIF
  !
  RETURN
  !
END SUBROUTINE add_vnlpsi_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE add_spsi_gpu( lda, outk, npw, nveck, &
                     nb1k, stx, ikblk, npol, spsi_d, nvecx, nset )
  !-----------------------------------------------------------------------
  !
  ! This routine adds the contribution of the non local part to  
  ! spsi. It assumes that pssk_d is already filled with the
  ! sum_n becp(n,ibnd) * q_mn.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : pssk_d, vkbk_d, nkb => nkb_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: lda
  !! input: the leading dimension of spsi_d
  INTEGER, INTENT(IN), VALUE :: nset, nvecx
  !! input: the number of sets
  !! input: the number of vectors in each set for psi
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: ikblk(nset)
  !! input: the index of the k point of this set in the block
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves for each set
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in hpsi_d and spsi_d
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  COMPLEX(DP), DEVICE :: spsi_d(lda*npol, nvecx*nset)
  !! in/out: the spsi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, ik1, i, j, j1, k, n, mb, nb1, stx_, sh
  !
  COMPLEX(DP) :: asum
  !
  ik=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  ik1=ikblk(ik)
  n = npw(ik)
  stx_=stx(ik)
  nb1=nb1k(ik)
  mb=nveck(ik) 

  i=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (i>n) RETURN
  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k = k0 + stx_ + nb1 - 1
  sh= nkb*(ik1-1)

  asum=(0.0_DP,0.0_DP)
  DO j=1, nkb
     j1=j+sh
     asum=asum+vkbk_d(i,j1)*pssk_d(j,1,k0,ik)
  ENDDO
  spsi_d(i,k)=spsi_d(i,k) + asum
  !
  IF (npol==2) THEN
     asum=(0.0_DP,0.0_DP)
     DO j=1, nkb
        j1=j+sh
        asum=asum+vkbk_d(i,j1)*pssk_d(j,2,k0,ik)
     ENDDO
     spsi_d(lda+i,k)=spsi_d(lda+i,k) + asum
  ENDIF
  !
  RETURN
  !
END SUBROUTINE add_spsi_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE compute_ps_gpu( outk, nveck, ikt, nset)
  !-----------------------------------------------------------------------
  !
  ! This routines fills psk_d and pssk_d with sum_n becp(n,ibnd) * D_mn 
  ! and sum_n becp(n,ibnd) * q_mn respectively.
  ! pssk_d is computed only when uspp is .TRUE..
  ! This routine is called only in the collinear case.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : psk_d, pssk_d, becpk_d, isk_d, nat=>nat_d, &
                             ntyp=>ntyp_d, ityp=>ityp_d, nh=>nh_d,      &
                             qq_at=>qq_at_d, deeq=>deeq_d, lsda => lsda_d, &
                             nkb => nkb_d, okvan => okvan_d

  USE uspp, ONLY : ofsbeta=>ofsbeta_d
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nset
  !! input: the number of sets
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)
  !! input: index of k in the global list of k points for this set
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, ik2, m, k, j, isk_, nt, na, nhnt
  !
  COMPLEX(DP) :: asum, bsum
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  IF (nkb == 0) RETURN
  m=nveck(ik) 
  ik2=ikt(ik)
  isk_=1
  IF (lsda) isk_=isk_d(ik2)

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>m) RETURN

  psk_d(:,1,k0,ik)=(0.0_DP,0.0_DP)
  pssk_d(:,1,k0,ik)=(0.0_DP,0.0_DP)

  DO nt = 1, ntyp
     !
     IF ( nh(nt) == 0 ) CYCLE
     !
     nhnt = nh(nt)
     !
     DO na = 1, nat
        !
        IF ( ityp(na) == nt ) THEN
           !
           DO k=1, nhnt
              asum=(0.0_DP,0.0_DP)
              bsum=(0.0_DP,0.0_DP)
              DO j=1, nhnt
                 asum=asum+deeq(k,j,na,isk_)*becpk_d(ofsbeta(na)+j,1,k0,ik)
                 IF (okvan) bsum=bsum+qq_at(k,j,na)* &
                                     becpk_d(ofsbeta(na)+j,1,k0,ik)
              ENDDO
              psk_d(ofsbeta(na)+k,1,k0,ik)=psk_d(ofsbeta(na)+k,1,k0,ik)+asum
              IF (okvan) pssk_d(ofsbeta(na)+k,1,k0,ik)=&
                        pssk_d(ofsbeta(na)+k,1,k0,ik)+bsum
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE compute_ps_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE compute_ps_nc_gpu(outk, nveck, nset, minus_b)
  !-----------------------------------------------------------------------
  !
  ! This routines fills psk_d and pssk_d with sum_n becp(n,ibnd) * D_mn 
  ! and sum_n becp(n,ibnd) * q_mn respectively.
  ! pssk_d is computed only when uspp is .TRUE..
  ! This routine is called only in the noncollinear case.
  !
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : psk_d, pssk_d, becpk_d, nat=>nat_d,        &
                             ntyp=>ntyp_d, ityp=>ityp_d, nh=>nh_d,      &
                             qq_at=>qq_at_d, deeq_nc=>deeq_nc_d,        &
                             qq_so=>qq_so_d, nkb => nkb_d,              &
                             lspinorb => lspinorb_d, okvan => okvan_d
  USE many_k_ph_mod,  ONLY:  deeq_nc_save=> deeq_nc_save_d
  USE uspp, ONLY : ofsbeta=>ofsbeta_d
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN), VALUE :: minus_b
  !! input: if true the second set of psi is calculate with minus B_xc
  INTEGER, INTENT(IN), VALUE :: nset
  !! input: the number of sets
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, mb, k, j, nt, na, nhnt, iset
  !
  COMPLEX(DP) :: asum, bsum
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  iset=1
  IF (minus_b.AND.(ik>nset/2)) iset=2
  IF (outk(ik)) RETURN
  IF (nkb == 0) RETURN
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN

  psk_d(:,:,k0,ik)=(0.0_DP,0.0_DP)
  pssk_d(:,:,k0,ik)=(0.0_DP,0.0_DP)

  DO nt = 1, ntyp
     !
     IF ( nh(nt) == 0 ) CYCLE
     !
     nhnt = nh(nt)
     !
     DO na = 1, nat
        !
        IF ( ityp(na) == nt ) THEN
           !
           DO k=1, nhnt
              asum=(0.0_DP,0.0_DP)
              IF (minus_b) THEN
                 DO j=1, nhnt
                    asum=asum+deeq_nc_save(k,j,na,1,iset)*      &
                              becpk_d(ofsbeta(na)+j,1,k0,ik)    &
                             +deeq_nc_save(k,j,na,2,iset)*      &
                              becpk_d(ofsbeta(na)+j,2,k0,ik)
                 ENDDO
              ELSE
                 DO j=1, nhnt
                    asum=asum+deeq_nc(k,j,na,1)*becpk_d(ofsbeta(na)+j,1,k0,ik)&
                             +deeq_nc(k,j,na,2)*becpk_d(ofsbeta(na)+j,2,k0,ik)
                 ENDDO
              ENDIF
              psk_d(ofsbeta(na)+k,1,k0,ik)=psk_d(ofsbeta(na)+k,1,k0,ik)+asum
              IF (okvan) THEN
                 bsum=(0.0_DP,0.0_DP)
                 IF (lspinorb) THEN
                    DO j=1, nhnt
                       bsum = bsum                                           &
                            + qq_so(k,j,1,nt)*becpk_d(ofsbeta(na)+j,1,k0,ik) &
                            + qq_so(k,j,2,nt)*becpk_d(ofsbeta(na)+j,2,k0,ik)
                    ENDDO
                 ELSE
                    DO j=1, nhnt
                       bsum=bsum+qq_at(k,j,na)*becpk_d(ofsbeta(na)+j,1,k0,ik) 
                    ENDDO
                 ENDIF
                 pssk_d(ofsbeta(na)+k,1,k0,ik)=pssk_d(ofsbeta(na)+k,1,k0,ik)&
                                               +bsum
              ENDIF
           ENDDO
           DO k=1, nhnt
              asum=(0.0_DP,0.0_DP)
              IF (minus_b) THEN
                 DO j=1, nhnt
                    asum=asum + deeq_nc_save(k,j,na,3,iset)*             &
                                becpk_d(ofsbeta(na)+j,1,k0,ik) &
                              + deeq_nc_save(k,j,na,4,iset)*             &
                                becpk_d(ofsbeta(na)+j,2,k0,ik)
                 ENDDO
              ELSE
                 DO j=1, nhnt
                    asum=asum + deeq_nc(k,j,na,3)*             &
                                becpk_d(ofsbeta(na)+j,1,k0,ik) &
                              + deeq_nc(k,j,na,4)*             &
                                becpk_d(ofsbeta(na)+j,2,k0,ik)
                 ENDDO
              ENDIF
              psk_d(ofsbeta(na)+k,2,k0,ik)=psk_d(ofsbeta(na)+k,2,k0,ik)+asum
              IF (okvan) THEN
                 bsum=(0.0_DP,0.0_DP)
                 IF (lspinorb) THEN
                    DO j=1, nhnt
                       bsum=bsum                                             &
                             +qq_so(k,j,3,nt)*becpk_d(ofsbeta(na)+j,1,k0,ik) &
                             +qq_so(k,j,4,nt)*becpk_d(ofsbeta(na)+j,2,k0,ik)
                    ENDDO
                 ELSE
                    DO j=1, nhnt
                       bsum=bsum+qq_at(k,j,na)*becpk_d(ofsbeta(na)+j,2,k0,ik) 
                    ENDDO
                 ENDIF
                 pssk_d(ofsbeta(na)+k,2,k0,ik)=pssk_d(ofsbeta(na)+k,2,k0,ik)&
                                               +bsum
              ENDIF
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE compute_ps_nc_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE compute_ps_s_gpu( outk, nveck, nset)
  !-----------------------------------------------------------------------
  !
  ! This routines fills pssk_d with sum_n becp(n,ibnd) * q_mn 
  ! respectively.
  ! This routine is called only in the collinear case.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : pssk_d, becpk_d, isk_d, nat=>nat_d, &
                             ntyp=>ntyp_d, ityp=>ityp_d, nh=>nh_d,      &
                             qq_at=>qq_at_d, deeq=>deeq_d, lsda => lsda_d, &
                             nkb => nkb_d, okvan => okvan_d

  USE uspp, ONLY : ofsbeta=>ofsbeta_d
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nset
  !! input: the number of sets
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, m, k, j, nt, na, nhnt
  !
  COMPLEX(DP) :: asum, bsum
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  IF (nkb == 0) RETURN
  m=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>m) RETURN

  pssk_d(:,1,k0,ik)=(0.0_DP,0.0_DP)

  DO nt = 1, ntyp
     !
     IF ( nh(nt) == 0 ) CYCLE
     !
     nhnt = nh(nt)
     !
     DO na = 1, nat
        !
        IF ( ityp(na) == nt ) THEN
           !
           DO k=1, nhnt
              bsum=(0.0_DP,0.0_DP)
              DO j=1, nhnt
                 bsum=bsum+qq_at(k,j,na)* &
                                     becpk_d(ofsbeta(na)+j,1,k0,ik)
              ENDDO
              pssk_d(ofsbeta(na)+k,1,k0,ik)=&
                        pssk_d(ofsbeta(na)+k,1,k0,ik)+bsum
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE compute_ps_s_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE compute_ps_s_nc_gpu(outk, nveck, nset)
  !-----------------------------------------------------------------------
  !
  ! This routines fills psk_d and pssk_d with sum_n becp(n,ibnd) * D_mn 
  ! and sum_n becp(n,ibnd) * q_mn respectively.
  ! pssk_d is computed only when uspp is .TRUE..
  ! This routine is called only in the noncollinear case.
  !
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : psk_d, pssk_d, becpk_d, nat=>nat_d,        &
                             ntyp=>ntyp_d, ityp=>ityp_d, nh=>nh_d,      &
                             qq_at=>qq_at_d, deeq_nc=>deeq_nc_d,        &
                             qq_so=>qq_so_d, nkb => nkb_d,              &
                             lspinorb => lspinorb_d, okvan => okvan_d
  USE many_k_ph_mod,  ONLY:  deeq_nc_save=> deeq_nc_save_d
  USE uspp, ONLY : ofsbeta=>ofsbeta_d
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nset
  !! input: the number of sets
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  !
  !  ... local variables
  !
  INTEGER :: k0, ik, mb, k, j, nt, na, nhnt, iset
  !
  COMPLEX(DP) :: bsum
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  IF (nkb == 0) RETURN
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN

  pssk_d(:,:,k0,ik)=(0.0_DP,0.0_DP)

  DO nt = 1, ntyp
     !
     IF ( nh(nt) == 0 ) CYCLE
     !
     nhnt = nh(nt)
     !
     DO na = 1, nat
        !
        IF ( ityp(na) == nt ) THEN
           !
           DO k=1, nhnt
              bsum=(0.0_DP,0.0_DP)
              IF (lspinorb) THEN
                 DO j=1, nhnt
                    bsum = bsum                                           &
                         + qq_so(k,j,1,nt)*becpk_d(ofsbeta(na)+j,1,k0,ik) &
                         + qq_so(k,j,2,nt)*becpk_d(ofsbeta(na)+j,2,k0,ik)
                 ENDDO
              ELSE
                 DO j=1, nhnt
                    bsum=bsum+qq_at(k,j,na)*becpk_d(ofsbeta(na)+j,1,k0,ik) 
                 ENDDO
              ENDIF
              pssk_d(ofsbeta(na)+k,1,k0,ik)=pssk_d(ofsbeta(na)+k,1,k0,ik)&
                                               +bsum
           ENDDO
           DO k=1, nhnt
              bsum=(0.0_DP,0.0_DP)
              IF (lspinorb) THEN
                 DO j=1, nhnt
                    bsum=bsum                                             &
                          +qq_so(k,j,3,nt)*becpk_d(ofsbeta(na)+j,1,k0,ik) &
                          +qq_so(k,j,4,nt)*becpk_d(ofsbeta(na)+j,2,k0,ik)
                 ENDDO
              ELSE
                 DO j=1, nhnt
                    bsum=bsum+qq_at(k,j,na)*becpk_d(ofsbeta(na)+j,2,k0,ik) 
                 ENDDO
              ENDIF
              pssk_d(ofsbeta(na)+k,2,k0,ik)=pssk_d(ofsbeta(na)+k,2,k0,ik)&
                                            +bsum
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE compute_ps_s_nc_gpu
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE vlocpsi_gpu_vp(outk, nveck, st, ikt, npol, &
                              psicr, nvec, nnr, nset, minus_b)
!----------------------------------------------------------------------- 
!  
  ! This routine applies the local potential to the wavefunctions 
  ! in parallel for all sets, bands and real space points.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE scf_gpum,       ONLY : vrs_d
  USE many_k_mod,     ONLY : isk => isk_d, lsda => lsda_d, domag => domag_d, &
                             noncolin => noncolin_d
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN), VALUE :: minus_b
  !! if .TRUE. in the noncollinear magnetic case applies -B_xc to the 
  !! second half of the sets
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nvec
  !! input: the number of sets
  !! input: the number of points of the smooth mesh
  !! input: the number of vectors for each set in psicr
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: start of each set in psicr
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)
  !! input: index of k in the global list of k points for this set
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnr, npol, nvec * nset)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: ik, ik2, k0, k2, i, mb, ipol, is_, st_

  REAL(DP) :: bsign

  COMPLEX(DP) :: aux1, aux2, sup, sdwn
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  bsign=1.0_DP
  IF (minus_b.AND.ik>nset/2) bsign=-1.0_DP
  st_=st(ik)
  ik2=ikt(ik)
  mb=nveck(ik) 
  is_=1
  IF (lsda) is_=isk(ik2)

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  i=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (i>nnr) RETURN

  k2 = k0 + st_ 

  IF (noncolin.AND.domag) THEN
     aux1 = CMPLX(psicr(1,i,1,k2),psicr(2,i,1,k2),KIND=DP)
     aux2 = CMPLX(psicr(1,i,2,k2),psicr(2,i,2,k2),KIND=DP)
     sup  = aux1 * (vrs_d(i,1)+bsign*vrs_d(i,4)) + &
            aux2 * (vrs_d(i,2)-(0.d0,1.d0)*vrs_d(i,3)) * bsign
     sdwn = aux2 * (vrs_d(i,1)-bsign*vrs_d(i,4)) + &
            aux1 * (vrs_d(i,2)+(0.d0,1.d0)*vrs_d(i,3)) * bsign
     psicr(1,i,1,k2) = DBLE(sup)
     psicr(2,i,1,k2) = AIMAG(sup)
     psicr(1,i,2,k2) = DBLE(sdwn)
     psicr(2,i,2,k2) = AIMAG(sdwn)
  ELSE
     DO ipol=1, npol
        psicr(1,i,ipol,k2) = psicr(1,i,ipol,k2) * vrs_d(i,is_)
        psicr(2,i,ipol,k2) = psicr(2,i,ipol,k2) * vrs_d(i,is_)
     ENDDO
  ENDIF
  RETURN
  !
END SUBROUTINE vlocpsi_gpu_vp
!
!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE vlocpsi_gpu_setpsi_zero( outk, nveck, st, &
                              npol, psicr, nvec, nnr, nset) 
!----------------------------------------------------------------------- 
!  
  !
  !  This routine sets to zero the psicr variable. It is
  !  parallelized on the GPU on the sets, bands, and fft mesh points.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nvec
  !! input: the number of sets
  !! input: the number of points of the smooth mesh
  !! input: the number of vectors for each set in psicr
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: start of each set in psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnr, npol, nvec*nset)
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, ik, i, mb, st_, ipol
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  st_=st(ik)
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  i=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (i>nnr) RETURN
!
!  zero all the psicr 
!
  DO ipol=1, npol
     psicr(1,i,ipol,k2)=0.0_DP
     psicr(2,i,ipol,k2)=0.0_DP
  ENDDO

  RETURN
  END SUBROUTINE vlocpsi_gpu_setpsi_zero

!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE put_psi_on_grid( lda, outk, npw, nveck, &
                nb1k, st, stx, ikt, npol, psi_d, psicr, nvec, &
                nvecx, nnr, nset) 
!----------------------------------------------------------------------- 
!  
  ! This routine distributes psi_d on the fft mesh psicr.
  ! It is parallelized for the GPU on the sets and on the bands.
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
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nvecx, nvec
  !! input: the number of sets
  !! input: the number of points of the smooth mesh
  !! input: the number of vectors for each set in psi_d
  !! input: the number of vectors for each set in psicr
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)
  !! input: index of k in the global list of k points for this set
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves. 
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute. 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: start of each set in psicr
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in psi_d
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each psi
  COMPLEX(DP), DEVICE, INTENT(IN) :: psi_d(lda*npol, nvecx * nset)
  !! input: the wavefunctions to distribute
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnr, npol, nvec * nset)
  !! output: the fft mesh with the distributed functions. psicr is assumed
  !!         to contain zeros everywhere.
  !
  !  ... local variables
  !
  INTEGER :: k0, k1, k2, ik, ik1, i, np, mb, nb1,  &
             stx_, st_, iv
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  np = npw(ik)
  stx_=stx(ik)
  st_=st(ik)
  nb1=nb1k(ik)
  mb=nveck(ik) 
  ik1=ikt(ik)

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k1 = k0 + stx_ + nb1 - 1
  k2 = k0 + st_ 
!
!  Distribute the psi_d of this thread in its FFT mesh
!
  DO i=1,np
     iv=nl_d(igk_k_d(i,ik1))
     psicr(1,iv,1,k2)= REAL(psi_d(i,k1))
     psicr(2,iv,1,k2)= AIMAG(psi_d(i,k1))
  ENDDO

  IF (npol==2) THEN
     DO i=1,np
        iv=nl_d(igk_k_d(i,ik1))
        psicr(1,iv,2,k2)= REAL(psi_d(lda+i,k1))
        psicr(2,iv,2,k2)= AIMAG(psi_d(lda+i,k1))
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE put_psi_on_grid
!
!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE add_grid_to_hpsi( lda, outk, npw, nveck, &
                   nb1k, st, stx, ikt, npol, hpsi_d, psicr, nvec,      &
                   nvecx, nnr, nset) 
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
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nvecx, nvec
  !! input: the number of sets
  !! input: the number of points of the smooth mesh
  !! input: the number of vectors for each set in hpsi_d
  !! input: the number of vectors for each set in psicr
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: if .TRUE. this set is not calculated
  INTEGER, INTENT(IN), DEVICE :: ikt(nset)
  !! input: index of k in the global list of k points for this set
  INTEGER, INTENT(IN), DEVICE :: npw(nset)
  !! input: the number of plane waves. 
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of vectors to compute (must be <= nvec). 
  INTEGER, INTENT(IN), DEVICE :: nb1k(nset)
  !! input: where to start computing for each set (referred to the start of
  !!                                               the set) 
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: start of each set in psicr
  INTEGER, INTENT(IN), DEVICE :: stx(nset)
  !! input: start of each set in hpsi_d
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: number of components of each hpsi_d
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: hpsi_d(lda*npol, nvecx * nset)
  !! inp/out: the vectors to which the grid is added
  REAL(DP), DEVICE, INTENT(IN) :: psicr(2, nnr, npol, nvec * nset)
  !! input: the grid to add
  !
  !  ... local variables
  !
  INTEGER :: k0, k1, k2, ik, ik2, i, np, mb, nb1,  &
             stx_, st_, iv
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  np = npw(ik)
  stx_=stx(ik)
  st_=st(ik)
  nb1=nb1k(ik)
  mb=nveck(ik) 
  ik2=ikt(ik)

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k1 = k0 + stx_ + nb1 - 1
  k2 = k0 + st_ 
!
!  Collect hpsi_d of this thread in its FFT mesh
!
  DO i=1,np
     iv=nl_d(igk_k_d(i,ik2))
     hpsi_d(i,k1)=hpsi_d(i,k1)+CMPLX(psicr(1,iv,1,k2), &
                                     psicr(2,iv,1,k2), KIND=DP) 
  ENDDO

  IF (npol==2) THEN
     DO i=1,np
        iv=nl_d(igk_k_d(i,ik2))
        hpsi_d(lda+i,k1)=hpsi_d(lda+i,k1)+CMPLX(psicr(1,iv,2,k2), &
                                 psicr(2,iv,2,k2), KIND=DP) 
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE add_grid_to_hpsi

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft1inv_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim) 
!----------------------------------------------------------------------- 
  !
  !  This routines makes a fft along the z direction for all the 
  !  points in the plane. It goes from reciprocal to real space.
  !  It makes the calculation on GPU parallelized on k points, on 
  !  the vectors to transform and on the x direction index.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : work=>work_d, lenwrk=>lenwrk_d, wsave=> wsave_d, &
                             lensav=>lensav_d, isindex_d
  !
  IMPLICIT NONE
#include<cfft1b_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                adim
  !! input: the number of sets 
  !! input: the fft mesh and max dimensions
  !! input: the number of vectors for each k
  !! input: the maximum between nr1, nr2, nr3
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: allow to skip some k
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of psi to tranform for each k
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: the starting point of each k in the list of psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnr, npol, nvec * nset)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, k3, ik, l, m, mb, st_, sh, nnrx, ier
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  st_=st(ik)
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  nnrx=nr1x*nr2x
  l=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (l>nr1) RETURN
  k3=(k2-1)*adim+l
!
!  Now perform the inverse fourier transform, direction z (sign +1)
!
  DO m=1,nr2
     sh=l + (m-1)*nr1x
     IF (isindex_d(sh)>0) &
        CALL cfft1b ( nr3, nnrx, psicr(1,sh,1,k2), nr3*nnrx-sh+1, &
                    wsave(1,3), lensav, work(1,k3), lenwrk, ier )
  ENDDO

  IF (npol==2) THEN
     DO m=1,nr2
        sh=l + (m-1)*nr1x
        IF (isindex_d(sh)>0) &
           CALL cfft1b ( nr3, nnrx, psicr(1,sh,2,k2), nr3*nnrx-sh+1, &
                    wsave(1,3), lensav, work(1,k3), lenwrk, ier )
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE fft1inv_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft2inv_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim) 
!----------------------------------------------------------------------- 
!  
  !  This routines makes a fft along the y direction for all the 
  !  points z and x. It goes from reciprocal to real space.
  !  It makes the calculation on GPU parallelized on k points, on vectors
  !  to transform and on the x direction index.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : work=>work_d, lenwrk=>lenwrk_d, wsave=>wsave_d, &
                             lensav=>lensav_d, iplane_d
  !
  IMPLICIT NONE
#include<cfft1b_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, nr2x, nvec, &
                                adim
  !! input: the number of sets
  !! input: the fft mesh and max dimensions
  !! input: the number of vectors for each k
  !! input: the maximum between nr1, nr2, nr3
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: allow to skip some k
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of psi to tranform for each k
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: the starting point of each k in the list of psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, k3, ik, k, l, mb, st_, sh, nnrx, ier
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  st_=st(ik)
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  nnrx=nr1x*nr2x
  l=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (l>nr1) RETURN
  IF (iplane_d(l)==0) RETURN
  k3=(k2-1)*adim+l
!
!  Now perform the Fourier transform, direction y sign +1
!
  DO k=1,nr3
     sh = l + (k-1) * nnrx
     CALL cfft1b ( nr2, nr1x, psicr(1,sh,1,k2), nnrx-l+1, &
                    wsave(1,2), lensav, work(1,k3), lenwrk, ier )
  ENDDO

  IF (npol==2) THEN
     DO k=1,nr3
        sh = l + (k-1) * nnrx
        CALL cfft1b ( nr2, nr1x, psicr(1,sh,2,k2), nnrx-l+1, &
                    wsave(1,2), lensav, work(1,k3), lenwrk, ier )
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE fft2inv_dev
!
!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE fft3inv_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim) 
!----------------------------------------------------------------------- 
!  
  !  This routines makes a fft along the z direction for all the 
  !  points x and y. It goes from reciprocal to real space.
  !  It makes the calculation on GPU parallelized on k points and
  !  vectors to transform.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : work=>work_d, lenwrk=>lenwrk_d, wsave=>wsave_d, &
                             lensav=>lensav_d
  !
  IMPLICIT NONE
#include<cfft1b_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, &
                                nr2x, nvec, adim
  !! input: the number of sets
  !! input: the fft mesh and max dimensions
  !! input: the number of vectors for each k
  !! input: the maximum between nr1, nr2, nr3
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: allow to skip some k
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of psi to tranform for each k
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: the starting point of each k in the list of psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, ik, j, k, mb, st_, sh, sh1, nnrx, ier
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  st_=st(ik)
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  nnrx=nr1x*nr2x
!
!  Now perform the fourier transform, direction x sign +1
!
  DO j=1,nr2
     sh1=(j-1) * nr1x 
     DO k=1,nr3
        sh=sh1 + (k-1) * nnrx
        CALL cfft1b ( nr1, 1, psicr(1,sh+1,1,k2), nr1x, &
                    wsave(1,1), lensav, work(1,k2), lenwrk, ier )
     ENDDO
  ENDDO

  IF (npol==2) THEN
     DO j=1,nr2
        sh1=(j-1) * nr1x 
        DO k=1,nr3
           sh=sh1 + (k-1) * nnrx
           CALL cfft1b ( nr1, 1, psicr(1,sh+1,2,k2), nr1x, &
                    wsave(1,1), lensav, work(1,k2), lenwrk, ier )
        ENDDO
     ENDDO
  ENDIF
  RETURN
  END SUBROUTINE fft3inv_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft1fwd_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim) 
!----------------------------------------------------------------------- 
!  
  !  This routines makes a fft along the z direction for all the 
  !  points x and y. It goes from real to reciprocal space.
  !  It makes the calculation on GPU parallelized on k points, 
  !  vectors to transform and on the x direction index.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : work=>work_d, lenwrk=>lenwrk_d, wsave=>wsave_d, &
                             lensav=>lensav_d, isindex_d
  !
  IMPLICIT NONE
#include<cfft1f_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, &
                                nr2x, nvec, adim
  !! input: the number of sets
  !! input: the fft mesh and max dimensions
  !! input: the number of vectors for each k
  !! input: the maximum between nr1, nr2, nr3
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: allow to skip some k
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of psi to tranform for each k
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: the starting point of each k in the list of psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components
  REAL(DP), INTENT(INOUT), DEVICE :: psicr(2, nnr, npol, nvec * nset)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, k3, ik, k, l, m, mb, st_, sh, nnrx, ier
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  st_=st(ik)
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  nnrx=nr1x*nr2x
  l=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (l>nr1) RETURN
  k3=(k2-1)*adim+l
!
!  Now perform the inverse Fourier transform, direction z (sign +1)
!
  DO m=1,nr2
     sh=l + (m-1)*nr1x
     IF (isindex_d(sh)>0) &
        CALL cfft1f ( nr3, nnrx, psicr(1,sh,1,k2), nr3*nnrx-sh+1, &
                    wsave(1,3), lensav, work(1,k3), lenwrk, ier )
  ENDDO

  IF (npol==2) THEN
     DO m=1,nr2
        sh=l + (m-1)*nr1x
        IF (isindex_d(sh)>0) &
           CALL cfft1f ( nr3, nnrx, psicr(1,sh,2,k2), nr3*nnrx-sh+1, &
                    wsave(1,3), lensav, work(1,k3), lenwrk, ier )
     ENDDO
  ENDIF
  RETURN
  END SUBROUTINE fft1fwd_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE fft2fwd_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim) 
!----------------------------------------------------------------------- 
!  
  !  This routines makes a fft along the y direction for all the
  !  points x and z. It goes from real to reciprocal space.
  !  It makes the calculation on GPU parallelized on k points, 
  !  on vectors to transform and on the x direction index.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : work=>work_d, lenwrk=>lenwrk_d, wsave=>wsave_d,  &
                             lensav=>lensav_d, iplane_d
  !
  IMPLICIT NONE
#include<cfft1f_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, &
                                nr2x, nvec, adim
  !! input: the number of sets
  !! input: the fft mesh and max dimensions
  !! input: the number of vectors for each k
  !! input: the maximum between nr1, nr2, nr3
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: allow to skip some k
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of psi to tranform for each k
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: the starting point of each k in the list of psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnr, npol, nvec * nset)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, k3, ik, j, k, l, mb, st_, sh, nnrx, ier
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  st_=st(ik)
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  nnrx=nr1x*nr2x
  l=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (l>nr1) RETURN
  IF (iplane_d(l)==0) RETURN
  k3=(k2-1)*adim+l
!
!  Now perform the Fourier transform, direction y sign +1
!
  DO k=1,nr3
     sh = l + (k-1) * nnrx
     CALL cfft1f ( nr2, nr1x, psicr(1,sh,1,k2), nnrx-l+1, &
                    wsave(1,2), lensav, work(1,k3), lenwrk, ier )
  ENDDO

  IF (npol==2) THEN
     DO k=1,nr3
        sh = l + (k-1) * nnrx
        CALL cfft1f ( nr2, nr1x, psicr(1,sh,2,k2), nnrx-l+1, &
                    wsave(1,2), lensav, work(1,k3), lenwrk, ier )
     ENDDO
  ENDIF 

  RETURN
  END SUBROUTINE fft2fwd_dev
!
!----------------------------------------------------------------------- 
ATTRIBUTES(GLOBAL) SUBROUTINE fft3fwd_dev(outk, nveck, st, npol, psicr, &
                   nvec, nr1, nr2, nr3, nr1x, nr2x, nnr, nset, adim) 
!----------------------------------------------------------------------- 
!  
  !  This routines makes a fft along the x direction for all the
  !  points y and z. It goes from real to reciprocal space.
  !  It makes the calculation on GPU parallelized on k points and
  !  vectors to transform.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : work=>work_d, lenwrk=>lenwrk_d, wsave=>wsave_d, &
                             lensav=>lensav_d
  !
  IMPLICIT NONE
#include<cfft1f_interf.f90>
  !
  INTEGER, INTENT(IN), VALUE :: nset, nnr, nr1, nr2, nr3, nr1x, &
                                nr2x, nvec, adim
  !! input: the number of sets
  !! input: the fft mesh and max dimensions
  !! input: the number of vectors for each k
  !! input: the maximum between nr1, nr2, nr3
  LOGICAL, INTENT(IN), DEVICE :: outk(nset)
  !! input: allow to skip some k
  INTEGER, INTENT(IN), DEVICE :: nveck(nset)
  !! input: the number of psi to tranform for each k
  INTEGER, INTENT(IN), DEVICE :: st(nset)
  !! input: the starting point of each k in the list of psicr
  INTEGER, INTENT(IN), VALUE :: npol
  !! input: the number of components
  REAL(DP), DEVICE, INTENT(INOUT) :: psicr(2, nnr, npol, nvec * nset)
  !! inp/out: the psi vector
  !
  !  ... local variables
  !
  INTEGER :: k0, k2, ik, j, k, mb, st_, sh, sh1, nnrx, ier
  !
  ik=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik>nset) RETURN
  IF (outk(ik)) RETURN
  st_=st(ik)
  mb=nveck(ik) 

  k0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (k0>mb) RETURN
  k2 = k0 + st_ 
  nnrx=nr1x*nr2x
!
!  Now perform the Fourier transform, direction x sign +1
!
  DO j=1,nr2
     sh1=(j-1) * nr1x
     DO k=1,nr3
        sh=sh1 + (k-1) * nnrx
        CALL cfft1f ( nr1, 1, psicr(1,sh+1,1,k2), nr1x, &
                   wsave(1,1), lensav, work(1,k2), lenwrk, ier )
     ENDDO
  ENDDO

  IF (npol==2) THEN
     DO j=1,nr2
        sh1=(j-1) * nr1x
        DO k=1,nr3
           sh=sh1 + (k-1) * nnrx
           CALL cfft1f ( nr1, 1, psicr(1,sh+1,2,k2), nr1x, &
                   wsave(1,1), lensav, work(1,k2), lenwrk, ier )
        ENDDO
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE fft3fwd_dev

#endif
