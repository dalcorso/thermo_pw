!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! Copyright (C) 2023-2024 Andrea Dal Corso (for many k, many perturbations 
!                                      generalization)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE cgsolve_all_many_k (ch_psi, cg_psi, e, d0psi, dpsi, h_diag,  &
     ndmx, ngk, ethr, kter, conv_root, anorm, nbnd, npol, nk, npe, nsolv, &
     ikks, ikqs, nbnd_occ, nks, nksq)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear systems (i=1,nbnd):
  !
  !                 ( H - e_i + Q ) * dpsi_i = d0psi_i                      (1)
  !
  !     where H is a complex hermitean matrix, e_i is a real scalar, Q is a
  !     projector on occupied states, dpsi_i and d0psi_ are complex vectors
  !
  !     on input:
  !                 ch_psi   EXTERNAL  name of a subroutine:
  !                          Calculates  (H-e+Q)*psi products.
  !                          Vectors psi and psip should be dimensioned
  !                          (ndmx,nbnd)
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 e        real     unperturbed eigenvalue.
  !
  !                 dpsi     contains an estimate of the solution
  !                          vector.
  !
  !                 d0psi    contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of dpsi, ecc.
  !
  !                 ndim     integer actual row dimension of dpsi
  !
  !                 ethr     real     convergence threshold. solution
  !                          improvement is stopped when the error in
  !                          eq (1), defined as l.h.s. - r.h.s., becomes
  !                          less than ethr in norm.
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   revised (extensively)       6 Apr 1997 by A. Dal Corso & F. Mauri
  !   revised (to reduce memory) 29 May 2004 by S. de Gironcoli
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,          ONLY : DP
  USE mp_bands,       ONLY : intra_bgrp_comm, inter_bgrp_comm, use_bgrp_in_hpsi
  USE mp,             ONLY : mp_sum, mp_barrier
  USE control_flags,  ONLY : gamma_only
  USE gvect,          ONLY : gstart
  USE uspp,           ONLY : okvan
  USE lsda_mod,       ONLY : lsda, current_spin, isk
  USE scf,            ONLY : vrs
  USE many_k_ph_mod,  ONLY : current_ikb_ph, startkb_ph
  USE uspp,           ONLY : okvan, deeq_nc
  USE fft_base,       ONLY : dffts
  USE uspp,           ONLY : nkb
  USE lr_nc_mag,      ONLY : deeq_nc_save
  USE control_lr,     ONLY : alpha_pv
  USE qpoint_aux,     ONLY : ikmks
  USE io_global,      ONLY : stdout
#if defined(__CUDA)
 USE cublas
#endif

  implicit none

#if defined(__CUDA)
#include<cgsolve_all_interf.f90>
#endif
  !
  !   first the I/O variables
  !
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter(nk*npe*nsolv), & ! output: counter on iterations
             nks,  & ! input: total number of k point of this pool
             nksq,  & ! input: total number of k point (without k+q)
             nbnd, & ! input: the maximum number of bands among k points
             npol, & ! input: number of components of the wavefunctions
             nk,   & ! input: the number of k points
             nsolv, & ! input: the number of linear systems to solve
             ikks(nksq), & ! input correspondence ik -> ikk
             ikqs(nksq), & ! input correspondence ik -> ikq
             ngk(nks),  & ! number of plane waves for each k
             nbnd_occ(nks), & ! number of occupated bands
             npe     ! input: the number of perturbations

  real(DP) :: &
             e(nbnd,nks), & ! input: the actual eigenvalue
             anorm(nk*npe*nsolv),  & ! output: the norm of the error 
                                     !         in the solution
             h_diag(ndmx*npol,nbnd*nk*npe*nsolv), & ! input: an estimate of 
                                              !         ( H - \epsilon )
             ethr       ! input: the required precision

  complex(DP) :: &
             dpsi (ndmx*npol, nbnd*nk*npe*nsolv), & ! output: the solution 
                                              !  of the linear system
             d0psi (ndmx*npol, nbnd*nk*npe*nsolv) ! input: the known term

  logical :: conv_root(nk*npe*nsolv) ! output: if true the root is converged
  external ch_psi      ! input: the routine computing ch_psi
  external cg_psi      ! input: the routine computing cg_psi

#if defined(__CUDA)
ATTRIBUTES(DEVICE) :: h_diag, d0psi, dpsi
#endif
  !
  !  here the local variables
  !
  integer, parameter :: maxter = 400
  ! the maximum number of iterations
  integer :: iter, ibnd, ibnd_, ik, ik1, ipert, st_, id, nnr, isolv
  ! counters on iteration, bands
  integer , allocatable :: conv (:)
  ! if 1 the root is converged

  complex(DP), allocatable :: g (:,:), t (:,:), h (:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  complex(DP), allocatable :: hpsi (:,:), spsi (:,:), ps(:,:)
  ! some space to save hpsi and spsi
  !
  complex(DP) ::  dcgamma, dclambda
  COMPLEX(DP), ALLOCATABLE :: dclambdak(:)
  COMPLEX(DP), ALLOCATABLE :: dcgammak(:)
  !  the ratio between rho
  !  step length
  !  the scalar product
  real(DP), allocatable :: rho (:), rhoold (:), eu (:), a(:), c(:)
  real(DP), ALLOCATABLE :: rho_h(:), edev(:,:)
  ! the residue
  ! auxiliary for h_diag
  real(DP), ALLOCATABLE :: kter_eff(:)
  ! account the number of iterations with b
  ! coefficient of quadratic form
  ! bgrp parallelization auxiliary variables
  INTEGER :: n_start, n_end, my_nbnd, ierr
  LOGICAL, ALLOCATABLE :: outk(:)
  ! if .TRUE. this k point and perturbations have converged
  INTEGER, ALLOCATABLE :: lbndk(:)
  !
  INTEGER, ALLOCATABLE :: st(:)
  ! starting point of each set of bands
  INTEGER, ALLOCATABLE :: lbnd(:)

  INTEGER, ALLOCATABLE :: kdimk(:)
  INTEGER, ALLOCATABLE :: npwk(:)
  INTEGER, ALLOCATABLE :: nbndk(:)
  INTEGER, ALLOCATABLE :: nb1k(:)
  INTEGER, ALLOCATABLE :: ikt(:)
  INTEGER, ALLOCATABLE :: ikblk(:)
#if defined(__CUDA)
  INTEGER, ALLOCATABLE, DEVICE :: kdimk_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: npwk_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: ikt_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: ikblk_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: nbndk_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: nb1k_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: conv_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: st_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: lbndk_d(:)
  LOGICAL, ALLOCATABLE, DEVICE :: outk_d(:)
  REAL(DP), ALLOCATABLE, DEVICE :: psicmr(:,:,:,:)
  ATTRIBUTES(DEVICE) :: hpsi, spsi, ps
  ATTRIBUTES(DEVICE) :: hold, h, g, t, eu, rho, edev
  COMPLEX(DP), ALLOCATABLE, DEVICE :: dclambdak_d(:)
  COMPLEX(DP), ALLOCATABLE, DEVICE :: dcgammak_d(:)
#endif

  INTEGER :: done_ik, ikk, ikq, i, j, ikmk, nbd
  logical :: lsave_use_bgrp_in_hpsi, minus_b
  !auxiliary for ddot
  !
  call start_clock ('cgsolve')

  call divide (inter_bgrp_comm,nbnd,n_start,n_end)
  my_nbnd = n_end - n_start + 1

  nnr=dffts%nnr
  minus_b=(nsolv==2)

  ALLOCATE(edev(nbnd,nks))
  edev=e
  ! allocate workspace (bgrp distributed)
  allocate ( conv(nbnd*nk*npe*nsolv) )
  allocate ( hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv))
  allocate ( spsi(ndmx*npol,my_nbnd*nk*npe*nsolv))
  allocate ( ps(nbnd,my_nbnd*nk*npe*nsolv))
  allocate ( g(ndmx*npol,my_nbnd*nk*npe*nsolv), &
             t(ndmx*npol,my_nbnd*nk*npe*nsolv), &
             h(ndmx*npol,my_nbnd*nk*npe*nsolv), &
             hold(ndmx*npol,my_nbnd*nk*npe*nsolv) )
  allocate ( a(my_nbnd*nk*npe*nsolv), c(my_nbnd*nk*npe*nsolv) )
  allocate ( rho(my_nbnd*nk*npe*nsolv), rhoold(my_nbnd*nk*npe*nsolv) )
  allocate ( rho_h(my_nbnd*nk*npe*nsolv) )
  allocate ( eu(my_nbnd*nk*npe*nsolv) )
  allocate ( dclambdak(my_nbnd*nk*npe*nsolv) )
  allocate ( dcgammak(my_nbnd*nk*npe*nsolv) )
  ALLOCATE (lbndk(my_nbnd*nk*npe*nsolv) )
  ALLOCATE (st(nk*npe*nsolv))
  ALLOCATE (lbnd(nk*npe*nsolv))
  ALLOCATE (kter_eff(nk*npe*nsolv))
  ALLOCATE (outk(nk*npe*nsolv))

  ALLOCATE(kdimk(nk*npe*nsolv))
  ALLOCATE(npwk(nk*npe*nsolv))
  ALLOCATE(nb1k(nk*npe*nsolv))
  ALLOCATE(nbndk(nk*npe*nsolv))
  ALLOCATE(ikt(nk*npe*nsolv))
  ALLOCATE(ikblk(nk*npe*nsolv))

#if defined(__CUDA)
  ALLOCATE(kdimk_d(nk*npe*nsolv))
  ALLOCATE(npwk_d(nk*npe*nsolv))
  ALLOCATE(nb1k_d(nk*npe*nsolv))
  ALLOCATE(nbndk_d(nk*npe*nsolv))
  ALLOCATE(ikt_d(nk*npe*nsolv))
  ALLOCATE(ikblk_d(nk*npe*nsolv))
  ALLOCATE(st_d(nk*npe*nsolv))
  ALLOCATE(outk_d(nk*npe*nsolv))
  ALLOCATE(psicmr(2,nnr,npol,my_nbnd*nk*npe*nsolv))
  ALLOCATE(lbndk_d(my_nbnd*nk*npe*nsolv))
  ALLOCATE(dclambdak_d(my_nbnd*nk*npe*nsolv))
  ALLOCATE(dcgammak_d(my_nbnd*nk*npe*nsolv))
  ALLOCATE(conv_d(nbnd*nk*npe*nsolv))
#endif
  WRITE(stdout,'(5x,"Memory allocated in cgsolve_all")')
  CALL print_gpu_memory()

  kter_eff(1:nk*npe*nsolv) = 0.d0 ; conv (1:nbnd*nk*npe*nsolv) = 0

  ! bgrp parallelization is done outside h_psi/s_psi. set use_bgrp_in_hpsi temporarily to false
  lsave_use_bgrp_in_hpsi = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
  !$acc enter data create(a(1:my_nbnd*nk*npe*nsolv),c(1:my_nbnd*nk*npe*nsolv))
  !copyin(e(1:nbnd,1:nks))
  !$acc kernels 
  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)
  !$acc end kernels 
  DO ik1=1,nk
     DO isolv=1,nsolv
        DO ipert=1,npe
           id=ik1 + (ipert-1)*nk + (isolv-1)*nk*npe
           ik=ik1+startkb_ph(current_ikb_ph)
           ikk=ikks(ik)
           ikq=ikqs(ik)
           IF (npol==2) THEN
              kdimk(id)=ndmx*npol
           ELSE
              kdimk(id)=ngk(ikq)
           ENDIF
           npwk(id)=ngk(ikq)
           ikt(id)=ikq
           ikblk(id)=ik1
           nb1k(id)=1
           nbndk(id)=nbnd_occ(ikk)
           st(id)=my_nbnd * (id-1) 
           outk(id)=.FALSE.
        ENDDO
     ENDDO
  ENDDO
#if defined(__CUDA)  
  kdimk_d=kdimk
  npwk_d=npwk
  nb1k_d=nb1k
  nbndk_d=nbndk
  ikt_d=ikt
  ikblk_d=ikblk
  st_d=st
  outk_d=outk
#endif
  DO ik1=1,nk
     ik=ik1+startkb_ph(current_ikb_ph)
     ikk=ikks(ik)
     nbd=nbnd_occ(ikk)
     DO isolv=1, nsolv
        IF (isolv==2) ikk=ikmks(ik)
        DO ipert=1,npe
           id=ik1+(ipert-1)*nk + (isolv-1) * nk * npe
           st_=st(id)
           !$acc parallel loop 
           DO ibnd = 1, nbd
              eu(st_+ibnd) = edev (ibnd,ikk)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

iterate:  do iter = 1, maxter
     !
     !    compute the gradient. can reuse information from previous step
     !
     IF (mod(iter,10)==0) WRITE(6,*) 'cgsolve_all_mk iteration:', iter
     if (iter == 1) then
#if defined(__CUDA)
        !$acc kernels 
         hpsi = (0.d0, 0.d0)
         spsi = (0.d0, 0.d0)
        !$acc end kernels
        IF (okvan) THEN
           CALL h_s_psik_dev(ndmx, outk_d, kdimk_d, npwk_d, nbndk_d, &
                    nb1k_d, st_d, st_d, ikt_d, ikblk_d, npol, dpsi, &
                    hpsi, spsi, psicmr, my_nbnd, my_nbnd, nnr,      &
                    nk*npe*nsolv, minus_b)
           ierr=cudaDeviceSynchronize()
        ELSE
           CALL h_psik_dev(ndmx, outk_d, kdimk_d, npwk_d, nbndk_d, nb1k_d, &
                   st_d, st_d, ikt_d, ikblk_d, npol, dpsi, hpsi, psicmr,  &
                   my_nbnd, my_nbnd, nnr, nk*npe*nsolv, minus_b)
           ierr=cudaDeviceSynchronize() 
           !$acc parallel loop 
           DO i=1,ndmx*npol
              spsi(i,1:my_nbnd*nk*npe*nsolv)=dpsi(i,1:my_nbnd*nk*npe*nsolv)
           ENDDO
        ENDIF
 
        call ch_psi_dev (ndmx, outk_d, kdimk_d, npwk_d, st_d, nbndk_d, &
              dpsi, g, hpsi, spsi, ps, eu, current_ikb_ph,             &
              npol, nk, npe, nsolv, nkb, nbnd, my_nbnd, alpha_pv)
#else
        DO ik1=1,nk
           ik=ik1+startkb_ph(current_ikb_ph)
           ikk=ikks(ik)
           ikq=ikqs(ik)
           ndim=ngk(ikq)
           IF (lsda) current_spin = isk (ikk)
           DO isolv=1, nsolv
              ikmk=ikk
              IF (isolv==2) THEN
                 ikmk=ikmks(ik)
                 vrs(:,2:4)=-vrs(:,2:4)
                 IF (okvan) deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,2)
              ENDIF
              DO ipert=1,npe
                 id=ik1 + (ipert-1)*nk + (isolv-1)*nk*npe
                 IF (outk(id)) CYCLE
                 st_=st(id)
                 CALL ch_psi (ndim, dpsi(1,st_+n_start), g(1,st_+1), &
                   hpsi(1,st_+1), spsi(1,st_+1), ps(1,st_+1),        &
                   e(n_start,ikmk), ik, nbnd_occ(ikk), id, isolv)
              ENDDO
              IF (isolv==2) THEN
                 vrs(:,2:4)=-vrs(:,2:4)
                 IF (okvan) deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,1)
              ENDIF
           ENDDO
        ENDDO
#endif        
        DO ik1=1,nk
           ik=ik1+startkb_ph(current_ikb_ph)
           ikk=ikks(ik)
           ikq=ikqs(ik)
           ndim=ngk(ikq)
           DO isolv=1, nsolv
              DO ipert=1,npe
                 id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
                 IF (outk(id)) CYCLE
                 st_=st(id)
                 do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                    IF (ibnd>nbnd_occ(ikk)) CYCLE
                    call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,st_+ibnd), 1, &
                                   g(1,st_+ibnd_), 1)
                 enddo
                 IF (npol==2) THEN
                    do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                       IF (ibnd>nbnd_occ(ikk)) CYCLE
                       call zaxpy (ndim, (-1.d0,0.d0), d0psi(ndmx+1,st_+ibnd),&
                                1, g(ndmx+1,st_+ibnd_), 1)
                    enddo
                 END IF
              ENDDO
           ENDDO
        ENDDO
     endif
     !
     !    compute preconditioned residual vector and convergence check
     !
#if defined(__CUDA)
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)

              lbnd(id) = 0
              do ibnd = n_start, n_end ;  ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    lbnd(id) = lbnd(id)+1
                    lbndk(st_+ibnd)=lbnd(id)
                 endif
              enddo
           ENDDO
        ENDDO
     ENDDO
     conv_d=conv
     lbndk_d=lbndk
     CALL cgsolve_all_loop2<<<dim3(nk,npe*nsolv,nbnd/4+1),dim3(1,1,4)>>>  &
             (ndmx, outk_d, st_d, conv_d, lbndk_d, h, g, h_diag, rho, &
              current_ikb_ph, npol, nk, npe, nsolv, nbnd, my_nbnd)
     ierr=cudaDeviceSynchronize() 
     rho_h=rho 
#else
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)

              lbnd(id) = 0
              do ibnd = n_start, n_end ;  ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    lbnd(id) = lbnd(id)+1
                    call zcopy (ndmx*npol, g (1, st_+ibnd_), 1, &
                                                      h (1, st_+ibnd_), 1)
                    call cg_psi(ndmx, ndim, 1, h(1,st_+ibnd_), &
                                                       h_diag(1,st_+ibnd) )
           
                    CALL MYDDOTV3 (2*ndmx*npol, h(1,st_+ibnd_), 1, &
                                      g(1,st_+ibnd_), 1, rho(st_+lbnd(id)))
                 endif
              enddo
           ENDDO
        ENDDO
     ENDDO
     rho_h=rho
#endif
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk + (isolv-1) * nk * npe
              IF (outk(id)) CYCLE
              st_=st(id)
              kter_eff(id) = kter_eff(id)+DBLE(lbnd(id))/DBLE (nbnd_occ(ikk))
              call mp_sum( rho_h(st_+1:st_+lbnd(id)), intra_bgrp_comm )
              do ibnd = n_end, n_start, -1 ; ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv(st_+ibnd).eq.0) then
                    rho_h(st_+ibnd_)=rho_h(st_+lbnd(id))
                    lbnd(id) = lbnd(id) -1
                    anorm(id) = sqrt (rho_h (st_+ibnd_) )
!                    WRITE(6,*) 'anorm', iter, id, anorm(id)
                    if (anorm(id).lt.ethr) conv (st_+ibnd) = 1
                 endif
              enddo
!
              conv_root(id) = .true.
              do ibnd = n_start, n_end
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 conv_root(id) = conv_root(id).and.(conv (st_+ibnd) .eq.1)
              enddo
              if (conv_root(id)) outk(id)=.TRUE.
           ENDDO
        ENDDO
     ENDDO
#if defined(__CUDA)
     outk_d=outk

     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk+(isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)
              !
              !   compute the step direction h. Conjugate it to previous step
              !
              lbnd(id) = 0
              do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    lbnd(id) = lbnd(id)+1
                    lbndk(st_+ibnd)=lbnd(id)
                    if (iter.ne.1) then
                       dcgammak(st_+ibnd) = rho_h (st_+ibnd_) / rhoold (st_+ibnd_)
                    endif
                 endif
              enddo
           ENDDO
        ENDDO
     ENDDO
     conv_d=conv
     lbndk_d=lbndk
     dcgammak_d=dcgammak
     CALL cgsolve_all_loop4<<<dim3(nk,npe*nsolv,nbnd),dim3(1,1,1)>>>(ndmx, &
         outk_d, st_d, conv_d, lbndk_d, h, hold, dcgammak_d, eu, edev,  &
         current_ikb_ph, npol, nk, npe, nsolv, nbnd, my_nbnd, iter, nks)
     ierr=cudaDeviceSynchronize()
#else
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           ikmk=ikk
           IF (isolv==2) ikmk=ikmks(ik)
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk+(isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)
              !
              !  compute the step direction h. Conjugate it to previous step
              !
              lbnd(id) = 0
              do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    call dscal (2 * ndmx * npol, - 1.d0, h (1, st_+ibnd_), 1)
                    if (iter.ne.1) then
                       dcgamma = rho_h (st_+ibnd_) / rhoold (st_+ibnd_)
                       call zaxpy (ndmx*npol, dcgamma, hold (1, st_+ibnd_), &
                                                1, h (1, st_+ibnd_), 1)
                    endif
!
! here hold is used as auxiliary vector in order to efficiently compute t = A*h
! it is later set to the current (becoming old) value of h
!
                    lbnd(id) = lbnd(id)+1
                    call zcopy (ndmx*npol, h (1, st_+ibnd_), 1, &
                                                  hold (1, st_+lbnd(id)), 1)
                    !$acc serial 
                    eu (st_+lbnd(id)) = e (ibnd,ikmk)
                    !$acc end serial
                 endif
              enddo
           ENDDO
        ENDDO
     ENDDO
#endif
#if defined(__CUDA)
    nbndk_d=lbnd
    !$acc kernels 
    hpsi (:,:) = (0.d0, 0.d0)
    spsi (:,:) = (0.d0, 0.d0)
    !$acc end kernels

    IF (okvan) THEN
       CALL h_s_psik_dev(ndmx, outk_d, kdimk_d, npwk_d, nbndk_d, nb1k_d, &
          st_d, st_d, ikt_d, ikblk_d, npol, hold, hpsi, spsi, psicmr,     &
          my_nbnd, my_nbnd, nnr, nk*npe*nsolv, minus_b)
       ierr=cudaDeviceSynchronize() 
    ELSE
       CALL h_psik_dev(ndmx, outk_d, kdimk_d, npwk_d, nbndk_d, nb1k_d,  &
          st_d, st_d, ikt_d, ikblk_d, npol, hold, hpsi, psicmr, my_nbnd, &
          my_nbnd, nnr, nk*npe*nsolv, minus_b)
       ierr=cudaDeviceSynchronize()
       !$acc parallel loop 
       DO i=1,ndmx*npol
          spsi(i,1:my_nbnd*nk*npe*nsolv)=hold(i,1:my_nbnd*nk*npe*nsolv)
       ENDDO
    ENDIF

    CALL ch_psi_dev (ndmx, outk_d, kdimk_d, npwk_d, st_d, nbndk_d, &
              hold, t, hpsi, spsi, ps, eu, current_ikb_ph,      &
              npol, nk, npe, nsolv, nkb, nbnd, my_nbnd, alpha_pv)
#else
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        IF (lsda) current_spin = isk (ikk)
        DO isolv=1, nsolv
           IF (isolv==2) THEN
              vrs(:,2:4)=-vrs(:,2:4)
              IF (okvan) deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,2)
           ENDIF
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk+(isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)
              !
              !        compute t = A*h
              !
              call ch_psi (ndim, hold(1,st_+1), t(1,st_+1), &
                        hpsi(1,st_+1), spsi(1,st_+1), ps(1,st_+1), &
                        eu(st_+1), ik, lbnd(id), id, isolv)
           ENDDO
           IF (isolv==2) THEN
              vrs(:,2:4)=-vrs(:,2:4)
              IF (okvan) deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,1)
           ENDIF
        ENDDO
     ENDDO
#endif     
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk+(isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)

              lbnd(id)=0
              do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    lbnd(id)=lbnd(id)+1
                    lbndk(st_+ibnd)=lbnd(id)
                 end if
              end do
           ENDDO
        ENDDO
     ENDDO

#if defined(__CUDA)    
     conv_d=conv
     lbndk_d=lbndk
     !$acc host_data use_device(a,c)
     CALL cgsolve_all_loop5<<<dim3(nk,npe*nsolv,nbnd/4+1),dim3(1,1,4)>>>&
             (ndmx, outk_d, st_d, conv_d, lbndk_d, g, h, t, a, c, &
             current_ikb_ph, npol, nk, npe, nsolv, nbnd, my_nbnd)
     ierr=cudaDeviceSynchronize()
     !$acc end host_data
     !$acc update host(a,c)
#else
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)

              lbnd(id)=0
              do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    lbnd(id)=lbnd(id)+1
                    !$acc host_data use_device(a,c)
                    CALL MYDDOTV3(2*ndmx*npol, h(1,st_+ibnd_), 1, &
                              g(1,st_+ibnd_), 1, a(st_+lbnd(id)))
                    CALL MYDDOTV3(2*ndmx*npol, h(1,st_+ibnd_), 1, &
                              t(1,st_+lbnd(id)), 1, c(st_+lbnd(id)))
                    !$acc end host_data
                 end if
              end do
           ENDDO
        ENDDO
     ENDDO
     !$acc update host(a,c)
#endif
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)

              call mp_sum(  a(st_+1:st_+lbnd(id)), intra_bgrp_comm )
              call mp_sum(  c(st_+1:st_+lbnd(id)), intra_bgrp_comm )
              lbnd(id)=0
              do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    lbnd(id)=lbnd(id)+1
                    lbndk(st_+ibnd)=lbnd(id)
                    dclambdak(st_+ lbnd(id)) = &
                          CMPLX( - a(st_+lbnd(id)) / c(st_+lbnd(id)), &
                                                           0.d0,kind=DP)
                    rhoold (st_+ibnd_) = rho_h (st_+ibnd_)
                 endif
              enddo
           ENDDO
        ENDDO
     ENDDO
#if defined(__CUDA)
 dclambdak_d=dclambdak
 conv_d=conv
 lbndk_d=lbndk
 CALL cgsolve_all_loop6<<<dim3(nk*npe*nsolv,nbnd,(ndmx*npol)/32+1),&
             dim3(1,1,32)>>> &
            (ndmx, outk_d, st_d, conv_d, lbndk_d, dpsi, g, h, hold, t, &
            dclambdak_d, current_ikb_ph, npol, nk, npe, nsolv, nbnd, my_nbnd)
 ierr=cudaDeviceSynchronize()
#else
     DO ik1=1,nk
        ik=ik1+startkb_ph(current_ikb_ph)
        ikk=ikks(ik)
        ikq=ikqs(ik)
        ndim=ngk(ikq)
        DO isolv=1, nsolv
           DO ipert=1,npe
              id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
              IF (outk(id)) CYCLE
              st_=st(id)

              do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
                 IF (ibnd>nbnd_occ(ikk)) CYCLE
                 if (conv (st_+ibnd) .eq.0) then
                    dclambda=dclambdak(st_+lbndk(st_+ibnd))
                    !
                    !    move to new position
                    !
                    call zaxpy (ndmx*npol, dclambda, h(1,st_+ibnd_), 1, &
                                                  dpsi(1,st_+ibnd), 1)
                    !!dpsi()=dpsi+dclambda*h
                    !
                    !    update to get the gradient
                    !
                    !g=g+lam
                    call zaxpy (ndmx*npol, dclambda, t(1,st_+lbndk(st_+ibnd)),&
                                            1, g(1,st_+ibnd_), 1)
                    !!g=g+dclambda*t
                    !
                    !    save current (now old) h and rho for later use
                    !
                    !!hold=h
                    call zcopy (ndmx*npol, h(1,st_+ibnd_), 1, &
                                                     hold(1,st_+ibnd_), 1)
                 endif
              enddo
           ENDDO
        ENDDO
     ENDDO
#endif
     done_ik = COUNT(outk(1:nk*npe*nsolv))
     IF (done_ik==nk*npe*nsolv) EXIT iterate
  enddo iterate
  !$acc exit data delete(a,c,e) 

  ! deallocate workspace not needed anymore
  deallocate (eu) ; deallocate (rho, rho_h, rhoold) ; deallocate (a,c) ; deallocate (g, t, h, hold) ; deallocate(hpsi,spsi,ps)

  ! wait for all bgrp to complete their task
  CALL mp_barrier( inter_bgrp_comm )

  ! check if all root converged across all bgrp
  call mp_sum( conv, inter_bgrp_comm )
  call mp_sum( kter_eff, inter_bgrp_comm )
  
  DO ik1=1,nk
     ik=ik1+startkb_ph(current_ikb_ph)
     ikk=ikks(ik)
     ikq=ikqs(ik)
     ndim=ngk(ikq)
     DO isolv=1, nsolv
        DO ipert=1,npe
           id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
           st_=st(id)
           conv_root(id) = .true.
           do ibnd = 1, nbnd
              IF (ibnd>nbnd_occ(ikk)) CYCLE
              conv_root(id) = conv_root(id).and. (conv (st_+ibnd) .eq.1)
           enddo

  ! collect the result
           if (n_start > 1 ) dpsi(:, st_+1:st_+n_start-1) = (0.d0,0.d0) 
           if (n_end < nbnd_occ(ikk)) dpsi(:, st_+n_end+1:st_+nbnd_occ(ikk)) &
                                                            = (0.d0,0.d0)
           kter(id) = kter_eff(id)
        ENDDO
     ENDDO
  ENDDO
  deallocate (conv)

  call mp_sum( dpsi, inter_bgrp_comm )

  ! restore the value of use_bgrp_in_hpsi to its saved value
  use_bgrp_in_hpsi = lsave_use_bgrp_in_hpsi

  DEALLOCATE(st)
  DEALLOCATE(lbnd)
  DEALLOCATE(kter_eff)
  DEALLOCATE(outk)

  DEALLOCATE(kdimk)
  DEALLOCATE(npwk)
  DEALLOCATE(nb1k)
  DEALLOCATE(nbndk)
  DEALLOCATE(ikt)
  DEALLOCATE(ikblk)

  deallocate ( dclambdak )
  deallocate ( dcgammak )
  deallocate ( lbndk )

#if defined(__CUDA)
  DEALLOCATE(kdimk_d)
  DEALLOCATE(npwk_d)
  DEALLOCATE(nb1k_d)
  DEALLOCATE(nbndk_d)
  DEALLOCATE(ikt_d)
  DEALLOCATE(ikblk_d)
  DEALLOCATE(st_d)
  DEALLOCATE(outk_d)
  DEALLOCATE(psicmr)
  DEALLOCATE(lbndk_d)
  DEALLOCATE(conv_d)
  DEALLOCATE(dclambdak_d)
  DEALLOCATE(dcgammak_d)
  DEALLOCATE(edev)
#endif

  call stop_clock ('cgsolve')
  return
end subroutine cgsolve_all_many_k

