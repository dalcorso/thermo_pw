!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! Copyright (C) Dic. 2022- Andrea Dal Corso (extension to many k with GPU
!                          mainly in CUDA fortran but some acc commands 
!                          (from I. Carmineo) still used to parallelize loops).
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cegterg_vk( h_psii, s_psii, uspp, g_psii, npw, npwx, nvec, &
                    nvecx, npol, evc, ethr, e, btype, notcnvk, lrot,  &
                    dav_iter, nhpsi, nk, nkb )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problems:
  !
  ! ... ( H (ik) - e (ik) S (ik) ) * evc (ik) = 0
  !
  ! ... where H (ik) are a set of hermitean operators, e(ik) is a real scalar,
  ! ... S(ik) are a set of overlap matrix, evc(ik) is a complex vector
  !
  ! see the original routine for an explanation of the input variables
  !
#if defined(__CUDA)
  use cudafor
  use cublas
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
                            nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum, mp_gather, mp_bcast, mp_size,&
                            mp_type_create_column_section, mp_type_free
  USE device_memcpy_m, ONLY : dev_memcpy, dev_memset
  USE io_global,      ONLY : stdout
  USE many_k_mod,     ONLY : current_ikb, startkb
  USE fft_base,       ONLY : dffts
  !
  IMPLICIT NONE

#if defined(__CUDA)
#include<cegterg_vk_interf.f90>
#include<g_psi_interf.f90>
#include<diago_interf.f90>
#endif
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npwx, nvec, nvecx, npol, nk, nkb
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
    ! number of Hamiltonians to solve
    ! number of nonlocal projectors
  INTEGER, INTENT(IN) :: npw(nk)
    ! number of plane waves for each k point
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol, nvec*nk)
    !  evc contains the refined estimates of the eigenvectors  
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec*nk)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec*nk)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter(nk), notcnvk(nk)
  ! integer number of iterations performed for each ik
  ! number of unconverged roots for each ik
  INTEGER, INTENT(OUT) :: nhpsi(nk)
  ! total number of indivitual hpsi for each ik
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, np, kdmx, done_ik, n, m, ipol, ik, nbn
    ! counter on iterations
    ! counter on the reduced basis vectors
    ! adapted npwx
    ! number of converged k points
    ! do-loop counters
  INTEGER, ALLOCATABLE :: nbasek(:), kdimk(:), nb1k(:)
  ! As below for each k
  INTEGER :: nbase, kdim, nb1, notcnv
    ! dimension of the reduced basis
    ! adapted npw
    ! dimension of the reduced basis plus one
    ! notcnv roots for each k
  INTEGER :: n_start, n_end, my_n
    !
    !   these variables are for band division. Probably not working anymore 
    !   here
    !
  INTEGER :: column_section_type
    ! defines a column section for communication
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc(:,:,:), sc(:,:,:), vc(:,:,:)
    ! Hamiltonians on the reduced basis
    ! S matrices on the reduced basis
    ! the eigenvectors of the Hamiltonians
  REAL(DP), ALLOCATABLE :: ew(:,:)
    ! eigenvalues of the reduced Hamiltonians
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr, time, scnds
    ! threshold for empty bands
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
    ! receive counts and memory offsets
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER, ALLOCATABLE :: numblock(:)
    ! chunking parameters
  INTEGER :: i, j, k 
    !
  INTEGER, ALLOCATABLE :: st(:), stx(:)
  ! starting point of each ik in the arrays of size nvec*nk.
  ! starting point of each ik in the arrays of size nvecx*nk.
  INTEGER, ALLOCATABLE :: ikt(:), ikblk(:)
  ! index in the total list of k of this ik
  ! index in the current block of this ik
  LOGICAL, ALLOCATABLE :: outk(:), enter(:)
  ! When outk(ik) is .TRUE. the ik point has finished its iterations
  ! When enter(ik) is .TRUE. the ik point needs to update the basis size
  REAL(DP), EXTERNAL :: MYDDOT_VECTOR_GPU
  !$acc routine(MYDDOT_VECTOR_GPU) vector
  ! auxiliary variables to contain the current value of the parameters
  INTEGER :: nnrs, st_, stx_, npw_, nkeff, lwork, lrwork, liwork
!
!  variables to diagonalize only a subset of the k points
!
  INTEGER, ALLOCATABLE :: times(:), start(:)
  INTEGER :: bs, nblock, nbs, rest, ib
  !
  !  declaration of device variables
  !
#if defined(__CUDA)

  !
  !   Several arrays of dimensions needed to pass the information
  !   to the global routines. See above for their use.
  !
  INTEGER, ALLOCATABLE, DEVICE :: nbasek_d(:), kdimk_d(:), nb1k_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: notcnvk_d(:)
  LOGICAL, ALLOCATABLE, DEVICE :: outk_d(:), enter_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: st_d(:), stx_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: npw_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: ikt_d(:), ikblk_d(:)
  LOGICAL, ALLOCATABLE, DEVICE :: conv_d(:)
  INTEGER, ALLOCATABLE, DEVICE :: dav_iter_d(:)
!
!  Auxiliary space to pass the eigenvalues (should not be needed...)
!
  REAL(DP), ALLOCATABLE, DEVICE :: e_d(:)
!
!  work space for the diago_dev is allocated on device here
!
  COMPLEX(DP), ALLOCATABLE, DEVICE :: work_d(:,:)
  REAL(DP), ALLOCATABLE, DEVICE :: rwork_d(:,:)
  REAL(DP), ALLOCATABLE, DEVICE :: hdiago_d(:,:)
  REAL(DP), ALLOCATABLE, DEVICE :: sdiago_d(:,:)
  INTEGER,  ALLOCATABLE, DEVICE :: ifail_d(:,:)
  INTEGER,  ALLOCATABLE, DEVICE :: iwork_d(:,:)
  INTEGER,  ALLOCATABLE, DEVICE :: m_d(:)
  INTEGER,  ALLOCATABLE, DEVICE :: info_d(:)
!
!  work space for the fft. We allocate a smooth grid for each k point.
!
  REAL(DP), ALLOCATABLE, DEVICE :: psicmr(:,:,:,:)

  ATTRIBUTES(DEVICE) :: psi, hpsi, spsi, vc, hc, sc, ew
#endif
  !
  EXTERNAL  h_psii,    s_psii,    g_psii
    ! h_psii(npwx,npw(ik),nvec,psi,hpsi,ik)
    !     calculates H|psi>
    ! s_psii(npwx,npw(ik),nvec,spsi,ik)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx*npol,nvecx*nk)
    ! g_psi(npwx,npw(ik),notcnv(ik),psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first stx(ik)+nvec columns contain the trial eigenvectors
  !
  CALL start_clock( 'cegterg' )!; write(*,*) 'start cegterg' ; FLUSH(6)

  !$acc data deviceptr(evc, e )
  nhpsi = 0
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdmx = npwx
     !
  ELSE
     !
     kdmx = npwx*npol
     !
  END IF
  !
#if ! defined(__CUDA)
  ! compute the number of chuncks
  ALLOCATE(numblock(nk))
  DO ik=1, nk
     numblock(ik)  = (npw(ik)+blocksize-1)/blocksize
  ENDDO
#endif
  !
  ALLOCATE(  psi( npwx*npol, nvecx * nk ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx*npol, nvecx * nk ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx*npol, nvecx * nk ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ALLOCATE( sc( nvecx, nvecx, nk ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc ', ABS(ierr) )
  ALLOCATE( hc( nvecx, nvecx, nk ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc ', ABS(ierr) )
  ALLOCATE( vc( nvecx, nvecx, nk ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc ', ABS(ierr) )
  ALLOCATE( ew( nvecx, nk ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew ', ABS(ierr) )
  ALLOCATE( conv( nvec * nk ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate conv ', ABS(ierr) )
  ALLOCATE( recv_counts(mp_size(inter_bgrp_comm))) 
  ALLOCATE( displs(mp_size(inter_bgrp_comm)) )

  nnrs=dffts%nnr

  ALLOCATE(nbasek(nk))
  ALLOCATE(kdimk(nk))
  ALLOCATE(nb1k(nk))
  ALLOCATE(outk(nk))
  ALLOCATE(enter(nk))
  ALLOCATE(st(nk))
  ALLOCATE(stx(nk))
  ALLOCATE(ikt(nk))
  ALLOCATE(ikblk(nk))

  CALL find_lwork(nvecx,lwork)
  lrwork=7 * nvecx  
  liwork=5 * nvecx  

#if defined(__CUDA)
  ALLOCATE(nbasek_d(nk))
  ALLOCATE(kdimk_d(nk))
  ALLOCATE(npw_d(nk))
  ALLOCATE(nb1k_d(nk))
  ALLOCATE(outk_d(nk))
  ALLOCATE(enter_d(nk))
  ALLOCATE(st_d(nk))
  ALLOCATE(stx_d(nk))
  ALLOCATE(ikt_d(nk))
  ALLOCATE(ikblk_d(nk))
  ALLOCATE(notcnvk_d(nk))
  ALLOCATE(conv_d(nvec*nk))

  ALLOCATE(work_d(lwork,nk))
  ALLOCATE(rwork_d(lrwork,nk))
  ALLOCATE(hdiago_d(nvecx,nk))
  ALLOCATE(sdiago_d(nvecx,nk))
  ALLOCATE(iwork_d(liwork,nk))
  ALLOCATE(ifail_d(nvecx,nk))
  ALLOCATE(m_d(nk))
  ALLOCATE(info_d(nk))

  ALLOCATE(dav_iter_d(nk))
  ALLOCATE(e_d(nvec*nk))
  ALLOCATE(psicmr(2,nnrs,npol,nvec*nk))

  bs=MAX(1,10000/nvecx)
  nblock=(nk-1)/bs+1
!  WRITE(stdout,'(5x,"GPU diagonalization with",i5," blocks")') nblock
  ALLOCATE(times(nblock))
  ALLOCATE(start(nblock))
  nbs= nk/nblock
  rest=nk-nbs*nblock
  DO ib=1,nblock
     times(ib)=nbs
     IF (ib<=rest) times(ib)=times(ib)+1
  ENDDO
  start(1)=1
  DO ib=2,nblock
     start(ib)=start(ib-1)+times(ib-1)
  ENDDO
#endif
WRITE(6,*) 'Allocated memory in cegter_vk'
CALL print_gpu_memory()
  !
  DO ik = 1, nk 
     nbasek(ik)  = nvec
     IF (npol==1) THEN
        kdimk(ik) = npw(ik)
     ELSE
        kdimk(ik) = npwx*npol
     ENDIF
     nb1k(ik)=1
     notcnvk(ik) = nvec
     outk(ik)=.FALSE.
     enter(ik)=.TRUE.
     st(ik) = (ik - 1)*nvec
     stx(ik) = (ik - 1)*nvecx
     ikt(ik) = ik + startkb(current_ikb)
     ikblk(ik) = ik 
  ENDDO
  conv   = .FALSE.
#if defined(__CUDA)
  nbasek_d=nbasek
  kdimk_d=kdimk
  npw_d=npw
  nb1k_d=nb1k
  notcnvk_d=notcnvk
  outk_d=outk
  enter_d=enter
  st_d=st
  stx_d=stx
  ikt_d=ikt
  ikblk_d=ikblk
#endif
  !
#if defined(__CUDA)  
  CALL copy_psi_gpu<<<dim3((npwx*npol)/4+1,nvec/4+1,nk/32+1),dim3(4,4,32)>>> &
       ( npwx, outk_d, enter_d, kdimk_d, stx_d, st_d, npol, psi, evc,  &
                                                       nvec, nvecx, nk )
#else
  DO ik=1, nk
     stx_=stx(ik)
     st_=st(ik)
     psi(:,stx_+1:stx_+nvec)=evc(:,st_+1:st_+nvec)
  ENDDO
#endif
  !
  ! ... hpsi contains h times the basis vectors
  !
#if defined(__CUDA)

  IF (uspp) THEN
     CALL h_s_psik_dev(npwx, outk_d, kdimk_d, npw_d, notcnvk_d, nb1k_d, &
             st_d, stx_d, ikt_d, ikblk_d, npol, psi, hpsi, spsi, psicmr, &
             nvec, nvecx, nnrs, nk, .FALSE.)
  ELSE
     CALL h_psik_dev(npwx, outk_d, kdimk_d, npw_d, notcnvk_d, nb1k_d, st_d, &
          stx_d, ikt_d, ikblk_d, npol, psi, hpsi, psicmr, nvec, nvecx, &
          nnrs, nk, .FALSE.)
  ENDIF
#endif

#if ! defined (__CUDA)
  DO ik=1, nk
     stx_=stx(ik)
     CALL h_psii( npwx, npw(ik), nvec, psi(1,stx_+1), &
                                      hpsi(1,stx_+1), ik )  
  !
  ! ... spsi contains s times the basis vectors
  !
     IF ( uspp ) CALL s_psii( npwx, npw(ik), nvec, psi(1,stx_+1), &
                                                  spsi(1,stx_+1), ik )
  ENDDO

#endif
  DO ik=1,nk
     nhpsi(ik) = nhpsi(ik) + nvec
  ENDDO

  CALL start_clock( 'cegterg:init' )
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
#if defined (__CUDA)
  IF (uspp) THEN
     CALL cegterg_init_us<<<dim3(nvec/4+1,nvecx/4+1,nk/32+1),       &
       dim3(4,4,32)>>>(outk_d, nbasek_d, kdimk_d, stx_d, hpsi,      &
                       spsi, psi, hc, sc, kdmx, nvecx, nk)
  ELSE
     CALL cegterg_init<<<dim3(nvec/4+1,nvecx/4+1,nk/32+1),          &
       dim3(4,4,32)>>>(outk_d, nbasek_d, kdimk_d, stx_d, hpsi, psi, &
                       hc, sc, kdmx, nvecx, nk)
  ENDIF
  ierr=cudaDeviceSynchronize()
#else
  DO ik = 1, nk
     nbase=nbasek(ik)
     kdim=kdimk(ik)
     stx_=stx(ik)
     CALL divide_all(inter_bgrp_comm,nbase,n_start,n_end,recv_counts,displs)
     CALL mp_type_create_column_section(sc(1,1,1), 0, nbase, nvecx, &
                                                          column_section_type)
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
   
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi(1,stx_+1), &
                kdmx, hpsi(1,stx_+n_start), kdmx, ZERO, &
                hc(1,n_start,ik), nvecx )

     IF (n_start .LE. n_end) &
        CALL mp_sum( hc(1:nbase, n_start:n_end, ik), intra_bgrp_comm )
     CALL mp_gather( hc(:,:,ik), column_section_type, recv_counts, displs, &
                  root_bgrp_id, inter_bgrp_comm )
     !
     IF ( uspp ) THEN
     !
        IF (n_start .LE. n_end) &
        CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi(1,stx_+1), &
                 kdmx, spsi(1,stx_+n_start), kdmx, &
                 ZERO, sc(1,n_start,ik), nvecx )
        !
     ELSE
        !
        IF (n_start .LE. n_end) &
        CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi(1,stx_+1), &
                 kdmx, psi(1,stx_+n_start), kdmx, &
                 ZERO, sc(1,n_start,ik), nvecx )
        !
     END IF
     IF (n_start .LE. n_end) &
         CALL mp_sum( sc(1:nbase, n_start:n_end, ik), intra_bgrp_comm )
     CALL mp_gather( sc(:,:,ik), column_section_type, recv_counts, displs, &
                                            root_bgrp_id, inter_bgrp_comm )

     CALL mp_type_free( column_section_type )
  ENDDO
#endif

#if defined (__CUDA)
  DO ik = 1, nk
     nbase=nbasek(ik)
     n_start=1
     n_end=nbase
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end

     IF (n_start .LE. n_end) & 
        CALL mp_sum( hc(:,:,ik), 1, nbase, n_start, n_end, intra_bgrp_comm )
   
     IF (n_start .LE. n_end) & 
         CALL mp_sum( sc(:,:,ik), 1, nbase, n_start, n_end, intra_bgrp_comm)
 
  ENDDO 
#endif
  !
  
  DO ik=1, nk
     !
     nbase=nbasek(ik)
     kdim=kdimk(ik)
     !$acc parallel vector_length(64) 
     !$acc loop gang 
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real
        !
        hc(n,n,ik) = CMPLX( REAL( hc(n,n,ik) ), 0.D0 ,kind=DP)
        sc(n,n,ik) = CMPLX( REAL( sc(n,n,ik) ), 0.D0 ,kind=DP)
        !
        !$acc loop vector 
        DO m = n + 1, nbase
           !
           hc(n,m,ik) = CONJG( hc(m,n,ik) )
           sc(n,m,ik) = CONJG( sc(m,n,ik) )
           !
        END DO
        !
     END DO
     !$acc end parallel
     !
  END DO
  !
  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
     DO ik = 1, nk
        nbase=nbasek(ik)
        CALL dev_memset(vc(:,:,ik), ZERO, (/1, nbase/), 1, &
                                          (/1, nbase/), 1)
     ENDDO
     !
     DO ik = 1, nk
        !
        nbase=nbasek(ik)
        st_=st(ik)
        !$acc parallel loop 
        DO n = 1, nbase
           !
           e(n+st_) = REAL( hc(n,n,ik) )
           !
           vc(n,n,ik) = ONE
           !
        END DO
        !$acc end parallel
        !
     END DO
     !
     CALL mp_bcast( e, root_bgrp_id, inter_bgrp_comm )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
#if defined(__CUDA)
     
     DO ib=1, nblock
        CALL diago_dev<<<times(ib),1>>>(times(ib),nvecx,              &
            nbasek_d(start(ib)),nvec,outk_d(start(ib)),hc(1,1,start(ib)),  &
            sc(1,1,start(ib)),ew(1,start(ib)),vc(1,1,start(ib)),      &
            work_d(1,start(ib)),lwork,hdiago_d(1,start(ib)),           &
            sdiago_d(1,start(ib)),rwork_d(1,start(ib)),lrwork,        &
            iwork_d(1,start(ib)), liwork, info_d(start(ib)),          &
            ifail_d(1,start(ib)), m_d(start(ib)) )
        ierr=cudaDeviceSynchronize()
     ENDDO
#else
     IF ( my_bgrp_id == root_bgrp_id ) THEN
        DO ik=1, nk
           nbase=nbasek(ik)
           CALL diaghg( nbase, nvec, hc(:,:,ik), sc(:,:,ik), nvecx, &
                     ew(:,ik), vc(:,:,ik), me_bgrp, root_bgrp, intra_bgrp_comm )
        ENDDO
     END IF
#endif
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
     DO ik=1, nk
        st_=st(ik)
        CALL dev_memcpy (e(st_+1:), ew(:,ik), (/ 1, nvec /), 1 )
     ENDDO
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
!     WRITE(stdout,*) 'cegterg iteration: ', kter
     CALL start_clock( 'cegterg:update' )
#if defined (__CUDA)
     conv_d=conv
     e_d=e
     CALL cegterg_upd0<<<nk,1>>>(outk_d, nbasek_d, st_d, dav_iter_d, conv_d, &
                           nb1k_d, vc, ew, e_d, nvecx, nvec, kter, nk)
     ierr=cudaDeviceSynchronize()

     nb1k=nb1k_d
     dav_iter=dav_iter_d
#else
     DO ik = 1, nk
        !
        IF (outk(ik)) CYCLE
        !
        nbase=nbasek(ik)
        st_=st(ik)
        !
        dav_iter(ik) = kter ; !write(*,*) kter, notcnv, conv
        !
        np = 0
        !
        DO n = 1, nvec
           !
           IF ( .NOT. conv(st_+n) ) THEN
              !
              ! ... this root not yet converged ... 
              !
              np = np + 1
              !
              ! ... reorder eigenvectors so that coefficients for unconverged
              ! ... roots come first. This allows to use quick matrix-matrix 
              ! ... multiplications to set a new basis vector (see below)
              !
              IF ( np /= n ) THEN
                 DO i = 1, nvecx
                    vc(i,np,ik) = vc(i,n,ik)
                 END DO
              END IF
              !
              ! ... for use in g_psi
              !
              ew(nbase+np,ik) = e(st_+n)
              !
           END IF
           !
        END DO
        !
        nb1k(ik) = nbase + 1
        !
    END DO
#endif
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
  !
     DO ik = 1, nk
        !
        IF (outk(ik)) CYCLE
        !
        nbase=nbasek(ik) 
        kdim=kdimk(ik) 
        nb1=nb1k(ik)
        notcnv=notcnvk(ik)
        stx_=stx(ik)
        !
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        !
        IF ( uspp ) THEN
           !
           IF (n_start .LE. n_end) &
              CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE,   &
              spsi(1,stx_+n_start), kdmx, vc(n_start,1,ik), nvecx, &
                    ZERO, psi(1,stx_+nb1), kdmx )
           !     
        ELSE
           !
           IF (n_start .LE. n_end) &
              CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, &
                    psi(1, stx_+n_start), kdmx, vc(n_start,1,ik), nvecx, &
                    ZERO, psi(1,stx_+nb1), kdmx )
           !
        END IF
        !
     END DO
     
! NB: must not call mp_sum over inter_bgrp_comm here because it is done later to the full correction
     !
#if defined(__CUDA)
     !$acc parallel loop gang copyin(outk,nbasek,notcnvk,stx)
     DO ik = 1, nk
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        notcnv=notcnvk(ik)
        npw_=npw(ik)
        stx_=stx(ik)
        !$acc loop collapse(3) 
        DO np = 1, notcnv
           DO ipol = 1, npol
              DO k=1,npw_
                 psi(k + (ipol-1)*npwx, stx_+nbase+np) = &
                    - ew(nbase+np,ik)* &
                    psi(k + (ipol-1)*npwx, stx_+nbase+np)
              END DO
           END DO
        END DO
     END DO
     !$acc end parallel 
#else
     DO ik=1, nk
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        !$omp parallel do collapse(3)
        DO n = 1, notcnvk(ik)
           DO ipol = 1, npol
              DO m = 1, numblock(ik)
                 psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw(ik), m*blocksize)+(ipol-1)*npwx, &
                            stx(ik)+nbase+n) = - ew(nbase+n,ik) * &
                 psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw(ik), m*blocksize)+(ipol-1)*npwx, &
                                stx(ik)+nbase+n)
              END DO
           END DO
        END DO
        !$omp end parallel do
     END DO
#endif
     !
#if defined(__CUDA)
     CALL cegterg_upd3<<<dim3(kdmx/32+1,nvec/32+1,nk),dim3(32,32,1)>>>&
         (outk_d, nbasek_d, kdimk_d, notcnvk_d, nb1k_d, stx_d, hpsi, &
          psi, vc, kdmx, nvecx, nk)
     ierr=cudaDeviceSynchronize()
#else
     DO ik=1, nk
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        kdim=kdimk(ik)
        nb1=nb1k(ik)
        notcnv=notcnvk(ik)
        stx_=stx(ik)
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end

        IF (n_start .LE. n_end) &
           CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, &
                 hpsi(1, stx_+n_start), kdmx, vc(n_start,1,ik), nvecx, &
                 ONE, psi(1,stx_+nb1), kdmx )
        CALL mp_sum( psi(:,stx_+nb1:stx(ik)+nbase+notcnv), inter_bgrp_comm )
     ENDDO
#endif
     CALL stop_clock( 'cegterg:update' )
        !
        ! clean up garbage if there is any
        !
     DO ik=1,nk
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        nb1=nb1k(ik)
        notcnv=notcnvk(ik)
        stx_=stx(ik)
        IF (npw(ik) < npwx) CALL dev_memset(psi, ZERO, [npw(ik)+1,npwx], 1, &
                  [stx_+nb1, stx_+nbase+notcnv])
        IF (npol == 2)  CALL dev_memset(psi, ZERO, [npwx+npw(ik)+1,2*npwx], 1,&
                  [stx_+nb1, stx_+nbase+notcnv])
        !
        ! ... approximate inverse iteration
        !
#if defined (__CUDA)
     ENDDO
     
     CALL g_psii<<<dim3(npwx,nvec,nk),dim3(1,1,1)>>>&
         ( npwx, outk_d, npw_d, notcnvk_d, nb1k_d, stx_d, npol, psi, &
                                                       ew, nvecx, nk )
     ierr=cudaDeviceSynchronize()
#else
        !
        ! ... approximate inverse iteration
        !
        CALL g_psii( npwx, npw(ik), notcnv, npol, psi(1,stx_+nb1), &
                     ew(nb1,ik), ik )
     ENDDO
#endif
     !
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
#if defined(__CUDA)
     CALL compute_dot_ew<<<dim3(nk,nvec,1),dim3(1,1,1)>>>(outk_d, npw_d, &
                nbasek_d, notcnvk_d, stx_d, psi, ew, npol, npwx, nvecx, nk)
     ierr=cudaDeviceSynchronize()
#else
     DO ik=1, nk
        !
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        notcnv=notcnvk(ik)
        stx_=stx(ik)
        !
        !$acc parallel vector_length(96) 
        !$acc loop gang private(nbn)
        DO n = 1, notcnv
           !
           nbn = stx_ + nbase + n
           !
           ew(n,ik) = MYDDOT_VECTOR_GPU( 2*npw(ik), psi(1,nbn), psi(1,nbn) )
           !
        END DO
        !$acc end parallel
        !
     END DO
     !
     IF (npol.ne.1) THEN 
        DO ik=1, nk
           !
           IF (outk(ik)) CYCLE
           !
           nbase=nbasek(ik)
           stx_=stx(ik)
           notcnv=notcnvk(ik)
           !
           !$acc loop gang private(nbn)
           DO n = 1, notcnv
              nbn = stx_ + nbase + n
              ew(n,ik) = ew(n,ik) + MYDDOT_VECTOR_GPU( 2*npw(ik), &
                                    psi(npwx+1,nbn), psi(npwx+1,nbn) ) 
           END DO
           !
        END DO 
        !
     END IF 
#endif
     !
     DO ik=1, nk
        IF (outk(ik)) CYCLE
        notcnv=notcnvk(ik)
        CALL mp_sum( ew( 1:notcnv, ik ), intra_bgrp_comm )
     ENDDO
     !
#if defined(__CUDA)
     !$acc parallel loop gang 
     DO ik=1, nk
        IF (outk_d(ik)) CYCLE
        nbase=nbasek_d(ik)
        notcnv=notcnvk_d(ik)
        npw_=npw(ik)
        stx_=stx_d(ik)
        !$acc loop vector collapse(3)
        DO i = 1, notcnv
           DO ipol = 1,npol
              DO k=1,npw_
                 psi(k + (ipol-1)*npwx,stx_+nbase+i) = &
                 psi(k+(ipol-1)*npwx,stx_+nbase+i)/SQRT( ew(i,ik) )
              END DO
           END DO
        END DO
     END DO
     !$acc end parallel
#else
     !$omp parallel do collapse(4)
     DO ik = 1, nk 
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        notcnv=notcnvk(ik)
        stx_=stx(ik)
        DO n = 1, notcnv
           DO ipol = 1, npol
              DO m = 1, numblock(ik)
                 psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw(ik), m*blocksize)+(ipol-1)*npwx, &
                                             stx_+nbase+n) = &
                 psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw(ik), m*blocksize)+(ipol-1)*npwx, &
                    stx_+nbase+n) / SQRT( ew(n,ik) )
              END DO
           END DO
        END DO
     END DO
     !$omp end parallel do
#endif
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
#if defined (__CUDA)
     IF (uspp) THEN
        CALL h_s_psik_dev(npwx, outk_d, kdimk_d, npw_d, notcnvk_d, nb1k_d, &
             st_d, stx_d, ikt_d, ikblk_d, npol, psi, hpsi, spsi, psicmr,   &
             nvec, nvecx, nnrs, nk, .FALSE.)
     ELSE
        CALL h_psik_dev(npwx, outk_d, kdimk_d, npw_d, notcnvk_d, nb1k_d, &
             st_d, stx_d, ikt_d, ikblk_d, npol, psi, hpsi, psicmr, nvec, &
             nvecx, nnrs, nk, .FALSE.)
     ENDIF
#else
     DO ik=1, nk
        !
        IF (outk(ik)) CYCLE
        nb1=nb1k(ik)
        notcnv=notcnvk(ik)
        stx_=stx(ik)
        !
        CALL h_psii( npwx, npw(ik), notcnv, psi(1,stx_+nb1), &
                   hpsi(1, stx_+nb1), ik ) 
        !
        IF ( uspp ) CALL s_psii( npwx, npw(ik), notcnv, &
                   psi(1,stx_ + nb1), spsi(1,stx_+nb1), ik )
        !
     ENDDO
#endif
     CALL start_clock( 'cegterg:overlap' )
     DO ik=1,nk
        IF (outk(ik)) CYCLE
        nhpsi(ik) = nhpsi(ik) + notcnvk(ik)
     ENDDO
#if defined(__CUDA)
     IF (uspp) THEN
        CALL cegterg_overlap_us<<<dim3(nvec/4+1,nvecx/4+1,nk/32+1),     &
          dim3(4,4,32)>>>(outk_d, nbasek_d, kdimk_d, notcnvk_d, nb1k_d, &
          stx_d, hpsi, spsi, psi, hc, sc, kdmx, nvecx, nk)
     ELSE
        CALL cegterg_overlap<<<dim3(nvec/4+1,nvecx/4+1,nk/32+1),        &
          dim3(4,4,32)>>>(outk_d, nbasek_d, kdimk_d, notcnvk_d, nb1k_d, &
          stx_d, hpsi, psi, hc, sc, kdmx, nvecx, nk)
     ENDIF
     ierr=cudaDeviceSynchronize()
#else
     DO ik=1, nk
        !
        IF (outk(ik)) CYCLE
        !
        nbase=nbasek(ik)
        kdim=kdimk(ik)
        nb1=nb1k(ik)
        notcnv=notcnvk(ik)
        stx_=stx(ik)
        !
        ! ... update the reduced hamiltonian
        !
        CALL divide_all(inter_bgrp_comm,nbase+notcnv,n_start,n_end,&
                        recv_counts,displs)
        CALL mp_type_create_column_section(sc(1,1,1), nbase, notcnv, &
                        nvecx, column_section_type)
        my_n = n_end - n_start + 1 !; write (*,*) nbase+notcnv,n_start,n_end
        CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, &
                  hpsi(1, stx_+nb1), kdmx, psi(1,stx_+n_start), &
                  kdmx, ZERO, hc(nb1,n_start,ik), nvecx )
        !
        IF (n_start .LE. n_end) &
           CALL mp_sum( hc(nb1:nbase+notcnv, n_start:n_end, ik), &
                                                         intra_bgrp_comm )
        CALL mp_gather( hc(:,:,ik), column_section_type, recv_counts, displs, &
                                           root_bgrp_id, inter_bgrp_comm )
        !
        CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)

        my_n = n_end - n_start + 1!; write (*,*) nbase+notcnv,n_start,n_end
        IF ( uspp ) THEN
           !
           CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, &
                   spsi(1,stx_+nb1), kdmx, psi(1, stx_+n_start), &
                   kdmx, ZERO, sc(nb1, n_start, ik), nvecx )
           !     
        ELSE
           !
           CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, &
                   psi(1,stx_+nb1), kdmx, psi(1,stx_+n_start), &
                   kdmx, ZERO, sc(nb1,n_start,ik), nvecx )
          !
        END IF
        !
        IF (n_start .LE. n_end) &
           CALL mp_sum( sc(nb1:nbase+notcnv, n_start:n_end, ik), &
                      intra_bgrp_comm )
        CALL mp_gather( sc(:,:,ik), column_section_type, recv_counts, displs, &
                                            root_bgrp_id, inter_bgrp_comm )
        CALL mp_type_free( column_section_type )
    END DO
#endif
#if defined(__CUDA)
     DO ik=1, 0
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        nb1=nb1k(ik)
        notcnv=notcnvk(ik)
        n_start=1
        n_end=nbase+notcnv
        my_n=n_end
        IF (n_start .LE. n_end) &
           CALL mp_sum( hc(:,:,ik), nb1, nbase+notcnv, n_start, &
                                                 n_end , intra_bgrp_comm )
        IF (n_start .LE. n_end) & 
           CALL mp_sum( sc(:,:,ik), nb1, nbase+notcnv, n_start, &
                n_end, intra_bgrp_comm )
     ENDDO
#endif
     !
     DO ik=1,nk
        IF (outk(ik)) CYCLE
        nbasek(ik) = nbasek(ik) + notcnvk(ik)
     ENDDO
     !
#if defined(__CUDA)
     nbasek_d=nbasek
     CALL cegterg_herm<<<dim3(nvecx,1,nk),dim3(1,1,1)>>>(outk_d, nbasek_d, &
                       nb1k_d, hc, sc, nvecx, nk)
     ierr=cudaDeviceSynchronize()
#else
     DO ik = 1, nk
        !
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        nb1=nb1k(ik)
        !
        !$acc parallel loop vector_length(64) 
        DO n = 1, nbase
           !
           ! ... the diagonal of hc and sc must be strictly real
           !
           IF ( n>=nb1 ) THEN
              hc(n,n,ik) = CMPLX( REAL( hc(n,n,ik) ), 0.D0 ,kind=DP)
              sc(n,n,ik) = CMPLX( REAL( sc(n,n,ik) ), 0.D0 ,kind=DP)
           ENDIF
           !
           !$acc loop vector
           DO m = MAX(n+1,nb1), nbase
              !
              hc(n,m,ik) = CONJG( hc(m,n,ik) )
              sc(n,m,ik) = CONJG( sc(m,n,ik) )
              !
           END DO
           !
        END DO
        !$acc end parallel 
        !
     END DO
#endif
     CALL stop_clock( 'cegterg:overlap' )
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:diag' )
#if defined(__CUDA)     
     DO ib=1, nblock
        CALL diago_dev<<<times(ib),1>>>(times(ib),nvecx,              &   
            nbasek_d(start(ib)),nvec,outk_d(start(ib)),hc(1,1,start(ib)),  &
            sc(1,1,start(ib)),ew(1,start(ib)),vc(1,1,start(ib)),      &
            work_d(1,start(ib)),lwork,hdiago_d(1,start(ib)),          &
            sdiago_d(1,start(ib)),rwork_d(1,start(ib)),lrwork,        &
            iwork_d(1,start(ib)),liwork,info_d(start(ib)),           & 
            ifail_d(1,start(ib)), m_d(start(ib)) )
        ierr=cudaDeviceSynchronize()
     ENDDO
#else
     DO ik = 1, nk
        !
        IF (outk(ik)) CYCLE
        nbase=nbasek(ik)
        !
        IF ( my_bgrp_id == root_bgrp_id ) THEN
           CALL diaghg( nbase, nvec, hc(:,:,ik), sc(:,:,ik), nvecx, &
                  ew(:,ik), vc(:,:,ik), me_bgrp, root_bgrp, intra_bgrp_comm )
!
!    This call is for testing purpose only. It uses the lapack routines
!    included with tpw, the same routines used by diago_dev.
!    To use it without CUDA info_d and rwork_d must be declared and allocated
!    outside the if defined(__CUDA)
!
!           CALL diago_tpw(nbase,nvecx,hc(:,:,ik),sc(:,:,ik),ew(:,ik), &
!                               vc(:,:,ik),rwork_d(:,ik),lrwork,info_d(ik))

        END IF
     END DO
#endif
     IF ( nbgrp > 1 ) THEN
        CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
     ! ... test for convergence
     !
     !$acc parallel loop gang copy(conv(1:nvec*nk)) copyin(btype(1:nvec*nk),outk(1:nk),st(1:nk))
     DO ik=1, nk
        IF (outk(ik)) CYCLE
        st_=st(ik)
        !$acc loop vector
        DO i = 1, nvec
           IF (btype(st_+i) == 1) THEN
              conv(st_+i) = ( ABS( ew(i,ik) - e(st_+i) ) < ethr ) 
           ELSE
              conv(st_+i) = ( ABS( ew(i,ik) - e(st_+i) ) < empty_ethr ) 
           END IF 
        ENDDO
     END DO 
     !
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     DO ik=1, nk
        IF (outk(ik)) CYCLE
        notcnvk(ik) = COUNT( .NOT. conv(st(ik)+1:st(ik)+nvec) )
     END DO
#if defined(__CUDA)
     notcnvk_d=notcnvk
#endif
     !
     DO ik=1, nk
        IF (outk(ik)) CYCLE
        st_=st(ik)
        CALL dev_memcpy (e(st_+1:st_+nvec), ew(1:nvec,ik), &
                                                (/ 1, nvec /), 1 )
     END DO
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     !
     CALL start_clock( 'cegterg:last' )
     DO ik=1, nk
        enter(ik)=.FALSE.
        IF (outk(ik)) CYCLE
        IF ( notcnvk(ik) == 0 .OR. &
           nbasek(ik)+notcnvk(ik) > nvecx .OR. dav_iter(ik) == maxter ) &
           enter(ik)=.TRUE.
        IF ( notcnvk(ik) == 0 .OR. dav_iter(ik) == maxter ) outk(ik)=.TRUE.
!        WRITE(6,*) 'outk', ik, outk(ik), notcnvk(ik), nbasek(ik)
     ENDDO
#if defined(__CUDA)
     outk_d=outk
     enter_d=enter
#endif
     DO ik=1, nk
        IF (.NOT.enter(ik)) CYCLE
        nbase=nbasek(ik)
        kdim=kdimk(ik)
        stx_=stx(ik)
        st_=st(ik)
        !
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, &
                     psi(1,stx_+n_start), kdmx, vc(n_start,1,ik), nvecx, &
                     ZERO, evc(1,st_+1), kdmx )

        IF (nbgrp>1) CALL mp_sum( evc(:,st_+1:st_+nvec), inter_bgrp_comm )
     ENDDO

#if defined(__CUDA)
     CALL copy_psi_gpu<<<dim3(kdmx,nvec,nk/32+1),dim3(1,1,32)>>>  &
          ( npwx, outk_d, enter_d, kdimk_d, stx_d, st_d, npol, psi, evc,  &
            nvec, nvecx, nk )
     ierr=cudaDeviceSynchronize()
#endif

     DO ik=1, nk
        IF (outk(ik).OR..NOT.enter(ik)) CYCLE
        !
        nbase=nbasek(ik)
        kdim=kdimk(ik)
        stx_=stx(ik)
        st_=st(ik)
        !
        ! ... refresh psi, H*psi and S*psi
        !
#if ! defined(__CUDA)
        psi(:,stx_+1:stx_+nvec)=evc(:,st_+1:st_+nvec)
#endif
!        CALL dev_memcpy(psi(:,stx_+1:stx_+nvec), &
!                              evc(:,st_+1:st_+nvec), &
!                                  (/ 1, npwx*npol /), 1,     &
!                                  (/1, nvec /), 1)
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        IF ( uspp ) THEN
        !
           CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, &
                  spsi(1,stx_+n_start), kdmx, vc(n_start,1,ik), nvecx, &
                       ZERO, psi(1,stx_+nvec+1), kdmx)
           CALL dev_memcpy(spsi(:,stx_+1:stx_+nvec),        &
                            psi(:,stx_+nvec+1:stx_+2*nvec), &
                                     (/1, npwx*npol/), 1, &
                                     (/1, nvec/), 1)
           CALL mp_sum( spsi(:,stx_+1:stx_+nvec), inter_bgrp_comm )
           !
        END IF
        !
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE,               &
              hpsi(1,stx_+n_start), kdmx, vc(n_start,1,ik), nvecx, &
                 ZERO, psi(1,stx_+nvec+1), kdmx )
        CALL dev_memcpy(hpsi(:,stx_+1:stx_+nvec), &
                         psi(:,stx_+nvec+1:stx_+2*nvec), &
                                     (/1, npwx*npol/), 1, &
                                     (/1, nvec/), 1)
        CALL mp_sum( hpsi(:,stx_+1:stx_+nvec), inter_bgrp_comm )
        !
     END DO
     !
     ! ... refresh the reduced hamiltonian 
     !
     !$acc kernels copyin(outk,enter,st) copy(nbasek)
     DO ik=1, nk
        IF (outk(ik).OR..NOT.enter(ik)) CYCLE
        nbasek(ik) = nvec
        nbase=nbasek(ik)
        st_=st(ik)
        !
        ! These variables are set to ZERO in the CUF Kernel below
        !hc(1:nbase,1:nbase) = ZERO
        !sc(1:nbase,1:nbase) = ZERO
        !vc(1:nbase,1:nbase) = ZERO
        !
        DO n = 1, nbase
           hc(n,n,ik) = CMPLX( e(st_+n), 0.0_DP ,kind=DP)
           sc(n,n,ik) = ONE
           vc(n,n,ik) = ONE
           DO j = n+1, nbase
              hc(j,n,ik) = ZERO
              hc(n,j,ik) = ZERO
              sc(j,n,ik) = ZERO
              sc(n,j,ik) = ZERO
              vc(j,n,ik) = ZERO
              vc(n,j,ik) = ZERO
           END DO
           !
        END DO
        !
     END DO
     !$acc end kernels
#if defined(__CUDA)
     nbasek_d=nbasek
#endif
     !
     CALL stop_clock( 'cegterg:last' )
!
!   now count how many ik are converged and exit if they are nk
!
     done_ik = COUNT(outk(1:nk))
     IF (done_ik==nk) EXIT iterate
     !
  END DO iterate
  !
  DEALLOCATE( recv_counts )
  DEALLOCATE( displs )

  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )

  DEALLOCATE(nbasek)
  DEALLOCATE(kdimk)
  DEALLOCATE(nb1k)
  DEALLOCATE(outk)
  DEALLOCATE(enter)
  DEALLOCATE(st)
  DEALLOCATE(stx)
  DEALLOCATE(ikt)
  DEALLOCATE(ikblk)

#if defined(__CUDA)
  DEALLOCATE(nbasek_d)
  DEALLOCATE(kdimk_d)
  DEALLOCATE(npw_d)
  DEALLOCATE(nb1k_d)
  DEALLOCATE(outk_d)
  DEALLOCATE(enter_d)
  DEALLOCATE(st_d)
  DEALLOCATE(stx_d)
  DEALLOCATE(ikt_d)
  DEALLOCATE(ikblk_d)
  DEALLOCATE(conv_d)
  DEALLOCATE(notcnvk_d)

  DEALLOCATE(work_d)
  DEALLOCATE(rwork_d)
  DEALLOCATE(iwork_d)
  DEALLOCATE(ifail_d)
  DEALLOCATE(m_d)
  DEALLOCATE(hdiago_d)
  DEALLOCATE(sdiago_d)
  DEALLOCATE(info_d)

  DEALLOCATE(dav_iter_d)
  DEALLOCATE(e_d)
  DEALLOCATE(psicmr)

  DEALLOCATE(times)
  DEALLOCATE(start)
#endif
  !
#if ! defined(__CUDA)
  ! compute the number of chuncks
  DEALLOCATE(numblock)
#endif
  !$acc end data 
  !
  CALL stop_clock( 'cegterg' )!; write(*,*) 'stop cegterg' ; FLUSH(6)
!  call print_clock( 'cegterg' )
!  call print_clock( 'cegterg:init' )
!  call print_clock( 'cegterg:diag' )
!  call print_clock( 'cegterg:update' )
!  call print_clock( 'cegterg:overlap' )
!  call print_clock( 'cegterg:last' )
  !
  RETURN
  !
END SUBROUTINE cegterg_vk
