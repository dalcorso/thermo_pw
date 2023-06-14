! Copyright (C) 2023 Dal Corso Andrea
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE many_k_ph_mod
!
!   This module contains the device variables of the phonon code
!   needed to compute simultaneously many k points on device.  
!   It contains also the variables needed to divide the k points
!   of the phonon code in blocks.
!
#if defined(__CUDA)
USE cudafor
#endif
USE kinds,     ONLY : DP
USE io_global, ONLY : stdout

IMPLICIT NONE
PRIVATE
!
!   Variables to divide the k points in blocks
!
INTEGER :: current_ikb_ph  ! the block that is currently running

INTEGER :: nkblocks_ph   ! number of k blocks
INTEGER :: nksbx_ph      ! maximum number of k per block
INTEGER, ALLOCATABLE :: nksb_ph(:)      ! number of k points per block
INTEGER, ALLOCATABLE :: startkb_ph(:)   ! last k of the previous block
INTEGER, ALLOCATABLE :: startkb_ph_d(:) ! the same on device
!
!  variables for putting on device many right hand sides and many
!  derivatives of the wavefunctions.
!
COMPLEX(DP), ALLOCATABLE :: dvpsik_d(:,:)
COMPLEX(DP), ALLOCATABLE :: dpsik_d(:,:)
COMPLEX(DP), ALLOCATABLE :: evqk_d(:,:)
REAL(DP),    ALLOCATABLE :: h_diagk_ph_d(:,:)
COMPLEX(DP), ALLOCATABLE :: becp1k_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: alphak_d(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: becptk_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: alphatk_d(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: dbecq_d(:,:,:)
COMPLEX(DP), ALLOCATABLE :: ps1k(:,:,:)
COMPLEX(DP), ALLOCATABLE :: ps1k_d(:,:,:)
COMPLEX(DP), ALLOCATABLE :: ps2k(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: ps2k_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: ps1k_nc_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: ps2k_nc_d(:,:,:,:,:)
!
!  Variables needed to put on device phonon variables on host
!
REAL(DP),    ALLOCATABLE :: deff_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: deff_nc_d(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: int1_d(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: int2_d(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: int2_so_d(:,:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: int1_nc_save_d(:,:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: int1_nc_d(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: deeq_nc_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: deeq_nc_save_d(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: dbecsum_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: u_d(:,:)
INTEGER,     ALLOCATABLE :: nbnd_occ_d(:) 
INTEGER,     ALLOCATABLE :: ikks_d(:) 
INTEGER,     ALLOCATABLE :: ikqs_d(:) 

#if defined(__CUDA)
ATTRIBUTES(DEVICE) :: startkb_ph_d
ATTRIBUTES(DEVICE) :: dvpsik_d
ATTRIBUTES(DEVICE) :: dpsik_d
ATTRIBUTES(DEVICE) :: evqk_d
ATTRIBUTES(DEVICE) :: h_diagk_ph_d
ATTRIBUTES(DEVICE) :: dbecq_d
ATTRIBUTES(DEVICE) :: becp1k_d
ATTRIBUTES(DEVICE) :: alphak_d
ATTRIBUTES(DEVICE) :: becptk_d
ATTRIBUTES(DEVICE) :: alphatk_d
ATTRIBUTES(DEVICE) :: ps1k_d
ATTRIBUTES(DEVICE) :: ps2k_d
ATTRIBUTES(DEVICE) :: ps1k_nc_d
ATTRIBUTES(DEVICE) :: ps2k_nc_d

ATTRIBUTES(DEVICE) :: deff_d
ATTRIBUTES(DEVICE) :: deff_nc_d
ATTRIBUTES(DEVICE) :: int1_d
ATTRIBUTES(DEVICE) :: int2_d
ATTRIBUTES(DEVICE) :: int1_nc_d
ATTRIBUTES(DEVICE) :: int2_so_d
ATTRIBUTES(DEVICE) :: int1_nc_save_d
ATTRIBUTES(DEVICE) :: deeq_nc_d
ATTRIBUTES(DEVICE) :: deeq_nc_save_d
ATTRIBUTES(DEVICE) :: dbecsum_d
ATTRIBUTES(DEVICE) :: u_d
ATTRIBUTES(DEVICE) :: nbnd_occ_d
ATTRIBUTES(DEVICE) :: ikks_d
ATTRIBUTES(DEVICE) :: ikqs_d
#endif
!
! variables needed to divide in blocks the phonon k points
!
PUBLIC nkblocks_ph, nksbx_ph, nksb_ph, startkb_ph, startkb_ph_d, &
       current_ikb_ph
!
! variables needed to put on device many of the quantities for many k points
!
PUBLIC dvpsik_d, dpsik_d, evqk_d, h_diagk_ph_d, becp1k_d, alphak_d, dbecq_d, &
       ps1k, ps1k_d, ps2k, ps2k_d, ps1k_nc_d, ps2k_nc_d, becptk_d, alphatk_d
!
!  Variables needed to copy on device the phonon variables
!
PUBLIC deff_d, deff_nc_d, int1_d, int2_d, int1_nc_d, int2_so_d, dbecsum_d, &
       u_d, nbnd_occ_d, ikks_d, ikqs_d, deeq_nc_d, deeq_nc_save_d,         &
       int1_nc_save_d
!
!  Routines of this module
!
PUBLIC allocate_many_k_ph, deallocate_many_k_ph, init_k_blocks_ph, &
       prepare_ph_device, copy_alphap_becp_device

CONTAINS
!
!------------------------------------------------------------------
SUBROUTINE init_k_blocks_ph(npwx,npol,nksq,nbnd,nspin,nhm,nkb,nat,nnrs,&
                            npe,nsolv)
!------------------------------------------------------------------
!
!   This routine receives the maximum number of k+q+G vectors, the number
!   of bands, the number of k points and the maximum memory available
!   on the GPU and divide the k points in blocks so that cgsolve_all 
!   does not run out of memory. The memory used by 
!   cgsolve_all is estimated in bytes while memgpu is given in Gbytes.
!
USE many_k_mod, ONLY : memgpu
USE gvect,      ONLY : ngm
USE ions_base,  ONLY : ntyp => nsp
USE noncollin_module, ONLY : nspin_mag

IMPLICIT NONE
INTEGER :: npwx, npol, nksq, nbnd, nnrs, nspin, nhm, nkb, nat, npe, nsolv
REAL(DP) :: totmem, fftmem, fixmem

INTEGER :: fact, i, rest

fact=npe*nsolv*3+1
fftmem=DBLE(nnrs) * DBLE(nbnd) * DBLE(nksq) * DBLE (npe) *DBLE(nsolv) * 16

fixmem= DBLE(nnrs) * DBLE(nbnd) * DBLE(npol) * 16 +   & ! psic
        DBLE(npwx) * DBLE(npol) * DBLE (nbnd) * 16 * 2 +  & ! evc + evq
        DBLE(npwx) * DBLE(npol) * DBLE (nbnd) * 16 * 2 +  & ! dpsi + dvpsi
        DBLE(npwx) * DBLE(nkb) * 16 +  &  ! vkb
        DBLE(npwx) * 16  +             &  ! g2kin
        DBLE(10 * ngm) * 8                ! g, gg, eigts
IF (npol==2) THEN
   fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(nspin) * DBLE(nat+ntyp) * 16
                                               ! deeq_nc and qq_so
   fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(nspin)*DBLE(nat+1) &
                 * DBLE(nat) * 3 *16 ! int1_nc_d + int2_so_d
ELSE
   fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(nspin+1) * DBLE(nat) * 8
                                               ! deeq and qq_at
ENDIF
fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(nspin+nat) * DBLE(nat) * 3 *16
                                               ! int1_d + int2_d

totmem=DBLE(npwx*npol) * DBLE(nbnd) * DBLE(nksq) * fact * 16      &
                          ! dvpsik_d, dpsik_d, h_diagk_ph_d, evqk_d
       + DBLE(nkb) * DBLE(nbnd) * DBLE(nksq) * (npe*nsolv + 4*npol + &
                          4 * npol * npe * nsolv )&
                          ! dbecq_d, beco1k_d, alphak_d, ps1k_nc_d, ps2k_nc_d
       + DBLE(nhm) * DBLE(nhm) * DBLE(nat) * DBLE(nbnd)  * DBLE(nksq) *       &
         DBLE(nspin) & !  deff_d  of deff_nc_d         
       + DBLE((nhm*(nhm+1))/2) * DBLE(nat)* DBLE(nspin_mag)* &
         DBLE(nbnd*nksq*npe*nsolv) & ! dbecsum_d
       + DBLE(npwx*npol) * DBLE(nbnd) * DBLE(nksq) * DBLE(npe*nsolv)*16*7 &
            ! dpsi, hpsi, spsi, g, t, h, hold
       + 2.0_DP*fftmem    ! psicmr + dpsicrm
WRITE(stdout,'(5x,"Estimated memory in cgsolve_all",f10.1," MB")') totmem/1.D6 
IF (fixmem > (memgpu * 1.D9)) THEN
   WRITE(stdout,'(5x,"memgpu must be greater than ",f8.2," to use many_k")') &
               5.0_DP * fixmem/1.D9
!
!  in this case put one k point per block
!
   nkblocks_ph=nksq
ELSE
   nkblocks_ph=MIN(FLOOR(totmem / (memgpu * 1.D9-fixmem)) + 1, nksq) 
ENDIF

ALLOCATE(nksb_ph(nkblocks_ph))
ALLOCATE(startkb_ph(nkblocks_ph))

rest=MOD(nksq,nkblocks_ph)
nksb_ph=nksq/nkblocks_ph
DO i=1,rest
   nksb_ph(i)=nksb_ph(i)+1
ENDDO
startkb_ph(1)=0
DO i=2,nkblocks_ph
   startkb_ph(i)=startkb_ph(i-1)+nksb_ph(i-1)
ENDDO

nksbx_ph=0
DO i=1,nkblocks_ph
   nksbx_ph=MAX(nksbx_ph,nksb_ph(i))
ENDDO

WRITE(stdout,'(/,5x,"Many-k conjugate gradient")')
WRITE(stdout,'(5x,i8," k-points divided into ",i5," blocks")') nksq, &
                                                      nkblocks_ph
WRITE(stdout,'(5x,"Maximum size of the blocks",i10)') nksbx_ph

WRITE(stdout,'(5x,"Block number",7x,"size",7x,"start",8x,"end")') 
DO i=1,nkblocks_ph-1
   WRITE(stdout,'(5x,i7,8x,i7,5x,i7,5x,i7)') i, nksb_ph(i), startkb_ph(i)+1, &
                                             startkb_ph(i+1)
ENDDO
WRITE(stdout,'(5x,i7,8x,i7,5x,i7,5x,i7,/)') i, nksb_ph(nkblocks_ph), &
                                             startkb_ph(nkblocks_ph)+1, nksq
RETURN
END SUBROUTINE init_k_blocks_ph
!
!---------------------------------------------------------------------
SUBROUTINE prepare_ph_device(nsolv)
!---------------------------------------------------------------------
USE modes,        ONLY : u
USE control_lr,   ONLY : nbnd_occ
USE qpoint,       ONLY : ikks, ikqs
USE phus,         ONLY : int1, int2, int1_nc, int2_so
USE uspp,         ONLY : deeq_nc
USE nc_mag_aux,   ONLY : int1_nc_save, deeq_nc_save
USE noncollin_module, ONLY : noncolin
USE uspp,         ONLY : okvan

IMPLICIT NONE
INTEGER :: nsolv

#if defined(__CUDA)
IF (okvan) THEN
   IF (noncolin) THEN
      IF (nsolv==2) THEN
         int1_nc_save_d=int1_nc_save
         deeq_nc_save_d=deeq_nc_save
      ELSE
         int1_nc_d=int1_nc
         deeq_nc_d=deeq_nc
      ENDIF
      int2_so_d=int2_so
      int2_d=int2
   ELSE
      int1_d=int1
      int2_d=int2
   ENDIF
ENDIF
u_d=u
nbnd_occ_d=nbnd_occ
startkb_ph_d=startkb_ph
ikqs_d=ikqs
ikks_d=ikks
#endif

RETURN
END SUBROUTINE prepare_ph_device
!
!------------------------------------------------------------------
SUBROUTINE allocate_many_k_ph(npe, nsolv)
!------------------------------------------------------------------
USE wvfct,            ONLY : nbnd, npwx
USE noncollin_module, ONLY : npol, nspin_mag, noncolin
USE ions_base,        ONLY : nat
USE uspp_param,       ONLY : nhm
USE lsda_mod,         ONLY : nspin
USE qpoint,           ONLY : nksq
USE klist,            ONLY : nks
USE uspp,             ONLY : nkb

IMPLICIT NONE
INTEGER, INTENT(IN) :: npe, nsolv

 ALLOCATE(dvpsik_d(npwx*npol,nbnd*nksbx_ph*npe*nsolv))
 ALLOCATE(dpsik_d(npwx*npol,nbnd*nksbx_ph*npe*nsolv))
 ALLOCATE(evqk_d(npwx*npol,nbnd*nksbx_ph*nsolv))
 ALLOCATE(h_diagk_ph_d(npwx*npol,nbnd*nksbx_ph*npe*nsolv))
 ALLOCATE(ps1k(nkb,nbnd,nksbx_ph*npe*nsolv))
 ALLOCATE(ps2k(nkb,nbnd,3,nksbx_ph*npe*nsolv))

#if defined(__CUDA) 
 ALLOCATE(startkb_ph_d(nkblocks_ph))

 ALLOCATE(dbecq_d(nkb,nbnd,nksbx_ph*npe*nsolv))
 ALLOCATE(becp1k_d(nkb,npol,nbnd,nksbx_ph))
 ALLOCATE(alphak_d(nkb,npol,nbnd,3,nksbx_ph))
 IF (nsolv==2) THEN
    ALLOCATE(becptk_d(nkb,npol,nbnd,nksbx_ph))
    ALLOCATE(alphatk_d(nkb,npol,nbnd,3,nksbx_ph))
 ENDIF
 ALLOCATE(ps1k_d(nkb,nbnd,nksbx_ph*npe*nsolv))
 ALLOCATE(ps2k_d(nkb,nbnd,3,nksbx_ph*npe*nsolv))
 ALLOCATE(ps1k_nc_d(nkb,npol,nbnd,nksbx_ph*npe*nsolv))
 ALLOCATE(ps2k_nc_d(nkb,npol,nbnd,3,nksbx_ph*npe*nsolv))
!
!  variables needed to copy the phonon variables
!
 IF (noncolin) THEN
    ALLOCATE(deff_nc_d(nhm, nhm, nat, nspin, nbnd*nksbx_ph*npe*nsolv))
    IF (nsolv==2) THEN
       ALLOCATE(int1_nc_save_d(nhm, nhm, 3, nat, nspin, nsolv))
    ELSE
       ALLOCATE(int1_nc_d(nhm, nhm, 3, nat, nspin))
    ENDIF
    ALLOCATE(int2_so_d(nhm, nhm, 3, nat, nat, nspin))
    ALLOCATE(deeq_nc_save_d( nhm, nhm, nat, nspin, nsolv))
 ELSE
    ALLOCATE(deff_d(nhm, nhm, nat, nbnd*nksbx_ph*npe*nsolv))
 ENDIF
 ALLOCATE(int1_d(nhm, nhm, 3, nat, nspin))
 ALLOCATE(int2_d(nhm, nhm, 3, nat, nat))
 ALLOCATE(dbecsum_d((nhm*(nhm+1))/2, nat, nspin_mag, nbnd*nksbx_ph*npe*nsolv))
 ALLOCATE(u_d(3*nat, 3*nat))
 ALLOCATE(nbnd_occ_d(nks))
 ALLOCATE(ikks_d(nksq))
 ALLOCATE(ikqs_d(nksq))
#endif

RETURN
END SUBROUTINE allocate_many_k_ph
!
!------------------------------------------------------------------
SUBROUTINE deallocate_many_k_ph()
!------------------------------------------------------------------
IMPLICIT NONE

IF (ALLOCATED(nksb_ph)) DEALLOCATE(nksb_ph)
IF (ALLOCATED(startkb_ph)) DEALLOCATE(startkb_ph)
IF (ALLOCATED(startkb_ph_d)) DEALLOCATE(startkb_ph_d)

IF (ALLOCATED(dvpsik_d)) DEALLOCATE(dvpsik_d)
IF (ALLOCATED(dpsik_d)) DEALLOCATE(dpsik_d)
IF (ALLOCATED(evqk_d)) DEALLOCATE(evqk_d)
IF (ALLOCATED(h_diagk_ph_d)) DEALLOCATE(h_diagk_ph_d)

IF (ALLOCATED(dbecq_d)) DEALLOCATE(dbecq_d)
IF (ALLOCATED(becp1k_d)) DEALLOCATE(becp1k_d)
IF (ALLOCATED(alphak_d)) DEALLOCATE(alphak_d)
IF (ALLOCATED(becptk_d)) DEALLOCATE(becptk_d)
IF (ALLOCATED(alphatk_d)) DEALLOCATE(alphatk_d)
IF (ALLOCATED(ps1k)) DEALLOCATE(ps1k)
IF (ALLOCATED(ps1k_d)) DEALLOCATE(ps1k_d)
IF (ALLOCATED(ps2k)) DEALLOCATE(ps2k)
IF (ALLOCATED(ps2k_d)) DEALLOCATE(ps2k_d)
IF (ALLOCATED(ps1k_nc_d)) DEALLOCATE(ps1k_nc_d)
IF (ALLOCATED(ps2k_nc_d)) DEALLOCATE(ps2k_nc_d)

IF (ALLOCATED(deff_d)) DEALLOCATE(deff_d)
IF (ALLOCATED(deff_nc_d)) DEALLOCATE(deff_nc_d)
IF (ALLOCATED(int1_d)) DEALLOCATE(int1_d)
IF (ALLOCATED(int2_d)) DEALLOCATE(int2_d)
IF (ALLOCATED(int1_nc_d)) DEALLOCATE(int1_nc_d)
IF (ALLOCATED(int2_so_d)) DEALLOCATE(int2_so_d)
IF (ALLOCATED(int1_nc_save_d)) DEALLOCATE(int1_nc_save_d)
IF (ALLOCATED(deeq_nc_save_d)) DEALLOCATE(deeq_nc_save_d)

IF (ALLOCATED(dbecsum_d)) DEALLOCATE(dbecsum_d)
IF (ALLOCATED(u_d)) DEALLOCATE(u_d)
IF (ALLOCATED(nbnd_occ_d)) DEALLOCATE(nbnd_occ_d)
IF (ALLOCATED(ikks_d)) DEALLOCATE(ikks_d)
IF (ALLOCATED(ikqs_d)) DEALLOCATE(ikqs_d)

RETURN
END SUBROUTINE deallocate_many_k_ph
!
!------------------------------------------------------------------
SUBROUTINE copy_alphap_becp_device(ikb, nsolv)
!------------------------------------------------------------------
USE phus,          ONLY : alphap
USE lrus,          ONLY : becp1
USE qpoint_aux,    ONLY : becpt, alphapt
USE uspp,          ONLY : nkb
USE wvfct,         ONLY : nbnd
USE noncollin_module, ONLY : noncolin, npol
IMPLICIT NONE
INTEGER :: ikb

INTEGER :: ik, ik1, ipol, nsolv
!
#if defined(__CUDA)
DO ik=startkb_ph(ikb)+1, startkb_ph(ikb)+nksb_ph(ikb)
   ik1=ik-startkb_ph(ikb)
   IF (noncolin) THEN
      becp1k_d(1:nkb,1:npol,1:nbnd,ik1)=becp1(ik)%nc(1:nkb,1:npol,1:nbnd)
      DO ipol=1,3
         alphak_d(1:nkb,1:npol,1:nbnd,ipol,ik1)=&
                             alphap(ipol,ik)%nc(1:nkb,1:npol,1:nbnd)
      ENDDO
      IF (nsolv==2) THEN
         becptk_d(1:nkb,1:npol,1:nbnd,ik1)=becpt(ik)%nc(1:nkb,1:npol,1:nbnd)
         DO ipol=1,3
            alphatk_d(1:nkb,1:npol,1:nbnd,ipol,ik1)=&
                             alphapt(ipol,ik)%nc(1:nkb,1:npol,1:nbnd)
         ENDDO
      ENDIF
   ELSE
      becp1k_d(1:nkb,1,1:nbnd,ik1)=becp1(ik)%k(1:nkb,1:nbnd)
      DO ipol=1,3
         alphak_d(1:nkb,1,1:nbnd,ipol,ik1)=alphap(ipol,ik)%k(1:nkb,1:nbnd)
      ENDDO
   ENDIF
ENDDO
#endif

RETURN
END SUBROUTINE copy_alphap_becp_device
!
END MODULE many_k_ph_mod
