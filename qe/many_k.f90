! Copyright (C) 2023 Dal Corso Andrea
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
MODULE many_k_mod
!-------------------------------------------------------------------------
!
!   This module declares the device quantities needed to compute
!   simultaneously many k points on the GPU. Moreover, it copies on device
!   variables of the pw.x code that are not already there. Finally it
!   allocates the variable needed to make the fft on device.
!   It has also variables to divide the k points in blocks so that
!   the amount of memory used in the GPU do not exceeds its actual
!   capacity.
!
#if defined(__CUDA)
USE cudafor
#endif
USE kinds, ONLY : DP
USE io_global, ONLY : stdout

IMPLICIT NONE
PRIVATE
!
REAL(DP) :: memgpu       ! available memory for cegterg (in GBytes)
                         ! given in input
!
!   Variables to divide the k points in blocks.
!
INTEGER :: nkblocks      ! number of k blocks
INTEGER :: nksbx         ! maximum number of k per block
INTEGER :: current_ikb   ! current_block
INTEGER, ALLOCATABLE :: nksb(:)   ! number of k points per block
INTEGER, ALLOCATABLE :: startkb(:) ! last k of the previous block

INTEGER, ALLOCATABLE :: startkb_d(:) ! last k of the previous block on device

LOGICAL :: alloc_many_k  ! If there is a single block do not recompute
!
!  fft on device variables
!
REAL(DP), ALLOCATABLE :: work_d(:,:)
REAL(DP), ALLOCATABLE :: wsave_d(:,:)
INTEGER,  ALLOCATABLE :: nl_d(:)      
INTEGER,  ALLOCATABLE :: isindex_d(:)  
INTEGER,  ALLOCATABLE :: iplane_d(:)
INTEGER :: lenwrk, lensav
INTEGER :: lenwrk_d, lensav_d
!
!  variables for many_k simultaneously on device
!
COMPLEX(DP), ALLOCATABLE :: becpk_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: psk_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: pssk_d(:,:,:,:)
REAL(DP),    ALLOCATABLE :: tab_beta_d(:,:,:)

COMPLEX(DP), ALLOCATABLE :: vkbk_d(:,:)
REAL(DP),    ALLOCATABLE :: g2kink_d(:,:)
COMPLEX(DP), ALLOCATABLE :: evck(:,:)
COMPLEX(DP), ALLOCATABLE :: evck_d(:,:)
REAL(DP),    ALLOCATABLE :: h_diagk_d(:,:,:)
REAL(DP),    ALLOCATABLE :: s_diagk_d(:,:,:)
!
!  device variables to copy pw.x variables
!
REAL(DP),    ALLOCATABLE :: tau_d(:,:)
INTEGER,     ALLOCATABLE :: ityp_d(:) 
INTEGER,     ALLOCATABLE :: nh_d(:) 
LOGICAL,     ALLOCATABLE :: type1ps(:)
LOGICAL,     ALLOCATABLE :: type1ps_d(:)
LOGICAL,     ALLOCATABLE :: tvanp(:)
LOGICAL,     ALLOCATABLE :: tvanp_d(:)
INTEGER,     ALLOCATABLE :: nbeta_d(:) 
INTEGER,     ALLOCATABLE :: mill_d(:,:)
INTEGER,     ALLOCATABLE :: ofsbeta_d(:)
INTEGER,     ALLOCATABLE :: nhtolm_d(:,:)
INTEGER,     ALLOCATABLE :: ijtoh_d(:,:,:)
INTEGER,     ALLOCATABLE :: nhtol_d(:,:)
INTEGER,     ALLOCATABLE :: indv_d(:,:)
REAL(DP),    ALLOCATABLE :: g_d(:,:)
REAL(DP),    ALLOCATABLE :: xk_d(:,:)
REAL(DP),    ALLOCATABLE :: wk_d(:)
INTEGER,     ALLOCATABLE :: isk_d(:) 
REAL(DP),    ALLOCATABLE :: qq_at_d(:,:,:)  
COMPLEX(DP), ALLOCATABLE :: qq_so_d(:,:,:,:)
REAL(DP),    ALLOCATABLE :: deeq_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: deeq_nc_d(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: eigts1_d(:,:)
COMPLEX(DP), ALLOCATABLE :: eigts2_d(:,:)
COMPLEX(DP), ALLOCATABLE :: eigts3_d(:,:)
REAL(DP),    ALLOCATABLE :: vrs_d(:,:)

INTEGER :: nat_d
INTEGER :: ntyp_d
INTEGER :: lmaxkb_d
INTEGER :: nhm_d 
INTEGER :: nkb_d 
INTEGER :: nqx_d
INTEGER :: npwx_d
INTEGER :: ngauss_d
REAL(DP) :: tpiba_d
REAL(DP) :: omega_d
REAL(DP) :: dq_d
REAL(DP) :: qcutz_d
REAL(DP) :: q2sigma_d
REAL(DP) :: ecfixed_d
REAL(DP) :: degauss_d
LOGICAL  :: ltetra_d
LOGICAL  :: lgauss_d
LOGICAL  :: okvan_d
LOGICAL :: noncolin_d
LOGICAL :: lspinorb_d
LOGICAL :: lsda_d
LOGICAL :: domag_d

#if defined(__CUDA)
   ATTRIBUTES(DEVICE) :: becpk_d
   ATTRIBUTES(DEVICE) :: psk_d
   ATTRIBUTES(DEVICE) :: pssk_d
   ATTRIBUTES(DEVICE) :: tab_beta_d
   ATTRIBUTES(DEVICE) :: vkbk_d
   ATTRIBUTES(DEVICE) :: g2kink_d
   ATTRIBUTES(DEVICE) :: evck_d
   ATTRIBUTES(DEVICE) :: h_diagk_d
   ATTRIBUTES(DEVICE) :: s_diagk_d

   ATTRIBUTES(DEVICE) :: startkb_d

   ATTRIBUTES(DEVICE) :: nl_d
   ATTRIBUTES(DEVICE) :: work_d
   ATTRIBUTES(DEVICE) :: vrs_d
   ATTRIBUTES(DEVICE) :: wsave_d
   ATTRIBUTES(DEVICE) :: isindex_d
   ATTRIBUTES(DEVICE) :: iplane_d
   ATTRIBUTES(DEVICE) :: lenwrk_d, lensav_d

   ATTRIBUTES(DEVICE) :: tau_d
   ATTRIBUTES(DEVICE) :: ityp_d
   ATTRIBUTES(DEVICE) :: nh_d
   ATTRIBUTES(DEVICE) :: type1ps_d
   ATTRIBUTES(DEVICE) :: tvanp_d
   ATTRIBUTES(DEVICE) :: nbeta_d
   ATTRIBUTES(DEVICE) :: g_d
   ATTRIBUTES(DEVICE) :: mill_d
   ATTRIBUTES(DEVICE) :: ofsbeta_d
   ATTRIBUTES(DEVICE) :: nhtolm_d
   ATTRIBUTES(DEVICE) :: ijtoh_d
   ATTRIBUTES(DEVICE) :: nhtol_d
   ATTRIBUTES(DEVICE) :: indv_d
   ATTRIBUTES(DEVICE) :: xk_d
   ATTRIBUTES(DEVICE) :: wk_d
   ATTRIBUTES(DEVICE) :: isk_d
   ATTRIBUTES(DEVICE) :: qq_at_d
   ATTRIBUTES(DEVICE) :: qq_so_d
   ATTRIBUTES(DEVICE) :: deeq_d
   ATTRIBUTES(DEVICE) :: deeq_nc_d
   ATTRIBUTES(DEVICE) :: eigts1_d
   ATTRIBUTES(DEVICE) :: eigts2_d
   ATTRIBUTES(DEVICE) :: eigts3_d

   ATTRIBUTES(DEVICE) :: nat_d
   ATTRIBUTES(DEVICE) :: ntyp_d
   ATTRIBUTES(DEVICE) :: lmaxkb_d
   ATTRIBUTES(DEVICE) :: nhm_d
   ATTRIBUTES(DEVICE) :: nkb_d
   ATTRIBUTES(DEVICE) :: nqx_d
   ATTRIBUTES(DEVICE) :: npwx_d
   ATTRIBUTES(DEVICE) :: ngauss_d
   ATTRIBUTES(DEVICE) :: tpiba_d
   ATTRIBUTES(DEVICE) :: omega_d
   ATTRIBUTES(DEVICE) :: dq_d
   ATTRIBUTES(DEVICE) :: qcutz_d
   ATTRIBUTES(DEVICE) :: q2sigma_d
   ATTRIBUTES(DEVICE) :: ecfixed_d
   ATTRIBUTES(DEVICE) :: degauss_d
   ATTRIBUTES(DEVICE) :: ltetra_d
   ATTRIBUTES(DEVICE) :: lgauss_d
   ATTRIBUTES(DEVICE) :: okvan_d
   ATTRIBUTES(DEVICE) :: noncolin_d
   ATTRIBUTES(DEVICE) :: lspinorb_d
   ATTRIBUTES(DEVICE) :: lsda_d
   ATTRIBUTES(DEVICE) :: domag_d
#endif
!
! variables for saving k-point dependent variables 
!
PUBLIC becpk_d, psk_d, pssk_d, vkbk_d, g2kink_d, evck_d, evck, h_diagk_d,  &
       s_diagk_d, tab_beta_d
!
! variables for dividing the k points in blocks
!
PUBLIC nkblocks, nksbx, nksb, startkb, startkb_d, memgpu, current_ikb,     &
       alloc_many_k
!
! variables for the fft on device
!
PUBLIC nl_d, work_d, wsave_d, lenwrk, lensav, isindex_d, iplane_d, lenwrk_d, &
       lensav_d 
!
! variable on device to copy pw.x variables
!
PUBLIC nat_d, ntyp_d, ityp_d, nh_d, isk_d, qq_at_d, deeq_d, lmaxkb_d,      &
       nhm_d, tpiba_d, xk_d, wk_d, nbeta_d, dq_d, qq_so_d, deeq_nc_d,      &
       tau_d, nqx_d, npwx_d, type1ps_d, noncolin_d, lspinorb_d, lsda_d,    &
       domag_d, tvanp_d, qcutz_d, ecfixed_d, q2sigma_d, nkb_d, okvan_d,    &
       eigts1_d, eigts2_d, eigts3_d, g_d, mill_d, lgauss_d, ltetra_d,      &
       ngauss_d, degauss_d, omega_d, ofsbeta_d, nhtolm_d, ijtoh_d, vrs_d,  &
       nhtol_d, indv_d
!
! host subroutines of the module
!
PUBLIC init_k_blocks, allocate_many_k, deallocate_many_k,                  &
       initialize_fft_factors, deallocate_fft_factors,                     &
       initialize_device_variables, allocate_becps_many_k,                 & 
       deallocate_becps_many_k
!
! global subroutines of the module
!
#if defined(__CUDA)
PUBLIC prepare_phase_factors
#endif

CONTAINS
!
!------------------------------------------------------------------
SUBROUTINE init_k_blocks(npwx,npol,nbndx,nks,nkstot,nbnd,nnr,nnrs,okvan)
!------------------------------------------------------------------
!
USE uspp,       ONLY : nkb
USE uspp_param, ONLY : nhm
USE lsda_mod,   ONLY : nspin
USE gvect,      ONLY : ngm
USE ions_base,  ONLY : nat, ntyp => nsp
USE noncollin_module, ONLY : lspinorb
!
!   This routine receives the maximum number of G vectors, the number
!   of bands, the number of k points and the maximum available memory
!   and divide the k points in blocks so that cegterg does not run
!   out of memory. It assumes that the memory used by cegterg is
!   in bytes npwx*npol*nbndx*nksbx*16 (*3 in us case, *2 in nc case).
!   memgpu is supposed to be in Gbytes
!
IMPLICIT NONE
INTEGER :: npwx, npol, nbndx, nks, nkstot, nbnd, nnr, nnrs
LOGICAL :: okvan
REAL(DP) :: totmem, fftmem, fixmem

INTEGER :: fact, i, rest

fact=2
IF (okvan) fact=3
!
!  This is an estimate of the memory not due to many_k already allocated,
!  either by cuda fortran or from the acc directives.
!
fixmem= DBLE(nnr) * DBLE(nbnd) * 16 +    & ! psic
        DBLE(nnr) * DBLE(nspin) * 16 +   & ! vrs
        DBLE(npwx) * DBLE(npol) * DBLE (nbnd) * 16 +  & ! evc
        DBLE(npwx) * DBLE(nks) * 8 * 2  +  &  ! igk_k (on acc), igk_k_d (CUDA)
        DBLE(npwx) * DBLE(nkb) * 16 +  &  ! vkb
        DBLE(npwx) * 16  +             &  ! g2kin
        DBLE(3*nbnd*nkstot) * 8 +      &  ! et, wg, btype
        DBLE(8 * ngm) * 8                 ! g, gg, mills, igtongl
IF (npol==2) THEN
   fixmem=fixmem + DBLE(nnr) * DBLE(nbnd) * DBLE(npol) * 16 !psic_nc
   fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(nspin+1) * nat * 16 
                                               ! deeq_nc + qq_at 
   fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(nspin)*DBLE(ntyp) * 16 
                                               ! dvan_so_d
   IF (lspinorb) &
      fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(ntyp) * 16 * 4
                                               ! fcoeff
ELSE
   fixmem=fixmem + DBLE(nhm) * DBLE(nhm) * DBLE(nspin+1) * DBLE(nat) * 8
                                               ! deeq and qq_at
ENDIF
fixmem=fixmem + DBLE(nhm) * DBLE(nhm+1) * DBLE(nspin) * DBLE(nat) * 8 !becsum_d
!
!  This is an estimate of the memory allocated by many_k in any case
!
fixmem=fixmem + DBLE(8 * ngm) * 8     ! g_d, gg_d, mills_d (another CUDA copy)
!
!  This is instead an estimate of the memory allocated by this modulus
!  or by c_bands or by cegterg if all k points where computed 
!  together
!
fftmem=DBLE(nnrs) * DBLE(nbnd) * DBLE(nks) * DBLE(npol) * 16  ! psicmr
totmem=DBLE(npwx*npol*nbndx*nks)*16*fact +   &  ! psi, hpsi, (spsi) 
       DBLE(nbndx*nbndx*nks)*16*3 +          &  ! hc, sc, vc
       DBLE(npwx*npol)*DBLE(nbnd)*DBLE(nks)*16 + &! evck_d
       DBLE(20*nbndx*nks)*8 + & ! work_d, rwork_d, hdiago_d, sdiago_d, ifail_d
       DBLE(nkb*npol)*DBLE(nbnd)*DBLE(nks)*16 * 3 + & ! becpk_d, psk_d, pssk_d
       fftmem 
WRITE(stdout,'(5x,"Estimated memory in cegterg",f10.1," GB")') &
                                                 (totmem+fixmem)/1.D9 
WRITE(stdout,'(5x,"Fixed memory in cegterg",f10.1," GB")')  fixmem/1.D9 
IF (fixmem > (memgpu * 1.D9)) THEN
   WRITE(stdout,'(5x,"memgpu must be greater than ",f8.2," to use many_k")') &
               5.0_DP * fixmem/1.D9
!
!  In this case put one k point per block.
!
   nkblocks=nks
ELSE
   nkblocks=MIN(FLOOR(totmem / (memgpu * 1.D9 - fixmem)) + 1, nks) 
ENDIF

ALLOCATE(nksb(nkblocks))
ALLOCATE(startkb(nkblocks))

rest=MOD(nks,nkblocks)
nksb=nks/nkblocks
DO i=1,rest
   nksb(i)=nksb(i)+1
ENDDO
startkb(1)=0
DO i=2,nkblocks
   startkb(i)=startkb(i-1)+nksb(i-1)
ENDDO

#if defined(__CUDA)
ALLOCATE(startkb_d(nkblocks))
startkb_d=startkb
#endif

nksbx=0
DO i=1,nkblocks
   nksbx=MAX(nksbx,nksb(i))
ENDDO
alloc_many_k=.FALSE.

WRITE(stdout,'(/,5x,"Many-k Davidson diagonalization")')
WRITE(stdout,'(5x,i8," k-points divided into ",i5," blocks")') nks, nkblocks
WRITE(stdout,'(5x,"Maximum size of the blocks",i10)') nksbx

WRITE(stdout,'(5x,"Block number",7x,"size",7x,"start",8x,"end")') 
DO i=1,nkblocks-1
   WRITE(stdout,'(5x,i7,8x,i7,5x,i7,5x,i7)') i, nksb(i), startkb(i)+1, &
                                             startkb(i+1)
ENDDO
WRITE(stdout,'(5x,i7,8x,i7,5x,i7,5x,i7,/)') i, nksb(nkblocks), &
                                             startkb(nkblocks)+1, nks
RETURN
END SUBROUTINE init_k_blocks
!
!------------------------------------------------------------------
SUBROUTINE allocate_many_k(nsolv)
!------------------------------------------------------------------
USE wvfct,            ONLY : nbnd, npwx
USE uspp,             ONLY : nkb
USE lsda_mod,         ONLY : nspin
USE noncollin_module, ONLY : npol, noncolin
USE beta_mod,         ONLY : nqx
USE uspp_param,       ONLY : nhm, nwfcm, nbetam
USE ions_base,        ONLY : nat, ntyp=>nsp
USE gvect,            ONLY : ngm
USE klist,            ONLY : nks
USE fft_base,         ONLY : dfftp

IMPLICIT NONE

INTEGER :: nsolv, nrxx

nrxx=dfftp%nr1x*dfftp%nr2x*dfftp%nr3
ALLOCATE(vkbk_d(npwx,nkb*nksbx))
ALLOCATE(g2kink_d(npwx,nksbx))
ALLOCATE(h_diagk_d(npwx,npol,nksbx))
ALLOCATE(s_diagk_d(npwx,npol,nksbx))

ALLOCATE(tau_d(3,nat))
ALLOCATE(ityp_d(nat))
ALLOCATE(nh_d(ntyp))
ALLOCATE(type1ps(ntyp))
ALLOCATE(type1ps_d(ntyp))
ALLOCATE(tvanp(ntyp))
ALLOCATE(tvanp_d(ntyp))
ALLOCATE(nbeta_d(ntyp))
ALLOCATE(xk_d(3,nks))
ALLOCATE(wk_d(nks))
ALLOCATE(isk_d(nks))
IF (noncolin) THEN
   ALLOCATE(deeq_nc_d(nhm,nhm,nat,nspin))
   ALLOCATE(qq_so_d(nhm,nhm,nspin,ntyp))
ELSE
   ALLOCATE(deeq_d(nhm,nhm,nat,nspin))
   ALLOCATE(qq_at_d(nhm,nhm,nat))
ENDIF
ALLOCATE(mill_d(3,ngm) )
ALLOCATE(ofsbeta_d(nat) )
ALLOCATE(nhtolm_d(nhm,ntyp))
ALLOCATE(ijtoh_d(nhm,nhm,ntyp))
ALLOCATE(nhtol_d(nhm,ntyp))
ALLOCATE(indv_d(nhm,ntyp))
ALLOCATE(g_d(3,ngm) )
ALLOCATE(eigts1_d(-dfftp%nr1:dfftp%nr1,nat) )
ALLOCATE(eigts2_d(-dfftp%nr2:dfftp%nr2,nat) )
ALLOCATE(eigts3_d(-dfftp%nr3:dfftp%nr3,nat) )
ALLOCATE(evck(npwx*npol,nbnd*nksbx*nsolv))
ALLOCATE(evck_d(npwx*npol,nbnd*nksbx*nsolv))
ALLOCATE(tab_beta_d(nqx,nbetam,ntyp))
ALLOCATE(vrs_d(dfftp%nnr,nspin))


RETURN
END SUBROUTINE allocate_many_k

!------------------------------------------------------------------
SUBROUTINE deallocate_many_k()
!------------------------------------------------------------------
IMPLICIT NONE

IF (ALLOCATED(nksb)) DEALLOCATE(nksb)
IF (ALLOCATED(startkb)) DEALLOCATE(startkb)
IF (ALLOCATED(startkb_d)) DEALLOCATE(startkb_d)

IF (ALLOCATED(vkbk_d)) DEALLOCATE(vkbk_d)
IF (ALLOCATED(g2kink_d)) DEALLOCATE(g2kink_d)
IF (ALLOCATED(h_diagk_d)) DEALLOCATE(h_diagk_d)
IF (ALLOCATED(s_diagk_d)) DEALLOCATE(s_diagk_d)

IF (ALLOCATED(tau_d)) DEALLOCATE(tau_d)
IF (ALLOCATED(ityp_d)) DEALLOCATE(ityp_d)
IF (ALLOCATED(nh_d)) DEALLOCATE(nh_d)
IF (ALLOCATED(vrs_d)) DEALLOCATE(vrs_d)
IF (ALLOCATED(type1ps)) DEALLOCATE(type1ps)
IF (ALLOCATED(type1ps_d)) DEALLOCATE(type1ps_d)
IF (ALLOCATED(tvanp)) DEALLOCATE(tvanp)
IF (ALLOCATED(tvanp_d)) DEALLOCATE(tvanp_d)
IF (ALLOCATED(nbeta_d)) DEALLOCATE(nbeta_d)
IF (ALLOCATED(xk_d)) DEALLOCATE(xk_d)
IF (ALLOCATED(wk_d)) DEALLOCATE(wk_d)
IF (ALLOCATED(isk_d)) DEALLOCATE(isk_d)
IF (ALLOCATED(deeq_d)) DEALLOCATE(deeq_d)
IF (ALLOCATED(deeq_nc_d)) DEALLOCATE(deeq_nc_d)
IF (ALLOCATED(qq_at_d)) DEALLOCATE(qq_at_d)
IF (ALLOCATED(qq_so_d)) DEALLOCATE(qq_so_d)
IF (ALLOCATED(evck_d))  DEALLOCATE(evck_d)
IF (ALLOCATED(evck))    DEALLOCATE(evck)
IF (ALLOCATED(tab_beta_d)) DEALLOCATE(tab_beta_d)
IF (ALLOCATED(ofsbeta_d))   DEALLOCATE(ofsbeta_d)
IF (ALLOCATED(nhtolm_d))   DEALLOCATE(nhtolm_d)
IF (ALLOCATED(ijtoh_d))  DEALLOCATE(ijtoh_d)
IF (ALLOCATED(nhtol_d))  DEALLOCATE(nhtol_d)
IF (ALLOCATED(indv_d))   DEALLOCATE(indv_d)
IF (ALLOCATED(g_d))      DEALLOCATE(g_d)
IF (ALLOCATED(mill_d))   DEALLOCATE(mill_d)
IF (ALLOCATED(eigts1_d)) DEALLOCATE(eigts1_d)
IF (ALLOCATED(eigts2_d)) DEALLOCATE(eigts2_d)
IF (ALLOCATED(eigts3_d)) DEALLOCATE(eigts3_d)

RETURN
END SUBROUTINE deallocate_many_k
!
!------------------------------------------------------------------
SUBROUTINE initialize_device_variables()
!------------------------------------------------------------------
USE klist,      ONLY : nks, xk, wk, igk_k, igk_k_d, lgauss, ltetra, &
                       ngauss, degauss
USE ions_base,  ONLY : nat, ntyp=>nsp, ityp, tau
USE cell_base,  ONLY : tpiba, omega
USE uspp,       ONLY : nkb, qq_at, deeq, qq_so, deeq_nc, okvan, ofsbeta, &
                       ijtoh, nhtolm, indv, nhtol
USE gvecw,      ONLY : ecfixed, qcutz, q2sigma
USE gvect,      ONLY : g, mill, eigts1, eigts2, eigts3
USE scf,        ONLY : vrs
USE beta_mod,   ONLY : tab_beta, dq, nqx
USE uspp_param, ONLY : upf, nh, nhm, lmaxkb
USE klist,      ONLY : xk
USE wvfct,      ONLY : npwx
USE lsda_mod,   ONLY : lsda, isk
USE noncollin_module, ONLY : noncolin, lspinorb, domag
IMPLICIT NONE

INTEGER :: nt, nbeta_cpu(ntyp)

#if defined(__CUDA)
igk_k_d=igk_k

tau_d=tau
ityp_d=ityp
nh_d=nh
xk_d(1:3,1:nks)=xk(1:3,1:nks)
wk_d(1:nks)=wk(1:nks)
isk_d(1:nks)=isk(1:nks)
deeq_d=deeq
qq_at_d=qq_at
IF (noncolin) THEN
   deeq_nc_d=deeq_nc
   IF (lspinorb) qq_so_d=qq_so
ENDIF
nat_d=nat
ntyp_d=ntyp
lmaxkb_d=lmaxkb
nhm_d=nhm
nqx_d=nqx
npwx_d=npwx
ngauss_d=ngauss
tpiba_d=tpiba
omega_d=omega
dq_d=dq
qcutz_d=qcutz
q2sigma_d=q2sigma
ecfixed_d=ecfixed
degauss_d=degauss
nkb_d=nkb
ltetra_d=ltetra
lgauss_d=lgauss
okvan_d=okvan
noncolin_d=noncolin
lspinorb_d=lspinorb
lsda_d=lsda
domag_d=domag

vrs_d=vrs
DO nt=1,ntyp
   nbeta_cpu(nt)=upf(nt)%nbeta
   type1ps(nt)=upf(nt)%tvanp .OR.upf(nt)%is_multiproj
   tvanp(nt)=upf(nt)%tvanp
ENDDO
mill_d=mill
g_d=g
eigts1_d=eigts1
eigts2_d=eigts2
eigts3_d=eigts3
nbeta_d=nbeta_cpu
type1ps_d=type1ps
tvanp_d=tvanp
tab_beta_d=tab_beta
ofsbeta_d=ofsbeta
nhtolm_d=nhtolm
ijtoh_d=ijtoh
nhtol_d=nhtol
indv_d=indv
#endif

RETURN
END SUBROUTINE initialize_device_variables
!
!------------------------------------------------------------------
SUBROUTINE allocate_becps_many_k(npe, nsolv)
!------------------------------------------------------------------
USE wvfct,            ONLY : nbnd, npwx
USE uspp,             ONLY : nkb
USE noncollin_module, ONLY : npol

IMPLICIT NONE
INTEGER :: npe, nsolv

ALLOCATE(becpk_d(nkb,npol,nbnd,nksbx*npe*nsolv))
ALLOCATE(psk_d(nkb,npol,nbnd,nksbx*npe*nsolv))
ALLOCATE(pssk_d(nkb,npol,nbnd,nksbx*npe*nsolv))

RETURN
END SUBROUTINE allocate_becps_many_k
!
!------------------------------------------------------------------
SUBROUTINE deallocate_becps_many_k()
!------------------------------------------------------------------
!
IMPLICIT NONE

IF (ALLOCATED(becpk_d)) DEALLOCATE(becpk_d)
IF (ALLOCATED(psk_d))   DEALLOCATE(psk_d)
IF (ALLOCATED(pssk_d))  DEALLOCATE(pssk_d)

RETURN
END SUBROUTINE deallocate_becps_many_k
!
!------------------------------------------------------------------
SUBROUTINE initialize_fft_factors(npe, nsolv)
!------------------------------------------------------------------
USE fft_base,   ONLY : dffts
USE constants,  ONLY : pi
USE wvfct,      ONLY : nbnd
IMPLICIT NONE
INTEGER :: npe, nsolv

INTEGER :: nr1, nr2, nr3, nr1x, nr2x, adim
REAL(DP) :: arg

nr1=dffts%nr1
nr2=dffts%nr2
nr3=dffts%nr3
nr1x=dffts%nr1x
nr2x=dffts%nr2x

ALLOCATE(nl_d(dffts%ngm))
   
#if defined(__CUDA)
nl_d=dffts%nl
adim=MAX(nr1,nr2,nr3)
lenwrk = 2 * adim
lensav = 2 * adim + INT ( LOG ( REAL ( adim, KIND = DP ) ) &
                                      / LOG ( 2.0_DP ) ) + 4
lenwrk_d=lenwrk
lensav_d=lensav
ALLOCATE(work_d(lenwrk,nbnd*nksbx*npe*nsolv*adim))
ALLOCATE(wsave_d(lensav,3))
ALLOCATE(isindex_d(nr1x*nr2x))
ALLOCATE(iplane_d(nr1x))
isindex_d=dffts%isind
iplane_d=dffts%iplw

CALL prepare_phase_factors<<<1,1>>>(wsave_d(1,1), lensav, nr1) 
CALL prepare_phase_factors<<<1,1>>>(wsave_d(1,2), lensav, nr2) 
CALL prepare_phase_factors<<<1,1>>>(wsave_d(1,3), lensav, nr3) 
#endif

RETURN
END SUBROUTINE initialize_fft_factors
!
!------------------------------------------------------------------
SUBROUTINE deallocate_fft_factors()
!------------------------------------------------------------------
!
IMPLICIT NONE

IF (ALLOCATED(nl_d)) DEALLOCATE(nl_d)
IF (ALLOCATED(wsave_d)) DEALLOCATE(wsave_d)
IF (ALLOCATED(work_d)) DEALLOCATE(work_d)
IF (ALLOCATED(isindex_d)) DEALLOCATE(isindex_d)
IF (ALLOCATED(iplane_d)) DEALLOCATE(iplane_d)

RETURN
END SUBROUTINE deallocate_fft_factors
!

#if defined(__CUDA)
!---------------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE prepare_phase_factors(wsave, lensav, n)
!---------------------------------------------------------------------------
IMPLICIT NONE
#include<cfft1i_interf.f90>
INTEGER, VALUE :: lensav, n
REAL(DP), DEVICE, INTENT(INOUT) :: wsave(lensav)

INTEGER :: i
INTEGER :: ier

i=(BlockIdx%x-1) * BlockDim%x + ThreadIdx%x
IF (i > 1) RETURN

wsave(:)=0.0_DP
CALL cfft1i ( n, wsave, lensav, ier )

RETURN
END SUBROUTINE prepare_phase_factors
#endif

END MODULE many_k_mod
