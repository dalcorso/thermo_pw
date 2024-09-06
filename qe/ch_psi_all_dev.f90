!
! Copyright (C) 2023 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
!
!----------------------------------------------------------------------
SUBROUTINE ch_psi_dev (ndmx, outk_d, kdimk_d, npwk_d, st_d, nbndk_d,  &
              dpsi, ah, hpsi, spsi, ps, eu, current_ikb_ph,           &
              npol, nk, npe, nsolv, nkb, nbnd, my_nbnd, alpha_pv)
!----------------------------------------------------------------------
!
!   This routine computes on GPU several terms that were computed by 
!   ch_psi.
!
!
USE cudafor
USE kinds,         ONLY : DP
USE many_k_ph_mod, ONLY : evqk_d
USE mp_bands,      ONLY : intra_bgrp_comm
USE uspp,          ONLY : okvan
USE mp,            ONLY : mp_sum
IMPLICIT NONE
#include<ch_psi_all_interf.f90>

INTEGER  :: ndmx, current_ikb_ph, npol, nk, npe, nsolv, nkb, nbnd, my_nbnd
REAL(DP) :: alpha_pv
LOGICAL, DEVICE :: outk_d(nk*npe*nsolv)
INTEGER, DEVICE :: npwk_d(nk*npe*nsolv)
INTEGER, DEVICE :: kdimk_d(nk*npe*nsolv)
INTEGER, DEVICE :: st_d(nk*npe*nsolv) 
INTEGER, DEVICE :: nbndk_d(nk*npe*nsolv)

REAL(DP), DEVICE :: eu(my_nbnd*nk*npe*nsolv)
COMPLEX(DP), DEVICE :: dpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
COMPLEX(DP), DEVICE :: ah(ndmx*npol,my_nbnd*nk*npe*nsolv)
COMPLEX(DP), DEVICE :: hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
COMPLEX(DP), DEVICE :: spsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
COMPLEX(DP), DEVICE :: ps(nbnd,my_nbnd*nk*npe*nsolv)

INTEGER :: ierr

!!$acc kernels present(ps)
!ps (:,:) = (0.d0, 0.d0)
!!$acc end kernels
CALL start_clock('chp_ps')
CALL ch_psi_computeps<<<dim3(nk*npe*nsolv,nbnd,nbnd/4+1),dim3(1,1,4)>>>(ndmx, &
     outk_d, kdimk_d, st_d, nbndk_d, evqk_d, spsi, ps, current_ikb_ph, &
     npol, nk, npe, nsolv, nbnd, my_nbnd, alpha_pv)
ierr=cudaDeviceSynchronize()
!CALL mp_sum ( ps, intra_bgrp_comm )
CALL stop_clock('chp_ps')
CALL start_clock('chp_ah')
CALL ch_psi_ah<<<dim3(nk,npe*nsolv,nbnd),dim3(1,1,1)>>>(ndmx, &
              outk_d, st_d, nbndk_d, ah, hpsi, spsi, eu, current_ikb_ph, &
              npol, nk, npe, nsolv, nbnd, my_nbnd)
ierr=cudaDeviceSynchronize()
CALL stop_clock('chp_ah')
CALL start_clock('chp_lo2')
CALL ch_psi_lo2<<<dim3(nk,npe*nsolv,nbnd/4+1),dim3(1,1,4)>>>&
            (ndmx, outk_d, st_d, nbndk_d, evqk_d, hpsi, ps, current_ikb_ph, &
             npol, nk, npe, nsolv, nbnd, my_nbnd)
ierr=cudaDeviceSynchronize()
CALL stop_clock('chp_lo2')
CALL start_clock('chp_calbec')
CALL ch_psi_calbec<<<dim3(nk*npe*nsolv,nkb,nbnd/4+1),dim3(1,1,4)>>>(ndmx,  &
    outk_d, st_d, nbndk_d, npwk_d, hpsi, current_ikb_ph, npol, nk, &
    npe, nsolv, nkb, nbnd, my_nbnd)
ierr=cudaDeviceSynchronize()
CALL stop_clock('chp_calbec')
IF (okvan) THEN
   CALL start_clock('chp_qqps')
   CALL ch_psi_qqps<<<dim3(nk*npe*nsolv,nbnd,1),dim3(1,1,1)>>>(outk_d, &
         nbndk_d, nkb, nk, npe, nsolv, npol )
   ierr=cudaDeviceSynchronize()
   CALL stop_clock('chp_qqps')
   CALL start_clock('chp_sqdpsi')
   CALL ch_psi_sqdpsi<<<dim3(nk*npe*nsolv,nbnd,ndmx/32+1),dim3(1,1,32)>>>(ndmx,&
        outk_d, npwk_d, st_d, nbndk_d, hpsi, ah, current_ikb_ph, npol, nk, &
        npe, nsolv, nkb, nbnd, my_nbnd)
   ierr=cudaDeviceSynchronize()
   CALL stop_clock('chp_sqdpsi')
ELSE
   CALL ch_psi_sdpsi<<<dim3(nk*npe*nsolv,ndmx,nbnd),dim3(1,1,1)>>>(ndmx, &
        outk_d, st_d, nbndk_d, hpsi, ah, current_ikb_ph, npol, nk, &
        npe, nsolv, nkb, nbnd, my_nbnd)
   ierr=cudaDeviceSynchronize()
ENDIF

RETURN
END SUBROUTINE ch_psi_dev

 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_computeps(ndmx, outk, kdimk, st, &
         nbndk, evqk, spsi, ps, current_ikb_ph, npol, nk, &
         npe, nsolv, nbnd, my_nbnd, alpha_pv)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_ph_mod, ONLY : ikks => ikks_d, ikqs => ikqs_d, &
                           startkb_ph => startkb_ph_d,     &
                           nbnd_occ => nbnd_occ_d
 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph
 REAL(DP), VALUE :: alpha_pv 

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)
 INTEGER, DEVICE :: kdimk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: evqk(ndmx*npol,my_nbnd*nk*nsolv)
 COMPLEX(DP), DEVICE :: spsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: ps(nbnd,my_nbnd*nk*npe*nsolv)


 INTEGER :: ik1, ik, ikk, ikq, isolv, ipert, isp, id, st_, ibnd, k, ig, &
            kstart, kdim

 COMPLEX(DP) :: asum

 id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (id>nk*npe*nsolv) RETURN
 IF (outk(id)) RETURN
 isp=(id-1)/nk+1
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1
 ik1=MOD(id-1,nk)+1
 kdim=kdimk(id)

 kstart=nbnd*(ik1-1) + (isolv - 1) * nk * nbnd
 ik=ik1+startkb_ph(current_ikb_ph)
 ikk=ikks(ik)
 ikq=ikqs(ik)

 st_=st(id)

 k=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (k>nbnd_occ(ikq)) RETURN
 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbndk(id)) RETURN

 asum=(0.0_DP,0.0_DP)
 DO ig=1, kdim
    asum=asum+CONJG(evqk(ig,kstart+k))*spsi(ig,st_+ibnd)
 ENDDO
 ps(k,st_+ibnd)=asum * alpha_pv

 RETURN

 END SUBROUTINE ch_psi_computeps
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_ah(ndmx, outk, st, &
         nbndk, ah, hpsi, spsi, eu, current_ikb_ph, npol, nk, &
         npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_ph_mod, ONLY : ikks => ikks_d, startkb_ph => startkb_ph_d,      &
                           nbnd_occ => nbnd_occ_d

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 REAL(DP), DEVICE :: eu(my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: ah(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: spsi(ndmx*npol,my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ik, ikk, isolv, ipert, isp, id, st_, ibnd, ig

 ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (ik1>nk) RETURN

 ik=ik1+startkb_ph(current_ikb_ph)
 ikk=ikks(ik)

 isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (isp>npe*nsolv) RETURN
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1

 id=ik1+(ipert-1)*nk + (isolv-1)*npe*nk
 IF (outk(id)) RETURN
 st_=st(id)

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbndk(id)) RETURN

 DO ig=1, ndmx*npol
    ah (ig, st_+ibnd)=hpsi(ig,st_+ibnd)-eu(st_+ibnd)*spsi(ig,st_+ibnd)
 ENDDO

 RETURN

 END SUBROUTINE ch_psi_ah
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_lo2(ndmx, outk, st, &
         nbndk, evqk, hpsi, ps, current_ikb_ph, npol, nk,      &
         npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_ph_mod, ONLY : ikqs => ikqs_d,                  &
                           startkb_ph => startkb_ph_d,      &
                           nbnd_occ => nbnd_occ_d

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: evqk(ndmx*npol,my_nbnd*nk*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: ps(nbnd,my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ik, ipert, id, st_, ibnd, ig, kstart, nbd, ikq, k, &
            isp, isolv
 COMPLEX(DP) :: asum

 ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (ik1>nk) RETURN
 isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (isp>npe*nsolv) RETURN
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1
 id=ik1 + (ipert-1)*nk + (isolv -1)*nk*npe
 IF (outk(id)) RETURN

 kstart=nbnd*(ik1-1) + (isolv-1) * nk * nbnd
 ik=ik1+startkb_ph(current_ikb_ph)
 ikq=ikqs(ik)

 st_=st(id)

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbndk(id)) RETURN

 nbd=nbnd_occ(ikq)

 DO ig=1,ndmx*npol
    asum=(0.0_DP, 0.0_DP)
    DO k=1, nbd
       asum=asum + evqk(ig,kstart+k) * ps(k, st_ + ibnd)
    ENDDO
    hpsi(ig,st_+ibnd)= asum
 ENDDO

 RETURN

 END SUBROUTINE ch_psi_lo2
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_calbec(ndmx, outk, st,  &
         nbndk, npw, hpsi, current_ikb_ph, npol, nk,   &
         npe, nsolv, nkb, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_mod,    ONLY : becpk_d, vkbk_d

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, nkb, &
                   current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: npw(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: hpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ik, ikk, ipert, id, st_, ibnd, ig, kstart, nbd, &
            k, j, n, i0, isp, isolv
 COMPLEX(DP) :: asum

 id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (id>nk*npe*nsolv) RETURN
 IF (outk(id)) RETURN
 isp=(id-1)/nk+1
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1
 ik1=MOD(id-1,nk)+1

 kstart=nkb*(ik1-1)
 n=npw(id)

 st_=st(id)

 i0=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (i0>nkb) RETURN

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbndk(id)) RETURN

 asum=(0.0_DP,0.0_DP)
 DO j=1, n
    asum=asum+CONJG(vkbk_d(j,kstart+i0))*hpsi(j,st_+ibnd)
 ENDDO
 becpk_d(i0,1,ibnd,id)=asum

 IF (npol==2) THEN
    asum=(0.0_DP,0.0_DP)
    DO j=1, n
       asum=asum+CONJG(vkbk_d(j,kstart+i0))*hpsi(ndmx+j,st_+ibnd)
    ENDDO
    becpk_d(i0,2,ibnd,id)=asum
 ENDIF

 RETURN

 END SUBROUTINE ch_psi_calbec
 !
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_sqdpsi(lda, outk, npw, st,  &
         nbndk, hpsi, ah, current_ikb_ph, npol, nk,   &
         npe, nsolv, nkb, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_mod,    ONLY : vkbk_d, pssk_d

 IMPLICIT NONE

 INTEGER, VALUE :: lda, nk, npe, nsolv, npol, nbnd, my_nbnd, nkb, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: npw(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: ah(lda*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(lda*npol,my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ipert, id, st_, ibnd, ig, kstart, nbd, ikq, n, &
            k, j, j1, i0, isp, isolv
 COMPLEX(DP) :: asum

 id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (id>nk*npe*nsolv) RETURN
 IF (outk(id)) RETURN
 isp=(id-1)/nk+1
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1
 ik1=MOD(id-1,nk)+1
 n=npw(id)

 kstart=nkb*(ik1-1)
 st_=st(id)


 ibnd=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (ibnd>nbndk(id)) RETURN
 k=st_+ibnd

 ig=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ig>n) RETURN

 asum=(0.0_DP,0.0_DP)
 DO j=1, nkb
    j1=j+kstart
    asum=asum+vkbk_d(ig,j1)*pssk_d(j,1,ibnd,id)
 ENDDO
 ah(ig,k) = ah(ig,k) + hpsi(ig,k) + asum

 IF (npol==2) THEN
    asum=(0.0_DP,0.0_DP)
    DO j=1, nkb
       j1=j+kstart
       asum=asum+vkbk_d(ig,j1)*pssk_d(j,2,ibnd,id)
    ENDDO
    ah(lda+ig,k) = ah(lda+ig,k) +hpsi(lda+ig,k) + asum
 ENDIF

 RETURN

 END SUBROUTINE ch_psi_sqdpsi
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_sdpsi(lda, outk, st,     &
         nbndk, hpsi, ah, current_ikb_ph, npol, nk, npe, nsolv, &
         nkb, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 IMPLICIT NONE

 INTEGER, VALUE :: lda, nk, npe, nsolv, npol, nbnd, my_nbnd, nkb, &
                   current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: nbndk(nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: ah(lda*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hpsi(lda*npol,my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ik, ikk, ipert, id, st_, ibnd, ig, nbd, ikq, &
            k, j, j1, i0, isolv, isp

 id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (id>nk*npe*nsolv) RETURN
 IF (outk(id)) RETURN
 isp=(id-1)/nk+1
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1
 ik1=MOD(id-1,nk)+1

 st_=st(id)

 ig=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (ig>lda) RETURN

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbndk(id)) RETURN
 k=st_+ibnd

 ah(ig,k)=ah(ig,k) + hpsi(ig,k) 

 IF (npol==2) THEN
    ah(lda+ig,k)=ah(lda+ig,k) + hpsi(lda+ig,k) 
 ENDIF

 RETURN

 END SUBROUTINE ch_psi_sdpsi

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ch_psi_qqps( outk, nbndk, nkb, nk, npe, &
                                            nsolv, npol )
  !-----------------------------------------------------------------------
  !
  ! This routine computes becp for all k points.
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE many_k_mod,     ONLY : pssk_d, becpk_d, isk_d, nat=>nat_d, &
                             ntyp=>ntyp_d, ityp=>ityp_d, nh=>nh_d,      &
                             qq_at=>qq_at_d, lspinorb => lspinorb_d,    &
                             tvanp => tvanp_d, qq_so => qq_so_d

  USE uspp,           ONLY : ofsbeta=>ofsbeta_d
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, nkb, npe, nsolv, npol
  !! input: the number of k points 
  !! input: the number of nonlocal projectors
  !! input: the number of perturbations
  !! input: the number of linear systems to solve
  !! input: the number of components of the wavefunctions
  LOGICAL, INTENT(IN), DEVICE :: outk(nk*npe*nsolv)
  INTEGER, INTENT(IN), DEVICE :: nbndk(nk*npe*nsolv)
  !
  !  ... local variables
  !
  INTEGER :: id, ibnd, m, k, j, nt, na, nhnt, ipol
  !
  COMPLEX(DP) :: asum
  !
  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe*nsolv) RETURN
  IF (outk(id)) RETURN
  IF (nkb == 0) RETURN

  m=nbndk(id) 

  ibnd=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (ibnd>m) RETURN

  DO ipol=1, npol
     pssk_d(:,ipol,ibnd,id)=(0.0_DP,0.0_DP)
  ENDDO

  IF (lspinorb) THEN

     DO nt = 1, ntyp
        !
        IF ( nh(nt) == 0 ) CYCLE
        !
        IF (.NOT.tvanp(nt)) CYCLE
        !
        nhnt = nh(nt)
        !
        DO na = 1, nat
           !
           IF ( ityp(na) == nt ) THEN
              !
              DO k=1, nhnt
                 !
                 asum=(0.0_DP,0.0_DP)
                 DO j=1, nhnt
                    asum=asum+qq_so(k,j,1,nt)*                          &
                                     becpk_d(ofsbeta(na)+j,1,ibnd,id) + &
                              qq_so(k,j,2,nt)*                          &
                                     becpk_d(ofsbeta(na)+j,2,ibnd,id) 
                 ENDDO
                 pssk_d(ofsbeta(na)+k,1,ibnd,id)=&
                                  pssk_d(ofsbeta(na)+k,1,ibnd,id)+asum

                 asum=(0.0_DP,0.0_DP)
                 DO j=1, nhnt
                    asum=asum+qq_so(k,j,3,nt)*                          &
                                     becpk_d(ofsbeta(na)+j,1,ibnd,id) + &
                              qq_so(k,j,4,nt)*                          &
                                     becpk_d(ofsbeta(na)+j,2,ibnd,id) 
                 ENDDO
                 pssk_d(ofsbeta(na)+k,2,ibnd,id)=&
                                  pssk_d(ofsbeta(na)+k,2,ibnd,id)+asum
                 !
              ENDDO
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
  ELSE
     DO nt = 1, ntyp
        !
        IF ( nh(nt) == 0 ) CYCLE
        !
        IF (.NOT.tvanp(nt)) CYCLE
        !
        nhnt = nh(nt)
        !
        DO na = 1, nat
           !
           IF ( ityp(na) == nt ) THEN
              !
              DO k=1, nhnt
                 !
                 DO ipol=1, npol
                    asum=(0.0_DP,0.0_DP)
                    DO j=1, nhnt
                       asum=asum+qq_at(k,j,na)*&
                                        becpk_d(ofsbeta(na)+j,ipol,ibnd,id)
                    ENDDO
                    pssk_d(ofsbeta(na)+k,ipol,ibnd,id)=&
                                     pssk_d(ofsbeta(na)+k,ipol,ibnd,id)+asum
                 ENDDO
                 !
              ENDDO
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO
     !
  ENDIF
     !
  RETURN
  !
END SUBROUTINE ch_psi_qqps
#endif
!
