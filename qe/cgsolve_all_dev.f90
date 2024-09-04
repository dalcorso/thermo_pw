!
! Copyright (C) 2023 Andrea Dal Corso  
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  This file contains some parts of the conjugate gradient algorithm
!  that can be done in parallel on the GPU.
!
#if defined(__CUDA)
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop2(ndmx, outk, st, conv, &
         lbndk, h, g, h_diag, rho, current_ikb_ph, npol, nk, npe, nsolv, &
         nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_ph_mod, ONLY : ikks => ikks_d, startkb_ph => startkb_ph_d,      &
                           nbnd_occ => nbnd_occ_d

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: g(ndmx*npol,my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: h_diag(ndmx*npol,nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: rho(my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ik, ikk, ipert, id, st_, ibnd, lbnd, ig, isp, isolv

 COMPLEX(DP) :: asum

 ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (ik1>nk) RETURN

 ik=ik1+startkb_ph(current_ikb_ph)
 ikk=ikks(ik)

 isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (isp>npe*nsolv) RETURN
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1

 id=ik1 + (ipert-1)*nk + (isolv-1)*nk*npe
 IF (outk(id)) RETURN
 st_=st(id)

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbnd_occ(ikk)) RETURN
 if (conv (st_+ibnd)==1) RETURN

 lbnd=lbndk(st_+ibnd)

 DO ig=1, ndmx*npol
    h (ig, st_+ibnd)=g (ig, st_+ibnd)
 ENDDO

 DO ig=1, ndmx
    h (ig, st_+ibnd)=h (ig, st_+ibnd) * h_diag (ig, st_+ibnd)
 ENDDO

 IF (npol==2) THEN
    DO ig=1, ndmx
       h (ig+ndmx, st_+ibnd)=h (ig+ndmx, st_+ibnd) * h_diag (ig+ndmx, st_+ibnd)
    ENDDO
 ENDIF

 asum=(0.0_DP,0.0_DP)
 DO ig=1, ndmx*npol
    asum=asum+CONJG(h(ig,st_+ibnd))*g(ig,st_+ibnd)
 ENDDO
 rho(st_+lbnd)=DBLE(asum)

 RETURN

 END SUBROUTINE cgsolve_all_loop2

 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop4(ndmx, outk, st, conv, &
         lbndk, h, hold, dcgammak, eu, e, current_ikb_ph, npol, nk, &
         npe, nsolv, nbnd, my_nbnd, iter, nks)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_ph_mod, ONLY : ikks => ikks_d, startkb_ph => startkb_ph_d,      &
                           nbnd_occ => nbnd_occ_d, ikmks => ikmks_d

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, &
                   current_ikb_ph, nks, iter

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hold(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: dcgammak(my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: eu(my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: e(nbnd,nks)

 INTEGER :: ik1, ik, ikmk, ikk, ipert, id, st_, ibnd, lbnd, ig, isp, isolv

 COMPLEX(DP) :: dcgamma

 ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (ik1>nk) RETURN

 ik=ik1+startkb_ph(current_ikb_ph)
 ikk=ikks(ik)
 ikmk=ikk

 isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (isp>npe*nsolv) RETURN
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1
 IF (isolv==2) ikmk=ikmks(ik)

 id=ik1+(ipert-1)*nk + (isolv-1)*nk*npe
 IF (outk(id)) RETURN
 st_=st(id)

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbnd_occ(ikk)) RETURN
 if (conv (st_+ibnd)==1) RETURN

 lbnd=lbndk(st_+ibnd)

 DO ig=1, ndmx*npol
    h(ig,st_+ibnd)=-1.0_DP*h(ig,st_+ibnd)
 ENDDO

 IF (iter>1) THEN
    dcgamma=dcgammak(st_+ibnd)
    DO ig=1, ndmx*npol
       h (ig, st_+ibnd) = h (ig, st_+ibnd) + dcgamma * hold (ig, st_+ibnd)
    ENDDO
 ENDIF

 DO ig=1, ndmx*npol
    hold (ig, st_+lbnd)= h (ig, st_+ibnd)
 ENDDO

 eu (st_+lbnd) = e (ibnd,ikmk)

 RETURN

 END SUBROUTINE cgsolve_all_loop4
 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop5(ndmx, outk, st, conv, &
         lbndk, g, h, t, a, c, current_ikb_ph, npol, nk, &
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
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 REAL(DP), DEVICE :: a(my_nbnd*nk*npe*nsolv)
 REAL(DP), DEVICE :: c(my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: g(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: t(ndmx*npol,my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ik, ikk, ipert, id, st_, ibnd, lbnd, ig, isp, isolv

 COMPLEX(DP) :: asum

 ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (ik1>nk) RETURN

 ik=ik1+startkb_ph(current_ikb_ph)
 ikk=ikks(ik)

 isp=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (isp>npe*nsolv) RETURN
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1

 id=ik1+(ipert-1)*nk+(isolv-1)*nk*npe
 IF (outk(id)) RETURN
 st_=st(id)

 ibnd=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ibnd>nbnd_occ(ikk)) RETURN
 if (conv (st_+ibnd)==1) RETURN

 lbnd=lbndk(st_+ibnd)

 asum=(0.0_DP,0.0_DP)
 DO ig=1, ndmx*npol
    asum=asum + CONJG(h(ig,st_+ibnd))* g(ig,st_+ibnd)
 ENDDO
 a(st_+lbnd)=DBLE(asum)

 asum=(0.0_DP,0.0_DP)
 DO ig=1, ndmx*npol
    asum=asum + CONJG(h(ig,st_+ibnd))* t(ig,st_+lbnd)
 ENDDO
 c(st_+lbnd)=DBLE(asum)

 RETURN

 END SUBROUTINE cgsolve_all_loop5

 !--------------------------------------------------------------------------
 ATTRIBUTES(GLOBAL) SUBROUTINE cgsolve_all_loop6(ndmx, outk, st, conv, &
         lbndk, dpsi, g, h, hold, t, dclambdak, current_ikb_ph, npol, nk, &
         npe, nsolv, nbnd, my_nbnd)
 !--------------------------------------------------------------------------
 USE cudafor
 USE util_param,    ONLY : DP

 USE many_k_ph_mod, ONLY : ikks => ikks_d, ikqs => ikqs_d,  &
                           startkb_ph => startkb_ph_d,      &
                           nbnd_occ => nbnd_occ_d

 IMPLICIT NONE

 INTEGER, VALUE :: ndmx, nk, npe, nsolv, npol, nbnd, my_nbnd, current_ikb_ph

 LOGICAL, DEVICE :: outk(nk*npe*nsolv)
 INTEGER, DEVICE :: st(nk*npe*nsolv)
 INTEGER, DEVICE :: lbndk(my_nbnd*nk*npe*nsolv)
 INTEGER, DEVICE :: conv(nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: dpsi(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: g(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: h(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: hold(ndmx*npol,my_nbnd*nk*npe*nsolv)
 COMPLEX(DP), DEVICE :: t(ndmx*npol,my_nbnd*nk*npe*nsolv)

 COMPLEX(DP), DEVICE :: dclambdak(my_nbnd*nk*npe*nsolv)

 INTEGER :: ik1, ik, ikk, ikq, ipert, id, st_, ibnd, lbnd, ig, isp, isolv

 COMPLEX(DP) :: dclambda

 id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
 IF (id>nk*npe*nsolv) RETURN
 isp=(id-1)/nk+1
 isolv=(isp-1)/npe+1
 ipert=MOD(isp-1,npe)+1
 ik1=MOD(id-1,nk)+1

 ik=ik1+startkb_ph(current_ikb_ph)
 ikk=ikks(ik)
 ikq=ikqs(ik)

 IF (outk(id)) RETURN
 st_=st(id)

 ibnd=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
 IF (ibnd>nbnd_occ(ikk)) RETURN
 IF (conv (st_+ibnd)==1) RETURN

 ig=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
 IF (ig > ndmx * npol) RETURN

 lbnd=lbndk(st_+ibnd)
 dclambda=dclambdak(st_+lbnd)

! DO ig=1, ndmx*npol
    dpsi(ig,st_+ibnd)=dpsi(ig,st_+ibnd) + dclambda * h(ig,st_+ibnd)
! ENDDO

! DO ig=1, ndmx*npol
    g(ig,st_+ibnd) = g(ig,st_+ibnd) + dclambda * t(ig,st_+lbnd)
! ENDDO

! DO ig=1,ndmx*npol
    hold(ig,st_+ibnd) = h(ig,st_+ibnd) 
! ENDDO

 RETURN
 END SUBROUTINE cgsolve_all_loop6
#endif
