!
! Copyright (C) 2024 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__CUDA)
!----------------------------------------------------------------------------
SUBROUTINE  set_hprec_dev( st_d, ikb, nk, g2kink_d, evqk_d, &
                        eprec_d, h_diag_ph_d, nbnd, npe, nsolv, npol, npwx)
!----------------------------------------------------------------------------                
USE cudafor
USE util_param,     ONLY : DP
IMPLICIT NONE
#include<ke_g2kin_interf.f90>
INTEGER :: ikb, nk, nbnd, npe, nsolv, npol, npwx

INTEGER, INTENT(IN), DEVICE :: st_d(nk*npe*nsolv)
REAL(DP), DEVICE, INTENT(INOUT) :: g2kink_d(npwx*npol, nk)
COMPLEX(DP), DEVICE, INTENT(IN) :: evqk_d(npwx*npol, nbnd*nk*nsolv)
REAL(DP), DEVICE, INTENT(INOUT) :: eprec_d(nbnd*nk*nsolv)
REAL(DP), DEVICE, INTENT(INOUT) :: h_diag_ph_d(npwx*npol, nbnd*nk*npe*nsolv)

INTEGER :: ierr

CALL start_clock('ke_g2kin')
CALL ke_g2kin<<<dim3(nk,npwx/32+1,1),dim3(1,32,1)>>>&
                                 ( ikb, nk, g2kink_d, npwx)
ierr=cudaDeviceSynchronize()
CALL stop_clock('ke_g2kin')

CALL start_clock('ke_eprec')
CALL ke_eprec<<<dim3(nk*nsolv,nbnd/32+1,1),dim3(1,32,1)>>>(ikb, nk, g2kink_d, &
                                  evqk_d, eprec_d, nbnd, npol, npwx, nsolv)
ierr=cudaDeviceSynchronize()
CALL start_clock('ke_eprec')

CALL start_clock('ke_hprec')
CALL ke_hprec<<<dim3(nk*npe*nsolv,nbnd,npwx/32+1),dim3(1,1,32)>>>( st_d, ikb, &
                nk, g2kink_d, h_diag_ph_d, eprec_d, nbnd, npe, nsolv, npol, &
                npwx)
ierr=cudaDeviceSynchronize()
CALL stop_clock('ke_hprec')
RETURN
END SUBROUTINE set_hprec_dev

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ke_g2kin( ikb, nk, g2kink_d, npwx)
!-----------------------------------------------------------------------
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d, igk_k => igk_k_d
  USE many_k_mod,     ONLY : qcutz => qcutz_d, ecfixed => ecfixed_d, &
                             q2sigma => q2sigma_d, tpiba => tpiba_d, &
                             xk => xk_d, g => g_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d, ikmks=>ikmks_d,    &
                             ikmkmqs=>ikmkmqs_d, startkb_ph => startkb_ph_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npwx, ikb
  REAL(DP), DEVICE, INTENT(INOUT) :: g2kink_d(npwx, nk)

  INTEGER :: ik1, ik, ikk, ikq, npwq, ig
  REAL(DP) :: xk1, xk2, xk3, tpiba2

  ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (ik1>nk) RETURN
  ik=ik1+startkb_ph(ikb) 

  ikk=ikks(ik)
  ikq=ikqs(ik)
  npwq = ngk(ikq)

  ig=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (ig>npwq) RETURN
  !
  tpiba2=tpiba**2
  xk1 = xk(1,ikq)
  xk2 = xk(2,ikq)
  xk3 = xk(3,ikq)
  !
  g2kink_d(ig,ik1) = ((xk1 + g(1,igk_k(ig,ikq)))*( xk1 + g(1,igk_k(ig,ikq))) + &
                ( xk2 + g(2,igk_k(ig,ikq)) )*( xk2 + g(2,igk_k(ig,ikq)) ) + &
                ( xk3 + g(3,igk_k(ig,ikq)) )*( xk3 + g(3,igk_k(ig,ikq)) ) ) &
                * tpiba2
  !
 IF ( qcutz > 0.D0 ) THEN
     !
     !
     g2kink_d(ig,ik1) = g2kink_d(ig,ik1) + qcutz * &
          ( 1.D0 + erf( ( g2kink_d(ig,ik1) - ecfixed ) / q2sigma ) )
     !
END IF

RETURN
END SUBROUTINE ke_g2kin
!
!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ke_eprec( ikb, nk, g2kink_d, evqk_d, &
                                        eprec_d, nbnd, npol, npwx, nsolv)
!-----------------------------------------------------------------------
  !
  !  This routine computes the array eprec for all k points
  !
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d
  USE many_k_mod,     ONLY : noncolin_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d,    &
                             startkb_ph => startkb_ph_d, &
                             nbnd_occ => nbnd_occ_d
  !
IMPLICIT NONE
!
INTEGER, INTENT(IN), VALUE :: nk, npwx, npol, ikb, nbnd, nsolv
REAL(DP), DEVICE, INTENT(IN) :: g2kink_d(npwx, nk)
COMPLEX(DP), DEVICE, INTENT(IN) :: evqk_d(npwx*npol, nbnd*nk*nsolv)
REAL(DP), DEVICE, INTENT(INOUT) :: eprec_d(nbnd*nk*nsolv)

INTEGER :: id, isolv, ik1, ik, ikk, ikq, npwq, ibnd, ikwf, k2, ig
COMPLEX(DP) :: aux

id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
IF (id>nk*nsolv) RETURN
isolv=(id-1)/nk+1
ik1=MOD(id-1,nk)+1
ik=ik1+startkb_ph(ikb)

ikk=ikks(ik)
ikq=ikqs(ik)
npwq = ngk(ikq)

ibnd=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
IF (ibnd>nbnd_occ(ikk)) RETURN

ikwf=ik1 + (isolv-1)*nk
k2= (ikwf-1)* nbnd + ibnd

aux=(0.d0,0.d0)
DO ig = 1, npwq
   aux = aux + g2kink_d (ig,ik1) * ABS(evqk_d (ig, k2))**2
END DO
IF (noncolin_d) THEN
   DO ig = 1, npwq
      aux = aux + g2kink_d (ig,ik1) * ABS(evqk_d (ig+npwx, k2))**2
   END DO
ENDIF

eprec_d(k2) = 1.35_DP * aux

RETURN
END SUBROUTINE ke_eprec

!-----------------------------------------------------------------------
ATTRIBUTES(GLOBAL) SUBROUTINE ke_hprec( st, ikb, nk, g2kink_d, h_diag_ph_d, &
                                        eprec_d, nbnd, npe, nsolv, npol, npwx)
!-----------------------------------------------------------------------
  USE cudafor
  USE util_param,     ONLY : DP
  USE klist,          ONLY : ngk => ngk_d 
  USE many_k_mod,     ONLY : noncolin_d
  USE many_k_ph_mod,  ONLY : ikks=>ikks_d, ikqs=>ikqs_d,   &
                             startkb_ph => startkb_ph_d,   &
                             nbnd_occ => nbnd_occ_d 
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN), VALUE :: nk, npwx, npol, ikb, nbnd, npe, nsolv
  INTEGER, INTENT(IN), DEVICE :: st(nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(IN) :: g2kink_d(npwx, nk)
  REAL(DP), DEVICE, INTENT(INOUT) :: h_diag_ph_d(npwx*npol, &
                                                     nbnd*nk*npe*nsolv)
  REAL(DP), DEVICE, INTENT(IN) :: eprec_d(nbnd*nk*nsolv)

  INTEGER :: id, isp, isolv, ipert, ik1, st_, ikk, ikq, ik, ibnd, ig, npwq, &
             ikwf, k2

  id=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
  IF (id>nk*npe*nsolv) RETURN
  isp=(id-1)/nk+1
  isolv=(isp-1)/npe+1
  ipert=MOD(isp-1,npe)+1
  ik1=MOD(id-1,nk)+1
  st_=st(id)
  ik=ik1+startkb_ph(ikb)

  ikk=ikks(ik)
  ikq=ikqs(ik)
  npwq=ngk(ikq)

  ibnd=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
  IF (ibnd>nbnd_occ(ikk)) RETURN

  ig=(BlockIdx%z-1)*BlockDim%z + ThreadIdx%z
  IF (ig>npwx) RETURN

  ikwf=ik1 + (isolv-1)*nk

  k2= (ikwf-1)* nbnd + ibnd

  IF (ig>npwq) THEN
     h_diag_ph_d(ig,st_+ibnd)=0.0_DP
  ELSE
     h_diag_ph_d(ig,st_+ibnd)=1.0_DP/MAX(1.0_DP,g2kink_d(ig,ik1)/eprec_d(k2))
  ENDIF
  IF (noncolin_d) THEN
     h_diag_ph_d(ig+npwx,st_+ibnd)=h_diag_ph_d(ig,st_+ibnd)
  ENDIF

RETURN  
END SUBROUTINE ke_hprec
#endif
