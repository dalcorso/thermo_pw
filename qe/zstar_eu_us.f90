!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
SUBROUTINE zstar_eu_us_tpw(dvscfin)
!----------===========-----------------------------------------
  !
  ! Calculate the part of the Born effective charges 
  ! due to the integral of the change of the exchange and correlation 
  ! potential and the change of the charge due to the orthogonality
  ! constraint and to the change of the augmentation functions.
  !
  ! The zstar calculated in the PAW case has an additional term 
  ! which is added here.
  ! The output of this routine is saved in zstareu0_rec
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp=>nsp, ityp
  USE cell_base,  ONLY : omega
  USE fft_base,   ONLY : dfftp
  USE uspp_param, ONLY :  upf, nh
  USE uspp,       ONLY : okvan
  USE paw_variables,    ONLY : okpaw
  USE noncollin_module, ONLY : nspin_mag, noncolin
  USE buffers,    ONLY : get_buffer
  USE phus,       ONLY : becsumort
  USE lrus,       ONLY : int3_paw
  USE zstar_add,  ONLY : zstareu0_rec
  USE modes,      ONLY : npert, nirr
  USE partial,    ONLY : done_irr, comp_irr
  USE units_ph,   ONLY : lrdrhous, iudrhous
  USE io_global,  ONLY : stdout

  USE mp_pools,             ONLY : inter_pool_comm, root_pool, my_pool_id
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER :: ih, jh, ijh, na, nt, is
  INTEGER :: jpol, imode0, mode0, irr, ipert, npe, imode, mode

  COMPLEX(DP), INTENT(IN) :: dvscfin( dfftp%nnr , nspin_mag, 3)
  COMPLEX(DP), ALLOCATABLE :: drhoscfh (:,:)
  COMPLEX(DP):: zstareu0_wrk (3,3*nat)
  COMPLEX(DP), EXTERNAL :: zdotc

  INTEGER :: icar, jcar

  IF (.NOT.okvan) RETURN

  CALL start_clock('zstar_eu_us')
!
!
! Calculate the integral of the perturbed Hartree and exchange and correlation 
! potential with the change of the charge density due to the change of the
! orthogonality constraint and to the change of the augmentation functions.  
!
!
  imode0 = 0
  ALLOCATE(drhoscfh(dfftp%nnr,nspin_mag)) 
  zstareu0_wrk=(0.0_DP, 0.0_DP)

  DO irr = 1, nirr
     npe = npert(irr)
     DO imode = 1, npe
        mode = imode0 + imode
!        IF (my_pool_id==0) &
           CALL get_buffer (drhoscfh, lrdrhous, iudrhous, mode)
!        CALL mp_bcast(drhoscfh, root_pool, inter_pool_comm)
        DO jpol = 1, 3
           zstareu0_wrk(jpol,mode) =  zstareu0_wrk(jpol,mode) -             &
               zdotc(dfftp%nnr*nspin_mag,dvscfin(1,1,jpol),1,drhoscfh,1) &
               * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3) 
        ENDDO
     ENDDO
     imode0 = imode0 + npe
  ENDDO
!
!  Collect the contribution of all g vectors
!
  CALL mp_sum ( zstareu0_wrk, intra_bgrp_comm )
  CALL mp_sum ( zstareu0_wrk, inter_pool_comm )
  zstareu0_rec= zstareu0_rec + zstareu0_wrk
!
!  Compute the PAW contribution. This term is present only in the PAW case.
!
  IF (okpaw) THEN
     mode0 = 0
     DO irr = 1, nirr
        DO ipert = 1, npert (irr)
           mode = mode0 + ipert
           DO nt=1,ntyp
              IF (upf(nt)%tpawp) THEN
                 ijh=0
                 DO ih=1,nh(nt)
                    DO jh=ih,nh(nt)
                       ijh=ijh+1
                       DO na=1,nat
                          IF (ityp(na)==nt) THEN
                             DO jpol = 1, 3
                                DO is=1,nspin_mag
                                 zstareu0_rec(jpol,mode)=&
                                            zstareu0_rec(jpol,mode)  &
                                    -int3_paw(ih,jh,na,is,jpol)* &
                                            becsumort(ijh,na,is,mode)
                                ENDDO
                             ENDDO
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
        mode0 = mode0 + npert (irr)
     ENDDO
  ENDIF

  DEALLOCATE ( drhoscfh )
  CALL stop_clock('zstar_eu_us')

  RETURN
END SUBROUTINE zstar_eu_us_tpw
