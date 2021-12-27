!
! Copyright (C) 2001 PWSCF group
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE local_dos1d_so (ik, kband, plan)
  !--------------------------------------------------------------------
  !
  !     calculates |psi|^2 for band kband at point ik.
  !     This routine is a generalization of the one contained in PP/src
  !     that calculates also the planar average of the contribution of each
  !     state to the magnetization density in the noncollinear case.
  !     Note that this contribution can be non zero also in a 
  !     non-magnetic system, so it is calculated even if domag is .FALSE.
  !
  USE kinds,     ONLY: dp
  USE cell_base, ONLY: omega
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp
  USE fft_base,  ONLY: dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvecs,   ONLY : doublegrid
  USE lsda_mod, ONLY: nspin, current_spin
  USE uspp, ONLY: becsum, okvan
  USE uspp_param, ONLY: upf, nh, nhm
  USE wvfct, ONLY: npwx, wg
  USE klist, ONLY : ngk, igk_k
  USE noncollin_module, ONLY: noncolin, npol, lspinorb
  USE upf_spinorb, ONLY : fcoef
  USE wavefunctions,  ONLY: evc, psic, psic_nc
  USE becmod, ONLY: bec_type, becp
  USE fft_interfaces,        ONLY : fft_interpolate

  IMPLICIT NONE
  !
  ! input variables
  !
  INTEGER :: ik, kband
  ! input: the k point
  ! input: the band

  REAL(DP) :: plan (dfftp%nr3,nspin)
  ! output: the planar average of this state
  !
  !    Additional local variables for Ultrasoft PP's
  !
  INTEGER :: ikb, jkb, ijkb0, ih, jh, na, ijh, ipol, np, npw
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for ijkb0
  ! counter on solid beta functions
  ! counter on solid beta functions
  ! counter on atoms
  ! counter on composite beta functions
  ! the pseudopotential
  !
  !    And here the local variables
  !
  INTEGER :: ir, ig, is, js, ispin, ibnd, is1, is2, kkb, kh
  ! counter on 3D r points
  ! counter on spin polarizations
  ! counter on g vectors
  ! counter on bands

  REAL(DP) :: w1
  ! the weight of one k point
  REAL(DP), ALLOCATABLE :: aux (:,:)
  ! auxiliary for rho

  COMPLEX(DP), ALLOCATABLE :: prho (:,:)
  ! complex charge for fft
  !
  COMPLEX(DP), ALLOCATABLE :: becsum_nc(:,:,:,:)

  ALLOCATE (prho(dfftp%nnr,nspin))
  ALLOCATE (aux(dfftp%nnr,nspin))

  aux(:,:) = 0.d0
  becsum(:,:,:) = 0.d0

  wg (kband, ik) = 1.d0
  npw=ngk(ik)
  !
  !
  !     First compute the square modulus of the state kband,ik on the smooth
  !     mesh and its contribution to the magnetization density
  !
  IF (noncolin) THEN
     psic_nc = (0.d0,0.d0)
     DO ig = 1, npw
        psic_nc (dffts%nl (igk_k (ig,ik) ), 1 ) = evc (ig     , kband)
        psic_nc (dffts%nl (igk_k (ig,ik) ), 2 ) = evc (ig+npwx, kband)
     ENDDO
     DO ipol=1,npol
        CALL invfft ('Wave', psic_nc(:,ipol), dffts)
     ENDDO

     w1 = wg (kband, ik) / omega
     DO ipol=1,npol
        DO ir = 1, dffts%nnr
           aux(ir,1) = aux(ir,1) + w1 * ( DBLE(psic_nc(ir,ipol))**2 + &
                                         AIMAG(psic_nc(ir,ipol))**2 )
        ENDDO
     ENDDO
!
!   magnetization density
!
     DO ir = 1, dffts%nnr
        aux(ir,2) = aux(ir,2) + w1 * 2.D0 * &
                       (DBLE(psic_nc(ir,1))* DBLE(psic_nc(ir,2)) + &
                       AIMAG(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)))
        aux(ir,3) = aux(ir,3) + w1 * 2.D0 * &
                       (DBLE(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)) - &
                        DBLE(psic_nc(ir,2))*AIMAG(psic_nc(ir,1)))
        aux(ir,4) = aux(ir,4) + w1 * &
                       (DBLE(psic_nc(ir,1))**2+AIMAG(psic_nc(ir,1))**2 &
                       -DBLE(psic_nc(ir,2))**2-AIMAG(psic_nc(ir,2))**2)
     END DO
  ELSE
     psic(1:dffts%nnr) = (0.d0,0.d0)
     DO ig = 1, npw
        psic (dffts%nl (igk_k (ig,ik) ) ) = evc (ig, kband)
     ENDDO
     CALL invfft ('Wave', psic, dffts)

     w1 = wg (kband, ik) / omega
     DO ir = 1, dffts%nnr
        aux(ir,1) = aux(ir,1) + w1 * (DBLE(psic(ir))**2 + AIMAG(psic(ir))**2)
     ENDDO
  ENDIF
  !
  !    If we have a US pseudopotential we compute here the becsum term
  !
  IF (okvan) THEN
     IF (noncolin) THEN
        ALLOCATE(becsum_nc(nhm*(nhm+1)/2,nat,npol,npol))
        becsum_nc=(0.d0, 0.d0)
     END IF
  ENDIF
  !
  ibnd = kband

  w1 = wg (ibnd, ik)
  ijkb0 = 0
  DO np = 1, ntyp
     IF ( upf(np)%tvanp ) THEN
        DO na = 1, nat
           IF (ityp(na)==np) THEN
              ijh = 1
              DO ih = 1, nh(np)
                 ikb = ijkb0 + ih
                 IF (noncolin) THEN
                    DO is=1,npol
                       DO js=1,npol
                          becsum_nc(ijh,na,is,js) =         &
                               becsum_nc(ijh,na,is,js)+w1 *  &
                               CONJG(becp%nc(ikb,is,ibnd)) * &
                                     becp%nc(ikb,js,ibnd)
                       END DO
                    END DO
                 ELSE
                    becsum(ijh,na,current_spin) = &
                          becsum(ijh,na,current_spin) + &
                          w1 * DBLE( CONJG( becp%k(ikb,ibnd) ) * &
                                            becp%k(ikb,ibnd) )
                 END IF
                 ijh = ijh + 1
                 DO jh = ( ih + 1 ), nh(np)
                    jkb = ijkb0 + jh
                    IF (noncolin) THEN
                       DO is=1,npol
                          DO js=1,npol
                             becsum_nc(ijh,na,is,js) =         &
                                 becsum_nc(ijh,na,is,js) + w1 * &
                                 CONJG(becp%nc(ikb,is,ibnd)) *  &
                                       becp%nc(jkb,js,ibnd)
                          END DO
                       END DO
                    ELSE
                       becsum(ijh,na,current_spin) = &
                           becsum(ijh,na,current_spin) + w1 * 2.D0 * &
                                     DBLE( CONJG( becp%k(ikb,ibnd) ) * &
                                                  becp%k(jkb,ibnd) )
                    ENDIF
                    ijh = ijh + 1
                 END DO
              END DO
              ijkb0 = ijkb0 + nh(np)
           END IF
        END DO
     ELSE
        DO na = 1, nat
           IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
        END DO
     END IF
  END DO
  IF (noncolin.AND.okvan) THEN
     DO np = 1, ntyp
        IF ( upf(np)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==np) THEN
                 IF (upf(np)%has_so) THEN
                    CALL transform_becsum_so(becsum_nc,becsum,na)
                 ELSE
                    CALL transform_becsum_nc(becsum_nc,becsum,na)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF
  !
  IF (okvan.AND.noncolin) DEALLOCATE( becsum_nc )
  !
  !    Interpolate on the thick mesh and pass to reciprocal space
  !
  DO ispin=1, nspin
     IF (doublegrid) THEN
        CALL fft_interpolate (dffts, aux(:,ispin), dfftp, aux(:,ispin))
     ENDIF
     DO ir = 1, dfftp%nnr
        prho (ir,ispin) = cmplx(aux (ir,ispin), 0.d0,kind=DP)
     ENDDO
     CALL fwfft ('Rho', prho(:,ispin), dfftp)
  ENDDO
  !
  !    Here we add the US contribution to the charge for the atoms which have
  !    it. Or compute the planar average in the NC case.
  !
  CALL addusdens1d_so (plan, prho)
  !
  DEALLOCATE (aux)
  DEALLOCATE (prho)
  !
  RETURN
END SUBROUTINE local_dos1d_so
