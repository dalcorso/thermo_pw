!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! Copyright (C) 2023 Andrea Dal Corso for the multi-k global version
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  !
  !----------------------------------------------------------------------
#if defined(__CUDA)
  ATTRIBUTES(GLOBAL) SUBROUTINE init_us_2_kernel(vkbk_d, ylm_d, vkb1_d, nhm, &
                              lmaxkb, nkb, npwx, nk, ikt )
    !----------------------------------------------------------------------
    !! Calculates beta functions (Kleinman-Bylander projectors) for nk
    !! k points in parallel on the GPU device. The parallelizzation
    !! is done on nkb and on the k points
    !! structure factor, for all atoms, in reciprocal space.
    !
    USE cudafor
    USE kinds,        ONLY : DP
    USE constants,    ONLY : tpi
    USE many_k_mod,   ONLY : nat=>nat_d, ityp=>ityp_d, tau=>tau_d, &
                             tpiba=>tpiba_d, ntyp=>ntyp_d, nh=> nh_d, &
                             nqx=>nqx_d, dq=>dq_d, &
                             xk=>xk_d, nbeta=>nbeta_d, &
                             eigts1=>eigts1_d, eigts2=>eigts2_d, &
                             eigts3=>eigts3_d, mill=>mill_d, g=>g_d, &
                             tab_at => tab_at_d
    USE klist,        ONLY : igk=>igk_k_d, ngk=>ngk_d
    USE uspp,         ONLY : nhtol=>nhtol_d, nhtolm=>nhtolm_d, indv=>indv_d
    USE uspp,         ONLY : ofsbeta=>ofsbeta_d
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN), VALUE :: npwx, lmaxkb, nhm, nkb, nk
    COMPLEX(DP), INTENT(OUT), DEVICE :: vkbk_d(npwx,nkb*nk)
    INTEGER, DEVICE :: ikt(nk)
    REAL(DP), INTENT(IN), DEVICE :: ylm_d((lmaxkb + 1) **2, npwx, nk)
    REAL(DP), INTENT(INOUT), DEVICE :: vkb1_d(nhm, npwx, nk)
    !
    !     Local variables
    !
    INTEGER :: i0, i1, i2, i3, ig, lm, na, nt, nb, ih, jkb, iv_d, &
               ik, ik1, startk, np, jkb1
    REAL(DP) :: px, ux, vx, wx, arg, q1, q2, q3, rv_d

    COMPLEX(DP) :: phase, pref, sk_d
    INTEGER :: iq
    REAL(DP) :: gk_d(3), qg_d, vq_d
    !
    IF (lmaxkb<0) RETURN
    !
    ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
    IF (ik1>nk) RETURN
    ik=ikt(ik1)
    startk=(ik1-1)*nkb
    np=ngk(ik)
    ig=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
    IF (ig>npwx) RETURN
    IF (ig>np) THEN
       DO jkb=1,nkb
          jkb1=startk+jkb
          vkbk_d(ig, jkb1) = CMPLX(0.0_DP,0.0_DP)
       ENDDO
       RETURN
    ENDIF

    q1 = xk(1,ik)
    q2 = xk(2,ik)
    q3 = xk(3,ik)

    iv_d = igk(ig,ik)
    gk_d (1) = q1 + g(1, iv_d )
    gk_d (2) = q2 + g(2, iv_d )
    gk_d (3) = q3 + g(3, iv_d )
    qg_d = gk_d(1)*gk_d(1) + &
           gk_d(2)*gk_d(2) + &
           gk_d(3)*gk_d(3)
    !
    ! set now qg=|q+G| in atomic units
    !
    qg_d = SQRT(qg_d)*tpiba

    ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
    jkb = 0
    DO nt = 1, ntyp
       DO nb = 1, nbeta(nt)
          rv_d = qg_d
          px = rv_d / dq - int (rv_d / dq)
          ux = 1.d0 - px
          vx = 2.d0 - px
          wx = 3.d0 - px
          i0 = INT( rv_d / dq ) + 1
          i1 = i0 + 1
          i2 = i0 + 2
          i3 = i0 + 3
          vq_d = ux*vx*(wx*tab_at(i0, nb, nt)+px*tab_at(i3, nb, nt)) / 6.d0 + &
              px*wx*(vx * tab_at(i1, nb, nt)-ux * tab_at(i2, nb, nt)) * 0.5d0

          ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
          DO ih = 1, nh (nt)
             IF (nb==indv (ih, nt) ) then
                !l = nhtol (ih, nt)
                lm =nhtolm (ih, nt)
                vkb1_d (ih,ig,ik1) = ylm_d (lm,ig,ik1) * vq_d 
             ENDIF
          ENDDO
       ENDDO
       !
       ! vkb1 contains all betas including angular part for type nt
       ! now add the structure factor and factor (-i)^l
       !
       DO na = 1, nat
          ! ordering: first all betas for atoms of type 1
          !           then  all betas for atoms of type 2  and so on
          IF (ityp (na)== nt) THEN
             arg = (q1 * tau (1, na) + &
                    q2 * tau (2, na) + &
                    q3 * tau (3, na) ) * tpi
             phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
             !
             sk_d = eigts1 (mill(1,igk(ig,ik)), na) * &
                    eigts2 (mill(2,igk(ig,ik)), na) * &
                    eigts3 (mill(3,igk(ig,ik)), na)
             !
             DO ih = 1, nh (nt)
                jkb = jkb + 1
                jkb1= startk + jkb
                pref = (0.d0, -1.d0) **nhtol (ih, nt) * phase
                vkbk_d(ig, jkb1) = vkb1_d (ih,ig,ik1) * sk_d * pref
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !
    RETURN
    END SUBROUTINE init_us_2_kernel
#endif
