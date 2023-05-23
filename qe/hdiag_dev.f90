! Copyright (C) 2023 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  !
  !----------------------------------------------------------------------
#if defined(__CUDA)
  ATTRIBUTES(GLOBAL) SUBROUTINE hdiag_kernel(g2kink, h_diagk, s_diagk, npwx, &
                                          v_of_0, npol, nkb, nk, current_ikb )
    !----------------------------------------------------------------------
    !
    ! Calculate an estimate of the diagonal elements of the Hamiltonian
    ! and of the overlap matrix S to be used in preconditioning.
    !
    USE cudafor
    USE kinds,        ONLY : DP
    USE constants,    ONLY : tpi 
    USE many_k_mod,   ONLY : tpiba=>tpiba_d, startkb=>startkb_d,      &
                             nat=>nat_d, ityp=>ityp_d, tau=>tau_d,    &
                             tpiba=>tpiba_d, ntyp=>ntyp_d, nh=> nh_d, &
                             vkbk_d, qq_at=>qq_at_d, deeq=>deeq_d,      &
                             type1ps=>type1ps_d, isk=>isk_d, xk=> xk_d, &
                             lsda => lsda_d, noncolin => noncolin_d,    &
                             lspinorb => lspinorb_d, deeq_nc => deeq_nc_d, &
                             qq_so => qq_so_d, ecfixed => ecfixed_d, &
                             q2sigma => q2sigma_d, qcutz => qcutz_d
    USE gvect,        ONLY : g=>g_d
    USE klist,        ONLY : igk=>igk_k_d, ngk=>ngk_d
    USE uspp,         ONLY : ofsbeta=>ofsbeta_d
    !
    IMPLICIT NONE
    !
    INTEGER, VALUE, INTENT(IN) :: npwx, nkb, nk, npol, current_ikb
    !! input: maximum number of plane waves
    !! input: maximum number of projectors
    !! input: number of k points
    !! input: number of components of the wavefunctions
    !! input: current block of k points
    REAL(DP), VALUE, INTENT(IN) :: v_of_0
    !! input the local potential at G=0.
    REAL(DP), DEVICE, INTENT(OUT) :: g2kink(npwx,nk)
    !! output: the kinetic energy
    REAL(DP), DEVICE, INTENT(OUT) :: h_diagk(npwx,npol,nk)
    !! output: the diagonal elements of the Hamiltonian
    REAL(DP), DEVICE, INTENT(OUT) :: s_diagk(npwx,npol,nk)
    !! output: the diagonal elements of the overlap matrix
    !
    !     Local variables
    !
    INTEGER :: ig, ik, ik1, iv_d, np, nh_, na, nt, ijk_start, ih, jh, &
               ikb, jkb, startk, ijkb_start, current_spin, ipol
    COMPLEX(DP) :: ar, cv
    REAL(DP) :: q1, q2, q3, gk_d(3), qg_d, sum_h, sum_s, sum_h1, sum_h4, &
                sum_s1, sum_s4
    !
    ik1=(BlockIdx%x-1)*BlockDim%x + ThreadIdx%x
    IF (ik1>nk) RETURN
    ik=ik1+startkb(current_ikb)
    startk=(ik1-1)*nkb
    current_spin=1
    IF (lsda) current_spin=isk(ik)
    np=ngk(ik)
    ig=(BlockIdx%y-1)*BlockDim%y + ThreadIdx%y
    IF (ig>np) RETURN

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
    qg_d = qg_d*tpiba*tpiba

    g2kink(ig,ik1)=qg_d
    DO ipol=1, npol
       h_diagk (ig,ipol,ik1)=qg_d + v_of_0
       s_diagk (ig,ipol,ik1)=1.0_DP
    ENDDO

    IF (lspinorb) THEN
       DO nt = 1, ntyp
          IF ( type1ps(nt) ) THEN
             DO na = 1, nat
                IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   sum_h1 = 0.d0
                   sum_h4 = 0.d0
                   sum_s1 = 0.d0
                   sum_s4 = 0.d0
                   DO ih = 1, nh_
                      DO jh = 1, nh_
                         ikb = ijkb_start + ih
                         cv = vkbk_d(ig,ikb+startk)
                         jkb = ijkb_start + jh
                         ar = cv*CONJG(vkbk_d(ig,jkb+startk))
                         sum_h1 = sum_h1 + DBLE(deeq_nc(ih,jh,na,1) * ar)
                         sum_h4 = sum_h4 + DBLE(deeq_nc(ih,jh,na,4) * ar)
                         sum_s1 = sum_s1 + DBLE(qq_so(ih,jh,1,nt) * ar)
                         sum_s4 = sum_s4 + DBLE(qq_so(ih,jh,4,nt) * ar)
                      ENDDO
                   ENDDO
                   h_diagk (ig,1,ik1) = h_diagk (ig,1,ik1) + sum_h1
                   h_diagk (ig,2,ik1) = h_diagk (ig,2,ik1) + sum_h4
                   s_diagk (ig,1,ik1) = s_diagk (ig,1,ik1) + sum_s1
                   s_diagk (ig,2,ik1) = s_diagk (ig,2,ik1) + sum_s4
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   sum_h1 = 0.d0
                   sum_h4 = 0.d0
                   sum_s1 = 0.d0
                   sum_s4 = 0.d0
                   DO ih = 1, nh_
                      ikb = ijkb_start + ih
                      ar = vkbk_d(ig,startk+ikb)*CONJG(vkbk_d(ig,startk+ikb))
                      sum_h1 = sum_h1 + DBLE(deeq_nc(ih,ih,na,1) * ar)
                      sum_h4 = sum_h4 + DBLE(deeq_nc(ih,ih,na,4) * ar)
                      sum_s1 = sum_s1 + DBLE(qq_so(ih,ih,1,nt) * ar)
                      sum_s4 = sum_s4 + DBLE(qq_so(ih,ih,4,nt) * ar)
                   END DO
                   h_diagk (ig,1,ik1) = h_diagk (ig,1,ik1) + sum_h1
                   h_diagk (ig,2,ik1) = h_diagk (ig,2,ik1) + sum_h4
                   s_diagk (ig,1,ik1) = s_diagk (ig,1,ik1) + sum_s1
                   s_diagk (ig,2,ik1) = s_diagk (ig,2,ik1) + sum_s4
                END IF
             END DO
          END IF
       END DO
    ELSE IF (noncolin) THEN
       DO nt = 1, ntyp
          IF ( type1ps(nt) ) THEN
             DO na = 1, nat
                IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   sum_h1 = 0.d0
                   sum_h4 = 0.d0
                   sum_s = 0.d0
                   DO ih = 1, nh_
                      DO jh = 1, nh_
                         ikb = ijkb_start + ih
                         cv = vkbk_d(ig,ikb+startk)
                         jkb = ijkb_start + jh
                         ar = cv*CONJG(vkbk_d(ig,jkb+startk))
                         sum_h1 = sum_h1 + DBLE(deeq_nc(ih,jh,na,1) * ar)
                         sum_h4 = sum_h4 + DBLE(deeq_nc(ih,jh,na,4) * ar)
                         sum_s = sum_s + DBLE(qq_at(ih,jh,na) * ar)
                      ENDDO
                   ENDDO
                   h_diagk (ig,1,ik1) = h_diagk (ig,1,ik1) + sum_h1
                   h_diagk (ig,2,ik1) = h_diagk (ig,2,ik1) + sum_h4
                   s_diagk (ig,1,ik1) = s_diagk (ig,1,ik1) + sum_s
                   s_diagk (ig,2,ik1) = s_diagk (ig,2,ik1) + sum_s
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   sum_h1 = 0.d0
                   sum_h4 = 0.d0
                   sum_s = 0.d0
                   DO ih = 1, nh_
                      ikb = ijkb_start + ih
                      ar = vkbk_d(ig,startk+ikb)*CONJG(vkbk_d(ig,startk+ikb))
                      sum_h1 = sum_h1 + DBLE(deeq_nc(ih,ih,na,1) * ar)
                      sum_h4 = sum_h4 + DBLE(deeq_nc(ih,ih,na,4) * ar)
                      sum_s = sum_s + DBLE(qq_at(ih,ih,na) * ar)
                   END DO
                   h_diagk (ig,1,ik1) = h_diagk (ig,1,ik1) + sum_h1
                   h_diagk (ig,2,ik1) = h_diagk (ig,2,ik1) + sum_h4
                   s_diagk (ig,1,ik1) = s_diagk (ig,1,ik1) + sum_s
                   s_diagk (ig,2,ik1) = s_diagk (ig,2,ik1) + sum_s
                END IF
             END DO
          END IF
       END DO
    ELSE
       DO nt = 1, ntyp
          IF ( type1ps(nt) ) THEN
             DO na = 1, nat
                IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   sum_h = 0.d0
                   sum_s = 0.d0
                   DO ih = 1, nh_
                      DO jh = 1, nh_
                         ikb = ijkb_start + ih
                         cv = vkbk_d(ig,ikb+startk)
                         jkb = ijkb_start + jh
                         ar = cv*CONJG(vkbk_d(ig,jkb+startk))
                         sum_h = sum_h + DBLE(deeq(ih,jh,na,current_spin) * ar)
                         sum_s = sum_s + DBLE(qq_at(ih,jh,na) * ar)
                      ENDDO
                   ENDDO
                   h_diagk (ig,1,ik1) = h_diagk (ig,1,ik1) + sum_h
                   s_diagk (ig,1,ik1) = s_diagk (ig,1,ik1) + sum_s
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF (ityp (na) == nt) THEN
                   ijkb_start = ofsbeta(na)
                   nh_ = nh(nt)
                   sum_h = 0.d0
                   sum_s = 0.d0
                   DO ih = 1, nh_
                      ikb = ijkb_start + ih
                      ar = vkbk_d(ig,startk+ikb)*CONJG(vkbk_d(ig,startk+ikb))
                      sum_h = sum_h + DBLE(deeq(ih,ih,na,current_spin) * ar)
                      sum_s = sum_s + DBLE(qq_at(ih,ih,na) * ar)
                   END DO
                   h_diagk (ig,1,ik1) = h_diagk (ig,1,ik1) + sum_h
                   s_diagk (ig,1,ik1) = s_diagk (ig,1,ik1) + sum_s
                END IF
             END DO
          END IF
       END DO
    ENDIF

    IF ( qcutz > 0.D0 ) THEN
    !
    !
       g2kink(ig,ik1) = g2kink(ig,ik1) + qcutz * &
            ( 1.D0 + erf( ( g2kink(ig,ik1) - ecfixed ) / q2sigma ) )
    !
    END IF
    !
    RETURN
    END SUBROUTINE hdiag_kernel
#endif
  !
