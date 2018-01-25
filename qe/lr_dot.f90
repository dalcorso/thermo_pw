!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION lr_dot(x,y)
  !---------------------------------------------------------------------
  !
  ! This subroutine calculates a dot product of the conjugate 
  ! of a complex vector x and a complex vector y 
  ! (sums over the bands and k-points).
  !
  !
  USE kinds,                ONLY : dp
  USE klist,                ONLY : nks, wk, ngk
  USE wvfct,                ONLY : npwx, nbnd, wg
  USE control_flags,        ONLY : gamma_only
  USE lr_global,            ONLY : rpert
  USE noncollin_module,     ONLY : noncolin, npol
  USE lsda_mod,             ONLY : nspin
  USE control_lr,           ONLY : nbnd_occ, lgamma
  USE qpoint,               ONLY : nksq
  !
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm

  IMPLICIT NONE
  !
  COMPLEX(kind=dp) :: x(npwx*npol,nbnd,nksq*rpert), &
                      y(npwx*npol,nbnd,nksq*rpert)

  COMPLEX(kind=dp) :: lr_dot
  COMPLEX(kind=dp) :: temp_k
  REAL(kind=dp) :: temp_gamma
  INTEGER :: ibnd, ik, ipert, ikp
  !
  CALL start_clock ('lr_dot')
  !
  lr_dot = (0.0d0,0.0d0)
  temp_gamma = 0.0d0
  temp_k = (0.0d0,0.0d0)
  !
  IF (gamma_only) THEN
     !
     CALL lr_dot_gamma()
     lr_dot = CMPLX(temp_gamma,0.0d0,DP)
     !
  ELSE
     !
     CALL lr_dot_k()
     lr_dot = temp_k
     !
  ENDIF
  IF (nspin==1) lr_dot=lr_dot/2.0_DP ! This factor is on the kpoint weights
                                       ! for spin, but should not be used here

  !
  CALL stop_clock ('lr_dot')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_dot_gamma()
    USE gvect,                ONLY : gstart

    IMPLICIT NONE
    REAL(kind=dp), EXTERNAL    :: DDOT
    INTEGER :: ibnd
    !
    ! Optical case: gamma_only
    ! Noncollinear case is not implemented
    !
    DO ibnd=1,nbnd
       !
       temp_gamma = temp_gamma + 2.D0*wg(ibnd,1)*DDOT(2*ngk(1),x(:,ibnd,1),1,&
                                                               y(:,ibnd,1),1)
       !
       ! G=0 has been accounted twice, so we subtract one contribution.
       !
       IF (gstart==2) temp_gamma = temp_gamma - wg(ibnd,1)*DBLE(x(1,ibnd,1))&
                                                          *DBLE(y(1,ibnd,1))
       !
    ENDDO
    !
    CALL mp_sum(temp_gamma, intra_bgrp_comm)
    !
    RETURN
    !
  END SUBROUTINE lr_dot_gamma
  !
  SUBROUTINE lr_dot_k
    !
    USE qpoint, ONLY : ikks, ikqs
    IMPLICIT NONE
    INTEGER :: ik, ikp, ikk, ikq, npwq, ipert, ibnd
    COMPLEX(kind=dp), EXTERNAL :: ZDOTC
    REAL(DP) :: weigth
    !
    DO ipert=1, rpert
       DO ik=1, nksq   
          ikp = ik + nksq * (ipert-1)
          ikk  = ikks(ik)
          ikq  = ikqs(ik)
          npwq = ngk(ikq)
          DO ibnd=1,nbnd_occ(ikk)
             !
             IF (lgamma) THEN
                weigth=wg(ibnd,ikk)
             ELSE
                weigth=wk(ikk)
             ENDIF
             temp_k = temp_k + weigth * ZDOTC(npwq,x(1,ibnd,ikp),1,&
                                                           y(1,ibnd,ikp),1)
             IF (noncolin) &
                temp_k = temp_k+weigth*ZDOTC(npwq,x(1+npwx,ibnd,ikp),1,&
                                                        y(1+npwx,ibnd,ikp),1)
          ENDDO
       ENDDO
    ENDDO
    !
    CALL mp_sum(temp_k, inter_pool_comm)
    CALL mp_sum(temp_k, intra_bgrp_comm)
    !
    RETURN
    !
  END SUBROUTINE lr_dot_k

END FUNCTION lr_dot

!-----------------------------------------------------------------------
FUNCTION lr_dot0(x,y)
  !---------------------------------------------------------------------
  !
  ! This subroutine calculates a dot product of the conjugate 
  ! of a complex vector x and a complex vector y 
  ! (sums over the bands and k-points).
  !
  !
  USE kinds,                ONLY : dp
  USE klist,                ONLY : nks, wk, ngk
  USE wvfct,                ONLY : npwx, nbnd, wg
  USE control_flags,        ONLY : gamma_only
  USE lr_global,            ONLY : rpert
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol
  USE control_lr,           ONLY : nbnd_occ, lgamma
  USE qpoint,               ONLY : nksq
  !
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm

  IMPLICIT NONE
  !
  COMPLEX(kind=dp) :: x(npwx*npol,nbnd,nksq), &
                      y(npwx*npol,nbnd,nksq)

  COMPLEX(kind=dp) :: lr_dot0
  COMPLEX(kind=dp) :: temp_k
  INTEGER :: ibnd, ik, ipert, ikp
  !
  CALL start_clock ('lr_dot')
  !
  lr_dot0 = (0.0d0,0.0d0)
  temp_k = (0.0d0,0.0d0)
  !
     !
  CALL lr_dot_k()
  lr_dot0 = temp_k
  IF (nspin==1) lr_dot0=lr_dot0/2.0_DP ! This factor is on the kpoint weights
                                       ! for spin, but should not be used here
     !
  !
  CALL stop_clock ('lr_dot')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_dot_k
    !
    USE qpoint, ONLY : ikks, ikqs
    IMPLICIT NONE
    INTEGER :: ik, ikk, ikq, npwq, ipert, ibnd
    COMPLEX(kind=dp), EXTERNAL :: ZDOTC
    REAL(DP) :: weigth
    !
    DO ik=1, nksq   
       ikk  = ikks(ik)
       ikq  = ikqs(ik)
       npwq = ngk(ikq)
       DO ibnd=1,nbnd_occ(ikk)
          !
          IF (lgamma) THEN
             weigth=wg(ibnd,ikk)
          ELSE
             weigth=wk(ikk)
          ENDIF
          temp_k = temp_k + weigth * ZDOTC(npwq,x(1,ibnd,ik),1,&
                                                        y(1,ibnd,ik),1)
          IF (noncolin) &
             temp_k = temp_k+weigth*ZDOTC(npwq,x(1+npwx,ibnd,ik),1,&
                                                     y(1+npwx,ibnd,ik),1)
       ENDDO
    ENDDO
    !
    CALL mp_sum(temp_k, inter_pool_comm)
    CALL mp_sum(temp_k, intra_bgrp_comm)
    !
    RETURN
    !
  END SUBROUTINE lr_dot_k

END FUNCTION lr_dot0

SUBROUTINE lr_prec(res, pres)

USE kinds,       ONLY : DP
USE lr_cg,       ONLY : prec_vec
USE klist,       ONLY : ngk
USE qpoint,      ONLY : ikks, ikqs
USE wvfct,       ONLY : npwx, nbnd
USE control_lr,  ONLY : nbnd_occ
USE noncollin_module, ONLY : noncolin, npol
USE lr_global,   ONLY : rpert, sevq0, evq0
USE qpoint,      ONLY : nksq

IMPLICIT NONE
COMPLEX(DP), INTENT(IN)  :: res(npwx*npol, nbnd, nksq*rpert)
COMPLEX(DP), INTENT(OUT) :: pres(npwx*npol, nbnd, nksq*rpert)

INTEGER :: ipert, ik, ikp, ikk, ikq, ig, npwq, ibnd

pres=(0.0_DP, 0.0_DP)
DO ipert=1, rpert
   DO ik=1, nksq   
      ikp = ik + nksq * (ipert-1)
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npwq = ngk(ikq)
      DO ibnd=1,nbnd_occ(ikk)
         pres(1:npwq,ibnd,ikp)=-prec_vec(1:npwq,ibnd,ik)*res(1:npwq,ibnd,ikp)
         IF (noncolin) &
            pres(npwx+1:npwx+npwq,ibnd,ikp)=&
              -prec_vec(npwx+1:npwx+npwq,ibnd,ik)*res(npwx+1:npwx+npwq,ibnd,ikp)
      ENDDO
      CALL orthogonalize(pres(1,1,ikp), evq0(1,1,ik), ikk, ikq, &
                                      sevq0(1,1,ik), npwq, .TRUE.)
   ENDDO
ENDDO

RETURN
END SUBROUTINE lr_prec
