!
! Copyright (C) 2017 Andrea Dal Corso
! This is extracted from the Quantum ESPRESSO routine (C) and adapted
! to thermo_pw.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE psh_lanczos_step(iter, flag)
!
!   This routine receives evc1 and evc1_new. Using these two vectors it
!   makes a Lanczos step and it gives in output evc1_old that contains 
!   the input evc1 and evc1 with the new evc1. On input in evc1_new(:,:,:,1)
!   there must be A evc1 or B evc1, where A is the linear operator to 
!   tridiagonalize.
!   At each iteration it computes beta and gamma. It receives also a
!   vector in d0psi (lgamma=.TRUE.) or in in d0psi2 (lgamma=.FALSE)
!   and computes the scalar product of the current Lanczos vector with 
!   this vector, needed later to compute the correct matrix element.
!   If flag=.TRUE. orthogonalize the new Lanczos vector to the occupied
!   valence states. Should not be needed, but it helps (presently not used).
!   
!   This version of the routine is only for the pseudo-hermitean case.
!
    USE kinds,        ONLY : DP
    USE io_global,    ONLY : ionode, stdout
    USE klist,        ONLY : nks, ngk
    USE qpoint,       ONLY : ikks, ikqs, nksq
    USE lr_lanczos,   ONLY : evc1, evc1_new, evc1_old, sevc1, beta_store, &
                             gamma_store, zeta_store, iulanczos
    USE lr_global,    ONLY : evq0, sevq0, d0psi, d0psi2, size_evc1, rpert
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iter
    LOGICAL, INTENT(IN) :: flag
    !
    ! Local variables
    !
    COMPLEX(kind=dp), EXTERNAL :: lr_dot, lr_dot0
    INTEGER :: ik, ikk, ikq, ikp, ip, ipert, jpert, ig
    REAL(kind=dp) :: beta, dgamma
    COMPLEX(kind=dp) :: zeta
    !
    CALL start_clock('lanczos_step')

    CALL lr_apply_s(evc1_new(:,:,:,1), sevc1(:,:,:,1))
    !
    ! Orthogonality requirement: <v|\bar{L}|v> = 1
    !
    beta = lr_dot(evc1(:,:,:,1), sevc1(:,:,:,1))
    !
    IF ( ABS(beta)<1.0D-12 ) THEN
       !
       WRITE(stdout,'(5x,"lr_lanczos: Left and right Lanczos vectors &
                      &are orthogonal, this is a violation of oblique &
                                                              &projection")')
       !
    ELSEIF ( beta<0.0D0 ) THEN
       !
       beta = SQRT(-beta)
       dgamma = -beta
       !
    ELSEIF ( beta>0.0D0 ) THEN 
       !
       beta = SQRT(beta)
       dgamma = beta
       !
    ENDIF
    !
    beta_store (iter) = beta
    gamma_store(iter) = dgamma
!    WRITE(stdout,'(5X,"beta (",i8.8,")=",f15.6)') iter, beta
    WRITE(stdout,'(5X,"gamma(",i8.8,")=",f15.6)') iter, dgamma
    IF (ionode) WRITE(iulanczos,'(i9,e25.15)') iter, dgamma
    !
    ! Renormalize q(i) and Lq(i)
    !
    CALL ZSCAL(size_evc1,CMPLX(1.0d0/beta,0.0d0,kind=dp),evc1(1,1,1,1),1)
    CALL ZSCAL(size_evc1,CMPLX(1.0d0/beta,0.0d0,kind=dp),evc1_new(1,1,1,1),1)
    !
    !
    ! Calculation of zeta coefficients.
    ! See Eq.(35) in Malcioglu et al., Comput. Phys. Commun. 182, 1744 (2011).
    !
    IF (mod(iter,2)==0) THEN
       !
       ! Optics: In the ultrasoft case, the S operator was already
       ! applied to d0psi, so we have <S*d0psi|evc1>.
       !
       IF (rpert>1) THEN
          DO ipert = 1, rpert
             DO jpert=1, rpert
                zeta = lr_dot0(d0psi(1,1,1,jpert),evc1(1,1,1+nks*(ipert-1),1))
                zeta_store (ipert,jpert,iter) = zeta
                WRITE(stdout,'(5x,"z1= ",1x,2i6,2(1x,e22.15))')jpert,ipert, &
                                                       real(zeta), aimag(zeta)
                IF (ionode) WRITE(iulanczos,'(2e25.15)') zeta
             ENDDO
          ENDDO
       ELSE
          zeta = lr_dot(d0psi2(:,:,:),evc1(:,:,:,1))
          zeta_store (1,1,iter) = zeta
          WRITE(stdout,'(5x,"z1= ",2(1x,e22.15))') real(zeta), aimag(zeta)
          IF (ionode) WRITE(iulanczos,'(2e25.15)') zeta
       ENDIF
       !
    ELSE
       !
       zeta = (0.0d0,0.0d0)
       DO ipert=1,rpert
          DO jpert=1,rpert
             zeta_store (ipert,jpert,iter) = zeta
          ENDDO
       ENDDO
!       WRITE(stdout,'(5x,"z1= ",2(1x,e22.15))') real(zeta),aimag(zeta)
       !
    ENDIF
    !
    ! X. Ge: q(i+1) = Lq(i) - beta(i)*q(i-1); 
    ! Renormalization will be done in the begining of the next iteration.
    !
    CALL ZAXPY(size_evc1,-CMPLX(dgamma,0.0d0,KIND=DP),evc1_old,1,evc1_new,1)
    !
    ! X. Ge: To increase the stability, apply lr_ortho.
    ! Presently not used in thermo_pw
    !
    IF (flag) THEN
       !
       DO ipert=1,rpert
          DO ik=1, nksq
             ikk=ikks(ik)
             ikq=ikqs(ik)
             ikp = ik + nksq * (ipert-1)
             CALL orthogonalize(evc1_new(:,:,ikp,1), evq0(:,:,ik), ikk, ikq, &
                                             sevq0(:,:,ik),ngk(ikq),.TRUE.)
          ENDDO
       ENDDO
!
!  orthogonalize changes the sign, so we have to change sign again
!
       evc1_new(:,:,:,:)=-evc1_new(:,:,:,:)
       !
    ENDIF
    !
    ! update the vectors setting evc1 in evc1_old and evc1_new in evc1
    !
    CALL ZCOPY(size_evc1,evc1,1,evc1_old,1) ! evc1_old = evc1
    CALL ZCOPY(size_evc1,evc1_new,1,evc1,1) ! evc1 = evc1_new
    !
    CALL stop_clock('lanczos_step')
    !
    RETURN
    !
END SUBROUTINE psh_lanczos_step
