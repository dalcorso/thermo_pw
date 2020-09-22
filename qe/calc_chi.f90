!
! Copyright (C) 2018 Andrea Dal Corso
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This routine is taken from TDDFPT/tools/calculate_spectrum.f90 and modified
! to interface with thermo_pw.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE extrapolate()
  !-----------------------------------------------------------------------
  !
  ! This subroutine applies the "extrapolation" scheme 
  ! for extrapolating the reduced matrix.
  !
  USE kinds,        ONLY : DP
  USE lr_lanczos,   ONLY : extrapolation, lanczos_steps, beta_store, &
                           gamma_store, beta_store_ext, gamma_store_ext, &
                           lanczos_steps_ext
  USE io_global,    ONLY : stdout
  IMPLICIT NONE

  INTEGER  :: ip, counter, i
  REAL(DP) :: average, av_amplitude
  LOGICAL  :: skip
  !
  !  Terminatore
  !
  skip = .FALSE.
  !
  IF (trim(extrapolation)/="no") THEN
     !
     average = 0.d0
     av_amplitude = 0.d0
     !
     beta_store_ext(1:lanczos_steps-1)=beta_store(2:lanczos_steps)
     gamma_store_ext(1:lanczos_steps-1)=gamma_store(2:lanczos_steps)
     beta_store_ext(lanczos_steps)=beta_store(lanczos_steps)
     gamma_store_ext(lanczos_steps)=gamma_store(lanczos_steps)
     !
     counter=0
     !
     DO i=151,lanczos_steps
        !
        IF (skip .EQV. .TRUE.) THEN
           skip=.FALSE.
           CYCLE
        ENDIF
        !
        IF (mod(i,2)==1) THEN
           !
           IF ( i/=151 .AND. ABS( beta_store_ext(i)-average/counter ) &
                                                                  > 2.d0 ) THEN
              !
              !if ( i.ne.151 .and. counter == 0) counter = 1
              skip=.TRUE.
              !
           ELSE
              !
              average=average+beta_store_ext(i)
              av_amplitude=av_amplitude+beta_store_ext(i)
              counter=counter+1
!              print *, "t1 ipol",i, ip,"conter", counter, "av_amp",av_amplitude(ip)
              !
           ENDIF
           !
        ELSE
           !
           IF ( i/=151 .AND. abs( beta_store_ext(i)-average/counter ) &
                                                                  > 2.d0 ) THEN
              !
              !if ( i.ne.151 .and. counter == 0) counter = 1
              skip=.TRUE.
              !
           ELSE
              !
              average=average+beta_store_ext(i)
              av_amplitude=av_amplitude-beta_store_ext(i)
              counter=counter+1
           ENDIF
!             print *, "t2 ipol",i, ip,"conter", counter, "av_amp",&
!                                                av_amplitude
           !
        ENDIF
        !
     ENDDO
     !
     IF (counter>0) THEN
        average=average/counter
        av_amplitude=av_amplitude/counter
     ENDIF
     !print *, "t3 ipol",ip,"av_amp",av_amplitude(ip)
     !
     WRITE(stdout,'(5x,"Lanczos coefficients:")')
     WRITE(stdout,'(5x,"Average =",3F15.8)') average
     WRITE(stdout,'(5x,"Average oscillation amplitude =",F15.8)') av_amplitude
     !
     IF (trim(extrapolation)=="constant") av_amplitude=0
     !
     !
     DO i=lanczos_steps,lanczos_steps_ext
        !
        IF (MOD(i,2)==1) THEN
           !
           beta_store_ext(i)=average+av_amplitude
           gamma_store_ext(i)=average+av_amplitude
           !
        ELSE
           !
           beta_store_ext(i)=average-av_amplitude
           gamma_store_ext(i)=average-av_amplitude
           !
        ENDIF
        !
     ENDDO
   !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE extrapolate

SUBROUTINE calc_chi(freq,broad,chi)
  !----------------------------------------------------------------------------
  !
  ! This subroutine calculates the susceptibility and the dielectric
  ! constant in the optical case.
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : fpi, e2
  USE lr_lanczos, ONLY : beta_store, beta_store_ext, gamma_store_ext, &
                         zeta_store, lanczos_steps, lanczos_steps_ext
  USE lr_global,  ONLY : rpert
  USE control_lr, ONLY : lgamma
  USE symme,      ONLY : symmatrix, crys_to_cart
  USE lsda_mod,   ONLY : nspin
  USE cell_base,  ONLY : omega
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: freq
  REAL(DP), INTENT(IN) :: broad
  COMPLEX(DP), INTENT(INOUT) :: chi(rpert,rpert)

  COMPLEX(DP) :: omeg_c, a(lanczos_steps_ext), b(lanczos_steps_ext), &
                 c(lanczos_steps_ext), r(lanczos_steps_ext)
  COMPLEX(DP) :: ZDOTC
  INTEGER     :: ip, ip2, i, info, ipert
  REAL(DP)    :: degspin, epsil(3,3), epsili(3,3), fact

  degspin=2.0_DP
  IF (nspin /=1) degspin=1.0_DP
  !
  omeg_c = CMPLX(freq,broad,DP)
  !
  DO ip =1, rpert
     !
     a(:) = omeg_c
     !
     DO i = 1,lanczos_steps_ext-1
        !
        b(i) = CMPLX(-beta_store_ext(i),0.0d0,DP)
        c(i) = CMPLX(-gamma_store_ext(i),0.0d0,DP)
        !
     ENDDO
     !
     r(:) = (0.0d0,0.0d0)
     r(1) = (1.0d0,0.0d0)
     !
     ! |w_t|=(w-L) |1,0,0,...,0|
     ! 
     CALL ZGTSV(lanczos_steps_ext,1,b,a,c,r(:),lanczos_steps_ext,info)
     !
     IF (info /= 0) CALL errore("calc_chi", "Unable to solve &
                                                      &tridiagonal system",1)
     !
     ! p=-div.rho'
     ! p= chi . E
     ! Thus, chi = - <zeta|w_t>
     !
     ! Notice that brodening has a positive sign, 
     ! thus the abs. coefficient is Im(tr(chi)) not -Im(Tr(chi)) as usual
     ! 
     DO ip2 = 1,rpert
         !
         chi(ip,ip2) = ZDOTC(lanczos_steps,zeta_store(ip,ip2,:),1,r(:),1)
         !
         ! Multiplication with a norm
         !
         chi(ip,ip2) = chi(ip,ip2) * CMPLX(beta_store(1),0.0d0,DP)
         !
         ! The response charge density is defined as 2*evc0*q, see Eq. (43) in
         ! JCP 128, 154105 (2008). 
         ! Therefore, the dipole is given by 2*degspin* zeta^T *
         ! (w-T^lanczos_steps)^-1 * e_1. See also Eq. (15) in that paper.
         ! Optics: The minus sign accounts for the negative electron charge
         ! (perturbation is -e E x, rather than E x)
         !
         ! 
         IF (lgamma) THEN
            chi(ip,ip2) = chi(ip,ip2) * CMPLX(-2.d0*degspin, 0.d0, DP)
         ELSE
            chi(ip,ip2) = chi(ip,ip2) * CMPLX( 2.d0*degspin, 0.d0, DP)
         ENDIF
         !
     ENDDO
     !
  ENDDO
  !
  !   In the optical case symmetrize the tensor. Real and imaginary parts
  !   are symmetrized independently. Finally add one on the diagonal.
  !
  IF (rpert==3) THEN
     fact=fpi*e2/omega
     epsil(:,:) = DBLE(chi(:,:)) 
     CALL crys_to_cart(epsil)
     CALL symmatrix(epsil)
     epsili(:,:) = AIMAG(chi(:,:)) 
     CALL crys_to_cart(epsili)
     CALL symmatrix(epsili)
     chi(:,:) = CMPLX(epsil(:,:), epsili(:,:))
     DO ipert=1,rpert
        chi(ipert,ipert)=chi(ipert,ipert)*fact+CMPLX(1.0_DP,0.0_DP)
     ENDDO
  ENDIF     

  RETURN
  !
END SUBROUTINE calc_chi

