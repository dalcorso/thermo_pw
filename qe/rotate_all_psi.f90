!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! The following routine is an alternative to the routines available 
! in Quantum ESPRESSO to rotate the wavefunctions. It works for scalar 
! relativistic one-component wavefunctions and for two-component spinor 
! wavefunctions. It requires less memory. Only one complete FFT mesh
! is necessary. This memory does not decrease by increasing the number of
! processors and is allocated by each processor.
!
SUBROUTINE rotate_all_psi_tpw(ik,psic_nc,evcr,s,ftau,d_spin,has_e,gk)
  !
  !  This subroutine rotates one-component or two-component
  !  wavefunctions according to the symmetry s (actually it applies s^-1). 
  !  d_spin contains the 2x2 rotation matrix in spin space (it must be
  !  the d_spin that corresponds to the inverse of s).
  !  has_e=-1 means that also a 360 degrees rotation is applied in spin space.
  !  For one-component wavefunctions d_spin and has_e are not used.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym
  USE fft_interfaces, ONLY : fwfft
  USE gvect,     ONLY : nl
  USE wvfct,     ONLY : nbnd, npwx
  USE klist,     ONLY : ngk, igk_k
  USE noncollin_module, ONLY : npol

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ik
  INTEGER, INTENT(IN) :: s(3,3), ftau(3), gk(3), has_e
  COMPLEX(DP), INTENT(IN) :: psic_nc(dfftp%nnr,npol,nbnd), d_spin(2,2)
  COMPLEX(DP), INTENT(INOUT) :: evcr(npwx*npol,nbnd)
  COMPLEX(DP), ALLOCATABLE :: psir(:), evcr_save(:)
  COMPLEX(DP), ALLOCATABLE :: phase1(:), phase2(:), phase3(:)
  COMPLEX(DP) :: phase, term
  REAL(DP) :: arg
  INTEGER :: i, j, k, ri, rj, rk, ir, ipol, jpol, ibnd, j0, k0, idx
  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, nr3x, nr12x, my_nrxx, nrxx, npw
  INTEGER :: ind1(2), ind2(2), ss(3,3)
  INTEGER, ALLOCATABLE :: iir(:), jir(:), kir(:)
  INTEGER, ALLOCATABLE :: rir(:)
  LOGICAL :: zone_border
  !
  COMPLEX (DP), ALLOCATABLE :: psir_collect(:)
!
!  Shorten the names of the fft sizes
!
  nr1=dfftp%nr1
  nr2=dfftp%nr2
  nr3=dfftp%nr3
  nr1x=dfftp%nr1x
  nr2x=dfftp%nr2x
  nr3x=dfftp%nr3x
  nr12x= nr1x * nr2x
  nrxx=dfftp%nnr
  my_nrxx=nr1x*dfftp%my_nr2p*dfftp%my_nr3p
  npw = ngk(ik)
!
!  find i, j, and k of each point that belong to this processor and
!  the ir of its rotated by {s^-1, ftau}
!
  ALLOCATE (iir(my_nrxx))
  ALLOCATE (jir(my_nrxx))
  ALLOCATE (kir(my_nrxx))
  ALLOCATE (rir(my_nrxx))

  ss (1, 1) = s (1, 1)
  ss (2, 1) = s (2, 1) * nr1 / nr2
  ss (3, 1) = s (3, 1) * nr1 / nr3
  ss (1, 2) = s (1, 2) * nr2 / nr1
  ss (2, 2) = s (2, 2)
  ss (3, 2) = s (3, 2) * nr2 / nr3
  ss (1, 3) = s (1, 3) * nr3 / nr1
  ss (2, 3) = s (2, 3) * nr3 / nr2
  ss (3, 3) = s (3, 3)

  j0 = dfftp%my_i0r2p
  k0 = dfftp%my_i0r3p
  rir=0
  DO ir = 1, my_nrxx
     idx = ir - 1
     k   = idx / (nr1x*dfftp%my_nr2p)
     idx = idx - (nr1x*dfftp%my_nr2p)*k
     k   = k + k0
     j   = idx / nr1x
     idx = idx - nr1x * j
     j   = j + j0
     i   = idx
     i=i + 1
     j=j + 1
     k=k + 1
     iir(ir)=i
     jir(ir)=j
     kir(ir)=k
     IF (i > nr1 .OR. j > nr2 .OR. k > nr3) CYCLE
     CALL ruotaijk_tpw (ss, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk ) 
     rir(ir)=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
  ENDDO
!
!  In the cases in which we are at zone border computes all the
!  necessary phases
!

  zone_border = (gk(1) /= 0 .OR. gk(2) /= 0 .OR. gk(3) /= 0 )

  IF (zone_border) THEN
     ALLOCATE(phase2(nr2))
     ALLOCATE(phase3(nr3))
     ALLOCATE(phase1(nr1))
     IF (gk(1) /= 0) THEN
        phase1(1)=(1.0_DP,0.0_DP)
        arg = tpi*gk(1)/DBLE(nr1)
        term= CMPLX(COS(arg), SIN(arg), KIND=DP)
        DO i = 2, nr1
           phase1(i) = phase1(i-1) * term
        END DO
     ELSE
        phase1(:)=(1.0_DP,0.0_DP)
     ENDIF
     IF (gk(2) /= 0) THEN
        phase2(1)=(1.0_DP,0.0_DP)
        arg = tpi*gk(2)/DBLE(nr2)
        term= CMPLX(COS(arg), SIN(arg), KIND=DP)
        DO j = 2, nr2
           phase2(j) = phase2(j-1) * term
        END DO
     ELSE
        phase2(:)=(1.0_DP,0.0_DP)
     ENDIF
     IF (gk(3) /= 0) THEN
        phase3(1)=(1.0_DP,0.0_DP)
        arg = tpi*gk(3)/DBLE(nr3)
        term= CMPLX(COS(arg), SIN(arg), KIND=DP)
        DO k = 2, nr3
           phase3(k) = phase3(k-1) * term
        END DO
     ELSE
        phase3(:)=(1.0_DP,0.0_DP)
     ENDIF
  ENDIF
  !
  !  Now we allocate a complete FFT mesh where we collect each wavefunction
  !  a distributed fft used to pass in reciprocal space
  !  and a vector evc used as an auxiliary space before applying d_spin the
  !  rotation in spin space
  !
  ALLOCATE(psir_collect(nr12x*nr3x))
  ALLOCATE(psir(nrxx))
  ALLOCATE(evcr_save(npwx*npol))

  npw = ngk(ik)
  DO ipol=1,npol
     ind1(ipol)=1 + (ipol-1)*npwx
     ind2(ipol)=npw + (ipol-1)*npwx
  ENDDO
  evcr=(0.d0,0.d0)
  DO ibnd=1,nbnd
     psir = ( 0.D0, 0.D0 )
     psir_collect=(0.d0,0.d0)
     DO ipol=1,npol
        CALL cgather_sym( dfftp, psic_nc(:,ipol,ibnd), psir_collect)
        IF (zone_border) THEN
           DO ir = 1, my_nrxx
              IF (rir(ir)==0) CYCLE
              i=iir(ir)
              j=jir(ir)
              k=kir(ir)
              phase=phase1(i)*phase2(j)*phase3(k)
              psir(ir)=psir_collect(rir(ir))*phase
           ENDDO
        ELSE
           DO ir = 1, my_nrxx
              IF (rir(ir)==0) CYCLE
              psir(ir)=psir_collect(rir(ir))
           ENDDO
        ENDIF
        CALL fwfft ('Dense', psir, dfftp)
        evcr(ind1(ipol):ind2(ipol),ibnd) = psir(nl(igk_k(1:npw,ik)))
     ENDDO
     IF (npol==2) THEN
        evcr_save(:)=evcr(:,ibnd)
        evcr(:,ibnd)=(0.0_DP,0.0_DP)
        DO ipol=1,npol
           DO jpol=1,npol
              evcr(ind1(ipol):ind2(ipol),ibnd)=&
                evcr(ind1(ipol):ind2(ipol),ibnd)+ &
                    d_spin(ipol,jpol)*evcr_save(ind1(jpol):ind2(jpol))
           ENDDO
        ENDDO
     ENDIF
  ENDDO

  IF (npol==2.AND.has_e==-1) evcr=-evcr
  !
  DEALLOCATE(iir)
  DEALLOCATE(jir)
  DEALLOCATE(kir)
  DEALLOCATE(rir)
  DEALLOCATE(evcr_save)
  DEALLOCATE(psir)
  DEALLOCATE(psir_collect)

  IF (zone_border) THEN
     DEALLOCATE(phase1)
     DEALLOCATE(phase2)
     DEALLOCATE(phase3)
  ENDIF

  RETURN
END SUBROUTINE rotate_all_psi_tpw

