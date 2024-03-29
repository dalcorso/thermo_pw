!
! Copyright (C) 2001-2003 PWSCF group
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE addusdens1d_so (plan, prho)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation. This is done only along the G_z direction in
  !  reciprocal space. The output of the routine is the planar average
  !  of the charge density. This routine generalizes the corresponding
  !  routine addusdens1d_so to the noncollinear case and it makes
  !  also the planar average of the three components of the magnetization
  !  densities
  !
  USE kinds,      ONLY: DP
  USE cell_base,  ONLY: alat, omega, celldm, tpiba
  USE ions_base,  ONLY: nat, ntyp => nsp, ityp
  USE fft_base,   ONLY: dfftp
  USE fft_scalar, ONLY: cft_1z
  USE gvect,      ONLY: eigts1, eigts2, eigts3, mill
  USE lsda_mod,   ONLY: nspin, current_spin
  USE uspp,       ONLY: becsum
  USE uspp_param, ONLY: upf, lmaxq, nh
  USE mp_global,  ONLY: intra_bgrp_comm
  USE mp,         ONLY: mp_sum

  !
  !     here the local variables
  !
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN) :: prho(dfftp%nnr,nspin)
  REAL(DP), INTENT(INOUT) :: plan(dfftp%nr3,nspin)
  INTEGER :: ig, na, nt, ih, jh, ijh, ispin, ngm1d, ig1dto3d (dfftp%nr3), &
       igtongl1d (dfftp%nr3), nl1d (dfftp%nr3)
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic types
  ! counter on beta functions
  ! counter on beta functions
  ! composite index ih jh
  ! the number of 1D G vectors on this processor
  ! correspondence 1D with 3D G vectors
  ! the correspondence 1D with the 3D shells
  ! correspondence 1D FFT mesh G with array G

  REAL(DP) :: dimz, g1d (3, dfftp%nr3), gg1d (dfftp%nr3), qmod (dfftp%nr3), &
              qgr(dfftp%nr3), qgi(dfftp%nr3), ylmk0(dfftp%nr3, lmaxq * lmaxq)
  !  the planar average
  !  dimension along z
  !  ngm1d 3D vectors with the 1D G of this proc
  !  ngm1d scalars with the modulus of 1D G
  ! the modulus of G
  ! real and
  ! imaginary part of qg
  ! the spherical harmonics

  COMPLEX(DP) :: skk, qg (dfftp%nr3x, nspin), qgout (dfftp%nr3)
  ! auxiliary variable
  ! auxiliary space for the charge
  ! auxiliary variable for FFT
  ! auxiliary variable for rho(G,nspin)
  COMPLEX(DP), ALLOCATABLE :: qgm(:), aux (:,:)


  CALL ggen1d (ngm1d, g1d, gg1d, ig1dto3d, nl1d, igtongl1d)
  ALLOCATE (qgm(ngm1d)) 
  ALLOCATE (aux(ngm1d,nspin))
  DO ig = 1, ngm1d
     qmod (ig) = sqrt (gg1d (ig) ) * tpiba
  ENDDO
  aux(:,:) = (0.d0, 0.d0)

  qg(:,:) = (0.d0, 0.d0)
  IF (ngm1d > 0) THEN
     CALL ylmr2 (lmaxq * lmaxq, ngm1d, g1d, gg1d, ylmk0)
     DO nt = 1, ntyp
        IF (upf(nt)%tvanp  ) THEN
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 CALL qvan2 (ngm1d, ih, jh, nt, qmod, qgm, ylmk0)
                 ijh = ijh + 1
                 DO na = 1, nat

                    IF (ityp (na) == nt) THEN
                       !
                       !  Multiply becsum and qg with the correct structure factor
                       !
                       DO ig = 1, ngm1d

                          skk = eigts1 (mill(1,ig1dto3d (ig) ), na) * &
                                eigts2 (mill(2,ig1dto3d (ig) ), na) * &
                                eigts3 (mill(3,ig1dto3d (ig) ), na)
                          DO ispin=1, nspin
                             aux (ig,ispin) = aux (ig,ispin) + &
                                qgm (ig) * skk * becsum (ijh, na, ispin)
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     !     adds to the charge density and converts to real space
     !
     DO ispin=1, nspin
        DO ig = 1, ngm1d
           qg (nl1d (ig), ispin ) = aux (ig, ispin) + &
                                    prho (dfftp%nl (ig1dto3d (ig) ), ispin)
        ENDDO
     END DO
  ENDIF
  CALL mp_sum(  qg, intra_bgrp_comm )
  dimz = alat * celldm (3)
  DO ispin=1,nspin
     CALL cft_1z (qg(:,ispin), 1, dfftp%nr3, dfftp%nr3x, 1, qgout)
     DO ig = 1, dfftp%nr3
        plan (ig,ispin) = DBLE(qgout (ig)) * omega / dimz
     ENDDO
  ENDDO
  DEALLOCATE (aux, qgm)

  RETURN
END SUBROUTINE addusdens1d_so
