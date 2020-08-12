!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE addusddenseq (drhoscf, dbecsum)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the change of the charge and magnetization
  !  densities due to an electric field perturbation
  !  the part due to the US augmentation.
  !  It assumes that the array dbecsum has already accumulated the
  !  change of the becsum term.
  !  The expression implemented is given in Eq. B32 of PRB 64, 235118
  !  (2001) with b=c=0.
  !

  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE cell_base, ONLY : tpiba
  USE fft_base,  ONLY: dfftp
  USE fft_interfaces, ONLY: invfft
  USE gvect,     ONLY : g, gg, ngm, eigts1, eigts2, eigts3, mill
  USE uspp,      ONLY: okvan
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm
  USE noncollin_module, ONLY : nspin_mag

  USE qpoint, ONLY : xq, eigqts
  IMPLICIT NONE
  !
  !   the dummy variables
  !

  ! input: if zero does not compute drho
  ! input: the number of perturbations

  COMPLEX(DP) :: drhoscf(dfftp%nnr,nspin_mag,1), &
                 dbecsum(nhm*(nhm+1)/2,nat,nspin_mag,1)

  ! inp/out: change of the charge density
  ! input: sum over kv of bec
  !
  !     here the local variables
  !

  INTEGER :: ig, na, nt, ih, jh, ijh, is

  ! counters

  REAL(DP), ALLOCATABLE  :: qmod(:), qpg(:,:), ylmk0(:,:)
  ! the modulus of q+G
  ! the spherical harmonics

  COMPLEX(DP) :: zsum
  COMPLEX(DP), ALLOCATABLE ::  sk (:), qg (:), qgm (:), aux (:,:)
  ! the structure factor
  ! work space

  IF (.NOT.okvan) RETURN
  CALL start_clock ('addusddenseq')
  ALLOCATE (aux(  ngm, nspin_mag))
  ALLOCATE (sk (  ngm))
  ALLOCATE (qg (  dfftp%nnr))
  ALLOCATE (ylmk0(ngm , lmaxq * lmaxq))
  ALLOCATE (qgm  (ngm))
  ALLOCATE (qmod (ngm))
  ALLOCATE (qpg(3, ngm))
  !
  !  And then we compute the additional charge in reciprocal space
  !
  CALL setqmod (ngm, xq, g, qmod, qpg)
  CALL ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
  DO ig = 1, ngm
     qmod (ig) = sqrt (qmod (ig) ) * tpiba
  ENDDO

  aux (:,:) = (0.d0, 0.d0)
  DO nt = 1, ntyp
     IF (upf(nt)%tvanp ) THEN
        ijh = 0
        DO ih = 1, nh (nt)
           DO jh = ih, nh (nt)
              CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
              DO na = 1, nat
                 IF (ityp (na) == nt) then
                    !
                    ! calculate the structure factor
                    !
                    DO ig = 1, ngm
                       sk(ig)=eigts1(mill(1,ig),na)*eigts2(mill(2,ig),na) &
                             *eigts3(mill(3,ig),na)*eigqts(na)*qgm(ig)
                    ENDDO
                    !
                    !  And qgmq and becp and dbecq
                    !
                    DO is=1,nspin_mag
                       zsum = dbecsum (ijh, na, is, 1)
                       CALL zaxpy(ngm,zsum,sk,1,aux(1,is),1)
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  !     convert aux to real space
  !
  DO is=1,nspin_mag
     qg (:) = (0.d0, 0.d0)
     qg (dfftp%nl (:) ) = aux (:, is)
     CALL invfft ('Rho', qg, dfftp)
     drhoscf(:,is,1) = drhoscf(:,is,1) + 2.d0*qg(:)
  ENDDO

  DEALLOCATE (qpg)
  DEALLOCATE (qmod)
  DEALLOCATE (qgm)
  DEALLOCATE (ylmk0)
  DEALLOCATE (qg)
  DEALLOCATE (sk)
  DEALLOCATE (aux)

  CALL stop_clock ('addusddenseq')
  RETURN
END SUBROUTINE addusddenseq
