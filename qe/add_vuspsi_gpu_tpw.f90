!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! Copyright (C) 2023 Andrea Dal Corso (for the multi k generalization)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE add_vuspsik_gpu( lda, n, m, hpsi_d, ik )
  !----------------------------------------------------------------------------
  !! This routine applies the Ultra-Soft Hamiltonian to a
  !! vector psi and puts the result in hpsi. 
  !! It requires the products of psi with all beta functions
  !! in array becp(nkb,m) (calculated by calbec).
  !
  USE kinds,           ONLY: DP
  USE ions_base,       ONLY: nat, ntyp => nsp, ityp
  USE lsda_mod,        ONLY: current_spin, isk
  USE many_k_mod,      ONLY: becpk_d, psk_d, vkbk_d
  USE control_flags,   ONLY: gamma_only
  USE noncollin_module
  USE uspp,            ONLY: ofsbeta, nkb, vkb, deeq, deeq_nc
  USE uspp_param,      ONLY: nh, nhm
                             
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  INTEGER, INTENT(IN) :: ik
  !! the k point
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(lda*npol,m)
  !! V_US|psi> is added to hpsi
  REAL(DP) :: time, scnds, tot_time
  !
  ! ... here the local variables
  !
#if defined(__CUDA)
  attributes(DEVICE) :: hpsi_d
  !
  !
  CALL start_clock_gpu( 'add_vuspsi' )  
  !
  IF ( gamma_only ) THEN
     !
     CALL errore('add_vuspsik_gpu','many_k and gamma incompatible',1)
     !
  ELSE IF ( noncolin) THEN
     !
     CALL add_vuspsik_nc_gpu ()
     !
  ELSE
     !
     CALL add_vuspsik_k_gpu()
     !
  END IF
  !
  CALL stop_clock_gpu( 'add_vuspsi' )  
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsik_k_gpu()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_gamma for comments
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#endif
       USE device_fbuff_m, ONLY : dev_buf
       !
       IMPLICIT NONE
       COMPLEX(DP), POINTER :: deeaux_d (:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: i, j, k, jkb, ikb, ih, jh, na, nt, ibnd, nhnt, isk_

#if defined(__CUDA)
       ATTRIBUTES( DEVICE ) :: deeaux_d
#endif
       !
       IF ( nkb == 0 ) RETURN
       !
       CALL dev_buf%lock_buffer(deeaux_d, (/ nhm, nhm /), ierr ) !ALLOCATE ( deeaux_d(nhm, nhm) )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate deeaux_d ', ABS( ierr ) )
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          !
          nhnt = nh(nt)
          !
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                ! deeq is real: copy it into a complex variable to perform
                ! a zgemm - simple but sub-optimal solution
                !
                !deeaux_d(:,:) = CMPLX(deeq(1:nh(nt),1:nh(nt),na,current_spin), 0.0_dp, KIND=dp )
                !
                isk_=isk(ik)
                !$acc parallel loop collapse(2) present(deeq)
                DO j = 1, nhnt
                   DO k = 1, nhnt
                      deeaux_d(k,j) = CMPLX(deeq(k,j,na,isk_), 0.0_dp, KIND=DP )
                   END DO
                END DO

                !
                call cpu_time(time)
                time=scnds-time
                tot_time=tot_time+time

                CALL ZGEMM('N','N', nhnt, m, nhnt, (1.0_dp,0.0_dp), &
                           deeaux_d, nhm, becpk_d(ofsbeta(na)+1,1,1,ik), nkb, &
                          (0.0_dp, 0.0_dp), psk_d(ofsbeta(na)+1,1,1,ik), nkb )

                call cpu_time(time)
                time=scnds-time
                tot_time=tot_time+time
                !
             END IF
             !
          END DO
          !
       END DO
       CALL dev_buf%release_buffer(deeaux_d, ierr) ! DEALLOCATE (deeaux_d)
       !
!!$acc data present(vkbk_d(:,:))
!!$acc host_data use_device(vkbk_d)
!       CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ) , &
!                   vkbk_d(:,nkb*(ik-1)+1:nkb*ik), &
!                   lda, ps_d, nkb, ( 1.D0, 0.D0 ) , hpsi_d, lda )
!!$acc end host_data
!!$acc end data
       !
!       CALL dev_buf%release_buffer(ps_d, ierr) !DEALLOCATE (ps_d)
       !
       RETURN
       !
     END SUBROUTINE add_vuspsik_k_gpu
     !  
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsik_nc_gpu()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_k for comments
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#endif
       USE device_fbuff_m,      ONLY : dev_buf
       IMPLICIT NONE
       COMPLEX(DP), POINTER :: ps_d (:,:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: na, nt, ibnd

#if defined(__CUDA)
       ATTRIBUTES( DEVICE ) :: ps_d
#endif
       !
       IF ( nkb == 0 ) RETURN
       !
       ! ALLOCATE (ps_d( nkb, npol, m), STAT=ierr )
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, npol, m /), ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating ps_d ', ABS( ierr ) )
       !
       !  OPTIMIZE HERE: possibly streamline
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                !$acc host_data use_device(deeq_nc)
                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,1), nhm, becpk_d(ofsbeta(na)+1,1,1,ik), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1,1), 2*nkb )

                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,2), nhm, becpk_d(ofsbeta(na)+1,2,1,ik), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1,1), 2*nkb )


                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,3), nhm, becpk_d(ofsbeta(na)+1,1,1,ik), 2*nkb, &
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,2,1), 2*nkb )


                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                           deeq_nc(1,1,na,4), nhm, becpk_d(ofsbeta(na)+1,2,1,ik), 2*nkb, &
                          (1.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,2,1), 2*nkb )
                !$acc end host_data
                !
!                DO ibnd = 1, m
!                   !
!                   DO jh = 1, nh(nt)
!                      !
!!$acc parallel loop present(deeq_nc)
!                      DO ih = 1, nh(nt)
!                         !
!                         ikb = ijkb0 + ih
!                         jkb = ijkb0 + jh   
!                         becpup_jkb = becp_nc_d(jkb,1,ibnd)
!                         becpdn_jkb = becp_nc_d(jkb,2,ibnd)
!                         !
!                         ps_d(ikb,1,ibnd) = ps_d(ikb,1,ibnd) +   & 
!                              deeq_nc(ih,jh,na,1)*becpup_jkb + & 
!                              deeq_nc(ih,jh,na,2)*becpdn_jkb
!                         ps_d(ikb,2,ibnd) = ps_d(ikb,2,ibnd) +   & 
!                              deeq_nc(ih,jh,na,3)*becpup_jkb + &
!                              deeq_nc(ih,jh,na,4)*becpdn_jkb
!                         !
!                      END DO
!                      !
!                   END DO
!                   !
!                END DO
                !
             END IF
             !
          END DO
          !
       END DO
       !
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
       call ZGEMM ('N', 'N', n, m*npol, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps_d, nkb, ( 1.D0, 0.D0 ) , hpsi_d, lda )
!$acc end host_data
!$acc end data
       !
       CALL dev_buf%release_buffer(ps_d, ierr ) ! DEALLOCATE (ps_d)
       !
       RETURN
       !
     END SUBROUTINE add_vuspsik_nc_gpu
#endif
     !
END SUBROUTINE add_vuspsik_gpu
!
