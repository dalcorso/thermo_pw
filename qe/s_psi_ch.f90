!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE s_psi_ch( lda, n, m, psi, spsi, ik, id )
  !--------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  !! \(\textit{Wrapper routine}\): performs bgrp parallelization on 
  !! non-distributed bands if suitable and required, calls old S\psi
  !! routine s_psi_ . See comments in h_psi.f90 about band 
  !! parallelization.
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE xc_lib,           ONLY : exx_is_active
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_allgather, mp_size, &
                               mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: ik
  !! index of the k points
  INTEGER, INTENT(IN) :: id
  !! composite index of k, perturbation, etc.
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !! S matrix dot wavefunctions psi
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
  !
  CALL start_clock( 's_psi_bgrp' )
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     ! use band parallelization here
     ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
     CALL divide_all( inter_bgrp_comm,m,m_start,m_end, recv_counts,displs )
     CALL mp_type_create_column_section( spsi(1,1), 0, lda*npol, lda*npol, column_type )
     !
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL s_psi_ch_( lda, n, m_end-m_start+1, psi(1,m_start), spsi(1,m_start), ik, id )
     CALL mp_allgather( spsi, column_type, recv_counts, displs, inter_bgrp_comm )
     !
     CALL mp_type_free( column_type )
     DEALLOCATE( recv_counts )
     DEALLOCATE( displs )
  ELSE
     ! don't use band parallelization here
     CALL s_psi_ch_( lda, n, m, psi, spsi, ik, id )
  ENDIF
  !
  CALL stop_clock( 's_psi_bgrp' )
  !
  !
  RETURN
  !
END SUBROUTINE s_psi_ch
!
!
!----------------------------------------------------------------------------
SUBROUTINE s_psi_ch_( lda, n, m, psi, spsi, ik, id)
  !----------------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  USE kinds,            ONLY: DP
  USE becmod,           ONLY: becp
  USE uspp,             ONLY: vkb, nkb, okvan, qq_at, qq_so, ofsbeta
  USE uspp_param,       ONLY: upf, nh, nhm
  USE ions_base,        ONLY: nat, nsp, ityp
  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol, noncolin, lspinorb
  USE realus,           ONLY: real_space, invfft_orbital_gamma,     &
                              fwfft_orbital_gamma, calbec_rs_gamma, &
                              s_psir_gamma, invfft_orbital_k,       &
                              fwfft_orbital_k, calbec_rs_k, s_psir_k
  USE many_k_mod,       ONLY: vkbk_d, becpk_d
  USE many_k_ph_mod,    ONLY: current_ikb_ph, startkb_ph
  USE wavefunctions,    ONLY: psic
  USE fft_base,         ONLY: dffts
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  INTEGER, INTENT(IN) :: ik
  !! index of the k points
  INTEGER, INTENT(IN) :: id
  !! composite index of k, perturbation, etc.
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !! S matrix dot wavefunctions psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd, ik1
  !
  !
  ! ... initialize  spsi
  !
  ik1=ik-startkb_ph(current_ikb_ph)
  CALL threaded_memcpy( spsi, psi, lda*npol*m*2 )
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
  !
  CALL start_clock( 's_psi' )  
  !
  ! ... The product with the beta functions
  !
  IF ( gamma_only ) THEN
     !
     CALL errore('s_psi_ch_','gamma only not included',1)
     !
  ELSEIF ( noncolin ) THEN
     !
     CALL s_psi_ch_nc()
     !
  ELSE 
     !
     IF ( real_space ) THEN
        !
        DO ibnd = 1, m
!SdG: the becp are already computed ! no need to invfft psi to real space.
!           CALL invfft_orbital_k( psi, ibnd, m )
!SdG: we just need to clean psic in real space ...
           CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
!SdG: ... before computing the us-only contribution ...
           CALL s_psir_k( ibnd, m )
!SdG: ... and add it to spsi (already containing psi).
           CALL fwfft_orbital_k( spsi, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        !
     ELSE
        !
        CALL s_psi_ch_k()
        !
     ENDIF    
     !
  ENDIF    
  !
  CALL stop_clock( 's_psi' )
  !
  !
  RETURN
  !
  CONTAINS
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_ch_k()
       !-----------------------------------------------------------------------
       !! k-points version of \(\textrm{s_psi}\) routine.
       !
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
       ! counters
       COMPLEX(DP), ALLOCATABLE :: ps(:,:), qqc(:,:)
       ! ps = product vkb and psi ; qqc = complex version of qq
       !
       ALLOCATE( ps( nkb, m ), STAT=ierr )
       !
       !
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_k ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             ! qq is real:  copy it into a complex variable to perform
             ! a zgemm - simple but sub-optimal solution
             ALLOCATE( qqc(nh(nt),nh(nt)) )
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   qqc(:,:) = CMPLX ( qq_at(1:nh(nt),1:nh(nt),na), 0.0_DP, KIND=DP )
                   CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_DP,0.0_DP), &
                        qqc, nh(nt), becpk_d(ofsbeta(na)+1,1,1,id), nkb, &
                        (0.0_DP,0.0_DP), ps(ofsbeta(na)+1,1), nkb )
                ENDIF
             ENDDO
             DEALLOCATE( qqc )
             !
          ELSE
             !
             IF (nh(nt)>0) THEN
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      ps(ofsbeta(na)+1:ofsbeta(na)+nh(nt),1:m) = (0.0_DP,0.0_DP)
                   ENDIF
                ENDDO
             ENDIF
             !
          ENDIF
          !
       ENDDO
       !
       IF ( m == 1 ) THEN
          !
          CALL ZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkbk_d(:,&
             (ik1-1)*nkb+1:ik1*nkb), lda, ps, 1, ( 1.D0, 0.D0 ), spsi, 1 )
          !
       ELSE
          !
          CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkbk_d(:,&
             (ik1-1)*nkb+1:ik1*nkb), lda, ps, nkb, ( 1.D0, 0.D0 ), spsi, lda )
          !
       ENDIF
       !
       DEALLOCATE( ps )
       !
       !
       RETURN
       !
     END SUBROUTINE s_psi_ch_k     
     !
     !
     !-----------------------------------------------------------------------
       SUBROUTINE s_psi_ch_nc ( )
       !-----------------------------------------------------------------------
       !! k-points noncolinear/spinorbit version of \(\textrm{s_psi}\) routine.
       !
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ipol, ierr
       ! counters
       COMPLEX (DP), ALLOCATABLE :: ps(:,:,:)
       ! the product vkb and psi
       !
       !
       ALLOCATE( ps(nkb,npol,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_nc ', ' cannot allocate memory (ps) ', ABS(ierr) )
       ps(:,:,:) = (0.D0,0.D0)
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             !
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ofsbeta(na) + ih
                      DO jh = 1, nh(nt)
                         jkb = ofsbeta(na) + jh
                         IF ( .NOT. lspinorb ) THEN
                            DO ipol = 1, npol
                               DO ibnd = 1, m
                                  ps(ikb,ipol,ibnd) = ps(ikb,ipol,ibnd) + &
                                    qq_at(ih,jh,na)*becpk_d(jkb,ipol,ibnd,id)
                               ENDDO
                            ENDDO
                         ELSE
                            DO ibnd = 1, m
                               ps(ikb,1,ibnd) = ps(ikb,1,ibnd) + &
                                    qq_so(ih,jh,1,nt)*becpk_d(jkb,1,ibnd,id)+ &
                                    qq_so(ih,jh,2,nt)*becpk_d(jkb,2,ibnd,id)
                               ps(ikb,2,ibnd) = ps(ikb,2,ibnd) + &
                                    qq_so(ih,jh,3,nt)*becpk_d(jkb,1,ibnd,id)+ &
                                    qq_so(ih,jh,4,nt)*becpk_d(jkb,2,ibnd,id)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             !
          ENDIF
          !
       ENDDO
       !
       CALL ZGEMM ( 'N', 'N', n, m*npol, nkb, (1.d0,0.d0) , vkbk_d(:,&
              (ik1-1)*nkb+1:ik1*nkb), lda, ps, nkb, (1.d0,0.d0), &
              spsi(1,1), lda )
       !
       DEALLOCATE( ps )
       !
       !
       RETURN
       !
    END SUBROUTINE s_psi_ch_nc
    !
END SUBROUTINE s_psi_ch_
