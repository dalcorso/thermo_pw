!
! Copyright (C) 2002-2022 Quantum ESPRESSO group
! Copyright (C) 2023 Andrea Dal Corso (generalization to many k)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psii( lda, n, m, psi, hpsi, ik )
  !----------------------------------------------------------------------------
  !! This routine computes the product of the Hamiltonian matrix with m 
  !! wavefunctions contained in psi.
  !
  !! \(\textit{Wrapper routine}\): performs bgrp parallelization on 
  !! non-distributed bands. If suitable and required, calls old H\psi 
  !! routine h_psi_ .
  !
  USE kinds,              ONLY: DP
  USE noncollin_module,   ONLY: npol
  USE xc_lib,             ONLY: exx_is_active
  USE mp_bands,           ONLY: use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,                 ONLY: mp_allgather, mp_size, &
                                mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  INTEGER, INTENT(IN) :: ik
  !! k point
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m) 
  !! the wavefunction
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
  !
  !
  CALL start_clock( 'h_psi_bgrp' ); !write (*,*) 'start h_psi_bgrp'; FLUSH(6)
  !
  ! band parallelization with non-distributed bands is performed if
  ! 1. enabled (variable use_bgrp_in_hpsi must be set to .T.)
  ! 2. exact exchange is not active (if it is, band parallelization is already
  !    used in exx routines called by Hpsi)
  ! 3. there is more than one band, otherwise there is nothing to parallelize
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
     !
     ! use band parallelization here
     ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
     CALL divide_all( inter_bgrp_comm, m, m_start, m_end, recv_counts,displs )
     CALL mp_type_create_column_section( hpsi(1,1), 0, lda*npol, lda*npol, column_type )
     !
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL h_psii_( lda, n, m_end-m_start+1, psi(1,m_start), hpsi(1,m_start), ik )
     CALL mp_allgather( hpsi, column_type, recv_counts, displs, inter_bgrp_comm )
     !
     CALL mp_type_free( column_type )
     DEALLOCATE( recv_counts )
     DEALLOCATE( displs )
     !
  ELSE
     ! don't use band parallelization here
     CALL h_psii_( lda, n, m, psi, hpsi, ik )
     !
  ENDIF
  !
  CALL stop_clock( 'h_psi_bgrp' )
  !
  RETURN
  !
END SUBROUTINE h_psii
!
!----------------------------------------------------------------------------
SUBROUTINE h_psii_( lda, n, m, psi, hpsi, ik )
  !----------------------------------------------------------------------------
  !! This routine computes the product of the Hamiltonian matrix with m 
  !! wavefunctions contained in psi.
  !
  USE kinds,                   ONLY: DP
  USE bp,                      ONLY: lelfield, l3dstring, gdir, efield, efield_cry
  USE becmod,                  ONLY: bec_type, becp, calbec
  USE lsda_mod,                ONLY: current_spin, isk
  USE scf,                     ONLY: vrs  
  USE many_k_mod,              ONLY: g2kink_d, vkbk_d
  USE uspp,                    ONLY: nkb
  USE ldaU,                    ONLY: lda_plus_u, Hubbard_projectors
  USE gvect,                   ONLY: gstart
  USE control_flags,           ONLY: gamma_only, scissor
  USE noncollin_module,        ONLY: npol, noncolin
  USE realus,                  ONLY: real_space, invfft_orbital_gamma, fwfft_orbital_gamma, &
                                     calbec_rs_gamma, add_vuspsir_gamma, invfft_orbital_k,  &
                                     fwfft_orbital_k, calbec_rs_k, add_vuspsir_k,           & 
                                     v_loc_psir_inplace
  USE fft_base,                ONLY: dffts
  USE exx,                     ONLY: use_ace, vexx, vexxace_gamma, vexxace_k
  USE xc_lib,                  ONLY: exx_is_active, xclib_dft_is
  USE sci_mod,                 ONLY: p_psi
  USE fft_helper_subroutines
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  INTEGER, INTENT(IN) :: ik
  !! the k point
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m) 
  !! the wavefunction
  COMPLEX(DP), INTENT(OUT) :: hpsi(lda*npol,m)
  !! Hamiltonian dot psi 
  !
  ! ... local variables
  !
#if ! defined (__CUDA)
  INTEGER :: ipol, ibnd
  REAL(DP) :: ee
  !
  !
  CALL start_clock( 'h_psi' ); !write (*,*) 'start h_psi';FLUSH(6)

  !
  ! ... Here we set the kinetic energy (k+G)^2 psi and clean up garbage
  !
  !$omp parallel do
  DO ibnd = 1, m
     hpsi(1:n,ibnd) = g2kink_d(1:n,ik) * psi(1:n,ibnd)
     IF (n<lda) hpsi(n+1:lda, ibnd) = (0.0_dp, 0.0_dp)
     IF ( noncolin ) THEN
        hpsi(lda+1:lda+n, ibnd) = g2kink_d(1:n,ik) * psi(lda+1:lda+n, ibnd)
        IF (n<lda) hpsi(lda+n+1:lda+lda, ibnd) = (0.0_dp, 0.0_dp)
     ENDIF
  ENDDO
  !$omp end parallel do

  CALL start_clock( 'h_psi:pot' ); !write (*,*) 'start h_psi:pot';FLUSH(6)
  !
  ! ... Here the product with the local potential V_loc psi
  !
  IF ( gamma_only ) THEN
     ! 
     CALL errore('h_psii_', 'multiple k and gamma_only not available',1)
     IF ( real_space .AND. nkb > 0  ) THEN
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%has_task_groups ) &
             CALL errore( 'h_psi', 'task_groups not implemented with real_space', 1 )

        DO ibnd = 1, m, 2
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_gamma( psi, ibnd, m )
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock( 'h_psi:calbec' ) 
           CALL calbec_rs_gamma( ibnd, m, becp%r )
     CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace( ibnd, m ) 
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_gamma( ibnd, m )
           ! ... transform psic back in reciprocal space and add it to hpsi
           CALL fwfft_orbital_gamma( hpsi, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        !
     ELSE
        ! ... usual reciprocal-space algorithm
        CALL vloc_psi_gamma( lda, n, m, psi, vrs(1,current_spin), hpsi ) 
        !
     ENDIF 
     !
  ELSEIF ( noncolin ) THEN 
     !
     CALL vloc_psik_nc( lda, n, m, psi, vrs, hpsi, ik )
     !
  ELSE  
     ! 
     IF ( real_space .AND. nkb > 0  ) THEN
        CALL errore('h_psii_', 'multiple k and realus not available',1)
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%has_task_groups ) &
             CALL errore( 'h_psi', 'task_groups not implemented with real_space', 1 )
        !
        DO ibnd = 1, m
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_k( psi, ibnd, m )
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock( 'h_psi:calbec' )
           CALL calbec_rs_k( ibnd, m )
     CALL stop_clock( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace( ibnd, m )
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_k( ibnd, m )
           ! ... transform psic back in reciprocal space and add it to hpsi
           CALL fwfft_orbital_k( hpsi, ibnd, m, add_to_orbital=.TRUE. )
           !
        ENDDO
        !
     ELSE
        !
        CALL vloc_psik_k( lda, n, m, psi, vrs(1,isk(ik)), hpsi, ik )
        !
     ENDIF
     !
  ENDIF  
  !
  ! ... Here the product with the non local potential V_NL psi
  ! ... (not in the real-space case: it is done together with V_loc)
  !
  IF ( nkb > 0 .AND. .NOT. real_space) THEN
     !
     !
     CALL start_clock( 'h_psi:calbec' )
     CALL calbec( n, vkbk_d(:,nkb*(ik-1)+1:nkb*ik), psi, becp, m )
     CALL stop_clock( 'h_psi:calbec' )
     CALL add_vuspsik( lda, n, m, hpsi, ik )
     !
  ENDIF
  !
  CALL stop_clock( 'h_psi:pot' ); !write (*,*) 'stop h_psi:pot';FLUSH(6)
  !  
!  IF (xclib_dft_is('meta')) CALL h_psi_meta( lda, n, m, psi, hpsi )
  IF (xclib_dft_is('meta')) CALL errore('h_psii_',&
                              'multiple k and meta not available',1)
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u .AND. Hubbard_projectors.NE."pseudo" ) THEN
     !
     CALL errore('h_psii_', 'multiple k and hubbard U not available',1)
     IF ( noncolin ) THEN
        CALL vhpsi_nc( lda, n, m, psi, hpsi )
     ELSE
        CALL vhpsi( lda, n, m, psi, hpsi )
     ENDIF
     !
  ENDIF
  !
  ! ... apply scissor operator
  !
  IF (scissor) CALL errore('h_psii_', 'multiple k and scissor not available',1)
  IF (scissor) call p_psi(lda,n,m,psi,hpsi) 
  !
  ! ... Here the exact-exchange term Vxx psi
  !
  IF ( exx_is_active() ) THEN
     CALL errore('h_psii_', 'multiple k and exx not available',1)
     IF ( use_ace ) THEN
        IF ( gamma_only ) THEN
           CALL vexxace_gamma( lda, m, psi, ee, hpsi )
        ELSE
           CALL vexxace_k( lda, m, psi, ee, hpsi )
        ENDIF
     ELSE
        CALL vexx( lda, n, m, psi, hpsi, becp )
     ENDIF
  ENDIF
  !
  ! ... electric enthalpy if required
  !
  IF ( lelfield ) THEN
     CALL errore('h_psii_', 'multiple k and lelfield not available',1)
     !
     IF ( .NOT.l3dstring ) THEN
        CALL h_epsi_her_apply( lda, n, m, psi, hpsi,gdir, efield )
     ELSE
        DO ipol = 1, 3
           CALL h_epsi_her_apply( lda, n, m, psi, hpsi,ipol,efield_cry(ipol) )
        ENDDO
     ENDIF
     !
  ENDIF
#if defined(__OSCDFT)
  IF ( use_oscdft ) THEN
     CALL oscdft_h_psi(oscdft_ctx, lda, n, m, psi, hpsi)
  END IF
#endif
  !
  ! ... With Gamma-only trick, Im(H*psi)(G=0) = 0 by definition,
  ! ... but it is convenient to explicitly set it to 0 to prevent trouble
  !
  IF ( gamma_only .AND. gstart == 2 ) &
      hpsi(1,1:m) = CMPLX( DBLE( hpsi(1,1:m) ), 0.D0, KIND=DP)
  !
  CALL stop_clock( 'h_psi' )
#endif
  !
  !
  RETURN
  !
END SUBROUTINE h_psii_
