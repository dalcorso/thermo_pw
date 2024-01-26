!
! Copyright (C) 2002-2022 Quantum ESPRESSO group
! Copyright (C) 2023 Andrea Dal Corso (generalization to many k and 
!                                 global subroutines)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psii_gpu( lda, n, m, psi_d, hpsi_d, ik )
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
  !! the k points
  COMPLEX(DP), INTENT(IN) :: psi_d(lda*npol,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(OUT) :: hpsi_d(lda*npol,m)
  !! Hamiltonian dot psi
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
  !
  !
  CALL start_clock_gpu( 'h_psi_bgrp' ); !write (*,*) 'start h_psi_bgrp'; FLUSH(6)
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
     CALL mp_type_create_column_section( hpsi_d(1,1), 0, lda*npol, lda*npol, column_type )
     !
     ! Check if there at least one band in this band group
     IF (m_end >= m_start) &
        CALL h_psii__gpu( lda, n, m_end-m_start+1, psi_d(1,m_start), hpsi_d(1,m_start), ik )
     CALL mp_allgather( hpsi_d, column_type, recv_counts, displs, inter_bgrp_comm)
     !
     CALL mp_type_free( column_type )
     DEALLOCATE( recv_counts )
     DEALLOCATE( displs )
     !
  ELSE
     ! don't use band parallelization here
     CALL h_psii__gpu( lda, n, m, psi_d, hpsi_d, ik )
     !
  ENDIF
  !
  CALL stop_clock_gpu( 'h_psi_bgrp' )
  !
  !
  RETURN
  !
END SUBROUTINE h_psii_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE h_psii__gpu( lda, n, m, psi_d, hpsi_d, ik )
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !! This routine computes the product of the Hamiltonian matrix with m 
  !! wavefunctions contained in psi.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                   ONLY: DP
  USE bp,                      ONLY: lelfield, l3dstring, gdir, efield, efield_cry
  USE becmod,                  ONLY: bec_type, becp, calbec
  USE lsda_mod,                ONLY: current_spin, isk
  USE scf_gpum,                ONLY: vrs_d, using_vrs_d
  USE uspp,                    ONLY: nkb, vkb
  USE many_k_mod,              ONLY: vkbk_d
  USE ldaU,                    ONLY: lda_plus_u, lda_plus_u_kind, Hubbard_projectors
  USE gvect,                   ONLY: gstart
  USE control_flags,           ONLY: gamma_only
  USE noncollin_module,        ONLY: npol, noncolin
  USE realus,                  ONLY: real_space, invfft_orbital_gamma, fwfft_orbital_gamma, &
                                     calbec_rs_gamma, add_vuspsir_gamma, invfft_orbital_k,  &
                                     fwfft_orbital_k, calbec_rs_k, add_vuspsir_k,           & 
                                     v_loc_psir_inplace
  USE fft_base,                ONLY: dffts
  USE exx,                     ONLY: use_ace, vexx, vexxace_gamma, vexxace_k
  USE xc_lib,                  ONLY: exx_is_active, xclib_dft_is
  USE fft_helper_subroutines
  USE device_memcpy_m,         ONLY: dev_memcpy, dev_memset
  !
  USE many_k_mod,              ONLY: g2kink_d, becpk_d
#if defined(__OSCDFT)
  USE plugin_flags,            ONLY : use_oscdft
  USE oscdft_base,             ONLY : oscdft_ctx
  USE oscdft_functions_gpu,    ONLY : oscdft_h_psi_gpu
#endif
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)     :: lda, n, m, ik
  COMPLEX(DP), INTENT(IN)  :: psi_d(lda*npol,m)
  COMPLEX(DP), INTENT(OUT) :: hpsi_d(lda*npol,m)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d
#endif
  !
  COMPLEX(DP), ALLOCATABLE :: psi_host(:,:)
  COMPLEX(DP), ALLOCATABLE :: hpsi_host(:,:)
#if defined(__CUDA)
  attributes(PINNED) :: psi_host, hpsi_host
#endif
  !
  INTEGER     :: ipol, ibnd, incr, i, ispin, ll, lm
  REAL(dp)    :: ee
  !
  LOGICAL     :: need_host_copy
  COMPLEX(DP) :: adata, aux(nkb,m)
  !
  CALL start_clock_gpu( 'h_psi' ); !write (*,*) 'start h_psi';FLUSH(6)
!  CALL using_vrs_d(0)
  !
  ! ... Here we add the kinetic energy (k+G)^2 psi and clean up garbage
  !
  need_host_copy = ( real_space .and. nkb > 0  ) .OR. &
                     xclib_dft_is('meta') .OR. &
                    (lda_plus_u .AND. Hubbard_projectors.NE."pseudo" ) .OR. &
                    exx_is_active() .OR. lelfield


  IF (need_host_copy) THEN
      ALLOCATE(psi_host(lda*npol,m) , hpsi_host(lda*npol,m) )
      CALL dev_memcpy(psi_host, psi_d)    ! psi_host = psi_d
  ENDIF

!   WRITE(6,*) 'apply kinetic energy'
  !!$acc parallel loop collapse(2) present(g2kink_d, hpsi_d, psi_d)
!  DO ibnd = 1, m
!     DO i=1, lda
!        IF (i <= n) THEN
!           hpsi_d (i, ibnd) = g2kink_d (i,ik) * psi_d (i, ibnd)
!           IF ( noncolin ) THEN
!              hpsi_d (lda+i, ibnd) = g2kink_d (i,ik) * psi_d (lda+i, ibnd)
!           END IF
!        ELSE
!           hpsi_d (i, ibnd) = (0.0_dp, 0.0_dp)
!           IF ( noncolin ) THEN
!              hpsi_d (lda+i, ibnd) = (0.0_dp, 0.0_dp)
!           END IF
!        END IF
!     END DO
!  END DO
!  WRITE(6,*) 'done apply kinetic energy'


  IF (need_host_copy) CALL dev_memcpy(hpsi_host, hpsi_d)    ! hpsi_host = hpsi_d

  CALL start_clock_gpu( 'h_psi:pot' ); !write (*,*) 'start h_pot';FLUSH(6)
  !
  ! ... Here the product with the local potential V_loc psi
  !
  IF ( gamma_only ) THEN
     ! 
     CALL errore('h_psii__gpu','multiple k and gamma_only not available',1) 
     IF ( real_space .AND. nkb > 0  ) THEN
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        IF ( dffts%has_task_groups ) &
             CALL errore( 'h_psi', 'task_groups not implemented with real_space', 1 )

        DO ibnd = 1, m, 2
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_gamma(psi_host, ibnd, m )
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock_gpu( 'h_psi:calbec' ) 
           CALL calbec_rs_gamma( ibnd, m, becp%r )
           !$acc update device(becp%r)
     CALL stop_clock_gpu( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace( ibnd, m ) 
           ! ... psic (hpsi) -> psic + vusp
           CALL  add_vuspsir_gamma( ibnd, m )
           ! ... transform psic back in reciprocal space and add it to hpsi
           CALL fwfft_orbital_gamma( hpsi_host, ibnd, m, add_to_orbital=.TRUE. )
        ENDDO
        CALL dev_memcpy(hpsi_d, hpsi_host) ! hpsi_d = hpsi_host
        !
     ELSE
        ! ... usual reciprocal-space algorithm
        CALL vloc_psi_gamma_gpu ( lda, n, m, psi_d, vrs_d(1,current_spin), hpsi_d )
        !
     ENDIF 
     !
  ELSE IF ( noncolin ) THEN 
     !
     CALL errore('h_psii__gpu','multiple k and noncolin not available',1) 
     CALL vloc_psi_nc_gpu ( lda, n, m, psi_d, vrs_d, hpsi_d )
     !
  ELSE  
     ! 
     IF ( real_space .and. nkb > 0  ) then 
        !
        ! ... real-space algorithm
        ! ... fixme: real_space without beta functions does not make sense
        !
        CALL errore('h_psii__gpu','multiple k and realus not available',1) 
        IF ( dffts%has_task_groups ) &
             CALL errore( 'h_psi', 'task_groups not implemented with real_space', 1 )
        !
        DO ibnd = 1, m
           ! ... transform psi to real space -> psic 
           CALL invfft_orbital_k(psi_host, ibnd, m )
           ! ... compute becp%r = < beta|psi> from psic in real space
     CALL start_clock_gpu( 'h_psi:calbec' )
           CALL calbec_rs_k( ibnd, m )
     CALL stop_clock_gpu( 'h_psi:calbec' )
           ! ... psic -> vrs * psic (psic overwritten will become hpsi)
           CALL v_loc_psir_inplace( ibnd, m )
           ! ... psic (hpsi) -> psic + vusp
           CALL add_vuspsir_k( ibnd, m )
           ! ... transform psic back in reciprocal space and add it to hpsi
           CALL fwfft_orbital_k( hpsi_host, ibnd, m, add_to_orbital=.TRUE. )
           !
        ENDDO
        IF (need_host_copy) CALL dev_memcpy(hpsi_d, hpsi_host) ! hpsi_d = hpsi_host
        !
     ELSE
        !
!  WRITE(6,*) 'call vloc_psi_gpu'
CALL start_clock( 'cegt:vloc' )
        ispin=isk(ik)
        CALL vloc_psii_k_gpu ( lda, n, m, psi_d, vrs_d(1,ispin), hpsi_d, ik )
CALL stop_clock( 'cegt:vloc' )
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
!     CALL add_vuspsik_gpu( lda, n, m, hpsi_d, ik )
     !
  END IF
  !  
  CALL stop_clock_gpu( 'h_psi:pot' )
  !
  IF (xclib_dft_is('meta')) THEN
     CALL errore('h_psii__gpu','multiple k and metagga not available',1) 
     CALL dev_memcpy(hpsi_host, hpsi_d) ! hpsi_host = hpsi_d
     call h_psi_meta (lda, n, m, psi_host, hpsi_host)
     CALL dev_memcpy(hpsi_d, hpsi_host) ! hpsi_d = hpsi_host
  end if
  !
  ! ... Here we add the Hubbard potential times psi
  !
  IF ( lda_plus_u .AND. Hubbard_projectors.NE."pseudo" ) THEN
     !
     CALL errore('h_psii__gpu','multiple k and lda+U not available',1) 
     CALL dev_memcpy(hpsi_host, hpsi_d )    ! hpsi_host = hpsi_d
     IF ( noncolin ) THEN
        CALL vhpsi_nc( lda, n, m, psi_host, hpsi_host )
        CALL dev_memcpy(hpsi_d, hpsi_host)  ! hpsi_d = hpsi_host
     ELSE
        IF ( lda_plus_u_kind.EQ.0 .OR. lda_plus_u_kind.EQ.1 ) THEN
          CALL vhpsi_gpu( lda, n, m, psi_d, hpsi_d )  ! DFT+U
        ELSEIF ( lda_plus_u_kind.EQ.2 ) THEN          ! DFT+U+V
          CALL vhpsi( lda, n, m, psi_host, hpsi_host )
          CALL dev_memcpy(hpsi_d, hpsi_host)
        ENDIF
     ENDIF
     !
  ENDIF
  !
  ! ... Here the exact-exchange term Vxx psi
  !
  IF ( exx_is_active() ) THEN
     CALL errore('h_psii__gpu','multiple k and exx not available',1) 
     CALL dev_memcpy(hpsi_host, hpsi_d ) ! hpsi_host = hpsi_d
     IF ( use_ace) THEN
        IF (gamma_only) THEN
           CALL vexxace_gamma(lda,m,psi_host,ee,hpsi_host)
        ELSE
           CALL vexxace_k(lda,m,psi_host,ee,hpsi_host) 
        END IF
     ELSE
        CALL vexx( lda, n, m, psi_host, hpsi_host, becp )
     END IF
     CALL dev_memcpy(hpsi_d, hpsi_host) ! hpsi_d = hpsi_host
  END IF
  !
  ! ... electric enthalpy if required
  !
  IF ( lelfield ) THEN
     !
     CALL errore('h_psii__gpu','multiple k and electric field not available',1) 
     CALL dev_memcpy(hpsi_host, hpsi_d ) ! hpsi_host = hpsi_d
     IF ( .NOT.l3dstring ) THEN
        CALL h_epsi_her_apply( lda, n, m, psi_host, hpsi_host,gdir, efield )
     ELSE
        DO ipol=1,3
           CALL h_epsi_her_apply( lda, n, m, psi_host, hpsi_host,ipol,efield_cry(ipol) )
        END DO
     END IF
     CALL dev_memcpy(hpsi_d, hpsi_host) ! hpsi_d = hpsi_host
     !
  END IF
#if defined(__OSCDFT)
  IF ( use_oscdft ) THEN
     CALL oscdft_h_psi_gpu(oscdft_ctx, lda, n, m, psi_d, hpsi_d)
  END IF
#endif
  !
  ! ... With Gamma-only trick, Im(H*psi)(G=0) = 0 by definition,
  ! ... but it is convenient to explicitly set it to 0 to prevent trouble
  !
  IF ( gamma_only .AND. gstart == 2 ) then
      !$cuf kernel do(1)
      do i=1,m
         hpsi_d(1,i) = CMPLX( DBLE( hpsi_d(1,i) ), 0.D0 ,kind=DP)
      end do
  end if
  !
  if (need_host_copy) then
      DEALLOCATE(psi_host , hpsi_host )
  end if
  CALL stop_clock_gpu( 'h_psi' )
  !
  RETURN
  !
END SUBROUTINE h_psii__gpu
!
