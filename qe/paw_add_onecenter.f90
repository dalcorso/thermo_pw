!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains some PAW functionalities not available in
! the paw_one_center.f90 module. 
!
! Since some of the routines of paw_onecenter.f90 needed by the present
! routines are private to the module, they had to be copied here.
! These routines will propably disapper from the next versions.
! The copyright for these routines is:
!
! Copyright (C) 2007-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE paw_add_onecenter
    !
    USE kinds,          ONLY : DP
    USE paw_variables,  ONLY : paw_info, rad, radial_grad_style, vs_rad
    USE paw_onecenter,  ONLY : paw_rho_lm
    USE mp_images,      ONLY : nproc_image, me_image, intra_image_comm
    USE mp,             ONLY : mp_sum
    !
    IMPLICIT NONE

    ! entry points:
    PUBLIC :: PAW_deqtranpotential ! change of the xc potential due
                             ! to transverse magnetic field
    PUBLIC :: paw_xc_potential
    !
    INTEGER, SAVE :: paw_comm, me_paw, nproc_paw
    !
    INTEGER, SAVE :: nx_loc, ix_s, ix_e  ! parallelization on the directions
    !
    PRIVATE

    REAL(DP), ALLOCATABLE   :: msmall_lm(:,:,:) ! magnetiz. due to small
    !                                             components expanded on Y_lm
    REAL(DP), ALLOCATABLE :: g_lm(:,:,:)  ! potential density as lm components
    !
    LOGICAL :: with_small_so = .FALSE.
    !
    ! the following global variable controls the use of several fine-grained clocks
    ! set it to .false. in order to disable them, set it to .true. to enable them.
    !
    LOGICAL, PARAMETER :: TIMING = .false.
    !
    INTEGER, EXTERNAL :: ldim_block
    INTEGER, EXTERNAL :: gind_block

 CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE PAW_xc_potential( i, rho_lm, rho_core, v_lm, energy )
    !--------------------------------------------------------------------
    !! Use the density produced by sum_rad_rho to compute xc potential
    !! and energy, as xc functional is not diagonal on angular momentum
    !! numerical integration is performed.
    !
    USE noncollin_module,       ONLY : nspin_mag
    USE constants,              ONLY : e2, eps12
    USE uspp_param,             ONLY : upf
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE funct,                  ONLY : dft_is_gradient
    USE xc_lda_lsda,            ONLY : xc
    USE constants,              ONLY : fpi ! REMOVE
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,i%l**2,nspin)
    !! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(i%m)
    !! core charge, radial and spherical
    REAL(DP), INTENT(OUT) :: v_lm(i%m,i%l**2,nspin)
    !! potential density as lm components
    REAL(DP),OPTIONAL,INTENT(OUT) :: energy
    !! XC energy (if required)
    !
    ! ... local variables
    !
    REAL(DP), ALLOCATABLE :: rho_loc(:,:)       ! local density (workspace), up and down
    REAL(DP) :: v_rad(i%m,rad(i%t)%nx,nspin)    ! radial potential (to be integrated)
    REAL(DP), ALLOCATABLE :: g_rad(:,:,:)       ! radial potential

    REAL(DP), ALLOCATABLE :: rho_rad(:,:)       ! workspace (only one radial slice of rho)
    !
    REAL(DP), ALLOCATABLE :: e_rad(:)           ! aux, used to store radial slices of energy
    REAL(DP), ALLOCATABLE :: e_of_tid(:)        ! aux, for openmp parallel reduce
    REAL(DP) :: e                               ! aux, used to integrate energy
    !
    INTEGER :: ix,k                             ! counters on directions and radial grid
    INTEGER :: lsd                              ! switch for local spin density
    REAL(DP) :: vs, amag
    INTEGER :: kpol 
    INTEGER :: mytid, ntids
    !
    REAL(DP), ALLOCATABLE :: arho(:,:)
    REAL(DP), ALLOCATABLE :: ex(:), ec(:)
    REAL(DP), ALLOCATABLE :: vx(:,:), vc(:,:)
    REAL(DP), PARAMETER   :: eps = 1.e-30_dp
    !
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
    !
    IF (TIMING) CALL start_clock( 'PAW_xc_pot' )
    !
    ! true if using spin
    lsd = 0
    IF (nspin == 2)  lsd = 1
    IF (with_small_so) THEN
       ALLOCATE( g_rad(i%m,rad(i%t)%nx,nspin) )
       g_rad = 0.0_DP
    ENDIF
    !
!$omp parallel default(private), &
!$omp shared(i,rad,v_lm,rho_lm,rho_core,v_rad,ix_s,ix_e,energy,e_of_tid, &
!$omp        nspin,g,lsd,nspin_mag,with_small_so,g_rad)
#if defined(_OPENMP)
    mytid = omp_get_thread_num()+1 ! take the thread ID
    ntids = omp_get_num_threads()  ! take the number of threads
#else
    mytid = 1
    ntids = 1
#endif
    ! This will hold the "true" charge density, without r**2 or other factors
    ALLOCATE( rho_loc(i%m,nspin_mag) ) 
    rho_loc = 0._DP
    !
    ALLOCATE( rho_rad(i%m,nspin_mag) ) 
    !
    ALLOCATE( arho(i%m,2) ) !^^^
    ALLOCATE( ex(i%m) )
    ALLOCATE( ec(i%m) )
    ALLOCATE( vx(i%m,2) )
    ALLOCATE( vc(i%m,2) )
    !
    IF (PRESENT(energy)) THEN
!$omp single
       energy = 0._DP
       ALLOCATE( e_of_tid(ntids) )
!$omp end single
       ALLOCATE( e_rad(i%m) )
       e_of_tid(mytid) = 0._DP
    ENDIF
!$omp workshare
    v_rad = 0.0_dp
!$omp end workshare
!$omp do
    DO ix = ix_s, ix_e
        !
        ! --- LDA (and LSDA) part (no gradient correction) ---
        ! convert _lm density to real density along ix
        !
        CALL PAW_lm2rad( i, ix, rho_lm, rho_rad, nspin_mag )
        !
        ! compute the potential along ix
        !
        IF ( nspin_mag==4 ) THEN
           IF (with_small_so .AND. i%ae==1) CALL add_small_mag( i, ix, rho_rad )
           !
           DO k = 1, i%m
              rho_loc(k,1:nspin) = rho_rad(k,1:nspin)*g(i%t)%rm2(k)
              rho_loc(k,1) = rho_loc(k,1) + rho_core(k)
           ENDDO
           !
           CALL xc( i%m, 4, 2, rho_loc, ex, ec, vx, vc )
           !
           DO k = 1, i%m
              IF (PRESENT(energy)) &
                  e_rad(k) = e2*(ex(k)+ec(k))*(rho_rad(k,1)+rho_core(k)*g(i%t)%r2(k))
              vs = e2*0.5D0*( vx(k,1) + vc(k,1) - vx(k,2) - vc(k,2) )
              v_rad(k,ix,1) = e2*(0.5D0*( vx(k,1) + vc(k,1) + vx(k,2) + vc(k,2)))
              amag = SQRT(rho_loc(k,2)**2+rho_loc(k,3)**2+rho_loc(k,4)**2)
              IF ( amag > eps12 ) THEN
                 v_rad(k,ix,2:4) =  vs * rho_loc(k,2:4) / amag
              ELSE
                 v_rad(k,ix,2:4)=0.0_DP
                 IF (PRESENT(energy)) e_rad(k)=0.0_DP
              ENDIF
           ENDDO
           !
           IF ( with_small_so ) CALL compute_g( i, ix, v_rad, g_rad )
        ELSEIF ( nspin==2 ) THEN
           DO k = 1, i%m
              rho_loc(k,1) = rho_rad(k,1)*g(i%t)%rm2(k)
              rho_loc(k,2) = rho_rad(k,2)*g(i%t)%rm2(k)
           ENDDO
        ELSE
           DO k = 1, i%m
              rho_loc(k,1) = rho_rad(k,1)*g(i%t)%rm2(k)
           ENDDO
        ENDIF
        !
        ! Integrate to obtain the energy
        !
        IF (nspin_mag <= 2 ) THEN
           !
           !
           IF ( lsd == 0 ) THEN
             !
             arho(:,1) = rho_loc(:,1) + rho_core
             !
             CALL xc( i%m, 1, 1, arho(:,1), ex, ec, vx(:,1), vc(:,1) )
             !
             v_rad(:,ix,1) = e2*( vx(:,1) + vc(:,1) )
             IF (PRESENT(energy)) e_rad = e2*( ex(:) + ec(:) )
             !
           ELSE
             !
             arho(:,1) = rho_loc(:,1) + rho_loc(:,2) + rho_core(:)
             arho(:,2) = rho_loc(:,1) - rho_loc(:,2)
             !
             CALL xc( i%m, 2, 2, arho, ex, ec, vx, vc )
             !
             v_rad(:,ix,:) = e2*( vx(:,:) + vc(:,:) )
             IF (PRESENT(energy)) e_rad(:) = e2*( ex(:) + ec(:) )
             !
           ENDIF
           !
           IF (PRESENT(energy)) THEN
              IF (nspin_mag < 2) THEN
                 e_rad = e_rad * ( rho_rad(:,1) + rho_core*g(i%t)%r2 )
              ELSEIF (nspin_mag == 2) THEN
                 e_rad = e_rad *(rho_rad(:,1)+rho_rad(:,2)+rho_core*g(i%t)%r2 )
              ENDIF
           ENDIF
           !
        ENDIF
        ! Integrate to obtain the energy
        IF (PRESENT(energy)) THEN
           CALL simpson( i%m, e_rad, g(i%t)%rab, e )
           e_of_tid(mytid) = e_of_tid(mytid) + e * rad(i%t)%ww(ix)
        ENDIF
        !
    ENDDO
!$omp end do nowait
    !
    IF ( PRESENT(energy) ) DEALLOCATE( e_rad )
    !
    DEALLOCATE( rho_rad )
    DEALLOCATE( rho_loc )
    !
    DEALLOCATE( arho )
    DEALLOCATE( ex )
    DEALLOCATE( ec )
    DEALLOCATE( vx )
    DEALLOCATE( vc )
    !
!$omp end parallel
    !
    !
    IF (PRESENT(energy)) THEN
       energy = SUM(e_of_tid)
       DEALLOCATE( e_of_tid )
       CALL mp_sum( energy, paw_comm )
    ENDIF
    !
    ! Recompose the sph. harm. expansion
    CALL PAW_rad2lm( i, v_rad, v_lm, i%l, nspin_mag )
    !
    IF ( with_small_so ) THEN
       CALL PAW_rad2lm( i, g_rad, g_lm, i%l, nspin_mag )
       DEALLOCATE( g_rad )
    ENDIF
    !
    ! Add gradient correction, if necessary
    IF ( dft_is_gradient() ) &
        CALL PAW_gcxc_potential( i, rho_lm, rho_core, v_lm, energy )
        !
    IF (TIMING) CALL stop_clock( 'PAW_xc_pot' )
    !
    RETURN
    !
  END SUBROUTINE PAW_xc_potential
  !
  !
  !------------------------------------------------------------------------------------
  SUBROUTINE PAW_gcxc_potential( i, rho_lm, rho_core, v_lm, energy )
    !---------------------------------------------------------------------------------
    !! Add gradient correction to v_xc, code mostly adapted from ../atomic/vxcgc.f90
    !! in order to support non-spherical charges (as Y_lm expansion).  
    !! Note that the first derivative in vxcgc becomes a gradient, while the second is
    !! a divergence.  
    !! We also have to temporarily store some additional Y_lm components in order not
    !! to loose precision during the calculation, even if only the ones up to 
    !! lmax_rho (the maximum in the density of charge) matter when computing \int v*rho.
    !
    USE lsda_mod,               ONLY : nspin
    USE noncollin_module,       ONLY : noncolin, nspin_mag, nspin_gga
    USE atom,                   ONLY : g => rgrid
    USE constants,              ONLY : sqrtpi, fpi,pi,e2
    USE funct,                  ONLY : igcc_is_lyp
    USE xc_gga,                 ONLY : xc_gcx
    USE mp,                     ONLY : mp_sum
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    REAL(DP), INTENT(IN) :: rho_lm(i%m,i%l**2,nspin)
    !! charge density as lm components
    REAL(DP), INTENT(IN) :: rho_core(i%m)
    !! core charge, radial and spherical
    REAL(DP), INTENT(INOUT) :: v_lm(i%m,i%l**2,nspin)
    !! potential to be updated
    REAL(DP), OPTIONAL, INTENT(INOUT) :: energy
    !! if present, add GC to energy
    !
    ! ... local variables
    !
    REAL(DP), PARAMETER :: epsr = 1.e-6_DP, epsg = 1.e-10_DP
    ! (as in PW/src/gradcorr.f90)
    !
    REAL(DP), ALLOCATABLE :: rho_rad(:,:) ! charge density sampled
    REAL(DP), ALLOCATABLE :: grad(:,:,:)  ! gradient
    REAL(DP), ALLOCATABLE :: gradx(:,:,:) ! gradient (swapped indexes)
    REAL(DP), ALLOCATABLE :: grad2(:,:)   ! square modulus of gradient
                                          ! (first of charge, than of hamiltonian)
    REAL(DP), ALLOCATABLE :: gc_rad(:,:,:)    ! GC correction to V (radial samples)
    REAL(DP), ALLOCATABLE :: gc_lm(:,:,:)     ! GC correction to V (Y_lm expansion)
    REAL(DP), ALLOCATABLE :: h_rad(:,:,:,:)   ! hamiltonian (vector field)
    REAL(DP), ALLOCATABLE :: h_lm(:,:,:,:)    ! hamiltonian (vector field)
                                    ! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !
    REAL(DP), ALLOCATABLE :: div_h(:,:,:)  ! div(hamiltonian)
    !
    REAL(DP), ALLOCATABLE :: rhoout_lm(:,:,:) ! charge density as lm components
    REAL(DP), ALLOCATABLE :: vout_lm(:,:,:)   ! potential as lm components
    REAL(DP), ALLOCATABLE :: segni_rad(:,:)   ! sign of the magnetization
    !
    REAL(DP), ALLOCATABLE :: arho(:,:), grad2_v(:)
    REAL(DP), ALLOCATABLE :: r_vec(:,:)
    !
    REAL(DP), DIMENSION(i%m,nspin_gga) :: v1x, v2x, v1c, v2c  !workspace
    REAL(DP), DIMENSION(i%m) :: sx, sc
    REAL(DP), ALLOCATABLE :: v2cud(:)
    !
    REAL(DP) :: vnull
    !
    REAL(DP), ALLOCATABLE :: e_rad(:)      ! aux, used to store energy
    REAL(DP) :: e, e_gcxc                  ! aux, used to integrate energy
    !
    INTEGER  :: k, ix, is, lm              ! counters on spin and mesh
    REAL(DP) :: sgn                        ! workspace
    REAL(DP) :: co2                        ! workspace
    !
    INTEGER :: mytid, ntids
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
    REAL(DP),ALLOCATABLE :: egcxc_of_tid(:)


    if(TIMING) CALL start_clock ('PAW_gcxc_v')
  
    e_gcxc = 0._dp

    ALLOCATE( gc_rad(i%m,rad(i%t)%nx,nspin_gga) )! GC correction to V (radial samples)
    ALLOCATE( gc_lm(i%m,i%l**2,nspin_gga)       )! GC correction to V (Y_lm expansion)
    ALLOCATE( h_rad(i%m,3,rad(i%t)%nx,nspin_gga))! hamiltonian (vector field)
    ALLOCATE( h_lm(i%m,3,(i%l+rad(i%t)%ladd)**2,nspin_gga) ) 
                                        ! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !
    ALLOCATE(div_h(i%m,i%l**2,nspin_gga))
    ALLOCATE(rhoout_lm(i%m,i%l**2,nspin_gga)) ! charge density as lm components
    ALLOCATE(vout_lm(i%m,i%l**2,nspin_gga))   ! potential as lm components
    ALLOCATE(segni_rad(i%m,rad(i%t)%nx))      ! charge density as lm components
    vout_lm=0.0_DP
    !
    IF ( nspin_mag == 2 .OR. nspin_mag == 4 ) THEN
       !   transform the noncollinear case into sigma-GGA case
       IF (noncolin) THEN
          CALL compute_rho_spin_lm(i, rho_lm, rhoout_lm, segni_rad)
       ELSE
          rhoout_lm=rho_lm
       ENDIF
    ENDIF
    !
!$omp parallel default(private), &
!$omp shared(i,g,nspin,nspin_gga,nspin_mag,rad,e_gcxc,egcxc_of_tid,gc_rad,h_rad,rho_lm,rhoout_lm,rho_core,energy,ix_s,ix_e)
    !
    mytid = 1
    ntids = 1
#if defined(_OPENMP)
    mytid = omp_get_thread_num()+1 ! take the thread ID
    ntids = omp_get_num_threads()  ! take the number of threads
#endif
    ALLOCATE( rho_rad(i%m,nspin_gga)) ! charge density sampled
    ALLOCATE( grad(i%m,3,nspin_gga) ) ! gradient
    ALLOCATE( grad2(i%m,nspin_gga)  ) ! square modulus of gradient
                                      ! (first of charge, than of hamiltonian)
!$omp workshare
    gc_rad = 0.0d0
    h_rad  = 0.0d0
!$omp end workshare nowait
    !
    IF (PRESENT(energy)) THEN
!$omp single
        ALLOCATE( egcxc_of_tid(ntids) )
!$omp end single
        egcxc_of_tid(mytid) = 0.0_dp
        ALLOCATE( e_rad(i%m) )
    ENDIF
    !
    spin:&
    !
    IF ( nspin_mag == 1 ) THEN
        !
        !     GGA case
        !
        ALLOCATE( arho(i%m,1), grad2_v(i%m) )
        ALLOCATE( gradx(3,i%m,1) )
        !
!$omp do
        DO ix = ix_s, ix_e
           !
           !  WARNING: the next 2 calls are duplicated for spin==2
           CALL PAW_lm2rad( i, ix, rho_lm, rho_rad, nspin_mag )
           CALL PAW_gradient( i, ix, rho_lm, rho_rad, rho_core, grad2, grad )
           !
           DO k = 1, i%m
              arho(k,1) = rho_rad(k,1)*g(i%t)%rm2(k) + rho_core(k)
              arho(k,1) = ABS(arho(k,1))
              gradx(:,k,1) = grad(k,:,1)
           ENDDO
           !
           CALL xc_gcx( i%m, 1, arho, gradx, sx, sc, v1x, v2x, v1c, v2c )
           !
           DO k = 1, i%m
              IF ( PRESENT(energy) ) &
                 e_rad(k)     = e2 * (sx(k)+sc(k)) * g(i%t)%r2(k)
              gc_rad(k,ix,1)  = (v1x(k,1)+v1c(k,1))  !*g(i%t)%rm2(k)
              h_rad(k,:,ix,1) = (v2x(k,1)+v2c(k,1))*grad(k,:,1)*g(i%t)%r2(k)
           ENDDO
           !
           ! integrate energy (if required)
           IF ( PRESENT(energy) ) THEN
               CALL simpson(i%m, e_rad, g(i%t)%rab, e)
               egcxc_of_tid(mytid) = egcxc_of_tid(mytid) + e*rad(i%t)%ww(ix)
           ENDIF
           !
        ENDDO
!$omp end do
        !
        DEALLOCATE( arho, grad2_v ) 
        DEALLOCATE( gradx )
        !
        !
    ELSEIF ( nspin_mag == 2 .OR. nspin_mag == 4 ) THEN
        !
        ALLOCATE( gradx(3,i%m,2) )
        ALLOCATE( r_vec(i%m,2) )
        ALLOCATE( v2cud(i%m) )
        !
        !   this is the \sigma-GGA case
        !
!$omp do
        DO ix = ix_s, ix_e
           !
           CALL PAW_lm2rad( i, ix, rhoout_lm, rho_rad, nspin_gga )
           CALL PAW_gradient( i, ix, rhoout_lm, rho_rad, rho_core,grad2, grad )
           !
           DO k = 1, i%m
               !
               ! Prepare the necessary quantities
               ! rho_core is considered half spin up and half spin down:
               co2 = rho_core(k)/2
               ! than I build the real charge dividing by r**2
               r_vec(k,1) = rho_rad(k,1)*g(i%t)%rm2(k) + co2
               r_vec(k,2) = rho_rad(k,2)*g(i%t)%rm2(k) + co2
               !
               !
               gradx(:,k,1) = grad(k,:,1)
               gradx(:,k,2) = grad(k,:,2)
           ENDDO
           !
           CALL xc_gcx( i%m, 2, r_vec, gradx, sx, sc, v1x, v2x, v1c, v2c, v2cud )
           !
           DO k = 1, i%m
              !
              IF ( PRESENT(energy) ) e_rad(k) = e2*(sx(k)+sc(k))*g(i%t)%r2(k)
              !
              ! first term of the gradient correction : D(rho*Exc)/D(rho)
              gc_rad(k,ix,1)  = (v1x(k,1)+v1c(k,1)) !*g(i%t)%rm2(k)
              gc_rad(k,ix,2)  = (v1x(k,2)+v1c(k,2)) !*g(i%t)%rm2(k)
              !
              ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
              ! h_rad(k,:,ix,1) =( (v2xup_vec(k)+v2c)*grad(k,:,1)+v2c*grad(k,:,2) )*g(i%t)%r2(k)
              ! h_rad(k,:,ix,2) =( (v2xdw_vec(k)+v2c)*grad(k,:,2)+v2c*grad(k,:,1) )*g(i%t)%r2(k)
              h_rad(k,:,ix,1) =( (v2x(k,1)+v2c(k,1))*grad(k,:,1) + &
                                  v2cud(k)*grad(k,:,2) )*g(i%t)%r2(k)
              h_rad(k,:,ix,2) =( (v2x(k,2)+v2c(k,2))*grad(k,:,2) + &
                                  v2cud(k)*grad(k,:,1) )*g(i%t)%r2(k)
              !
           ENDDO
           !
           ! integrate energy (if required)
           ! NOTE: this integration is duplicated for every spin, FIXME!
           IF (PRESENT(energy)) THEN
               CALL simpson( i%m, e_rad, g(i%t)%rab, e )
               egcxc_of_tid(mytid) = egcxc_of_tid(mytid) + e * rad(i%t)%ww(ix)
           ENDIF
           !
        ENDDO ! ix
!$omp end do nowait
        !
        DEALLOCATE( gradx )
        DEALLOCATE( r_vec )
        DEALLOCATE( v2cud )
        !
    ELSE spin
    !
!$omp master
        CALL errore( 'PAW_gcxc_v', 'unknown spin number', 2 )
!$omp end master
    ENDIF spin
    !
    IF ( PRESENT(energy) ) THEN
       DEALLOCATE( e_rad )
    ENDIF
    !
    DEALLOCATE( rho_rad )
    DEALLOCATE( grad  )
    DEALLOCATE( grad2 )
!$omp end parallel
    !
    IF ( PRESENT(energy) ) THEN
       e_gcxc = SUM(egcxc_of_tid)
       CALL mp_sum( e_gcxc, paw_comm )
       energy = energy + e_gcxc
    ENDIF
    !
    IF ( PRESENT(energy) ) THEN
       DEALLOCATE( egcxc_of_tid )
    ENDIF
    !
    ! convert the first part of the GC correction back to spherical harmonics
    CALL PAW_rad2lm( i, gc_rad, gc_lm, i%l, nspin_gga )
    !
    ! Note that the expansion into spherical harmonics of the derivative 
    ! with respect to theta of the spherical harmonics, is very slow to
    ! converge and would require a huge angular momentum ladd.
    ! This derivative divided by sin_th is much faster to converge, so
    ! we divide here before calculating h_lm and keep into account for
    ! this factor sin_th in the expression of the divergence.
    !
    ! ADC 30/04/2009.
    ! 
    DO ix = ix_s, ix_e
       h_rad(1:i%m,3,ix,1:nspin_gga) = h_rad(1:i%m,3,ix,1:nspin_gga) / &
                                       rad(i%t)%sin_th(ix)
    ENDDO
    ! We need the gradient of H to calculate the last part of the exchange
    ! and correlation potential. First we have to convert H to its Y_lm expansion
    CALL PAW_rad2lm3( i, h_rad, h_lm, i%l+rad(i%t)%ladd, nspin_gga )
    !
    ! Compute div(H)
    CALL PAW_divergence( i, h_lm, div_h, i%l+rad(i%t)%ladd, i%l )
    !                       input max lm --^  output max lm-^
    !
    ! Finally sum it back into v_xc
    DO is = 1,nspin_gga
      DO lm = 1,i%l**2
         vout_lm(1:i%m,lm,is) = vout_lm(1:i%m,lm,is) + e2*(gc_lm(1:i%m,lm,is)-div_h(1:i%m,lm,is))
      ENDDO
    ENDDO
    !
    IF (nspin_mag == 4 ) THEN
       CALL compute_pot_nonc( i, vout_lm, v_lm, segni_rad, rho_lm )
    ELSE
       v_lm(:,:,1:nspin_mag) = v_lm(:,:,1:nspin_mag)+vout_lm(:,:,1:nspin_mag)
    ENDIF
    !
    DEALLOCATE( gc_rad )
    DEALLOCATE( gc_lm  )
    DEALLOCATE( h_rad  )
    DEALLOCATE( h_lm   )
    DEALLOCATE( div_h  )
    DEALLOCATE( rhoout_lm )
    DEALLOCATE( vout_lm   )
    DEALLOCATE( segni_rad )
    !
    ! if(PRESENT(energy)) write(*,*) "gcxc -->", e_gcxc
    IF (TIMING) CALL stop_clock( 'PAW_gcxc_v' )
    !
  END SUBROUTINE PAW_gcxc_potential
  !
  !
  !-------------------------------------------------------------------------------
  SUBROUTINE PAW_divergence( i, F_lm, div_F_lm, lmaxq_in, lmaxq_out )
    !---------------------------------------------------------------------------
    !! Compute divergence of a vector field (actually the hamiltonian).  
    !! It is assumed that:  
    !! 1. the input function is multiplied by \(r^2\);  
    !! 2. the output function is multiplied by \(r^2\) too.
    !
    USE constants,              ONLY : sqrtpi, fpi, e2
    USE noncollin_module,       ONLY : nspin_gga
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    INTEGER, INTENT(IN) :: lmaxq_in
    !! max angular momentum to derive (divergence is reliable up to lmaxq_in-2)
    INTEGER, INTENT(IN) :: lmaxq_out
    !! max angular momentum to reconstruct for output
    REAL(DP), INTENT(IN) :: F_lm(i%m,3,lmaxq_in**2,nspin_gga)
    !! Y_lm expansion of F
    REAL(DP), INTENT(OUT):: div_F_lm(i%m,lmaxq_out**2,nspin_gga)
    !! div(F) 
    !
    ! ... local variables
    !
    REAL(DP) :: div_F_rad(i%m,rad(i%t)%nx,nspin_gga) ! div(F) on rad. grid
    REAL(DP) :: aux(i%m)!,aux2(i%m)                  ! workspace
    ! counters on: spin, angular momentum, radial grid point:
    INTEGER :: is, lm, ix
    !
    IF (TIMING) CALL start_clock( 'PAW_div' )
    !
    ! This is the divergence in spherical coordinates:
    !     {1 \over r^2}{\partial ( r^2 A_r ) \over \partial r} 
    !   + {1 \over r\sin\theta}{\partial \over \partial \theta} (  A_\theta\sin\theta )
    !   + {1 \over r\sin\theta}{\partial A_\phi \over \partial \phi}
    !
    ! The derivative sum_LM d(Y_LM sin(theta) )/dtheta will be expanded as:
    ! sum_LM ( Y_lm cos(theta) + sin(theta) dY_lm/dtheta )
    !
    ! The radial component of the divergence is computed last, for practical reasons
    !
    !     CALL errore('PAW_divergence', 'More angular momentum components are needed (in input)'//&
    !                 ' to provide the number you have requested (in output)', lmaxq_out-lmaxq_in+2)
    !
    ! phi component
    !
    div_F_rad = 0.0_DP
    !
    DO is = 1, nspin_gga
      DO ix = ix_s, ix_e
        aux(:) = 0._DP
        ! this derivative has no spherical component, so lm starts from 2
        DO lm = 2, lmaxq_in**2
            aux(1:i%m) = aux(1:i%m) + rad(i%t)%dylmp(ix,lm)* (F_lm(1:i%m,2,lm,is))! &
                                    !* g(i%t)%rm1(1:i%m) !/sin_th(ix) 
        ! as for PAW_gradient this is already present in dylmp --^
        ENDDO
        div_F_rad(1:i%m,ix,is) = aux(1:i%m)
      ENDDO
    ENDDO
    !
    ! theta component
    DO is = 1, nspin_gga
      DO ix = ix_s, ix_e
        aux(:) = 0._DP
        ! this derivative has a spherical component too!
        DO lm = 1, lmaxq_in**2
           aux(1:i%m) = aux(1:i%m) + F_lm(1:i%m,3,lm,is) &
                        * (rad(i%t)%dylmt(ix,lm)*rad(i%t)%sin_th(ix)&
                        + 2.0_DP*rad(i%t)%ylm(ix,lm)*rad(i%t)%cos_th(ix))
                      ! * (rad(i%t)%dylmt(ix,lm)  & 
                      ! + rad(i%t)%ylm(ix,lm)*rad(i%t)%cotg_th(ix) )
        ENDDO
        div_F_rad(1:i%m,ix,is) = div_F_rad(1:i%m,ix,is) + aux(1:i%m)
      ENDDO
    ENDDO
    !
    ! Convert what I have done so far to Y_lm
    CALL PAW_rad2lm( i, div_F_rad, div_F_lm, lmaxq_out, nspin_gga )
    ! Multiply by 1/r**3: 1/r is for theta and phi componente only
    ! 1/r**2 is common to all the three components.
    DO is = 1, nspin_gga
      DO lm = 1, lmaxq_out**2
        div_F_lm(1:i%m,lm,is) = div_F_lm(1:i%m,lm,is) * g(i%t)%rm3(1:i%m)
      ENDDO
    ENDDO
    !
    ! Compute partial radial derivative d/dr
    DO is = 1, nspin_gga
      DO lm = 1, lmaxq_out**2
        ! Derive along \hat{r} (F already contains a r**2 factor, otherwise
        ! it may be better to expand (1/r**2) d(A*r**2)/dr = dA/dr + 2A/r)
        CALL radial_gradient( F_lm(1:i%m,1,lm,is), aux, g(i%t)%r, i%m, radial_grad_style )
        ! Sum it in the divergence: it is already in the right Y_lm form
        aux(1:i%m) = aux(1:i%m)*g(i%t)%rm2(1:i%m)
        !
        div_F_lm(1:i%m,lm,is) = div_F_lm(1:i%m,lm,is) + aux(1:i%m)
      ENDDO
    ENDDO
    !
    IF (TIMING) CALL stop_clock( 'PAW_div' )
    !
  END SUBROUTINE PAW_divergence
  !
  !
  !---------------------------------------------------------------------------------------
  SUBROUTINE PAW_gradient( i, ix, rho_lm, rho_rad, rho_core, grho_rad2, grho_rad )
    !--------------------------------------------------------------------------------------
    !! Build gradient of radial charge distribution from its spherical harmonics expansion
    !
    USE constants,              ONLY : fpi
    USE noncollin_module,       ONLY : nspin_gga
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    !
    INTEGER, INTENT(IN)  :: ix
    !! line of the dylm2 matrix to use actually it is  one of the nx spherical
    !! integration directions
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    REAL(DP), INTENT(IN) :: rho_lm(i%m,i%l**2,nspin_gga)
    !! Y_lm expansion of rho
    REAL(DP), INTENT(IN) :: rho_rad(i%m,nspin_gga)
    !! radial density along direction ix
    REAL(DP), INTENT(IN) :: rho_core(i%m)
    !! core density
    REAL(DP), INTENT(OUT) :: grho_rad2(i%m,nspin_gga)
    !! |grad(rho)|^2 on rad. grid
    REAL(DP), OPTIONAL,INTENT(OUT):: grho_rad(i%m,3,nspin_gga)
    !! vector gradient (only for gcxc)
    ! r, theta and phi components ---^
    !
    REAL(DP) :: aux(i%m), aux2(i%m), fact  ! workspace
    INTEGER :: is, lm                      ! counters on: spin, angular momentum
    !
    IF (TIMING) CALL start_clock( 'PAW_grad' )
    !
    ! 1. build real charge density = rho/r**2 + rho_core
    ! 2. compute the partial derivative of rho_rad
    fact = 1.0_DP/DBLE(nspin_gga)
    grho_rad2(:,:) = 0._DP
    !
    DO is = 1, nspin_gga
        ! build real charge density
        aux(1:i%m) = rho_rad(1:i%m,is)*g(i%t)%rm2(1:i%m) &
                          + rho_core(1:i%m)*fact
        CALL radial_gradient( aux, aux2, g(i%t)%r, i%m, radial_grad_style )
        ! compute the square
        grho_rad2(:,is) = aux2(:)**2
        ! store in vector gradient, if present:
        IF (PRESENT(grho_rad)) grho_rad(:,1,is) = aux2(:)
    ENDDO
    !
    spin: &
    DO is = 1, nspin_gga
        aux(:)  = 0._DP
        aux2(:) = 0._DP
        ! Spherical (lm=1) component (that would also include core correction) can be omitted
        ! as its contribution to non-radial derivative is zero
        DO lm = 2, i%l**2
            ! 5. [ \sum_{lm} rho(r) (dY_{lm}/dphi /cos(theta))  ]**2
            aux(1:i%m) = aux(1:i%m) + rad(i%t)%dylmp(ix,lm) * rho_lm(1:i%m,lm,is)
            ! 6. [ \sum_{lm} rho(r) (dY_{lm}/dtheta)  ]**2
            aux2(1:i%m) = aux2(1:i%m) + rad(i%t)%dylmt(ix,lm) * rho_lm(1:i%m,lm,is)
        ENDDO
        ! Square and sum up these 2 components, the (1/r**2)**3 factor come from:
        !  a. 1/r**2 from the derivative in spherical coordinates
        !  b. (1/r**2)**2 from rho_lm being multiplied by r**2 
        !     (as the derivative is orthogonal to r you can multiply after deriving)
        grho_rad2(1:i%m,is) = grho_rad2(1:i%m,is) &
                                + (aux(1:i%m)**2 + aux2(1:i%m)**2) &
                                    * g(i%t)%rm2(1:i%m)**3
        ! Store vector components:
        IF (PRESENT(grho_rad)) THEN
            grho_rad(1:i%m,2,is) = aux(1:i%m)  * g(i%t)%rm3(1:i%m) ! phi 
            grho_rad(1:i%m,3,is) = aux2(1:i%m) * g(i%t)%rm3(1:i%m) ! theta
        ENDIF
    ENDDO spin
    !
    IF (TIMING) CALL stop_clock( 'PAW_grad' )
    !
  END SUBROUTINE PAW_gradient
  !
  !-----------------------------------------------------------------------------------
  SUBROUTINE PAW_lm2rad( i, ix, F_lm, F_rad, nspin )
    !---------------------------------------------------------------------------------
    !! Build radial charge distribution from its spherical harmonics expansion.
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    INTEGER :: ix
    !! line of the ylm matrix to use
    !! actually it is one of the nx directions
    INTEGER, INTENT(IN) :: nspin
    !! number of spin components
    REAL(DP), INTENT(IN) :: F_lm(i%m,i%l**2,nspin)
    !! Y_lm expansion of rho
    REAL(DP), INTENT(OUT) :: F_rad(i%m,nspin)
    !! charge density on rad. grid
    !
    ! ... local variables
    !
    INTEGER :: ispin, lm ! counters on angmom and spin
    !
    IF (TIMING) CALL start_clock( 'PAW_lm2rad' )
    !
    F_rad(:,:) = 0._DP
    ! cycling on spin is a bit less general...
    spins: DO ispin = 1,nspin
        DO lm = 1, i%l**2
            F_rad(:,ispin) = F_rad(:,ispin) + &
                    rad(i%t)%ylm(ix,lm)*F_lm(:,lm,ispin)
        ENDDO ! lm
    ENDDO spins
    !
    IF (TIMING) CALL stop_clock( 'PAW_lm2rad' )
    !
  END SUBROUTINE PAW_lm2rad
  !
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE PAW_rad2lm( i, F_rad, F_lm, lmax_loc, nspin )
    !------------------------------------------------------------------------------
    !! Computes:
    !! \[ F_{lm}(r) = \int d \Omega\ F(r,\text{th},\text{ph})\ Y_{lm}(\text{th},
    !! \text{ph}) \]
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    INTEGER, INTENT(IN) :: nspin
    !! spin configuration label
    INTEGER,  INTENT(IN) :: lmax_loc
    !! In some cases I have to keep higher angular components
    !! than the default ones (=lmaxq =the ones present in rho)
    REAL(DP), INTENT(OUT):: F_lm(i%m, lmax_loc**2, nspin)
    !! lm component of F up to lmax_loc
    REAL(DP), INTENT(IN) :: F_rad(i%m, rad(i%t)%nx, nspin)
    !! radial samples of F
    !
    ! ... local variables
    !
    INTEGER :: ix    ! counter for integration
    INTEGER :: lm    ! counter for angmom
    INTEGER :: ispin ! counter for spin
    INTEGER :: j
    !
    IF (TIMING) CALL start_clock( 'PAW_rad2lm' )
    !
!$omp parallel default(shared), private(ispin,lm,ix,j)
    DO ispin = 1, nspin
!$omp do
      DO lm = 1, lmax_loc**2
        F_lm(:,lm,ispin) = 0._dp
        DO ix = ix_s, ix_e
          DO j  = 1, i%m
             F_lm(j, lm, ispin) = F_lm(j, lm, ispin) + F_rad(j,ix,ispin)* rad(i%t)%wwylm(ix,lm)
          ENDDO
        ENDDO
      ENDDO
!$omp end do
    ENDDO
!$omp end parallel
    !
    ! This routine recollects the result within the paw communicator
    !
    CALL mp_sum( F_lm, paw_comm )
    !
    IF (TIMING) CALL stop_clock( 'PAW_rad2lm' )
    !
  END SUBROUTINE PAW_rad2lm
  !
  !
  !--------------------------------------------------------------------------------------
  SUBROUTINE PAW_rad2lm3( i, F_rad, F_lm, lmax_loc, nspin )
    !-----------------------------------------------------------------------------------
    !! Computes:
    !! \[ F_{lm}(r) = \int d \Omega\ F(r,\text{th},\text{ph})\ Y_{lm}(\text{th},
    !! \text{ph}) \] 
    !! Duplicated version to work on vector fields, necessary for performance reasons.
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    INTEGER, INTENT(IN) :: lmax_loc
    !! in some cases I have to keep higher angular components.
    !! than the default ones (=lmaxq =the ones present in rho).
    INTEGER, INTENT(IN)  :: nspin
    !! spin configuration label
    REAL(DP), INTENT(OUT):: F_lm(i%m,3,lmax_loc**2,nspin)
    !! lm component of F up to lmax_loc
    REAL(DP), INTENT(IN) :: F_rad(i%m,3,rad(i%t)%nx,nspin)
    !! radial samples of F
    !
    ! ... local variables
    !
    REAL(DP) :: aux(i%m) ! optimization
    !
    INTEGER :: ix    ! counter for integration
    INTEGER :: lm    ! counter for angmom
    INTEGER :: ispin ! counter for spin
    !
    IF (TIMING) CALL start_clock( 'PAW_rad2lm3' )
    !
    ! Third try: 50% faster than blind implementation (60% with prefetch)
    DO ispin = 1,nspin
      DO lm = 1,lmax_loc**2
        !
        aux(:) = 0._DP
        DO ix = ix_s, ix_e
           aux(1:i%m) = aux(1:i%m) + F_rad(1:i%m,1,ix,ispin) * rad(i%t)%wwylm(ix,lm)
           !CALL MM_PREFETCH( F_rad(1:i%m,1,MIN(ix+1,rad(i%t)%nx),ispin), 1 )
        ENDDO
        F_lm(1:i%m, 1, lm, ispin) = aux(1:i%m)
        !
        aux(:) = 0._DP
        DO ix = ix_s, ix_e
           aux(1:i%m) = aux(1:i%m) + F_rad(1:i%m,2,ix,ispin) * rad(i%t)%wwylm(ix,lm)
           !CALL MM_PREFETCH( F_rad(1:i%m,2,MIN(ix+1,rad(i%t)%nx),ispin), 1 )
        ENDDO
        F_lm(1:i%m, 2, lm, ispin) = aux(1:i%m)
        !
        aux(:) = 0._DP
        DO ix = ix_s, ix_e
           aux(1:i%m) = aux(1:i%m) + F_rad(1:i%m,3,ix,ispin) * rad(i%t)%wwylm(ix,lm)
           !CALL MM_PREFETCH( F_rad(1:i%m,3,MIN(ix+1,rad(i%t)%nx),ispin), 1 )
        ENDDO
        F_lm(1:i%m, 3, lm, ispin) = aux(1:i%m)
        !
      ENDDO
    ENDDO
    !
    ! NB: this routine collects the result among the paw communicator 
    !
    CALL mp_sum( F_lm, paw_comm )
    !
    IF (TIMING) CALL stop_clock( 'PAW_rad2lm3' )
    !
  END SUBROUTINE PAW_rad2lm3
  !
  !

SUBROUTINE PAW_deqtranpotential(dbecsum, becsum, int3)
   USE atom,              ONLY : g => rgrid
   USE ions_base,         ONLY : nat, ityp
   USE mp,                ONLY : mp_comm_split, mp_comm_free, mp_size, mp_rank
   USE noncollin_module,  ONLY : nspin_mag
   USE lsda_mod,          ONLY : nspin
   USE uspp_param,        ONLY : nh, nhm, upf

   REAL(DP), INTENT(IN) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag) ! cross band 
                                                           ! occupations 
   COMPLEX(DP), INTENT(IN) :: dbecsum(nhm*(nhm+1)/2,nat,nspin_mag)! 
   
   COMPLEX(DP), INTENT(OUT) :: int3(nhm,nhm,1,nat,nspin_mag) ! change of 
                                           !descreening coefficients (AE - PS)
   INTEGER, PARAMETER      :: AE = 1, PS = 2,&      ! All-Electron and Pseudo
                              XC = 1, H  = 2        ! XC and Hartree
   REAL(DP), POINTER       :: rho_core(:)           ! pointer to AE/PS core charge density 
   TYPE(paw_info)          :: i                     ! minimal info on atoms
   INTEGER                 :: i_what                ! counter on AE and PS
   INTEGER                 :: is                    ! spin index
   INTEGER                 :: lm                    ! counters on angmom and radial grid
   INTEGER                 :: nb, mb, nmb           ! augfun indexes
   INTEGER                 :: ia,mykey,ia_s,ia_e    ! atoms counters and indexes
   !
   REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:) ! density expanded on Y_lm
   REAL(DP), ALLOCATABLE   :: v_lm(:,:,:)   ! workspace: potential
   REAL(DP), ALLOCATABLE   :: drhor_lm(:,:,:,:) ! change of density expanded 
                                              ! on Y_lm (real part)
   REAL(DP), ALLOCATABLE   :: drhoi_lm(:,:,:,:) ! change of density expanded 
                                              ! on Y_lm (imaginary part)
   REAL(DP), ALLOCATABLE   :: savedvr_lm(:,:,:,:)   ! workspace: potential
   REAL(DP), ALLOCATABLE   :: savedvi_lm(:,:,:,:)   ! workspace: potential
   REAL(DP), ALLOCATABLE   :: aux_lm(:) ! auxiliary radial function
   ! fake cross band occupations to select only one pfunc at a time:
   REAL(DP)                :: becfake(nhm*(nhm+1)/2,nat,nspin_mag)
   REAL(DP)                :: energy
   REAL(DP)                :: integral_r           ! workspace
   REAL(DP)                :: integral_i           ! workspace
   REAL(DP)                :: sgn                ! +1 for AE -1 for PS
   INTEGER  :: ipert, k

   CALL start_clock('PAW_dpot')
   ! Some initialization
   becfake(:,:,:) = 0._dp
   int3 = (0.0_DP, 0.0_DP)
   !
   ! Parallel: divide tasks among all the processor for this image
   ! (i.e. all the processors except for NEB and similar)
   CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )
   ! build the group of all the procs associated with the same atom
   !
   CALL mp_comm_split( intra_image_comm, ia_s - 1, me_image, paw_comm )
   !
   me_paw    = mp_rank( paw_comm )
   nproc_paw = mp_size( paw_comm )
   !

   atoms: DO ia = ia_s, ia_e
      !
      i%a = ia                      ! atom's index
      i%t = ityp(ia)                ! type of atom ia
      i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
      i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
      i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia
      !
      ifpaw: IF (upf(i%t)%tpawp) THEN
!
!    Initialize parallelization over the directions
!
         nx_loc = ldim_block( rad(i%t)%nx, nproc_paw, me_paw )
         ix_s   = gind_block( 1, rad(i%t)%nx, nproc_paw, me_paw )
         ix_e   = ix_s + nx_loc - 1
         !
         ! Arrays are allocated inside the cycle to allow reduced
         ! memory usage as differnt atoms have different meshes
         !
         ALLOCATE(v_lm(i%m,i%l**2,nspin))
         ALLOCATE(savedvr_lm(i%m,i%l**2,nspin_mag,1))
         ALLOCATE(savedvi_lm(i%m,i%l**2,nspin_mag,1))
         ALLOCATE(rho_lm(i%m,i%l**2,nspin_mag))
         ALLOCATE(drhor_lm(i%m,i%l**2,nspin_mag,1))
         ALLOCATE(drhoi_lm(i%m,i%l**2,nspin_mag,1))
         ALLOCATE(aux_lm(i%m))
         !
         whattodo: DO i_what = AE, PS
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
               CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%pfunc, rho_lm)
               rho_core => upf(i%t)%paw%ae_rho_atc
               sgn = +1._dp
            ELSE
               CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%qfuncl)
               rho_core => upf(i%t)%rho_atc 
               sgn = -1._dp                 
            ENDIF
!
!           Compute the change of the charge density. Complex because the
!           displacements might be complex
!
            DO ipert=1,1
               IF (i_what == AE) THEN
                  becfake(:,ia,:)=DBLE(dbecsum(:,ia,:))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, drhor_lm(1,1,1,ipert))
                  becfake(:,ia,:)=AIMAG(dbecsum(:,ia,:))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, drhoi_lm(1,1,1,ipert))
               ELSE
                  becfake(:,ia,:)=DBLE(dbecsum(:,ia,:))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, drhor_lm(1,1,1,ipert), upf(i%t)%qfuncl)
                  becfake(:,ia,:)=AIMAG(dbecsum(:,ia,:))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, drhoi_lm(1,1,1,ipert), upf(i%t)%qfuncl)
               END IF
            END DO

            savedvr_lm(:,:,:,:) = 0._dp
            savedvi_lm(:,:,:,:) = 0._dp

            DO ipert=1,1
               !
               ! Change of Exchange-correlation potential
               !
               CALL PAW_xc_potential(i, rho_lm, rho_core, v_lm, energy)

               DO is=1,nspin_mag
                  DO lm = 1,i%l**2
                     DO k=1,i%m
                        IF (ABS(rho_lm(k,lm,1)-rho_lm(k,lm,2))>1.D-11) THEN
                           savedvr_lm(k,lm,is,1)=0.5_DP*(v_lm(k,lm,1)- &
                                                       v_lm(k,lm,2)) &
                                          /(rho_lm(k,lm,1)-rho_lm(k,lm,2)) * &
                                          drhor_lm(k,lm,is,1)
                           savedvi_lm(k,lm,is,1)=0.5_DP*(v_lm(k,lm,1)- &
                                                       v_lm(k,lm,2)) &
                                          /(rho_lm(k,lm,1)-rho_lm(k,lm,2)) * &
                                          drhoi_lm(k,lm,is,1)
                        ELSE
                           savedvr_lm(k,lm,is,1)=0.0_DP
                           savedvi_lm(k,lm,is,1)=0.0_DP
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            END DO
            !
            spins: DO is = 1, 1
               nmb = 0
               ! loop on all pfunc for this kind of pseudo
               becfake=0.0_DP
               DO nb = 1, nh(i%t)
                  DO mb = nb, nh(i%t)
                     nmb = nmb+1 
                     becfake(nmb,ia,is) = 1._dp
                     IF (i_what == AE) THEN
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, rho_lm)
                     ELSE
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, &
                                           rho_lm, upf(i%t)%qfuncl)
                     ENDIF
!
!                 Integrate the change of xc potential and the partial waves
!                 to find the change of the D coefficients: D^1-~D^1
!
                     DO ipert=1,1
                        DO lm = 1,i%l**2
                           aux_lm(1:i%m)=rho_lm(1:i%m,lm,is)* &
                                               savedvr_lm(1:i%m,lm,is,ipert) 
                           CALL simpson (upf(i%t)%kkbeta,aux_lm, &
                                                      g(i%t)%rab,integral_r)
                           aux_lm(1:i%m)=rho_lm(1:i%m,lm,is)* &
                                            savedvi_lm(1:i%m,lm,is,ipert) 
                           CALL simpson (upf(i%t)%kkbeta,aux_lm, &
                                                      g(i%t)%rab,integral_i)
                           int3(nb,mb,ipert,i%a,is) = &
                                        int3(nb,mb,ipert,i%a,is) &
                                       + sgn * CMPLX(integral_r, integral_i,kind=DP)
                        ENDDO
                        IF (nb /= mb)  int3(mb,nb,ipert,i%a,is) = &
                                                    int3(nb,mb,ipert,i%a,is) 
                     ENDDO
                     becfake(nmb,ia,is) = 0._dp
                  ENDDO ! mb
               ENDDO ! nb
            ENDDO spins
         ENDDO whattodo
         ! cleanup
         DEALLOCATE(rho_lm)
         DEALLOCATE(drhor_lm)
         DEALLOCATE(drhoi_lm)
         DEALLOCATE(savedvr_lm)
         DEALLOCATE(savedvi_lm)
         DEALLOCATE(v_lm)
         DEALLOCATE(aux_lm)
         !
      ENDIF ifpaw
   ENDDO atoms

#if defined(__MPI)
    IF( mykey /= 0 ) int3 = 0.0_dp
    CALL mp_sum(int3, intra_image_comm)
#endif

   CALL mp_comm_free( paw_comm )

   CALL stop_clock('PAW_dpot')

END SUBROUTINE PAW_deqtranpotential

SUBROUTINE compute_rho_spin_lm(i,rho_lm,rhoout_lm,segni_rad)
!
!   This subroutine diagonalizes the spin density matrix and gives 
!   the spin-up and spin-down components of the charge. In input
!   the spin_density is decomposed into the lm components and in
!   output the spin-up and spin-down densities are decomposed into 
!   the lm components. segni_rad is an output variable with the sign
!   of the direction of the magnetization in each point.
!
USE kinds, ONLY : dp
USE constants, ONLY: eps12
USE lsda_mod,  ONLY : nspin
USE noncollin_module, ONLY : ux, nspin_gga, nspin_mag
USE uspp_param,  ONLY : upf
USE atom,       ONLY : g => rgrid
USE io_global,  ONLY :  stdout

TYPE(paw_info), INTENT(IN) :: i
         
REAL(DP), INTENT(IN) ::  rho_lm(i%m, i%l**2, nspin)
             ! input: the four components of the charge 
REAL(DP), INTENT(OUT) :: rhoout_lm(i%m, i%l**2, nspin_gga)
             ! output: the spin up and spin down charge
REAL(DP), INTENT(OUT) :: segni_rad(i%m, rad(i%t)%nx)
             ! output: keep track of the spin direction

REAL(DP) :: rho_rad(i%m, nspin)    ! auxiliary: the charge+mag along a line
REAL(DP) :: rhoout_rad(i%m, rad(i%t)%nx, nspin_gga) ! auxiliary: rho up and down along a line
REAL(DP) :: mag             ! modulus of the magnetization
REAL(DP) :: m(3)

INTEGER :: ix, k, ipol, kpol      ! counter on mesh points

IF (nspin /= 4) CALL errore('compute_rho_spin_lm','called in the wrong case',1)

segni_rad=0.0_DP

DO ix = ix_s, ix_e
   CALL PAW_lm2rad(i, ix, rho_lm, rho_rad, nspin)
   IF (with_small_so) CALL add_small_mag(i,ix,rho_rad)
   DO k=1, i%m
      rho_rad(k, 1:nspin) = rho_rad(k, 1:nspin)*g(i%t)%rm2(k)
      mag = sqrt( rho_rad(k,2)**2 + rho_rad(k,3)**2 + rho_rad(k,4)**2 )
!
! Choose rhoup and rhodw depending on the projection of the magnetization
! on the chosen direction
!
      IF (mag.LT.eps12) THEN
         segni_rad(k,ix)=1.0_DP
      ELSE
         DO ipol=1,3
            m(ipol)=rho_rad(k,1+ipol)/mag
         ENDDO
!
!  The axis ux is chosen in the corresponding routine in real space.
!
         segni_rad(k,ix)=SIGN(1.0_DP, m(1)*ux(1)+m(2)*ux(2)+m(3)*ux(3))
      ENDIF
      rhoout_rad(k, ix, 1)= 0.5d0*( rho_rad(k,1) + segni_rad(k,ix)*mag )* &
                                   g(i%t)%r2(k)
      rhoout_rad(k, ix, 2)= 0.5d0*( rho_rad(k,1) - segni_rad(k,ix)*mag )* &
                                   g(i%t)%r2(k)
   ENDDO
ENDDO   
CALL PAW_rad2lm(i, rhoout_rad, rhoout_lm, i%l, nspin_gga)

#if defined(__MPI)
CALL mp_sum( segni_rad, paw_comm )
#endif

RETURN
END SUBROUTINE compute_rho_spin_lm
!
SUBROUTINE compute_pot_nonc(i,vout_lm,v_lm,segni_rad,rho_lm)
!
!   This subroutine receives the GGA potential for spin up and
!   spin down and calculates the exchange and correlation potential and 
!   magnetic field.
!
USE kinds, ONLY : dp
USE constants, ONLY: eps12
USE lsda_mod,  ONLY : nspin
USE noncollin_module, ONLY : nspin_gga, nspin_mag
USE uspp_param, ONLY : upf
USE atom,       ONLY : g => rgrid
USE io_global,  ONLY :  stdout

TYPE(paw_info), INTENT(IN) :: i
         
REAL(DP), INTENT(IN) ::  rho_lm(i%m, i%l**2, nspin)
             ! input: the charge and magnetization densities
REAL(DP), INTENT(IN) :: vout_lm(i%m, i%l**2, nspin_gga)
             ! input: the spin up and spin down charges
REAL(DP), INTENT(IN) :: segni_rad(i%m, rad(i%t)%nx)
             ! input: keep track of the direction of the magnetization
REAL(DP), INTENT(INOUT) :: v_lm(i%m, i%l**2, nspin)
             ! output: the xc potential and magnetic field

REAL(DP) :: vsave_lm(i%m, i%l**2, nspin) ! auxiliary: v_lm is updated
REAL(DP) :: gsave_lm(i%m, i%l**2, nspin) ! auxiliary: g_lm is updated

REAL(DP) :: vout_rad(i%m, nspin_gga)  ! auxiliary: the potential along a line

REAL(DP) :: rho_rad(i%m, nspin)       ! auxiliary: the charge+mag along a line

REAL(DP) :: v_rad(i%m, rad(i%t)%nx, nspin) ! auxiliary: rho up and down along a line
REAL(DP) :: g_rad(i%m, rad(i%t)%nx, nspin) ! auxiliary: rho up and down along a line
REAL(DP) :: mag            ! modulus of the magnetization
integer :: ix, k, ipol, kpol     ! counter on mesh points

IF (nspin /= 4) CALL errore('compute_pot_nonc','called in the wrong case',1)

v_rad=0.0_DP
IF (upf(i%t)%has_so.and.i%ae==1) g_rad=0.0_DP

DO ix = ix_s, ix_e
   CALL PAW_lm2rad(i, ix, vout_lm, vout_rad, nspin_gga)
   CALL PAW_lm2rad(i, ix, rho_lm, rho_rad, nspin_mag)
   IF (with_small_so) CALL add_small_mag(i,ix,rho_rad)
   DO k=1, i%m
      rho_rad(k, 1:nspin) = rho_rad(k, 1:nspin) * g(i%t)%rm2(k)
      mag = sqrt( rho_rad(k,2)**2 + rho_rad(k,3)**2 + rho_rad(k,4)**2 )
      v_rad(k, ix, 1) = 0.5_DP * ( vout_rad(k,1) + vout_rad(k,2) )
      vs_rad(k,ix,i%a) = 0.5_DP * ( vout_rad(k,1) - vout_rad(k,2) )
!
! Choose rhoup and rhodw depending on the projection of the magnetization
! on the chosen direction
!
      IF (mag.GT.eps12) THEN
         DO ipol=2,4
            v_rad(k, ix, ipol) = vs_rad(k,ix,i%a) * segni_rad(k,ix) * & 
                                                    rho_rad(k,ipol) / mag
         ENDDO
      ENDIF
   ENDDO
   IF (with_small_so) CALL compute_g(i,ix,v_rad,g_rad)
ENDDO   

CALL PAW_rad2lm(i, v_rad, vsave_lm, i%l, nspin)

v_lm=v_lm+vsave_lm

IF (with_small_so) THEN
   CALL PAW_rad2lm(i, g_rad, gsave_lm, i%l, nspin)
   g_lm=g_lm+gsave_lm
ENDIF

RETURN
END SUBROUTINE compute_pot_nonc
!
SUBROUTINE add_small_mag(i, ix, rho_rad)

USE noncollin_module, ONLY : nspin_mag
!
!  This subroutine computes the contribution of the small component to the
!  magnetization in the noncollinear case and adds its to rho_rad.
!  The calculation is done along the radial line ix.
!
!  NB: Both the input and the output magnetizations are multiplied by
!      r^2.
!
TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
INTEGER, INTENT(IN) :: ix ! the line
REAL(DP), INTENT(INOUT)  :: rho_rad(i%m,nspin_mag)  ! the magnetization 

REAL(DP) :: msmall_rad(i%m, nspin_mag) ! auxiliary: the mag of the small 
                                          ! components along a line
REAL(DP) :: hatr(3)
INTEGER  :: k, ipol, kpol
 
CALL PAW_lm2rad(i, ix, msmall_lm, msmall_rad, nspin_mag)
hatr(1)=rad(i%t)%sin_th(ix)*rad(i%t)%cos_phi(ix)
hatr(2)=rad(i%t)%sin_th(ix)*rad(i%t)%sin_phi(ix)
hatr(3)=rad(i%t)%cos_th(ix)

DO k=1,i%m
   DO ipol=1,3
      DO kpol=1,3
         rho_rad(k,ipol+1) = rho_rad(k,ipol+1) - &
                 msmall_rad(k,kpol+1) * hatr(ipol) * hatr(kpol) * 2.0_DP
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE add_small_mag
!
SUBROUTINE compute_g(i, ix, v_rad, g_rad)
!
!   This routine receives as input B_{xc} and calculates the function G
!   described in Phys. Rev. B 82, 075116 (2010). The same routine can 
!   be used when v_rad contains the induced B_{xc}. In this case the 
!   output is the change of G.
!
   USE noncollin_module, ONLY : nspin_mag

   TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
   INTEGER, INTENT(IN) :: ix         ! the line
   REAL(DP), INTENT(IN)    :: v_rad(i%m,rad(i%t)%nx,nspin_mag) ! radial pot 
   REAL(DP), INTENT(INOUT) :: g_rad(i%m,rad(i%t)%nx,nspin_mag) 
                                                ! radial potential (small comp)

   REAL(DP) :: hatr(3)

   INTEGER :: k, ipol, kpol

   hatr(1)=rad(i%t)%sin_th(ix)*rad(i%t)%cos_phi(ix)
   hatr(2)=rad(i%t)%sin_th(ix)*rad(i%t)%sin_phi(ix)
   hatr(3)=rad(i%t)%cos_th(ix)

   DO k=1, i%m
      DO ipol=1,3
         DO kpol=1,3
!
!    v_rad contains -B_{xc} with the notation of the papers
!
            g_rad(k,ix,ipol+1)=g_rad(k,ix,ipol+1) - &
                               v_rad(k,ix,kpol+1)*hatr(kpol)*hatr(ipol)*2.0_DP
         ENDDO
      ENDDO
   ENDDO

   RETURN
END SUBROUTINE compute_g
!
!
END MODULE paw_add_onecenter
