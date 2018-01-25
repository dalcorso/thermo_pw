!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables that can be used by thermo_pw to control
!     internal working of QE routines. These variables are used
!     only in the routines of this directory and in those of thermo_pw
!     that sets them
!
MODULE control_qe
  USE kinds,  ONLY : DP
  !
  SAVE

  LOGICAL  :: tcollect_all=.FALSE.
  LOGICAL  :: force_band_calculation=.FALSE.
  LOGICAL  :: use_ph_images=.FALSE.

END MODULE control_qe

MODULE optical
  USE kinds,  ONLY : DP
  !
  SAVE

  REAL(DP), ALLOCATABLE :: fru(:)       ! real part of the frequency
  COMPLEX(DP) :: current_w              ! current frequency

  LOGICAL :: lcfreq,        &    ! complex frequency calculation
             freq_line,     &    ! in input frequencies are given
                                 ! along lines
             lcharge,       &    ! if .true. computes the charge-charge susc.
             lchimag,       &    ! if .true. computes the mag_z-mag_z susc.
             lmagnon,       &    ! if .true. computes \chi_+-
             lall_tensor         ! if .true. computes \chi_-+

  COMPLEX (DP), ALLOCATABLE ::      &
                  intq(:,:,:),      &! nhm, nhm, nat),    integral of e^iqr Q 
                  intq_nc(:,:,:,:)   ! nhm, nhm, nat, nspin), integral of 
                                     ! e^iqr Q in the noncollinear case

  REAL(DP), ALLOCATABLE :: dmuxc_tran(:)  ! contains 1\|m| d B_xc / d|m|


  COMPLEX(DP), ALLOCATABLE :: chirr(:), &  ! charge-charge \chi
                              chirz(:), &  ! charge-mag_z \chi
                              chizr(:), &  ! mag_z-charge \chi
                              chizz(:), &  ! mag_z-mag_z \chi
                              chipm(:), &  ! \chi_+-
                              chimp(:), &  ! \chi_-+
                              chixx(:), &  ! \chi_xx
                              chixy(:), &  ! \chi_xy
                              epsm1(:)     ! epsm1

  COMPLEX(DP), ALLOCATABLE :: polarc(:,:,:)   ! polarizability (computed
                                              ! via Clausius-Mossotti relation)
                                              ! for molecules only
  COMPLEX(DP), ALLOCATABLE :: epsilonc(:,:,:) ! complex dielectric constant
  COMPLEX(DP), ALLOCATABLE :: epsilonm1c(:,:,:) ! inverse of the complex
                                              ! dielectric constant

  INTEGER :: iu1dwf, lr1dwf     ! unit for response wavefunctions at -w
 
  INTEGER :: start_freq,    &   ! initial frequency in the job
             last_freq          ! last frequency in the job

  LOGICAL :: linear_im_freq     ! if .TRUE. the imaginary part of the
                                ! frequency grows linearly
  LOGICAL :: lfreq_ev           ! when .TRUE. the frequencies are given in eV.
END MODULE optical

MODULE images_omega
  USE kinds,  ONLY : DP

  SAVE
  LOGICAL, ALLOCATABLE :: comp_f(:)

  INTEGER :: omega_group

END MODULE images_omega

MODULE zstar_add
  USE kinds,  ONLY : DP

  SAVE
  COMPLEX (DP), ALLOCATABLE ::          &
                zstareu0_rec(:,:)        ! 3, 3 * nat)
 
  LOGICAL :: done_start_zstar=.FALSE.

END MODULE zstar_add

MODULE band_computation
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npk

  SAVE
  LOGICAL :: diago_bands(npk)    ! If .TRUE. this band is not available and
                                 ! must be recomputed by diagonalization
  INTEGER :: isym_bands(npk)     ! Symmetry operation to use to rotate the
                                 ! wavefunctions of the original k
  INTEGER :: ik_origin(npk)      ! Index of the original k use to generate
                                 ! the wavefunctions of this k.
  INTEGER :: nks0                ! number of k point after reduction with 
                                 ! the point group of the solid.
                                 ! In principle this is the total number of
                                 ! k point in which we have to compute the
                                 ! bands if the mesh of k+q coincides with 
                                 ! that of k. 
  LOGICAL :: sym_for_diago=.FALSE. ! if .TRUE. when possible the bands are
                                 ! calculated by symmetrization instead of
                                 ! diagonalization

END MODULE band_computation

MODULE lr_global
  USE kinds,      ONLY : DP
  SAVE

  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       evc0(:,:,:)           ! the ground state wavefunctions at k

  COMPLEX(DP), POINTER :: &
       evq0(:,:,:)           ! the ground state wavefunctions at k+q

  COMPLEX(KIND=DP), POINTER :: &
       sevq0(:,:,:)          ! S * ground state wavefunctions at k+q

  COMPLEX(KIND=DP), ALLOCATABLE ::     &
       d0psi(:,:,:,:),    &  ! for saving the original starting vectors
       d0psi2(:,:,:)         ! for saving the original starting vectors 

  INTEGER :: size_evc1    ! size of the global vectors
  INTEGER :: rpert        ! number of perturbations that have to be computed
                          ! together in order to symmetrize
  LOGICAL :: pseudo_hermitian ! if .TRUE. use a pseudo-hermitean algorithm

END MODULE lr_global

MODULE lr_lanczos
  USE kinds,      ONLY : DP
  SAVE

  LOGICAL :: llanczos  ! if .TRUE. the lanczos algorithm is used

  COMPLEX(KIND=DP), ALLOCATABLE :: &

       evc1_old(:,:,:,:), &    ! response wavefunctions in the pw basis (last
                               ! index 1: q' using rotated SBR 2: p')
       evc1(:,:,:,:),     &    !  "    "
       evc1_new(:,:,:,:), &    !  "    "

       sevc1(:,:,:,:),    &    ! S * "    "
       sevc1_new(:,:,:,:)      ! S * "    "


  REAL(KIND=DP), ALLOCATABLE ::      &  
       beta_store(:),      &  ! coefficients of Lanczos chain
       gamma_store(:),     &  ! coefficients of Lanczos chain
       beta_store_ext(:),  &  ! extrapolated coefficients
       gamma_store_ext(:)     ! extrapolated coefficients

  COMPLEX(kind=DP), ALLOCATABLE :: zeta_store(:,:,:) ! perturbation projected 

  INTEGER :: &
       lanczos_steps, &  ! steps of the Lanczos chain
       lanczos_steps_ext ! steps of the extrapolated lanczos chain

  CHARACTER(LEN=256) :: extrapolation ! extrapolation type


  COMPLEX(KIND=DP), ALLOCATABLE :: bbk(:,:,:)  ! coefficients of S^{-1}
  COMPLEX(KIND=DP), ALLOCATABLE :: bbnc(:,:,:) ! coefficients of the inverse
                                               ! of S
  INTEGER :: iulanczos          ! iunit where the Lanczos chain is saved

  LOGICAL :: only_spectrum      ! if .TRUE. assumes that the Lanczos chain
                                ! coefficients are on file

END MODULE lr_lanczos

MODULE lr_cg
  USE kinds,      ONLY : DP
  SAVE

  LOGICAL :: lcg  ! if .TRUE. the conjugate gradient algorithm is used where
                  ! available

  COMPLEX(KIND=DP), ALLOCATABLE ::  &
      evc1(:,:,:,:),                &  ! current position
      res(:,:,:,:),                 &  ! current residual
      pres(:,:,:,:),                &  ! current preconditioned residual
      dir(:,:,:,:),                 &  ! current direction
      dir_new(:,:,:,:)                 ! A applied to the current direction

  REAL(DP), ALLOCATABLE :: &
      prec_vec(:,:,:)                  ! The preconditioning vector

END MODULE lr_cg
