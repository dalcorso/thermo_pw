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

END MODULE optical

MODULE images_omega
  USE kinds,  ONLY : DP

  LOGICAL, ALLOCATABLE :: comp_f(:)

  INTEGER :: omega_group

END MODULE images_omega

MODULE zstar_add
  USE kinds,  ONLY : DP

  COMPLEX (DP), ALLOCATABLE ::          &
                zstareu0_rec(:,:)        ! 3, 3 * nat)
 
  LOGICAL :: done_start_zstar=.FALSE.

END MODULE zstar_add
