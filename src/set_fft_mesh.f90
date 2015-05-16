!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE set_fft_mesh
  !----------------------------------------------------------------------------
  !
  !  This routine sets the FFT mesh after knowing the geometry.
  !  It must be used in thermo_pw because setup is not able to deal
  !  with fft_fact. 
  !
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : pi
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2
  USE gvect,              ONLY : gcutm
  USE fft_base,           ONLY : dfftp, dffts
  USE grid_subroutines,   ONLY : realspace_grid_init
  USE thermo_sym,         ONLY : fft_fact
  USE gvecs,              ONLY : doublegrid, gcutms, dual
  USE wvfct,              ONLY : ecutwfc
  !
  IMPLICIT NONE
  !
  tpiba  = 2.D0 * pi / alat
  tpiba2 = tpiba**2
  !
  ! ... Compute the cut-off of the G vectors
  !
  doublegrid = ( dual > 4.D0 )
  gcutm = dual * ecutwfc / tpiba2
  !
  IF ( doublegrid ) THEN
     !
     gcutms = 4.D0 * ecutwfc / tpiba2
     !
  ELSE
     !
     gcutms = gcutm
     !
  END IF
  !
  ! ... calculate dimensions of the FFT grid
  !
  CALL clean_dfft()
  CALL realspace_grid_init ( dfftp, at, bg, gcutm, fft_fact )
  CALL realspace_grid_init ( dffts, at, bg, gcutms, fft_fact )
  !
  RETURN
  END SUBROUTINE set_fft_mesh
