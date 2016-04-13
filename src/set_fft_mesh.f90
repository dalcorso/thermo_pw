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
  USE gvecw,              ONLY : ecutwfc
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
  IF ( gcutms == gcutm ) THEN
     IF ( dffts%nr1 ==0 .AND. dffts%nr2==0 .AND. dffts%nr3==0) THEN
          dffts%nr1 = dfftp%nr1     
          dffts%nr2 = dfftp%nr2     
          dffts%nr3 = dfftp%nr3
          dffts%nr1x= dfftp%nr1x
          dffts%nr2x= dfftp%nr2x     
          dffts%nr3x= dfftp%nr3x
     END IF
  END IF
  CALL realspace_grid_init ( dffts, at, bg, gcutms, fft_fact )
  !
  RETURN
  END SUBROUTINE set_fft_mesh
