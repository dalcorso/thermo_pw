!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE set_fft_mesh()
  !----------------------------------------------------------------------------
  !
  !  This routine sets the FFT mesh after knowing the geometry.
  !  It must be used in thermo_pw because setup is not able to deal
  !  with fft_fact. 
  !
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : pi
  USE cell_base,          ONLY : ibrav, at, bg, alat, tpiba, tpiba2
  USE gvect,              ONLY : gcutm
  USE fft_base,           ONLY : dfftp, dffts
  USE thermo_sym,         ONLY : fft_fact
  USE gvecs,              ONLY : doublegrid, gcutms, dual
  USE gvecw,              ONLY : ecutwfc
  USE mp_bands,           ONLY : intra_bgrp_comm

  !
  IMPLICIT NONE
  INTEGER :: nmax
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
  ENDIF
  !
  ! ... calculate dimensions of the FFT grid
  !
  CALL clean_dfft()
  dfftp%comm=intra_bgrp_comm

  CALL realspace_grid_init_tpw ( dfftp, at, bg, gcutm, fft_fact )
  IF (ibrav==10) THEN
!
!  The face-centered orthorombic lattice needs the three values of nr1, nr2,
!  nr3 equal to exploit all the symmetry
!
     nmax=MAX(dfftp%nr1, dfftp%nr2, dfftp%nr3)
     dfftp%nr1=nmax
     dfftp%nr2=nmax
     dfftp%nr3=nmax
  ENDIF

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
  dffts%comm=intra_bgrp_comm
  CALL realspace_grid_init_tpw ( dffts, at, bg, gcutms, fft_fact )
  IF (ibrav==10) THEN
!
!  The face-centered orthorombic lattice needs the three values of nr1, nr2,
!  nr3 equal to exploit all the symmetry
!
     nmax=MAX(dffts%nr1, dffts%nr2, dffts%nr3)
     dffts%nr1=nmax
     dffts%nr2=nmax
     dffts%nr3=nmax
  ENDIF
  !
  RETURN
  END SUBROUTINE set_fft_mesh
