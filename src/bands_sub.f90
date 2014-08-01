!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE bands_sub()
  !-----------------------------------------------------------------------
  !
  USE control_thermo,   ONLY : filband, spin_component
  USE control_paths,    ONLY : q2d
  USE control_flags,    ONLY : gamma_only
  USE control_bands,    ONLY : lsym
  USE klist,            ONLY : two_fermi_energies
  USE noncollin_module, ONLY : i_cons
  !
  IMPLICIT NONE
  !
  IF (q2d) lsym=.FALSE.

  CALL start_clock('bands')
  CALL clean_pw(.TRUE.)
  CALL read_file()
  !
  IF (gamma_only) CALL errore('bands_sub','gamma_only case not implemented',1)
  !
  IF (two_fermi_energies .OR. i_cons /= 0) &
     CALL errore('bands_sub',&
     'The bands code with constrained magnetization has not been tested',1)
  !
  CALL openfil_pp()
  !
  IF (q2d) THEN
     CALL punch_band_2d(filband, spin_component)
  ELSE
     CALL write_bands(filband, spin_component)
     IF (lsym) CALL sym_band_sub(filband,spin_component)
  END IF
  CALL clean_pw(.FALSE.)

  CALL print_clock('bands')
  CALL stop_clock('bands')
  !
  RETURN
END SUBROUTINE bands_sub

SUBROUTINE punch_band_2d(filband,spin_component)
!
!  This routine opens a file for each band and writes on output 
!  kx, ky, energy, 
!  kx, ky, energy
!  .., .., ..
!  where kx and ky are proportional to the length
!  of the vectors k_1 and k_2 specified in the input of the 2d plot.
!
!  The k points are supposed to be in the form
!  xk(i,j) = xk_0 + dkx *(i-1) + dky * (j-1)      1<i<n1, 1<j<n2
!
!  kx(i,j) = (i-1) |dkx|
!  ky(i,j) = (j-1) |dky|
!
   USE kinds,     ONLY : DP
   USE constants, ONLY : eps8, rytoev
   USE lsda_mod,  ONLY : nspin
   USE klist,     ONLY : xk, nkstot, nks
   USE wvfct,     ONLY : et, nbnd
   USE io_files,  ONLY : iuntmp
   USE io_global, ONLY : ionode, ionode_id
   USE mp,        ONLY : mp_bcast
   USE mp_images, ONLY : intra_image_comm

   IMPLICIT NONE
   CHARACTER(LEN=256),INTENT(IN) :: filband
   INTEGER, INTENT(IN) :: spin_component
   REAL(DP) :: xk0(3), xk1(3), xk2(3), dkx(3), dky(3), xkdum(3), mdkx, mdky
   INTEGER :: n1, n2
   INTEGER :: ik, i, i1, i2, ibnd, ijk, start_k, last_k, nks_eff, j, ios
   CHARACTER(LEN=256) :: filename
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   REAL(DP), ALLOCATABLE :: xk_collect(:,:), et_collect(:,:)
   
   ALLOCATE(xk_collect(3,nkstot))
   ALLOCATE(et_collect(nbnd,nkstot))
   CALL xk_et_collect( xk_collect, et_collect, xk, et, nkstot, nks, nbnd )

   start_k=1
   last_k=nkstot
   nks_eff=nkstot
   IF (nspin==2) THEN
      nks_eff=nkstot/2
      IF (spin_component==1) THEN
         start_k=1
         last_k=nks_eff
      ELSE
         start_k=nks_eff+1
         last_k=nkstot
      ENDIF
   ENDIF
!
!  Determine xk0
!
   xk0(:)=xk_collect(:,start_k)
!
! Determine dkx
!
   dky(:)=xk_collect(:,start_k+1)-xk0(:)
!
! Determine n2 and dky
!

loop_k:  DO j=start_k+2, nkstot
     xkdum(:)=xk0(:)+(j-1)*dky(:)
     IF (ABS(xk_collect(1,j)-xkdum(1))>eps8.OR.   &
         ABS(xk_collect(2,j)-xkdum(2))>eps8.OR.   &
         ABS(xk_collect(3,j)-xkdum(3))>eps8) THEN    
         n2=j-1
         dkx(:)=xk_collect(:,j)-xk0(:)
         EXIT loop_k
     ENDIF
  ENDDO  loop_k
  n1=nks_eff/n2
  IF (n1*n2 /= nks_eff) CALL errore('punch_band_2d',&
                                    'Problems with k points',1)
  mdkx = sqrt( dkx(1)**2 + dkx(2)**2 + dkx(3)**2 )
  mdky = sqrt( dky(1)**2 + dky(2)**2 + dky(3)**2 )
!   
!  write the output, a band per file
!
  DO ibnd=1,nbnd
     filename=TRIM(filband) // '.' // TRIM(int_to_char(ibnd))
     IF (ionode) &
     open(unit=iuntmp,file=filename,status='unknown', err=100, iostat=ios)
     CALL mp_bcast(ios,ionode_id, intra_image_comm)
100  CALL errore('punch_band_2d','Problem opening outputfile',ios)
     ijk=0
     DO i1=1,n1
        DO i2=1,n2
           ijk=ijk+1
           IF (ionode) &
           WRITE(iuntmp,'(3f16.6)') mdkx*(i1-1), mdky*(i2-1), &
                                    et_collect(ibnd,ijk)*rytoev
        ENDDO 
     ENDDO
     IF (ionode) CLOSE(unit=iuntmp,status='KEEP')
  ENDDO

  DEALLOCATE(xk_collect)
  DEALLOCATE(et_collect)

  RETURN
  END SUBROUTINE punch_band_2d
!
!-----------------------------------------------------------------------
SUBROUTINE write_bands (filband, spin_component)
  !-----------------------------------------------------------------------
  !
  !    This routine writes the band energies on a file. 
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE lsda_mod,             ONLY : nspin
  USE klist,                ONLY : xk, nks, nkstot
  USE wvfct,                ONLY : et, nbnd
  USE io_global,            ONLY : ionode, ionode_id
  USE mp,                   ONLY : mp_bcast
  USE mp_images,            ONLY : intra_image_comm

  IMPLICIT NONE
  CHARACTER (LEN=*), INTENT(IN) :: filband
  INTEGER, INTENT(IN) :: spin_component
  REAL(DP), ALLOCATABLE :: xk_collect(:,:), et_collect(:,:)
  INTEGER :: iunpun, ios, ibnd, ik

  IF (filband == ' ') RETURN

  iunpun = 18
  !
  IF ( ionode ) THEN
     !
     OPEN (UNIT = iunpun, FILE = TRIM(filband), STATUS = 'unknown', FORM = &
          'formatted', IOSTAT = ios)
     REWIND (iunpun)
     !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
  IF ( ios /= 0 ) &
     CALL errore ('punch_band', 'Opening filband file', abs(ios) )
  !
  !
  IF ( spin_component/=1 .AND. nspin/=2 ) &
     CALL errore('write_bands','incorrect spin_component',1)
  IF (spin_component<1.or.spin_component>2) &
     CALL errore('write_bands','incorrect lsda spin_component',1)

  ALLOCATE(xk_collect(3,nkstot))
  ALLOCATE(et_collect(nbnd,nkstot))
  CALL xk_et_collect( xk_collect, et_collect, xk, et, nkstot, nks, nbnd )
  !
  IF ( ionode ) THEN
     !
     DO ik=1,nkstot
        IF (ik == 1) THEN
           WRITE (iunpun, '(" &plot nbnd=",i4,", nks=",i6," /")') &
             nbnd, nkstot
        ENDIF
        WRITE (iunpun, '(10x,3f10.6)') xk_collect(:,ik)
        WRITE (iunpun, '(10f8.3)') (et_collect(ibnd, ik)*rytoev,ibnd = 1, nbnd)
        !
     ENDDO
     CLOSE(UNIT=iunpun, STATUS='KEEP')
  ENDIF
  !
  DEALLOCATE(xk_collect)
  DEALLOCATE(et_collect)
  !
  RETURN
  !
END SUBROUTINE write_bands

