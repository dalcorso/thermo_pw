MODULE io_bands
!
!   This module is used to write and read the band structure or the
!   phonon frequencies and the representations files. The file that
!   contains the bands has the form 
!   
!   a namelist called plot that contains nks and nbnd
!   nks records of the form
!   xk(1,i)  xk(2,i)  xk(3, i)
!   e(1,i), e(2,i), e(3,i) ... e(nbnd,i)
!
!   
!   The file that contains the representations has the form
!   a namelist called plot_rap that contains nks and nbnd
!   nks records of the form
!   xk(1,i), xk(2,i), xk(3,i), high_symmetry(i), gcodek(i), aux_ind(i), 
!                              gcodek_ext(i), ptypek(1,i), ptypek(2,i), 
!                              ptypek(3,i), lprojk(i), nsym(i), same_next(i)
!   if lprojk(i) is 1 thereis another line with real numbers
!      gauge(1,i), gauge(2,i), ... gauge(nsym,i) 
!   where nsym is the number of symmetry of the point group gcodek(i)
!   Finally there is a list of nbnd integers
!   rap(1,i), rap(2,i), ... rap(nbnd,i)
!   
!   xk are real,
!   high_symmetry is a logical
!   gcodek, aux_ind, gcode_ext, ptype, lprojk are integers
!   same_next is a logical and force a continous path between the two
!   points even if they are distant in reciprocal space.
!

  USE kinds,      ONLY : DP
  !
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER :: nks_=0, nbnd_=0

  PUBLIC  read_parameters, read_bands, write_bands, read_representations,  &
          write_representations

CONTAINS

  SUBROUTINE write_bands(nks, nbnd, xk, et, fact, filedata)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nks, nbnd
  REAL(DP), INTENT(IN) :: xk(3,nks), et(nbnd,nks), fact
  CHARACTER(LEN=256) :: filedata

  INTEGER :: ios, ik, ibnd, iunpun
  INTEGER :: find_free_unit

  IF ( ionode ) THEN
     iunpun=find_free_unit()
     OPEN (UNIT = iunpun, FILE = TRIM(filedata), STATUS = 'unknown', &
                           FORM = 'formatted', IOSTAT = ios)
  ENDIF
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
  IF ( ios /= 0 ) &
     CALL errore ('write_bands', 'Opening filband file', ABS(ios) )

  IF (ionode) THEN
     WRITE (iunpun, '(" &plot nbnd=",i4,", nks=",i6," /")') nbnd, nks
     DO ik=1,nks
        !
!        WRITE (iunpun, '(10x,3f10.6)') xk(:,ik)
!        WRITE (iunpun, '(10f9.3)') (et(ibnd, ik), ibnd = 1, nbnd)
        WRITE (iunpun, '(10x,3f15.10)') xk(:,ik)
        WRITE (iunpun, '(8f15.8)') (et(ibnd, ik), ibnd = 1, nbnd)
        !
     ENDDO
     CLOSE( UNIT=iunpun, STATUS='KEEP' )
  ENDIF

  RETURN
  END SUBROUTINE write_bands

  SUBROUTINE write_representations(nks, nbnd, xk, rap, high_symmetry,      &
                          gcodek, aux_ind, gcodek_ext, ptypek, lprojk, & 
                          same_next, gaugek, filedata, ik0, nks_eff) 

  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  USE point_group, ONLY : nsym_group

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nks, nbnd, ik0, nks_eff
  INTEGER, INTENT(IN) :: rap(nbnd, nks), gcodek(nks), gcodek_ext(nks), &
                         aux_ind(nks), ptypek(3,nks), lprojk(nks)
  LOGICAL, INTENT(IN) :: high_symmetry(nks), same_next(nks)
  REAL(DP), INTENT(IN) :: xk(3,nks), gaugek(48,nks)
  CHARACTER(LEN=256) :: filedata

  INTEGER :: iunpun, ios, ik, ibnd, isym, irap
  INTEGER :: find_free_unit

  IF ( ionode ) THEN
     iunpun=find_free_unit()
     OPEN (UNIT = iunpun, FILE = TRIM(filedata), STATUS = 'unknown', &
                           FORM = 'formatted', IOSTAT = ios)
  ENDIF
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
  IF ( ios /= 0 ) &
     CALL errore ('write_representations', 'Opening representation file', &
                                                                   ABS(ios) )
  IF (ionode) THEN
     WRITE (iunpun, '(" &plot_rap nbnd_rap=",i4,", nks_rap=",i6," /")') &
                                                       nbnd, nks_eff-ik0
     DO ik=ik0+1, nks_eff
        WRITE (iunpun, '(10x,3f10.6,l5,7i4,l3)') xk(1:3,ik),  &
               high_symmetry(ik), gcodek(ik), aux_ind(ik), &
               gcodek_ext(ik), ptypek(1:3,ik), lprojk(ik), same_next(ik)
        IF (lprojk(ik)==1) &
           WRITE (iunpun, '(5f16.11)') (gaugek(isym,ik)/pi, isym=1,&
                                            nsym_group(gcodek(ik)))
        WRITE (iunpun, '(10i8)') (rap(ibnd,ik), ibnd=1,nbnd)
     ENDDO
     CLOSE( UNIT=iunpun, STATUS='KEEP' )
  ENDIF

  RETURN
  END SUBROUTINE write_representations

  SUBROUTINE read_parameters(nks, nbnd, filedata)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: nks, nbnd
  CHARACTER(LEN=256) :: filedata
  
  INTEGER :: iunpun, ios
  INTEGER :: find_free_unit

  NAMELIST /plot/ nks, nbnd

  IF (ionode) THEN
     iunpun=find_free_unit()
     OPEN(UNIT=iunpun,FILE=TRIM(filedata),FORM='formatted',STATUS='OLD',ERR=10,&
                                                  IOSTAT=ios)
  ENDIF
10  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('read_parameters','opening band file',ABS(ios))

  IF (ionode) READ (iunpun, plot, ERR=20, IOSTAT=ios)
20  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('read_parameters','reading plot namelist',ABS(ios))
  CALL mp_bcast(nks, ionode_id, intra_image_comm)
  CALL mp_bcast(nbnd, ionode_id, intra_image_comm)
  !
  nks_=nks
  nbnd_=nbnd

  IF (ionode) CLOSE( UNIT=iunpun, STATUS='KEEP' )
  RETURN
  END SUBROUTINE read_parameters

  SUBROUTINE read_bands(nks, nbnd, xk, et, filedata)

  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: nks, nbnd
  REAL(DP), INTENT(INOUT) :: xk(3,nks), et(nbnd,nks)
  CHARACTER(LEN=256) :: filedata

  INTEGER :: iunpun, ik, ipol, ibnd, ios
  INTEGER :: find_free_unit

  NAMELIST /plot/ nks, nbnd

  IF ((nks_==0).OR.(nbnd_==0)) CALL errore('read_band',&
                      'Parameters not initialized',1)
  IF ((nks /= nks_) .OR. (nbnd/=nbnd_)) CALL errore('read_band',&
                      'Wrong parameters',1)
  
  IF (ionode) THEN
     iunpun=find_free_unit()
     OPEN(UNIT=iunpun,FILE=TRIM(filedata),FORM='formatted',STATUS='OLD',ERR=10,&
                                                  IOSTAT=ios)
  END IF
10  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('read_bands','opening band file',ABS(ios))

  IF (ionode) READ (iunpun, plot, ERR=20, IOSTAT=ios)
20  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('read_bands','reading plot namelist',ABS(ios))
  CALL mp_bcast(nks, ionode_id, intra_image_comm)
  CALL mp_bcast(nbnd, ionode_id, intra_image_comm)
  !
  IF (nks /= nks_ .or. nbnd /= nks_) THEN
     CALL errore('read_bands','incorrect dimensions',ABS(ios))
  ELSE
     WRITE(stdout, '(/,5x,"Reading ",i4," bands at ",i6," k-points")') nbnd, nks
  ENDIF

  IF (ionode) THEN
     DO ik=1,nks
        READ(iunpun,*,END=100,ERR=100,IOSTAT=ios) (xk(ipol,ik), ipol=1,3 )
        READ(iunpun,*,END=100,ERR=100,IOSTAT=ios) (et(ibnd,ik),ibnd=1,nbnd)
     END DO
     CLOSE( UNIT=iunpun, STATUS='KEEP' )
  ENDIF
100 CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('read_bands','reading kpoint or bands',ABS(ios))

  CALL mp_bcast(xk, ionode_id, intra_image_comm)
  CALL mp_bcast(et, ionode_id, intra_image_comm)

  RETURN
  END SUBROUTINE read_bands

  SUBROUTINE read_representations(nks, nbnd, xk, rap, high_symmetry,      &
                          gcodek, aux_ind, gcodek_ext, ptypek, lprojk, & 
                          same_next, gaugek, exist_rap, filedata) 

  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  USE point_group, ONLY : nsym_group

  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: nks, nbnd
  INTEGER, INTENT(INOUT) :: rap(nbnd, nks), gcodek(nks), gcodek_ext(nks), &
                         aux_ind(nks), ptypek(3,nks), lprojk(nks)
  LOGICAL, INTENT(INOUT) :: high_symmetry(nks), same_next(nks)
  REAL(DP), INTENT(INOUT) :: xk(3,nks), gaugek(48,nks)
  LOGICAL, INTENT(OUT) :: exist_rap
  CHARACTER(LEN=256) :: filedata

  INTEGER :: iunpun, ik, ipol, isym, ibnd, irap, ios, nks_rap, nbnd_rap
  INTEGER :: find_free_unit

  NAMELIST /plot_rap/ nks_rap, nbnd_rap

  IF ((nks_==0).OR.(nbnd_==0)) CALL errore('read_representations',&
                      'Parameters not initialized',1)
  IF ((nks /= nks_) .OR. (nbnd/=nbnd_)) CALL errore('read_representations',&
                      'Wrong parameters',1)
  exist_rap=.TRUE.

  IF (ionode) THEN
     iunpun=find_free_unit()
     OPEN(UNIT=iunpun,FILE=TRIM(filedata),FORM='formatted',STATUS='OLD',ERR=10,&
                                                  IOSTAT=ios)
  ENDIF
10  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) THEN
     exist_rap=.FALSE.
     RETURN
  ENDIF

  IF (ionode) READ (iunpun, plot_rap, ERR=110, IOSTAT=ios)
110  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios == 0 ) THEN
     CALL mp_bcast(nks_rap, ionode_id, intra_image_comm)
     CALL mp_bcast(nbnd_rap, ionode_id, intra_image_comm)
     IF (nks_rap/=nks.or.nbnd_rap/=nbnd) THEN
        WRITE(stdout,'("file with representations not compatible &
                       & with bands")')
        exist_rap=.FALSE.
        RETURN
     ENDIF
  ELSE
     WRITE(stdout,'("Problem reading representation file")')
     exist_rap=.FALSE.
     RETURN
  ENDIF

  IF (ionode) THEN
     gaugek=0.0_DP
     DO ik=1, nks
        READ(iunpun,*,END=100,ERR=100,IOSTAT=ios) (xk(ipol,ik),ipol=1,3), &
                       high_symmetry(ik), gcodek(ik), aux_ind(ik),   &
                       gcodek_ext(ik), ptypek(1:3,ik), lprojk(ik),   &
                       same_next(ik)
        IF (lprojk(ik)==1) THEN
           READ(iunpun,*,END=100,ERR=100,IOSTAT=ios) &
                             (gaugek(isym,ik),isym=1,nsym_group(gcodek(ik)))
           gaugek(:,ik)=gaugek(:,ik) * pi
        ENDIF
        READ(iunpun,*,END=100,ERR=100,IOSTAT=ios) (rap(ibnd,ik),ibnd=1,nbnd)
     ENDDO
     CLOSE( UNIT=iunpun, STATUS='KEEP' )
  ENDIF
100 CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) THEN
     exist_rap=.FALSE.
     RETURN
  ENDIF

  CALL mp_bcast(xk, ionode_id, intra_image_comm)
  CALL mp_bcast(rap, ionode_id, intra_image_comm)
  CALL mp_bcast(gcodek, ionode_id, intra_image_comm)
  CALL mp_bcast(aux_ind, ionode_id, intra_image_comm)
  CALL mp_bcast(high_symmetry, ionode_id, intra_image_comm)
  CALL mp_bcast(gcodek_ext, ionode_id, intra_image_comm)
  CALL mp_bcast(ptypek, ionode_id, intra_image_comm)
  CALL mp_bcast(lprojk, ionode_id, intra_image_comm)
  CALL mp_bcast(gaugek, ionode_id, intra_image_comm)
  CALL mp_bcast(same_next, ionode_id, intra_image_comm)

  RETURN
  END SUBROUTINE read_representations

  SUBROUTINE clean_band_reader()

  IMPLICIT NONE

  nks_=0
  nbnd_=0
  
  RETURN
  END SUBROUTINE clean_band_reader

END MODULE io_bands

