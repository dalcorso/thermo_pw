!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_gruneisen_band(file_disp)
  !
  ! read data files produced by "bands_sub" for ngeo geometries, writes
  ! a file with the gruneisen parameters: the derivatives of the phonon 
  ! dispersions with respect to the volume (gruneisen parameters)
  ! 
  USE kinds,          ONLY : DP
  USE control_thermo, ONLY : flgrun
  USE thermo_mod,     ONLY : tot_ngeo, omega_geo
  USE ph_freq_anharmonic, ONLY : vminf_t
  USE mp,             ONLY : mp_bcast
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE mp_images,      ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE

  CHARACTER(LEN=256), INTENT(IN) :: file_disp

  INTEGER, PARAMETER :: m1 = 3   ! number of polynomial coefficients
  REAL(DP) :: alpha(m1)          ! the polynomial coefficients

  REAL(DP), ALLOCATABLE :: freq_geo(:,:,:), k(:,:)
  REAL(DP), ALLOCATABLE :: k_rap(:,:), frequences(:), gruneisen(:,:)
  REAL(DP) :: k1(3), k2(3), ps, omega, freq
  REAL(DP) :: eps=1.d-4
  INTEGER, ALLOCATABLE :: rap_geo(:,:,:)
  INTEGER, ALLOCATABLE :: level(:,:)
  INTEGER :: nks = 0, nbnd = 0
  INTEGER :: nks_rap = 0, nbnd_rap = 0
  INTEGER :: ibnd, jbnd, irap, ios, i, n, ierr, igeo, geo, use_geo
  INTEGER :: iu_grun
  LOGICAL, ALLOCATABLE :: high_symmetry(:), is_gamma(:)
  LOGICAL :: copy_before
  CHARACTER(LEN=256) :: filename, filedata
  CHARACTER(LEN=6), EXTERNAL :: int_to_char

  NAMELIST /plot/ nks, nbnd
  NAMELIST /plot_rap/ nks_rap, nbnd_rap

  IF ( my_image_id /= root_image ) RETURN

  IF (flgrun == ' ') RETURN

  IF ( .NOT. ALLOCATED( level ) ) ALLOCATE (level(12,tot_ngeo))
  IF ( .NOT. ALLOCATED( frequences ) ) ALLOCATE (frequences(tot_ngeo))

  DO igeo = 1, tot_ngeo

     filedata = TRIM(file_disp)//'.g'//TRIM(int_to_char(igeo))

     IF (ionode) &
        OPEN(UNIT=1,FILE=TRIM(filedata),FORM='formatted',STATUS='OLD',ERR=101,&
                                                  IOSTAT=ios)
101  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band','opening dispersion file',ABS(ios))

     IF (ionode) READ (1, plot, IOSTAT=ios)
     CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band','reading plot namelist',ABS(ios))
     CALL mp_bcast(nks, ionode_id, intra_image_comm)
     CALL mp_bcast(nbnd, ionode_id, intra_image_comm)
     !
     IF (nks <= 0 .or. nbnd <= 0) THEN
        CALL errore('write_gruneisen_band','reading plot namelist',ABS(ios))
     ELSE
        WRITE(stdout, '("Reading ",i4," dispersions at ",i6," k-points for&
                       & geometry",i4)') nbnd, nks, igeo
     ENDIF

     filename=TRIM(filedata)//".rap"
     IF (ionode) OPEN(UNIT=21, FILE=TRIM(filename), FORM='formatted', &
                   STATUS='old', ERR=100, IOSTAT=ios)
100  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band','representations are needed',ABS(ios))

     IF (ionode) READ (21, plot_rap, ERR=110, IOSTAT=ios)
110  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band','problem reading &
                                           &representations',ABS(ios))
     CALL mp_bcast(nks_rap, ionode_id, intra_image_comm)
     CALL mp_bcast(nbnd_rap, ionode_id, intra_image_comm)
     IF ( nks_rap/=nks .OR. nbnd_rap/=nbnd ) &
        CALL errore('write_gruneisen_band','("file with representations &
                       & not compatible with bands")')
     !
     IF ( .NOT. ALLOCATED( freq_geo ) )  ALLOCATE (freq_geo(nbnd,nks,tot_ngeo))
     IF ( .NOT. ALLOCATED( rap_geo ) )   ALLOCATE (rap_geo(nbnd,nks,tot_ngeo))
     IF ( .NOT. ALLOCATED( k ) )         ALLOCATE (k(3,nks)) 
     IF ( .NOT. ALLOCATED( k_rap ) )     ALLOCATE (k_rap(3,nks))
     IF ( .NOT. ALLOCATED( high_symmetry ) ) ALLOCATE (high_symmetry(nks))
     IF ( .NOT. ALLOCATED( is_gamma ) ) ALLOCATE (is_gamma(nks))
     IF ( .NOT. ALLOCATED( gruneisen ) ) ALLOCATE (gruneisen(nbnd,nks))

     IF (ionode) THEN
        ierr=0
        DO n=1,nks
           READ(1,*,end=220,err=220)  (k(i,n), i=1,3 )
           READ(1,*,end=220,err=220)  (freq_geo(i,n,igeo),i=1,nbnd)
           READ(21,*,end=220,err=220) (k_rap(i,n),i=1,3), high_symmetry(n)
           READ(21,*,end=220,err=220) (rap_geo(i,n,igeo),i=1,nbnd)
           IF (abs(k(1,n)-k_rap(1,n))+abs(k(2,n)-k_rap(2,n))+  &
               abs(k(3,n)-k_rap(3,n))  > eps ) &
                  CALL errore('write_gruneisen_band',&
                       '("Incompatible k points in rap file")',1)
           is_gamma(n) = (( k(1,n)**2 + k(2,n)**2 + k(3,n)**2) < 1.d-12)
        ENDDO
        CLOSE(UNIT=1, STATUS='KEEP')
        CLOSE(UNIT=21, STATUS='KEEP')
        GOTO 222
220     ierr=1
        GOTO 222
     ENDIF
222  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
     IF (ierr==1) CALL errore('plotband_sub','problem reading data',1)
  ENDDO
  CALL mp_bcast(k, ionode_id, intra_image_comm)
  CALL mp_bcast(k_rap, ionode_id, intra_image_comm)
  CALL mp_bcast(is_gamma, ionode_id, intra_image_comm)
  CALL mp_bcast(high_symmetry, ionode_id, intra_image_comm)
  CALL mp_bcast(freq_geo, ionode_id, intra_image_comm)
  CALL mp_bcast(rap_geo, ionode_id, intra_image_comm)
!
!  Part two: Compute the Gruneisen parameters
!
  copy_before=.FALSE.
  IF (mod(tot_ngeo,2)==0) THEN
     use_geo=4
     geo=tot_ngeo/2-2
  ELSE
     use_geo=3
     geo=tot_ngeo/2-1
  ENDIF
  DO n = 1,nks
     IF (is_gamma(n)) THEN
!
!    In the gamma point the Gruneisen parameters are not defined.
!    In order to have a continuous plot we take the same parameters
!    of the previous point, if this point exists and is not gamma.
!    Otherwise at the next point we copy the parameters in the present one
!
        copy_before=.FALSE.
        IF (n==1) THEN
           copy_before=.TRUE.
        ELSEIF (is_gamma(n-1)) THEN
           copy_before=.TRUE.
        ELSE
           DO ibnd=1,nbnd
              gruneisen(ibnd,n)=gruneisen(ibnd,n-1)
           ENDDO
        ENDIF
     ELSE
        level=1
        DO ibnd=1,nbnd
!
!   there are several representation files, we choose to order the
!   Gruneisen parameters on file as those of the central geometry
!
           irap=rap_geo(ibnd,n,geo+2)
            
           IF (irap == -1) THEN
              DO igeo=1,use_geo
                 frequences(igeo)=freq_geo(ibnd,n,geo+igeo)
              ENDDO
           ELSE
              DO igeo=1,use_geo
                 DO jbnd=level(irap,igeo),nbnd
                    IF (rap_geo(jbnd,n,geo+igeo)==irap) THEN
                       frequences(igeo)=freq_geo(jbnd,n,geo+igeo)
                       level(irap,igeo)=jbnd+1
                       GOTO 20
                    ENDIF
                 ENDDO
                 CALL errore('write_gruneisen_band','representation not found',1)
20               CONTINUE
              ENDDO
           ENDIF
!
!    Use the volume and frequencies of the central geometry if there is 
!    an odd number of geometries, or the average of the two central
!    geometries if there is an even number
!
           IF (mod(tot_ngeo,2)==0) THEN
              omega = 0.5_DP * (omega_geo(tot_ngeo/2) + omega_geo(tot_ngeo/2+1))
              freq  = 0.5_DP*(freq_geo(ibnd,n,tot_ngeo/2)+freq_geo(ibnd,n,&
                                              tot_ngeo/2+1))
           ELSE
              omega = omega_geo(tot_ngeo/2+1)
              freq = freq_geo(ibnd,n,tot_ngeo/2+1)
           ENDIF

           CALL polifit( omega_geo(geo+1), frequences, use_geo, alpha, m1 )
           gruneisen(ibnd,n) = - (alpha(2) + 2.0_DP * alpha(3) * &
                                              vminf_t(1)) * omega / freq
           IF (copy_before) gruneisen(ibnd,n-1)=gruneisen(ibnd,n)
        ENDDO
        copy_before=.FALSE.
     ENDIF
  ENDDO
!
!  Third part: writes Gruneisen parameters on file
!
iu_grun=2
IF (ionode) &
   OPEN(UNIT=iu_grun, FILE=TRIM(flgrun), FORM='formatted', STATUS='UNKNOWN', &
             ERR=10, IOSTAT=ios)
10 CALL mp_bcast(ios, ionode_id, intra_image_comm)
   CALL errore('write_gruneisen_band','opening dispersion file',ABS(ios))

   IF (ionode) THEN
      DO n=1,nks
         IF (n == 1 ) &
            WRITE (iu_grun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
                  nbnd, nks
         WRITE (iu_grun, '(10x,3f10.6)') k(1,n),k(2,n),k(3,n)
         WRITE (iu_grun, '(10f8.3)') (gruneisen(ibnd,n), ibnd = 1, nbnd)
      ENDDO
      CLOSE(iu_grun)
   ENDIF

   DEALLOCATE ( gruneisen )
   DEALLOCATE ( is_gamma )
   DEALLOCATE ( high_symmetry )
   DEALLOCATE ( k_rap )
   DEALLOCATE ( k ) 
   DEALLOCATE ( rap_geo )
   DEALLOCATE ( freq_geo )
   DEALLOCATE ( level )
   DEALLOCATE ( frequences )

   RETURN
END SUBROUTINE write_gruneisen_band
