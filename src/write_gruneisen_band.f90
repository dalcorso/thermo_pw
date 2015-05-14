!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_gruneisen_band(file_disp, file_vec)
  !
  ! read data files produced by "bands_sub" for ngeo geometries, writes
  ! a file with the gruneisen parameters: the derivatives of the phonon 
  ! dispersions with respect to the volume (gruneisen parameters)
  ! This is calculated at the volume given in input or at the volume
  ! that corresponds to the temperature given in input.
  ! 
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ntyp => nsp, ityp, amass, ityp
  USE control_thermo, ONLY : with_eigen
  USE data_files,     ONLY : flgrun
  USE thermo_mod,     ONLY : ngeo, omega_geo
  USE ph_freq_anharmonic, ONLY : vminf_t
  USE control_grun,   ONLY : temp_ph, volume_ph
  USE temperature,    ONLY : temp, ntemp
  USE mp,             ONLY : mp_bcast
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE mp_images,      ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE

  CHARACTER(LEN=256), INTENT(IN) :: file_disp

  REAL(DP) :: eps=1.d-4
  REAL(DP), ALLOCATABLE :: freq_geo(:,:,:), k(:,:), k_rap(:,:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:,:)
  INTEGER, ALLOCATABLE :: rap_geo(:,:,:)
  INTEGER :: nks, nbnd, nks_rap, nbnd_rap 
  INTEGER :: ibnd, jbnd, irap, ios, i, n, ierr, igeo
  INTEGER :: iu_grun, iumode
  INTEGER :: poly_order
  REAL(DP), ALLOCATABLE :: poly_grun(:,:), frequency(:,:), gruneisen(:,:)
  REAL(DP) :: vm
  LOGICAL, ALLOCATABLE :: high_symmetry(:), is_gamma(:)
  LOGICAL :: copy_before
  CHARACTER(LEN=256) :: filename, filedata, file_vec
  CHARACTER(LEN=6), EXTERNAL :: int_to_char

  NAMELIST /plot/ nks, nbnd
  NAMELIST /plot_rap/ nks_rap, nbnd_rap

  IF ( my_image_id /= root_image ) RETURN

  IF (flgrun == ' ') RETURN

  iumode=23
  WRITE(stdout,*)
  DO igeo = 1, ngeo(1)

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
        WRITE(stdout, '(5x,"Reading ",i4," dispersions at ",i6," k-points for&
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

     IF (with_eigen) THEN
        filename=TRIM(file_vec)//".g"//TRIM(int_to_char(igeo))
        IF (ionode) OPEN(UNIT=iumode, FILE=TRIM(filename), FORM='formatted', &
                   STATUS='old', ERR=210, IOSTAT=ios)
210     CALL mp_bcast(ios, ionode_id, intra_image_comm)
        CALL errore('write_gruneisen_band','modes are needed',ABS(ios))
     END IF
     !
     IF ( .NOT. ALLOCATED( freq_geo ) )  ALLOCATE (freq_geo(nbnd,ngeo(1),nks))
     IF ( .NOT. ALLOCATED( displa_geo ) .AND. with_eigen )  &
                                  ALLOCATE (displa_geo(nbnd,nbnd,ngeo(1),nks))
     IF ( .NOT. ALLOCATED( rap_geo ) )   ALLOCATE (rap_geo(nbnd,ngeo(1),nks))
     IF ( .NOT. ALLOCATED( k ) )         ALLOCATE (k(3,nks)) 
     IF ( .NOT. ALLOCATED( k_rap ) )     ALLOCATE (k_rap(3,nks))
     IF ( .NOT. ALLOCATED( high_symmetry ) ) ALLOCATE (high_symmetry(nks))
     IF ( .NOT. ALLOCATED( is_gamma ) ) ALLOCATE (is_gamma(nks))

     IF (ionode) THEN
        ierr=0
        DO n=1,nks
           READ(1,*,end=220,err=220)  (k(i,n), i=1,3 )
           READ(1,*,end=220,err=220)  (freq_geo(i,igeo,n),i=1,nbnd)
           READ(21,*,end=220,err=220) (k_rap(i,n),i=1,3), high_symmetry(n)
           READ(21,*,end=220,err=220) (rap_geo(i,igeo,n),i=1,nbnd)
           IF (abs(k(1,n)-k_rap(1,n))+abs(k(2,n)-k_rap(2,n))+  &
               abs(k(3,n)-k_rap(3,n))  > eps ) &
                  CALL errore('write_gruneisen_band',&
                       '("Incompatible k points in rap file")',1)
           is_gamma(n) = (( k(1,n)**2 + k(2,n)**2 + k(3,n)**2) < 1.d-12)
        ENDDO
        IF (with_eigen) &
           CALL readmodes(nat,nks,k,displa_geo,ngeo(1),igeo,ntyp,ityp,  &
                                                                 amass,iumode)

        CLOSE(UNIT=1, STATUS='KEEP')
        CLOSE(UNIT=21, STATUS='KEEP')
        CLOSE(UNIT=iumode, STATUS='KEEP')
        GOTO 222
220     ierr=1
        GOTO 222
     ENDIF
222  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
     IF (ierr==1) CALL errore('write_gruneisen_band','problem reading data',1)
  ENDDO
  CALL mp_bcast(k, ionode_id, intra_image_comm)
  CALL mp_bcast(k_rap, ionode_id, intra_image_comm)
  CALL mp_bcast(is_gamma, ionode_id, intra_image_comm)
  CALL mp_bcast(high_symmetry, ionode_id, intra_image_comm)
  CALL mp_bcast(freq_geo, ionode_id, intra_image_comm)
  IF (with_eigen) CALL mp_bcast(displa_geo, ionode_id, intra_image_comm)
  CALL mp_bcast(rap_geo, ionode_id, intra_image_comm)
!
!  Part two: Compute the Gruneisen parameters
!
  copy_before=.FALSE.
  poly_order=3
  ALLOCATE(poly_grun(poly_order,nbnd))
  ALLOCATE(frequency(nbnd,nks))
  ALLOCATE(gruneisen(nbnd,nks))
  frequency(:,:)= 0.0_DP
  gruneisen(:,:)= 0.0_DP
  IF (volume_ph==0.0_DP) THEN
     CALL evaluate_vm(temp_ph, vm, ntemp, temp, vminf_t)
  ELSE
     vm=volume_ph
  ENDIF

  WRITE(stdout,'(/,5x,"Plotting Gruneisen parameters at volume",f17.8)') vm
  IF (volume_ph==0.0_DP) &
            WRITE(stdout,'(5x,"Corresponding to T=",f17.8)') temp_ph

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
              frequency(ibnd,n)=frequency(ibnd,n-1)
           ENDDO
!
!  At the gamma point the first three frequencies vanishes
!
           DO ibnd=1,3
              frequency(ibnd,n)=0.0_DP
           ENDDO
        ENDIF
     ELSE
        IF (with_eigen) THEN
           CALL compute_freq_derivative_eigen(ngeo(1),freq_geo(1,1,n),   &
                        omega_geo, displa_geo(1,1,1,n),poly_order,poly_grun)
        ELSE 
           CALL compute_freq_derivative(ngeo,freq_geo(1,1,n),rap_geo(1,1,n), &
                                   omega_geo,poly_order,poly_grun)
        ENDIF

!
!  frequencies and gruneisen parameters are calculated at the chosen
!  volume
!
        DO ibnd=1,nbnd
           DO i=1,poly_order
              frequency(ibnd,n) = frequency(ibnd,n) + &
                    poly_grun(i,ibnd) * vm**(i-1)
              gruneisen(ibnd,n) = gruneisen(ibnd,n) - &
                    poly_grun(i,ibnd) * vm**(i-1) * (i-1.0_DP)
           END DO
           IF (frequency(ibnd,n) > 0.0_DP ) THEN
              gruneisen(ibnd,n) = gruneisen(ibnd,n) / frequency(ibnd,n)
           ELSE
              gruneisen(ibnd,n) = 0.0_DP
           ENDIF
           IF (copy_before) THEN
              gruneisen(ibnd,n-1) = gruneisen(ibnd,n)
              frequency(ibnd,n-1) = frequency(ibnd,n)
           ENDIF 
        ENDDO
        copy_before=.FALSE.
!
!  At the gamma point the first three frequencies vanishes
!
        IF (n>1.AND.is_gamma(n-1)) THEN
           DO ibnd=1,3
              frequency(ibnd,n-1)=0.0_DP
           ENDDO
        ENDIF
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
         WRITE (iu_grun, '(10x,3f13.7)') k(1,n),k(2,n),k(3,n)
         WRITE (iu_grun, '(6f13.7)') (gruneisen(ibnd,n), ibnd = 1, nbnd)
      ENDDO
      CLOSE(iu_grun)
   ENDIF
!
!  writes frequencies at the chosen volume on file
!
iu_grun=2
IF (ionode) &
   OPEN(UNIT=iu_grun, FILE=TRIM(flgrun)//'_freq', FORM='formatted', &
                  STATUS='UNKNOWN', ERR=20, IOSTAT=ios)
20 CALL mp_bcast(ios, ionode_id, intra_image_comm)
   CALL errore('write_gruneisen_band','opening dispersion file 1',ABS(ios))

   IF (ionode) THEN
      DO n=1,nks
         IF (n == 1 ) &
            WRITE (iu_grun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
                  nbnd, nks
         WRITE (iu_grun, '(10x,3f13.7)') k(1,n),k(2,n),k(3,n)
         WRITE (iu_grun, '(6f13.7)') (frequency(ibnd,n), ibnd = 1, nbnd)
      ENDDO
      CLOSE(iu_grun)
   ENDIF

   DEALLOCATE ( is_gamma )
   DEALLOCATE ( high_symmetry )
   DEALLOCATE ( k_rap )
   DEALLOCATE ( k ) 
   DEALLOCATE ( rap_geo )
   IF (with_eigen) DEALLOCATE ( displa_geo )
   DEALLOCATE ( freq_geo )
   DEALLOCATE ( poly_grun )
   DEALLOCATE ( frequency )
   DEALLOCATE ( gruneisen )

   RETURN
END SUBROUTINE write_gruneisen_band


SUBROUTINE evaluate_vm(temp_ph, vm, ntemp, temp, vminf_t)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), vminf_t(ntemp)
REAL(DP), INTENT(INOUT) :: temp_ph, vm

INTEGER :: itemp0, itemp1, itemp

itemp0=1
DO itemp=1,ntemp
   IF (temp(itemp) < temp_ph) itemp0=itemp
ENDDO

IF (itemp0 == ntemp) THEN
   WRITE(stdout,'(5x,"temp_ph too large setting to",f15.8 )') temp(ntemp-1)
   temp_ph=temp(ntemp-1)
   vm=vminf_t(ntemp-1)
   RETURN
ENDIF

IF (itemp0 == 1) THEN
   WRITE(stdout,'(5x,"temp_ph too small setting to",f15.8 )') temp(2)
   temp_ph=temp(2)
   vm=vminf_t(2)
   RETURN
ENDIF

itemp1=itemp0+1

vm = vminf_t(itemp0) + (temp_ph - temp(itemp0)) *          &
                       (vminf_t(itemp1)-vminf_t(itemp0)) / &
                       (temp(itemp1)-temp(itemp0))

RETURN
END SUBROUTINE evaluate_vm
