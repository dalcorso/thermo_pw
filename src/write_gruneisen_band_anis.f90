!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_gruneisen_band_anis(file_disp, file_vec)
  !
  ! read data files produced by "bands_sub" for ngeo geometries, writes
  ! several files with the Gruneisen parameters: the derivatives of the 
  ! phonon dispersions with respect to the celldm parameters.
  ! This is calculated at the celldm given in input or at the celldm
  ! that corresponds to the temperature given in input.
  ! 
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ntyp => nsp
  USE cell_base,      ONLY : celldm
  USE data_files,     ONLY : flgrun
  USE control_paths,  ONLY : nqaux, disp_q, disp_nqs
  USE thermo_mod,     ONLY : ngeo, omega_geo, celldm_geo, no_ph
  USE anharmonic, ONLY : celldm_t
  USE control_grun,   ONLY : temp_ph, volume_ph, celldm_ph
  USE control_pwrun,  ONLY : ibrav_save, amass_save, ityp_save
  USE temperature,    ONLY : temp, ntemp
  USE quadratic_surfaces, ONLY : evaluate_fit_quadratic, &
                                 evaluate_fit_grad_quadratic
  USE mp,             ONLY : mp_bcast
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE mp_images,      ONLY : intra_image_comm, root_image, my_image_id

  IMPLICIT NONE

  CHARACTER(LEN=256), INTENT(IN) :: file_disp, file_vec

  REAL(DP) :: eps=1.d-4
  REAL(DP), ALLOCATABLE :: freq_geo(:,:,:), k(:,:), k_rap(:,:)
  COMPLEX(DP), ALLOCATABLE :: displa_geo(:,:,:,:)
  INTEGER, ALLOCATABLE :: rap_geo(:,:,:)
  INTEGER :: nks, nbnd, nks_rap, nbnd_rap 
  INTEGER :: ibnd, jbnd, irap, ios, i, n, ierr, igeo, nwork
  INTEGER :: iu_grun, iufreq, iurap, iumode
  INTEGER :: nvar, degree, icrys, compute_nwork
  REAL(DP), ALLOCATABLE :: poly_grun(:,:), frequency(:,:), gruneisen(:,:,:)
  REAL(DP) :: vm, cm(6), f
  REAL(DP), ALLOCATABLE :: grad(:), x(:)
  LOGICAL, ALLOCATABLE :: high_symmetry(:), is_gamma(:)
  LOGICAL :: copy_before, exst_rap
  CHARACTER(LEN=256) :: filename, filedata
  CHARACTER(LEN=6), EXTERNAL :: int_to_char

  NAMELIST /plot/ nks, nbnd
  NAMELIST /plot_rap/ nks_rap, nbnd_rap

  IF ( my_image_id /= root_image ) RETURN

  IF (flgrun == ' ') RETURN

  iufreq=1
  iurap=21
  iumode=22
  exst_rap=.TRUE.

  nwork=compute_nwork()
  DO igeo = 1, nwork
     
     IF (no_ph(igeo)) CYCLE
     filedata = TRIM(file_disp)//'.g'//TRIM(int_to_char(igeo))

     IF (ionode) &
        OPEN(UNIT=iufreq,FILE=TRIM(filedata),FORM='formatted',&
                                             STATUS='OLD', ERR=101, IOSTAT=ios)
101  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band_anis','opening dispersion file',ABS(ios))

     IF (ionode) READ (iufreq, plot, IOSTAT=ios)
     CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band_anis','reading plot namelist',ABS(ios))
     CALL mp_bcast(nks, ionode_id, intra_image_comm)
     CALL mp_bcast(nbnd, ionode_id, intra_image_comm)
     !
     IF (nks <= 0 .or. nbnd <= 0) THEN
        CALL errore('write_gruneisen_band_anis','reading plot namelist',&
                                                                 ABS(ios))
     ELSE
        WRITE(stdout, '("Reading ",i4," dispersions at ",i6," k-points for&
                       & geometry",i4)') nbnd, nks, igeo
     ENDIF

     filename=TRIM(filedata)//".rap"
     IF (ionode) OPEN(UNIT=iurap, FILE=TRIM(filename), FORM='formatted', &
                   STATUS='old', ERR=100, IOSTAT=ios)
100  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     IF (ios /= 0) exst_rap=.FALSE.

     IF (exst_rap) THEN
        IF (ionode) READ (iurap, plot_rap, ERR=110, IOSTAT=ios)
110     CALL mp_bcast(ios, ionode_id, intra_image_comm)
        CALL errore('write_gruneisen_band_anis','problem reading &
                                           &representations',ABS(ios))
        CALL mp_bcast(nks_rap, ionode_id, intra_image_comm)
        CALL mp_bcast(nbnd_rap, ionode_id, intra_image_comm)
        IF ( nks_rap/=nks .OR. nbnd_rap/=nbnd ) &
        CALL errore('write_gruneisen_band_anis','("file with representations &
                       & not compatible with bands")')
     ELSE
        nks_rap=nks
        nbnd_rap=nbnd
     ENDIF

     filename=TRIM(file_vec)//".g"//TRIM(int_to_char(igeo))
     IF (ionode) OPEN(UNIT=iumode, FILE=TRIM(filename), FORM='formatted', &
                   STATUS='old', ERR=200, IOSTAT=ios)
200  CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band_anis','modes are needed',ABS(ios))

     !
     IF ( .NOT. ALLOCATED( freq_geo ) )  ALLOCATE (freq_geo(nbnd,nwork,nks))
     IF ( .NOT. ALLOCATED( displa_geo ) )  &
                                  ALLOCATE (displa_geo(nbnd,nbnd,nwork,nks))
     IF (exst_rap) THEN
        IF ( .NOT. ALLOCATED( k_rap ) )     ALLOCATE (k_rap(3,nks))
        IF ( .NOT. ALLOCATED( rap_geo ) )   ALLOCATE (rap_geo(nbnd,nwork,nks))
     ENDIF
     IF ( .NOT. ALLOCATED( k ) )         ALLOCATE (k(3,nks)) 
     IF ( .NOT. ALLOCATED( high_symmetry ) ) ALLOCATE (high_symmetry(nks))
     IF ( .NOT. ALLOCATED( is_gamma ) ) ALLOCATE (is_gamma(nks))

     IF (ionode) THEN
        ierr=0
        DO n=1,nks
           READ(iufreq,*,end=220,err=220)  (k(i,n), i=1,3 )
           READ(iufreq,*,end=220,err=220)  (freq_geo(i,igeo,n),i=1,nbnd)
           IF (exst_rap) THEN
              READ(iurap,*,end=220,err=220) (k_rap(i,n),i=1,3), high_symmetry(n)
              READ(iurap,*,end=220,err=220) (rap_geo(i,igeo,n),i=1,nbnd)
              IF (abs(k(1,n)-k_rap(1,n))+abs(k(2,n)-k_rap(2,n))+  &
                  abs(k(3,n)-k_rap(3,n))  > eps ) &
                     CALL errore('write_gruneisen_band_anis',&
                       '("Incompatible k points in rap file")',1)
           ENDIF
           is_gamma(n) = (( k(1,n)**2 + k(2,n)**2 + k(3,n)**2) < 1.d-12)
        ENDDO
!
!  readmodes reads a file with the displacements, but writes in displa_geo
!  the normalized eigenvectors of the dynamical matrix
!
        CALL readmodes(nat,nks,k,displa_geo,nwork,igeo,ntyp,ityp_save, &
                                       amass_save, iumode)
        CLOSE(UNIT=iumode, STATUS='KEEP')
        CLOSE(UNIT=iufreq, STATUS='KEEP')
        CLOSE(UNIT=iurap, STATUS='KEEP')
        GOTO 222
220     ierr=1
        GOTO 222
     ENDIF
222  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
     IF (ierr==1) CALL errore('write_gruneisen_band_anis','problem reading data',1)
  ENDDO
  CALL mp_bcast(k, ionode_id, intra_image_comm)
  CALL mp_bcast(is_gamma, ionode_id, intra_image_comm)
  CALL mp_bcast(high_symmetry, ionode_id, intra_image_comm)
  CALL mp_bcast(freq_geo, ionode_id, intra_image_comm)
  CALL mp_bcast(displa_geo, ionode_id, intra_image_comm)
  IF (exst_rap) THEN
     CALL mp_bcast(rap_geo, ionode_id, intra_image_comm)
     CALL mp_bcast(k_rap, ionode_id, intra_image_comm)
  ENDIF
!
!  Part two: Compute the Gruneisen parameters
!
  copy_before=.FALSE.
  CALL compute_degree(ibrav_save,degree,nvar)
  
  ALLOCATE(poly_grun(nvar,nbnd))
  ALLOCATE(frequency(nbnd,nks))
  ALLOCATE(gruneisen(degree,nbnd,nks))
  ALLOCATE(x(degree))
  ALLOCATE(grad(degree))
  frequency(:,:)= 0.0_DP
  gruneisen(:,:,:)= 0.0_DP
  IF (celldm_ph(1)==0.0_DP) THEN
     CALL evaluate_celldm(temp_ph, cm, ntemp, temp, celldm_t)
  ELSE
     cm=celldm_ph(:)
  ENDIF
  CALL compute_x(cm,x,degree,ibrav_save)
!
!  calculate the BZ path that corresponds to the cm parameters
!
  celldm(:)=cm(:)
  CALL set_bz_path()
  IF (nqaux > 0) CALL set_paths_disp()
  IF (disp_nqs /= nks) &
     CALL errore('write_gruneisen_band_anis','Problem with the path',1)

  WRITE(stdout,'(/,5x,"Plotting Gruneisen parameters at celldm:")') 
  WRITE(stdout,'(5x,6f15.7)') cm 
  IF (celldm_ph(1)==0.0_DP) &
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
              gruneisen(:,ibnd,n)=gruneisen(:,ibnd,n-1)
              frequency(ibnd,n)=frequency(ibnd,n-1)
           END DO
!
!  At the gamma point the first three frequencies vanishes
!
           DO ibnd=1,3
              frequency(ibnd,n)=0.0_DP
           END DO
        ENDIF
     ELSE
        CALL compute_freq_derivative_anis_eigen(nwork, freq_geo(1,1,n),&
              celldm_geo, displa_geo(1,1,1,n), degree, nvar, ibrav_save, &
                                                       no_ph, poly_grun)
!
!  frequencies and gruneisen parameters are calculated at the chosen
!  volume
!
        DO ibnd=1,nbnd
           CALL evaluate_fit_quadratic(degree,nvar,x,f,poly_grun(1,ibnd))
           CALL evaluate_fit_grad_quadratic(degree,nvar,x,grad,&
                                                          poly_grun(1,ibnd))
           frequency(ibnd,n) = f 
           gruneisen(:,ibnd,n) = -grad(:)
           IF (frequency(ibnd,n) > 0.0_DP ) THEN
              DO i=1,degree
                 gruneisen(i,ibnd,n) = x(i) * gruneisen(i,ibnd,n) /  &
                                              frequency(ibnd,n)
              END DO
           ELSE
              gruneisen(:,ibnd,n) = 0.0_DP
           ENDIF
           IF (copy_before) THEN
              gruneisen(:,ibnd,n-1) = gruneisen(:,ibnd,n)
              frequency(ibnd,n-1) = frequency(ibnd,n)
           ENDIF 
        ENDDO
        copy_before=.FALSE.
!
!  At the gamma point the first three frequencies vanish
!
        IF (n>1.AND.is_gamma(n-1)) THEN
           DO ibnd=1,3
              frequency(ibnd,n-1)=0.0_DP
           ENDDO
        ENDIF
     ENDIF
  ENDDO
!
!  Writes Gruneisen parameters on file
!
  DO icrys=1, degree
     iu_grun=2
     IF (ionode) &
        OPEN(UNIT=iu_grun, FILE=TRIM(flgrun)//'_'//TRIM(INT_TO_CHAR(icrys)), &
             FORM='formatted', STATUS='UNKNOWN', ERR=20, IOSTAT=ios)
20   CALL mp_bcast(ios, ionode_id, intra_image_comm)
     CALL errore('write_gruneisen_band_anis','opening c-gruneisen file', &
                                              ABS(ios))

     IF (ionode) THEN
        DO n=1,nks
           IF (n == 1 ) &
              WRITE (iu_grun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
                    nbnd, nks
           WRITE (iu_grun, '(10x,3f13.7)') disp_q(1,n),disp_q(2,n),disp_q(3,n)
           WRITE (iu_grun, '(6f13.7)') (gruneisen(icrys,ibnd,n), &
                                                     ibnd = 1, nbnd)
        END DO
        CLOSE(iu_grun)
     END IF
  END DO
!
!  writes frequencies at the chosen geometry on file
!
iu_grun=2
IF (ionode) &
   OPEN(UNIT=iu_grun, FILE=TRIM(flgrun)//'_freq', FORM='formatted', &
                  STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, ionode_id, intra_image_comm)
   CALL errore('write_gruneisen_band_anis','opening dispersion file 1',ABS(ios))

IF (ionode) THEN
   DO n=1,nks
      IF (n == 1 ) &
          WRITE (iu_grun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
                  nbnd, nks
      WRITE (iu_grun, '(10x,3f13.7)') disp_q(1,n),disp_q(2,n),disp_q(3,n)
      WRITE (iu_grun, '(6f13.7)') (frequency(ibnd,n), ibnd = 1, nbnd)
   ENDDO
   CLOSE(iu_grun)
ENDIF

DEALLOCATE ( is_gamma )
DEALLOCATE ( high_symmetry )
DEALLOCATE ( k ) 
DEALLOCATE ( freq_geo )
DEALLOCATE ( displa_geo )
DEALLOCATE ( poly_grun )
DEALLOCATE ( frequency )
DEALLOCATE ( gruneisen )
DEALLOCATE ( x )
DEALLOCATE ( grad )
IF (exst_rap) THEN
   DEALLOCATE ( k_rap )
   DEALLOCATE ( rap_geo )
ENDIF

RETURN
END SUBROUTINE write_gruneisen_band_anis


SUBROUTINE evaluate_celldm(temp_ph, cm, ntemp, temp, celldmf_t)

USE kinds, ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: temp(ntemp), celldmf_t(6,ntemp)
REAL(DP), INTENT(INOUT) :: temp_ph, cm(6)

INTEGER :: itemp0, itemp1, itemp

itemp0=1
DO itemp=1,ntemp
   IF (temp(itemp) < temp_ph) itemp0=itemp
ENDDO

IF (itemp0 == ntemp) THEN
   WRITE(stdout,'(5x,"temp_ph too large setting to",f15.8 )') temp(ntemp-1)
   temp_ph=temp(ntemp-1)
   cm(:)=celldmf_t(:,ntemp-1)
   RETURN
ENDIF

IF (itemp0 == 1) THEN
   WRITE(stdout,'(5x,"temp_ph too small setting to",f15.8 )') temp(2)
   temp_ph=temp(2)
   cm(:)=celldmf_t(:,2)
   RETURN
ENDIF

itemp1=itemp0+1

cm(:) = celldmf_t(:,itemp0) + (temp_ph - temp(itemp0)) *           &
                       (celldmf_t(:,itemp1)-celldmf_t(:,itemp0)) / &
                       (temp(itemp1)-temp(itemp0))

RETURN
END SUBROUTINE evaluate_celldm

SUBROUTINE compute_x(cm,x,degree,ibrav)

USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: degree, ibrav
REAL(DP), INTENT(IN) :: cm(6)
REAL(DP), INTENT(INOUT) :: x(degree)

SELECT CASE (ibrav)
   CASE(1,2,3)
      x(1) = cm(1)
   CASE(4,5,6,7)
      x(1) = cm(1)
      x(2) = cm(3)
      IF (ibrav==5) x(2) = ACOS(cm(4))
   CASE(8,9,91,10,11)
      x(1) = cm(1)
      x(2) = cm(2)
      x(3) = cm(3)
   CASE(12,-12,13,-13)
      x(1) = cm(1)
      x(2) = cm(2)
      x(3) = cm(3)
      IF (ibrav>0) THEN
!
!   c unique
!
         x(4) = ACOS(cm(4))
      ELSE
!
!   b unique
!
         x(4) = ACOS(cm(5))
      ENDIF
   CASE DEFAULT
      x(1) = cm(1)
      x(2) = cm(2)
      x(3) = cm(3)
      x(4) = ACOS( cm(4) )
      x(5) = ACOS( cm(5) )
      x(6) = ACOS( cm(6) )
END SELECT

RETURN
END SUBROUTINE compute_x
