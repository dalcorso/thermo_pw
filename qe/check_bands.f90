!
! Copyright (C) 2015-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__OLDXML)
!----------------------------------------------------------------------------
LOGICAL FUNCTION check_bands( outdir, xq, iq )
  !----------------------------------------------------------------------------
  ! ...
  ! ... This routine checks if the bands for the q point iq are available 
  ! ... on disk and returns .true. if it finds them. 
  ! ... In input it requires the name of the outdir directory where the
  ! ... bands have to be searched. iq must be the current_iq
  !
  USE kinds,      ONLY : DP
  USE disp,       ONLY : lgamma_iq
  USE io_files,   ONLY : tmp_dir, prefix
  USE control_ph, ONLY : lqdir, newgrid
  USE io_global,  ONLY : ionode, ionode_id
  USE qexml_module, ONLY : qexml_init, qexml_openfile, qexml_closefile, &
                         qexml_read_bz
  USE iotk_module, ONLY : iotk_free_unit
  USE mp,         ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  REAL(DP), INTENT(IN) :: xq(3)
  CHARACTER(LEN=*) :: outdir
  CHARACTER(LEN=256)  :: dirname, filename, dir_phq, tmp_dir_save
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst, exst_restart, exst_recover
  INTEGER :: iunaux, nks, ierr
  REAL(DP), ALLOCATABLE :: xk(:,:), wk(:)
  !
  check_bands=.FALSE.
  tmp_dir_save=tmp_dir
  !
  ! We check if the file data-file.xml is present 
  ! in the directory where it should be. 
  ! If the file is present and there is a restart file, the bands are not
  ! done yet.
  ! For the gamma point done_bands might be false only with newgrid. 
  !
  IF (lqdir.AND..NOT. lgamma_iq(iq)) THEN
     dir_phq= TRIM (outdir) // TRIM(prefix) //  &
                    & '.q_' // TRIM(int_to_char(iq)) // '/'
     dirname= TRIM (dir_phq) //TRIM(prefix)//'.save'
     tmp_dir=dir_phq
  ELSE
     dirname = TRIM( outdir ) // TRIM( prefix ) // '.save'
     tmp_dir=outdir
  ENDIF
  !
  filename=TRIM(dirname) // '/data-file.xml'
  !
  IF (ionode) inquire (file = TRIM(filename), exist = exst)
  !
  CALL mp_bcast( exst, ionode_id, intra_image_comm )
  !
  exst_restart=.FALSE.
  IF (exst) CALL check_restart_recover(exst_recover, exst_restart)
  !
  IF (exst .AND. .NOT. exst_restart) THEN
     IF ( ionode ) THEN
        CALL iotk_free_unit( iunaux )
        CALL qexml_init(iunaux)
        CALL qexml_openfile( filename, 'read', BINARY = .FALSE., IERR = ierr )
        CALL qexml_read_bz( NUM_K_POINTS= nks, IERR =  ierr ) 
     END IF 
     CALL mp_bcast( nks, ionode_id, intra_image_comm )
     IF (nks > 1) THEN
!
!   check that this file has the correct set of k points. Allocate
!   enough space for the lsda casa
!
        ALLOCATE(xk(3,2*nks))
        ALLOCATE(wk(2*nks))
        IF (ionode) THEN
           CALL qexml_read_bz(XK=xk, Wk=wk, IERR = ierr) 
           CALL qexml_closefile('read', ierr)
        ENDIF
        CALL mp_bcast( xk, ionode_id, intra_image_comm )
        CALL mp_bcast( wk, ionode_id, intra_image_comm )
        IF ( ABS(xq(1)+xk(1,1)-xk(1,2)) + ABS(xq(2)+xk(2,1)-xk(2,2)) + &
             ABS(xq(3)+xk(3,1)-xk(3,2)) < 1.D-10 .AND. wk(2) == 0.0_DP) &
                                                       check_bands=.TRUE.  
        IF (lgamma_iq(iq).AND.wk(2)/=0.0_DP) check_bands=.TRUE.
        DEALLOCATE(wk)
        DEALLOCATE(xk)
     END IF
  END IF   
  !
  tmp_dir=tmp_dir_save
  !
  RETURN
  END FUNCTION check_bands

#else
!
!----------------------------------------------------------------------------
LOGICAL FUNCTION check_bands( outdir, xq, iq )
  !----------------------------------------------------------------------------
  ! ...
  ! ... This routine checks if the bands for the q point iq are available 
  ! ... on disk and returns .true. if it finds them. 
  ! ... In input it requires the name of the outdir directory where the
  ! ... bands have to be searched. iq must be the current_iq
  !
  USE kinds,      ONLY : DP
  USE disp,       ONLY : lgamma_iq
  USE io_files,   ONLY : tmp_dir, prefix
  USE control_ph, ONLY : lqdir, newgrid
  USE io_global,  ONLY : ionode, ionode_id
  USE pw_restart_new,   ONLY : pw_readschema_file
  USE qes_types_module, ONLY : output_type
  USE bcast_qes_types_module, ONLY : bcast_output_type
  USE qes_libs_module,      ONLY :  qes_reset_output
  USE io_files,   ONLY : xmlpun_schema
  USE iotk_module, ONLY : iotk_free_unit
  USE mp,         ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  REAL(DP), INTENT(IN) :: xq(3)
  CHARACTER(LEN=*) :: outdir
  CHARACTER(LEN=256)  :: dirname, filename, dir_phq, tmp_dir_save
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst, exst_restart, exst_recover
  INTEGER :: ik, iunaux, nks, ierr
  REAL(DP) :: xk(3,2), wk(2)

  TYPE(output_type) :: output_obj
  !
  check_bands=.FALSE.
  tmp_dir_save=tmp_dir
  !
  ! We check if the file data-file.xml is present 
  ! in the directory where it should be. 
  ! If the file is present and there is a restart file, the bands are not
  ! done yet.
  ! For the gamma point done_bands might be false only with newgrid. 
  !
  IF (lqdir.AND..NOT. lgamma_iq(iq)) THEN
     dir_phq= TRIM (outdir) // TRIM(prefix) //  &
                    & '.q_' // TRIM(int_to_char(iq)) // '/'
     dirname= TRIM (dir_phq) //TRIM(prefix)//'.save'
     tmp_dir=dir_phq
  ELSE
     dirname = TRIM( outdir ) // TRIM( prefix ) // '.save'
     tmp_dir=outdir
  ENDIF
  !
  filename=TRIM(dirname) // '/' // TRIM(xmlpun_schema)
  !
  IF (ionode) inquire (file = TRIM(filename), exist = exst)
  !
  CALL mp_bcast( exst, ionode_id, intra_image_comm )
  !
  exst_restart=.FALSE.
  exst_recover=.FALSE.
  IF (exst.AND.ionode) CALL check_restart_recover(exst_recover, exst_restart)
  CALL mp_bcast( exst_restart, ionode_id, intra_image_comm ) 
  CALL mp_bcast( exst_recover, ionode_id, intra_image_comm ) 
  
  !
  IF (exst .AND. .NOT. exst_restart) THEN
     IF (ionode) &
        CALL pw_readschema_file(ierr, output_obj)
     CALL mp_bcast( ierr, ionode_id, intra_image_comm )
     CALL errore('check_bands','xmlfile not readable',ABS(ierr))
     CALL bcast_output_type(output_obj, ionode_id, intra_image_comm )
     nks = output_obj%band_structure%nks 
     IF (nks > 1) THEN
!
!   check that this file has the correct set of k points. 
!
        DO ik=1,2
           xk(:,ik)=output_obj%band_structure%ks_energies(ik)%k_point%k_point(:)
           wk(ik)=output_obj%band_structure%ks_energies(ik)%k_point%weight
        ENDDO
        IF ( ABS(xq(1)+xk(1,1)-xk(1,2)) + ABS(xq(2)+xk(2,1)-xk(2,2)) + &
             ABS(xq(3)+xk(3,1)-xk(3,2)) < 1.D-10 .AND. wk(2) == 0.0_DP) &
                                                       check_bands=.TRUE.  
        IF (lgamma_iq(iq).AND.wk(2)/=0.0_DP) check_bands=.TRUE.
     END IF
     CALL qes_reset_output ( output_obj )
  END IF   
  !
  tmp_dir=tmp_dir_save
  !
  RETURN
  END FUNCTION check_bands
#endif
