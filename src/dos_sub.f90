!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE dos_sub()
  !--------------------------------------------------------------------
  !
  ! Calculates the Density of States (DOS),
  ! separated into up and down components for LSDA
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : rytoev
  USE klist,         ONLY : wk, degauss, ngauss, lgauss, ltetra, &
                            nks, nelec
  USE ktetra,        ONLY : tetra, tetra_dos_t
  USE wvfct,         ONLY : nbnd, et
  USE lsda_mod,      ONLY : nspin
  USE ener,          ONLY : ef
  USE control_eldos, ONLY : legauss, sigmae, deltae, ndose, save_ndos
  USE control_bands, ONLY : emax_input
  USE data_files,    ONLY : fleldos
  USE pw_restart_new, ONLY : read_xml_file
  USE io_global,     ONLY : stdout, ionode
  USE mp,            ONLY : mp_sum
  USE mp_images,     ONLY : my_image_id, root_image, intra_image_comm
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: filedos
  REAL(DP) :: emin, emax
  REAL(DP), ALLOCATABLE :: e(:), dosofe(:,:), dosint(:)
  REAL(DP) :: save_degauss, ef1
  INTEGER  :: save_ngauss
  LOGICAL  :: save_lgauss, save_ltetra, wfc_is_collected
  INTEGER  :: n, ndos, nstart, nlast, iu_dos
  INTEGER  :: find_free_unit
  !
  IF ( my_image_id /= root_image ) RETURN
  !
  CALL clean_pw(.TRUE.)
  CALL read_xml_file( wfc_is_collected )
  !
  save_degauss=degauss
  save_ngauss=ngauss
  save_ltetra=ltetra
  save_lgauss=lgauss
  IF (sigmae/=0.d0) THEN
     degauss=sigmae
     IF (legauss) ngauss = 0
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
             &        "ngauss,degauss=",i4,f12.6," Ry", /)') ngauss, degauss
     ltetra=.false.
     lgauss=.true.
  ELSEIF (ltetra) THEN
     WRITE( stdout,'(/5x,"Tetrahedra used"/)')
  ELSEIF (legauss) THEN
     ngauss=0
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6, " Ry",/)') ngauss,degauss
     ltetra=.FALSE.
     lgauss=.TRUE.
  ELSE
     ngauss=0
     degauss=0.01_DP
     WRITE( stdout,'(/5x,"Broadening (default values): ",&
             &        "ngauss,degauss=",i4,f12.6," Ry"/)') ngauss,degauss
  ENDIF
  !
  ! find min and max energy for plot (band extrema if not set)
  ! bands are in Ry as emin and emax
  !
  emin = MINVAL ( et(1, 1:nks) )
  IF ( degauss > 0.0_dp ) emin = emin - 5.0_dp * degauss

  IF ( emax_input  == 0.0_DP ) THEN
     emax = MAXVAL ( et(nbnd, 1:nks) )
     IF ( degauss > 0.0_dp ) emax = emax + 3.0_dp * degauss
  ELSE 
!
!   The input emax_input given in eV is converted in Ry
!
     emax = emax_input/rytoev
  END IF
  !
  IF (deltae > 0.0_DP) THEN
     ndos = NINT ( (emax - emin) / deltae)
  ELSEIF (ndose > 0 ) THEN
     ndos = ndose
     deltae = (emax-emin) / (ndos - 1)
  ELSE
     deltae=0.01_DP
     ndos = NINT ( (emax - emin) / deltae)
  ENDIF
  WRITE(stdout,'(/,5x,"Delta e=", f15.8, " Ry, ndos= ", i9)') deltae, ndos

  save_ndos=ndos
  ALLOCATE(e(ndos))
  ALLOCATE(dosofe(2,ndos))
  ALLOCATE(dosint(ndos))

  dosofe=0.0_DP
  e=0.0_DP
  !
  WRITE(stdout,*)

  CALL divide (intra_image_comm, ndos, nstart, nlast)
  DO n= nstart, nlast
     IF (MOD(n,5000)==0) WRITE(stdout,&
                     '(5x,"Computing point number ", i10," /",i10  &
                           & ," Total ", i10)') n, nlast-nstart+1, ndos
     e(n) = emin + (n - 1) * deltae
     IF (ltetra) THEN
        CALL tetra_dos_t(et,nspin,nbnd,nks,e(n),dosofe(1,n))
     ELSE
!
!   note that here even with pools nks is the total number of k points
!   and et and wk contains the bands and the weights of all k points
!   This is because the pw variables have been cleaned and reread at
!   the beginning.
!
        CALL dos_g(et,nspin,nbnd,nks,wk,degauss,ngauss,e(n),dosofe(1,n))
     ENDIF
  ENDDO
!
!  sum over all processors of this image
!
  CALL mp_sum(e, intra_image_comm)
  CALL mp_sum(dosofe, intra_image_comm)
!
!  and integrate the dos to get the number of states 
!
  IF (nspin==2) THEN
     dosint(1) = (dosofe (1,1) + dosofe (2,1)) * deltae
     DO n=2, ndos
        dosint(n)=dosint(n-1) + (dosofe(1,n) + dosofe(2,n)) * deltae
     ENDDO
  ELSE
     dosint(1) = dosofe (1,1) * deltae
     DO n=2, ndos
        dosint(n)=dosint(n-1) + dosofe(1,n) * deltae
     ENDDO
  ENDIF

  IF (ltetra .OR. lgauss) THEN
!
!  Recalculate the Fermi level using the integrated dos
!
     WRITE(stdout,'(/,5x, "Fermi energy from nscf run", f19.8, " eV")')  &
                                                               ef * rytoev
     ef1=0.0_DP
     DO n=1,ndos
        IF (dosint(n)<=nelec) ef1=e(n)
     ENDDO

     WRITE(stdout,'(/,5x, "Fermi energy from integrated dos", &
              & f13.8, " eV +/- ",f13.8," eV")')  ef1 * rytoev, deltae*rytoev
  ENDIF
  !
  !  Write the results on file
  !
  IF ( ionode ) THEN

     iu_dos=find_free_unit()
     filedos='therm_files/'//TRIM(fleldos)
     OPEN (unit=iu_dos, file=TRIM(filedos), status='unknown', form='formatted')

     IF (nspin==1 .OR. nspin==4) THEN
        WRITE(iu_dos,'("#  E (eV)   dos(E)     Int dos(E)")')
     ELSE
        WRITE(iu_dos,'("#  E (eV)   dosup(E)     dosdw(E)   Int dos(E)")')
     ENDIF

     DO n= 1, ndos
        IF (nspin==1 .OR. nspin==4) THEN
           WRITE (iu_dos, '(f15.9,2e16.8)') e(n)*rytoev, dosofe(1,n)/rytoev,&
                                            dosint(n)
        ELSE
           WRITE (iu_dos, '(f15.9,3e16.8)') e(n)*rytoev, dosofe(:,n)/rytoev,&
                                            dosint(n)
        ENDIF
     ENDDO

     CLOSE (unit = iu_dos)
     !
  ENDIF

  degauss=save_degauss
  ngauss=save_ngauss
  ltetra=save_ltetra
  lgauss=save_lgauss

  DEALLOCATE(dosint)
  DEALLOCATE(dosofe)
  DEALLOCATE(e)
  !
  RETURN
END SUBROUTINE dos_sub
