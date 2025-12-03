!
! Copyright (C) 2009-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE phescf_tpw()
  !-----------------------------------------------------------------------
  !
  !! This is the main driver for the calculation of the
  !! response to an electric field and related quantities.
  !
  USE io_global,       ONLY : stdout
  USE paw_variables,   ONLY : okpaw
  USE klist,           ONLY : lgauss
  USE uspp,            ONLY : okvan
  USE uspp_param,      ONLY : nhm
  USE fft_base,        ONLY : dffts
  USE ions_base,       ONLY : nat
  USE noncollin_module,ONLY : noncolin, nspin_mag, npol, domag
  USE lsda_mod,        ONLY : nspin
  USE control_ph,      ONLY : zeu, lnoloc, done_epsil, done_zeu, epsil
  USE eqv,             ONLY : drhos
  USE io_files,        ONLY : tmp_dir
  USE wvfct,           ONLY : nbnd, npwx
  USE control_flags,   ONLY : io_level
  USE output,          ONLY : fildrho
  USE ph_restart,      ONLY : ph_writefile
  USE freq_ph
  USE optical,         ONLY : current_w, fru, polarc, epsilonc, epsilonm1c, &
                              lr1dwf, iu1dwf, lcfreq, start_freq, last_freq
  USE lr_lanczos,      ONLY : llanczos, iulanczos, only_spectrum
  USE lr_global,       ONLY : pseudo_hermitian
  USE lr_cg,           ONLY : lcg
  USE partial,         ONLY : comp_irr
  USE images_omega,    ONLY : comp_f
  USE ramanm,          ONLY : ramtns, lraman, elop, done_lraman, done_elop
  USE lrus,            ONLY : int3, int3_nc, int3_paw
  USE control_lr,      ONLY : convt, lrpa, rec_code, rec_code_read, where_rec
  USE buffers,         ONLY : close_buffer, open_buffer
  USE io_global,       ONLY : ionode, ionode_id
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_bcast
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax
  USE ldaU_lr,         ONLY : dnsscf
  USE ldaU_ph,         ONLY : dnsscf_all_modes
  USE control_flags,   ONLY : iverbosity
  USE magnetic_charges, ONLY : alpha_me
  USE control_lr,       ONLY : lgamma
  USE write_hub

  !
  IMPLICIT NONE
  !
  INTEGER :: iu, ipol, jpol, iu_epsil, ierr, ios
  !
  INTEGER :: find_free_unit
  !
  INTEGER, ALLOCATABLE :: computed(:)
  !
  LOGICAL :: exst_mem, exst
  !
  IF ( .NOT. comp_irr(0)  ) RETURN

  IF ( rec_code_read >  1 ) THEN
     IF (done_epsil) call summarize_epsilon()
     IF (done_zeu) call summarize_zeu()
     IF (done_elop) call summarize_elopt()
     IF (done_lraman) call write_ramtns(6,ramtns)
     RETURN
  ENDIF
  !
  IF (okvan) THEN
     ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, 3))
     IF (okpaw) ALLOCATE (int3_paw ( nhm, nhm, nat, nspin_mag, 3))
     IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, nat, nspin, 3))
  ENDIF
  ! Set symmetry representation in lr_symm_base
  !
  CALL ph_set_upert_e()
  !
  ALLOCATE (drhos( dffts%nnr, nspin_mag, 3))
  !
  ! DFPT+U: dnsscf in the electric field calculation
  ! is the scf change of atomic occupations ns induced by the electric field.
  ! dnsscf_all_modes = dnsscf because nirr=1, number of perturbations = 3.
  !
  IF (lda_plus_u) THEN
     ALLOCATE (dnsscf(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,3))
     ALLOCATE (dnsscf_all_modes(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,3))
     dnsscf = (0.d0, 0.d0)
     dnsscf_all_modes = (0.d0, 0.d0)
  ENDIF
  !
  IF (fpol) THEN    ! calculate freq. dependent polarizability
     !
     WRITE( stdout, '(/,5X,"Frequency Dependent TD-DFPT Dielectric constants",/)' )
     IF (lnoloc) &
        WRITE( stdout, '(/,5X,"lnoloc=.true.: Independent electrons",/)' )

     IF (lrpa) &
        WRITE( stdout, '(/,5X,"lrpa=.true.: RPA with local fields",/)' )
     !
     ALLOCATE(epsilonc(3,3,nfs))
     ALLOCATE(epsilonm1c(3,3,nfs))
     ALLOCATE(polarc(3,3,nfs))
     ALLOCATE(computed(nfs))
     IF (lcfreq) THEN
        iu1dwf = 42
        lr1dwf =  nbnd * npwx * npol
        CALL open_buffer (iu1dwf, 'mwf', lr1dwf, io_level, exst_mem, &
                                                           exst, tmp_dir)
     ENDIF
     !
     computed=0
     IF (llanczos) THEN
        CALL allocate_lanczos()
        IF (ionode) THEN
           OPEN(UNIT=99, FILE='buffer') 
           iulanczos=find_free_unit()
           CLOSE(UNIT=99, STATUS='delete') 
           OPEN(UNIT=iulanczos, FILE='dynamical_matrices/save_chain',&
            STATUS='unknown', POSITION='append', ERR=100, IOSTAT=ios)
        ENDIF
100     CALL mp_bcast(ios,ionode_id,intra_image_comm)
        CALL errore('phescf','opening save_chain file',ABS(ios))
        IF (only_spectrum) THEN
           CALL read_lanczos_chain()
        ELSE
           IF (pseudo_hermitian) THEN
              CALL do_lanczos_psh()
           ELSE
              CALL do_lanczos()
           ENDIF
        ENDIF
        CALL extrapolate()
        IF (ionode) THEN
           epsilonm1c = (0.0_DP,0.0_DP)
           DO iu=start_freq, last_freq
              computed(iu)=1
              CALL calc_chi(fru(iu),fiu(iu),epsilonc(1,1,iu))
              DO ipol=1,3
                 DO jpol=1,3
                    IF (ABS(epsilonc(ipol,jpol,iu))>1.D-8) THEN
                       epsilonm1c(ipol,jpol,iu)=CMPLX(1.0_DP,0.0_DP) / &
                        epsilonc(ipol,jpol,iu)
                    END IF
                 END DO
              END DO
              CALL write_epsilon_on_disk(iu)
           ENDDO
           CLOSE(UNIT=iulanczos, STATUS='keep')
        ENDIF
        CALL deallocate_lanczos()
     ELSE
        epsilonm1c = (0.0_DP,0.0_DP)
        DO iu = start_freq, last_freq
           !
           IF (.NOT.comp_f(iu)) CYCLE
           computed(iu)=1
           !
           current_w=CMPLX(fru(iu), fiu(iu))
           WRITE( stdout, '(/,5x,70("-"))') 
           WRITE( stdout, '(6x,"Dielectric constant at &
                &frequency",f9.4," +",f9.4," i Ry", i6," /",i6)') current_w, &
                iu, nfs
           !
           CALL solve_e_fpolc( iu )

           IF ( convt ) CALL polarizc (iu)
           !
           !  Save also the inverse of the macroscopic dielectric constant
           !  and its square root, the refractive_index
           !
           DO ipol=1,3
              DO jpol=1,3
                 IF (ABS(epsilonc(ipol,jpol,iu))>1.D-8) THEN
                    epsilonm1c(ipol,jpol,iu)=CMPLX(1.0_DP,0.0_DP) / &
                        epsilonc(ipol,jpol,iu)
                 ENDIF
              ENDDO
           ENDDO

           CALL write_epsilon_on_disk(iu)

        ENDDO 
     ENDIF
     CALL write_optical_on_disk(epsilonc, fru, fiu, computed, nfs)

     IF (lcfreq) CALL close_buffer(iu1dwf, 'delete')

     DEALLOCATE(polarc)
     DEALLOCATE(epsilonc)
     DEALLOCATE(epsilonm1c)
     DEALLOCATE(computed)
     !
     WRITE( stdout, '(/,5X,"End of Frequency Dependent Polarizability Calculation")' )
     !
  ELSE
  !
     IF (((epsil.AND..NOT.done_epsil).OR.(zeu.AND..NOT.done_zeu).OR.  &
         (lraman.AND..NOT.done_lraman).OR.(elop.AND..NOT.done_elop))  &
                                         .AND..NOT.lgauss) THEN

        WRITE( stdout, '(/,5X,"Electric Fields Calculation")' )
        !
        IF (lcg) THEN
           CALL allocate_cg(3,1)
           CALL do_cg_e(drhos)
           CALL deallocate_cg()
        ELSE
           CALL solve_e_tpw(drhos)
        ENDIF
        !
        WRITE( stdout, '(/,5X,"End of electric fields calculation")' )
        !
        IF ( convt ) THEN
           !
           ! ... calculate the dielectric tensor epsilon
           !
           IF (.NOT. done_epsil) THEN
              CALL dielec_tpw()
           ELSE
              CALL summarize_epsilon()
           ENDIF
           !
           ! ... calculate the effective charges Z(E,Us) (E=scf,Us=bare)
           !
           IF (.NOT.(lrpa.OR.lnoloc).AND.(zeu.AND..NOT.done_zeu)) THEN
              CALL zstar_eu_tpw(drhos)
           ELSEIF (done_zeu) THEN
              CALL summarize_zeu()
           ENDIF
           ! 
           ! ... write magnetoelectric tensor 
           !
           IF (noncolin.AND.domag.AND.lgamma) CALL summarize_alpha()
           !
           IF ( fildrho /= ' ' ) CALL punch_plot_e()
           !
        ELSE
           !
           CALL stop_ph( .FALSE. )
           !
        END IF
        !
        IF ( (lraman.AND..NOT.done_lraman) .OR. (elop.AND..NOT.done_elop) &
                     .AND..NOT.noncolin) CALL raman()
        !
        where_rec='after_diel'
        rec_code=2
        CALL ph_writefile('status_ph',0,0,ierr)
     ELSE
        IF (done_epsil) call summarize_epsilon()
        IF (done_zeu) call summarize_zeu()
        IF (done_elop) call summarize_elopt()
        IF (done_lraman) call write_ramtns(6,ramtns)
     ENDIF
  ENDIF
  !
  IF (okvan) THEN
     DEALLOCATE (int3)
     IF (okpaw) DEALLOCATE (int3_paw)
     IF (noncolin) DEALLOCATE(int3_nc)
  ENDIF

  DEALLOCATE (drhos)
  !
  CALL ph_deallocate_upert()
  !
  ! DFPT+U
  !
  IF (lda_plus_u) THEN
     !
     ! Write dnsscf_all_modes in the cartesian coordinates
     ! to the standard output
     !
     IF (iverbosity==1) CALL write_dnsscf_e()
     !
     DEALLOCATE (dnsscf)
     DEALLOCATE (dnsscf_all_modes)
     !
  ENDIF

  RETURN
  !
END SUBROUTINE phescf_tpw

!-----------------------------------------------------------------
SUBROUTINE write_epsilon_on_disk(iu)
!-----------------------------------------------------------------
USE kinds,    ONLY : DP
USE io_global, ONLY : ionode
USE klist,    ONLY : nkstot
USE lsda_mod, ONLY : lsda
USE freq_ph,  ONLY : fiu
USE optical,  ONLY : fru, epsilonc, epsilonm1c, polarc
USE mp_images, ONLY : my_image_id

IMPLICIT NONE
INTEGER, INTENT(IN) :: iu
INTEGER :: iu_epsil
LOGICAL :: exst
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6)   :: int_to_char

IF (ionode) THEN

   iu_epsil=2
   filename='dynamical_matrices/epsilon_re'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   INQUIRE(FILE=TRIM(filename), exist=exst)
   IF (exst.AND.iu>1) THEN
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old',&
                              POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                FORM='formatted')
      WRITE(iu_epsil,'("#  Re(w)     Im(w)    e11            e22&
              &            e33            e12            e13            e23")')
   END IF

   WRITE(iu_epsil,'(2f10.5,6e15.7)') fru(iu), fiu(iu),  &
                DREAL(epsilonc(1,1,iu)), DREAL(epsilonc(2,2,iu)), &
                DREAL(epsilonc(3,3,iu)), DREAL(epsilonc(1,2,iu)), &
                DREAL(epsilonc(1,3,iu)), DREAL(epsilonc(2,3,iu) )

   CLOSE(iu_epsil)

   iu_epsil=2
   filename='dynamical_matrices/epsilon_im'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   INQUIRE(FILE=TRIM(filename), exist=exst)
   IF (exst.AND.iu>1) THEN
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old',&
                             POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                FORM='formatted')
      WRITE(iu_epsil,'("#  Re(w)     Im(w)    e11            e22&
              &            e33            e12            e13            e23")')
   END IF

   WRITE(iu_epsil,'(2f10.5,6e15.7)') fru(iu), fiu(iu), &
                DIMAG(epsilonc(1,1,iu)), DIMAG(epsilonc(2,2,iu)), &
                DIMAG(epsilonc(3,3,iu)), DIMAG(epsilonc(1,2,iu)), &
                DIMAG(epsilonc(1,3,iu)), DIMAG(epsilonc(2,3,iu))

   CLOSE(iu_epsil)

   iu_epsil=2
   filename='dynamical_matrices/epsilonm1_re'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   INQUIRE(FILE=TRIM(filename), exist=exst)
   IF (exst.AND.iu>1) THEN
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                       POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                          FORM='formatted')
      WRITE(iu_epsil,'("#  Re(w)     Im(w)  em111          em122&
              &          em133          em112          em113          em123")')
   END IF
   WRITE(iu_epsil,'(2f10.5,6e15.7)') fru(iu), fiu(iu), &
                DREAL(epsilonm1c(1,1,iu)), DREAL(epsilonm1c(2,2,iu)), &
                DREAL(epsilonm1c(3,3,iu)), DREAL(epsilonm1c(1,2,iu)), &
                DREAL(epsilonm1c(1,3,iu)), DREAL(epsilonm1c(2,3,iu))
   CLOSE(iu_epsil)

   iu_epsil=2
   filename='dynamical_matrices/epsilonm1_im'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   INQUIRE(FILE=TRIM(filename), exist=exst)
   IF (exst.AND.iu>1) THEN
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                      POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                FORM='formatted')
      WRITE(iu_epsil,'("#  Re(w)     Im(w)  em111          em122&
              &          em133          em112          em113          em123")')
   END IF

   WRITE(iu_epsil,'(2f10.5,6e15.7)') fru(iu), fiu(iu), &
              DIMAG(epsilonm1c(1,1,iu)),DIMAG(epsilonm1c(2,2,iu)),&
              DIMAG(epsilonm1c(3,3,iu)),DIMAG(epsilonm1c(1,2,iu)),&
              DIMAG(epsilonm1c(1,3,iu)),DIMAG(epsilonm1c(2,3,iu))

   CLOSE(iu_epsil)

   IF (nkstot==1 .OR. (nkstot==2.AND.lsda)) THEN
!
!   In the molecular case write also the polarization
!
      iu_epsil=2
      filename='dynamical_matrices/polariz_re'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst.AND.iu>1) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                              POSITION='append', FORM='formatted')
      ELSE
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                   FORM='formatted')
         WRITE(iu_epsil,'("#  Re(w)     Im(w)    a11            a22&
              &            a33            a12            a13            a23")')
      END IF

      WRITE(iu_epsil,'(2f10.5,6e15.7)') fru(iu), fiu(iu),  &
                   DREAL(polarc(1,1,iu)), DREAL(polarc(2,2,iu)), &
                   DREAL(polarc(3,3,iu)), DREAL(polarc(1,2,iu)), &
                   DREAL(polarc(1,3,iu)), DREAL(polarc(2,3,iu) )

      CLOSE(iu_epsil)

      iu_epsil=2
      filename='dynamical_matrices/polariz_im'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst.AND.iu>1) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                              POSITION='append', FORM='formatted')
      ELSE
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                   FORM='formatted')
         WRITE(iu_epsil,'("#  Re(w)     Im(w)    a11            a22&
              &            a33            a12            a13            a23")')
      END IF

      WRITE(iu_epsil,'(2f10.5,6e15.7)') fru(iu), fiu(iu), &
                   DIMAG(polarc(1,1,iu)), DIMAG(polarc(2,2,iu)), &
                   DIMAG(polarc(3,3,iu)), DIMAG(polarc(1,2,iu)), &
                   DIMAG(polarc(1,3,iu)), DIMAG(polarc(2,3,iu))

      CLOSE(iu_epsil)
   END IF
END IF

RETURN
END SUBROUTINE write_epsilon_on_disk

!-----------------------------------------------------------
SUBROUTINE collect_all_epsilon()
!-----------------------------------------------------------
!
!   This subroutine reads from disk all the dielectric constants written
!   by the different images and writes them on a single file
!
USE kinds,     ONLY : DP
USE constants, ONLY : fpi
USE klist,     ONLY : nkstot
USE cell_base, ONLY : omega
USE lsda_mod,  ONLY : lsda
USE freq_ph,   ONLY : fiu, nfs
USE optical,   ONLY : fru, epsilonc, epsilonm1c, polarc

USE io_global, ONLY : meta_ionode, ionode, ionode_id
USE mp_images, ONLY : my_image_id, intra_image_comm, inter_image_comm, nimage
USE mp,        ONLY : mp_bcast, mp_sum

IMPLICIT NONE
INTEGER  :: iu_epsil, ipol, jpol, iuf, iu
LOGICAL  :: exst
REAL(DP) :: buffer(3,3), alpha(3,3), epsilonr(3,3,nfs), epsiloni(3,3,nfs),&
            ro, ri
INTEGER, ALLOCATABLE :: computed(:)
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6)   :: int_to_char

ALLOCATE(epsilonc(3,3,nfs))
ALLOCATE(epsilonm1c(3,3,nfs))
ALLOCATE(polarc(3,3,nfs))
ALLOCATE(computed(nfs))


computed=0
IF (ionode) THEN

   iu_epsil=2
   filename='dynamical_matrices/epsilon_re'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   INQUIRE(FILE=TRIM(filename), exist=exst)
   IF (exst) THEN
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      READ(iu_epsil, *) 
   END IF

   epsilonr=0.0_DP
   DO iuf=1,nfs
      buffer=0.0_DP
      READ(iu_epsil,'(2f10.5,6e15.7)',END=100) ro, ri, buffer(1,1), &
           buffer(2,2), buffer(3,3), buffer(1,2), buffer(1,3), buffer(2,3)
      !
      !   Now search which frequency is this
      !
      DO iu=1,nfs
         IF ((ABS(ro-fru(iu))+ABS(ri-fiu(iu)))<1.D-4) THEN
            epsilonr(:,:,iu)=buffer(:,:)
            computed(iu)=1
            EXIT
         ENDIF
      ENDDO       
   ENDDO
100 CONTINUE

   CLOSE(iu_epsil)

   iu_epsil=2
   filename='dynamical_matrices/epsilon_im'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   INQUIRE(FILE=TRIM(filename), exist=exst)
   IF (exst) THEN
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      READ(iu_epsil, *) 
   END IF

   epsiloni=0.0_DP
   DO iuf=1,nfs
      READ(iu_epsil,'(2f10.5,6e15.7)',END=200) ro, ri, buffer(1,1), &
                buffer(2,2), buffer(3,3), buffer(1,2), buffer(1,3), buffer(2,3)
      !
      !   Now search which frequency is this
      !
      DO iu=1,nfs
         IF ((ABS(ro-fru(iu))+ABS(ri-fiu(iu)))<1.D-4) THEN
            epsiloni(:,:,iu)=buffer(:,:)
            computed(iu)=1
            EXIT
         ENDIF
      ENDDO
   ENDDO
200 CONTINUE
   CLOSE(iu_epsil)
!
!   Finally set the dynamical matrix with all frequencies
!
   epsilonc(:,:,:)=CMPLX(epsilonr(:,:,:), epsiloni(:,:,:))
   epsilonc(2,1,:)=epsilonc(1,2,:)
   epsilonc(3,1,:)=epsilonc(1,3,:)
   epsilonc(3,2,:)=epsilonc(2,3,:)
ENDIF
!
! Now transfer the dielectric constant to all nodes. First collect
! all frequencies and then send to all the nodes of each image
!
CALL mp_sum(computed, inter_image_comm) 
CALL mp_bcast(computed, ionode_id, intra_image_comm)
CALL mp_sum(epsilonc, inter_image_comm) 
CALL mp_bcast(epsilonc, ionode_id, intra_image_comm)
DO iu=1,nfs
   IF (computed(iu)>1) epsilonc(:,:,iu) = epsilonc(:,:,iu) / computed(iu)
END DO
!
!  compute epsilonm1 and polarc
!
DO iu=1,nfs
   DO ipol=1,3
      DO jpol=1,3
         IF (ABS(epsilonc(ipol,jpol,iu))>1.D-8) THEN
             epsilonm1c(ipol,jpol,iu)=CMPLX(1.0_DP,0.0_DP) / &
                     epsilonc(ipol,jpol,iu)
         END IF
      END DO
   END DO
END DO

polarc=(0.0_DP,0.0_DP)
IF (nkstot==1 .OR. (nkstot==2.AND.lsda)) THEN
!
!   In the molecular case write also the polarization
!
   DO iu=1,nfs
      alpha=(0.0_DP, 0.0_DP)
      DO ipol = 1, 3
         DO jpol = 1, 3
            IF (ABS(epsilonc(ipol,jpol,iu)) > 1.D-4) &
                  alpha(ipol, jpol)=(3.d0*omega/fpi)*(epsilonc(ipol,jpol,iu) &
                            - 1.0_DP )/(epsilonc(ipol,jpol,iu)+2.0_DP)
         ENDDO
      ENDDO
      polarc(:,:,iu)=alpha(:,:)
   ENDDO
ENDIF
!
iu_epsil=2
!
!  First each image deletes the files with the dielectric constant that it
!  has produced
!
IF (ionode) THEN
   filename='dynamical_matrices/epsilon_re'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
   CLOSE (UNIT=iu_epsil, STATUS='DELETE')

   filename='dynamical_matrices/epsilon_im'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
   CLOSE (UNIT=iu_epsil, STATUS='DELETE')

   filename='dynamical_matrices/epsilonm1_re'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
   CLOSE (UNIT=iu_epsil, STATUS='DELETE')

   filename='dynamical_matrices/epsilonm1_im'
   IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
   OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
   CLOSE (UNIT=iu_epsil, STATUS='DELETE')

   IF (nkstot==1 .OR. (nkstot==2.AND.lsda)) THEN
      filename='dynamical_matrices/polariz_re'
      IF (my_image_id>0) filename=TRIM(filename)//'_'// &
                                                     int_to_char(my_image_id)
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      CLOSE (UNIT=iu_epsil, STATUS='DELETE')

      filename='dynamical_matrices/polariz_im'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      CLOSE (UNIT=iu_epsil, STATUS='DELETE')
   ENDIF
ENDIF
!
!  now rewrite on a single file all the collected dielectric constants.
!  Only one processor writes on disk
!
IF (meta_ionode) THEN
   DO iu=1,nfs
      IF (computed(iu)>0) CALL write_epsilon_on_disk(iu)
   END DO
   CALL write_optical_on_disk(epsilonc, fru, fiu, computed, nfs)
ENDIF

DEALLOCATE(polarc)
DEALLOCATE(epsilonc)
DEALLOCATE(epsilonm1c)
DEALLOCATE(computed)

RETURN
END SUBROUTINE collect_all_epsilon
!
!-------------------------------------------------------------------
SUBROUTINE write_optical_on_disk(epsilonc, fru, fiu, computed, nfs)
!-------------------------------------------------------------------
!
!  This routine writes on disk the optical quantities:
!  epsilon1, epsilon2  real and imaginary part of the dielectric constant
!  n, k                real and imaginary part of the refractive index
!  For cubic solids it writes also
!  R, alpha            the reflectivity for normal incidence and the &
!                      absorption coefficient
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE optical_module, ONLY : compute_refractive_index, compute_lambda, &
                           compute_reflectivity, compute_alpha
IMPLICIT NONE
INTEGER, INTENT(IN)     :: nfs
COMPLEX(DP), INTENT(IN) :: epsilonc(3,3,nfs)
REAL(DP), INTENT(IN)    :: fru(nfs), fiu(nfs)
INTEGER, INTENT(IN)     :: computed(nfs)

CHARACTER(LEN=256) :: filename
INTEGER  :: find_free_unit
INTEGER  :: iu_optic, iu
REAL(DP) :: epsilon1, epsilon2, enne, kappa, ref, lambda, alpha
LOGICAL  :: cubic

iu_optic=find_free_unit()
filename='dynamical_matrices/optical_xx'
OPEN (UNIT=iu_optic, FILE=TRIM(filename), STATUS='unknown', FORM='formatted')
WRITE(iu_optic,'("# Im(w)=",e15.7," Ry")') fiu(1)
cubic=(ABS(ibrav)>0.AND.ibrav<4)
IF (cubic) THEN
   WRITE(iu_optic,'("# Re(w) (Ry)      lambda         e1             e2&
         &          n              k              R         alpha (cm^-1)")')
ELSE
   WRITE(iu_optic,'("# Re(w) (Ry)      lambda         e1             e2&
         &          n              k")')
ENDIF
DO iu=1, nfs
   IF (computed(iu)>0) THEN
      epsilon1=DBLE(epsilonc(1,1,iu))
      epsilon2=AIMAG(epsilonc(1,1,iu))
      CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
      CALL compute_lambda(fru(iu), lambda)
      IF (ABS(ibrav)>0.AND.ibrav<4) THEN
         CALL compute_alpha(fru(iu), kappa, alpha)
         CALL compute_reflectivity(enne, kappa, ref)
         WRITE(iu_optic,'(8e15.7)') fru(iu), lambda, epsilon1, epsilon2, &
                                             enne, kappa, ref, alpha
      ELSE
         WRITE(iu_optic,'(6e15.7)') fru(iu), lambda, epsilon1, epsilon2, &
                                             enne, kappa
      ENDIF
   ENDIF
ENDDO
CLOSE(iu_optic)

IF (ibrav>3.AND.ABS(ibrav)<12) THEN 
   filename='dynamical_matrices/optical_zz'
   OPEN (UNIT=iu_optic, FILE=TRIM(filename), STATUS='unknown', &
                                                     FORM='formatted')
   WRITE(iu_optic,'("# Im(w)=",e15.7," Ry")') fiu(1)
   WRITE(iu_optic,'("# Re(w) (Ry)      lambda         e1             e2&
         &          n              k   ")')
   DO iu=1, nfs
      IF (computed(iu)>0) THEN
         epsilon1=DBLE(epsilonc(3,3,iu))
         epsilon2=AIMAG(epsilonc(3,3,iu))
         CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
         CALL compute_lambda(fru(iu), lambda)
         WRITE(iu_optic,'(6e15.7)') fru(iu), lambda, epsilon1, epsilon2, &
                                             enne, kappa
      ENDIF
   ENDDO
   CLOSE(iu_optic)
ENDIF

IF (ibrav>7.AND.ABS(ibrav)<12) THEN 
   filename='dynamical_matrices/optical_yy'
   OPEN (UNIT=iu_optic, FILE=TRIM(filename), STATUS='unknown', &
                                                     FORM='formatted')
   WRITE(iu_optic,'("# Im(w)=",e15.7," Ry")') fiu(1)
   WRITE(iu_optic,'("# Re(w) (Ry)      lambda         e1             e2&
         &          n              k          ")')
   DO iu=1, nfs
      IF (computed(iu)>0) THEN
         epsilon1=DBLE(epsilonc(2,2,iu))
         epsilon2=AIMAG(epsilonc(2,2,iu))
         CALL compute_refractive_index(epsilon1, epsilon2, enne, kappa)
         CALL compute_lambda(fru(iu), lambda)
         WRITE(iu_optic,'(8e15.7)') fru(iu), lambda, epsilon1, epsilon2, &
                                             enne, kappa
      ENDIF
   ENDDO
   CLOSE(iu_optic)
ENDIF

RETURN
END SUBROUTINE write_optical_on_disk
