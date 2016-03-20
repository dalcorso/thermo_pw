!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE phescf_tpw()
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver for the calculation of the
  ! ... response to an electric field and related quantities.
  !
  USE io_global,       ONLY : stdout
  USE paw_variables,   ONLY : okpaw
  USE klist,           ONLY : lgauss
  USE uspp,            ONLY : okvan
  USE uspp_param,      ONLY : nhm
  USE ions_base,       ONLY : nat
  USE noncollin_module,ONLY : noncolin, nspin_mag, npol
  USE lsda_mod,        ONLY : nspin
  USE control_ph,      ONLY : convt, zeu, rec_code, rec_code_read, lnoloc, &
                              lrpa, where_rec, done_epsil, done_zeu, epsil
  USE io_files,        ONLY : tmp_dir
  USE wvfct,           ONLY : nbnd, npwx
  USE control_flags,   ONLY : io_level
  USE output,          ONLY : fildrho
  USE ph_restart,      ONLY : ph_writefile
  USE phus,            ONLY : int3, int3_nc, int3_paw
  USE freq_ph
  USE optical,         ONLY : current_w, fru, polarc, epsilonc, epsilonm1c, &
                              lr1dwf, iu1dwf, lcfreq
  USE ramanm,          ONLY : ramtns, lraman, elop, done_lraman, done_elop
  USE buffers,         ONLY : close_buffer, open_buffer
  !
  IMPLICIT NONE
  !
  INTEGER :: iu, ipol, jpol, iu_epsil, ierr
  !
  LOGICAL :: exst_mem, exst
  !
  IF ( rec_code_read >  1 ) THEN
     IF (done_epsil) call summarize_epsilon()
     IF (done_zeu) call summarize_zeu()
     IF (done_elop) call summarize_elopt()
     IF (done_lraman) call write_ramtns(6,ramtns)
     RETURN
  ENDIF
  !
  IF (okvan) THEN
     ALLOCATE (int3 ( nhm, nhm, 3, nat, nspin_mag))
     IF (okpaw) ALLOCATE (int3_paw ( nhm, nhm, 3, nat, nspin_mag))
     IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, 3, nat, nspin))
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
     IF (lcfreq) THEN
        iu1dwf = 42
        lr1dwf =  nbnd * npwx * npol
        CALL open_buffer (iu1dwf, 'mwf', lr1dwf, io_level, exst_mem, &
                                                           exst, tmp_dir)
     ENDIF
     !
     epsilonm1c = (0.0_DP,0.0_DP)
     DO iu = 1, nfs
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
              END IF
           END DO
        END DO

        CALL write_epsilon_on_disk(iu)

     END DO 

     IF (lcfreq) CALL close_buffer(iu1dwf, 'delete')

     DEALLOCATE(polarc)
     DEALLOCATE(epsilonc)
     DEALLOCATE(epsilonm1c)
     !
     WRITE( stdout, '(/,5X,"End of Frequency Dependent Polarizability Calculation")' )
     !
  ENDIF
  !
  IF (((epsil.AND..NOT.done_epsil).OR.(zeu.AND..NOT.done_zeu).OR.  &
      (lraman.AND..NOT.done_lraman).OR.(elop.AND..NOT.done_elop))  &
                                         .AND..NOT.lgauss) THEN

     WRITE( stdout, '(/,5X,"Electric Fields Calculation")' )
     !

     CALL solve_e()
     !
     WRITE( stdout, '(/,5X,"End of electric fields calculation")' )
     !
     IF ( convt ) THEN
        !
        ! ... calculate the dielectric tensor epsilon
        !
        IF (.NOT. done_epsil) THEN
           CALL dielec()
        ELSE
           CALL summarize_epsilon()
        ENDIF
       !
       ! ... calculate the effective charges Z(E,Us) (E=scf,Us=bare)
       !
       IF (.NOT.(lrpa.OR.lnoloc).AND.(zeu.AND..NOT.done_zeu)) THEN
           CALL zstar_eu()
        ELSEIF (done_zeu) THEN
           CALL summarize_zeu()
        ENDIF
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
  !
  IF (okvan) THEN
     DEALLOCATE (int3)
     IF (okpaw) DEALLOCATE (int3_paw)
     IF (noncolin) DEALLOCATE(int3_nc)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE phescf_tpw


SUBROUTINE write_epsilon_on_disk(iu)
USE kinds, ONLY : DP
USE io_global, ONLY : meta_ionode
USE klist, ONLY : nkstot
USE lsda_mod, ONLY : lsda
USE freq_ph,  ONLY : fiu
USE optical, ONLY : fru, epsilonc, epsilonm1c, polarc

IMPLICIT NONE
INTEGER, INTENT(IN) :: iu
INTEGER :: iu_epsil
LOGICAL :: exst

IF (meta_ionode) THEN

   iu_epsil=2
   INQUIRE(FILE="epsilon_re", exist=exst)
   IF (exst) THEN
      OPEN (UNIT=iu_epsil, FILE='epsilon_re', STATUS='old', &
                              POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE='epsilon_re', STATUS='unknown', &
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
   INQUIRE(FILE="epsilon_im", exist=exst)
   IF (exst) THEN
      OPEN (UNIT=iu_epsil, FILE='epsilon_im', STATUS='old', &
                              POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE='epsilon_im', STATUS='unknown', &
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
   INQUIRE(FILE="epsilonm1_re", exist=exst)
   IF (exst) THEN
      OPEN (UNIT=iu_epsil, FILE='epsilonm1_re', STATUS='old', &
                              POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE='epsilonm1_re', STATUS='unknown', &
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
   INQUIRE(FILE="epsilonm1_im", exist=exst)
   IF (exst) THEN
      OPEN (UNIT=iu_epsil, FILE='epsilonm1_im', STATUS='old', &
                           POSITION='append', FORM='formatted')
   ELSE
      OPEN (UNIT=iu_epsil, FILE='epsilonm1_im', STATUS='unknown', &
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
      INQUIRE(FILE="polariz_re", exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE='polariz_re', STATUS='old', &
                                 POSITION='append', FORM='formatted')
      ELSE
         OPEN (UNIT=iu_epsil, FILE='polariz_re', STATUS='unknown', &
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
      INQUIRE(FILE="polariz_im", exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE='polariz_im', STATUS='old', &
                                 POSITION='append', FORM='formatted')
      ELSE
         OPEN (UNIT=iu_epsil, FILE='polariz_im', STATUS='unknown', &
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
