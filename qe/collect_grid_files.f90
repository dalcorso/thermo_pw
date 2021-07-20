!
! Copyright (C) 2012-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
   SUBROUTINE collect_grid_files_tpw()
!-----------------------------------------------------------------------------
   !
   !  This subroutine collects all the xml files contained in different
   !  directories and created by the diffent images in the phsave directory
   !  of the image 0
   !
   USE io_files,  ONLY : tmp_dir, prefix
   USE control_ph, ONLY : tmp_dir_ph
   USE save_ph,   ONLY : tmp_dir_save
   USE disp,      ONLY : nqs, lgamma_iq
   USE grid_irr_iq,  ONLY : comp_irr_iq, irr_iq
   USE control_ph,  ONLY : ldisp, epsil, zue, zeu
   USE klist,       ONLY : lgauss, ltetra
   USE el_phon,     ONLY : elph
   USE clib_wrappers,  ONLY : f_copy
   USE mp,        ONLY : mp_barrier
   USE mp_images, ONLY : my_image_id, nimage, intra_image_comm, nimage
   USE io_global, ONLY : stdout, ionode

   IMPLICIT NONE

   INTEGER :: iq, irr, ios, ima
   LOGICAL :: exst
   CHARACTER(LEN=256) :: file_input, file_output
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   CHARACTER(LEN=8) :: phpostfix

   CALL mp_barrier(intra_image_comm)
   IF (nimage==1) RETURN

#if defined (_WIN32)
   phpostfix='.phsave\'
#else
   phpostfix='.phsave/'
#endif

   DO iq=1,nqs
      DO irr=0, irr_iq(iq)
         IF (comp_irr_iq(irr,iq).and.ionode) THEN
            file_input=TRIM( tmp_dir_ph ) // &
                    & TRIM( prefix ) // phpostfix //'dynmat.'  &
                    &  // TRIM(int_to_char(iq))&
                    &  // '.' // TRIM(int_to_char(irr)) // '.xml'
!
!   copy all that has been calculated in all the images, so we can restart
!   without computing what is already available
!
            INQUIRE (FILE = TRIM(file_input), EXIST = exst)
            IF (exst) THEN
               DO ima=0, nimage-1
!                 ima=MOD(iq, nimage)
                  file_output=TRIM( tmp_dir_save ) // '/_ph'//&
                       &      TRIM(int_to_char(ima)) &
                       &    //'/'// TRIM( prefix ) // phpostfix //'dynmat.' &
                       &    // TRIM(int_to_char(iq))  &
                       &    // '.' // TRIM(int_to_char(irr)) // '.xml'
                  IF (ima /= my_image_id) ios = f_copy(file_input, file_output)
               ENDDO
            ENDIF
            IF ( elph .AND. irr>0 ) THEN

               file_input=TRIM( tmp_dir_ph ) // &
                    & TRIM( prefix ) // phpostfix //'elph.'  &
                    &  // TRIM(int_to_char(iq))&
                    &  // '.' // TRIM(int_to_char(irr)) // '.xml'

               file_output=TRIM( tmp_dir_save ) // '/_ph0/' // &
                    &   TRIM( prefix ) // phpostfix //'elph.' &
                    &    // TRIM(int_to_char(iq))  &
                    &    // '.' // TRIM(int_to_char(irr)) // '.xml'

               INQUIRE (FILE = TRIM(file_input), EXIST = exst)
               IF (exst) ios = f_copy(file_input, file_output)
            ENDIF
         ENDIF
      ENDDO
      IF ((ldisp.AND..NOT. (lgauss .OR. ltetra)).OR.(epsil.OR.zeu.OR.zue)) THEN
         IF (lgamma_iq(iq).AND.comp_irr_iq(0,iq).AND.ionode) THEN
            file_input=TRIM( tmp_dir_ph ) // &
                      TRIM( prefix ) // phpostfix //'tensors.xml'
            INQUIRE (FILE = TRIM(file_input), EXIST = exst)
            IF (exst) THEN
               DO ima=0, nimage-1
                  file_output=TRIM( tmp_dir_save ) // '/_ph'//&
                              TRIM(int_to_char(ima)) &
                      //'/'// TRIM( prefix ) // phpostfix //'tensors.xml'

                  IF (ima/=my_image_id) ios = f_copy(file_input, file_output)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   RETURN
   END SUBROUTINE collect_grid_files_tpw
