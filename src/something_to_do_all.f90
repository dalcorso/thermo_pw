!
! Copyright (C) 2013-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION something_to_do_all(iwork, igeom, iqw, irrw)

USE initial_conf, ONLY : collect_info_save
USE control_qe,   ONLY : use_ph_images
USE mp_images,    ONLY : nimage

IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork, igeom, iqw, irrw

LOGICAL :: std
INTEGER :: iq, irr, image

std=.FALSE.
IF (use_ph_images) THEN
   image=MOD(iwork-1, nimage) + 1 
   DO iq=1, collect_info_save(igeom)%nqs
      DO irr=0, collect_info_save(igeom)%irr_iq(iq)
         IF ((collect_info_save(igeom)%comp_irr_iq(irr,iq,image)==1).AND.&
             (collect_info_save(igeom)%done_irr_iq(irr,iq,image)==0)) std=.TRUE.
      ENDDO
   ENDDO
ELSE
   IF (collect_info_save(igeom)%done_irr_iq(irrw,iqw,1)==0) std=.TRUE.
ENDIF
something_to_do_all=std

RETURN
END FUNCTION something_to_do_all
