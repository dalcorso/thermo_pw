!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE set_2d_bz_path()
!--------------------------------------------------------------------
!
!  This routine sets a default path for each of the 2d Brillouin zones
!  Its sets also the number of points for each line of the path.
!
USE kinds,            ONLY : DP
USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                             trd_ht, rd_ht, cell_units
USE cell_base,        ONLY : at, celldm_cb => celldm, cell_base_init
USE control_paths,    ONLY : xqaux, wqaux, npk_label, letter,     &
                             label_list, nqaux, point_label_type, &
                             label_disp_q, letter_path, nrap_plot_in, &
                             rap_plot_in
USE bz_2d_form,       ONLY : find_ibrav_2d
IMPLICIT NONE

INTEGER :: ibrav_2d, i
REAL(DP) :: celldm_2d(3)
!
!  first select the 2d Brillouin zone type
!
CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                      trd_ht, rd_ht, cell_units )

CALL find_ibrav_2d(at, ibrav_2d, celldm_cb, celldm_2d)
!
!  Then set the number of path in the bz
!
npk_label=0
SELECT CASE (ibrav_2d) 

CASE (1) 
     npk_label=4
CASE (2)
     npk_label=6
CASE (3) 
     npk_label=6
CASE (4)
     npk_label=4
CASE (5) 
     npk_label=4
END SELECT
IF (npk_label==0) CALL errore('set_2d_bz_path','bz_type not supported',1)
nqaux=npk_label
!
!  Allocate the points and label necessary to specify this path
!
ALLOCATE(xqaux(3,nqaux))
ALLOCATE(wqaux(nqaux))
ALLOCATE(letter(nqaux))
ALLOCATE(letter_path(nqaux))
ALLOCATE(label_list(nqaux))
ALLOCATE(label_disp_q(nqaux))
ALLOCATE(nrap_plot_in(nqaux))
ALLOCATE(rap_plot_in(12,nqaux))
!
! and set the default path.
!
SELECT CASE (ibrav_2d) 

   CASE (1) 
!
!  Oblique bz
!
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG' /)
      wqaux(1:npk_label) = (/  30,   30,   30,   1   /)
      label_list(1:npk_label)=(/ ( i, i=1,npk_label) /)
      letter_path=letter
   CASE (2) 
!
!  rectangular bz
!
      letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG', 'S ' /)
      wqaux(1:npk_label) = (/ 30,   30,   30,    30,  30, 1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      letter_path=letter
   CASE (3) 
!
!   centered rectangular bz
!
      IF (celldm_2d(2) < 1.0_DP) THEN
         letter(1:npk_label)= (/ 'gG ', 'X  ', 'K  ', 'M  ', 'K1 ', 'gG ' /)  
         wqaux(1:npk_label)=  (/  30,     30,    30,    30,     30,   1   /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         letter(1:npk_label)= (/ 'gG ', 'K  ', 'M  ', 'K1 ', 'Y  ', 'gG ' /)  
         wqaux(1:npk_label)=  (/  30,     30,    30,     30,    30,   1   /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ENDIF
      letter_path=letter
   CASE (4)
!
! square bz
!
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG' /)  
      wqaux(1:npk_label) = (/  30,   30,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      letter_path=letter
CASE (5) 
!
!   hexagonal
!
      letter(1:npk_label)= (/ 'gG', 'K ', 'M ', 'gG' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      letter_path=letter
END SELECT
nrap_plot_in=0
rap_plot_in=0

RETURN
END SUBROUTINE set_2d_bz_path
