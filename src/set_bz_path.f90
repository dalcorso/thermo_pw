!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_bz_path()
USE kinds,            ONLY : DP
USE input_parameters, ONLY : ibrav, celldm
USE control_paths,    ONLY : xqaux, wqaux, npk_label, letter,     &
                             label_list, nqaux, point_label_type, &
                             label_disp_q, letter_path
USE bz_form,          ONLY : find_bz_type
USE thermo_mod,       ONLY : what
IMPLICIT NONE

INTEGER :: bzt, i
!
!  first select the Brillouin zone type
!
CALL find_bz_type(ibrav, celldm, bzt)
!
!  Then set the number of path in the bz
!
SELECT CASE (bzt) 

CASE (1, 2, 3) 
     nqaux=6
     IF (what=='mur_lc_t') nqaux=7
CASE (4) 
     nqaux=12
     IF (what=='mur_lc_t') nqaux=13
CASE (5) 
     nqaux=11
     IF (what=='mur_lc_t') nqaux=12
CASE (6) 
     nqaux=10
CASE (7) 
     nqaux=12
     IF (what=='mur_lc_t') nqaux=13
CASE (8) 
     nqaux=0
CASE (9) 
     nqaux=0
CASE (10) 
     nqaux=0
CASE (11) 
     nqaux=0
CASE (12) 
     nqaux=0
CASE (13) 
     nqaux=12
     IF (what=='mur_lc_t') nqaux=13
CASE (14) 
     nqaux=12
CASE (15) 
     nqaux=7
END SELECT
IF (nqaux==0) CALL errore('set_bz_type','bz_type not supported',1)
!
!  Allocate the points and label necessary to specify this path
!
ALLOCATE(xqaux(3,nqaux))
ALLOCATE(wqaux(nqaux))
ALLOCATE(letter(nqaux))
ALLOCATE(letter_path(nqaux))
ALLOCATE(label_list(nqaux))
ALLOCATE(label_disp_q(nqaux))
!
! and set the default path.
!
SELECT CASE (bzt) 

CASE (1) 
!
!  Simple cubic bz
!
   IF (what=='mur_lc_t') THEN
      npk_label=7
      letter(1:7)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'R ', 'X ' /)
      wqaux(1:7) = (/  30,   30,   45,    0,   50,   45,   1   /)
      label_list(1:7)=(/ ( i, i=1,7) /)
   ELSE
      npk_label=6
      letter(1:6)= (/ 'gG', 'X ', 'M ', 'gG', 'R ', 'X ' /)
      wqaux(1:6) = (/  30,   30,   45,   50,   45,   1   /)
      label_list(1:6)=(/ ( i, i=1,6) /)
   ENDIF
   letter_path=letter
CASE (2) 
!
!  fcc bz
!
   IF (what=='mur_lc_t') THEN
      npk_label=7
      letter(1:7)= (/ 'gG', 'X ', 'W ', 'M ', 'gG', 'gG', 'L ' /)
      wqaux(1:7) = (/ 40,   20,   20,    55,    0,   40,  1 /)
      label_list(1:7) =(/ (i, i=1, 7) /)
   ELSE
      npk_label=6
      letter(1:6)= (/ 'gG', 'X ', 'W ', 'M ', 'gG', 'L ' /)
      wqaux(1:6) = (/ 40,   20,   20,    55,    40,  1 /)
      label_list(1:6) =(/ (i, i=1, 6) /)
   ENDIF
   letter_path=letter
   letter_path(4)='X'
   point_label_type='BI'
CASE (3) 
!
!   bcc bz
!
   IF (what=='mur_lc_t') THEN
      npk_label=7
      letter(1:7)= (/ 'gG', 'H ', 'N ', 'gG', 'gG', 'P ', 'H ' /)  
      wqaux(1:7)=  (/  50,  40,    40,    0,    50,  50,   1   /)
      label_list(1:7) =(/ (i, i=1, 7) /)
   ELSE
      npk_label=6
      letter(1:6)= (/ 'gG', 'H ', 'N ', 'gG', 'P ', 'H ' /)  
      wqaux(1:6)=  (/  50,  40,    40,   50,   50,   1   /)
      label_list(1:6) =(/ (i, i=1, 6) /)
   ENDIF
   letter_path=letter
CASE (4)
!
! simple tetragonal lattice
!
   IF (what=='mur_lc_t') THEN
      npk_label=13
      letter(1:13)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'Z ', 'R ',  &
                       'A ', 'Z ', 'A ', 'M ', 'X ', 'R ' /)  
      wqaux(1:13) =  (/  30,   30,   45,  0,  40,   30,  30, &
                         45,   45,   40,    30,   40,  1 /)
      label_list(1:13) =(/ (i, i=1, 13) /)
   ELSE
      npk_label=12
      letter(1:12)= (/ 'gG', 'X ', 'M ', 'gG', 'Z ', 'R ',  &
                       'A ', 'Z ', 'A ', 'M ', 'X ', 'R ' /)  
      wqaux(1:12) =  (/  30,   30,   45,    40,   30,  30, &
                         45,   45,   40,    30,   40,  1 /)
      label_list(1:12) =(/ (i, i=1, 12) /)
   ENDIF
CASE (5) 
!
!   bct c < a
!
   IF (what=='mur_lc_t') THEN
      npk_label=12
      letter(1:12)= (/ 'gG', 'X ', 'L ', 'gG', 'gG', 'W ', 'H ',  &
                       'H1', 'L ', 'X ', 'P ', 'V ' /)  
      wqaux(1:12) =  (/  30,   30,   30,  0,  30,   30,  30, &
                         30,   30,   30,    30,   1 /)
      label_list(1:12) =(/ (i, i=1, 12) /)
   ELSE
      npk_label=11
      letter(1:11)= (/ 'gG', 'X ', 'L ', 'gG', 'W ', 'H ',  &
                       'H1', 'L ', 'X ', 'P ', 'V ' /)  
      wqaux(1:11) =  (/  30,   30,   30,    30,   30,  30, &
                         30,   30,   30,    30,   1 /)
      label_list(1:11) =(/ (i, i=1, 11) /)
   ENDIF
   letter_path=letter
CASE (6) 
!
!   bct c > a
!
   npk_label=10
   letter(1:10)= (/ 'gG ', 'gS ', 'Y  ', 'X  ', 'P  ', 'Y1 ',  &
                    'Z  ', 'gS1', 'N  ', 'gG ' /)  
   wqaux(1:10) =  (/  30,   30,   30,    30,   30,  30, &
                      30,   30,   30,     1 /)
   label_list(1:10) =(/ (i, i=1, 10) /)
   letter_path=letter
CASE (7) 
!
!  Simple orthorombic lattice
!
   IF (what=='mur_lc_t') THEN
      npk_label=13
      letter(1:13)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG', 'gG', 'Z ',  &
                       'U ', 'R ', 'T ', 'Z ', 'R ', 'S ' /)  
      wqaux(1:13) =  (/  30,   30,   30,    30,   0,  30,  30, &
                         30,   30,   30,    45,   30,  1 /)
      label_list(1:13) =(/ (i, i=1, 13) /)
   ELSE
      npk_label=12
      letter(1:12)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG', 'Z ',  &
                       'U ', 'R ', 'T ', 'Z ', 'R ', 'S ' /)  
      wqaux(1:12) =  (/  30,   30,   30,    30,   30,  30, &
                         30,   30,   30,    45,   30,  1 /)
      label_list(1:12) =(/ (i, i=1, 12) /)
   ENDIF
   letter_path=letter
CASE (8) 
CASE (9) 
CASE (10) 
CASE (11) 
CASE (12)
CASE (13) 
!
!  simple haxagonal bz
!
   IF (what=='mur_lc_t') THEN
      npk_label=13
      letter(1:13) = (/ 'gG', 'K ', 'M ', 'gG', 'gG', 'A ', 'H ', &
                        'L ', 'A ', 'H ', 'K ', 'M ', 'L ' /)
      wqaux(1:13) =  (/  30,   30,   30,   0,    30,   30,  30, &
                         30,   30,   40,   30,   30,  1 /)
      label_list(1:13) =(/ (i, i=1, 13) /)
   ELSE
      npk_label=12
      letter(1:12) = (/ 'gG', 'K ', 'M ', 'gG', 'A ', 'H ', &
                        'L ', 'A ', 'H ', 'K ', 'M ', 'L ' /)
      wqaux(1:12) =  (/  30,   30,   30,    30,   30,  30, &
                         30,   30,   40,    30,   30,  1 /)
      label_list(1:12) =(/ (i, i=1, 12) /)
   ENDIF
CASE (14) 
!
!  rombohedral (or trigonal) lattice  alpha < 90
!
   npk_label=12
   letter(1:12) = (/ 'gG', 'L ', 'X ', 'L1', 'P2', 'F ', &
                     'P1', 'Z ', 'B ', 'Q ', 'B1', 'gG' /)
   wqaux(1:12) =  (/  30,   30,   30,    30,   30,  30, &
                      30,   30,   40,    30,   30,  1 /)
   label_list(1:12) =(/ (i, i=1, 12) /)
   letter_path=letter
CASE (15) 
!
!  rombohedral (or trigonal) lattice  alpha > 90
!
   npk_label=7
   letter(1:7) = (/ 'gG', 'P ', 'P1', 'Q1', 'Q ', 'Z ', 'gG' /)
   wqaux(1:7) =  (/  30,   30,   30,    30,   30,  30, 1 /)
   label_list(1:7) =(/ (i, i=1, 7) /)
   letter_path=letter
END SELECT

RETURN
END SUBROUTINE set_bz_path
