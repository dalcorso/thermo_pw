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

CASE (1, 3) 
     npk_label=6
     IF (what=='mur_lc_t') npk_label=7
CASE (2)
     npk_label=7
     IF (what=='mur_lc_t') npk_label=8
CASE (4, 7, 12, 13) 
     npk_label=12
     IF (what=='mur_lc_t') npk_label=13
CASE (5,10) 
     npk_label=10
     IF (what=='mur_lc_t') npk_label=11
CASE (6) 
     npk_label=10
CASE (8) 
     npk_label=12
     IF (what=='mur_lc_t') npk_label=13
CASE (9) 
     npk_label=14
     IF (what=='mur_lc_t') npk_label=15
CASE (11) 
     npk_label=15
     IF (what=='mur_lc_t') npk_label=16
CASE (14) 
     npk_label=12
CASE (15) 
     npk_label=7
END SELECT
IF (npk_label==0) CALL errore('set_bz_type','bz_type not supported',1)
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
!
! and set the default path.
!
SELECT CASE (bzt) 

CASE (1) 
!
!  Simple cubic bz
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'R ', 'X ' /)
      wqaux(1:npk_label) = (/  30,   30,   45,    0,   50,   45,   1   /)
      label_list(1:npk_label)=(/ ( i, i=1,npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'R ', 'X ' /)
      wqaux(1:npk_label) = (/  30,   30,   45,   50,   45,   1   /)
      label_list(1:npk_label)=(/ ( i, i=1,npk_label) /)
   ENDIF
   letter_path=letter
CASE (2) 
!
!  fcc bz
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'X ', 'W ', 'M ', 'K ', 'gG', 'gG', 'L ' /)
      wqaux(1:npk_label) = (/ 40,   20,   20,    15,  40,  0,   40,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'W ', 'M ', 'K ', 'gG', 'L ' /)
      wqaux(1:npk_label) = (/ 40,   20,   20,    15,  40,  40,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
   letter_path(4)='X'
   point_label_type='BI'
CASE (3) 
!
!   bcc bz
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'H ', 'N ', 'gG', 'gG', 'P ', 'H ' /)  
      wqaux(1:npk_label)=  (/  50,  40,    40,    0,    50,  50,   1   /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'H ', 'N ', 'gG', 'P ', 'H ' /)  
      wqaux(1:npk_label)=  (/  50,  40,    40,   50,   50,   1   /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (4)
!
! simple tetragonal lattice
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'Z ', 'R ',  &
                              'A ', 'Z ', 'A ', 'M ', 'X ', 'R ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   45,  0,  40,   30,  30, &
                                45,   45,   40,    30,   40,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'Z ', 'R ',  &
                       'A ', 'Z ', 'A ', 'M ', 'X ', 'R ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   45,    40,   30,  30, &
                         45,   45,   40,    30,   40,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (5) 
!
!   bct c < a
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'Z ', 'P ',  &
                              'N ', 'Z1', 'M ', 'P ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,  0,  30,   30,  30, &
                                30,   30,   30,  1   /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'Z ', 'P ',  &
                              'N ', 'Z1', 'M ', 'P ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  &
                                30,   30,   30,    1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (6) 
!
!   bct c > a
!
   letter(1:npk_label)= (/ 'gG ', 'gS ', 'Y  ', 'X  ', 'P  ', 'Y1 ',  &
                    'Z  ', 'gS1', 'N  ', 'gG ' /)  
   wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, &
                      30,   30,   30,     1 /)
   label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   letter_path=letter
CASE (7) 
!
!  Simple orthorombic lattice
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG', 'gG', 'Z ',  &
                       'U ', 'R ', 'T ', 'Z ', 'R ', 'S ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,    30,   0,  30,  30, &
                         30,   30,   30,    45,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG', 'Z ',  &
                       'U ', 'R ', 'T ', 'Z ', 'R ', 'S ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, &
                         30,   30,   30,    45,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (8) 
!
!  Face centered Orthorombic case 1
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'X ', 'A1', 'Y ', 'T ', 'Z ', 'gG',  &
                       'gG', 'Y ', 'X1', 'L ', 'A ', 'Z ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,    0, &
                         30,   30,    30,   30,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'A1', 'Y ', 'T ', 'Z ',  &
                       'gG', 'Y ', 'X1', 'L ', 'A ', 'Z ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, &
                         30,   30,   30, 30,  30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (9) 
!
!  Face centered orthorombic case 2
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'Y ', 'C ', 'D ', 'X ', 'gG', 'gG',  &
                       'Z ', 'D1', 'C1', 'H1', 'D ', 'D1', 'Z ', 'gG' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,    0,   30,  &
                         30,   30,   30,   45,   30,   30,   30,   1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'Y ', 'C ', 'D ', 'X ', 'gG',        &
                       'Z ', 'D1', 'C1', 'H1', 'D ', 'D1', 'Z ', 'gG' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,        &
                         30,   30,   30,   45,   30,   30,   30,   1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (10) 
!
!  Face centered orthorombic case 3
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label)= (/ 'gG', 'Y ', 'T ', 'Z ', 'gG', 'gG',  &
                       'X ', 'A1', 'L ', 'A ', 'Z ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   0,   30,  &
                         30,   30,   30,   30,     1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'Y ', 'T ', 'Z ', 'gG',        &
                       'X ', 'A1', 'L ', 'A ', 'Z ' /)  
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,        &
                         30,   30,   30,   30,    1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (11) 
!
!  Body centered orthorombic
!
   IF (what=='mur_lc_t') THEN
      IF (celldm(2) < 1.0_DP .AND. celldm(3) < 1.0_DP ) THEN
         letter(1:npk_label)= (/ 'gG', 'Y ', 'T ', 'W ', 'S ', 'L2', 'Z ', &
                           'gG', 'gG', 'X ', 'Y1', 'L ', 'Z1', 'R ', 'W ', 'L '/) 
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,  &
                              0,   30,   30,   30,   30,   30,  30,  30, 1  /)
      ELSEIF (celldm(2) > 1.0_DP .AND. celldm(2) > celldm(3) ) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'T ', 'W ', 'R ', 'L2', 'Z ', &
                           'gG', 'gG', 'Y ', 'X1', 'L1', 'Z1', 'S ', 'W ', 'L '/) 
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,  &
                              0,   30,   30,   30,   30,   30,  30,  30, 1  /)
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'L ', 'T ', 'W ', 'R ', 'X1', &
                           'Z ', 'gG', 'gG', 'Y ', 'L1', 'W ', 'S ', 'Y1', 'Z '/) 
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,  &
                         30,     0,   30,   30,   30,   30,   30,  30,  1  /)
      ENDIF
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      IF (celldm(2) < 1.0_DP .AND. celldm(3) < 1.0_DP ) THEN
         letter(1:npk_label)= (/ 'gG', 'Y ', 'T ', 'W ', 'S ', 'L2', 'Z ',   &
                           'gG', 'X ', 'Y1', 'L ', 'Z1', 'R ', 'W ', 'L ' /)  
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,    &
                                   30,   30,   30,   30,   30,  30,  30,  1 /)
      ELSEIF (celldm(2) > 1.0_DP .AND. celldm(2) > celldm(3) ) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'T ', 'W ', 'R ', 'L2', 'Z ',   &
                           'gG', 'Y ', 'X1', 'L1', 'Z1', 'S ', 'W ', 'L ' /)  
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,    &
                                   30,   30,   30,   30,   30,  30,  30,  1 /)
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'L ', 'T ', 'W ', 'R ', 'X1',   &
                           'Z ', 'gG', 'Y ', 'L1', 'W ', 'S ', 'Y1', 'Z ' /)  
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,    &
                                   30,   30,   30,   30,   30,  30,  30,  1 /)
      ENDIF
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (12)
!
!  One phase centered orthorombic
!
   IF (what=='mur_lc_t') THEN
      IF (celldm(2) > 1.0_DP) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'R ', 'A ', 'Z ', 'gG', &
                                 'gG', 'Y ', 'X1', 'A1', 'T ', 'Y ' /)  
      ELSE
         letter(1:npk_label)= (/ 'gG', 'Y ', 'S ', 'R ', 'T ', 'Z ', 'gG',   &
                                 'gG', 'X ', 'Y1', 'A1', 'A ', 'X ' /)  
      ENDIF
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,  &
                                 0,   30,   30,   30,   30,    1  /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      IF (celldm(2) > 1.0_DP) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'R ', 'A ', 'Z ',    &
                                 'gG', 'Y ', 'X1', 'A1', 'T ', 'Y ' /)  
      ELSE
         letter(1:npk_label)= (/ 'gG', 'Y ', 'S ', 'R ', 'T ', 'Z ',    &
                                 'gG', 'X ', 'Y1', 'A1', 'A ', 'X ' /)  
      ENDIF
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,    &
                                30,   30,   30,   30,   30,    1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (13) 
!
!  simple hexagonal bz
!
   IF (what=='mur_lc_t') THEN
      letter(1:npk_label) = (/ 'gG', 'K ', 'M ', 'gG', 'gG', 'A ', 'H ', &
                        'L ', 'A ', 'H ', 'K ', 'M ', 'L ' /)
      wqaux(1:npk_label) =  (/  30,   30,   30,   0,    30,   30,  30, &
                         30,   30,   40,   30,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label) = (/ 'gG', 'K ', 'M ', 'gG', 'A ', 'H ', &
                        'L ', 'A ', 'H ', 'K ', 'M ', 'L ' /)
      wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, &
                         30,   30,   40,    30,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
CASE (14) 
!
!  rombohedral (or trigonal) lattice  alpha < 90
!
   letter(1:npk_label) = (/ 'gG', 'L ', 'X ', 'L1', 'P2', 'F ', &
                            'P1', 'Z ', 'B ', 'Q ', 'B1', 'gG' /)
   wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, &
                      30,   30,   40,    30,   30,  1 /)
   label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   letter_path=letter
CASE (15) 
!
!  rombohedral (or trigonal) lattice  alpha > 90
!
   letter(1:npk_label) = (/ 'gG', 'P ', 'P1', 'Q1', 'Q ', 'Z ', 'gG' /)
   wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, 1 /)
   label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   letter_path=letter
END SELECT

RETURN
END SUBROUTINE set_bz_path
