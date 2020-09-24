!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
SUBROUTINE set_bz_path()
!-------------------------------------------------------------------
!
!   This routine sets a default path for each three dimensional Brillouin
!   zone.
!
USE kinds,            ONLY : DP
USE cell_base,        ONLY : ibrav, celldm
USE control_paths,    ONLY : xqaux, wqaux, npk_label, letter,     &
                             label_list, nqaux, point_label_type, &
                             label_disp_q, letter_path, long_path, old_path, &
                             path_fact
USE bz_form,          ONLY : find_bz_type
USE thermo_sym,       ONLY : sg_number
USE thermo_mod,       ONLY : what
IMPLICIT NONE

INTEGER :: bzt, i
LOGICAL :: double_gamma
!
!  skip unsupported ibrav
!
nqaux=0
IF (ibrav==91 .OR. ibrav==-13 .OR. ibrav==13 .OR. ibrav==14 .OR. ibrav==0) &
                                               RETURN
!
!  cases in which we need to double gamma
!
double_gamma=(what=='mur_lc_t'.OR.what=='elastic_constants_t'.OR.&
              what=='scf_disp'.OR.what=='mur_lc_disp')
!
!  first select the Brillouin zone type
!
CALL find_bz_type(ibrav, celldm, bzt)
!
!  Then set the number of path in the bz
!
SELECT CASE (bzt) 

CASE (1) 
     IF (long_path) THEN
        npk_label=8
        IF (double_gamma) npk_label=9
     ELSE
        npk_label=4
     ENDIF
CASE (2)
     IF (long_path) THEN
        npk_label=11
        IF (double_gamma) npk_label=12
     ELSE
        npk_label=6
        IF (double_gamma) npk_label=7
     ENDIF
CASE(3)
     npk_label=8
     IF (double_gamma) npk_label=9
CASE(4)
     IF (long_path) THEN
        npk_label=12
        IF (double_gamma) npk_label=13
     ELSE
        npk_label=4
     ENDIF
CASE(5) 
    IF (long_path) THEN
       npk_label=11
       IF (double_gamma) npk_label=12
    ELSE
       npk_label=5
       IF (double_gamma) npk_label=6
    ENDIF
CASE(7)
     IF (long_path) THEN
        npk_label=16
        IF (double_gamma) npk_label=17
     ELSE
        npk_label=5
     ENDIF
CASE (6) 
     npk_label=13
     IF (double_gamma) npk_label=15
CASE (8) 
     npk_label=12
     IF (double_gamma) npk_label=13
CASE (9) 
     npk_label=14
     IF (double_gamma) npk_label=15
CASE (10)
     npk_label=10
     IF (double_gamma) npk_label=11
CASE (11) 
     npk_label=15
     IF (double_gamma) npk_label=16
CASE (12) 
     IF (long_path) THEN
        npk_label=12
        IF (double_gamma) npk_label=13
     ELSE
        npk_label=6
     ENDIF
CASE (13) 
     IF (long_path.OR.old_path) THEN
        npk_label=12
        IF (double_gamma) npk_label=13
     ELSE
        npk_label=4
     ENDIF
CASE (14) 
     npk_label=12
CASE (15) 
     npk_label=7
CASE (16) 
   IF (long_path.OR.ibrav==-12) THEN
      npk_label=7
   ELSE
      npk_label=3
   ENDIF
END SELECT
IF (npk_label==0) CALL errore('set_bz_path','bz_type not supported',1)
nqaux=npk_label
!
!  Allocate the points and label necessary to specify this path
!
IF (ALLOCATED(xqaux)) DEALLOCATE(xqaux)
IF (ALLOCATED(wqaux)) DEALLOCATE(wqaux)
IF (ALLOCATED(letter)) DEALLOCATE(letter)
IF (ALLOCATED(letter_path)) DEALLOCATE(letter_path)
IF (ALLOCATED(label_list)) DEALLOCATE(label_list)
IF (ALLOCATED(label_disp_q)) DEALLOCATE(label_disp_q)

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
   IF (long_path) THEN
      IF (double_gamma) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'R ', 'X ', &
                                 'R ', 'M ' /)
         wqaux(1:npk_label) = (/  30,   30,   45,    0,   50,   45,     0, & 
                                  30,    1  /)
         label_list(1:npk_label)=(/ ( i, i=1,npk_label) /)
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'R ', 'X ', &
                                 'R ', 'M '  /)
         wqaux(1:npk_label) = (/  30,   30,   45,   50,   45,    0,  &
                                  30,    1   /)
         label_list(1:npk_label)=(/ ( i, i=1,npk_label) /)
      ENDIF
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG' /)
      wqaux(1:npk_label) = (/  30,   30,   45,   1   /)
      label_list(1:npk_label)=(/ ( i, i=1,npk_label) /)
   ENDIF
   letter_path=letter
   point_label_type='SC'
CASE (2) 
!
!  fcc bz
!
   IF (long_path) THEN
      IF (double_gamma) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'K ', 'gG', &
                                 'gG', 'L ', 'W ', 'X ', 'L ', 'K ', 'W ' /)
         wqaux(1:npk_label) = (/ 40,  0,  15,  40,  0,   40,  40, &
                                 20,  0, 40, 30, 1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'K ', 'gG', 'L ',&
                                 'W ', 'X ', 'L ', 'K ', 'W ' /)
         wqaux(1:npk_label) = (/ 40,  0,  15,  30,  40,  40, 20, 0, 40, 30, 1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ENDIF
      letter_path=letter
      letter_path(3)='X'
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'K ', 'gG', 'L '/)
      wqaux(1:npk_label) = (/ 30,      0,   15,   30,   30,  1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      letter_path=letter
      letter_path(3)='X'
   ENDIF
   point_label_type='BI'
CASE (3) 
!
!   bcc bz
!
   IF (double_gamma) THEN
      letter(1:npk_label)= (/ 'gG', 'H ', 'N ', 'gG', 'gG', 'P ', 'H ', &
                              'P ', 'N ' /)  
      wqaux(1:npk_label)=  (/  40,  30,    30,    0,    30,   30,    0, &
                                30,   1   /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG', 'H ', 'N ', 'gG', 'P ', 'H ', &
                              'P ', 'N ' /)  
      wqaux(1:npk_label)=  (/  40,  30,    30,   40,   40,     0, &
                               30,   1   /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
   point_label_type='SC'
CASE (4)
!
! simple tetragonal lattice
!
   IF (long_path) THEN
      IF (double_gamma) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'Z ', 'R ',  &
                                 'A ', 'Z ', 'X ', 'R ', 'M ', 'A '  /)  
         wqaux(1:npk_label) =  (/  30,   30,   45,  0,  40,   30,  30, &
                                45,    0,   30,  0,  30,  1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'Z ', 'R ',  &
                          'A ', 'Z ', 'X ', 'R ', 'M ', 'A ' /)  
         wqaux(1:npk_label) =  (/  30,   30,   45,    40,   30,  30, &
                         45,   0,  30,  0,  30,  1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ENDIF
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG' /)  
      wqaux(1:npk_label) = (/  30,   30,   45,    1  /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
   point_label_type='SC'
CASE (5) 
!
!   bct c < a
!
   IF (long_path) THEN
      IF (double_gamma) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'Z ', 'Z1', &
                                 'M ', 'X ', 'P ', 'N ', 'gG' /)  
         wqaux(1:npk_label) =  (/  30,   30,   30,  0,  30,   0,   &
                                   30,    0,   30,  30,  30,  1   /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
         letter_path=letter
         letter_path(7)='Z '
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'Z ',   &
                                 'Z1', 'M ', 'X ', 'P ', 'N ', 'gG' /)  
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,    0,    &
                                30,    0,   30,   30,   30,    1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
         letter_path=letter
         letter_path(6)='Z '
      ENDIF
   ELSE
      IF (double_gamma) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'gG', 'Z ' /)
         wqaux(1:npk_label) =  (/  30,   30,   30,  0,  30,   1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'M ', 'gG', 'Z ' /)
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ENDIF
      letter_path=letter
   ENDIF
   point_label_type='SC'
CASE (6) 
!
!   bct c > a
!
   IF (double_gamma) THEN
      letter(1:npk_label)= (/ 'gG ', 'X  ', 'P  ', 'N  ', 'gG ', 'gG ', &
                              'M  ', 'S  ', 'S0 ', 'gG ', 'X  ', 'R  ', &
                              'gG ', 'gG ', 'M  ' /)  
      wqaux(1:npk_label) =  (/   30,    30,    30,    30,     0,    30,  &  
                                 30,     0,    30,     0,    30,     0,  &  
                                 0,     30,     1  /)

      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ELSE
      letter(1:npk_label)= (/ 'gG ', 'X  ', 'P  ', 'N  ', 'gG ', 'M  ', 'S  ', &
                              'S0 ', 'gG ', 'X  ', 'R  ', 'gG ', 'M  ' /)  
      wqaux(1:npk_label) =  (/   30,    30,    30,    30,    30,    30,    0, &
                                 30,     0,    30,     0,    30,      1  /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
   point_label_type='BI'

CASE (7) 
!
!  Simple orthorhombic lattice
!
   IF (long_path) THEN
      IF (double_gamma) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG', 'gG', 'Z ', &
                                 'U ', 'R ', 'T ', 'Z ', 'X ', 'U ', 'Y ', &
                                 'T ', 'S ', 'R '  /)  
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,    0,   30,   30, &
                                   30,   30,   30,    0,   30,    0,   30, & 
                                    0,   30,    1   /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG', 'Z ',       &
                                 'U ', 'R ', 'T ', 'Z ', 'X ', 'U ', 'Y ', & 
                                 'T ', 'S ', 'R '  /)  
         wqaux(1:npk_label) =  (/  30,   30,   30,  30,   30,  30,        &
                                   30,   30,   30,   0,   30,    0,   30, &  
                                    0,   30,    1   /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ENDIF
   ELSE
      letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'Y ', 'gG' /)  
      wqaux(1:npk_label) = (/  30,    30,   30,   30,   1  /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF

   letter_path=letter
   point_label_type='SC'
CASE (8) 
!
!  Face centered Orthorombic case 1
!
   IF (double_gamma) THEN
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
   point_label_type='SC'
CASE (9) 
!
!  Face centered orthorhombic case 2
!
   IF (double_gamma) THEN
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
   point_label_type='SC'
CASE (10) 
!
!  Face centered orthorhombic case 3
!
   IF (double_gamma) THEN
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
   point_label_type='SC'
CASE (11) 
!
!  Body centered orthorhombic
!
   IF (double_gamma) THEN
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
   point_label_type='SC'
CASE (12)
!
!  Base centered orthorhombic
!
   IF (long_path) THEN
      IF (double_gamma) THEN
         IF (celldm(2) >= 1.0_DP) THEN
            letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'R ', 'A ', 'Z ', 'gG', &
                                    'gG', 'Y ', 'X1', 'A1', 'T ', 'Y ' /)  
         ELSE
            letter(1:npk_label)= (/ 'gG', 'Y ', 'S ', 'R ', 'T ', 'Z ', 'gG',  &
                                    'gG', 'X ', 'Y1', 'A1', 'A ', 'X ' /)  
         ENDIF
         wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,   30,  30,  &
                                    0,   30,   30,   30,   30,    1  /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         IF (celldm(2) >= 1.0_DP) THEN
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
   ELSE
      IF (celldm(2) >= 1.0_DP) THEN
         letter(1:npk_label)= (/ 'gG', 'X ', 'S ', 'X1', 'Y ', 'gG' /)  
      ELSE
         letter(1:npk_label)= (/ 'gG', 'X ', 'Y1', 'S ', 'Y ', 'gG' /) 
      ENDIF
      wqaux(1:npk_label) =  (/  30,   30,   30,   30,   30,    1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   ENDIF
   letter_path=letter
   point_label_type='SC'
CASE (13) 
!
!  simple hexagonal bz
!
   IF (old_path) THEN
      IF (double_gamma) THEN
         letter(1:npk_label) = (/ 'gG', 'K ', 'M ', 'gG', 'gG', 'A ', 'H ', &
                                  'L ', 'A ', 'M ', 'L ', 'K ', 'H ' /)
         wqaux(1:npk_label) =  (/  30,   30,   30,   0,    30,   30,  30, &
                                   30,   30,   30,   30,   30,    1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         letter(1:npk_label) = (/ 'gG', 'K ', 'M ', 'gG', 'A ', 'H ', &
                           'L ', 'A ', 'M ', 'L ', 'K ', 'H ' /)
         wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, &
                                   30,   30,   30,    30,   30,  1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ENDIF
   ELSEIF (long_path) THEN
      IF (double_gamma) THEN
         letter(1:npk_label) = (/ 'gG', 'M ', 'K ', 'gG', 'gG', 'A ', 'L ', &
                                  'H ', 'A ', 'L ', 'M ', 'H ', 'K ' /)
         wqaux(1:npk_label) =  (/  30,   30,   30,   0,    30,   30,  30, &
                                   30,    0,   30,    0,   30,    1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ELSE
         letter(1:npk_label) = (/ 'gG', 'M ', 'K ', 'gG', 'A ', 'L ', &
                           'H ', 'A ', 'L ', 'M ', 'H ', 'K ' /)
         wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, &
                                   30,    0,   30,     0,   30,  1 /)
         label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
      ENDIF
   ELSE
      letter(1:npk_label) = (/ 'gG', 'M ', 'K ', 'gG' /)
      wqaux(1:npk_label) =  (/  30,   30,   30,   1 /)
      label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   END IF
   letter_path=letter
   point_label_type='SC'
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
   point_label_type='SC'
CASE (15) 
!
!  rombohedral (or trigonal) lattice  alpha > 90
!
   letter(1:npk_label) = (/ 'gG', 'P ', 'P1', 'Q1', 'Q ', 'Z ', 'gG' /)
   wqaux(1:npk_label) =  (/  30,   30,   30,    30,   30,  30, 1 /)
   label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   letter_path=letter
   point_label_type='SC'
CASE (16) 
!
!  simple monoclinic lattice
!
   IF (long_path.OR.ibrav==-12) THEN
      IF (ibrav==12) THEN
         letter(1:npk_label) = (/ 'gG', 'X ', 'A ', 'Z ', 'D ', 'Y ', 'gG' /)
         wqaux(1:npk_label) =  (/  30,    30,   30,   30,   30,  30,    1 /)
      ELSE
         letter(1:npk_label) = (/ 'gG', 'X ', 'A ', 'Y ', 'D ', 'Z ', 'gG' /)
         wqaux(1:npk_label) =  (/  30,    30,   30,   30,   30,  30,    1 /)
      ENDIF
   ELSE
      IF (ibrav==12) THEN
         letter(1:npk_label) = (/ 'Y ', 'gG', 'X '  /)
         wqaux(1:npk_label) =  (/  30,    30,   1   /)
      ENDIF
   ENDIF
   label_list(1:npk_label) =(/ (i, i=1, npk_label) /)
   letter_path=letter
   point_label_type='SC'
END SELECT
wqaux(1:npk_label)=wqaux(1:npk_label)*path_fact

RETURN
END SUBROUTINE set_bz_path
