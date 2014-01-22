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
IMPLICIT NONE

INTEGER :: bzt
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
CASE (4) 
     nqaux=12
CASE (5) 
     nqaux=0
CASE (6) 
     nqaux=0
CASE (7) 
     nqaux=0
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
     nqaux=0
CASE (14) 
     nqaux=0
CASE (15) 
     nqaux=0
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
   npk_label=6
   letter(1)='gG'
   label_list(1)=1
   wqaux(1)=30
   letter(2)='X'
   label_list(2)=2
   wqaux(2)=30
   letter(3)='M'
   label_list(3)=3
   wqaux(3)=45
   letter(4)='gG'
   label_list(4)=4
   wqaux(4)=50
   letter(5)='R'
   label_list(5)=5
   wqaux(5)=45
   letter(6)='X'
   label_list(6)=6
   wqaux(6)=1
   letter_path=letter
CASE (2) 
!
!  fcc bz
!
   npk_label=6
   letter(1)='gG'
   label_list(1)=1
   wqaux(1)=40
   letter(2)='X'
   label_list(2)=2
   wqaux(2)=20
   letter(3)='W'
   label_list(3)=3
   wqaux(3)=20
   letter(4)='M'
   label_list(4)=4
   wqaux(4)=55
   letter(5)='gG'
   label_list(5)=5
   wqaux(5)=40
   letter(6)='L'
   label_list(6)=6
   wqaux(6)=1
   letter_path=letter
   letter_path(4)='X'
   point_label_type='BI'
CASE (3) 
   npk_label=6
   letter(1)='gG'
   label_list(1)=1
   wqaux(1)=50
   letter(2)='H'
   label_list(2)=2
   wqaux(2)=40
   letter(3)='H'
   label_list(3)=3
   wqaux(3)=40
   letter(4)='gG'
   label_list(4)=4
   wqaux(4)=50
   letter(5)='P'
   label_list(5)=5
   wqaux(5)=50
   letter(6)='H'
   label_list(6)=6
   wqaux(6)=1
   letter_path=letter
CASE (4) 
   npk_label=12
   letter(1)='gG'
   label_list(1)=1
   wqaux(1)=30
   letter(2)='K'
   label_list(2)=2
   wqaux(2)=30
   letter(3)='M'
   label_list(3)=3
   wqaux(3)=30
   letter(4)='gG'
   label_list(4)=4
   wqaux(4)=30
   letter(5)='A'
   label_list(5)=5
   wqaux(5)=30
   letter(6)='H'
   label_list(6)=6
   wqaux(6)=30
   letter(7)='L'
   label_list(7)=7
   wqaux(7)=30
   letter(8)='A'
   label_list(8)=8
   wqaux(8)=30
   letter(9)='H'
   label_list(9)=9
   wqaux(9)=40
   letter(10)='K'
   label_list(10)=10
   wqaux(10)=40
   letter(11)='M'
   label_list(11)=11
   wqaux(11)=30
   letter(12)='L'
   label_list(12)=12
   wqaux(12)=30
   letter_path=letter
CASE (5) 
CASE (6) 
CASE (7) 
CASE (8) 
CASE (9) 
CASE (10) 
CASE (11) 

END SELECT

RETURN
END SUBROUTINE set_bz_path
