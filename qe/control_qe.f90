!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables that can be used by thermo_pw to control
!     internal working of QE routines. These variables are used
!     only in the routines of this directory and in those of thermo_pw
!     that sets them
!
MODULE control_qe
  USE kinds,  ONLY : DP
  !
  !
  SAVE

  LOGICAL  :: tcollect_all=.FALSE.

END MODULE control_qe

