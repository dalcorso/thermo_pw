!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE save_ph_variables()
!--------------------------------------------------------------------------
!
USE control_ph,   ONLY : epsil, zeu, zue, start_q, last_q
USE initial_conf, ONLY : epsil_save, zeu_save, zue_save, start_q_save, &
                         last_q_save
IMPLICIT NONE

epsil_save=epsil
zeu_save=zeu
zue_save=zue
start_q_save=start_q
last_q_save=last_q

RETURN
END SUBROUTINE save_ph_variables
