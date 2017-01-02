!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_xrdp(ext)

USE kinds,            ONLY : DP

USE cell_base,        ONLY : at, bg
USE ions_base,        ONLY : nat, atm, tau, nsp, ityp
USE control_pressure, ONLY : pressure_kb
USE initial_conf,     ONLY : ibrav_save
USE control_xrdp,     ONLY : lxrdp, lambda, flxrdp, lcm
USE equilibrium_conf, ONLY : celldm0

USE xrdp_module,      ONLY : compute_xrdp

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: ext

CHARACTER(LEN=256) :: filename
CHARACTER(LEN=8)   :: float_to_char

IF (ext /= ' ') THEN
   filename=TRIM(flxrdp)//TRIM(ext)
ELSE
   filename=TRIM(flxrdp)
ENDIF
IF (pressure_kb /= 0.0_DP) &
     filename=TRIM(flxrdp)//'.'//TRIM(float_to_char(pressure_kb,1))

CALL compute_xrdp(at,bg,celldm0(1),nat,tau,nsp,ityp,atm, &
                                lambda,ibrav_save,lcm,filename)
CALL plot_xrdp(ext)

RETURN
END SUBROUTINE manage_xrdp

