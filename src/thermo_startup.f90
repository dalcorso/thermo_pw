!
! Copyright (C) 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------
SUBROUTINE thermo_startup(code)
!-------------------------------------------------
!
! Initialize MPI, clocks, pools, bands, images, etc. print initial messages
!
USE mp_global,        ONLY : mp_startup
USE environment,      ONLY : environment_start
USE mp_world,         ONLY : world_comm
USE mp_pools,         ONLY : intra_pool_comm
USE mp_bands,         ONLY : intra_bgrp_comm, inter_bgrp_comm
USE command_line_options,  ONLY : ndiag_

IMPLICIT NONE
INCLUDE 'laxlib.fh'
CHARACTER(LEN=*) :: code
  !
  CALL mp_startup( start_images=.TRUE. )
  !
  CALL laxlib_start ( ndiag_, world_comm, intra_bgrp_comm, &
                         do_distr_diag_inside_bgrp_ = .TRUE. )
  !
  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
       inter_bgrp_comm )
  !
  CALL environment_start ( code )
  !
  CALL start_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE thermo_startup
!
!-------------------------------------------------
SUBROUTINE thermo_end(code)
!-------------------------------------------------
!
!  This subroutine deallocates all the variables of thermo_pw, close
!  the files and closes MPI
!
USE mp_global,        ONLY : mp_global_end
USE environment,      ONLY : environment_end

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: code
   !
   CALL deallocate_thermo()
   !
   CALL laxlib_end()
   !
   CALL environment_end( code )
   !
   CALL mp_global_end ()
   !
   RETURN
END SUBROUTINE thermo_end
