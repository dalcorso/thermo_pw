!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE read_lanczos_chain()
!-----------------------------------------------------------------------

USE kinds,      ONLY : DP
USE cell_base,  ONLY : omega
USE lr_lanczos, ONLY : iulanczos, lanczos_steps, beta_store, gamma_store, &
                       zeta_store
USE lr_global,  ONLY : rpert
USE io_global,  ONLY : ionode, ionode_id
USE mp,         ONLY : mp_bcast
USE mp_images,  ONLY : intra_image_comm

IMPLICIT NONE
INTEGER :: ipert, jpert, iter, idum, ios

IF (ionode) THEN
   REWIND(iulanczos)
   DO iter = 1, lanczos_steps
      READ(iulanczos, *, ERR=100, END=100, IOSTAT=ios) idum, gamma_store(iter)
      DO ipert = 1, rpert
         DO jpert = 1, rpert
            IF (MOD(iter,2)==0) THEN
               READ(iulanczos,'(2e25.15)',ERR=100,END=100,IOSTAT=ios) &
                                        zeta_store(ipert, jpert, iter)

            ELSE
               zeta_store(ipert, jpert, iter)=(0.0_DP,0.0_DP)
            ENDIF
         ENDDO
      ENDDO
      beta_store(iter)=ABS(gamma_store(iter))
   ENDDO
ENDIF
100 CALL mp_bcast(ios, ionode_id, intra_image_comm)
CALL errore('read_lanczos_chain', 'The save file has not enough iteration', &
                                                         ABS(ios))

CALL mp_bcast(beta_store, ionode_id, intra_image_comm)
CALL mp_bcast(gamma_store, ionode_id, intra_image_comm)
CALL mp_bcast(zeta_store, ionode_id, intra_image_comm)

RETURN
END SUBROUTINE read_lanczos_chain
