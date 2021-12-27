!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE symmetrize_drho(drhoscf, dbecsum, irr, npe, code)
!
!  This subroutine symmetrizes drhoscf and, in the PAW case dbecsum, 
!  calling the appropriate routines according to code
!  
!  code
!   1   : phonon perturbation
!   2   : electric field perturbation q=0
!   3   : electric field perturbation q/=0 (for eels) 
!   4   : electric field perturbation q=0 (complex quantity for optical)
!
USE kinds,            ONLY : DP
USE ions_base,        ONLY : nat
USE fft_base,         ONLY : dfftp
USE uspp_param,       ONLY : nhm
USE modes,            ONLY : npertx, u, t, tmq
USE lr_symm_base,     ONLY : irotmq, minus_q, nsymq, rtau
USE paw_symmetry,     ONLY : paw_dusymmetrize, paw_dumqsymmetrize, &
                             paw_desymmetrize
USE paw_add_symmetry, ONLY : paw_deqsymmetrize
USE noncollin_module, ONLY : noncolin, nspin_mag, domag
USE uspp,             ONLY : okvan
USE paw_variables,    ONLY : okpaw
USE control_ph,       ONLY : lgamma_gamma
USE qpoint,           ONLY : xq
!USE lr_sym_mod,       ONLY : psyme_tpw, psymeq_tpw
USE lr_sym_mod,       ONLY : psymdvscf_tpw, psyme_tpw, psyme_fpol_tpw, &
                             psymeq_tpw

IMPLICIT NONE

INTEGER :: npe, irr, code
COMPLEX(DP) :: drhoscf(dfftp%nnr, nspin_mag, npe)
COMPLEX(DP) :: dbecsum((nhm * (nhm + 1))/2, nat, nspin_mag, npe)

IF (lgamma_gamma) RETURN

SELECT CASE (code) 
   CASE (1)
      CALL psymdvscf_tpw (npe, irr, drhoscf)
      IF (okpaw) THEN
         IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,npe,irr,npertx,irotmq,&
                                                           rtau,xq,tmq)
         CALL PAW_dusymmetrize(dbecsum,npe,irr,npertx,nsymq,rtau,xq,t)
      END IF
   CASE (2)
      CALL psyme_tpw (drhoscf)
!      IF (okpaw) CALL PAW_desymmetrize(dbecsum)
   CASE (3)
      CALL psymeq_tpw (drhoscf)
!      IF (okpaw) CALL PAW_deqsymmetrize(dbecsum)
   CASE (4)
      CALL psyme_fpol_tpw (drhoscf)
!      IF (okpaw) CALL PAW_desymmetrize(dbecsum)
   CASE DEFAULT
      CALL errore('symmetrize_drho','case not programmed',1)
END SELECT

RETURN
END SUBROUTINE symmetrize_drho
