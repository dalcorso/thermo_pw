!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE recover_mod_tpw

  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE

  INTEGER ::  iunrec=99

  PUBLIC :: read_rec_tpw

CONTAINS

  SUBROUTINE read_rec_tpw(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh, dbecsum)
    !
    !  General restart reading routine
    !
    USE kinds, ONLY : DP
    USE ions_base, ONLY : nat
    USE uspp_param, ONLY : nhm
    USE gvecs, ONLY : doublegrid
    USE fft_base, ONLY : dfftp, dffts
    USE fft_interfaces, ONLY : fft_interpolate
    USE uspp,  ONLY : okvan, nlcc_any
    USE lsda_mod, ONLY : nspin
    USE noncollin_module, ONLY : noncolin, nspin_mag, domag
    USE units_ph, ONLY : this_pcxpsi_is_on_file
    USE control_lr, ONLY : ext_recover, convt
    USE efield_mod, ONLY : zstareu0, zstarue0
    USE phus, ONLY : int1, int2
    USE io_files, ONLY : seqopn

    USE lrus, ONLY : int3
    USE eqv,  ONLY : drhos

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter0
    INTEGER, INTENT(IN)  :: npe
    REAL(DP), INTENT(OUT) :: dr2
    COMPLEX(DP), INTENT(OUT) :: dvscfin (dfftp%nnr, nspin_mag, npe)
    COMPLEX(DP), INTENT(OUT) :: dvscfins (dffts%nnr, nspin_mag, npe)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: drhoscfh (dfftp%nnr, nspin_mag, npe)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: dbecsum((nhm*(nhm+1))/2,nat,nspin_mag,npe)

    INTEGER :: is, ipol
    LOGICAL :: exst

    CALL start_clock ('read_rec')
    CALL seqopn (iunrec, 'recover', 'unformatted', exst)
    READ (iunrec) iter0, dr2, convt
    READ (iunrec) this_pcxpsi_is_on_file
    READ (iunrec) zstareu0, zstarue0
    READ (iunrec) dvscfin
    IF (convt.AND.nlcc_any) READ(iunrec) drhoscfh
    IF (convt.AND.ALLOCATED(drhos)) READ(iunrec) drhos
    IF (PRESENT(dbecsum)) READ(iunrec) dbecsum
    IF (okvan) THEN
       READ (iunrec) int1, int2, int3
       IF (noncolin) THEN
          IF (domag) THEN
             CALL set_int12_nc(0)
             CALL compute_int3_coeff(dvscfin, dbecsum, npe)
          ELSE
             CALL set_int12_nc(0)
             CALL set_int3_nc(npe)
          END IF
       END IF
    END IF
    CLOSE (UNIT = iunrec, STATUS = 'keep')
    IF (doublegrid) THEN
       DO is=1,nspin_mag
          DO ipol=1,npe
             CALL fft_interpolate (dfftp, dvscfin(:,is,ipol), dffts, dvscfins(:,is,ipol))
          END DO
       END DO
    END IF
    ext_recover=.FALSE.
    CALL stop_clock ('read_rec')

    RETURN
  END SUBROUTINE read_rec_tpw

END MODULE recover_mod_tpw
