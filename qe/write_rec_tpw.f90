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

  SUBROUTINE read_rec_tpw(dr2, iter0, dfpt_data)
    !
    !  General restart reading routine
    !
    USE kinds, ONLY : DP
    USE ions_base, ONLY : nat
    USE uspp_param, ONLY : nhm
    USE gvecs, ONLY : doublegrid
    USE fft_base, ONLY : dfftp, dffts
    USE fft_interfaces, ONLY : fft_interpolate
    USE paw_variables,  ONLY : okpaw
    USE uspp,  ONLY : okvan, nlcc_any
    USE lsda_mod, ONLY : nspin
    USE noncollin_module, ONLY : noncolin, nspin_mag, domag
    USE units_ph, ONLY : this_pcxpsi_is_on_file
    USE control_lr, ONLY : ext_recover, convt
    USE efield_mod, ONLY : zstareu0, zstarue0
    USE phus, ONLY : int1, int2
    USE io_files, ONLY : seqopn
    USE dfpt_type, ONLY : dfpt_data_type

    USE lrus, ONLY : int3

    IMPLICIT NONE
    TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
    !! Output: Data that describes linear response quantities

    INTEGER, INTENT(OUT) :: iter0
    REAL(DP), INTENT(OUT) :: dr2

    INTEGER :: is, ipol, npe
    LOGICAL :: exst

    CALL start_clock ('read_rec')
    npe=dfpt_data%npert
    CALL seqopn (iunrec, 'recover', 'unformatted', exst)
    READ (iunrec) iter0, dr2, convt
    READ (iunrec) this_pcxpsi_is_on_file
    READ (iunrec) zstareu0, zstarue0
    READ (iunrec) dfpt_data%dvscfp
    IF (convt.AND.nlcc_any) READ(iunrec) dfpt_data%drhop
    IF (convt.AND.ALLOCATED(dfpt_data%drhos)) READ(iunrec) dfpt_data%drhos
    IF (okpaw) READ(iunrec) dfpt_data%dbecsum
    IF (okvan) THEN
       READ (iunrec) int1, int2, int3
       IF (noncolin) THEN
          IF (domag) THEN
             CALL set_int12_nc(0)
             CALL compute_int3_coeff(dfpt_data%dvscfp, dfpt_data%dbecsum, npe)
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
             CALL fft_interpolate (dfftp, dfpt_data%dvscfp(:,is,ipol), &
                                   dffts, dfpt_data%dvscfs(:,is,ipol))
          END DO
       END DO
    END IF
    ext_recover=.FALSE.
    CALL stop_clock ('read_rec')

    RETURN
  END SUBROUTINE read_rec_tpw

END MODULE recover_mod_tpw
