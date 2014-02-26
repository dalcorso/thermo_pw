!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE bcast_thermo_input()
  !-----------------------------------------------------------------------
  !
  !  This routine broadcasts to all the images the input of thermo_pw.
  !
  USE thermo_mod,      ONLY : what, ngeo, step_ngeo, ntry
  USE control_mur,     ONLY : vmin_input, vmax_input, deltav, nvol
  USE control_thermo,  ONLY : outdir_thermo, flevdat,    &
                              flfrc, flfrq, fldos, fltherm, flanhar, &
                              filband, flkeconv, flnkconv, flgrun
  USE temperature,     ONLY : tmin, tmax, deltat, ntemp
  USE ifc,             ONLY : nq1_d, nq2_d, nq3_d, ndos_input, deltafreq, zasr, &
                              freqmin_input, freqmax_input
  USE control_paths,   ONLY : q_in_band_form, q_in_cryst_coord, q2d, &
                              point_label_type
  USE control_bands,   ONLY : flpband, emin_input, emax_input, nbnd_bands, lsym
  USE control_grun,    ONLY : flpgrun
  USE control_gnuplot, ONLY : flgnuplot, flpsband, flpsdisp, &
                              flpsdisp, flpsdos, flpstherm, &
                              flpsanhar, flpsmur, flpskeconv, flpsnkconv, &
                              flpsgrun, &
                              lgnuplot, gnuplot_command
  USE control_conv,    ONLY : nke, deltake, nkeden, deltakeden, &
                              nnk, deltank, nsigma, deltasigma
  USE mp_world,        ONLY : world_comm
  USE mp,              ONLY : mp_bcast
  USE io_global,       ONLY : meta_ionode_id
  !
  IMPLICIT NONE
  !
  CALL mp_bcast( what, meta_ionode_id, world_comm )
  CALL mp_bcast( ngeo, meta_ionode_id, world_comm )
  CALL mp_bcast( ntry, meta_ionode_id, world_comm )
  CALL mp_bcast( step_ngeo, meta_ionode_id, world_comm )
  CALL mp_bcast( zasr, meta_ionode_id, world_comm )
  CALL mp_bcast( flfrc, meta_ionode_id, world_comm )
  CALL mp_bcast( flfrq, meta_ionode_id, world_comm )
  CALL mp_bcast( fldos, meta_ionode_id, world_comm )
  CALL mp_bcast( fltherm, meta_ionode_id, world_comm )
  CALL mp_bcast( flanhar, meta_ionode_id, world_comm )
  CALL mp_bcast( flkeconv, meta_ionode_id, world_comm )
  CALL mp_bcast( flnkconv, meta_ionode_id, world_comm )
  CALL mp_bcast( flgrun, meta_ionode_id, world_comm )
  CALL mp_bcast( flpgrun, meta_ionode_id, world_comm )
  CALL mp_bcast( nq1_d, meta_ionode_id, world_comm )
  CALL mp_bcast( nq2_d, meta_ionode_id, world_comm )
  CALL mp_bcast( nq3_d, meta_ionode_id, world_comm )
  CALL mp_bcast( freqmin_input, meta_ionode_id, world_comm )
  CALL mp_bcast( freqmax_input, meta_ionode_id, world_comm )
  CALL mp_bcast( ndos_input, meta_ionode_id, world_comm )
  CALL mp_bcast( deltafreq, meta_ionode_id, world_comm )
  CALL mp_bcast( nbnd_bands, meta_ionode_id, world_comm )
  CALL mp_bcast( lsym, meta_ionode_id, world_comm )
  CALL mp_bcast( vmin_input, meta_ionode_id, world_comm )
  CALL mp_bcast( vmax_input, meta_ionode_id, world_comm )
  CALL mp_bcast( deltav, meta_ionode_id, world_comm )
  CALL mp_bcast( nvol, meta_ionode_id, world_comm )
  CALL mp_bcast( nke, meta_ionode_id, world_comm )
  CALL mp_bcast( deltake, meta_ionode_id, world_comm )
  CALL mp_bcast( nkeden, meta_ionode_id, world_comm )
  CALL mp_bcast( deltakeden, meta_ionode_id, world_comm )
  CALL mp_bcast( nnk, meta_ionode_id, world_comm )
  CALL mp_bcast( deltank, meta_ionode_id, world_comm )
  CALL mp_bcast( nsigma, meta_ionode_id, world_comm )
  CALL mp_bcast( deltasigma, meta_ionode_id, world_comm )
  CALL mp_bcast( tmin, meta_ionode_id, world_comm )
  CALL mp_bcast( tmax, meta_ionode_id, world_comm )
  CALL mp_bcast( deltat, meta_ionode_id, world_comm )
  CALL mp_bcast( ntemp, meta_ionode_id, world_comm )
  CALL mp_bcast( flpband, meta_ionode_id, world_comm )
  CALL mp_bcast( filband, meta_ionode_id, world_comm )
  CALL mp_bcast( flgnuplot, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsmur, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsband, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsdos, meta_ionode_id, world_comm )
  CALL mp_bcast( flpstherm, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsanhar, meta_ionode_id, world_comm )
  CALL mp_bcast( flpskeconv, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsnkconv, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsgrun, meta_ionode_id, world_comm )
  CALL mp_bcast( emin_input, meta_ionode_id, world_comm )
  CALL mp_bcast( emax_input, meta_ionode_id, world_comm )
  CALL mp_bcast( flevdat, meta_ionode_id, world_comm )
  CALL mp_bcast( q_in_band_form, meta_ionode_id, world_comm )
  CALL mp_bcast( q_in_cryst_coord, meta_ionode_id, world_comm )
  CALL mp_bcast( point_label_type, meta_ionode_id, world_comm )
  CALL mp_bcast( q2d, meta_ionode_id, world_comm )
  CALL mp_bcast( lgnuplot, meta_ionode_id, world_comm )
  CALL mp_bcast( gnuplot_command, meta_ionode_id, world_comm )

  RETURN
END SUBROUTINE bcast_thermo_input
