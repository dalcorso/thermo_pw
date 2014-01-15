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
  USE thermo_mod,      ONLY : what
  USE control_thermo,  ONLY : outdir_thermo, file_ev_dat,    &
                              flfrc, flfrq, fldos, fltherm, flanhar, &
                              filband, spin_component
  USE thermodynamics,  ONLY : ngeo, tmin, tmax, deltat, ntemp
  USE ifc,             ONLY : nq1_d, nq2_d, nq3_d, ndos, deltae, zasr
  USE control_thermo,  ONLY : file_ev_dat
  USE control_paths,   ONLY : q_in_band_form, q_in_cryst_coord, q2d, &
                              point_label_type, nbnd_bands
  USE control_bands,   ONLY : flpband, emin_input, emax_input
  USE control_gnuplot, ONLY : flgnuplot, flpsband, flpsdisp, &
                              flpsdisp, flpsdos, flpstherm, &
                              flpsanhar 
  USE mp_world,        ONLY : world_comm
  USE mp,              ONLY : mp_bcast
  USE io_global,       ONLY : meta_ionode_id
  !
  IMPLICIT NONE
  !
  CALL mp_bcast( what, meta_ionode_id, world_comm )
  CALL mp_bcast( ngeo, meta_ionode_id, world_comm )
  CALL mp_bcast( zasr, meta_ionode_id, world_comm )
  CALL mp_bcast( flfrc, meta_ionode_id, world_comm )
  CALL mp_bcast( flfrq, meta_ionode_id, world_comm )
  CALL mp_bcast( fldos, meta_ionode_id, world_comm )
  CALL mp_bcast( fltherm, meta_ionode_id, world_comm )
  CALL mp_bcast( flanhar, meta_ionode_id, world_comm )
  CALL mp_bcast( nq1_d, meta_ionode_id, world_comm )
  CALL mp_bcast( nq2_d, meta_ionode_id, world_comm )
  CALL mp_bcast( nq3_d, meta_ionode_id, world_comm )
  CALL mp_bcast( ndos, meta_ionode_id, world_comm )
  CALL mp_bcast( deltae, meta_ionode_id, world_comm )
  CALL mp_bcast( nbnd_bands, meta_ionode_id, world_comm )
  CALL mp_bcast( tmin, meta_ionode_id, world_comm )
  CALL mp_bcast( tmax, meta_ionode_id, world_comm )
  CALL mp_bcast( deltat, meta_ionode_id, world_comm )
  CALL mp_bcast( ntemp, meta_ionode_id, world_comm )
  CALL mp_bcast( flpband, meta_ionode_id, world_comm )
  CALL mp_bcast( filband, meta_ionode_id, world_comm )
  CALL mp_bcast( flgnuplot, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsband, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsdos, meta_ionode_id, world_comm )
  CALL mp_bcast( flpstherm, meta_ionode_id, world_comm )
  CALL mp_bcast( flpsanhar, meta_ionode_id, world_comm )
  CALL mp_bcast( spin_component, meta_ionode_id, world_comm )
  CALL mp_bcast( emin_input, meta_ionode_id, world_comm )
  CALL mp_bcast( emax_input, meta_ionode_id, world_comm )
  CALL mp_bcast( file_ev_dat, meta_ionode_id, world_comm )
  CALL mp_bcast( q_in_band_form, meta_ionode_id, world_comm )
  CALL mp_bcast( q_in_cryst_coord, meta_ionode_id, world_comm )
  CALL mp_bcast( point_label_type, meta_ionode_id, world_comm )
  CALL mp_bcast( q2d, meta_ionode_id, world_comm )

  RETURN
END SUBROUTINE bcast_thermo_input
