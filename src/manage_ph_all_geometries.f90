!
! Copyright (C) 2017-2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE manage_ph_all_geometries()
!--------------------------------------------------------------------------
!
!  This routine controls the calculation of the phonon dispersions 
!  when geometries are made in parallel. It first reads the input files
!  of all the pw.x calculations and computes all the tasks to do.
!  Then it runs asynchronously all the tasks and finally collects the 
!  results and computes the thermodynamic properties for all geometries.
!
  USE thermo_mod,            ONLY : no_ph, start_geometry, last_geometry
  USE control_phrun,         ONLY : auxdyn
  USE control_thermo,        ONLY : after_disp, lq2r
  USE distribute_collection, ONLY : me_igeom
!
!   pw.x and ph.x variables
!
  USE io_global,             ONLY : stdout
  USE control_ph,            ONLY : always_run
  USE output,                ONLY : fildyn
!
!  parallelization variables
!
  USE mp_asyn,               ONLY : stop_signal_activated
  USE mp,                    ONLY : mp_barrier
  USE mp_world,              ONLY : world_comm

IMPLICIT NONE

INTEGER           :: part, nwork, igeom
CHARACTER(LEN=80) :: message
LOGICAL           :: something_todo

always_run=.TRUE.
CALL start_clock( 'PHONON' )
!
!   First check if for some geometry the dynamical matrices are already
!   available. This check works only if fildyn has been set in the input
!   thermo_control.
!
CALL check_dynmat_all_geo_on_file()
!
IF (after_disp) GOTO 50   ! Skip all the phonon calculation if after_disp 
                          ! is true.
!
!  Initialize the work of all the geometries and analyze what is on disk
!  This routine must be called by all processors
!
CALL check_initial_status_all_geo(something_todo)
!
IF (something_todo) THEN
!
!  Initialize the asynchronous work
!
   part=2
   CALL initialize_thermo_work(nwork, part)
!
!  Asynchronous work starts here. No communication is
!  allowed except through the master/slave mechanism
!  
   CALL run_thermo_asynchronously(nwork, part, 1)

   IF (stop_signal_activated) GOTO 100
ENDIF
!
!  Now all calculations are done, we collect the results. 
!  Each image acts independently. The total number of collection tasks
!  are divided among images. The processors must be resynchronized here
!  otherwise some partial dynamical matrix could be missing.
!
CALL mp_barrier(world_comm)
CALL divide_all_collection_work()
!
!  loop over all geometries and skip those not done by this image
!
DO igeom=start_geometry, last_geometry
   IF (.NOT.me_igeom(igeom)) CYCLE
   WRITE(message,'(5x,"Collecting geometry ", i5)') igeom
   CALL decorated_write(message)
   CALL set_outdir_name(igeom)
   !
   ! reads the pw.x output and the phonon input. With the flag 1 check
   ! initial geometry does not recheck what is on disk. Here all parts 
   ! of the dynamical matrices must have been computed by some image.
   !
   CALL fast_phq_readin(.TRUE., igeom)
   CALL set_fildyn_name(igeom)
   CALL check_initial_geometry(auxdyn,1)
   !
   ! and collect all the results
   !
   CALL manage_collection(auxdyn, igeom)
   !
   ! deallocate the phonon variables and close pw and ph
   !
   CALL close_ph_geometry(.TRUE.)
ENDDO
!
!  deallocate the variables needed to divide the collection work among
!  images.
!
CALL clean_collection_work()
50 CONTINUE
!
!  resynchronize all processors, otherwise some dynamical matrix could be
!  missing. The calculation of the thermodynamical properties is
!  parallelized over all processors, so the following routines must
!  be called by all images.
!
CALL mp_barrier(world_comm)
!
!  a second loop over the geometries. We compute here the thermodynamic
!  properties
!
DO igeom=start_geometry, last_geometry
   IF (no_ph(igeom)) CYCLE
   WRITE(message,'(5x,"Computing thermodynamic properties", i5)') igeom
   CALL decorated_write(message)

   CALL set_files_names(igeom)
   auxdyn=fildyn
!
!  Compute the phonon dispersions and the thermodynamic properties
!
   IF (lq2r) CALL manage_ph_postproc(igeom)
ENDDO

100 CONTINUE
CALL restore_files_names()

RETURN
END SUBROUTINE manage_ph_all_geometries

