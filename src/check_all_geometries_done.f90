!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE check_all_geometry_done(all_geometry_done)

USE thermo_mod, ONLY : tot_ngeo
USE output, ONLY : fildyn

IMPLICIT NONE

LOGICAL, INTENT(OUT) :: all_geometry_done

INTEGER :: igeom
LOGICAL  :: check_dyn_file_exists
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: auxdyn

all_geometry_done=.TRUE.
DO igeom=1,tot_ngeo
   auxdyn=TRIM(fildyn)//'.g'//TRIM(int_to_char(igeom))//'.'
   IF (all_geometry_done) all_geometry_done=all_geometry_done.AND. &
        check_dyn_file_exists(auxdyn)
   IF (.NOT.all_geometry_done) RETURN
ENDDO
!
!  When start_geometry and last_geometry are used and we arrive here the
!  dynamical matrices for the missing geometries are on file and
!  we read the thermal properties of the geometries not computed in this run
!
CALL collect_all_geometries()

RETURN
END SUBROUTINE check_all_geometry_done

SUBROUTINE collect_all_geometries()

USE thermo_mod,            ONLY : tot_ngeo, ibrav_geo, celldm_geo, no_ph, &
                                  max_geometries, start_geometry, last_geometry
USE cell_base,             ONLY : ibrav, celldm
USE control_ph,            ONLY : ldisp
USE control_thermo,        ONLY : set_internal_path, lq2r
USE output,                ONLY : fildyn
USE io_global,             ONLY : stdout

IMPLICIT NONE
INTEGER :: igeom
CHARACTER(LEN=256) :: auxdyn

DO igeom=1, tot_ngeo
   IF (no_ph(igeom)) CYCLE
   IF ((igeom<start_geometry.OR.igeom>last_geometry.OR.&
                                 max_geometries<tot_ngeo).AND.lq2r) THEN
      WRITE(stdout,'(/,5x,40("%"))')
      WRITE(stdout,'(5x,"Recomputing geometry ", i5)') igeom
      WRITE(stdout,'(5x,40("%"),/)')
      ldisp=.TRUE.
      ibrav=ibrav_geo(igeom)
      celldm(:)=celldm_geo(:,igeom)
      IF (set_internal_path) CALL set_bz_path()
      CALL set_paths_disp()
      CALL set_files_names(igeom)
      auxdyn=fildyn
      CALL manage_ph_dispersions(auxdyn, igeom)
   ENDIF
ENDDO   

RETURN
END SUBROUTINE collect_all_geometries

