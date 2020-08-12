!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
SUBROUTINE collect_all_geometries()
!--------------------------------------------------------
!
!  This routine should be called only when the dynamical matrices
!  for all geometries are found on disk. It checks if the
!  files with the harmonic thermodynamic quantities are on disk and if
!  it does not find them it recomputes the missing quantities.
!
USE thermo_mod,            ONLY : tot_ngeo, ibrav_geo, celldm_geo, no_ph, &
                                  max_geometries, start_geometry, last_geometry
USE control_elastic_constants, ONLY : start_geometry_qha, last_geometry_qha, &
                                  ngeom
USE cell_base,             ONLY : ibrav, celldm
USE control_thermo,        ONLY : set_internal_path, lq2r
USE output,                ONLY : fildyn
USE io_global,             ONLY : stdout

IMPLICIT NONE
INTEGER :: igeom, igeom_qha, iwork, work_base
CHARACTER(LEN=256) :: auxdyn

work_base=tot_ngeo/ngeom
DO igeom_qha=start_geometry_qha, last_geometry_qha
   DO iwork=1,work_base
      igeom=(igeom_qha-1)*work_base+iwork
      IF (no_ph(igeom)) CYCLE
      IF ((igeom<start_geometry.OR.igeom>last_geometry.OR.&
                                 max_geometries<tot_ngeo).AND.lq2r) THEN
         WRITE(stdout,'(/,5x,40("%"))')
         WRITE(stdout,'(5x,"Recomputing geometry ", i5)') igeom
         WRITE(stdout,'(5x,40("%"),/)')
         ibrav=ibrav_geo(igeom)
         celldm(:)=celldm_geo(:,igeom)
         IF (set_internal_path) CALL set_bz_path()
         CALL set_files_names(igeom)
         auxdyn=fildyn
         CALL manage_ph_postproc(auxdyn, igeom)
      ENDIF
   ENDDO   
ENDDO

RETURN
END SUBROUTINE collect_all_geometries

