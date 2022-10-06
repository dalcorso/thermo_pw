! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  This module provides the routines to read and write 
!  a file that contains a certain number of strained
!  configurations of a solid. The format of the file is
!  
!  ngeo           ! the number of configurations
!  celldm(.,1)    ! the 6 celldm of the first configuration
!  celldm(.,2)    ! the 6 celldm of the second configuration
!  ...
!  celldm(.,ngeo) ! the 6 celldm of the ngeo configuration
!
!  Atomic coordinates are not saved. They are obtained by 
!  straining uniformly those read from the input file.
!
MODULE geometry_file

USE kinds, ONLY : DP

IMPLICIT NONE
SAVE
PRIVATE

REAL(DP), ALLOCATABLE :: celldm_geo_file(:,:)
REAL(DP), ALLOCATABLE :: press_file(:)
INTEGER :: ngeo_file

PUBLIC read_geometry_file, write_geometry_file, set_celldm_geo_from_file, &
       ngeo_file, celldm_geo_file, deallocate_geometry_file, &
       write_geometry_output

CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE read_geometry_file(ngeo)
!-----------------------------------------------------------------------
!
USE data_files, ONLY : flgeom
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
USE io_global,  ONLY : meta_ionode, meta_ionode_id
IMPLICIT NONE
INTEGER, INTENT(OUT) :: ngeo(6)
CHARACTER(LEN=256) :: filename
REAL(DP) :: rdum
INTEGER :: iu_geom, igeo, ios
INTEGER :: find_free_unit

filename='./'//TRIM(flgeom)//'.dat'
IF (meta_ionode) THEN
!
!  meta_ionode reads the data and broadcast them to all nodes
!
   iu_geom=find_free_unit()
   OPEN (UNIT=iu_geom, FILE=TRIM(filename), STATUS='old',&
                                         FORM='formatted', ERR=100, IOSTAT=ios)
   READ(iu_geom,*) ngeo_file
   ALLOCATE(celldm_geo_file(6,ngeo_file))
   DO igeo=1,ngeo_file
      READ(iu_geom,*) rdum, celldm_geo_file(:,igeo)
   ENDDO
   CLOSE(UNIT=iu_geom, STATUS='KEEP')
   CALL mp_bcast(ngeo_file,meta_ionode_id,world_comm)
   CALL mp_bcast(celldm_geo_file,meta_ionode_id,world_comm)
ELSE
!
!  the other processors receive the data
!
   CALL mp_bcast(ngeo_file,meta_ionode_id,world_comm)
   ALLOCATE(celldm_geo_file(6,ngeo_file))
   CALL mp_bcast(celldm_geo_file,meta_ionode_id,world_comm)
ENDIF
100 CALL mp_bcast(ios,meta_ionode_id,world_comm)
CALL errore('read_geometry_file','Problem with the geometry file',ios)
!
!  all nodes set here the number of geometries to pass to the main code
!  the actual geometries are not copied here in celldm_geo since
!  this array is not yet allocated.
!
ngeo=1
ngeo(1)=ngeo_file

RETURN
END SUBROUTINE read_geometry_file
!
!----------------------------------------------------------------------
SUBROUTINE write_geometry_file()
!----------------------------------------------------------------------
!
!  This routine writes the file with the geometries. It assumes that the data
!  are in the variables of this module ngeo_file and celldm_geo_file,
!  so the calling routine must set these data. Only meta_ionode writes
!  the file
!
USE data_files,  ONLY : flgeom
USE io_global,   ONLY : meta_ionode
IMPLICIT NONE
CHARACTER(LEN=256) :: filename
INTEGER :: iu_geom, igeo
INTEGER :: find_free_unit

filename='./'//TRIM(flgeom)//'.dat'
IF (meta_ionode) THEN
   iu_geom=find_free_unit()
   OPEN (UNIT=iu_geom, FILE=TRIM(filename), STATUS='unknown',&
                                                     FORM='formatted')
   WRITE(iu_geom,'(i8)') ngeo_file
   DO igeo=1,ngeo_file
      WRITE(iu_geom,'(7f13.8)') press_file(igeo), celldm_geo_file(:,igeo)
   ENDDO
   CLOSE(UNIT=iu_geom, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_geometry_file
!
!----------------------------------------------------------------------
SUBROUTINE set_celldm_geo_from_file(celldm_geo, ngeo)
!----------------------------------------------------------------------
!
!  celldm_geo is copied here and not read directly because we
!  need to allocate it outside this routine before setting it.
!
USE data_files, ONLY : flgeom
USE io_global,  ONLY : meta_ionode
IMPLICIT NONE
INTEGER, INTENT(IN) :: ngeo
REAL(DP), INTENT(INOUT) :: celldm_geo(6,ngeo)

INTEGER :: igeo

DO igeo=1,ngeo
   celldm_geo(:,igeo)=celldm_geo_file(:,igeo)
ENDDO

RETURN
END SUBROUTINE set_celldm_geo_from_file
!
!-----------------------------------------------------------------------
SUBROUTINE deallocate_geometry_file()
!-----------------------------------------------------------------------
IMPLICIT NONE

IF (ALLOCATED(celldm_geo_file)) DEALLOCATE(celldm_geo_file)
  
END SUBROUTINE deallocate_geometry_file
!
!-----------------------------------------------------------------------
SUBROUTINE write_geometry_output(npress,press,celldmp)
!-----------------------------------------------------------------------
!
!   The geometries written on output have the same celldm(1) of the
!   grid points and the other celldm(2-6) optimized to minimize the
!   energy. Note that in output we put the points of the mesh only if
!   their pressure is within the range of the probed pressures.
!
USE thermo_mod, ONLY : ngeo, celldm_geo
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), celldmp(6,npress)

INTEGER :: ipress, igeo, i, ip1, ip2
REAL(DP) :: mind

ngeo_file=ngeo(1)
ALLOCATE(celldm_geo_file(6,ngeo_file))
ALLOCATE(press_file(ngeo_file))
ngeo_file=0
DO igeo=1, ngeo(1)
!
!  find the two points of the mesh with the celldmp(1) closer to
!  celldm_geo(1,igeo)
!
   ip1=0
   mind=1.D8
   DO ipress=1,npress
      IF ((celldmp(1,ipress)<celldm_geo(1,igeo)).AND. &
          ABS(celldmp(1,ipress)-celldm_geo(1,igeo))<mind) THEN
          ip1=ipress 
          mind=ABS(celldmp(1,ipress)-celldm_geo(1,igeo))
      ENDIF
   ENDDO
!
!  If they exist set celldm_geo_file(1) equal to celldm_geo(1,igeo),
!  and interpolate linearly celldm_geo_file(2-6) and the pressure
!
   IF (ip1>0.AND.ip1<npress) THEN
      ip2=ip1+1
      ngeo_file=ngeo_file+1
      celldm_geo_file(1,ngeo_file)=celldm_geo(1,igeo)
      DO i=2,6
         celldm_geo_file(i,ngeo_file)= celldmp(i,ip1) +       &
                 (celldmp(i,ip2)-celldmp(i,ip1)) *            & 
                 (celldm_geo(1,igeo) - celldmp(1,ip1)) /      &
                 (celldmp(1,ip2) - celldmp(1,ip1)) 
      ENDDO
      press_file(ngeo_file)=press(ip1) +                      &
                 (press(ip2)-press(ip1)) *                    &
                 (celldm_geo(1,igeo) - celldmp(1,ip1)) /      &
                 (celldmp(1,ip2) - celldmp(1,ip1))
   ENDIF
ENDDO

CALL write_geometry_file()
DEALLOCATE(press_file)
DEALLOCATE(celldm_geo_file)

RETURN
END SUBROUTINE write_geometry_output

END MODULE geometry_file
