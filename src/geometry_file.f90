! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  This module provides the routines to read and write 
!  a file that contains ngeo strained configurations of a solid 
!  at some pressures. The format of the file is
!  
!  ngeo           ! the number of configurations
!  press(1), energy(1), celldm(.,1), omega(1)  ! the pressure, the energy, 
!                 ! the 6 celldm and the volume of the first configuration
!  press(2), energy(2), celldm(.,2), omega(2)  ! the pressure, the energy, 
!                 ! the 6 celldm and the volume of the second configuration
!  ...
!  press(ngeo), energy(ngeo), celldm(.,ngeo) ! the pressure, the energy, 
!                 ! the 6 celldm and the volume of the ngeo configuration
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
REAL(DP), ALLOCATABLE :: energy_geo_file(:,:)
REAL(DP), ALLOCATABLE :: omega_file(:)
REAL(DP), ALLOCATABLE :: press_file(:)
REAL(DP), ALLOCATABLE :: energy_file(:)
INTEGER :: ngeo_file

REAL(DP) :: vmin_file
REAL(DP) :: b0_file
REAL(DP) :: b01_file
REAL(DP) :: b02_file
REAL(DP) :: emin_file

PUBLIC read_geometry_file, write_geometry_file, set_celldm_geo_from_file, &
       ngeo_file, celldm_geo_file, deallocate_geometry_file,              &
       write_geometry_output, compute_celldm_geo_file, energy_geo_file,   &
       press_file, energy_file, vmin_file, b0_file, b01_file, b02_file,   &
       emin_file, omega_file, do_ev_geometry

CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE read_geometry_file(ngeo)
!-----------------------------------------------------------------------
!
USE data_files, ONLY : flgeom
USE initial_conf, ONLY : ibrav_save
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
USE io_global,  ONLY : meta_ionode, meta_ionode_id
IMPLICIT NONE
INTEGER, INTENT(OUT) :: ngeo(6)
CHARACTER(LEN=256) :: filename
REAL(DP) :: rdum
REAL(DP) :: compute_omega_geo
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
   ALLOCATE(press_file(ngeo_file))
   ALLOCATE(energy_file(ngeo_file))
   ALLOCATE(omega_file(ngeo_file))
   DO igeo=1,ngeo_file
      READ(iu_geom,*) press_file(igeo), energy_file(igeo), &
                      celldm_geo_file(:,igeo), omega_file(igeo)
   ENDDO
   CLOSE(UNIT=iu_geom, STATUS='KEEP')
   CALL mp_bcast(ngeo_file,meta_ionode_id,world_comm)
   CALL mp_bcast(celldm_geo_file,meta_ionode_id,world_comm)
   CALL mp_bcast(press_file,meta_ionode_id,world_comm)
   CALL mp_bcast(energy_file,meta_ionode_id,world_comm)
   CALL mp_bcast(omega_file,meta_ionode_id,world_comm)
ELSE
!
!  the other processors receive the data
!
   CALL mp_bcast(ngeo_file,meta_ionode_id,world_comm)
   ALLOCATE(celldm_geo_file(6,ngeo_file))
   ALLOCATE(press_file(ngeo_file))
   ALLOCATE(energy_file(ngeo_file))
   ALLOCATE(omega_file(ngeo_file))
   CALL mp_bcast(celldm_geo_file,meta_ionode_id,world_comm)
   CALL mp_bcast(press_file,meta_ionode_id,world_comm)
   CALL mp_bcast(energy_file,meta_ionode_id,world_comm)
   CALL mp_bcast(omega_file,meta_ionode_id,world_comm)
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
DO igeo=1,ngeo_file
   omega_file(igeo)=compute_omega_geo(ibrav_save,celldm_geo_file(1,igeo))
ENDDO

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
      WRITE(iu_geom,'(9f13.8)') press_file(igeo), energy_file(igeo), &
                                celldm_geo_file(:,igeo), omega_file(igeo)
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
IF (ALLOCATED(press_file)) DEALLOCATE(press_file)
IF (ALLOCATED(omega_file)) DEALLOCATE(omega_file)
IF (ALLOCATED(energy_file)) DEALLOCATE(energy_file)
  
RETURN
END SUBROUTINE deallocate_geometry_file
!
!-----------------------------------------------------------------------
SUBROUTINE write_geometry_output(npress,press,celldmp,energyp)
!-----------------------------------------------------------------------
!
!   The geometries written on output have the same celldm(1) of the
!   grid points and the other celldm(2-6) optimized so that in the
!   solid there is a uniform pressure. Note that in output we put 
!   the points of the mesh only if their pressure is within the range 
!   of the probed pressures.
!
USE thermo_mod, ONLY : ngeo, celldm_geo
USE initial_conf, ONLY : ibrav_save
IMPLICIT NONE
INTEGER, INTENT(IN) :: npress
REAL(DP), INTENT(IN) :: press(npress), celldmp(6,npress), energyp(npress)

INTEGER :: ipress, igeo, i, ip1, ip2
REAL(DP) :: mind, celldm1_max, celldm1_min, distance
REAL(DP) :: compute_omega_geo

celldm1_max=0.0_DP
celldm1_min=1.0D10
DO ipress=1, npress
   IF (celldmp(1,ipress)>celldm1_max) celldm1_max=celldmp(1,ipress) 
   IF (celldmp(1,ipress)<celldm1_min) celldm1_min=celldmp(1,ipress) 
ENDDO

ngeo_file=ngeo(1)
ALLOCATE(celldm_geo_file(6,ngeo_file))
ALLOCATE(omega_file(ngeo_file))
ALLOCATE(press_file(ngeo_file))
ALLOCATE(energy_file(ngeo_file))
ngeo_file=0
DO igeo=1, ngeo(1)
!
!   Check that celldm_geo(1,igeo) is within the limits of celldm1
!
    IF (celldm_geo(1,igeo)< celldm1_min.OR.celldm_geo(1,igeo)>celldm1_max) &
       CYCLE
!
!  find the point of the mesh with the celldmp(1) closer to
!  celldm_geo(1,igeo)
!
   ip1=0
   mind=1.D8
   DO ipress=1,npress
      distance=ABS(celldmp(1,ipress)-celldm_geo(1,igeo))
      IF (distance<mind) THEN
          ip1=ipress 
          mind=distance
      ENDIF
   ENDDO
!
!  find the second point on the mesh. Here we assume that celldmp(1) 
!  decreases with pressure
!
   IF (celldm_geo(1,igeo) > celldmp(1,ip1)) THEN
      IF (ip1==1) THEN
         ip2=2
      ELSE
         ip2=ip1-1
      ENDIF
   ELSE
      IF (ip1==npress) THEN
         ip2=npress-1
      ELSE
         ip2=ip1+1
      ENDIF
   ENDIF
!   WRITE(6,*) 'found indices ', ip1, ip2
!   WRITE(6,*) 'celldm_geo igeo', igeo, celldm_geo(1, igeo)
!   WRITE(6,*) 'celldmp(ip1), celldmp(ip2) ', celldmp(1,ip1), celldmp(1,ip2)
!
!  Now set celldm_geo_file(1) equal to celldm_geo(1,igeo),
!  and interpolate linearly the pressure
!
   ngeo_file=ngeo_file+1
   celldm_geo_file(1,ngeo_file)=celldm_geo(1,igeo)
   press_file(ngeo_file)=press(ip1) +                      &
                 (press(ip2)-press(ip1)) *                    &
                 (celldm_geo(1,igeo) - celldmp(1,ip1)) /      &
                 (celldmp(1,ip2) - celldmp(1,ip1))
   DO i=2,6
      celldm_geo_file(i,ngeo_file)= celldmp(i,ip1) +       &
                 (celldmp(i,ip2)-celldmp(i,ip1)) *         & 
                 (press_file(ngeo_file) - press(ip1)) /    &
                  (press(ip2)- press(ip1))
   ENDDO
   energy_file(ngeo_file)= energyp(ip1) +       &
                 (energyp(ip2)-energyp(ip1)) *         & 
                 (press_file(ngeo_file) - press(ip1)) /    &
                  (press(ip2)- press(ip1))
   omega_file(ngeo_file)=compute_omega_geo(ibrav_save,&
                                         celldm_geo_file(1,ngeo_file))
ENDDO

CALL do_ev_geometry()
CALL write_geometry_file()
CALL deallocate_geometry_file()

RETURN
END SUBROUTINE write_geometry_output
!
!------------------------------------------------------------------------
SUBROUTINE compute_celldm_geo_file(vmin,celldm0,target_press)
!------------------------------------------------------------------------
!
USE initial_conf, ONLY : ibrav_save
!
IMPLICIT NONE
REAL(DP) :: vmin, celldm0(6)
REAL(DP) :: target_press
REAL(DP) :: compute_omega_geo
REAL(DP) :: volume, pmin_file
INTEGER  :: ip1, ip2, i, ip
!
!  Find the celldm(2-6) at pressure target_press and use those as celldm at
!  vmin. Start by finding the two pressures closest to zero
!
ip1=0
pmin_file=1.D50
DO ip=1,ngeo_file
   IF (ABS (press_file(ip)-target_press)< pmin_file) THEN
      ip1=ip
      pmin_file=ABS(press_file(ip)-target_press)
   ENDIF
ENDDO

IF (press_file(ip) > target_press) THEN
    IF (ip1==ngeo_file) THEN
       ip2=ngeo_file-1
    ELSE
       ip2=ip1+1
    ENDIF
ELSE
   IF (ip1==1) THEN
      ip2=2
   ELSE
      ip2=ip1-1
   ENDIF
ENDIF
!
!  Find the celldm(2-6) from the condition that the pressure is zero
!
DO i=2,6
   celldm0(i)=celldm_geo_file(i,ip1)+(target_press-press_file(ip1)) *  &
             (celldm_geo_file(i,ip2)-celldm_geo_file(i,ip1)) /         &
             (press_file(ip2)-press_file(ip1))
END DO

celldm0(1)=1.0_DP
volume=compute_omega_geo(ibrav_save,celldm0)

celldm0(1)=(vmin/volume)**(1.0_DP/3.0_DP)

RETURN
END SUBROUTINE compute_celldm_geo_file

!------------------------------------------------------------------------
SUBROUTINE do_ev_geometry()
!------------------------------------------------------------------------

USE kinds,  ONLY : DP
USE control_ev,  ONLY : npt, v0, e0, ieos
USE io_global, ONLY : stdout, meta_ionode_id
USE mp,        ONLY : mp_bcast
USE mp_world,  ONLY : world_comm

IMPLICIT NONE
INTEGER :: ipt
CHARACTER(LEN=80) :: eos_label(4)

npt=ngeo_file
ALLOCATE(v0(npt))
ALLOCATE(e0(npt))
DO ipt=1, npt
   v0(ipt)=omega_file(ipt)
   e0(ipt)=energy_file(ipt) 
ENDDO
CALL ev_sub_nodisk(vmin_file, b0_file, b01_file, b02_file, emin_file )
DEALLOCATE(e0)
DEALLOCATE(v0)

CALL mp_bcast(vmin_file, meta_ionode_id, world_comm)
CALL mp_bcast(b0_file, meta_ionode_id, world_comm)
CALL mp_bcast(b01_file, meta_ionode_id, world_comm)
CALL mp_bcast(b02_file, meta_ionode_id, world_comm)
CALL mp_bcast(emin_file, meta_ionode_id, world_comm)

eos_label(1)="Birch-Murnaghan third-order interpolation"
eos_label(2)="Birch-Murnaghan fourth-order interpolation"
eos_label(4)="Murnaghan interpolation"

WRITE(stdout, '(/,1x, 76("-"))') 
WRITE(stdout, '(/,5x, a," of geometry_file data:")') TRIM(eos_label(ieos))

WRITE(stdout, '(/,5x," Equilibrium Volume: ", f13.4, " (a.u.)^3")') vmin_file
WRITE(stdout, '(5x, " Bulk modulus: ", 6x, f13.4, " kbar")') b0_file
WRITE(stdout, '(5x, " dB_0/dp: ", 11x, f13.4)') b01_file
IF (ieos==2) WRITE(stdout, '(5x, " d^2B_0/dp^2: ", 7x, f13.4, " kbar^-1")') &
                                           b02_file
WRITE(stdout, '(5x, " E_min: ", 17x, f13.8," Ry")') emin_file

WRITE(stdout, '(1x, 76("-"))') 

RETURN
END SUBROUTINE do_ev_geometry

END MODULE geometry_file
