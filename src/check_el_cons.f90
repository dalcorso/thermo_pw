! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE check_el_cons()
  !-----------------------------------------------------------------------
  !
  !  This routine tries to read the elastic constants from file for all
  !  the geometries considered in the calculation.
  !  If it finds one geometry it sets el_cons_available to true.
  !  If it finds all geometries it sets also el_cons_geo_available to true.
  !  In the variables el_cons and el_compliances and macro_el it sets
  !  the elastic properties of the central geometry if available or of
  !  the last geometry read if not available.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout
  USE cell_base,         ONLY : ibrav
  USE elastic_constants, ONLY : read_elastic, el_con, el_compliances, &
                                compute_elastic_compliances,  &
                                print_macro_elasticity,  &
                                print_elastic_constants, &
                                print_elastic_compliances
  USE control_elastic_constants, ONLY : el_cons_available, frozen_ions, &
                                        el_cons_geo_available, el_con_geo, &
                                        el_con_d_geo
  USE control_macro_elasticity,  ONLY : macro_el
  USE thermo_mod,        ONLY : tot_ngeo, no_ph
  USE data_files,        ONLY : fl_el_cons
  !
  IMPLICIT NONE
  CHARACTER(LEN=256) :: filelastic
  INTEGER            :: igeo, central_geo
  LOGICAL            :: exst, found(tot_ngeo)
  !
  !   If this is already allocated the elastic constants are already 
  !   available and we do not reread them
  !
  IF (ALLOCATED(el_con_geo)) RETURN

  ALLOCATE( el_con_geo(6,6,tot_ngeo) )
!  ALLOCATE( el_con_d_geo(6,6,tot_ngeo) )
  el_con_geo=0.0_DP
  found=.FALSE.
  el_cons_available=.FALSE.
  el_cons_geo_available=.FALSE.
  DO igeo = 1, tot_ngeo
     CALL add_geometry_number('elastic_constants/', fl_el_cons, &
                                filelastic, igeo)
     CALL read_elastic(filelastic, exst)
     IF (.NOT.exst) CYCLE
     found(igeo)=.TRUE.
     el_con_geo(:,:,igeo)=el_con(:,:)
     WRITE(stdout,'(5x,"Geometry number",i5," elastic constants found")') igeo
     el_cons_available=.TRUE.
  ENDDO
  IF (.NOT.el_cons_available) THEN
     DEALLOCATE(el_con_geo)
     RETURN
  ENDIF
!
!  If there are the elastic constants of the central geometry we
!  set the el_cons and el_compliances of that geometry. Otherwise
!  they are those of the last geometry read.
!
  CALL find_central_geo(tot_ngeo,no_ph,central_geo) 
  IF (found(central_geo)) THEN
     el_con(:,:)=el_con_geo(:,:,central_geo)
     WRITE(stdout,'(/,5x,"Central geometry is number",i5)') central_geo
  ENDIF
  CALL compute_elastic_compliances(el_con,el_compliances)
  CALL print_elastic_constants(el_con, frozen_ions)
  CALL print_elastic_compliances(el_compliances, frozen_ions)
  CALL print_macro_elasticity(ibrav,el_con,el_compliances,macro_el,.TRUE.)
!
!  If the code arrives here we check if the elastic constants have been
!  found for all geometries and in that case set the appropriate flag.
! 
  el_cons_geo_available=.TRUE.
  DO igeo=1,tot_ngeo
     el_cons_geo_available=el_cons_geo_available.AND.found(igeo)
  ENDDO     

  RETURN
END SUBROUTINE check_el_cons

!---------------------------------------------------------------------
SUBROUTINE check_piezo_tensor()
  !-----------------------------------------------------------------------
  !
  !  This routine tries to read the piezoelectric tensor from file for all
  !  the geometries considered in the calculation.
  !  If it finds one geometry it sets piezo_available to true.
  !  If it finds all geometries it sets also piezo_geo_available to true.
  !  It fills the variables g_piezo_tensor_geo, eg_piezo_tensor_geo, 
  !  e_piezo_tensor_geo, d_piezo_tensor_geo
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout
  USE cell_base,         ONLY : ibrav
  USE piezoelectric_tensor, ONLY : read_piezo_tensor, &
                                   g_piezo_tensor,    &
                                   eg_piezo_tensor,   &
                                   e_piezo_tensor,    &
                                   d_piezo_tensor

  USE control_piezoelectric_tensor, ONLY : piezo_available,  &
                                        piezo_geo_available, &
                                        g_piezo_tensor_geo,  &
                                        eg_piezo_tensor_geo, &
                                        e_piezo_tensor_geo,  &
                                        d_piezo_tensor_geo,  &
                                        polar0_geo

  USE thermo_mod,        ONLY : tot_ngeo
  USE data_files,        ONLY : fl_piezo
  !
  IMPLICIT NONE
  CHARACTER(LEN=256) :: filepiezo
  INTEGER            :: igeo
  LOGICAL            :: exst, found(tot_ngeo)
  !
  !   If the piezoelectric tensor is already allocated we do not reread them
  !
  IF (ALLOCATED(g_piezo_tensor_geo)) RETURN

  ALLOCATE( g_piezo_tensor_geo(3,6,tot_ngeo) )
  ALLOCATE( eg_piezo_tensor_geo(3,6,tot_ngeo) )
  ALLOCATE( e_piezo_tensor_geo(3,6,tot_ngeo) )
  ALLOCATE( d_piezo_tensor_geo(3,6,tot_ngeo) )
  ALLOCATE( polar0_geo(3,tot_ngeo) )
  g_piezo_tensor_geo=0.0_DP
  eg_piezo_tensor_geo=0.0_DP
  e_piezo_tensor_geo=0.0_DP
  d_piezo_tensor_geo=0.0_DP
  polar0_geo=0.0_DP
  found=.FALSE.
  piezo_available=.FALSE.
  piezo_geo_available=.FALSE.
  DO igeo = 1, tot_ngeo
     CALL add_geometry_number('elastic_constants/', fl_piezo, &
                                filepiezo, igeo)
     CALL read_piezo_tensor(filepiezo, polar0_geo(1,igeo), exst)
     IF (.NOT.exst) CYCLE
     found(igeo)=.TRUE.
     g_piezo_tensor_geo(:,:,igeo)=g_piezo_tensor(:,:)
     eg_piezo_tensor_geo(:,:,igeo)=eg_piezo_tensor(:,:)
     e_piezo_tensor_geo(:,:,igeo)=e_piezo_tensor(:,:)
     d_piezo_tensor_geo(:,:,igeo)=d_piezo_tensor(:,:)
     WRITE(stdout,'(5x,"Geometry number",i5," piezoelectric tensor found")') &
                                                                          igeo
     piezo_available=.TRUE.
  ENDDO
  IF (.NOT.piezo_available) THEN
     DEALLOCATE(g_piezo_tensor_geo)
     DEALLOCATE(eg_piezo_tensor_geo)
     DEALLOCATE(e_piezo_tensor_geo)
     DEALLOCATE(d_piezo_tensor_geo)
     DEALLOCATE(polar0_geo)
     RETURN
  ENDIF
!
!  If the code arrives here we check if the piezoelectric tensors have been
!  found for all geometries and in that case set the appropriate flag.
! 
  piezo_geo_available=.TRUE.
  DO igeo=1,tot_ngeo
     piezo_geo_available=piezo_geo_available.AND.found(igeo)
  ENDDO     

  RETURN
END SUBROUTINE check_piezo_tensor

