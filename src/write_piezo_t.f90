!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE write_piezo_t( )
!----------------------------------------------------------------------
!
!  This routine writes on file the piezoelectric tensor as a function of
!  temperature intepolating them at the volume that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE control_piezoelectric_tensor, ONLY : e_piezo_tensor_geo, lpiezo, lpiezof,&
                                         lpiezo_d, lpiezof_d
USE piezoelectric_tensor, ONLY : write_piezo_tensor_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE control_mur,    ONLY : lmurn
USE anharmonic,     ONLY : celldm_t, e_piezo_tensor_t, d_piezo_tensor_t, &
                           el_comp_t
USE ph_freq_anharmonic, ONLY : celldmf_t, e_piezo_tensorf_t, &
                               d_piezo_tensorf_t, el_compf_t
USE rap_point_group, ONLY : code_group
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filepiezo
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: pt_p1(:,:)
TYPE(poly2), ALLOCATABLE :: pt_p2(:,:)
TYPE(poly3), ALLOCATABLE :: pt_p3(:,:)
TYPE(poly4), ALLOCATABLE :: pt_p4(:,:)
INTEGER :: i, j, idata, itemp
INTEGER :: compute_nwork

ibrav=ibrav_geo(1)
nvar=crystal_parameters(ibrav)
IF (lmurn) nvar=1

ndata=compute_nwork()

ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(pt_p1(3,6))
ALLOCATE(pt_p2(3,6))
ALLOCATE(pt_p3(3,6))
ALLOCATE(pt_p4(3,6))

DO i=1,3
   DO j=1,6
      CALL init_poly(nvar,pt_p1(i,j))
      CALL init_poly(nvar,pt_p2(i,j))
      CALL init_poly(nvar,pt_p3(i,j))
      CALL init_poly(nvar,pt_p4(i,j))
   ENDDO
ENDDO
!
!  Part 1 evaluation of the polynomial coefficients
!
IF (lmurn) THEN
   IF (lgeo_from_file) THEN
      DO idata=1,ndata
         x(1,idata)=omega_geo(idata)
      ENDDO
   ELSE
      DO idata=1,ndata
         x(1,idata)=celldm_geo(1,idata)
      ENDDO
   ENDIF
ELSE
   CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo_eos)
ENDIF
DO i=1,3
   DO j=1,6
      IF (ABS(e_piezo_tensor_geo(i,j,1))>1D-7) THEN
         WRITE(stdout,'(/,5x,"Fitting piezoelectric tensor e(",i4,",",i4,&
                              &")")') i, j
         DO idata=1,ndata
            f(idata)=e_piezo_tensor_geo(i,j,idata)
         END DO

         IF (poly_degree_elc==4) THEN
            CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,pt_p4(i,j)) 
         ELSEIF (poly_degree_elc==3) THEN
            CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,pt_p3(i,j)) 
         ELSEIF (poly_degree_elc==2) THEN
            CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,pt_p2(i,j))      
         ELSEIF (poly_degree_elc==1) THEN
            CALL fit_multi_linear(ndata,nvar,lsolve,x,f,pt_p1(i,j))
         ELSE
            CALL errore('write_piezo_t','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the piezoelectric tensor at the temperature 
!          dependent geometry and computation of the compliances and
!          bulk modulus
!
IF (ltherm_dos) THEN
   CALL interpolate_piezo_tensor(celldm_t, nvar, ibrav, pt_p1, pt_p2, pt_p3, &
               pt_p4, poly_degree_elc, e_piezo_tensor_t)
   lpiezo=.TRUE.
   filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo'
   CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
                               e_piezo_tensor_t, filepiezo, 0, 1)
   CALL compute_d_piezo_tensor_t(e_piezo_tensor_t,el_comp_t,d_piezo_tensor_t,&
                                ntemp)
   lpiezo_d=.TRUE.
   filepiezo='anhar_files/'//TRIM(flanhar)//'.d_piezo'
   CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
                               d_piezo_tensor_t, filepiezo, 0, 2)
ENDIF

IF (ltherm_freq) THEN
   CALL interpolate_piezo_tensor(celldmf_t, nvar, ibrav, pt_p1, pt_p2, pt_p3, &
               pt_p4, poly_degree_elc, e_piezo_tensorf_t)
   lpiezof=.TRUE.
   filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_ph'
   CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
                                    e_piezo_tensorf_t, filepiezo, 0, 1)
   CALL compute_d_piezo_tensor_t(e_piezo_tensorf_t,el_compf_t,&
                                             d_piezo_tensorf_t,ntemp)
   lpiezof_d=.TRUE.
   filepiezo='anhar_files/'//TRIM(flanhar)//'.d_piezo_ph'
   CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
                               d_piezo_tensorf_t, filepiezo, 0, 2)
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DO i=1,3
   DO j=1,6
      CALL clean_poly(pt_p1(i,j))
      CALL clean_poly(pt_p2(i,j))
      CALL clean_poly(pt_p3(i,j))
      CALL clean_poly(pt_p4(i,j))
   ENDDO
ENDDO
DEALLOCATE(pt_p1)
DEALLOCATE(pt_p2)
DEALLOCATE(pt_p3)
DEALLOCATE(pt_p4)

RETURN
END SUBROUTINE write_piezo_t
!
!--------------------------------------------------------------------------
SUBROUTINE write_piezo_p()
!--------------------------------------------------------------------------
!
!  This routine writes on file the piezoelectric tensor as a function of
!  pressure intepolating them at the crystal parameters found at each
!  pressure.
!  This is done only if npress>0, the crystal parameters have 
!  been calculated and the piezoelectric tensor for all geometries are 
!  available on file.
!  This routine is used when lmurn=.FALSE..
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo_eos
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces,    ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces,     ONLY : fit_multi_cubic
USE quartic_surfaces,   ONLY : fit_multi_quartic

USE initial_conf,       ONLY : ibrav_save
USE elastic_constants,  ONLY : write_el_cons_on_file, write_macro_el_on_file, &
                               write_sound_on_file, print_sound_velocities
USE control_elastic_constants, ONLY : el_con_geo, lelastic_p, &
                               el_cons_geo_available
USE lattices,           ONLY : crystal_parameters
USE control_thermo,     ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE control_macro_elasticity, ONLY: macro_el
USE polynomial,         ONLY : poly1, poly2, poly3, poly4, init_poly, &
                               clean_poly
USE data_files,         ONLY : fl_el_cons
USE uniform_pressure,   ONLY : celldm_p, el_cons_p, el_comp_p, b0ec_p, &
                               macro_el_p, density_p, v_p
USE control_pressure,   ONLY : npress, press

IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP),    ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: ec_p1(:,:)
TYPE(poly2), ALLOCATABLE :: ec_p2(:,:)
TYPE(poly3), ALLOCATABLE :: ec_p3(:,:)
TYPE(poly4), ALLOCATABLE :: ec_p4(:,:)
INTEGER :: i, j, idata, itemp, ipress
INTEGER :: compute_nwork

lelastic_p=.FALSE.
IF (npress==0) RETURN
CALL check_el_cons()
IF (.NOT.el_cons_geo_available) RETURN

ibrav=ibrav_geo(1)
nvar=crystal_parameters(ibrav)

ndata=compute_nwork()
ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(ec_p1(6,6))
ALLOCATE(ec_p2(6,6))
ALLOCATE(ec_p3(6,6))
ALLOCATE(ec_p4(6,6))
ALLOCATE(el_cons_p(6,6,npress))
ALLOCATE(el_comp_p(6,6,npress))
ALLOCATE(macro_el_p(8,npress))
ALLOCATE(v_p(3,npress))
ALLOCATE(b0ec_p(npress))

DO i=1,6
   DO j=1,6
      CALL init_poly(nvar,ec_p1(i,j))
      CALL init_poly(nvar,ec_p2(i,j))
      CALL init_poly(nvar,ec_p3(i,j))
      CALL init_poly(nvar,ec_p4(i,j))
   ENDDO
ENDDO
!
!  Part 1 evaluation of the polynomial coefficients
!
CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo_eos)

DO i=1,6
   DO j=1,6
      IF (el_con_geo(i,j,1)>0.1_DP) THEN
         WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")')&
                       i, j
         DO idata=1,ndata
            f(idata)=el_con_geo(i,j,idata)
         END DO

         IF (poly_degree_elc==4) THEN
            CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,ec_p4(i,j)) 
         ELSEIF (poly_degree_elc==3) THEN
            CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,ec_p3(i,j)) 
         ELSEIF (poly_degree_elc==2) THEN
            CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,ec_p2(i,j))      
         ELSEIF (poly_degree_elc==1) THEN
            CALL fit_multi_linear(ndata,nvar,lsolve,x,f,ec_p1(i,j))
         ELSE
            CALL errore('write_elastic_t','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the elastic constants at the pressure
!          dependent geometry and computation of the compliances, 
!          the bulk, young, and shear moduli 
!
CALL interpolate_el_cons_p(celldm_p, nvar, ibrav, ec_p1, ec_p2, ec_p3, &
            ec_p4, poly_degree_elc, el_cons_p, el_comp_p, macro_el_p)

DO ipress=1,npress
   CALL print_sound_velocities(ibrav_save, el_cons_p(:,:,ipress), &
           el_comp_p(:,:,ipress), density_p(ipress), v_p(1,ipress), &
           v_p(2,ipress), v_p(3,ipress),.FALSE.)
ENDDO

b0ec_p(:)= (macro_el_p(1,:)+macro_el_p(5,:)) * 0.5_DP
lelastic_p=.TRUE.
filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.el_cons_p'
CALL write_el_cons_on_file(press, npress, ibrav, laue, el_cons_p, b0ec_p, &
                                                    filelastic, 2)

filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.el_comp_p'
CALL write_el_cons_on_file(press, npress, ibrav, laue, el_comp_p, b0ec_p, & 
                                                       filelastic, 3)

filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.macro_el_p'
CALL write_macro_el_on_file(press, npress, macro_el_p, filelastic, 1)

filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.sound_vel_p'
CALL write_sound_on_file(press, npress, v_p, filelastic, 1)

DEALLOCATE(x)
DEALLOCATE(f)
DO i=1,6
   DO j=1,6
      CALL clean_poly(ec_p1(i,j))
      CALL clean_poly(ec_p2(i,j))
      CALL clean_poly(ec_p3(i,j))
      CALL clean_poly(ec_p4(i,j))
   ENDDO
ENDDO
DEALLOCATE(b0ec_p)
DEALLOCATE(v_p)
DEALLOCATE(macro_el_p)
DEALLOCATE(el_comp_p) 
DEALLOCATE(el_cons_p) 
DEALLOCATE(ec_p1)
DEALLOCATE(ec_p2)
DEALLOCATE(ec_p3)
DEALLOCATE(ec_p4)

RETURN
END SUBROUTINE write_piezo_p
!
!----------------------------------------------------------------------
SUBROUTINE write_piezo_pt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the piezoelectric tensor as a function of
!  temperature interpolating them at the structure that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T. This routine
!  does this for selected pressures.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos
USE rap_point_group, ONLY : code_group
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE piezoelectric_tensor, ONLY : write_piezo_tensor_on_file
USE control_piezoelectric_tensor, ONLY : e_piezo_tensor_geo, lpiezo_pt, &
                           lpiezof_pt, lpiezo_d_pt, lpiezof_d_pt
USE lattices,       ONLY : crystal_parameters
USE control_mur,    ONLY : lmurn
USE control_pressure, ONLY : npress_plot, ipress_plot, press
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE anharmonic_pt,     ONLY : celldm_pt, e_piezo_tensor_pt, &
                                         d_piezo_tensor_pt, el_comp_pt
USE ph_freq_anharmonic_pt,     ONLY : celldmf_pt, e_piezo_tensorf_pt, &
                                      d_piezo_tensorf_pt, el_compf_pt
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filepiezo
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: pt_p1(:,:)
TYPE(poly2), ALLOCATABLE :: pt_p2(:,:)
TYPE(poly3), ALLOCATABLE :: pt_p3(:,:)
TYPE(poly4), ALLOCATABLE :: pt_p4(:,:)
INTEGER :: i, j, idata, itemp, ipressp, ipress
INTEGER :: compute_nwork

IF (npress_plot==0) RETURN

ibrav=ibrav_geo(1)
IF (lmurn) THEN
   nvar=1
ELSE
   nvar=crystal_parameters(ibrav)
ENDIF

ndata=compute_nwork()
ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(pt_p1(3,6))
ALLOCATE(pt_p2(3,6))
ALLOCATE(pt_p3(3,6))
ALLOCATE(pt_p4(3,6))

DO i=1,3
   DO j=1,6
      CALL init_poly(nvar,pt_p1(i,j))
      CALL init_poly(nvar,pt_p2(i,j))
      CALL init_poly(nvar,pt_p3(i,j))
      CALL init_poly(nvar,pt_p4(i,j))
   ENDDO
ENDDO
!
!  Part 1 evaluation of the polynomial coefficients
!
IF (lmurn) THEN
   IF (lgeo_from_file) THEN
      DO idata=1, ndata
         x(1,idata)=omega_geo(idata)
      ENDDO
   ELSE
      DO idata=1, ndata
         x(1,idata)=celldm_geo(1,idata)
      ENDDO
   ENDIF
ELSE
   CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo_eos)
ENDIF

DO i=1,3
   DO j=1,6
      IF (ABS(e_piezo_tensor_geo(i,j,1))>1D-6) THEN
         WRITE(stdout,'(/,5x,"Fitting piezoelectric tensor e(",i4,",",i4,")")')&
                       i, j
         DO idata=1,ndata
            f(idata)=e_piezo_tensor_geo(i,j,idata)
         END DO

         IF (poly_degree_elc==4) THEN
            CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,pt_p4(i,j)) 
         ELSEIF (poly_degree_elc==3) THEN
            CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,pt_p3(i,j)) 
         ELSEIF (poly_degree_elc==2) THEN
            CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,pt_p2(i,j))      
         ELSEIF (poly_degree_elc==1) THEN
            CALL fit_multi_linear(ndata,nvar,lsolve,x,f,pt_p1(i,j))
         ELSE
            CALL errore('write_piezo_pt','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the elastic constants at the temperature 
!          dependent geometry and computation of the compliances and
!          bulk modulus
!
IF (ltherm_dos) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_piezo_tensor(celldm_pt(:,:,ipressp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           e_piezo_tensor_pt(:,:,:,ipressp))
      lpiezo_pt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_press'
      CALL add_value(filepiezo, press(ipress))
      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
           e_piezo_tensor_pt(:,:,:,ipressp), filepiezo, 0, 1)
      CALL compute_d_piezo_tensor_t(e_piezo_tensor_pt(:,:,:,ipressp), &
           el_comp_pt(:,:,:,ipressp),d_piezo_tensor_pt(:,:,:,ipressp), ntemp)
      lpiezo_d_pt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.d_piezo_press'
      CALL add_value(filepiezo, press(ipress))

      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
           d_piezo_tensor_pt(:,:,:,ipressp), filepiezo, 0, 2)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_piezo_tensor(celldmf_pt(:,:,ipressp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           e_piezo_tensorf_pt(:,:,:,ipressp))
      lpiezof_pt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_ph_press'
      CALL add_value(filepiezo, press(ipress))
      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
           e_piezo_tensorf_pt(:,:,:,ipressp), filepiezo, 0, 1)
      CALL compute_d_piezo_tensor_t(e_piezo_tensorf_pt(:,:,:,ipressp), &
           el_compf_pt(:,:,:,ipressp),d_piezo_tensorf_pt(:,:,:,ipressp), ntemp)
      lpiezof_d_pt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.d_piezo_ph_press'
      CALL add_value(filepiezo, press(ipress))
      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav, code_group, &
           d_piezo_tensorf_pt(:,:,:,ipressp), filepiezo, 0, 2)
   ENDDO
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DO i=1,3
   DO j=1,6
      CALL clean_poly(pt_p1(i,j))
      CALL clean_poly(pt_p2(i,j))
      CALL clean_poly(pt_p3(i,j))
      CALL clean_poly(pt_p4(i,j))
   ENDDO
ENDDO
DEALLOCATE(pt_p1)
DEALLOCATE(pt_p2)
DEALLOCATE(pt_p3)
DEALLOCATE(pt_p4)

RETURN
END SUBROUTINE write_piezo_pt
!
!----------------------------------------------------------------------
SUBROUTINE write_piezo_ptt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the piezoelectric tensor as a function of
!  pressure intepolating it at the crystal structure that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T. It does
!  this at selected temperatures.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos

USE rap_point_group, ONLY : code_group
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE piezoelectric_tensor, ONLY : write_piezo_tensor_on_file
USE control_piezoelectric_tensor, ONLY : e_piezo_tensor_geo, lpiezo_ptt, &
                                         lpiezof_ptt, lpiezo_d_ptt,      &
                                         lpiezof_d_ptt
USE lattices,       ONLY : crystal_parameters
USE temperature,    ONLY : ntemp_plot, itemp_plot
USE control_pressure,  ONLY : npress, press
USE control_thermo,    ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE anharmonic_ptt,    ONLY : celldm_ptt, e_piezo_tensor_ptt, &
                              el_comp_ptt, d_piezo_tensor_ptt
USE ph_freq_anharmonic_ptt,  ONLY : celldmf_ptt, e_piezo_tensorf_ptt, &
                              el_compf_ptt, d_piezo_tensorf_ptt
USE control_mur, ONLY : lmurn
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filepiezo
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: pt_p1(:,:)
TYPE(poly2), ALLOCATABLE :: pt_p2(:,:)
TYPE(poly3), ALLOCATABLE :: pt_p3(:,:)
TYPE(poly4), ALLOCATABLE :: pt_p4(:,:)
INTEGER :: i, j, idata, itemp, itempp
INTEGER :: compute_nwork

IF (ntemp_plot==0) RETURN

ibrav=ibrav_geo(1)
IF (lmurn) THEN
   nvar=1
ELSE
   nvar=crystal_parameters(ibrav)
ENDIF

ndata=compute_nwork()
ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(pt_p1(3,6))
ALLOCATE(pt_p2(3,6))
ALLOCATE(pt_p3(3,6))
ALLOCATE(pt_p4(3,6))

DO i=1,3
   DO j=1,6
      CALL init_poly(nvar,pt_p1(i,j))
      CALL init_poly(nvar,pt_p2(i,j))
      CALL init_poly(nvar,pt_p3(i,j))
      CALL init_poly(nvar,pt_p4(i,j))
   ENDDO
ENDDO
!
!  Part 1 evaluation of the polynomial coefficients
!
IF (lmurn) THEN
   IF (lgeo_from_file) THEN
      DO idata=1, ndata
         x(1,idata)=omega_geo(idata)
      ENDDO
   ELSE
      DO idata=1, ndata
         x(1,idata)=celldm_geo(1,idata)
      ENDDO
   ENDIF
ELSE
   CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo_eos)
ENDIF

DO i=1,3
   DO j=1,6
      IF (ABS(e_piezo_tensor_geo(i,j,1))>1.D-6) THEN
         WRITE(stdout,'(/,5x,"Fitting the piezoelectric tensor &
                     & e(",i4,",",i4,")")') i, j
         DO idata=1,ndata
            f(idata)=e_piezo_tensor_geo(i,j,idata)
         END DO

         IF (poly_degree_elc==4) THEN
            CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,pt_p4(i,j)) 
         ELSEIF (poly_degree_elc==3) THEN
            CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,pt_p3(i,j)) 
         ELSEIF (poly_degree_elc==2) THEN
            CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,pt_p2(i,j))      
         ELSEIF (poly_degree_elc==1) THEN
            CALL fit_multi_linear(ndata,nvar,lsolve,x,f,pt_p1(i,j))
         ELSE
            CALL errore('write_piezo_ptt','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the elastic constants at the temperature 
!          dependent geometry and computation of the compliances and
!          bulk modulus
!
IF (ltherm_dos) THEN
   DO itempp=1, ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_piezo_tensor_p(celldm_ptt(:,:,itempp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           e_piezo_tensor_ptt(:,:,:,itempp))
      lpiezo_ptt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_temp'
      CALL add_value(filepiezo, temp(itemp))
      CALL write_piezo_tensor_on_file(press, npress, ibrav, code_group, &
           e_piezo_tensor_ptt(:,:,:,itempp), filepiezo, 2, 1)
      CALL compute_d_piezo_tensor_t(e_piezo_tensor_ptt(:,:,:,itempp), &
           el_comp_ptt(:,:,:,itempp),d_piezo_tensor_ptt(:,:,:,itempp), npress)
      lpiezo_d_ptt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.d_piezo_temp'
      CALL add_value(filepiezo, temp(itemp))
      CALL write_piezo_tensor_on_file(press, npress, ibrav, code_group, &
           d_piezo_tensor_ptt(:,:,:,itempp), filepiezo, 2, 2)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO itempp=1, ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_piezo_tensor_p(celldmf_ptt(:,:,itempp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           e_piezo_tensorf_ptt(:,:,:,itempp))
      lpiezof_ptt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_ph_temp'
      CALL add_value(filepiezo, temp(itemp))
      CALL write_piezo_tensor_on_file(press, npress, ibrav, code_group, &
           e_piezo_tensorf_ptt(:,:,:,itempp), filepiezo, 2, 1)
      CALL compute_d_piezo_tensor_t(e_piezo_tensorf_ptt(:,:,:,itempp), &
           el_compf_ptt(:,:,:,itempp),d_piezo_tensorf_ptt(:,:,:,itempp), npress)
      lpiezof_d_ptt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.d_piezo_ph_temp'
      CALL add_value(filepiezo, temp(itemp))
      CALL write_piezo_tensor_on_file(press, npress, ibrav, code_group, &
           d_piezo_tensorf_ptt(:,:,:,itempp), filepiezo, 2, 2)
   ENDDO
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)

DO i=1,3
   DO j=1,6
      CALL clean_poly(pt_p1(i,j))
      CALL clean_poly(pt_p2(i,j))
      CALL clean_poly(pt_p3(i,j))
      CALL clean_poly(pt_p4(i,j))
   ENDDO
ENDDO
DEALLOCATE(pt_p1)
DEALLOCATE(pt_p2)
DEALLOCATE(pt_p3)
DEALLOCATE(pt_p4)

RETURN
END SUBROUTINE write_piezo_ptt
