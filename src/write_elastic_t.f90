!
! Copyright (C) 2016-2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_elastic_t( )
!
!  This routine writes on file the elastic constants as a function of
!  temperature intepolating them at the volume that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE elastic_constants, ONLY : write_el_cons_on_file
USE control_elastic_constants, ONLY : el_con_geo, lelastic, lelasticf
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_macro_elasticity, ONLY: macro_el
USE anharmonic,     ONLY : celldm_t, el_cons_t, el_comp_t, b0_t
USE ph_freq_anharmonic, ONLY : celldmf_t, el_consf_t, el_compf_t, b0f_t
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: ec_p1(:,:)
TYPE(poly2), ALLOCATABLE :: ec_p2(:,:)
TYPE(poly3), ALLOCATABLE :: ec_p3(:,:)
TYPE(poly4), ALLOCATABLE :: ec_p4(:,:)
INTEGER :: i, j, idata, itemp
INTEGER :: compute_nwork

ibrav=ibrav_geo(1)
nvar=crystal_parameters(ibrav)

ndata=compute_nwork()
ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(ec_p1(6,6))
ALLOCATE(ec_p2(6,6))
ALLOCATE(ec_p3(6,6))
ALLOCATE(ec_p4(6,6))

DO i=1,6
   DO j=i,6
      CALL init_poly(nvar,ec_p1(i,j))
      CALL init_poly(nvar,ec_p2(i,j))
      CALL init_poly(nvar,ec_p3(i,j))
      CALL init_poly(nvar,ec_p4(i,j))
   ENDDO
ENDDO
!
!  Part 1 evaluation of the polynomial coefficients
!
CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)

DO i=1,6
   DO j=i,6
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
!  Part 2: interpolation of the elastic constants at the temperature 
!          dependent geometry and computation of the compliances and
!          bulk modulus
!
IF (ltherm_dos) THEN
   CALL interpolate_el_cons(celldm_t, nvar, ibrav, ec_p1, ec_p2, ec_p3, &
               ec_p4, poly_degree_elc, el_cons_t, el_comp_t, b0_t)
   lelastic=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_cons_t, b0_t, &
                                                       filelastic, 0)
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_comp_t, b0_t, & 
                                                       filelastic, 1)
   
ENDIF

IF (ltherm_freq) THEN
   CALL interpolate_el_cons(celldmf_t, nvar, ibrav, ec_p1, ec_p2, ec_p3, &
               ec_p4, poly_degree_elc, el_consf_t, el_compf_t, b0f_t)
   lelasticf=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_ph'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_consf_t, b0f_t, &
                                                          filelastic, 0)

   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp_ph'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_compf_t, b0f_t, &
                                                           filelastic,1)
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DO i=1,6
   DO j=i,6
      CALL clean_poly(ec_p1(i,j))
      CALL clean_poly(ec_p2(i,j))
      CALL clean_poly(ec_p3(i,j))
      CALL clean_poly(ec_p4(i,j))
   ENDDO
ENDDO
DEALLOCATE(ec_p1)
DEALLOCATE(ec_p2)
DEALLOCATE(ec_p3)
DEALLOCATE(ec_p4)

RETURN
END SUBROUTINE write_elastic_t
!
!--------------------------------------------------------------------------
SUBROUTINE write_elastic_p()
!--------------------------------------------------------------------------
!
!  This routine writes on file the elastic constants as a function of
!  pressure intepolating them at the crystal parameters found at each
!  pressure.
!  This is done only if npress>0 and the crystal parameters have 
!  been calculated and the elastic constant for all geometries are 
!  available on file.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces,    ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces,     ONLY : fit_multi_cubic
USE quartic_surfaces,   ONLY : fit_multi_quartic

USE elastic_constants,  ONLY : write_el_cons_on_file
USE control_elastic_constants, ONLY : el_con_geo, lelastic_p, &
                               el_cons_t_available
USE lattices,           ONLY : crystal_parameters
USE control_thermo,     ONLY : ltherm_dos, ltherm_freq
USE control_macro_elasticity, ONLY: macro_el
USE polynomial,         ONLY : poly1, poly2, poly3, poly4, init_poly, &
                               clean_poly
USE data_files,         ONLY : fl_el_cons
USE uniform_pressure,   ONLY : celldm_p, el_cons_p, el_comp_p, b0ec_p
USE control_pressure,   ONLY : npress, press

IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP),    ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: ec_p1(:,:)
TYPE(poly2), ALLOCATABLE :: ec_p2(:,:)
TYPE(poly3), ALLOCATABLE :: ec_p3(:,:)
TYPE(poly4), ALLOCATABLE :: ec_p4(:,:)
INTEGER :: i, j, idata, itemp
INTEGER :: compute_nwork

lelastic_p=.FALSE.
IF (npress==0) RETURN
CALL check_el_cons()
IF (.NOT.el_cons_t_available) RETURN

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
ALLOCATE(b0ec_p(npress))

DO i=1,6
   DO j=i,6
      CALL init_poly(nvar,ec_p1(i,j))
      CALL init_poly(nvar,ec_p2(i,j))
      CALL init_poly(nvar,ec_p3(i,j))
      CALL init_poly(nvar,ec_p4(i,j))
   ENDDO
ENDDO
!
!  Part 1 evaluation of the polynomial coefficients
!
CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)

DO i=1,6
   DO j=i,6
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
!          dependent geometry and computation of the compliances and
!          the bulk modulus
!
CALL interpolate_el_cons_p(celldm_p, nvar, ibrav, ec_p1, ec_p2, ec_p3, &
            ec_p4, poly_degree_elc, el_cons_p, el_comp_p, b0ec_p)

lelastic_p=.TRUE.
filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.el_cons_p'
CALL write_el_cons_on_file(press, npress, ibrav, laue, el_cons_p, b0ec_p, &
                                                    filelastic, 2)

filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.el_comp_p'
CALL write_el_cons_on_file(press, npress, ibrav, laue, el_comp_p, b0ec_p, & 
                                                       filelastic, 3)
DEALLOCATE(x)
DEALLOCATE(f)
DO i=1,6
   DO j=i,6
      CALL clean_poly(ec_p1(i,j))
      CALL clean_poly(ec_p2(i,j))
      CALL clean_poly(ec_p3(i,j))
      CALL clean_poly(ec_p4(i,j))
   ENDDO
ENDDO
DEALLOCATE(b0ec_p)
DEALLOCATE(el_comp_p) 
DEALLOCATE(el_cons_p) 
DEALLOCATE(ec_p1)
DEALLOCATE(ec_p2)
DEALLOCATE(ec_p3)
DEALLOCATE(ec_p4)

RETURN
END SUBROUTINE write_elastic_p
