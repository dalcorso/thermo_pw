!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE write_epsilon_infty_t( )
!----------------------------------------------------------------------
!
!  This routine writes on file the dielectric constant as a function of
!  temperature intepolating it at the geometry that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       epsilon_infty_geo
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE control_epsilon_infty, ONLY : lepsilon_infty, lepsilon_inftyf
USE dielectric_constant, ONLY : write_epsilon_infty_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE control_mur,    ONLY : lmurn
USE anharmonic,     ONLY : celldm_t, epsilon_infty_t
USE ph_freq_anharmonic, ONLY : celldmf_t, epsilon_inftyf_t
USE rap_point_group, ONLY : code_group
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: fileepsilon
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
ALLOCATE(pt_p1(3,3))
ALLOCATE(pt_p2(3,3))
ALLOCATE(pt_p3(3,3))
ALLOCATE(pt_p4(3,3))

DO i=1,3
   DO j=1,3
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
   DO j=1,3
      IF (ABS(epsilon_infty_geo(i,j,1))>1D-7) THEN
         WRITE(stdout,'(/,5x,"Fitting dielectric constants e(",i4,",",i4,&
                              &")")') i, j
         DO idata=1,ndata
            f(idata)=epsilon_infty_geo(i,j,idata)
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
            CALL errore('write_epsilon_infty_t','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the dielectric constant at the temperature 
!          dependent geometry 
!
IF (ltherm_dos) THEN
   CALL interpolate_epsilon_infty(celldm_t, nvar, ibrav, pt_p1, pt_p2, pt_p3, &
               pt_p4, poly_degree_elc, epsilon_infty_t)
   lepsilon_infty=.TRUE.
   fileepsilon='anhar_files/'//TRIM(flanhar)//'.dielec'
   CALL write_epsilon_infty_on_file(temp, ntemp, ibrav, code_group, &
                               epsilon_infty_t, fileepsilon, 0)
ENDIF

IF (ltherm_freq) THEN
   CALL interpolate_epsilon_infty(celldmf_t, nvar, ibrav, pt_p1, pt_p2, pt_p3,&
               pt_p4, poly_degree_elc, epsilon_inftyf_t)
   lepsilon_inftyf=.TRUE.
   fileepsilon='anhar_files/'//TRIM(flanhar)//'.dielec_ph'
   CALL write_epsilon_infty_on_file(temp, ntemp, ibrav, code_group, &
                                    epsilon_inftyf_t, fileepsilon, 0)
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DO i=1,3
   DO j=1,3
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
END SUBROUTINE write_epsilon_infty_t
!
!----------------------------------------------------------------------
SUBROUTINE write_epsilon_infty_pt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the dielectric constant as a function of
!  temperature interpolating them at the structure that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T. This routine
!  does this for selected pressures.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       epsilon_infty_geo
USE rap_point_group, ONLY : code_group
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE dielectric_constant, ONLY : write_epsilon_infty_on_file
USE control_epsilon_infty, ONLY : lepsilon_infty_pt, lepsilon_inftyf_pt
USE lattices,       ONLY : crystal_parameters
USE control_mur,    ONLY : lmurn
USE control_pressure, ONLY : npress_plot, ipress_plot, press
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE anharmonic_pt,     ONLY : celldm_pt, epsilon_infty_pt
USE ph_freq_anharmonic_pt,     ONLY : celldmf_pt, epsilon_inftyf_pt
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp

IMPLICIT NONE
CHARACTER(LEN=256) :: fileepsilon
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
ALLOCATE(pt_p1(3,3))
ALLOCATE(pt_p2(3,3))
ALLOCATE(pt_p3(3,3))
ALLOCATE(pt_p4(3,3))

DO i=1,3
   DO j=1,3
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
   DO j=1,3
      IF (ABS(epsilon_infty_geo(i,j,1))>1D-6) THEN
         WRITE(stdout,'(/,5x,"Fitting dielectric constant e(",i4,",",i4,")")')&
                       i, j
         DO idata=1,ndata
            f(idata)=epsilon_infty_geo(i,j,idata)
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
            CALL errore('write_epsilon_infty_pt','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the dielectric constant at the temperature 
!          dependent geometry
!
IF (ltherm_dos) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_epsilon_infty(celldm_pt(:,:,ipressp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           epsilon_infty_pt(:,:,:,ipressp))
      lepsilon_infty_pt=.TRUE.
      fileepsilon='anhar_files/'//TRIM(flanhar)//'.dielec_press'
      CALL add_value(fileepsilon, press(ipress))
      CALL write_epsilon_infty_on_file(temp, ntemp, ibrav, code_group, &
           epsilon_infty_pt(:,:,:,ipressp), fileepsilon, 0)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_epsilon_infty(celldmf_pt(:,:,ipressp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           epsilon_inftyf_pt(:,:,:,ipressp))
      lepsilon_inftyf_pt=.TRUE.
      fileepsilon='anhar_files/'//TRIM(flanhar)//'.dielec_ph_press'
      CALL add_value(fileepsilon, press(ipress))
      CALL write_epsilon_infty_on_file(temp, ntemp, ibrav, code_group, &
           epsilon_inftyf_pt(:,:,:,ipressp), fileepsilon, 0)
   ENDDO
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DO i=1,3
   DO j=1,3
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
END SUBROUTINE write_epsilon_infty_pt
!
!----------------------------------------------------------------------
SUBROUTINE write_epsilon_infty_ptt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the dielectric constant as a function of
!  pressure intepolating it at the crystal structure that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T. It does
!  this at selected temperatures.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       epsilon_infty_geo

USE rap_point_group, ONLY : code_group
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE dielectric_constant, ONLY : write_epsilon_infty_on_file
USE control_epsilon_infty, ONLY : lepsilon_infty_ptt, lepsilon_inftyf_ptt
USE lattices,       ONLY : crystal_parameters
USE temperature,    ONLY : ntemp_plot, itemp_plot
USE control_pressure,  ONLY : npress, press
USE control_thermo,    ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE anharmonic_ptt,    ONLY : celldm_ptt, epsilon_infty_ptt
USE ph_freq_anharmonic_ptt,  ONLY : celldmf_ptt, epsilon_inftyf_ptt 
USE control_mur, ONLY : lmurn
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: fileepsilon
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
ALLOCATE(pt_p1(3,3))
ALLOCATE(pt_p2(3,3))
ALLOCATE(pt_p3(3,3))
ALLOCATE(pt_p4(3,3))

DO i=1,3
   DO j=1,3
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
   DO j=1,3
      IF (ABS(epsilon_infty_geo(i,j,1))>1.D-6) THEN
         WRITE(stdout,'(/,5x,"Fitting the dielectric constant &
                     & e(",i4,",",i4,")")') i, j
         DO idata=1,ndata
            f(idata)=epsilon_infty_geo(i,j,idata)
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
            CALL errore('write_epsilon_infty_ptt','wrong poly_degree_elc',1)
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
      CALL interpolate_epsilon_infty_p(celldm_ptt(:,:,itempp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           epsilon_infty_ptt(:,:,:,itempp))
      lepsilon_infty_ptt=.TRUE.
      fileepsilon='anhar_files/'//TRIM(flanhar)//'.dielec_temp'
      CALL add_value(fileepsilon, temp(itemp))
      CALL write_epsilon_infty_on_file(press, npress, ibrav, code_group, &
           epsilon_infty_ptt(:,:,:,itempp), fileepsilon, 2)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO itempp=1, ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_epsilon_infty_p(celldmf_ptt(:,:,itempp), nvar, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           epsilon_inftyf_ptt(:,:,:,itempp))
      lepsilon_inftyf_ptt=.TRUE.
      fileepsilon='anhar_files/'//TRIM(flanhar)//'.dielec_ph_temp'
      CALL add_value(fileepsilon, temp(itemp))
      CALL write_epsilon_infty_on_file(press, npress, ibrav, code_group, &
           epsilon_inftyf_ptt(:,:,:,itempp), fileepsilon, 2)
   ENDDO
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)

DO i=1,3
   DO j=1,3
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
END SUBROUTINE write_epsilon_infty_ptt
