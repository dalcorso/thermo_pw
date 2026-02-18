!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE write_epsilonm1_t( )
!----------------------------------------------------------------------
!
!  This routine writes on file the inverse of the static dielectric constant 
!  as a function of temperature intepolating them at the volume that 
!  minimizes the Helmholtz (or Gibbs) free energy at temperature T.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       epsilon_zerom1_geo
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE control_dielectric_constant, ONLY : lepsilon_zerom1, lepsilon_zerom1f
USE dielectric_constant, ONLY : write_epsilon_infty_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE control_mur,    ONLY : lmurn
USE anharmonic,     ONLY : celldm_t, epsilon_zerom1_t
USE ph_freq_anharmonic, ONLY : celldmf_t, epsilon_zerom1f_t
USE rap_point_group, ONLY : code_group
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp

IMPLICIT NONE
CHARACTER(LEN=80) :: astring
CHARACTER(LEN=256) :: filepsilon
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
      IF (ABS(epsilon_zerom1_geo(i,j,1))>1D-5) THEN
         WRITE(stdout,'(/,5x,"Fitting epsilon_0^(-1) (",i4,",",i4,&
                              &")")') i, j
         DO idata=1,ndata
            f(idata)=epsilon_zerom1_geo(i,j,idata)
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
            CALL errore('write_epsilonm1_t','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the inverse of the dielectric constants
!           at the temperature dependent geometry 
!
IF (ltherm_dos) THEN
   CALL interpolate_epsilon_infty(celldm_t, nvar, ibrav, pt_p1, pt_p2, pt_p3, &
               pt_p4, poly_degree_elc, epsilon_zerom1_t)
   lepsilon_zerom1=.TRUE.
   filepsilon='anhar_files/'//TRIM(flanhar)//'.epsilon_zerom1'
   astring="#     epsilon_zero^{-1} (adimensional)"
   CALL write_epsilon_infty_on_file(temp, ntemp, epsilon_zerom1_t, &
                                                 astring, filepsilon, 0)
ENDIF

IF (ltherm_freq) THEN
   CALL interpolate_epsilon_infty(celldmf_t, nvar, ibrav, pt_p1, pt_p2, &
               pt_p3, pt_p4, poly_degree_elc, epsilon_zerom1f_t)
   lepsilon_zerom1f=.TRUE.
   filepsilon='anhar_files/'//TRIM(flanhar)//'.epsilon_zerom1_ph'
   astring="#     epsilon_zero^{-1} (adimensional)"
   CALL write_epsilon_infty_on_file(temp, ntemp, epsilon_zerom1f_t, &
                                                 astring, filepsilon, 0)
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
END SUBROUTINE write_epsilonm1_t

!----------------------------------------------------------------------
SUBROUTINE write_epsilonm1_pt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the inverse of the static dielectric constant 
!  as a function of temperature intepolating them at the volume that 
!  minimizes the Helmholtz (or Gibbs) free energy at temperature T.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       epsilon_zerom1_geo
USE thermo_sym, ONLY : laue
USE control_pressure, ONLY : press, npress_plot, ipress_plot, npress
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE control_dielectric_constant, ONLY : lepsilon_zerom1_pt, lepsilon_zerom1f_pt
USE dielectric_constant, ONLY : write_epsilon_infty_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE control_mur,    ONLY : lmurn
USE anharmonic_pt,     ONLY : celldm_pt, epsilon_zerom1_pt
USE ph_freq_anharmonic_pt, ONLY : celldmf_pt, epsilon_zerom1f_pt
USE rap_point_group, ONLY : code_group
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp

IMPLICIT NONE
CHARACTER(LEN=80) :: astring
CHARACTER(LEN=256) :: filepsilon
INTEGER :: igeo, ibrav, nvar, ndata, ipressp, ipress
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
      IF (ABS(epsilon_zerom1_geo(i,j,1))>1D-5) THEN
         WRITE(stdout,'(/,5x,"Fitting epsilon_0^(-1) (",i4,",",i4,&
                              &")")') i, j
         DO idata=1,ndata
            f(idata)=epsilon_zerom1_geo(i,j,idata)
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
            CALL errore('write_epsilonm1_t','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the inverse of the dielectric constants
!           at the temperature dependent geometry 
!
IF (ltherm_dos) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_epsilon_infty(celldm_pt(1,1,ipressp), nvar, ibrav, &
               pt_p1, pt_p2, pt_p3,  pt_p4, poly_degree_elc, &
               epsilon_zerom1_pt(1,1,1,ipressp))
      lepsilon_zerom1_pt=.TRUE.
      filepsilon='anhar_files/'//TRIM(flanhar)//'.epsilon_zerom1_press'
      CALL add_value(filepsilon, press(ipress))
      astring="#     epsilon_zero^{-1} (adimensional)"
      CALL write_epsilon_infty_on_file(temp, ntemp,                     &
                    epsilon_zerom1_pt(1,1,1,ipressp), astring, filepsilon, 0)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_epsilon_infty(celldmf_pt(1,1,ipressp), nvar, ibrav, &
               pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
               epsilon_zerom1f_pt(1,1,1,ipressp))
      lepsilon_zerom1f_pt=.TRUE.
      filepsilon='anhar_files/'//TRIM(flanhar)//'.epsilon_zerom1_ph_press'
      CALL add_value(filepsilon, press(ipress))
      astring="#     epsilon_zero^{-1} (adimensional)"
      CALL write_epsilon_infty_on_file(temp, ntemp, &
           epsilon_zerom1f_pt(1,1,1,ipressp), astring, filepsilon, 0)
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
END SUBROUTINE write_epsilonm1_pt

!----------------------------------------------------------------------
SUBROUTINE write_epsilonm1_ptt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the inverse of the static dielectric constant 
!  as a function of pressure for several temperatures
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       epsilon_zerom1_geo
USE thermo_sym, ONLY : laue
USE temperature, ONLY : ntemp_plot, itemp_plot, temp
USE control_pressure, ONLY : press, npress
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE control_dielectric_constant, ONLY : lepsilon_zerom1_ptt, &
                                        lepsilon_zerom1f_ptt
USE dielectric_constant, ONLY : write_epsilon_infty_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE control_mur,    ONLY : lmurn
USE anharmonic_ptt, ONLY : celldm_ptt, epsilon_zerom1_ptt
USE ph_freq_anharmonic_ptt, ONLY : celldmf_ptt, epsilon_zerom1f_ptt
USE rap_point_group, ONLY : code_group
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar

IMPLICIT NONE
CHARACTER(LEN=80) :: astring
CHARACTER(LEN=256) :: filepsilon
INTEGER :: igeo, ibrav, nvar, ndata, itempp, ipress
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
      IF (ABS(epsilon_zerom1_geo(i,j,1))>1D-5) THEN
         WRITE(stdout,'(/,5x,"Fitting epsilon_0^(-1) (",i4,",",i4,&
                              &")")') i, j
         DO idata=1,ndata
            f(idata)=epsilon_zerom1_geo(i,j,idata)
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
            CALL errore('write_epsilonm1_t','wrong poly_degree_elc',1)
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
!  Part 2: interpolation of the inverse of the dielectric constants
!           at the temperature dependent geometry 
!
IF (ltherm_dos) THEN
   DO itempp=1, ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_epsilon_infty_p(celldm_ptt(1,1,itempp), nvar, ibrav, &
               pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc,             &
               epsilon_zerom1_ptt(1,1,1,itempp))
      lepsilon_zerom1_ptt=.TRUE.
      filepsilon='anhar_files/'//TRIM(flanhar)//'.epsilon_zerom1_temp'
      CALL add_value(filepsilon, temp(itemp))
      astring="#     epsilon_zero^{-1} (adimensional)"
      CALL write_epsilon_infty_on_file(press, npress,     &
                epsilon_zerom1_ptt(1,1,1,itempp), astring, filepsilon, 2)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO itempp=1, ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_epsilon_infty_p(celldmf_ptt(1,1,itempp), nvar, ibrav, &
               pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc,              &
               epsilon_zerom1f_ptt(1,1,1,itempp))
      lepsilon_zerom1f_ptt=.TRUE.
      filepsilon='anhar_files/'//TRIM(flanhar)//'.epsilon_zerom1_ph_temp'
      CALL add_value(filepsilon, temp(itemp))
      astring="#     epsilon_zero^{-1} (adimensional)"
      CALL write_epsilon_infty_on_file(press, npress, &
                epsilon_zerom1f_ptt(1,1,1,itempp), astring, filepsilon, 2)
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
END SUBROUTINE write_epsilonm1_ptt
