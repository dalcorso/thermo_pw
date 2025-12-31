!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE write_zeu_t( )
!----------------------------------------------------------------------
!
!  This routine writes on file the Born effective charges as a function of
!  temperature intepolating them at the geometry that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE ions_base,  ONLY : nat
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       zeu_geo
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE control_epsilon_infty, ONLY : lzeu, lzeuf
USE dielectric_constant, ONLY : write_zeu_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE control_mur,    ONLY : lmurn
USE anharmonic,     ONLY : celldm_t, zeu_t
USE ph_freq_anharmonic, ONLY : celldmf_t, zeuf_t
USE rap_point_group, ONLY : code_group
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filezeu
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: pt_p1(:,:,:)
TYPE(poly2), ALLOCATABLE :: pt_p2(:,:,:)
TYPE(poly3), ALLOCATABLE :: pt_p3(:,:,:)
TYPE(poly4), ALLOCATABLE :: pt_p4(:,:,:)
INTEGER :: i, j, na, idata, itemp
INTEGER :: compute_nwork

ibrav=ibrav_geo(1)
nvar=crystal_parameters(ibrav)
IF (lmurn) nvar=1

ndata=compute_nwork()

ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(pt_p1(3,3,nat))
ALLOCATE(pt_p2(3,3,nat))
ALLOCATE(pt_p3(3,3,nat))
ALLOCATE(pt_p4(3,3,nat))

DO na=1, nat
   DO i=1,3
      DO j=1,3
         CALL init_poly(nvar,pt_p1(i,j,na))
         CALL init_poly(nvar,pt_p2(i,j,na))
         CALL init_poly(nvar,pt_p3(i,j,na))
         CALL init_poly(nvar,pt_p4(i,j,na))
      ENDDO
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

DO na=1, nat
   DO i=1,3
      DO j=1,3
         IF (ABS(zeu_geo(i,j,na,1))>1D-7) THEN
            WRITE(stdout,'(/,5x,"Fitting Z^*(",i4,",",i4,&
                              &",i4)")') i, j, na
            DO idata=1,ndata
               f(idata)=zeu_geo(i,j,na,idata)
            END DO

            IF (poly_degree_elc==4) THEN
               CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,pt_p4(i,j,na)) 
            ELSEIF (poly_degree_elc==3) THEN
               CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,pt_p3(i,j,na)) 
            ELSEIF (poly_degree_elc==2) THEN
               CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,pt_p2(i,j,na))      
            ELSEIF (poly_degree_elc==1) THEN
               CALL fit_multi_linear(ndata,nvar,lsolve,x,f,pt_p1(i,j,na))
            ELSE
               CALL errore('write_zeu_t','wrong poly_degree_elc',1)
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!  Part 2: interpolation of the dielectric constant at the temperature 
!          dependent geometry 
!
IF (ltherm_dos) THEN
   CALL interpolate_zeu(celldm_t, nvar, nat, ibrav, pt_p1, pt_p2, pt_p3, &
               pt_p4, poly_degree_elc, zeu_t)
   lzeu=.TRUE.
   filezeu='anhar_files/'//TRIM(flanhar)//'.zeu'
   CALL write_zeu_on_file(temp, ntemp, nat, ibrav, code_group, zeu_t, filezeu, 0)
ENDIF

IF (ltherm_freq) THEN
   CALL interpolate_zeu(celldmf_t, nvar, nat, ibrav, pt_p1, pt_p2, pt_p3,&
                        pt_p4, poly_degree_elc, zeuf_t)
   lzeuf=.TRUE.
   filezeu='anhar_files/'//TRIM(flanhar)//'.zeu_ph'
   CALL write_zeu_on_file(temp, ntemp, nat, ibrav, code_group, &
                                                   zeuf_t, filezeu, 0)
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DO na=1, nat
   DO i=1,3
      DO j=1,3
         CALL clean_poly(pt_p1(i,j,na))
         CALL clean_poly(pt_p2(i,j,na))
         CALL clean_poly(pt_p3(i,j,na))
         CALL clean_poly(pt_p4(i,j,na))
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(pt_p1)
DEALLOCATE(pt_p2)
DEALLOCATE(pt_p3)
DEALLOCATE(pt_p4)

RETURN
END SUBROUTINE write_zeu_t
!
!----------------------------------------------------------------------
SUBROUTINE write_zeu_pt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the Born effective charges as a function of
!  temperature interpolating them at the structure that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T. This routine
!  does this for selected pressures.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE ions_base,  ONLY : nat
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       zeu_geo
USE rap_point_group, ONLY : code_group
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE dielectric_constant, ONLY : write_zeu_on_file
USE control_epsilon_infty, ONLY : lzeu_pt, lzeuf_pt
USE lattices,       ONLY : crystal_parameters
USE control_mur,    ONLY : lmurn
USE control_pressure, ONLY : npress_plot, ipress_plot, press
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE anharmonic_pt,     ONLY : celldm_pt, zeu_pt
USE ph_freq_anharmonic_pt,     ONLY : celldmf_pt, zeuf_pt
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp

IMPLICIT NONE
CHARACTER(LEN=256) :: filezeu
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: pt_p1(:,:,:)
TYPE(poly2), ALLOCATABLE :: pt_p2(:,:,:)
TYPE(poly3), ALLOCATABLE :: pt_p3(:,:,:)
TYPE(poly4), ALLOCATABLE :: pt_p4(:,:,:)
INTEGER :: i, j, na, idata, itemp, ipressp, ipress
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
ALLOCATE(pt_p1(3,3,nat))
ALLOCATE(pt_p2(3,3,nat))
ALLOCATE(pt_p3(3,3,nat))
ALLOCATE(pt_p4(3,3,nat))

DO na=1, nat
   DO i=1,3
      DO j=1,3
         CALL init_poly(nvar,pt_p1(i,j,na))
         CALL init_poly(nvar,pt_p2(i,j,na))
         CALL init_poly(nvar,pt_p3(i,j,na))
         CALL init_poly(nvar,pt_p4(i,j,na))
      ENDDO
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

DO na=1,nat
   DO i=1,3
      DO j=1,3
         IF (ABS(zeu_geo(i,j,na,1))>1D-6) THEN
            WRITE(stdout,'(/,5x,"Fitting dielectric constant &
                &e(",i4,",",i4,",",i4,")")') i, j, na
            DO idata=1,ndata
               f(idata)=zeu_geo(i,j,na,idata)
            END DO

            IF (poly_degree_elc==4) THEN
               CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,pt_p4(i,j,na)) 
            ELSEIF (poly_degree_elc==3) THEN
               CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,pt_p3(i,j,na)) 
            ELSEIF (poly_degree_elc==2) THEN
               CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,pt_p2(i,j,na))   
            ELSEIF (poly_degree_elc==1) THEN
               CALL fit_multi_linear(ndata,nvar,lsolve,x,f,pt_p1(i,j,na))
            ELSE
               CALL errore('write_epsilon_infty_pt','wrong poly_degree_elc',1)
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!  Part 2: interpolation of the dielectric constant at the temperature 
!          dependent geometry
!
IF (ltherm_dos) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_zeu(celldm_pt(:,:,ipressp), nvar, nat, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           zeu_pt(:,:,:,:,ipressp))
      lzeu_pt=.TRUE.
      filezeu='anhar_files/'//TRIM(flanhar)//'.zeu_press'
      CALL add_value(filezeu, press(ipress))
      CALL write_zeu_on_file(temp, ntemp, nat, ibrav, code_group, &
           zeu_pt(:,:,:,:,ipressp), filezeu, 0)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO ipressp=1, npress_plot
      ipress=ipress_plot(ipressp)
      CALL interpolate_zeu(celldmf_pt(:,:,ipressp), nvar, nat, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           zeuf_pt(:,:,:,:,ipressp))
      lzeuf_pt=.TRUE.
      filezeu='anhar_files/'//TRIM(flanhar)//'.zeu_ph_press'
      CALL add_value(filezeu, press(ipress))
      CALL write_zeu_on_file(temp, ntemp, nat, ibrav, code_group, &
           zeuf_pt(:,:,:,:,ipressp), filezeu, 0)
   ENDDO
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)

DO na=1, nat
   DO i=1,3
      DO j=1,3
         CALL clean_poly(pt_p1(i,j,na))
         CALL clean_poly(pt_p2(i,j,na))
         CALL clean_poly(pt_p3(i,j,na))
         CALL clean_poly(pt_p4(i,j,na))
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(pt_p1)
DEALLOCATE(pt_p2)
DEALLOCATE(pt_p3)
DEALLOCATE(pt_p4)

RETURN
END SUBROUTINE write_zeu_pt
!
!----------------------------------------------------------------------
SUBROUTINE write_zeu_ptt( )
!----------------------------------------------------------------------
!
!  This routine writes on file the Born effective charges as a function of
!  pressure intepolating it at the crystal structure that minimizes the
!  Helmholtz (or Gibbs) free energy at temperature T. It does
!  this at selected temperatures.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE ions_base,  ONLY : nat
USE thermo_mod, ONLY : ibrav_geo, celldm_geo, omega_geo, celldm_geo_eos, &
                       zeu_geo

USE rap_point_group, ONLY : code_group
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc

USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic

USE dielectric_constant, ONLY : write_zeu_on_file
USE control_epsilon_infty, ONLY : lzeu_ptt, lzeuf_ptt
USE lattices,       ONLY : crystal_parameters
USE temperature,    ONLY : ntemp_plot, itemp_plot
USE control_pressure,  ONLY : npress, press
USE control_thermo,    ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE anharmonic_ptt,    ONLY : celldm_ptt, zeu_ptt
USE ph_freq_anharmonic_ptt,  ONLY : celldmf_ptt, zeuf_ptt 
USE control_mur, ONLY : lmurn
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filezeu
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:)
TYPE(poly1), ALLOCATABLE :: pt_p1(:,:,:)
TYPE(poly2), ALLOCATABLE :: pt_p2(:,:,:)
TYPE(poly3), ALLOCATABLE :: pt_p3(:,:,:)
TYPE(poly4), ALLOCATABLE :: pt_p4(:,:,:)
INTEGER :: i, j, na, idata, itemp, itempp
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
ALLOCATE(pt_p1(3,3,nat))
ALLOCATE(pt_p2(3,3,nat))
ALLOCATE(pt_p3(3,3,nat))
ALLOCATE(pt_p4(3,3,nat))

DO na=1,nat
   DO i=1,3
      DO j=1,3
         CALL init_poly(nvar,pt_p1(i,j,na))
         CALL init_poly(nvar,pt_p2(i,j,na))
         CALL init_poly(nvar,pt_p3(i,j,na))
         CALL init_poly(nvar,pt_p4(i,j,na))
      ENDDO
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

DO na=1,nat
   DO i=1,3
      DO j=1,3
         IF (ABS(zeu_geo(i,j,na,1))>1.D-6) THEN
            WRITE(stdout,'(/,5x,"Fitting the dielectric constant &
                     & e(",i4,",",i4,",",i4,")")') i, j, na
            DO idata=1,ndata
               f(idata)=zeu_geo(i,j,na,idata)
            END DO

            IF (poly_degree_elc==4) THEN
               CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,pt_p4(i,j,na)) 
            ELSEIF (poly_degree_elc==3) THEN
               CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,pt_p3(i,j,na)) 
            ELSEIF (poly_degree_elc==2) THEN
               CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,pt_p2(i,j,na))   
            ELSEIF (poly_degree_elc==1) THEN
               CALL fit_multi_linear(ndata,nvar,lsolve,x,f,pt_p1(i,j,na))
            ELSE
               CALL errore('write_zeu_ptt','wrong poly_degree_elc',1)
            ENDIF
         ENDIF
      ENDDO
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
      CALL interpolate_zeu_p(celldm_ptt(:,:,itempp), nvar, nat, ibrav, &
         pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, zeu_ptt(:,:,:,:,itempp))
      lzeu_ptt=.TRUE.
      filezeu='anhar_files/'//TRIM(flanhar)//'.zeu_temp'
      CALL add_value(filezeu, temp(itemp))
      CALL write_zeu_on_file(press, npress, nat, ibrav, code_group, &
           zeu_ptt(:,:,:,:,itempp), filezeu, 2)
   ENDDO
ENDIF

IF (ltherm_freq) THEN
   DO itempp=1, ntemp_plot
      itemp=itemp_plot(itempp)
      CALL interpolate_zeu_p(celldmf_ptt(:,:,itempp), nvar, nat, ibrav, &
           pt_p1, pt_p2, pt_p3, pt_p4, poly_degree_elc, &
           zeuf_ptt(:,:,:,:,itempp))
      lzeuf_ptt=.TRUE.
      filezeu='anhar_files/'//TRIM(flanhar)//'.zeu_ph_temp'
      CALL add_value(filezeu, temp(itemp))
      CALL write_zeu_on_file(press, npress, nat, ibrav, code_group, &
                        zeuf_ptt(:,:,:,:,itempp), filezeu, 2)
   ENDDO
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)

DO na=1,nat
   DO i=1,3
      DO j=1,3
         CALL clean_poly(pt_p1(i,j,na))
         CALL clean_poly(pt_p2(i,j,na))
         CALL clean_poly(pt_p3(i,j,na))
         CALL clean_poly(pt_p4(i,j,na))
      ENDDO
   ENDDO
ENDDO

DEALLOCATE(pt_p1)
DEALLOCATE(pt_p2)
DEALLOCATE(pt_p3)
DEALLOCATE(pt_p4)

RETURN
END SUBROUTINE write_zeu_ptt
