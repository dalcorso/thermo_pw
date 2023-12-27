!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit()
  !-----------------------------------------------------------------------
  !
  !   This routine receives the total energy for several values of 
  !   celldm and fits them with a quadratic polynomial of dimension 
  !   equal to the number of indipendent parameters in celldm. 
  !   It finds the minimum of the quadratic function.
  !   If lquartic is true it fits the data with a quartic polynomial 
  !   and starting from the quadratic minimum it finds the minimum of
  !   the quartic polynomial.
  !   
  !   the output of this routine is celldm0(:) and emin the energy at
  !   the mimimum.
  !
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : ibrav
  USE control_mur,  ONLY : emin, vmin
  USE equilibrium_conf, ONLY : celldm0
  USE initial_conf, ONLY : ibrav_save
  USE thermo_mod,   ONLY : celldm_geo, omega_geo, energy_geo
  USE control_pressure, ONLY : pressure, pressure_kb
  USE control_quadratic_energy, ONLY : hessian_v, hessian_e, x_pos_min, &
                                       p2, nvar
  USE control_quartic_energy, ONLY : p4, x_min_4, lquartic, lsolve,     &
                                     hessian4_v, hessian4_e
  USE lattices,     ONLY : expand_celldm, crystal_parameters
  USE quadratic_surfaces, ONLY : fit_multi_quadratic,  &
                          find_quadratic_extremum,     &
                          write_quadratic_hessian,     &
                          print_quadratic_polynomial,  &
                          summarize_fitting_data,      &
                          compare_quadratic_fit,       &
                          introduce_quadratic_fit, print_chisq_quadratic
  USE quartic_surfaces, ONLY : fit_multi_quartic, &
                          find_quartic_extremum, print_quartic_polynomial, &
                          print_chisq_quartic, introduce_quartic_fit,      &
                          write_quartic_hessian
  USE polynomial,   ONLY : init_poly
  USE vector_mod,       ONLY : write_vector
  USE io_global,    ONLY : stdout
  IMPLICIT NONE

  INTEGER  :: tncoeff4, ndata, ncoeff
  REAL(DP), ALLOCATABLE :: x(:,:), y(:), f(:)
  REAL(DP) :: ymin, ymin4
  INTEGER  :: compute_nwork
  REAL(DP) :: compute_omega_geo
  !
  ! Only the first image does the calculation
  !
  celldm0(:)=0.0_DP
  emin=0.0_DP
  !
  WRITE(stdout,'(/,5x,71("-"))')
  !
  nvar=crystal_parameters(ibrav)
  !
  ndata = compute_nwork()

  IF (pressure>0.0_DP) THEN
     WRITE(stdout,'(/,5x,"Fitting the enthalpy with a function of the crystal &
                                         &parameters.")') 
     WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
  ELSE
     WRITE(stdout,'(/,5x,"Fitting the energy with a function of the crystal &
                                         &parameters.")') 
  ENDIF

  ALLOCATE(x(nvar,ndata))
  ALLOCATE(x_pos_min(nvar))
  ALLOCATE(hessian_e(nvar))
  ALLOCATE(hessian_v(nvar, nvar))
  ALLOCATE(f(ndata))

  CALL init_poly(nvar,p2)
  ncoeff=1+nvar+p2%ncoeff2
  CALL introduce_quadratic_fit(nvar, ncoeff, ndata)

  f(:)=energy_geo(:) + pressure * omega_geo(:)

  CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)
  !
!  CALL summarize_fitting_data(nvar, ndata, x, f)

  CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f,p2)
 
  CALL print_quadratic_polynomial(nvar, p2)

  IF (pressure>0.0_DP) THEN
     WRITE(stdout,'(7x,"Enthalpy (1)      Fitted enthalpy (2)   &
                                                    &DeltaE (1)-(2)")') 
  ELSE
     WRITE(stdout,'(7x,"Energy (1)      Fitted energy (2)   &
                                                    &DeltaE (1)-(2)")') 
  ENDIF
  CALL compare_quadratic_fit(ndata, nvar, x, f, p2)
  CALL print_chisq_quadratic(ndata, nvar, x, f, p2)

  CALL find_quadratic_extremum(nvar,x_pos_min,ymin,p2)
  !
  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')

  CALL write_vector(nvar,x_pos_min)
  CALL print_energy(ymin)
  !
  CALL write_quadratic_hessian(nvar,p2,hessian_v,hessian_e)
  !
  CALL expand_celldm(celldm0, x_pos_min, nvar, ibrav)
  emin=ymin
  !
  IF (lquartic) THEN
     ALLOCATE(x_min_4(nvar))
     ALLOCATE(hessian4_e(nvar))
     ALLOCATE(hessian4_v(nvar, nvar))
     CALL init_poly(nvar,p4)
     tncoeff4=1+nvar+p4%ncoeff2+p4%ncoeff3+p4%ncoeff4
     CALL introduce_quartic_fit(nvar,tncoeff4,ndata)

     CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,p4)
     !
     CALL print_quartic_polynomial(nvar,p4) 
!    WRITE(stdout,'(/,7x,"Energy (1)    Fitted energy (2)   DeltaE (1)-(2)")') 
     CALL print_chisq_quartic(ndata, nvar, x, f, p4)
!
!   searching the minimum starting from the minimum of the quadratic
!
     x_min_4=x_pos_min
     CALL find_quartic_extremum(nvar,x_min_4,ymin4,p4)

     WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     CALL write_vector(nvar,x_min_4)
     CALL print_energy(ymin4)
     CALL write_quartic_hessian(nvar, x_min_4, p4, hessian4_v, hessian4_e)
!
     CALL expand_celldm(celldm0, x_min_4, nvar, ibrav)
     emin=ymin4
  ENDIF
  vmin=compute_omega_geo(ibrav_save, celldm0)

  DEALLOCATE(f)
  DEALLOCATE(x)
  !
  RETURN
END SUBROUTINE quadratic_fit

!-------------------------------------------------------------------
SUBROUTINE quadratic_fit_t_run()
!-------------------------------------------------------------------
USE kinds, ONLY : DP
USE temperature, ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE anharmonic, ONLY : celldm_t, free_e_min_t, p1t_t, p2t_t, p3t_t, p4t_t
USE control_eldos, ONLY : lel_free_energy
USE io_global, ONLY : stdout
USE mp_world, ONLY : world_comm
USE mp,       ONLY : mp_sum

IMPLICIT NONE
INTEGER :: itempp, itemp, startt, lastt

CALL divide(world_comm, ntemp, startt, lastt)

WRITE(stdout,'(/,5x,70("-"))')
IF (pressure_kb > 0.0_DP) THEN
   WRITE(stdout,'(5x, "Gibbs energy from phdos, at ", i5, " Temperatures")')&
                       ntemp 
   WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
ELSE
   WRITE(stdout,'(5x, "Helmholtz free energy from phdos, at ",i5, &
                      &" Temperatures" )') ntemp
ENDIF
IF (lel_free_energy) WRITE(stdout,'(5x,"Adding electron free energy")') 

celldm_t=0.0_DP
free_e_min_t=0.0_DP
DO itemp = startt, lastt
   CALL quadratic_fit_t(itemp, celldm_t(:,itemp), free_e_min_t(itemp), &
                        p1t_t(itemp), p2t_t(itemp), p3t_t(itemp),      &
                        p4t_t(itemp))
ENDDO
CALL mp_sum(celldm_t, world_comm)
CALL mp_sum(free_e_min_t, world_comm)

RETURN
END SUBROUTINE quadratic_fit_t_run

!-------------------------------------------------------------------
SUBROUTINE quadratic_fit_t_noe_run()
!-------------------------------------------------------------------
USE kinds, ONLY : DP
USE temperature, ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE anharmonic, ONLY : celldm_noe_t, free_e_min_noe_t, p1t_noe_t, &
                       p2t_noe_t, p3t_noe_t, p4t_noe_t
USE control_eldos, ONLY : lel_free_energy
USE io_global, ONLY : stdout
USE mp_world, ONLY : world_comm
USE mp,       ONLY : mp_sum

IMPLICIT NONE
INTEGER :: itempp, itemp, idata, startt, lastt

IF (.NOT.lel_free_energy) RETURN

CALL divide(world_comm, ntemp, startt, lastt)

WRITE(stdout,'(/,5x,70("-"))')
IF (pressure_kb > 0.0_DP) THEN
   WRITE(stdout,'(5x, "Gibbs energy from phdos, at ", i5, " Temperatures")')&
                       ntemp 
   WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
ELSE
   WRITE(stdout,'(5x, "Helmholtz free energy from phdos, at ",i5, &
                      &" Temperatures" )') ntemp
ENDIF
WRITE(stdout,'(5x,"Electron free energy not added")') 

celldm_noe_t=0.0_DP
free_e_min_noe_t=0.0_DP
DO itemp = startt, lastt
   CALL quadratic_fit_t(itemp, celldm_noe_t(:,itemp), free_e_min_noe_t(itemp),&
                        p1t_noe_t(itemp), p2t_noe_t(itemp), p3t_noe_t(itemp), &
                        p4t_noe_t(itemp))
ENDDO
CALL mp_sum(celldm_noe_t, world_comm)
CALL mp_sum(free_e_min_noe_t, world_comm)

RETURN
END SUBROUTINE quadratic_fit_t_noe_run
!
!-------------------------------------------------------------------
SUBROUTINE quadratic_fitf_t_run()
!-------------------------------------------------------------------
USE kinds, ONLY : DP
USE temperature, ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE ph_freq_anharmonic, ONLY : celldmf_t, free_e_minf_t, p1tf_t, &
                        p2tf_t, p3tf_t, p4tf_t
USE control_eldos, ONLY : lel_free_energy
USE io_global, ONLY : stdout
USE mp_world, ONLY : world_comm
USE mp,       ONLY : mp_sum

IMPLICIT NONE
INTEGER :: itempp, itemp, startt, lastt

CALL divide(world_comm, ntemp, startt, lastt)

WRITE(stdout,'(/,5x,70("-"))')
IF (pressure_kb > 0.0_DP) THEN
   WRITE(stdout,'(5x, "Gibbs energy from BZ int., at ", i5, " Temperatures")')&
                       ntemp 
   WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
ELSE
   WRITE(stdout,'(5x, "Helmholtz free energy from BZ int., at ",i5, &
                      & " Temperatures" )') ntemp
ENDIF
IF (lel_free_energy) WRITE(stdout,'(5x,"Adding electron free energy")') 

celldmf_t=0.0_DP
free_e_minf_t=0.0_DP
DO itemp = startt, lastt
   CALL quadratic_fit_t(itemp, celldmf_t(:,itemp), free_e_minf_t(itemp), &
                        p1tf_t(itemp), p2tf_t(itemp), p3tf_t(itemp),      &
                        p4tf_t(itemp))
ENDDO
CALL mp_sum(celldmf_t, world_comm)
CALL mp_sum(free_e_minf_t, world_comm)

RETURN
END SUBROUTINE quadratic_fitf_t_run

!-------------------------------------------------------------------
SUBROUTINE quadratic_fitf_t_noe_run()
!-------------------------------------------------------------------
USE kinds, ONLY : DP
USE temperature, ONLY : ntemp, temp
USE control_pressure, ONLY : pressure_kb
USE ph_freq_anharmonic, ONLY : celldmf_noe_t, free_e_minf_noe_t, p1tf_noe_t, &
                       p2tf_noe_t, p3tf_noe_t, p4tf_noe_t
USE control_eldos, ONLY : lel_free_energy
USE io_global, ONLY : stdout
USE mp_world, ONLY : world_comm
USE mp,       ONLY : mp_sum

IMPLICIT NONE
INTEGER :: itempp, itemp, idata, startt, lastt

IF (.NOT.lel_free_energy) RETURN

CALL divide(world_comm, ntemp, startt, lastt)

WRITE(stdout,'(/,5x,70("-"))')
IF (pressure_kb > 0.0_DP) THEN
   WRITE(stdout,'(5x, "Gibbs energy from BZ int., at ", i5, " Temperatures")')&
                       ntemp 
   WRITE(stdout,'(5x, "Pressure is :",f12.6)') pressure_kb
ELSE
   WRITE(stdout,'(5x, "Helmholtz free energy from BZ int., at ",i5, &
                      &" Temperatures" )') ntemp
ENDIF
WRITE(stdout,'(5x,"Electron free energy not added")') 

celldmf_noe_t=0.0_DP
free_e_minf_noe_t=0.0_DP
DO itemp = startt, lastt
   CALL quadratic_fit_t(itemp, celldmf_noe_t(:,itemp),      &
                        free_e_minf_noe_t(itemp),           &
                        p1tf_noe_t(itemp), p2tf_noe_t(itemp), &
                        p3tf_noe_t(itemp), p4tf_noe_t(itemp))
ENDDO
CALL mp_sum(celldmf_noe_t, world_comm)
CALL mp_sum(free_e_minf_noe_t, world_comm)

RETURN
END SUBROUTINE quadratic_fitf_t_noe_run

!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_t(itemp, celldm_t, free_e_min_t, pt1, pt2, pt3, pt4 )
  !-----------------------------------------------------------------------
  !
  !   This routine receives the polynomial interpolation of the 
  !   free energy at temperature itemp
  !   Then it adds the quadratic polynomial that fits the energy
  !   (or enthalphy) and finds the minimum.
  !   If lquartic is true it interpolates the free energy with a polynomial 
  !   of degree poly_degree_ph, adds the quartic polynomial that fits 
  !   the energy (or enthalpy) and finds its minimum starting from
  !   the minimum of the quadratic polynomial.
  !
  !   The output of this routine is celldm_t at the given temperature itemp
  !   and the free energy at the minimum free_e_min_t
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE thermo_mod,  ONLY : celldm_geo, no_ph
  USE control_quadratic_energy, ONLY : nvar, enthalpy_p2 => p2
  USE control_quartic_energy, ONLY :  p4, lquartic, poly_degree_ph, lsolve
  USE lattices,    ONLY : compress_celldm, expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : fit_multi_quadratic, &
                      find_two_quadratic_extremum,    &
                      print_quadratic_polynomial,     &
                      summarize_fitting_data,         &
                      introduce_quadratic_fit, print_chisq_quadratic

  USE cubic_surfaces, ONLY : fit_multi_cubic, print_chisq_cubic
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      fit_multi_quartic, find_two_quartic_extremum,     &
                      find_quartic_cubic_extremum, print_chisq_quartic, &
                      find_quartic_linear_extremum

  USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  INTEGER  :: itemp
  REAL(DP) :: celldm_t(6), free_e_min_t
  REAL(DP), ALLOCATABLE ::  x_pos_min(:)
  TYPE(poly1) :: pt1
  TYPE(poly2) :: pt2
  TYPE(poly3) :: pt3             
  TYPE(poly4) :: pt4             
  REAL(DP) :: ymin
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, enthalpy_p2, pt2)
!  WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!  CALL write_vector(nvar,x_pos_min)
!  CALL print_genergy(ymin)

  IF (lquartic) THEN
     IF (poly_degree_ph==4) THEN
!        WRITE(stdout,'(/,5x, "Fit improved with a fourth order polynomial")') 
!        CALL init_poly(nvar,pt4)
!        CALL fit_multi_quartic(ndata, nvar, lsolve, x, f, pt4)
!        CALL print_chisq_quartic(ndata, nvar, x, f, pt4)
        CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, p4, pt4) 
!        CALL clean_poly(pt4)
!        WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
     ELSEIF (poly_degree_ph==3) THEN
!        WRITE(stdout,'(/,5x, "Fit improved with a third order polynomial")') 
!        CALL init_poly(nvar,pt3)
!        CALL fit_multi_cubic(ndata, nvar, lsolve, x, f, pt3)
!        CALL print_chisq_cubic(ndata, nvar, x, f, pt3)
        CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, p4, pt3)
!        CALL clean_poly(pt3)
!        WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
     ELSEIF (poly_degree_ph==1) THEN
!        WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
!        CALL init_poly(nvar,pt1)
!        CALL fit_multi_linear(ndata, nvar, lsolve, x, f, pt1)
!        CALL print_chisq_linear(ndata, nvar, x, f, pt1)
        CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, p4, pt1)
!        CALL clean_poly(pt1)
!        WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
     ELSE
!        WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
        CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, p4, pt2)
     ENDIF
!     CALL write_vector(nvar,x_pos_min)
!     CALL print_genergy(ymin)
  ENDIF

  free_e_min_t=ymin
  CALL expand_celldm(celldm_t, x_pos_min, nvar, ibrav)

  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_t_pm()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at pressure plus and minus dp.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine are celldm_t_p1  and celldm_t_m1
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE temperature, ONLY : ntemp
  USE anharmonic,  ONLY : p1t_t, p2t_t, p3t_t, p4t_t, celldm_t_p1, celldm_t_m1
  USE uniform_pressure, ONLY : p2_p_p1, p2_p_m1, p4_p_p1, p4_p_m1

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                            p2_p_p1, p2t_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_p1, p4t_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_p1, p3t_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_p1, p1t_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_p1, p2t_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldm_t_p1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p_m1, p2t_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_m1, p4t_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_m1, p3t_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_m1, p1t_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_m1, p2t_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldm_t_m1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_t_pm
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fitf_t_pm()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at pressure plus and minus dp.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine are celldmf_t_p1  and celldmf_t_m1
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE temperature, ONLY : ntemp
  USE ph_freq_anharmonic,  ONLY : p1tf_t, p2tf_t, p3tf_t, p4tf_t, celldmf_t_p1, celldmf_t_m1
  USE uniform_pressure, ONLY : p2_p_p1, p2_p_m1, p4_p_p1, p4_p_m1

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                            p2_p_p1, p2tf_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_p1, p4tf_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_p1, p3tf_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_p1, p1tf_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_p1, p2tf_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldmf_t_p1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p_m1, p2tf_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_m1, p4tf_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_m1, p3tf_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_m1, p1tf_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_m1, p2tf_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldmf_t_m1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fitf_t_pm
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_noe_t_pm()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at pressure plus and minus dp.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine are celldm_noe_t_p1  and celldm_noe_t_m1
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE temperature, ONLY : ntemp
  USE anharmonic,  ONLY : p1t_noe_t, p2t_noe_t, p3t_noe_t, p4t_noe_t, &
                          celldm_noe_t_p1, celldm_noe_t_m1
  USE uniform_pressure, ONLY : p2_p_p1, p2_p_m1, p4_p_p1, p4_p_m1

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE control_eldos,  ONLY : lel_free_energy

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  IF (.NOT.lel_free_energy) RETURN

  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                            p2_p_p1, p2t_noe_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_p1, p4t_noe_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_p1, p3t_noe_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_p1, p1t_noe_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_p1, p2t_noe_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldm_noe_t_p1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p_m1, p2t_noe_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_m1, p4t_noe_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_m1, p3t_noe_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_m1, p1t_noe_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_m1, p2t_noe_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldm_noe_t_m1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_noe_t_pm
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fitf_noe_t_pm()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at pressure plus and minus dp.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine are celldmf_noe_t_p1  and celldmf_noe_t_m1
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE temperature, ONLY : ntemp
  USE ph_freq_anharmonic,  ONLY : p1tf_noe_t, p2tf_noe_t, p3tf_noe_t, &
                       p4tf_noe_t, celldmf_noe_t_p1, celldmf_noe_t_m1
  USE uniform_pressure, ONLY : p2_p_p1, p2_p_m1, p4_p_p1, p4_p_m1

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE control_eldos,  ONLY : lel_free_energy

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  IF (.NOT.lel_free_energy) RETURN
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                            p2_p_p1, p2tf_noe_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_p1, p4tf_noe_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_p1, p3tf_noe_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_p1, p1tf_noe_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_p1, p2tf_noe_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldmf_noe_t_p1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO

  DO itemp=1, ntemp
     CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p_m1, p2tf_noe_t(itemp))
     IF (lquartic) THEN
        IF (poly_degree_ph==4) THEN
           CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p_m1, p4tf_noe_t(itemp)) 
        ELSEIF (poly_degree_ph==3) THEN
           CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p_m1, p3tf_noe_t(itemp))
        ELSEIF (poly_degree_ph==1) THEN
           CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p_m1, p1tf_noe_t(itemp))
        ELSE
           CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p_m1, p2tf_noe_t(itemp))
        ENDIF
     ENDIF
     CALL expand_celldm(celldmf_noe_t_m1(:,itemp), x_pos_min, nvar, ibrav)
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fitf_noe_t_pm
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_pt()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at all temperatures for the set of pressures specified by press_plot.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine is celldm_pt and the Gibbs energy at 
  !   the minimum emin_pt.
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE control_pressure, ONLY : npress_plot, ipress_plot
  USE temperature, ONLY : ntemp, ntemp_plot, temp_plot
  USE anharmonic,  ONLY : p1t_t, p2t_t, p3t_t, p4t_t
  USE anharmonic_pt, ONLY : celldm_pt, emin_pt
  USE uniform_pressure, ONLY : p2_p, p4_p

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp, ipress, ipressp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO ipressp=1, npress_plot
     ipress=ipress_plot(ipressp)
     DO itemp=1, ntemp

        IF (itemp==1.OR..NOT.lquartic) THEN
           CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p(ipress), p2t_t(itemp))
        ENDIF
!        WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!        CALL write_vector(nvar,x_pos_min)
!        CALL print_genergy(ymin)

        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!           WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                          p4_p(ipress), p4t_t(itemp)) 
!           WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!           WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                      &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                    p4_p(ipress), p3t_t(itemp))
!           WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!           WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                   p4_p(ipress), p1t_t(itemp))
!           WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!           WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
             CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p(ipress), p2t_t(itemp))
           ENDIF
!          CALL write_vector(nvar,x_pos_min)
!          CALL print_genergy(ymin)
        ENDIF

        emin_pt(itemp,ipressp)=ymin
        CALL expand_celldm(celldm_pt(:,itemp,ipressp), x_pos_min, nvar, ibrav)

     ENDDO
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_pt
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fitf_pt()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at all temperatures for the set of pressures specified by press_plot.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine is celldmf_pt and the Gibbs energy at 
  !   the minimum eminf_pt.
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE control_pressure, ONLY : npress_plot, ipress_plot
  USE temperature, ONLY : ntemp, ntemp_plot, temp_plot
  USE ph_freq_anharmonic,  ONLY : p1tf_t, p2tf_t, p3tf_t, p4tf_t
  USE ph_freq_anharmonic_pt, ONLY : celldmf_pt, eminf_pt
  USE uniform_pressure, ONLY : p2_p, p4_p

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp, ipress, ipressp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO ipressp=1, npress_plot
     ipress=ipress_plot(ipressp)
     DO itemp=1, ntemp

        IF (itemp==1.OR..NOT.lquartic) THEN
           CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p(ipress), p2tf_t(itemp))
!        WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!        CALL write_vector(nvar,x_pos_min)
!        CALL print_genergy(ymin)

        ENDIF
        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                   &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p(ipress), p4tf_t(itemp)) 
!              WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                         &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p(ipress), p3tf_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!              WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p(ipress), p1tf_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!              WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p(ipress), p2tf_t(itemp))
           ENDIF
!           CALL write_vector(nvar,x_pos_min)
!           CALL print_genergy(ymin)
        ENDIF

        eminf_pt(itemp,ipressp)=ymin
        CALL expand_celldm(celldmf_pt(:,itemp,ipressp), x_pos_min, nvar, ibrav)

     ENDDO
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fitf_pt
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_pt_pm()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at all temperatures for the set of pressures specified by press_plot
  !   plus and minus dp.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine is celldm_pt_p1 and celldm_pt_m1. 
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE control_pressure, ONLY : npress, npress_plot, ipress_plot
  USE temperature, ONLY : ntemp, ntemp_plot, temp_plot
  USE anharmonic,  ONLY : p1t_t, p2t_t, p3t_t, p4t_t
  USE anharmonic_pt, ONLY : celldm_pt_p1, celldm_pt_m1
  USE uniform_pressure, ONLY : p2_p, p4_p

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp, ipress, ipressp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO ipressp=1, npress_plot
     ipress=ipress_plot(ipressp)
     IF (ipress==1.OR.ipress==npress) &
                  CALL errore('quadratic_fit_pt_pm',&
                            'increase pmax and/or decrease pmin',1)
     DO itemp=1, ntemp
        CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p(ipress+1), p2t_t(itemp))
!        WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!        CALL write_vector(nvar,x_pos_min)
!        CALL print_genergy(ymin)

        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                   &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p(ipress+1), p4t_t(itemp)) 
!              WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                         &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p(ipress+1), p3t_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!              WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p(ipress+1), p1t_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!              WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p(ipress+1), p2t_t(itemp))
           ENDIF
!           CALL write_vector(nvar,x_pos_min)
!           CALL print_genergy(ymin)
        ENDIF

        CALL expand_celldm(celldm_pt_p1(:,itemp,ipressp), x_pos_min, nvar, ibrav)
     ENDDO
     DO itemp=1, ntemp
        CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p(ipress-1), p2t_t(itemp))
!        WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!        CALL write_vector(nvar,x_pos_min)
!        CALL print_genergy(ymin)

        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                   &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p(ipress-1), p4t_t(itemp)) 
!              WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                         &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p(ipress-1), p3t_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!              WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p(ipress-1), p1t_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!              WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p(ipress-1), p2t_t(itemp))
           ENDIF
        ENDIF

        CALL expand_celldm(celldm_pt_m1(:,itemp,ipressp), x_pos_min, nvar, ibrav)
     ENDDO
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fit_pt_pm
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fitf_pt_pm()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the entalphy to calculate the minimum of the Gibbs energy 
  !   at all temperatures for the set of pressures specified by press_plot
  !   plus and minus dp.
  !   It provides the crystal parameters as a function of temperature
  !   for the given pressures.
  !
  !   The output of this routine is celldmf_pt_p1 and celldmf_pt_m1
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE control_pressure, ONLY : npress, npress_plot, ipress_plot
  USE temperature, ONLY : ntemp, ntemp_plot, temp_plot
  USE ph_freq_anharmonic,  ONLY : p1tf_t, p2tf_t, p3tf_t, p4tf_t
  USE ph_freq_anharmonic_pt, ONLY : celldmf_pt_p1, celldmf_pt_m1
  USE uniform_pressure, ONLY : p2_p, p4_p

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:)
  REAL(DP) :: ymin
  INTEGER  :: itemp, ipress, ipressp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))

  DO ipressp=1, npress_plot
     ipress=ipress_plot(ipressp)
     IF (ipress==1.OR.ipress==npress) &
                  CALL errore('quadratic_fit_pt_pm',&
                            'increase pmax and/or decrease pmin',1)
     DO itemp=1, ntemp
        CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p(ipress+1), p2tf_t(itemp))
!        WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!        CALL write_vector(nvar,x_pos_min)
!        CALL print_genergy(ymin)

        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                   &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p(ipress+1), p4tf_t(itemp)) 
!              WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                         &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p(ipress+1), p3tf_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!              WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p(ipress+1), p1tf_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!              WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p(ipress+1), p2tf_t(itemp))
           ENDIF
!           CALL write_vector(nvar,x_pos_min)
!           CALL print_genergy(ymin)
        ENDIF

        CALL expand_celldm(celldmf_pt_p1(:,itemp,ipressp), x_pos_min, nvar, ibrav)
     ENDDO
     DO itemp=1, ntemp
        CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p(ipress-1), p2tf_t(itemp))
!        WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!        CALL write_vector(nvar,x_pos_min)
!        CALL print_genergy(ymin)

        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                   &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p(ipress-1), p4tf_t(itemp)) 
!              WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!              WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                         &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                       p4_p(ipress-1), p3tf_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!              WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                      p4_p(ipress-1), p1tf_t(itemp))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!              WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                        p4_p(ipress-1), p2tf_t(itemp))
           ENDIF
        ENDIF
        CALL expand_celldm(celldmf_pt_m1(:,itemp,ipressp), x_pos_min, nvar, ibrav)
     ENDDO
  ENDDO
  DEALLOCATE(x_pos_min)
  !
  RETURN
  !
END SUBROUTINE quadratic_fitf_pt_pm
!
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fit_ptt()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational 
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the enthalpy to calculate the minimum of the Gibbs energy 
  !   at all pressures for the set of temperatures specified by temp_plot.
  !   It provides the crystal parameters as a function of pressure
  !   for the given temperatures.
  !
  !   The output of this routine is celldm_pt and the gibbs energy at 
  !   the minimum emin_pt.
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE temperature, ONLY : ntemp_plot, itemp_plot
  USE control_pressure, ONLY : npress
  USE anharmonic,  ONLY : p1t_t, p2t_t, p3t_t, p4t_t
  USE anharmonic_ptt,  ONLY : celldm_ptt, emin_ptt, celldm_ptt_p1, &
                              celldm_ptt_m1, emin_ptt_p1, emin_ptt_m1
  USE uniform_pressure, ONLY : p2_p, p4_p

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:), x_pos_min_p1(:), x_pos_min_m1(:)
  REAL(DP) :: ymin, ymin_p1, ymin_m1
  INTEGER  :: ipress, itempp, itemp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))
  ALLOCATE(x_pos_min_p1(nvar))
  ALLOCATE(x_pos_min_m1(nvar))

  DO itempp=1,ntemp_plot
     itemp=itemp_plot(itempp)
     DO ipress=1, npress
        IF (ipress==1.OR..NOT.lquartic) THEN
           CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                     p2_p(ipress), p2t_t(itemp))
           CALL find_two_quadratic_extremum(nvar, x_pos_min_p1, ymin_p1, &
                                     p2_p(ipress), p2t_t(itemp+1))
           CALL find_two_quadratic_extremum(nvar, x_pos_min_m1, ymin_m1, &
                                     p2_p(ipress), p2t_t(itemp-1))
!        WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!        CALL write_vector(nvar,x_pos_min)
!        CALL print_genergy(ymin)
        ENDIF

        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!             WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                             p4_p(ipress), p4t_t(itemp)) 
              CALL find_two_quartic_extremum(nvar, x_pos_min_p1, ymin_p1, &
                                             p4_p(ipress), p4t_t(itemp+1)) 
              CALL find_two_quartic_extremum(nvar, x_pos_min_m1, ymin_m1, &
                                             p4_p(ipress), p4t_t(itemp-1)) 
!             WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!             WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                      &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                    p4_p(ipress), p3t_t(itemp))
              CALL find_quartic_cubic_extremum(nvar, x_pos_min_p1, ymin_p1, &
                    p4_p(ipress), p3t_t(itemp+1))
              CALL find_quartic_cubic_extremum(nvar, x_pos_min_m1, ymin_m1, &
                    p4_p(ipress), p3t_t(itemp-1))
!             WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!             WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                   p4_p(ipress), p1t_t(itemp))
              CALL find_quartic_linear_extremum(nvar, x_pos_min_p1, ymin_p1, &
                                   p4_p(ipress), p1t_t(itemp+1))
              CALL find_quartic_linear_extremum(nvar, x_pos_min_m1, ymin_m1, &
                                   p4_p(ipress), p1t_t(itemp-1))
!              WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!             WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                     p4_p(ipress), p2t_t(itemp))
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min_p1, &
                     ymin_p1, p4_p(ipress), p2t_t(itemp+1))
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min_m1, &
                     ymin_m1, p4_p(ipress), p2t_t(itemp-1))
           ENDIF
!          CALL write_vector(nvar,x_pos_min)
!          CALL print_genergy(ymin)
        ENDIF
        emin_ptt(ipress,itempp)=ymin
        CALL expand_celldm(celldm_ptt(:,ipress,itempp), x_pos_min, nvar, ibrav)
        emin_ptt_p1(ipress,itempp)=ymin_p1
        CALL expand_celldm(celldm_ptt_p1(:,ipress,itempp), x_pos_min_p1, &
                                                           nvar, ibrav)
        emin_ptt_m1(ipress,itempp)=ymin_m1
        CALL expand_celldm(celldm_ptt_m1(:,ipress,itempp), x_pos_min_m1, &
                                                           nvar, ibrav)
     ENDDO
  ENDDO
  DEALLOCATE(x_pos_min)
  DEALLOCATE(x_pos_min_p1)
  DEALLOCATE(x_pos_min_m1)
  !
  RETURN
END SUBROUTINE quadratic_fit_ptt
!-----------------------------------------------------------------------
SUBROUTINE quadratic_fitf_ptt()
  !-----------------------------------------------------------------------
  !
  !   This uses the polynomial interpolation of the vibrational 
  !   (plus electronic) free energy and the polynomial interpolation 
  !   of the enthalpy to calculate the minimum of the Gibbs energy 
  !   at all pressures for the set of temperatures specified by temp_plot.
  !   It provides the crystal parameters as a function of pressure
  !   for the given temperatures.
  !
  !   The output of this routine is celldmf_pt and the Gibbs energy at 
  !   the minimum eminf_pt.
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : expand_celldm, crystal_parameters
  USE io_global,   ONLY : stdout
  USE temperature, ONLY : ntemp_plot, itemp_plot
  USE control_pressure, ONLY : npress
  USE ph_freq_anharmonic,  ONLY : p1tf_t, p2tf_t, p3tf_t, p4tf_t
  USE ph_freq_anharmonic_ptt, ONLY : celldmf_ptt, eminf_ptt, celldmf_ptt_p1, &
                                     eminf_ptt_p1, celldmf_ptt_m1,           &
                                     eminf_ptt_m1
  USE uniform_pressure, ONLY : p2_p, p4_p

  USE linear_surfaces,  ONLY : fit_multi_linear, print_chisq_linear
  USE quadratic_surfaces, ONLY : find_two_quadratic_extremum
  USE quartic_surfaces, ONLY : find_quartic_quadratic_extremum,     &
                      find_two_quartic_extremum, find_quartic_cubic_extremum,&
                      find_quartic_linear_extremum
  USE vector_mod, ONLY : write_vector

  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: x_pos_min(:), x_pos_min_p1(:), x_pos_min_m1(:)
  REAL(DP) :: ymin, ymin_p1, ymin_m1
  INTEGER  :: ipress, itempp, itemp
  INTEGER  :: compute_nwork, compute_nwork_ph
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos_min(nvar))
  ALLOCATE(x_pos_min_p1(nvar))
  ALLOCATE(x_pos_min_m1(nvar))

  DO itempp=1, ntemp_plot
     itemp=itemp_plot(itempp)
     DO ipress=1, npress
        IF (ipress==1.OR..NOT.lquartic) THEN
           CALL find_two_quadratic_extremum(nvar, x_pos_min, ymin, &
                                        p2_p(ipress), p2tf_t(itemp))
           CALL find_two_quadratic_extremum(nvar, x_pos_min_p1, ymin_p1, &
                                        p2_p(ipress), p2tf_t(itemp+1))
           CALL find_two_quadratic_extremum(nvar, x_pos_min_m1, ymin_m1, &
                                        p2_p(ipress), p2tf_t(itemp-1))
!          WRITE(stdout,'(/,5x,"Extremum of the quadratic found at:")')
!          CALL write_vector(nvar,x_pos_min)
!          CALL print_genergy(ymin)
        ENDIF

        IF (lquartic) THEN
           IF (poly_degree_ph==4) THEN
!           WRITE(stdout,'(/,5x, "Fit improved with a fourth &
!                                                &order polynomial")') 
              CALL find_two_quartic_extremum(nvar, x_pos_min, ymin, &
                                          p4_p(ipress), p4tf_t(itemp)) 
              CALL find_two_quartic_extremum(nvar, x_pos_min_p1, ymin_p1, &
                                          p4_p(ipress), p4tf_t(itemp+1)) 
              CALL find_two_quartic_extremum(nvar, x_pos_min_m1, ymin_m1, &
                                          p4_p(ipress), p4tf_t(itemp-1)) 
!           WRITE(stdout,'(/,5x,"Extremum of the quartic found at:")')
           ELSEIF (poly_degree_ph==3) THEN
!           WRITE(stdout,'(/,5x, "Fit improved with a third order &
!                                                      &polynomial")') 
              CALL find_quartic_cubic_extremum(nvar, x_pos_min, ymin, &
                    p4_p(ipress), p3tf_t(itemp))
              CALL find_quartic_cubic_extremum(nvar, x_pos_min_p1, ymin_p1, &
                    p4_p(ipress), p3tf_t(itemp+1))
              CALL find_quartic_cubic_extremum(nvar, x_pos_min_m1, ymin_m1, &
                    p4_p(ipress), p3tf_t(itemp-1))
!           WRITE(stdout,'(/,5x,"Extremum of the quartic+cubic found at:")')
           ELSEIF (poly_degree_ph==1) THEN
!           WRITE(stdout,'(/,5x, "Fit with a fist order polynomial")') 
              CALL find_quartic_linear_extremum(nvar, x_pos_min, ymin, &
                                   p4_p(ipress), p1tf_t(itemp))
              CALL find_quartic_linear_extremum(nvar, x_pos_min_p1, ymin_p1, &
                                   p4_p(ipress), p1tf_t(itemp+1))
              CALL find_quartic_linear_extremum(nvar, x_pos_min_m1, ymin_m1, &
                                   p4_p(ipress), p1tf_t(itemp-1))
!           WRITE(stdout,'(/,5x,"Extremum of the quartic+linear found at:")')
           ELSE
!           WRITE(stdout,'(/,5x,"Quartic fit used only a T=0:")')
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min, ymin, &
                     p4_p(ipress), p2tf_t(itemp))
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min_p1, &
                                   ymin_p1, p4_p(ipress), p2tf_t(itemp+1))
              CALL find_quartic_quadratic_extremum(nvar, x_pos_min_m1, &
                                   ymin_m1, p4_p(ipress), p2tf_t(itemp-1))
           ENDIF
!           CALL write_vector(nvar,x_pos_min)
!           CALL print_genergy(ymin)
        ENDIF

        eminf_ptt(ipress,itempp)=ymin
        CALL expand_celldm(celldmf_ptt(:,ipress,itempp), x_pos_min, nvar, ibrav)

        eminf_ptt_p1(ipress,itempp)=ymin_p1
        CALL expand_celldm(celldmf_ptt_p1(:,ipress,itempp), x_pos_min_p1, &
                                                                 nvar, ibrav)
        eminf_ptt_m1(ipress,itempp)=ymin_m1
        CALL expand_celldm(celldmf_ptt_m1(:,ipress,itempp), x_pos_min_m1, &
                                                                 nvar, ibrav)
     ENDDO
  ENDDO
  DEALLOCATE(x_pos_min)
  DEALLOCATE(x_pos_min_p1)
  DEALLOCATE(x_pos_min_m1)
  !
  RETURN
END SUBROUTINE quadratic_fitf_ptt
!
!-------------------------------------------------------------------
SUBROUTINE set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)
!-------------------------------------------------------------------
!
!  this rouotine receives and array of values of celldm celldm_geo(6,ndata)
!  and transform it in a compact array x(nvar,ndata), where nvar depends
!  on the Bravais lattice
!
USE kinds, ONLY : DP
USE lattices, ONLY : compress_celldm

IMPLICIT NONE
INTEGER, INTENT(IN)     :: ibrav, nvar, ndata
REAL(DP), INTENT(IN)    :: celldm_geo(6,ndata)
REAL(DP), INTENT(INOUT) :: x(nvar,ndata)

INTEGER :: idata

DO idata=1,ndata
   CALL compress_celldm(celldm_geo(1,idata), x(1,idata), nvar, ibrav)
ENDDO

RETURN
END SUBROUTINE set_x_from_celldm

!-------------------------------------------------------------------
SUBROUTINE print_genergy(ymin)
!-------------------------------------------------------------------
USE kinds, ONLY : DP
USE control_pressure, ONLY : pressure
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) :: ymin

IF (pressure > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Gibbs energy at the extremum:",f22.12," Ry")') ymin
ELSE
   WRITE(stdout,'(5x,"Free energy at the extremum:",f22.12," Ry")') ymin
ENDIF

RETURN
END SUBROUTINE print_genergy

!-------------------------------------------------------------------
SUBROUTINE print_energy(ymin)
!-------------------------------------------------------------------

USE kinds, ONLY : DP
USE control_pressure, ONLY : pressure
USE io_global, ONLY : stdout

IMPLICIT NONE
REAL(DP) :: ymin

IF (pressure > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Enthalpy at the extremum:",f21.12," Ry")') ymin
ELSE
   WRITE(stdout,'(5x,"Energy at the extremum:",f21.12," Ry")') ymin
ENDIF

RETURN
END SUBROUTINE print_energy
!
!-----------------------------------------------------------------------
SUBROUTINE compute_anhar_poly(celldm_t, free_e_min_t, pt1, pt2, pt3, pt4 )
  !-----------------------------------------------------------------------
  !
  !   This routine receives the polynomial interpolation of the 
  !   free energy and computes the free energy at the celldm given in input.
  !
  !   The output of this routine is free_e_min_t at the given temperature 
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : ibrav
  USE control_quadratic_energy, ONLY : nvar
  USE control_quartic_energy, ONLY :  lquartic, poly_degree_ph
  USE lattices,    ONLY : compress_celldm, crystal_parameters
  USE io_global,   ONLY : stdout

  USE linear_surfaces, ONLY : evaluate_fit_linear
  USE quadratic_surfaces, ONLY : evaluate_fit_quadratic
  USE cubic_surfaces, ONLY : evaluate_fit_cubic
  USE quartic_surfaces, ONLY : evaluate_fit_quartic

  USE polynomial, ONLY : poly1, poly2, poly3, poly4

  IMPLICIT NONE
  REAL(DP) :: celldm_t(6), free_e_min_t
  REAL(DP), ALLOCATABLE ::  x_pos(:)
  TYPE(poly1) :: pt1
  TYPE(poly2) :: pt2
  TYPE(poly3) :: pt3             
  TYPE(poly4) :: pt4             
  REAL(DP) :: ymin
  !
  nvar=crystal_parameters(ibrav)
  !
  ALLOCATE(x_pos(nvar))
  CALL compress_celldm(celldm_t, x_pos, nvar, ibrav)
  CALL evaluate_fit_quadratic(nvar,x_pos,ymin,pt2)
  
  IF (lquartic) THEN
     IF (poly_degree_ph==4) THEN
        CALL evaluate_fit_quartic(nvar,x_pos,ymin,pt4)
     ELSEIF (poly_degree_ph==3) THEN
        CALL evaluate_fit_cubic(nvar,x_pos,ymin,pt3)
     ELSEIF (poly_degree_ph==1) THEN
        CALL evaluate_fit_linear(nvar,x_pos,ymin,pt1)
     ENDIF
  ENDIF

  free_e_min_t=ymin

  DEALLOCATE(x_pos)
  !
  RETURN
  !
END SUBROUTINE compute_anhar_poly
