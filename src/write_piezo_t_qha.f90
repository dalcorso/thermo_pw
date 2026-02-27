!
! Copyright (C) 2026 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE write_piezo_t_qha()
!---------------------------------------------------------------------------
!
!  This routine interpolates the temperature dependent FFEM
!  piezoelectric tensor at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at temperature T.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : omega_geo_eos, celldm_geo_eos
USE thermo_sym, ONLY : code_group_save
USE initial_conf, ONLY : ibrav_save
USE control_mur, ONLY : lmurn
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE thermodynamics, ONLY : e_piezo_tensor_eos_t
USE ph_freq_thermodynamics, ONLY : e_piezo_tensorf_eos_t
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_piezoelectric_tensor, ONLY : piezo_qha_geo_available, &
                             piezof_qha_geo_available, lpiezo, lpiezof, &
                             found_ph_pt, found_dos_pt, fnph_pt, fndos_pt, &
                             f_geodos_pt, f_geoph_pt
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE piezoelectric_tensor, ONLY : write_piezo_tensor_on_file
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE anharmonic,     ONLY : celldm_t, e_piezo_tensor_t
USE ph_freq_anharmonic, ONLY : celldmf_t, e_piezo_tensorf_t
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filepiezo
INTEGER :: igeo, nvar, ndata, jdata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)
REAL(DP) :: compute_omega_geo
TYPE(poly1) :: pt_p1
TYPE(poly2) :: pt_p2
TYPE(poly3) :: pt_p3
TYPE(poly4) :: pt_p4

INTEGER :: i, j, idata, itemp, startt, lastt, pdelc_dos, pdelc_ph
INTEGER :: compute_nwork 

nvar=crystal_parameters(ibrav_save)
IF (lmurn) nvar=1

ndata=compute_nwork()

ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(x1(nvar,ndata))
!
!  Part 1 evaluation of the polynomial coefficients
!
ALLOCATE(xfit(nvar))
CALL divide(world_comm, ntemp, startt, lastt)

IF (lmurn) THEN
   IF (lgeo_from_file) THEN
      DO idata=1,ndata
         x(1,idata)=omega_geo_eos(idata)
      ENDDO
   ELSE
      DO idata=1,ndata
         x(1,idata)=celldm_geo_eos(1,idata)
      ENDDO
   ENDIF
ELSE
   CALL set_x_from_celldm(ibrav_save, nvar, ndata, x, celldm_geo_eos)
ENDIF

e_piezo_tensor_t=0.0_DP
e_piezo_tensorf_t=0.0_DP
pdelc_dos=MIN(poly_degree_elc,fndos_pt-1)
pdelc_ph=MIN(poly_degree_elc,fnph_pt-1)
IF (pdelc_dos/=poly_degree_elc) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for dos")') pdelc_dos
IF (pdelc_ph/=poly_degree_elc) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for ph")') pdelc_ph

DO itemp=startt,lastt
   IF (itemp==1.OR.itemp==ntemp) CYCLE
   
   IF (piezo_qha_geo_available.AND.ltherm_dos) THEN
      IF (lgeo_from_file) THEN
         xfit(1)=compute_omega_geo(ibrav_save, celldm_t(:,itemp))
      ELSE
         CALL compress_celldm(celldm_t(:,itemp),xfit,nvar,ibrav_save)
      ENDIF
      f=0.0_DP
      DO i=1,3
         DO j=1,6
            IF (ABS(e_piezo_tensor_eos_t(i,j,itemp,f_geodos_pt))>1.D-7) THEN
               WRITE(stdout,'(/,5x,"Fitting the piezoelectric tensor &
                                     &e(",i4,",",i4,")")') i,j
               WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

               jdata=0
               DO idata=1,ndata
                  IF (found_dos_pt(idata)) THEN
                     jdata=jdata+1
                     x1(:,jdata)=x(:,idata)
                     f(jdata)=e_piezo_tensor_eos_t(i,j,itemp,idata)
                  ENDIF
               END DO

               IF (pdelc_dos==4) THEN
                  CALL init_poly(nvar,pt_p4)
                  CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,pt_p4) 
                  CALL evaluate_fit_quartic(nvar, xfit, &
                                    e_piezo_tensor_t(i,j,itemp),pt_p4)
                  CALL clean_poly(pt_p4)
               ELSEIF (pdelc_dos==3) THEN
                  CALL init_poly(nvar,pt_p3)
                  CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,pt_p3)
                  CALL evaluate_fit_cubic(nvar,xfit, &
                                e_piezo_tensor_t(i,j,itemp), pt_p3)
                  CALL clean_poly(pt_p3)
               ELSEIF (pdelc_dos==2) THEN
                  CALL init_poly(nvar,pt_p2)
                  CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f, pt_p2)
                  CALL evaluate_fit_quadratic(nvar,xfit,&
                                   e_piezo_tensor_t(i,j,itemp), pt_p2)
                  CALL clean_poly(pt_p2)
               ELSEIF (pdelc_dos==1) THEN
                  CALL init_poly(nvar,pt_p1)
                  CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,pt_p1)
                  CALL evaluate_fit_linear(nvar, xfit, &
                                       e_piezo_tensor_t(i,j,itemp), pt_p1)
                  CALL clean_poly(pt_p1)
               ELSE
                  CALL errore('write_piezo_t_qha','wrong poly_degree_elc',1)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   END IF

   IF (piezof_qha_geo_available.AND.ltherm_freq) THEN
      IF (lgeo_from_file) THEN
         xfit(1)=compute_omega_geo(ibrav_save, celldmf_t(:,itemp))
      ELSE
         CALL compress_celldm(celldmf_t(:,itemp),xfit,nvar,ibrav_save)
      ENDIF
      f=0.0_DP
      DO i=1,3
         DO j=1,6
            IF (ABS(e_piezo_tensorf_eos_t(i,j,itemp,f_geoph_pt))>1D-7) THEN
               WRITE(stdout,'(/,5x,"Fitting piezoelectric tensor &
                                         &e(",i4,",",i4,")")') i,j
               WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

               jdata=0
               DO idata=1,ndata
                  IF (found_ph_pt(idata)) THEN
                     jdata=jdata+1
                     x1(:,jdata)=x(:,idata)
                     f(jdata)=e_piezo_tensorf_eos_t(i,j,itemp,idata)
                  ENDIF
               END DO

               IF (pdelc_ph==4) THEN
                  CALL init_poly(nvar,pt_p4)
                  CALL fit_multi_quartic(jdata, nvar, lsolve, x1, f, pt_p4)
                  CALL evaluate_fit_quartic(nvar,xfit,&
                                       e_piezo_tensorf_t(i,j,itemp), pt_p4) 
                  CALL clean_poly(pt_p4)
               ELSEIF (pdelc_ph==3) THEN
                  CALL init_poly(nvar,pt_p3)
                  CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,pt_p3)
                  CALL evaluate_fit_cubic(nvar,xfit, &
                                     e_piezo_tensorf_t(i,j,itemp), pt_p3)
                  CALL clean_poly(pt_p3)
               ELSEIF (pdelc_ph==2) THEN
                  CALL init_poly(nvar,pt_p2)
                  CALL fit_multi_quadratic(jdata, nvar, lsolve, x1, f, pt_p2)
                  CALL evaluate_fit_quadratic(nvar,xfit, &
                                   e_piezo_tensorf_t(i,j,itemp), pt_p2)
                  CALL clean_poly(pt_p2)
               ELSEIF (pdelc_ph==1) THEN
                  CALL init_poly(nvar,pt_p1)
                  CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,pt_p1)
                  CALL evaluate_fit_linear(nvar,xfit,&
                                    e_piezo_tensorf_t(i,j,itemp), pt_p1)
                  CALL clean_poly(pt_p1)
               ELSE
                  CALL errore('write_piezo_t_qha','wrong poly_degree_elc',1)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDIF
ENDDO
!

IF (ltherm_dos) THEN
   CALL mp_sum(e_piezo_tensor_t, world_comm)
   lpiezo=.TRUE.
   filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo'
   CALL write_piezo_tensor_on_file(temp, ntemp, ibrav_save,          &
                       code_group_save, e_piezo_tensor_t(:,:,:),     &
                                              filepiezo, 0, 1)
ENDIF
!
IF (ltherm_freq) THEN
   CALL mp_sum(e_piezo_tensorf_t, world_comm)
   lpiezof=.TRUE.
   filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_ph'
   CALL write_piezo_tensor_on_file(temp, ntemp, ibrav_save,          &
                       code_group_save, e_piezo_tensorf_t(:,:,:),    &
                                              filepiezo, 0, 1)
ENDIF

DEALLOCATE(x1)
DEALLOCATE(f)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_piezo_t_qha

!---------------------------------------------------------------------------
SUBROUTINE write_piezo_pt_qha()
!---------------------------------------------------------------------------
!
!  This routine interpolates the temperature dependent FFEM
!  piezoelectric tensors at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at temperature T for several pressures.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : omega_geo_eos, celldm_geo_eos
USE thermo_sym, ONLY : code_group_save
USE initial_conf, ONLY : ibrav_save
USE control_pressure, ONLY : npress_plot, ipress_plot, press
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE control_mur, ONLY : lmurn
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_piezoelectric_tensor, ONLY : piezo_qha_geo_available, &
                             piezof_qha_geo_available, lpiezo_pt, &
                             lpiezof, lpiezof_pt, &
                             found_ph_pt, found_dos_pt, fnph_pt, fndos_pt, &
                             f_geodos_pt, f_geoph_pt
USE lattices,       ONLY : compress_celldm
USE piezoelectric_tensor, ONLY : write_piezo_tensor_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE thermodynamics,     ONLY : e_piezo_tensor_eos_t
USE ph_freq_thermodynamics, ONLY : e_piezo_tensorf_eos_t
USE anharmonic_pt,  ONLY : celldm_pt, e_piezo_tensor_pt
USE ph_freq_anharmonic_pt,  ONLY : celldmf_pt, e_piezo_tensorf_pt
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filepiezo
INTEGER :: igeo, nvar, ndata, jdata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)
REAL(DP) :: compute_omega_geo
TYPE(poly1) :: pt_p1
TYPE(poly2) :: pt_p2
TYPE(poly3) :: pt_p3
TYPE(poly4) :: pt_p4

INTEGER :: i, j, idata, itemp, ipressp, ipress, startt, lastt, &
                 pdelc_dos, pdelc_ph
INTEGER :: compute_nwork 

IF (npress_plot==0) RETURN

nvar=crystal_parameters(ibrav_save)
IF (lmurn) nvar=1

ndata=compute_nwork()

ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(x1(nvar,ndata))
!
!  Part 1 evaluation of the polynomial coefficients
!
ALLOCATE(xfit(nvar))
CALL divide(world_comm, ntemp, startt, lastt)

IF (lmurn) THEN
   IF (lgeo_from_file) THEN
      DO idata=1,ndata
         x(1,idata)=omega_geo_eos(idata)
      ENDDO
   ELSE
      DO idata=1,ndata
         x(1,idata)=celldm_geo_eos(1,idata)
      ENDDO
   ENDIF
ELSE
   CALL set_x_from_celldm(ibrav_save, nvar, ndata, x, celldm_geo_eos)
ENDIF

e_piezo_tensor_pt=0.0_DP

e_piezo_tensorf_pt=0.0_DP

pdelc_dos=MIN(poly_degree_elc,fndos_pt-1)
pdelc_ph=MIN(poly_degree_elc,fnph_pt-1)
IF (pdelc_dos/=poly_degree_elc.AND.ltherm_dos) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for dos")') pdelc_dos
IF (pdelc_ph/=poly_degree_elc.AND.ltherm_freq) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for ph")') pdelc_ph

DO ipressp=1,npress_plot
   ipress=ipress_plot(ipressp)
   IF (ltherm_dos.AND.piezo_qha_geo_available) THEN
      DO itemp=startt,lastt
         IF (itemp==1.OR.itemp==ntemp) CYCLE
   
         IF (lgeo_from_file) THEN
            xfit(1)=compute_omega_geo(ibrav_save, celldm_pt(:,itemp,ipressp))
         ELSE
            CALL compress_celldm(celldm_pt(:,itemp,ipressp),xfit,nvar,&
                                                                   ibrav_save)
         ENDIF
         f=0.0_DP
         DO i=1,3
            DO j=1,6
               IF (ABS(e_piezo_tensor_eos_t(i,j,itemp,f_geodos_pt))>1.D-7) THEN
                  WRITE(stdout,'(/,5x,"Fitting piezoelectric tensor &
                                  &e(",i4,",",i4,")")') i,j
                  WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

                  jdata=0
                  DO idata=1,ndata
                     IF (found_dos_pt(idata)) THEN
                        jdata=jdata+1
                        x1(:,jdata)=x(:,idata)
                        f(jdata)=e_piezo_tensor_eos_t(i,j,itemp,idata)
                     ENDIF
                  END DO

                  IF (pdelc_dos==4) THEN
                     CALL init_poly(nvar,pt_p4)
                     CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,pt_p4) 
                     CALL evaluate_fit_quartic(nvar, xfit, &
                                   e_piezo_tensor_pt(i,j,itemp,ipressp),pt_p4)
                     CALL clean_poly(pt_p4)
                  ELSEIF (pdelc_dos==3) THEN
                     CALL init_poly(nvar,pt_p3)
                     CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,pt_p3)
                     CALL evaluate_fit_cubic(nvar,xfit, &
                       e_piezo_tensor_pt(i,j,itemp,ipressp), pt_p3)
                     CALL clean_poly(pt_p3)
                  ELSEIF (pdelc_dos==2) THEN
                     CALL init_poly(nvar,pt_p2)
                     CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f,pt_p2)
                     CALL evaluate_fit_quadratic(nvar,xfit, &
                                e_piezo_tensor_pt(i,j,itemp,ipressp), pt_p2)
                     CALL clean_poly(pt_p2)
                  ELSEIF (pdelc_dos==1) THEN
                     CALL init_poly(nvar,pt_p1)
                     CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,pt_p1)
                     CALL evaluate_fit_linear(nvar, xfit, &
                         e_piezo_tensor_pt(i,j,itemp,ipressp), pt_p1)
                     CALL clean_poly(pt_p1)
                  ELSE
                     CALL errore('write_piezo_pt_qha',&
                                              'wrong poly_degree_elc',1)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      CALL mp_sum(e_piezo_tensor_pt(:,:,:,ipressp), world_comm)
      lpiezo_pt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_press'
      CALL add_value(filepiezo, press(ipress))
      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav_save,          &
           code_group_save, e_piezo_tensor_pt(:,:,:,ipressp), filepiezo, 0, 1)
   ENDIF

   IF (ltherm_freq.AND.piezof_qha_geo_available) THEN
      DO itemp=startt,lastt
         IF (itemp==1.OR.itemp==ntemp) CYCLE
   
         IF (lgeo_from_file) THEN
            xfit(1)=compute_omega_geo(ibrav_save, celldmf_pt(:,itemp,ipressp))
         ELSE
            CALL compress_celldm(celldmf_pt(:,itemp,ipressp),xfit,&
                                                        nvar,ibrav_save)
         ENDIF
         f=0.0_DP
         DO i=1,3
            DO j=1,6
               IF (ABS(e_piezo_tensorf_eos_t(i,j,itemp,f_geoph_pt))>1.D-7) THEN
                  WRITE(stdout,'(/,5x,"Fitting piezoelectric tensor &
                              &e(",i4,",",i4,")")') i,j
                  WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

                  jdata=0
                  DO idata=1,ndata
                     IF (found_ph_pt(idata)) THEN
                        jdata=jdata+1
                        x1(:,jdata)=x(:,idata)
                        f(jdata)=e_piezo_tensorf_eos_t(i,j,itemp,idata)
                     ENDIF
                  END DO

                  IF (pdelc_ph==4) THEN
                     CALL init_poly(nvar,pt_p4)
                     CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,pt_p4) 
                     CALL evaluate_fit_quartic(nvar, xfit, &
                                   e_piezo_tensorf_pt(i,j,itemp,ipressp),pt_p4)
                     CALL clean_poly(pt_p4)
                  ELSEIF (pdelc_ph==3) THEN
                     CALL init_poly(nvar,pt_p3)
                     CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,pt_p3)
                     CALL evaluate_fit_cubic(nvar,xfit, &
                       e_piezo_tensorf_pt(i,j,itemp,ipressp), pt_p3)
                     CALL clean_poly(pt_p3)
                  ELSEIF (pdelc_ph==2) THEN
                     CALL init_poly(nvar,pt_p2)
                     CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f,pt_p2)
                     CALL evaluate_fit_quadratic(nvar,xfit, &
                                 e_piezo_tensorf_pt(i,j,itemp,ipressp), pt_p2)
                     CALL clean_poly(pt_p2)
                  ELSEIF (pdelc_ph==1) THEN
                     CALL init_poly(nvar,pt_p1)
                     CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,pt_p1)
                     CALL evaluate_fit_linear(nvar, xfit, &
                         e_piezo_tensorf_pt(i,j,itemp,ipressp), pt_p1)
                     CALL clean_poly(pt_p1)
                  ELSE
                     CALL errore('write_piezo_pt_qha',&
                                           'wrong poly_degree_elc',1)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      CALL mp_sum(e_piezo_tensorf_pt(:,:,:,ipressp), world_comm)
      lpiezof_pt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_ph_press'
      CALL add_value(filepiezo, press(ipress))
      CALL write_piezo_tensor_on_file(temp, ntemp, ibrav_save,          &
           code_group_save, e_piezo_tensorf_pt(:,:,:,ipressp), filepiezo, 0, 1)
   ENDIF
ENDDO

DEALLOCATE(x1)
DEALLOCATE(f)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_piezo_pt_qha

!---------------------------------------------------------------------------
SUBROUTINE write_piezo_ptt_qha()
!---------------------------------------------------------------------------
!
!  This routine interpolates the temperature dependent FFEM 
!  piezoelectric tensor at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at pressure p for several temperatures.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : omega_geo_eos, celldm_geo_eos
USE thermo_sym, ONLY : code_group_save
USE initial_conf, ONLY : ibrav_save
USE control_mur, ONLY : lmurn
USE control_pressure, ONLY : npress, press
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_piezoelectric_tensor, ONLY : piezo_qha_geo_available, &
                             piezof_qha_geo_available, lpiezo_ptt, &
                             lpiezof_ptt, &
                             found_ph_pt, found_dos_pt, fnph_pt, fndos_pt, &
                             f_geodos_pt, f_geoph_pt
USE lattices,       ONLY : compress_celldm
USE piezoelectric_tensor, ONLY : write_piezo_tensor_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq, lgeo_from_file
USE thermodynamics,     ONLY : e_piezo_tensor_eos_t
USE ph_freq_thermodynamics,     ONLY : e_piezo_tensorf_eos_t
USE anharmonic_ptt,  ONLY : celldm_ptt, e_piezo_tensor_ptt
USE ph_freq_anharmonic_ptt,  ONLY : celldmf_ptt, e_piezo_tensorf_ptt 
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp, ntemp_plot, itemp_plot
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filepiezo
INTEGER :: igeo, nvar, ndata, jdata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)
REAL(DP) :: compute_omega_geo
TYPE(poly1) :: pt_p1
TYPE(poly2) :: pt_p2
TYPE(poly3) :: pt_p3
TYPE(poly4) :: pt_p4

INTEGER :: i, j, idata, itemp, itempp, ipress, startp, lastp, &
                 pdelc_dos, pdelc_ph
INTEGER :: compute_nwork 

IF (ntemp_plot==0) RETURN

nvar=crystal_parameters(ibrav_save)
IF (lmurn) nvar=1

ndata=compute_nwork()

ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(x1(nvar,ndata))
!
ALLOCATE(xfit(nvar))
CALL divide(world_comm, npress, startp, lastp)

IF (lmurn) THEN
   IF (lgeo_from_file) THEN
      DO idata=1,ndata
         x(1,idata)=omega_geo_eos(idata)
      ENDDO
   ELSE
      DO idata=1,ndata
         x(1,idata)=celldm_geo_eos(1,idata)
      ENDDO
   ENDIF
ELSE
   CALL set_x_from_celldm(ibrav_save, nvar, ndata, x, celldm_geo_eos)
ENDIF

e_piezo_tensor_ptt=0.0_DP
e_piezo_tensorf_ptt=0.0_DP

pdelc_dos=MIN(poly_degree_elc,fndos_pt-1)
pdelc_ph=MIN(poly_degree_elc,fnph_pt-1)
IF (pdelc_dos/=poly_degree_elc.AND.ltherm_dos) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for dos")') pdelc_dos
IF (pdelc_ph/=poly_degree_elc.AND.ltherm_freq) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for ph")') pdelc_ph

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   IF (ltherm_dos.AND.piezo_qha_geo_available) THEN
      DO ipress=startp,lastp
         IF (ipress==1.OR.ipress==npress) CYCLE
         IF (lgeo_from_file) THEN
            xfit(1)=compute_omega_geo(ibrav_save, celldm_ptt(:,ipress,itempp))
         ELSE
            CALL compress_celldm(celldm_ptt(:,ipress,itempp),xfit,&
                                                          nvar,ibrav_save)
         ENDIF
         f=0.0_DP
         DO i=1,3
            DO j=1,6
               IF (ABS(e_piezo_tensor_eos_t(i,j,itemp,f_geodos_pt))>1.D-7) THEN
                  jdata=0
                  DO idata=1,ndata
                     IF (found_dos_pt(idata)) THEN
                        jdata=jdata+1
                        x1(:,jdata)=x(:,idata)
                        f(jdata)=e_piezo_tensor_eos_t(i,j,itemp,idata)
                     ENDIF
                  END DO

                  IF (pdelc_dos==4) THEN
                     CALL init_poly(nvar,pt_p4)
                     CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,pt_p4) 
                     CALL evaluate_fit_quartic(nvar, xfit, &
                            e_piezo_tensor_ptt(i,j,ipress,itempp), pt_p4)
                     CALL clean_poly(pt_p4)
                  ELSEIF (pdelc_dos==3) THEN
                     CALL init_poly(nvar,pt_p3)
                     CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,pt_p3)
                     CALL evaluate_fit_cubic(nvar,xfit, &
                         e_piezo_tensor_ptt(i,j,ipress,itempp), pt_p3)
                     CALL clean_poly(pt_p3)
                  ELSEIF (pdelc_dos==2) THEN
                     CALL init_poly(nvar,pt_p2)
                     CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f,pt_p2)
                     CALL evaluate_fit_quadratic(nvar,xfit, &
                               e_piezo_tensor_ptt(i,j,ipress,itempp), pt_p2)
                     CALL clean_poly(pt_p2)
                  ELSEIF (pdelc_dos==1) THEN
                     CALL init_poly(nvar,pt_p1)
                     CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,pt_p1)
                     CALL evaluate_fit_linear(nvar, xfit, &
                         e_piezo_tensor_ptt(i,j,ipress,itempp), pt_p1)
                     CALL clean_poly(pt_p1)
                  ELSE
                     CALL errore('write_piezo_ptt_qha',&
                                            'wrong poly_degree_elc',1)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      CALL mp_sum(e_piezo_tensor_ptt(:,:,:,itempp), world_comm)
      lpiezo_ptt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_temp'
      CALL add_value(filepiezo, temp(itemp))
      CALL write_piezo_tensor_on_file(press, npress, ibrav_save,          &
           code_group_save, e_piezo_tensor_ptt(:,:,:,itempp), filepiezo, 2, 1)
   ENDIF

   IF (ltherm_freq.AND.piezof_qha_geo_available) THEN
      DO ipress=startp,lastp
         IF (ipress==1.OR.ipress==npress) CYCLE
         IF (lgeo_from_file) THEN
            xfit(1)=compute_omega_geo(ibrav_save, celldmf_ptt(:,ipress,itempp))
         ELSE
            CALL compress_celldm(celldmf_ptt(:,ipress,itempp),xfit,&
                                                         nvar,ibrav_save)
         ENDIF
         f=0.0_DP
         DO i=1,3
            DO j=1,6
               IF (ABS(e_piezo_tensorf_eos_t(i,j,itemp,f_geoph_pt))>1.D-7) THEN
                  jdata=0
                  DO idata=1,ndata
                     IF (found_ph_pt(idata)) THEN
                        jdata=jdata+1
                        x1(:,jdata)=x(:,idata)
                        f(jdata)=e_piezo_tensorf_eos_t(i,j,itemp,idata)
                     ENDIF
                  END DO

                  IF (pdelc_ph==4) THEN
                     CALL init_poly(nvar,pt_p4)
                     CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,pt_p4) 
                     CALL evaluate_fit_quartic(nvar, xfit, &
                            e_piezo_tensorf_ptt(i,j,ipress,itempp), pt_p4)
                     CALL clean_poly(pt_p4)
                  ELSEIF (pdelc_ph==3) THEN
                     CALL init_poly(nvar,pt_p3)
                     CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,pt_p3)
                     CALL evaluate_fit_cubic(nvar,xfit, &
                         e_piezo_tensorf_ptt(i,j,ipress,itempp), pt_p3)
                     CALL clean_poly(pt_p3)
                  ELSEIF (pdelc_ph==2) THEN
                     CALL init_poly(nvar,pt_p2)
                     CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f,pt_p2)
                     CALL evaluate_fit_quadratic(nvar,xfit, &
                               e_piezo_tensorf_ptt(i,j,ipress,itempp), pt_p2)
                     CALL clean_poly(pt_p2)
                  ELSEIF (pdelc_dos==1) THEN
                     CALL init_poly(nvar,pt_p1)
                     CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,pt_p1)
                     CALL evaluate_fit_linear(nvar, xfit, &
                         e_piezo_tensorf_ptt(i,j,ipress,itempp), pt_p1)
                     CALL clean_poly(pt_p1)
                  ELSE
                     CALL errore('write_piezo_ptt_qha',&
                                              'wrong poly_degree_elc',1)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      CALL mp_sum(e_piezo_tensorf_ptt(:,:,:,itempp), world_comm)
      lpiezof_ptt=.TRUE.
      filepiezo='anhar_files/'//TRIM(flanhar)//'.e_piezo_ph_temp'
      CALL add_value(filepiezo, temp(itemp))
      CALL write_piezo_tensor_on_file(press, npress, ibrav_save,          &
           code_group_save, e_piezo_tensorf_ptt(:,:,:,itempp), filepiezo, 2, 1)
   ENDIF
ENDDO

DEALLOCATE(x1)
DEALLOCATE(f)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_piezo_ptt_qha

