!
! Copyright (C) 2019 C. Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE write_elastic_t_qha()
!---------------------------------------------------------------------------
!
!  This routine interpolates the temperature dependent quasi harmonic 
!  elastic constants at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at temperature T.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE thermo_sym, ONLY : laue
USE control_mur, ONLY : lmurn
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_elastic_constants, ONLY : el_cons_qha_geo_available, &
                             el_consf_qha_geo_available, lelastic, lelasticf, &
                             found_ph_ec, found_dos_ec, fnph_ec, fndos_ec, &
                             f_geodos_ec, f_geoph_ec
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE elastic_constants, ONLY : el_con, el_compliances, write_el_cons_on_file, &
                              compute_elastic_compliances, &
                              print_macro_elasticity
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE anharmonic,     ONLY : celldm_t, el_cons_t, el_comp_t, b0_ec_t, &
                           el_con_geo_t, macro_el_t
USE ph_freq_anharmonic, ONLY : celldmf_t, el_consf_t, el_compf_t, b0f_ec_t, &
                           el_conf_geo_t, macro_elf_t
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata, jdata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)
TYPE(poly1) :: ec_p1
TYPE(poly2) :: ec_p2
TYPE(poly3) :: ec_p3
TYPE(poly4) :: ec_p4

INTEGER :: i, j, idata, itemp, startt, lastt, pdelc_dos, pdelc_ph
INTEGER :: compute_nwork 

ibrav=ibrav_geo(1)
IF (lmurn) THEN
   nvar=1
ELSE
   nvar=crystal_parameters(ibrav)
ENDIF

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
   DO idata=1,ndata
      x(1,idata)=celldm_geo(1,idata)
   ENDDO
ELSE
   CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)
ENDIF

el_cons_t=0.0_DP
el_consf_t=0.0_DP
el_comp_t=0.0_DP
el_compf_t=0.0_DP
macro_el_t=0.0_DP
macro_elf_t=0.0_DP
b0_ec_t=0.0_DP
b0f_ec_t=0.0_DP
pdelc_dos=MIN(poly_degree_elc,fndos_ec-1)
pdelc_ph=MIN(poly_degree_elc,fnph_ec-1)
IF (pdelc_dos/=poly_degree_elc) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for dos")') pdelc_dos
IF (pdelc_ph/=poly_degree_elc) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for ph")') pdelc_ph

DO itemp=startt,lastt
   IF (itemp==1.OR.itemp==ntemp) CYCLE
   
   IF (el_cons_qha_geo_available.AND.ltherm_dos) THEN
      CALL compress_celldm(celldm_t(:,itemp),xfit,nvar,ibrav)
      f=0.0_DP
      DO i=1,6
         DO j=1,6
            IF (el_con_geo_t(i,j,itemp,f_geodos_ec)>0.1_DP) THEN
               WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")') i,j
               WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

               jdata=0
               DO idata=1,ndata
                  IF (found_dos_ec(idata)) THEN
                     jdata=jdata+1
                     x1(:,jdata)=x(:,idata)
                     f(jdata)=el_con_geo_t(i,j,itemp,idata)
                  ENDIF
               END DO

               IF (pdelc_dos==4) THEN
                  CALL init_poly(nvar,ec_p4)
                  CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,ec_p4) 
                  CALL evaluate_fit_quartic(nvar, xfit, el_cons_t(i,j,itemp),&
                                         ec_p4)
                  CALL clean_poly(ec_p4)
               ELSEIF (pdelc_dos==3) THEN
                  CALL init_poly(nvar,ec_p3)
                  CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
                  CALL evaluate_fit_cubic(nvar,xfit, el_cons_t(i,j,itemp), &
                       ec_p3)
                  CALL clean_poly(ec_p3)
               ELSEIF (pdelc_dos==2) THEN
                  CALL init_poly(nvar,ec_p2)
                  CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f, ec_p2)
                  CALL evaluate_fit_quadratic(nvar,xfit,el_cons_t(i,j,itemp), &
                                      ec_p2)
                  CALL clean_poly(ec_p2)
               ELSEIF (pdelc_dos==1) THEN
                  CALL init_poly(nvar,ec_p1)
                  CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
                  CALL evaluate_fit_linear(nvar, xfit, el_cons_t(i,j,itemp), &
                                                       ec_p1)
                  CALL clean_poly(ec_p1)
               ELSE
                  CALL errore('write_elastic_t_qha','wrong poly_degree_elc',1)
               ENDIF
            ENDIF
!            IF (i/=j) el_cons_t(j,i,itemp)=el_cons_t(i,j,itemp)
         ENDDO
      ENDDO
      CALL compute_elastic_compliances(el_cons_t(:,:,itemp), &
                                                      el_comp_t(:,:,itemp))
      CALL print_macro_elasticity(ibrav,el_cons_t(:,:,itemp), &
                         el_comp_t(:,:,itemp), macro_el_t(:,itemp), .FALSE.)
      b0_ec_t(itemp)=(macro_el_t(1,itemp) + macro_el_t(5,itemp)) * 0.5_DP
   END IF

   IF (el_consf_qha_geo_available.AND.ltherm_freq) THEN
      CALL compress_celldm(celldmf_t(:,itemp),xfit,nvar,ibrav)
      f=0.0_DP
      DO i=1,6
         DO j=1,6
            IF (el_conf_geo_t(i,j,itemp,f_geoph_ec)>0.1_DP) THEN
               WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")') i,j
               WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

               jdata=0
               DO idata=1,ndata
                  IF (found_ph_ec(idata)) THEN
                     jdata=jdata+1
                     x1(:,jdata)=x(:,idata)
                     f(jdata)=el_conf_geo_t(i,j,itemp,idata)
                  ENDIF
               END DO

               IF (pdelc_ph==4) THEN
                  CALL init_poly(nvar,ec_p4)
                  CALL fit_multi_quartic(jdata, nvar, lsolve, x1, f, ec_p4)
                  CALL evaluate_fit_quartic(nvar,xfit,el_consf_t(i,j,itemp), &
                                                                      ec_p4) 
                  CALL clean_poly(ec_p4)
               ELSEIF (pdelc_ph==3) THEN
                  CALL init_poly(nvar,ec_p3)
                  CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
                  CALL evaluate_fit_cubic(nvar,xfit, el_consf_t(i,j,itemp), &
                      ec_p3)
                  CALL clean_poly(ec_p3)
               ELSEIF (pdelc_ph==2) THEN
                  CALL init_poly(nvar,ec_p2)
                  CALL fit_multi_quadratic(jdata, nvar, lsolve, x1, f, ec_p2)
                  CALL evaluate_fit_quadratic(nvar,xfit, &
                                   el_consf_t(i,j,itemp), ec_p2)
                  CALL clean_poly(ec_p2)
               ELSEIF (pdelc_ph==1) THEN
                  CALL init_poly(nvar,ec_p1)
                  CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
                  CALL evaluate_fit_linear(nvar,xfit,el_consf_t(i,j,itemp), &
                                 ec_p1)
                  CALL clean_poly(ec_p1)
               ELSE
                  CALL errore('write_elastic_t_qha','wrong poly_degree_elc',1)
               ENDIF
            ENDIF
!            IF (i/=j) el_consf_t(j,i,itemp)=el_consf_t(i,j,itemp)
         ENDDO
      ENDDO
      CALL compute_elastic_compliances(el_consf_t(:,:,itemp), &
                                                      el_compf_t(:,:,itemp))
      CALL print_macro_elasticity(ibrav,el_consf_t(:,:,itemp), &
                          el_compf_t(:,:,itemp), macro_elf_t(:,itemp),.FALSE.)
      b0f_ec_t(itemp)=(macro_elf_t(1,itemp) + macro_elf_t(5,itemp))*0.5_DP
   END IF
ENDDO
!

IF (ltherm_dos) THEN
   CALL mp_sum(el_cons_t, world_comm)
   CALL mp_sum(el_comp_t, world_comm)
   CALL mp_sum(b0_ec_t, world_comm)
   CALL mp_sum(macro_el_t, world_comm)
   lelastic=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_cons_t, b0_ec_t, &
                                                       filelastic, 0)
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_comp_t, b0_ec_t, & 
                                                       filelastic, 1)
ENDIF
!
IF (ltherm_freq) THEN
   CALL mp_sum(el_consf_t, world_comm)
   CALL mp_sum(el_compf_t, world_comm)
   CALL mp_sum(b0f_ec_t, world_comm)
   CALL mp_sum(macro_elf_t, world_comm)
   lelasticf=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_ph'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_consf_t, b0f_ec_t, &
                                                          filelastic, 0)

   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp_ph'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_compf_t, b0f_ec_t, &
                                                           filelastic,1)
ENDIF

DEALLOCATE(x1)
DEALLOCATE(f)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_elastic_t_qha
!
! Copyright (C) 2023 Andrea Dal Corso
!
! These two routines are similar to the previous one, but make the
! interpolation for a few pressures, or for a few temperatures as a 
! function of pressure
!
!---------------------------------------------------------------------------
SUBROUTINE write_elastic_pt_qha()
!---------------------------------------------------------------------------
!
!  This routine interpolates the temperature dependent quasi harmonic 
!  elastic constants at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at temperature T for several pressures.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE thermo_sym, ONLY : laue
USE control_pressure, ONLY : npress_plot, ipress_plot, press
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE control_mur, ONLY : lmurn
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_elastic_constants, ONLY : el_cons_qha_geo_available, &
                             el_consf_qha_geo_available, lelastic_pt, &
                             lelasticf, lelasticf_pt, &
                             found_ph_ec, found_dos_ec, fnph_ec, fndos_ec, &
                             f_geodos_ec, f_geoph_ec
USE lattices,       ONLY : compress_celldm
USE elastic_constants, ONLY : el_con, el_compliances, write_el_cons_on_file, &
                              compute_elastic_compliances, &
                              print_macro_elasticity
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE anharmonic,     ONLY : el_con_geo_t
USE ph_freq_anharmonic,     ONLY : el_conf_geo_t
USE anharmonic_pt,  ONLY : celldm_pt, el_cons_pt, el_comp_pt, b0_ec_pt, &
                           macro_el_pt
USE ph_freq_anharmonic_pt,  ONLY : celldmf_pt, el_consf_pt, el_compf_pt, & 
                                   b0f_ec_pt, macro_elf_pt
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata, jdata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)
TYPE(poly1) :: ec_p1
TYPE(poly2) :: ec_p2
TYPE(poly3) :: ec_p3
TYPE(poly4) :: ec_p4

INTEGER :: i, j, idata, itemp, ipressp, ipress, startt, lastt, &
                 pdelc_dos, pdelc_ph
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
ALLOCATE(x1(nvar,ndata))
!
!  Part 1 evaluation of the polynomial coefficients
!
ALLOCATE(xfit(nvar))
CALL divide(world_comm, ntemp, startt, lastt)

IF (lmurn) THEN
   DO idata=1,ndata
      x(1,idata)=celldm_geo(1,idata)
   ENDDO
ELSE
   CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)
ENDIF

el_cons_pt=0.0_DP
el_comp_pt=0.0_DP
macro_el_pt=0.0_DP
b0_ec_pt=0.0_DP

el_consf_pt=0.0_DP
el_compf_pt=0.0_DP
macro_elf_pt=0.0_DP
b0f_ec_pt=0.0_DP

pdelc_dos=MIN(poly_degree_elc,fndos_ec-1)
pdelc_ph=MIN(poly_degree_elc,fnph_ec-1)
IF (pdelc_dos/=poly_degree_elc.AND.ltherm_dos) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for dos")') pdelc_dos
IF (pdelc_ph/=poly_degree_elc.AND.ltherm_freq) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for ph")') pdelc_ph

DO ipressp=1,npress_plot
   ipress=ipress_plot(ipressp)
   IF (ltherm_dos) THEN
      DO itemp=startt,lastt
         IF (itemp==1.OR.itemp==ntemp) CYCLE
   
         IF (el_cons_qha_geo_available) THEN
            CALL compress_celldm(celldm_pt(:,itemp,ipressp),xfit,nvar,ibrav)
            f=0.0_DP
            DO i=1,6
               DO j=1,6
                  IF (el_con_geo_t(i,j,itemp,f_geodos_ec)>0.1_DP) THEN
                     WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")') i,j
                     WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

                     jdata=0
                     DO idata=1,ndata
                        IF (found_dos_ec(idata)) THEN
                           jdata=jdata+1
                           x1(:,jdata)=x(:,idata)
                           f(jdata)=el_con_geo_t(i,j,itemp,idata)
                        ENDIF
                     END DO

                     IF (pdelc_dos==4) THEN
                        CALL init_poly(nvar,ec_p4)
                        CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,ec_p4) 
                        CALL evaluate_fit_quartic(nvar, xfit, &
                                      el_cons_pt(i,j,itemp,ipressp),&
                                         ec_p4)
                        CALL clean_poly(ec_p4)
                     ELSEIF (pdelc_dos==3) THEN
                        CALL init_poly(nvar,ec_p3)
                        CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
                        CALL evaluate_fit_cubic(nvar,xfit, &
                          el_cons_pt(i,j,itemp,ipressp), ec_p3)
                        CALL clean_poly(ec_p3)
                     ELSEIF (pdelc_dos==2) THEN
                        CALL init_poly(nvar,ec_p2)
                        CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f, ec_p2)
                        CALL evaluate_fit_quadratic(nvar,xfit, &
                                       el_cons_pt(i,j,itemp,ipressp), ec_p2)
                        CALL clean_poly(ec_p2)
                     ELSEIF (pdelc_dos==1) THEN
                        CALL init_poly(nvar,ec_p1)
                        CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
                        CALL evaluate_fit_linear(nvar, xfit, &
                            el_cons_pt(i,j,itemp,ipressp), ec_p1)
                        CALL clean_poly(ec_p1)
                     ELSE
                        CALL errore('write_elastic_pt_qha',&
                                                 'wrong poly_degree_elc',1)
                     ENDIF
                  ENDIF
!                  IF (i/=j) el_cons_pt(j,i,itemp,ipressp)=&
!                                               el_cons_pt(i,j,itemp,ipressp)
               ENDDO
            ENDDO
         ENDIF
         CALL compute_elastic_compliances(el_cons_pt(:,:,itemp,ipressp),       &
                                              el_comp_pt(:,:,itemp,ipressp))
         CALL print_macro_elasticity(ibrav,el_cons_pt(:,:,itemp,ipressp),      &
               el_comp_pt(:,:,itemp,ipressp), macro_el_pt(:,itemp,ipressp), &
                                                                     .FALSE.)
         b0_ec_pt(itemp,ipressp)=(macro_el_pt(1,itemp,ipressp) + &
                               macro_el_pt(5,itemp,ipressp) ) * 0.5_DP
      ENDDO

      CALL mp_sum(el_cons_pt(:,:,:,ipressp), world_comm)
      CALL mp_sum(el_comp_pt(:,:,:,ipressp), world_comm)
      CALL mp_sum(b0_ec_pt(:,ipressp), world_comm)
      CALL mp_sum(macro_el_pt(:,:,ipressp), world_comm)
      lelastic_pt=.TRUE.
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_press'
      CALL add_value(filelastic, press(ipress))
      CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, &
              el_cons_pt(1,1,1,ipressp), b0_ec_pt(1,ipressp), filelastic, 0)
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp_press'
      CALL add_value(filelastic, press(ipress))
      CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, &
              el_comp_pt(1,1,1,ipressp), b0_ec_pt(1,ipressp), filelastic, 1)
   ENDIF

   IF (ltherm_freq) THEN
      DO itemp=startt,lastt
         IF (itemp==1.OR.itemp==ntemp) CYCLE
   
         IF (el_consf_qha_geo_available) THEN
            CALL compress_celldm(celldmf_pt(:,itemp,ipressp),xfit,nvar,ibrav)
            f=0.0_DP
            DO i=1,6
               DO j=1,6
                  IF (el_conf_geo_t(i,j,itemp,f_geoph_ec)>0.1_DP) THEN
                     WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")') i,j
                     WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

                     jdata=0
                     DO idata=1,ndata
                        IF (found_ph_ec(idata)) THEN
                           jdata=jdata+1
                           x1(:,jdata)=x(:,idata)
                           f(jdata)=el_conf_geo_t(i,j,itemp,idata)
                        ENDIF
                     END DO

                     IF (pdelc_ph==4) THEN
                        CALL init_poly(nvar,ec_p4)
                        CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,ec_p4) 
                        CALL evaluate_fit_quartic(nvar, xfit, &
                                      el_consf_pt(i,j,itemp,ipressp),&
                                         ec_p4)
                        CALL clean_poly(ec_p4)
                     ELSEIF (pdelc_ph==3) THEN
                        CALL init_poly(nvar,ec_p3)
                        CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
                        CALL evaluate_fit_cubic(nvar,xfit, &
                          el_consf_pt(i,j,itemp,ipressp), ec_p3)
                        CALL clean_poly(ec_p3)
                     ELSEIF (pdelc_ph==2) THEN
                        CALL init_poly(nvar,ec_p2)
                        CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f, ec_p2)
                        CALL evaluate_fit_quadratic(nvar,xfit, &
                                    el_consf_pt(i,j,itemp,ipressp), ec_p2)
                        CALL clean_poly(ec_p2)
                     ELSEIF (pdelc_ph==1) THEN
                        CALL init_poly(nvar,ec_p1)
                        CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
                        CALL evaluate_fit_linear(nvar, xfit, &
                            el_consf_pt(i,j,itemp,ipressp), ec_p1)
                        CALL clean_poly(ec_p1)
                     ELSE
                        CALL errore('write_elastic_pt_qha',&
                                                 'wrong poly_degree_elc',1)
                     ENDIF
                  ENDIF
!                  IF (i/=j) el_consf_pt(j,i,itemp,ipressp)=&
!                                               el_consf_pt(i,j,itemp,ipressp)
               ENDDO
            ENDDO
         ENDIF
         CALL compute_elastic_compliances(el_consf_pt(:,:,itemp,ipressp), &
                                           el_compf_pt(:,:,itemp,ipressp))
         CALL print_macro_elasticity(ibrav,el_consf_pt(:,:,itemp,ipressp),&
               el_compf_pt(:,:,itemp,ipressp), macro_elf_pt(:,itemp,ipressp), &
                                                                     .FALSE.)
         b0f_ec_pt(itemp,ipressp)=(macro_elf_pt(1,itemp,ipressp) + &
                                macro_elf_pt(5,itemp,ipressp) ) * 0.5_DP
      ENDDO
      CALL mp_sum(el_consf_pt(:,:,:,ipressp), world_comm)
      CALL mp_sum(el_compf_pt(:,:,:,ipressp), world_comm)
      CALL mp_sum(b0f_ec_pt(:,ipressp), world_comm)
      CALL mp_sum(macro_elf_pt(:,:,ipressp), world_comm)
      lelasticf_pt=.TRUE.
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_ph_press'
      CALL add_value(filelastic, press(ipress))
      CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, &
              el_consf_pt(1,1,1,ipressp), b0f_ec_pt(1,ipressp), filelastic, 0)
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp_ph_press'
      CALL add_value(filelastic, press(ipress))
      CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, &
              el_compf_pt(1,1,1,ipressp), b0f_ec_pt(1,ipressp), filelastic, 1)
   ENDIF
ENDDO

DEALLOCATE(x1)
DEALLOCATE(f)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_elastic_pt_qha

!---------------------------------------------------------------------------
SUBROUTINE write_elastic_ptt_qha()
!---------------------------------------------------------------------------
!
!  This routine interpolates the temperature dependent quasi harmonic 
!  elastic constants at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at pressure p for several temperatures.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE thermo_sym, ONLY : laue
USE control_mur, ONLY : lmurn
USE control_pressure, ONLY : npress, press
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_elastic_constants, ONLY : el_cons_qha_geo_available, &
                             el_consf_qha_geo_available, lelastic_ptt, &
                             lelasticf, lelasticf_ptt, &
                             found_ph_ec, found_dos_ec, fnph_ec, fndos_ec, &
                             f_geodos_ec, f_geoph_ec
USE lattices,       ONLY : compress_celldm
USE elastic_constants, ONLY : el_con, el_compliances, write_el_cons_on_file, &
                              compute_elastic_compliances, &
                              print_macro_elasticity
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE anharmonic,     ONLY : el_con_geo_t
USE ph_freq_anharmonic,     ONLY : el_conf_geo_t
USE anharmonic_ptt,  ONLY : celldm_ptt, el_cons_ptt, el_comp_ptt, b0_ec_ptt, &
                            macro_el_ptt
USE ph_freq_anharmonic_ptt,  ONLY : celldmf_ptt, el_consf_ptt, el_compf_ptt, &
                                    b0f_ec_ptt, macro_elf_ptt
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp, ntemp_plot, itemp_plot
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata, jdata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)
TYPE(poly1) :: ec_p1
TYPE(poly2) :: ec_p2
TYPE(poly3) :: ec_p3
TYPE(poly4) :: ec_p4

INTEGER :: i, j, idata, itemp, itempp, ipress, startp, lastp, &
                 pdelc_dos, pdelc_ph
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
ALLOCATE(x1(nvar,ndata))
!
ALLOCATE(xfit(nvar))
CALL divide(world_comm, npress, startp, lastp)

IF (lmurn) THEN
   DO idata=1,ndata
      x(1,idata)=celldm_geo(1,idata)
   ENDDO
ELSE
   CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)
ENDIF

el_cons_ptt=0.0_DP
el_comp_ptt=0.0_DP
macro_el_ptt=0.0_DP
b0_ec_ptt=0.0_DP

el_consf_ptt=0.0_DP
el_compf_ptt=0.0_DP
macro_elf_ptt=0.0_DP
b0f_ec_ptt=0.0_DP

pdelc_dos=MIN(poly_degree_elc,fndos_ec-1)
pdelc_ph=MIN(poly_degree_elc,fnph_ec-1)
IF (pdelc_dos/=poly_degree_elc.AND.ltherm_dos) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for dos")') pdelc_dos
IF (pdelc_ph/=poly_degree_elc.AND.ltherm_freq) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for ph")') pdelc_ph

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   IF (ltherm_dos) THEN
      DO ipress=startp,lastp
         IF (el_cons_qha_geo_available) THEN
            CALL compress_celldm(celldm_ptt(:,ipress,itempp),xfit,nvar,ibrav)
            f=0.0_DP
            DO i=1,6
               DO j=1,6
                  IF (el_con_geo_t(i,j,itemp,f_geodos_ec)>0.1_DP) THEN
                     jdata=0
                     DO idata=1,ndata
                        IF (found_dos_ec(idata)) THEN
                           jdata=jdata+1
                           x1(:,jdata)=x(:,idata)
                           f(jdata)=el_con_geo_t(i,j,itemp,idata)
                        ENDIF
                     END DO

                     IF (pdelc_dos==4) THEN
                        CALL init_poly(nvar,ec_p4)
                        CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,ec_p4) 
                        CALL evaluate_fit_quartic(nvar, xfit, &
                               el_cons_ptt(i,j,ipress,itempp), ec_p4)
                        CALL clean_poly(ec_p4)
                     ELSEIF (pdelc_dos==3) THEN
                        CALL init_poly(nvar,ec_p3)
                        CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
                        CALL evaluate_fit_cubic(nvar,xfit, &
                            el_cons_ptt(i,j,ipress,itempp), ec_p3)
                        CALL clean_poly(ec_p3)
                     ELSEIF (pdelc_dos==2) THEN
                        CALL init_poly(nvar,ec_p2)
                        CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f, ec_p2)
                        CALL evaluate_fit_quadratic(nvar,xfit, &
                                  el_cons_ptt(i,j,ipress,itempp), ec_p2)
                        CALL clean_poly(ec_p2)
                     ELSEIF (pdelc_dos==1) THEN
                        CALL init_poly(nvar,ec_p1)
                        CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
                        CALL evaluate_fit_linear(nvar, xfit, &
                            el_cons_ptt(i,j,ipress,itempp), ec_p1)
                        CALL clean_poly(ec_p1)
                     ELSE
                        CALL errore('write_elastic_ptt_qha',&
                                                 'wrong poly_degree_elc',1)
                     ENDIF
                  ENDIF
!                  IF (i/=j) el_cons_ptt(j,i,ipress,itempp)=&
!                                           el_cons_ptt(i,j,ipress,itempp)
               ENDDO
            ENDDO
         ENDIF
         CALL compute_elastic_compliances(el_cons_ptt(:,:,ipress,itempp),     &
                                          el_comp_ptt(:,:,ipress,itempp))
         CALL print_macro_elasticity(ibrav,el_cons_ptt(:,:,ipress,itempp),   &
           el_comp_ptt(:,:,ipress,itempp), macro_el_ptt(:,ipress,itempp),  &
                                                                     .FALSE.)
         b0_ec_ptt(ipress,itempp)=(macro_el_ptt(1,ipress,itempp) + &
                                macro_el_ptt(5,ipress,itempp) ) * 0.5_DP
      ENDDO
!
      CALL mp_sum(el_cons_ptt(:,:,:,itempp), world_comm)
      CALL mp_sum(el_comp_ptt(:,:,:,itempp), world_comm)
      CALL mp_sum(b0_ec_ptt(:,itempp), world_comm)
      CALL mp_sum(macro_el_ptt(:,:,itempp), world_comm)
      lelastic_ptt=.TRUE.
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_temp'
      CALL add_value(filelastic, temp(itemp))
      CALL write_el_cons_on_file(press, npress, ibrav, laue, &
              el_cons_ptt(1,1,1,itempp), b0_ec_ptt(1,itempp), filelastic, 2)
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp_temp'
      CALL add_value(filelastic, temp(itemp))
      CALL write_el_cons_on_file(press, npress, ibrav, laue, &
              el_comp_ptt(1,1,1,itempp), b0_ec_ptt(1,itempp), filelastic, 3)
   ENDIF

   IF (ltherm_freq) THEN
      DO ipress=startp,lastp
         IF (el_consf_qha_geo_available) THEN
            CALL compress_celldm(celldmf_ptt(:,ipress,itempp),xfit,nvar,ibrav)
            f=0.0_DP
            DO i=1,6
               DO j=1,6
                  IF (el_conf_geo_t(i,j,itemp,f_geoph_ec)>0.1_DP) THEN
                     jdata=0
                     DO idata=1,ndata
                        IF (found_ph_ec(idata)) THEN
                           jdata=jdata+1
                           x1(:,jdata)=x(:,idata)
                           f(jdata)=el_conf_geo_t(i,j,itemp,idata)
                        ENDIF
                     END DO

                     IF (pdelc_ph==4) THEN
                        CALL init_poly(nvar,ec_p4)
                        CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,ec_p4) 
                        CALL evaluate_fit_quartic(nvar, xfit, &
                               el_consf_ptt(i,j,ipress,itempp), ec_p4)
                        CALL clean_poly(ec_p4)
                     ELSEIF (pdelc_ph==3) THEN
                        CALL init_poly(nvar,ec_p3)
                        CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
                        CALL evaluate_fit_cubic(nvar,xfit, &
                            el_consf_ptt(i,j,ipress,itempp), ec_p3)
                        CALL clean_poly(ec_p3)
                     ELSEIF (pdelc_ph==2) THEN
                        CALL init_poly(nvar,ec_p2)
                        CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f, ec_p2)
                        CALL evaluate_fit_quadratic(nvar,xfit, &
                                  el_consf_ptt(i,j,ipress,itempp), ec_p2)
                        CALL clean_poly(ec_p2)
                     ELSEIF (pdelc_dos==1) THEN
                        CALL init_poly(nvar,ec_p1)
                        CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
                        CALL evaluate_fit_linear(nvar, xfit, &
                            el_consf_ptt(i,j,ipress,itempp), ec_p1)
                        CALL clean_poly(ec_p1)
                     ELSE
                        CALL errore('write_elastic_ptt_qha',&
                                                 'wrong poly_degree_elc',1)
                     ENDIF
                  ENDIF
!                  IF (i/=j) el_consf_ptt(j,i,ipress,itempp)=&
!                                           el_consf_ptt(i,j,ipress,itempp)
               ENDDO
            ENDDO
         ENDIF
         CALL compute_elastic_compliances(el_consf_ptt(:,:,ipress,itempp),     &
                                          el_compf_ptt(:,:,ipress,itempp))
         CALL print_macro_elasticity(ibrav,el_consf_ptt(:,:,ipress,itempp),   &
           el_compf_ptt(:,:,ipress,itempp), macro_elf_ptt(:,ipress,itempp),  &
                                                                     .FALSE.)
         b0f_ec_ptt(ipress,itempp)=(macro_elf_ptt(1,ipress,itempp)+ &
                                 macro_elf_ptt(5,ipress,itempp) ) * 0.5_DP
      ENDDO
!
      CALL mp_sum(el_consf_ptt(:,:,:,itempp), world_comm)
      CALL mp_sum(el_compf_ptt(:,:,:,itempp), world_comm)
      CALL mp_sum(b0f_ec_ptt(:,itempp), world_comm)
      CALL mp_sum(macro_elf_ptt(:,:,itempp), world_comm)
      lelasticf_ptt=.TRUE.
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_ph_temp'
      CALL add_value(filelastic, temp(itemp))
      CALL write_el_cons_on_file(press, npress, ibrav, laue, &
              el_consf_ptt(1,1,1,itempp), b0f_ec_ptt(1,itempp), filelastic, 2)
      filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp_ph_temp'
      CALL add_value(filelastic, temp(itemp))
      CALL write_el_cons_on_file(press, npress, ibrav, laue, &
              el_compf_ptt(1,1,1,itempp), b0f_ec_ptt(1,itempp), filelastic, 3)
   ENDIF

ENDDO

DEALLOCATE(x1)
DEALLOCATE(f)
DEALLOCATE(x)

RETURN
END SUBROUTINE write_elastic_ptt_qha
!
!---------------------------------------------------------------------------
SUBROUTINE write_dyde_t_qha(istep)
!---------------------------------------------------------------------------
!
!  This routine interpolates the temperature dependent derivative of
!  the crystal parameters with respect to strain at the crystal parameters 
!  that minimize the Helmholtz (or Gibbs) free energy at temperature T.
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_elastic_constants, ONLY : el_cons_qha_geo_available, &
                             el_consf_qha_geo_available, lelastic, lelasticf, &
                             found_ph_ec, found_dos_ec, fnph_ec, fndos_ec, &
                             f_geodos_ec, f_geoph_ec, dyde
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE anharmonic,     ONLY : celldm_t, dyde_t, vmin_t
USE ph_freq_anharmonic, ONLY : celldmf_t, dydef_t, vminf_t
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
INTEGER :: istep
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata, jdata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:), x1(:,:)
CHARACTER(LEN=6) :: int_to_char
TYPE(poly1) :: ec_p1
TYPE(poly2) :: ec_p2
TYPE(poly3) :: ec_p3
TYPE(poly4) :: ec_p4

INTEGER :: i, j, idata, itemp, startt, lastt, pdelc_dos, pdelc_ph
INTEGER :: compute_nwork 

ibrav=ibrav_geo(1)
nvar=crystal_parameters(ibrav)
ndata=compute_nwork()

ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
ALLOCATE(x1(nvar,ndata))
!
!  Part 1 evaluation of the polynomial coefficients
!
ALLOCATE(xfit(nvar))
CALL divide(world_comm, ntemp, startt, lastt)

CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)

dyde_t(istep,:)=0.0_DP
dydef_t(istep,:)=0.0_DP
pdelc_dos=MIN(poly_degree_elc,fndos_ec-1)
pdelc_ph=MIN(poly_degree_elc,fnph_ec-1)
IF (pdelc_dos/=poly_degree_elc) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for dos")') pdelc_dos
IF (pdelc_ph/=poly_degree_elc) &
   WRITE(stdout,'(5x,"Poly_degree_elc decreased to ",i3," for ph")') pdelc_ph

DO itemp=startt,lastt
   IF (itemp==1.OR.itemp==ntemp) CYCLE
   
   IF (el_cons_qha_geo_available.AND.ltherm_dos) THEN
      CALL compress_celldm(celldm_t(:,itemp),xfit,nvar,ibrav)
      f=0.0_DP
      WRITE(stdout,'(/,5x,"Fitting dos dyde for step", i5)') istep
      WRITE(stdout,'(/,5x,"at Temperature",f15.5," K")') temp(itemp)
      jdata=0
      DO idata=1,ndata
         IF (found_dos_ec(idata)) THEN
            jdata=jdata+1
            x1(:,jdata)=x(:,idata)
            f(jdata)=dyde(istep,idata,itemp)
            IF (itemp==5) WRITE(stdout,*) idata, jdata, x1(:,jdata), &
                                          dyde(istep,idata,itemp)
         ENDIF
      END DO

      IF (pdelc_dos==4) THEN
         CALL init_poly(nvar,ec_p4)
         CALL fit_multi_quartic(jdata,nvar,lsolve,x1,f,ec_p4) 
         CALL evaluate_fit_quartic(nvar, xfit, dyde_t(istep,itemp),ec_p4)
         CALL clean_poly(ec_p4)
      ELSEIF (pdelc_dos==3) THEN
         CALL init_poly(nvar,ec_p3)
         CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
         CALL evaluate_fit_cubic(nvar,xfit, dyde_t(istep,itemp),ec_p3)
         CALL clean_poly(ec_p3)
      ELSEIF (pdelc_dos==2) THEN
         CALL init_poly(nvar,ec_p2)
         CALL fit_multi_quadratic(jdata,nvar,lsolve,x1,f,ec_p2)
         CALL evaluate_fit_quadratic(nvar,xfit,dyde_t(istep,itemp),ec_p2)
         CALL clean_poly(ec_p2)
      ELSEIF (pdelc_dos==1) THEN
         CALL init_poly(nvar,ec_p1)
         CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
         CALL evaluate_fit_linear(nvar,xfit,dyde_t(istep,itemp),ec_p1)
         CALL clean_poly(ec_p1)
      ELSE
         CALL errore('write_elastic_t_qha','wrong poly_degree_elc',1)
      ENDIF
   ENDIF

   IF (el_consf_qha_geo_available.AND.ltherm_freq) THEN
      CALL compress_celldm(celldmf_t(:,itemp),xfit,nvar,ibrav)
      f=0.0_DP
      WRITE(stdout,'(/,5x,"Fitting freq dyde for step",i5)') istep
      WRITE(stdout,'(/,5x,"at Temperature",f15.5," K")') temp(itemp)
      jdata=0
      DO idata=1,ndata
         IF (found_ph_ec(idata)) THEN
            jdata=jdata+1
            x1(:,jdata)=x(:,idata)
            f(jdata)=dyde(istep,idata,itemp)
         ENDIF
      ENDDO

      IF (pdelc_ph==4) THEN
         CALL init_poly(nvar,ec_p4)
         CALL fit_multi_quartic(jdata, nvar, lsolve, x1, f, ec_p4)
         CALL evaluate_fit_quartic(nvar,xfit,dydef_t(istep,itemp), ec_p4) 
         CALL clean_poly(ec_p4)
      ELSEIF (pdelc_ph==3) THEN
         CALL init_poly(nvar,ec_p3)
         CALL fit_multi_cubic(jdata,nvar,lsolve,x1,f,ec_p3)
         CALL evaluate_fit_cubic(nvar,xfit, dydef_t(istep,itemp), ec_p3)
         CALL clean_poly(ec_p3)
      ELSEIF (pdelc_ph==2) THEN
         CALL init_poly(nvar,ec_p2)
         CALL fit_multi_quadratic(jdata, nvar, lsolve, x1, f, ec_p2)
         CALL evaluate_fit_quadratic(nvar,xfit, dydef_t(istep,itemp), ec_p2)
         CALL clean_poly(ec_p2)
      ELSEIF (pdelc_ph==1) THEN
         CALL init_poly(nvar,ec_p1)
         CALL fit_multi_linear(jdata,nvar,lsolve,x1,f,ec_p1)
         CALL evaluate_fit_linear(nvar,xfit,dydef_t(istep,itemp), ec_p1)
         CALL clean_poly(ec_p1)
      ELSE
         CALL errore('write_elastic_t_qha','wrong poly_degree_elc',1)
      ENDIF
   ENDIF
ENDDO
!
IF (ltherm_dos) THEN
   CALL mp_sum(dyde_t(istep,:), world_comm)
   lelastic=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.int_rel.'//int_to_char(istep)
   CALL write_dyde_on_file(temp, ntemp, dyde_t, vmin_t, istep, filelastic)
ENDIF
!
IF (ltherm_freq) THEN
   CALL mp_sum(dydef_t(istep,:), world_comm)
   lelasticf=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.int_relf.'//int_to_char(istep)
   CALL write_dyde_on_file(temp, ntemp, dydef_t, vminf_t, istep, filelastic)
ENDIF

DEALLOCATE(f)
DEALLOCATE(x)
DEALLOCATE(x1)

RETURN
END SUBROUTINE write_dyde_t_qha
!
!-------------------------------------------------------------------------
SUBROUTINE write_dyde_on_file(temp, ntemp, dyde, omega, istep, filename)
!-------------------------------------------------------------------------
!
USE kinds, ONLY : DP
USE io_global,      ONLY : meta_ionode

IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, istep
REAL(DP), INTENT(IN) :: temp(ntemp), dyde(21,ntemp), omega(ntemp)
CHARACTER(LEN=256), INTENT(IN) :: filename

INTEGER :: iu_therm, itemp
INTEGER :: find_free_unit

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')

   WRITE(iu_therm,'("# Temp (K)        dyde (a.u.)      volume (a.u.)**3")')
   DO itemp = 2, ntemp-1
      WRITE(iu_therm, '(e12.5,2e22.13)') temp(itemp), dyde(istep,itemp), &
                                         omega(itemp)
   ENDDO
   CLOSE(iu_therm, STATUS='KEEP')
ENDIF

RETURN
END SUBROUTINE write_dyde_on_file
!
