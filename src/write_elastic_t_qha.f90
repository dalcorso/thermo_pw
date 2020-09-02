!
! Copyright (C) 2019 C. Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_elastic_t_qha()
!
!  This routine interpolate the temperature dependent quasi harmonic 
!  elastic constants at the crystal parameters that minimize the Helmholtz 
!  (or Gibbs) free energy at temperature T.
!
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : stdout
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE linear_surfaces,    ONLY : fit_multi_linear, evaluate_fit_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic
USE control_elastic_constants, ONLY : el_cons_qha_geo_available, &
                             el_consf_qha_geo_available  
USE lattices,       ONLY : compress_celldm
USE elastic_constants, ONLY : el_con, el_compliances, write_el_cons_on_file
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE anharmonic,     ONLY : celldm_t, el_cons_t, el_comp_t, b0_t, lelastic, &
                           el_con_geo_t
USE ph_freq_anharmonic, ONLY : celldmf_t, el_consf_t, el_compf_t, b0f_t, &
                           lelasticf, el_conf_geo_t
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
USE mp,          ONLY : mp_sum
USE mp_world,    ONLY : world_comm

IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ndata
REAL(DP), ALLOCATABLE :: x(:,:), f(:), xfit(:)
TYPE(poly1) :: ec_p1
TYPE(poly2) :: ec_p2
TYPE(poly3) :: ec_p3
TYPE(poly4) :: ec_p4

INTEGER :: i, j, idata, itemp, startt, lastt
INTEGER :: compute_nwork 

ibrav=ibrav_geo(1)
nvar=crystal_parameters(ibrav)
ndata=compute_nwork()

ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))
!
!  Part 1 evaluation of the polynomial coefficients
!
ALLOCATE(xfit(nvar))
CALL divide(world_comm, ntemp, startt, lastt)

CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)

el_cons_t=0.0_DP
el_consf_t=0.0_DP

DO itemp=startt,lastt
   IF (itemp==1.OR.itemp==ntemp) CYCLE
   
   IF (el_cons_qha_geo_available) THEN
      CALL compress_celldm(celldm_t(:,itemp),xfit,nvar,ibrav)
      f=0.0_DP
      DO i=1,6
         DO j=i,6
            IF (el_con_geo_t(i,j,itemp,1)>0.1_DP) THEN
               WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")') i,j
               WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

               DO idata=1,ndata
                  f(idata)=el_con_geo_t(i,j,itemp,idata)
               END DO

               IF (poly_degree_elc==4) THEN
                  CALL init_poly(nvar,ec_p4)
                  CALL fit_multi_quartic(ndata,nvar,lsolve,x,f,ec_p4) 
                  CALL evaluate_fit_quartic(nvar, xfit, el_cons_t(i,j,itemp),&
                                         ec_p4)
                  CALL clean_poly(ec_p4)
               ELSEIF (poly_degree_elc==3) THEN
                  CALL init_poly(nvar,ec_p3)
                  CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,ec_p3)
                  CALL evaluate_fit_cubic(nvar,xfit, el_cons_t(i,j,itemp), &
                       ec_p3)
                  CALL clean_poly(ec_p3)
               ELSEIF (poly_degree_elc==2) THEN
                  CALL init_poly(nvar,ec_p2)
                  CALL fit_multi_quadratic(ndata,nvar,lsolve,x,f, ec_p2)
                  CALL evaluate_fit_quadratic(nvar,xfit,el_cons_t(i,j,itemp), &
                                      ec_p2)
                  CALL clean_poly(ec_p2)
               ELSEIF (poly_degree_elc==1) THEN
                  CALL init_poly(nvar,ec_p1)
                  CALL fit_multi_linear(ndata,nvar,lsolve,x,f,ec_p1)
                  CALL evaluate_fit_linear(nvar, xfit, el_cons_t(i,j,itemp), &
                                                       ec_p1)
                  CALL clean_poly(ec_p1)
               ELSE
                  CALL errore('write_elastic_t_qha','wrong poly_degree_elc',1)
               ENDIF
            ENDIF
            IF (i/=j) el_cons_t(j,i,itemp)=el_cons_t(i,j,itemp)
         ENDDO
      ENDDO
   END IF

   IF (el_consf_qha_geo_available) THEN
      CALL compress_celldm(celldmf_t(:,itemp),xfit,nvar,ibrav)
      f=0.0_DP
      DO i=1,6
         DO j=i,6
            IF (el_conf_geo_t(i,j,itemp,1)>0.1_DP) THEN
               WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")') i,j
               WRITE(stdout,'(/,5x,"at Temperature",f15.5,"K")') temp(itemp)

               DO idata=1,ndata
                  f(idata)=el_conf_geo_t(i,j,itemp,idata)
               END DO

               IF (poly_degree_elc==4) THEN
                  CALL init_poly(nvar,ec_p4)
                  CALL fit_multi_quartic(ndata,nvar, lsolve, x, f, ec_p4)
                  CALL evaluate_fit_quartic(nvar,xfit,el_consf_t(i,j,itemp), &
                                                                      ec_p4) 
                  CALL clean_poly(ec_p4)
               ELSEIF (poly_degree_elc==3) THEN
                  CALL init_poly(nvar,ec_p3)
                  CALL fit_multi_cubic(ndata,nvar,lsolve,x,f,ec_p3)
                  CALL evaluate_fit_cubic(nvar,xfit, el_consf_t(i,j,itemp), &
                      ec_p3)
                  CALL clean_poly(ec_p3)
               ELSEIF (poly_degree_elc==2) THEN
                  CALL init_poly(nvar,ec_p2)
                  CALL fit_multi_quadratic(ndata, nvar, lsolve, x, f, ec_p2)
                  CALL evaluate_fit_quadratic(nvar,xfit, &
                                   el_consf_t(i,j,itemp), ec_p2)
                  CALL clean_poly(ec_p2)
               ELSEIF (poly_degree_elc==1) THEN
                  CALL init_poly(nvar,ec_p1)
                  CALL fit_multi_linear(ndata,nvar,lsolve,x,f,ec_p1)
                  CALL evaluate_fit_linear(nvar,xfit,el_consf_t(i,j,itemp), &
                                 ec_p1)
                  CALL clean_poly(ec_p1)
               ELSE
                  CALL errore('write_elastic_t_qha','wrong poly_degree_elc',1)
               ENDIF
            ENDIF
            IF (i/=j) el_consf_t(j,i,itemp)=el_consf_t(i,j,itemp)
         ENDDO
      ENDDO
   END IF
ENDDO
!
IF (ltherm_dos) THEN
   CALL mp_sum(el_cons_t, world_comm)
   CALL compute_el_comp_t(el_cons_t,el_comp_t,b0_t)
   lelastic=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_cons_t, b0_t, &
                                                       filelastic, 0)
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_comp_t, b0_t, & 
                                                       filelastic, 1)
ENDIF

IF (ltherm_freq) THEN
   CALL mp_sum(el_consf_t, world_comm)
   CALL compute_el_comp_t(el_consf_t,el_compf_t,b0f_t)
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

RETURN
END SUBROUTINE write_elastic_t_qha
