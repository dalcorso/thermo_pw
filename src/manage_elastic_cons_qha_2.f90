!
! Copyright (C) 2019 Cristiano Malica and Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_elastic_cons_qha_2()

USE kinds,             ONLY : DP
USE cell_base,      ONLY : ibrav
USE initial_conf,      ONLY : ibrav_save
USE thermo_mod,     ONLY : energy_geo, tot_ngeo
USE thermo_sym,        ONLY : laue
USE elastic_constants, ONLY : epsilon_geo, el_con,           &
                              compute_elastic_constants_ene, &
                              write_el_cons_on_file
USE control_elastic_constants, ONLY : ngeo_strain, elcpvar, ngeom, &
                           work_base, el_con_omega_geo, epsil_geo
USE anharmonic,        ONLY : el_cons_t, el_comp_t, b0_t, lelastic
USE control_quartic_energy, ONLY : poly_degree_elc
USE linear_surfaces, ONLY : evaluate_fit_linear
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic
USE cubic_surfaces,   ONLY : evaluate_fit_cubic
USE quartic_surfaces, ONLY : evaluate_fit_quartic
USE temperature,    ONLY : ntemp, temp
USE lattices,       ONLY : expand_celldm, crystal_parameters
USE polynomial,     ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE data_files,     ONLY : flanhar
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE

INTEGER :: itemp, startt, lastt, igroup, ngroup, istrain, ndata, ind, nvar
CHARACTER(LEN=256)  :: filename, filelastic
REAL(DP), ALLOCATABLE :: free_ener(:), epsilon_geo_loc(:,:,:), x_t(:,:)
REAL(DP), ALLOCATABLE :: y(:)
REAL(DP) :: celldm_t(6), vmin_t(ntemp), compute_omega_geo

TYPE(poly1), ALLOCATABLE :: p1(:,:)
TYPE(poly2), ALLOCATABLE :: p2(:,:)
TYPE(poly3), ALLOCATABLE :: p3(:,:)
TYPE(poly4), ALLOCATABLE :: p4(:,:)

nvar=crystal_parameters(ibrav_save)
ngroup = work_base/ngeo_strain
ALLOCATE(x_t(nvar,ntemp))

filename='anhar_files/'//TRIM(flanhar)//'.celldm'
CALL read_alpha_anis(ibrav_save, nvar, x_t, temp, ntemp, filename)

CALL divide(world_comm, ntemp, startt, lastt)

ALLOCATE(p1(startt:lastt,ngroup))
ALLOCATE(p2(startt:lastt,ngroup))
ALLOCATE(p3(startt:lastt,ngroup))
ALLOCATE(p4(startt:lastt,ngroup))

ALLOCATE(y(nvar+1))
ALLOCATE(free_ener(work_base))
ALLOCATE(epsilon_geo_loc(3,3,work_base))


CALL interpolate_free_ener_strain(p1,p2,p3,p4,startt,lastt,ngroup)

DO itemp = startt, lastt
   CALL expand_celldm(celldm_t, x_t(1,itemp), nvar, ibrav_save)
   vmin_t(itemp)=compute_omega_geo(ibrav_save, celldm_t)
   DO igroup=1, ngroup
      DO istrain=1, ngeo_strain
         ind = istrain + (igroup-1)*ngeo_strain
            
         epsilon_geo_loc(:,:,ind) = epsilon_geo(:,:,ind)
         
         y(1:nvar) = x_t(1:nvar, itemp)
         y(nvar+1) = epsil_geo(ind) 
         
         IF (poly_degree_elc==4) THEN
            CALL evaluate_fit_quartic(nvar+1, y, & 
                                free_ener(ind), p4(itemp,igroup))
        
         ELSEIF (poly_degree_elc==3) THEN
            CALL evaluate_fit_cubic(nvar+1, y, & 
                                free_ener(ind), p3(itemp,igroup))
         ELSEIF (poly_degree_elc==2) THEN
            CALL evaluate_fit_quadratic(nvar+1, y, & 
                                free_ener(ind), p2(itemp,igroup))

         ELSEIF (poly_degree_elc==1) THEN 
            CALL evaluate_fit_linear(nvar+1, y, & 
                                free_ener(ind), p1(itemp,igroup))
         
         ELSE 
            CALL errore('manage_elastic_constants_qha_2',&
                                             'wrong poly_degree_elc',1)
         ENDIF
      ENDDO
      IF (poly_degree_elc==4) THEN
            CALL clean_poly(p4(itemp,igroup))
         ELSEIF (poly_degree_elc==3) THEN
            CALL clean_poly(p3(itemp,igroup))
         ELSEIF (poly_degree_elc==2) THEN
            CALL clean_poly(p2(itemp,igroup))
         ELSEIF (poly_degree_elc==1) THEN
            CALL clean_poly(p1(itemp,igroup))
      ENDIF
   ENDDO

   !ibrav_save

   CALL compute_elastic_constants_ene(free_ener, epsilon_geo_loc,  &
                               work_base, ngeo_strain, ibrav_save, laue,  &
                               vmin_t(itemp), elcpvar)

   el_cons_t(:,:,itemp) = el_con(:,:)
ENDDO  

CALL mp_sum(el_cons_t, world_comm)
CALL compute_el_comp_t(el_cons_t,el_comp_t,b0_t)

filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue, el_cons_t, b0_t, &
                                                          filelastic, 0)
filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp'
CALL write_el_cons_on_file(temp, ntemp, ibrav_save, laue, el_comp_t, b0_t, &
                                                          filelastic, 1)
DEALLOCATE(p1)
DEALLOCATE(p2)
DEALLOCATE(p3)
DEALLOCATE(p4)

DEALLOCATE(free_ener)
DEALLOCATE(epsilon_geo_loc)

lelastic=.TRUE.
ngeom=1
CALL plot_elastic_t(0,.FALSE.)
CALL plot_elastic_t(1,.FALSE.)

RETURN
END SUBROUTINE manage_elastic_cons_qha_2

SUBROUTINE read_alpha_anis(ibrav, nvar, x_t, temp, ntemp, filename)
USE kinds, ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibrav, ntemp, nvar
REAL(DP), INTENT(IN) :: temp(ntemp)
REAL(DP), INTENT(OUT) :: x_t(nvar,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER :: itemp, iu_therm, ios
INTEGER :: find_free_unit

REAL :: rdum_temp, rdum

x_t=0.0_DP

IF (meta_ionode) THEN
   iu_therm=find_free_unit()
   OPEN(UNIT=iu_therm, FILE=TRIM(filename), FORM='FORMATTED', & 
                                   STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
ENDIF
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('read_alpha_anis','opening celldm (T) file', ABS(ios)) 
                                     
IF (meta_ionode) THEN
   READ(iu_therm,*)
   IF (ibrav==1 .OR. ibrav==2 .OR. ibrav==3 ) THEN
      DO itemp = 1, ntemp-1
         READ(iu_therm, '(e12.5,2e20.9)') rdum_temp, x_t(1,itemp), rdum
         IF (ABS(rdum_temp-temp(itemp))>1D-5) &
            CALL errore('read_alpha_anis','uncorrect temperature', 1)
   END DO

   ELSEIF (ibrav==4 .OR. ibrav==6 .OR. ibrav==7 ) THEN
      DO itemp = 1, ntemp-1
         READ(iu_therm, '(e12.5,4e20.9)') rdum_temp, x_t(1,itemp), &
                                                     x_t(2,itemp), &
                                                     rdum, rdum
         IF (ABS(rdum_temp-temp(itemp))>1D-5) &
            CALL errore('read_alpha_anis','uncorrect temperature', 1)
      END DO

   ELSEIF ( ibrav==5 ) THEN
      DO itemp = 1, ntemp-1
         READ(iu_therm, '(e12.5,4e20.9)') rdum_temp, x_t(1,itemp), &
                                                     x_t(2,itemp), &
                                                     rdum, rdum
         IF (ABS(rdum_temp-temp(itemp))>1D-5) &
            CALL errore('read_alpha_anis','uncorrect temperature', 1)
      END DO

   ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
      DO itemp = 1, ntemp-1
         READ(iu_therm, '(e12.5,6e20.9)') rdum_temp, x_t(1,itemp), &
                                                     x_t(2,itemp), &
                                                     x_t(3,itemp), &
                                                     rdum, rdum, rdum
         IF (ABS(rdum_temp-temp(itemp))>1D-5) &
            CALL errore('read_alpha_anis','uncorrect temperature', 1)
      END DO

   ELSEIF (ibrav==12 .OR. ibrav==13) THEN
      DO itemp = 1, ntemp-1
         READ(iu_therm, '(e12.5,4e17.9)') rdum_temp, x_t(1,itemp), &
                                                     x_t(2,itemp), &
                                                     x_t(3,itemp), &
                                                     x_t(4,itemp)
         IF (ABS(rdum_temp-temp(itemp))>1D-5) &
            CALL errore('read_alpha_anis','uncorrect temperature', 1)
      END DO
   ELSEIF (ibrav==-12 .OR. ibrav==-13) THEN
      DO itemp = 1, ntemp-1
         READ(iu_therm, '(e12.5,4e17.9)') rdum_temp, x_t(1,itemp), &
                                                     x_t(2,itemp), &
                                                     x_t(3,itemp), &
                                                     x_t(4,itemp)
         IF (ABS(rdum_temp-temp(itemp))>1D-5) &
             CALL errore('read_alpha_anis','uncorrect temperature', 1)
      END DO
   ELSEIF (ibrav==14) THEN
      DO itemp = 1, ntemp-1
         READ(iu_therm, '(e12.5,6e15.7)') rdum_temp, x_t(1,itemp), &
                                                     x_t(2,itemp), &
                                                     x_t(3,itemp), &
                                                     x_t(4,itemp), &
                                                     x_t(5,itemp), &
                                                     x_t(6,itemp)
         IF (ABS(rdum_temp-temp(itemp))>1D-5) &
              CALL errore('read_alpha_anis','uncorrect temperature', 1)
      END DO
   ELSE IF (ibrav==0) THEN
!
!  In this case we write nothing but do not stop
!
   ELSE
      CALL errore('read_alpha_anis','ibrav not programmed',1)
   END IF
   x_t(:,ntemp)=x_t(:,ntemp-1)

   CLOSE(iu_therm)
ENDIF

CALL mp_bcast(x_t, meta_ionode_id, world_comm)

RETURN
END SUBROUTINE read_alpha_anis
