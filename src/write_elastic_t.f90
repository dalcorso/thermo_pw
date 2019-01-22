!
! Copyright (C) 2016 Andrea Dal Corso
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
USE thermo_mod, ONLY : ngeo, ibrav_geo, celldm_geo
USE control_quartic_energy, ONLY : lquartic_elc, lsolve
USE quadratic_surfaces, ONLY : print_chisq_quadratic, fit_multi_quadratic,    &
                             introduce_quadratic_fit, summarize_fitting_data, &
                             print_quadratic_polynomial, quadratic_var
USE quartic_surfaces, ONLY : quartic_var, print_chisq_quartic,         &
                             fit_multi_quartic, introduce_quartic_fit, &
                             print_quartic_polynomial
USE control_elastic_constants, ONLY : el_con_geo,  el_cons_available, &
                                      el_cons_t_available
USE control_grun,   ONLY : lb0_t
USE control_mur,    ONLY : b0
USE elastic_constants, ONLY : el_con, el_compliances
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_macro_elasticity, ONLY: macro_el
USE anharmonic,     ONLY : celldm_t, el_cons_t, el_comp_t, b0_t, lelastic
USE ph_freq_anharmonic, ONLY : celldmf_t, el_consf_t, el_compf_t, b0f_t, &
                           lelasticf
USE grun_anharmonic, ONLY : el_cons_grun_t, el_comp_grun_t, lelastic_grun
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, degree, nvar, ndata
REAL(DP), ALLOCATABLE :: el_cons_coeff(:,:,:), x(:,:), f(:)
INTEGER :: i, j, idata, itemp
INTEGER :: compute_nwork

IF (.NOT.lb0_t.AND.el_cons_available) THEN
   b0=macro_el(5)
   b0_t=macro_el(5)
   b0f_t=macro_el(5)
   DO itemp=1,ntemp
      el_cons_t(:,:,itemp)=el_con(:,:)
      el_consf_t(:,:,itemp)=el_con(:,:)
      el_cons_grun_t(:,:,itemp)=el_con(:,:)
      el_comp_grun_t(:,:,itemp)=el_compliances(:,:)
   END DO
   lelastic=.TRUE.
   lelasticf=.TRUE.
   lelastic_grun=.TRUE.
   RETURN
END IF

IF (.NOT.el_cons_t_available) RETURN

ibrav=ibrav_geo(1)
degree=crystal_parameters(ibrav)

ndata=compute_nwork()
ALLOCATE(x(degree,ndata))
ALLOCATE(f(ndata))

IF (lquartic_elc) THEN
   nvar=quartic_var(degree)
ELSE
   nvar=quadratic_var(degree)
ENDIF
ALLOCATE(el_cons_coeff(nvar,6,6))
!
!  Part 1 evaluation of the polynomial coefficients
!
CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_geo)

DO i=1,6
   DO j=i,6
      IF (el_con_geo(i,j,1)>0.1_DP) THEN
         WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")')&
                       i, j
         DO idata=1,ndata
            f(idata)=el_con_geo(i,j,idata)
         END DO

         IF (lquartic_elc) THEN
!            CALL introduce_quartic_fit(degree, nvar, ndata)
            CALL fit_multi_quartic(ndata,degree,nvar,lsolve,x,f,&
                                                   el_cons_coeff(:,i,j))
!            CALL print_quartic_polynomial(degree, nvar, el_cons_coeff(:,i,j))
       !    WRITE(stdout,'(/,7x,"Energy (1)    Fitted energy (2)   &
       !                                                   &DeltaE (1)-(2)")') 
!            CALL print_chisq_quartic(ndata, degree, nvar, x, f, &
!                                                      el_cons_coeff(:,i,j))
         ELSE
!            CALL introduce_quadratic_fit(degree, nvar, ndata)
!            CALL summarize_fitting_data(degree, ndata, x, f)
            CALL fit_multi_quadratic(ndata,degree,nvar,x,f,el_cons_coeff(:,i,j))
!            CALL print_quadratic_polynomial(degree, nvar, el_cons_coeff(:,i,j))
!            CALL print_chisq_quadratic(ndata, degree, nvar, x, f, &
!                                                        el_cons_coeff(:,i,j))
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
   CALL interpolate_el_cons(celldm_t, nvar, degree, ibrav, el_cons_coeff,&
                            lquartic_elc, el_cons_t, el_comp_t, b0_t)
   el_cons_grun_t=el_cons_t
   el_comp_grun_t=el_comp_t
   lelastic=.TRUE.
   lelastic_grun=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, el_cons_t, b0_t, &
                                                       filelastic, 0)
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, el_comp_t, b0_t, & 
                                                       filelastic, 1)
   
ENDIF

IF (ltherm_freq) THEN
   CALL interpolate_el_cons(celldmf_t, nvar, degree, ibrav, el_cons_coeff,&
                                  lquartic_elc, el_consf_t, el_compf_t, b0f_t)
   el_cons_grun_t=el_consf_t
   el_comp_grun_t=el_compf_t
   lelasticf=.TRUE.
   lelastic_grun=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_ph'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, el_consf_t, b0f_t, &
                                                          filelastic, 0)

   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp_ph'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, el_compf_t, b0f_t, &
                                                           filelastic,1)
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DEALLOCATE(el_cons_coeff)

RETURN
END SUBROUTINE write_elastic_t

SUBROUTINE write_el_cons_on_file(temp, ntemp, ibrav, el_cons_t, b0, &
                                                             filename, iflag)
!
!  iflag=0 writes the elastic constants
!  iflag=1 writes the elastic compliances
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE thermo_sym, ONLY : laue
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, ibrav, iflag
REAL(DP), INTENT(IN) :: temp(ntemp), el_cons_t(6,6,ntemp), b0(ntemp)
REAL(DP) :: b0_t(ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_el_cons, ios
INTEGER :: find_free_unit

iu_el_cons=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_el_cons, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_el_cons_on_file','opening elastic constants (T) file',&
                                                             ABS(ios))
!
! Plot the compressibility togheter with the elastic complinances
!

IF (iflag/=0) THEN
   b0_t=1.0_DP/b0
ELSE
   b0_t=b0
ENDIF

IF (meta_ionode) THEN
   SELECT CASE (laue)
      CASE(29,32)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
                  & 13x, "     C_12 ", 13x, "     C_44 ")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
                  & 13x, "     S_12 ", 13x, "     S_44 ")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,4e20.12)') temp(itemp), b0_t(itemp), &
                 el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                 el_cons_t(4,4,itemp)
         ENDDO
      CASE(25)
!
!     D_3d
!            
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
                  & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, &
                       &" C_44 ", 13x, " C_14")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
                  & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, &
                       &" S_44 ", 13x, " S_14")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,7e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(1,4,itemp)
         ENDDO
      CASE(27)
!
!     S_6
!            
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", " C_11 ", 13x, &
                  &" C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44", 13x, &
                  &" C_14", 13x, "C_25" )')
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", " S_11 ", 13x, &
                  &" S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44", 13x, &
                  &" S_14", 13x, "S_25" )')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,8e20.12)')  temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(1,4,itemp), &
                  el_cons_t(2,5,itemp)
         ENDDO
      CASE(19,23)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", 13x,&
                  &" C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44" )')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", 13x,&
                  &" S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44" )')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,6e20.12)') temp(itemp),  b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp)
         ENDDO
      CASE(22)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", 13x,&
                  &" C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44", 13x, & 
                  &" C_66 " )')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", 13x,&
                  &" S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44", 13x, & 
                  &" S_66 " )')
         ENDIF

         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,7e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp)
         ENDDO
      CASE(20)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B", 13x, " C_11 ", 13x,&
                  &" C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", 13x,&
                  &" C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K", 13x, " S_11 ", 13x,&
                  &" S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", 13x,&
                  &" S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,10e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                  el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                  el_cons_t(6,6,itemp)
         ENDDO
      CASE(18)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
                  & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_33 ", 13x, "C_44", &
                  & 13x, " C_66 ", 13x, " C_16 ")')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
                  & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_33 ", 13x, "S_44", &
                  & 13x, " S_66 ", 13x, " S_16 ")')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,8e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp), &
                  el_cons_t(1,6,itemp)
         ENDDO
      CASE(16)
         IF (ibrav > 0) THEN
            !
            !  b unique
            !
            IF (iflag==0) THEN
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
               & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", &
               & 13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ", &
               & 13x, " C_15 ", 13x, " C_25 ", 13x, " C_35 ", 13x, " C_46 ")')
            ELSE
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
               & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", &
               & 13x, " S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ", &
               & 13x, " S_15 ", 13x, " S_25 ", 13x, " S_35 ", 13x, " S_46 ")')
            ENDIF
            DO itemp=2,ntemp-1
               WRITE(iu_el_cons,'(e16.8,14e20.12)') temp(itemp), b0_t(itemp),&
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,5,itemp), &
                     el_cons_t(2,5,itemp), el_cons_t(3,5,itemp), &
                     el_cons_t(4,6,itemp)
            ENDDO
         ELSE
            !
            !  c unique
            !
            IF (iflag==0) THEN
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ",13x," C_11 ",  &
               & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", &
               & 13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ", &
               & 13x, " C_16 ", 13x, " C_26 ", 13x, " C_36 ", 13x, " C_45 ")')
            ELSE
               WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ",13x," S_11 ",  &
               & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", &
               & 13x, " S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ", &
               & 13x, " S_16 ", 13x, " S_26 ", 13x, " S_36 ", 13x, " S_45 ")')
            ENDIF
            DO itemp=2,ntemp-1
               WRITE(iu_el_cons,'(e16.8,14e20.12)') temp(itemp), b0_t(itemp),&
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,6,itemp), &
                     el_cons_t(2,6,itemp), el_cons_t(3,6,itemp), &
                     el_cons_t(4,5,itemp)
            ENDDO
         ENDIF
      CASE(2)
         IF (iflag==0) THEN
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " B ", 13x, " C_11 ", &
               & 13x, " C_12 ", 13x, " C_13 ", 13x, " C_22 ", 13x, " C_23 ", &
               & 13x, " C_33 ", 13x, " C_44 ", 13x, " C_55 ", 13x, " C_66 ", &
               & 13x, " C_14 ", 13x, " C_15 ", 13x, " C_16 ", 13x, " C_24 ", &
               & 13x, " C_25 ", 13x, " C_26 ", 13x, " C_34 ", 13x, " C_35 ", &
               & 13x, " C_36 ", 13x, " C_45 ", 13x, " C_46 ", 13x, " C_55 ", &
               & 13x, " C_56 ", 13x, " C_66 " )')
         ELSE
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " K ", 13x, " S_11 ", &
               & 13x, " S_12 ", 13x, " S_13 ", 13x, " S_22 ", 13x, " S_23 ", &
               & 13x, " S_33 ", 13x, " S_44 ", 13x, " S_55 ", 13x, " S_66 ", &
               & 13x, " S_14 ", 13x, " S_15 ", 13x, " S_16 ", 13x, " S_24 ", &
               & 13x, " S_25 ", 13x, " S_26 ", 13x, " S_34 ", 13x, " S_35 ", &
               & 13x, " S_36 ", 13x, " S_45 ", 13x, " S_46 ", 13x, " S_55 ", &
               & 13x, " S_56 ", 13x, " S_66 " )')
         ENDIF
         DO itemp=2,ntemp-1
            WRITE(iu_el_cons,'(e16.8,24e20.12)') temp(itemp), b0_t(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                  el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                  el_cons_t(6,6,itemp), el_cons_t(1,4,itemp), &
                  el_cons_t(1,5,itemp), el_cons_t(1,6,itemp), &
                  el_cons_t(2,4,itemp), el_cons_t(2,5,itemp), &
                  el_cons_t(2,6,itemp), el_cons_t(3,4,itemp), &
                  el_cons_t(3,5,itemp), el_cons_t(3,6,itemp), &
                  el_cons_t(4,5,itemp), el_cons_t(4,6,itemp), &
                  el_cons_t(5,5,itemp), el_cons_t(5,6,itemp), &
                  el_cons_t(6,6,itemp)
         ENDDO
   END SELECT
   CLOSE(iu_el_cons)
ENDIF

RETURN
END SUBROUTINE write_el_cons_on_file
