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
!  Helmholtz free energy at temperature T.
!
USE kinds, ONLY : DP
USE io_global, ONLY : stdout
USE mp_images, ONLY : my_image_id, root_image
USE thermo_mod, ONLY : tot_ngeo
USE control_quartic_energy, ONLY :  nvar4, lquartic, lsolve
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic, fit_multi_quadratic, &
                             introduce_quadratic_fit, summarize_fitting_data, &
                             print_quadratic_polynomial
USE quartic_surfaces, ONLY : compute_quartic_var, evaluate_fit_quartic, &
                             fit_multi_quartic, introduce_quartic_fit, &
                             print_quartic_polynomial
USE control_elastic_constants, ONLY : el_con_geo, el_con_celldm_geo, &
                                      el_con_ibrav_geo, el_cons_available, &
                                      el_cons_t_available
USE elastic_constants, ONLY : print_macro_elasticity, &
                              compute_elastic_compliances, el_con
USE control_grun,   ONLY : lb0_t
USE equilibrium_conf,    ONLY : celldm0
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_macro_elasticity, ONLY: macro_el
USE anharmonic,     ONLY : celldm_t, el_cons_t, el_comp_t, macro_el_t, b0_t, &
                           lelastic
USE ph_freq_anharmonic, ONLY : celldmf_t, el_consf_t, el_compf_t, macro_elf_t, &
                           b0f_t, lelasticf
USE thermo_sym, ONLY : laue

USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, degree, nvar, ndata
REAL(DP), ALLOCATABLE :: el_cons_coeff(:,:,:), x(:,:), f(:), xfit(:), &
                         el_cons_coeff4(:,:,:)
REAL(DP) :: chisq, perc, aux
INTEGER :: iu_el_cons, i, j, idata, itemp, ios

IF ( my_image_id /= root_image ) RETURN

IF (.NOT.lb0_t.AND.el_cons_available) THEN
   IF (.NOT. (ltherm_dos.OR.ltherm_freq)) THEN
      DO itemp=1,ntemp
         celldm_t(:,itemp)=celldm0(:)
      ENDDO
   ENDIF
   b0_t=macro_el(5)
   b0f_t=macro_el(5)
   DO itemp=1,ntemp
      el_cons_t(:,:,itemp)=el_con(:,:)
      el_consf_t(:,:,itemp)=el_con(:,:)
   END DO
   lelastic=.TRUE.
   IF (ltherm_freq) lelasticf=.TRUE.
   RETURN
END IF

IF (.NOT.el_cons_t_available) RETURN

ibrav=el_con_ibrav_geo(1)
CALL compute_degree(ibrav,degree,nvar)

ndata=tot_ngeo
ALLOCATE(x(degree,ndata))
ALLOCATE(f(ndata))
ALLOCATE(el_cons_coeff(nvar,6,6))
ALLOCATE(xfit(degree))
IF (lquartic) THEN
   nvar4=compute_quartic_var(degree)
   ALLOCATE(el_cons_coeff4(nvar4,6,6))
ENDIF

CALL introduce_quadratic_fit(degree, nvar, ndata)
!
!  Part 1 evaluation of the polynomial coefficients
!
CALL set_x_from_celldm(ibrav, degree, ndata, x, el_con_celldm_geo)

DO i=1,6
   DO j=i,6
      IF (el_con_geo(i,j,1)>0.1_DP) THEN
         WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")')&
                       i, j
         DO idata=1,ndata
            f(idata)=el_con_geo(i,j,idata)
         END DO

         CALL summarize_fitting_data(degree, ndata, x, f)
         CALL fit_multi_quadratic(ndata,degree,nvar,x,f,el_cons_coeff(:,i,j))
         CALL print_quadratic_polynomial(degree, nvar, el_cons_coeff(:,i,j))

         chisq=0.0_DP
         perc=0.0_DP
         DO idata=1,ndata
            CALL evaluate_fit_quadratic(degree,nvar,x(1,idata),aux,&
                                                   el_cons_coeff(:,i,j))
            WRITE(stdout,'(3f19.12)') f(idata), aux, f(idata)-aux
            chisq = chisq + (aux - f(idata))**2 
            perc= perc + ABS((f(idata)-aux) / f(idata)) 
         ENDDO
         WRITE(stdout,'(5x,"chi square=",e18.5," relative error",e18.5,&
                                     &" %",/)') chisq / ndata, perc / ndata

         IF (lquartic) THEN

            CALL introduce_quartic_fit(degree, nvar4, ndata)

            CALL fit_multi_quartic(ndata,degree,nvar4,lsolve,x,f,&
                                                   el_cons_coeff4(:,i,j))
            CALL print_quartic_polynomial(degree, nvar4, el_cons_coeff4(:,i,j))
       !    WRITE(stdout,'(/,7x,"Energy (1)    Fitted energy (2)   &
       !                                                   &DeltaE (1)-(2)")') 
            chisq=0.0_DP
            perc=0.0_DP
            DO idata=1,ndata
               CALL evaluate_fit_quartic(degree,nvar4,x(1,idata),aux,&
                                                         el_cons_coeff4(:,i,j))
!
!              WRITE(stdout,'(2f10.4,2f19.12,e19.12)') x(1,idata), x(2,idata),
!                                   f(idata), aux, f(idata)-aux
               chisq = chisq + (aux - f(idata))**2
               perc= perc + ABS((f(idata)-aux) / f(idata)) 
            ENDDO
            WRITE(stdout,'(5x,"chi square=",e18.5," relative error",e18.5,&
                           & " %",/)') chisq / ndata, perc / ndata
         ENDIF
      ENDIF
   ENDDO
ENDDO

IF (ltherm_dos) THEN
   DO itemp=1,ntemp
      CALL compute_x(celldm_t(:,itemp),xfit,degree,ibrav)
      DO i=1,6
         DO j=i,6
            IF (el_con_geo(i,j,1)>0.1_DP) THEN
               IF (lquartic) THEN
                  CALL evaluate_fit_quartic(degree,nvar4,xfit,aux,&
                                                     el_cons_coeff4(:,i,j))
               ELSE
                  CALL evaluate_fit_quadratic(degree,nvar,xfit,aux,&
                                                   el_cons_coeff(:,i,j))
               ENDIF
               el_cons_t(i,j,itemp)=aux
               el_cons_t(j,i,itemp)=aux
            ENDIF
         END DO
      END DO
      CALL compute_elastic_compliances(el_cons_t(:,:,itemp),&
                                            el_comp_t(:,:,itemp))
      CALL print_macro_elasticity(ibrav,el_cons_t(:,:,itemp), &
                          el_comp_t(:,:,itemp),macro_el_t(:,itemp),.FALSE.)
      b0_t(itemp)=macro_el_t(5,itemp)
   END DO
   lelastic=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
   CALL write_el_cons_on_file(temp, ntemp, laue, ibrav, el_cons_t, b0_t, &
                                                              filelastic)
ENDIF

IF (ltherm_freq) THEN
   DO itemp=1,ntemp
      CALL compute_x(celldmf_t(:,itemp),xfit,degree,ibrav)
      DO i=1,6
         DO j=i,6
            IF (el_con_geo(i,j,1)>0.1_DP) THEN
               IF (lquartic) THEN
                  CALL evaluate_fit_quartic(degree,nvar4,xfit,aux,&
                                                     el_cons_coeff4(:,i,j))
               ELSE
                  CALL evaluate_fit_quadratic(degree,nvar,xfit,aux,&
                                                   el_cons_coeff(:,i,j))
               ENDIF
               el_consf_t(i,j,itemp)=aux
               el_consf_t(j,i,itemp)=aux
            ENDIF
         END DO
      END DO
      CALL compute_elastic_compliances(el_consf_t(:,:,itemp),&
                                            el_compf_t(:,:,itemp))
      CALL print_macro_elasticity(ibrav,el_consf_t(:,:,itemp), &
                          el_compf_t(:,:,itemp),macro_elf_t(:,itemp),.FALSE.)
      b0f_t(itemp)=macro_elf_t(5,itemp)
   END DO
   lelasticf=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons_ph'
   CALL write_el_cons_on_file(temp, ntemp, laue, ibrav, el_consf_t, b0f_t, &
                                                               filelastic)
ENDIF

DEALLOCATE(x)
DEALLOCATE(f)
DEALLOCATE(el_cons_coeff)
DEALLOCATE(xfit)
IF (lquartic) DEALLOCATE(el_cons_coeff4)

RETURN
END SUBROUTINE write_elastic_t

SUBROUTINE write_el_cons_on_file(temp, ntemp, laue, ibrav, el_cons_t, &
                                                           b0_t, filename)

USE kinds,     ONLY : DP
USE io_global, ONLY : ionode, ionode_id, stdout
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, laue, ibrav
REAL(DP), INTENT(IN) :: temp(ntemp), el_cons_t(6,6,ntemp), b0_t(ntemp)
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_el_cons, ios

iu_el_cons=2
IF (ionode) &
   OPEN(UNIT=iu_el_cons, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, ionode_id, intra_image_comm)
   CALL errore('write_el_cons_on_file','opening elastic constants (T) file',&
                                                             ABS(ios))

IF (ionode) THEN
   SELECT CASE (laue)
      CASE(29,32)
         WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, "  C_12 ",&
                  & 9x, " C_44 ", 9x, " B " )')
         DO itemp=1,ntemp
            WRITE(iu_el_cons,'(e16.8,4e20.12)') temp(itemp), &
                 el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                 el_cons_t(4,4,itemp), b0_t(itemp)
         ENDDO
      CASE(19,23)
         WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, " C_12 ",&
                  & 9x, " C_13 ", 9x, " C_33 ", 9x, "C_44", 9x, " B " )')
         DO itemp=1,ntemp
            WRITE(iu_el_cons,'(e16.8,6e20.12)') temp(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), b0_t(itemp)
         ENDDO
      CASE(22)
         WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, " C_12 ",&
                  & 9x, " C_13 ", 9x, " C_33 ", 9x, "C_44", 9x, & 
                  &     " C_66 ", 9x, " B " )')
         DO itemp=1,ntemp
            WRITE(iu_el_cons,'(e16.8,7e20.12)') temp(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp), b0_t(itemp)
         ENDDO
      CASE(20)
         WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, " C_12 ", &
                  & 9x, " C_13 ", 9x, " C_22 ", 9x, " C_23 ", 9x, " C_33 ",  &
                  & 9x, " C_44 ", 9x, " C_55 ", 9x, " C_66 ", 9x, " B " )')
         DO itemp=1,ntemp
            WRITE(iu_el_cons,'(e16.8,10e20.12)') temp(itemp), &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                  el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                  el_cons_t(6,6,itemp), b0_t(itemp)
         ENDDO
      CASE(18)
         WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, " C_12 ",&
                  & 9x, " C_13 ", 9x, " C_33 ", 9x, "C_44", 9x, &
                  &     " C_66 ", 9x, " C_16 ", 9x, " B " )')

         DO itemp=1,ntemp
            WRITE(iu_el_cons,'(e16.8,8e20.12)') temp(itemp),  &
                  el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                  el_cons_t(1,3,itemp), el_cons_t(3,3,itemp), &
                  el_cons_t(4,4,itemp), el_cons_t(6,6,itemp), &
                  el_cons_t(1,6,itemp), b0_t(itemp)
         ENDDO
      CASE(16)
         IF (ibrav > 0) THEN
            !
            !  b unique
            !
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, " C_12 ", &
                  & 9x, " C_13 ", 9x, " C_22 ", 9x, " C_23 ", 9x, " C_33 ",  &
                  & 9x, " C_44 ", 9x, " C_55 ", 9x, " C_66 ", 9x, " C_15 ",  &
                  & 9x, " C_25 ", 9x, " C_35 ", 9x, " C_46 ", 9x, " B " )')
            DO itemp=1,ntemp
               WRITE(iu_el_cons,'(e16.8,14e20.12)') temp(itemp), &
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,5,itemp), &
                     el_cons_t(2,5,itemp), el_cons_t(3,5,itemp), &
                     el_cons_t(4,6,itemp), b0_t(itemp)
            ENDDO
         ELSE
            !
            !  c unique
            !
            WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, " C_12 ", &
                  & 9x, " C_13 ", 9x, " C_22 ", 9x, " C_23 ", 9x, " C_33 ",  &
                  & 9x, " C_44 ", 9x, " C_55 ", 9x, " C_66 ", 9x, " C_16 ",  &
                  & 9x, " C_26 ", 9x, " C_36 ", 9x, " C_45 ", 9x, " B " )')
            DO itemp=1,ntemp
               WRITE(iu_el_cons,'(e16.8,14e20.12)') temp(itemp), &
                     el_cons_t(1,1,itemp), el_cons_t(1,2,itemp), &
                     el_cons_t(1,3,itemp), el_cons_t(2,2,itemp), &
                     el_cons_t(2,3,itemp), el_cons_t(3,3,itemp), &
                     el_cons_t(4,4,itemp), el_cons_t(5,5,itemp), &
                     el_cons_t(6,6,itemp), el_cons_t(1,6,itemp), &
                     el_cons_t(2,6,itemp), el_cons_t(3,6,itemp), &
                     el_cons_t(4,5,itemp), b0_t(itemp)
            ENDDO
         ENDIF
      CASE(2)
         WRITE(iu_el_cons,'("#",5x,"   T  ", 10x, " C_11 ", 9x, " C_12 ",   &
                  & 9x, " C_13 ", 9x, " C_22 ", 9x, " C_23 ", 9x, " C_33 ", &
                  & 9x, " C_44 ", 9x, " C_55 ", 9x, " C_66 ", 9x, &
                  & 9x, " C_14 ", 9x, " C_15 ", 9x, " C_16 ", 9x, &
                  & 9x, " C_24 ", 9x, " C_25 ", 9x, " C_26 ", 9x, &
                  & 9x, " C_34 ", 9x, " C_35 ", 9x, " C_36 ", 9x, &
                  & 9x, " C_45 ", 9x, " C_46 ", 9x, " C_55 ", 9x, &
                  & 9x, " C_56 ", 9x, " C_66 ", 9x, " B " )')
         DO itemp=1,ntemp
            WRITE(iu_el_cons,'(e16.8,24e20.12)') temp(itemp), &
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
                  el_cons_t(6,6,itemp), b0_t(itemp)
         ENDDO

   END SELECT
   CLOSE(iu_el_cons)
ENDIF

RETURN
END SUBROUTINE write_el_cons_on_file
