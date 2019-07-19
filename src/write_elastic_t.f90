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
USE thermo_sym, ONLY : laue
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE quadratic_surfaces, ONLY : quadratic_ncoeff, fit_multi_quadratic
USE quartic_surfaces, ONLY : quartic_ncoeff, fit_multi_quartic
USE cubic_surfaces,   ONLY : cubic_ncoeff, fit_multi_cubic
USE elastic_constants, ONLY : write_el_cons_on_file
USE control_elastic_constants, ONLY : el_con_geo
USE lattices,       ONLY : crystal_parameters
USE control_thermo, ONLY : ltherm_dos, ltherm_freq
USE control_macro_elasticity, ONLY: macro_el
USE anharmonic,     ONLY : celldm_t, el_cons_t, el_comp_t, b0_t, lelastic
USE ph_freq_anharmonic, ONLY : celldmf_t, el_consf_t, el_compf_t, b0f_t, &
                           lelasticf
USE data_files, ONLY : flanhar
USE temperature, ONLY : ntemp, temp
IMPLICIT NONE
CHARACTER(LEN=256) :: filelastic
INTEGER :: igeo, ibrav, nvar, ncoeff, ndata
REAL(DP), ALLOCATABLE :: el_cons_coeff(:,:,:), x(:,:), f(:)
INTEGER :: i, j, idata, itemp
INTEGER :: compute_nwork

ibrav=ibrav_geo(1)
nvar=crystal_parameters(ibrav)

ndata=compute_nwork()
ALLOCATE(x(nvar,ndata))
ALLOCATE(f(ndata))

IF (poly_degree_elc==4) THEN
   ncoeff=quartic_ncoeff(nvar)
ELSEIF(poly_degree_elc==3) THEN
   ncoeff=cubic_ncoeff(nvar)
ELSE
   ncoeff=quadratic_ncoeff(nvar)
ENDIF
ALLOCATE(el_cons_coeff(ncoeff,6,6))
!
!  Part 1 evaluation of the polynomial coefficients
!
CALL set_x_from_celldm(ibrav, nvar, ndata, x, celldm_geo)

DO i=1,6
   DO j=i,6
      IF (el_con_geo(i,j,1)>0.1_DP) THEN
         WRITE(stdout,'(/,5x,"Fitting elastic constants C(",i4,",",i4,")")')&
                       i, j
         DO idata=1,ndata
            f(idata)=el_con_geo(i,j,idata)
         END DO

         IF (poly_degree_elc==4) THEN
            CALL fit_multi_quartic(ndata,nvar,ncoeff,lsolve,x,f,&
                                                   el_cons_coeff(:,i,j))
         ELSEIF (poly_degree_elc==3) THEN
            CALL fit_multi_cubic(ndata,nvar,ncoeff,lsolve,x,f,&
                                                   el_cons_coeff(:,i,j))
         ELSE
            CALL fit_multi_quadratic(ndata,nvar,ncoeff,x,f,el_cons_coeff(:,i,j))
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
   CALL interpolate_el_cons(celldm_t, ncoeff, nvar, ibrav, el_cons_coeff,&
                            poly_degree_elc, el_cons_t, el_comp_t, b0_t)
   lelastic=.TRUE.
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_cons'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_cons_t, b0_t, &
                                                       filelastic, 0)
   filelastic='anhar_files/'//TRIM(flanhar)//'.el_comp'
   CALL write_el_cons_on_file(temp, ntemp, ibrav, laue, el_comp_t, b0_t, & 
                                                       filelastic, 1)
   
ENDIF

IF (ltherm_freq) THEN
   CALL interpolate_el_cons(celldmf_t, ncoeff, nvar, ibrav, el_cons_coeff,&
                               poly_degree_elc, el_consf_t, el_compf_t, b0f_t)
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
DEALLOCATE(el_cons_coeff)

RETURN
END SUBROUTINE write_elastic_t

