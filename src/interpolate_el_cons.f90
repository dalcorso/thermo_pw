!
! Copyright (C) 2018 - 2020 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE interpolate_el_cons(celldm_t, nvar, ibrav, ec_p1, ec_p2, &
                    ec_p3, ec_p4, poly_degree_elc, el_cons_t, el_comp_t, b0_t)
!----------------------------------------------------------------------------
!
! This routine receives as input the coeffients of polynomials which 
! interpolate the elastic constants on a grid of geometries,  
! the celldm that correspond to each temperature and gives as output
! the elastic constants, the elastic compliances and the bulk modulus
! interpolated at the celldm that corresponds to each temperature.
!  
USE kinds,              ONLY : DP
USE temperature,        ONLY : ntemp
USE elastic_constants,  ONLY : print_macro_elasticity, &
                              compute_elastic_compliances, el_con
USE control_elastic_constants, ONLY : el_con_geo
USE linear_surfaces,    ONLY : evaluate_fit_linear
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic
USE cubic_surfaces,     ONLY : evaluate_fit_cubic
USE quartic_surfaces,   ONLY : evaluate_fit_quartic
USE lattices,           ONLY : compress_celldm
USE polynomial,         ONLY : poly1, poly2, poly3, poly4
USE mp_world,           ONLY : world_comm
USE mp,                 ONLY : mp_sum

IMPLICIT NONE
INTEGER :: ibrav, nvar
INTEGER :: poly_degree_elc
REAL(DP) :: celldm_t(6, ntemp), el_cons_t(6,6,ntemp), el_comp_t(6,6,ntemp), &
            b0_t(ntemp)
TYPE(poly1) :: ec_p1(6,6)
TYPE(poly2) :: ec_p2(6,6)
TYPE(poly3) :: ec_p3(6,6)
TYPE(poly4) :: ec_p4(6,6)

REAL(DP), ALLOCATABLE :: xfit(:)

INTEGER  :: i, j, itemp, startt, lastt
REAL(DP) :: aux, macro_el(8)

ALLOCATE(xfit(nvar))
CALL divide(world_comm, ntemp, startt, lastt)
el_cons_t=0.0_DP
el_comp_t=0.0_DP
b0_t=0.0_DP
DO itemp=startt,lastt
   CALL compress_celldm(celldm_t(:,itemp),xfit,nvar,ibrav)
   DO i=1,6
      DO j=1,6
         IF (el_con_geo(i,j,1)>0.1_DP) THEN
            IF (poly_degree_elc==4) THEN
               CALL evaluate_fit_quartic(nvar,xfit,aux,ec_p4(i,j)) 
            ELSEIF (poly_degree_elc==3) THEN
               CALL evaluate_fit_cubic(nvar,xfit,aux,ec_p3(i,j))
            ELSEIF (poly_degree_elc==2) THEN
               CALL evaluate_fit_quadratic(nvar,xfit,aux,ec_p2(i,j))
            ELSEIF (poly_degree_elc==1) THEN
               CALL evaluate_fit_linear(nvar,xfit,aux,ec_p1(i,j))
            ELSE
               CALL errore('interpolate_el_cons','wrong poly_degree_elc',1)
            ENDIF
            el_cons_t(i,j,itemp)=aux
         ENDIF
      ENDDO
   ENDDO
   CALL compute_elastic_compliances(el_cons_t(:,:,itemp), el_comp_t(:,:,itemp))
   CALL print_macro_elasticity(ibrav,el_cons_t(:,:,itemp), &
                          el_comp_t(:,:,itemp),macro_el,.FALSE.)
   b0_t(itemp)=macro_el(5)
ENDDO
CALL mp_sum(el_cons_t, world_comm)
CALL mp_sum(el_comp_t, world_comm)
CALL mp_sum(b0_t, world_comm)

DEALLOCATE(xfit)

RETURN
END SUBROUTINE interpolate_el_cons
!
!----------------------------------------------------------------------------
SUBROUTINE interpolate_el_cons_p(celldm_p, nvar, ibrav, ec_p1, ec_p2, &
           ec_p3, ec_p4, poly_degree_elc, el_cons_p, el_comp_p, macro_el_p)
!----------------------------------------------------------------------------
!
! This routine receives as input the coeffients of polynomials which 
! interpolate the elastic constants on a grid of geometries,  
! the celldm_p that correspond to each pressure and gives as output
! the elastic constants, the elastic compliances and the bulk modulus
! interpolated at the celldm_p that corresponds to each pressure.
!  
USE kinds,              ONLY : DP
USE control_pressure,        ONLY : npress
USE elastic_constants,  ONLY : print_macro_elasticity, &
                              compute_elastic_compliances, el_con
USE control_elastic_constants, ONLY : el_con_geo
USE linear_surfaces,    ONLY : evaluate_fit_linear
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic
USE cubic_surfaces,     ONLY : evaluate_fit_cubic
USE quartic_surfaces,   ONLY : evaluate_fit_quartic
USE lattices,           ONLY : compress_celldm
USE polynomial,         ONLY : poly1, poly2, poly3, poly4
USE mp_world,           ONLY : world_comm
USE mp,                 ONLY : mp_sum

IMPLICIT NONE
INTEGER :: ibrav, nvar
INTEGER :: poly_degree_elc
REAL(DP) :: celldm_p(6, npress), el_cons_p(6,6,npress), &
            el_comp_p(6,6,npress), macro_el_p(8,npress)
TYPE(poly1) :: ec_p1(6,6)
TYPE(poly2) :: ec_p2(6,6)
TYPE(poly3) :: ec_p3(6,6)
TYPE(poly4) :: ec_p4(6,6)

REAL(DP), ALLOCATABLE :: xfit(:)

INTEGER  :: i, j, ipress, startp, lastp
REAL(DP) :: aux

ALLOCATE(xfit(nvar))
CALL divide(world_comm, npress, startp, lastp)
el_cons_p=0.0_DP
el_comp_p=0.0_DP
macro_el_p=0.0_DP
DO ipress=startp,lastp
   CALL compress_celldm(celldm_p(:,ipress),xfit,nvar,ibrav)
   DO i=1,6
      DO j=1,6
         IF (el_con_geo(i,j,1)>0.1_DP) THEN
            IF (poly_degree_elc==4) THEN
               CALL evaluate_fit_quartic(nvar,xfit,aux,ec_p4(i,j)) 
            ELSEIF (poly_degree_elc==3) THEN
               CALL evaluate_fit_cubic(nvar,xfit,aux,ec_p3(i,j))
            ELSEIF (poly_degree_elc==2) THEN
               CALL evaluate_fit_quadratic(nvar,xfit,aux,ec_p2(i,j))
            ELSEIF (poly_degree_elc==1) THEN
               CALL evaluate_fit_linear(nvar,xfit,aux,ec_p1(i,j))
            ELSE
               CALL errore('interpolate_el_cons','wrong poly_degree_elc',1)
            ENDIF
            el_cons_p(i,j,ipress)=aux
         ENDIF
      ENDDO
   ENDDO
   CALL compute_elastic_compliances(el_cons_p(:,:,ipress), &
                                    el_comp_p(:,:,ipress))
   CALL print_macro_elasticity(ibrav,el_cons_p(:,:,ipress), &
              el_comp_p(:,:,ipress),macro_el_p(:,ipress),.FALSE.)
ENDDO
CALL mp_sum(el_cons_p, world_comm)
CALL mp_sum(el_comp_p, world_comm)
CALL mp_sum(macro_el_p, world_comm)

DEALLOCATE(xfit)

RETURN
END SUBROUTINE interpolate_el_cons_p
