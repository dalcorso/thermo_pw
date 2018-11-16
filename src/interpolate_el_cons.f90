!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE interpolate_el_cons(celldm_t, nvar, degree, ibrav, el_cons_coeff,&
                                     lquartic, el_cons_t, el_comp_t, b0_t)

USE kinds,              ONLY : DP
USE temperature,        ONLY : ntemp
USE elastic_constants,  ONLY : print_macro_elasticity, &
                              compute_elastic_compliances, el_con
USE control_elastic_constants, ONLY : el_con_geo
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic
USE quartic_surfaces,   ONLY : evaluate_fit_quartic
USE lattices,           ONLY : compress_celldm
USE mp_world,           ONLY : world_comm
USE mp,                 ONLY : mp_sum

IMPLICIT NONE
INTEGER :: nvar, ibrav, degree
LOGICAL :: lquartic
REAL(DP) :: celldm_t(6, ntemp), el_cons_coeff(nvar,6,6), el_cons_t(6,6,ntemp),&
            el_comp_t(6,6,ntemp), b0_t(ntemp)

REAL(DP), ALLOCATABLE :: xfit(:)

INTEGER  :: i, j, itemp, startt, lastt
REAL(DP) :: aux, macro_el(8)

ALLOCATE(xfit(degree))
CALL divide(world_comm, ntemp, startt, lastt)
el_cons_t=0.0_DP
el_comp_t=0.0_DP
b0_t=0.0_DP
DO itemp=startt,lastt
   CALL compress_celldm(celldm_t(:,itemp),xfit,degree,ibrav)
   DO i=1,6
      DO j=i,6
         IF (el_con_geo(i,j,1)>0.1_DP) THEN
            IF (lquartic) THEN
               CALL evaluate_fit_quartic(degree,nvar,xfit,aux,&
                                                  el_cons_coeff(:,i,j))
            ELSE
               CALL evaluate_fit_quadratic(degree,nvar,xfit,aux,&
                                                el_cons_coeff(:,i,j))
            ENDIF
            el_cons_t(i,j,itemp)=aux
            el_cons_t(j,i,itemp)=aux
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
