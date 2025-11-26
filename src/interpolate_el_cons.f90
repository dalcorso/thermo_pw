!
! Copyright (C) 2018 - 2025 Andrea Dal Corso
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
USE control_thermo,     ONLY : lgeo_from_file
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
REAL(DP) :: compute_omega_geo
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
   IF (itemp==1.OR.itemp==ntemp) CYCLE
   IF (lgeo_from_file) THEN
      xfit(1)=compute_omega_geo(ibrav, celldm_t(:,itemp))
   ELSE
      CALL compress_celldm(celldm_t(:,itemp),xfit,nvar,ibrav)
   ENDIF
   DO i=1,6
      DO j=1,6
         IF (ABS(el_con_geo(i,j,1))>0.1_DP) THEN
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
USE control_thermo,     ONLY : lgeo_from_file
USE lattices,           ONLY : compress_celldm
USE polynomial,         ONLY : poly1, poly2, poly3, poly4
USE mp_world,           ONLY : world_comm
USE mp,                 ONLY : mp_sum

IMPLICIT NONE
INTEGER :: ibrav, nvar
INTEGER :: poly_degree_elc
REAL(DP) :: celldm_p(6, npress), el_cons_p(6,6,npress), &
            el_comp_p(6,6,npress), macro_el_p(8,npress)
REAL(DP) :: compute_omega_geo
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
   IF (ipress==1.OR.ipress==npress) CYCLE
   IF (lgeo_from_file) THEN
      xfit(1)=compute_omega_geo(ibrav, celldm_p(:,ipress))
   ELSE
      CALL compress_celldm(celldm_p(:,ipress),xfit,nvar,ibrav)
   ENDIF
   DO i=1,6
      DO j=1,6
         IF (ABS(el_con_geo(i,j,1))>0.1_DP) THEN
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
!
!----------------------------------------------------------------------------
SUBROUTINE interpolate_piezo_tensor(celldm_t, nvar, ibrav, pt_p1, pt_p2, &
                    pt_p3, pt_p4, poly_degree_elc, e_piezo_tensor_t)
!----------------------------------------------------------------------------
!
! This routine receives as input the coeffients of polynomials which 
! interpolate the piezoelectric tensor on a grid of geometries,  
! the celldm that correspond to each temperature and gives as output
! the piezoelectric tensor interpolated at the celldm that corresponds 
! to each temperature.
!  
USE kinds,              ONLY : DP
USE temperature,        ONLY : ntemp
USE control_piezoelectric_tensor, ONLY : e_piezo_tensor_geo
USE control_thermo,     ONLY : lgeo_from_file
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
REAL(DP) :: celldm_t(6, ntemp), e_piezo_tensor_t(3,6,ntemp)
REAL(DP) :: compute_omega_geo
TYPE(poly1) :: pt_p1(3,6)
TYPE(poly2) :: pt_p2(3,6)
TYPE(poly3) :: pt_p3(3,6)
TYPE(poly4) :: pt_p4(3,6)

REAL(DP), ALLOCATABLE :: xfit(:)

INTEGER  :: i, j, itemp, startt, lastt
REAL(DP) :: aux, macro_el(8)

ALLOCATE(xfit(nvar))
CALL divide(world_comm, ntemp, startt, lastt)
e_piezo_tensor_t=0.0_DP
DO itemp=startt,lastt
   IF (itemp==1.OR.itemp==ntemp) CYCLE
   IF (lgeo_from_file) THEN
      xfit(1)=compute_omega_geo(ibrav, celldm_t(:,itemp))
   ELSE
      CALL compress_celldm(celldm_t(:,itemp),xfit,nvar,ibrav)
   ENDIF
   DO i=1,3
      DO j=1,6
         IF (ABS(e_piezo_tensor_geo(i,j,1))>1.D-6) THEN
            IF (poly_degree_elc==4) THEN
               CALL evaluate_fit_quartic(nvar,xfit,aux,pt_p4(i,j)) 
            ELSEIF (poly_degree_elc==3) THEN
               CALL evaluate_fit_cubic(nvar,xfit,aux,pt_p3(i,j))
            ELSEIF (poly_degree_elc==2) THEN
               CALL evaluate_fit_quadratic(nvar,xfit,aux,pt_p2(i,j))
            ELSEIF (poly_degree_elc==1) THEN
               CALL evaluate_fit_linear(nvar,xfit,aux,pt_p1(i,j))
            ELSE
               CALL errore('interpolate_piezo_tensor',&
                                           'wrong poly_degree_elc',1)
            ENDIF
            e_piezo_tensor_t(i,j,itemp)=aux
         ENDIF
      ENDDO
   ENDDO
ENDDO
CALL mp_sum(e_piezo_tensor_t, world_comm)

DEALLOCATE(xfit)

RETURN
END SUBROUTINE interpolate_piezo_tensor
!
!----------------------------------------------------------------------------
SUBROUTINE interpolate_piezo_tensor_p(celldm_p, nvar, ibrav, pt_p1, pt_p2, &
           pt_p3, pt_p4, poly_degree_elc, e_piezo_tensor_p )
!----------------------------------------------------------------------------
!
! This routine receives as input the coeffients of polynomials which 
! interpolate the piezoelectric tensor on a grid of geometries,  
! the celldm_p that corresponds to each pressure and gives as output
! the piezoelectric tensor interpolated at the celldm_p that corresponds 
! to each pressure.
!  
USE kinds,              ONLY : DP
USE control_pressure,   ONLY : npress
USE control_piezoelectric_tensor, ONLY : e_piezo_tensor_geo

USE linear_surfaces,    ONLY : evaluate_fit_linear
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic
USE cubic_surfaces,     ONLY : evaluate_fit_cubic
USE quartic_surfaces,   ONLY : evaluate_fit_quartic

USE control_thermo,     ONLY : lgeo_from_file
USE lattices,           ONLY : compress_celldm
USE polynomial,         ONLY : poly1, poly2, poly3, poly4
USE mp_world,           ONLY : world_comm
USE mp,                 ONLY : mp_sum

IMPLICIT NONE
INTEGER :: ibrav, nvar
INTEGER :: poly_degree_elc
REAL(DP) :: celldm_p(6, npress), e_piezo_tensor_p(3,6,npress)
REAL(DP) :: compute_omega_geo
TYPE(poly1) :: pt_p1(3,6)
TYPE(poly2) :: pt_p2(3,6)
TYPE(poly3) :: pt_p3(3,6)
TYPE(poly4) :: pt_p4(3,6)

REAL(DP), ALLOCATABLE :: xfit(:)

INTEGER  :: i, j, ipress, startp, lastp
REAL(DP) :: aux

ALLOCATE(xfit(nvar))
CALL divide(world_comm, npress, startp, lastp)
e_piezo_tensor_p=0.0_DP
DO ipress=startp,lastp
   IF (ipress==1.OR.ipress==npress) CYCLE
   IF (lgeo_from_file) THEN
      xfit(1)=compute_omega_geo(ibrav, celldm_p(:,ipress))
   ELSE
      CALL compress_celldm(celldm_p(:,ipress),xfit,nvar,ibrav)
   ENDIF
   DO i=1,3
      DO j=1,6
         IF (ABS(e_piezo_tensor_geo(i,j,1))>1.D-6) THEN
            IF (poly_degree_elc==4) THEN
               CALL evaluate_fit_quartic(nvar,xfit,aux,pt_p4(i,j)) 
            ELSEIF (poly_degree_elc==3) THEN
               CALL evaluate_fit_cubic(nvar,xfit,aux,pt_p3(i,j))
            ELSEIF (poly_degree_elc==2) THEN
               CALL evaluate_fit_quadratic(nvar,xfit,aux,pt_p2(i,j))
            ELSEIF (poly_degree_elc==1) THEN
               CALL evaluate_fit_linear(nvar,xfit,aux,pt_p1(i,j))
            ELSE
               CALL errore('interpolate_piezo_tensor_p',&
                                          'wrong poly_degree_elc',1)
            ENDIF
            e_piezo_tensor_p(i,j,ipress)=aux
         ENDIF
      ENDDO
   ENDDO
ENDDO
CALL mp_sum(e_piezo_tensor_p, world_comm)

DEALLOCATE(xfit)

RETURN
END SUBROUTINE interpolate_piezo_tensor_p

!----------------------------------------------------------------------------
SUBROUTINE compute_d_piezo_tensor_t(e_piezo_tensor_t,el_comp_t,&
                                        d_piezo_tensor_t,ntemp)
!----------------------------------------------------------------------------

USE kinds,              ONLY : DP
USE piezoelectric_tensor, ONLY : e_piezo_tensor, d_piezo_tensor, &
                                 compute_d_piezo_tensor

IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: e_piezo_tensor_t(3,6,ntemp), el_comp_t(6,6,ntemp)
REAL(DP), INTENT(OUT) :: d_piezo_tensor_t(3,6,ntemp)

INTEGER :: itemp

DO itemp=1,ntemp
   e_piezo_tensor(:,:)=e_piezo_tensor_t(:,:,itemp)
   CALL compute_d_piezo_tensor(el_comp_t(1,1,itemp))
   d_piezo_tensor_t(:,:,itemp)=d_piezo_tensor(:,:)
ENDDO

RETURN
END SUBROUTINE compute_d_piezo_tensor_t
