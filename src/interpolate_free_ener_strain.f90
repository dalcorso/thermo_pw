! Copyright (C) 2019 Cristiano Malica Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE interpolate_free_ener_strain(p1,p2,p3,p4,startt,lastt,ngroup)
!------------------------------------------------------------------------

USE kinds,             ONLY : DP
USE initial_conf,      ONLY : ibrav_save
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE thermo_mod,        ONLY : energy_geo, tot_ngeo
USE control_elastic_constants, ONLY : ngeo_strain, elcpvar, ngeom, &
                              work_base, el_con_celldm_geo, epsil_geo
USE thermodynamics,    ONLY : ph_free_ener
USE temperature,       ONLY : ntemp
USE control_quartic_energy, ONLY : lsolve, poly_degree_elc
USE linear_surfaces, ONLY : fit_multi_linear
USE quadratic_surfaces, ONLY : fit_multi_quadratic
USE cubic_surfaces, ONLY : fit_multi_cubic
USE quartic_surfaces, ONLY : fit_multi_quartic
USE polynomial, ONLY : poly1, poly2, poly3, poly4, init_poly, clean_poly
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
INTEGER :: startt, lastt, ngroup

TYPE(poly1) :: p1(startt:lastt,ngroup)
TYPE(poly2) :: p2(startt:lastt,ngroup)
TYPE(poly3) :: p3(startt:lastt,ngroup)
TYPE(poly4) :: p4(startt:lastt,ngroup)

INTEGER :: itemp, igroup, ndata, igeom, ind, nvar, istrain 

REAL(DP), ALLOCATABLE :: x(:,:), y(:,:), g(:)


nvar=crystal_parameters(ibrav_save)

ngroup = work_base/ngeo_strain

ALLOCATE(x(nvar,ngeom))
ALLOCATE(y(nvar+1, ngeom*ngeo_strain))
ALLOCATE(g(ngeom*ngeo_strain))

DO igeom=1, ngeom
   CALL compress_celldm(el_con_celldm_geo(1,igeom), x(1,igeom), nvar, ibrav_save)
ENDDO

DO itemp = startt, lastt
   DO igroup=1, ngroup
      ndata=0
      DO igeom=1, ngeom
         DO istrain=1, ngeo_strain
            ndata=ndata+1
            ind = istrain + (igeom-1)*work_base + &
                  (igroup-1)*ngeo_strain
            y(1:nvar,ndata) = x(1:nvar,igeom)
            y(nvar+1,ndata) = epsil_geo(ind)
            g(ndata) = energy_geo(ind) + ph_free_ener(itemp,ind)
         ENDDO
      ENDDO
      IF (poly_degree_elc==4) THEN
         CALL init_poly(nvar+1,p4(itemp,igroup))
         CALL fit_multi_quartic(ndata, nvar+1, lsolve, y, g, p4(itemp,igroup))
      ELSEIF (poly_degree_elc==3) THEN
         CALL init_poly(nvar+1,p3(itemp,igroup))
         CALL fit_multi_cubic(ndata, nvar+1, lsolve, y, g, p3(itemp,igroup))
      ELSEIF (poly_degree_elc==2) THEN
         CALL init_poly(nvar+1,p2(itemp,igroup))
         CALL fit_multi_quadratic(ndata, nvar+1, lsolve, y, g, p2(itemp,igroup)) 
      ELSEIF (poly_degree_elc==1) THEN
         CALL init_poly(nvar+1,p1(itemp,igroup))
         CALL fit_multi_linear(ndata, nvar+1, lsolve, y, g, p1(itemp,igroup))
      ELSE
         CALL errore('interpolate_free_ener_strain','wrong poly_degree',1)
      END IF
   ENDDO
ENDDO 
!
RETURN
END SUBROUTINE interpolate_free_ener_strain
