!
! Copyright (C) 2018 Cristiano Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE interpolate_b_fact(vmin_t, ph_b_fact, bfact_t)
!
!  This subroutine receives the 6 independent components of
!  B factor (BF) matrix for each atom in the unit cell
!  for several volumes and the equilibrium volume as a function
!  of temperature. Then it interpolates each component of BF
!  matrices at the volume that corresponds to a given temperature.
!
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph, tot_ngeo
USE temperature,    ONLY : ntemp, temp
USE io_global,      ONLY : stdout
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic, &
                               introduce_quadratic_fit
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic, &
                             introduce_quartic_fit, compute_quartic_var


IMPLICIT NONE
REAL(DP), INTENT(IN)  :: vmin_t(ntemp), ph_b_fact(3,3,nat,ntemp,tot_ngeo)
REAL(DP), INTENT(INOUT) :: bfact_t(6,nat,ntemp)
REAL(DP), ALLOCATABLE :: b(:), coeff(:), coeff4(:)
REAL(DP) :: asum, aux, chisq
INTEGER :: itemp, startt, lastt
INTEGER :: na, ipol, jpol, ijpol, idata, ndata, degree, nvar, nvar4, lsolve
INTEGER :: compute_nwork

CALL compute_degree(ibrav, degree, nvar)

nvar4=compute_quartic_var(degree)

ndata = compute_nwork()

lsolve=2

ALLOCATE(b(ndata))
ALLOCATE(coeff(nvar))
ALLOCATE(coeff4(nvar4))

WRITE(stdout,'(/,5x,"B FACTOR INTERPOLATION")')

bfact_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp=startt,lastt
   
   DO na=1, nat
      ijpol=0
      DO ipol=1,3
         DO jpol=ipol,3
            ijpol=ijpol+1
               
               asum=0.0_DP
               DO idata=1, ndata
                  b(idata)=ph_b_fact(ipol,jpol,na,itemp,idata)
                  asum=asum+b(idata)**2
               END DO
               IF (asum<1.D-10) CYCLE
           
               WRITE(stdout,'(/,5x,"Temperature:",f12.3," K,   atom n.",i4,",   B component: ",i1)') temp(itemp), na, ijpol     
           
               CALL introduce_quadratic_fit(degree, nvar, ndata)
               CALL fit_multi_quadratic(ndata, degree, nvar, omega_geo, b, coeff)

               chisq=0.0_DP
               DO idata=1,ndata
                  CALL evaluate_fit_quadratic(degree, nvar, omega_geo(idata), aux, coeff)
                  chisq = chisq + (aux - b(idata))**2
               ENDDO
               WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq
               
               CALL evaluate_fit_quadratic(degree, nvar, vmin_t(itemp), bfact_t(ijpol,na,itemp), coeff)

               !CALL introduce_quartic_fit(degree, nvar4, ndata)

               CALL fit_multi_quartic(ndata, degree, nvar4, lsolve, omega_geo, b, coeff4)

               chisq=0.0_DP
               DO idata=1,ndata
                  CALL evaluate_fit_quartic(degree, nvar4, omega_geo(idata), aux, coeff4)
                  chisq = chisq + (aux - b(idata))**2
               ENDDO
                  
               WRITE(stdout,'(/,5x,"chi square of the quartic fit =",e18.5,/)') chisq
         END DO
      END DO
   END DO

END DO

CALL mp_sum(bfact_t, world_comm)

DEALLOCATE(b)
DEALLOCATE(coeff)
DEALLOCATE(coeff4)


RETURN

END SUBROUTINE interpolate_b_fact

SUBROUTINE interpolate_b_fact_anis(celldm_t, ph_b_fact, bfact_t)
!
!  This subroutine receives the 6 independent components of
!  B factor (BF) matrix for each atom in the unit cell
!  for several volumes and the equilibrium lattice parameters 
!  as a function of temperature. Then it interpolates each component 
!  of BF matrices at the lattice parameters that corresponds to 
!  a given temperature.
!
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : ngeo, omega_geo, no_ph, tot_ngeo, celldm_geo
USE temperature,    ONLY : ntemp, temp
USE io_global,      ONLY : stdout
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic, &
                               introduce_quadratic_fit
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic, &
                             introduce_quartic_fit, compute_quartic_var

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: celldm_t(6,ntemp), ph_b_fact(3,3,nat,ntemp,tot_ngeo)
REAL(DP), INTENT(OUT) :: bfact_t(6,nat,ntemp)
REAL(DP), ALLOCATABLE :: b(:), coeff(:), coeff4(:), x(:,:), y(:,:)
REAL(DP) :: asum, aux, chisq
INTEGER :: itemp, startt, lastt
INTEGER :: na, ipol, jpol, ijpol, idata, ndata, degree, nvar, nvar4, lsolve
INTEGER :: compute_nwork

CALL compute_degree(ibrav, degree, nvar)

nvar4=compute_quartic_var(degree)

ndata = compute_nwork()

lsolve=2

ALLOCATE(b(ndata))
ALLOCATE(coeff(nvar))
ALLOCATE(coeff4(nvar4))
ALLOCATE(x(degree,ndata))
ALLOCATE(y(degree,ndata))

WRITE(stdout,'(/,5x,"B FACTOR INTERPOLATION")')

CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_geo)

bfact_t=0.0_DP

CALL divide(world_comm, ntemp, startt, lastt)

DO itemp=startt,lastt
   CALL set_x_from_celldm(ibrav, degree, 1, y, celldm_t(1,itemp))

   DO na=1, nat
      ijpol=0
      DO ipol=1,3
         DO jpol=ipol,3
            ijpol=ijpol+1

               asum=0.0_DP
               DO idata=1, ndata
                  b(idata)=ph_b_fact(ipol,jpol,na,itemp,idata)
                  asum=asum+b(idata)**2
               END DO
               IF (asum<1.D-10) CYCLE
               
               WRITE(stdout,'(/,5x,"Temperature:",f12.3," K, atom n.",i4,",   B component: ",i1)') temp(itemp), na, ijpol
               
               CALL introduce_quadratic_fit(degree, nvar, ndata)
               CALL fit_multi_quadratic(ndata, degree, nvar, x, b, coeff)

               chisq=0.0_DP
               DO idata=1,ndata
                  CALL evaluate_fit_quadratic(degree, nvar, x(1,idata), aux, coeff)
                  chisq = chisq + (aux - b(idata))**2
               ENDDO
               WRITE(stdout,'(5x,"chi square=",e18.5,/)') chisq

               CALL evaluate_fit_quadratic(degree, nvar, y, bfact_t(ijpol,na,itemp), coeff)

               !CALL introduce_quartic_fit(degree, nvar4, ndata)
               CALL fit_multi_quartic(ndata, degree, nvar4, lsolve, x, b, coeff4)

               chisq=0.0_DP
               DO idata=1,ndata
                  CALL evaluate_fit_quartic(degree, nvar4, x(1,idata), aux, coeff4)
                  chisq = chisq + (aux - b(idata))**2
               ENDDO

               WRITE(stdout,'(/,5x,"chi square of the quartic fit =",e18.5,/)') chisq

         END DO
      END DO
   END DO

END DO

CALL mp_sum(bfact_t, world_comm)

DEALLOCATE(b)
DEALLOCATE(coeff)
DEALLOCATE(coeff4)
DEALLOCATE(x)
DEALLOCATE(y)

RETURN

END SUBROUTINE interpolate_b_fact_anis
