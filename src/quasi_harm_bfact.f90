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
!  matrices at the equilibrium volume at each temperature.
!
USE kinds,          ONLY : DP
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : ngeo, omega_geo, tot_ngeo
USE temperature,    ONLY : ntemp, temp
USE io_global,      ONLY : stdout
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic,  &
                        introduce_quadratic_fit, print_quadratic_polynomial, &
                        summarize_fitting_data, print_chisq_quadratic
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic,   &
                        introduce_quartic_fit, compute_quartic_var,     &
                        print_quartic_polynomial, evaluate_fit_quartic, &
                        print_chisq_quartic
IMPLICIT NONE
REAL(DP), INTENT(IN)  :: vmin_t(ntemp), ph_b_fact(3,3,nat,ntemp,tot_ngeo)
REAL(DP), INTENT(OUT) :: bfact_t(6,nat,ntemp)
REAL(DP), ALLOCATABLE :: b(:), coeff(:), coeff4(:)
REAL(DP) :: asum
INTEGER :: itemp, startt, lastt
INTEGER :: na, ipol, jpol, ijpol, idata, ndata, degree, nvar, nvar4, lsolve
INTEGER :: compute_nwork
!
!  in this routine the B factor is given as a function of volume so the
!  polynomial has one variable.
!
degree=1
nvar=3
nvar4=compute_quartic_var(degree)
ndata = compute_nwork()
lsolve=2

ALLOCATE(b(ndata))
ALLOCATE(coeff(nvar))
ALLOCATE(coeff4(nvar4))

WRITE(stdout,'(5x,"B FACTOR INTERPOLATION")')
CALL divide(world_comm, ntemp, startt, lastt)
bfact_t=0.0_DP
DO itemp=startt,lastt
   WRITE(stdout,'(5x,70("_"))')
   WRITE(stdout,'(/,5x,"Temperature:",f12.3," K")') temp(itemp)
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
            WRITE(stdout,'(5x,"atom n.",i4,",   B component: ",i1)') na, ijpol

            IF (itemp==startt) CALL introduce_quadratic_fit(degree, nvar, ndata)
!            CALL summarize_fitting_data(degree, ndata, omega_geo, b)
            CALL fit_multi_quadratic(ndata, degree, nvar, omega_geo, b, coeff)
!            CALL print_quadratic_polynomial(degree, nvar, coeff)
            CALL print_chisq_quadratic(ndata, degree, nvar, omega_geo, b, coeff)
!            CALL evaluate_fit_quadratic(degree, nvar, vmin_t(itemp), &
!                                        bfact_t(ijpol,na,itemp), coeff)
            IF (itemp==startt) CALL introduce_quartic_fit(degree, nvar4, ndata)
            CALL fit_multi_quartic(ndata, degree, nvar4, lsolve, &
                                                  omega_geo, b, coeff4)
            !CALL print_quartic_polynomial(degree, nvar4, coeff4)
            CALL evaluate_fit_quartic(degree, nvar4, vmin_t(itemp), &
                             bfact_t(ijpol,na,itemp), coeff4)
            CALL print_chisq_quartic(ndata, degree, nvar4, omega_geo, &
                                                            b, coeff4)
         ENDDO
      ENDDO
   ENDDO
ENDDO

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
!  for several geometries and the equilibrium lattice parameters 
!  as a function of temperature. Then it interpolates each component 
!  of BF matrices at the equilibrium  lattice parameters at each
!  temperature.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : ngeo, tot_ngeo, celldm_geo
USE temperature,    ONLY : ntemp, temp
USE io_global,      ONLY : stdout
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic, &
                    introduce_quadratic_fit, summarize_fitting_data,        &
                    print_quadratic_polynomial, print_chisq_quadratic, &
                    quadratic_var
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic, &
                    introduce_quartic_fit, compute_quartic_var,       &
                    print_quartic_polynomial, print_chisq_quartic 

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: celldm_t(6,ntemp), ph_b_fact(3,3,nat,ntemp,tot_ngeo)
REAL(DP), INTENT(OUT) :: bfact_t(6,nat,ntemp)
REAL(DP), ALLOCATABLE :: b(:), coeff(:), coeff4(:), x(:,:), y(:,:)
REAL(DP) :: asum
INTEGER :: itemp, startt, lastt
INTEGER :: na, ipol, jpol, ijpol, idata, ndata, degree, nvar, nvar4, lsolve
INTEGER :: compute_nwork

degree=crystal_parameters(ibrav)
nvar=quadratic_var(degree)
nvar4=compute_quartic_var(degree)
ndata = compute_nwork()
lsolve=2

ALLOCATE(b(ndata))
ALLOCATE(coeff(nvar))
ALLOCATE(coeff4(nvar4))
ALLOCATE(x(degree,ndata))
ALLOCATE(y(degree,1))

WRITE(stdout,'(5x,"B FACTOR INTERPOLATION")')

CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_geo)

CALL divide(world_comm, ntemp, startt, lastt)
bfact_t=0.0_DP
DO itemp=startt,lastt
   WRITE(stdout,'(5x,70("_"))')
   WRITE(stdout,'(/,5x,"Temperature:",f12.3," K")') temp(itemp)
   CALL compress_celldm(celldm_t(1,itemp), y, degree, ibrav)
   DO na=1, nat
      ijpol=0
      DO ipol=1,3
         DO jpol=ipol,3
            ijpol=ijpol+1
            asum=0.0_DP
            DO idata=1, ndata
               b(idata)=ph_b_fact(ipol,jpol,na,itemp,idata)
               asum=asum+b(idata)**2
            ENDDO
            IF (asum<1.D-10) CYCLE
            WRITE(stdout,'(5x,"Atom n.",i4,",   B component: ",i1)') na, ijpol
            
            IF (itemp==startt) CALL introduce_quadratic_fit(degree, nvar, ndata)
!            CALL summarize_fitting_data(degree, ndata, x, b)
            CALL fit_multi_quadratic(ndata, degree, nvar, x, b, coeff)
            CALL print_quadratic_polynomial(degree, nvar, coeff)
            CALL print_chisq_quadratic(ndata, degree, nvar, x, b, coeff)

!            CALL evaluate_fit_quadratic(degree, nvar, y, &
!                                           bfact_t(ijpol,na,itemp), coeff)
            IF (itemp==startt) CALL introduce_quartic_fit(degree, nvar4, ndata)
            CALL fit_multi_quartic(ndata, degree, nvar4, lsolve, &
                                                            x, b, coeff4)
            CALL print_quartic_polynomial(degree, nvar4, coeff4)
            CALL evaluate_fit_quartic(degree, nvar4, y, &
                                         bfact_t(ijpol,na,itemp), coeff4)
            CALL print_chisq_quadratic(ndata, degree, nvar4, x, b, coeff)
         ENDDO
      ENDDO
   ENDDO
ENDDO

CALL mp_sum(bfact_t, world_comm)

DEALLOCATE(b)
DEALLOCATE(coeff)
DEALLOCATE(coeff4)
DEALLOCATE(x)
DEALLOCATE(y)

RETURN
END SUBROUTINE interpolate_b_fact_anis
