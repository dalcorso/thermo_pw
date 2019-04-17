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
USE thermo_mod,     ONLY : omega_geo, tot_ngeo, no_ph
USE temperature,    ONLY : ntemp, temp
USE io_global,      ONLY : stdout
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic, &
                        introduce_quadratic_fit, quadratic_var
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic,   &
                        introduce_quartic_fit, quartic_var,             &
                        print_chisq_quartic
USE cubic_surfaces, ONLY : fit_multi_cubic, evaluate_fit_cubic, cubic_var, &
                        introduce_cubic_fit
USE control_quartic_energy, ONLY : lsolve, poly_degree_bfact
IMPLICIT NONE
REAL(DP), INTENT(IN)  :: vmin_t(ntemp), ph_b_fact(3,3,nat,ntemp,tot_ngeo)
REAL(DP), INTENT(OUT) :: bfact_t(6,nat,ntemp)
REAL(DP), ALLOCATABLE :: x(:), b(:), coeff(:)
REAL(DP) :: asum
INTEGER :: itemp, startt, lastt
INTEGER :: na, ipol, jpol, ijpol, igeo, ndata, degree, nvar
INTEGER :: compute_nwork_ph
!
!  in this routine the B factor is given as a function of volume so the
!  polynomial has one variable.
!
degree=1
ndata = compute_nwork_ph(no_ph,tot_ngeo)
IF (poly_degree_bfact==4) THEN
   nvar=quartic_var(degree)
ELSEIF (poly_degree_bfact==3) THEN
   nvar=cubic_var(degree)
ELSE
   nvar=quadratic_var(degree)
ENDIF
ALLOCATE(coeff(nvar))
ALLOCATE(x(ndata))
ALLOCATE(b(ndata))

WRITE(stdout,'(5x,"B FACTOR INTERPOLATION")')
!
!  collect the volumes of the geometries for which phonons have been 
!  calculated
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   ndata=ndata+1
   x(ndata)=omega_geo(igeo)
ENDDO

CALL divide(world_comm, ntemp, startt, lastt)
bfact_t=0.0_DP
DO itemp=startt,lastt
!   WRITE(stdout,'(5x,70("_"))')
!   WRITE(stdout,'(/,5x,"Temperature:",f12.3," K")') temp(itemp)
   DO na=1, nat
      ijpol=0
      DO ipol=1,3
         DO jpol=ipol,3
            ijpol=ijpol+1
            asum=0.0_DP
            ndata=0
            DO igeo=1, tot_ngeo
               IF (no_ph(igeo)) CYCLE
               ndata=ndata+1
               b(ndata)=ph_b_fact(ipol,jpol,na,itemp,igeo)
               asum=asum+b(ndata)**2
            ENDDO
            IF (asum<1.D-10) CYCLE
!            WRITE(stdout,'(5x,"atom n.",i4,",   B component: ",i1)') na, ijpol

            IF (poly_degree_bfact==4) THEN
               IF (itemp==startt) CALL introduce_quartic_fit(degree, nvar, &
                                                              ndata)
               CALL fit_multi_quartic(ndata, degree, nvar, lsolve, x, b,   &
                                                                     coeff)
               CALL evaluate_fit_quartic(degree, nvar, vmin_t(itemp),      &
                                            bfact_t(ijpol,na,itemp), coeff)
            ELSEIF(poly_degree_bfact==3) THEN
               IF (itemp==startt) CALL introduce_cubic_fit(degree, nvar, &
                                                              ndata)
               CALL fit_multi_cubic(ndata, degree, nvar, lsolve, x, b, coeff)
               CALL evaluate_fit_cubic(degree, nvar, vmin_t(itemp),      &
                                            bfact_t(ijpol,na,itemp), coeff)
            ELSE
               IF (itemp==startt) CALL introduce_quadratic_fit(degree, nvar, &
                                                              ndata)
               CALL fit_multi_quadratic(ndata, degree, nvar, x, b, coeff)
               CALL evaluate_fit_quadratic(degree, nvar, vmin_t(itemp),      &
                                            bfact_t(ijpol,na,itemp), coeff)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

CALL mp_sum(bfact_t, world_comm)

DEALLOCATE(x)
DEALLOCATE(b)
DEALLOCATE(coeff)

RETURN
END SUBROUTINE interpolate_b_fact

SUBROUTINE interpolate_b_fact_anis(celldm_t, ph_b_fact, bfact_t)
!
!  This subroutine receives the 6 independent components of
!  B factor (BF) matrix for each atom in the unit cell
!  for several geometries and the equilibrium lattice parameters 
!  as a function of temperature. Then it interpolates each component 
!  of the BF matrices at the equilibrium  lattice parameters at each
!  temperature.
!
USE kinds,          ONLY : DP
USE cell_base,      ONLY : ibrav
USE ions_base,      ONLY : nat
USE thermo_mod,     ONLY : tot_ngeo, celldm_geo, no_ph
USE temperature,    ONLY : ntemp
USE io_global,      ONLY : stdout
USE lattices,       ONLY : compress_celldm, crystal_parameters
USE quadratic_surfaces, ONLY : fit_multi_quadratic, evaluate_fit_quadratic, &
                           introduce_quadratic_fit, quadratic_var
USE quartic_surfaces, ONLY : fit_multi_quartic, evaluate_fit_quartic, &
                           introduce_quartic_fit, quartic_var
USE cubic_surfaces,   ONLY : fit_multi_cubic, evaluate_fit_cubic, &
                             introduce_cubic_fit, cubic_var
USE control_quartic_energy, ONLY : lsolve, poly_degree_bfact
USE mp_world,       ONLY : world_comm
USE mp,             ONLY : mp_sum

IMPLICIT NONE
REAL(DP), INTENT(IN)  :: celldm_t(6,ntemp), ph_b_fact(3,3,nat,ntemp,tot_ngeo)
REAL(DP), INTENT(OUT) :: bfact_t(6,nat,ntemp)
REAL(DP), ALLOCATABLE :: b(:), coeff(:), x(:,:), y(:,:)
REAL(DP) :: asum
INTEGER :: itemp, startt, lastt
INTEGER :: na, ipol, jpol, ijpol, igeo, ndata, degree, nvar
INTEGER :: compute_nwork_ph

degree=crystal_parameters(ibrav)
ndata = compute_nwork_ph(no_ph,tot_ngeo)
IF (poly_degree_bfact==4) THEN
   nvar=quartic_var(degree)
ELSEIF (poly_degree_bfact==3) THEN
   nvar=cubic_var(degree)
ELSE
   nvar=quadratic_var(degree)
ENDIF

ALLOCATE(coeff(nvar))
ALLOCATE(b(ndata))
ALLOCATE(x(degree,ndata))
ALLOCATE(y(degree,1))

WRITE(stdout,'(5x,"B FACTOR INTERPOLATION")')
!
!  collect the geometric information of the geometries for which
!  phonons have been calculated
!
ndata=0
DO igeo=1, tot_ngeo
   IF (no_ph(igeo)) CYCLE
   ndata=ndata+1
   CALL compress_celldm(celldm_geo(1,igeo), x(1,ndata), degree, ibrav)
ENDDO

CALL divide(world_comm, ntemp, startt, lastt)
bfact_t=0.0_DP
DO itemp=startt,lastt
!   WRITE(stdout,'(5x,70("_"))')
!   WRITE(stdout,'(/,5x,"Temperature:",f12.3," K")') temp(itemp)
   CALL compress_celldm(celldm_t(1,itemp), y, degree, ibrav)
   DO na=1, nat
      ijpol=0
      DO ipol=1,3
         DO jpol=ipol,3
            ijpol=ijpol+1
            asum=0.0_DP
            ndata=0
            DO igeo=1, tot_ngeo
               IF (no_ph(igeo)) CYCLE
               ndata=ndata+1
               b(ndata)=ph_b_fact(ipol,jpol,na,itemp,igeo)
               asum=asum+b(ndata)**2
            ENDDO
            IF (asum<1.D-10) CYCLE
!            WRITE(stdout,'(5x,"Atom n.",i4,",   B component: ",i1)') na, ijpol
            
            IF (poly_degree_bfact==4) THEN
               IF (itemp==startt) CALL introduce_quartic_fit(degree, nvar, &
                                                                     ndata)
               CALL fit_multi_quartic(ndata, degree, nvar, lsolve, x, b, coeff)
               CALL evaluate_fit_quartic(degree, nvar, y, &
                                         bfact_t(ijpol,na,itemp), coeff)

            ELSEIF (poly_degree_bfact==3) THEN
               IF (itemp==startt) CALL introduce_cubic_fit(degree, nvar, ndata)
               CALL fit_multi_cubic(ndata, degree, nvar, lsolve, x, b, coeff)
               CALL evaluate_fit_cubic(degree, nvar, y, &
                                         bfact_t(ijpol,na,itemp), coeff)
            ELSE
               IF (itemp==startt) CALL introduce_quadratic_fit(degree, nvar, &
                                                                     ndata)
               CALL fit_multi_quadratic(ndata, degree, nvar, x, b, coeff)
               CALL evaluate_fit_quadratic(degree, nvar, y, &
                                         bfact_t(ijpol,na,itemp), coeff)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

CALL mp_sum(bfact_t, world_comm)

DEALLOCATE(b)
DEALLOCATE(coeff)
DEALLOCATE(x)
DEALLOCATE(y)

RETURN
END SUBROUTINE interpolate_b_fact_anis
