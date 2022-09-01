!
! Copyright (C) 2021 Andrea Dal Corso and Xuejun Gong
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
PROGRAM test_eos
!--------------------------------------------------------------------
!
!  This program tests the module eos.f90 of the library.
!  It reads in input the parameters of an equation of state and
!  produces on output the energy, the pressure and bulk modulus
!  and the first and possibly the second derivative of the bulk
!  modulus with respect to pressure as a function of volume.
!  It checks these quantities with quantities computed by 
!  numerical finite differences.
!  It is also possible to add a polynomial of arbitrary degree to
!  the equation of state by givin in input its coefficients.
!  The input file has the following form
!  ieos       ! choose the equation of state 1 - Birch-Murnaghan  3 order
!                                            2 - Birch-Murnaghan  4 order
!                                            4 - Murnaghan
!  v0         ! the equilibrium volume in (a.u.)^3.
!  b0         ! the bulk modulus (in kbar).
!  b01        ! derivative of the bulk modulus with respect to pressure.
!  b02        ! second derivative of the bulk modulus with respect to 
!             ! pressure, set it to 0.0 for ioes 1 or 4. In (1/kbar).
!  vmax       ! maximum volume of the output functions in (a.u.)^3.
!  vmin       ! minimum volume of the output functions in (a.u.)^3.
!  npt        ! number of points in the grid of volumes.
!  'Yes'      ! to use a polynomial. Anything else do not use it.
!  degree     ! degree of the polynomial.
!  a(1)       ! coefficient of the polynomial for volume^0
!  a(2)       ! coefficients of the polynomial for volume^1
!  ...        !
!  a(degree+1)! coefficient of the polynomial for volume^degree
!
USE kinds, ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE constants,        ONLY : ry_kbar
USE environment,      ONLY : environment_start, environment_end
USE eos,              ONLY : eos_energy, eos_bulk, eos_press, eos_dpress, &
                             eos_energy_pol, eos_bulk_pol, eos_press_pol, &
                             eos_dpress_pol
USE io_global,        ONLY : stdin, stdout

IMPLICIT NONE
CHARACTER(LEN=9) :: code='EOS'
CHARACTER(LEN=3) :: polyout
REAL(DP) :: v0, b0, b01, b02, deltav, vmax, vmin
INTEGER :: i, j, ieos, degree, m1
INTEGER :: npt
LOGICAL :: lpoly
REAL(DP), ALLOCATABLE :: vm(:), eout(:), pout(:), dpout(:),     &
                         bulkout(:), dbulkout(:), d2bulkout(:), &
                         p_num(:), dp_num(:), bulk_num(:),      &
                         dbulk_num(:), d2bulk_num(:),a(:)

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'(/,5x,"Choose equation of state: (1: BM-3 2: BM-4 4: Murnaghan)")')
READ(stdin,*) ieos

WRITE(stdout,'(/,5x,"V0 (unit: (a.u.)^3)?")')
READ(stdin,*) v0

WRITE(stdout,'(/,5x,"B0 (unit: kbar)?")')
READ(stdin,*) b0

WRITE(stdout,'(/,5x,"dB/dP?")')
READ(stdin,*) b01

WRITE(stdout,'(/,5x,"d^2B/dP^2 (unit: 1/kbar)?")')
READ(stdin,*) b02

WRITE(stdout,'(/,5x,"Maximum volume (unit: (a.u.)^3)?")')
READ(stdin,*) vmax

WRITE(stdout,'(/,5x,"Minimum volume (unit: (a.u.)^3)?")')
READ(stdin,*) vmin

WRITE(stdout,'(/,5x,"Number of points between maximum and minimum volume ?")')
READ(stdin,*) npt

ALLOCATE(vm(npt))
ALLOCATE(eout(npt))
ALLOCATE(pout(npt))
ALLOCATE(dpout(npt))
ALLOCATE(bulkout(npt))
ALLOCATE(dbulkout(npt))
ALLOCATE(d2bulkout(npt))
ALLOCATE(p_num(npt))
ALLOCATE(dp_num(npt))
ALLOCATE(bulk_num(npt))
ALLOCATE(dbulk_num(npt))
ALLOCATE(d2bulk_num(npt))

WRITE(stdout,'(/,5x,"polynomial? (Yes/No)")')
READ(stdin,*) polyout

lpoly = .FALSE.

IF ((TRIM(polyout)=='Yes') .OR. (TRIM(polyout)=='YES') .OR. (TRIM(polyout)=='yes')) THEN
    lpoly=.TRUE.
    WRITE(stdout,'(/,5x,"The degree of the polynominal? ")')
    READ(stdin,*) degree
    m1 = degree+1
    ALLOCATE(a(m1))
    WRITE(stdout,'(/,5x,"Coefficients of the polynomial: ")')
    DO i=1,degree+1
        READ(stdin,*) a(i)
    ENDDO
ENDIF

WRITE(stdout,*) lpoly
deltav= (vmax-vmin)/npt

IF (lpoly) THEN
   DO j=1,npt
      vm(j)=vmin+(j-1)*deltav
      CALL eos_energy_pol(ieos, vm(j), eout(j), v0, b0/ry_kbar, b01,  &
                                                   b02*ry_kbar, a, m1)
      CALL eos_press_pol(ieos, vm(j), pout(j), v0, b0/ry_kbar, b01,   &
                                                   b02*ry_kbar, a, m1)
      CALL eos_dpress_pol(ieos, vm(j), dpout(j), v0, b0/ry_kbar, b01, &
                                                   b02*ry_kbar, a, m1)
      CALL eos_bulk_pol(ieos, vm(j), bulkout(j), dbulkout(j),         &
                        d2bulkout(j), v0, b0/ry_kbar, b01,            &
                                                   b02*ry_kbar, a, m1)
   ENDDO
ELSE
   DO j=1,npt
      vm(j)=vmin+(j-1)*deltav
      CALL eos_energy(ieos, vm(j), eout(j), v0, b0/ry_kbar, b01, b02*ry_kbar)
      CALL eos_press(ieos, vm(j), pout(j), v0, b0/ry_kbar, b01, b02*ry_kbar)
      CALL eos_dpress(ieos, vm(j), dpout(j), v0, b0/ry_kbar, b01, b02*ry_kbar)
      CALL eos_bulk(ieos, vm(j), bulkout(j), dbulkout(j), d2bulkout(j), &
                                    v0, b0/ry_kbar, b01, b02*ry_kbar)
   ENDDO
ENDIF

OPEN(UNIT=28, FILE='outdata', STATUS='UNKNOWN')

WRITE(28,'("#",2x,"Volume(a.u.)^3",5x,"Energy (Ry)",5x,"Press (kbar)",3x,&
                     "dPress",11x,"B0 (kbar)",9x,"B01",9x,"B02 (1/kbar)")')
DO j=1,npt
    WRITE(28,'(7e16.6)') vm(j), eout(j), pout(j)*ry_kbar, dpout(j)*ry_kbar,&
                         bulkout(j)*ry_kbar, dbulkout(j), d2bulkout(j)/ry_kbar
ENDDO

CLOSE(UNIT=28, STATUS='KEEP')

OPEN(UNIT=28, FILE='numerical_results', STATUS='UNKNOWN')

DO j=2,npt-1
   p_num(j)= -(eout(j+1)-eout(j-1))/(vm(j+1)-vm(j-1))
   dp_num(j)= -(pout(j+1)-pout(j-1))/(vm(j+1)-vm(j-1))
   bulk_num(j)= dp_num(j)*vm(j)
   dbulk_num(j)= (bulkout(j+1)-bulkout(j-1))/(pout(j+1)-pout(j-1))
   d2bulk_num(j)= (dbulkout(j+1)-dbulkout(j-1))/(pout(j+1)-pout(j-1))
ENDDO

IF (ieos==2) THEN
   WRITE(28,'("#",2x,"Volume(a.u.)^3",5x,"dP/P",11x,"d dP/dP",&
                    &10x,"dB0/B0",8x," dB01/B01",8x,"dB02/B02")')
!  DO j=2,npt-1
!     WRITE(28,'(6e16.6)') vm(j), pout(j)*ry_kbar-p_num(j)*ry_kbar, &
!                          dpout(j)*ry_kbar-dp_num(j)*ry_kbar, &
!     bulkout(j)*ry_kbar-bulk_num(j)*ry_kbar, dbulkout(j)-dbulk_num(j), &
!     d2bulkout(j)/ry_kbar-d2bulk_num(j)/ry_kbar
!  ENDDO

   DO j=2,npt-1
      WRITE(28,'(6e16.6)') vm(j), (pout(j)-p_num(j))/pout(j),         &
                           (dpout(j)-dp_num(j))/dpout(j),             &
                           (bulkout(j)-bulk_num(j))/bulkout(j),       &
                           (dbulkout(j)-dbulk_num(j))/dbulkout(j),    &
                           (d2bulkout(j)-d2bulk_num(j))/d2bulk_num(j)
   ENDDO
ELSE
   WRITE(28,'("#",2x,"Volume(a.u.)^3",5x,"dP/P",11x,"d dP/dP",&
                    &10x,"dB0/B0",8x," dB01/B01")')
!  DO j=2,npt-1
!     WRITE(28,'(5e16.6)') vm(j), pout(j)*ry_kbar-p_num(j)*ry_kbar, &
!                          dpout(j)*ry_kbar-dp_num(j)*ry_kbar, &
!     bulkout(j)*ry_kbar-bulk_num(j)*ry_kbar, dbulkout(j)-dbulk_num(j)
!  ENDDO

  DO j=2,npt-1
      WRITE(28,'(5e16.6)') vm(j), (pout(j)-p_num(j))/pout(j),    &
                           (dpout(j)-dp_num(j))/dpout(j),        &
                           (bulkout(j)-bulk_num(j))/bulkout(j),  &
                           (dbulkout(j)-dbulk_num(j))/dbulkout(j)
   ENDDO
ENDIF

CLOSE(UNIT=28, STATUS='KEEP')

DEALLOCATE(vm)
DEALLOCATE(eout)
DEALLOCATE(pout)
DEALLOCATE(dpout)
DEALLOCATE(bulkout)
DEALLOCATE(dbulkout)
DEALLOCATE(d2bulkout) 
DEALLOCATE(p_num)
DEALLOCATE(dp_num)
DEALLOCATE(bulk_num)
DEALLOCATE(dbulk_num)
DEALLOCATE(d2bulk_num)

IF (lpoly) DEALLOCATE(a)

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM test_eos

