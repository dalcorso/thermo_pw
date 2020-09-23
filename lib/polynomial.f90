!
! Copyright (C) 2019 A. Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE polynomial
!
!  This module contains the low level definitions of polynomials.
!  It is used by the linear_surfaces, quadratic_surfaces, cubic_surfaces,
!  quartic_surfaces modules.
!  It can deal with polynomial of any number of variables but it is
!  currently limited to polynomial of up to 4 degree.
!
!  It defines the following types:
!  poly1   A polynomial of degree 1
!  poly2   A polynomial of degree 2
!  poly3   A polynomial of degree 3
!  poly4   A polynomial of degree 4
!
!  init_poly Is an interface that receives a polynomial of degree 1, 2, 3, or
!            4 and calls the routine needed to initialize the polynomial.
!
!  clean_poly Is an interface that receives a polynomial of degree 1, 2, 3, or
!            4 and deallocates its variables.
!  
!  poly_ncoeff : receives the number of variables of the polynomial and its
!                degree d and gives the number of coefficients of a polynomial
!                that has all terms up to those of degree d
!
!  poly_ncoeff_sep : receives the number of variables of the polynomial and
!                 the degree d and gives the number of coeffiecients of a
!                 polynomial that has only the terms of degree d.
!
USE kinds, ONLY : DP
IMPLICIT NONE
PRIVATE
SAVE

TYPE poly1
   INTEGER  :: nvar                  ! the number of variables
   REAL(DP) :: a0                    ! the zero degree coefficient
   REAL(DP), ALLOCATABLE :: phi1(:)  ! coefficients of the linear polynomial
END TYPE poly1

TYPE poly2
   INTEGER  :: nvar                  ! number of variables
   INTEGER  :: ncoeff2               ! number of coefficients of the quadratic 
                                     ! polynomial
   REAL(DP) :: a0                    ! the zero degree coefficient
   REAL(DP), ALLOCATABLE :: phi1(:)  ! coefficients of the linear polynomial
   REAL(DP), ALLOCATABLE :: phi2(:)  ! coefficients of the quadratic polynomial
END TYPE poly2

TYPE poly3
   INTEGER  :: nvar                  ! number of variables
   INTEGER  :: ncoeff2               ! number of coefficients of the quadratic
                                     ! polynomial
   INTEGER  :: ncoeff3               ! number of coefficients of the cubic
                                     ! polynomial 
   REAL(DP) :: a0                    ! the zero degree coefficient
   REAL(DP), ALLOCATABLE :: phi1(:)  ! coefficients of the linear polynomial
   REAL(DP), ALLOCATABLE :: phi2(:)  ! coefficients of the quadratic polynomial
   REAL(DP), ALLOCATABLE :: phi3(:)  ! coefficients of the cubic polynomial
END TYPE poly3

TYPE poly4
   INTEGER  :: nvar                  ! number of variables
   INTEGER  :: ncoeff2               ! number of coefficients of the quadratic
                                     ! polynomial
   INTEGER  :: ncoeff3               ! number of coefficients of the cubic
                                     ! polynomial
   INTEGER  :: ncoeff4               ! number of coefficients of the quartic
                                     ! polynomial
   REAL(DP) :: a0                    ! the zero degree coefficient
   REAL(DP), ALLOCATABLE :: phi1(:)  ! coefficients of the linear polynomial
   REAL(DP), ALLOCATABLE :: phi2(:)  ! coefficients of the quadratic polynomial
   REAL(DP), ALLOCATABLE :: phi3(:)  ! coefficients of the cubic polynomial
   REAL(DP), ALLOCATABLE :: phi4(:)  ! coefficients of the quartic polynomial
END TYPE poly4

PUBLIC :: poly1, poly2, poly3, poly4, init_poly, clean_poly, poly_ncoeff, &
          poly_ncoeff_sep

INTERFACE init_poly

   MODULE PROCEDURE init_poly1, &
                    init_poly2, &
                    init_poly3, &
                    init_poly4  
END INTERFACE init_poly

INTERFACE clean_poly

   MODULE PROCEDURE clean_poly1, &
                    clean_poly2, &
                    clean_poly3, &
                    clean_poly4
END INTERFACE clean_poly

CONTAINS
!
!   The following four routines allocate the variables of a polynomial of
!   degree 1, 2, 3, or 4.
!
!--------------------------------------------------------------------
   SUBROUTINE init_poly1(nvar,poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nvar
   TYPE(poly1), INTENT(INOUT) :: poly_in

   poly_in%nvar=nvar
   ALLOCATE(poly_in%phi1(nvar))

   RETURN
   END SUBROUTINE init_poly1
!
! --------------------------------------------------------------------------
!
   SUBROUTINE init_poly2(nvar,poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nvar
   TYPE(poly2), INTENT(INOUT) :: poly_in

   poly_in%nvar=nvar
   CALL poly_ncoeff_sep(nvar,2,poly_in%ncoeff2)
   ALLOCATE(poly_in%phi1(nvar))
   ALLOCATE(poly_in%phi2(poly_in%ncoeff2))

   RETURN
   END SUBROUTINE init_poly2
!
! --------------------------------------------------------------------------
!
   SUBROUTINE init_poly3(nvar,poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nvar
   TYPE(poly3), INTENT(INOUT) :: poly_in

   poly_in%nvar=nvar
   CALL poly_ncoeff_sep(nvar,2,poly_in%ncoeff2)
   CALL poly_ncoeff_sep(nvar,3,poly_in%ncoeff3)
   ALLOCATE(poly_in%phi1(poly_in%nvar))
   ALLOCATE(poly_in%phi2(poly_in%ncoeff2))
   ALLOCATE(poly_in%phi3(poly_in%ncoeff3))

   RETURN
   END SUBROUTINE init_poly3
!
! --------------------------------------------------------------------------
!
   SUBROUTINE init_poly4(nvar,poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nvar
   TYPE(poly4), INTENT(INOUT) :: poly_in

   poly_in%nvar=nvar
   CALL poly_ncoeff_sep(nvar,2,poly_in%ncoeff2)
   CALL poly_ncoeff_sep(nvar,3,poly_in%ncoeff3)
   CALL poly_ncoeff_sep(nvar,4,poly_in%ncoeff4)
   ALLOCATE(poly_in%phi1(nvar))
   ALLOCATE(poly_in%phi2(poly_in%ncoeff2))
   ALLOCATE(poly_in%phi3(poly_in%ncoeff3))
   ALLOCATE(poly_in%phi4(poly_in%ncoeff4))

   RETURN
   END SUBROUTINE init_poly4
!
! --------------------------------------------------------------------------
!
!   The following four routines deallocate the variables of a polynomial of
!   degree 1, 2, 3, or 4.
!
   SUBROUTINE clean_poly1(poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(poly1), INTENT(INOUT) :: poly_in

   IF (ALLOCATED(poly_in%phi1)) DEALLOCATE(poly_in%phi1)

   RETURN
   END SUBROUTINE clean_poly1
!
! --------------------------------------------------------------------------
!
   SUBROUTINE clean_poly2(poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(poly2), INTENT(INOUT) :: poly_in

   IF (ALLOCATED(poly_in%phi1)) DEALLOCATE(poly_in%phi1)
   IF (ALLOCATED(poly_in%phi2)) DEALLOCATE(poly_in%phi2)

   RETURN
   END SUBROUTINE clean_poly2
!
! --------------------------------------------------------------------------
!
   SUBROUTINE clean_poly3(poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(poly3), INTENT(INOUT) :: poly_in

   IF (ALLOCATED(poly_in%phi1)) DEALLOCATE(poly_in%phi1)
   IF (ALLOCATED(poly_in%phi2)) DEALLOCATE(poly_in%phi2)
   IF (ALLOCATED(poly_in%phi3)) DEALLOCATE(poly_in%phi3)

   RETURN
   END SUBROUTINE clean_poly3
!
! --------------------------------------------------------------------------
!
   SUBROUTINE clean_poly4(poly_in)
!--------------------------------------------------------------------
   IMPLICIT NONE
   TYPE(poly4), INTENT(INOUT) :: poly_in

   IF (ALLOCATED(poly_in%phi1)) DEALLOCATE(poly_in%phi1)
   IF (ALLOCATED(poly_in%phi2)) DEALLOCATE(poly_in%phi2)
   IF (ALLOCATED(poly_in%phi3)) DEALLOCATE(poly_in%phi3)
   IF (ALLOCATED(poly_in%phi4)) DEALLOCATE(poly_in%phi4)

   RETURN
   END SUBROUTINE clean_poly4
!
! --------------------------------------------------------------------------
!
   SUBROUTINE poly_ncoeff_sep(n,d,ncoeff)
!--------------------------------------------------------------------
!
!  This routine gives the number of coefficients of a polynomial of n
!  variables of degree d (only the terms of degree d are accounted for).
!
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: n, d
   INTEGER, INTENT(OUT) :: ncoeff
   INTEGER, ALLOCATABLE :: nmat(:,:)
   INTEGER :: ni, di, dj

   ALLOCATE(nmat(n,0:d))

   nmat(:,:)=0
   nmat(1,:)=1
   DO ni=2, n
      DO di=0, d
         DO dj=0, di
            nmat(ni,di)=nmat(ni,di) + nmat(ni-1,dj)
         ENDDO
      ENDDO
   ENDDO

   ncoeff=nmat(n,d)
   DEALLOCATE(nmat)

   RETURN
   END SUBROUTINE poly_ncoeff_sep
!
! --------------------------------------------------------------------------
!
   SUBROUTINE poly_ncoeff(n,d,ncoeff)
!--------------------------------------------------------------------
!
!  This routine gives the total number of coefficients of a polynomial of
!  n variables of degree d (all terms of degree up to d are accounted for).
!
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: n, d
   INTEGER, INTENT(OUT) :: ncoeff
   INTEGER :: ni, di, nsum

   nsum=0
   DO di=0,d
      CALL poly_ncoeff_sep(n,di,ni)
      nsum=nsum+ni
   ENDDO

   ncoeff=nsum

   RETURN
   END SUBROUTINE poly_ncoeff

END MODULE polynomial
