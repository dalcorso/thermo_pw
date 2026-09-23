!
! Copyright (C) 2025 - present Afrasyiab Ahmed and Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE numerical_derivatives_module
!---------------------------------------------------------------------------
USE kinds, ONLY : DP
IMPLICIT NONE
PRIVATE
SAVE
PUBLIC numerical_derivatives_nd, numerical_derivatives_nd_3,    &
       numerical_first_derivatives_nd, build_poly2_from_taylor, &
       build_poly1_from_taylor, build_stencil_19

CONTAINS

!==============================================================
! 19‑POINT STENCIL BUILDER 
!==============================================================
SUBROUTINE build_stencil_19(nvar, njump, stencil)
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: nvar, njump
   INTEGER, INTENT(OUT) :: stencil(nvar,19)
   INTEGER :: i, j, k, n
   n = 0
   DO k = -1, 1
      DO j = -1, 1
         DO i = -1, 1
            IF (ABS(i) + ABS(j) + ABS(k) <= 2) THEN
               n = n + 1
               stencil(:,n) = njump * (/ i, j, k /)
            END IF
         END DO
      END DO
   END DO

END SUBROUTINE build_stencil_19
!==============================================================
! NEW get_f (19‑point stencil)
!==============================================================
REAL(DP) FUNCTION get_f(f_vals, idx, nvar, stencil)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
INTEGER, INTENT(IN) :: idx(nvar)
INTEGER, INTENT(IN) :: stencil(:, :)
REAL(DP), INTENT(IN) :: f_vals(:)
INTEGER :: i, j
LOGICAL :: match

!DO i = 1, SIZE(f_vals)
DO i = 1, MIN(SIZE(f_vals), SIZE(stencil,2))
   match = .TRUE.
   DO j = 1, nvar
      IF (stencil(j,i) /= idx(j)) THEN
         match = .FALSE.
         EXIT
      ENDIF
   END DO
   IF (match) THEN
      get_f = f_vals(i)
      RETURN
   ENDIF
END DO

CALL errore('get_f','point not found in 19‑point stencil',1)

END FUNCTION get_f


!==============================================================
! NEW numerical_derivatives_nd (1+2*n*n‑point stencil)
!==============================================================
!----------------------------------------------------------------------------
SUBROUTINE numerical_derivatives_nd(f_vals, x, nvar, xi0, f0, njump, dfdx, &
                                    d2fdxdy)
!----------------------------------------------------------------------------
  IMPLICIT NONE

  ! ====== INTENT declarations ======
  INTEGER, INTENT(IN) :: nvar, njump
  REAL(DP), INTENT(IN) :: f_vals(:)
  REAL(DP), INTENT(IN) :: x(:,:)
  REAL(DP), INTENT(OUT) :: xi0(nvar)
  REAL(DP), INTENT(OUT) :: f0
  REAL(DP), INTENT(OUT) :: dfdx(nvar)
  REAL(DP), INTENT(OUT) :: d2fdxdy(nvar*(nvar+1)/2)

  ! ====== Local variables ======
  INTEGER :: ivar, jvar, indi
  INTEGER :: idx(nvar)
  REAL(DP) :: h(nvar)
  REAL(DP) :: f_plus, f_minus, f_pp, f_pm, f_mp, f_mm
  REAL(DP) :: d2fdx2(nvar)
  REAL(DP) :: x_plus, x_minus
  REAL(DP), PARAMETER :: EPS_H = 1.0D-14
  REAL(DP) :: denom
  INTEGER, ALLOCATABLE :: stencil(:,:)

  ! ====== Build stencil internally ======
  ALLOCATE(stencil(nvar, 1 + 2*nvar*nvar))
  CALL build_stencil_19(nvar, njump, stencil)

  ! ====== Compute derivatives ======
  idx(:) = 0
  f0 = get_f(f_vals, idx, nvar, stencil)
  !xi0(:) = 0.0_DP
! new lines
  DO ivar = 1, nvar
     xi0(ivar) = get_f(x(ivar,:), idx, nvar, stencil)
  END DO
!new above lines
  DO ivar = 1, nvar

     idx(:) = 0
     idx(ivar) = njump
     f_plus = get_f(f_vals, idx, nvar, stencil)
     x_plus = get_f(x(ivar,:), idx, nvar, stencil)

    idx(:) = 0
    idx(ivar) = -njump
    f_minus = get_f(f_vals, idx, nvar, stencil)
    x_minus = get_f(x(ivar,:), idx, nvar, stencil)

    ! positive step size
    h(ivar) = ABS(x_plus - x_minus) * 0.5_DP
     ! guard against zero/too-small step
    IF (h(ivar) <= EPS_H) &
       CALL errore('numerical_derivatives_nd', &
               'zero or too small step size h', ivar)
    dfdx(ivar) = (f_plus - f_minus) / (2.0_DP * h(ivar))
    d2fdx2(ivar) = (f_plus - 2.0_DP*f0 + f_minus) / (h(ivar)**2)

  END DO

  indi = 1
  DO ivar = 1, nvar
     DO jvar = ivar, nvar

        IF (ivar == jvar) THEN
           d2fdxdy(indi) = d2fdx2(ivar)
        ELSE
           idx(:) = 0
           idx(ivar) = njump
           idx(jvar) = njump
           f_pp = get_f(f_vals, idx, nvar, stencil)

           idx(:) = 0
           idx(ivar) = njump
           idx(jvar) = -njump
           f_pm = get_f(f_vals, idx, nvar, stencil)

           idx(:) = 0
           idx(ivar) = -njump
           idx(jvar) = njump
           f_mp = get_f(f_vals, idx, nvar, stencil)

           idx(:) = 0
           idx(ivar) = -njump
           idx(jvar) = -njump
           f_mm = get_f(f_vals, idx, nvar, stencil)

           denom = 4.0_DP * h(ivar) * h(jvar)
           IF (ABS(denom) <= EPS_H) &
              CALL errore('numerical_derivatives_nd',                   &
                          'zero denominator in mixed derivative', indi)

           d2fdxdy(indi) = (f_pp - f_pm - f_mp + f_mm) / denom

        END IF

        indi = indi + 1
     END DO
  END DO

  DEALLOCATE(stencil)

END SUBROUTINE numerical_derivatives_nd
!----------------------------------------------------------------------
SUBROUTINE numerical_derivatives_nd_3(f_vals, x, nvar, xi0, f0, &
                                               dfdx, d2fdxdy, ind_gen3)
!----------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, INTENT(IN)   :: nvar
REAL(DP), INTENT(IN)  :: f_vals((nvar*(nvar+3))/2+1)
REAL(DP), INTENT(IN)  :: x(nvar, (nvar*(nvar+3))/2+1)
INTEGER, INTENT(IN)   :: ind_gen3(nvar, (nvar*(nvar+3))/2+1)
REAL(DP), INTENT(OUT) :: xi0(nvar)
REAL(DP), INTENT(OUT) :: f0
REAL(DP), INTENT(OUT) :: dfdx(nvar)
REAL(DP), INTENT(OUT) :: d2fdxdy(nvar*(nvar+1)/2)

INTEGER :: ivar, jvar, indi
INTEGER :: idx(nvar)
REAL(DP) :: f_plus, f_minus, f_pp, f_pm, f_mp, f_mm, x_plus, x_minus
REAL(DP) :: h(nvar), d2fdx2(nvar), xdata((nvar*(nvar+3))/2+1)


idx(:) = 0
f0 = get_f3(f_vals, idx, ind_gen3, nvar)
!
! First derivatives
!
do ivar = 1, nvar
    xdata(1:(nvar*(nvar+3))/2+1)=x(ivar,1:(nvar*(nvar+3))/2+1)
    idx(:) = 0
    x_plus = get_f3(xdata, idx, ind_gen3, nvar)
    xi0(ivar) = x_plus
    idx(ivar) = 1
    f_plus = get_f3(f_vals, idx, ind_gen3, nvar)   ! f(x_i + h_i)
    x_plus = get_f3(xdata, idx, ind_gen3, nvar)
    idx(:) = 0
    idx(ivar) = -1
    f_minus = get_f3(f_vals, idx, ind_gen3,  nvar)  ! f(x_i - h_i)
    x_minus = get_f3(xdata, idx, ind_gen3, nvar)
    h(ivar) = (x_plus - x_minus) * 0.5_DP
    dfdx(ivar) = (f_plus - f_minus) / (2.0_DP * h(ivar))
end do

indi=1
DO ivar = 1, nvar
!
! Second derivatives - diagonal terms
!
    idx(:) = 0
    idx(ivar) = 1
    f_plus = get_f3(f_vals, idx, ind_gen3, nvar)   ! f(x_i + h_i)

    idx(:) = 0
    idx(ivar) = -1
    f_minus = get_f3(f_vals, idx, ind_gen3, nvar)  ! f(x_i - h_i)

    d2fdxdy(indi) = (f_plus - 2.0_DP*f0 + f_minus) / (h(ivar)**2)
    indi=indi+1
!
! Mixed second derivatives (off-diagonal)
!
    DO jvar = ivar+1, nvar
        idx(:) = 0
        idx(ivar) = 1
        idx(jvar) = 1
        f_pp = get_f3(f_vals, idx, ind_gen3, nvar)  ! f(x_i + h_i, x_j + h_j)

        ! f_pm = f(x_i + h_i, x_j)
        idx(:) = 0
        idx(ivar) = 1
        f_pm = get_f3(f_vals, idx, ind_gen3, nvar)

        ! f_mp = f(x_i, x_j + h_j)
        idx(:) = 0
        idx(jvar) = 1
        f_mp = get_f3(f_vals, idx, ind_gen3, nvar)

        ! f_mm = f(x) = f0
        d2fdxdy(indi) = (f_pp - f_pm - f_mp + f0) / (h(ivar)*h(jvar))

        indi=indi+1
    ENDDO
 ENDDO

RETURN
END SUBROUTINE numerical_derivatives_nd_3
!
!----------------------------------------------------------------------
SUBROUTINE numerical_first_derivatives_nd(f_vals, x, nvar, xi0, f0, &
                                                     dfdx )
!----------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nvar
REAL(DP), INTENT(IN) :: f_vals(2*nvar+1)
REAL(DP), INTENT(IN) :: x(nvar, 2*nvar+1)
REAL(DP), INTENT(OUT) :: xi0(nvar)
REAL(DP), INTENT(OUT) :: f0
REAL(DP), INTENT(OUT) :: dfdx(nvar)

INTEGER :: ivar
INTEGER :: offset
REAL(DP) :: f_plus, f_minus, x_plus, x_minus
REAL(DP) :: h(nvar), xdata(2*nvar+1)

f0=get_f_1d(f_vals, 0, nvar)
!
! First partial derivatives
!
DO ivar = 1, nvar
    xdata(1:2*nvar+1)=x(ivar,1:2*nvar+1)
    offset=0
    x_plus = get_f_1d(xdata, offset, nvar)
    xi0(ivar) = x_plus
    offset = ivar
    f_plus = get_f_1d(f_vals, offset, nvar)
    f_minus = get_f_1d(f_vals,- offset, nvar)
    x_plus = get_f_1d(xdata, offset, nvar)
    x_minus = get_f_1d(xdata, - offset, nvar)
    h(ivar) = (x_plus - x_minus) * 0.5_DP
    dfdx(ivar) = (f_plus - f_minus) / (2.0_DP * h(ivar))
END DO
RETURN
END SUBROUTINE numerical_first_derivatives_nd
!
!----------------------------------------------------------------------
 SUBROUTINE build_poly2_from_taylor(poly, xi0, f0, dfdx, d2fdxdy, nvar)
!----------------------------------------------------------------------
!
! Build poly2 from Taylor coefficients
!
  USE polynomial, ONLY : poly2
  USE quadratic_surfaces, ONLY : print_quadratic_polynomial
  IMPLICIT NONE
  TYPE(poly2), INTENT(INOUT) :: poly
  INTEGER, INTENT(IN)   :: nvar
  REAL(DP), INTENT(IN)  :: f0, dfdx(nvar), d2fdxdy((nvar*(nvar + 1)) / 2), &
                              xi0(nvar)
  INTEGER :: ivar, jvar, idx

  poly%nvar    = nvar
  poly%ncoeff2 = (nvar * (nvar + 1)) / 2
!
!  zero degree terms
!
  poly%a0      = f0
!
!  first degree terms
!
  DO ivar = 1, nvar
     poly%a0=poly%a0-dfdx(ivar)*xi0(ivar)
     poly%phi1(ivar)= dfdx(ivar)
  ENDDO
!
!  second degree terms
!
  idx = 1
  DO ivar = 1, nvar
     DO jvar = ivar, nvar
        IF (ivar == jvar) THEN
           poly%a0=poly%a0 + 0.5_DP * d2fdxdy(idx) * xi0(ivar) * xi0(ivar)
           poly%phi1(ivar) = poly%phi1(ivar) - d2fdxdy(idx) * xi0(jvar)
           poly%phi2(idx) = 0.5_DP * d2fdxdy(idx)  
        ELSE    
           poly%a0=poly%a0 + d2fdxdy(idx) * xi0(ivar) * xi0(jvar)
           poly%phi1(ivar) = poly%phi1(ivar) - d2fdxdy(idx) * xi0(jvar)
           poly%phi1(jvar) = poly%phi1(jvar) - d2fdxdy(idx) * xi0(ivar)
           poly%phi2(idx) = d2fdxdy(idx)  
        ENDIF
        idx = idx + 1
      ENDDO
   ENDDO
   CALL print_quadratic_polynomial(nvar, poly)
   RETURN
END SUBROUTINE build_poly2_from_taylor
!
!----------------------------------------------------------------------
 SUBROUTINE build_poly1_from_taylor(poly, xi0, f0, dfdx, nvar)
!----------------------------------------------------------------------
!
! Build poly1 from Taylor coefficients
!
  USE polynomial, ONLY : poly1
  IMPLICIT NONE
  TYPE(poly1), INTENT(INOUT) :: poly
  INTEGER, INTENT(IN)   :: nvar
  REAL(DP), INTENT(IN)  :: f0, dfdx(nvar), xi0(nvar)
  INTEGER :: ivar

  poly%nvar    = nvar
!
!  zero degree terms
!
  poly%a0      = f0
!
!  first degree terms
!
  DO ivar = 1, nvar
     poly%a0=poly%a0-dfdx(ivar)*xi0(ivar)
     poly%phi1(ivar)= dfdx(ivar)
  ENDDO
  RETURN
END SUBROUTINE build_poly1_from_taylor
!
!----------------------------------------------------------------------
REAL(DP) FUNCTION get_f3(f_vals, idx, ind_gen3, nvar) 
!----------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: idx(nvar), nvar
INTEGER, INTENT(IN) :: ind_gen3(nvar,(nvar*(nvar+3))/2+1)
REAL(DP), INTENT(IN) :: f_vals((nvar*(nvar+3))/2+1)

INTEGER :: flat_idx, iaux, ivar, jvar

flat_idx=0
DO ivar=1, (nvar*(nvar+3))/2+1
   iaux=0
   DO jvar=1,nvar
      IF (iaux/=0) CYCLE
      iaux= iaux+ ABS(idx(jvar)-ind_gen3(jvar,ivar))
   ENDDO
   IF (iaux==0) flat_idx=ivar
   IF (flat_idx /= 0) EXIT
ENDDO

get_f3 = f_vals(flat_idx)
RETURN
END FUNCTION get_f3

REAL(DP) FUNCTION get_f_1d(f_vals, offset, nvar) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: offset, nvar
  REAL(DP), INTENT(IN) :: f_vals(2*nvar + 1)
  INTEGER :: pos

  IF (offset < -nvar .OR. offset > nvar) CALL errore('get_f_1d',&
                                                     'wrong offset',1)

  !
  ! Index mapping: offset = -nvar -> 1, 0 -> nvar+1, +nvar -> 2 nvar+1
  !
  pos = offset + nvar + 1
  get_f_1d = f_vals(pos)

END FUNCTION get_f_1d


END MODULE numerical_derivatives_module




