!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE h_pcg_step(conv, thresh, dr2)
!
!   This routine receives the vectors evc1, res, pres, dir that contain the
!   present position, residual, preconditioned residual, and search direction 
!   and dir_new that contains the action of the matrix A on the vector dir. 
!   It substitutes these vectors with the corresponding evc1, res, pres,
!   and dir at the following step.
!   If the scalar product of the new residual with its preconditioned form
!   is lower than thresh it sets conv to true, otherwise to false.
!   This version of the routine is only for hermitean matrices.
!   The routine requires two auxiliary routines that must be called
!   lr_dot, to make the scalar product between vectors and lr_prec, to 
!   precondition the residual vector.
!
    USE kinds,        ONLY : DP
    USE lr_global,    ONLY : size_evc1
    USE lr_cg,        ONLY : evc1, res, pres, dir, dir_new
    !
    IMPLICIT NONE
    LOGICAL,  INTENT(OUT) :: conv
    REAL(DP), INTENT(IN)  :: thresh
    REAL(DP), INTENT(OUT) :: dr2
    !
    ! Local variables
    !
    COMPLEX(DP), EXTERNAL :: lr_dot
    COMPLEX(DP) :: alpha, beta, den, num, num1
    !
    CALL start_clock('pcg_step')
    !
    !   new position and new residual
    !
    num = lr_dot(res, pres)
    den = lr_dot(dir, dir_new)
    alpha = num / den
    CALL ZAXPY(size_evc1,alpha,dir,1,evc1,1)
    CALL ZAXPY(size_evc1,-alpha,dir_new,1,res,1)
    !
    !  apply preconditioning
    !
    CALL lr_prec(res, pres)
    !
    !   check for convergence
    !
    num1 = lr_dot(res, pres)
    dr2=ABS(num1)
    conv = (dr2 < thresh)
    IF (conv) RETURN
    !
    !  computes the conjugate search direction
    !
    beta = num1 / num
    !
    dir_new = pres
    CALL ZAXPY(size_evc1,beta,dir,1,dir_new,1)
    dir=dir_new
    !
    CALL stop_clock('pcg_step')
    !
    RETURN
    !
END SUBROUTINE h_pcg_step
