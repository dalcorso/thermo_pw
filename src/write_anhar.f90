!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_anharmonic()
USE kinds,          ONLY : DP
USE thermodynamics, ONLY : ntemp, temp
USE thermo_mod,     ONLY : vmin_t, b0_t
USE anharmonic,     ONLY : alpha_t, beta_t, cp_t, cv_t, b0_s
USE control_thermo, ONLY : flanhar
USE io_global,      ONLY : ionode
USE mp_images,      ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER :: itemp, iu_therm

IF (my_image_id /= root_image) RETURN

IF ( .NOT. ALLOCATED(alpha_t) ) ALLOCATE( alpha_t(ntemp) )
IF ( .NOT. ALLOCATED(beta_t) )  ALLOCATE( beta_t(ntemp) )

DO itemp = 1, ntemp
   IF (itemp==1) THEN
      beta_t(1)=( vmin_t(itemp+1) - vmin_t(itemp) ) / &
             ( temp(itemp+1,1) - temp(itemp,1) ) / vmin_t(itemp)
   ELSEIF (itemp==ntemp) THEN
      beta_t(itemp) = (vmin_t(itemp)-vmin_t(itemp-1))/                &
                       (temp(itemp,1)-temp(itemp-1,1)) / vmin_t(itemp)
   ELSE
      beta_t(itemp) = (vmin_t(itemp+1)-vmin_t(itemp-1)) / &
                       (temp(itemp+1,1)-temp(itemp-1,1)) / vmin_t(itemp)
   ENDIF
END DO
alpha_t = beta_t / 3.0_DP

CALL compute_cp()

IF (ionode) THEN
   iu_therm=2
   OPEN(UNIT=iu_therm, FILE=TRIM(flanhar), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_therm,'("# alpha is the linear thermal expansion ")')
   WRITE(iu_therm,'("#   T (K)     V(T) (a.u.)^3   alpha (10^(-6)) &
               &  B (T) (kbar)   C_p(T) (Ry/cell)  (C_p - C_v)(T) &
               &  (B_S - B_T) (T) (kbar) " )' )
   DO itemp = 1, ntemp
      WRITE(iu_therm, '(7e15.7)') temp(itemp,1), vmin_t(itemp),   &
                   alpha_t(itemp)*1.D6, b0_t(itemp), cp_t(itemp), &
                                     cp_t(itemp) - cv_t(itemp),   &
                                     b0_s(itemp) - b0_t(itemp)
   END DO
   CLOSE(iu_therm)
END IF

RETURN
END SUBROUTINE write_anharmonic
