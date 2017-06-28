!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE isoentropic
!
!   this module contains the support routines for the calculation
!   of the isoentropic quantities given the corresponding isothermal
!   quantities and the parameters that link the two.
!   It can calculate:
!   The isoentropic elastic constants given the isothermal elastic constants,
!   the isoentropic bulk modulus given the isothermal bulk modulus,
!   the isobaric specific heat given the isochoric one,
!   the average gruneisen parameter.
!   The routine does not receive the isothermal quantities and
!   provide only the difference between isoentropic and isothermal
!   quantities. 
!
!
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast 
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC  isobaric_heat_capacity, isoentropic_bulk_modulus,  &
          isostress_heat_capacity,                           &
          average_gruneisen, thermal_stress,                 &
          isoentropic_elastic_constants, gen_average_gruneisen

CONTAINS
!
SUBROUTINE isobaric_heat_capacity(volume,b0_t,beta_t,temp,cpmcv,ntemp)
!
!  This routine receives the volume, the bulk modulus, the volume thermal
!  expansion at ntemp temperatures and provides the difference between
!  isobaric and isochoric heat capacity.
!
USE kinds, ONLY : DP
USE constants, ONLY : ry_kbar
IMPLICIT NONE

INTEGER :: ntemp
REAL(DP), INTENT(IN) :: volume(ntemp), b0_t(ntemp), beta_t(ntemp), &
                        temp(ntemp)
REAL(DP), INTENT(INOUT) :: cpmcv(ntemp)

INTEGER :: itemp

cpmcv=0.0_DP
DO itemp=2, ntemp-1
   cpmcv(itemp) = temp(itemp) * volume(itemp) * &
                           beta_t(itemp)**2 * b0_t(itemp) / ry_kbar
END DO

RETURN
END SUBROUTINE isobaric_heat_capacity

SUBROUTINE isoentropic_bulk_modulus(volume,b0_t,beta_t,cp_t,temp,bsmbt,ntemp)
!
!  This routine receives the volume, the bulk modulus, the volume thermal
!  expansion and the isobaric heat capacity at ntemp temperatures and 
!  provides the difference between isoentropic and isothermal bulk modulus
!

USE kinds, ONLY : DP
USE constants, ONLY : ry_kbar
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: volume(ntemp), b0_t(ntemp), beta_t(ntemp), &
                        cp_t(ntemp), temp(ntemp)
REAL(DP), INTENT(INOUT) :: bsmbt(ntemp)

INTEGER :: itemp

bsmbt=0.0_DP
DO itemp=2,ntemp-1
   bsmbt(itemp) =  temp(itemp) * volume(itemp) * beta_t(itemp)**2 * &
                   b0_t(itemp)**2 / cp_t(itemp) / ry_kbar
END DO

RETURN
END SUBROUTINE isoentropic_bulk_modulus

SUBROUTINE average_gruneisen(volume,b0_t,beta_t,cv_t,temp,gamma_t,ntemp)
!
!  This routine receives the volume, the bulk modulus, the volume thermal
!  expansion and the isochoric heat capacity at ntemp temperatures and 
!  provides the average Gruneisen parameter
!
USE kinds, ONLY : DP
USE constants, ONLY : ry_kbar
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: volume(ntemp), b0_t(ntemp), beta_t(ntemp), &
                        cv_t(ntemp), temp(ntemp)
REAL(DP), INTENT(OUT) :: gamma_t(ntemp)

INTEGER :: itemp

gamma_t=0.0_DP
DO itemp=2,ntemp-1
   gamma_t(itemp)=beta_t(itemp) * b0_t(itemp) * volume(itemp) / cv_t(itemp) &
                                              / ry_kbar
END DO

RETURN
END SUBROUTINE average_gruneisen

SUBROUTINE thermal_stress(el_con_t,alpha_t,bths,ntemp)
!
!  This routine receives the isothermal elastic constants and
!  the thermal expansion tensor and gives as output the thermal stress as a
!  function of temperature. 
!  Elastic constants and thermal expansion in input are in Voigt 
!  notation, but thermal stresses are given in standard notation as a
!  3x3 matrix.
!  On output thermal stresses are in the same units of the elastic
!  constants
!

USE kinds,      ONLY : DP
USE constants,  ONLY : ry_kbar
USE voigt,      ONLY : to_voigt4
USE strain_mod, ONLY : trans_epsilon
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: el_con_t(6,6,ntemp), alpha_t(6,ntemp) 
REAL(DP), INTENT(OUT) :: bths(3,3,ntemp)

REAL(DP) :: elcon(3,3,3,3), alp(3,3), aux(6,6)
INTEGER :: itemp, i,j,m,n

bths=0.0_DP
DO itemp=2,ntemp-1
   aux(:,:)=el_con_t(:,:,itemp)
   CALL to_voigt4(aux, elcon, .FALSE.)
   CALL trans_epsilon(alpha_t(1,itemp), alp, 1)
   DO i=1,3 
      DO j=1,3
         DO m=1,3
            DO n=1,3
               bths(i,j,itemp) = bths(i,j,itemp) - elcon(i,j,m,n) * alp(m,n)
            END DO
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE thermal_stress

SUBROUTINE isostress_heat_capacity(volume,el_con_t,alpha_t,temp,cpmcv,ntemp)
!
!  This routine receives the isothermal elastic constants and
!  the thermal expansion tensor and gives as output the difference
!  between the constant stress and the constant strain heat capacity.
!

USE kinds,      ONLY : DP
USE constants,  ONLY : ry_kbar
USE voigt,      ONLY : to_voigt4
USE strain_mod, ONLY : trans_epsilon
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: el_con_t(6,6,ntemp), alpha_t(6,ntemp), temp(ntemp), &
                        volume(ntemp)
REAL(DP), INTENT(OUT) :: cpmcv(ntemp)

REAL(DP) :: elcon(3,3,3,3), alp(3,3), aux(6,6)
INTEGER :: itemp, i, j, m, n

cpmcv=0.0_DP
DO itemp=2,ntemp-1
   aux(:,:)=el_con_t(:,:,itemp)
   CALL to_voigt4(aux, elcon, .FALSE.)
   CALL trans_epsilon(alpha_t(1,itemp), alp, 1)
   DO i=1,3 
      DO j=1,3
         DO m=1,3
            DO n=1,3
               cpmcv(itemp) = cpmcv(itemp) + alp(i,j)*elcon(i,j,m,n)*alp(m,n)
            END DO
         END DO
      END DO
   END DO
   cpmcv(itemp)=cpmcv(itemp) * temp(itemp) * volume(itemp) / ry_kbar
END DO

RETURN
END SUBROUTINE isostress_heat_capacity

SUBROUTINE isoentropic_elastic_constants(volume,bths,cv_t,temp,csmct,ntemp)
!
!  This routine receives the thermal stress, in kbar, the specific heat
!  in Ry/cell and gives as output the correction to the elastic constants
!  in kbar
!
USE kinds,     ONLY : DP
USE constants, ONLY : ry_kbar
USE voigt,     ONLY : to_voigt4
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: volume(ntemp), bths(3,3,ntemp), cv_t(ntemp), &
                        temp(ntemp)
REAL(DP), INTENT(INOUT) :: csmct(6,6,ntemp)

REAL(DP) :: aux(3,3,3,3)
INTEGER :: itemp, i, j, m, n

csmct=0.0_DP
DO itemp=2,ntemp-1
   DO i=1,3
      DO j=1,3
         DO m=1,3
            DO n=1,3
               aux(i,j,m,n)= - temp(itemp) * volume(itemp) * bths(i,j,itemp)* &
                               bths(m,n,itemp) / cv_t(itemp) / ry_kbar
            END DO
         END DO
      END DO
   END DO 
!
!  transform to voigt indeces and copy aux in csmct
!
   CALL to_voigt4(csmct(1,1,itemp), aux, .TRUE.) 
END DO

RETURN
END SUBROUTINE isoentropic_elastic_constants

SUBROUTINE gen_average_gruneisen(volume,bths,cv_t,temp,ggamma_t,ntemp)
!
!  This routine receives the thermal stress, the constant strain specific
!  heat and gives the generalized average Gruneisen parameters.
!  Note that we use here the definition reported in 
!  D.C. Wallace, Thermodynamics of Crystals.
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp
REAL(DP), INTENT(IN) :: volume(ntemp), bths(3,3,ntemp), cv_t(ntemp), &
                        temp(ntemp)
REAL(DP), INTENT(INOUT) :: ggamma_t(3,3,ntemp)

INTEGER :: itemp, i, j

ggamma_t=0.0_DP
DO itemp = 2,ntemp-1
   DO i = 1,3
      DO j = 1,3
         ggamma_t(i,j,itemp) = volume(itemp) * bths(i,j,itemp) / cv_t(itemp)
      END DO
   END DO 
END DO

RETURN
END SUBROUTINE gen_average_gruneisen

END MODULE isoentropic
