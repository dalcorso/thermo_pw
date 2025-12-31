!
! Copyright (C) 2025 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_tau_acc(celldm_geo, tau_acc, uint_geo, nwork,            &
                                                   nat, nint_var, it)
!-----------------------------------------------------------------------
!
!   This routine sets the grid of values on uint_geo.
!   Data are arranged so that the ninternal displacements for the
!   first set of external parameters is written first. Then there 
!   are the ninternal displacements for the second set of external
!   parameters, etc.
!
USE kinds,         ONLY : DP
USE thermo_mod,    ONLY : ibrav_geo, at_geo
USE control_elastic_constants, ONLY : ninternal_ec, int_ngeo_ec, &
                       nint_var_ec, int_step_ngeo_ec, iconstr_internal_ec
USE control_atomic_pos, ONLY : max_nint_var
USE cell_base,     ONLY : celldm
USE ions_base,     ONLY : tau

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, nat, nint_var, it
REAL(DP), INTENT(IN) :: celldm_geo(6,nwork)
REAL(DP), INTENT(INOUT) :: tau_acc(3,nat,nwork)
REAL(DP), INTENT(INOUT) :: uint_geo(max_nint_var,nwork)

INTEGER  :: igeo(max_nint_var), id(max_nint_var)
INTEGER  :: iwork, total_work, iw, nexternal, iexternal, iinternal, ivar
REAL(DP) :: delta(nint_var)
!
!  then generate a mesh of dimensions nint_var of values of u centered
!  in 0.0
!
delta=0.0_DP
DO iw=1,nint_var
   IF (MOD(int_ngeo_ec(iw,it),2)==0) delta(iw)=int_step_ngeo_ec(iw,it)/2.0_DP
   id(iw)=int_ngeo_ec(iw,it)/2+1 
ENDDO
!
!  compute the displacement with respect to the strained atomic
!  positions
!
uint_geo=0.0_DP
total_work=0
nexternal=nwork/ninternal_ec(it)
DO iexternal=1,nexternal
   DO iinternal=1,ninternal_ec(it)
      CALL find_ipoint(iinternal, nint_var_ec(it), int_ngeo_ec(:,it), igeo)
      total_work=total_work+1
      DO ivar=1,nint_var_ec(it) 
         uint_geo(ivar,total_work)=(igeo(ivar)-id(ivar))*&
                   int_step_ngeo_ec(ivar,it) + delta(ivar) 
      ENDDO
   ENDDO
ENDDO
!
!  and determine for each value of u the atomic coodinates
!
DO iwork=1,nwork
   CALL internal_to_tau(ibrav_geo(iwork), celldm_geo(1,iwork), &
                       tau_acc(1,1,iwork), at_geo(1,1,iwork),  &
                       uint_geo(1,iwork), nat, nint_var_ec(it),  &
                       iconstr_internal_ec(it), 1)
ENDDO

RETURN
END SUBROUTINE set_tau_acc
!
!-----------------------------------------------------------------------
SUBROUTINE internal_to_tau(ibrav_geo, celldm_geo, tau_geo, at_geo, uint_geo, &
                               nat, nint_var, iconstr_internal, iflag)
!-----------------------------------------------------------------------
!
!   This routine transforms the values of uint_geo into atomic 
!   coordinates tau_geo (iflag=1) or from the values of tau_geo
!   it finds u_geo (iflag=2)
!
USE kinds,         ONLY : DP
USE control_atomic_pos, ONLY : ninternal
USE io_global, ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: iflag    ! 1 use uint_geo to obtain tau_geo
                                ! 2 use tau_geo to obtain u
INTEGER, INTENT(IN) :: nat, nint_var
INTEGER, INTENT(IN) :: iconstr_internal, ibrav_geo
REAL(DP), INTENT(IN) :: celldm_geo(6), at_geo(3,3)
REAL(DP), INTENT(INOUT) :: tau_geo(3,nat)
REAL(DP), INTENT(INOUT) :: uint_geo(nint_var)

INTEGER  :: na
REAL(DP) :: delta(nint_var), tau_aux(3,nat), a, c_a, fact

IF ( iconstr_internal==1) THEN
!
!  In this constraint u is one dimensional and it is the parameter that
!  determines the atomic positions of the wurtzite structure. Atomic
!  coordinates on output are cartesian in units of celldm(1)
!
   IF (iflag==1) THEN
      tau_geo(:,:)=0.0_DP
      tau_geo(1,2)=0.5_DP
      tau_geo(2,2)=-SQRT(3.0_DP) / 6.0_DP
      tau_geo(3,2)= celldm_geo(3)/2.0_DP
      tau_geo(3,3)= celldm_geo(3) * uint_geo(1)
      tau_geo(1,4)=0.5_DP
      tau_geo(2,4)=-SQRT(3.0_DP) / 6.0_DP
      tau_geo(3,4)= celldm_geo(3)*(1.0_DP/2.0_DP+uint_geo(1))
   ELSEIF (iflag==2) THEN
!
!     This routine works only if atom 1 and 2 are of the same type
!     and are in the origin and in the middle of the cell (0,0,c/2a).
!     Atom 3 and 4 are assumed to be in (1/2, -\sqrt(3)/6, u c/a) and 
!     (1/2, -\sqrt(3)/6, (u+1/2) c/a). A global translation of all atoms
!     is also allowed
!     
!     Adjust tau if ibrav is zero
      IF (ibrav_geo==0) THEN
         a=at_geo(1,1)*celldm_geo(1)
         c_a= at_geo(3,3) / at_geo(1,1)
         fact = celldm_geo(1) / a
      ELSE
         fact=1.0_DP
         a=celldm_geo(1)
         c_a=celldm_geo(3)
      ENDIF

      tau_aux=tau_geo
      DO na=1,nat
         tau_aux(:,na)=tau_aux(:,na)-tau_geo(:,1)
      ENDDO

      tau_aux=tau_aux * fact
      uint_geo(1)=0.0_DP
      DO na=1,nat
         IF (uint_geo(1)==0.0_DP.AND.ABS(tau_aux(3,na)-c_a*0.5_DP)>1.D-2) THEN
             uint_geo(1)=tau_aux(3,na) / c_a
             IF (uint_geo(1)>0.5_DP) uint_geo(1)=uint_geo(1)-0.5_DP
         ENDIF
      ENDDO


!      WRITE(6,*) 'a= ', a
!      WRITE(6,*) 'c/a= ', c_a
!      WRITE(6,*) 'u= ', uint_geo(1)
!      DO na=1,nat
!         WRITE(stdout,'(3f20.12)') tau_aux(1,na), tau_aux(2,na), tau_aux(3,na)
!      ENDDO
   ENDIF
ELSEIF ( iconstr_internal==2) THEN
!
!  In this constraint the hcp structure is distorted with a strain
! (\epsilon,0,0,0,0,0). In this case the routine receives the displacement
! u along y (one dimensional) and gives as output the displacements of
! the coordinates of the two atoms of the hcp structure
!
    IF (iflag==1) THEN
       tau_geo(:,:) = 0.0_DP
       tau_geo(2,1) = - uint_geo(1) / 2.0_DP
       tau_geo(2,2) =   uint_geo(1) / 2.0_DP
    ELSEIF (iflag==2) THEN
       uint_geo(1) = tau_geo(2,2) * 2.0_DP
    ENDIF
ELSEIF ( iconstr_internal==3) THEN
!
!  This constrain is the monoclinic distortion of a wurtzite structure
!  needed to compute the piezoelectric tensor e_15 or the elastic constant
!  C_44. In this case the internal degree of freedom is the x coordinate
!  of the first atom (the atom in c/2 also moves of x, while the other
!  two atoms moves at -x).
!
!   Adjust tau if ibrav is zero

   IF (iflag==1) THEN
      tau_geo=0.0_DP
      tau_geo(1,3)= uint_geo(1)
      tau_geo(1,4)= uint_geo(1)

   ELSEIF (iflag==2) THEN
!
!  bring the first atom in the origin
!
!      tau_aux=tau_geo/ fact
!      DO na=1,nat
!         tau_aux(:,na)=tau_aux(:,na)-tau_geo(:,1)
!      ENDDO
      uint_geo(1)=tau_geo(1,1)
   ENDIF
ELSEIF ( iconstr_internal==4) THEN
!
!  This constrain is the orthorombic distortion of a wurtzite structure
!  needed to compute the elastic constant C_11. In this case the internal 
!  degree of freedom is the y coordinate of the atom at heigth c/2 and
!  at height (1/2+u)c, and the parameter u of the wurtzite. Here however
!  we have to put the difference between u and u_0 (the equilibrium value
!  of u)
   IF (ibrav_geo==0) THEN
      a=at_geo(1,1)*celldm_geo(1)
      c_a= at_geo(3,3) / at_geo(1,1)
      fact = at_geo(1,1) 
   ELSE
      fact=1.0_DP
      a=celldm_geo(1)
      c_a=celldm_geo(3)
   ENDIF
!  
    IF (iflag==1) THEN
      tau_geo=0.0_DP
      tau_geo(2,2)= uint_geo(1)
      tau_geo(2,4)= uint_geo(1)
      tau_geo(3,3)= uint_geo(2) * c_a
      tau_geo(3,4)= uint_geo(2) * c_a
    ELSEIF (iflag==2) THEN
      uint_geo(1)= tau_geo(2,2) 
      uint_geo(2)= tau_geo(3,3) / c_a
    ENDIF
ELSEIF ( iconstr_internal==5) THEN
!
!  This constrain is the distortion of a wurtzite structure
!  that keeps the hexagonal lattice. In this case only the parameter
!  u changes.
!  Here put the difference between u and u_0 (the equilibrium value
!  of u)
!  

    IF (iflag==1) THEN
      tau_geo=0.0_DP
      tau_geo(3,3)= uint_geo(1) * celldm_geo(3)
      tau_geo(3,4)= uint_geo(1) * celldm_geo(3)
    ELSEIF (iflag==2) THEN
      uint_geo(1)= tau_geo(3,3) / celldm_geo(3)
    ENDIF
ELSEIF ( iconstr_internal==6) THEN
!
!  This constrain is the orthorombic distortion of a wurtzite structure
!  needed to compute the elastic constant C_11 using 3 parameters. 
!  In this case the internal 
!  degree of freedom is the y coordinate of the atom at heigth c/2 and
!  at height (1/2+u)c, and the parameter u of the wurtzite. Here however
!  we have to put the difference between u and u_0 (the equilibrium value
!  of u)
!  

    IF (iflag==1) THEN
      tau_geo=0.0_DP
      tau_geo(2,2)= uint_geo(1)
      tau_geo(2,3)= uint_geo(2)
      tau_geo(2,4)= uint_geo(1) - uint_geo(2)
      tau_geo(3,3)= uint_geo(3) * celldm_geo(3)
      tau_geo(3,4)= uint_geo(3) * celldm_geo(3)
    ELSEIF (iflag==2) THEN
      uint_geo(1)= tau_geo(2,2)
      uint_geo(2)= tau_geo(2,3)
      uint_geo(3)= tau_geo(3,3) / celldm_geo(3)
    ENDIF
ELSEIF ( iconstr_internal==7) THEN
!
!  This constrain is the distortion of a zincblend structure
!  with a rhombohedral straini e_5. The atoms moves along the x direction.
!  Here put the difference between u and u_0 (the equilibrium value
!  of u)
!  
    IF (iflag==1) THEN
      tau_geo=0.0_DP
      tau_geo(1,2)= uint_geo(1) 
    ELSEIF (iflag==2) THEN
      uint_geo(1)= tau_geo(1,2) - 0.25_DP 
    ENDIF
ENDIF

RETURN
END SUBROUTINE internal_to_tau
!
