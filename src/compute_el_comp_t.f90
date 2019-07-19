! Copyright (C) 2019 C. Malica
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_el_comp_t(el_cons_t, el_comp_t, b0_t)
!
!   This routine computes the elastic compliances and the bulk modulus
!   at all temperatures given the elastic constants at all temperatures
!
 USE kinds,              ONLY : DP
 USE cell_base,          ONLY : ibrav
 USE temperature,        ONLY : ntemp
 USE elastic_constants,  ONLY : compute_elastic_compliances, &
                                print_macro_elasticity
 USE mp_world,           ONLY : world_comm
 USE mp,                 ONLY : mp_sum

 IMPLICIT NONE

 REAL(DP) :: el_cons_t(6,6,ntemp), el_comp_t(6,6,ntemp), b0_t(ntemp)
 REAL(DP) :: macro_el(8)
 INTEGER  :: itemp, startt, lastt

 CALL divide(world_comm, ntemp, startt, lastt)
 el_comp_t=0.0_DP
 b0_t=0.0_DP
 DO itemp=startt,lastt
    IF (itemp==1.OR.itemp==ntemp) CYCLE
    CALL compute_elastic_compliances(el_cons_t(:,:,itemp),el_comp_t(:,:,itemp))
    CALL print_macro_elasticity(ibrav,el_cons_t(:,:,itemp), &
                          el_comp_t(:,:,itemp),macro_el,.FALSE.)
    b0_t(itemp)=macro_el(5)
 ENDDO

 CALL mp_sum(el_comp_t, world_comm)
 CALL mp_sum(b0_t, world_comm)

 RETURN
END SUBROUTINE compute_el_comp_t
