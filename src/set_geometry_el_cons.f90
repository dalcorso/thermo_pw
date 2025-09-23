!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE set_geometry_el_cons(iwork)
!----------------------------------------------------------------------
!
!  This routine receives as input the current unperturbed geometry
!  and set the variables of pw with this unperturbed geometry.
!
USE kinds, ONLY : DP
USE cell_base, ONLY : ibrav, celldm
USE ions_base, ONLY : tau, nat, atm, ityp
USE control_elastic_constants, ONLY : el_con_ibrav_geo, el_con_celldm_geo, &
                                el_con_tau_crys_geo, work_base
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork

REAL(DP) :: omega, at(3,3)
INTEGER  :: igeom, ipol, na

igeom= (iwork-1)/work_base + 1
ibrav=el_con_ibrav_geo(igeom)
celldm(:)=el_con_celldm_geo(:,igeom)
CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
at=at/celldm(1)
tau(:,:)=el_con_tau_crys_geo(:,:,igeom)
CALL cryst_to_cart( nat, tau, at, 1 )
CALL set_equilibrium_conf( celldm, tau, at, omega )

WRITE(stdout, '(/,80("*"))')
WRITE(stdout, '(5x,"Computing geometry ", i5)') igeom
WRITE(stdout, '(5x,"Computing the elastic constant for celldm")')
WRITE(stdout, '(5x,6f12.5)') celldm(:)
WRITE(stdout, '(5x,"Unit cell volume",f12.5)') omega
WRITE(stdout, '(5x,"Cartesian axes")')
WRITE(stdout, '(5x,"site n.     atom                  positions (alat units)")')
WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
     (na, atm(ityp(na)), na, (tau(ipol,na), ipol=1,3), na=1,nat)
WRITE(stdout,'(80("*"),/)')

RETURN
END SUBROUTINE set_geometry_el_cons
