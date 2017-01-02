!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE init_elastic_constants_t( )
!
!  This routine copies in its internal variables the lattice parameters
!  of each geometry decided by the energy minimizer. These are all the
!  geometries that must be computed by the elastic constants routine.
!  All these geometries are the unperturbed geometries of the elastic 
!  constants calculation.
!
USE kinds,      ONLY : DP
USE thermo_mod, ONLY : tot_ngeo
USE ions_base,  ONLY : tau, nat
USE control_elastic_constants, ONLY : el_con_geo, el_con_ibrav_geo,  &
                                      el_con_celldm_geo, el_con_tau_geo
USE thermo_mod, ONLY : ibrav_geo, celldm_geo
IMPLICIT NONE
INTEGER :: igeo

ALLOCATE(el_con_geo(6,6,tot_ngeo))
ALLOCATE(el_con_ibrav_geo(tot_ngeo))
ALLOCATE(el_con_celldm_geo(6,tot_ngeo))
ALLOCATE(el_con_tau_geo(3,nat,tot_ngeo))

DO igeo=1,tot_ngeo
   el_con_ibrav_geo(igeo)=ibrav_geo(igeo)
   el_con_celldm_geo(:,igeo)=celldm_geo(:,igeo)
   el_con_tau_geo(:,:,igeo)=tau(:,:)
ENDDO

RETURN
END SUBROUTINE init_elastic_constants_t

SUBROUTINE set_geometry_el_cons(igeom)
!
!  This routine receives as input the current unperturbed geometry
!  and set the variables of pw with this unperturbed geometry.
!
USE kinds, ONLY : DP
USE cell_base, ONLY : ibrav, celldm
USE ions_base, ONLY : tau, nat, atm, ityp
USE control_elastic_constants, ONLY : el_con_ibrav_geo, el_con_celldm_geo, &
                                      el_con_tau_geo
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

REAL(DP) :: omega, at(3,3)
INTEGER  :: ipol, na

ibrav=el_con_ibrav_geo(igeom)
celldm(:)=el_con_celldm_geo(:,igeom)
tau(:,:)=el_con_tau_geo(:,:,igeom)
CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
at=at/celldm(1)
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
