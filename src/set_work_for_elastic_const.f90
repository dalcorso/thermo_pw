!
! Copyright (C) 2013-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE set_work_for_elastic_const(iwork)
!-----------------------------------------------------------------------
USE kinds,            ONLY : DP
!
!  variables of thermo_pw
!
USE control_thermo,   ONLY : outdir_thermo
USE equilibrium_conf, ONLY : celldm0, at0, tau0_crys
USE thermo_mod,       ONLY : celldm_geo, ibrav_geo, tau_geo
USE control_elastic_constants, ONLY : frozen_ions, elastic_algorithm,  &
                             rot_mat, ngeom, tau_acc
USE control_thermo,   ONLY : ltau_from_file
!
!  library routines
!
USE elastic_constants, ONLY : epsilon_geo
USE strain_mod,       ONLY : apply_strain, print_strain
USE rotate,           ONLY : rotate_vect
!
!  variables of qe
!
USE cell_base,        ONLY : cell_base_init, at
USE ions_base,        ONLY : tau, nat
USE io_global,        ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: iwork

INTEGER  :: ivec, na, ipol, jpol, ibrav, irot
REAL(DP), ALLOCATABLE :: tau_ocoord(:,:)
REAL(DP) :: rd_ht(3,3), zero, celldm_(6)
LOGICAL  :: trd_ht
CHARACTER(LEN=256) :: outdir
CHARACTER(LEN=10)  :: cell_units
CHARACTER(LEN=6)   :: int_to_char

WRITE(stdout,'(/,2x,76("-"))')
CALL print_strain(epsilon_geo(:,:,iwork))

IF (ngeom>1) CALL set_geometry_el_cons(iwork)
!
!  entering here we have:
!  
!  at0 that contains the unstrained at vectors in units of celldm0(1)
!
!  tau0_crys that contains the crystal coordinates of the atoms in the
!  basis of the at0. In a uniform strain these coordinates do not change.
!
!  first strain the at0. 
!
DO ivec=1, 3
   CALL apply_strain(at0(1,ivec), at(1,ivec), epsilon_geo(1,1,iwork))
ENDDO
!
!   Now find tau strained in cartesian coordinates using the strained at
!
tau=tau0_crys
CALL cryst_to_cart( nat, tau, at, 1 )
!
!  here there is the possibility to add a fixed delta tau to the strained
!  coordinates. This has to be requested by the user.
!
tau(:,:)=tau(:,:)+tau_acc(:,:,iwork)
!
zero=0.0_DP
IF (elastic_algorithm=='standard'.OR.elastic_algorithm=='energy_std') THEN
   ibrav=0
   rd_ht = TRANSPOSE( at )
   trd_ht=.TRUE.
   cell_units='alat'
   CALL cell_base_init ( ibrav, celldm0, zero, zero, zero, zero, &
                     zero, zero, trd_ht, rd_ht, cell_units )
!
!  the atomic coordinates are strained uniformely and are not modified here
!  except when they are read from file
!
   IF (ltau_from_file) tau(:,:)=tau_geo(:,:,iwork)
ELSEIF (elastic_algorithm=='advanced' .OR. &
                                 elastic_algorithm=='energy') THEN
!
!  compute the at on the basis of ibrav_geo and celldm_geo
!
   ibrav = ibrav_geo(iwork)
   celldm_(:)=celldm_geo(:,iwork)
   trd_ht=.FALSE.
   rd_ht=0.0_DP
   CALL cell_base_init ( ibrav, celldm_, zero, zero, zero, zero, &
                         zero, zero, .FALSE., rd_ht, ' ' )
!
!   In this scheme sometimes the cartesian axes of the strained 
!   and unstrained cells are different. We rotate all the atomic positions
!   already strained to the new axis.
!
   ALLOCATE(tau_ocoord(3,nat))
   tau_ocoord=tau
   CALL rotate_vect(rot_mat(1,1,iwork), nat, tau_ocoord, tau, 1)
   DEALLOCATE(tau_ocoord)
!
!  bring the tau in the correct units of the new alat
!
   tau=tau * celldm0(1) / celldm_(1)
!
!  Superseed the previous calculation if atomic positions have been
!  read from file
!
   IF (ltau_from_file) tau(:,:)=tau_geo(:,:,iwork)
!
!  find the optimal fft mesh
!
   CALL find_fft_fact()
ENDIF
CALL set_fft_mesh()
!
!  Set the tmp_dir directory for this geometry
!
outdir=TRIM(outdir_thermo)//'/g'//TRIM(int_to_char(iwork))//'/'
CALL set_tmp_dir( outdir )
IF (.NOT.frozen_ions) CALL clean_bfgs_history()

RETURN
END SUBROUTINE set_work_for_elastic_const
