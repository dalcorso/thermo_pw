!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE manage_elastic_cons(nwork, igeom)

USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : energy_geo, density
USE control_elastic_constants, ONLY : ngeo_strain, frozen_ions,            &
                              elastic_algorithm, rot_mat, elcpvar
USE initial_conf,      ONLY : ibrav_save
USE thermo_sym,        ONLY : laue
USE elastic_constants, ONLY : print_elastic_constants,                     &
                              compute_elastic_constants, epsilon_geo,      &
                              sigma_geo, el_con, el_compliances,           &
                              compute_elastic_compliances,                 &
                              print_elastic_compliances, read_elastic,     &
                              write_elastic, print_macro_elasticity,       &
                              compute_elastic_constants_adv,               &
                              compute_elastic_constants_ene,               &
                              print_sound_velocities
USE equilibrium_conf,  ONLY : omega0
USE control_macro_elasticity, ONLY : macro_el, vp, vb, vg, approx_debye_t
USE debye_module,      ONLY : compute_debye_temperature,                   &
                              compute_debye_temperature_poisson
USE control_debye,     ONLY : debye_t
USE data_files,        ONLY : fl_el_cons
USE ions_base,         ONLY : nat
USE cell_base,         ONLY : omega
USE mp_images,         ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: nwork, igeom

REAL(DP) :: poisson, bulkm
LOGICAL  :: exst
CHARACTER(LEN=6)    :: int_to_char
CHARACTER(LEN=256)  :: filelastic
!
!  the elastic constants are calculated here
!
IF (elastic_algorithm=='standard') THEN
   CALL compute_elastic_constants(sigma_geo, epsilon_geo, nwork, &
                          ngeo_strain, ibrav_save, laue, elcpvar)
ELSE IF (elastic_algorithm=='advanced') THEN
   CALL compute_elastic_constants_adv(sigma_geo, epsilon_geo, &
                               nwork, ngeo_strain, ibrav_save, laue, rot_mat, &
                                                   elcpvar)
ELSE IF (elastic_algorithm=='energy') THEN
   CALL compute_elastic_constants_ene(energy_geo, epsilon_geo, &
                               nwork, ngeo_strain, ibrav_save, laue, omega0, &
                                                               elcpvar)
END IF
CALL print_elastic_constants(el_con, frozen_ions)
!
!  now compute the elastic compliances and prints them
!
CALL compute_elastic_compliances(el_con,el_compliances)
CALL print_elastic_compliances(el_compliances, frozen_ions)
CALL print_macro_elasticity(ibrav_save,el_con,el_compliances,macro_el, .TRUE.)
!
!  here compute the sound velocities, using the density of the solid and
!  the elastic constants
!
CALL print_sound_velocities( ibrav_save, el_con, el_compliances, &
                                       density, vp, vb, vg )
!
!  here we compute the Debye temperature approximatively from the
!  poisson ratio and the bulk modulus
!
poisson=(macro_el(4)+macro_el(8) ) * 0.5_DP
bulkm=(macro_el(1)+macro_el(5) ) * 0.5_DP
CALL compute_debye_temperature_poisson(poisson, bulkm, &
                            density, nat, omega, approx_debye_t)
!
!  compute the Debye temperature and the thermodynamic quantities
!  within the Debye model
!
CALL compute_debye_temperature(el_con, density, nat, omega, debye_t)
CALL write_thermo_debye(igeom)
CALL plot_thermo_debye(igeom)
!
!  save elastic constants and compliances on file
!
filelastic='elastic_constants/'//TRIM(fl_el_cons)//'.g'//&
                                                     TRIM(int_to_char(igeom))
IF (my_image_id==root_image) CALL write_elastic(filelastic)
RETURN
END SUBROUTINE manage_elastic_cons

