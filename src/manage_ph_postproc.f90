!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE manage_ph_postproc(igeom)
!-----------------------------------------------------------------------
!
!  This driver computes the interatomic force constants and
!  interpolates the phonon frequencies to make a phonon dispersion
!  and then computes them on a thick mesh to compute the harmonic
!  thermodynamic quantities.
!
USE kinds,            ONLY : DP
USE constants,        ONLY : amu_ry, ry_to_cmm1
USE thermo_mod,       ONLY : z_geo, freq_geo, omega_geo, epsilon_zero_geo, what
USE control_thermo,   ONLY : ltherm, ltherm_dos, ltherm_freq, set_internal_path
USE ifc,              ONLY : has_zstar, epsil_ifc, zeu
USE ions_base,        ONLY : nat, ityp, amass
USE control_epsilon_infty, ONLY : lepsilon_infty_geo, lzeu_geo
USE ph_freq_thermodynamics, ONLY : ph_freq_save
USE thermo_mod,       ONLY : epsilon_infty_geo, zeu_geo
USE control_paths,    ONLY : disp_nqs
USE control_phrun,    ONLY : auxdyn
USE control_atomic_pos, ONLY : linterpolate_tau
USE data_files,       ONLY : fl_dielectric
USE dielectric_constant, ONLY : write_dielectric_properties_to_file, &
                                polar_mode_permittivity_tpw
USE mp_images,          ONLY : my_image_id, root_image

IMPLICIT NONE
INTEGER, INTENT(IN) :: igeom

CHARACTER(LEN=256) :: filedata, filerap, fileout, gnu_filename, filenameps, &
                      filepbs, filename
CHARACTER(LEN=6) :: int_to_char
COMPLEX(DP), ALLOCATABLE :: z(:,:)
REAL(DP), ALLOCATABLE :: w2(:)
INTEGER :: imode, na, nta, ipol, jpol
!
!   Compute the interatomic force constants from the dynamical matrices
!   written on file
!
CALL q2r_sub(auxdyn) 
!
IF (linterpolate_tau) CALL save_tau_on_tpw(igeom)
!
!    compute interpolated dispersions
!
IF (set_internal_path) CALL set_bz_path()
CALL write_ph_dispersions()
CALL set_files_for_plot(2, ' ', filedata, filerap, fileout, &
                                            gnu_filename, filenameps, filepbs)
IF (disp_nqs>0) CALL plotband_sub(2, filedata, filerap, fileout, &
                                            gnu_filename, filenameps, filepbs)
!
!   Compute the harmonic thermodynamic quantities
!
IF (ltherm) THEN
!
!    the frequencies on a uniform mesh are interpolated or read from disk 
!    if available and saved in ph_freq_save 
!
   CALL write_ph_gamma(igeom)
   CALL write_ph_freq(igeom)
   IF (ltherm_freq.OR..NOT.ltherm_dos) CALL write_thermo_ph(igeom)
 
   IF (ltherm_dos) THEN
      CALL write_phdos(igeom)
      CALL plot_phdos()
      CALL write_thermo(igeom)
   ENDIF
   CALL plot_thermo()
   IF (has_zstar) THEN
      epsilon_infty_geo(:,:,igeom)=epsil_ifc(:,:)
      zeu_geo(:,:,:,igeom)=zeu(:,:,:)
      lepsilon_infty_geo(igeom)=.TRUE.
      lzeu_geo(igeom)=.TRUE.
!
!   Compute epsilon_zero. The routine requires the displacement
!
      ALLOCATE(w2(3*nat))
      ALLOCATE(z(3*nat,3*nat))
      DO imode = 1, 3*nat
         DO na = 1,nat
            nta = ityp(na)
            DO ipol = 1,3
               z((na-1)*3+ipol,imode) = z_geo((na-1)*3+ipol,imode,igeom)/ &
                                    SQRT(amu_ry*amass(nta))
            ENDDO
         ENDDO
         w2(imode)=SIGN((freq_geo(imode,igeom) /  RY_TO_CMM1)**2, &
                                               freq_geo(imode,igeom) )
      ENDDO
      CALL polar_mode_permittivity_tpw( nat, epsilon_infty_geo(1,1,igeom), &
            z, zeu_geo(1,1,1,igeom), w2,         &
            omega_geo(igeom), epsilon_zero_geo(1,1,igeom))
      DEALLOCATE(z)
      DEALLOCATE(w2)
      filename='elastic_constants/'//TRIM(fl_dielectric)//'.g'//&
                                              TRIM(int_to_char(igeom))
      IF (my_image_id==root_image) &
         CALL write_dielectric_properties_to_file(filename, nat, &
             epsilon_infty_geo(1,1,igeom), zeu_geo(1,1,1,igeom) )
   ENDIF

   IF (what/='mur_lc_t') DEALLOCATE(ph_freq_save(igeom)%nu)
ENDIF

CALL clean_ifc_variables()

RETURN
END SUBROUTINE manage_ph_postproc

