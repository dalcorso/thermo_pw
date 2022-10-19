!
! Copyright (C) 2022 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE manage_hot_electrons()
!---------------------------------------------------------------------
!
!  This routine reads from files contained in a directory
!  restart2, restart3, ... restart#nsigma the energies of the 
!  ngeo geometries calculated at several temperatures (i.e. smearing sigma) 
!  with the Fermi-Dirac smearing and fits these energies with a polinomial 
!  of temperature. From the derivatives with respect to
!  temperature it gets the contribution of the excited electrons
!  to the energy, free energy, entropy and isochoric heat capacity at
!  any temperature.
!  These quantities are then added to the free energy in other routines.
!
USE kinds,         ONLY : DP
USE constants,     ONLY : rydberg_si, k_boltzmann_si
USE thermo_mod,    ONLY : ngeo, energy_geo
USE el_thermodynamics, ONLY : el_free_ener, el_ener, el_entr, el_ce
USE polyfit_mod,   ONLY : polyfit, compute_poly, compute_poly_deriv, &
                          compute_poly_deriv2
USE temperature,   ONLY : temp, ntemp, temp_sigma, sigma_ry
USE control_conv,  ONLY : nsigma

USE mp_images,     ONLY : intra_image_comm
USE mp,            ONLY : mp_bcast
USE io_global,     ONLY : ionode, ionode_id, stdout

IMPLICIT NONE

CHARACTER(LEN=256) :: direc, filename
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: exst
INTEGER :: iu_ene, igeom, isigma, itemp, ios, degree

INTEGER :: find_free_unit

REAL(DP) :: a(nsigma), y(nsigma)

ALLOCATE(temp_sigma(nsigma))

degree=nsigma-1
DO isigma=1,nsigma
   temp_sigma(isigma)= sigma_ry(isigma) * rydberg_si / k_boltzmann_si
!   WRITE(6,*) 'isigma, temp', isigma, temp_sigma(isigma)
ENDDO

IF (ionode) THEN
   iu_ene=find_free_unit()
   DO igeom=1,ngeo(1)
      DO isigma=2,nsigma
         direc='restart'//TRIM(int_to_char(isigma))
         filename=TRIM(direc)//'/e_work_part.'//TRIM(int_to_char(igeom))//'.1'
         INQUIRE(FILE=TRIM(filename),EXIST=exst)
         IF (exst) THEN
            WRITE(stdout,'(5x,"Data found on file")')
            OPEN(UNIT=iu_ene, FILE=TRIM(filename), STATUS='OLD', &
                 FORM='FORMATTED', ERR=20, IOSTAT=ios)
            READ(iu_ene,*) y(isigma)
            WRITE(stdout,'(5x,"Sigma=",f20.8," Total energy = ",f20.8," Ry")')&
                                  sigma_ry(isigma), y(isigma)
            CLOSE(UNIT=iu_ene, STATUS='KEEP')
         ELSE
            CALL errore('manage_hot_electrons','missing energy',1)
         ENDIF
20       CALL errore('manage_hot_electrons','problem opening file',ios)
      ENDDO
      y(1)=energy_geo(igeom)
!
!   Interpolate with a polynomial the free energy as a function of sigma
!   and computes the contribution to the energy, free_energy, entropy
!   and isochoric heat capacity.
!
      CALL polyfit(temp_sigma,y,nsigma,a,degree)
!      WRITE(6,*) 'a', a(1), a(2), a(3)
      DO itemp=1,ntemp
         CALL compute_poly(temp(itemp), degree, a, el_free_ener(itemp,igeom))
         CALL compute_poly_deriv(temp(itemp), degree, a, &
                                                        el_entr(itemp,igeom))
         CALL compute_poly_deriv2(temp(itemp), degree, a, el_ce(itemp,igeom))
         el_entr(itemp,igeom)=-el_entr(itemp,igeom)
         el_ener(itemp,igeom)=el_free_ener(itemp,igeom) + &
                                            temp(itemp) * el_entr(itemp,igeom)
         el_ce(itemp,igeom)= -temp(itemp) * el_ce(itemp,igeom)
         el_free_ener(itemp,igeom)=el_free_ener(itemp,igeom)-energy_geo(igeom)
         el_ener(itemp,igeom)=el_ener(itemp,igeom)-energy_geo(igeom)
!         WRITE(6,*) temp(itemp), el_free_ener(itemp,igeom), &
!          el_ener(itemp,igeom), el_entr(itemp,igeom), el_ce(itemp,igeom)
      ENDDO
   ENDDO
ENDIF
!
!   Broadcast the results to all processors within the image
!
CALL mp_bcast(el_free_ener,ionode_id,intra_image_comm)
CALL mp_bcast(el_ener,ionode_id,intra_image_comm)
CALL mp_bcast(el_entr,ionode_id,intra_image_comm)
CALL mp_bcast(el_ce,ionode_id,intra_image_comm)

RETURN
END SUBROUTINE manage_hot_electrons

