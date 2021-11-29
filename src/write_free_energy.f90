!
! Copyright (C) 2021 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE write_free_energy()
!----------------------------------------------------------------------------
!
!  This subroutines saves on file the vibrational free energy and
!  (if calculated) the electronic free energy at the ntemp_plot 
!  temperatures.
!
USE kinds,             ONLY : DP
USE thermo_mod,        ONLY : omega_geo, energy_geo, ngeo, no_ph
USE temperature,       ONLY : ntemp, temp, ntemp_plot, itemp_plot
USE thermodynamics,    ONLY : ph_free_ener
USE el_thermodynamics, ONLY : el_free_ener
USE control_eldos,     ONLY : lel_free_energy
USE data_files,        ONLY : flevdat
USE io_global,         ONLY : ionode

IMPLICIT NONE
CHARACTER(LEN=256) :: filedata
CHARACTER(LEN=8) :: float_to_char
INTEGER :: itempp, itemp, idata, iu_free
INTEGER :: find_free_unit

DO itempp=1,ntemp_plot
   itemp=itemp_plot(itempp)
   filedata="anhar_files/"//TRIM(flevdat)//"_free."//&
                                      TRIM(float_to_char(temp(itemp),1))
   CALL add_pressure(filedata)

   IF (ionode) THEN
      iu_free=find_free_unit()
      OPEN(UNIT=iu_free, FILE=TRIM(filedata), STATUS='UNKNOWN', &
                                                           FORM='FORMATTED')
      IF (lel_free_energy) THEN
         WRITE(iu_free,'("# Volume (a.u)^3,",10x,"energy (Ry)",10x,&
                          &"ph free-energy (Ry)",10x,"el free-energy (Ry)")')
      ELSE
         WRITE(iu_free,'("# Volume (a.u)^3,",17x,"energy (Ry)",10x,&
                                                   &"ph free-energy (Ry)")')
      ENDIF
      DO idata=1,ngeo(1)
         IF (no_ph(idata)) CYCLE
         IF (lel_free_energy) THEN
            WRITE(iu_free,'(4f20.12)') omega_geo(idata), energy_geo(idata), &
                        ph_free_ener(itemp,idata), el_free_ener(itemp,idata)
         ELSE
            WRITE(iu_free,'(3f25.12)') omega_geo(idata), energy_geo(idata),&
                                       ph_free_ener(itemp,idata)
         ENDIF
      ENDDO
      CLOSE(UNIT=iu_free, STATUS='KEEP')
   ENDIF
ENDDO

RETURN
END SUBROUTINE write_free_energy
