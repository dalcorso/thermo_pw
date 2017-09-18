!
! Copyright (C) 2013 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE deallocate_q2r()
  !-----------------------------------------------------------------------
  !
  !  This routine deallocates the variables allocated in 
  !  manage_ph_dispersions
  !
  USE thermo_mod,     ONLY : tot_ngeo
  USE phdos_module,   ONLY : destroy_phdos
  USE ph_freq_module, ONLY : destroy_ph_freq
  USE ph_freq_thermodynamics, ONLY : ph_freq_save
  USE thermodynamics, ONLY : phdos_save
  USE control_dosq,   ONLY : dos_q, dos_wq
  USE control_thermo, ONLY : lph, ltherm, ltherm_dos

  IMPLICIT NONE
  INTEGER :: igeom
  !

  IF ( ALLOCATED (dos_q) )           DEALLOCATE(dos_q)
  IF ( ALLOCATED (dos_wq) )          DEALLOCATE(dos_wq)

  IF (lph) THEN
     IF (ltherm.AND.ltherm_dos) THEN
        IF (ALLOCATED(phdos_save)) THEN
           DO igeom=1,tot_ngeo
              CALL destroy_phdos(phdos_save(igeom))
           ENDDO
           DEALLOCATE(phdos_save)
        ENDIF
        IF (ALLOCATED(ph_freq_save)) THEN
           DO igeom=1,tot_ngeo
              CALL destroy_ph_freq(ph_freq_save(igeom))
           ENDDO
           DEALLOCATE(ph_freq_save)
        ENDIF
     ENDIF
  ENDIF
  ! 
  RETURN
  !
END SUBROUTINE deallocate_q2r
