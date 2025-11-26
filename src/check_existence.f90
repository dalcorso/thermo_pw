!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
SUBROUTINE check_existence(iwork, part, run)
!--------------------------------------------------------
!
!   This routine checks if on file there are the information of
!   the self-consitent run corresponding to iwork, part.
!   If it finds them, it reads the energy and the stress if necessary
!   and sets run to .FALSE.. If it does not find them it sets
!   run to .TRUE. so that pw.x is called.
!
!   This routine can be called separately by each image
!
  USE kinds,           ONLY : DP
  USE constants,       ONLY : ry_kbar

  USE control_thermo,  ONLY : lstress, lberry

  USE thermo_mod,      ONLY : energy_geo
  USE elastic_constants, ONLY : sigma_geo
  USE polarization_vector, ONLY : mod_tot
  USE piezoelectric_tensor, ONLY : polar_strain, tot_b_phase

  USE io_global,       ONLY : ionode, ionode_id, stdout
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iwork, part
  LOGICAL, INTENT(OUT) :: run

  INTEGER  :: iu_ene, ios, ipol, jpol
  INTEGER  :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename
  REAL(DP) :: energy, stress(3,3), polar(3), tot_b_ph(3)
  LOGICAL :: exst
  run=.TRUE.
  
  IF (ionode) THEN
     iu_ene=find_free_unit()
     filename='restart/e_work_part.'//TRIM(int_to_char(iwork))//'.'//&
                              TRIM(int_to_char(part))
     INQUIRE(FILE=TRIM(filename),EXIST=exst)
     IF (exst) THEN
        WRITE(stdout,'(5x,"Data found on file")')
        OPEN(UNIT=iu_ene, FILE=TRIM(filename), STATUS='OLD', FORM='FORMATTED', &
            ERR=20, IOSTAT=ios)

        READ(iu_ene,*) energy
        WRITE(stdout,'(5x,"Total energy = ",f20.8," Ry")') energy
        IF (iwork > 0) THEN
           IF (lstress(iwork)) THEN
              WRITE(stdout,'(5x,"Stress (kbar)")')
              DO ipol=1,3
                 READ(iu_ene,*) (stress(ipol,jpol),jpol=1,3)
                 WRITE(stdout,'(3f15.7)') (stress(ipol,jpol)*ry_kbar,jpol=1,3)
              ENDDO
           END IF
           IF (lberry(iwork)) THEN
              WRITE(stdout,'(5x,"Polarization (e/bohr**2)")')
              READ(iu_ene,*) (polar(ipol),ipol=1,3)
              WRITE(stdout,'(3f15.7)') (polar(ipol),ipol=1,3)
              WRITE(stdout,'(5x,"Berry phase")')
              READ(iu_ene,*) (tot_b_ph(ipol),ipol=1,3)
              WRITE(stdout,'(3f15.7)') (tot_b_ph(ipol),ipol=1,3)
              IF (mod_tot==2) then
                 DO ipol=1,3
                    tot_b_ph(ipol)=tot_b_ph(ipol)-2.d0*nint(tot_b_ph(ipol)/2.d0)
                 ENDDO
              ELSE
                 DO ipol=1,3
                    tot_b_ph(ipol)=tot_b_ph(ipol)-nint(tot_b_ph(ipol))
                 ENDDO
              ENDIF
           END IF
        ELSE
           WRITE(stdout,'(5x,"Stress (kbar)")')
           DO ipol=1,3
              READ(iu_ene,*) (stress(ipol,jpol),jpol=1,3)
              WRITE(stdout,'(3f15.7)') (stress(ipol,jpol)*ry_kbar,jpol=1,3)
           ENDDO
        END IF
        run=.FALSE.
        CLOSE(iu_ene)
     ENDIF
!
!  If the file cannot be opened it is considered unexistent and the program
!  continue
!
20   CONTINUE
  ENDIF
  CALL mp_bcast(run, ionode_id, intra_image_comm)
  IF (.NOT.run .AND. iwork > 0) THEN
     CALL mp_bcast(energy, ionode_id, intra_image_comm)
     energy_geo(iwork)=energy
     IF (lstress(iwork)) THEN
        CALL mp_bcast(stress, ionode_id, intra_image_comm)
        sigma_geo(:,:,iwork)=stress(:,:)
     END IF
     IF (lberry(iwork)) THEN
        CALL mp_bcast(polar, ionode_id, intra_image_comm)
        polar_strain(:,iwork)=polar(:)
        CALL mp_bcast(tot_b_ph, ionode_id, intra_image_comm)
        tot_b_phase(:,iwork)=tot_b_ph(:)
     END IF
  END IF

  RETURN
END SUBROUTINE check_existence

!--------------------------------------------------------
SUBROUTINE save_existence(iwork, part)
!--------------------------------------------------------
!
!  This routine saves on file the total energy, and
!  possibly the stress of a self-consistent run. 
!
  USE control_thermo,  ONLY : lstress, lberry

  USE thermo_mod,      ONLY : energy_geo, what
  USE elastic_constants, ONLY : sigma_geo
  USE piezoelectric_tensor, ONLY : polar_strain, tot_b_phase

  USE ener,            ONLY : etot
  USE force_mod,       ONLY : sigma
  USE ions_base,       ONLY : tau, nat
  USE control_elastic_constants, ONLY : tau_save_ec

  USE io_global,       ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iwork, part
  INTEGER :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename
  INTEGER :: iu_ene, ios, ipol, jpol

  IF (ionode) THEN
     iu_ene=find_free_unit()
     filename='restart/e_work_part.'//TRIM(int_to_char(iwork))//'.'//&
                              TRIM(int_to_char(part))
     OPEN(UNIT=iu_ene, FILE=TRIM(filename), STATUS='UNKNOWN', &
                       FORM='FORMATTED', ERR=20, IOSTAT=ios)

     IF (iwork > 0) THEN
        WRITE(iu_ene,*) energy_geo(iwork)
        IF (lstress(iwork)) THEN
           DO ipol=1,3
              WRITE(iu_ene,*) (sigma_geo(ipol,jpol,iwork),jpol=1,3)
           ENDDO
        ENDIF
        IF (lberry(iwork)) THEN
           WRITE(iu_ene,*) (polar_strain(ipol,iwork), ipol=1,3)
           WRITE(iu_ene,*) (tot_b_phase(ipol,iwork), ipol=1,3)
        ENDIF
     ELSE
        WRITE(iu_ene,*) etot
        DO ipol=1,3
           WRITE(iu_ene,*) (sigma(ipol,jpol),jpol=1,3)
        ENDDO
     ENDIF
     CLOSE(iu_ene)
!
!  If the file cannot be written we do nothing
!
20   CONTINUE
  END IF

  IF (what=='elastic_constants_geo') THEN
     tau_save_ec(1:3,1:nat,iwork)=tau(1:3,1:nat)
  ENDIF
  
  RETURN
END SUBROUTINE save_existence
!
!--------------------------------------------------------
SUBROUTINE save_geometry(iwork, part, iwho)
!--------------------------------------------------------
!
!  This routine saves on file the crystal parameters, and
!  the relaxed atomic coordinates. 
!  The tau are saved as they are in cartesian coordinates 
!  in units of alat.
!
  USE cell_base,       ONLY : ibrav, celldm, at
  USE ions_base,       ONLY : atm, tau, nat, ityp

  USE io_global,       ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iwork, part, iwho
  INTEGER :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename, label
  INTEGER :: iu_geo, ios, na, ipol

  IF (iwho==0) RETURN
  IF (ionode) THEN
     iu_geo=find_free_unit()
     IF (iwho==1) label='restart/geo_work_part.mlc.'
     IF (iwho==2) label='restart/geo_work_part.ecg.'
     filename=TRIM(label)//TRIM(int_to_char(iwork))//'.'//&
                              TRIM(int_to_char(part))
     OPEN(UNIT=iu_geo, FILE=TRIM(filename), STATUS='UNKNOWN', &
                       FORM='FORMATTED', ERR=20, IOSTAT=ios)

     IF (iwork > 0) THEN
        WRITE(iu_geo,*) ibrav
        WRITE(iu_geo,*) celldm(:)
        WRITE(iu_geo,*) 
        DO ipol=1,3
           WRITE(iu_geo, *) at(:,ipol)
        ENDDO
        WRITE(iu_geo,*) nat
        WRITE(iu_geo,*) (ityp(na), na=1,nat)
        DO na=1,nat
           WRITE(iu_geo,'(a6,3e23.14)') atm(ityp(na)), tau(1,na), tau(2,na), &
                                                                  tau(3,na)
        ENDDO
     ENDIF
     CLOSE(iu_geo)
!
!  If the file cannot be written we do nothing
!
20   CONTINUE
  END IF

  RETURN
END SUBROUTINE save_geometry
!
!--------------------------------------------------------
SUBROUTINE check_geometry_exist(iwork, part, iwho)
!--------------------------------------------------------
!
!  This routine reads from file the crystal parameters, and
!  the relaxed atomic coordinates. 
!  The tau are supposed to be in cartesian coordinates 
!  in units of alat.
!
  USE thermo_mod,      ONLY : ibrav_geo, celldm_geo, tau_geo, at_geo
  USE initial_conf,    ONLY : atm_save
  USE ions_base,       ONLY : nat, ityp

  USE io_global,       ONLY : ionode, ionode_id
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iwork, part, iwho
  INTEGER :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename, label
  INTEGER :: iu_geo, ios, na, ipol

  IF (iwho==0) RETURN
  IF (ionode) THEN
     ios=0
     iu_geo=find_free_unit()
     IF (iwho==1) label='restart/geo_work_part.mlc.'
     IF (iwho==2) label='restart/geo_work_part.ecg.'
     filename=TRIM(label)//TRIM(int_to_char(iwork))//'.'//&
                              TRIM(int_to_char(part))
     OPEN(UNIT=iu_geo, FILE=TRIM(filename), STATUS='OLD', &
                       FORM='FORMATTED', ERR=20, IOSTAT=ios)

     IF (iwork > 0) THEN
        READ(iu_geo, *, ERR=20, IOSTAT=ios) ibrav_geo(iwork)
        READ(iu_geo, *, ERR=20, IOSTAT=ios) celldm_geo(:,iwork)
        READ(iu_geo, *, ERR=20, IOSTAT=ios) 
        DO ipol=1,3
           READ(iu_geo, *, ERR=20, IOSTAT=ios) at_geo(:,ipol,iwork)
        ENDDO
        READ(iu_geo, *, ERR=20, IOSTAT=ios) nat
        READ(iu_geo, *, ERR=20, IOSTAT=ios) (ityp(na), na=1,nat)
        DO na=1,nat
           READ(iu_geo,*, ERR=20, IOSTAT=ios) &
                          atm_save(ityp(na)), tau_geo(1,na,iwork), &
                          tau_geo(2,na,iwork), tau_geo(3,na,iwork)
        ENDDO
     ENDIF
     CLOSE(iu_geo)
!
!  If the file cannot be read we stop. This routine is called only 
!  if explicitely required by the user that must check that the file
!  is present. The file is produced after a previous self-consistent
!  calculation on this geometry.
!
  END IF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
20  CALL errore('check_geometry_exist',&
         'geometry file unreadable or inexistent',ABS(ios))
  CALL mp_bcast(ibrav_geo(iwork), ionode_id, intra_image_comm)
  CALL mp_bcast(celldm_geo(:,iwork), ionode_id, intra_image_comm)
  CALL mp_bcast(nat, ionode_id, intra_image_comm)
  CALL mp_bcast(ityp, ionode_id, intra_image_comm)
  CALL mp_bcast(atm_save, ionode_id, intra_image_comm)
  CALL mp_bcast(at_geo(:,:,iwork), ionode_id, intra_image_comm)
  CALL mp_bcast(tau_geo(:,:,iwork), ionode_id, intra_image_comm)

  RETURN
END SUBROUTINE check_geometry_exist

!--------------------------------------------------------
SUBROUTINE check_geometry_el_cons_exist(iwork, part)
!--------------------------------------------------------
!
!  This routine reads from file the crystal parameters, and
!  the relaxed atomic coordinates. 
!  The tau are supposed to be in cartesian coordinates 
!  in units of alat. It sets them in the unperturbed geometries
!  for elastic_constants_geo calculation
!
  USE kinds,  ONLY : DP
  USE control_elastic_constants,  ONLY : el_con_ibrav_geo, el_con_celldm_geo, &
                                    el_con_tau_crys_geo, el_con_omega_geo,    &
                                    el_con_tau_geo, el_con_at_geo
  USE initial_conf,    ONLY : atm_save
  USE ions_base,       ONLY : nat, ityp

  USE io_global,       ONLY : ionode, ionode_id
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iwork, part
  INTEGER :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(LEN=256) :: filename
  REAL(DP) :: tau_(3,nat), celldm_(6), at_(3,3), bg_(3,3), omega_, &
              compute_omega_geo
  INTEGER :: iu_geo, ios, na, ibrav_, ipol

  IF (ionode) THEN
     ios=0
     iu_geo=find_free_unit()
     filename='restart/geo_work_part.mlc.'//TRIM(int_to_char(iwork))//'.'//&
                              TRIM(int_to_char(part))
     OPEN(UNIT=iu_geo, FILE=TRIM(filename), STATUS='OLD', &
                       FORM='FORMATTED', ERR=20, IOSTAT=ios)

     IF (iwork > 0) THEN
        READ(iu_geo,*, ERR=20, IOSTAT=ios) el_con_ibrav_geo(iwork)
        READ(iu_geo,*, ERR=20, IOSTAT=ios) el_con_celldm_geo(:,iwork)
        READ(iu_geo,*, ERR=20, IOSTAT=ios) 
        DO ipol=1,3
           READ(iu_geo, *, ERR=20, IOSTAT=ios) el_con_at_geo(:,ipol,iwork)
        ENDDO
        READ(iu_geo,*, ERR=20, IOSTAT=ios) nat
        READ(iu_geo,*, ERR=20, IOSTAT=ios) (ityp(na), na=1,nat)
        DO na=1,nat
           READ(iu_geo,*, ERR=20, IOSTAT=ios) &
                          atm_save(ityp(na)), tau_(1,na), &
                          tau_(2,na), tau_(3,na)
        ENDDO
        el_con_tau_geo(:,:,iwork)=tau_(:,:)
        el_con_tau_crys_geo(:,:,iwork)=tau_(:,:)
        celldm_(1:6)=el_con_celldm_geo(1:6,iwork)
        ibrav_=el_con_ibrav_geo(iwork)
        at_(:,:)=el_con_at_geo(:,:,iwork)
        CALL recips(at_(1,1),at_(1,2),at_(1,3),bg_(1,1),bg_(1,2),bg_(1,3))
!
!   put atomic coordinates in the crystal basis
!
        CALL cryst_to_cart( nat, el_con_tau_crys_geo(:,:,iwork), bg_, -1 )
        el_con_omega_geo(iwork)=compute_omega_geo(ibrav_,celldm_)
     ENDIF
     CLOSE(iu_geo)
!
!  If the file cannot be read we stop. This routine is called only 
!  if explicitely required by the user that must check that the file
!  is present. The file is produced after a previous self-consistent
!  calculation on this geometry .
!
  END IF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
20  CALL errore('check_geometry_el_cons_exist',&
         'geometry file unreadable or inexistent',ABS(ios))
  CALL mp_bcast(el_con_ibrav_geo(iwork), ionode_id, intra_image_comm)
  CALL mp_bcast(el_con_celldm_geo(:,iwork), ionode_id, intra_image_comm)
  CALL mp_bcast(el_con_omega_geo(iwork), ionode_id, intra_image_comm)
  CALL mp_bcast(el_con_tau_crys_geo(:,:,iwork), ionode_id, intra_image_comm)
  CALL mp_bcast(el_con_tau_geo(:,:,iwork), ionode_id, intra_image_comm)

  RETURN
END SUBROUTINE check_geometry_el_cons_exist

