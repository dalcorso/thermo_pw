!
! Copyright (C) 2016-present Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE write_el_thermo(igeom)
!------------------------------------------------------------------------
!
!  This routine writes on file the electronic thermodynamical quantities
!
USE kinds,          ONLY : DP
USE constants,      ONLY : rytoev, electronvolt_si, rydberg_si, avogadro
USE eldos_module,   ONLY : eldos_type, el_free_energy, el_energy, el_entropy, &
                           el_specific_heat_cv, destroy_eldos, &
                           el_chem_pot, read_eldos_data
USE el_thermodynamics, ONLY : el_ener, el_free_ener, el_entr, el_mu, &
                           el_ce
USE ener,           ONLY : ef
USE lsda_mod,       ONLY : nspin
USE klist,          ONLY : degauss, nelec, ltetra
USE lsda_mod,       ONLY : lsda
USE temperature,    ONLY : ntemp, temp
USE mp_images,      ONLY : root_image, my_image_id, intra_image_comm
USE mp,             ONLY : mp_bcast, mp_sum
USE io_global,      ONLY : ionode, ionode_id, stdout
USE data_files,     ONLY : fleltherm, fleldos

IMPLICIT NONE
REAL(DP), PARAMETER :: caltoj=4.184_DP

INTEGER, INTENT(IN) :: igeom
INTEGER          :: itemp
INTEGER          :: iu_therm
TYPE(eldos_type) :: eldos

CHARACTER(LEN=256) :: fileeltherm, fileeldos, message1, message2

REAL(DP) :: ene0, mu0
REAL(DP) ::  dos1, dos2, dosef, udosef, ddosef, ddos1, ddos2, &
             uddosde, dddosde, ddosde
INTEGER :: n1, n2, n, ndos, idum, nstart, nlast
INTEGER :: find_free_unit
LOGICAL :: check_file_exists, do_read
!
!  check if the data are already on file
!
do_read=.FALSE.
fileeltherm="therm_files/"//TRIM(fleltherm)
IF (check_file_exists(fileeltherm)) do_read=.TRUE.
!
IF (my_image_id /= root_image) RETURN
!
!  Electrons contribute to the thermodynamical properties of metals only.
!  For insulators return.
!
IF (degauss==0.0_DP.AND..NOT.ltetra) RETURN

IF (do_read) THEN
   IF (ionode) THEN
      iu_therm=find_free_unit()
      fileeltherm='therm_files/'//TRIM(fleltherm)
      OPEN (UNIT=iu_therm, FILE=TRIM(fileeltherm), STATUS='old', &
                                                           FORM='formatted')
      DO idum=1,12
         READ(iu_therm,*)
      ENDDO
      DO itemp = 1, ntemp
         READ(iu_therm, '(e16.8,5e20.12)') temp(itemp), &
                    el_ener(itemp,igeom), el_free_ener(itemp,igeom), &
                    el_entr(itemp,igeom), el_ce(itemp,igeom),    &
                    el_mu(itemp,igeom)
      END DO
      CLOSE(iu_therm)
   END IF
   CALL mp_bcast(el_ener(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(el_free_ener(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(el_entr(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(el_ce(:,igeom), ionode_id, intra_image_comm)
   CALL mp_bcast(el_mu(:,igeom), ionode_id, intra_image_comm)
   RETURN
ENDIF

message1="     Computing the thermodynamic properties from electron dos"
WRITE(message2,'(5x,"Writing on file ",a)') TRIM(fleltherm)
CALL decorated1_write2(message1, message2)

fileeldos='therm_files/'//TRIM(fleldos)
CALL read_eldos_data(eldos, lsda, fileeldos)

!
!  The energy at very low temperature is considered as the zero of the energy
!  The default low temperature is 4. K or the lowest temperature
!  required in input.
!
CALL el_chem_pot(eldos, MIN(4.0_DP,temp(1)), nelec, mu0)
CALL el_energy(eldos, MIN(4.0_DP,temp(1)), mu0, ene0)
WRITE(stdout,'(/,5x, "Chemical potential (mu) at T=", f13.5," K  =", &
                   &f13.8, "  eV")') min(4.0_DP,temp(1)),  mu0 * rytoev

ndos=eldos%number_of_points
n1=1
n2=1
DO n=2,ndos
   IF (eldos%e(n) < mu0) n1=n
   IF (eldos%e(n) > mu0 .AND. n2==1) n2=n
END DO
IF (n1>=n2) CALL errore('write_el_thermo','some problems with dos or ef',1) 
dos1=eldos%dos(n1)
dos2=eldos%dos(n2)
dosef=dos1 + (dos2 - dos1) * ( mu0 - eldos%e(n1) ) / &
                                          ( eldos%e(n2) - eldos%e(n1) )
ddosde= (dos2 - dos1) / (eldos%e(n2) - eldos%e(n1))
IF (nspin==2) THEN
   ddos1 = eldos%ddos(n1)
   ddos2 = eldos%ddos(n2)
   udosef = dosef
   ddosef = ddos1 + (ddos2 - ddos1) * ( mu0 - eldos%e(n1) ) /  &
                                      ( eldos%e(n2) - eldos%e(n1) )
   uddosde = ddosde
   dddosde = (ddos2 - ddos1) / (eldos%e(n2) - eldos%e(n1))
   dosef= udosef + ddosef
   ddosde= uddosde + dddosde

   WRITE(stdout,'(/,5x,"Density of up states at mu =     ",&
                      &f15.8,"  states/(eV spin cell)")')  udosef / rytoev
   WRITE(stdout,'(5x,  "Density of down states at mu =   ",&
                      &f15.8,"  states/(eV spin cell)")')  ddosef / rytoev

   WRITE(stdout,'(/,5x,"Derivative of the up dos at mu =   ",&
                      &f13.8,"  states/(eV^2 spin cell)")')  uddosde/rytoev**2
   WRITE(stdout,'(5x,  "Derivative of the down dos at mu = ",&
                      &f13.8,"  states/(eV^2 spin cell)")')  dddosde/rytoev**2
END IF
WRITE(stdout,'(/,5x,"Density of states (g) at mu =",f19.8,&
                           &"  states/(eV cell)")')  dosef / rytoev
WRITE(stdout,'(/,5x,"Derivative of the dos (g'') at mu =",f14.8,&
                           &"  states/(eV^2 cell)")')  ddosde / rytoev**2
WRITE(stdout,'(/,5x,"g''/g at mu =", 21x, f15.8,&
                           &"  eV^(-1)")')  ddosde / dosef /rytoev

WRITE(stdout,'(/,5x,"Computing the electronic thermodynamic properties",/)')

el_mu(:,igeom)=0.0_DP
el_free_ener(:,igeom)=0.0_DP
el_ener(:,igeom)=0.0_DP
el_entr(:,igeom)=0.0_DP
el_ce(:,igeom)=0.0_DP

CALL divide (intra_image_comm, ntemp, nstart, nlast)
DO itemp = nstart, nlast
   IF (MOD(itemp,50)==0.OR.itemp==ntemp) &
                  WRITE(stdout,'(5x, "Computing temperature ", i8,&
                       &" /",i8, " Total", i8)') itemp, nlast-nstart+1, ntemp
   CALL el_chem_pot(eldos, temp(itemp), nelec, el_mu(itemp,igeom))
   CALL el_free_energy(eldos, temp(itemp), el_mu(itemp,igeom),  &
                                                el_free_ener(itemp,igeom))
   CALL el_energy(eldos, temp(itemp), el_mu(itemp,igeom), el_ener(itemp,igeom))
   CALL el_entropy(eldos, temp(itemp), el_mu(itemp,igeom), &
                                                     el_entr(itemp,igeom))
   CALL el_specific_heat_cv(eldos, temp(itemp), el_mu(itemp,igeom), &
                                                     el_ce(itemp,igeom))
END DO

CALL mp_sum(el_mu(:,igeom), intra_image_comm)
CALL mp_sum(el_free_ener(:,igeom), intra_image_comm)
CALL mp_sum(el_ener(:,igeom), intra_image_comm)
CALL mp_sum(el_entr(:,igeom), intra_image_comm)
CALL mp_sum(el_ce(:,igeom), intra_image_comm)
!
!   The zero of the energy is taken when all the valence bands are occupied
!
el_ener(:,igeom)=el_ener(:,igeom)-ene0
el_free_ener(:,igeom)=el_free_ener(:,igeom)-ene0

IF (ionode) THEN
   iu_therm=find_free_unit()
   fileeltherm='therm_files/'//TRIM(fleltherm)
   OPEN (UNIT=iu_therm, FILE=TRIM(fileeltherm), STATUS='unknown',&
                                                     FORM='formatted')
   WRITE(iu_therm,'("#")')  
   WRITE(iu_therm,'("# Temperature T in K, ")')
   WRITE(iu_therm,'("# Chemical potential in Ry")') 
   WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
   WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
   WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
   WRITE(iu_therm,'("# Multiply by ",f7.4," to have energies in &
                       &eV/cell etc..")') rytoev
   WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to have &
                  &energies in cal/(N mol).")') rytoev, electronvolt_si &
                         * avogadro / caltoj, rydberg_si*avogadro/caltoj
   WRITE(iu_therm,'("# Multiply by ",f7.4," x ",f8.2," = ",f9.1," to &
                  &have energies in J/(N mol).")') rytoev, electronvolt_si&
                         * avogadro, rydberg_si*avogadro
   WRITE(iu_therm,'("# N is the number of formula units per cell.")')
   WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                   &energies in cal/mol etc. ")')
   WRITE(iu_therm,'("#",5x,"   T  ", 10x, " energy ", 9x, "  free energy ",&
                  & 9x, " entropy ", 12x, " Cv ", 9x, "  chemical pot" )')

   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,5e20.12)') temp(itemp), &
             el_ener(itemp,igeom), el_free_ener(itemp,igeom), &
             el_entr(itemp,igeom), el_ce(itemp,igeom), el_mu(itemp,igeom)
   END DO

   CLOSE(iu_therm)
END IF

CALL destroy_eldos(eldos)

RETURN
END SUBROUTINE write_el_thermo
