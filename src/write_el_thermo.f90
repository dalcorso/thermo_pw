!
! Copyright (C) 2016-present Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_el_thermo()
!
!  This routine writes on file the electronic thermodynamical quantities
!
USE kinds,          ONLY : DP
USE constants,      ONLY : rytoev
USE eldos_module,   ONLY : eldos_type, el_free_energy, el_energy, el_entropy, &
                           el_specific_heat_cv, destroy_eldos, &
                           el_chem_pot, read_eldos_data
USE ener,           ONLY : ef
USE lsda_mod,       ONLY : nspin
USE klist,          ONLY : degauss, nelec
USE ktetra,         ONLY : ltetra
USE lsda_mod,       ONLY : lsda
USE temperature,    ONLY : tmin, tmax, deltat, ntemp, temp
USE mp_images,      ONLY : root_image, my_image_id, intra_image_comm
USE mp,             ONLY : mp_bcast
USE io_global,      ONLY : ionode, ionode_id, stdout
USE data_files,     ONLY : fleltherm, fleldos

IMPLICIT NONE

INTEGER          :: itemp
INTEGER          :: iu_therm
TYPE(eldos_type) :: eldos

CHARACTER(LEN=256) :: fileeltherm, fileeldos
REAL(DP), ALLOCATABLE :: el_mu(:)
REAL(DP), ALLOCATABLE :: el_free_ener(:)
REAL(DP), ALLOCATABLE :: el_ener(:)
REAL(DP), ALLOCATABLE :: el_entr(:)
REAL(DP), ALLOCATABLE :: el_cv(:)
REAL(DP) :: ene0, mu0
REAL(DP) ::  dos1, dos2, dosef, udosef, ddosef, ddos1, ddos2, &
             uddosde, dddosde, ddosde
INTEGER :: n1, n2, n, ndos
!
IF (my_image_id /= root_image) RETURN
!
!  Electrons contribute to the thermodynamical properties of metals only.
!  For insulators return.
!
IF (degauss==0.0_DP.AND..NOT.ltetra) RETURN

WRITE(stdout,'(/,2x,76("+"))')
WRITE(stdout,'(5x,"Computing the thermodynamic properties from electron dos")')
WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(fleltherm)
WRITE(stdout,'(2x,76("+"),/)')

fileeldos='therm_files/'//TRIM(fleldos)
CALL read_eldos_data(eldos, lsda, fileeldos)

ALLOCATE(el_mu(ntemp))
ALLOCATE(el_free_ener(ntemp))
ALLOCATE(el_ener(ntemp))
ALLOCATE(el_entr(ntemp))
ALLOCATE(el_cv(ntemp))
!
!  The energy at very low temperature is consider as the zero of the energy
!  The default low temperature is 4. K or the lowest temperature
!  required in input.
!
CALL el_chem_pot(eldos, min(4.0_DP,temp(1)), nelec, mu0)
CALL el_energy(eldos, min(4.0_DP,temp(1)), mu0, ene0)
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

DO itemp = 1, ntemp
   CALL el_chem_pot(eldos, temp(itemp), nelec, el_mu(itemp))
   CALL el_free_energy(eldos, temp(itemp), el_mu(itemp), el_free_ener(itemp))
   CALL el_energy(eldos, temp(itemp), el_mu(itemp), el_ener(itemp))
   CALL el_entropy(eldos, temp(itemp), el_mu(itemp), el_entr(itemp))
   CALL el_specific_heat_cv(eldos, temp(itemp), el_mu(itemp), el_cv(itemp))
END DO

IF (ionode) THEN
   iu_therm=2
   fileeltherm='therm_files/'//TRIM(fleltherm)
   OPEN (UNIT=iu_therm, FILE=TRIM(fileeltherm), STATUS='unknown',&
                                                     FORM='formatted')
   WRITE(iu_therm,'("#")')  
   WRITE(iu_therm,'("# Temperature T in K, ")')
   WRITE(iu_therm,'("# Chemical potential in Ry")') 
   WRITE(iu_therm,'("# Energy and free energy in Ry/cell,")')
   WRITE(iu_therm,'("# Entropy in Ry/cell/K,")')
   WRITE(iu_therm,'("# Heat capacity Cv in Ry/cell/K.")')
   WRITE(iu_therm,'("# Multiply by 13.6058 to have energies in &
                       &eV/cell etc..")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  &energies in cal/(N mol).")')
   WRITE(iu_therm,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  &have energies in J/(N mol).")')
   WRITE(iu_therm,'("# N is the number of formula units per cell.")')
   WRITE(iu_therm,'("# For instance in silicon N=2. Divide by N to have &
                   &energies in cal/mol etc. ")')
   WRITE(iu_therm,'("#",5x,"   T  ", 10x, " energy ", 9x, "  free energy ",&
                  & 9x, " entropy ", 12x, " Cv ", 9x, "  chemical pot" )')

   DO itemp = 1, ntemp
      WRITE(iu_therm, '(e16.8,5e20.12)') temp(itemp), &
                    el_ener(itemp)-ene0, el_free_ener(itemp)-ene0, &
                    el_entr(itemp), el_cv(itemp), el_mu(itemp)
   END DO

   CLOSE(iu_therm)
END IF

CALL destroy_eldos(eldos)
DEALLOCATE( el_mu )
DEALLOCATE( el_free_ener )
DEALLOCATE( el_ener )
DEALLOCATE( el_entr )
DEALLOCATE( el_cv )

RETURN
END SUBROUTINE write_el_thermo
