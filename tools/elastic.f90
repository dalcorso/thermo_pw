!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM elastic
!
!  This program reads in input the elastic constants of a solid and
!  its Bravais lattice. Sets the elastic constants matrix and computes
!  a few auxiliary quantities such as the bulk modulus and some 
!  polycrystalline averages using the routines of the thermo_pw library
!
!  The input variables are:
!
!  ibrav : the Bravais lattice index
!
!  If necessary:
!  laue :  the Laue class
!  
!  c_ij : the elastic constants as requested and dependent on the laue
!         class.
!
!  optionally you can give also the density in kg/m^3. In this case the
!  code computes also the Voigt-Reuss-Hill averaged sound velocities.
!  If you do not know the density set it to 0.0.
!
!  Note 1 g/cm^3 -> 10^3 kg/m^3. If you have the density in g/cm^3 multiply
!  by 10^3 to have it in kg/m^3
!
!  If you give also the number of atoms per cell and the cell volume
!  (in (a.u.)^3) the code will compute the debye averaged sound velocity 
!  and the debye temperature. It produces also a file that contains the
!  Debye vibrational energy, free energy, entropy and heat capacity
!  in the same units used by thermo_pw. The temperature if from 1. to 900.
!  degrees in steps of 1. degree. These parameters cannot be changed from
!  input. You have to recompile the code to change them. These 
!  quantities are written in a file called thermo_debye.dat.
!  If you do not know set nat to 0.
!
USE kinds, ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end
USE elastic_constants, ONLY : print_elastic_constants,     &
                              compute_elastic_compliances, &
                              print_elastic_compliances,   &
                              print_macro_elasticity,      &
                              print_sound_velocities
USE debye_module, ONLY : compute_debye_temperature, debye_cv, &
                         debye_vib_energy, debye_free_energy, debye_entropy, &
                         compute_debye_temperature_poisson, debye_e0
USE io_global, ONLY : stdout, ionode

IMPLICIT NONE

INTEGER :: ibrav, laue, nat
REAL(DP) :: el_con(6,6)          ! the elastic constants
REAL(DP) :: el_compliances(6,6)  ! the elastic constants
REAL(DP) :: density, omega, debye_t, approx_debye_t, poisson, bulkm, deltat
REAL(DP) :: macro_el(8), vp, vb, vg, deb_e0
INTEGER :: i, ntemp, ios
REAL(DP), ALLOCATABLE :: temp(:), deb_cv(:), deb_energy(:), &
                         deb_free_energy(:), deb_entropy(:)
CHARACTER(LEN=9) :: code='elastic'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

WRITE(stdout,'(5x,"Bravais lattice index")') 
READ(5,*) ibrav

el_con=0.0_DP
IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
!
!  cubic
!
   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)
   el_con(2,2)=el_con(1,1)
   el_con(3,3)=el_con(1,1)
   el_con(5,5)=el_con(4,4)
   el_con(6,6)=el_con(4,4)
   el_con(1,3)=el_con(1,2)
   el_con(2,3)=el_con(1,2)
   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)
   el_con(5,5)=el_con(4,4)
   el_con(6,6)=el_con(4,4)
ELSEIF (ibrav==4) THEN
!
! hexagonal
!
   laue=23

   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(5,*) el_con(1,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(5,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)
   
   el_con(2,2)=el_con(1,1)
   el_con(5,5)=el_con(4,4)
   el_con(6,6)=0.5_DP * ( el_con(1,1) - el_con(1,2) )
   el_con(2,3)=el_con(1,3)
   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)
ELSEIF (ibrav==5) THEN
!
!  trigonal
!
   WRITE(stdout,'(5x,"laue class? (25 (D_3d) or 27 (S_6)) ")')
   READ(5,*) laue
   IF (laue /= 25 .AND. laue /= 27) CALL errore('elastic','Wrong Laue class',1)
   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(5,*) el_con(1,3)
   WRITE(stdout,'(5x,"C14?")')
   READ(5,*) el_con(1,4)
   WRITE(stdout,'(5x,"C33?")')
   READ(5,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)

   el_con(2,2)=el_con(1,1)
   el_con(5,5)=el_con(4,4)
   el_con(6,6)=0.5_DP * ( el_con(1,1) - el_con(1,2) )
   el_con(2,3)=el_con(1,3)
   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)

   el_con(2,4)=-el_con(1,4) 
   el_con(5,6)=el_con(1,4) 
   el_con(4,1)=el_con(1,4)
   el_con(4,2)=el_con(2,4)
   el_con(6,5)=el_con(5,6)

   IF (laue==27) THEN
      WRITE(stdout,'(5x,"C15?")')
      READ(5,*) el_con(1,5)
      el_con(2,5)=-el_con(1,5)
      el_con(4,6)=-el_con(1,5)
      el_con(5,1)=el_con(1,5)
      el_con(5,2)=el_con(2,5)
      el_con(6,4)=el_con(4,6)
   END IF
ELSEIF (ibrav==6 .OR. ibrav==7) THEN
!
! tetragonal
!
   WRITE(stdout,'(5x,"laue class? (18 (C_4h) or 22 (D_2d)) ")')
   READ(5,*) laue
   IF (laue /= 18 .AND. laue /= 22) CALL errore('elastic','Wrong Laue class',1)
   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(5,*) el_con(1,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(5,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)
   WRITE(stdout,'(5x,"C66?")')
   READ(5,*) el_con(6,6)

   el_con(2,2)=el_con(1,1)
   el_con(5,5)=el_con(4,4)
   el_con(2,3)=el_con(1,3)
   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)

   IF (laue==22) THEN
      WRITE(stdout,'(5x,"C16?")')
      READ(5,*) el_con(1,6)
      el_con(2,6)=-el_con(1,6)
      el_con(6,1)=el_con(1,6) 
      el_con(6,2)=el_con(2,6) 
   ENDIF
ELSEIF (ibrav==8 .OR. ibrav==9 .OR. ibrav==10 .OR. ibrav==11) THEN
!
!  Orthorombic, base-centered orthorhombic, face-centered orthorhombic,
!  body-centered orthorhombic
!
   laue=20
   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(5,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(5,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(5,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(5,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(5,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(5,*) el_con(6,6)

   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)
ELSEIF (ibrav==12 .OR. ibrav==13) THEN
!
! c-unique monoclinic or base-centered monoclinic
!
   laue=16
   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(5,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(5,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(5,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(5,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(5,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(5,*) el_con(6,6)
   WRITE(stdout,'(5x,"C16?")')
   READ(5,*) el_con(1,6)
   WRITE(stdout,'(5x,"C26?")')
   READ(5,*) el_con(2,6)
   WRITE(stdout,'(5x,"C36?")')
   READ(5,*) el_con(3,6)
   WRITE(stdout,'(5x,"C45?")')
   READ(5,*) el_con(4,5)

   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)

   el_con(6,1)=el_con(1,6)
   el_con(6,2)=el_con(2,6)
   el_con(6,3)=el_con(3,6)

   el_con(5,4)=el_con(4,5)

ELSEIF (ibrav==-12 .OR. ibrav==-13) THEN
!
!  b-unique monoclinic or base-centered monoclinic
!
   laue=16
   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(5,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(5,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(5,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(5,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(5,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(5,*) el_con(6,6)
   WRITE(stdout,'(5x,"C15?")')
   READ(5,*) el_con(1,5)
   WRITE(stdout,'(5x,"C25?")')
   READ(5,*) el_con(2,5)
   WRITE(stdout,'(5x,"C35?")')
   READ(5,*) el_con(3,5)
   WRITE(stdout,'(5x,"C46?")')
   READ(5,*) el_con(4,6)

   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)

   el_con(5,1)=el_con(1,5)
   el_con(5,2)=el_con(2,5)
   el_con(5,3)=el_con(3,5)

   el_con(6,4)=el_con(4,6)
ELSEIF (ibrav==14) THEN
!
!  triclinic
!
   laue=2
   WRITE(stdout,'(5x,"C11?")')
   READ(5,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(5,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(5,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(5,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(5,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(5,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(5,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(5,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(5,*) el_con(6,6)
   WRITE(stdout,'(5x,"C14?")')
   READ(5,*) el_con(1,4)
   WRITE(stdout,'(5x,"C15?")')
   READ(5,*) el_con(1,5)
   WRITE(stdout,'(5x,"C16?")')
   READ(5,*) el_con(1,6)
   WRITE(stdout,'(5x,"C24?")')
   READ(5,*) el_con(2,4)
   WRITE(stdout,'(5x,"C25?")')
   READ(5,*) el_con(2,5)
   WRITE(stdout,'(5x,"C26?")')
   READ(5,*) el_con(2,6)
   WRITE(stdout,'(5x,"C34?")')
   READ(5,*) el_con(3,4)
   WRITE(stdout,'(5x,"C35?")')
   READ(5,*) el_con(3,5)
   WRITE(stdout,'(5x,"C36?")')
   READ(5,*) el_con(3,6)
   WRITE(stdout,'(5x,"C45?")')
   READ(5,*) el_con(4,5)
   WRITE(stdout,'(5x,"C46?")')
   READ(5,*) el_con(4,6)
   WRITE(stdout,'(5x,"C56?")')
   READ(5,*) el_con(5,6)

   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(4,1)=el_con(1,4)
   el_con(5,1)=el_con(1,5)
   el_con(6,1)=el_con(1,6)

   el_con(3,2)=el_con(2,3)
   el_con(4,2)=el_con(2,4)
   el_con(5,2)=el_con(2,5)
   el_con(6,2)=el_con(2,6)

   el_con(4,3)=el_con(3,4)
   el_con(5,3)=el_con(3,5)
   el_con(6,3)=el_con(3,6)

   el_con(5,4)=el_con(4,5)
   el_con(6,4)=el_con(4,6)

   el_con(6,5)=el_con(5,6)

ELSE
   CALL errore('elastic','Bravais lattice not programmed',1)
ENDIF

WRITE(stdout,'(5x,"Density (kg/m^3)?")')
READ(5,*) density
IF (density > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Number of atoms per cell?")')
   READ(5,*) nat
   IF (nat>0) THEN
      WRITE(stdout,'(5x,"Volume per cell?")')
      READ(5,*) omega
   END IF
END IF

CALL print_elastic_constants(el_con, .FALSE.)
!
!  now compute the elastic compliances and prints them
!
CALL compute_elastic_compliances(el_con,el_compliances)
CALL print_elastic_compliances(el_compliances, .FALSE.)
!
!  now compute the macro elasticity quantities
!
CALL print_macro_elasticity(ibrav, el_con, el_compliances, macro_el, .TRUE.)

IF (density > 0.0_DP) THEN
!
!  print the Voigt-Reuss-Hill average sound velocities
!
   CALL print_sound_velocities( ibrav, el_con, el_compliances, &
                                               density, vp, vb, vg )
!
!  compute the Debye temperature
!
   IF (nat > 0 .AND. omega > 0.0_DP) THEN
!
!   first the approximate debye temperature
!
      poisson=(macro_el(4) + macro_el(8)) * 0.5_DP
      bulkm=(macro_el(1) + macro_el(5)) * 0.5_DP
      CALL compute_debye_temperature_poisson(poisson, bulkm, &
                               density, nat, omega, approx_debye_t)
!
!  then the real one
!
      CALL compute_debye_temperature(el_con, density, nat, omega, debye_t)
      ntemp=900
      deltat=1._DP
      ALLOCATE(temp(ntemp))
      ALLOCATE(deb_cv(ntemp))
      ALLOCATE(deb_energy(ntemp))
      ALLOCATE(deb_entropy(ntemp))
      ALLOCATE(deb_free_energy(ntemp))
      DO i=1,ntemp
         temp(i)=deltat*i
      ENDDO
      CALL debye_e0 (debye_t, nat, deb_e0)
      CALL debye_cv (debye_t, temp, ntemp, nat, deb_cv)
      CALL debye_vib_energy (debye_t, temp, ntemp, nat, deb_energy)
      CALL debye_free_energy (debye_t, temp, ntemp, nat, deb_free_energy)
      CALL debye_entropy(debye_t, temp, ntemp, nat, deb_entropy)

      IF (ionode) THEN
         OPEN(UNIT=25, FILE='thermo_debye.dat', STATUS='unknown', &
                       FORM='formatted', IOSTAT=ios)
         WRITE(25,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  &energies in cal/(N mol).")')
         WRITE(25,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  &have energies in J/(N mol).")')
         WRITE(25,'("# N is the number of formula units per cell.")')
         WRITE(25,'("#",5x,"   T  ", 10x, " energy ", 9x, "  free energy ",&
                  & 9x, " entropy ", 12x, " Cv ")')

         DO i=1,ntemp
            WRITE(25,'(5E16.8)') temp(i), deb_energy(i)+deb_e0,  &
                              deb_free_energy(i)+deb_e0, deb_entropy(i), &
                              deb_cv(i)
         ENDDO
         CLOSE(25)
      END IF
      DEALLOCATE(temp)
      DEALLOCATE(deb_cv)
      DEALLOCATE(deb_energy)
      DEALLOCATE(deb_entropy)
      DEALLOCATE(deb_free_energy)
   END IF
END IF

CALL environment_end( code )
!
CALL mp_global_end ()

END PROGRAM elastic
