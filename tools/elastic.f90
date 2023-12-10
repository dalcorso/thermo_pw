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
USE constants, ONLY : pi
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end
USE elastic_constants, ONLY : print_elastic_constants,     &
                              compute_elastic_compliances, &
                              print_elastic_compliances,   &
                              print_macro_elasticity,      &
                              print_sound_velocities,      &
                              compute_sound

USE voigt,             ONLY : to_voigt4 
USE gnuplot,           ONLY : gnuplot_start, gnuplot_end,  &
                              gnuplot_write_header,        &
                              gnuplot_ylabel,              &
                              gnuplot_xlabel,              &
                              gnuplot_write_command,       &
                              gnuplot_write_file_mul_data

USE debye_module, ONLY : compute_debye_temperature, debye_cv, &
                         debye_vib_energy, debye_free_energy, debye_entropy, &
                         compute_debye_temperature_poisson, debye_e0, &
                         compute_debye_temperature_macro_el
USE io_global, ONLY : stdout, ionode

IMPLICIT NONE

INTEGER :: ibrav, laue, nat
REAL(DP) :: el_con(6,6)          ! the elastic constants
REAL(DP) :: el_compliances(6,6)  ! the elastic constants
REAL(DP) :: density, omega, debye_t, approx_debye_t, poisson, bulkm, deltat
REAL(DP) :: macro_el(8), vp, vb, vg, deb_e0, debye_macro_el
INTEGER :: i, ntemp, iundeb, ios, iunsound
INTEGER :: find_free_unit, stdin
REAL(DP), ALLOCATABLE :: temp(:), deb_cv(:), deb_energy(:), &
                         deb_free_energy(:), deb_entropy(:)
CHARACTER(LEN=9) :: code='elastic'

LOGICAL :: flag
INTEGER :: l
INTEGER :: ichoice
INTEGER, PARAMETER :: lmax=101
REAL(DP) :: q1, q2, sound_speeds(3,lmax), sound_disp(3,3), theta, factor, &
            qvec(3), el_con3(3,3,3,3)

CHARACTER(LEN=256) :: filename, filename1, gnu_filename
INTEGER :: ierr, system
CHARACTER(LEN=8) :: flext
CHARACTER(LEN=80) :: gnuplot_command

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

stdin=5
WRITE(stdout,'(5x,"Bravais lattice index")') 
READ(stdin,*) ibrav
WRITE(stdout,'(5x,"Input elastic constants in kbar")')

el_con=0.0_DP
IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
!
!  cubic
!
   WRITE(stdout,'(5x,"C11?")')
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)
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
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(stdin,*) el_con(1,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(stdin,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)
   
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
   READ(stdin,*) laue
   IF (laue /= 25 .AND. laue /= 27) CALL errore('elastic','Wrong Laue class',1)
   WRITE(stdout,'(5x,"C11?")')
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(stdin,*) el_con(1,3)
   WRITE(stdout,'(5x,"C14?")')
   READ(stdin,*) el_con(1,4)
   WRITE(stdout,'(5x,"C33?")')
   READ(stdin,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)

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
      READ(stdin,*) el_con(1,5)
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
   READ(stdin,*) laue
   IF (laue /= 18 .AND. laue /= 22) CALL errore('elastic','Wrong Laue class',1)
   WRITE(stdout,'(5x,"C11?")')
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(stdin,*) el_con(1,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(stdin,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)
   WRITE(stdout,'(5x,"C66?")')
   READ(stdin,*) el_con(6,6)

   el_con(2,2)=el_con(1,1)
   el_con(5,5)=el_con(4,4)
   el_con(2,3)=el_con(1,3)
   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)

   IF (laue==22) THEN
      WRITE(stdout,'(5x,"C16?")')
      READ(stdin,*) el_con(1,6)
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
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(stdin,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(stdin,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(stdin,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(stdin,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(stdin,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(stdin,*) el_con(6,6)

   el_con(2,1)=el_con(1,2)
   el_con(3,1)=el_con(1,3)
   el_con(3,2)=el_con(2,3)
ELSEIF (ibrav==12 .OR. ibrav==13) THEN
!
! c-unique monoclinic or base-centered monoclinic
!
   laue=16
   WRITE(stdout,'(5x,"C11?")')
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(stdin,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(stdin,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(stdin,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(stdin,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(stdin,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(stdin,*) el_con(6,6)
   WRITE(stdout,'(5x,"C16?")')
   READ(stdin,*) el_con(1,6)
   WRITE(stdout,'(5x,"C26?")')
   READ(stdin,*) el_con(2,6)
   WRITE(stdout,'(5x,"C36?")')
   READ(stdin,*) el_con(3,6)
   WRITE(stdout,'(5x,"C45?")')
   READ(stdin,*) el_con(4,5)

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
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(stdin,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(stdin,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(stdin,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(stdin,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(stdin,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(stdin,*) el_con(6,6)
   WRITE(stdout,'(5x,"C15?")')
   READ(stdin,*) el_con(1,5)
   WRITE(stdout,'(5x,"C25?")')
   READ(stdin,*) el_con(2,5)
   WRITE(stdout,'(5x,"C35?")')
   READ(stdin,*) el_con(3,5)
   WRITE(stdout,'(5x,"C46?")')
   READ(stdin,*) el_con(4,6)

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
   READ(stdin,*) el_con(1,1)
   WRITE(stdout,'(5x,"C12?")')
   READ(stdin,*) el_con(1,2)
   WRITE(stdout,'(5x,"C13?")')
   READ(stdin,*) el_con(1,3)
   WRITE(stdout,'(5x,"C22?")')
   READ(stdin,*) el_con(2,2)
   WRITE(stdout,'(5x,"C23?")')
   READ(stdin,*) el_con(2,3)
   WRITE(stdout,'(5x,"C33?")')
   READ(stdin,*) el_con(3,3)
   WRITE(stdout,'(5x,"C44?")')
   READ(stdin,*) el_con(4,4)
   WRITE(stdout,'(5x,"C55?")')
   READ(stdin,*) el_con(5,5)
   WRITE(stdout,'(5x,"C66?")')
   READ(stdin,*) el_con(6,6)
   WRITE(stdout,'(5x,"C14?")')
   READ(stdin,*) el_con(1,4)
   WRITE(stdout,'(5x,"C15?")')
   READ(stdin,*) el_con(1,5)
   WRITE(stdout,'(5x,"C16?")')
   READ(stdin,*) el_con(1,6)
   WRITE(stdout,'(5x,"C24?")')
   READ(stdin,*) el_con(2,4)
   WRITE(stdout,'(5x,"C25?")')
   READ(stdin,*) el_con(2,5)
   WRITE(stdout,'(5x,"C26?")')
   READ(stdin,*) el_con(2,6)
   WRITE(stdout,'(5x,"C34?")')
   READ(stdin,*) el_con(3,4)
   WRITE(stdout,'(5x,"C35?")')
   READ(stdin,*) el_con(3,5)
   WRITE(stdout,'(5x,"C36?")')
   READ(stdin,*) el_con(3,6)
   WRITE(stdout,'(5x,"C45?")')
   READ(stdin,*) el_con(4,5)
   WRITE(stdout,'(5x,"C46?")')
   READ(stdin,*) el_con(4,6)
   WRITE(stdout,'(5x,"C56?")')
   READ(stdin,*) el_con(5,6)

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
READ(stdin,*) density
IF (density > 0.0_DP) THEN
   WRITE(stdout,'(5x,"Number of atoms per cell?")')
   READ(stdin,*) nat
   IF (nat>0) THEN
      WRITE(stdout,'(5x,"Volume per cell (a.u.)^3?")')
      READ(stdin,*) omega
   END IF

   WRITE(stdout,'(5x,"Plane to plot? (1: yz, 2: xz, 3: xy)")')
   READ(stdin,*) ichoice
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
                                               density, vp, vb, vg, .TRUE. )
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
!   the debye temperature with macro elasticity
!
      CALL compute_debye_temperature_macro_el(vp, vg, density, nat, &
                                                      omega, debye_macro_el)

      WRITE(stdout,'(/,5x,"Debye temperature from speed of sound",f15.7,&
                           " K")') debye_macro_el
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
         iundeb=find_free_unit()
         OPEN(UNIT=iundeb, FILE='thermo_debye.dat', STATUS='unknown', &
                       FORM='formatted', IOSTAT=ios)
         WRITE(iundeb,'("# Multiply by 13.6058 x 23060.35 = 313 754.5 to have &
                  &energies in cal/(N mol).")')
         WRITE(iundeb,'("# Multiply by 13.6058 x 96526.0 = 1 313 313 to &
                  &have energies in J/(N mol).")')
         WRITE(iundeb,'("# N is the number of formula units per cell.")')
         WRITE(iundeb,'("#",5x,"   T  ", 10x, " energy ", 9x, "  free energy ",&
                  & 9x, " entropy ", 12x, " Cv ")')

         DO i=1,ntemp
            WRITE(iundeb,'(5E16.8)') temp(i), deb_energy(i)+deb_e0,  &
                              deb_free_energy(i)+deb_e0, deb_entropy(i), &
                              deb_cv(i)
         ENDDO
         CLOSE(iundeb)
      END IF
      DEALLOCATE(temp)
      DEALLOCATE(deb_cv)
      DEALLOCATE(deb_energy)
      DEALLOCATE(deb_entropy)
      DEALLOCATE(deb_free_energy)
   END IF
!
!  Contributed by X. Gong
!
   flag = .false. ! from 6*6 to 3*3*3*3
   call to_voigt4(el_con, el_con3, flag)

   DO l=1,lmax,1
      theta = 2.0_DP*pi*(l-1)/(lmax-1)
      q1 = cos(theta)
      q2 = sin(theta)
      ! z1z2 surface q = (q1,q2,0)
      IF (ichoice == 3) THEN 
         qvec=(/q1,q2,0.0_DP/) 
      ELSEIF (ichoice == 2) THEN 
      ! z2z3 surface q = (0,q2,q3)
         qvec=(/q1,0.0_DP,q2/)
      ELSEIF (ichoice == 1) THEN 
         qvec=(/0.0_DP,q1,q2/)
      ELSE 
         CALL errore("elastic", "ichoice must be 1, 2 or 3", 1)
      ENDIF
      CALL compute_sound(el_con3, qvec, density, sound_speeds(1,l), sound_disp)
   ENDDO

   filename1="velocity_surfaces.dat"

   IF (ionode) THEN
      iunsound=find_free_unit()
      OPEN(UNIT=iunsound, FILE=filename1, STATUS='unknown', &
                          FORM='formatted', IOSTAT=ios)
      factor = 1.D-3
      DO l=1,lmax,1
         theta = 2.0_DP*pi*(l-1)/(lmax-1)
         q1 = cos(theta)
         q2 = sin(theta)
         WRITE(iunsound,*) theta, sound_speeds(1,l)*factor,        &
         sound_speeds(2,l)*factor, sound_speeds(3,l)*factor,       &
         sound_speeds(1,l)*factor*q1, sound_speeds(1,l)*factor*q2, &
         sound_speeds(2,l)*factor*q1, sound_speeds(2,l)*factor*q2, &
         sound_speeds(3,l)*factor*q1, sound_speeds(3,l)*factor*q2
      ENDDO
      CLOSE(UNIT=iunsound, STATUS='KEEP')
   ENDIF

   gnu_filename="plot_velocity_surfaces.gnu"
   filename="velocity_surfaces.ps"
   CALL gnuplot_start(gnu_filename)

   flext = ".ps"
   CALL gnuplot_write_header(filename, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                            1.0_DP, flext ) 
   CALL gnuplot_write_command("set size ratio     -1",.FALSE.)
   CALL gnuplot_write_command("set noborder",.FALSE.)
   CALL gnuplot_write_command("set xzeroaxis",.FALSE.)
   CALL gnuplot_write_command("set yzeroaxis",.FALSE.)
   CALL gnuplot_write_command("set xtics axis",.FALSE.)     
   CALL gnuplot_write_command("set ytics axis",.FALSE.)
   CALL gnuplot_write_command("set xtics add ("""" 0)",.FALSE.)  
   CALL gnuplot_write_command("set ytics add ("""" 0)",.FALSE.)
   CALL gnuplot_write_command("set arrow 1 from 0,0 to graph 1, &
                                             &first 0 filled head",.FALSE.)
   CALL gnuplot_write_command("set arrow 2 from 0,0 to first 0, &
                                             &graph 1 filled head",.FALSE.)
   CALL gnuplot_write_command("set arrow 3 from 0,0 to graph 0, &
                                             &first 0 filled nohead",.FALSE.)
   CALL gnuplot_write_command("set arrow 4 from 0,0 to first 0, &
                                             &graph 0 filled nohead",.FALSE.)
   CALL gnuplot_write_file_mul_data(filename1,5,6,'color_red', & 
                                              .TRUE.,.FALSE., .FALSE.)
   CALL gnuplot_write_file_mul_data(filename1,7,8,'color_green',&
                                              .FALSE.,.FALSE., .FALSE.)
   CALL gnuplot_write_file_mul_data(filename1,9,10,'color_blue',&
                                              .FALSE.,.TRUE., .FALSE.)
   CALL gnuplot_end()

   gnuplot_command="gnuplot"

   IF (ionode) &
      ierr=system(TRIM(gnuplot_command)//' '//TRIM(gnu_filename))

   !IF (lgnuplot.AND.ionode) &
   !   CALL EXECUTE_COMMAND_LINE(TRIM(gnuplot_command)//' '&
   !                                       //TRIM(gnu_filename), WAIT=.FALSE.)
END IF

CALL environment_end( code )
!
CALL mp_global_end ()

END PROGRAM elastic
