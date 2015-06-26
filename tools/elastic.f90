!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM elastic
!
!  This small program reads in input the elastic constants of a solid and
!  its Bravais lattice. Sets the elastic constants matrix and computes
!  a few auxiliary quantities such as the bulk modulus 
!  and some microcrystalline averages using the routines of the
!  thermo_pw library
!
!  The input variables are:
!
!  ibrav : the Bravais lattice index
!
!  If necessary index:
!  laue :  the Laue class
!  
!  c_ij : the elastic constants as requested and dependent on the laue
!         class.
!
!  Limitation: still limited to cubic and hexagonal crystals.
!
!
USE kinds, ONLY : DP
USE elastic_constants, ONLY : print_elastic_constants,     &
                              compute_elastic_compliances, &
                              print_elastic_compliances,   &
                              print_macro_elasticity
USE io_global, ONLY : stdout

IMPLICIT NONE

INTEGER :: ibrav, laue
REAL(DP) :: el_con(6,6)          ! the elastic constants
REAL(DP) :: el_compliances(6,6)  ! the elastic constants

WRITE(stdout,'(5x,"Bravais lattice index")') 
READ(5,*) ibrav

el_con=0.0_DP
IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
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
ELSE
   CALL errore('elastic','Bravais lattice not programmed',1)
ENDIF

CALL print_elastic_constants(el_con, .FALSE.)
!
!  now compute the elastic compliances and prints them
!
CALL compute_elastic_compliances(el_con,el_compliances)
CALL print_elastic_compliances(el_compliances, .FALSE.)
!
!  now compute the macro elasticity quantities
!
CALL print_macro_elasticity(ibrav, el_con, el_compliances)

END PROGRAM elastic
