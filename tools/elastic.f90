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
!  If necessary:
!  laue :  the Laue class
!  
!  c_ij : the elastic constants as requested and dependent on the laue
!         class.
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
!  Orthorombic, base-centered orthorombic, face-centered orthorombic,
!  body-centered orthorombic
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
