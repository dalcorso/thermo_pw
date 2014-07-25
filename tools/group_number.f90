!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM group_number

USE space_groups, ONLY : find_space_group_number, find_space_group_names, &
                         set_add_info
IMPLICIT NONE
INTEGER :: g_number
CHARACTER(LEN=11) :: group_name

WRITE(6,'(5x,"Space group name?")') 
READ(5,'(a)') group_name

CALL set_add_info()
CALL find_space_group_number(group_name,g_number)

IF (g_number > 0 ) THEN
   WRITE(6,'(5x,"The space group number of ",a," is ",i5)') TRIM(group_name), &
                                                            g_number
   WRITE(6,'(5x,"Other names of the same group:",/)')
   CALL find_space_group_names(g_number)
   WRITE(6,'(/,5x,"The first name is in the ITA tables. ")')
   WRITE(6,'(5x,"s short name, S Schoenflies symbol, * name used in 1935 edition of IT.")')
   IF (g_number > 15 .AND. g_number < 75) THEN
      WRITE(6,'(5x,"Orthorombic group.")')
      WRITE(6,'(5x,"Apply the second transformation to the coordinates and")')
      WRITE(6,'(5x,"to celldm to pass to the ITA group, and the first &
                    &transformation")')
      WRITE(6,'(5x,"to return to the given group.")')
   ENDIF
ELSE
   WRITE(6,'(5x,"This group name is unknown" )') 
ENDIF


END PROGRAM group_number
