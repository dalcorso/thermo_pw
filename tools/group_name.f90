PROGRAM group_name

USE space_groups, ONLY : find_space_group_names, set_add_info
IMPLICIT NONE
INTEGER :: g_number

WRITE(6,'(5x,"Space group number?")') 
READ(5,*) g_number

CALL set_add_info()

IF (g_number > 0.AND. g_number <= 230 ) THEN
   WRITE(6,'(5x,"The name(s) of the space group ",i5, " is (are):")') g_number
   CALL find_space_group_names(g_number)
   WRITE(6,'(/,5x,"The first name is in the ITA tables. ")')
   WRITE(6,'(5x,"s short name, S Schoenflies symbol, * name used in 1935 edition of the IT.")')
   IF (g_number > 15 .AND. g_number < 75) THEN
      WRITE(6,'(5x,"Orthorombic group.")')
      WRITE(6,'(5x,"Apply the second transformation to the coordinates and")')
      WRITE(6,'(5x,"to celldm to pass to the ITA group, and the first  &
                    &transformation")')
      WRITE(6,'(5x,"to return to the given group.")')
   ENDIF
ELSE
   WRITE(6,'(5x,"Group number must be between 1 and 230" )') 
ENDIF

END PROGRAM group_name
