!
! Copyright (C) 2021 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM change_e_name
!
! This program changes the geometry number of all the dynamical matrices
! files indicated in input.
! It reads a file with the format:
! efile         ! name of the files to change
! n             ! number of geometries to change
! orig1  new1   ! original and new number of the first change
! orig2  new2 
! ...    ...
! orign  newn   ! original and new number of the n change
!
! Note that all changes are simultaneous, so orig# must be all different
! and the order of the changes is not relevant.
! If new# is zero the file orig# is removed.
!
IMPLICIT NONE
CHARACTER(LEN=256) :: filename, file_in, file_out
INTEGER :: n
INTEGER, ALLOCATABLE :: orig(:), new(:)
CHARACTER(LEN=256), ALLOCATABLE :: save_name(:)
CHARACTER(LEN=6) :: int_to_char
LOGICAL :: xmldyn, has_xml
INTEGER :: i, len1, nfiles, ios
INTEGER :: ierr, system

READ(5,*) filename
READ(5,*) n

ALLOCATE(orig(n))
ALLOCATE(new(n))
ALLOCATE(save_name(n))

DO i=1, n
   READ(5,*) orig(i), new(i)
ENDDO

nfiles=0
DO i=1,n
   file_in=TRIM(filename)//'.g'//TRIM(int_to_char(orig(i)))
   IF (new(i)==0) THEN
      ierr=system('rm '//TRIM(file_in))
   ELSE
      file_out=TRIM(filename)//'.g'//TRIM(int_to_char(new(i)))//'.b'
      ierr=system('mv '//TRIM(file_in)//' '//TRIM(file_out))
      nfiles=nfiles+1
      save_name(nfiles)=file_out
   ENDIF
ENDDO

DO i=1,nfiles
   file_in=save_name(i)
   len1=LEN(TRIM(file_in))
   file_out=TRIM(file_in(1:len1-2))
   ierr=system('mv '//TRIM(file_in)//' '//TRIM(file_out))
ENDDO

DEALLOCATE(save_name)
DEALLOCATE(new)
DEALLOCATE(orig)

END PROGRAM change_e_name
