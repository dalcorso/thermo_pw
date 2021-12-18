!
! Copyright (C) 2021 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM change_dynmat_name
!
! This program changes the geometry number of all the dynamical matrices
! files indicated in input.
! It reads a file with the format:
! fildyn        ! name of the dynamical matrix file
! n             ! number of geometries to change
! orig1  new1   ! original and new number of the first change
! orig2  new2 
! ...    ...
! orign  newn   ! original and new number of the n change
!
! Note that all changes are simultaneous, so orig# must be all different
! and the order of the changes is not relevant.
! If new# is zero the dynamical matrix orig# is removed.
!
USE io_global, ONLY : stdout
IMPLICIT NONE
CHARACTER(LEN=256) :: filename, fildyn_in, fildyn_out
INTEGER :: n
INTEGER, ALLOCATABLE :: orig(:), new(:), ndyn(:)
CHARACTER(LEN=256), ALLOCATABLE :: save_name(:)
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=4) :: ext
LOGICAL :: xmldyn, has_xml
INTEGER :: i, j, leng, len1, nfiles, nfiles_tot, ios
INTEGER :: ierr, system

READ(5,*) filename

leng=LEN(filename)
xmldyn=has_xml(filename)
ext=''
IF (xmldyn) THEN
   filename=filename(1:leng-4)
   ext='.xml'
ENDIF
READ(5,*) n

ALLOCATE(orig(n))
ALLOCATE(new(n))
ALLOCATE(ndyn(n))

DO i=1, n
   READ(5,*) orig(i), new(i)
ENDDO

DO i=1,n
   fildyn_in=TRIM(filename)//'.g'//TRIM(int_to_char(orig(i)))//'.0'
   WRITE(stdout, '("Opening ",a)') TRIM(fildyn_in)
   OPEN(UNIT=28, FILE=TRIM(fildyn_in), STATUS='old', ERR=100, IOSTAT=ios)
   READ(28,*)
   READ(28,*) ndyn(i)
   CLOSE(28)
ENDDO
100 CALL errore('change_name','opening fildyn',ABS(ios))

nfiles=0
DO i=1, n
   nfiles=nfiles + ndyn(i) + 1
ENDDO

nfiles_tot=nfiles
ALLOCATE(save_name(nfiles_tot))

nfiles=0
DO i=1,n
   fildyn_in=TRIM(filename)//'.g'//TRIM(int_to_char(orig(i)))//'.0'
   IF (new(i)==0) THEN
      ierr=system('rm '//TRIM(fildyn_in))
   ELSE   
      fildyn_out=TRIM(filename)//'.g'//TRIM(int_to_char(new(i)))//'.0_b'
      ierr=system('mv '//TRIM(fildyn_in)//' '//TRIM(fildyn_out))
      nfiles=nfiles+1
      save_name(nfiles)=fildyn_out
   ENDIF
   DO j=1,ndyn(i)
      fildyn_in=TRIM(filename)//'.g'//TRIM(int_to_char(orig(i)))//'.'//&
                                &TRIM(int_to_char(j))//TRIM(ext)
      IF (new(i)==0) THEN
         ierr=system('rm '//TRIM(fildyn_in))
      ELSE
         fildyn_out=TRIM(filename)//'.g'//TRIM(int_to_char(new(i)))//'.'&
                                 &//TRIM(int_to_char(j))//TRIM(ext)//'_b'
         ierr=system('mv '//TRIM(fildyn_in)//' '//TRIM(fildyn_out))
         nfiles=nfiles+1
         IF (nfiles>nfiles_tot) CALL errore('change_name',&
                                                'wrong number of files',1)
         save_name(nfiles)=fildyn_out
      ENDIF
   ENDDO
ENDDO

DO i=1,nfiles
   fildyn_in=save_name(i)
   len1=LEN(TRIM(fildyn_in))
   fildyn_out=TRIM(fildyn_in(1:len1-2))
   ierr=system('mv '//TRIM(fildyn_in)//' '//TRIM(fildyn_out))
ENDDO

DEALLOCATE(save_name)
DEALLOCATE(ndyn)
DEALLOCATE(new)
DEALLOCATE(orig)

END PROGRAM change_dynmat_name
