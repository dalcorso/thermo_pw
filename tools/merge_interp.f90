!
! Copyright (C) 2021 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
PROGRAM merge_interp
!--------------------------------------------------------------------
!
!  This program is a driver for the routine dosplineint_1D.
!  The input is the following
!  file1    the name of file 1 without any ' ' or " "
!  ncol     number of columns of file 1
!  nx1, ny1 the column with x and f(x) for file1
!  fact1    the column nx1 and ny1 are multiplied by fact1
!  file2    the name of file 2
!  ncol     number of columns of file 2
!  nx2, ny2 the column with x and f(x) for file2
!  fact2    the column nx2 and ny2 are multiplied by fact2
!  fileout  the name of the output file.
!
!  on output it provides a file with three columns
!  nx2*fact2 ny1*fact1,  ny2*fact2 the column nx2, the column ny1 interpolated 
!                      on the mesh nx2 and a copy of the column ny2
!                      fact1 and fact2 can be used if the quantities refer
!                      to different number of atoms. Otherwise set them to 1.0
!
!  Note that the spline is not good to extrapolate, so 
!  only the points of nx2 in the range of column nx1 are written on output 
!
USE kinds, ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end

USE splinelib,        ONLY : dosplineint

USE io_global,        ONLY : stdout, meta_ionode

IMPLICIT NONE
CHARACTER(LEN=9) :: code='MERGE_INT'

CHARACTER(LEN=256) :: line_read, file1, file2, fileout
REAL(DP), ALLOCATABLE :: field1(:,:), field2(:,:), new_field(:)
REAL(DP) :: xmin, xmax, fact1, fact2
INTEGER, PARAMETER :: maxrow=100000
INTEGER :: ios, nline1, nline2, i
INTEGER :: ncol1, ncol2, nx1, nx2, ny1, ny2, iu_file

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

IF (meta_ionode) THEN
   READ(5,'(a)') file1
   READ(5,*) ncol1
   READ(5,*) nx1, ny1
   READ(5,*) fact1
   READ(5,'(a)') file2
   READ(5,*) ncol2
   READ(5,*) nx2, ny2
   READ(5,*) fact2
   READ(5,*) fileout

   WRITE(stdout,'(/,5x,"Reading from ",a," and ",a)') TRIM(file1), TRIM(file2)
   WRITE(stdout,'(5x,"Writing on   ",a)') TRIM(fileout)
   WRITE(stdout,'(/,5x,"File1: columns ",i5,",",i5," of ",i5)') nx1, ny1, ncol1
   WRITE(stdout,'(5x,"Multiplication factor:",f15.5)') fact1
   WRITE(stdout,'(5x,"File2: columns ",i5,",",i5," of ",i5)') nx2, ny2, ncol2
   WRITE(stdout,'(5x,"Multiplication factor:",f15.5)') fact2

   IF (nx1<1.OR.nx1>ncol1) CALL errore('merge_interp','error in nx1',1)
   IF (ny1<1.OR.ny1>ncol1) CALL errore('merge_interp','error in ny1',1)
   IF (nx2<1.OR.nx2>ncol2) CALL errore('merge_interp','error in nx2',1)
   IF (ny1<1.OR.ny2>ncol2) CALL errore('merge_interp','error in ny2',1)

   ALLOCATE(field1(maxrow,ncol1))
   ALLOCATE(field2(maxrow,ncol2))
   ALLOCATE(new_field(maxrow))

   iu_file=25
   OPEN(UNIT=iu_file,FILE=TRIM(file1), STATUS='OLD',ERR=1000, IOSTAT=ios)
   nline1=0
   DO 
      READ(iu_file,'(a)',ERR=100,END=100,IOSTAT=ios) line_read
      IF (line_read(1:1)=='#') CYCLE
      nline1=nline1+1
      READ(line_read,*) (field1(nline1,i),i=1,ncol1)
   ENDDO
100 CONTINUE
   CLOSE(iu_file,STATUS='KEEP')
   DO i=1,nline1
      field1(i,nx1)=field1(i,nx1)*fact1
      field1(i,ny1)=field1(i,ny1)*fact1
   ENDDO

   OPEN(UNIT=iu_file,FILE=TRIM(file2), STATUS='OLD',ERR=1000, IOSTAT=ios)
   nline2=0
   DO 
      READ(iu_file,'(a)',ERR=200,END=200,IOSTAT=ios) line_read
      IF (line_read(1:1)=='#') CYCLE
      nline2=nline2+1
      READ(line_read,*) (field2(nline2,i),i=1,ncol2)
   ENDDO
200 CONTINUE
   CLOSE(iu_file,STATUS='KEEP')

   DO i=1,nline2
      field2(i,nx2)=field2(i,nx2)*fact2
      field2(i,ny2)=field2(i,ny2)*fact2
   ENDDO

   CALL dosplineint(field1(1:nline1,nx1), field1(1:nline1,ny1),  &
                       field2(1:nline2,nx2), new_field(1:nline2)) 

   xmax=-1.D50
   xmin=1.D50
   DO i=1,nline1
      IF (field1(i,nx1)< xmin) xmin=field1(i,nx1)
      IF (field1(i,nx1)> xmax) xmax=field1(i,nx1)
   ENDDO

   OPEN(UNIT=iu_file,FILE=TRIM(fileout), STATUS='UNKNOWN',ERR=1000, IOSTAT=ios)
   DO i=1,nline2
      IF (field2(i,nx2)>=xmin.AND.field2(i,nx2)<=xmax) &
         WRITE(iu_file,'(3e20.12)') field2(i,nx2), new_field(i), field2(i,ny2)
   ENDDO
   CLOSE(iu_file,STATUS='KEEP')

   DEALLOCATE(field1)
   DEALLOCATE(field2)
   DEALLOCATE(new_field)
ENDIF
1000 CONTINUE

CALL environment_end( code )
CALL mp_global_end ()

END PROGRAM merge_interp
