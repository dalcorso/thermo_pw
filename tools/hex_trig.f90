PROGRAM hex_trig
!
!  This program reads the parameters a and c of an hexagonal lattice
!  (in Angstrom) and gives as output the parameters a and cos(gamma)
!  of a trigonal lattice equivalent to the hexagonal one.
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP)  :: a, c, at, cgammat, csua
INTEGER   :: stdin, stdout

stdin=5
stdout=6
WRITE(stdout,'(5x,"insert a and c of hexagonal lattice in A")')  
READ(stdin,*) a, c  

at=SQRT(3.0_DP*a**2+c**2) / 3.0_DP
csua=c/a
cgammat = (csua**2 -1.5_DP) / (csua**2 + 3.0_DP)

WRITE(stdout,'(5x,"at=",f12.7," gamma=",f12.7, " cos(gamma)=",f12.7)') at, ACOS(cgammat), &
                                                                  cgammat
WRITE(stdout,'(5x,"celldm(1)=",f12.7)') at / 0.529177_DP
WRITE(stdout,'(5x,"celldm(4)=",f12.7)') cgammat

END PROGRAM hex_trig
