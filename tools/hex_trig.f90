PROGRAM hex_trig
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP)  :: a, c, at, cgammat, csua

WRITE(6,'(5x,"insert a and c of hexagonal lattice in A")')  
READ(5,*) a, c  

at=SQRT(3.0_DP*a**2+c**2) / 3.0_DP
csua=c/a
cgammat = (csua**2 -1.5_DP) / (csua**2 + 3.0_DP)

WRITE(6,'(5x,"at=",f12.7," gamma=",f12.7, " cos(gamma)=",f12.7)') at, ACOS(cgammat), &
                                                                  cgammat
WRITE(6,'(5x,"celldm(1)=",f12.7)') at / 0.529177_DP
WRITE(6,'(5x,"celldm(4)=",f12.7)') cgammat

END PROGRAM hex_trig
