!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE readmodes (nat,nks,q,displa,ngeo,igeo,ntyp,ityp,amass,iflag,inunit)
  !-----------------------------------------------------------------------
  !
  !   read modes from a file and transform them into eigenvectors
  !   of the dynamical matrix if iflag==1, otherwise assumes to
  !   read the eigenvectors and copy them in displa
  !
  USE kinds, ONLY: dp
  USE constants, ONLY : amu_ry
  IMPLICIT NONE
  ! input
  INTEGER, INTENT(IN) :: nat, nks, igeo, ngeo, ntyp, iflag, inunit
  INTEGER, INTENT(IN) :: ityp(nat)
  REAL(DP), INTENT(IN) :: q(3,nks), amass(ntyp)
  COMPLEX(DP), INTENT(INOUT) :: displa(3*nat,3*nat,ngeo,nks)
  ! local
  INTEGER nat3, na, ipol, iq, nta, sna, i, j
  REAL(DP):: freq(3*nat), inq(3)
  REAL(DP):: znorm
  REAL(DP) :: eps=1.d-3
  !
  nat3=3*nat
  !
  !  read the normalised displacements
  !
  DO iq=1,nks
     READ(inunit,*)
     READ(inunit,*)
     READ(inunit,'(5x,3f12.4)') inq(1:3)
     IF (ABS(inq(1)-q(1,iq))+ABS(inq(2)-q(2,iq))+ABS(inq(3)-q(3,iq)) > eps ) &
        CALL errore('readmodes','Incompatible q points in modes file',1)
     READ(inunit,*)
     DO i = 1, nat3
      
        READ (inunit,*)
        DO na = 1,nat
           nta = ityp(na)
           sna=(na-1)*3
           READ (inunit,9020) (displa(sna+ipol,i,igeo,iq),ipol=1,3)
           IF (iflag==1) THEN
              DO ipol=1,3
                 displa(sna+ipol,i,igeo,iq)=displa(sna+ipol,i,igeo,iq)* &
                                         SQRT(amu_ry*amass(nta))
              END DO
           END IF
        END DO
        IF (iflag==1) THEN
           znorm = 0.0d0
           DO j=1,nat3
              znorm=znorm+ABS(displa(j,i,igeo,iq))**2
           END DO
           znorm = SQRT(znorm)
           displa(:,i,igeo,iq)=displa(:,i,igeo,iq)/znorm
        END IF
     END DO   
     READ(inunit,*)
  END DO
  !
  RETURN
  !
9020 FORMAT (2x,3(f10.6,1x,f10.6,3x))
  !
END SUBROUTINE readmodes
