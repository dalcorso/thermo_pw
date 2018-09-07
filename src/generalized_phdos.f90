!
! Copyright (C) 2018 C. Malica and A. Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
SUBROUTINE generalized_phdos(et, nbnd, nks, wk, degauss, ngauss, e, &
                             displa, nat, dosg)
!------------------------------------------------------------------------------
!
  USE kinds,            ONLY : DP
  
  IMPLICIT NONE
  INTEGER :: nks, nbnd, ngauss

  REAL(DP) :: wk(nks), et(nbnd, nks), degauss, e, dosg(6,nat)
  REAL(DP) :: w0gauss, weight
  INTEGER :: n, ns, nk0, nk, ik
  INTEGER :: nspin0
  COMPLEX(DP) :: displa(3*nat,3*nat,nks)
  INTEGER :: ipol, jpol, ijpol, na, nat, indi, indj
  COMPLEX(DP) :: u1, u2, ufact
  EXTERNAL w0gauss
  !
  nk = nks
  !
  dosg = 0.0_DP
  DO ik = 1, nk
     DO n = 1, nbnd
        weight = w0gauss ( (E-et(n,ik) ) / Degauss, ngauss)
        DO na=1, nat
           ijpol=0
           DO ipol=1, 3
              indi=3*(na-1)+ipol
              DO jpol=ipol, 3
                 indj=3*(na-1)+jpol
                 ijpol=ijpol+1
                 u1=displa(indi,n,ik)
                 u2=displa(indj,n,ik)
                 ufact = DREAL(u1*CONJG(u2))
                 dosg(ijpol,na)=dosg(ijpol,na) + wk(ik)*ufact*weight
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  dosg = dosg / degauss
  !
  RETURN
END SUBROUTINE generalized_phdos
