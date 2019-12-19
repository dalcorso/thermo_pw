!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dyndiag_tpw (nat,ntyp,amass,ityp,dyn,w2,z,flag)
  !-----------------------------------------------------------------------
  !
  !   diagonalise the dynamical matrix
  !   On input:  amass = masses, in amu
  !   On output: w2 = energies, z = displacements (flag=1)
  !                             z = eigenvectors  (flag/=1)
  !
  use kinds, only: dp
  use constants, only: amu_ry
  implicit none
  ! input
  integer nat, ntyp, ityp(nat)
  complex(DP) dyn(3,3,nat,nat)
  real(DP) amass(ntyp)
  integer :: flag
  ! output
  real(DP) w2(3*nat)
  complex(DP) z(3*nat,3*nat)
  ! local
  real(DP) diff, dif1, difrel
  integer nat3, na, nta, ntb, nb, ipol, jpol, i, j
  complex(DP), allocatable :: dyn2(:,:)
  !
  !  fill the two-indices dynamical matrix
  !
  nat3 = 3*nat
  allocate(dyn2 (nat3, nat3))
  !
  do na = 1,nat
     do nb = 1,nat
        do ipol = 1,3
           do jpol = 1,3
              dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = dyn(ipol,jpol,na,nb)
           end do
        end do
     end do
  end do
  !
  !  impose hermiticity
  !
  diff = 0.d0
  difrel=0.d0
  do i = 1,nat3
     dyn2(i,i) = CMPLX( DBLE(dyn2(i,i)),0.d0,kind=DP)
     do j = 1,i - 1
        dif1 = abs(dyn2(i,j)-CONJG(dyn2(j,i)))
        if ( dif1 > diff .and. &
             max ( abs(dyn2(i,j)), abs(dyn2(j,i))) > 1.0d-6) then
           diff = dif1
           difrel=diff / min ( abs(dyn2(i,j)), abs(dyn2(j,i)))
        end if
        dyn2(i,j) = 0.5d0* (dyn2(i,j)+CONJG(dyn2(j,i)))
        dyn2(j,i) = CONJG(dyn2(i,j))
     end do
  end do
  if ( diff > 1.d-6 ) write (6,'(5x,"Max |d(i,j)-d*(j,i)| = ",f9.6,/,5x, &
       & "Max |d(i,j)-d*(j,i)|/|d(i,j)|: ",f8.4,"%")') diff, difrel*100
  !
  !  divide by the square root of masses
  !
  do na = 1,nat
     nta = ityp(na)
     do nb = 1,nat
        ntb = ityp(nb)
        do ipol = 1,3
           do jpol = 1,3
             dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = &
                  dyn2((na-1)*3+ipol, (nb-1)*3+jpol) / &
                  (amu_ry*sqrt(amass(nta)*amass(ntb)))
          end do
       end do
    end do
 end do
 !
 !  diagonalisation
 !
 call cdiagh2_tpw(nat3,dyn2,nat3,w2,z)
 !
 deallocate(dyn2)
 !
 !  displacements are eigenvectors divided by sqrt(amass)
 !  if flag=1 they are given as output of the routine, otherwise
 !  the output are the eigenvectors of the dynamical matrix
 !
 IF (flag==1) THEN
    do i = 1,nat3
       do na = 1,nat
          nta = ityp(na)
          do ipol = 1,3
             z((na-1)*3+ipol,i) = z((na-1)*3+ipol,i)/ sqrt(amu_ry*amass(nta))
          end do
       end do
    end do
 ENDIF
 !
 return
end subroutine dyndiag_tpw
!
!-----------------------------------------------------------------------
subroutine cdiagh2_tpw (n,h,ldh,e,v)
  !-----------------------------------------------------------------------
  !
  !   calculates all the eigenvalues and eigenvectors of a complex
  !   hermitean matrix H . On output, the matrix is unchanged
  !
  use kinds, only: dp
  implicit none
  !
  ! on INPUT
  integer          n,       &! dimension of the matrix to be diagonalized
       &           ldh       ! leading dimension of h, as declared
  ! in the calling pgm unit
  complex(DP)  h(ldh,n)  ! matrix to be diagonalized
  !
  ! on OUTPUT
  real(DP)     e(n)      ! eigenvalues
  complex(DP)  v(ldh,n)  ! eigenvectors (column-wise)
  !
  ! LOCAL variables (LAPACK version)
  !
  integer          lwork,   &! aux. var.
       &           ILAENV,  &! function which gives block size
       &           nb,      &! block size
       &           info      ! flag saying if the exec. of libr. routines was ok
  !
  real(DP), allocatable::   rwork(:)
  complex(DP), allocatable:: work(:)
  !
  !     check for the block size
  !
  nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
  if (nb.lt.1) nb=max(1,n)
  if (nb.eq.1.or.nb.ge.n) then
     lwork=2*n-1
  else
     lwork = (nb+1)*n
  endif
  !
  ! allocate workspace
  !
  call zcopy(n*ldh,h,1,v,1)
  allocate(work (lwork))
  allocate(rwork (3*n-2))
  call ZHEEV('V','U',n,v,ldh,e,work,lwork,rwork,info)
  call errore ('cdiagh2_tpw','info =/= 0',abs(info))
  ! deallocate workspace
  deallocate(rwork)
  deallocate(work)
  !
  return
end subroutine cdiagh2_tpw

