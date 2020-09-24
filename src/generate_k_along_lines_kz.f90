!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE generate_k_along_lines_kz(nkaux, xkaux, wkaux, xk, wk, nkstot, nkz)
!---------------------------------------------------------------------------
!
!  This routine recieves as input a set of k point (xkaux) and integer weights
!  (wkaux) and generates a set of k points along the lines 
!  xkaux(:,i+1)-xkaux(:,i). Each line contains wkaux(i) points.
!  The weights of each k point wk(i) is the lenght of the path from xk(:,1)
!  to xk(i). Points with wkaux=0 do not increase the path lenght.
!  In this routines the points xkaux are contained in the xy plane.
!  The kz coordinate is chosen so that there are nkz copies of the path,
!  one at a different kz. This routine is used for the projected band
!  structure.
!  The total number of output points must be nkstot, and xk and wk must
!  be an array of lenght nkstot.
!
USE kinds, ONLY : DP
USE cell_base, ONLY : bg
USE thermo_sym, ONLY : code_group_save
USE point_group, ONLY : has_sigma_h
USE control_paths, ONLY : nrap_plot_in, rap_plot_in, nrap_plot, rap_plot
IMPLICIT NONE

INTEGER, INTENT(IN) :: nkaux, nkstot, wkaux(nkaux), nkz
REAL(DP), INTENT(IN) :: xkaux(3,nkaux)
REAL(DP), INTENT(OUT) :: xk(3,nkstot), wk(nkstot)

INTEGER :: nkstot_, i, j, ikz
REAL(DP) :: delta, xkmod, deltakz(3)
LOGICAL :: type1

!
!  If the point group of the slab has sigma_h we choose only positive
!  k_z, otherwise we have to sample both positive and negative k_z
!
nkstot_= 0
IF (nkz > 1) THEN
   type1=has_sigma_h(code_group_save)
   IF (type1) THEN
      deltakz(:) = bg(:,3) * 0.5_DP / (nkz - 1)
   ELSE
      deltakz(:) = bg(:,3) / nkz
   ENDIF
ELSE
   type1=.TRUE.
   deltakz(:) = 0.0_DP 
ENDIF
DO ikz=1,nkz
   nkstot_=nkstot_+1
   wk(nkstot_)=0.0_DP
   IF (type1) THEN
      xk(:,nkstot_)=xkaux(:,1) + deltakz(:) * ( ikz - 1 )
   ELSE
      xk(:,nkstot_)=xkaux(:,1) + deltakz(:) * ikz  - bg(:,3) * 0.5_DP
   ENDIF
   nrap_plot(nkstot_)=nrap_plot_in(1)
   IF (nrap_plot(nkstot_)>0) &
      rap_plot(1:nrap_plot(nkstot_),nkstot_)=rap_plot_in(1:nrap_plot(nkstot_),1)
   DO i=2,nkaux
      IF (wkaux(i-1)>0) THEN
         delta=1.0_DP/wkaux(i-1)
         DO j=1,wkaux(i-1)
            nkstot_=nkstot_+1
            IF (nkstot_ > nkstot) CALL errore ('generate_k_along_lines_kz', &
                                           'internal error 1: wrong nkstot',i)
            IF (type1) THEN
               xk(:,nkstot_)= xkaux(:,i-1)+delta*j*(xkaux(:,i)-xkaux(:,i-1)) &
                         + deltakz(:) * ( ikz - 1 )
            ELSE
               xk(:,nkstot_)= xkaux(:,i-1)+delta*j*(xkaux(:,i)-xkaux(:,i-1)) &
                         + deltakz(:) * ikz - bg(:,3) * 0.5_DP
            ENDIF
            xkmod=SQRT( (xk(1,nkstot_)-xk(1,nkstot_-1))**2 +   &
                        (xk(2,nkstot_)-xk(2,nkstot_-1))**2 +   &
                        (xk(3,nkstot_)-xk(3,nkstot_-1))**2 )
            wk(nkstot_)=wk(nkstot_-1) + xkmod
            nrap_plot(nkstot_)=nrap_plot_in(i-1)
            IF (nrap_plot(nkstot_)>0) &
               rap_plot(1:nrap_plot(nkstot_),nkstot_)=&
                         rap_plot_in(1:nrap_plot(nkstot_),i-1)
         ENDDO
      ELSEIF (wkaux(i-1)==0) THEN
         nkstot_=nkstot_+1
         IF (nkstot_ > nkstot) CALL errore ('generate_k_along_lines_kz', &
                                           'internal error 2: wrong nkstot',i)
         IF (nkstot_ ==1 ) CALL errore ('generate_k_along_lines_kz', &
                                            'problems with weights',i)
         IF (type1) THEN
            xk(:,nkstot_)=xkaux(:,i) + deltakz(:) * ( ikz - 1 )
         ELSE
            xk(:,nkstot_)=xkaux(:,i) + deltakz(:) * ikz - bg(:,3) * 0.5_DP
         ENDIF
         wk(nkstot_)=wk(nkstot_-1) 
         nrap_plot(nkstot_)=nrap_plot_in(i-1)
         IF (nrap_plot(nkstot_)>0) &
            rap_plot(1:nrap_plot(nkstot_),nkstot_)=&
                                   rap_plot_in(1:nrap_plot(nkstot_),i-1)
      ELSE
         CALL errore ('generate_k_along_lines_kz', 'wrong number of points',i)
      ENDIF
   ENDDO
END DO
IF (nkstot_ /= nkstot) CALL errore ('generate_k_along_lines_kz', &
                                    'internal error 3: wrong nkstot',nkstot_)

RETURN
END SUBROUTINE generate_k_along_lines_kz


