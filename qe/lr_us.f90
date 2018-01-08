!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE lr_apply_s(vect, svect)
    !--------------------------------------------------------------------------
    !   
    ! This subroutine applies the S operator to vect and puts the result in 
    ! svect.
    !
    USE kinds,              ONLY : dp
    USE io_global,          ONLY : stdout
    USE uspp,               ONLY : okvan, vkb, nkb
    USE wvfct,              ONLY : npwx, nbnd
    USE klist,              ONLY : nks, xk, ngk, igk_k
    USE becmod,             ONLY : becp, calbec
    USE noncollin_module,   ONLY : npol
    USE lr_lanczos,         ONLY : rpert
    USE qpoint,             ONLY : nksq, ikks, ikqs
    USE control_lr,         ONLY : nbnd_occ
    
    IMPLICIT NONE
    !
    COMPLEX(dp), INTENT(IN)  ::  vect(npwx*npol,nbnd,nksq*rpert)
    COMPLEX(dp), INTENT(OUT) :: svect(npwx*npol,nbnd,nksq*rpert)
    !
    ! Local variables
    !
    INTEGER :: ik,    &
               ikk,   & ! index of the point k
               ikq,   & ! index of the point k+q
               npwq,  & ! number of the plane-waves at point k+q
               ibnd,  & ! index of bands
               ipert, & ! index on perturbations
               ikp      ! index on k and the perturbation
    
   IF (okvan) THEN
      DO ik = 1, nksq
         !
         ikk  = ikks(ik)
         ikq  = ikqs(ik)
         npwq = ngk(ikq)
         !
         ! Calculate beta-functions vkb at point k+q
         !
         CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb)
         !
         ! Calculate the product of beta-functions vkb with vect:
         ! becp%k = <vkb|vect>
         !
         DO ipert=1,rpert
            ikp = ik + nksq * (ipert-1)
            CALL calbec(npwq, vkb, vect(:,:,ikp), becp, nbnd_occ(ikk))
            !
            ! Apply the S operator
            !
            CALL s_psi(npwx, npwq, nbnd_occ(ikk), vect(:,:,ikp), svect(:,:,ikp))
         ENDDO
         !
      ENDDO
   ELSE
     svect=vect
   ENDIF
   !
   RETURN
   !
END SUBROUTINE lr_apply_s

