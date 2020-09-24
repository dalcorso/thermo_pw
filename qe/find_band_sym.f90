!
! Copyright (C) 2006-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE find_band_sym_tpw (ik,evc,nsym,s,ftau,gk,invs,rap_et,times,&
                              ngroup,istart,accuracy)
!------------------------------------------------------------------------
  !
  !   This subroutine finds the irreducible representations which give
  !   the transformation properties of the wavefunctions.
  !   The output is in rap_et. rap_et contains the number of the 
  !   representation or zero if the latter could not be found.
  !   This routine does not support points at zone border if the space group of
  !   the crystal has fractionary translations (non-symmorphic space groups).
  !
  !
  USE io_global,       ONLY : stdout
  USE kinds,           ONLY : DP
  USE rap_point_group, ONLY : nclass, nelem, elem, which_irr, char_mat
  USE wvfct,           ONLY : nbnd, npwx
  USE klist,           ONLY : ngk, igk_k
  USE uspp,            ONLY : vkb, nkb, okvan
  USE becmod,          ONLY : bec_type, becp, calbec, allocate_bec_type, &
                              deallocate_bec_type
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : invfft
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum

  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: ik
  REAL(DP), INTENT(IN) :: accuracy

  INTEGER ::                  &
       nsym,             &
       rap_et(nbnd),     &
       ftau(3,48),       &
       gk(3,48),         &
       s(3,3,48),        &
       invs(48),         &
       ngroup,           &  ! number of different frequencies groups
       istart(nbnd+1)       ! position of the first band of each group

  COMPLEX(DP) ::         &
       d_spin_dum(2, 2), &
       times(nbnd,24),   &
       evc(npwx, nbnd)

  INTEGER ::      &
       ibnd,      &
       igroup,    &
       dim_rap,   &
       irot,      &
       irap,      &
       iclass,    &
       shift,     &
       mult,      &
       i,         &
       dimen,     &
       nrxx, npw, &
       has_e_dum

  COMPLEX(DP) :: zdotc

  REAL(DP) :: sumt
  COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), trace(:,:), psic(:,:), w(:,:)
  !
  !    Divide the bands on the basis of the band degeneracy.
  !
  nrxx=dfftp%nnr
  ALLOCATE(evcr(npwx,nbnd))
  ALLOCATE(psic(nrxx,nbnd))
  ALLOCATE(trace(48,nbnd))
  ALLOCATE(w(48,nbnd))
  IF (okvan) CALL allocate_bec_type ( nkb, nbnd, becp )

  rap_et=0
  has_e_dum=1
  d_spin_dum=(0.0_DP,0.0_DP)
!
!   bring all the bands in real space
!
  npw = ngk(ik)
  psic=(0.0_DP,0.0_DP)
  DO ibnd=1,nbnd
     psic(dfftp%nl(igk_k(1:npw,ik)),ibnd) = evc(1:npw,ibnd)
     CALL invfft ('Rho', psic(:,ibnd), dfftp)
  ENDDO
  !
  !  Find the character of one symmetry operation per class
  !
  DO iclass=1,nclass
     irot=elem(1,iclass)
     !
     !   Rotate all the bands together.
     !   NB: rotate_psi assume that s is in the small group of k. 
     !
     !
     IF (irot==1) THEN
        evcr=evc
     ELSE
        CALL rotate_all_psi_tpw(ik,psic,evcr,s(1,1,invs(irot)), &
                  ftau(1,invs(irot)),d_spin_dum,has_e_dum,gk(1,invs(irot)))
     ENDIF
     !
     !   and apply S if necessary
     !
     IF ( okvan ) THEN
        CALL calbec( npw, vkb, evcr, becp )
        CALL s_psi( npwx, npw, nbnd, evcr, evcr )
     ENDIF
     !
     !  find the diagonal element of the representation of dimension nbnd
     !
     DO ibnd=1,nbnd
        w(iclass,ibnd) = ZDOTC(npw,evc(1,ibnd),1,evcr(1,ibnd),1)
     ENDDO
  ENDDO
  !
  CALL mp_sum( w, intra_bgrp_comm )
  !
  !  Divide the bands in groups of degenerate eigenvalues. 
  !  Computes the trace for each group of degenerate modes.
  !  We continue to add diagonal elements to the trace, until we
  !  find a set of traces whose sum of square moduli over all elements
  !  of the group is an integer multiple of the group order. 
  !
  ! 
  ngroup=1
  istart(ngroup)=1
  trace=(0.d0,0.d0)
  DO ibnd=1, nbnd
     DO iclass=1,nclass
        trace(iclass,ngroup)=trace(iclass,ngroup) + w(iclass, ibnd)
     ENDDO
     sumt=0.0_DP
     DO iclass=1,nclass
        sumt=sumt+ABS(trace(iclass,ngroup))**2*nelem(iclass)
     ENDDO
     sumt=sumt/nsym
!
!    If sumt is an integer we have found an irreducible representation or
!    an integer number of irreducible representations. We can start to 
!    identify a new group of modes.
!
     IF (ABS(NINT(sumt)-sumt) < 1.d-5) THEN
        ngroup=ngroup+1
        istart(ngroup)=ibnd+1
     ENDIF
  ENDDO
  ngroup=ngroup-1

  !DO iclass=1,nclass
  !   write(6,'(i5,3(2f11.8,1x))') iclass,trace(iclass,4),trace(iclass,5), &
  !                                       trace(iclass,6)
  !ENDDO
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of bands
  !
  DO igroup=1,ngroup
     dim_rap=istart(igroup+1)-istart(igroup)
     shift=0
     DO irap=1,nclass
        times(igroup,irap)=(0.d0,0.d0)
        DO iclass=1,nclass
           times(igroup,irap)=times(igroup,irap) &
                +trace(iclass,igroup)*CONJG(char_mat(irap,which_irr(iclass)))&
                *nelem(iclass)
        ENDDO
        times(igroup,irap)=times(igroup,irap)/nsym
        mult=NINT(DBLE(times(igroup,irap)))
        IF ((ABS(mult-DBLE(times(igroup,irap))) > accuracy).OR. &
            (ABS(AIMAG(times(igroup,irap))) > accuracy) ) THEN
           GOTO 300
        ELSEIF (ABS(times(igroup,irap)) > accuracy) THEN
           ibnd=istart(igroup)+shift
           dimen=NINT(DBLE(char_mat(irap,1)))
           DO i=1,dimen*mult
              ibnd=istart(igroup)+shift+i-1
              rap_et(ibnd)=irap
           ENDDO
           shift=shift+dimen*mult
        ENDIF
     ENDDO
300  CONTINUE
  ENDDO

  DEALLOCATE(trace)
  DEALLOCATE(w)
  DEALLOCATE(evcr)
  DEALLOCATE(psic)
  IF (okvan) CALL deallocate_bec_type (becp)

  RETURN
END SUBROUTINE find_band_sym_tpw

!------------------------------------------------------------------------
SUBROUTINE find_band_sym_so_tpw (ik,evc,nsym,s,ftau,d_spin,gk, &
                                 invs,rap_et,times,ngroup,istart,accuracy)
!------------------------------------------------------------------------
  !
  !   This subroutine finds the irreducible representations of the
  !   double group which give the transformation properties of the
  !   spinor wavefunctions evc.
  !   The output is in rap_et. rap_et contains the number of the 
  !   representation or zero if the latter could not be found.
  !   This routine does not support points at zone border if the space group of
  !   the crystal has fractionary translations (non-symmorphic space groups).
  !
  USE io_global,          ONLY : stdout
  USE kinds,              ONLY : DP
  USE rap_point_group,    ONLY : nclass
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, which_irr_so, &
                                 char_mat_so
  USE wvfct,              ONLY : nbnd, npwx
  USE klist,              ONLY : ngk, igk_k
  USE fft_base,           ONLY : dfftp
  USE fft_interfaces,     ONLY : invfft
  USE uspp,               ONLY : vkb, nkb, okvan
  USE noncollin_module,   ONLY : npol
  USE becmod,             ONLY : bec_type, becp, calbec, allocate_bec_type, &
                                 deallocate_bec_type
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp,                 ONLY : mp_sum

  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: ik
  REAL(DP), INTENT(IN) :: accuracy

  INTEGER ::             &
       nsym,             &
       ngroup,           &
       istart(nbnd+1),   &
       rap_et(nbnd),     &
       ftau(3,48),       &
       gk(3,48),         &
       invs(48),         &
       s(3,3,48)

  COMPLEX(DP) ::         &
       times(nbnd,24),   &
       d_spin(2,2,48),   &
       evc(npwx*npol, nbnd)

  INTEGER ::      &
       ibnd,      &
       igroup,    &
       dim_rap,   &  ! counters
       irot,      &
       irap,      &
       shift,     &
       iclass,    &
       mult,      &
       i,         &
       ipol,      & 
       dimen,     &
       npw,       &
       ind1, ind2

  COMPLEX(DP) :: zdotc   ! scalar product routine

  REAL(DP) :: sumt

  COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), & ! the rotated of each wave function
       psic(:,:,:),   &! the wavefunctions in real space
       w(:,:),        &! the diagonals of the big representation
       trace(:,:)      ! the trace of the symmetry matrix
  !
  ALLOCATE(evcr(npwx*npol,nbnd))
  ALLOCATE(psic(dfftp%nnr,npol,nbnd))
  ALLOCATE(trace(48,nbnd))
  ALLOCATE(w(48,nbnd))
  IF (okvan) CALL allocate_bec_type ( nkb, nbnd, becp )
  !
  !  Bring each wavefunction in real space on the thick mesh
  !
  rap_et=0
  npw = ngk(ik)
  psic=(0.0_DP,0.0_DP)
  DO ipol=1, npol
     ind1 = 1 + (ipol-1)*npwx
     ind2 = npw + (ipol-1)*npwx
     DO ibnd=1,nbnd
        psic(dfftp%nl(igk_k(1:npw,ik)),ipol,ibnd) = evc(ind1:ind2,ibnd)
        CALL invfft ('Rho', psic(:,ipol,ibnd), dfftp)
     ENDDO
  ENDDO
  !
  !  Rotate each band with a symmetry element per class and find the
  !  diagonal element of the matrix representation of dimension nbnd
  !
  DO iclass=1,nclass
     irot=elem_so(1,iclass)
     !
     !   Rotate all the bands together.
     !   NB: rotate_psi assumes that s is in the small group of k. 
     !
      evcr=(0.0_DP,0.0_DP)
      CALL rotate_all_psi_tpw(ik,psic,evcr,s(1,1,invs(irot)),        &
           ftau(1,invs(irot)),d_spin(1,1,irot),has_e(1,iclass),gk(1,invs(irot)))
     !
     !   and apply S in the US case.
     !
     IF ( okvan ) THEN
        CALL calbec( npw, vkb, evcr, becp )
        CALL s_psi( npwx, npw, nbnd, evcr, evcr )
     ENDIF
     !
     !  find the diagonal element of the representation of dimension nbnd
     !
     DO ibnd=1,nbnd
        w(iclass,ibnd) = ZDOTC(npwx*npol,evc(1,ibnd),1,evcr(1,ibnd),1)
     ENDDO
  ENDDO
  !
  CALL mp_sum(w,intra_bgrp_comm)
  !
  !  Divide the bands in groups of degenerate eigenvalues. 
  !  Computes the trace for each group of degenerate modes.
  !  We continue to add diagonal elements to the trace, until we
  !  find a set of traces whose sum of square moduli over all elements
  !  of the group is an integer multiple of the group order. 
  !
  ngroup=1
  istart(ngroup)=1
  trace=(0.d0,0.d0)
  DO ibnd=1, nbnd
     DO iclass=1,nclass
        trace(iclass,ngroup)=trace(iclass,ngroup) + w(iclass, ibnd)
     ENDDO
     sumt=0.0_DP
     DO iclass=1,nclass
        sumt=sumt+ABS(trace(iclass,ngroup))**2*DBLE(nelem_so(iclass))
     ENDDO
     sumt=sumt/(2*nsym)
!
!    If sumt is an integer we have found an irreducible representation or
!    an integer number of irreducible representations. We can start to 
!    identify a new group of modes.
!
     IF (ABS(NINT(sumt)-sumt) < 1.d-5) THEN
        ngroup=ngroup+1
        istart(ngroup)=ibnd+1
     ENDIF
  ENDDO
  ngroup=ngroup-1
  !
!  DO iclass=1,nclass
!     write(6,'(i5,3(2f11.8,1x))') iclass,trace(iclass,1),trace(iclass,2), &
!                                         trace(iclass,3)
!  ENDDO
!
!  Use the character tables of the irreducible representation to identify
!  the symmetry of the representation for each group of degenerate bands
!
  DO igroup=1,ngroup
     dim_rap=istart(igroup+1)-istart(igroup)
     shift=0
     DO irap=1,nrap
        times(igroup,irap)=(0.d0,0.d0)
        DO iclass=1,nclass
           times(igroup,irap)=times(igroup,irap) &
                +trace(iclass,igroup)*CONJG(char_mat_so(irap, &
                which_irr_so(iclass)))*DBLE(nelem_so(iclass))
        ENDDO
        times(igroup,irap)=times(igroup,irap)/2/nsym
        mult=NINT(DBLE(times(igroup,irap)))
        IF ((ABS(mult-DBLE(times(igroup,irap)))> accuracy).OR. &
            (ABS(AIMAG(times(igroup,irap))) > accuracy) ) THEN
           GOTO 300
        ENDIF
        IF (ABS(times(igroup,irap)) > accuracy) THEN
           dimen=NINT(DBLE(char_mat_so(irap,1)))
           ibnd=istart(igroup) + shift
           DO i=1,dimen*mult
              ibnd=istart(igroup)+shift+i-1
              rap_et(ibnd)=irap
           ENDDO
           shift=shift+dimen*mult
        ENDIF
     ENDDO
300  CONTINUE
  ENDDO

  DEALLOCATE(trace)
  DEALLOCATE(w)
  DEALLOCATE(psic)
  DEALLOCATE(evcr)
  IF (okvan) CALL deallocate_bec_type ( becp )
  RETURN
END SUBROUTINE find_band_sym_so_tpw
