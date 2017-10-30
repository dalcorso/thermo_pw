!
! Copyright (C) 2006-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_band_sym_tpw (ik,evc,et,nsym,s,ftau,gk,invs,rap_et,times,&
                          ngroup,istart,accuracy)
  !
  !   This subroutine finds the irreducible representations which give
  !   the transformation properties of the wavefunctions.
  !   Presently it does NOT work at zone border if the space group of
  !   the crystal has fractionary translations (non-symmorphic space groups).
  !
  !
  USE io_global,       ONLY : stdout
  USE kinds,           ONLY : DP
  USE constants,       ONLY : rytoev
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat
  USE gvect,           ONLY : nl
  USE wvfct,           ONLY : nbnd, npwx
  USE klist,           ONLY : ngk, igk_k
  USE uspp,            ONLY : vkb, nkb, okvan
  USE becmod,          ONLY : bec_type, becp, calbec, &
       allocate_bec_type, deallocate_bec_type
  USE fft_base,        ONLY : dfftp
  USE fft_interfaces,  ONLY : invfft
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ik
  REAL(DP), INTENT(in) :: accuracy

  INTEGER ::                  &
       nsym,             &
       rap_et(nbnd),     &
       ftau(3,48),       &
       gk(3,48),         &
       s(3,3,48),        &
       invs(48),         &
       ngroup,           &  ! number of different frequencies groups
       istart(nbnd+1)

  REAL(DP) ::                 &
       et(nbnd)

  COMPLEX(DP) ::  &
       d_spin_dum(2, 2), &
       times(nbnd,24), &
       evc(npwx, nbnd)

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::      &
       ibnd,      &
       igroup,    &
       dim_rap,   &
       irot,      &
       irap,      &
       iclass,    &
       shift,     &
       na, i, j,  &
       ig, dimen, &
       nrxx, npw, &
       has_e_dum

  COMPLEX(DP) :: zdotc

  REAL(DP), ALLOCATABLE ::  w1(:)
  COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), trace(:,:), psic(:,:)
  !
  !    Divide the bands on the basis of the band degeneracy.
  !
  nrxx=dfftp%nnr
  ALLOCATE(w1(nbnd))
  ALLOCATE(evcr(npwx,nbnd))
  ALLOCATE(psic(nrxx,nbnd))
  ALLOCATE(trace(48,nbnd))
  IF (okvan) CALL allocate_bec_type ( nkb, nbnd, becp )

  rap_et=-1
  w1=et*rytoev

  ngroup=1
  istart(ngroup)=1
  DO ibnd=2,nbnd
     IF (abs(w1(ibnd)-w1(ibnd-1)) > 0.0001d0) THEN
        ngroup=ngroup+1
        istart(ngroup)=ibnd
     ENDIF
  ENDDO
  istart(ngroup+1)=nbnd+1
!
!   bring all the bands in real space
!
  npw = ngk(ik)
  psic=(0.0_DP,0.0_DP)
  DO ibnd=1,nbnd
     psic(nl(igk_k(1:npw,ik)),ibnd) = evc(1:npw,ibnd)
     CALL invfft ('Dense', psic(:,ibnd), dfftp)
  ENDDO
  !
  !  Find the character of one symmetry operation per class
  !
  DO iclass=1,nclass
     irot=elem(1,iclass)
     !
     !   Rotate all the bands together.
     !   NB: rotate_psi assume that s is in the small group of k. It does not
     !       rotate the k point.
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
     !  Compute the trace of the representation for each group of bands
     !
     DO igroup=1,ngroup
        dim_rap=istart(igroup+1)-istart(igroup)
        trace(iclass,igroup)=(0.d0,0.d0)
        DO i=1,dim_rap
           ibnd=istart(igroup)+i-1
           trace(iclass,igroup)=trace(iclass,igroup) + &
                zdotc(npw,evc(1,ibnd),1,evcr(1,ibnd),1)
        ENDDO
        !      write(6,*) igroup, iclass, trace(iclass,igroup)
     ENDDO
  ENDDO
  !
  CALL mp_sum( trace, intra_bgrp_comm )

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
        IF ((abs(nint(dble(times(igroup,irap)))-dble(times(igroup,irap))) &
             > accuracy).OR. (abs(aimag(times(igroup,irap))) > eps) ) THEN
           ibnd=istart(igroup)
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dim_rap
                 ibnd=istart(igroup)+i-1
                 rap_et(ibnd)=0
              ENDDO
           ENDIF
           GOTO 300
        ELSEIF (abs(times(igroup,irap)) > accuracy) THEN
           ibnd=istart(igroup)+shift
           dimen=nint(dble(char_mat(irap,1)))
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dimen*nint(dble(times(igroup,irap)))
                 ibnd=istart(igroup)+shift+i-1
                 rap_et(ibnd)=irap
              ENDDO
              shift=shift+dimen*nint(dble(times(igroup,irap)))
           ENDIF
        ENDIF
     ENDDO
300  CONTINUE
  ENDDO

  DEALLOCATE(trace)
  DEALLOCATE(w1)
  DEALLOCATE(evcr)
  DEALLOCATE(psic)
  IF (okvan) CALL deallocate_bec_type (becp)

  RETURN
END SUBROUTINE find_band_sym_tpw

SUBROUTINE find_band_sym_so_tpw (ik,evc,et,nsym,s,ftau,d_spin,gk, &
     invs,rap_et,times,ngroup,istart,accuracy)
  !
  !   This subroutine finds the irreducible representations of the
  !   double group which give the transformation properties of the
  !   spinor wavefunctions evc.
  !   Presently it does NOT work at zone border if the space group of
  !   the crystal has fractionary translations (non-symmorphic space groups).
  !
  !
  USE io_global,          ONLY : stdout
  USE kinds,              ONLY : DP
  USE constants,          ONLY : rytoev
  USE rap_point_group,    ONLY : code_group, nclass
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, which_irr_so, &
       char_mat_so
  USE gvect,              ONLY : nl
  USE wvfct,              ONLY : nbnd, npwx
  USE klist,              ONLY : ngk, igk_k
  USE fft_base,           ONLY : dfftp
  USE fft_interfaces,     ONLY : invfft
  USE spin_orb,           ONLY : domag
  USE uspp,               ONLY : vkb, nkb, okvan
  USE noncollin_module,   ONLY : npol
  USE becmod,             ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp,                 ONLY : mp_sum

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ik
  REAL(DP), INTENT(in) :: accuracy

  INTEGER ::                  &
       nsym,             &
       ngroup,           &
       istart(nbnd+1),   &
       rap_et(nbnd),     &
       ftau(3,48),       &
       gk(3,48),         &
       invs(48),         &
       s(3,3,48)

  REAL(DP) ::                 &
       et(nbnd)

  COMPLEX(DP) ::  &
       times(nbnd,24), &
       d_spin(2,2,48), &
       evc(npwx*npol, nbnd)

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::         &
       ibnd,      &
       igroup,    &
       dim_rap,   &  ! counters
       irot,      &
       irap,      &
       shift,     &
       iclass,    &
       na, i, j, ig, ipol, jpol, jrap, dimen, npw, ind1, ind2

  COMPLEX(DP) :: zdotc          ! moltiplication factors

  REAL(DP), ALLOCATABLE ::  w1(:)      ! list of energy eigenvalues in eV
  COMPLEX(DP), ALLOCATABLE ::  evcr(:,:), & ! the rotated of each wave function
       psic(:,:,:),   &! the wavefunctions in real space
       trace(:,:)   ! the trace of the symmetry matrix
  ! within a given group
  !
  !    Divide the bands on the basis of the band degeneracy.
  !
  ALLOCATE(w1(nbnd))
  ALLOCATE(evcr(npwx*npol,nbnd))
  ALLOCATE(psic(dfftp%nnr,npol,nbnd))
  ALLOCATE(trace(48,nbnd))
  IF (okvan) CALL allocate_bec_type ( nkb, nbnd, becp )

  rap_et=-1
  w1=et*rytoev
  npw = ngk(ik)
  psic=(0.0_DP,0.0_DP)
  DO ipol=1, npol
     ind1 = 1 + (ipol-1)*npwx
     ind2 = npw + (ipol-1)*npwx
     DO ibnd=1,nbnd
        psic(nl(igk_k(1:npw,ik)),ipol,ibnd) = evc(ind1:ind2,ibnd)
        CALL invfft ('Dense', psic(:,ipol,ibnd), dfftp)
     ENDDO
  ENDDO
  !
  !  divide the energies in groups of degenerate eigenvalues. Two eigenvalues
  !  are assumed to be degenerate if their difference is less than 0.0001 eV.
  !
  ngroup=1
  istart(ngroup)=1
  DO ibnd=2,nbnd
     IF (abs(w1(ibnd)-w1(ibnd-1)) > 0.0001d0) THEN
        ngroup=ngroup+1
        istart(ngroup)=ibnd
     ENDIF
  ENDDO
  istart(ngroup+1)=nbnd+1

  trace=(0.d0,0.d0)
  DO iclass=1,nclass
     irot=elem_so(1,iclass)
     !
     !   Rotate all the bands together.
     !   NB: rotate_psi assumes that s is in the small group of k. It does not
     !       rotate the k point.
     !
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
     !  Compute the trace of the representation for each group of bands
     !
     DO igroup=1,ngroup
        dim_rap=istart(igroup+1)-istart(igroup)
        DO i=1,dim_rap
           ibnd=istart(igroup)+i-1
           trace(iclass,igroup)=trace(iclass,igroup) +            &
                zdotc(2*npwx,evc(1,ibnd),1,evcr(1,ibnd),1)
        ENDDO
        !      write(6,*) igroup, iclass, dim_rap, trace(iclass,igroup)
     ENDDO
  ENDDO
  !
  CALL mp_sum(trace,intra_bgrp_comm)
  !
!  DO iclass=1,nclass
!     write(6,'(i5,3(2f11.8,1x))') iclass,trace(iclass,1),trace(iclass,2), &
!                                         trace(iclass,3)
!  ENDDO

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

        IF ((abs(nint(dble(times(igroup,irap)))-dble(times(igroup,irap)))&
             > accuracy).or. (abs(aimag(times(igroup,irap))) > accuracy) ) THEN
           ibnd=istart(igroup)
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dim_rap
                 ibnd=istart(igroup)+i-1
                 rap_et(ibnd)=0
              ENDDO
           ENDIF
           GOTO 300
        ENDIF
        IF (abs(times(igroup,irap)) > accuracy) THEN
           dimen=nint(dble(char_mat_so(irap,1)))
           ibnd=istart(igroup) + shift
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dimen*nint(dble(times(igroup,irap)))
                 ibnd=istart(igroup)+shift+i-1
                 rap_et(ibnd)=irap
              ENDDO
              shift=shift+dimen*nint(dble(times(igroup,irap)))
           ENDIF
        ENDIF
     ENDDO
300  CONTINUE
  ENDDO

  DEALLOCATE(trace)
  DEALLOCATE(w1)
  DEALLOCATE(psic)
  DEALLOCATE(evcr)
  IF (okvan) CALL deallocate_bec_type ( becp )
  RETURN
END SUBROUTINE find_band_sym_so_tpw
