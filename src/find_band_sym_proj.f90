!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE band_symmetry
  !
  !    this module contains routines to rotate the wavefunctions and to
  !    check which projective or vector representations is contained in the 
  !    representation provided by the degenerate groups of wavefunctions.
  !    
  !
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC find_band_sym_proj, rotate_all_psi_tpw

CONTAINS

SUBROUTINE find_band_sym_proj (ik,evc,et,nsym,s,ftau,d_spin,gk,invs, &
     rap_et,times,ngroup,istart,accuracy,cge)
  !
  !   This subroutine finds the irreducible representations 
  !   (vector or projective) which give the transformation 
  !   properties of the wavefunctions evc.
  !   This routine is used when the group has non symmorphic operations
  !   and the k point is at zone border. It works for both one or two
  !   components spinors.
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : rytoev, tpi
  USE gvect,              ONLY : nl
  USE klist,              ONLY : ngk, igk_k
  USE wvfct,              ONLY : nbnd, npwx
  USE fft_base,           ONLY : dfftp
  USE fft_interfaces,     ONLY : invfft
  USE uspp,               ONLY : vkb, nkb, okvan
  USE noncollin_module,   ONLY : npol
  USE becmod,             ONLY : bec_type, becp, calbec, &
                                 allocate_bec_type, deallocate_bec_type
  USE proj_rap_point_group,   ONLY : which_elem, char_mat_proj, nrap_proj
  USE point_group,        ONLY : find_factor_system
  USE io_global,          ONLY : stdout
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp,                 ONLY : mp_sum

  IMPLICIT NONE

  INTEGER, INTENT(IN) ::   &
       ik,               &   ! the k point
       nsym,             &   ! the number of symmetries
       s(3,3,48),        &   ! the rotation matrices  S_isym
       ftau(3,48),       &   ! the fractional translation 
       invs(48),         &   ! the inverse of each operation
       gk(3,48),         &   ! the gk in crystal coordinates G_S^{-1}_isym
       cge                   ! the extended group code

  INTEGER, INTENT(OUT) ::  &
       ngroup,             & ! the number of groups of degenerate bands
       istart(nbnd+1),     & ! the starting point of each group of bands
       rap_et(nbnd)          ! the representation of each point    

  REAL(DP), INTENT(IN) ::  &
       et(nbnd),           & ! the eigenvalues
       accuracy              ! accepted error in the noninteger part of times

  COMPLEX(DP), INTENT(IN) ::  &
       d_spin(2,2,48),        &   ! the spin rotation matrices
       evc(npwx*npol, nbnd)       ! the wavefunctions

  COMPLEX(DP), INTENT(OUT) ::   &
       times(nbnd,24)       ! the number of times each representation appears
                            ! in each group

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::         &
       ibnd,      &    ! counter on bands
       jbnd,      &    ! counter on bands
       igroup,    &    ! counter on groups of bands
       dim_rap,   &    ! the dimension of the representation
       dim_rap_t, &    ! the dimension of the test representation
       i,j,       &    ! counter on the representation dimensions
       irot,      &    ! counter on symmetries
       jrot,      &    ! counter on symmetries
       irap,      &    ! counte on representations
       shift,     &    ! used when some representation is reducible
       ipol,      &    ! counter on spin components
       nr1, nr2, nr3, & ! fft mesh dimension. ftau is in these units
       npw,          & ! the number of plane waves
       ind1, ind2,   & ! the indeces within the plane waves
       igroup_t,     & ! choice of the representation
       has_e           ! auxiliary, here we do not use -E

  COMPLEX(DP) :: &
       zdotc,    &     ! the scalar product routines     
       factor(48,48), &! the factor system of the representation
       pha             ! the phase to be applied at zone border

  REAL(DP) :: &
       arg,    &       ! the argument of the phase
       ft(3,48)        ! the fractional translations in crystal coordinates

  REAL(DP), ALLOCATABLE ::  &
       w1(:)           ! list of energy eigenvalues in eV

  COMPLEX(DP), ALLOCATABLE ::  &
       psic_nc(:,:,:),    & ! wavefunctions in real space
       evcr(:,:),         & ! the rotated of each wave function in rec. space
       trace(:,:),        & ! the trace of the symmetry matrix
       sym_mat(:,:,:)       ! the matrices of one representation

!
!  Allocate the necessary quantities
!
  ALLOCATE(w1(nbnd))
  ALLOCATE(evcr(npwx*npol,nbnd))
  ALLOCATE(psic_nc(dfftp%nnr,npol,nbnd))
  ALLOCATE(trace(96,nbnd))
  IF (okvan) CALL allocate_bec_type ( nkb, nbnd, becp )
  !
  !  A few initializations
  !
  rap_et=-1
  w1=et*rytoev
  has_e=1
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
!
!  Compute the fractional translations 
!
  nr1=dfftp%nr1
  nr2=dfftp%nr2
  nr3=dfftp%nr3
  DO irot=1,nsym
     ft(1,irot) = FLOAT(ftau(1,irot)) / FLOAT(nr1)
     ft(2,irot) = FLOAT(ftau(2,irot)) / FLOAT(nr2)
     ft(3,irot) = FLOAT(ftau(3,irot)) / FLOAT(nr3)
  ENDDO
!
!  Bring the wavefunctions in real space
!
  npw = ngk(ik)
  psic_nc=(0.0_DP,0.0_DP)
  DO ipol=1, npol
     ind1 = 1 + (ipol-1)*npwx
     ind2 = npw + (ipol-1)*npwx
     DO ibnd=1,nbnd
        psic_nc(nl(igk_k(1:npw,ik)),ipol,ibnd) = evc(ind1:ind2,ibnd)
        CALL invfft ('Dense', psic_nc(:,ipol,ibnd), dfftp)
     ENDDO
  ENDDO

!  igroup_t=1
!  dim_rap_t=istart(igroup_t+1)-istart(igroup_t)
!  ALLOCATE(sym_mat(dim_rap_t,dim_rap_t,nsym))
!
!  Rotate the wavefunction an compute the trace block by block of the
!  matrices that represent the rotation operators
!
  trace=(0.d0,0.d0)
  DO irot=1,nsym
     !
     !   Rotate all the bands together.
     !   NB: rotate_psi assumes that s is in the small group of k. It does not
     !       rotate the k point.
     !   With this call evcr contains the wavefunctions rotated with S_irot
     !   (not invs(irot))
     !
     CALL rotate_all_psi_tpw(ik,psic_nc,evcr,s(1,1,invs(irot)),        &
          ftau(1,invs(irot)),d_spin(1,1,irot),has_e,gk(1,invs(irot)))
!
!  At zone border there is an additional phase to add.
!  ft already is minus the fractionary translation. gk and ft both in crystal
!  coordinates of the reciprocal and direct lattice respectively.
!
     arg = tpi * ( gk(1,invs(irot))*ft(1,irot) +  &
                   gk(2,invs(irot))*ft(2,irot) +  &
                   gk(3,invs(irot))*ft(3,irot) )
     pha=CMPLX(COS(arg), SIN(arg),kind=DP)
     evcr(:,1:nbnd)=evcr(:,1:nbnd)*pha
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
           trace(irot,igroup)=trace(irot,igroup) +            &
                zdotc(npol*npwx,evc(1,ibnd),1,evcr(1,ibnd),1)
        ENDDO
        !      write(6,*) igroup, irot, dim_rap, trace(irot,igroup)
     ENDDO
!
! For testing purposes we can compute all the matrix of the 
! representation for group igroup. Set igroup_t with the group of interest
! The matrices sym_mat are saved already in the correct position
! using the stanrdar order of the point groups
!
!     DO i=1,dim_rap_t
!        ibnd=istart(igroup_t)+i-1
!        DO j=1,dim_rap_t
!           jbnd=istart(igroup_t)+j-1
!           sym_mat(i,j,which_elem(irot))= zdotc(npol*npwx,evc(1,ibnd),1,&
!                                                evcr(1,jbnd),1)
!        ENDDO
!     ENDDO
  ENDDO
  !
  CALL mp_sum(trace,intra_bgrp_comm)
  !
!  CALL mp_sum(sym_mat,intra_bgrp_comm)
!  WRITE(stdout,'(/,5x,"The factor system of the band representation",/)')
!  CALL find_factor_system(sym_mat, dim_rap_t, nsym, cge, factor, .TRUE.)
!  DEALLOCATE(sym_mat)

!  DO igroup=1,ngroup
!     DO irot=1,nsym
!        DO jrot=1,nsym
!           IF (which_elem(jrot)==irot) &
!              WRITE(stdout,'(2i5,2f11.8)') igroup,jrot,trace(jrot,igroup) 
!        ENDDO
!     ENDDO
!     WRITE(stdout,*)
!  ENDDO
!
! Now for each group of band decompose the representation in the irreducible
! representations that it contains
!
  DO igroup=1,ngroup
     dim_rap=istart(igroup+1)-istart(igroup)
     shift=0
     DO irap=1,nrap_proj
        times(igroup,irap)=(0.d0,0.d0)
        DO irot=1,nsym
           times(igroup,irap)=times(igroup,irap) &
                +trace(irot,igroup)*CONJG(char_mat_proj(irap,which_elem(irot)))
        ENDDO
        times(igroup,irap)=times(igroup,irap)/nsym
!       WRITE(6,*) 'group representation ', igroup, irap, times(igroup,irap)

        IF ((ABS(NINT(DBLE(times(igroup,irap)))-DBLE(times(igroup,irap))) &
             > accuracy).OR. (ABS(AIMAG(times(igroup,irap))) > eps) ) THEN
!
!   If times is not an integer number or its imaginary part is not zero
!   this representation is not the correct one. 
!   We set all this group of bands to representation 0 to signal that
!   something went wrong
!
           ibnd=istart(igroup)+shift
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dim_rap
                 ibnd=istart(igroup)+i-1
                 rap_et(ibnd)=0
              ENDDO
              shift = shift + dim_rap
           ENDIF
           GOTO 300
        ELSEIF (ABS(times(igroup,irap)) > accuracy) THEN
!
!   In this case the group of bands belong to the representation irap
!   and we can set rap_et for these bands. shift is used because 
!   a group of bands can be the basis for a reducible representation
!   in case of accidental degeneracy. In this case we do not distinguish
!   which band belongs to each degenerate representation.
!
           ibnd=istart(igroup)+shift
           dim_rap=NINT(DBLE(char_mat_proj(irap,1)))
           IF (rap_et(ibnd)==-1) THEN
              DO i=1,dim_rap*NINT(DBLE(times(igroup,irap)))
                 ibnd=istart(igroup)+shift+i-1
                 rap_et(ibnd)=irap
              ENDDO
              shift = shift + dim_rap*NINT(DBLE(times(igroup,irap)))
           ENDIF
        ENDIF
     ENDDO
300  CONTINUE
  ENDDO

  !WRITE( stdout, '(/,1x,74("*"))')

  DEALLOCATE(trace)
  DEALLOCATE(w1)
  DEALLOCATE(evcr)
  DEALLOCATE(psic_nc)
  IF (okvan) CALL deallocate_bec_type ( becp )

  RETURN
END SUBROUTINE find_band_sym_proj

SUBROUTINE rotate_all_psi_tpw(ik,psic_nc,evcr,s,ftau,d_spin,has_e,gk)
  !
  !  This subroutine rotates a one component or a two component
  !  wavefunction according to the symmetry s (actually it applies s^-1). 
  !  d_spin contains the 2x2 rotation matrix in the spin space (it must be
  !  the d_spin that corresponds to the inverse of s).
  !  has_e=-1 means that also a 360 degrees rotation is applied in spin space.
  !  For one component wavefunctions d_spin and has_e are not used.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym_many, cscatter_sym_many
  USE fft_interfaces, ONLY : fwfft
  USE gvect,     ONLY : nl
  USE wvfct,     ONLY : nbnd, npwx
  USE klist,     ONLY : ngk, igk_k
  USE noncollin_module, ONLY : npol
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ik
  INTEGER, INTENT(IN) :: s(3,3), ftau(3), gk(3), has_e
  COMPLEX(DP) :: psic_nc(dfftp%nnr,npol,nbnd), evcr(npwx*npol,nbnd), &
                 d_spin(2,2)
  COMPLEX(DP), ALLOCATABLE :: psic(:,:), psir(:), evcr_save(:,:)
  COMPLEX(DP), ALLOCATABLE :: phase1(:), phase2(:), phase3(:)
  COMPLEX(DP) :: phase
  REAL(DP) :: arg
  INTEGER :: i, j, k, ri, rj, rk, ir, rir, ipol, jpol, ibnd
  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, npw
  INTEGER :: ind1, ind2, ind3, ind4
  INTEGER :: start_band, last_band, my_nbnd_proc
  INTEGER :: start_band_proc(dfftp%nproc), nbnd_proc(dfftp%nproc)
  !
  COMPLEX (DP), ALLOCATABLE :: psir_collect(:)
  COMPLEX (DP), ALLOCATABLE :: psic_collect(:,:)

  nr1=dfftp%nr1
  nr2=dfftp%nr2
  nr3=dfftp%nr3
  nr1x=dfftp%nr1x
  nr2x=dfftp%nr2x
  nr3x=dfftp%nr3x
  nrxx=dfftp%nnr

  call divide (intra_bgrp_comm, nbnd, start_band, last_band)
  start_band_proc=0
  start_band_proc(dfftp%mype+1)=start_band
  nbnd_proc=0
  my_nbnd_proc=last_band-start_band+1
  nbnd_proc(dfftp%mype+1)=my_nbnd_proc
  CALL mp_sum(start_band_proc, intra_bgrp_comm)
  CALL mp_sum(nbnd_proc, intra_bgrp_comm)

  ALLOCATE (psic_collect(nr1x*nr2x*nr3x,my_nbnd_proc))
  ALLOCATE (psir_collect(nr1x*nr2x*nr3x))
  ALLOCATE(psic(nrxx,nbnd))
  ALLOCATE(psir(nrxx))
  ALLOCATE(evcr_save(npwx*npol,nbnd))
  ALLOCATE(phase1(nr1))
  ALLOCATE(phase2(nr2))
  ALLOCATE(phase3(nr3))
  evcr_save=(0.0_DP,0.0_DP)
  DO i = 1, nr1
     arg = tpi * (gk(1)*(i-1))/DBLE(nr1)
     phase1(i) = CMPLX(COS(arg), SIN(arg), KIND=DP)
  END DO
  DO j = 1, nr2
     arg = tpi * (gk(2)*(j-1))/DBLE(nr2)
     phase2(j) = CMPLX(COS(arg), SIN(arg), KIND=DP)
  END DO
  DO k = 1, nr3
     arg = tpi * (gk(3)*(k-1))/DBLE(nr3)
     phase3(k) = CMPLX(COS(arg), SIN(arg), KIND=DP)
  END DO
  !
  npw = ngk(ik)
  DO ipol=1,npol
     !
     ind1=1 + (ipol-1)*npwx
     ind2=npw + (ipol-1)*npwx
     !
     psic = ( 0.D0, 0.D0 )
     psir = ( 0.D0, 0.D0 )
     !
     DO ibnd=1,nbnd
        psic(:,ibnd) = psic_nc(:,ipol,ibnd)
     ENDDO
     !
     CALL cgather_sym_many( dfftp, psic, psic_collect, nbnd, nbnd_proc, &
                                                       start_band_proc  )
     psir_collect=(0.d0,0.d0)
     DO ibnd=1,my_nbnd_proc
        !
        DO k = 1, nr3
           DO j = 1, nr2
              DO i = 1, nr1
                 CALL ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk )
                 ir=i+(j-1)*nr1x+(k-1)*nr1x*nr2x
                 rir=ri+(rj-1)*nr1x+(rk-1)*nr1x*nr2x
                 phase=phase1(i)*phase2(j)*phase3(k)
                 psir_collect(ir)=psic_collect(rir,ibnd)*phase
              ENDDO
           ENDDO
        ENDDO
        psic_collect(:,ibnd) = psir_collect(:)
        !
     ENDDO
     !
     DO ibnd=1,nbnd
        !
        CALL cscatter_sym_many(dfftp, psic_collect, psir, ibnd, nbnd, &
                               nbnd_proc, start_band_proc)

        CALL fwfft ('Dense', psir, dfftp)
        !
        evcr_save(ind1:ind2,ibnd) = psir(nl(igk_k(1:npw,ik)))
     ENDDO
     !
  ENDDO

  IF (npol==2) THEN
     evcr=(0.d0,0.d0)
     DO ibnd=1,nbnd 
        DO ipol=1,npol
           ind1 = 1 + (ipol-1)*npwx
           ind2 = npw + (ipol-1)*npwx
           DO jpol=1,npol
              ind3 = 1 + (jpol-1)*npwx
              ind4 = npw + (jpol-1)*npwx
              evcr(ind1:ind2,ibnd)=evcr(ind1:ind2,ibnd)+ &
                             d_spin(ipol,jpol)*evcr_save(ind3:ind4,ibnd)
           ENDDO
        ENDDO
     ENDDO
     IF (has_e==-1) evcr=-evcr
  ELSE
     evcr=evcr_save
  ENDIF
  !
  DEALLOCATE(phase1)
  DEALLOCATE(phase2)
  DEALLOCATE(phase3)
  DEALLOCATE(evcr_save)
  DEALLOCATE(psic)
  DEALLOCATE(psir)
  DEALLOCATE(psic_collect)
  DEALLOCATE(psir_collect)

  RETURN
END SUBROUTINE rotate_all_psi_tpw

END MODULE band_symmetry
