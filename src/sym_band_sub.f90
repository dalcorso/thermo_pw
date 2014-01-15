!
! Copyright (C) 2006-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE sym_band_sub(filband, spin_component)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba2, at, bg, ibrav
  USE constants,            ONLY : rytoev
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, nl, g
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : et, nbnd, npwx, npw, igk, g2kin, ecutwfc
  USE klist,                ONLY : xk, nks, nkstot
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE buffers,              ONLY : get_buffer
  USE symm_base,            ONLY : s, ftau, nsym, t_rev, sname
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so,   ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, &
       name_class_so, d_spin, name_class_so1
  USE rap_point_group_is,   ONLY : nsym_is, sr_is, ftau_is, gname_is, &
       sname_is, code_group_is
  USE uspp,                 ONLY : nkb, vkb
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_images,             ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j, irot, iclass, ig, ibnd
  INTEGER :: spin_component, nks1, nks2, firstk, lastk
  INTEGER :: iunout, igroup, irap, dim_rap, ios
  INTEGER :: sk(3,3,48), ftauk(3,48), gk(3,48), sk_is(3,3,48), &
       gk_is(3,48), t_revk(48), nsymk, isym, ipol, jpol
  LOGICAL :: is_complex, is_complex_so, is_symmorphic, is_high_sym, search_sym
  REAL(DP), PARAMETER :: accuracy=1.d-4
  COMPLEX(DP) :: d_spink(2,2,48), d_spin_is(2,2,48), zdotc
  COMPLEX(DP),ALLOCATABLE :: times(:,:,:)
  REAL(DP) :: dxk(3), dkmod, dkmod_save
  INTEGER, ALLOCATABLE :: rap_et(:,:), code_group_k(:)
  INTEGER, ALLOCATABLE :: ngroup(:), istart(:,:)
  CHARACTER(len=11) :: group_name
  CHARACTER(len=45) :: snamek(48)
  CHARACTER (len=256) :: filband, namefile
  !
  IF (spin_component/=1.and.nspin/=2) &
       CALL errore('sym_band_sub','incorrect spin_component',1)
  IF (spin_component<1.or.spin_component>2) &
       CALL errore('sym_band_sub','incorrect lsda spin_component',1)

  ALLOCATE(rap_et(nbnd,nkstot))
  ALLOCATE(code_group_k(nkstot))
  ALLOCATE(times(nbnd,24,nkstot))
  ALLOCATE(ngroup(nkstot))
  ALLOCATE(istart(nbnd+1,nkstot))

  code_group_k=0
  rap_et=-1
  times=(0.0_DP,0.0_DP)

  ios=0
  IF ( ionode ) THEN
     iunout=58
     namefile=trim(filband)//".rap"
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', iostat = ios)
     REWIND (iunout)
  ENDIF

  CALL mp_bcast ( ios, ionode_id, intra_image_comm )
  IF ( ios /= 0) CALL errore ('sym_band_suc', 'Opening filband file', abs (ios) )

  DO ik = 1, nks
     !
     !    prepare the indices of this k point
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     CALL init_us_2 (npw, igk, xk (1, ik), vkb)
     !
     !   read eigenfunctions
     !
     CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
     !
     ! Find the small group of k
     !
     CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, nsym, sk, ftauk, &
          gk, t_revk, snamek, nsymk)
     !
     !  character of the irreducible representations
     !
     CALL find_info_group(nsymk,sk,t_revk,ftauk,d_spink,gk,snamek,&
          sk_is,d_spin_is,gk_is, &
          is_symmorphic,search_sym)
     code_group_k(ik)=code_group
     !
     IF (.not.search_sym) THEN
        rap_et(:,ik)=-1
        GOTO 100
     ENDIF
     !
     !  Find the symmetry of each state
     !
     IF (noncolin) THEN
        IF (domag) THEN
           CALL find_band_sym_so(evc,et(1,ik),at,nbnd,npw,nsym_is, &
                ngm,sk_is,ftau_is,d_spin_is,gk_is,xk(1,ik),igk,nl,dfftp%nr1,dfftp%nr2,&
                dfftp%nr3,dfftp%nr1x,dfftp%nr2x,dfftp%nr3x,dfftp%nnr,npwx,rap_et(1,ik),times(1,1,ik), &
                ngroup(ik),istart(1,ik),accuracy)
        ELSE
           CALL find_band_sym_so(evc,et(1,ik),at,nbnd,npw,nsymk,ngm, &
                sk,ftauk,d_spink,gk,xk(1,ik),igk,nl,dfftp%nr1,dfftp%nr2,dfftp%nr3,dfftp%nr1x, &
                dfftp%nr2x,dfftp%nr3x,dfftp%nnr,npwx,rap_et(1,ik),times(1,1,ik),ngroup(ik),&
                istart(1,ik),accuracy)
        ENDIF
     ELSE
        CALL find_band_sym (evc, et(1,ik), at, nbnd, npw, nsymk, ngm, &
             sk, ftauk, gk, xk(1,ik), igk, nl, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, &
             dfftp%nr2x, dfftp%nr3x, dfftp%nnr, npwx, rap_et(1,ik), times(1,1,ik), ngroup(ik),&
             istart(1,ik),accuracy)
     ENDIF

100  CONTINUE
  ENDDO

#ifdef __MPI
  !
  !  Only the symmetry of a set of k points is calculated by this
  !  processor with pool. Here we collect the results into ionode
  !
  CALL ipoolrecover(code_group_k,1,nkstot,nks)
  CALL ipoolrecover(rap_et,nbnd,nkstot,nks)
  CALL poolrecover(times,2*24*nbnd,nkstot,nks)
  CALL ipoolrecover(ngroup,1,nkstot,nks)
  CALL ipoolrecover(istart,nbnd+1,nkstot,nks)
#endif
  IF (ionode) THEN
     is_high_sym=.false.
     DO ik=1, nkstot
        CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, &
             nsym, sk, ftauk, gk, t_revk, snamek, nsymk)
        CALL find_info_group(nsymk,sk,t_revk,ftauk,d_spink,gk,snamek,&
             sk_is,d_spin_is,gk_is, &
             is_symmorphic,search_sym)
        IF (code_group_k(ik) /= code_group) &
             CALL errore('sym_band','problem with code_group',1)
        WRITE(stdout, '(/,1x,74("*"))')
        WRITE(stdout, '(/,20x,"xk=(",2(f10.5,","),f10.5,"  )")') &
             xk(1,ik), xk(2,ik), xk(3,ik)
        IF (.not.search_sym) THEN
           WRITE(stdout,'(/,5x,"zone border point and non-symmorphic group ")')
           WRITE(stdout,'(5x,"symmetry decomposition not available")')
           WRITE(stdout, '(/,1x,74("*"))')
        ENDIF
        IF (ik == 1) THEN
           WRITE (iunout, '(" &plot_rap nbnd_rap=",i4,", nks_rap=",i4," /")') &
                nbnd, nkstot
           IF (search_sym) CALL write_group_info(.true.)
           is_high_sym=.true.
           dxk(:) = xk(:,2) - xk(:,1)
           dkmod_save = sqrt( dxk(1)**2 + dxk(2)**2 + dxk(3)**2 )
        ELSE
           IF (code_group_k(ik)/=code_group_k(ik-1).and.search_sym) &
                CALL write_group_info(.true.)
!
!    When the symmetry changes the point must be considered a high
!    symmetry point. If the previous point was also high_symmetry, there
!    are two possibilities. The two points are distant and in this case
!    both of them must be considered high symmetry. If they are close only
!    the first point is a high symmetry point. First compute the distance
!
           dxk(:) = xk(:,ik) - xk(:,ik-1)
           dkmod= sqrt( dxk(1)**2 + dxk(2)**2 + dxk(3)**2 )
           IF (dkmod < 5.0_DP * dkmod_save) THEN
!
!    In this case the two points are considered close
!
              is_high_sym= ((code_group_k(ik)/=code_group_k(ik-1)) &
                 .and..not.is_high_sym) 
              IF (dkmod > 1.d-3) dkmod_save=dkmod
           ELSE
!
!    Points are distant. They are all high symmetry
!
              is_high_sym= .TRUE.
           ENDIF
        ENDIF
        WRITE (iunout, '(10x,3f10.6,l5)') xk(1,ik),xk(2,ik),xk(3,ik), &
             is_high_sym
        WRITE (iunout, '(10i8)') (rap_et(ibnd,ik), ibnd=1,nbnd)
        IF (.not.search_sym) CYCLE
        IF (noncolin) THEN
           IF (domag) THEN
              WRITE(stdout,'(/,5x,"Band symmetry, ",a11," [",a11, &
                   & "] magnetic double point group,")') gname, gname_is
              WRITE(stdout,'(5x,"using ",a11,/)') gname_is
           ELSE
              WRITE(stdout,'(/,5x,"Band symmetry, ",a11,&
                   & " double point group:",/)') gname
           ENDIF
        ELSE
           WRITE(stdout,'(/,5x,"Band symmetry, ",a11," point group:",/)') &
                group_name(code_group_k(ik))
        ENDIF

        DO igroup=1,ngroup(ik)
           dim_rap=istart(igroup+1,ik)-istart(igroup,ik)
           DO irap=1,nclass
              IF (noncolin) THEN
                 IF ((abs(nint(dble(times(igroup,irap,ik)))-  &
                      dble(times(igroup,irap,ik))) > accuracy).or. &
                      (abs(aimag(times(igroup,irap,ik))) > accuracy) ) THEN
                    WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,&
                         &"eV",3x,i3,3x, "-->   ?")') &
                         istart(igroup,ik), istart(igroup+1,ik)-1, &
                         et(istart(igroup,ik),ik)*rytoev, dim_rap
                    exit
                 ELSEIF (abs(times(igroup,irap,ik)) > accuracy) THEN
                    IF (abs(nint(dble(times(igroup,irap,ik))-1.d0)) < &
                         accuracy) THEN
                       WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, &
                            dim_rap, name_rap_so(irap)
                    ELSE
                       WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",i3," ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, dim_rap, &
                            nint(dble(times(igroup,irap,ik))), name_rap_so(irap)
                    ENDIF
                 ENDIF
              ELSE
                 IF ((abs(nint(dble(times(igroup,irap,ik)))-  &
                      dble(times(igroup,irap,ik))) > accuracy).or. &
                      (abs(aimag(times(igroup,irap,ik))) > accuracy) ) THEN
                    WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,&
                         &"eV",3x,i3,3x, "-->   ?")') &
                         istart(igroup,ik), istart(igroup+1,ik)-1, &
                         et(istart(igroup,ik),ik)*rytoev, dim_rap
                    exit
                 ELSEIF (abs(times(igroup,irap,ik)) > accuracy) THEN
                    IF (abs(nint(dble(times(igroup,irap,ik))-1.d0)) < &
                         accuracy) THEN
                       WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, &
                            dim_rap, name_rap(irap)
                    ELSE
                       WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",&
                            &f12.5,2x,"eV",3x,i3,3x,"--> ",i3," ",a15)') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, dim_rap, &
                            nint(dble(times(igroup,irap,ik))), name_rap(irap)
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
        WRITE( stdout, '(/,1x,74("*"))')
     ENDDO
     CLOSE(iunout)
  ENDIF

  !
  DEALLOCATE(times)
  DEALLOCATE(code_group_k)
  DEALLOCATE(rap_et)
  DEALLOCATE(ngroup)
  DEALLOCATE(istart)
  !
  RETURN
END SUBROUTINE sym_band_sub
!
