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
  USE wvfct,                ONLY : et, nbnd, npwx, g2kin
  USE klist,                ONLY : ngk, igk_k
  USE cell_base,            ONLY : celldm
  USE gvecw,                ONLY : ecutwfc
  USE klist,                ONLY : xk, nks, nkstot
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE buffers,              ONLY : get_buffer
  USE symm_base,            ONLY : s, ftau, nsym, t_rev, sname
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram, elem_name
  USE rap_point_group_so,   ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, &
       name_class_so, d_spin, name_class_so1, elem_name_so
  USE rap_point_group_is,   ONLY : nsym_is, sr_is, ftau_is, gname_is, &
       sname_is, code_group_is
  USE control_2d_bands,     ONLY : nkz, averag, vacuum, aux_ind_sur,  &
                                   sym_divide, nlayers, identify_sur, &
                                   surface1, surface2
  USE control_bands,        ONLY : lsym
  USE control_paths,        ONLY : high_sym_path, disp_nqs
  USE data_files,           ONLY : flprojlayer
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
  INTEGER :: ik, i, j, irot, iclass, ig, ibnd, ilayer, ishift
  INTEGER :: npw, nat_, nbnd_, nkstot_, idum, iun
  INTEGER :: spin_component, nks1, nks2, firstk, lastk
  INTEGER :: nks1tot, nks2tot
  INTEGER :: iunout, igroup, irap, dim_rap, ik2, ike, ikz, nks_, ios
  INTEGER :: sk(3,3,48), ftauk(3,48), gk(3,48), sk_is(3,3,48), &
       gk_is(3,48), t_revk(48), nsymk, isym, ipol, jpol
  LOGICAL :: is_complex, is_complex_so, is_symmorphic, search_sym
  LOGICAL, ALLOCATABLE :: high_symmetry(:)
  LOGICAL :: exst
  REAL(DP), PARAMETER :: accuracy=1.d-4
  COMPLEX(DP) :: d_spink(2,2,48), d_spink1(2,2,48),d_spin_is(2,2,48)
  COMPLEX(DP),ALLOCATABLE :: times(:,:,:)
  REAL(DP) :: dxk(3), dkmod, dkmod_save, k1(3), k2(3), modk1, modk2, ps
  INTEGER, ALLOCATABLE :: rap_et(:,:), code_group_k(:), aux_ind(:)
  INTEGER, ALLOCATABLE :: ngroup(:), istart(:,:)
  REAL(DP) :: factor_dx, sizeb, sizec
  CHARACTER(LEN=11) :: group_name
  CHARACTER(LEN=45) :: snamek(48)
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER (LEN=256) :: filband, namefile
  !
  IF (.NOT.lsym) GOTO 450
  IF (spin_component/=1.and.nspin/=2) &
       CALL errore('sym_band_sub','incorrect spin_component',1)
  IF (spin_component<1.or.spin_component>2) &
       CALL errore('sym_band_sub','incorrect lsda spin_component',1)

  ALLOCATE(rap_et(nbnd,nkstot))
  ALLOCATE(code_group_k(nkstot))
  ALLOCATE(aux_ind(nkstot))
  ALLOCATE(times(nbnd,24,nkstot))
  ALLOCATE(ngroup(nkstot))
  ALLOCATE(istart(nbnd+1,nkstot))
  ALLOCATE(high_symmetry(nkstot))

  code_group_k=0
  rap_et=-1
  times=(0.0_DP,0.0_DP)

  CALL find_nks1nks2(1,nkstot,nks1tot,nks1,nks2tot,nks2,spin_component)

  ios=0
  IF ( ionode ) THEN
     iunout=58
     namefile="band_files/"//TRIM(filband)//".rap"
     IF (nspin==2) &
        namefile="band_files/"//TRIM(filband)//"."// &
                                TRIM(int_to_char(spin_component))//".rap"
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', iostat = ios)
     REWIND (iunout)
  ENDIF

  CALL mp_bcast ( ios, ionode_id, intra_image_comm )
  IF ( ios /= 0) CALL errore ('sym_band_sub', 'Opening filband file', abs (ios) )

  DO ik = nks1, nks2
     !
     !    prepare the indices of this k point
     !
     npw = ngk(ik)
     !
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
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
          sk_is,d_spin_is,gk_is,is_symmorphic,search_sym)
     code_group_k(ik)=code_group
     IF (noncolin.AND.domag) code_group_k(ik)=code_group_is
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
           CALL find_band_sym_so(ik, evc,et(1,ik),nsym_is, &
                sk_is,ftau_is,d_spin_is,gk_is,&
                rap_et(1,ik),times(1,1,ik), &
                ngroup(ik),istart(1,ik),accuracy)
        ELSE
           CALL find_band_sym_so(ik, evc,et(1,ik),nsymk, &
                sk,ftauk,d_spink,gk, &
                rap_et(1,ik),times(1,1,ik),ngroup(ik),&
                istart(1,ik),accuracy)
        ENDIF
     ELSE
        CALL find_band_sym (ik, evc, et(1,ik),  nsymk, sk, ftauk, gk, &
             rap_et(1,ik), times(1,1,ik), ngroup(ik), &
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
     IF (disp_nqs /= nks2tot - nks1tot + 1) &
        CALL errore('sym_band_sub','problem with number of k points',1)
!
!   avoid that cells too long in one direction confuse the
!   plotting program
!
     IF (celldm(2)>0.0_DP) THEN
        sizeb=celldm(2)
     ELSE
        sizeb=1.0_DP
     ENDIF

     IF (celldm(3)>0.0_DP) THEN
        sizec=celldm(3)
     ELSE
        sizec=1.0_DP
     ENDIF
     factor_dx = MAX(5.0_DP, 2.0_DP * sizeb, 2.0_DP * sizec, &
                                       2.0_DP/sizeb, 2.0_DP/sizec)

     high_symmetry(nks1tot:nks2tot)=high_sym_path(1:disp_nqs)
     DO ik=nks1tot, nks2tot
        CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, &
             nsym, sk, ftauk, gk, t_revk, snamek, nsymk)
        CALL find_info_group(nsymk,sk,t_revk,ftauk,d_spink,gk,snamek,&
             sk_is,d_spin_is,gk_is, &
             is_symmorphic,search_sym)
        IF (noncolin) THEN
           CALL set_class_el_name_so(nsymk,snamek,has_e,nclass,nelem_so, &
                                     elem_so,elem_name_so)
        ELSE
           CALL set_class_el_name(nsymk,snamek,nclass,nelem,elem,elem_name)
        ENDIF
        IF (noncolin.AND.domag) THEN
           IF (code_group_k(ik) /= code_group_is) &
             CALL errore('sym_band_sub','problem with code_group',1)
        ELSE
           IF (code_group_k(ik) /= code_group) &
             CALL errore('sym_band_sub','problem with code_group',1)
        ENDIF
        WRITE(stdout, '(/,1x,74("*"))')
        WRITE(stdout, '(/,20x,"xk=(",2(f10.5,","),f10.5,"  )")') &
             xk(1,ik), xk(2,ik), xk(3,ik)
        IF (.not.search_sym) THEN
           WRITE(stdout,'(/,5x,"zone border point and non-symmorphic group ")')
           WRITE(stdout,'(5x,"symmetry decomposition not available")')
           WRITE(stdout, '(/,1x,74("*"))')
        ENDIF
        IF (ik == nks1tot) THEN
           WRITE (iunout, '(" &plot_rap nbnd_rap=",i8,", nks_rap=",i8," /")') &
                nbnd, nks2tot - nks1tot +1
           IF (search_sym) CALL write_group_info(.true.)
           dxk(:) = xk(:,2) - xk(:,1)
           dkmod_save = sqrt( dxk(1)**2 + dxk(2)**2 + dxk(3)**2 )
        ELSE
           IF (code_group_k(ik)/=code_group_k(ik-1).and.search_sym) &
                CALL write_group_info(.true.)
!
!    When the symmetry changes the point must be considered a high
!    symmetry point. If the previous point was also high_symmetry, there
!    are three possibilities. The two points are distant and in this case
!    both of them must be considered high symmetry. The two points coincide
!    and in this case they have both the same symmetry. If they are close only
!    the first point is a high symmetry point. First compute the distance
!
           dxk(:) = xk(:,ik) - xk(:,ik-1)
           dkmod= sqrt( dxk(1)**2 + dxk(2)**2 + dxk(3)**2 )
           IF (dkmod < 1.D-6) THEN
              !
              !   In this case is_high_sym does not change because the point
              !   is the same
              high_symmetry(ik)=high_symmetry(ik-1)
              !
           ELSE IF (dkmod < factor_dx * dkmod_save) THEN
!
!    In this case the two points are considered close
!
              IF (.NOT. high_symmetry(ik-1)) &
                 high_symmetry(ik) = code_group_k(ik) /= code_group_k(ik-1) &
                                  .OR. high_symmetry(ik)
!
!   avoid that a single too small dkmod makes all the other points appear
!   distant. 
!
              dkmod_save= MAX( 0.5_DP * dkmod_save, dkmod)
           ELSE
!
!    Points are distant. They are all high symmetry
!
              high_symmetry(ik)= .TRUE.
           ENDIF
        ENDIF

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

     aux_ind=0
     IF (nks2tot - nks1tot > 1) THEN
        CALL find_aux_ind_xk(xk(1,nks1tot), xk(1,nks1tot+1), aux_ind(nks1tot+1))
        CALL find_aux_ind_xk(xk(1,nks2tot), xk(1,nks2tot-1), aux_ind(nks2tot-1))
     END IF
     DO ik=nks1tot+1,nks2tot-1
        IF (high_symmetry(ik).AND..NOT.high_symmetry(ik+1)) &
           CALL find_aux_ind_xk(xk(1,ik), xk(1,ik+1), aux_ind(ik+1))
        IF (high_symmetry(ik).AND..NOT.high_symmetry(ik-1)) &
           CALL find_aux_ind_xk(xk(1,ik), xk(1,ik-1), aux_ind(ik-1))
     ENDDO

     DO ik=nks1tot,nks2tot
        WRITE (iunout, '(10x,3f10.6,l5,2i5)') xk(1,ik), xk(2,ik), xk(3,ik), &
             high_symmetry(ik), code_group_k(ik), aux_ind(ik)
        WRITE (iunout, '(10i8)') (rap_et(ibnd,ik), ibnd=1,nbnd)
     ENDDO
     CLOSE(iunout)
!
!  In a pbs calculation we need also the aux_ind for the symmetry descent
!  of the different layers
!
     IF (nkz>1.AND.sym_divide) THEN
        nks_ = (nks2tot - nks1tot + 1) / nkz
        ALLOCATE(aux_ind_sur(nks_,nkz))
        IF (MOD(nkz,2)==0) THEN
           ishift=0
        ELSE
           ishift=nks_
        ENDIF
        aux_ind_sur=0
        DO ik = 1, nks_
           ik2 = ik + ishift
           DO ikz = nks1tot, nks1tot + nkz - 1
              ike = ik + nks_ * (ikz - 1)
              CALL find_aux_ind_xk(xk(1,ike),xk(1,ik2),aux_ind_sur(ik,ikz))
           END DO
        END DO
     END IF
  ELSE
     IF (sym_divide.AND.nkz>1) THEN
        nks_ = (nks2tot - nks1tot + 1) / nkz
        ALLOCATE(aux_ind_sur(nks_,nkz))
     ENDIF
  ENDIF

  IF (sym_divide.AND.nkz>1) &
     CALL mp_bcast(aux_ind_sur,ionode_id,intra_image_comm)
  DEALLOCATE(times)
  DEALLOCATE(code_group_k)
  DEALLOCATE(aux_ind)
  DEALLOCATE(rap_et)
  DEALLOCATE(ngroup)
  DEALLOCATE(high_symmetry)
  DEALLOCATE(istart)
  !
450 CONTINUE
  IF (identify_sur) THEN
     IF (ionode) &
        INQUIRE( FILE = TRIM(flprojlayer), EXIST = exst )
     CALL mp_bcast(exst,ionode_id,intra_image_comm)
!
!   the file with the projections is created here if it does not exist,
!   otherwise we assume that it has been already calculated in a previous run
!
     IF (exst) GOTO 500
     ALLOCATE(averag(nat, nspin, nbnd, nkstot))
     ALLOCATE(vacuum(nspin, nbnd, nkstot))
     CALL plan_avg_sub(averag, vacuum, nat, nbnd, nkstot, nlayers, &
                               surface1, surface2)
     IF (ionode) THEN
        iun=39
        OPEN(UNIT=iun,FILE=TRIM(flprojlayer),STATUS='unknown',ERR=400,&
                                                           IOSTAT=ios)
        WRITE(iun, '(5i8)') nat, nlayers, nbnd, nkstot, nspin     
        WRITE(iun, '(4i8)') surface1, surface2    
        DO ik=nks1tot,nks2tot
           DO ibnd=1, nbnd
              WRITE(iun,'(2i8)') ik, ibnd
              DO ilayer=1,nlayers
                WRITE(iun,'(i8,4f17.12)') ilayer, averag(ilayer, 1:nspin, &
                                                  ibnd, ik)
              ENDDO
              WRITE(iun,'(4f20.12)') vacuum(1:nspin,ibnd, ik)
           ENDDO
        ENDDO
        CLOSE(iun)
     ENDIF
400  CALL mp_bcast(ios,ionode_id,intra_image_comm)
     IF (ios /= 0) CALL errore('sym_band_sub','problems with flprojlayer',1)
     DEALLOCATE (vacuum)
     DEALLOCATE (averag)
  END IF
500 CONTINUE
  !
  RETURN
END SUBROUTINE sym_band_sub
!
SUBROUTINE find_aux_ind_xk(xk1, xk2, aux_ind)
USE kinds, ONLY : DP
USE cell_base, ONLY : at, bg
USE symm_base, ONLY : s, ftau, t_rev, sname, nsym
USE noncollin_module, ONLY : noncolin
USE spin_orb,        ONLY : domag
USE rap_point_group, ONLY : code_group
USE rap_point_group_so, ONLY : d_spin
USE rap_point_group_is, ONLY : code_group_is, nsym_is
USE point_group,        ONLY : find_aux_ind_two_groups
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xk1(3), xk2(3)
INTEGER, INTENT(OUT) :: aux_ind
INTEGER :: sk(3,3,48), ftauk(3,48), gk(3,48), sk_is(3,3,48), &
           gk_is(3,48), t_revk(48), nsymk
INTEGER :: sk1(3,3,48), ftauk1(3,48), gk1(3,48), sk1_is(3,3,48), &
           t_revk1(48), nsymk1
INTEGER :: group_a, group_b, nsym_isa, nsym_isb
LOGICAL :: is_symmorphic, search_sym
CHARACTER(len=45) :: snamek(48), snamek1(48)
COMPLEX(DP) :: d_spink(2,2,48), d_spink1(2,2,48), d_spin_is(2,2,48)

CALL smallgk (xk1, at, bg, s, ftau, t_rev, sname, &
              nsym, sk, ftauk, gk, t_revk, snamek, nsymk)
CALL find_info_group(nsymk, sk, t_revk, ftauk, d_spink, gk, snamek,&
             sk_is, d_spin_is, gk_is, is_symmorphic, search_sym)
group_a=code_group
IF (noncolin.AND.domag) THEN
   group_a=code_group_is
   nsym_isa=nsym_is
ENDIF
CALL smallgk (xk2, at, bg, s, ftau, t_rev, sname, &
            nsym, sk1, ftauk1, gk1, t_revk1, snamek1, nsymk1)
CALL find_info_group(nsymk1, sk1, t_revk1, ftauk1, d_spink1, gk1,&
            snamek1, sk1_is, d_spin_is, gk_is,is_symmorphic,search_sym)
group_b=code_group
IF (noncolin.AND.domag) THEN
   group_b=code_group_is
   nsym_isb=nsym_is
ENDIF

IF (group_b /= group_a) THEN
   IF (noncolin.AND.domag) THEN
      CALL find_aux_ind_two_groups(nsym_isa, nsym_isb, sk_is, sk1_is, at, bg, &
                                group_a, group_b, aux_ind)
   ELSE
      CALL find_aux_ind_two_groups(nsymk, nsymk1, sk, sk1, at, bg, &
                                group_a, group_b, aux_ind)
   ENDIF
ELSE
   aux_ind=0
ENDIF

RETURN
END SUBROUTINE find_aux_ind_xk
