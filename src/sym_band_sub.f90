!
! Copyright (C) 2006-2007 Quantum ESPRESSO group
! Copyright (C) 2016 Andrea Dal Corso 
!    (added symmetry analysis at zone border point in nonsymmorphic space
!     groups)
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE sym_band_sub(filband, spin_component)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at, bg
  USE constants,            ONLY : rytoev, pi
  USE fft_base,             ONLY : dfftp
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : et, nbnd, npwx
  USE klist,                ONLY : ngk, igk_k
  USE cell_base,            ONLY : celldm
  USE klist,                ONLY : xk, nks, nkstot
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE buffers,              ONLY : get_buffer
  USE symm_base,            ONLY : s, ftau, nsym, t_rev, sname
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem,  &
       name_rap, name_class, gname, elem_name
  USE rap_point_group_so,   ONLY : nrap, nelem_so, elem_so, has_e, &
       name_rap_so, d_spin, elem_name_so
  USE rap_point_group_is,   ONLY : nsym_is, sr_is, ftau_is, gname_is, &
       sname_is, code_group_is
  USE proj_rap_point_group, ONLY : which_elem, char_mat_proj, nrap_proj, &
                                   name_rap_proj, group_desc, nsym_proj
  USE point_group,          ONLY : find_group_info_ext, find_irr_proj, &
                                   find_projection_type, nsym_group,  &
                                   has_sigma_h
  USE thermo_sym,           ONLY : code_group_save
  USE control_2d_bands,     ONLY : nkz, averag, vacuum, aux_ind_sur,  &
                                   sym_divide, nlayers, identify_sur, &
                                   surface1, surface2
  USE band_symmetry,        ONLY : find_band_sym_proj
  USE lattices,             ONLY : zone_border, same_star
  USE control_bands,        ONLY : lsym
  USE control_paths,        ONLY : high_sym_path, disp_nqs, nrap_plot, &
                                   rap_plot, dkmod_save
  USE data_files,           ONLY : flprojlayer
  USE uspp,                 ONLY : nkb, vkb
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : evc
  USE io_bands,             ONLY : write_representations
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_images,            ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j, irot, iclass, ig, ibnd, ilayer, ishift
  INTEGER :: npw, nat_, nbnd_, nkstot_, iun
  INTEGER :: spin_component, nks1, nks2
  INTEGER :: nks1tot, nks2tot
  INTEGER :: igroup, irap, dim_rap, ik2, ike, ikz, nks_, ios
  INTEGER :: sk(3,3,48), ftauk(3,48), gk(3,48), invsk(48), t_revk(48), &
             nsymk, isym
  INTEGER :: sk_is(3,3,48), gk_is(3,48), invs_is(48)
  INTEGER :: sk_in(3,3,48), ftau_in(3,48), gk_in(3,48), invs_in(48)
  LOGICAL :: is_symmorphic, search_sym, type1
  LOGICAL, ALLOCATABLE :: high_symmetry(:)
  LOGICAL :: exst, lprinted
  REAL(DP), PARAMETER :: accuracy=1.d-4
  COMPLEX(DP) :: d_spink(2,2,48), d_spin_is(2,2,48), d_spin_in(2,2,48)
  COMPLEX(DP), ALLOCATABLE :: times(:,:,:)
  REAL(DP) :: dxk(3), dkmod, k1(3), k2(3), ps, &
              gauge(48), argument(48,48), srk(3,3,48)
  INTEGER, ALLOCATABLE :: rap_et(:,:), code_group_k(:), aux_ind(:), &
                          code_group_ext_k(:), lprojk(:), nrapk(:), &
                          ptypek(:,:), ngroup(:), istart(:,:)
  INTEGER :: code_group1, code_group_ext, code_group_in, nsym_in, ptype(3), &
             nr1, nr2, nr3, sg_number
  INTEGER :: find_free_unit
  LOGICAL :: lwrite
  LOGICAL, ALLOCATABLE :: same_next(:)
  REAL(DP) :: factor_dx, sizeb, sizec, ft_in(3,48), rot(3,3), ur(3,3), &
              s01(3), s02(3)
  REAL(DP), ALLOCATABLE :: gaugek(:,:)
  CHARACTER(LEN=11) :: group_name
  CHARACTER(LEN=45) :: snamek(48), name_out(12)
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
  ALLOCATE(code_group_ext_k(nkstot))
  ALLOCATE(gaugek(48,nkstot))
  ALLOCATE(lprojk(nkstot))
  ALLOCATE(ptypek(3,nkstot))
  ALLOCATE(nrapk(nkstot))
  ALLOCATE(aux_ind(nkstot))
  ALLOCATE(times(nbnd,24,nkstot))
  ALLOCATE(ngroup(nkstot))
  ALLOCATE(istart(nbnd+1,nkstot))
  ALLOCATE(high_symmetry(nkstot))
  ALLOCATE(same_next(nkstot))

  nr1=dfftp%nr1
  nr2=dfftp%nr2
  nr3=dfftp%nr3
  code_group_k=0
  code_group_ext_k=0
  gaugek=0.0_DP
  lprojk=0
  ptypek=1
  IF (noncolin) ptypek(1,:)=-1
  nrapk=0
  rap_et=-1
  times=(0.0_DP,0.0_DP)
  d_spin_in=(0.0_DP,0.0_DP)
!
!   Find which k points must be done by this pool
!
  CALL find_nks1nks2(1,nkstot,nks1tot,nks1,nks2tot,nks2,spin_component)
!
!  Make the symmetry analysis first
!
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
     CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, nsym, &
          sk, ftauk, gk, t_revk, invsk, snamek, nsymk)
     !
     !  character of the irreducible representations
     !
     CALL find_info_group(nsymk,sk,t_revk,ftauk,d_spink,gk,snamek,&
          sk_is,d_spin_is,gk_is,invs_is,is_symmorphic,search_sym)
!
!    select the symmetries to use and set the fractional translations
!    in crystal coordinates
!
     IF (noncolin.AND. domag) THEN
        nsym_in=nsym_is
        code_group_in=code_group_is
        DO isym=1, nsym_in
           sk_in(:,:,isym)=sk_is(:,:,isym)
           d_spin_in(:,:,isym)=d_spin_is(:,:,isym)
           ftau_in(:,isym)=ftau_is(:,isym)
           gk_in(:,isym)=gk_is(:,isym)
           invs_in(isym)=invs_is(isym)
           ft_in(1,isym) = DBLE( ftau_is(1,isym) ) / DBLE(nr1)
           ft_in(2,isym) = DBLE( ftau_is(2,isym) ) / DBLE(nr2)
           ft_in(3,isym) = DBLE( ftau_is(3,isym) ) / DBLE(nr3)
        END DO
        code_group_k(ik)=code_group_is
        nrapk(ik)=nrap
     ELSE
        nsym_in=nsymk
        code_group_in=code_group
        DO isym=1, nsym_in
           sk_in(:,:,isym)=sk(:,:,isym)
           d_spin_in(:,:,isym)=d_spink(:,:,isym)
           ftau_in(:,isym)=ftauk(:,isym)
           gk_in(:,isym)=gk(:,isym)
           invs_in(isym)=invsk(isym)
           ft_in(1,isym) = DBLE( ftauk(1,isym) ) / DBLE(nr1)
           ft_in(2,isym) = DBLE( ftauk(2,isym) ) / DBLE(nr2)
           ft_in(3,isym) = DBLE( ftauk(3,isym) ) / DBLE(nr3)
        END DO
        code_group_k(ik)=code_group
        nrapk(ik)=nclass
        IF (noncolin) nrapk(ik)=nrap
     END IF
     !
     !   find additional information on the point group.
     !
     DO isym=1,nsym_in
        CALL s_axis_to_cart (sk_in(1,1,isym), srk(1,1,isym), at, bg)
     END DO

     CALL find_group_info_ext(nsym_in, srk, code_group1, code_group_ext, &
                                                    which_elem, group_desc)
     code_group_ext_k(ik)=code_group_ext
     !
     !  Find the symmetry of each state
     !
     IF (search_sym) THEN
!
!    This part uses the standard routines of QE
!
        IF (noncolin) THEN
           CALL find_band_sym_so(ik,evc,et(1,ik),nsym_in, &
                   sk_in,ftau_in,d_spin_in,gk_in,invs_in,&
                   rap_et(1,ik),times(1,1,ik), &
                   ngroup(ik),istart(1,ik),accuracy)
        ELSE
           CALL find_band_sym (ik, evc, et(1,ik), nsymk, sk, ftauk, gk, &
                invsk, rap_et(1,ik), times(1,1,ik), ngroup(ik), &
                istart(1,ik),accuracy)
        ENDIF
        IF (zone_border(xk(1,ik),at,bg,-1)) lprojk(ik)=2
           
     ELSE
!
!  This part uses the routines of thermo_pw
!
        lprojk(ik)=1
        !
        !  Set the factor system of the current representation
        !
        CALL set_factor_system(argument, sk_in, ft_in, gk_in, nsym_in, &
                               which_elem,.FALSE.,noncolin,code_group_ext)
        !
        !  find the p-equivalence and the gauge
        !
        CALL find_projection_type(code_group_in, code_group_ext, argument, &
                                  ptype, gauge, .FALSE.)
        gaugek(:,ik)=gauge(:)
        ptypek(:,ik)=ptype(:)
        !
        !  set the irreducible representations corresponding to the
        !  p-equivalence and the gauge
        !
        CALL find_irr_proj(code_group_ext,char_mat_proj,name_rap_proj, &
                      nrap_proj, nsym_in, ptype, gauge, .FALSE.)

        nrapk(ik)=nrap_proj
        !
        !  finally use these representations to classify the band structure
        !
        CALL find_band_sym_proj (ik, evc, et(1,ik),  nsym_in,  &
                         sk_in, ftau_in, d_spin_in, &
                         gk_in, invs_in, rap_et(1,ik), times(1,1,ik), &
                         ngroup(ik), istart(1,ik), accuracy, code_group_ext)
     END IF
  END DO
  !
  !  Only the symmetry of a set of k points is calculated by this
  !  processor within its pool. Here we collect the results into ionode
  !
  CALL ipoolrecover(code_group_k,1,nkstot,nks)
  CALL ipoolrecover(code_group_ext_k,1,nkstot,nks)
  CALL ipoolrecover(lprojk,1,nkstot,nks)
  CALL ipoolrecover(ptypek,3,nkstot,nks)
  CALL ipoolrecover(nrapk,1,nkstot,nks)
  CALL ipoolrecover(rap_et,nbnd,nkstot,nks)
  CALL ipoolrecover(ngroup,1,nkstot,nks)
  CALL ipoolrecover(istart,nbnd+1,nkstot,nks)

  CALL poolrecover(times,2*24*nbnd,nkstot,nks)
  CALL poolrecover(gaugek,48,nkstot,nks)

  ios=0
  IF ( ionode ) THEN
     namefile="band_files/"//TRIM(filband)//".rap"
     IF (nspin==2) &
        namefile="band_files/"//TRIM(filband)//"."// &
                                TRIM(int_to_char(spin_component))//".rap"
  ENDIF

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
     factor_dx = MAX(6.0_DP, 2.0_DP * sizeb, 2.0_DP * sizec, &
                                       2.0_DP/sizeb, 2.0_DP/sizec)

     high_symmetry(nks1tot:nks2tot)=high_sym_path(1:disp_nqs)
!
!   Now repeat the loop over all the k points and write the results on
!   output
!
     DO ik=nks1tot, nks2tot

        CALL smallgk (xk(1,ik), at, bg, s, ftau, t_rev, sname, &
             nsym, sk, ftauk, gk, t_revk, invsk, snamek, nsymk)
        CALL find_info_group(nsymk,sk,t_revk,ftauk,d_spink,gk,snamek,&
             sk_is,d_spin_is,gk_is,invs_is,is_symmorphic,search_sym)

        IF (noncolin.AND.domag) THEN
           nsym_in=nsym_is
           code_group_in=code_group_is
           DO isym=1, nsym_is
              sk_in(:,:,isym)=sk_is(:,:,isym)
              ftau_in(:,isym)=ftau_is(:,isym)
              gk_in(:,isym)=gk_is(:,isym)
              ft_in(1,isym) = DBLE( ftau_is(1,isym) ) / DBLE(nr1)
              ft_in(2,isym) = DBLE( ftau_is(2,isym) ) / DBLE(nr2)
              ft_in(3,isym) = DBLE( ftau_is(3,isym) ) / DBLE(nr3)
           END DO
        ELSE
           nsym_in=nsymk
           code_group_in=code_group
           DO isym=1, nsymk
              sk_in(:,:,isym)=sk(:,:,isym)
              ftau_in(:,isym)=ftauk(:,isym)
              gk_in(:,isym)=gk(:,isym)
              ft_in(1,isym) = DBLE( ftauk(1,isym) ) / DBLE(nr1)
              ft_in(2,isym) = DBLE( ftauk(2,isym) ) / DBLE(nr2)
              ft_in(3,isym) = DBLE( ftauk(3,isym) ) / DBLE(nr3)
           END DO
        ENDIF

        IF (noncolin) THEN
           IF (domag) THEN
              CALL set_class_el_name_so(nsym_is,sname_is,has_e,nclass,&
                                       nelem_so,elem_so,elem_name_so)
           ELSE
              CALL set_class_el_name_so(nsymk,snamek,has_e,nclass,nelem_so, &
                                     elem_so,elem_name_so)
           ENDIF
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

        IF (.NOT.search_sym) THEN
           WRITE(stdout,'(/,5x,"zone border point and non-symmorphic &
                                                         &operations ")')
           lwrite=.FALSE.
           IF (code_group_ext_k(ik)/=code_group_ext_k(ik-1).OR. &
               ptypek(1,ik)/=ptypek(1,ik-1).OR. &
               ptypek(2,ik)/=ptypek(2,ik-1).OR. &
               ptypek(3,ik)/=ptypek(3,ik-1)) lwrite=.TRUE.

           DO isym=1,nsym_in
              CALL s_axis_to_cart (sk_in(1,1,isym), srk(1,1,isym), at, bg)
           END DO

           CALL find_group_info_ext(nsym_in, srk, code_group1, code_group_ext, &
                                                    which_elem, group_desc)

           CALL set_factor_system(argument, sk_in, ft_in, gk_in, nsym_in, &
                           which_elem,lwrite,noncolin,code_group_ext_k(ik))

           CALL find_projection_type(code_group_in, code_group_ext_k(ik), &
                                      argument, ptype, gauge, .FALSE.)

           CALL find_irr_proj(code_group_ext_k(ik), char_mat_proj, &
                       name_rap_proj, nrap_proj, nsym_in, ptype, gauge, lwrite)
        ENDIF
        IF (ik == nks1tot) THEN
           IF (search_sym) CALL write_group_info(.true.)
           dxk(:) = xk(:,2) - xk(:,1)
        ELSE
           IF (code_group_k(ik)/=code_group_k(ik-1).AND.search_sym) &
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
           ELSE
!
!    Points are distant. They are all high symmetry
!
              high_symmetry(ik)= .TRUE.
           ENDIF
        ENDIF

        IF (noncolin) THEN
           IF (domag) THEN
              WRITE(stdout,'(/,5x,"Band symmetry, ",a11," [",a11, &
                   & "] magnetic double point group,")') gname, gname_is
              IF (ptypek(1,ik)==1) &
                 WRITE(stdout,'(/,5x,"Switching to the point group")')
              WRITE(stdout,'(5x,"using ",a11,/)') gname_is
           ELSE
              IF (ptypek(1,ik)==-1) THEN
                  WRITE(stdout,'(/,5x,"Band symmetry, ",a11,&
                         & " double point group:",/)') gname
              ELSE
                 WRITE(stdout,'(/,5x,"Switching to the point group ")')
                 WRITE(stdout,'(/,5x,"Band symmetry, ",a11,&
                        &" point group:",/)') group_name(code_group_k(ik))
              ENDIF
           ENDIF
        ELSE
           IF (ptypek(1,ik)==1) THEN
              WRITE(stdout,'(/,5x,"Band symmetry, ",a11," point group:",/)') &
                    group_name(code_group_k(ik))
           ELSE
              WRITE(stdout,'(/,5x,"Switching to the double point group")')
              WRITE(stdout,'(5x,"Band symmetry, ",a11,&
                      &" double point group:",/)') group_name(code_group_k(ik))
           ENDIF
        ENDIF

        DO igroup=1,ngroup(ik)
           dim_rap=istart(igroup+1,ik)-istart(igroup,ik)
           IF (search_sym) THEN
              IF (noncolin) THEN
                 DO irap=1,nrapk(ik)
                    name_out(irap)=name_rap_so(irap)
                 ENDDO
              ELSE
                 DO irap=1,nrapk(ik)
                    name_out(irap)=name_rap(irap)
                 ENDDO
              END IF
           ELSE
              DO irap=1,nrapk(ik)
                 name_out(irap)=name_rap_proj(irap)
              ENDDO
           ENDIF

           DO irap=1,nrapk(ik)
              IF ((ABS(NINT(DBLE(times(igroup,irap,ik)))-  &
                   DBLE(times(igroup,irap,ik))) > accuracy).or. &
                  (ABS(AIMAG(times(igroup,irap,ik))) > accuracy)) THEN
                       WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",f12.5,2x,&
                            &"eV",3x,i3,3x, "-->   ?")') &
                            istart(igroup,ik), istart(igroup+1,ik)-1, &
                            et(istart(igroup,ik),ik)*rytoev, dim_rap
                 EXIT
              ELSEIF (ABS(times(igroup,irap,ik)) > accuracy) THEN
                 IF (ABS(NINT(DBLE(times(igroup,irap,ik))-1.d0)) < accuracy) &
                                                       THEN
                    WRITE(stdout,'(5x, "e(",i3," -",i3,") = ",&
                          &f12.5,2x,"eV",3x,i3,3x,"--> ",a31)') &
                           istart(igroup,ik), istart(igroup+1,ik)-1, &
                           et(istart(igroup,ik),ik)*rytoev, &
                           dim_rap, name_out(irap)
                 ELSE
                    WRITE(stdout,'(5x,"e(",i3," -",i3,") = ",&
                          &f12.5,2x,"eV",3x,i3,3x,"--> ",i3," ",a28)') &
                           istart(igroup,ik), istart(igroup+1,ik)-1, &
                           et(istart(igroup,ik),ik)*rytoev, dim_rap, &
                           NINT(DBLE(times(igroup,irap,ik))), &
                           name_out(irap)
                 END IF
              END IF
           END DO
        END DO
        WRITE( stdout, '(/,1x,74("*"))')
     ENDDO
!
!    Here we find which subgroup we have when the symmetry of the point group 
!    changes. 
!
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
!
!    Here write the symmetry information on the output file.
!
     DO ik=nks1tot,nks2tot
        IF (ik==nks2tot) THEN
           same_next(ik)=.FALSE.
        ELSE
           same_next(ik)= same_star(nsym, s, xk(1,ik), xk(1,ik+1), at) 
        ENDIF
     ENDDO
!  In a pbs calculation we need also the aux_ind for the symmetry descent
!  of the different layers
!
     IF (nkz>1.AND.sym_divide) THEN
        nks_ = (nks2tot - nks1tot + 1) / nkz
        ALLOCATE(aux_ind_sur(nks_,nkz))
        type1=has_sigma_h(code_group_save)
        IF (type1) THEN
           ishift=nks_
        ELSE
           ishift=0
        ENDIF
!        IF (MOD(nkz,2)==0) THEN
!           ishift=0
!        ELSE
!           ishift=nks_
!        ENDIF
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

  CALL write_representations(nkstot, nbnd, xk, rap_et, high_symmetry,      &
                        code_group_k, aux_ind, code_group_ext_k, ptypek, &
                        lprojk, same_next, gaugek, namefile)

  IF (sym_divide.AND.nkz>1) &
     CALL mp_bcast(aux_ind_sur,ionode_id,intra_image_comm)

  DEALLOCATE(times)
  DEALLOCATE(code_group_k)
  DEALLOCATE(code_group_ext_k)
  DEALLOCATE(gaugek)
  DEALLOCATE(lprojk)
  DEALLOCATE(ptypek)
  DEALLOCATE(nrapk)
  DEALLOCATE(aux_ind)
  DEALLOCATE(rap_et)
  DEALLOCATE(ngroup)
  DEALLOCATE(high_symmetry)
  DEALLOCATE(same_next)
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
        iun=find_free_unit()
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
           gk_is(3,48), t_revk(48), invsk(48), invs_is(48), nsymk
INTEGER :: sk1(3,3,48), ftauk1(3,48), gk1(3,48), sk1_is(3,3,48), &
           t_revk1(48), invsk1(48), nsymk1
INTEGER :: group_a, group_b, nsym_isa, nsym_isb
LOGICAL :: is_symmorphic, search_sym
CHARACTER(len=45) :: snamek(48), snamek1(48)
COMPLEX(DP) :: d_spink(2,2,48), d_spink1(2,2,48), d_spin_is(2,2,48)

CALL smallgk (xk1, at, bg, s, ftau, t_rev, sname, &
              nsym, sk, ftauk, gk, t_revk, invsk, snamek, nsymk)
CALL find_info_group(nsymk, sk, t_revk, ftauk, d_spink, gk, snamek,&
             sk_is, d_spin_is, gk_is, invs_is, is_symmorphic, search_sym)
group_a=code_group
IF (noncolin.AND.domag) THEN
   group_a=code_group_is
   nsym_isa=nsym_is
ENDIF
CALL smallgk (xk2, at, bg, s, ftau, t_rev, sname, &
            nsym, sk1, ftauk1, gk1, t_revk1, invsk1, snamek1, nsymk1)
CALL find_info_group(nsymk1, sk1, t_revk1, ftauk1, d_spink1, gk1, &
            snamek1, sk1_is, d_spin_is, gk_is, invs_is, &
            is_symmorphic,search_sym)
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

SUBROUTINE set_factor_system(argument, s, ft_in, gk, nsym, which_elem, &
                             verbose, lso, cge)
!
!  This routine sets the factor system of a given projective representation
!  of the small point group of a k vector that depends on the fractional
!  translations and on the gk vectors. If lso is true the routine multiplies the
!  phase factors by the +-1 signs that give the projective representation
!  induced by the double point group (to be used to classify two-component
!  spinor wavefunctions).
!
USE kinds, ONLY : DP
USE constants, ONLY : tpi, pi
USE proj_rap_point_group, ONLY : group_desc
USE point_group, ONLY : sym_label, find_double_product_table
USE io_global, ONLY : stdout  

IMPLICIT NONE
INTEGER, INTENT(IN) :: nsym, cge
INTEGER, INTENT(IN) :: s(3,3,48), gk(3,48), which_elem(48)
REAL(DP), INTENT(IN) :: ft_in(3,48)
LOGICAL, INTENT(IN) :: verbose, lso
REAL(DP), INTENT(OUT) :: argument(48,48)
 
REAL(DP) :: arg

INTEGER :: isym, jsym, gk_inp(3,48)
INTEGER :: epos(48,48), prd(48,48)
REAL(DP) :: ft(3,48)

DO isym = 1, nsym
   DO jsym = 1, nsym
      IF (which_elem(jsym) == isym) THEN 
         ft(:,isym) = ft_in(:,jsym)
         gk_inp(:,isym)=gk(:,jsym)
      END IF
   END DO
END DO

CALL find_double_product_table(prd, epos, cge)
IF (.NOT.lso) epos=1

IF (verbose) THEN
   WRITE(stdout,'(/,3x, "Symmetry",12x,"Fractional translation",12x,&
                              &"G_k (S^-1 k=k+G_k) ")')
   DO isym=1,nsym
      WRITE(stdout,'(a8,3x,3f13.4,5x,3i4)') &
                  TRIM(sym_label(group_desc(isym))), &
                  ft(1,isym), ft(2,isym), ft(3,isym), gk_inp(:,isym)
   END DO
   WRITE(stdout,*)
ENDIF

DO isym=1,nsym
   DO jsym=1,nsym
!
!  note that ft contains minus the fractionary translation and gki(:,isym) the
!  the G that correspond to S^-1 k = k + G where S is the symmetry isym
!
      arg = tpi * ( DBLE(gk_inp(1,isym)) * ft(1,jsym) + &
                    DBLE(gk_inp(2,isym)) * ft(2,jsym) + &
                    DBLE(gk_inp(3,isym)) * ft(3,jsym) )
      argument(isym,jsym)=arg
      IF (epos(isym,jsym)==-1) argument(isym,jsym)=argument(isym,jsym)+pi
   ENDDO
ENDDO

RETURN
END SUBROUTINE set_factor_system
