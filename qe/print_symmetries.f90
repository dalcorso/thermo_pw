!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE print_symmetries_tpw ( iverbosity, noncolin, domag )
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : dp
  USE io_global,       ONLY : stdout 
  USE symm_base,       ONLY : nsym, nsym_ns, nsym_na, invsym, s, sr, &
                              t_rev, ft, sname
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, &
       which_irr, char_mat, name_rap, name_class, gname, ir_ram, elem_name
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, name_class_so, d_spin, &
       name_class_so1, elem_name_so
  USE rap_point_group_is, ONLY : nsym_is, sr_is, ft_is, d_spin_is, &
       gname_is, sname_is, code_group_is
  USE cell_base,       ONLY : at, ibrav
  USE fft_base, ONLY : dfftp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iverbosity
  LOGICAL, INTENT(IN) :: noncolin, domag
  !
  INTEGER :: nclass_ref   ! The number of classes of the point group
  INTEGER :: isym, ipol, ftau(3,48)
  REAL (dp) :: ft1, ft2, ft3
  !
  !
  IF (nsym <= 1) THEN
     WRITE( stdout, '(/5x,"No symmetry found")')
  ELSE
     ftau(1,1:nsym) = NINT ( ft(1,1:nsym)*dfftp%nr1 )
     ftau(2,1:nsym) = NINT ( ft(2,1:nsym)*dfftp%nr2 )
     ftau(3,1:nsym) = NINT ( ft(3,1:nsym)*dfftp%nr3 )
     IF (invsym) THEN
        IF ( nsym_ns > 0 ) THEN
           WRITE( stdout, '(/5x,i2," Sym. Ops., with inversion, found ", &
                    &  "(",i2," have fractional translation)")' ) nsym, nsym_ns
        ELSE 
           WRITE( stdout, '(/5x,i2," Sym. Ops., with inversion, found")' )&
                         nsym
        END IF
     ELSE
        IF ( nsym_ns > 0 ) THEN
           WRITE( stdout, '(/5x,i2," Sym. Ops. (no inversion) found ",&
                    &  "(",i2," have fractional translation)")' ) nsym, nsym_ns
        ELSE
           WRITE( stdout,'(/5x,i2," Sym. Ops. (no inversion) found")' ) nsym
        END IF
     ENDIF
  ENDIF
  IF ( nsym_na > 0 ) THEN 
      WRITE( stdout, '(10x,"(note: ",i2," additional sym.ops. were found ", &
                   &   "but ignored",/,10x," their fractional translations ",&
                   &   "are incommensurate with FFT grid)",/)') nsym_na
  ELSE
      WRITE( stdout, '(/)' )
  END IF
  IF ( iverbosity > 0 ) THEN
     IF (noncolin) THEN
        WRITE( stdout, '(26x,"s",18x,"frac. trans.",4x,"U(a,b)")')
     ELSE
        WRITE( stdout, '(26x,"s",18x,"frac. trans.")')
     ENDIF
     nsym_is=0
     DO isym = 1, nsym
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
        IF (noncolin) THEN
           IF (domag) THEN
              WRITE(stdout,*) 'Time Reversal ', t_rev(isym)
              IF (t_rev(isym)==0) THEN
                 nsym_is=nsym_is+1
                 sr_is(:,:,nsym_is) = sr(:,:,isym)
                 CALL find_u(sr_is(1,1,nsym_is), d_spin_is(1,1,nsym_is))
                 ft_is(:,nsym_is)=ft(:,isym)
                 sname_is(nsym_is)=sname(isym)
              ENDIF
           ENDIF
           CALL find_u(sr(1,1,isym),d_spin(1,1,isym))
        END IF
        IF (noncolin) THEN

           IF ( ftau(1,isym).NE.0 .OR. ftau(2,isym).NE.0 .OR. &
                ftau(3,isym).NE.0) THEN
              ft1 = at(1,1)*ftau(1,isym)/dfftp%nr1 + &
                    at(1,2)*ftau(2,isym)/dfftp%nr2 + &
                    at(1,3)*ftau(3,isym)/dfftp%nr3
              ft2 = at(2,1)*ftau(1,isym)/dfftp%nr1 + &
                    at(2,2)*ftau(2,isym)/dfftp%nr2 + &
                    at(2,3)*ftau(3,isym)/dfftp%nr3
              ft3 = at(3,1)*ftau(1,isym)/dfftp%nr1 + &
                    at(3,2)*ftau(2,isym)/dfftp%nr2 + &
                    at(3,3)*ftau(3,isym)/dfftp%nr3
              WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i3,4x), &
                &        " )    f =( ",f8.4," )")') &
                isym, (s(1,ipol,isym),ipol=1,3), &
                         DBLE(ftau(1,isym))/DBLE(dfftp%nr1)
              WRITE( stdout, '(17x," (",3(i3,4x), " )       ( ",f8.4," )")') &
                          (s(2,ipol,isym),ipol=1,3), &
                          DBLE(ftau(2,isym))/DBLE(dfftp%nr2)
              WRITE( stdout, '(17x," (",3(i3,4x), " )       ( ",f8.4," )"/)') &
                   (s(3,ipol,isym),ipol=1,3), DBLE(ftau(3,isym))/DBLE(dfftp%nr3)
              WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f7.3, &
                   &        " )    f =( ",f8.4," ) a=",2f8.4)') &
                   isym, (sr(1,ipol,isym),ipol=1,3), ft1, d_spin(1,1,isym)
              WRITE( stdout, '(17x," (",3f7.3, " )       ( ",f8.4," ) &
                         &b=",2f8.4)') &
                   (sr(2,ipol,isym),ipol=1,3), ft2, d_spin(1,2,isym)
              WRITE( stdout, '(17x," (",3f7.3, " )       ( ",f8.4," )"/)') &
                   (sr(3,ipol,isym),ipol=1,3), ft3
           ELSE
              WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = &
                                                 &(",3(i3,4x), " )")') &
                   isym,  (s (1, ipol, isym) , ipol = 1,3)
              WRITE( stdout, '(17x," (",3(i3,4x)," )")')  &
                            (s(2,ipol,isym), ipol=1,3)
              WRITE( stdout, '(17x," (",3(i3,4x)," )"/)') &
                            (s(3,ipol,isym), ipol=1,3)
              WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f7.3," )",&
                            &18x,"a=",2f8.4)') &
                     isym,  (sr (1, ipol,isym) , ipol = 1, 3), d_spin(1,1,isym)
              WRITE( stdout, '(17x," (",3f7.3," )",18x,"b=",2f8.4)')  &
                          (sr (2, ipol,isym) , ipol = 1, 3), d_spin(1,2,isym)
              WRITE( stdout, '(17x," (",3f7.3," )"/)') &
                                  (sr (3, ipol,isym) , ipol = 1, 3)
           END IF
        ELSE
           IF ( ftau(1,isym).NE.0 .OR. ftau(2,isym).NE.0 .OR. &
                ftau(3,isym).NE.0) THEN
              ft1 = at(1,1)*ftau(1,isym)/dfftp%nr1 &
                  + at(1,2)*ftau(2,isym)/dfftp%nr2 + &
                    at(1,3)*ftau(3,isym)/dfftp%nr3
              ft2 = at(2,1)*ftau(1,isym)/dfftp%nr1 + &
                    at(2,2)*ftau(2,isym)/dfftp%nr2 + &
                    at(2,3)*ftau(3,isym)/dfftp%nr3
              ft3 = at(3,1)*ftau(1,isym)/dfftp%nr1 + &
                    at(3,2)*ftau(2,isym)/dfftp%nr2 + &
                    at(3,3)*ftau(3,isym)/dfftp%nr3
              WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i3,4x), &
                &        " )    f =( ",f8.4," )")') &
                isym, (s(1,ipol,isym),ipol=1,3), &
                                    DBLE(ftau(1,isym))/DBLE(dfftp%nr1)
              WRITE( stdout, '(17x," (",3(i3,4x), " )       ( ",f8.4," )")') &
                   (s(2,ipol,isym),ipol=1,3), &
                                    DBLE(ftau(2,isym))/DBLE(dfftp%nr2)
              WRITE( stdout, '(17x," (",3(i3,4x), " )       ( ",f8.4," )"/)')&
                (s(3,ipol,isym),ipol=1,3), &
                                    DBLE(ftau(3,isym))/DBLE(dfftp%nr3)
              WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f7.3, &
                &        " )    f =( ",f8.4," )")') &
                isym, (sr(1,ipol,isym),ipol=1,3), ft1
              WRITE( stdout, '(17x," (",3f7.3, " )       ( ",f8.4," )")') &
                   (sr(2,ipol,isym),ipol=1,3), ft2
              WRITE( stdout, '(17x," (",3f7.3, " )       ( ",f8.4," )"/)') &
                   (sr(3,ipol,isym),ipol=1,3), ft3
           ELSE
              WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = &
                            &(",3(i3,2x), " )")') &
                   isym,  (s (1, ipol, isym) , ipol = 1,3)
              WRITE( stdout, '(17x," (",3(i3,2x)," )")')  &
                                          (s(2,ipol,isym), ipol=1,3)
              WRITE( stdout, '(17x," (",3(i3,2x)," )"/)') &
                                          (s(3,ipol,isym), ipol=1,3)
              WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f7.3," )")') &
                              isym,  (sr (1, ipol,isym) , ipol = 1, 3)
              WRITE( stdout, '(17x," (",3f7.3," )")')  &
                                     (sr (2, ipol,isym) , ipol = 1, 3)
              WRITE( stdout, '(17x," (",3f7.3," )"/)') &
                                     (sr (3, ipol,isym) , ipol = 1, 3)
           END IF
        END IF
     END DO
     !
     CALL find_group(nsym,sr,gname,code_group)
     !
     ! ... Do not attempt calculation of classes if the lattice is provided
     ! ... in input as primitive vectors and not computed from parameters:
     ! ... the resulting vectors may not be accurate enough for the algorithm 
     !
     IF ( ibrav == 0 ) RETURN
     !
     IF (noncolin.AND.domag) THEN
        CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
        CALL set_irr_rap_so(code_group_is,nclass_ref,nrap,char_mat_so, &
             name_rap_so,name_class_so,name_class_so1)
        CALL divide_class_so_tpw(code_group_is,nsym_is,sr_is,d_spin_is, &
             has_e,nclass,nelem_so,elem_so,which_irr_so)
        IF (nclass.ne.nclass_ref) CALL errore('summary', &
             'point double group ?',1)
        CALL set_class_el_name_so(nsym_is,sname_is,has_e,nclass,nelem_so, &
                                  elem_so,elem_name_so)
     ELSE
        IF (noncolin) THEN
           CALL set_irr_rap_so(code_group,nclass_ref,nrap,char_mat_so, &
                name_rap_so,name_class_so,name_class_so1)
           CALL divide_class_so_tpw(code_group,nsym,sr,d_spin,has_e,nclass,  &
                nelem_so, elem_so,which_irr_so)
           IF (nclass.ne.nclass_ref) CALL errore('summary', &
                'point double group ?',1)
           CALL set_class_el_name_so(nsym,sname,has_e,nclass,nelem_so, &
                                     elem_so,elem_name_so)
        ELSE
           CALL set_irr_rap(code_group,nclass_ref,char_mat,name_rap, &
                name_class,ir_ram)
           CALL divide_class_tpw(code_group,nsym,sr,nclass,nelem,elem,which_irr)
           IF (nclass.ne.nclass_ref) CALL errore('summary','point group ?',1)
           CALL set_class_el_name(nsym,sname,nclass,nelem,elem,elem_name)
        ENDIF
     ENDIF
     CALL write_group_info(.true.)
     !
  END IF
  !
END SUBROUTINE print_symmetries_tpw
