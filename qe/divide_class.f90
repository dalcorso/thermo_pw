!
! Copyright (C) 2006-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------------
SUBROUTINE divide_class_tpw( code_group, nrot, smat, nclass, nelem, elem, which_irr )
!-----------------------------------------------------------------------------
!! This subroutine receives as input a set of \(\text{nrot}\) 3x3 matrices 
!! \(\text{smat}\), which are assumed to be the operations of the point group given
!! by \(\texttt{code_group}\). The matrices \(\text{smat}\) are in cartesian 
!! coordinates.  
!! This routine divides the group in classes and finds:
!
!! * \(\text{nclass}\): the number of classes of the group;
!! * \(\text{nelem(iclass)}\): for each class, the number of elements of the class
!! * \(\text{elem(i,iclass)}\): \(1<i<\text{nelem}(\text{iclass})\) for each class tells
!!                              which matrices \(\text{smat}\) belong to that class;
!! * \(\text{which_irr(iclass)}\): for each class gives the position of that class in the
!!                                 character table associated with the group and provided
!!                                 by the routine \(\texttt{set_irr_rap}\).
!! NB: changing the order of the elements in the character table must be reflected in 
!! a change to \(\text{which_irr}\). Presently the character tables are those given 
!! by P.W. Atkins, M.S. Child, and C.S.G. Phillips, "Tables for group theory".  
!! Several equivalent names for the irreducible representation are given. D, G, L, S 
!! are used for Delta, Gamma, Lambda and Sigma.
!
USE kinds, ONLY : DP
USE point_group, ONLY : angle_rot_s_tpw
!
IMPLICIT NONE
!
INTEGER :: code_group
!! The code of the point group
INTEGER :: nrot
!! The number of symmetry operation
INTEGER :: nclass
!! The number of classes
INTEGER :: nelem(12)
!! The elements of each class 
INTEGER :: elem(8,12)
!! Which elements in the smat list for each class
INTEGER :: which_irr(12)
!! See the subroutine description 
REAL(DP) :: smat(3,3,nrot)
!! The set of nrot 3x3 matrices
!
! ... local variables
!
REAL(DP) :: cmat(3,3), ax(3), bx(3), cx(3), ars
!
INTEGER :: done(48), irot, jrot, krot, iclass, i, other, other1
INTEGER :: tipo_sym, ipol, axis, axis1, axis2, ts, iax, ibx, icx, aclass, &
           bclass, cclass, imax, imbx, imcx, amclass, bmclass, cmclass, ind2(3)
REAL(DP), PARAMETER :: eps = 1.d-7
REAL(DP) :: angle_rot, angle_rot_s, angle_vectors, ax_save(3,2:4)
LOGICAL :: compare_mat, is_axis, is_parallel, isok, isok1, done_ax(6)
!
! Divide the group in classes.
!
nclass=0
nelem=0
done=0
DO irot=1,nrot
   IF (done(irot)==0) THEN
      nclass=nclass+1
      DO jrot=1,nrot
         CALL coniug_mat(smat(1,1,jrot),smat(1,1,irot),cmat)
         DO krot=1,nrot
            IF (compare_mat(cmat,smat(1,1,krot)).AND.done(krot)==0) THEN
               nelem(nclass)=nelem(nclass)+1
               elem(nelem(nclass),nclass)=krot
               done(krot)=1 
            ENDIF 
         ENDDO
      ENDDO
   ENDIF
ENDDO
!
!  For each class we should now decide which_irr. This depends on the group
!  and on the tables of characters of the irreducible representations,
!  so we must make different things for different groups.
!
which_irr(1)=1
IF (code_group==1) THEN
   IF (nclass /= 1) CALL errore('divide_class','Wrong classes for C_1',1)
   !
   !  C_1 
   !
ELSEIF (code_group==2.OR.code_group==3.OR.code_group==4) THEN
   !
   !  C_i, C_s, C_2
   !
   IF (nclass /= 2) &
        CALL errore('divide_class','Wrong classes for C_i, C_s or C_2',1)
   which_irr(2)=2
ELSEIF (code_group==5) THEN
   !
   !  C_3
   !
   ! The function angle_rot(smat) provides the rotation angle of the matrix smat
   !
   IF (nclass /= 3) CALL errore('divide_class','Wrong classes for C_3',1)
   DO iclass=2,nclass
      IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
         which_irr(iclass)=2
      ELSE
         which_irr(iclass)=3
      ENDIF
   ENDDO
   !
ELSEIF (code_group==6) THEN
   !
   !  C_4
   !
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for C_4',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-90.d0)<eps) THEN
            which_irr(iclass)=2
         ELSEIF (ABS(ars-270.d0)<eps) THEN
            which_irr(iclass)=4
         ELSE
            CALL errore('divide_class','wrong angle',1)
         ENDIF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==7) THEN
   !
   !  C_6
   !
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for C_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=4
      ELSEIF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=2
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=3
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=5
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=6
         ELSE
            CALL errore('divide_class','wrong angle',1)
         ENDIF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
   !
ELSEIF (code_group==8) THEN
   !
   !  D_2  
   !
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for D_2',1)
   !
   CALL versor(smat(1,1,elem(1,2)),ax)
   CALL which_c2(ax,iax)
   CALL versor(smat(1,1,elem(1,3)),bx)
   CALL which_c2(bx,ibx)
   CALL versor(smat(1,1,elem(1,4)),cx)
   CALL which_c2(cx,icx)
   !
   CALL is_d2(iax, ibx, icx, ind2)
   !
   which_irr(2)=ind2(1) + 1
   which_irr(3)=ind2(2) + 1
   which_irr(4)=ind2(3) + 1
   !
ELSEIF (code_group==9) THEN
   !
   !  D_3
   !
   IF (nclass /= 3) CALL errore('divide_class','Wrong classes for D_3',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==3) THEN
         which_irr(iclass)=2
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==10) THEN
   !
   !  D_4
   !
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for D_4',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         axis=0
         DO ipol=1,3
            IF (is_axis(ax,ipol)) axis=ipol
         ENDDO 
         axis1=MOD(ipol,3)+1
         axis2=MOD(ipol+1,3)+1
         IF (axis==0) call errore('divide_class','unknown D_4 axis ',1)
      ENDIF
   END DO
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=3
         ELSEIF (is_axis(ax,axis1).or.is_axis(ax,axis2)) THEN
            which_irr(iclass)=4
         ELSE
            which_irr(iclass)=5
         END IF
      ELSEIF (ts.ne.3) THEN
         CALL errore('divide_class','wrong sym_type',1)
      END IF
   END DO
ELSEIF (code_group==11) THEN
   !
   !  D_6
   !
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for D_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps) ) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         ENDIF
      ELSEIF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,3)) THEN
            which_irr(iclass)=4
         ELSE 
            CALL which_c2(ax, iax)
            IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
               which_irr(iclass)=5
            ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
               which_irr(iclass)=6
            ELSE
               CALL errore('divide_class','D_6 C_2 axis not recognized',1)
            END IF
         END IF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      END IF
   END DO
ELSEIF (code_group==12) THEN
   !
   !  C_2v
   !
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for C_2v',1)
   iax=0
   ibx=0
   icx=0
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2( ax, iax) 
         which_irr(iclass)=2
      ELSEIF (ts==5) THEN
         IF (ibx==0) THEN
            CALL mirror_axis(smat(1,1,elem(1,iclass)), bx)
            CALL which_c2( bx, ibx) 
            bclass=iclass
         ELSE
            CALL mirror_axis(smat(1,1,elem(1,iclass)), bx)
            CALL which_c2( bx, icx) 
            cclass=iclass
         ENDIF
      ENDIF
   ENDDO
   CALL is_c2v(iax, ibx, icx, isok)
   IF (isok) THEN
      which_irr(bclass)=3
      which_irr(cclass)=4
   ELSE
      CALL is_c2v(iax, icx, ibx, isok1)
      IF (.NOT.isok1) CALL errore('divide_class','problem with C_2v',1)
      which_irr(bclass)=4
      which_irr(cclass)=3
   ENDIF
ELSEIF (code_group==13) THEN
   !
   !  C_3v
   !
   IF (nclass /= 3) CALL errore('divide_class','Wrong classes for C_3v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSEIF (ts==5) THEN
         which_irr(iclass)=3
      ELSE
         CALL errore('divide_class','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==14) THEN
   !
   !  C_4v
   !
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for C_4v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSEIF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2(ax, iax)
         IF (iax < 4) THEN
            !
            !   x, y, or z
            ! 
            which_irr(iclass)=4
         ELSE
            !
            !  all other cases
            !
            which_irr(iclass)=5
         ENDIF
      ENDIF
   ENDDO

ELSEIF (code_group==15) THEN
   !
   !  C_6v
   !
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for C_6v',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         ENDIF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=4
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)), ax)
         CALL which_c2(ax, iax)
         IF (iax==2 .OR. iax==12 .OR. iax==13) THEN
            which_irr(iclass)=6
         ELSEIF (iax==1 .OR. iax==10 .OR. iax==11) THEN
            which_irr(iclass)=5
         ELSE
            CALL errore('divide_class','C_6v mirror not recognized',1)
         ENDIF
      ENDIF
   ENDDO
ELSEIF (code_group==16) THEN
   !
   !  C_2h
   !
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for C_2h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=2
      ELSEIF (ts==2) THEN
         which_irr(iclass)=3
      ELSEIF (ts==5) THEN
         which_irr(iclass)=4
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==17) THEN
   !
   !  C_3h
   !
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for C_3h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         END IF
      ELSEIF (ts==5) THEN
         which_irr(iclass)=4
      ELSEIF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=5
         ELSE
            which_irr(iclass)=6
         END IF
      ELSE
         CALL errore('divide_class','wrong sym_type',1)
      ENDIF
   ENDDO
ELSEIF (code_group==18) THEN
   !
   !  C_4h
   !
   IF (nclass /= 8) CALL errore('divide_class','Wrong classes for C_4h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (angle_rot(smat(1,1,elem(1,iclass)))-90.d0<eps) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=4
         END IF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=3
      ELSEIF (ts==2) THEN
         which_irr(iclass)=5
      ELSEIF (ts==5) THEN
         which_irr(iclass)=7
      ELSEIF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-90.d0)<eps) THEN
            which_irr(iclass)=8
         ELSE
            which_irr(iclass)=6
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==19) THEN
!
!  C_6h
!
   IF (nclass /= 12) CALL errore('divide_class','Wrong classes for C_6h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=2
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=3
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=5
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=6
         END IF
      ELSEIF (ts==4) THEN
         which_irr(iclass)=4
      ELSEIF (ts==2) THEN
         which_irr(iclass)=7
      ELSEIF (ts==5) THEN
         which_irr(iclass)=10
      ELSEIF (ts==6) THEN
         ars=angle_rot_s(smat(1,1,elem(1,iclass)))
         IF (ABS(ars-60.d0)<eps) THEN
            which_irr(iclass)=11
         ELSEIF (ABS(ars-120.d0)<eps) THEN
            which_irr(iclass)=12
         ELSEIF (ABS(ars-240.d0)<eps) THEN
            which_irr(iclass)=8
         ELSEIF (ABS(ars-300.d0)<eps) THEN
            which_irr(iclass)=9
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      ENDIF
   ENDDO
ELSEIF (code_group==20) THEN
   !
   !  D_2h
   !
   !  mirror_axis gives the normal to the mirror plane
   !
   IF (nclass /= 8) CALL errore('divide_class','Wrong classes for D_2h',1)
   iax=0
   ibx=0
   icx=0
   imax=0
   imbx=0
   imcx=0
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (iax==0) THEN
            CALL which_c2(ax, iax)
            aclass=iclass
         ELSEIF (ibx==0) THEN
            CALL which_c2(ax, ibx)
            bclass=iclass
         ELSEIF (icx==0) THEN
            CALL which_c2(ax, icx)
            cclass=iclass
         ELSE
            CALL errore('divide_class','D_2h too many C_2 axis',1)
         ENDIF 
      ELSEIF (ts==2) THEN
         which_irr(iclass)=5
      ELSEIF (ts==5) THEN
         CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
         IF (imax==0) THEN
            CALL which_c2(ax, imax)
            amclass=iclass
         ELSEIF (imbx==0) THEN
            CALL which_c2(ax, imbx)
            bmclass=iclass
         ELSEIF (imcx==0) THEN
            CALL which_c2(ax, imcx)
            cmclass=iclass
         ELSE
            CALL errore('divide_class','D_2h too many mirrors',1)
         ENDIF 
      ELSE
         CALL errore('divide_class','D_2h operation not recognized',1)
      ENDIF
   ENDDO
   !
   CALL is_d2(iax, ibx, icx, ind2)
   !
   which_irr(aclass)=ind2(1) + 1 
   which_irr(bclass)=ind2(2) + 1 
   which_irr(cclass)=ind2(3) + 1  
   !
   IF (imax==iax) which_irr(amclass) = which_irr(aclass) + 4
   IF (imax==ibx) which_irr(amclass) = which_irr(bclass) + 4
   IF (imax==icx) which_irr(amclass) = which_irr(cclass) + 4
   IF (imbx==iax) which_irr(bmclass) = which_irr(aclass) + 4
   IF (imbx==ibx) which_irr(bmclass) = which_irr(bclass) + 4
   IF (imbx==icx) which_irr(bmclass) = which_irr(cclass) + 4
   IF (imcx==iax) which_irr(cmclass) = which_irr(aclass) + 4
   IF (imcx==ibx) which_irr(cmclass) = which_irr(bclass) + 4
   IF (imcx==icx) which_irr(cmclass) = which_irr(cclass) + 4
   !
ELSEIF (code_group==21) THEN
   !
   !  D_3h
   !
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for D_3h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==5) THEN
         IF (nelem(iclass)>1) THEN
            which_irr(iclass)=6
         ELSE 
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==6) THEN
         which_irr(iclass)=5
      END IF
   END DO
ELSEIF (code_group==22) THEN
   !
   !  D_4h
   !
   !
   !  First search the order 4 axis
   !
   IF (nclass /= 10) CALL errore('divide_class','Wrong classes for D_4h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         axis=0
         DO ipol=1,3
            IF (is_axis(ax,ipol)) axis=ipol
         ENDDO 
         IF (axis==0) call errore('divide_class','unknown D_4h axis ',1)
      ENDIF
   END DO
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         which_irr(iclass)=0
         CALL versor(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=3
         ELSE
            DO ipol=1,3
               IF (is_axis(ax,ipol)) which_irr(iclass)=4
            ENDDO
            IF (which_irr(iclass)==0) which_irr(iclass)=5
         END IF
      ELSEIF (ts==2) THEN
         which_irr(iclass)=6
      ELSEIF (ts==5) THEN
         which_irr(iclass)=0
         CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
         IF (is_axis(ax,axis)) THEN
            which_irr(iclass)=8
         ELSE 
            DO ipol=1,3
               IF (is_axis(ax,ipol)) which_irr(iclass)=9
            ENDDO
            IF (which_irr(iclass)==0) which_irr(iclass)=10
         END IF
      ELSEIF (ts==6) THEN
         which_irr(iclass)=7
      END IF
   END DO
ELSEIF (code_group==23) THEN
   !
   !  D_6h
   !
   IF (nclass /= 12) CALL errore('divide_class','Wrong classes for D_6h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         ars=angle_rot(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         END IF
      ELSE IF (ts==4) THEN
         IF (nelem(iclass)==1) THEN
            which_irr(iclass)=4
         ELSE
            CALL versor(smat(1,1,elem(1,iclass)),ax)
            CALL which_c2(ax, iax)
            IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
               which_irr(iclass)=5
            ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
               which_irr(iclass)=6
            ELSE
               CALL errore('divide_class','Problem with C_2 of D_6h',1)
            ENDIF
         END IF
      ELSE IF (ts==2) THEN
          which_irr(iclass)=7
      ELSE IF (ts==5) THEN
          IF (nelem(iclass)==1) THEN
             which_irr(iclass)=10
          ELSE 
             CALL mirror_axis(smat(1,1,elem(1,iclass)),ax)
             CALL which_c2(ax, iax)
             IF (iax==1 .OR. iax==10 .OR. iax==11) THEN
                which_irr(iclass)=11
             ELSEIF (iax==2 .OR. iax==12 .OR. iax==13) THEN
                which_irr(iclass)=12
             ELSE
                CALL errore('divide_class','Problem with mirror of D_6h',1)
             ENDIF
          END IF
      ELSE IF (ts==6) THEN
         ars=angle_rot_s(smat(1,1,elem(1,iclass)))
         IF ((ABS(ars-60.d0)<eps).OR.(ABS(ars-300.d0)<eps)) THEN
             which_irr(iclass)=9
         ELSE
             which_irr(iclass)=8
         END IF
      END IF
   END DO
ELSEIF (code_group==24) THEN
   !
   !  D_2d
   !
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for D_2d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==6) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         IF (nelem(iclass)==1) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==5) THEN
         which_irr(iclass)=5
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==25) THEN
   !
   !  D_3d
   !
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for D_3d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==2) THEN
         which_irr(iclass)=4
      ELSE IF (ts==6) THEN
         which_irr(iclass)=5
      ELSE IF (ts==5) THEN
         which_irr(iclass)=6
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==26) THEN
    !
    !  S_4
    !
    IF (nclass /= 4) CALL errore('divide_class','Wrong classes for S_4',1)
    DO iclass=2,nclass
       ts=tipo_sym(smat(1,1,elem(1,iclass)))
       IF (ts==4) THEN
          which_irr(iclass)=3
       ELSE IF (ts==6) THEN
          IF (ABS(angle_rot_s_tpw(smat(1,1,elem(1,iclass)))-90.d0)<eps) THEN
             which_irr(iclass)=4
          ELSE
             which_irr(iclass)=2
          END IF
       ELSE
          CALL errore('divide_class','wrong operation',1)
       END IF
    END DO
ELSE IF (code_group==27) THEN
   !
   !  S_6
   !
   IF (nclass /= 6) CALL errore('divide_class','Wrong classes for S_6',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=3
         END IF
      ELSE IF (ts==2) THEN
         which_irr(iclass)=4
      ELSE IF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-60.d0)<eps) THEN
            which_irr(iclass)=6
         ELSE
            which_irr(iclass)=5
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==28) THEN
   !
   !  T
   !
   IF (nclass /= 4) CALL errore('divide_class','Wrong classes for T',1)
   DO iclass=2,nclass
      ars=angle_rot(smat(1,1,elem(1,iclass)))
      IF (ABS(ars-120.d0)<eps) THEN
         which_irr(iclass)=3
      ELSE IF (ABS(ars-240.d0)<eps) THEN
         which_irr(iclass)=4
      ELSE IF (ABS(ars-180.d0)<eps) THEN
         which_irr(iclass)=2
      ELSE
         CALL errore('divide_class','wrong angle',1)
      END IF
   END DO
ELSE IF (code_group==29) THEN
   !
   !  T_h
   !
   IF (nclass /= 8) CALL errore('divide_class','Wrong classes for T_h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         IF (ABS(angle_rot(smat(1,1,elem(1,iclass)))-120.d0)<eps) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==4) THEN
         which_irr(iclass)=2
      ELSE IF (ts==2) THEN
         which_irr(iclass)=5
      ELSE IF (ts==6) THEN
         IF (ABS(angle_rot_s(smat(1,1,elem(1,iclass)))-60.d0)<eps) THEN
            which_irr(iclass)=8
         ELSE
            which_irr(iclass)=7
         END IF
      ELSE IF (ts==5) THEN
         which_irr(iclass)=6
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==30) THEN
   !
   !  T_d
   !
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for T_d',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==3) THEN
         which_irr(iclass)=2
      ELSE IF (ts==4) THEN
         which_irr(iclass)=3
      ELSE IF (ts==6) THEN
         which_irr(iclass)=4
      ELSE IF (ts==5) THEN
         which_irr(iclass)=5
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==31) THEN
   !
   !  O
   !
   IF (nclass /= 5) CALL errore('divide_class','Wrong classes for O',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         IF (nelem(iclass)==3) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=5
         END IF
      ELSE IF (ts==3) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE
         CALL errore('divide_class','wrong operation',1)
      END IF
   END DO
ELSEIF (code_group==32) THEN
   !
   !  O_h
   !
   IF (nclass /= 10) CALL errore('divide_class','Wrong classes for O_h',1)
   DO iclass=2,nclass
      ts=tipo_sym(smat(1,1,elem(1,iclass)))
      IF (ts==4) THEN
         IF (nelem(iclass)==3) THEN
            which_irr(iclass)=3
         ELSE
            which_irr(iclass)=5
         END IF
      ELSE IF (ts==3) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=2
         ELSE
            which_irr(iclass)=4
         END IF
      ELSE IF (ts==2) THEN
         which_irr(iclass)=6
      ELSE IF (ts==5) THEN
         IF (nelem(iclass)==6) THEN
            which_irr(iclass)=10
         ELSE
            which_irr(iclass)=8
         END IF
      ELSE IF (ts==6) THEN
         IF (nelem(iclass)==8) THEN
            which_irr(iclass)=7
         ELSE
            which_irr(iclass)=9
         END IF
      ELSE
         CALL errore('divide_class_so','wrong operation',1)
      END IF
   ENDDO
ELSE
 CALL errore('divide_class','code_group not correct',1)
ENDIF
!
RETURN
!
END SUBROUTINE divide_class_tpw

