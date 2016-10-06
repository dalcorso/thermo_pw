!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_mode_sym_proj (u, w2, tau, nat, nsym, s, sr, ft, gk, invs, &
                    irt, xq, rtau, amass, ntyp, ityp, flag, lmolecule, &
                    lstop, num_rap_mode, ierr)
  !
  !   This subroutine finds the irreducible representations which give
  !   the transformation properties of eigenvectors of the dynamical
  !   matrix. It is used at zone border for non symmorphic space groups.
  !   if flag=1 the true displacements are given in input, otherwise the
  !   eigenvectors of the dynamical matrix are given.
  !   The output of this routine is only num_rap_mode, the number of
  !   the irreducible representation for each mode.
  !   error conditions:
  !   num_rap_mode(i)= 0   ! the routine could not determine mode symmetry
  !
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry, RY_TO_CMM1, tpi
  USE rap_point_group, ONLY : code_group
  USE point_group,  ONLY : find_factor_system
  USE proj_rap_point_group,   ONLY : which_elem, char_mat_proj, nrap_proj, &
                                     code_groupq_ext
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::             &
       nat,         &     ! number of atoms
       nsym,        &     ! number of symmetries
       flag,        &     ! if 1 u are displacements, if 0 u are eigenvectors
       ntyp,        &     ! number of atomic types
       ityp(nat),   &     ! the type of each atom
       gk(3,48),    &     ! the g vector associated to each S
       s(3,3,48),   &     ! the symmetry matrices
       invs(48),    &     ! the inverse of each matrix
       irt(48,nat)        ! the rotated of each atom

  INTEGER, INTENT(OUT) :: num_rap_mode ( 3 * nat )

  INTEGER, INTENT(OUT) :: ierr ! 0 if the routine determined mode symmetry

  REAL(DP), INTENT(IN) ::   &
       xq(3),          &  ! the q vector of the modes
       tau(3,nat),     &  ! the atomic coordinates
       rtau(3,48,nat), &  ! the R vector for each rotated atom
       amass(ntyp),    &  ! the mass of the atoms
       ft(3,48),       &  ! fractional translations in crystal coordinates
       w2(3*nat),      &  ! the square of the frequencies
       sr(3,3,48)         ! the rotation matrices in real space.

  COMPLEX(DP), INTENT(IN) ::  &
       u(3*nat, 3*nat)       ! The eigenvectors or the displacement pattern

  LOGICAL, INTENT(IN) :: lmolecule, & ! if .true. these are eigenvalues of an
                                   ! isolated system and do not find the
                                   ! symmetry of the first six eigenvectors,
                                   ! or five for a linear molecule.
                         lstop     ! if .true. the routine stops if it
                                   ! does not understand the symmetry of a 
                                   ! mode

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::      &
       ngroup,    &   ! number of different frequencies groups
       nmodes,    &   ! number of modes
       imode,     &   ! counter on modes
       igroup,    &   ! counter on groups
       igroup_t,  &   ! index of the test groups
       dim_rap_t, &   ! dimension of the test representation
       nu_i, nu_j, mu,  &   ! counters on modes
       irot,      &   ! select a rotation
       irap,      &   ! counter on representations
       na,        &   ! counter on atoms
       i,j            ! generic counters

  INTEGER, ALLOCATABLE :: istart(:), dim_rap(:)

  COMPLEX(DP) :: times, phase              

  COMPLEX(DP), EXTERNAL :: zdotc
  REAL(DP), ALLOCATABLE :: w1(:)
  REAL(DP) :: arg
  COMPLEX(DP), ALLOCATABLE ::  rmode(:,:), trace(:,:), z(:,:), sym_mat(:,:,:)
  COMPLEX(DP) :: factor(48,48)
  LOGICAL :: is_linear
  INTEGER :: counter, counter_s, isym, jsym
!
!   A few initializations
!
  ierr=0
  num_rap_mode=0
  nmodes=3*nat

  ALLOCATE(istart(nmodes+1))
  ALLOCATE(dim_rap(nmodes))
  ALLOCATE(z(nmodes,nmodes))
  ALLOCATE(w1(nmodes))
  ALLOCATE(rmode(nmodes,nmodes))
  ALLOCATE(trace(48,nmodes))

  IF (flag==1) THEN
     !
     !  Find the eigenvalues of the dynmaical matrix
     !  Note that amass is in amu; amu_ry converts it to Ry au
     !
     DO nu_i = 1, nmodes
        DO mu = 1, nmodes
           na = (mu - 1) / 3 + 1
           z (mu, nu_i) = u (mu, nu_i) * SQRT (amu_ry*amass (ityp (na) ) )
        END DO
     END DO
  ELSE
     z=u
  ENDIF
!
!  Compute the mode frequency in cm-1. Two modes are considered degenerate
!  if their frequency is lower 0.05 cm-1
! 
  w1(:)=SIGN(SQRT(ABS(w2(:)))*RY_TO_CMM1,w2(:))

  ngroup=1
  istart(ngroup)=1
!
!  The symmetry of these modes is not computed
!
  IF (lmolecule) THEN
     istart(1)=7
     IF(is_linear(nat,tau)) istart(1)=6
  ENDIF
!
! The other modes are divided into groups of degenerate modes
!
  DO imode=istart(1)+1,nmodes
     IF (ABS(w1(imode)-w1(imode-1)) > 5.0d-2) THEN
        ngroup=ngroup+1
        istart(ngroup)=imode
     END IF
  END DO
  istart(ngroup+1)=nmodes+1
  !
  !  Find the character of all the symmetry operations 
  !
!  igroup_t=1
!  dim_rap_t=istart(igroup_t+1) - istart(igroup_t)
!  ALLOCATE(sym_mat(dim_rap_t,dim_rap_t,48))

  DO igroup=1,ngroup
     dim_rap(igroup)=istart(igroup+1)-istart(igroup)
     DO irot=1,nsym
!
!    rotate all modes together
!
        CALL rotate_mod(z,rmode,sr(1,1,irot),irt,rtau,xq,nat,invs(irot))

        arg=tpi * ( gk(1,invs(irot))*ft(1,irot) +  &
                    gk(2,invs(irot))*ft(2,irot) +  &
                    gk(3,invs(irot))*ft(3,irot) )
        phase=cmplx(cos(arg),sin(arg),kind=DP)
        rmode=rmode*phase

        trace(irot,igroup)=(0.d0,0.d0)
        DO i=1,dim_rap(igroup)
           nu_i=istart(igroup)+i-1
           trace(irot,igroup)=trace(irot,igroup) + &
                zdotc(3*nat,z(1,nu_i),1,rmode(1,nu_i),1)
        END DO
!
!       write(6,*) 'group,class',igroup, irot, trace(irot,igroup)
!
!   For testing purposes it is possible to compute the complete matrix 
!   of the representation and to study the factor system
!
!
!        DO i=1,dim_rap_t
!           nu_i=istart(igroup)+i-1
!           DO j=1,dim_rap_t
!              nu_j=istart(igroup)+j-1
!              sym_mat(i,j,which_elem(irot)) = zdotc(3*nat,z(1,nu_i),1, &
!                                                          rmode(1,nu_j),1)
!           END DO
!        END DO

     END DO
  END DO

!  WRITE(stdout,'(/,5x,"The factor system of the phonon representation",/)')
!  CALL find_factor_system(sym_mat, dim_rap_t, nsym, code_groupq_ext, &
!                          factor, .TRUE.)
!  DEALLOCATE(sym_mat)
!
!  DO igroup=1,ngroup
!     DO isym=1,nsym
!        DO jsym=1,nsym
!           IF (which_elem(jsym)==isym) &
!              WRITE(stdout,'(2i5,4f11.8)') igroup,jsym,trace(jsym,igroup) 
!        ENDDO
!     ENDDO
!     WRITE(stdout,*)
!  ENDDO
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of modes
  !
  DO igroup=1,ngroup
     counter=istart(igroup)
!
!   If the frequency is so small probably it has not been calculated.
!   This value 
!
     IF (ABS(w1(counter))<1.d-3) CYCLE
     DO irap=1,nrap_proj
        times=(0.d0,0.d0)
        DO irot=1,nsym
           times=times+trace(irot,igroup)*CONJG(char_mat_proj(irap, &
                                               which_elem(irot)))
        ENDDO
        times=times/nsym
!        WRITE(6,*) 'igroup, irap, times', igroup, irap, times
!
!   times must be a positive integer or zero, otherwise some error occured
!   somewhere
!
        IF ((ABS(NINT(ABS(DBLE(times)))-DBLE(times)) > 1.d-4).OR. &
             (ABS(AIMAG(times)) > eps) ) THEN 
           IF (lstop) THEN
              CALL errore('find_mode_sym','unknown mode symmetry',1)
           ELSE
              counter=counter + dim_rap(igroup)-1
              ierr=1
           ENDIF
        ELSE
!
!    If the code arrives here, no error occured and we can set the mode
!    symmetry for all the modes of the group
!
           IF (ABS(times) > eps) THEN
              IF (ABS(NINT(DBLE(times))-DBLE(times)) < 1.d-4) THEN
                 counter_s=counter
                 DO imode=counter_s, counter_s+NINT(DBLE(times))*&
                                       NINT(DBLE(char_mat_proj(irap,1)))-1
                    num_rap_mode(imode) = irap
                    counter=counter+1
                 ENDDO
              END IF
           END IF
        END IF
     END DO
  END DO

100 CONTINUE

  DEALLOCATE(trace)
  DEALLOCATE(z)
  DEALLOCATE(w1)
  DEALLOCATE(rmode)
  DEALLOCATE(dim_rap)
  DEALLOCATE(istart)

  RETURN
END SUBROUTINE find_mode_sym_proj

SUBROUTINE print_mode_sym_proj(w2, num_rap_mode, ptype)
!
!  This routine prints the vibrational frequencies and the 
!  symmetry of the eigenvectors of the dynamical matrix. It is used 
!  at zone border for space groups that contains nonsymmorphic operations
!
USE kinds, ONLY : DP
USE constants, ONLY : ry_to_cmm1
USE noncollin_module, ONLY : nspin_mag
USE ions_base, ONLY : nat
USE io_global, ONLY : stdout
USE proj_rap_point_group, ONLY : char_mat_proj, name_rap_proj, nrap_proj, &
                                 code_groupq_ext
USE point_group, ONLY : print_ptype_info
USE rap_point_group, ONLY : gname
USE rap_point_group_is, ONLY : gname_is

IMPLICIT NONE
REAL(DP), INTENT(IN) :: w2( 3*nat )
INTEGER, INTENT(IN) :: num_rap_mode( 3*nat ), ptype(3)

REAL(DP) :: w1( 3*nat )
INTEGER :: next, irap, imode
!
!  Transform the frequencies to cm^-1
!
w1(:)=SIGN(SQRT(ABS(w2(:)))*ry_to_cmm1,w2(:))
!
!  prints the name of the point group 
!
CALL print_ptype_info(ptype, code_groupq_ext)
IF ( nspin_mag == 4 ) THEN
   WRITE(stdout,  &
          '(5x,"Mode symmetry, ",a11," [",a11,"] magnetic group:",/)') gname, &
                                                                       gname_is
ELSE
   WRITE(stdout,'(/,5x,"Mode symmetry:",/)') 
END IF
!
! for each mode, or group of degenerate modes, writes the name of the
! irreducible representation
!
next=0
DO imode = 1, 3 * nat
   IF ( imode < next .OR. ABS(w1(imode)) < 1.d-3 ) CYCLE
   IF (num_rap_mode(imode) == 0)  THEN
      WRITE(stdout,'(5x,"freq (",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x, "-->   ?")') imode, imode, w1(imode)
   ELSE
      irap=num_rap_mode(imode)
      next = imode + NINT(DBLE(char_mat_proj(irap,1)))
      WRITE(stdout,'(5x,"freq (",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,&
                    &"--> ",a)') &
           imode, next-1, w1(imode), TRIM(name_rap_proj(irap))
   ENDIF
ENDDO

RETURN
END SUBROUTINE print_mode_sym_proj

SUBROUTINE prepare_sym_analysis_proj(nsymq,s,sr,ft,gii,ptype,gcode_old)
!
!  set gcode_old to -1 to avoid any output writing from this routine
!
USE kinds,         ONLY : DP
USE point_group,   ONLY : find_group_info_ext, find_projection_type, &
                          find_irr_proj
USE proj_rap_point_group, ONLY : which_elem, group_desc, char_mat_proj, & 
                          name_rap_proj, nrap_proj, code_groupq_ext, & 
                          qptype, qgauge 
USE rap_point_group, ONLY : code_group, gname

IMPLICIT NONE

INTEGER :: nsymq, gcode_old
INTEGER :: s(3,3,48), gii(3,48), ptype(3)
REAL(DP) :: sr(3,3,48), ft(3,48), argument(48,48), gauge(48)
INTEGER :: code_group_ext
LOGICAL :: lwrite

CALL find_group(nsymq,sr,gname,code_group)
lwrite=.FALSE.
IF (code_group/=gcode_old .AND. gcode_old /= -1 ) lwrite=.TRUE.

CALL find_group_info_ext(nsymq,sr,code_group,code_group_ext, &
                                                    which_elem, group_desc)
CALL set_factor_system(argument, s, ft, gii, nsymq, &
                            which_elem, .FALSE., .FALSE., code_group_ext)

CALL find_projection_type(code_group, code_group_ext, argument, &
                                      ptype, gauge, .FALSE.)
CALL find_irr_proj(code_group_ext,char_mat_proj,name_rap_proj, &
                      nrap_proj, nsymq, ptype, gauge, lwrite)
qptype(:)=ptype(:)
qgauge(:)=gauge(:)
code_groupq_ext=code_group_ext

RETURN
END SUBROUTINE prepare_sym_analysis_proj
