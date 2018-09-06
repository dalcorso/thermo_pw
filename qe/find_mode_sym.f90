!
! Copyright (C) 2006-2011 Quantum ESPRESSO group
! Copyright (C) 2017 Andrea Dal Corso 
! Removed the necessity to separate the modes into groups using
! a threshold on the frequencies.
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_mode_sym_tpw (u, w2, tau, nat, nsym, s, sr, irt, xq,    &
     rtau, amass, ntyp, ityp, flag, lmolecule, lstop, num_rap_mode, ierr)
  !
  !   This subroutine finds the irreducible representations which give
  !   the transformation properties of eigenvectors of the dynamical
  !   matrix. It does NOT work at zone border in non symmorphic space groups.
  !   if flag=1 the displacements are given in input, otherwise the
  !   eigenvectors of the dynamical matrix.
  !   The output of this routine is only num_rap_mode, the number of
  !   the irreducible representation for each mode.
  !   error conditions:
  !   num_rap_mode(i)= 0   ! the routine could not determine mode symmetry
  !
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry, ry_to_cmm1
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_is, ONLY : gname_is
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::             &
       nat,         &     ! number of atoms
       nsym,        &     ! number of symmetries
       flag,        &     ! if 1 u are displacements, if 0 u are eigenvectors
       ntyp,        &     ! number of atomic types
       ityp(nat),   &     ! the type of each atom
       irt(48,nat)        ! the rotated of each atom
  INTEGER, INTENT(OUT) :: num_rap_mode ( 3 * nat )

  INTEGER, INTENT(OUT) :: ierr ! 0 if the routine determined mode symmetry

  REAL(DP), INTENT(IN) ::   &
       xq(3),          &  ! the q vector of the modes
       tau(3,nat),     &  ! the atomic coordinates
       rtau(3,48,nat), &  ! the R vector for each rotated atom
       amass(ntyp),    &  ! the mass of the atoms
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
       s(3,3,48), &   ! rotation matrices
       nmodes,    &   ! number of modes
       imode,     &   ! counter on modes
       igroup,    &   ! counter on groups
       nu_i, mu,  &   ! counters on modes
       irot,      &   ! select a rotation
       irap,      &   ! counter on representations
       iclass,    &   ! counter on classes
       mult,      &   ! multiplicity of the representation
       na,        &   ! counter on atoms
       i              ! generic counter

  INTEGER, ALLOCATABLE :: istart(:)

  COMPLEX(DP) :: times              ! safe dimension
  ! in case of accidental degeneracy
  COMPLEX(DP), EXTERNAL :: zdotc
  REAL(DP) :: sumt
  COMPLEX(DP), ALLOCATABLE ::  rmode(:,:), trace(:,:), z(:,:), w(:,:)
  LOGICAL :: is_linear
  INTEGER :: counter, counter_s, dim_rap
  LOGICAL :: found
  INTEGER :: invs(48), ss(3,3), isym, jsym
  !
  !    Divide the modes on the basis of the mode degeneracy.
  !
  ierr=0
  num_rap_mode=0
  nmodes=3*nat

  ALLOCATE(istart(nmodes+1))
  ALLOCATE(z(nmodes,nmodes))
  ALLOCATE(rmode(nmodes,nmodes))
  ALLOCATE(w(12,nmodes))
  ALLOCATE(trace(12,nmodes))

  IF (flag==1) THEN
     !
     !  Find the eigenvectors of the dynamical matrix
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

  DO isym = 1, nsym
     found = .false.
     DO jsym = 1, nsym
        !
        ss = matmul (s(:,:,jsym),s(:,:,isym))
        ! s(:,:,1) is the identity
        IF ( all ( s(:,:,1) == ss(:,:) ) ) THEN
           invs (isym) = jsym
           found = .true.
        ENDIF
     ENDDO
     IF ( .NOT.found) CALL errore ('inverse_s', ' Not a group', 1)
  ENDDO
!
!  The symmetry of these modes is not computed
!
  IF (lmolecule) THEN
     istart(1)=7
     IF(is_linear(nat,tau)) istart(1)=6
  ELSE
     istart(1)=1
  ENDIF
  !
  !  Find the character of one symmetry operation per class on all the
  !  modes
  !
  DO iclass=1,nclass
     irot=elem(1,iclass)
!
!   rotate all modes together
!
     CALL rotate_mod(z,rmode,sr(1,1,irot),irt,rtau,xq,nat,invs(irot))
     DO nu_i=istart(1),3*nat
        w(iclass,nu_i)=zdotc(3*nat,z(1,nu_i),1,rmode(1,nu_i),1)
     ENDDO
  END DO
  !
  !  And now computes the trace for each group of degenerate modes.
  !  We continue to add diagonal elements to the trace, until we
  !  find a set of traces whose sum of square moduli is an integer
  !  multiple of the group order. 
  ! 
  trace=(0.d0,0.d0)
  ngroup=1
  DO nu_i=istart(1), nmodes
     DO iclass=1,nclass
        trace(iclass,ngroup)=trace(iclass,ngroup) + w(iclass, nu_i)
     ENDDO
     sumt=0.0_DP
     DO iclass=1,nclass
        sumt=sumt+ABS(trace(iclass,ngroup))**2*nelem(iclass)
     ENDDO
     sumt=sumt/nsym
!
!    If sumt is an integer we have found an irreducible representation or
!    an integer number of irreducible representations.
!    We can start to identify a new group of modes.
!
     IF (ABS(NINT(sumt)-sumt) < 1.d-5) THEN
        ngroup=ngroup+1
        istart(ngroup)=nu_i+1
     ENDIF
  ENDDO
  ngroup=ngroup-1
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of modes
  !
  DO igroup=1,ngroup
     counter=istart(igroup)
     dim_rap=istart(igroup+1)-istart(igroup)
!
!   If the frequency is so small probably it has not been calculated.
!
     IF (SQRT(ABS(w2(counter)))*ry_to_cmm1<1.d-3) CYCLE
     DO irap=1,nclass
        times=(0.d0,0.d0)
        DO iclass=1,nclass
           times=times+trace(iclass,igroup)*CONJG(char_mat(irap, &
                which_irr(iclass)))*nelem(iclass)
           !         write(6,*) igroup, irap, iclass, which_irr(iclass)
        ENDDO
        times=times/nsym
        mult=NINT(ABS(DBLE(times)))
!
!   times must be a positive integer or zero, otherwise some error occured
!   somewhere
!
        IF ((ABS(mult-DBLE(times)) > 1.d-4).OR. &
             (ABS(AIMAG(times)) > eps) ) THEN 
           IF (lstop) THEN
              CALL errore('find_mode_sym','unknown mode symmetry',1)
           ELSE
              counter=counter + dim_rap - 1
              ierr=1
           ENDIF
        ELSE
!
!    If the code arrives here, no error occured and we can set the mode
!    symmetry for all the modes of the group
!
           IF (ABS(times) > eps) THEN
              IF (ABS(mult-DBLE(times)) < 1.d-4) THEN
                 counter_s=counter
                 DO imode=counter_s, counter_s+mult*&
                                              NINT(DBLE(char_mat(irap,1)))-1
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
  DEALLOCATE(w)
  DEALLOCATE(rmode)
  DEALLOCATE(istart)

  RETURN
END SUBROUTINE find_mode_sym_tpw
