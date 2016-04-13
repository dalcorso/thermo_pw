!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE find_symmetry(fft_fact)
  !----------------------------------------------------------------------------
  !
  ! This routine is a stripped version of setup that calls only the
  ! routines necessary to set the fft dimension and set the symmetry
  ! matrices. It can have, as input, the factors that must be contained
  ! in the fft mesh to support the fractionary translations.
  ! ...
  ! ...    finds actual crystal symmetry:
  ! ...    s         symmetry matrices in the direct lattice vectors basis
  ! ...    nsym      number of crystal symmetry operations
  ! ...    nrot      number of lattice symmetry operations
  ! ...    ft        fractionary translations
  ! ...    irt       for each atom gives the corresponding symmetric
  ! ...    invsym    if true the system has inversion symmetry
  !
  !
  USE kinds,              ONLY : DP
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp
  USE gvect,              ONLY : gcutm
  USE fft_base,           ONLY : dfftp, dffts
  USE grid_subroutines,   ONLY : realspace_grid_init
  USE gvecs,              ONLY : doublegrid, gcutms, dual
  USE gvecw,              ONLY : ecutwfc
  USE symm_base,          ONLY : s, t_rev, irt, nrot, nsym, invsym, nosym, &
                                 set_sym_bl, find_sym
  USE noncollin_module, ONLY : m_loc, noncolin, i_cons, npol, angle1, angle2
  USE lsda_mod,         ONLY : starting_magnetization, nspin
  USE spin_orb,         ONLY : domag, lspinorb
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: fft_fact(3)
  INTEGER :: na
  LOGICAL :: magnetic_sym
  !

  magnetic_sym=noncolin.AND.domag
  ALLOCATE( m_loc( 3, nat ) )
  m_loc=0.0_DP
  ! time reversal operation is set up to 0 by default
  t_rev = 0
  IF ( noncolin ) THEN
     !
     ! ... wavefunctions are spinors with 2 components
     !
     npol = 2
     !
     ! ... Set the domag variable to make a spin-orbit calculation with zero
     ! ... magnetization
     !
     IF ( lspinorb ) THEN
        !
        domag = ANY ( ABS( starting_magnetization(1:ntyp) ) > 1.D-6 )
        !
     ELSE
        !
        domag = .TRUE.
        !
     END IF
     !
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     END DO
     !
  ELSE
     !
     ! ... wavefunctions are scalars
     !
     IF (lspinorb)  CALL errore( 'find_symmetry ',  &
         'spin orbit requires a non collinear calculation', 1 )
     npol = 1
     !
     !
     IF ( i_cons == 1) then
        do na=1,nat
           m_loc(1,na) = starting_magnetization(ityp(na))
        end do
     end if
     IF ( i_cons /= 0 .AND. nspin ==1) &
        CALL errore( 'find_symmetry', 'this i_cons requires a magnetic calculation ', 1 )
     IF ( i_cons /= 0 .AND. i_cons /= 1 ) &
        CALL errore( 'find_symmetry', 'this i_cons requires a non colinear run', 1 )
  END IF
  !
  ! ... Set the units in real and reciprocal space
  !
  tpiba  = 2.D0 * pi / alat
  tpiba2 = tpiba**2
  !
  ! ... Compute the cut-off of the G vectors
  !
  doublegrid = ( dual > 4.D0 )
  gcutm = dual * ecutwfc / tpiba2
  !
  IF ( doublegrid ) THEN
     !
     gcutms = 4.D0 * ecutwfc / tpiba2
     !
  ELSE
     !
     gcutms = gcutm
     !
  END IF
  !
  ! ... calculate dimensions of the FFT grid
  !
  CALL realspace_grid_init ( dfftp, at, bg, gcutm, fft_fact )
  IF ( gcutms == gcutm ) THEN
     IF ( dffts%nr1 ==0 .AND. dffts%nr2==0 .AND. dffts%nr3==0) THEN
          dffts%nr1 = dfftp%nr1     
          dffts%nr2 = dfftp%nr2     
          dffts%nr3 = dfftp%nr3
          dffts%nr1x= dfftp%nr1x
          dffts%nr2x= dfftp%nr2x     
          dffts%nr3x= dfftp%nr3x
     ENDIF
  END IF
  CALL realspace_grid_init ( dffts, at, bg, gcutms, fft_fact )
  !
  !  ... generate transformation matrices for the crystal point group
  !  ... First we generate all the symmetry matrices of the Bravais lattice
  !
  call set_sym_bl ( )
  !
  !
  ! ... eliminate rotations that are not symmetry operations
  !
  CALL find_sym ( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                 magnetic_sym, m_loc )

  DEALLOCATE(m_loc)
  IF (ALLOCATED(irt)) DEALLOCATE(irt)
  !
  RETURN
  !
END SUBROUTINE find_symmetry
