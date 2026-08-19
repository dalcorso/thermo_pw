!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE setup_nscf_tpw ( newgrid, xq, elph_mat )
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes variables for the non-scf calculations at k
  ! ... and k+q required by the linear response calculation at finite q.
  ! ... In particular: finds the symmetry group of the crystal that leaves
  ! ... the phonon q-vector (xq) or the single atomic displacement (modenum)
  ! ... unchanged; determines the k- and k+q points in the irreducible BZ
  ! ... Needed on input (read from data file):
  ! ... "nsym" crystal symmetries s, t_rev, "nrot" lattice symetries "s"
  ! ... "nkstot" k-points in the irreducible BZ wrt lattice symmetry
  ! ... Produced on output:
  ! ... symmetries ordered with the "nsymq" phonon symmetries first
  ! ... "nkstot" k- and k+q-points in the IBZ calculated for the phonon sym.)
  ! ... Misc. data needed for running the non-scf calculation
  !
  USE kinds,              ONLY : DP
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, tpiba
  USE ions_base,          ONLY : nat, tau, ityp, zv
  USE force_mod,          ONLY : force
  USE basis,              ONLY : natomwfc
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 ltetra, nkstot, qnorm
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
  USE symm_base,          ONLY : s, t_rev, nrot, nsym, time_reversal
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : ethr, isolve, david, max_cg_iter, &
                                 noinv, use_para_diag
  USE mp_pools,           ONLY : kunit
  USE control_lr,         ONLY : lgamma
  USE noncollin_module,   ONLY : noncolin, domag
  USE start_k,            ONLY : nks_start, xk_start, wk_start, &
                                 nk1, nk2, nk3, k1, k2, k3
  USE optical,            ONLY : lmagnon
  USE paw_variables,      ONLY : okpaw
  USE lr_symm_base,       ONLY : nsymq, invsymq, minus_q
  USE upf_ions,           ONLY : n_atom_wfc
  USE ktetra,             ONLY : tetra, tetra_type, opt_tetra_init
  USE band_computation,   ONLY : diago_bands, isym_bands, ik_origin, &
                                 sym_for_diago, nks0
 
  !
  IMPLICIT NONE
  !
  REAL (DP), INTENT(IN) :: xq(3)
  LOGICAL, INTENT (IN) :: newgrid
  LOGICAL, INTENT (IN) :: elph_mat
  !
  REAL (DP), ALLOCATABLE :: rtau (:,:,:)
  INTEGER  :: t_rev_eff(48), ik
  LOGICAL  :: magnetic_sym, sym(48)
  LOGICAL  :: skip_equivalence
  LOGICAL, EXTERNAL  :: check_para_diag
  !
  IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
  !
  ! ... threshold for diagonalization ethr - should be good for all cases
  !
  ethr= 1.0D-9 / nelec
  !
  ! ... variables for iterative diagonalization
  ! ... Davidson: isolve=0, david=4 ; CG: isolve=1, david=1
  IF (isolve == 0) THEN
     david = 4
  ELSE IF (isolve == 1) THEN
     david = 1
  ELSE
     call errore('setup_nscf_tpw','erroneous value for diagonalization method. Should be isolve=0 (david) or 1 (cg)',1)
  END IF

  nbndx = david*nbnd
  max_cg_iter=20
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  !
  CALL set_para_diag( nbnd, use_para_diag )
  !
  ! ... Symmetry and k-point section
  !
  ! ... time_reversal = use q=>-q symmetry for k-point generation
  !
  magnetic_sym = noncolin .AND. domag
  !
  ! ... smallg_q flags in symmetry operations of the crystal
  ! ... that are not symmetry operations of the small group of q
  !
  CALL set_small_group_of_q(nsymq,invsymq,minus_q)
  !
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  !
  if( nks_start > 0 .AND. .NOT. newgrid ) then
     !
     !  In this case I keep the same points of the Charge density
     !  calculations
     !
     nkstot = nks_start
     xk(:,1:nkstot) = xk_start(:,1:nkstot)
     wk(1:nkstot)   = wk_start(1:nkstot)
  else
     !
     ! In this case I generate a new set of k-points
     !
     ! In the case of electron-phonon matrix element with wannier functions 
     ! (and possibly in other cases as well) the k-points should not be reduced
     !
     skip_equivalence = elph_mat
     t_rev_eff=0
     CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev_eff, &
                      bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
  endif

  !
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "irreducible_BZ" computes the missing k-points.
  !
  IF(.NOT.elph_mat) THEN
    IF (sym_for_diago) THEN
!
!      Reduce the k points with the point group of the solid, all these points
!      must be explicitely diagonalized
!
       CALL irreducible_BZ(nrot, s, nsym, minus_q, magnetic_sym, &
                       at, bg, npk, nkstot, xk, wk, t_rev)
       diago_bands(1:nkstot)=.TRUE.
       isym_bands(1:nkstot)=1
       DO ik=1,nkstot
          ik_origin(ik)=ik
       ENDDO
       nks0=nkstot
!
!     Now reduce the k points with the small group of q. The k points added
!     here can be calculated using the symmetries of the point group.
!
       CALL irreducible_BZ_tpw (nsym, s, nsymq, minus_q, magnetic_sym, &
                       at, bg, npk, nkstot, xk, wk, t_rev)
       CALL distribute_diago()
    ELSE
       CALL irreducible_BZ(nrot, s, nsymq, minus_q, magnetic_sym, &
                       at, bg, npk, nkstot, xk, wk, t_rev)
       diago_bands(1:nkstot)=.TRUE.
    ENDIF
  ENDIF
  !
  ! ... add k+q to the list of k
  !
  IF (sym_for_diago) THEN
     IF (noncolin.AND.domag) THEN
        CALL set_kplusq_nc_tpw( xk, wk, xq, nkstot, npk, diago_bands, &
                                                isym_bands, ik_origin)
     ELSE
        CALL set_kplusq_tpw( xk, wk, xq, nkstot, npk, diago_bands, isym_bands, &
                                                ik_origin )
     ENDIF
!
!   diagonalize all the bands and forget what has been found by the
!   the last two variables are not used, but we initialize them.
!
!     DO ik=1,nkstot
!        WRITE(6,*) 'ik, idiago_bands, ik_origin', ik, diago_bands(ik), &
!                                                      ik_origin(ik)
!     ENDDO
  ELSE
     IF (noncolin.AND.domag) THEN
        CALL set_kplusq_nc( xk, wk, xq, nkstot, npk)
     ELSE
        CALL set_kplusq( xk, wk, xq, nkstot, npk)
     ENDIF
     diago_bands(1:nkstot)=.TRUE.
     nks0=nkstot
  ENDIF
  !
  ! ... set the granularity for k-point distribution
  !
  IF ( lgamma  ) THEN
     !
     kunit = 1
     IF (noncolin.AND.domag) kunit = 2
     !
  ELSE
     !
     kunit = 2
     IF (noncolin.AND.domag) kunit = 4
     !
  ENDIF
  !
  ! ... Map each k point in the irr.-BZ into tetrahedra
  !
  IF ( ltetra .AND. (tetra_type /= 0) ) THEN
     IF (ALLOCATED(tetra)) DEALLOCATE(tetra)
     CALL opt_tetra_init(nsymq, s, time_reversal .AND. minus_q, t_rev, at, bg,&
          npk, k1, k2, k3, nk1, nk2, nk3, nkstot, xk, kunit)
  END IF
  !
!
!   use in any case set_kplusq_tpw. Needed for magnons.
!
  IF ( lsda ) THEN
     !
     ! ... LSDA case: two different spin polarizations,
     ! ...            each with its own kpoints
     !
     if (nspin /= 2) call errore ('setup_nscf','nspin should be 2; check iosys',1)
     !
     IF (lmagnon) THEN
        CALL set_kup_and_kdw_lm_tpw( xk, wk, isk, nkstot, npk, lgamma, &
                                         diago_bands, isym_bands, ik_origin )
     ELSEIF (sym_for_diago) THEN
        CALL set_kup_and_kdw_tpw( xk, wk, isk, nkstot, npk, lgamma, &
                                         diago_bands, isym_bands, ik_origin )
     ELSE
        CALL set_kup_and_kdw( xk, wk, isk, nkstot, npk)
        diago_bands(1:nkstot)=.TRUE.
     ENDIF
     !
  ELSE IF ( noncolin ) THEN
     !
     ! ... noncolinear magnetism: potential and charge have dimension 4 (1+3)
     !
     if (nspin /= 4) call errore ('setup_nscf','nspin should be 4; check iosys',1)
     current_spin = 1
     !
  ELSE
     !
     ! ... LDA case: the two spin polarizations are identical
     !
     wk(1:nkstot)    = wk(1:nkstot) * degspin
     current_spin = 1
     !
     IF ( nspin /= 1 ) &
        CALL errore( 'setup_nscf', 'nspin should be 1; check iosys', 1 )
     !
  END IF
  !
  IF ( nkstot > npk ) CALL errore( 'setup_nscf', 'too many k points', nkstot )
  !
  ! ...notice: qnorm is used by allocate_nlpot to determine
  ! the correct size of the interpolation table "qrad"
  !
  qnorm = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2) * tpiba
  !
#if defined(__MPI)
  !
  ! ... distribute k-points (and their weights and spin indices)
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  IF (nks==0) CALL errore('setup_nscf_tpw','some pools have no k point',1)
  !
#else
  !
  nks = nkstot
  !
#endif
  !
  RETURN
  !
END SUBROUTINE setup_nscf_tpw

!-------------------------------------------------------------------------
SUBROUTINE distribute_diago()
!-------------------------------------------------------------------------
!
!   This routine distributes the k points that must be diagonalized as
!   uniformely as possible among the k points, so that if the k points
!   are divided in pools each pool receives a similar number of k points
!   to diagonalize
!
USE kinds, ONLY : DP
USE klist, ONLY : nkstot, xk, wk
USE band_computation,  ONLY : diago_bands, isym_bands, ik_origin, nks0

IMPLICIT NONE

INTEGER :: ik, ik1, nstep, nstep1, resto, iktarget
REAL(DP) :: buffer(3), wb
LOGICAL :: lb
INTEGER :: ib

nstep = nkstot / nks0
nstep1= nstep+1
resto = MOD(nkstot,nks0)

DO ik=nks0, 1, -1
   IF (ik<=resto) THEN
      iktarget= (ik-1)*nstep1 + 1
   ELSE
      iktarget= resto*nstep1 + (ik-resto-1)*nstep + 1
   ENDIF 
   IF (ik /= iktarget) THEN
      buffer(:)=xk(:,iktarget)
      xk(:,iktarget)=xk(:,ik)
      xk(:,ik)=buffer(:)

      wb=wk(iktarget)
      wk(iktarget)=wk(ik)
      wk(ik)=wb

      lb=diago_bands(iktarget)
      diago_bands(iktarget)=diago_bands(ik)
      diago_bands(ik)=lb

      ib=ik_origin(iktarget)
      ik_origin(iktarget)=ik_origin(ik)
      ik_origin(ik)=ib

      ib=isym_bands(iktarget)     
      isym_bands(iktarget)=isym_bands(ik)
      isym_bands(ik)=ib

      DO ik1=1, nkstot
         IF (ik_origin(ik1)==ik) THEN
            ik_origin(ik1)=iktarget
         ELSEIF (ik_origin(ik1)==iktarget) THEN
            ik_origin(ik1)=ik
         ENDIF
      ENDDO
   ENDIF
ENDDO

RETURN
END SUBROUTINE distribute_diago
