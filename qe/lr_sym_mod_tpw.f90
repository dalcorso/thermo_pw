! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains routines that rotate the charge density and the
! magnetization induced by a phonon or an electric field perturbation.
!
MODULE lr_sym_mod_tpw

IMPLICIT NONE
SAVE
PRIVATE

PUBLIC psyme_tpw, psyme_fpol_tpw, psymdvscf_tpw

CONTAINS
!
!---------------------------------------------------------------------
SUBROUTINE psyme_fpol_tpw (dvsym)
!---------------------------------------------------------------------
!
!     This routine symmetrizes the change of the charge/magnetization or 
!     potential/magnetic field due to an electric field perturbation. It 
!     is assumed that the perturbations are on the basis of the crystal.
!     The routine works also for complex input fields as needed in frequency
!     dependent calculations.
!
USE kinds,     ONLY : DP
USE cell_base, ONLY : at, bg
USE fft_base,  ONLY : dfftp
USE symm_base, ONLY : nsym, s, sname, invs, t_rev
USE noncollin_module, ONLY : noncolin, nspin_lsda, nspin_mag, domag
USE scatter_mod, ONLY : cgather_sym
USE lr_sym_mod, ONLY : rotate_mesh, ccryst_to_cart_t
IMPLICIT NONE

COMPLEX(DP), INTENT(INOUT) :: dvsym(dfftp%nnr, nspin_mag, 3)
  ! the potential to symmetrize
COMPLEX(DP), ALLOCATABLE ::  aux(:,:), aux_nc(:,:,:)
  ! auxiliary quantities
INTEGER :: is, isym, ipol, jpol, kpol, my_nrxx, nr1x, ir
  ! counter on spin polarization
  ! counter on rotations
  ! 3 counters on polarizations
  ! number of points on this processor
  ! nrx1 size of the FFT mesh array along direction 1
INTEGER, ALLOCATABLE :: rir(:,:)
  ! The rotated of each point for each symmetry
COMPLEX(DP) :: dmags(3), mag(3)
  ! auxiliary to save the magnetization
INTEGER :: stilde(3,3,nsym)
  ! The matrices that rotate the magnetization
  !
IF (nsym == 1) RETURN

nr1x=dfftp%nr1x
my_nrxx=nr1x*dfftp%my_nr2p*dfftp%my_nr3p

ALLOCATE(rir(my_nrxx, nsym))
ALLOCATE (aux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 3))

CALL rotate_mesh(my_nrxx, nsym, rir)

DO is = 1, nspin_lsda
   !
   !  collect all the quantity to symmetrize in all processors
   !
   DO ipol = 1, 3
#if defined(__MPI)
      CALL cgather_sym( dfftp, dvsym(:, is, ipol), aux(:, ipol))
#else
      aux(:,ipol)= dvsym(:, is, ipol)
#endif
      dvsym(:,is,ipol) = (0.0_DP, 0.0_DP)
   ENDDO
   !
   !  symmmetrize. Each processor symmetrizes only the points that it has
   !
   DO isym = 1, nsym
      DO ir = 1, my_nrxx
         IF (rir(ir,isym)==0) CYCLE
         DO ipol = 1, 3
            DO jpol = 1, 3
               dvsym(ir,is,ipol) = dvsym(ir,is,ipol) + &
                          s(ipol,jpol,isym) * aux(rir(ir,isym),jpol)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

DEALLOCATE (aux)
!
!  Rotate also the magnetization in the noncollinear magnetic case
!
IF (noncolin.AND.domag) THEN
   ALLOCATE (aux_nc(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 3, 3))
!
!  set the symmetry matrices that rotate the magnetization
!
   CALL set_stilde(s, stilde, sname, t_rev, invs, nsym)  
!
!  Bring the magnetization in crystal coordinates and
!  collect all the quantity to symmetrize in all processors
!
   DO ipol = 1, 3
      CALL ccryst_to_cart_t (dfftp%nnr, dvsym(:,2:4,ipol), bg, -1)
      DO is = 2, nspin_mag
#if defined(__MPI)
         CALL cgather_sym( dfftp, dvsym(:, is, ipol), aux_nc(:, is-1, ipol))
#else
         aux_nc(:, is-1, ipol)=dvsym(:, is, ipol)
#endif
         dvsym(:,is,ipol) = (0.0_DP, 0.0_DP)
      ENDDO
   ENDDO
   !
   DO isym=1, nsym
      DO ir = 1, my_nrxx
         IF (rir(ir,isym)==0) CYCLE
         DO ipol = 1, 3
            dmags=(0.0_DP,0.0_DP)
            DO is = 1, 3
               DO jpol = 1, 3
                  dmags(is)=dmags(is) + &
                           s(ipol,jpol,isym) * aux_nc(rir(ir,isym),is,jpol)
               ENDDO
            ENDDO
!
! rotate the magnetic moment
!
            DO kpol = 1, 3
               mag(kpol) = stilde(1,kpol,isym)*dmags(1) + &
                           stilde(2,kpol,isym)*dmags(2) + &
                           stilde(3,kpol,isym)*dmags(3)
            ENDDO
!
!   and add to the complete field
!
            IF (t_rev(isym)==1) THEN
               DO kpol = 1, 3
                  dvsym(ir,kpol+1,ipol) = dvsym(ir,kpol+1,ipol) + &
                                                           CONJG(mag(kpol))
               ENDDO
            ELSE
               DO kpol = 1, 3
                  dvsym(ir,kpol+1,ipol) = dvsym(ir,kpol+1,ipol) + mag(kpol)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE (aux_nc)
!
! go back to cartesian coordinates
!
   DO ipol=1,3
      CALL ccryst_to_cart_t (dfftp%nnr, dvsym(:,2:4,ipol), at, 1)
   ENDDO
ENDIF

dvsym = dvsym / DBLE(nsym)

DEALLOCATE(rir)

RETURN
END SUBROUTINE psyme_fpol_tpw

!---------------------------------------------------------------------
SUBROUTINE psymdvscf_tpw (npe, irr, dvsym)
!---------------------------------------------------------------------
!
!     This routine symmetrizes the change of the charge/magnetization 
!     or potential/magnetic field due to npe phonon perturbations that 
!     transform among themselves according to the representation 
!     t(:,:,irot,irr). It works also if irot needs t_rev but t must be
!     correctly defined.
!     The three components of the magnetization or magnetic field are 
!     in cartesian coodinates.
!
USE kinds,     ONLY : DP
USE cell_base, ONLY : at, bg
USE fft_base,  ONLY : dfftp
USE symm_base, ONLY : nsym, s, ft, invs, sname, t_rev
USE noncollin_module, ONLY : noncolin, nspin_lsda, nspin_mag, domag
USE modes,     ONLY : t, tmq
USE lr_symm_base, ONLY : nsymq, gi, minus_q, irotmq, gimq
USE scatter_mod,  ONLY : cgather_sym
USE lr_sym_mod,   ONLY : find_mesh_ijk, rotate_mesh, rotate_mesh_1s, &
                         compute_phase, ccryst_to_cart_t

IMPLICIT NONE
INTEGER, INTENT(IN) :: npe, irr

COMPLEX(DP), INTENT(INOUT) :: dvsym(dfftp%nnr, nspin_mag, npe)
  ! the potential to symmetrize
COMPLEX(DP), ALLOCATABLE ::  aux(:,:), aux1(:,:)
  ! auxiliary quantity
COMPLEX(DP), ALLOCATABLE ::  aux_nc(:,:,:)
  ! auxiliary quantity
COMPLEX(DP), ALLOCATABLE :: phase1(:,:), phase2(:,:), phase3(:,:)
  ! auxiliary quantity for the phases
LOGICAL, ALLOCATABLE :: zb(:)
  ! if true the symmetry needs a G vector
INTEGER :: is, isym, ipert, jpert, kpol, my_nrxx, nr1x, i, j, k, ir, &
           nr1, nr2, nr3
  ! counter on spin polarization
  ! counter on rotations
  ! 2 counters on perturbations
  ! counter on polarization
  ! number of points on this processor
  ! nrx1 size of the FFT mesh array along direction 1
  ! counters on the mesh
  ! counters on the FFT mesh
INTEGER, ALLOCATABLE :: iir(:), jir(:), kir(:), rir(:,:)
  ! the indices of the mesh and the rotated points
COMPLEX(DP) :: dmags(3,npe), mag(3,npe), phase
  ! the magnetization component and the phase
INTEGER :: stilde(3,3,nsym), ftau(3,48)
  ! the symmetry matrices that rotate the magnetization
  ! the fractional translations in the FFT mesh

IF (nsymq == 1) RETURN

nr1=dfftp%nr1
nr2=dfftp%nr2
nr3=dfftp%nr3
nr1x=dfftp%nr1x
my_nrxx=nr1x*dfftp%my_nr2p*dfftp%my_nr3p

ftau(1,1:nsym) = NINT ( ft(1,1:nsym)*dfftp%nr1 ) 
ftau(2,1:nsym) = NINT ( ft(2,1:nsym)*dfftp%nr2 ) 
ftau(3,1:nsym) = NINT ( ft(3,1:nsym)*dfftp%nr3 )

ALLOCATE(iir(my_nrxx))
ALLOCATE(jir(my_nrxx))
ALLOCATE(kir(my_nrxx))
ALLOCATE(rir(my_nrxx, nsymq))

ALLOCATE(phase1(nr1, nsymq))
ALLOCATE(phase2(nr2, nsymq))
ALLOCATE(phase3(nr3, nsymq))
ALLOCATE(zb(nsymq))
ALLOCATE (aux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, npe))
ALLOCATE (aux1(dfftp%nnr, npe))

CALL find_mesh_ijk(my_nrxx, iir, jir, kir)
IF (minus_q) THEN
   CALL rotate_mesh_1s(my_nrxx, s(1,1,irotmq), ftau(1,irotmq), rir) 
   CALL compute_phase(phase1, phase2, phase3, nr1, nr2, nr3, 1, gimq, zb)

   DO is=1, nspin_lsda
      !
      !  collect all the quantity to symmetrize in all processors
      !
      DO ipert = 1, npe
#if defined(__MPI)
         CALL cgather_sym( dfftp, dvsym(:, is, ipert), aux(:, ipert))
#else
         aux(:, ipert)= dvsym(:, is, ipert)
#endif
      ENDDO
      !
      !  symmmetrize. Each processor symmetrizes only the points that it has
      !
      aux1=(0.0_DP, 0.0_DP)
      IF (zb(1)) THEN
         DO ir = 1, my_nrxx
            IF (rir(ir,1)==0) CYCLE
            i=iir(ir)
            j=jir(ir)
            k=kir(ir)
            phase=phase1(i,1)*phase2(j,1)*phase3(k,1)
            DO ipert = 1, npe
               DO jpert = 1, npe
                  aux1(ir, ipert) = aux1(ir, ipert) + tmq(jpert, ipert, irr) * &
                           aux (rir(ir,1), jpert)*phase
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO ir = 1, my_nrxx
            IF (rir(ir,1)==0) CYCLE
            DO ipert = 1, npe
               DO jpert = 1, npe
                  aux1(ir, ipert) = aux1(ir, ipert) + tmq(jpert, ipert, irr) * &
                            aux (rir(ir,1), jpert) 
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      dvsym(:,is,:) = (dvsym(:,is,:) + CONJG(aux1(:,:)))*0.5_DP
   ENDDO
ENDIF

CALL rotate_mesh(my_nrxx, nsymq, rir)
CALL compute_phase(phase1, phase2, phase3, nr1, nr2, nr3, nsymq, gi, zb)

DO is=1, nspin_lsda
   !
   !  collect all the quantities to symmetrize in all processors
   !
   DO ipert = 1, npe
#if defined(__MPI)
      CALL cgather_sym( dfftp, dvsym(:, is, ipert), aux(:, ipert))
#else
      aux(:, ipert)=dvsym(:, is, ipert)
#endif
   ENDDO
   !
   !  symmmetrize. Each processor symmetrizes only the points that it has
   !
   dvsym(:,is,:)=(0.0_DP,0.0_DP)
   DO isym = 1, nsymq
      aux1=(0.0_DP, 0.0_DP)
      IF (zb(isym)) THEN
         DO ir = 1, my_nrxx
            IF (rir(ir,isym)==0) CYCLE
            i=iir(ir)
            j=jir(ir)
            k=kir(ir)
            phase=phase1(i,isym)*phase2(j,isym)*phase3(k,isym)
            DO ipert = 1, npe
               DO jpert = 1, npe
                  aux1 (ir, ipert) = aux1 (ir, ipert) +  &
                           t(jpert, ipert, isym, irr) *  &
                           aux (rir(ir,isym), jpert)  * phase
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO ir = 1, my_nrxx
            IF (rir(ir,isym)==0) CYCLE
            DO ipert = 1, npe
               DO jpert = 1, npe
                  aux1(ir, ipert) = aux1(ir, ipert) +      &
                            t (jpert, ipert, isym, irr) *  &
                            aux (rir(ir,isym), jpert) 
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (t_rev(isym)==1) THEN
         dvsym(:,is,:) = dvsym(:,is,:) + CONJG(aux1(:,:))
      ELSE
         dvsym(:,is,:) = dvsym(:,is,:) + aux1(:,:)
      ENDIF
   ENDDO
ENDDO

DEALLOCATE (aux1)
DEALLOCATE (aux)
!
!  Rotate also the magnetization in the noncollinear magnetic case
!
IF (noncolin.AND.domag) THEN
   ALLOCATE (aux_nc(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, 3, npe))
   !
   !  set the symmetry matrices that rotate the magnetization
   !
   CALL set_stilde(s, stilde, sname, t_rev, invs, nsymq)  
   !
   !  bring the magnetization in crystal coordinates
   !  collect all the quantity to symmetrize in all processors 
   !
   DO ipert = 1, npe
      CALL ccryst_to_cart_t (dfftp%nnr, dvsym(:,2:4,ipert), bg, -1)
      DO is = 2, nspin_mag
#if defined(__MPI)
         CALL cgather_sym( dfftp, dvsym(:, is, ipert), aux_nc(:, is-1, ipert))
#else
         aux_nc(:, is-1, ipert)=dvsym(:, is, ipert)
#endif
         dvsym(:,is,ipert) = (0.0_DP, 0.0_DP)
      ENDDO
   ENDDO
   !
   DO isym = 1, nsymq
      DO ir = 1, my_nrxx
         IF (rir(ir,isym)==0) CYCLE
         dmags=(0.0_DP,0.0_DP)
         IF (zb(isym)) THEN
            i=iir(ir)
            j=jir(ir)
            k=kir(ir)
            phase=phase1(i,isym)*phase2(j,isym)*phase3(k,isym)
            DO ipert = 1, npe
               DO jpert = 1, npe
                  DO is=1,3
                     dmags(is,ipert)=dmags(is,ipert) +     &
                          t(jpert, ipert, isym, irr) *     &
                          aux_nc (rir(ir,isym), is, jpert) * phase
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO ipert = 1, npe
               DO jpert = 1, npe
                  DO is=1,3
                     dmags(is,ipert)=dmags(is,ipert) +     &
                          t(jpert, ipert, isym, irr) *     &
                          aux_nc (rir(ir,isym), is, jpert) 
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
!
! rotate the magnetic moment
!
         DO kpol = 1, 3
            mag(kpol,:) = stilde(1,kpol,isym)*dmags(1,:) + &
                          stilde(2,kpol,isym)*dmags(2,:) + &
                          stilde(3,kpol,isym)*dmags(3,:)
         ENDDO

         IF (t_rev(isym)==1) THEN
            DO kpol=1,3
               dvsym(ir,kpol+1,:) = dvsym(ir,kpol+1,:) + CONJG(mag(kpol,:))
            ENDDO
         ELSE
            DO kpol=1,3
               dvsym(ir,kpol+1,:) = dvsym(ir,kpol+1,:) + mag(kpol,:)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   DEALLOCATE (aux_nc)
!
! go back to cartesian coordinates
!
   DO ipert=1,npe
      CALL ccryst_to_cart_t(dfftp%nnr, dvsym(:,2:4,ipert), at, 1)
   ENDDO
ENDIF

dvsym = dvsym / DBLE(nsymq)

DEALLOCATE(iir)
DEALLOCATE(jir)
dEALLOCATE(kir)
DEALLOCATE(rir)
DEALLOCATE(phase1)
DEALLOCATE(phase2)
DEALLOCATE(phase3)
DEALLOCATE(zb)

RETURN
END SUBROUTINE psymdvscf_tpw

!---------------------------------------------------------------------
SUBROUTINE psyme_tpw (dvsym)
!---------------------------------------------------------------------
!
!  This routine is similar to syme_fpol, but forces the quantity to
!  symmetrize to be real.
!
USE kinds,    ONLY : DP
USE fft_base, ONLY : dfftp
USE noncollin_module, ONLY : nspin_mag

IMPLICIT NONE
COMPLEX(DP) :: dvsym(dfftp%nnr, nspin_mag, 3)

dvsym(:,:,:) = CMPLX(DBLE(dvsym(:,:,:)), 0.0_DP, KIND=DP)

CALL psyme_fpol_tpw(dvsym)

RETURN
END SUBROUTINE psyme_tpw

SUBROUTINE set_stilde(s, stilde, sname, t_rev, invs, nsym)
!
!  This routine sets the matrices needed to rotate the magnetization.
!  These matrices contain the proper part of the inverse of s 
!  (when t_rev==0) or minus the proper part of the inverse of s 
!  (when t_rev==1).
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nsym
INTEGER, INTENT(IN) :: s(3,3,nsym), t_rev(nsym), invs(nsym)
INTEGER, INTENT(INOUT) :: stilde(3,3,nsym)
CHARACTER(LEN=45), INTENT(IN) :: sname(nsym)

INTEGER :: isym

DO isym = 1, nsym
   stilde(:,:,isym)=s(:,:,invs(isym))
   IF (sname(isym)(1:3)=='inv') stilde(:,:,isym)=-stilde(:,:,isym)
   IF (t_rev(isym)==1) stilde(:,:,isym)=-stilde(:,:,isym)
ENDDO

RETURN
END SUBROUTINE set_stilde

END MODULE lr_sym_mod_tpw
