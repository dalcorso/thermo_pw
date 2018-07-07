!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE plan_avg_sub(averag, vacuum, nat_, nbnd_, nks_, ninter, &
           surface1, surface2 )
!-----------------------------------------------------------------------
  !
  ! calculate planar averages of each wavefunction
  ! on output this routine sets the array averag with the planar average
  ! integrated on each layer of each band at each k point
  ! the array vacuum sets the integral of the state on vacuum.
  ! In the noncollinear case avegar contains also the planar average
  ! in each layer of the magnetization densities, in the spin 2:4
  !
  !
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE klist,     ONLY : nkstot
  USE cell_base, ONLY : ibrav, celldm
  USE ions_base, ONLY : nat, ityp, ntyp => nsp, tau, atm
  USE wvfct,     ONLY : nbnd
  USE lsda_mod,  ONLY : nspin
  USE control_flags, ONLY : gamma_only
  USE control_2d_bands, ONLY : dump_states
  USE becmod,    ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE uspp,      ONLY : nkb
  USE klist,     ONLY : nks, xk, nkstot
  USE wrappers,  ONLY : f_mkdir_safe
  USE io_global, ONLY : ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp_bands,  ONLY : nbgrp, me_bgrp, root_bgrp, intra_bgrp_comm
  USE mp_pools,  ONLY : me_pool, root_pool, intra_pool_comm
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !

  INTEGER, INTENT(IN) :: nat_, nbnd_, nks_
  REAL(DP), INTENT(INOUT) :: averag( nat_, nspin, nbnd_, nks_)
  REAL(DP), INTENT(INOUT) :: vacuum( nspin, nbnd_, nks_)
  INTEGER, INTENT(OUT) :: ninter, surface1, surface2
  !
  INTEGER :: ik_index, i1(nat), vacuum1, vacuum2
  REAL(DP) :: zdim
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=6) :: int_to_char
  !
  REAL(DP), ALLOCATABLE :: plan (:,:,:)
  INTEGER :: ir, na, ibnd, ik, iun, ios
  INTEGER :: find_free_unit
  !
  !   Now clean completely pw, allocate space for pwscf variables, 
  !   read and check them.
  !
  IF (nks_ /= nkstot) CALL errore('plan_avg_sub','inconsistent k points',1)
  CALL close_files(.TRUE.)
  CALL clean_pw(.TRUE.)
  CALL read_file ( )
  IF (nat_ /= nat .OR. nbnd_ /= nbnd) &
                     CALL errore('plan_avg_sub','some problems',1)
  !
  IF (gamma_only) CALL errore ('plan_avg_sub', &
       ' planar average with gamma tricks not yet implemented',2)
  !
  CALL openfil_pp ( )
  !
  ALLOCATE(plan(dfftp%nr3, nspin, nbnd))
  !
  CALL prepare_plan_avg(ninter, zdim, i1, vacuum1, vacuum2, surface1, surface2)
  averag(:,:,:,:) = 0.d0
  vacuum(:,:,:) = 0.d0
  ios=0
  IF (dump_states.AND.ionode) ios = f_mkdir_safe( 'dump' )
  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  IF (ios > 0) CALL errore('plan_avg_sub','problems opening dump directory',1)
  !
  IF (dump_states) THEN
     IF (ionode) THEN
        iun=find_free_unit()
        OPEN(UNIT=iun,FILE='dump/dump_states',STATUS='unknown',ERR=400,&
                                                              IOSTAT=ios)
        WRITE(iun,'(i5)') ibrav
        WRITE(iun,'(6f15.8)') celldm 
        WRITE(iun,'(i8)') nat
        DO na=1, nat
           WRITE(iun,'(a3,3f20.10)') atm(ityp(na)), tau(:,na)
        END DO
        WRITE(iun, '(4i8)') dfftp%nr3, nbnd, nkstot, nspin
        CLOSE(iun) 
     END IF
  END IF
400  CALL mp_bcast(ios,ionode_id,intra_image_comm)
     IF (ios /= 0) CALL errore('plan_avg_sub','problems with dump of states',1)

  CALL allocate_bec_type ( nkb, nbnd, becp )
  ! create the directory dump with all states
  IF (nbgrp > 1.AND.dump_states) &
      CALL errore('plan_avg_sub','dump state not implemented with band &
                                  &parallization',1)
  DO ik=1,nks
     plan(:,:,:) = 0.d0
     CALL do_plan_avg (ik, averag, vacuum, plan, zdim, ninter, i1, &
                                       vacuum1, vacuum2, surface1, surface2 )
!
!   to be checked if there is band parallelization
!
     IF (dump_states.AND.me_pool==root_pool) THEN
        iun=find_free_unit()
        ik_index = find_global_ik(ik)
        filename='dump/state_k_'//TRIM(int_to_char(ik_index))
        OPEN(UNIT=iun,FILE=TRIM(filename),STATUS='unknown',ERR=300,&
                                                              IOSTAT=ios)
        DO ibnd=1, nbnd
           WRITE(iun,'(2i8)') ik_index, ibnd
           DO ir=1, dfftp%nr3
              WRITE(iun,'(i8,4f12.7)') ir, plan(ir, 1:nspin, ibnd)
           ENDDO
        ENDDO
        CLOSE(iun)
     ENDIF
300     CALL mp_bcast(ios,root_pool,intra_pool_comm)
        CALL mp_bcast(ios,root_bgrp,intra_bgrp_comm)
     IF (ios /= 0) CALL errore('plan_avg_sub','problems with dump of states',1)
  ENDDO

  CALL deallocate_bec_type (becp)
  DEALLOCATE(plan)

  CALL poolrecover (averag, nat * nbnd * nspin, nkstot, nks)
  CALL poolrecover (vacuum, nbnd * nspin, nkstot, nks)
  CALL poolrecover (xk, 3, nkstot, nks)

  RETURN

CONTAINS
!
SUBROUTINE prepare_plan_avg (ninter, zdim, i1, vacuum1, vacuum2, surface1, surface2)
  !
  !  This routine prepare the calculation of the planar average by computing
  !  1) ninter   the number of layers of the surface
  !  2) i1       the starting point of each layer on the fft mesh 
  !  3) vacuum1, vacuum2 ! the initial and final point of vacuum on the fft mesh
  !  4) surface1, surface2 ! the layer in which the surface starts
  !
  !
  USE cell_base, ONLY: celldm, alat, tpiba2
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp, tau
  USE control_2d_bands, ONLY : sp_min
  USE io_global, ONLY : stdout

  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ninter, i1(nat), vacuum1, vacuum2, surface1, surface2
  ! output: see above
  REAL(DP), INTENT(OUT) :: zdim
  ! output: lenght of the cell in a.u.
  !
  !      Local variables
  !
  INTEGER :: ik, ibnd, iin, na, ir, ij, ind, ntau (nat + 1), ind1, &
             ind2, i3, max_dist
  ! counter on k points
  ! counter on bands
  ! counter on planes
  ! counter on atoms
  ! counter on points
  ! counter on coordinates and planes
  ! starting point of each plane
  ! the number of tau per plane
  INTEGER, ALLOCATABLE :: distance(:)

  REAL(DP) :: avg (nat), z1 (nat)
  ! the average position of each plane
  ! auxiliary for coordinates

  IF ( celldm(3) == 0.d0 ) celldm(3) = celldm(1)
  zdim = alat * celldm (3)
  !
  !     Compute the number of planes and the coordinates on the mesh of the
  !     points which define each plane
  !
  avg(:) = 0.d0
  ninter = 1
  z1 (ninter) = tau (3, 1)
  avg (ninter) = tau (3, 1)
  ntau (ninter) = 1
  DO na = 2, nat
     DO iin = 1, ninter
        IF (abs (mod (z1(iin)-tau(3,na), celldm(3)) ) < sp_min) THEN
           avg (iin) = avg (iin) + tau (3, na)
           ntau (iin) = ntau (iin) + 1
           GOTO 100
        ENDIF
     ENDDO
     ninter = ninter + 1
     z1 (ninter) = tau (3, na)
     avg (ninter) = tau (3, na)
     ntau (ninter) = 1
100  CONTINUE
  ENDDO
  !
  !     for each plane compute the average position of the central plane
  !     and first point in the fft mesh
  !
  DO iin = 1, ninter
     z1 (iin) = mod (avg (iin), celldm (3) ) / ntau (iin)
     ind = (z1 (iin) / celldm (3) ) * dfftp%nr3 + 1
     IF (ind<=0) ind = ind + dfftp%nr3
     i1 (iin) = ind
  ENDDO
  !
  !    order the points
  !
  DO iin = 1, ninter
     ntau (iin) = i1 (iin)
     DO ik = iin + 1, ninter
        IF (i1 (ik) <ntau (iin) ) THEN
           ij = ntau (iin)
           ntau (iin) = i1 (ik)
           i1 (ik) = ij
        ENDIF
     ENDDO
  ENDDO
  ntau (ninter + 1) = ntau (1) + dfftp%nr3
  !
  !    and compute the point associated to each layer
  !
  DO iin = 1, ninter
     i1 (iin) = (ntau (iin) + ntau (iin + 1) ) / 2
  ENDDO
  !
  !  Identify the two surface layers: they are the two consecutive layers
  !  that have the largest lenght. Distances are the lenght of two consecutive
  !  layers, calculated using periodic boundary conditions
  !
  ALLOCATE(distance(ninter))
  DO iin=2, ninter-1
     distance(iin)= i1(iin+1)-i1(iin-1)
  ENDDO
  distance(1)= i1(2) + dfftp%nr3 - i1(ninter)
  distance(ninter)= i1(1) + dfftp%nr3 - i1(ninter-1) 
  max_dist=0
  DO iin=1, ninter
     IF (distance(iin) > max_dist) THEN
        max_dist=distance(iin)
        surface1=iin
     ENDIF
  ENDDO
  surface2=surface1+1
  IF (surface2 > ninter) surface2=1
  !
  !  Vacuum starts at i1(surface1-1) + the lenght of the layer surface1-1 and
  !  ends at i1(surface2+1) + the lenght of the layer surface2+1
  ! 
  ind1=surface1-1
  IF (ind1<1) ind1=ind1+ninter
  ind2=surface1-2
  IF (ind2<1) ind2=ind2+ninter
  i3 = i1(ind1)-i1(ind2)
  IF (i3 < 1) i3=i3 + dfftp%nr3
  vacuum1=i1(ind1) + 2*i3
  IF (vacuum1 > dfftp%nr3 ) vacuum1=vacuum1-dfftp%nr3
  ind1=surface2+1
  IF (ind1>ninter) ind1=ind1-ninter
  i3 = i1(ind1) - i1(surface2)
  IF (i3 < 1) i3=i3 + dfftp%nr3
  vacuum2=i1(surface2) -  2*i3
  IF (vacuum2 < 1) vacuum2=vacuum2 + dfftp%nr3
  
  WRITE(stdout,'(/,5x,"Computing the average charge of each state on each &
                               &layer")')
  WRITE(stdout,'(5x,"Found ",i5," layers. FFT nr3 is", i5)') ninter,dfftp%nr3  
  WRITE(stdout,'(5x, "Layer number    starts    ends")') 
  iin=1
  WRITE(stdout,'(5x,3i9)') iin, i1(ninter), i1(iin)-1
  DO iin=1, ninter-1
     WRITE(stdout,'(5x,3i9)') iin+1, i1(iin), i1(iin+1)-1
  ENDDO
  WRITE(stdout,'(5x,"vacuum",2i9)') vacuum1, vacuum2
  WRITE(stdout,*)

  WRITE(stdout,'(5x, "The points of the fft mesh closest &
                                          &to each atom have i3 equal to")') 
  WRITE(stdout,'(5x, "Atom number     i3       tau_z")')
  DO na=1,nat
     ind = (tau (3,na) / celldm (3) ) * dfftp%nr3 + 1
     IF (ind<=0) ind = ind+dfftp%nr3
     WRITE(stdout,'(5x,2i9,f15.7)') na, ind, tau(3,na)
  END DO
  WRITE(stdout,*)

  DEALLOCATE(distance)

  RETURN
END SUBROUTINE prepare_plan_avg

SUBROUTINE do_plan_avg (ik, averag, vacuum, plan, zdim, ninter, i1, &
                                    vacuum1, vacuum2, surface1, surface2)
  !
  !    This routine computes the planar average on the xy plane
  !    for the charge density of all state of the system at ik.
  !    The routine should work on parallel machines.
  !    The index ik is the index inside the pool, and in principle each pool call
  !    this routine for its own k points and then all the results are collected
  !    In the US case the augmentation part is added only in one
  !    dimension, so that no overload with respect to the NC case
  !    is expected.
  !
  !    Furthermore the amount of charge contained in each plane is
  !    evaluated and given as output. The number of planes is
  !    computed starting from the atomic positions
  !
  USE cell_base, ONLY: tpiba2
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp, tau
  USE gvect,  ONLY : g, ngm
  USE klist, ONLY: nks, nkstot, xk
  USE lsda_mod, ONLY: nspin, lsda, current_spin, isk
  USE uspp, ONLY: vkb, nkb
  USE wvfct, ONLY: npwx, nbnd, wg, g2kin
  USE gvecw, ONLY : ecutwfc
  USE klist, ONLY : ngk, igk_k
  USE wavefunctions_module,  ONLY: evc
  USE noncollin_module, ONLY : noncolin, npol
  USE io_files, ONLY: iunwfc, nwordwfc
  USE io_global, ONLY : stdout
  USE becmod, ONLY: becp, calbec

  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ninter, surface1, surface2
  ! output: the number of planes
  REAL(DP), INTENT(INOUT) :: averag (nat, nspin, nbnd, nkstot), &
                             vacuum (nspin, nbnd, nkstot), &
                             plan (dfftp%nr3, nspin, nbnd)

  ! output: the average charge on ea
  ! output: the planar average
  !
  !      Local variables
  !
  INTEGER :: ik, ibnd, iin, na, ir, ij, ind, i1 (nat), ntau (nat + 1), ind1, &
                       ind2, i3, vacuum1, vacuum2, max_dist, ispin, npw
  INTEGER, ALLOCATABLE :: distance(:)
  ! counter on k points
  ! counter on bands
  ! counter on planes
  ! counter on atoms
  ! counter on points
  ! counter on coordinates and planes
  ! starting point of each plane
  ! the number of tau per plane

  REAL(DP) :: sp_min, avg (nat), z1 (nat), zdim
  ! minimum plane distance
  ! the average position of each plane
  ! auxiliary for coordinates
  ! length in a.u. of the cell along z
  !
  !     for each state compute the planar average
  !
  IF (lsda) current_spin = isk (ik)
  npw=ngk(ik)
!  CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
  CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)

  CALL calbec ( npw, vkb, evc, becp)

  DO ibnd = 1, nbnd
     CALL local_dos1d_so (ik, ibnd, plan (1, 1, ibnd) )
     !
     !     compute the integrals of the charge
     !
     DO ispin=1,nspin
        DO ir = 1, i1 (1) - 1
           averag (1, ispin, ibnd, ik) = averag (1, ispin, ibnd, ik) &
                                       + plan (ir, ispin, ibnd)
        ENDDO
        DO ir = i1 (ninter), dfftp%nr3
           averag (1, ispin, ibnd, ik) = averag (1, ispin, ibnd, ik) &
                                       + plan (ir, ispin, ibnd)
        ENDDO
        averag (1, ispin, ibnd, ik) = averag (1, ispin, ibnd, ik) * zdim &
                                                                  / dfftp%nr3
        DO iin = 2, ninter
           DO ir = i1 (iin - 1), i1 (iin) - 1
              averag(iin,ispin,ibnd,ik) = averag(iin,ispin,ibnd,ik) + &
                         plan(ir,ispin,ibnd)
           ENDDO
           averag (iin, ispin, ibnd, ik) = averag (iin, ispin, ibnd, ik) &
                                                 * zdim / dfftp%nr3
        ENDDO
        IF (vacuum1 < vacuum2) THEN
!
!   standard case the vacuum is in the center of the cell 
!
           DO ir=vacuum1, vacuum2
              vacuum(ispin,ibnd,ik) = vacuum(ispin,ibnd,ik) + &
                                      plan(ir,ispin,ibnd)
           ENDDO
        ELSE
!
!   opposite case the slab is in the center of the cell
!
           DO ir=vacuum1, dfftp%nr3
              vacuum(ispin,ibnd,ik) = vacuum(ispin,ibnd,ik) + &
                                      plan(ir,ispin,ibnd)
           ENDDO
           DO ir=1,vacuum2
              vacuum(ispin,ibnd,ik) = vacuum(ispin,ibnd,ik) + &
                                      plan(ir,ispin,ibnd)
           ENDDO
        ENDIF
        vacuum(ispin,ibnd,ik) = vacuum(ispin,ibnd,ik) * zdim / dfftp%nr3
     END DO
  END DO
  !
  RETURN
  !
END SUBROUTINE do_plan_avg


INTEGER FUNCTION find_global_ik(ik)
  USE klist,  ONLY : nkstot
  USE mp_pools, ONLY : my_pool_id, npool, kunit
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik
  INTEGER :: nks, rest, nbase

  nks  = kunit * ( nkstot / kunit / npool )
  rest = ( nkstot - nks * npool ) / kunit
  IF ( ( my_pool_id + 1 ) <= rest ) nks = nks + kunit
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit

  find_global_ik = nbase + ik

  RETURN
END FUNCTION find_global_ik

END SUBROUTINE plan_avg_sub

