! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  This is a rewriting of the epsilon.x routine of the QE distribution.
!  
! Copyright (C) 2004-2009 Andrea Benassi and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  The formulas are the same, but the present implementation supports
!  NC, US, and PAW pseudopotentials. The commutator of the r operator with
!  the nonlocal PP term is calculated. Moreover it allows the use of 
!  spinors and uses point group symmetry to reduce the number of k points 
!  needed to make the Brillouin zone integrations. It makes sums over empty 
!  states and at convergence with the number of empty bands it should give 
!  the same results of thermo_pw with the flag lnoloc=.TRUE.. 
!  The imaginary part of the inverse of the dielectric constant is written 
!  on file. The calculation of the joint density of states is also supported,
!  but at variance with the original epsilon.f90, this routine does not 
!  support metals. 
!  It can be launched with several images, but only one image is doing
!  the calculations. It is parallelized on G vectors and on pools, but
!  must be launched with the same options used in the thermo_pw.x 
!  calculation. 
!  A typical run is as the epsilon.f90 routine:
!  A run of thermo_pw.x (or pw.x) with some input and what='scf' or 
!  what='scf_dos'.
!  A run of epsilon_tpw.x with the same prefix and the same outdir
!  (with an /g1 extension in the case in which thermo_pw.x has been used).
!  Input variables in the input_epsilon namelist:
!  
!  Mandatory:
!
!   prefix,      ! as in thermo_pw.x (or pw.x) input
!   outdir,      ! the scratch directory (the same as pw.x or outdir/g1 if
                 ! thermo_pw.x is used).
!   
!  The following variables have a default value. You can specify them only
!  to change the plot.
!
!   calculation, ! epsilon or jdos                  Default : epsilon
!   intersmear,  ! the linewidth (in Ry)            Default : 0.01
!   nfs,         ! the number of frequencies        Default : 600
!   wmin,        ! the minimum frequency (in Ry)    Default : 0
!   wmax,        ! the maximum frequency (in Ry)    Default : 2.5
!   nbndmin,     ! the minimum band                 Default : 1
!   nbndmax,     ! the maximum band                 Default : nbnd
!   shift        ! an energy shift of the condution Default : 0.0
                 ! bands (scissor operator) in eV
!
!------------------------------
PROGRAM epsilon_tpw
!------------------------------
  !
  ! Compute the complex macroscopic dielectric function epsilon,
  ! neglecting local field effects.
  !
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : rytoev
  USE io_global,   ONLY : stdout, ionode, ionode_id, meta_ionode, &
                          meta_ionode_id
  USE io_files,    ONLY : tmp_dir, prefix, outdir, iunwfc
  USE ions_base,   ONLY : ntyp => nsp
  USE uspp_param,  ONLY : nhm
  USE uspp,        ONLY : okvan
  USE spin_orb,    ONLY : lspinorb
  USE lsda_mod,    ONLY : nspin
  USE lrus,        ONLY : dpqq, dpqq_so
  USE klist,       ONLY : lgauss
  USE ktetra,      ONLY : ltetra
  USE wvfct,       ONLY : nbnd

  USE mp_world,    ONLY : world_comm
  USE mp,          ONLY : mp_bcast
  USE mp_global,   ONLY : mp_startup, mp_global_end
  USE environment, ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=256) :: fileps
  !
  ! input variables
  !
  INTEGER                 :: nfs, nbndmin, nbndmax
  REAL(DP)                :: intersmear, wmax, wmin, shift
  CHARACTER(10)           :: calculation
  CHARACTER(LEN=9) :: code='epsil_tpw'
  !
  NAMELIST / input_epsilon / prefix,      & ! as in pw.x input
                             outdir,      & ! the scratch directory
                             intersmear,  & ! the linewidth (in Ry)
                             calculation, & ! epsilon or jdos
                             nfs,         & ! the number of frequencies
                             wmin,        & ! the minimum frequency
                             wmax,        & ! the maximum frequency
                             nbndmin,     & ! the minimum band
                             nbndmax,     & ! the maximum band
                             shift,       & ! an energy shift of the condution 
                                            ! bands (scissor operator) in eV
                             fileps         ! the name of the postcript file
                                            ! with the results
  !
  ! local variables
  !
  INTEGER :: ios
  !
  ! initialise environment
  !
  CALL mp_startup ( start_images=.true. )
  CALL environment_start ( code )
  !
  ! Set default values for variables in namelist
  !
  prefix       = 'pwscf'
  outdir       = './'
  calculation  = 'epsilon'
  intersmear   = 0.01_DP
  nfs          = 600
  wmin         = 0.0_DP
  wmax         = 2.5_DP
  nbndmin      = 1
  nbndmax      = 0
  shift        = 0.0_DP
  fileps       = 'output_epsilon.ps'
  !
  ! read input file
  !
  IF (ionode) WRITE( stdout, "( 2/, 5x, 'Reading input file...' ) " )
  ios = 0
  !
  IF (meta_ionode) READ( 5, input_epsilon, IOSTAT = ios )
  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore( 'epsilon', 'reading input_epsilon namelist', ABS( ios ) )
  !
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( outdir, ionode_id, world_comm )
  CALL mp_bcast( intersmear, ionode_id, world_comm)
  CALL mp_bcast( calculation, ionode_id, world_comm)
  CALL mp_bcast( nfs, ionode_id, world_comm )
  CALL mp_bcast( wmax, ionode_id, world_comm )
  CALL mp_bcast( wmin, ionode_id, world_comm )
  CALL mp_bcast( nbndmin, ionode_id, world_comm )
  CALL mp_bcast( nbndmax, ionode_id, world_comm )
  CALL mp_bcast( shift, ionode_id, world_comm )
  CALL mp_bcast( fileps, ionode_id, world_comm )

  tmp_dir = trimcheck(outdir)
  shift = shift / rytoev
  !
  ! read PW simulation parameters from prefix.save/data-file.xml
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )
  !
  CALL read_file()
  CALL openfil_pp()
  !
  IF (lgauss .or. ltetra) THEN
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a metal...' ) " )
      CALL errore('epsilon_tpw','This routine is not working for metals',1)
  ELSE
      IF (ionode) WRITE( stdout, "( 5x, 'The system is an insulator ...' ) " )
  ENDIF

  IF (nbndmax == 0) nbndmax = nbnd
  IF (okvan) THEN
!
!   Determine the dipole of the augmentation functions
!
     ALLOCATE (dpqq( nhm, nhm, 3, ntyp))
     CALL compute_qdipol(dpqq)
     IF (lspinorb) THEN
        ALLOCATE(dpqq_so( nhm, nhm, nspin, 3, ntyp))
        CALL compute_qdipol_so(dpqq, dpqq_so)
     ENDIF
     CALL qdipol_cryst()
  ENDIF
  !
  CALL eps_calc (intersmear, nfs, wmax, wmin, nbndmin, nbndmax, shift, &
                 calculation, fileps)
  !
  IF ( ionode ) WRITE( stdout , "(/)" )
  CALL print_clock( 'eps_calc' )
  CALL print_clock( 'dipole_calc' )
  IF ( ionode ) WRITE( stdout, *  )
  !
  CLOSE ( iunwfc, status='KEEP' )

  CALL environment_end ( code )
  CALL mp_global_end ()

END PROGRAM epsilon_tpw
!
!-----------------------------------------------------------------------------
SUBROUTINE eps_calc (intersmear, nw, wmax, wmin, nbndmin, nbndmax, shift, &
                     calculation, fileps)
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : omega
  USE wvfct,                ONLY : nbnd, et, wg
  USE klist,                ONLY : nks, nkstot, degauss, wk, nelec
  USE symme,                ONLY : symmatrix, crys_to_cart
  USE lsda_mod,             ONLY : current_spin, isk
  !
  USE io_global,            ONLY : ionode, stdout, meta_ionode
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! input variables
  !
  INTEGER,         INTENT(IN)    :: nw, nbndmin, nbndmax
  REAL(DP),        INTENT(IN)    :: wmax, wmin, intersmear, shift
  CHARACTER(LEN=10), INTENT(IN)  :: calculation
  CHARACTER(LEN=256), INTENT(IN) :: fileps
  !
  ! local variables
  !
  INTEGER       :: i, ik, iband1, iband2, ipol, jpol
  INTEGER       :: iw, ierr, iu_epsil
  REAL(DP)      :: etrans, const, w, sumweight
  !
  REAL(DP), ALLOCATABLE    :: epsr(:,:,:), epsi(:,:,:) 
  REAL(DP), ALLOCATABLE    :: eels(:,:,:), jdos(:)
  REAL(DP), ALLOCATABLE    :: wgrid(:)
  COMPLEX(DP),ALLOCATABLE  :: dipole(:,:,:,:), dipole_aux(:,:,:)
  COMPLEX(DP) :: den, aux(3,3)
  REAL(DP) :: deltaw, intjdos
  LOGICAL :: exst
  !
  !  Generate the frequency grid
  !
  CALL start_clock('eps_calc')
  ALLOCATE(wgrid(nw))
  deltaw = (wmax-wmin) / (nw-1)
  DO i=1,nw
     wgrid(i) = wmin + deltaw * (i-1)
  END DO

  IF (calculation=='jdos') THEN
     ALLOCATE(jdos(nw)) 
     jdos(:) = 0.0_DP
  ELSE
     !
     ! allocate main spectral and auxiliary quantities
     !
     ALLOCATE( dipole(3, 3, nbnd, nbnd), STAT=ierr )
     IF (ierr/=0) CALL errore('epsilon_tpw','allocating dipole', abs(ierr) )
     !
     ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
     IF (ierr/=0) CALL errore('epsilon_tpw','allocating dipole_aux', abs(ierr) )
     !
     ALLOCATE(epsr(3, 3, nw)) 
     ALLOCATE(epsi(3, 3, nw))
     ALLOCATE(eels(3, 3, nw)) 
     !
     ! initialize response functions
     !
     epsr(:,:,:) = 0.0_DP
     epsi(:,:,:) = 0.0_DP
     eels(:,:,:) = 0.0_DP
  ENDIF
  !
  ! main kpt loop
  !
  sumweight=0.0_DP
  DO ik = 1, nks
     !
     ! compute the dipole matrix elements
     !
     current_spin=isk(ik)
     IF (calculation == 'jdos' ) THEN
     ELSE
        CALL dipole_calc( ik, dipole_aux, nbndmin, nbndmax, shift)
        !
        DO ipol=1,3
           dipole(ipol,ipol,:,:)= dipole_aux(ipol,:,:) * &
                                      CONJG( dipole_aux(ipol,:,:))
        END DO

        DO ipol=1,3
           DO jpol=ipol+1,3
              dipole(ipol,jpol,:,:)= dipole_aux(ipol,:,:) * &
                                             CONJG( dipole_aux(jpol,:,:))
              dipole(jpol,ipol,:,:) = CONJG(dipole(ipol,jpol,:,:))
           END DO
        END DO
     END IF
     ! 
     ! and compute the dielectric constants or the jdos
     !
     DO iband2 = nbndmin,nbndmax

        DO iband1 = nbndmin,nbndmax
           !
           IF (iband1==iband2) CYCLE
!
!          iband1 is a valence band
!
           IF (ABS(wg(iband1,ik)) >= 1.d-4 * wk(ik) ) THEN
              !
              IF (abs(wg(iband2,ik)-wg(iband1,ik))< 1e-3 * wk(ik) ) CYCLE
              !
              ! transition energy
              !
              sumweight = sumweight + wg(iband1,ik)
              etrans = ( et(iband2,ik) -et(iband1,ik) ) + shift 
              !
              ! loop over frequencies
              !
              DO iw = 1, nw
                 !
                 w = wgrid(iw)
                 !
                 IF (calculation=='jdos') THEN
                    den = etrans - CMPLX(w, intersmear)
                    jdos(iw) = jdos(iw) + AIMAG(wg(iband1,ik)/den)
                 ELSE
                    den = etrans**2 - (CMPLX(w, intersmear))**2 
                    aux(:,:)=wg(iband1,ik)*dipole(:,:,iband1,iband2)*etrans/den
                    epsr(:,:,iw) = epsr(:,:,iw) + 2.0_DP * REAL(aux(:,:))       
                    epsi(:,:,iw) = epsi(:,:,iw) + 2.0_DP * AIMAG(aux(:,:))
                 ENDIF
              END DO
           END IF
        END DO
     END DO
  END DO
  !
  ! recover over kpt parallelization (inter_pool)
  !
  IF (calculation=='jdos') THEN

     CALL mp_sum( jdos, inter_pool_comm )
     const = 1.0_DP / pi / sumweight
     jdos(:) = jdos(:) * const
     intjdos = deltaw * SUM( jdos(:) )
     WRITE(stdout,'(/,5x, "Integral of jdos...",f15.9)') intjdos
     WRITE(stdout,'(/,5x, "Weights...", f15.9, "; ideal...", f15.9  )') &
                                      sumweight, (nbnd -nelec/2.0_DP)* nelec
                     

     IF ( meta_ionode ) THEN
        !
        WRITE(stdout,"(/,5x, 'Writing output on file...' )")
        !
        ! write results on data files
        !
        iu_epsil=2
        INQUIRE(FILE="jdos", exist=exst)
        IF (exst) THEN
           OPEN (UNIT=iu_epsil, FILE='jdos', STATUS='old', &
                                 POSITION='append', FORM='formatted')
        ELSE
           OPEN (UNIT=iu_epsil, FILE='jdos', STATUS='unknown', &
                                                    FORM='formatted')
           WRITE(iu_epsil,'("#  Re(w)       Im(w)          jdos        ")') 
           WRITE(iu_epsil,'("# Frequency in Ry")')
        END IF

        DO iw =1, nw
           WRITE(iu_epsil,'(2f10.5,e15.7)') wgrid(iw), intersmear, jdos(iw)
        END DO
        CLOSE(iu_epsil)
     END IF
     CALL plot_jdos(fileps, nw, wgrid, jdos)
  ELSE
     CALL mp_sum( epsr, inter_pool_comm )
     CALL mp_sum( epsi, inter_pool_comm )
     !
     ! impose the correct normalization
     !
     const = 8.0_DP * pi / omega 
     epsr(:,:,:) = epsr(:,:,:) * const
     epsi(:,:,:) = epsi(:,:,:) * const
     !
     !  Symmetrize
     !
     DO iw = 1, nw
        CALL crys_to_cart ( epsr(:,:,iw) )
        CALL symmatrix ( epsr(:,:,iw) )
        CALL crys_to_cart ( epsi(:,:,iw) )
        CALL symmatrix ( epsi(:,:,iw) )
     END DO
     !
     ! Add one on the diagonal   
     !
     DO ipol=1,3
        epsr(ipol,ipol,:) = 1.0_DP + epsr(ipol,ipol,:)
     END DO
     !
     ! Calculation of eels spectrum
     !
     DO iw = 1, nw
        !
        DO ipol=1,3
           DO jpol=ipol,3
              IF (ABS(epsi(ipol,jpol,iw))> 1.D-10) &
                 eels(ipol,jpol,iw)=epsi(ipol,jpol,iw)/&
                        (epsr(ipol,jpol,iw)**2+epsi(ipol,jpol,iw)**2)
           END DO
        END DO
        !
     END DO
     !
     IF ( meta_ionode ) THEN
        !
        WRITE(stdout,"(/,5x, 'Writing output on file...' )")
        !
        ! write results on data files
        !
        iu_epsil=2
        INQUIRE(FILE="epsilon_re", exist=exst)
        IF (exst) THEN
           OPEN (UNIT=iu_epsil, FILE='epsilon_re', STATUS='old', &
                                 POSITION='append', FORM='formatted')
        ELSE
           OPEN (UNIT=iu_epsil, FILE='epsilon_re', STATUS='unknown', &
                                                    FORM='formatted')
           WRITE(iu_epsil,'("#  Re(w)       Im(w)          e11            e22&
              &            e33            e12            e13            e23")')
           WRITE(iu_epsil,'("# Frequency in Ry")')
        END IF

        DO iw =1, nw
           WRITE(iu_epsil,'(2f10.5,6e15.7)') wgrid(iw), intersmear,       &
                   epsr(1,1,iw), epsr(2,2,iw), epsr(3,3,iw), &
                   epsr(1,2,iw), epsr(1,3,iw), epsr(2,3,iw)
        END DO
        CLOSE(iu_epsil)
 
        iu_epsil=2
        INQUIRE(FILE="epsilon_im", exist=exst)
        IF (exst) THEN
           OPEN (UNIT=iu_epsil, FILE='epsilon_im', STATUS='old', &
                                 POSITION='append', FORM='formatted')
        ELSE
           OPEN (UNIT=iu_epsil, FILE='epsilon_im', STATUS='unknown', &
                                                    FORM='formatted')
           WRITE(iu_epsil,'("#  Re(w)       Im(w)          e11            e22&
                &            e33            e12            e13            e23")')
           WRITE(iu_epsil,'("# Frequency in Ry")')
        END IF

        DO iw =1, nw
           WRITE(iu_epsil,'(2f10.5,6e15.7)') wgrid(iw), intersmear,   &
                  epsi(1,1,iw), epsi(2,2,iw), epsi(3,3,iw),  &
                  epsi(1,2,iw), epsi(1,3,iw), epsi(2,3,iw)
        END DO
        CLOSE(iu_epsil)

        iu_epsil=2
        INQUIRE(FILE="epsilonm1_im", exist=exst)
        IF (exst) THEN
           OPEN (UNIT=iu_epsil, FILE='epsilonm1_im', STATUS='old', &
                                 POSITION='append', FORM='formatted')
        ELSE
           OPEN (UNIT=iu_epsil, FILE='epsilonm1_im', STATUS='unknown', &
                                                    FORM='formatted')
           WRITE(iu_epsil,'("#  Re(w)       Im(w)        em111          em122&
                &          em133          em112          em113          em123")')
           WRITE(iu_epsil,'("# Frequency in Ry")')
        END IF

        DO iw =1, nw
           WRITE(iu_epsil,'(2f10.5,6e15.7)') wgrid(iw), intersmear,   &
                  eels(1,1,iw), eels(2,2,iw), eels(3,3,iw),  &
                  eels(1,2,iw), eels(1,3,iw), eels(2,3,iw)
        END DO
        CLOSE(iu_epsil)
        !
     ENDIF
     CALL plot_epsilon(fileps, nw, wgrid, epsr, epsi, eels)

     DEALLOCATE (epsr) 
     DEALLOCATE (epsi) 
     DEALLOCATE (eels) 
     !
     DEALLOCATE (dipole) 
     DEALLOCATE (dipole_aux)
  ENDIF
  DEALLOCATE (wgrid)

  CALL stop_clock('eps_calc')
  RETURN
END SUBROUTINE eps_calc

!--------------------------------------------------------------------
SUBROUTINE dipole_calc( ik, dipole_aux, nbndmin, nbndmax, shift )
  !------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi, degspin
  USE wvfct,                ONLY : wg, npw, nbnd, igk, g2kin, npwx, et
  USE gvecw,                ONLY : ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE ener,                 ONLY : ef
  USE klist,                ONLY : xk, wk, nelec, ngk, igk_k
  USE noncollin_module,     ONLY : noncolin, npol
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE becmod,               ONLY : bec_type, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum

  IMPLICIT NONE
  !
  ! global variables
  INTEGER, INTENT(IN)        :: ik, nbndmin, nbndmax
  COMPLEX(DP), INTENT(INOUT) :: dipole_aux(3,nbnd,nbnd)
  REAL(DP), INTENT(IN)       :: shift
  !
  ! local variables
  !
  TYPE(bec_type) :: becp1, becp2  ! the scalar products between wavefunctions
                                  ! and projectors
  INTEGER :: iband1, iband2, ig, nbnd_occ, ipol, ibnd
  !
  COMPLEX(DP)   :: ZDOTC
  COMPLEX(DP), ALLOCATABLE :: dpsi(:,:), dvpsi(:,:)
  REAL(DP) :: small, xmax, fac, targete, etrans
  !
  CALL start_clock( 'dipole_calc' )
  !
  ALLOCATE ( dpsi ( npwx*npol, nbnd))
  !
  IF (okvan) ALLOCATE ( dvpsi ( npwx*npol, nbnd))
  !
  ! setup k+G grids for each kpt
  !
  CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  !
  !  should be in readfile but it is not there
  !
  ngk(ik) = npw
  !
  igk_k(1:npw,ik) = igk(1:npw)
  !
  ! read wfc for the given kpt
  !
  CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
  !
  CALL init_us_2 (npw, igk, xk (1, ik), vkb)
  !
  CALL allocate_bec_type ( nkb, nbnd, becp1)
  !
  CALL calbec (npw, vkb, evc, becp1 )
  !
  CALL allocate_bec_type ( nkb, nbnd, becp2)
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)

  IF (noncolin) THEN
     nbnd_occ = NINT (nelec) 
  ELSE
     nbnd_occ = NINT (nelec) / degspin
  ENDIF

  DO ipol=1,3 
     !
     CALL commutator_Hx_psi (ik, nbnd_occ, becp1, becp2, ipol, dpsi )
     IF (okvan) THEN
        dvpsi=(0.0_DP,0.0_DP)
        CALL adddvepsi_us_tpw( nbnd_occ, becp1, becp2, ipol, ik, dvpsi)
     ENDIF
     !
     DO iband2 = nbndmin, nbndmax
        !
        !  sum on the conduction bands and a few valence in metals
        !
        IF ( wg(iband2,ik) < wk(ik)  ) THEN
           DO iband1 = nbndmin, nbndmax
              !
              ! sum over valence band
              !
              IF ( ABS(wg(iband1,ik)) >= 1.d-4 ) THEN
                 dipole_aux(ipol,iband1,iband2)= &
                              ZDOTC(npw, evc(1,iband2),1,dpsi(1,iband1),1)
                 IF (noncolin) &
                    dipole_aux(ipol,iband1,iband2)= &
                              dipole_aux(ipol,iband1,iband2) +  &
                              ZDOTC(npw, evc(1+npwx,iband2),1,  &
                                                     dpsi(1+npwx,iband1),1)
              END IF
              etrans = ( et(iband2,ik) -et(iband1,ik) ) + shift 
              dipole_aux(ipol,iband1,iband2)=dipole_aux(ipol,iband1,iband2) &
                                 / etrans
              IF (okvan) THEN
                 dipole_aux(ipol,iband1,iband2)= &
                              dipole_aux(ipol,iband1,iband2) + &
                              ZDOTC(npw, evc(1,iband2),1,dvpsi(1,iband1),1)
                 IF (noncolin) &
                    dipole_aux(ipol,iband1,iband2)= &
                              dipole_aux(ipol,iband1,iband2) +  &
                              ZDOTC(npw, evc(1+npwx,iband2),1,  &
                                                     dvpsi(1+npwx,iband1),1)
              ENDIF
           END DO
        END IF
     END DO
  END DO
  !
  ! recover over G parallelization (intra_bgrp)
  !
  CALL mp_sum( dipole_aux, intra_bgrp_comm )
  !
  CALL deallocate_bec_type (becp1)
  CALL deallocate_bec_type (becp2)
  DEALLOCATE(dpsi)
  IF (okvan) DEALLOCATE(dvpsi)
  !
  CALL stop_clock( 'dipole_calc' )
  !
END SUBROUTINE dipole_calc

SUBROUTINE plot_jdos(fileps, nw, wgrid, jdos )
!
USE kinds, ONLY : DP
USE constants, ONLY : rytoev
USE gnuplot,   ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                      gnuplot_xlabel, gnuplot_ylabel, gnuplot_write_command, &
                      gnuplot_write_file_mul_data
USE io_global, ONLY : ionode
USE mp_images, ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: fileps

INTEGER, INTENT(IN) :: nw
REAL(DP) :: wgrid(nw), jdos(nw)
REAL(DP) :: ymin, ymax
INTEGER  :: iw, ierr

CHARACTER(LEN=256) :: gnu_filename, filename, ylabel, xlabel

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_tmp_epsilon'
CALL gnuplot_start(gnu_filename)

ymin=1.D10
ymax=0.0_DP
DO iw=1,nw
   IF (jdos(iw) > ymax) ymax=jdos(iw)
   IF (jdos(iw) < ymin) ymin=jdos(iw)
END DO
ymax=ymax*1.1_DP

IF (TRIM(fileps)=='output_epsilon.ps') THEN
   filename='output_jdos.ps'
ELSE
   filename=TRIM(fileps)
ENDIF

CALL gnuplot_write_header(filename, wgrid(1), wgrid(nw), ymin, ymax, rytoev )

xlabel='{/Symbol w} (eV)'
ylabel='joint-dos (states/ spin / (Ry N_c N_v ))'
CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
CALL gnuplot_xlabel(TRIM(xlabel), .FALSE.)
CALL gnuplot_write_command('plot_width=2',.FALSE.)

CALL gnuplot_write_file_mul_data('jdos',1,3,'color_red',.TRUE.,.TRUE., .FALSE.)

CALL gnuplot_end()

IF (ionode) &
   ierr=system('gnuplot '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_jdos

SUBROUTINE plot_epsilon(fileps, nw, wgrid, epsr, epsi, eels )
!
USE kinds, ONLY : DP
USE constants, ONLY : rytoev
USE gnuplot,   ONLY : gnuplot_start, gnuplot_end, gnuplot_write_header, &
                      gnuplot_xlabel, gnuplot_ylabel, gnuplot_write_command, &
                      gnuplot_write_file_mul_data
USE io_global, ONLY : ionode
USE mp_images, ONLY : my_image_id, root_image

IMPLICIT NONE
CHARACTER(LEN=256), INTENT(IN) :: fileps

INTEGER, INTENT(IN) :: nw
REAL(DP) :: wgrid(nw), epsr(3,3,nw), epsi(3,3,nw), eels(3,3,nw)
REAL(DP) :: ymin, ymax
INTEGER  :: iw, ierr

CHARACTER(LEN=256) :: gnu_filename, filename, ylabel, xlabel

IF ( my_image_id /= root_image ) RETURN

gnu_filename='gnuplot_tmp_epsilon'
CALL gnuplot_start(gnu_filename)

filename=TRIM(fileps)

CALL gnuplot_write_header(filename,wgrid(1),wgrid(nw),0.0_DP,0.0_DP,rytoev)

xlabel='{/Symbol w} (eV)'
ylabel='{/Symbol e}_1 ({/Symbol w})'
CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
CALL gnuplot_xlabel(TRIM(xlabel), .FALSE.)
CALL gnuplot_write_command('plot_width=2',.FALSE.)

CALL gnuplot_write_file_mul_data('epsilon_re',1,3,'color_blue',.TRUE.,.FALSE., .FALSE.)
CALL gnuplot_write_file_mul_data('epsilon_re',1,4,'color_green',.FALSE.,.FALSE., .FALSE.)
CALL gnuplot_write_file_mul_data('epsilon_re',1,5,'color_red',.FALSE.,.TRUE., .FALSE.)

ylabel='{/Symbol e}_2 ({/Symbol w})'
CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
CALL gnuplot_write_command('plot_width=2',.FALSE.)

CALL gnuplot_write_file_mul_data('epsilon_im',1,3,'color_blue',.TRUE.,.FALSE., .FALSE.)
CALL gnuplot_write_file_mul_data('epsilon_im',1,4,'color_green',.FALSE.,.FALSE., .FALSE.)
CALL gnuplot_write_file_mul_data('epsilon_im',1,5,'color_red',.FALSE.,.TRUE., .FALSE.)

ylabel='Im 1/{/Symbol e} ({/Symbol w})'
CALL gnuplot_ylabel(TRIM(ylabel), .FALSE.)
CALL gnuplot_write_command('plot_width=2',.FALSE.)

CALL gnuplot_write_file_mul_data('epsilonm1_im',1,3,'color_blue',.TRUE.,.FALSE., .FALSE.)
CALL gnuplot_write_file_mul_data('epsilonm1_im',1,4,'color_green',.FALSE.,.FALSE., .FALSE.)
CALL gnuplot_write_file_mul_data('epsilonm1_im',1,5,'color_red',.FALSE.,.TRUE., .FALSE.)

CALL gnuplot_end()

IF (ionode) &
   ierr=system('gnuplot '//TRIM(gnu_filename))

RETURN
END SUBROUTINE plot_epsilon
