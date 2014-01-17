! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE matdyn_sub(do_dos, igeom)
  !-----------------------------------------------------------------------
  !  this program calculates the phonon frequencies for a list of generic
  !  q vectors starting from the interatomic force constants generated
  !  from the dynamical matrices as written by DFPT phonon code through
  !  the companion program q2r
  !
  !  matdyn can generate a supercell of the original cell for mass
  !  approximation calculation. If supercell data are not specified
  !  in input, the unit cell, lattice vectors, atom types and positions
  !  are read from the force constant file
  !
  !  Input cards: namelist &input
  !     flfrc     file produced by q2r containing force constants (needed)
  !               It is the same as in the input of q2r.x (+ the .xml extension
  !               if the dynamical matrices produced by ph.x were in xml
  !               format). No default value: must be specified.
  !      asr      (character) indicates the type of Acoustic Sum Rule imposed
  !               - 'no': no Acoustic Sum Rules imposed (default)
  !               - 'simple':  previous implementation of the asr used
  !                  (3 translational asr imposed by correction of
  !                  the diagonal elements of the force constants matrix)
  !               - 'crystal': 3 translational asr imposed by optimized
  !                  correction of the force constants (projection).
  !               - 'one-dim': 3 translational asr + 1 rotational asr
  !                  imposed by optimized correction of the force constants
  !                  (the rotation axis is the direction of periodicity;
  !                   it will work only if this axis considered is one of
  !                   the cartesian axis).
  !               - 'zero-dim': 3 translational asr + 3 rotational asr
  !                  imposed by optimized correction of the force constants
  !               Note that in certain cases, not all the rotational asr
  !               can be applied (e.g. if there are only 2 atoms in a
  !               molecule or if all the atoms are aligned, etc.).
  !               In these cases the supplementary asr are cancelled
  !               during the orthonormalization procedure (see below).
  !     flvec     output file for normalized phonon displacements 
  !               (default: 'matdyn.modes'). The normalized phonon displacements
  !               are the eigenvectors divided by the mass and then normalized.
  !               As such they are not orthogonal.
  !              
  !     at        supercell lattice vectors - must form a superlattice of the
  !               original lattice (default: use original cell)
  !     ntyp      number of atom types in the supercell (default: ntyp of the
  !               original cell)
  !     amass     masses of atoms in the supercell (a.m.u.), one per atom type
  !               (default: use masses read from file flfrc)
  !  If q = 0, the direction qhat (q=>0) for the non-analytic part
  !  is extracted from the sequence of q-points as follows:
  !     qhat = q(n) - q(n-1)  or   qhat = q(n) - q(n+1)
  !  depending on which one is available and nonzero.
  !  For low-symmetry crystals, specify twice q = 0 in the list
  !  if you want to have q = 0 results for two different directions
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm, my_image_id, root_image
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_ifc_param, read_ifc
  USE cell_base,  ONLY : at, bg, celldm
  USE constants,  ONLY : RY_TO_THZ, RY_TO_CMM1, amu_ry
  USE symm_base,  ONLY : set_sym
  USE rap_point_group,  ONLY : code_group

  USE ifc,        ONLY : frc, atm, zeu, tau_blk, ityp_blk, m_loc, &
                         disp_q, disp_nqs, nq1_d, nq2_d, &
                         nq3_d, deltafreq, freqmin, freqmax, ndos, zasr, &
                         freqmin_input, freqmax_input
 
  USE thermodynamics, ONLY : phdos_save
  USE control_thermo, ONLY : flfrc, flfrq, fldos, ldos
  USE ions_base, ONLY : amass
  USE phdos_module, ONLY : set_phdos
  !
  IMPLICIT NONE
  !
  ! variables *_blk refer to the original cell, other variables
  ! to the (super)cell (which may coincide with the original cell)
  !
  LOGICAL, INTENT(IN) :: do_dos
  INTEGER, INTENT(IN) :: igeom
  INTEGER :: nax, nax_blk
  INTEGER, PARAMETER:: ntypx=10, nrwsx=200
  REAL(DP), PARAMETER :: eps=1.0d-6
  INTEGER :: nr1, nr2, nr3, nsc, ntetra, ibrav
  CHARACTER(LEN=256) :: flvec, filename
  CHARACTER(LEN=10)  :: asr
  LOGICAL :: has_zstar
  COMPLEX(DP), ALLOCATABLE :: dyn(:,:,:,:), dyn_blk(:,:,:,:), frc_ifc(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: z(:,:)
  REAL(DP), ALLOCATABLE:: tau(:,:), q(:,:), w2(:,:), freq(:,:)
  INTEGER, ALLOCATABLE:: tetra(:,:), ityp(:), itau_blk(:)
  REAL(DP) ::     omega,alat, &! cell parameters and volume
                  at_blk(3,3), bg_blk(3,3),  &! original cell
                  omega_blk,                 &! original cell volume
                  epsil(3,3),                &! dielectric tensor
                  amass_blk(ntypx),          &! original atomic masses
                  atws(3,3),      &! lattice vector for WS initialization
                  rws(0:3,nrwsx)   ! nearest neighbor list, rws(0,*) = norm^2
  !
  INTEGER :: nat, nat_blk, ntyp, ntyp_blk, &
             nrws,                         & ! number of nearest neighbor
             code_group_old

  INTEGER :: nspin_mag, nqs, ios
  !
  LOGICAL :: xmlifc, lo_to_split, na_ifc, fd
  !
  REAL(DP) :: qhat(3), qh, e, emin, emax, dosofe(2), qq
  REAL(DP) :: delta, pathL
  REAL(DP), ALLOCATABLE :: xqaux(:,:)
  INTEGER, ALLOCATABLE :: nqb(:)
  INTEGER :: n, i, j, it, nq, nqx, na, nb, iout, nqtot
  LOGICAL, EXTERNAL :: has_xml
  CHARACTER(LEN=15), ALLOCATABLE :: name_rap_mode(:)
  INTEGER, ALLOCATABLE :: num_rap_mode(:,:)
  LOGICAL, ALLOCATABLE :: high_sym(:)
  LOGICAL :: q_in_band_form
  ! .... variables for band plotting based on similarity of eigenvalues
  COMPLEX(DP), ALLOCATABLE :: f_of_q(:,:,:,:)
  INTEGER :: location(1), isig
  CHARACTER(LEN=6) :: int_to_char
  INTEGER            :: npk_label, nch
  CHARACTER(LEN=3), ALLOCATABLE :: letter(:)
  INTEGER, ALLOCATABLE :: label_list(:)
  LOGICAL :: tend, terr, dos
  !
!  NAMELIST /input/ flfrc, amass,  at,   &
!       &           ntyp,  &
!       &           na_ifc
  !
  IF ( my_image_id /= root_image ) RETURN

  dos = do_dos .AND. ldos
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x,"Interpolating the dynamical matrices")')
  IF (dos) THEN
     WRITE(stdout,'(5x,"Dos written on file ",a)') TRIM(fldos)
  ELSE
     WRITE(stdout,'(5x,"Frequencies written on file ",a)') TRIM(flfrq)
  ENDIF
  WRITE(stdout,'(2x,76("+"),/)')
  !
  ! ... all calculations are done by the first cpu
  !
  ! set namelist default
  !
  asr  = zasr
  flvec='matdyn.modes'
  amass_blk(:) =0.d0
  at(:,:) = 0.d0
  ntyp = 0
  na_ifc=.FALSE.
  fd=.FALSE.
  !
  ! read force constants
  !
  ntyp_blk = ntypx ! avoids fake out-of-bound error
  xmlifc=has_xml(TRIM(flfrc))
  IF (xmlifc) THEN
     CALL read_dyn_mat_param(flfrc,ntyp_blk,nat_blk)
     ALLOCATE (m_loc(3,nat_blk))
     ALLOCATE (tau_blk(3,nat_blk))
     ALLOCATE (ityp_blk(nat_blk))
     ALLOCATE (atm(ntyp_blk))
     ALLOCATE (zeu(3,3,nat_blk))
     CALL read_dyn_mat_header(ntyp_blk, nat_blk, ibrav, nspin_mag, &
               celldm, at_blk, bg_blk, omega_blk, atm, amass_blk, &
               tau_blk, ityp_blk,  m_loc, nqs, has_zstar, epsil, zeu )
     alat=celldm(1)
     call volume(alat,at_blk(1,1),at_blk(1,2),at_blk(1,3),omega_blk)
     CALL read_ifc_param(nr1,nr2,nr3)
     ALLOCATE(frc(nr1,nr2,nr3,3,3,nat_blk,nat_blk))
     CALL read_ifc(nr1,nr2,nr3,nat_blk,frc)
  ELSE
     CALL readfc ( flfrc, nr1, nr2, nr3, epsil, nat_blk, &
          ibrav, alat, at_blk, ntyp_blk, &
          amass_blk, omega_blk, has_zstar)
  ENDIF
  !
  CALL recips ( at_blk(1,1),at_blk(1,2),at_blk(1,3),  &
       bg_blk(1,1),bg_blk(1,2),bg_blk(1,3) )
  !
  ! set up (super)cell
  !
  if (ntyp < 0) then
     call errore ('matdyn','wrong ntyp ', abs(ntyp))
  else if (ntyp == 0) then
     ntyp=ntyp_blk
  end if
  !
  ! masses (for mass approximation)
  !
  DO it=1,ntyp
     IF (amass(it) < 0.d0) THEN
        CALL errore ('matdyn','wrong mass in the namelist',it)
     ELSE IF (amass(it) == 0.d0) THEN
        IF (it.LE.ntyp_blk) THEN
           WRITE (stdout,'(a,i3,a,a)') ' mass for atomic type ',it,      &
                &                     ' not given; uses mass from file ',flfrc
           amass(it) = amass_blk(it)
        ELSE
           CALL errore ('matdyn','missing mass in the namelist',it)
        END IF
     END IF
  END DO
  !
  ! lattice vectors
  !
  IF (SUM(ABS(at(:,:))) == 0.d0) THEN
     at(:,1) = at_blk(:,1)
     at(:,2) = at_blk(:,2)
     at(:,3) = at_blk(:,3)
  END IF
  !
  CALL check_at(at,bg_blk,alat,omega)
  !
  ! the supercell contains "nsc" times the original unit cell
  !
  nsc = NINT(omega/omega_blk)
  IF (ABS(omega/omega_blk-nsc) > eps) &
        CALL errore ('matdyn', 'volume ratio not integer', 1)
  !
  ! read/generate atomic positions of the (super)cell
  !
  nat = nat_blk * nsc
  nax = nat
  nax_blk = nat_blk
  !
  ALLOCATE ( tau (3, nat), ityp(nat), itau_blk(nat) )
  !
  CALL set_tau  &
       (nat, nat_blk, at, at_blk, tau, tau_blk, ityp, ityp_blk, itau_blk)
  !
  !
  ! reciprocal lattice vectors
  !
  CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  ! build the WS cell corresponding to the force constant grid
  !
  atws(:,1) = at_blk(:,1)*DBLE(nr1)
  atws(:,2) = at_blk(:,2)*DBLE(nr2)
  atws(:,3) = at_blk(:,3)*DBLE(nr3)
  ! initialize WS r-vectors
  CALL wsinit(rws,nrwsx,nrws,atws)
  !
  ! end of (super)cell setup
  !
  IF (dos) THEN
     ntetra = 6 * nq1_d * nq2_d * nq3_d
     nqx = nq1_d * nq2_d * nq3_d
     ALLOCATE ( tetra(4,ntetra), q(3,nqx) )
     CALL gen_qpoints (ibrav, at, bg, nat, tau, ityp, nq1_d, nq2_d, nq3_d, &
             ntetra, nqx, nq, q, tetra)
  ELSE
     !
     ! read q-point list
     !
     nqtot=disp_nqs
     ALLOCATE( q(3,nqtot) )
     ALLOCATE( tetra(1,1) )
     q(:,:)=disp_q(:,:)
     nq=nqtot
     ! 
  END IF
  !
  IF (asr /= 'no') THEN
     CALL set_asr (asr, nr1, nr2, nr3, frc, zeu, nat_blk, ibrav, tau_blk)
  END IF
  !
  IF (flvec.EQ.' '.OR. dos) THEN
     iout=0
  ELSE
     iout=4
     IF (ionode) OPEN (unit=iout,file=flvec,status='unknown',form='formatted')
  END IF


  ALLOCATE ( dyn(3,3,nat,nat), dyn_blk(3,3,nat_blk,nat_blk) )
  ALLOCATE ( z(3*nat,3*nat), w2(3*nat,nq), f_of_q(3,3,nat,nat) )

  IF (xmlifc) CALL set_sym(nat, tau, ityp, nspin_mag, m_loc, 6, 6, 6 )

  ALLOCATE(num_rap_mode(3*nat,nq))
  ALLOCATE(high_sym(nq))
  num_rap_mode=-1
  high_sym=.TRUE.

  DO n=1, nq
     dyn(:,:,:,:) = (0.d0, 0.d0)

     lo_to_split=.FALSE.
     f_of_q(:,:,:,:)=CMPLX(0.d0,0.d0)

     IF (na_ifc) THEN

        qq=sqrt(q(1,n)**2+q(2,n)**2+q(3,n)**3)
        if (qq == 0.0_DP) qq=1.0_DP
        qhat(1)=q(1,n)/qq
        qhat(2)=q(2,n)/qq
        qhat(3)=q(3,n)/qq

        CALL nonanal_ifc (nat,nat_blk,itau_blk,epsil,qhat,zeu,omega,dyn, &
                           nr1, nr2, nr3, f_of_q)
     END IF

     CALL setupmat (q(1,n), dyn, nat, at, bg, tau, itau_blk, nsc, alat, &
          dyn_blk, nat_blk, at_blk, bg_blk, tau_blk, omega_blk,  &
          epsil, zeu, frc, nr1,nr2,nr3, has_zstar, rws, nrws, na_ifc,f_of_q,fd)

     qhat(1) = q(1,n)*at(1,1)+q(2,n)*at(2,1)+q(3,n)*at(3,1)
     qhat(2) = q(1,n)*at(1,2)+q(2,n)*at(2,2)+q(3,n)*at(3,2)
     qhat(3) = q(1,n)*at(1,3)+q(2,n)*at(2,3)+q(3,n)*at(3,3)
     IF ( ABS( qhat(1) - NINT (qhat(1) ) ) <= eps .AND. &
          ABS( qhat(2) - NINT (qhat(2) ) ) <= eps .AND. &
          ABS( qhat(3) - NINT (qhat(3) ) ) <= eps ) THEN
        !
        ! q = 0 : we need the direction q => 0 for the non-analytic part
        !
        IF ( n == 1 ) THEN
           ! if q is the first point in the list
           IF ( nq > 1 ) THEN
              ! one more point
              qhat(:) = q(:,n) - q(:,n+1)
           ELSE
              ! no more points
              qhat(:) = 0.d0
           END IF
        ELSE IF ( n > 1 ) THEN
           ! if q is not the first point in the list
           IF ( q(1,n-1)==0.d0 .AND. &
                q(2,n-1)==0.d0 .AND. &
                q(3,n-1)==0.d0 .AND. n < nq ) THEN
               ! if the preceding q is also 0 :
              qhat(:) = q(:,n) - q(:,n+1)
           ELSE
              ! if the preceding q is npt 0 :
              qhat(:) = q(:,n) - q(:,n-1)
           END IF
        END IF
        qh = SQRT(qhat(1)**2+qhat(2)**2+qhat(3)**2)
        IF (qh /= 0.d0) qhat(:) = qhat(:) / qh
        IF (qh /= 0.d0 .AND. .NOT. has_zstar) THEN
           CALL infomsg  &
             ('matdyn','Z* not found in file '//TRIM(flfrc)// &
                          ', TO-LO splitting at q=0 will be absent!')
        ELSE
           lo_to_split=.TRUE.
        ENDIF
        !
        CALL nonanal (nat, nat_blk, itau_blk, epsil, qhat, zeu, omega, dyn)
        !
     END IF

     CALL dyndiag(nat,ntyp,amass,ityp,dyn,w2(1,n),z)
     !
     ! Cannot use the small group of \Gamma to analize the symmetry
     ! of the mode if there is an electric field.
     !
     IF (xmlifc.AND..NOT.lo_to_split) THEN
        ALLOCATE(name_rap_mode(3*nat))
        WRITE(stdout,'(10x,"xq=",3F8.4)') q(:,n)
        CALL find_representations_mode_q(nat,ntyp,q(:,n), &
                    w2(:,n),z,tau,ityp,amass,name_rap_mode, &
                    num_rap_mode(:,n), nspin_mag)
        IF (code_group==code_group_old.OR.high_sym(n-1)) high_sym(n)=.FALSE.
        code_group_old=code_group
        DEALLOCATE(name_rap_mode)
     ENDIF

     IF (ionode.AND.iout.NE.0) CALL writemodes(nat,q(1,n),w2(1,n),z,iout)
  END DO  !nq
  !
  IF(iout .NE. 0.and.ionode) CLOSE(unit=iout)
  !
  ALLOCATE (freq(3*nat, nq))
  DO n=1,nq
     ! freq(i,n) = frequencies in cm^(-1), with negative sign if omega^2 is negative
     DO i=1,3*nat
        freq(i,n)= SQRT(ABS(w2(i,n))) * RY_TO_CMM1
        IF (w2(i,n) < 0.0d0) freq(i,n) = -freq(i,n)
     END DO
  END DO
  !
  IF (flfrq.NE.' '.AND.ionode.AND..NOT.dos) THEN
     OPEN (unit=2,file=flfrq ,status='unknown',form='formatted')
     WRITE(2, '(" &plot nbnd=",i4,", nks=",i4," /")') 3*nat, nq
     DO n=1, nq
        WRITE(2, '(10x,3f10.6)')  q(1,n), q(2,n), q(3,n)
        WRITE(2,'(6f10.4)') (freq(i,n), i=1,3*nat)
     END DO
     CLOSE(unit=2)
  END IF
  !
  !  If the force constants are in the xml format we write also
  !  the file with the representations of each mode
  !
  IF (flfrq.NE.' '.AND.xmlifc.AND.ionode.AND..NOT.dos) THEN
     filename=TRIM(flfrq)//'.rap'
     OPEN (unit=2,file=filename ,status='unknown',form='formatted')
     WRITE(2, '(" &plot_rap nbnd_rap=",i4,", nks_rap=",i4," /")') 3*nat, nq
     DO n=1, nq
        WRITE(2,'(10x,3f10.6,l6)')  q(1,n), q(2,n), q(3,n), high_sym(n)
        WRITE(2,'(6i10)') (num_rap_mode(i,n), i=1,3*nat)
     END DO
     CLOSE(unit=2)
  END IF
  !
  !  write the dos file
  !
  IF (dos) THEN
     emin = 0.0d0
     emax = 0.0d0
     DO n=1,nq
        DO i=1, 3*nat
           emin = MIN (emin, freq(i,n))
           emax = MAX (emax, freq(i,n))
        END DO
     END DO
     IF (freqmin_input > 0.0_DP) THEN
        emin=freqmin_input
        freqmin=emin
     ELSE
        freqmin=emin
     ENDIF
     IF (freqmax_input > 0.0_DP) THEN
        emax=freqmax_input
        freqmax=emax
     ELSE
        freqmax=emax
     ENDIF
     !
     IF (ndos > 1) THEN
        deltafreq = (emax - emin)/(ndos-1)
     ELSE
        ndos = NINT ( (emax - emin) / deltafreq + 1.51d0 )
     END IF

     CALL set_phdos(phdos_save(igeom),ndos,deltafreq)

     IF (ionode) OPEN (unit=2,file=fldos,status='unknown',form='formatted')
     DO n= 1, ndos
        e = emin + (n - 1) * deltafreq
        CALL dos_t(freq, 1, 3*nat, nq, ntetra, tetra, e, dosofe)
        !
        ! The factor 0.5 corrects for the factor 2 in dos_t,
        ! that accounts for the spin in the electron DOS.
        !
        !WRITE (2, '(F15.10,F15.2,F15.6,F20.5)') &
        !     E, E*RY_TO_CMM1, E*RY_TO_THZ, 0.5d0*DOSofE(1)
        IF (ionode) WRITE (2, '(ES12.4,ES12.4)') e, 0.5d0*dosofe(1)
        phdos_save(igeom)%nu(n) = e
        phdos_save(igeom)%phdos(n) = 0.5d0*dosofe(1)
     END DO
     IF (ionode) CLOSE(unit=2)
  END IF  !dos
  !
  DEALLOCATE (z) 
  DEALLOCATE (w2) 
  DEALLOCATE (dyn) 
  DEALLOCATE (dyn_blk)
  DEALLOCATE (freq)
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (high_sym)
  DEALLOCATE (tau)
  DEALLOCATE (ityp)
  DEALLOCATE (itau_blk)
  DEALLOCATE(q)
  DEALLOCATE(tetra)
  DEALLOCATE(tau_blk)
  DEALLOCATE(ityp_blk)
  DEALLOCATE(zeu)
  DEALLOCATE(frc)
  IF (xmlifc) THEN
     DEALLOCATE(m_loc)
     DEALLOCATE(atm)
  ENDIF
  !
  RETURN
END SUBROUTINE matdyn_sub
!
!-----------------------------------------------------------------------
SUBROUTINE readfc ( flfrc, nr1, nr2, nr3, epsil, nat,    &
                    ibrav, alat, at, ntyp, amass, omega, has_zstar )
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE ifc,        ONLY : tau => tau_blk, ityp => ityp_blk, frc, zeu
  USE cell_base,  ONLY : celldm
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE mp,         ONLY : mp_bcast 
  USE mp_images,  ONLY : intra_image_comm 
  USE constants,  ONLY : amu_ry
  !
  IMPLICIT NONE
  ! I/O variable
  CHARACTER(LEN=256) :: flfrc
  INTEGER :: ibrav, nr1,nr2,nr3,nat, ntyp
  REAL(DP) :: alat, at(3,3), epsil(3,3)
  LOGICAL :: has_zstar
  ! local variables
  INTEGER :: i, j, na, nb, m1,m2,m3
  INTEGER :: ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  REAL(DP) :: amass(ntyp), amass_from_file, omega
  INTEGER :: nt
  CHARACTER(LEN=3) :: atm
  !
  !
  IF (ionode) OPEN (unit=1,file=TRIM(flfrc),status='old',form='formatted')
  !
  !  read cell data
  !
  IF (ionode)THEN
     READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
     if (ibrav==0) then
        read(1,*) ((at(i,j),i=1,3),j=1,3)
     end if
  ENDIF
  CALL mp_bcast(ntyp, ionode_id, intra_image_comm)
  CALL mp_bcast(nat, ionode_id, intra_image_comm)
  CALL mp_bcast(ibrav, ionode_id, intra_image_comm)
  CALL mp_bcast(celldm, ionode_id, intra_image_comm)
  IF (ibrav==0) CALL mp_bcast(at, ionode_id, intra_image_comm)
  !
  CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  alat = celldm(1)
  at = at / alat !  bring at in units of alat
  CALL volume(alat,at(1,1),at(1,2),at(1,3),omega)
  !
  !  read atomic types, positions and masses
  !
  DO nt = 1,ntyp
     IF (ionode) READ(1,*) i,atm,amass_from_file
     CALL mp_bcast(i,ionode_id, intra_image_comm)
     CALL mp_bcast(atm,ionode_id, intra_image_comm)
     CALL mp_bcast(amass_from_file,ionode_id, intra_image_comm)
     IF (i.NE.nt) CALL errore ('readfc','wrong data read',nt)
     IF (amass(nt).EQ.0.d0) THEN
        amass(nt) = amass_from_file/amu_ry
     ELSE
        WRITE(stdout,*) 'for atomic type',nt,' mass from file not used'
     END IF
  END DO
  !
  ALLOCATE (tau(3,nat), ityp(nat), zeu(3,3,nat))
  !
  DO na=1,nat
     IF (ionode) READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
     CALL mp_bcast(i,ionode_id, intra_image_comm)
     IF (i.NE.na) CALL errore ('readfc','wrong data read',na)
  END DO
  CALL mp_bcast(ityp,ionode_id, intra_image_comm)
  CALL mp_bcast(tau,ionode_id, intra_image_comm)
  !
  !  read macroscopic variable
  !
  IF (ionode) READ (1,*) has_zstar
  CALL mp_bcast(has_zstar,ionode_id,intra_image_comm)
  IF (has_zstar) THEN
     IF (ionode) READ(1,*) ((epsil(i,j),j=1,3),i=1,3)
     CALL mp_bcast(epsil,ionode_id,intra_image_comm)
     IF (ionode) THEN
        DO na=1,nat
           READ(1,*)
           READ(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
        END DO
     ENDIF
     CALL mp_bcast(zeu,ionode_id, intra_image_comm)
  ELSE
     zeu  (:,:,:) = 0.d0
     epsil(:,:) = 0.d0
  END IF
  !
  IF (ionode) READ (1,*) nr1,nr2,nr3
  CALL mp_bcast(nr1,ionode_id, intra_image_comm)
  CALL mp_bcast(nr2,ionode_id, intra_image_comm)
  CALL mp_bcast(nr3,ionode_id, intra_image_comm)
  !
  !  read real-space interatomic force constants
  !
  ALLOCATE ( frc(nr1,nr2,nr3,3,3,nat,nat) )
  frc(:,:,:,:,:,:,:) = 0.d0
  DO i=1,3
     DO j=1,3
        DO na=1,nat
           DO nb=1,nat
              IF (ionode) READ (1,*) ibid, jbid, nabid, nbbid
              CALL mp_bcast(ibid,ionode_id, intra_image_comm)
              CALL mp_bcast(jbid,ionode_id, intra_image_comm)
              CALL mp_bcast(nabid,ionode_id, intra_image_comm)
              CALL mp_bcast(nbbid,ionode_id, intra_image_comm)
              IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
                 na.NE.nabid .OR. nb.NE.nbbid)                      &
                 CALL errore  ('readfc','error in reading',1)
              IF (ionode) READ (1,*) (((m1bid, m2bid, m3bid,        &
                          frc(m1,m2,m3,i,j,na,nb),                  &
                           m1=1,nr1),m2=1,nr2),m3=1,nr3)
               
              CALL mp_bcast(frc(:,:,:,i,j,na,nb),ionode_id, intra_image_comm)
           END DO
        END DO
     END DO
  END DO
  !
  IF (ionode) CLOSE(unit=1)
  !
  RETURN
END SUBROUTINE readfc
!
!-----------------------------------------------------------------------
SUBROUTINE frc_blk(dyn,q,tau,nat,nr1,nr2,nr3,frc,at,bg,rws,nrws,f_of_q,fd)
  !-----------------------------------------------------------------------
  ! calculates the dynamical matrix at q from the (short-range part of the)
  ! force constants
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  INTEGER nr1, nr2, nr3, nat, n1, n2, n3, &
          ipol, jpol, na, nb, m1, m2, m3, nint, i,j, nrws, nax
  COMPLEX(DP) dyn(3,3,nat,nat), f_of_q(3,3,nat,nat)
  REAL(DP) frc(nr1,nr2,nr3,3,3,nat,nat), tau(3,nat), q(3), arg, &
               at(3,3), bg(3,3), r(3), weight, r_ws(3),  &
               total_weight, rws(0:3,nrws), alat
  REAL(DP), EXTERNAL :: wsweight
  REAL(DP),SAVE,ALLOCATABLE :: wscache(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: ttt(:,:,:,:,:), tttx(:,:)
  LOGICAL,SAVE :: first=.true.
  LOGICAL :: fd
  !
  FIRST_TIME : IF (first) THEN
    first=.false.
    ALLOCATE( wscache(-2*nr3:2*nr3, -2*nr2:2*nr2, -2*nr1:2*nr1, nat,nat) )
    DO na=1, nat
       DO nb=1, nat
          total_weight=0.0d0
          !
          DO n1=-2*nr1,2*nr1
             DO n2=-2*nr2,2*nr2
                DO n3=-2*nr3,2*nr3
                   DO i=1, 3
                      r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                      r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                      if (fd) r_ws(i) = r(i) + tau(i,nb)-tau(i,na)
                   END DO
                   wscache(n3,n2,n1,nb,na) = wsweight(r_ws,rws,nrws)
                ENDDO
             ENDDO
          ENDDO
      ENDDO
    ENDDO
  ENDIF FIRST_TIME
  !
  ALLOCATE(ttt(3,nat,nr1,nr2,nr3))
  ALLOCATE(tttx(3,nat*nr1*nr2*nr3))
  ttt(:,:,:,:,:)=0.d0

  DO na=1, nat
     DO nb=1, nat
        total_weight=0.0d0
        DO n1=-2*nr1,2*nr1
           DO n2=-2*nr2,2*nr2
              DO n3=-2*nr3,2*nr3
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 DO i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                 END DO

                 weight = wscache(n3,n2,n1,nb,na) 
                 IF (weight .GT. 0.0d0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1+1,nr1)
                    IF(m1.LE.0) m1=m1+nr1
                    m2 = MOD(n2+1,nr2)
                    IF(m2.LE.0) m2=m2+nr2
                    m3 = MOD(n3+1,nr3)
                    IF(m3.LE.0) m3=m3+nr3
                 !   write(*,'(6i4)') n1,n2,n3,m1,m2,m3
                    !
                    ! FOURIER TRANSFORM
                    !
                    do i=1,3
                       ttt(i,na,m1,m2,m3)=tau(i,na)+m1*at(i,1)+m2*at(i,2)+m3*at(i,3)
                       ttt(i,nb,m1,m2,m3)=tau(i,nb)+m1*at(i,1)+m2*at(i,2)+m3*at(i,3)
                    end do

                    arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                    DO ipol=1, 3
                       DO jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               (frc(m1,m2,m3,ipol,jpol,na,nb)+f_of_q(ipol,jpol,na,nb))     &
                               *CMPLX(COS(arg),-SIN(arg),kind=DP)*weight
                       END DO
                    END DO
                 END IF
                 total_weight=total_weight + weight
              END DO
           END DO
        END DO
        IF (ABS(total_weight-nr1*nr2*nr3).GT.1.0d-8) THEN
           WRITE(stdout,*) total_weight
           CALL errore ('frc_blk','wrong total_weight',1)
        END IF
     END DO
  END DO
  !
!  alat=10.2
!  nax=0
!  DO n1=1,nr1
!     DO n2=1,nr2
!        DO n3=1,nr3
!           do na=1,nat
!              nax=nax+1
!              do i=1,3
!                 tttx(i,nax)=ttt(i,na,n1,n2,n3)*alat*0.529177
!              end do
!           end do
!        end do
!     end do
!  end do
!
!  do nb=1,nat
!     write(6,'(3(f15.9,1x))') tau(1,nb),tau(2,nb),tau(3,nb)
!  enddo
!  print*, '========='
!  do nb=1,nat*nr1*nr2*nr3
!     write(6,'(3(f15.9,1x))') tttx(1,nb),tttx(2,nb),tttx(3,nb)
!  enddo
! 
  RETURN
END SUBROUTINE frc_blk
!
!-----------------------------------------------------------------------
SUBROUTINE setupmat (q,dyn,nat,at,bg,tau,itau_blk,nsc,alat, &
     &         dyn_blk,nat_blk,at_blk,bg_blk,tau_blk,omega_blk, &
     &                 epsil,zeu,frc,nr1,nr2,nr3,has_zstar,rws,nrws,na_ifc,f_of_q,fd)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER:: nr1, nr2, nr3, nat, nat_blk, nsc, nrws, itau_blk(nat)
  REAL(DP) :: q(3), tau(3,nat), at(3,3), bg(3,3), alat,      &
                  epsil(3,3), zeu(3,3,nat_blk), rws(0:3,nrws),   &
                  frc(nr1,nr2,nr3,3,3,nat_blk,nat_blk)
  REAL(DP) :: tau_blk(3,nat_blk), at_blk(3,3), bg_blk(3,3), omega_blk
  COMPLEX(DP) dyn_blk(3,3,nat_blk,nat_blk), f_of_q(3,3,nat,nat)
  COMPLEX(DP) ::  dyn(3,3,nat,nat)
  LOGICAL :: has_zstar, na_ifc, fd
  !
  ! local variables
  !
  REAL(DP) :: arg
  COMPLEX(DP) :: cfac(nat)
  INTEGER :: i,j,k, na,nb, na_blk, nb_blk, iq
  REAL(DP) :: qp(3), qbid(3,nsc) ! automatic array
  !
  !
  CALL q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !
  DO iq=1,nsc
     !
     DO k=1,3
        qp(k)= q(k) + qbid(k,iq)
     END DO
     !
     dyn_blk(:,:,:,:) = (0.d0,0.d0)
     CALL frc_blk (dyn_blk,qp,tau_blk,nat_blk,              &
          &              nr1,nr2,nr3,frc,at_blk,bg_blk,rws,nrws,f_of_q,fd)
      IF (has_zstar .and. .not.na_ifc) &
           CALL rgd_blk(nr1,nr2,nr3,nat_blk,dyn_blk,qp,tau_blk,   &
                        epsil,zeu,bg_blk,omega_blk,+1.d0)
     !
     DO na=1,nat
        na_blk = itau_blk(na)
        DO nb=1,nat
           nb_blk = itau_blk(nb)
           !
           arg=tpi* ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
           !
           cfac(nb) = CMPLX(COS(arg),SIN(arg),kind=DP)/nsc
           !
        END DO ! nb
        !
        DO i=1,3
           DO j=1,3
              !
              DO nb=1,nat
                 nb_blk = itau_blk(nb)
                 dyn(i,j,na,nb) = dyn(i,j,na,nb) + cfac(nb) * &
                      dyn_blk(i,j,na_blk,nb_blk)
              END DO ! nb
              !
           END DO ! j
        END DO ! i
     END DO ! na
     !
  END DO ! iq
  !
  RETURN
END SUBROUTINE setupmat
!
!
!----------------------------------------------------------------------
SUBROUTINE set_asr (asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  CHARACTER (LEN=10), intent(in) :: asr
  INTEGER, intent(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(DP), intent(in) :: tau(3,nat)
  REAL(DP), intent(inout) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1
  REAL(DP) :: zeu_new(3,3,nat)
  REAL(DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
  type vector
     real(DP),pointer :: vec(:,:,:,:,:,:,:)
  end type vector
  !
  type (vector) u(6*3*nat)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  integer :: u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:,:,:)
  real(DP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
  ! should be only 2 of them) and the value of that element. We do so in order
  ! to limit the amount of memory used.
  !
  real(DP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
  ! temporary vectors and parameters
  real(DP) :: scal,norm2, sum
  !
  real(DP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (if asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim')) then
     call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  endif
  !
  if(asr.eq.'simple') then
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
              sum = sum + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           end do
        end do
     end do
     !
     ! Simple Acoustic Sum Rule on force constants in real space
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              sum=0.0d0
               do nb=1,nat
                  do n1=1,nr1
                     do n2=1,nr2
                        do n3=1,nr3
                           sum=sum+frc(n1,n2,n3,i,j,na,nb)
                        end do
                     end do
                  end do
               end do
               frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
               !               write(6,*) ' na, i, j, sum = ',na,i,j,sum
            end do
         end do
      end do
      !
      return
      !
   end if

  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('set_asr','too many directions of &
             & periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'asr: rotational axis may be wrong'
     endif
     write(stdout,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(asr.eq.'zero-dim') n=6
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:,:,:,:)=0.0d0
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu_new(i,j,na)=zeu(i,j,na)
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p=p+1
        zeu_u(p,i,j,:)=1.0d0
        !
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p=p+1
        do na=1,nat
           zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
        enddo
        !
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p=p+1
           do na=1,nat
              zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
              zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
           enddo
           !
        enddo
     enddo
  endif
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less=0
  do k=1,p
     zeu_w(:,:,:)=zeu_u(k,:,:,:)
     zeu_x(:,:,:)=zeu_u(k,:,:,:)
     do q=1,k-1
        r=1
        do izeu_less=1,nzeu_less
           if (zeu_less(izeu_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
        endif
     enddo
     call sp_zeu(zeu_w,zeu_w,nat,norm2)
     if (norm2.gt.1.0d-16) then
        zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
     else
        nzeu_less=nzeu_less+1
        zeu_less(nzeu_less)=k
     endif
  enddo
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:,:,:)=0.0d0
  do k=1,p
     r=1
     do izeu_less=1,nzeu_less
        if (zeu_less(izeu_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        call sp_zeu(zeu_x,zeu_new,nat,scal)
        zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
     endif
  enddo
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
  call sp_zeu(zeu_w,zeu_w,nat,norm2)
  write(stdout,'("Norm of the difference between old and new effective ", &
       & "charges: ",F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection of zeu")')
  !do k=1,p
  !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
  !  call sp_zeu(zeu_x,zeu_new,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu(i,j,na)=zeu_new(i,j,na)
        enddo
     enddo
  enddo
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  do k=1,18*nat
     allocate(u(k) % vec(nr1,nr2,nr3,3,3,nat,nat))
     u(k) % vec (:,:,:,:,:,:,:)=0.0d0
  enddo
  ALLOCATE (frc_new(nr1,nr2,nr3,3,3,nat,nat))
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc_new(n1,n2,n3,i,j,na,nb)=frc(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        do na=1,nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           u(p) % vec (:,:,:,i,j,na,:)=1.0d0
           !
        enddo
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        do na=1,nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           do nb=1,nat
              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
           enddo
           !
        enddo
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              do nb=1,nat
                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
              enddo
              !
           enddo
        enddo
     enddo
  endif
  !
  allocate (ind_v(9*nat*nat*nr1*nr2*nr3,2,7), v(9*nat*nat*nr1*nr2*nr3,2) )
  m=0
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       ! These are the vectors associated with the symmetry constraints
                       q=1
                       l=1
                       do while((l.le.m).and.(q.ne.0))
                          if ((ind_v(l,1,1).eq.n1).and.(ind_v(l,1,2).eq.n2).and. &
                               (ind_v(l,1,3).eq.n3).and.(ind_v(l,1,4).eq.i).and. &
                               (ind_v(l,1,5).eq.j).and.(ind_v(l,1,6).eq.na).and. &
                               (ind_v(l,1,7).eq.nb)) q=0
                          if ((ind_v(l,2,1).eq.n1).and.(ind_v(l,2,2).eq.n2).and. &
                               (ind_v(l,2,3).eq.n3).and.(ind_v(l,2,4).eq.i).and. &
                               (ind_v(l,2,5).eq.j).and.(ind_v(l,2,6).eq.na).and. &
                               (ind_v(l,2,7).eq.nb)) q=0
                          l=l+1
                       enddo
                       if ((n1.eq.MOD(nr1+1-n1,nr1)+1).and.(n2.eq.MOD(nr2+1-n2,nr2)+1) &
                            .and.(n3.eq.MOD(nr3+1-n3,nr3)+1).and.(i.eq.j).and.(na.eq.nb)) q=0
                       if (q.ne.0) then
                          m=m+1
                          ind_v(m,1,1)=n1
                          ind_v(m,1,2)=n2
                          ind_v(m,1,3)=n3
                          ind_v(m,1,4)=i
                          ind_v(m,1,5)=j
                          ind_v(m,1,6)=na
                          ind_v(m,1,7)=nb
                          v(m,1)=1.0d0/DSQRT(2.0d0)
                          ind_v(m,2,1)=MOD(nr1+1-n1,nr1)+1
                          ind_v(m,2,2)=MOD(nr2+1-n2,nr2)+1
                          ind_v(m,2,3)=MOD(nr3+1-n3,nr3)+1
                          ind_v(m,2,4)=j
                          ind_v(m,2,5)=i
                          ind_v(m,2,6)=nb
                          ind_v(m,2,7)=na
                          v(m,2)=-1.0d0/DSQRT(2.0d0)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less=0
  allocate (w(nr1,nr2,nr3,3,3,nat,nat), x(nr1,nr2,nr3,3,3,nat,nat))
  do k=1,p
     w(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     do l=1,m
        !
        call sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
        do r=1,2
           n1=ind_v(l,r,1)
           n2=ind_v(l,r,2)
           n3=ind_v(l,r,3)
           i=ind_v(l,r,4)
           j=ind_v(l,r,5)
           na=ind_v(l,r,6)
           nb=ind_v(l,r,7)
           w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)-scal*v(l,r)
        enddo
     enddo
     if (k.le.(9*nat)) then
        na1=MOD(k,nat)
        if (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat,3)+1
        i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
     else
        q=k-9*nat
        if (n.eq.4) then
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           i1=MOD((q-na1)/nat,3)+1
        else
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           j1=MOD((q-na1)/nat,3)+1
           i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
        endif
     endif
     do q=1,k-1
        r=1
        do i_less=1,n_less
           if (u_less(i_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp3(x,u(q) % vec (:,:,:,:,:,:,:), i1,na1,nr1,nr2,nr3,nat,scal)
           w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) - scal* u(q) % vec (:,:,:,:,:,:,:)
        endif
     enddo
     call sp1(w,w,nr1,nr2,nr3,nat,norm2)
     if (norm2.gt.1.0d-16) then
        u(k) % vec (:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / DSQRT(norm2)
     else
        n_less=n_less+1
        u_less(n_less)=k
     endif
  enddo
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:,:,:,:,:,:,:)=0.0d0
  do l=1,m
     call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
     do r=1,2
        n1=ind_v(l,r,1)
        n2=ind_v(l,r,2)
        n3=ind_v(l,r,3)
        i=ind_v(l,r,4)
        j=ind_v(l,r,5)
        na=ind_v(l,r,6)
        nb=ind_v(l,r,7)
        w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)+scal*v(l,r)
     enddo
  enddo
  do k=1,p
     r=1
     do i_less=1,n_less
        if (u_less(i_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
        call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
        w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) + scal*u(k)%vec(:,:,:,:,:,:,:)
     endif
     deallocate(u(k) % vec)
  enddo
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - w(:,:,:,:,:,:,:)
  call sp1(w,w,nr1,nr2,nr3,nat,norm2)
  write(stdout,'("Norm of the difference between old and new force-constants:",&
       &     F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection IFC")')
  !do l=1,m
  !  call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("l= ",I8," frc_new|v(l)= ",F15.10)') l,scal
  !enddo
  !do k=1,p
  !  x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
  !  call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," frc_new|u(k)= ",F15.10)') k,scal
  !  deallocate(u(k) % vec)
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc(n1,n2,n3,i,j,na,nb)=frc_new(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate (x, w)
  deallocate (v, ind_v)
  deallocate (frc_new)
  !
  return
end subroutine set_asr
!
!----------------------------------------------------------------------
subroutine sp1(u,v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        do nb=1,nat
          do n1=1,nr1
            do n2=1,nr2
              do n3=1,nr3
                scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v,ind_v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
  ! but v is coded as explained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  integer ind_v(2,7)
  real(DP) v(2)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),ind_v(i,5),ind_v(i,6), &
         ind_v(i,7))*v(i)
  enddo
  !
  return
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v,i,na,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do j=1,3
    do nb=1,nat
      do n1=1,nr1
        do n2=1,nr2
          do n3=1,nr3
            scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp3
!
!-----------------------------------------------------------------------
SUBROUTINE q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !-----------------------------------------------------------------------
  ! generate list of q (qbid) that are G-vectors of the supercell
  ! but not of the bulk
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: nsc
  REAL(DP) qbid(3,nsc), at_blk(3,3), bg_blk(3,3), at(3,3), bg(3,3)
  !
  INTEGER, PARAMETER:: nr1=4, nr2=4, nr3=4, &
                       nrm=(2*nr1+1)*(2*nr2+1)*(2*nr3+1)
  REAL(DP), PARAMETER:: eps=1.0d-7
  INTEGER :: i, j, k,i1, i2, i3, idum(nrm), iq
  REAL(DP) :: qnorm(nrm), qbd(3,nrm) ,qwork(3), delta
  LOGICAL lbho
  !
  i = 0
  DO i1=-nr1,nr1
     DO i2=-nr2,nr2
        DO i3=-nr3,nr3
           i = i + 1
           DO j=1,3
              qwork(j) = i1*bg(j,1) + i2*bg(j,2) + i3*bg(j,3)
           END DO ! j
           !
           qnorm(i)  = qwork(1)**2 + qwork(2)**2 + qwork(3)**2
           !
           DO j=1,3
              !
              qbd(j,i) = at_blk(1,j)*qwork(1) + &
                         at_blk(2,j)*qwork(2) + &
                         at_blk(3,j)*qwork(3)
           END DO ! j
           !
           idum(i) = 1
           !
        END DO ! i3
     END DO ! i2
  END DO ! i1
  !
  DO i=1,nrm-1
     IF (idum(i).EQ.1) THEN
        DO j=i+1,nrm
           IF (idum(j).EQ.1) THEN
              lbho=.TRUE.
              DO k=1,3
                 delta = qbd(k,i)-qbd(k,j)
                 lbho = lbho.AND. (ABS(NINT(delta)-delta).LT.eps)
              END DO ! k
              IF (lbho) THEN
                 IF(qnorm(i).GT.qnorm(j)) THEN
                    qbd(1,i) = qbd(1,j)
                    qbd(2,i) = qbd(2,j)
                    qbd(3,i) = qbd(3,j)
                    qnorm(i) = qnorm(j)
                 END IF
                 idum(j) = 0
              END IF
           END IF
        END DO ! j
     END IF
  END DO ! i
  !
  iq = 0
  DO i=1,nrm
     IF (idum(i).EQ.1) THEN
        iq=iq+1
        qbid(1,iq)= bg_blk(1,1)*qbd(1,i) +  &
                    bg_blk(1,2)*qbd(2,i) +  &
                    bg_blk(1,3)*qbd(3,i)
        qbid(2,iq)= bg_blk(2,1)*qbd(1,i) +  &
                    bg_blk(2,2)*qbd(2,i) +  &
                    bg_blk(2,3)*qbd(3,i)
        qbid(3,iq)= bg_blk(3,1)*qbd(1,i) +  &
                    bg_blk(3,2)*qbd(2,i) +  &
                    bg_blk(3,3)*qbd(3,i)
     END IF
  END DO ! i
  !
  IF (iq.NE.nsc) CALL errore('q_gen',' probably nr1,nr2,nr3 too small ', iq)
  RETURN
END SUBROUTINE q_gen
!
!-----------------------------------------------------------------------
SUBROUTINE check_at(at,bg_blk,alat,omega)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP) :: at(3,3), bg_blk(3,3), alat, omega
  REAL(DP) :: work(3,3)
  INTEGER :: i,j
  REAL(DP), PARAMETER :: small=1.d-6
  !
  work(:,:) = at(:,:)
  CALL cryst_to_cart(3,work,bg_blk,-1)
  !
  DO j=1,3
     DO i =1,3
        IF ( ABS(work(i,j)-NINT(work(i,j))) > small) THEN
           WRITE (stdout,'(3f9.4)') work(:,:)
           CALL errore ('check_at','at not multiple of at_blk',1)
        END IF
     END DO
  END DO
  !
  omega =alat**3 * ABS(at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))- &
                       at(1,2)*(at(2,1)*at(3,3)-at(2,3)*at(3,1))+ &
                       at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))
  !
  RETURN
END SUBROUTINE check_at
!
!-----------------------------------------------------------------------
SUBROUTINE set_tau (nat, nat_blk, at, at_blk, tau, tau_blk, &
     ityp, ityp_blk, itau_blk)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER nat, nat_blk,ityp(nat),ityp_blk(nat_blk), itau_blk(nat)
  REAL(DP) at(3,3),at_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  REAL(DP) bg(3,3), r(3) ! work vectors
  INTEGER i,i1,i2,i3,na,na_blk
  REAL(DP) small
  INTEGER NN1,NN2,NN3
  PARAMETER (NN1=8, NN2=8, NN3=8, small=1.d-8)
  !
  CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  na = 0
  !
  DO i1 = -NN1,NN1
     DO i2 = -NN2,NN2
        DO i3 = -NN3,NN3
           r(1) = i1*at_blk(1,1) + i2*at_blk(1,2) + i3*at_blk(1,3)
           r(2) = i1*at_blk(2,1) + i2*at_blk(2,2) + i3*at_blk(2,3)
           r(3) = i1*at_blk(3,1) + i2*at_blk(3,2) + i3*at_blk(3,3)
           CALL cryst_to_cart(1,r,bg,-1)
           !
           IF ( r(1).GT.-small .AND. r(1).LT.1.d0-small .AND.          &
                r(2).GT.-small .AND. r(2).LT.1.d0-small .AND.          &
                r(3).GT.-small .AND. r(3).LT.1.d0-small ) THEN
              CALL cryst_to_cart(1,r,at,+1)
              !
              DO na_blk=1, nat_blk
                 na = na + 1
                 IF (na.GT.nat) CALL errore('set_tau','too many atoms',na)
                 tau(1,na)    = tau_blk(1,na_blk) + r(1)
                 tau(2,na)    = tau_blk(2,na_blk) + r(2)
                 tau(3,na)    = tau_blk(3,na_blk) + r(3)
                 ityp(na)     = ityp_blk(na_blk)
                 itau_blk(na) = na_blk
              END DO
              !
           END IF
           !
        END DO
     END DO
  END DO
  !
  IF (na.NE.nat) CALL errore('set_tau','too few atoms: increase NNs',na)
  !
  RETURN
END SUBROUTINE set_tau
!
!-----------------------------------------------------------------------
SUBROUTINE gen_qpoints (ibrav, at_, bg_, nat, tau, ityp, nk1, nk2, nk3, &
     ntetra, nqx, nq, q, tetra)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : set_sym_bl, find_sym, s, irt, nsym, &
                         nrot, t_rev, time_reversal,  sname
  !
  IMPLICIT NONE
  ! input
  INTEGER :: ibrav, nat, nk1, nk2, nk3, ntetra, ityp(*)
  REAL(DP) :: at_(3,3), bg_(3,3), tau(3,nat)
  ! output
  INTEGER :: nqx, nq, tetra(4,ntetra)
  REAL(DP) :: q(3,nqx)
  ! local
  REAL(DP) :: xqq(3), wk(nqx), mdum(3,nat)
  LOGICAL :: magnetic_sym=.FALSE., skip_equivalence=.FALSE.
  !
  time_reversal = .true.
  t_rev(:) = 0
  xqq (:) =0.d0
  at = at_
  bg = bg_
  CALL set_sym_bl ( )
  !
  CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, bg, nqx, &
                           0,0,0, nk1,nk2,nk3, nq, q, wk)
  !
  CALL find_sym ( nat, tau, ityp, 6, 6, 6, .not.time_reversal, mdum )
  !
  CALL irreducible_BZ (nrot, s, nsym, time_reversal, magnetic_sym, &
                       at, bg, nqx, nq, q, wk, t_rev)
  !
  IF (ntetra /= 6 * nk1 * nk2 * nk3) &
       CALL errore ('gen_qpoints','inconsistent ntetra',1)
  !
  CALL tetrahedra (nsym, s, time_reversal, t_rev, at, bg, nqx, 0, 0, 0, &
       nk1, nk2, nk3, nq, q, wk, ntetra, tetra)
  !
  RETURN
END SUBROUTINE gen_qpoints
!
!
SUBROUTINE find_representations_mode_q ( nat, ntyp, xq, w2, u, tau, ityp, &
                  amass, name_rap_mode, num_rap_mode, nspin_mag )

  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : find_sym, s, sr, ftau, irt, nsym, &
                         nrot, t_rev, time_reversal, sname, copy_sym, &
                         s_axis_to_cart

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, ntyp, nspin_mag
  REAL(DP), INTENT(IN) :: xq(3), amass(ntyp), tau(3,nat)
  REAL(DP), INTENT(IN) :: w2(3*nat)
  INTEGER, INTENT(IN) :: ityp(nat)
  COMPLEX(DP), INTENT(IN) :: u(3*nat,3*nat)
  CHARACTER(15), INTENT(OUT) :: name_rap_mode(3*nat)
  INTEGER, INTENT(OUT) :: num_rap_mode(3*nat)
  REAL(DP) :: gi (3, 48), gimq (3), sr_is(3,3,48), rtau(3,48,nat)
  INTEGER :: irotmq, nsymq, nsym_is, isym, i, ierr
  LOGICAL :: minus_q, search_sym, sym(48), magnetic_sym
!
!  find the small group of q
!
  time_reversal=.TRUE.
  IF (.NOT.time_reversal) minus_q=.FALSE.

  sym(1:nsym)=.true.
  call smallg_q (xq, 0, at, bg, nsym, s, ftau, sym, minus_q)
  nsymq=copy_sym(nsym,sym )
  call s_axis_to_cart ()
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
!
!  if the small group of q is non symmorphic,
!  search the symmetries only if there are no G such that Sq -> q+G
!
  search_sym=.TRUE.
  IF ( ANY ( ftau(:,1:nsymq) /= 0 ) ) THEN
     DO isym=1,nsymq
        search_sym=( search_sym.and.(abs(gi(1,isym))<1.d-8).and.  &
                                    (abs(gi(2,isym))<1.d-8).and.  &
                                    (abs(gi(3,isym))<1.d-8) )
     END DO
  END IF
!
!  Set the representations tables of the small group of q and
!  find the mode symmetry
!
  IF (search_sym) THEN
     magnetic_sym=(nspin_mag==4)
     CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)
     sym (1:nsym) = .TRUE.
     CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
     CALL find_mode_sym_new (u, w2, tau, nat, nsymq, sr, irt, xq,    &
             rtau, amass, ntyp, ityp, 1, .FALSE., .FALSE., num_rap_mode, ierr)

  ENDIF
  RETURN
  END SUBROUTINE find_representations_mode_q
