!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE q2r_sub(fildyn)
  !----------------------------------------------------------------------------
  ! 
  ! Reads force constant matrices C(q) produced by the phonon code
  ! for a grid of q-points, calculates the corresponding set of
  ! interatomic force constants (IFC), C(R)
  !
  ! This is a version of the program q2r in the form of a subroutine.
  ! It can be called by the phonon and uses its variables.
  ! Two input variables of this code have been added to the input of
  ! the phonon code. The others are those used by the phonon code.
  !
  !
  !  Input data that are supposed to be read from the phonon code.
  !          
  !     flfrc      :  output file containing the IFC in real space
  !                   (character, must be specified)
  !     zasr       :  Indicates type of Acoustic Sum Rules used for the Born
  !                   effective charges (character):
  !                   - 'no': no Acoustic Sum Rules imposed (default)
  !                   - 'simple':  previous implementation of the asr used
  !                     (3 translational asr imposed by correction of
  !                     the diagonal elements of the force-constants matrix)
  !                   - 'crystal': 3 translational asr imposed by optimized
  !                      correction of the IFC (projection).
  !                   - 'one-dim': 3 translational asr + 1 rotational asr
  !                     imposed by optimized correction of the IFC (the
  !                     rotation axis is the direction of periodicity; it
  !                     will work only if this axis considered is one of
  !                     the cartesian axis).
  !                   - 'zero-dim': 3 translational asr + 3 rotational asr
  !                     imposed by optimized correction of the IFC.
  !                   Note that in certain cases, not all the rotational asr
  !                   can be applied (e.g. if there are only 2 atoms in a
  !                   molecule or if all the atoms are aligned, etc.).
  !                   In these cases the supplementary asr are cancelled
  !                   during the orthonormalization procedure (see below).
  !
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_bcast
  USE dynamicalq, ONLY : phiq, tau, ityp, zeu
  USE ifc,        ONLY : zasr
  USE fft_scalar, ONLY : cfft3d
  USE io_global,  ONLY : stdout, meta_ionode, meta_ionode_id
  USE mp_images,  ONLY : my_image_id
  USE mp_world,   ONLY : world_comm
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail, &
                         write_dyn_mat_header, write_ifc
  USE data_files, ONLY : flfrc
  USE control_ph, ONLY : ldisp, xmldyn
  USE rigid,      ONLY : rgd_blk
  !
  IMPLICIT NONE
  !
  INTEGER,  PARAMETER  :: ntypx = 10
  REAL(DP), PARAMETER  :: eps=1.D-5, eps12=1.d-12
  INTEGER              :: nr1, nr2, nr3, nr(3)
  !     dimensions of the FFT grid formed by the q-point grid
  !
  CHARACTER(len=256) :: fildyn, filin, filefrc
  CHARACTER(len=3)   :: atm(ntypx)
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  LOGICAL :: lq, lrigid, lrigid1
  INTEGER :: m1, m2, m3, m(3), i, j, j1, j2, na1, na2, ipol, nn
  INTEGER :: nat, nq, ntyp, iq, icar, nfile, ifile, nqs, nq_log
  INTEGER :: na, nt
  !
  INTEGER :: ibrav, ierr, nspin_mag, iundyn, iunfrc, ios
  INTEGER :: find_free_unit
  !
  INTEGER,     ALLOCATABLE ::  nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: phid(:,:,:,:,:)
  REAL(DP),    ALLOCATABLE :: m_loc(:,:)
  !
  REAL(DP) :: celldm(6), at(3,3), bg(3,3)
  REAL(DP) :: q(3,48), omega, xq, amass(ntypx), resi
  REAL(DP) :: epsil(3,3), smat(3,3), angle_rot
  !
  ! Only one image run this routine, but the results are broadcasted 
  ! to all images
  !
  IF (flfrc == ' '.OR. .NOT. ldisp) RETURN
  !
  filefrc="phdisp_files/"//TRIM(flfrc)
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x,"Computing the interatomic force constants")')
  WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(filefrc)
  WRITE(stdout,'(2x,76("+"),/)')
  !
  IF (meta_ionode) THEN
     iundyn=find_free_unit()
     OPEN (UNIT=iundyn, FILE=TRIM(fildyn)//'0', STATUS='old', &
                                            FORM='formatted', IOSTAT=ierr)
  ENDIF
  CALL mp_bcast(ierr, meta_ionode_id, world_comm)
  IF (ierr /= 0) CALL errore('q2r_sub','No grid information on file',1)
  IF (meta_ionode) THEN
     WRITE (stdout,'(/,5x,"Reading q grid from file ")') 
     WRITE (stdout,'(5x,a)') TRIM(fildyn)//'0'
     READ (iundyn, *) nr1, nr2, nr3
     READ (iundyn, *) nfile
     CLOSE (UNIT=iundyn, STATUS='KEEP')
  ENDIF
  CALL mp_bcast(nr1, meta_ionode_id, world_comm)
  CALL mp_bcast(nr2, meta_ionode_id, world_comm)
  CALL mp_bcast(nr3, meta_ionode_id, world_comm)
  CALL mp_bcast(nfile, meta_ionode_id, world_comm)
     !
  IF (nr1 < 1 .OR. nr1 > 1024) CALL errore ('q2r_sub',' nr1 wrong or missing',1)
  IF (nr2 < 1 .OR. nr2 > 1024) CALL errore ('q2r_sub',' nr2 wrong or missing',1)
  IF (nr3 < 1 .OR. nr2 > 1024) CALL errore ('q2r_sub',' nr3 wrong or missing',1)
  IF (nfile < 1 .OR. nfile > 1024) &
     CALL errore ('q2r_sub','too few or too many file',MAX(1,nfile))
     !
     ! copy nrX -> nr(X)
     !
  nr(1) = nr1
  nr(2) = nr2
  nr(3) = nr3
  !
  ! D matrix (analytical part)
  !
  ntyp = ntypx ! avoids spurious out-of-bound errors
  !
  ALLOCATE ( nc(nr1,nr2,nr3) )
  nc = 0
  !
  ! Force constants in reciprocal space read from file
  !
  DO ifile=1,nfile
     filin = TRIM(fildyn) // TRIM( int_to_char( ifile ) )
     WRITE (stdout,'(/,5x,"Reading force constants from file")')
     WRITE (stdout,'(5x,a)') TRIM(filin)
     IF (xmldyn) THEN
        IF (my_image_id==0) CALL read_dyn_mat_param(filin,ntyp,nat)
        CALL mp_bcast(ntyp,meta_ionode_id,world_comm)
        CALL mp_bcast(nat,meta_ionode_id,world_comm)

        IF (ifile==1) THEN
           ALLOCATE (m_loc(3,nat))
           ALLOCATE (tau(3,nat))
           ALLOCATE (ityp(nat))
           ALLOCATE (zeu(3,3,nat))
        ENDIF
        IF (my_image_id==0) THEN
           IF (ifile==1) THEN
              CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
                  celldm, at, bg, omega, atm, amass, tau, ityp, &
                  m_loc, nqs, lrigid, epsil, zeu )
           ELSE
              CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
                 celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqs)
           ENDIF
        ENDIF

        IF (ifile==1) THEN
           CALL mp_bcast(nspin_mag, meta_ionode_id, world_comm)
           CALL mp_bcast(lrigid, meta_ionode_id, world_comm)
           CALL mp_bcast(zeu, meta_ionode_id, world_comm)
        ENDIF
        CALL mp_bcast(nqs, meta_ionode_id, world_comm)
        ALLOCATE (phiq(3,3,nat,nat,nqs) )
        IF (my_image_id==0) THEN
           DO iq=1,nqs
              CALL read_dyn_mat(nat,iq,q(:,iq),phiq(:,:,:,:,iq))
           ENDDO
           CALL read_dyn_mat_tail(nat)
        ENDIF
        CALL mp_bcast(q, meta_ionode_id, world_comm)
        CALL mp_bcast(phiq, meta_ionode_id, world_comm)
     ELSE
        IF (meta_ionode) THEN
           iundyn=find_free_unit()
           OPEN (UNIT=iundyn, FILE=TRIM(filin), STATUS='old', &
                                                FORM='formatted', IOSTAT=ierr)
        ENDIF
        CALL mp_bcast(ierr, meta_ionode_id, world_comm)
        IF (ierr /= 0) CALL errore('q2r_sub','file '//TRIM(filin)&
                                                           //' missing!',1)
        IF (my_image_id==0) THEN
           CALL read_dyn_from_file_tpw (nqs, q, epsil, lrigid,  &
                    ntyp, nat, ibrav, celldm, at, atm, amass, ifile, iundyn)
        ENDIF
        CALL mp_bcast(nat, meta_ionode_id, world_comm)
        CALL mp_bcast(ntyp, meta_ionode_id, world_comm)
        CALL mp_bcast(nqs, meta_ionode_id, world_comm)
        WRITE(6,*) nat
        IF (my_image_id/=0.AND.ifile==1) THEN
           ALLOCATE (tau(3,nat))
           ALLOCATE (ityp(nat))
           ALLOCATE (zeu(3,3,nat))
           ALLOCATE (phiq(3,3,nat,nat,48))
        ENDIF
        CALL mp_bcast(phiq, meta_ionode_id, world_comm)
        CALL mp_bcast(q, meta_ionode_id, world_comm)
        IF (ifile==1) ALLOCATE (m_loc(3,nat))
        IF (meta_ionode) CLOSE(unit=iundyn)
     ENDIF
     IF (ifile==1) THEN
        CALL mp_bcast(ibrav, meta_ionode_id, world_comm)
        CALL mp_bcast(celldm, meta_ionode_id, world_comm)
        CALL mp_bcast(at, meta_ionode_id, world_comm)
        CALL mp_bcast(atm, meta_ionode_id, world_comm)
        CALL mp_bcast(amass, meta_ionode_id, world_comm)
        CALL mp_bcast(tau, meta_ionode_id, world_comm)
        CALL mp_bcast(ityp, meta_ionode_id, world_comm)
        CALL mp_bcast(m_loc, meta_ionode_id, world_comm)
        CALL mp_bcast(lrigid, meta_ionode_id, world_comm)
        CALL mp_bcast(epsil, meta_ionode_id, world_comm)
        ! it must be allocated here because nat is read from file
        ALLOCATE (phid(nr1*nr2*nr3,3,3,nat,nat) )
        !
        lrigid1=lrigid

        CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
           at = at / celldm(1)  !  bring at in units of alat

        CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
        CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
        IF (lrigid .AND. (zasr.NE.'no')) &
           CALL set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
     END IF
     IF (lrigid.AND..NOT.lrigid1) CALL errore('q2r_sub', &
           & 'file with dyn.mat. at q=0 should be first of the list',ifile)
     !
     WRITE(stdout,*)
     DO nq = 1,nqs
        WRITE(stdout,'(5x,"q= ",3f12.8)') (q(i,nq),i=1,3)
        lq = .TRUE.
        DO ipol=1,3
           xq = 0.0d0
           DO icar=1,3
              xq = xq + at(icar,ipol) * q(icar,nq) * nr(ipol)
           END DO
           lq = lq .AND. (ABS(NINT(xq) - xq) .LT. eps)
           iq = NINT(xq)
           !
           m(ipol)= MOD(iq,nr(ipol)) + 1
           IF (m(ipol) .LT. 1) m(ipol) = m(ipol) + nr(ipol)
        END DO
        IF (.NOT.lq) CALL errore('q2r_sub','q not allowed',1)

        IF(nc(m(1),m(2),m(3)).EQ.0) THEN
           nc(m(1),m(2),m(3))=1
           IF (lrigid) THEN
              CALL rgd_blk (nr1,nr2,nr3,nat,phiq(1,1,1,1,nq),q(1,nq), &
                  tau,epsil,zeu,bg,omega,celldm(1), .false., -1.d0)
           END IF
           CALL trasl ( phid, phiq, nq, nr1,nr2,nr3, nat, m(1),m(2),m(3))
        ELSE
           WRITE (stdout,'(3i4)') (m(i),i=1,3)
           CALL errore('q2r_sub',' nc already filled: wrong q grid or wrong nr',1)
        END IF
     END DO
     IF (xmldyn) DEALLOCATE(phiq)
  END DO
  !
  ! Check grid dimension
  !
  nq_log = SUM (nc)
  IF (nq_log == nr1*nr2*nr3) THEN
     WRITE (stdout,'(/5x,"q grid ok,  number points =",i5)') nq_log
  ELSE
     CALL errore('q2r_sub',' missing q-point(s)!',1)
  END IF
  !
  ! dyn.mat. FFT (use serial version)
  !
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              CALL cfft3d ( phid (:,j1,j2,na1,na2), &
                      nr1,nr2,nr3, nr1,nr2,nr3, 1, 1 )
              phid(:,j1,j2,na1,na2) = &
                   phid(:,j1,j2,na1,na2) / DBLE(nr1*nr2*nr3)
           END DO
        END DO
     END DO
  END DO
  !
  ! Real space force constants written to file (analytical part)
  !
  IF (xmldyn) THEN
     IF (my_image_id==0) THEN
        IF (lrigid) THEN
           CALL write_dyn_mat_header( filefrc, ntyp, nat, ibrav, nspin_mag,  &
                celldm, at, bg, omega, atm, amass, tau, ityp,   &
                m_loc, nqs, epsil, zeu)
        ELSE
           CALL write_dyn_mat_header( filefrc, ntyp, nat, ibrav, nspin_mag,  &
                celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqs)
        ENDIF
        CALL write_ifc(nr1,nr2,nr3,nat,phid)
     ENDIF
  ELSE 
     IF (meta_ionode) THEN
        iunfrc=find_free_unit()
        OPEN(unit=iunfrc,file=filefrc,status='unknown',form='formatted')
        WRITE(iunfrc,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
        IF (ibrav==0) WRITE (2,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
        DO nt = 1,ntyp
           WRITE(iunfrc,*) nt," '",atm(nt),"' ",amass(nt)
        END DO
        DO na=1,nat
           WRITE(iunfrc,'(2i5,3f18.10)') na,ityp(na),(tau(j,na),j=1,3)
        END DO
        WRITE (iunfrc,*) lrigid
        IF (lrigid) THEN
           WRITE(iunfrc,'(3f15.7)') ((epsil(i,j),j=1,3),i=1,3)
           DO na=1,nat
              WRITE(iunfrc,'(i5)') na
              WRITE(iunfrc,'(3f15.7)') ((zeu(i,j,na),j=1,3),i=1,3)
           END DO
        END IF
        WRITE (iunfrc,'(4i4)') nr1, nr2, nr3
        DO j1=1,3
           DO j2=1,3
              DO na1=1,nat
                 DO na2=1,nat
                    WRITE (iunfrc,'(4i4)') j1,j2,na1,na2
                    nn=0
                    DO m3=1,nr3
                       DO m2=1,nr2
                          DO m1=1,nr1
                             nn=nn+1
                             WRITE (iunfrc,'(3i4,2x,1pe18.11)')   &
                                m1,m2,m3, DBLE(phid(nn,j1,j2,na1,na2))
                           END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
        CLOSE(iunfrc)
     END IF
  END IF
  resi = SUM ( ABS (AIMAG ( phid ) ) )
  IF (resi > eps12) THEN
     WRITE (stdout,"(/,5x,'fft-check warning: sum of imaginary terms = ',&
                                               &es12.6)") resi
  ELSE
     WRITE (stdout,"(/,5x,'fft-check success (sum of imaginary terms < &
                                                          &10^-12)')")
  END IF
  !

  CALL interface_with_tpw(phid, nr1, nr2, nr3, nat, ntyp, lrigid, zeu, &
                            epsil, atm, m_loc, tau, ityp, at, bg, omega)
  !
  DEALLOCATE(nc)
  DEALLOCATE(phid) 
  DEALLOCATE(zeu) 
  DEALLOCATE(tau)
  DEALLOCATE(ityp)
  DEALLOCATE(m_loc)
  IF (.NOT. xmldyn) DEALLOCATE(phiq)
  !
  RETURN
  !
END SUBROUTINE q2r_sub
!
!----------------------------------------------------------------------------
SUBROUTINE trasl( phid, phiq, nq, nr1, nr2, nr3, nat, m1, m2, m3 )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) ::  nr1, nr2, nr3, m1, m2, m3, nat, nq
  COMPLEX(DP), INTENT(IN) :: phiq(3,3,nat,nat,48)
  COMPLEX(DP), INTENT(OUT) :: phid(nr1,nr2,nr3,3,3,nat,nat)
  !
  INTEGER :: j1, j2,  na1, na2
  !
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              phid(m1,m2,m3,j1,j2,na1,na2) = &
                   0.5d0 * (      phiq(j1,j2,na1,na2,nq) +  &
                          CONJG(phiq(j2,j1,na2,na1,nq)))
           END DO
        END DO
     END DO
  END DO
  !
  RETURN
END SUBROUTINE trasl
!----------------------------------------------------------------------
subroutine set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
  !-----------------------------------------------------------------------
  !
  ! Impose ASR - refined version by Nicolas Mounet
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  implicit none
  character(len=10) :: zasr
  integer :: ibrav,nr1,nr2,nr3,nr,m,p,k,l,q,r
  integer :: n,i,j,n1,n2,n3,na,nb,nat,axis,i1,j1,na1
  !
  real(DP) :: sum, zeu(3,3,nat)
  real(DP) :: tau(3,nat), zeu_new(3,3,nat)
  !
  real(DP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat),scal,norm2
  ! temporary vectors and parameters

  ! Initialization.
  ! n is the number of sum rules to be considered (if zasr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2'
  ! and (Oz) if axis='3')
  !
  if((zasr.ne.'simple').and.(zasr.ne.'crystal').and.(zasr.ne.'one-dim') &
                       .and.(zasr.ne.'zero-dim')) then
      call errore('set_zasr','invalid Acoustic Sum Rulei for Z*:' // zasr, 1)
  endif
  if(zasr.eq.'crystal') n=3
  if(zasr.eq.'one-dim') then
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
        call errore('set_zasr','too many directions of &
             &   periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'zasr: rotational axis may be wrong'
     endif
     write(stdout,'("zasr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(zasr.eq.'zero-dim') n=6

  ! Acoustic Sum Rule on effective charges
  !
  if(zasr.eq.'simple') then
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
   else
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
           &  "charges: " , F25.20)') SQRT(norm2)
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
   endif
   !
   !
   return
 end subroutine set_zasr

SUBROUTINE interface_with_tpw(frc_, nr1, nr2, nr3, nat_, ntyp_, has_zstar_, &
                zeu_, epsil_, atm_, m_loc_, tau_, ityp_, at_, bg_, omega_ )
!
!  This routine is used to copy the variables produced by q2r in the variables
!  of the thermo_pw code, avoiding to read the file on disk.
!
USE kinds,  ONLY : DP
USE ifc,    ONLY : frc, atm, zeu, m_loc, epsil_ifc, has_zstar
USE ions_base, ONLY : nat, ntyp=>nsp, tau, ityp
USE cell_base, ONLY : at, bg, omega
USE disp,   ONLY : nq1, nq2, nq3

IMPLICIT NONE
INTEGER,     INTENT(IN) :: nr1, nr2, nr3, nat_, ntyp_, ityp_(nat_)
COMPLEX(DP),    INTENT(IN) :: frc_(nr1*nr2*nr3,3,3,nat_,nat_)
REAL(DP),    INTENT(IN) :: zeu_(3,3,nat_), m_loc_(3,nat_), epsil_(3,3), &
                           at_(3,3), bg_(3,3), omega_, tau_(3,nat_) 
LOGICAL,     INTENT(IN) :: has_zstar_
CHARACTER(LEN=3), INTENT(IN) :: atm_(ntyp_)
INTEGER :: i,j,k,ijk

ALLOCATE (frc(nr1,nr2,nr3,3,3,nat_,nat_))
ALLOCATE (zeu(3,3,nat_))
ALLOCATE (atm(ntyp_))
ALLOCATE (m_loc(3,nat_))

IF (nq1==0 .AND. nq2==0 .AND. nq3==0) THEN
   nq1=nr1
   nq2=nr2
   nq3=nr3
ELSEIF (nq1/=nr1 .OR. nq2/=nr2 .OR. nq3/=nr3) THEN
   CALL errore('interface_with_tpw','nq1, nq2, or nq3 different from &
                          &nr1, nr2, or nr3',1)
ENDIF

DO k=1,nr3
   DO j=1,nr2
      DO i=1,nr1
         ijk=i+(j-1)*nr1+(k-1)*nr1*nr2
         frc(i,j,k,:,:,:,:)=DBLE(frc_(ijk,:,:,:,:))
      ENDDO
   ENDDO
ENDDO
zeu=zeu_
atm=atm_
m_loc=m_loc_
epsil_ifc=epsil_
has_zstar=has_zstar_
!
!  initialize atomic positions and the cell. This is necessary when using 
!  after_disp=.TRUE.
!
at=at_
bg=bg_
nat=nat_
ntyp=ntyp_
omega=omega_
tau=tau_
ityp=ityp_

RETURN
END SUBROUTINE interface_with_tpw

SUBROUTINE clean_ifc_variables()

USE kinds,          ONLY : DP
USE ifc,            ONLY : frc, atm, zeu, m_loc
USE phonon_save,    ONLY : freq_save, z_save
IMPLICIT NONE

IF (ALLOCATED(frc)) DEALLOCATE(frc)
IF (ALLOCATED(atm)) DEALLOCATE(atm)
IF (ALLOCATED(zeu)) DEALLOCATE(zeu)
IF (ALLOCATED(m_loc)) DEALLOCATE(m_loc)
IF (ALLOCATED(freq_save)) DEALLOCATE (freq_save)
IF (ALLOCATED(z_save)) DEALLOCATE (z_save)

RETURN
END SUBROUTINE clean_ifc_variables
