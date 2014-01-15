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
  USE io_global,  ONLY : ionode_id, ionode, stdout
  USE mp_images,  ONLY : intra_image_comm, root_image, my_image_id
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail, &
                         write_dyn_mat_header, write_ifc
  USE control_thermo, ONLY : ldos, flfrc
  USE control_ph, ONLY : ldisp
  !
  IMPLICIT NONE
  !
  INTEGER,  PARAMETER  :: ntypx = 10
  REAL(DP), PARAMETER  :: eps=1.D-5, eps12=1.d-12
  INTEGER              :: nr1, nr2, nr3, nr(3)
  !     dimensions of the FFT grid formed by the q-point grid
  !
  CHARACTER(len=256) :: fildyn, filin
  CHARACTER(len=3)   :: atm(ntypx)
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  LOGICAL :: lq, lrigid, lrigid1, xmldyn
  INTEGER :: m1, m2, m3, m(3), i, j, j1, j2, na1, na2, ipol, nn
  INTEGER :: nat, nq, ntyp, iq, icar, nfile, ifile, nqs, nq_log
  INTEGER :: na, nt
  !
  INTEGER :: ibrav, ierr, nspin_mag, ios
  !
  INTEGER,     ALLOCATABLE ::  nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: phid(:,:,:,:,:)
  REAL(DP),    ALLOCATABLE :: m_loc(:,:)
  !
  REAL(DP) :: celldm(6), at(3,3), bg(3,3)
  REAL(DP) :: q(3,48), omega, xq, amass(ntypx), resi
  REAL(DP) :: epsil(3,3)
  !
  LOGICAL, EXTERNAL :: has_xml
  !
  ! Only one image run this routine 
  !
  IF ( my_image_id /= root_image ) RETURN
  IF (flfrc == ' '.OR. .NOT. ldisp) RETURN
  !
  WRITE(stdout,'(/,2x,76("+"))')
  WRITE(stdout,'(5x,"Computing the interatomic force constants")')
  WRITE(stdout,'(5x,"Writing on file ",a)') TRIM(flfrc)
  WRITE(stdout,'(2x,76("+"),/)')
  !
  xmldyn=has_xml(fildyn)

  IF (ionode) OPEN (unit=1, file=TRIM(fildyn)//'0', status='old', &
                    form='formatted', iostat=ierr)
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF (ionode) THEN
     IF (ierr /= 0) CALL errore('q2r_sub','No grid information on file',1)
     WRITE (stdout,'(/,4x," reading grid info from file ",a)') &
                                                        TRIM(fildyn)//'0'
     READ (1, *) nr1, nr2, nr3
     READ (1, *) nfile
     CLOSE (UNIT=1, STATUS='KEEP')
  ENDIF
  CALL mp_bcast(nr1, ionode_id, intra_image_comm)
  CALL mp_bcast(nr2, ionode_id, intra_image_comm)
  CALL mp_bcast(nr3, ionode_id, intra_image_comm)
  CALL mp_bcast(nfile, ionode_id, intra_image_comm)
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
     WRITE (stdout,*) ' reading force constants from file ',TRIM(filin)

     IF (xmldyn) THEN
        CALL read_dyn_mat_param(filin,ntyp,nat)
        IF (ifile==1) THEN
           ALLOCATE (m_loc(3,nat))
           ALLOCATE (tau(3,nat))
           ALLOCATE (ityp(nat))
           ALLOCATE (zeu(3,3,nat))
        ENDIF
        IF (ifile==1) THEN
           CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
              celldm, at, bg, omega, atm, amass, tau, ityp, &
              m_loc, nqs, lrigid, epsil, zeu )
        ELSE
           CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
              celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqs)
        ENDIF
        ALLOCATE (phiq(3,3,nat,nat,nqs) )
        DO iq=1,nqs
           CALL read_dyn_mat(nat,iq,q(:,iq),phiq(:,:,:,:,iq))
        ENDDO
        CALL read_dyn_mat_tail(nat)
     ELSE
        IF (ionode) &
           OPEN (unit=1,file=filin,status='old',form='formatted',iostat=ierr)
        CALL mp_bcast(ierr, ionode_id, intra_image_comm)
        IF (ierr /= 0) CALL errore('q2r_sub','file '//TRIM(filin)//' missing!',1)
        CALL read_dyn_from_file (nqs, q, epsil, lrigid,  &
                ntyp, nat, ibrav, celldm, at, atm, amass)
        IF (ionode) CLOSE(unit=1)
     ENDIF
     IF (ifile == 1) THEN
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
     DO nq = 1,nqs
        WRITE(stdout,'(a,3f12.8)') ' q= ',(q(i,nq),i=1,3)
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
                  tau,epsil,zeu,bg,omega,-1.d0)
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
     WRITE (stdout,'(/5x,a,i4)') ' q-space grid ok, #points = ',nq_log
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
                      nr1,nr2,nr3, nr1,nr2,nr3, 1 )
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
     IF (lrigid) THEN
        CALL write_dyn_mat_header( flfrc, ntyp, nat, ibrav, nspin_mag,  &
             celldm, at, bg, omega, atm, amass, tau, ityp,   &
             m_loc, nqs, epsil, zeu)
     ELSE
        CALL write_dyn_mat_header( flfrc, ntyp, nat, ibrav, nspin_mag,  &
             celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqs)
     ENDIF
     CALL write_ifc(nr1,nr2,nr3,nat,phid)
  ELSE 
     IF (ionode) THEN
        OPEN(unit=2,file=flfrc,status='unknown',form='formatted')
        WRITE(2,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
        IF (ibrav==0) WRITE (2,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
        DO nt = 1,ntyp
           WRITE(2,*) nt," '",atm(nt),"' ",amass(nt)
        END DO
        DO na=1,nat
           WRITE(2,'(2i5,3f18.10)') na,ityp(na),(tau(j,na),j=1,3)
        END DO
        WRITE (2,*) lrigid
        IF (lrigid) THEN
           WRITE(2,'(3f15.7)') ((epsil(i,j),j=1,3),i=1,3)
           DO na=1,nat
              WRITE(2,'(i5)') na
              WRITE(2,'(3f15.7)') ((zeu(i,j,na),j=1,3),i=1,3)
           END DO
        END IF
        WRITE (2,'(4i4)') nr1, nr2, nr3
        DO j1=1,3
           DO j2=1,3
              DO na1=1,nat
                 DO na2=1,nat
                    WRITE (2,'(4i4)') j1,j2,na1,na2
                    nn=0
                    DO m3=1,nr3
                       DO m2=1,nr2
                          DO m1=1,nr1
                             nn=nn+1
                             WRITE (2,'(3i4,2x,1pe18.11)')   &
                                m1,m2,m3, DBLE(phid(nn,j1,j2,na1,na2))
                           END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
        CLOSE(2)
     END IF
  END IF
  resi = SUM ( ABS (AIMAG ( phid ) ) )
  IF (resi > eps12) THEN
     WRITE (stdout,"(/5x,' fft-check warning: sum of imaginary terms = ',es12.6)") resi
  ELSE
     WRITE (stdout,"(/5x,' fft-check success (sum of imaginary terms < 10^-12)')")
  END IF
  !
  DEALLOCATE(nc)
  DEALLOCATE(phid) 
  DEALLOCATE(zeu) 
  DEALLOCATE(tau)
  DEALLOCATE(ityp)
  IF (xmldyn) THEN
     DEALLOCATE(m_loc)
  ELSE 
     DEALLOCATE(phiq)
  ENDIF
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
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: i,j,na,nat
  REAL(DP) :: zeu_u(3,3,nat)
  REAL(DP) :: zeu_v(3,3,nat)
  REAL(DP) :: scal
  !
  !
  scal=0.0d0
  DO i=1,3
    DO j=1,3
      DO na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      ENDDO
    ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE sp_zeu
