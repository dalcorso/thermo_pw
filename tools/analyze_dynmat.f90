!
! Copyright (C) 2025 Dal Corso Andrea
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
PROGRAM analyze_dynmat
  !----------------------------------------------------------------------------
  !! A small utility that reads two dynamical matrix files (either
  !! xml or plain text), the first at the Gamma point of the Brillouin zone
  !! the second at a finite q, applies the acoustic sum rule to both, and 
  !! diagonalize the two matrices.
  !
  !! Syntax:  
  !!   \(\texttt{analyze_dynamt.x}\) < filein 
  !
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : amu_ry
  USE parameters,         ONLY : ntypx
  USE mp,                 ONLY : mp_bcast
  USE mp_global,          ONLY : mp_startup, mp_global_end
  USE mp_world,           ONLY : world_comm
  USE io_global,          ONLY : ionode_id, ionode, stdout, stdin
  USE environment,        ONLY : environment_start, environment_end
  ! for reading the dyn.mat.
  USE cell_base,          ONLY : at, bg, celldm, ibrav, omega
  USE ions_base,          ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
  ! as above, unused here
  USE control_ph,         ONLY : xmldyn
  USE noncollin_module,   ONLY : m_loc, nspin_mag
  !
  ! for non-xml file only:
  USE dynamicalq,         ONLY : dq_phiq => phiq, dq_tau => tau, &
                                 dq_ityp => ityp, dq_zeu => zeu 
  ! fox xml files only
  USE io_dyn_mat,         ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                                 read_dyn_mat, read_dyn_mat_tail, &
                                 write_dyn_mat_header
  USE rigid,  ONLY : nonanal
  IMPLICIT NONE
  !
  CHARACTER(len=7),PARAMETER :: code="analyze_dynmat"
  CHARACTER(len=256) :: fildyng, fildynq, filout
  INTEGER :: ierr
  !
  INTEGER       :: nqs, imq, nqq
  REAL(DP)      :: xq(3), xqs(3,48), xqs_(3,48), epsilon(3,3), xqhat(3)
  !
  LOGICAL :: xmldyng, xmldynq
  LOGICAL, EXTERNAL :: has_xml
  !
  COMPLEX(DP),ALLOCATABLE :: phig(:,:,:,:), phiq(:,:,:,:), phg(:,:), phq(:,:)
  REAL(DP), ALLOCATABLE :: w2g(:), w2q(:), m_loc_(:,:), tau_(:,:), &
                           zeu_(:,:,:), zeu(:,:,:)
  INTEGER, ALLOCATABLE :: ityp_(:), itau(:)
  COMPLEX(DP) :: asum
  LOGICAL :: lrigid, lrigid_
  CHARACTER(LEN=10) :: zasr
  INTEGER :: i, j, icar, jcar, na, nb
  INTEGER :: nat_, ntyp_, nqs_, ibrav_, nspin_mag_
  REAL(DP) :: celldm_(6), at_(3,3), bg_(3,3), omega_, epsilon_(3,3), &
              amass_(ntypx), xqg(3), xmod
  !
  CALL mp_startup()
  CALL environment_start(code)
  !
  WRITE(stdout,'(5x,"Dynamical matrix at gamma?")') 
  READ(stdin,'(a)') fildyng
  WRITE(stdout,'(a)') TRIM(fildyng)
  WRITE(stdout,'(5x,"Dynamical matrix at q? ")') 
  READ(stdin,'(a)') fildynq
  WRITE(stdout,'(a)') TRIM(fildynq)
  WRITE(stdout,'(5x,"zasr? ")') 
  READ(stdin,'(a)') zasr
  WRITE(stdout,'(a)') TRIM(zasr)
  !
  ! check input
  IF (fildyng == ' ')  CALL errore (code,' bad fildyn at gamma',1)
  IF (fildynq == ' ')  CALL errore (code,' bad fildyn at q',1)
  xmldyng=has_xml(fildyng)
  xmldynq=has_xml(fildynq)
  IF (xmldyng .NEQV. xmldynq) CALL errore(code,&
                            'matrices must have the same format',1)
  xmldyn=xmldyng
  imq=0
!
!  Reading the dynamical matrix at gamma
!
  XML_FORMAT_READ : &
  IF (xmldyng) THEN
    ! read params
    CALL read_dyn_mat_param(fildyng,ntyp,nat)
    ALLOCATE(m_loc(3,nat))
    ALLOCATE(tau(3,nat))
    ALLOCATE(ityp(nat))
    ALLOCATE(zeu(3,3,nat))
    ALLOCATE(phig(3,3,nat,nat))
    ! read system information
    CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
                             celldm, at, bg, omega, atm, amass, tau, ityp, &
                             m_loc, nqs, lrigid, epsilon, zeu )
    ! read dyn.mat.
    CALL read_dyn_mat(nat,1,xqg,phig)
    ! close file
    CALL read_dyn_mat_tail(nat)
    !
  ELSE XML_FORMAT_READ
    ! open file
    IF (ionode) OPEN (unit=1, file=fildyng,status='old',form='formatted',&
                                                            iostat=ierr)
    CALL mp_bcast(ierr, ionode_id,world_comm)
    IF (ierr /= 0) CALL errore(code,'file '//TRIM(fildyng)//' missing!',1)
    ! read everything, this use global variables
    ntyp = ntypx
    CALL read_dyn_from_file (nqs, xqs, epsilon, lrigid,  &
        ntyp, nat, ibrav, celldm, at, atm, amass)
    !
    IF (ionode) CLOSE(unit=1)
    !
    xqg = xqs(:,1)
    ALLOCATE(phig(3,3,nat,nat))
    ALLOCATE(tau(3,nat))
    ALLOCATE(ityp(nat))
    ALLOCATE(zeu_(3,3,nat))
    phig  = dq_phiq(:,:,:,:,1)
    tau =  dq_tau
    ityp = dq_ityp
    zeu = dq_zeu
    amass = amass/amu_ry
!
!   These are allocated by read_dyn_from_file
!
    DEALLOCATE(dq_phiq) 
    DEALLOCATE(dq_tau) 
    DEALLOCATE(dq_ityp)
    DEALLOCATE(dq_zeu)
    !
  ENDIF XML_FORMAT_READ
  !
  ! regenerate the lattice
  !
  CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  at = at / celldm(1)  !  bring at in units of alat
  CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
  CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
!
!  Reading the dynamical matrix at gamma
!

  XML_FORMAT_READ_Q : &
  IF (xmldyng) THEN
    ! read params
    CALL read_dyn_mat_param(fildynq,ntyp_,nat_)
    IF ((nat_/=nat).OR.(ntyp_/=ntyp)) CALL errore(code,'different systems',1)
    ALLOCATE(m_loc_(3,nat))
    ALLOCATE(tau_(3,nat))
    ALLOCATE(ityp_(nat))
    ALLOCATE(zeu_(3,3,nat))
    ALLOCATE(phiq(3,3,nat,nat))
    ! read system information
    CALL read_dyn_mat_header(ntyp, nat, ibrav_, nspin_mag_, &
                             celldm_, at_, bg_, omega_, atm, amass, tau_, &
                             ityp_, m_loc_, nqs_, lrigid_, epsilon_, zeu_ )
    IF (ABS(omega_-omega)>1.D-9) CALL errore(code,'different omega',1)
    ! read dyn.mat.
    CALL read_dyn_mat(nat,1,xq,phiq)
    ! close file
    CALL read_dyn_mat_tail(nat)
    !
  ELSE XML_FORMAT_READ_Q
    ! open file
    IF (ionode) OPEN(unit=1, file=fildynq,status='old',form='formatted',&
                                                                iostat=ierr)
    CALL mp_bcast(ierr, ionode_id,world_comm)
    IF (ierr /= 0) CALL errore(code,'file '//TRIM(fildynq)//' missing!',1)
    ! read everything, this use global variables
    ntyp = ntypx
    CALL read_dyn_from_file (nqs_, xqs_, epsilon_, lrigid_,  &
        ntyp_, nat_, ibrav_, celldm_, at_, atm, amass)
    IF (ibrav/=ibrav_) CALL errore(code,'wrong ibrav',1)
    IF ((ABS(celldm(1)-celldm_(1))<1.D-7).OR. &
        (ABS(celldm(2)-celldm_(2))<1.D-7).OR. &
        (ABS(celldm(3)-celldm_(3))<1.D-7).OR. &
        (ABS(celldm(4)-celldm_(4))<1.D-7).OR. &
        (ABS(celldm(5)-celldm_(5))<1.D-7).OR. &
        (ABS(celldm(6)-celldm_(6))<1.D-7 )) CALL errore(code,'wrong celldm',1)
    !
    IF (ionode) CLOSE(unit=1)
    !
    xq = xqs(:,1)
    ALLOCATE(phiq(3,3,nat,nat))
    ALLOCATE(tau_(3,nat))
    ALLOCATE(ityp_(nat))
    ALLOCATE(zeu_(3,3,nat))
    phiq  = dq_phiq(:,:,:,:,1)
    tau_ =  dq_tau
    ityp_ = dq_ityp
    zeu_=dq_zeu
    amass_ = amass/amu_ry
!
!   These are allocated by read_dyn_from_file
!
    DEALLOCATE(dq_phiq) 
    DEALLOCATE(dq_tau) 
    DEALLOCATE(dq_ityp)
    DEALLOCATE(dq_zeu)
    !
  ENDIF XML_FORMAT_READ_Q

  filout='analyze_dynmat.g.out'
  XML_FORMAT_WRITE : &
  IF (xmldyn) THEN
     nqq=nqs
     IF (imq==0) nqq=2*nqs
     CALL write_dyn_mat_header( filout, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqq)
  ELSE XML_FORMAT_WRITE
      OPEN (unit=1, file=filout,status='unknown',form='formatted',iostat=ierr)
      IF (ierr /= 0) CALL errore(code,'opening output file',1)
      CALL write_old_dyn_mat_head(1)
  ENDIF XML_FORMAT_WRITE
!
!  If an insulator add the nonanalytic part to the matrix at gamma  
!
  IF (lrigid) THEN
     WRITE(stdout,'(5x,"Adding the nonanalytic part")')
     ALLOCATE (itau(nat))
     DO na=1,nat
        itau(na)=na
     END DO
     xmod=SQRT(xq(1)**2+xq(2)**2+xq(3)**2)
     xqhat=xq/xmod
     WRITE(6,'(5x,"qhat=",3f20.10)') xqhat(:)
     CALL nonanal ( nat, nat, itau, epsilon, xqhat, zeu, omega, phig )
     DEALLOCATE (itau)
  ENDIF

!
!   Apply the acoustic sum rule to both matrices
!  
  IF (zasr=='simple') THEN
     DO icar=1,3 
        DO jcar=1,3 
           DO na=1, nat
              asum=(0.0_DP,0.0_DP)      
              DO nb=1,nat
                 asum=asum+phig(icar,jcar,na,nb)
              ENDDO
              phig(icar,jcar,na,na)= phig(icar,jcar,na,na) - asum
              phiq(icar,jcar,na,na)= phiq(icar,jcar,na,na) - asum
           ENDDO
        ENDDO
     ENDDO
  ELSEIF (zasr=='no') THEN
  ELSE
     CALL errore(code,'unknown zasr',1)
  ENDIF
!
!  recompact them
!
  ALLOCATE(phg(3*nat, 3*nat))
  ALLOCATE(phq(3*nat, 3*nat))
  DO i = 1, 3 * nat
    na = (i - 1) / 3 + 1
    icar = i - 3 * (na - 1)
    DO j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcar = j - 3 * (nb - 1)
        phg (i, j) = phig(icar, jcar, na, nb)
        phq (i, j) = phiq(icar, jcar, na, nb)
    ENDDO
  ENDDO
!
!  Diagonalize and write on output the frequencies
!   
  ALLOCATE(w2g(3*nat))
  CALL dyndia (xqg, 3*nat, nat, ntyp, ityp, amass, 1, phg, w2g)
  IF (.NOT.xmldyn) THEN
     WRITE(1, '(/,3a,/)') "File generated with q2qstar.x from '", &
          TRIM(fildyng), "'" ! <-- to prevent crash with old versions of q2r.x
     CLOSE(1)
  ENDIF
  filout='analyze_dynmat.q.out'
  XML_FORMAT_WRITE_Q : &
  IF (xmldyn) THEN
     nqq=nqs_
     IF (imq==0) nqq=2*nqs_
     CALL write_dyn_mat_header( filout, ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, nqq)
  ELSE XML_FORMAT_WRITE_Q
      OPEN (unit=1, file=filout,status='unknown',form='formatted',iostat=ierr)
      IF (ierr /= 0) CALL errore(code,'opening output file',1)
      CALL write_old_dyn_mat_head(1)
  ENDIF XML_FORMAT_WRITE_Q

  ALLOCATE(w2q(3*nat))
  CALL dyndia (xq, 3*nat, nat, ntyp, ityp, amass, 1, phq, w2q)
  IF (.NOT.xmldyn) THEN
     WRITE(1, '(/,3a,/)') "File generated with q2qstar.x from '", &
          TRIM(fildynq), "'" ! <-- to prevent crash with old versions of q2r.x
     CLOSE(1)
  ENDIF
  !
  DEALLOCATE(phig)
  DEALLOCATE(phg)
  DEALLOCATE(w2g)
  DEALLOCATE(tau) 
  DEALLOCATE(ityp)
  DEALLOCATE(zeu)
  IF (xmldyn) THEN
     DEALLOCATE(m_loc)
     DEALLOCATE(m_loc_)
  ENDIF

  DEALLOCATE(phiq) 
  DEALLOCATE(phq) 
  DEALLOCATE(w2q)
  DEALLOCATE(tau_)
  DEALLOCATE(ityp_)
  DEALLOCATE(zeu_) 

  CALL environment_end( code )
  CALL mp_global_end ()

END PROGRAM analyze_dynmat
