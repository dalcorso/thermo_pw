!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE ev_sub(vmin,b0,b01,b02,emin_out,inputfile)
!-----------------------------------------------------------------------
!
!      fit of E(v) or H(V) at finite pressure to an equation of state (EOS)
!
!      Interactive input:
!         au or Ang
!         structure
!         equation of state
!         input data file
!         output data file
!
!      Input data file format for cubic systems:
!         a0(1)  Etot(1)
!         ...
!         a0(n)  Etot(n)
!      where a0 is the lattice parameter (a.u. or Ang)
!      Input data file format for noncubic (e.g. hexagonal) systems:
!         V0(1)  Etot(1)
!         ...
!         V0(n)  Etot(n)
!      where V0 is the unit-cell volume (a.u.^3 or Ang^3)
!      e.g. for an hexagonal cell,
!         V0(i)  = sqrt(3)/2 * a^2 * c    unit-cell volume
!         Etot(i)= min Etot(c)   for the given volume V0(i)
!      Etot in atomic (Rydberg) units
!
!      Output data file format  for cubic systems:
!      # a0=... a.u., K0=... kbar, dk0=..., d2k0=... kbar^-1, Emin=... Ry
!      # a0=... Ang,  K0=... GPa , V0=... (a.u.)^3, V0 = Ang^3
!         a0(1)  Etot(1) Efit(1)  Etot(1)-Efit(1)  Pfit(1)  Enth(1)
!         ...
!         a0(n)  Etot(n) Efit(n)  Etot(n)-Efit(n)  Pfit(n)  Enth(n)
!      Output data file format  for noncubic systems:
!      # V0=...(a.u.)^3, K0=... kbar, dk0=..., d2k0=... kbar^-1, Emin=... Ry
!      # V0=...Ang^3,  K0=... GPa
!         V0(1)  Etot(1) Efit(1)  Etot(1)-Efit(1)  Pfit(1)  Enth(1)
!         ...
!         V0(n)  Etot(n) Efit(n)  Etot(n)-Efit(n)  Pfit(n)  Enth(n)
!      where
!            a0(i), V0(i), Etot(i) as in input
!            Efit(i) is the fitted value from the EOS
!            Pfit(i) is the corresponding pressure from the EOS (GPa)
!            Enth(i)=Efit(i)+Pfit(i)*V0(i) is the enthalpy (Ry)
!!
   USE kinds, ONLY: DP
   USE constants, ONLY: bohr_radius_angs
   USE ev_xml,    ONLY : write_evdata_xml
   USE ev_mod,    ONLY : initialize_data_ev, find_minimum, emin
   USE control_pressure, ONLY : pressure_kb
   USE mp,        ONLY : mp_bcast
   USE io_global, ONLY : ionode, ionode_id, stdout
   USE mp_images, ONLY : my_image_id, root_image, intra_image_comm

   IMPLICIT NONE
   REAL(DP), INTENT(OUT)  :: vmin, b0, b01, b02, emin_out
   CHARACTER(LEN=*) :: inputfile
   INTEGER, PARAMETER:: nmaxpar=4, nmaxpt=100
   INTEGER :: npar,npt,ieos, ierr
   CHARACTER(LEN=3) :: bravais, au_unit
   REAL(DP) :: par(nmaxpar), v0(nmaxpt), etot(nmaxpt), &
               efit(nmaxpt), fac, chisq, a
   LOGICAL :: in_angstrom
   INTEGER :: iu_ev
   INTEGER :: find_free_unit
   CHARACTER(LEN=256) :: filin, fileout
  !
  IF (my_image_id /= root_image) RETURN

  IF ( ionode ) THEN
     iu_ev=find_free_unit()
     OPEN(UNIT=iu_ev, FILE=TRIM(inputfile), STATUS='OLD', FORM='FORMATTED')

     READ(iu_ev,'(a)') au_unit
     in_angstrom = au_unit=='Ang' .OR. au_unit=='ANG' .OR. &
                   au_unit=='ang'
     READ(iu_ev, '(a)') bravais
!
     IF (bravais=='fcc'.OR.bravais=='FCC') THEN
        fac = 0.25d0
     ELSEIF(bravais=='bcc'.OR.bravais=='BCC') THEN
        fac = 0.50d0
     ELSEIF(bravais=='sc'.OR.bravais=='SC') THEN
        fac = 1.0d0
     ELSEIF(bravais=='noncubic'.OR.bravais=='NONCUBIC'.OR.  &
            bravais=='hex'.OR.bravais=='HEX' ) THEN
         fac = 0.0_DP ! not used
     ELSE
        CALL errore('ev_sub','ev: unexpected lattice '//TRIM(bravais), 1)
     ENDIF
!
     READ (iu_ev,*) ieos
     IF (ieos==1 .OR. ieos==4) THEN
        npar=3
     ELSEIF(ieos==2) THEN
        npar=4
     ELSE
        CALL errore('ev_sub', 'Unexpected eq. of state', ieos)
     ENDIF
     READ(iu_ev, '(a)') filin
     READ(iu_ev, '(a)') fileout

     CLOSE(iu_ev)
!
!  reading the data
!
     OPEN(UNIT=iu_ev,FILE=TRIM(filin),STATUS='old',FORM='formatted',IOSTAT=ierr)
     IF (ierr/=0) THEN
        ierr= 1 
        GO TO 99
     END IF
  10 CONTINUE
     emin=1.d10
     DO npt=1,nmaxpt
        IF (bravais=='noncubic' .OR. bravais=='NONCUBIC' .OR. &
            bravais=='hex' .OR. bravais=='HEX' ) THEN
           READ(iu_ev,*,ERR=10,END=20) v0(npt), etot(npt)
           IF (in_angstrom) v0(npt)=v0(npt)/bohr_radius_angs**3
        ELSE
           READ(iu_ev,*,ERR=10,END=20) a, etot(npt)
           IF (in_angstrom) a = a/bohr_radius_angs
           v0  (npt) = fac*a**3
        ENDIF
        IF (etot(npt)<emin) THEN
           par(1) = v0(npt)
           emin = etot(npt)
        ENDIF
     ENDDO
     npt = nmaxpt+1
  20 CLOSE(iu_ev)
     npt = npt-1
!
! par(1) = V, Volume of the unit cell in (a.u.^3)
! par(2) = B, Bulk Modulus (in KBar)
! par(3) = dB/dP (adimensional)
! par(4) = d^2B/dP^2 (in KBar^(-1), used only by 2nd order formulae)
!
     par(2) = 500.0d0
     par(3) = 5.0d0
     par(4) = -0.01d0
!
     CALL initialize_data_ev(v0, etot, npt, ieos)
!
     CALL find_minimum (npar,par,chisq)
!
     CALL write_results &
          (npt,in_angstrom,fac,v0,etot,efit,ieos,par,npar,emin,chisq,fileout)
!
     CALL write_evdata_xml  &
          (npt,fac,v0,etot,efit,ieos,par,npar,emin,pressure_kb,&
                                                        chisq,fileout,ierr)
     IF (ierr /= 0) GO TO 99
  ENDIF
99  CALL mp_bcast ( ierr, ionode_id, intra_image_comm )

  IF ( ierr == 1) THEN
     CALL errore( 'ev_sub', 'file '//trim(filin)//' cannot be opened', ierr )
  ELSEIF ( ierr == 2 ) THEN
     CALL errore( 'ev_sub', 'file '//trim(fileout)//' cannot be opened', ierr )
  ELSEIF ( ierr == 11 ) THEN
     CALL errore( 'ev_sub', 'no free units to write ', ierr )
  ELSEIF ( ierr == 12 ) THEN
     CALL errore( 'ev_sub', 'error opening the xml file ', ierr )
  ENDIF
  CALL mp_bcast(par, ionode_id, intra_image_comm)
  CALL mp_bcast(emin, ionode_id, intra_image_comm)
    
  vmin=par(1)
  b0=par(2)
  b01=par(3)
  b02=par(4)
  emin_out=emin

  RETURN
  END SUBROUTINE ev_sub
!
!-----------------------------------------------------------------------
SUBROUTINE ev_sub_nodisk(vmin,b0,b01,b02,emin_out)
!-----------------------------------------------------------------------
!
!  This routine is similar to ev_sub, but it receives the input data
!  directly from the shared variables without the need to write on
!  disk. It does not write any output. It can be called by any processor.
!  This routine is not parallel. The CPU that calls it receive the result.
!  
!  Before calling this routine the user must set in the module control_ev
!  ieos : the equation of state to use
!         1 - Birch-Murnaghan first order
!         2 - Birch-Murnaghan third order
!         3 - Keane
!         4 - Murnaghan
!  npt : the number of points
!  v0(npt)  : the volume for each point
!  e0(npt)  : the energy or enthalpy for each point
!
!  v0 and e0 must be allocated and deallocated by the user of this routine
!
USE kinds, ONLY: DP
USE control_ev, ONLY : ieos, npt, v0, etot => e0
USE ev_mod,    ONLY : initialize_data_ev, find_minimum, emin

IMPLICIT NONE
REAL(DP), INTENT(OUT)  :: vmin, b0, b01, b02, emin_out
INTEGER, PARAMETER:: nmaxpar=4
INTEGER :: npar,ipt,ierr
REAL(DP) :: par(nmaxpar), chisq

IF (ieos==1 .OR. ieos==4) THEN
   npar=3
ELSEIF(ieos==2) THEN
   npar=4
ELSE
   CALL errore('ev_sub_nodisk', 'Unexpected eq. of state', ieos)
ENDIF
!
!  find emin and the initial guess for the volume
!
emin=1.d10
DO ipt=1,npt
   IF (etot(ipt)<emin) THEN
      par(1) = v0(ipt)
      emin = etot(ipt)
   ENDIF
ENDDO
!
! par(1) = V, Volume of the unit cell in (a.u.^3)
! par(2) = B, Bulk Modulus (in KBar)
! par(3) = dB/dP (adimensional)
! par(4) = d^2B/dP^2 (in KBar^(-1), used only by 2nd order formulae)
!
par(2) = 500.0d0
par(3) = 5.0d0
par(4) = -0.01d0

CALL initialize_data_ev(v0, etot, npt, ieos)
   !
CALL find_minimum (npar,par,chisq)
    
vmin=par(1)
b0=par(2)
b01=par(3)
b02=par(4)
emin_out=emin

RETURN
END SUBROUTINE ev_sub_nodisk
!
!-----------------------------------------------------------------------
SUBROUTINE write_results &
      (npt,in_angstrom,fac,v0,etot,efit,istat,par,npar,emin,chisq,filout)
!-----------------------------------------------------------------------
!
    USE kinds, ONLY: DP
    USE constants, ONLY : bohr_radius_angs, ry_kbar 
    USE control_pressure, ONLY : pressure_kb, pressure
    USE eos, ONLY : eos_press
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npt, istat, npar
    REAL(DP), INTENT(IN):: v0(npt), etot(npt), efit(npt), emin, chisq, fac
    REAL(DP), INTENT(INOUT):: par(npar)
    REAL(DP), PARAMETER :: gpa_kbar = 10.0_dp
    LOGICAL, INTENT(IN) :: in_angstrom
    CHARACTER(LEN=256), INTENT(IN) :: filout
    !
    REAL(DP) :: p(npt), epv(npt)
    REAL(DP) :: k0, dk0, d2k0, vol0
    INTEGER :: i, iun, ierr
    INTEGER :: find_free_unit
    LOGICAL :: exst


    IF (filout/=' ') THEN
       iun=find_free_unit()
       INQUIRE(FILE=TRIM(filout),EXIST=exst)
       IF (exst) PRINT '(5x,"Beware: file ",A," will be overwritten")',&
                  TRIM(filout)
       OPEN(UNIT=iun,FILE=TRIM(filout),FORM='formatted',STATUS='unknown', &
              IOSTAT=ierr)
       IF (ierr/=0) THEN
          ierr= 2 
          GO TO 99
       END IF
    ELSE
       iun=6
    ENDIF

    IF (istat==1) THEN
       WRITE(iun,'("# equation of state: birch 3st order.  chisq = ", &
                   & d10.4)') chisq
    ELSEIF(istat==2) THEN
       WRITE(iun,'("# equation of state: birch 4rd order.  chisq = ", &
                   & d10.4)') chisq
    ELSEIF(istat==4) THEN
       WRITE(iun,'("# equation of state: murnaghan.        chisq = ", &
                   & d10.4)') chisq
    ENDIF
!
    vol0 = par(1)
    k0   = par(2)/ry_kbar ! converts k0 to Ry atomic units...
    dk0  = par(3)
    d2k0 = par(4)*ry_kbar ! and d2k0/dp2 to (Ry a.u.)^(-1)
!
    DO i=1,npt
       CALL eos_press(istat, v0(i), p(i), vol0, k0, dk0, d2k0)
    ENDDO

    DO i=1,npt
       epv(i) = etot(i) + p(i)*v0(i) 
    ENDDO
!
!   from now on the pressure is in kbar
!
    p(1:npt)=p(1:npt)*ry_kbar

    IF ( fac /= 0.0_dp ) THEN
! cubic case
       IF (pressure_kb /= 0.0_DP) THEN
          WRITE(iun,'("# a0 =",f8.4," a.u., k0 =",i5," kbar, dk0 =", &
                    &f6.2," d2k0 =",f7.3," Hmin =",f11.5)') &
                  (par(1)/fac)**(1d0/3d0), int(par(2)), par(3), par(4), emin

       ELSE
          WRITE(iun,'("# a0 =",f8.4," a.u., k0 =",i5," kbar, dk0 =", &
                    &f6.2," d2k0 =",f7.3," emin =",f11.5)') &
                  (par(1)/fac)**(1d0/3d0), int(par(2)), par(3), par(4), emin
       ENDIF
       WRITE(iun,'("# a0 =",f9.5," Ang, k0 =", f6.1," GPa,  V0 = ", &
                  & f7.3," (a.u.)^3,  V0 =", f7.3," A^3 ",/)') &
           & (par(1)/fac)**(1d0/3d0)*bohr_radius_angs, par(2)/gpa_kbar, &
             par(1), par(1)*bohr_radius_angs**3

       WRITE(iun,'(73("#"))')
       IF (pressure_kb /= 0.0_DP) THEN
          WRITE(iun,'("# Lat.Par", 4x, "(E+pV)_calc", 2x, "(E+pV)_fit", 3x, &
            & "(E+pV)_diff", 2x, "Pressure", 6x, "Enthalpy")')
       ELSE
          WRITE(iun,'("# Lat.Par", 7x, "E_calc", 8x, "E_fit", 7x, &
             & "E_diff", 4x, "Pressure", 6x, "Enthalpy")')
       ENDIF
       IF (in_angstrom) THEN
          WRITE(iun,'("# Ang", 13x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
          WRITE(iun,'(73("#"))')
          WRITE(iun,'(f9.5,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
                ( (v0(i)/fac)**(1d0/3d0)*bohr_radius_angs, etot(i), efit(i),  &
                etot(i)-efit(i), (p(i)+pressure_kb)/gpa_kbar, &
                                epv(i), i=1,npt )
       ELSE
          WRITE(iun,'("# a.u.",12x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
          WRITE(iun,'(73("#"))')
          WRITE(iun,'(f9.5,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
               ( (v0(i)/fac)**(1d0/3d0), etot(i), efit(i),  &
               etot(i)-efit(i), (p(i)+pressure_kb)/gpa_kbar,  &
                                 epv(i), i=1,npt )
       ENDIF

    ELSE
! noncubic case
       IF (pressure_kb /= 0.0_DP) THEN
          WRITE(iun,'("# V0 =",f8.2," a.u.^3,  k0 =",i5," kbar,  dk0 =", &
                    & f6.2,"  d2k0 =",f7.3,"  Hmin =",f11.5)') &
                    & par(1), int(par(2)), par(3), par(4), emin
       ELSE
          WRITE(iun,'("# V0 =",f8.2," a.u.^3,  k0 =",i5," kbar,  dk0 =", &
                    & f6.2,"  d2k0 =",f7.3,"  emin =",f11.5)') &
                    & par(1), int(par(2)), par(3), par(4), emin
       ENDIF

       WRITE(iun,'("# V0 =",f8.2,"  Ang^3,  k0 =",f6.1," GPa"/)') &
                    & par(1)*bohr_radius_angs**3, par(2)/gpa_kbar

       WRITE(iun,'(74("#"))')
       IF (pressure_kb /= 0.0_DP) THEN
          WRITE(iun,'("# Vol.", 6x, "(E+pV)_calc", 2x, "(E+pV)_fit", 4x, &
            & "(E+pV)_diff", 2x, "Pressure", 6x, "Enthalpy")')
       ELSE
          WRITE(iun,'("# Vol.", 8x, "E_calc", 8x, "E_fit", 7x, &
            & "E_diff", 4x, "Pressure", 6x, "Enthalpy")')
       ENDIF
       IF (in_angstrom) THEN
          WRITE(iun,'("# Ang^3", 9x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
          WRITE(iun,'(74("#"))')
          WRITE(iun,'(f8.2,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              ( v0(i)*bohr_radius_angs**3, etot(i), efit(i),  &
               etot(i)-efit(i), (p(i)+pressure_kb)/gpa_kbar, epv(i), i=1,npt )
       ELSE
          WRITE(iun,'("# a.u.^3",8x, "Ry", 11x, "Ry", 12x, &
             & "Ry", 8x, "GPa", 11x, "Ry")')
          WRITE(iun,'(74("#"))')
          WRITE(iun,'(f8.2,2x,f12.5, 2x,f12.5, f12.5, 3x, f8.2, 3x,f12.5)') &
              ( v0(i), etot(i), efit(i),  &
               etot(i)-efit(i), (p(i)+pressure_kb)/gpa_kbar, epv(i), i=1,npt )
       ENDIF

     ENDIF
     IF(filout/=' ') CLOSE(UNIT=iun,STATUS='KEEP')
 99  RETURN
 END SUBROUTINE write_results
