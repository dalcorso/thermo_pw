!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ev_xml
!
!   This module contains routines to write the information obtained by the
!   ev.x program in an xml file.
!
USE kinds,     ONLY : DP

IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: write_evdata_xml, birch, keane

  INTEGER :: iunout
  !

  CONTAINS
!-----------------------------------------------------------------------
    SUBROUTINE write_evdata_xml &
        (npt,fac,v0,etot,efit,istat,par,npar,emin,pressure_kb,&
                                                     chisq,filout, ierr)
!-----------------------------------------------------------------------
!
  USE kinds,     ONLY : dp
  USE constants, ONLY : ry_kbar, bohr_radius_angs
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npt, istat, npar
  REAL(DP), INTENT(IN):: v0(npt), etot(npt), efit(npt), emin, chisq, fac
  REAL(DP), INTENT(In):: par(npar), pressure_kb
  CHARACTER(LEN=256), INTENT(IN) :: filout
  INTEGER, INTENT(OUT) :: ierr
  !
  INTEGER :: iunout = 11
  REAL(DP) :: p(npt), volume(2), a0(2), alldata(6,npt)
  INTEGER :: i, iun
  CHARACTER(len=256) :: filename

  IF (filout/=' ') THEN
     filename = TRIM(filout) // '.xml'
  ELSE
     filename = 'ev.xml'
  ENDIF
  !
  ! ... open XML descriptor
  !
  OPEN ( UNIT=iunout, FILE = TRIM( filename ), FORM='formatted', IOSTAT = ierr )
  IF ( ierr /= 0 ) THEN
     WRITE (6,*) 'Failed opening file ' // TRIM(filename)
     RETURN
  END IF

  WRITE (iunout,'("<xml>")')
  WRITE (iunout,'("<EQUATIONS_OF_STATE>")')
  WRITE (iunout,'("<EQUATION_TYPE>")')
  IF (istat==1) THEN
     WRITE (iunout,'("Birch 1st order")')
  ELSEIF (istat==2) THEN
     WRITE (iunout,'("Birch 2nd order")')
  ELSEIF (istat==3) THEN
     WRITE (iunout,'("Keane")')
  ELSEIF (istat==4) THEN
     WRITE (iunout,'("Murnaghan")')
  ENDIF
  WRITE (iunout,'("</EQUATION_TYPE>")')
  WRITE (iunout,'("<CHI_SQUARE>")')
  WRITE (iunout,'(1pe25.12)') chisq
  WRITE (iunout,'("</CHI_SQUARE>")')
  WRITE (iunout,'("</EQUATIONS_OF_STATE>")')

  IF (istat==1 .or. istat==2) THEN
     DO i=1,npt
        p(i)=birch(v0(i)/par(1),par(2),par(3),par(4))
     ENDDO
  ELSE
     DO i=1,npt
        p(i)=keane(v0(i)/par(1),par(2),par(3),par(4))
     ENDDO
  ENDIF

  DO i=1,npt
     alldata (1,i) = v0(i)
     alldata (2,i) = etot(i) 
     alldata (3,i) = efit(i)
     alldata (4,i) = etot(i) - efit(i)
     alldata (5,i) = p(i) 
     alldata (6,i) = etot(i) + p(i) * v0(i) / ry_kbar
  ENDDO

  WRITE (iunout,'("<EQUATIONS_PARAMETERS>")')

  volume(1)=par(1)
  volume(2)=par(1)*bohr_radius_angs**3
  WRITE (iunout, '("<EQUILIBRIUM_VOLUME_AU_A>")')
  WRITE (iunout, '(2(1pe25.15))') volume(:)
  WRITE (iunout, '("</EQUILIBRIUM_VOLUME_AU_A>")')
  WRITE (iunout, '("<BULK_MODULUS_KBAR>")')
  WRITE (iunout, '(1pe25.15)') par(2)
  WRITE (iunout, '("</BULK_MODULUS_KBAR>")')
  WRITE (iunout, '("<DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '(1pe25.15)') par(3)
  WRITE (iunout, '("</DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '("<SECOND_DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '(1pe25.15)') par(4)
  WRITE (iunout, '("</SECOND_DERIVATIVE_BULK_MODULUS>")')
  WRITE (iunout, '("<MINIMUM_ENERGY_RY>")')
  WRITE (iunout, '(1pe25.15)') emin
  WRITE (iunout, '("</MINIMUM_ENERGY_RY>")')
  WRITE (iunout, '("<CELL_FACTOR>")')
  WRITE (iunout, '(1pe25.15)') fac
  WRITE (iunout, '("</CELL_FACTOR>")')
  IF (fac /= 0.0_DP) THEN
     a0(1) = (par(1)/fac)**(1d0/3d0)
     a0(2) = (par(1)/fac)**(1d0/3d0) * bohr_radius_angs
     WRITE (iunout, '("<CELL_PARAMETER_AU_A>")')
     WRITE (iunout, '(2(1pe25.15))') a0
     WRITE (iunout, '("</CELL_PARAMETER_AU_A>")')
  ENDIF
  WRITE (iunout,'("</EQUATIONS_PARAMETERS>")')

  WRITE (iunout,'("<FIT_CHECK>")')
  WRITE (iunout,'("<NUMBER_OF_DATA>")')
  WRITE (iunout,'(i8)') npt
  WRITE (iunout,'("</NUMBER_OF_DATA>")')
  WRITE (iunout,'("<VOL_ENE_EFIT_DELTA_P_GIBBS>")')
  WRITE (iunout,'(6(1pe25.15))') alldata(:,:)
  WRITE (iunout,'("</VOL_ENE_EFIT_DELTA_P_GIBBS>")')

  WRITE (iunout,'("</FIT_CHECK>")')
  CLOSE (unit=iunout, status='keep')
  
  RETURN
END SUBROUTINE write_evdata_xml
!
!---------------------------------------------------------------------------
    FUNCTION birch(x,k0,dk0,d2k0)
!---------------------------------------------------------------------------
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP) birch, x, k0,dk0, d2k0
    REAL(DP) c0, c1

    IF(d2k0/=0.d0) THEN
       c0 = (9.d0*k0*d2k0 + 9.d0*dk0**2 - 63.d0*dk0 + 143.d0 )/48.d0
    ELSE
       c0 = 0.0d0
    ENDIF
    c1 = 3.d0*(dk0-4.d0)/8.d0
    birch = 3.d0*k0*( (-0.5d0+  c1-  c0)*x**( -5.d0/3.d0) &
         +( 0.5d0-2.d0*c1+3.0d0*c0)*x**( -7.d0/3.d0) &
         +(       c1-3.d0*c0)*x**( -9.0d0/3d0) &
         +(            c0)*x**(-11.0d0/3d0) )
    RETURN
    END FUNCTION birch
!
!---------------------------------------------------------------------------
    FUNCTION keane(x,k0,dk0,d2k0)
!---------------------------------------------------------------------------
!
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP) keane, x, k0, dk0, d2k0, ddk

    ddk = dk0 + k0*d2k0/dk0
    keane = k0*dk0/ddk**2*( x**(-ddk) - 1d0 ) + (dk0-ddk)/ddk*log(x)

    RETURN
    END FUNCTION keane

!
END MODULE ev_xml

