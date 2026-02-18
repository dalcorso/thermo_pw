!
! Copyright (C) 2025 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE polarization_vector
!---------------------------------------------------------------------------
!
!   this module contains the support routines for the calculation
!   of the polarization. Presently it contains only routines to read
!   and write the polarization on file.
!   It contains also a routine that computes the contribution
!   to the pyroelectric tensor due to the piezoelectric tensor
!
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  PRIVATE

  INTEGER :: mod_tot
!
!   Some array to simplify dealing with pyroelectric tensor
!
  INTEGER, PARAMETER :: py_elements=3
  
  CHARACTER(LEN=6) :: py_names(py_elements)

  DATA  py_names / 'p_{1} ', 'p_{2} ', 'p_{3} ' /

  INTEGER, PARAMETER :: py_types=10

  INTEGER :: py_code_group(py_types)  ! code of the point group for each type
  DATA  py_code_group /  1, 3, 4, 5, 6, 7, 12, 13, 14, 15 /

  INTEGER  :: py_present(py_elements, py_types)

  DATA py_present / &
       1,2,3, & ! 1  C_1
       1,0,3, & ! 3  C_s
       0,2,0, & ! 4  C_2
       0,0,3, & ! 5  C_3
       0,0,3, & ! 6  C_4
       0,0,3, & ! 7  C_6
       0,0,3, & ! 12 C_2v
       0,0,3, & ! 13 C_3v
       0,0,3, & ! 14 C_4v
       0,0,3  / ! 15 C_6v


  PUBLIC write_polarization, read_polarization, mod_tot, &
         write_pyro_on_file, piezo_pyro, compute_pyro,   &
         py_names, py_types, py_present, py_code_group,  &
         get_py_type, py_elements
        

CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE write_polarization(filename,polar,berry_phase)
!-------------------------------------------------------------------------
!
!  This routine writes the Berry phase and the polarization on file.
!  It must be called after computing the polarization
!  It saves: 
!  the spontaneous polarization of the current geometry
!
USE io_global, ONLY : ionode, ionode_id
USE constants, ONLY : electron_si, bohr_radius_si
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(IN) :: polar(3), berry_phase(3)
REAL(DP) :: fact
INTEGER :: find_free_unit
INTEGER :: outunit, i, ios

IF (ionode) THEN
   outunit=find_free_unit()
   OPEN(UNIT=outunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('write_polarization','ploblem opening output file', ABS(ios))

fact= electron_si / (bohr_radius_si)**2
IF (ionode) THEN
   WRITE(outunit,'("Berry phase calculated along the three &
                                          &reciprocal vectors")')
   WRITE(outunit,'(3e20.10)') (berry_phase(i), i=1,3)
   WRITE(outunit,*)
   WRITE(outunit,'("Spontaneous polarization (cartesian &
                                       &coordinates) (e/bohr**2)")')
   WRITE(outunit,'(3e20.10)') (polar(i), i=1,3)
   WRITE(outunit,*)
   WRITE(outunit,'("Spontaneous polarization (cartesian &
                                       &coordinates) (C/m**2)")')
   WRITE(outunit,'(3e20.10)') (polar(i)*fact, i=1,3)
   CLOSE(outunit)
ENDIF

RETURN
END SUBROUTINE write_polarization

!-------------------------------------------------------------------------
SUBROUTINE read_polarization(filename,polar,berry_phase)
!-------------------------------------------------------------------------
!
!  This routine writes the Berry phase and the polarization on file.
!  It must be called after computing the polarization
!  It saves: 
!  the spontaneous polarization of the current geometry
!
USE io_global, ONLY : ionode, ionode_id
USE mp_images, ONLY : intra_image_comm
USE mp,        ONLY : mp_bcast 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(DP), INTENT(OUT) :: polar(3), berry_phase(3)
REAL(DP) :: polar_cm2(3)
INTEGER :: find_free_unit
INTEGER :: inunit, i, ios

IF (ionode) THEN
   inunit=find_free_unit()
   OPEN(UNIT=inunit, FILE=TRIM(filename), STATUS='unknown', FORM='formatted', &
        ERR=100, IOSTAT=ios)
ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('read_polarization','ploblem opening input file', ABS(ios))

IF (ionode) THEN
   READ(inunit,*)
   READ(inunit,'(3e20.10)') (berry_phase(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   READ(inunit,'(3e20.10)') (polar(i), i=1,3)
   READ(inunit,*)
   READ(inunit,*)
   READ(inunit,'(3e20.10)') (polar_cm2(i), i=1,3)
   CLOSE(inunit)
ENDIF
CALL mp_bcast(berry_phase,ionode_id,intra_image_comm)
CALL mp_bcast(polar,ionode_id,intra_image_comm)

RETURN
END SUBROUTINE read_polarization

!-------------------------------------------------------------------------
SUBROUTINE compute_pyro(zeu, dtau_duint,  &
          alpha_int, pyro, nat, max_nint_var, nint_var)  
!-------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: nat, max_nint_var, nint_var
REAL(DP), INTENT(IN) :: zeu(3,3,nat), dtau_duint(3,nat,max_nint_var)
REAL(DP), INTENT(IN) :: alpha_int(max_nint_var)
REAL(DP), INTENT(OUT) :: pyro(3)

REAL(DP) :: zeu_int(3,nint_var) 
INTEGER  :: ipol, na, ivar
!
!  Compute the Born effective charges for the internal parameter
!
pyro=0.0_DP
zeu_int=0.0_DP
DO na=1,nat
   DO ipol=1,3
      DO ivar=1,nint_var
         zeu_int(:,ivar)= zeu_int(:,ivar) + zeu(:,ipol,na) * &
                                        dtau_duint(ipol,na,ivar)
      ENDDO
   ENDDO
ENDDO
!
! Computes the pyroelectric coefficient
!
DO ivar=1, nint_var
   pyro(:) = pyro(:) - zeu_int(:, ivar) * alpha_int(ivar)
ENDDO

RETURN
END SUBROUTINE compute_pyro

!-------------------------------------------------------------------------
SUBROUTINE piezo_pyro(e_piezo_tensor, thermal_exp, pyro)
!-------------------------------------------------------------------------
!
!  This routine receives as input the stress piezoelectric tensor 
!  (in e/(a.u.)**2) and the thermal expansion tensor (adimensional)
!  and gives as output the contribution of the piezoelectric effect
!  to pyroelectricity. 
!  NB: the first index of the piezoelectric tensor refer to the polarization,
!  the other two (that in input are condensed in one index in Voigt form) 
!  to the strain.
!
USE kinds, ONLY : DP
USE voigt, ONLY : to_voigt3

IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: e_piezo_tensor(3,6), thermal_exp(3,3)
REAL(DP), INTENT(OUT) :: pyro(3)

INTEGER :: ipol, jpol, kpol, lpol, mpol, npol
REAL(DP) :: fact
REAL(DP) :: e_piezo_aux(3,3,3)

CALL to_voigt3(e_piezo_tensor, e_piezo_aux, 1.0_DP, .FALSE.)
!
! In output the polarization is in e/(a.u)^2 units as the input
! piezoelectric tensor
!
pyro(:)=0.0_DP
DO ipol=1,3
   DO kpol=1,3
      DO lpol=1,3
         pyro(ipol)=pyro(ipol)- e_piezo_aux(ipol,kpol,lpol) * &
                                thermal_exp(kpol,lpol)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE piezo_pyro
!
!-------------------------------------------------------------------------
SUBROUTINE write_pyro_on_file(temp, ntemp, pyro_t, piezo_pyro_t, &
                                       astring, filename, iflag)
!-------------------------------------------------------------------------
!
!  iflag=0 writes the pyroelectricity as a function of temperature
!  iflag=2 writes the pyroelectricity as a function of pressure
!
USE kinds,      ONLY : DP
USE io_global,  ONLY : meta_ionode, meta_ionode_id, stdout
USE mp_world,   ONLY : world_comm
USE mp,         ONLY : mp_bcast
IMPLICIT NONE
INTEGER, INTENT(IN) :: ntemp, iflag
REAL(DP), INTENT(IN) :: temp(ntemp), pyro_t(3,ntemp), piezo_pyro_t(3,ntemp)
CHARACTER(LEN=*), INTENT(IN) :: astring
CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER :: itemp, iu_pyro, ios
INTEGER :: find_free_unit
CHARACTER(LEN=7) :: label

iu_pyro=find_free_unit()
IF (meta_ionode) &
   OPEN(UNIT=iu_pyro, FILE=TRIM(filename), FORM='formatted', &
                                       STATUS='UNKNOWN', ERR=30, IOSTAT=ios)
30 CALL mp_bcast(ios, meta_ionode_id, world_comm)
   CALL errore('write_pyro_on_file', 'opening output file', ABS(ios))
!
!  Choose if to plot as a function of temperature or pressure
!
IF (iflag<2) THEN
   label='  T(K) '
ELSE
   label='p(kbar)'
ENDIF

IF (meta_ionode) THEN

   WRITE(iu_pyro,'(a)') astring
   WRITE(iu_pyro,'("#    multiply by 0.57214766E+02 to have it in C/m^2")')
   WRITE(iu_pyro,'(a)') astring
   WRITE(iu_pyro,'("#",5x, a7, 13x, " p_1 ", 13x, " p_2",          &
                      &13x, " p_3", 13x, " pp_1", 13x, " pp_2",    &
                      &13x, " pp_3")') label

   DO itemp=2,ntemp-1
      WRITE(iu_pyro,'(e16.8,6e20.12)') temp(itemp),                &
           pyro_t(1,itemp), pyro_t(2,itemp), pyro_t(3,itemp),      &
           piezo_pyro_t(1,itemp), piezo_pyro_t(2,itemp), piezo_pyro_t(3,itemp)
   ENDDO
   CLOSE(iu_pyro)
ENDIF

RETURN
END SUBROUTINE write_pyro_on_file
!
!-----------------------------------------------------------------------
FUNCTION get_py_type(code_group)
!-----------------------------------------------------------------------
INTEGER :: get_py_type
INTEGER, INTENT(IN) :: code_group

INTEGER :: itype, aux_type

aux_type=0
DO itype=1,py_types
   IF (py_code_group(itype)==code_group) aux_type=itype
ENDDO
IF (aux_type==0) CALL errore('get_py_type','code_group not available',1)

get_py_type=aux_type
RETURN
END FUNCTION get_py_type

END MODULE polarization_vector
