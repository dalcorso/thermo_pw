!
! Copyright (C) 2016-present Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE eldos_module
!
!  This module provides methods to read electron dos files and to calculate 
!  their contribution to the free energy assuming independent electrons. 
!  It defines a type eldos that contains the electron dos as a function 
!  of energy. For metals it contains also the Fermi level.
!
USE kinds, ONLY : DP
USE constants, ONLY :  k_boltzmann_ry, rytoev
IMPLICIT NONE
SAVE
PRIVATE

REAL(DP), PARAMETER :: kb=k_boltzmann_ry ! Boltzmann constant in Ry/K
REAL(DP), PARAMETER :: kb1=1.0_DP/kb     ! inverse Boltzmann 
                                         ! constant in Ry/K

TYPE eldos_type
   INTEGER :: number_of_points        ! the number of points
   LOGICAL :: lsda                    ! if true dos for spin up and down
   REAL(DP) :: de                     ! interval of the mesh of  (ev)
   REAL(DP), ALLOCATABLE :: e(:)      ! the energies (ev)
   REAL(DP), ALLOCATABLE :: dos(:)    ! the electron dos (states/ ev)
   REAL(DP), ALLOCATABLE :: ddos(:)   ! the electron dos (states/ ev) for lsda
   REAL(DP), ALLOCATABLE :: intdos(:) ! the integrated electron dos (states/ ev)
END TYPE eldos_type

PUBLIC :: eldos_type, read_eldos_data, el_free_energy, el_energy, el_entropy, &
          el_specific_heat_cv, el_chem_pot, set_eldos, destroy_eldos
          
CONTAINS

!--------------------------------------------------------------------
SUBROUTINE set_eldos(eldos,ndiv,lsda,deltae)
!--------------------------------------------------------------------
IMPLICIT NONE
TYPE(eldos_type), INTENT(INOUT) :: eldos
INTEGER, INTENT(IN) :: ndiv
LOGICAL, INTENT(IN) :: lsda
REAL(DP), INTENT(IN) :: deltae

eldos%number_of_points=ndiv
eldos%de=deltae
eldos%lsda=lsda
ALLOCATE(eldos%e(ndiv))
ALLOCATE(eldos%dos(ndiv))
IF (lsda) ALLOCATE(eldos%ddos(ndiv))
ALLOCATE(eldos%intdos(ndiv))

RETURN
END SUBROUTINE set_eldos

!--------------------------------------------------------------------
SUBROUTINE destroy_eldos(eldos)
!--------------------------------------------------------------------
!
IMPLICIT NONE
TYPE(eldos_type), INTENT(INOUT) :: eldos

IF (ALLOCATED(eldos%e)) DEALLOCATE(eldos%e)
IF (ALLOCATED(eldos%dos)) DEALLOCATE(eldos%dos)
IF (ALLOCATED(eldos%ddos)) DEALLOCATE(eldos%ddos)
IF (ALLOCATED(eldos%intdos)) DEALLOCATE(eldos%intdos)

RETURN
END SUBROUTINE destroy_eldos

!--------------------------------------------------------------------
SUBROUTINE read_eldos_data(eldos, lsda, filename)
!--------------------------------------------------------------------
!
!  This subroutine reads the eldos from a file. It allocates space,
!  opens and closes the eldos file.
!
USE constants, ONLY : rytoev
USE mp_images, ONLY : intra_image_comm
USE io_global, ONLY : ionode_id, ionode, stdout
USE mp,        ONLY : mp_bcast

IMPLICIT NONE
TYPE(eldos_type), INTENT(INOUT) :: eldos
CHARACTER(LEN=256), INTENT(IN) :: filename
LOGICAL, INTENT(IN) :: lsda
INTEGER :: iunit, ios
INTEGER, PARAMETER :: ndivx=10000000
REAL(DP), ALLOCATABLE :: e(:), dos(:,:), intdos(:)
REAL(DP) :: de, de_, deltae
INTEGER :: i, ndiv
INTEGER :: find_free_unit

IF (ionode) THEN
   iunit=find_free_unit()
   OPEN(file=TRIM(filename), unit=iunit, status='old', &
                            form='formatted', err=100, iostat=ios)
ENDIF
100 CALL mp_bcast(ios, ionode_id, intra_image_comm)
IF (ios /= 0) CALL errore('read_eldos_data', &
                          'opening file'//TRIM(filename), ABS(ios))

ALLOCATE(e(ndivx))
IF (lsda) THEN
   ALLOCATE(dos(ndivx,2))
ELSE
   ALLOCATE(dos(ndivx,1))
ENDIF
ALLOCATE(intdos(ndivx))
de = 0d0
IF (ionode) THEN
   READ(iunit, *, ERR=10, IOSTAT=ios) 
   DO i=1,ndivx
      IF (lsda) THEN
         READ(iunit, *, END=20, ERR=10, IOSTAT=ios) e(i), dos(i,1), dos(i,2), &
                                                          intdos(i)
      ELSE
         READ(iunit, *, END=20, ERR=10, IOSTAT=ios) e(i), dos(i,1), intdos(i)
      ENDIF
      IF ( i ==2 ) de_ = e(2) - e(1)
      IF (i > 2) THEN
         de = e(i) - e(i-1)
         IF ( ABS(de - de_) > 1.0d-4 ) &
            CALL errore('read_eldos_data','nonuniform grid',1)
      END IF
      ndiv=i
   ENDDO
10 IF (ios /= 0 ) CALL errore('read_eldos_data', 'problem reading eldos', 1)
   IF (ndiv == ndivx ) CALL errore('read_eldos_data', 'increase ndivx', 1)
20 CONTINUE
ENDIF
CALL mp_bcast(ndiv,ionode_id,intra_image_comm)
CALL mp_bcast(de,ionode_id,intra_image_comm)
CALL mp_bcast(e,ionode_id,intra_image_comm)
CALL mp_bcast(dos,ionode_id,intra_image_comm)
CALL mp_bcast(intdos,ionode_id,intra_image_comm)

deltae=de / rytoev   ! internally energies are in Ry
CALL set_eldos(eldos,ndiv,lsda,deltae)

eldos%e(:) = e(1:ndiv) / rytoev            ! internally energies are in Ry
eldos%dos(:) = dos(1:ndiv,1) * rytoev        ! internally dos in states/Ry
IF (lsda) eldos%ddos(:) = dos(1:ndiv,2) * rytoev ! internally dos in states/Ry
eldos%intdos(:) = intdos(1:ndiv)

DEALLOCATE(e)
DEALLOCATE(dos)
DEALLOCATE(intdos)
IF (ionode) CLOSE(iunit)

RETURN
END SUBROUTINE read_eldos_data

!--------------------------------------------------------------------
SUBROUTINE el_entropy(eldos, temp, mu, entropy)
!--------------------------------------------------------------------
!
!  This routine receives as input an eldos, a temperature and a chemical
!  potential and gives as output the electron entropy at that temperature. 
!  temp in Kelvin, mu in Ry, entropy in Ry/Kelvin
!
TYPE(eldos_type), INTENT(IN) :: eldos
REAL(DP), INTENT(IN) :: temp  ! the temperature
REAL(DP), INTENT(IN) :: mu    ! the chemical potential
REAL(DP), INTENT(OUT) :: entropy

INTEGER :: ndiv, i
REAL(DP) :: arg, temp1, earg, gauss, e
REAL(DP) :: wgauss

entropy=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=eldos%number_of_points

DO i=1,ndiv
   e=eldos%e(i)
   arg= kb1 * (mu-e) * temp1
   gauss = wgauss( arg, -99 )
   IF ( (gauss > 0.0_DP) .AND. (gauss < 1.0_DP) ) THEN
      earg = gauss * log( gauss ) + ( 1.0_DP - gauss ) * log( 1.0_DP - gauss )
      entropy = entropy - eldos%dos(i)* earg
      IF (eldos%lsda) entropy = entropy - eldos%ddos(i)* earg
   END IF
END DO

entropy = entropy * eldos%de * kb

RETURN
END SUBROUTINE el_entropy

!--------------------------------------------------------------------
SUBROUTINE el_energy(eldos, temp, mu, ener)
!--------------------------------------------------------------------
!
!  This routine receives as input an eldos and a temperature and gives as 
!  output the indipendent electron energy at that temperature. 
!  The temperature is in kelvin, mu and the energy in Ry.
!
TYPE(eldos_type), INTENT(IN) :: eldos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(IN) :: mu
REAL(DP), INTENT(OUT) :: ener

INTEGER :: ndiv, i
REAL(DP) :: nu, temp1, arg, earg, e, ener0
REAL(DP) :: wgauss

ener=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv=eldos%number_of_points
DO i=1,ndiv
   e=eldos%e(i)
   arg= kb1 * (mu-e) * temp1
   earg = wgauss( arg, -99 )
   ener = ener + eldos%dos(i)* e * earg 
   IF (eldos%lsda) ener = ener + eldos%ddos(i)* e * earg 
ENDDO
ener = ener * eldos%de 

RETURN
END SUBROUTINE el_energy

!--------------------------------------------------------------------
SUBROUTINE el_free_energy(eldos, temp, mu, free_energy)
!--------------------------------------------------------------------
!
!  This routine receives as input an eldos, a temperature and a
!  chemical potential and gives as output the electron free_energy 
!  at that temperature. 
!
TYPE(eldos_type), INTENT(IN) :: eldos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(IN) :: mu
REAL(DP), INTENT(OUT) :: free_energy
REAL(DP) :: ener, entr

CALL el_entropy (eldos, temp, mu, entr)
CALL el_energy (eldos, temp, mu, ener)

IF (temp > 0.0_DP) THEN
   free_energy = ener - temp * entr
ELSE
   free_energy = 0.0_DP
ENDIF

RETURN
END SUBROUTINE el_free_energy

!--------------------------------------------------------------------
SUBROUTINE el_specific_heat_cv(eldos, temp, mu, cv)
!--------------------------------------------------------------------
!
!  This routine receives as input an electron dos, a temperature and
!  a chemical potential and gives as output the constant volume specific 
!  heat at that temperature and chemical potential. 
!  The output cv is in Ry / K.
!
TYPE(eldos_type), INTENT(IN) :: eldos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(IN) :: mu
REAL(DP), INTENT(OUT) :: cv

INTEGER :: ndiv, i
REAL(DP) :: temp1, arg, earg, e
REAL(DP) :: w0gauss

cv=0.0_DP
IF (temp <= 1.E-9_DP) RETURN
temp1 = 1.0_DP / temp
ndiv = eldos%number_of_points
DO i = 1, ndiv
   e = eldos%e(i)
   arg = kb1 * ( mu - e ) * temp1
   earg = w0gauss( arg, -99 )
   cv = cv + eldos%dos(i) * earg * ( e - mu )**2
   IF (eldos%lsda) cv = cv + eldos%ddos(i) * earg * ( e - mu )**2
ENDDO
cv = cv * eldos%de * temp1 ** 2 / kb

RETURN
END SUBROUTINE el_specific_heat_cv

!--------------------------------------------------------------------
SUBROUTINE el_chem_pot(eldos, temp, nelec, mu)
!--------------------------------------------------------------------
!
!  This routine receives as input an electron dos, a temperature, and 
!  the number of electrons and gives as output the chemical potential 
!  at that temperature. The output mu is in Ry.
!
USE io_global, ONLY : stdout
TYPE(eldos_type), INTENT(IN) :: eldos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(IN) :: nelec
REAL(DP), INTENT(OUT) :: mu

INTEGER :: maxter=100
INTEGER :: iter, ndiv
REAL(DP) :: mux, mu1, mu2, f1, f2, fx

ndiv = eldos%number_of_points
mu1 = eldos%e(1)
mu2 = eldos%e(ndiv)
f1 = integrated_dos(eldos, temp, mu1) - nelec
f2 = integrated_dos(eldos, temp, mu2) - nelec


IF ( f1 * f2 > 0.0_DP) THEN
   WRITE(stdout,'(4f20.8)') mu1, mu2, f1, f2
   CALL errore('el_chem_pot','Problem finding mu',1)
ENDIF

DO iter=1, maxter
   mux = mu1 + (mu2 - mu1) * 0.5_DP
   fx = integrated_dos(eldos, temp, mux) - nelec
   IF (ABS(fx) < 1.D-8) THEN
      mu=mux
      EXIT
   END IF
   IF (fx*f1 > 0.0_DP) THEN
      mu1 = mux
      f1=fx
   ELSE
      mu2 = mux
      f2  = fx
   ENDIF
END DO
IF (iter==maxter) CALL errore('el_chem_pot','Difficulty with bisection',1)

RETURN
END SUBROUTINE el_chem_pot

!--------------------------------------------------------------------
FUNCTION integrated_dos (eldos, temp, mu)
!--------------------------------------------------------------------
!
!  This routine receives as input an electron dos, a temperature and
!  a chemical potential and gives as output the number of electrons 
!  that correnspond to that chemical potential at that temperature
!
REAL(DP) :: integrated_dos
TYPE(eldos_type), INTENT(IN) :: eldos
REAL(DP), INTENT(IN) :: temp
REAL(DP), INTENT(IN) :: mu

INTEGER :: ndiv, i
REAL(DP) :: temp1, arg, earg, e, aux
REAL(DP) :: wgauss

temp1 = 1.0_DP / temp
ndiv = eldos%number_of_points
aux=0.0_DP
DO i = 1, ndiv
   e = eldos%e(i)
   arg = kb1 * ( mu - e ) * temp1
   earg = wgauss( arg, -99 )
   aux = aux + eldos%dos(i) * earg 
   IF (eldos%lsda) aux = aux + eldos%ddos(i) * earg 
ENDDO
integrated_dos = aux * eldos%de 

RETURN
END FUNCTION integrated_dos

END MODULE eldos_module
