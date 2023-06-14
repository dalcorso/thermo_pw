!
! Copyright (C) 2023 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#if defined(__CUDA)
!-------------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE compute_deff_dev( nhm, nat, current_spin, okvan, deff, et )
  !-----------------------------------------------------------------------
  !! This routine computes the effective value of the D-eS coefficients
  !! which appear often in many expressions in the US or PAW case. 
  !! This routine is for the collinear case.
  !
  USE cudafor
  USE kinds,       ONLY: DP
  USE many_k_mod,  ONLY: qq_at => qq_at_d, deeq =>deeq_d
                         
  !
  IMPLICIT NONE
  !
  INTEGER, VALUE :: current_spin, nhm, nat
  LOGICAL, VALUE :: okvan
  !
  REAL(DP), INTENT(IN), VALUE :: et
  !! The eigenvalues of the hamiltonian
  REAL(DP), INTENT(OUT), DEVICE :: deff(nhm,nhm,nat)
  !! Effective values of the D-eS coefficients
  !
  IF (.NOT. okvan) THEN
     !
     deff(1:nhm,1:nhm,1:nat) = deeq(1:nhm,1:nhm,1:nat,current_spin)
     !
  ELSE
     !
     deff(1:nhm,1:nhm,1:nat) = deeq(1:nhm,1:nhm,1:nat,current_spin) - &
                                          et*qq_at(1:nhm,1:nhm,1:nat)
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE compute_deff_dev
!-------------------------------------------------------------------------
ATTRIBUTES(DEVICE) SUBROUTINE compute_deff_nc_dev( nhm, nat, nspin, isolv, &
                                         npol, okvan, nsolv, deff_nc, et )
  !-----------------------------------------------------------------------
  !! This routine computes the effective value of the D-eS coefficients
  !! which appear often in many expressions in the US or PAW case. 
  !! This routine is for the noncollinear case.
  !
  USE cudafor
  USE kinds,       ONLY: DP
  USE many_k_mod,  ONLY: qq_at => qq_at_d, &
                         qq_so => qq_so_d, lspinorb => lspinorb_d, &
                         ntyp => ntyp_d, ityp => ityp_d
  USE many_k_ph_mod, ONLY : deeq_nc_s =>deeq_nc_save_d, deeq_nc => deeq_nc_d 
                         
  IMPLICIT NONE
  !
  INTEGER, VALUE :: nspin, nhm, nat, npol, isolv, nsolv
  LOGICAL, VALUE :: okvan
  !
  REAL(DP), INTENT(IN), VALUE :: et
  !! The eigenvalues of the hamiltonian
  COMPLEX(DP), INTENT(OUT), DEVICE :: deff_nc(nhm,nhm,nat,nspin)
  !! Effective values of the D-eS coefficients

  INTEGER :: i, j, nt, ia, is, ijs
  !
  IF (nsolv==2) THEN
     deff_nc(:,:,:,:) = deeq_nc_s(:,:,:,:,isolv)
  ELSE
     deff_nc(:,:,:,:) = deeq_nc(:,:,:,:)
  ENDIF
  IF (okvan) THEN
     !
     IF (lspinorb) THEN
        DO nt=1, ntyp
           DO ia = 1, nat
              IF (ityp(ia)==nt) THEN
                 DO i = 1, nhm
                    DO j = 1, nhm
                       IF (nsolv==2) THEN
                          deff_nc(i,j,ia,1:nspin) = &
                                   deeq_nc_s(i,j,ia,1:nspin,isolv) - &
                                   et*qq_so(i,j,1:nspin,nt)
                       ELSE
                          deff_nc(i,j,ia,1:nspin) = &
                                   deeq_nc(i,j,ia,1:nspin) - &
                                   et*qq_so(i,j,1:nspin,nt)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ELSE
        DO nt=1, ntyp
           DO ia = 1, nat
              IF (ityp(ia)==nt) THEN
                 DO i = 1, nhm
                    DO j = 1, nhm
                       DO is = 1, npol
                          ijs = (is-1)*npol + is
                          IF (nsolv==2) THEN
                             deff_nc(i,j,ia,ijs) = deeq_nc_s(i,j,ia,ijs,isolv)-&
                                                 et*qq_at(i,j,ia)
                          ELSE
                             deff_nc(i,j,ia,ijs) = deeq_nc(i,j,ia,ijs) - &
                                                 et*qq_at(i,j,ia)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDIF
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE compute_deff_nc_dev
#endif
!
