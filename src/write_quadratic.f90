!
! Copyright (C) 2015 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE write_quadratic()
!---------------------------------------------------------------------
!
! This routine writes the energy as a function of the celldm(1) for
! the polynomial fit of cubic systems.
! In input emin in Ry, omega0 in (a.u.)**3, b0in in kbar,
! in output omega in (a.u.)**3, p in kbar, e in Ry
!
USE kinds,          ONLY : DP
USE data_files,     ONLY : flevdat
USE thermo_mod,     ONLY : ngeo, celldm_geo_eos
USE cell_base,      ONLY : ibrav
USE control_vol,    ONLY : nvol
USE control_quadratic_energy, ONLY : p2
USE control_quartic_energy, ONLY : p4, lquartic
USE control_pressure,   ONLY : pressure_kb
USE quadratic_surfaces, ONLY : evaluate_fit_quadratic, evaluate_quadratic_grad
USE quartic_surfaces, ONLY : evaluate_fit_quartic, evaluate_quartic_grad
USE lattices,       ONLY : crystal_parameters
USE mp_images,      ONLY : root_image, my_image_id
USE constants,      ONLY : ry_kbar
USE io_global,      ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename
REAL(DP) :: a(1), e, p(1), fact, xmax, xmin, deltaa, x2max(2), x2min(2), &
                 x2(2), delta2(2)
INTEGER :: i, j, iu_mur, iwork, nwork, nvar
INTEGER :: compute_nwork, find_free_unit

IF (my_image_id /= root_image) RETURN

IF (ibrav<1.OR.ibrav>7.OR.ibrav==5) RETURN
IF (ibrav==1) fact=1.0_DP
IF (ibrav==2) fact=0.25_DP
IF (ibrav==3) fact=0.5_DP

nwork=compute_nwork()
nvar=crystal_parameters(ibrav)

filename='energy_files/'//TRIM(flevdat)//'_quadratic'
CALL add_pressure(filename)

IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
   xmax=0.0_DP
   xmin=100000.0_DP
   DO iwork = 1, ngeo(1)
      IF (celldm_geo_eos(1,iwork) > xmax) xmax=celldm_geo_eos(1,iwork)
      IF (celldm_geo_eos(1,iwork) < xmin) xmin=celldm_geo_eos(1,iwork)
   ENDDO
   xmin=xmin*0.99_DP
   xmax=xmax*1.01_DP
   IF (nvol <= 1) nvol=50
   deltaa = (xmax - xmin)/(nvol-1)
ELSEIF (ibrav==4.OR.ibrav==5.OR.ibrav==6.OR.ibrav==7) THEN
   x2max=0.0_DP
   x2min=100000.0_DP
   DO iwork = 1, nwork
      IF (celldm_geo_eos(1,iwork) > x2max(1)) x2max(1)=celldm_geo_eos(1,iwork)
      IF (celldm_geo_eos(1,iwork) < x2min(1)) x2min(1)=celldm_geo_eos(1,iwork)
      IF (ibrav==5) THEN
         IF (celldm_geo_eos(4,iwork) > x2max(2)) x2max(2)=ACOS(celldm_geo_eos(4,iwork))
         IF (celldm_geo_eos(4,iwork) < x2min(2)) x2min(2)=ACOS(celldm_geo_eos(4,iwork))
      ELSE
         IF (celldm_geo_eos(3,iwork) > x2max(2)) x2max(2)=celldm_geo_eos(3,iwork)
         IF (celldm_geo_eos(3,iwork) < x2min(2)) x2min(2)=celldm_geo_eos(3,iwork)
      ENDIF
   ENDDO
   IF (nvol <= 1) nvol=25
   delta2(1) = (x2max(1) - x2min(1))/(nvol-1)
   delta2(2) = (x2max(2) - x2min(2))/(nvol-1)
ENDIF

IF (ionode) THEN
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (ibrav==1.OR.ibrav==2.OR.ibrav==3) THEN
      IF (pressure_kb /= 0.0_DP) THEN
         WRITE(iu_mur,'( "# ",7x,"a (a.u.)",10x,"enthalpy (Ry)",6x,&
                                                       &"pressure (kbar)" )')
      ELSE
         WRITE(iu_mur,'( "# ",7x,"a (a.u.)",12x,"energy (Ry)",6x,&
                                                       &"pressure (kbar)" )')
      END IF

      DO i=1,nvol
         a= xmin + deltaa * (i-1)
         IF (lquartic) THEN
            CALL evaluate_fit_quartic(1,a,e,p4)
            CALL evaluate_quartic_grad(1,a,p,p4)
         ELSE
            CALL evaluate_fit_quadratic(1,a,e,p2)
            CALL evaluate_quadratic_grad(1,a,p,p2)
         ENDIF
         p= - p * ry_kbar / (fact*3.0_DP * a(1)**2)
         WRITE(iu_mur,'(3f20.10)') a, e, p + pressure_kb
      ENDDO 
   ELSEIF (ibrav==4.OR.ibrav==5.OR.ibrav==6.OR.ibrav==7) THEN
      DO i=1,nvol
         x2(1)=x2min(1)+delta2(1)*(i-1)
         DO j=1,nvol
            x2(2)=x2min(2)+delta2(2)*(j-1)
            IF (ibrav==5) x2(2)=COS(x2(2))
            IF (lquartic) THEN
               CALL evaluate_fit_quartic(2,x2,e,p4)
            ELSE
               CALL evaluate_fit_quadratic(2,x2,e,p2)
            ENDIF
            WRITE(iu_mur,'(3f20.10)') x2(1), x2(2), e
         END DO
      END DO  
   ENDIF
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
END IF

RETURN
END SUBROUTINE write_quadratic
