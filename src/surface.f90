!
! Copyright (C) 2016 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE convert_rap_surface(nbnd, nks, nkz, high_symmetry, gcodek, aux_ind, &
                               gcodek_ext, ptypek, rap, gaugek)
!
!   In a projected band structure calculation, we change the representations
!   of the points groups of higher symmetry that might occur for particular 
!   values of k_z and bring the representations to those of the subgroup
!   of the surface. It is assumed that the second k_z if nzk is odd, or the 
!   first k_z if nkz is even has no more symmetry than the surface
!     

  USE kinds, ONLY : DP
  USE thermo_sym,    ONLY : code_group_save
  USE control_2d_bands, ONLY : aux_ind_sur
  USE point_group,   ONLY : has_sigma_h, convert_rap_new

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nks, nkz, nbnd
  INTEGER, INTENT(INOUT) ::  gcodek(nks), aux_ind(nks), gcodek_ext(nks), &
                             ptypek(3,nks), rap(nbnd, nks)

  LOGICAL, INTENT(INOUT) :: high_symmetry(nks)
  REAL(DP), INTENT(INOUT) :: gaugek(48,nks)

  LOGICAL :: type1
  INTEGER :: ishift, nks_, ikz, ike, ik2, ik
  INTEGER, ALLOCATABLE :: rapin(:)

  ALLOCATE (rapin(nbnd))
  nks_= nks / nkz
  type1=has_sigma_h(code_group_save)
  IF (type1) THEN
     ishift=nks_
  ELSE
     ishift=0
  ENDIF
  DO ikz=1,nkz
     DO ik = 1, nks_
        ike = ik + nks_ * ( ikz - 1 )
        ik2 = ik + ishift
        IF (gcodek(ike) /= gcodek(ik2)) THEN
           rapin(:)=rap(:,ike)
           CALL convert_rap_new(nbnd,rapin,rap(1,ike),&
                     gcodek_ext(ike),&
                     gcodek_ext(ik2), aux_ind_sur(ik,ikz),&
                     ptypek(1,ike),ptypek(1,ik2),&
                     gaugek(1,ike),gaugek(1,ik2))
           gaugek(:,ike)=gaugek(:,ik2)
           ptypek(:,ike)=ptypek(:,ik2)
           gcodek(ike)=gcodek(ik2)
           aux_ind(ike) = aux_ind(ik2)

!   a point must be high symmetry in all planes.
!
           high_symmetry(ike)=high_symmetry(ik2)
        ENDIF
     END DO
  END DO
  
  DEALLOCATE(rapin)
  RETURN
  END SUBROUTINE convert_rap_surface


SUBROUTINE identify_surface_states(nat, nbnd, nkstot, e, rap)
!
!  This routine searches, among the bands the surface states using the
!  information given in input:
!  sur_layers is the number of surface atoms,
!  sur_thr is the percentage (from 0 to 1) of charge that a state must
!          have on the surface atoms to be considered as a surface state
!  averag  is the projection on each layer of each state
!
USE kinds,            ONLY : DP
USE control_2d_bands, ONLY : averag, vacuum, nlayers, sur_thr, sur_layers, &
                             sym_divide, surface1, surface2, lsurface_state, &
                             subtract_vacuum
USE io_global,        ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: nat, nbnd, nkstot
REAL(DP),INTENT(IN) :: e(nbnd,nkstot)
INTEGER, INTENT(IN) :: rap(nbnd,nkstot)

REAL(DP) :: suggested_sur_thr, maxsum
REAL(DP), ALLOCATABLE :: sumna(:,:)
LOGICAL, ALLOCATABLE :: plot(:)
INTEGER :: na, ibnd, ik, iter, ilayers, npoints, surface_layer

CALL read_state_densities()
ALLOCATE(lsurface_state(nbnd,nkstot))
ALLOCATE(sumna(nbnd,nkstot))
ALLOCATE(plot(nlayers))
!
!  surface1 and surface2 are the surface layers
!  
!
! assume always two equal surfaces
!
IF (sur_layers * 2 > nat ) &
   CALL errore('identify_surface_states','too many surface layers',1)
plot=.FALSE.
DO na=1, sur_layers
   plot(surface1-na+1)=.TRUE.
   plot(surface2+na-1)=.TRUE.
ENDDO
WRITE(stdout,'(/,5x, "Identifing surface states using charge &
                                &density per layer" )') 
WRITE(stdout,'(5x, "with nlayers=",i5," layers per surface",/)') &
         sur_layers
DO ilayers=1,nlayers
   IF (plot(ilayers)) WRITE(stdout,'(5x, "Surface layer", i5)') ilayers
ENDDO

WRITE(stdout,'(/,5x,"Layers are sets of FFT mesh planes perpendicular to the z")')
WRITE(stdout,'(5x, "direction. The first layer contains the origin,")')
WRITE(stdout,'(5x, "the other layers continue with positive z up to &
                   &alat*celldm(3)")')
                               
!
!  plot the bands that have more than plot_thr percentage of the charge
!  on the selected atoms
!
sumna=0.0_DP
maxsum=-1.D20
DO ik=1,nkstot
   DO ibnd=1,nbnd
      DO ilayers=1,nlayers
         IF (plot(ilayers)) sumna(ibnd,ik)=sumna(ibnd,ik)+averag(ilayers,1,ibnd, ik)
      ENDDO
      IF (subtract_vacuum) sumna(ibnd,ik)=sumna(ibnd,ik)-vacuum(1,ibnd,ik)
      IF (sumna(ibnd,ik)>maxsum) maxsum=sumna(ibnd,ik)
   ENDDO
ENDDO

DO iter=1, 6
   suggested_sur_thr = MIN(0.6_DP, MAX(0.2_DP, maxsum-0.15_DP))
   npoints = 0
   DO ik=1,nkstot
      DO ibnd=1,nbnd
         IF (sumna(ibnd,ik)> suggested_sur_thr) npoints=npoints+1
      ENDDO
   ENDDO
   IF (npoints > 10) EXIT
   IF (suggested_sur_thr == 0.2_DP) EXIT
END DO
!
! the suggested sur_thr is not larger than 0.6 and not smaller than 0.2
!
WRITE(stdout, '(/,5x,"Maximum density on the chosen layers", f15.3)') maxsum 
WRITE(stdout, '(5x,"Number of layers", i5)') nlayers

IF (sur_thr == 0.0_DP) THEN
   sur_thr=suggested_sur_thr
   WRITE(stdout,'(5x,"Using sur_thr =",f15.3)') sur_thr
ELSE
   WRITE(stdout,'(5x,"Suggested sur_thr =",f15.3)') suggested_sur_thr
   WRITE(stdout,'(5x,"Using sur_thr =",f15.3,/)') sur_thr
ENDIF
WRITE(stdout,'(25x,30("-"),/)') 

lsurface_state=.FALSE.
WRITE(stdout,'(5x,"Searching surface states for ",i6," k-points and ",&
                            &i6, " bands:")') nkstot, nbnd
IF (subtract_vacuum) &
WRITE(stdout,'(5x,"Vacuum charge has been subtracted")')
WRITE(stdout,'(5x,"ik,    ibnd,   charge on surface layers vacuum charge &
                                                 & energy rap")')
DO ik=1,nkstot
   DO ibnd=1,nbnd
      IF (sumna(ibnd,ik)> sur_thr) THEN
         WRITE(stdout,'(5x,2i8,3f13.7,i5)') ik, ibnd, sumna(ibnd,ik), &
                                              vacuum(1,ibnd,ik), e(ibnd,ik), &
                                              rap(ibnd,ik)
         lsurface_state(ibnd,ik)=.TRUE.
      ENDIF
   ENDDO
ENDDO

DEALLOCATE(plot)
DEALLOCATE(sumna)

RETURN
END SUBROUTINE identify_surface_states

SUBROUTINE plot_surface_states(nbnd, nks, nlines, kx, e_rap, emin, emax, eref, &
                  nrap, nbnd_rapk, start_rapk, start_point, last_point, &
                  nrap_plot, rap_plot )
!
!  This routine writes on the gnuplot scripts the command to
!  plot the surface states. The surface states must have been already
!  identified along each line 
!
USE kinds,            ONLY : DP
USE control_2d_bands, ONLY : sym_divide, lsurface_state_rap
USE control_bands,    ONLY : lsym
USE gnuplot,          ONLY : gnuplot_line, gnuplot_circle, gnuplot_write_command

IMPLICIT NONE
INTEGER, INTENT(IN)  :: nbnd, nks, nlines
INTEGER, INTENT(IN)  :: nrap(nlines), nbnd_rapk(12,nks), start_rapk(12,nks)
INTEGER, INTENT(IN)  :: start_point(nlines), last_point(nlines)
INTEGER, INTENT(IN)  :: nrap_plot(nks), rap_plot(12,nks)
REAL(DP), INTENT(IN) :: kx(nks), e_rap(nbnd, nks)
REAL(DP), INTENT(IN) :: emin, emax, eref
REAL(DP) :: x(2), y(2), ys(2), delta
INTEGER :: ilines, ibnd, ibnd1, jbnd, ik, irap, ir, nrapp
LOGICAL :: dorap

CALL gnuplot_write_command('surface_state_width=2',.FALSE.)
CALL gnuplot_write_command('surface_state_color="blue"',.FALSE.)
CALL gnuplot_write_command('surface_state_radius=0.002*xscale',.FALSE.)

DO ilines=1, nlines
   DO ik=start_point(ilines),last_point(ilines)
      IF (lsym) THEN
         DO irap=1,nrap(ilines)
            dorap=.TRUE.
            IF (sym_divide) THEN
               dorap=.FALSE.
               nrapp= nrap_plot(ik)
               DO ir=1,nrapp
                  dorap=dorap.OR.(rap_plot(ir,ik)==irap)
               END DO
               IF (nrapp==0) dorap=.TRUE.
            END IF
            IF (dorap) THEN
               DO jbnd=1, nbnd_rapk(irap,ik)
                  ibnd=start_rapk(irap,ik)+jbnd-1
                  IF (lsurface_state_rap(ibnd,ik)) THEN
                     x(1)=kx(ik)
                     ys(1)=e_rap(ibnd,ik) - eref
                     y(1)=MIN(emax, ys(1))
                     y(1)=MAX(emin, y(1))
                     IF (.NOT.(y(1)==emax .OR. y(1)==emin) ) THEN
                        CALL gnuplot_circle(x(1),y(1),'surface_state_radius',&
                                       '2', &
                                       'surface_state_color')
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO ibnd=1, nbnd
            IF (lsurface_state_rap(ibnd,ik)) THEN
               x(1)=kx(ik)
               ys(1)=e_rap(ibnd,ik) - eref
               y(1)=MIN(emax, ys(1))
               y(1)=MAX(emin, y(1))
               IF (.NOT.(y(1)==emax .OR. y(1)==emin) ) THEN
                   CALL gnuplot_circle(x(1),y(1),'surface_state_radius',&
                                       '2', &
                                       'surface_state_color')
               ENDIF
            ENDIF
         ENDDO
      END IF
   ENDDO
ENDDO

RETURN
END SUBROUTINE plot_surface_states
