!
! Copyright (C) 1997 Fabio Favot 
! Copyright (C) 2012-2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------
 SUBROUTINE proj_band_structure(kx, e, nks_, nbnd_, emin, emax, eref_, &
                 e_rap, nrap, nbnd_rapk, start_rapk, nlines_, start_point_, &
                 last_point_, nrap_plot, rap_plot )
!------------------------------------------------------------
!
!     This routine uses gnuplot library to plot 
!     the bulk Projected Band Structure  
!     using the method of the superlattice. The bands of
!     a periodic bulk slab taken in the direction of the surface normal  
!     are used to identify the projected bulk bands.
!
!     To calculate the projected bands several paths are given in 
!     input. Each path contains all the k_\\ to the surface and
!     a different k_z (z is the normal to the surface). There are nkz
!     such path.
!     Each one of these paths must have the same number of k_\\ and
!     the same number of bands. If available the file with the 
!     representation number (produced by bands.x with the option 
!     lsym=.true.) of each band can be used to plot separately
!     projected band structure of different symmetry.
!
!     This subroutine adds to the a gnuplot script produced by the
!     plotband_sub subroutine that must have already initialized the
!     plot and written the bands labels and the framework the commands
!     that plot the projected band structure.
!
!     History:
!
!     Fabio Favot  18-02-97
!     created a stand alone version of this program interfaced with the pgplot
!     library. 
!
!     ADC (Feb. 2006) Added the possibility to read the projection file 
!     to select surface states.
!
!     R. Mazzarello used this program to plot the relativistic band
!     structure of Au(111) published in Surf. Sci. 602, 893 (2008) and
!     added some support for reading the representation files. 
!     This part has been lost.
!
!     ADC (Nov. 2012) removed the pgplot part and substituted it with
!     gnuplot commands.  
!
!     ADC (Aug. 2014)  A large part of the original program contained
!     functionalities that are already available in the thermo_pw code.
!     So only this routine of the original program  has been reutilized 
!     and it has been interfaced with the thermo_pw program.
!     The input now is given through the thermo_pw input and the gnuplot
!     initialization as well as the plot of the slab bands is 
!     done by other routines. This routine writes only the commands
!     to plot the projected band structure in the gnuplot script.
!     Support for plotting in different panels the different representations
!     has been reintroduced. The reading of the representation files 
!     is done by other thermo_pw routines.
!
!
USE kinds, ONLY : DP
USE control_2d_bands, ONLY : nkz, gap_thr, sym_divide
USE data_files, ONLY : flpbs
USE gnuplot,  ONLY : gnuplot_polygon, gnuplot_line, gnuplot_write_command
USE point_group, ONLY : convert_rap
USE spin_orb, ONLY : lspinorb
USE io_global, ONLY : ionode, ionode_id, stdout
USE mp_images, ONLY : intra_image_comm
USE mp, ONLY : mp_bcast

IMPLICIT NONE 

INTEGER, INTENT(IN) :: nks_, nbnd_, nlines_
INTEGER, INTENT(IN) :: nrap(nlines_)
INTEGER, INTENT(IN) :: start_point_(nlines_), last_point_(nlines_), &
                       nbnd_rapk(12,nks_), start_rapk(12,nks_)
INTEGER, INTENT(IN) :: nrap_plot(nks_), rap_plot(12,nks_)
REAL(DP), INTENT(IN) :: emin, emax, eref_
REAL(DP) :: e(nbnd_,nks_), kx(nks_), e_rap(nbnd_, nks_)

INTEGER :: nks   ! this is the number of k point of a single path
INTEGER :: ik, ibnd, i, ike, ikz, ik2, it, nbnd, irap, spe, lpe, nbnd_ilines, iun
INTEGER :: ilines, nlines, jbnd, jrap, code_group_line, ios, nrapp, rapp(12)
INTEGER :: find_free_unit
REAL(DP), ALLOCATABLE :: et1(:,:), et2(:,:), eth(:,:)
REAL(DP), ALLOCATABLE :: x(:), y(:)
INTEGER, ALLOCATABLE :: nbnd_plot(:), start_point(:), last_point(:)
REAL(DP) :: q1, q2, m1, m2, xc, yc 
REAL(DP) :: eref
LOGICAL :: doif, plot

LOGICAL, ALLOCATABLE :: lth(:,:)

nks = nks_ / nkz
nbnd = nbnd_
eref = eref_
nlines = nlines_ / nkz
!
!  First convert all the layers 
!
ALLOCATE(start_point(nlines_))
ALLOCATE(last_point(nlines_))
start_point=start_point_
last_point=last_point_

IF ( nkz > 1) THEN
!
   ALLOCATE(et1(nbnd,nks))
   ALLOCATE(et2(nbnd,nks))
   ALLOCATE(nbnd_plot(nks_))
   et1=1.0D+30
   et2=-1.0D+30
   nbnd_plot=0
   DO ilines=1, nlines
! et1,et2 contain the maximum and minimum (along k_z) 
! for each state (ibnd,ik)

!
!  in this run we compute the PBS and plot it
!
      DO ikz = 1, nkz
         DO ik = start_point(ilines), last_point(ilines)
            ike = ik + nks * (ikz - 1) 
            nrapp=nrap_plot(ike)
            nbnd_plot(ike)=0
            IF (nrapp>0) THEN
               rapp(1:nrapp)=rap_plot(1:nrapp,ike)
               DO jrap=1,nrapp 
                  DO irap=1,nrap(ilines)
                     DO jbnd=1,nbnd_rapk(irap,ike)
                        ibnd=start_rapk(irap,ike)+jbnd-1
                        IF (rapp(jrap)==irap.OR..NOT.sym_divide) THEN 
                           nbnd_plot(ike)=nbnd_plot(ike)+1
                           et1(nbnd_plot(ike),ik)=MIN(et1(nbnd_plot(ike),ik), &
                                                     e_rap(ibnd,ike)-eref)
                           et2(nbnd_plot(ike),ik)=MAX(et2(nbnd_plot(ike),ik), &
                                                     e_rap(ibnd,ike)-eref)
                        END IF
                     END DO
                  END DO
               ENDDO
            ELSE
               DO ibnd=1,nbnd
                  nbnd_plot(ike)=nbnd_plot(ike)+1
                  et1(nbnd_plot(ike),ik)=MIN(et1(nbnd_plot(ike),ik), &
                                            e(ibnd,ike)-eref)
                  et2(nbnd_plot(ike),ik)=MAX(et2(nbnd_plot(ike),ik), &
                                            e(ibnd,ike)-eref)
               END DO
            END IF
         END DO
      END DO
   END DO
!
!   For each k use the minimum number of bands available for all k_z
!
   DO ik=1,nks
      DO ikz=1,nkz
         ike=ik + nks * ( ikz - 1 )
         nbnd_plot(ik)=MIN(nbnd_plot(ik), nbnd_plot(ike))
      ENDDO
   ENDDO
!
!  Put the same number of bands along each line
!
   DO ilines=1,nlines
      nbnd_ilines=10000
      DO ik=start_point(ilines),last_point(ilines)
         nbnd_ilines=MIN(nbnd_plot(ik), nbnd_ilines)  
      ENDDO
      DO ik=start_point(ilines),last_point(ilines)
         nbnd_plot(ik)=nbnd_ilines
      ENDDO
   ENDDO

   IF (ionode) THEN
      iun=find_free_unit()
      OPEN(UNIT=iun,FILE=TRIM(flpbs),STATUS='unknown',ERR=100,IOSTAT=ios)
      WRITE(iun, '(3i5,f12.6)') nbnd, nks, nlines, eref
      DO ilines=1, nlines
         WRITE(iun,'(2i5)') start_point(ilines), last_point(ilines)
         DO ik=start_point(ilines), last_point(ilines)
            WRITE(iun, '(i5,f12.7)', ERR=100, IOSTAT=ios) nbnd_plot(ik),&
                                                            kx(ik)
            WRITE(iun, '(8f9.4)', ERR=100, IOSTAT=ios) (et1(ibnd,ik),&
                                 ibnd=1,nbnd_plot(ik))
            WRITE(iun, '(8f9.4)', ERR=100, IOSTAT=ios) (et2(ibnd,ik),&
                                 ibnd=1,nbnd_plot(ik))
         END DO
      END DO
      CLOSE(iun)
   ENDIF
100 CALL mp_bcast(ios,ionode_id,intra_image_comm)
    CALL errore('proj_band_structure','problem writing PBS on file',ABS(ios))
ELSE 
!
!  In this run the PBS information is read from file
!
   IF (ionode) THEN
      OPEN(UNIT=iun, FILE=TRIM(flpbs), STATUS='old', ERR=200, IOSTAT=ios)     
      READ(iun, '(3i5,f12.6)') nbnd, nks, nlines, eref
   ENDIF
200 CALL mp_bcast(ios,ionode_id,intra_image_comm)
   IF (ios /=0) CALL errore('proj_band_structure','opening file',ABS(ios))
   CALL mp_bcast(nbnd,ionode_id,intra_image_comm)
   CALL mp_bcast(nks,ionode_id,intra_image_comm)
   CALL mp_bcast(nlines,ionode_id,intra_image_comm)
   CALL mp_bcast(eref,ionode_id,intra_image_comm)
   ALLOCATE(et1(nbnd,nks))
   ALLOCATE(et2(nbnd,nks))
   ALLOCATE(nbnd_plot(nks))
   IF (ionode) THEN
      DO ilines=1,nlines
         READ(iun,'(2i5)') start_point(ilines), last_point(ilines)
         DO ik=start_point(ilines), last_point(ilines)
            READ(iun,'(i5,f12.7)',ERR=200,END=200,IOSTAT=ios) &
                           nbnd_plot(ik), kx(ik)
            READ(iun,'(8f9.4)',ERR=200,END=200,IOSTAT=ios) (et1(ibnd,ik),&
                                  ibnd=1, nbnd_plot(ik))
            READ(iun,'(8f9.4)',ERR=200,END=200,IOSTAT=ios) (et2(ibnd,ik),&
                                  ibnd=1,nbnd_plot(ik))
         END DO
      END DO
      CLOSE(iun)
   ENDIF
300 CALL mp_bcast(ios,ionode_id,intra_image_comm)
   IF (ios /= 0) THEN
      WRITE(stdout,'("Problems reading from file; no PBS plotted")')
      IF (ALLOCATED(et1)) DEALLOCATE(et1)
      IF (ALLOCATED(et2)) DEALLOCATE(et2)
      IF (ALLOCATED(nbnd_plot)) DEALLOCATE(nbnd_plot)
      RETURN
   ELSE
      CALL mp_bcast(start_point,ionode_id,intra_image_comm)
      CALL mp_bcast(last_point,ionode_id,intra_image_comm)
      CALL mp_bcast(nbnd_plot,ionode_id,intra_image_comm)
      CALL mp_bcast(kx,ionode_id,intra_image_comm)
      CALL mp_bcast(et1,ionode_id,intra_image_comm)
      CALL mp_bcast(et2,ionode_id,intra_image_comm)
   ENDIF
ENDIF
!     
! et1 and et2 (min/max/min/max...) is in a single variable eth(2*nbnd,nks)
!
ALLOCATE(eth(2*nbnd,nks))
ALLOCATE(lth(2*nbnd,nks))
ALLOCATE(x(nks))
ALLOCATE(y(nks))


DO ik=1, nks
   i=1
   DO ibnd=1,nbnd_plot(ik)
      eth(i,ik)=et1(ibnd,ik)
      eth(i+1,ik)=et2(ibnd,ik)
      i=i+2
   END DO
END DO
!
! the extremes of the PBS to true
!
lth=.FALSE.
DO ik=1,nks
   lth(1,ik)=.TRUE.
   lth(2*nbnd_plot(ik),ik)=.TRUE.
   DO ibnd=2,2*nbnd_plot(ik)-2,2
      IF ((eth(ibnd+1,ik)-eth(ibnd,ik)) > gap_thr) THEN
         lth(ibnd,ik)=.TRUE.
         lth(ibnd+1,ik)=.TRUE.
      END IF  
   END DO
END DO
!
! a few definition. Change these commands inside the gnu script to change colors,
! do not change here.
!
CALL gnuplot_write_command('pbs_color="yellow"',.FALSE.)
CALL gnuplot_write_command('pbs_opacity=0.3',.FALSE.)
CALL gnuplot_write_command('pbs_border_color="black"',.FALSE.)
CALL gnuplot_write_command('pbs_border_width=2',.FALSE.)
!
! pbs as a filled area
!
DO ilines=1,nlines
   DO ik=start_point(ilines),last_point(ilines)-1
      DO ibnd=1,2*nbnd_plot(ik)-1
         IF (ik==last_point(ilines)-1) THEN
            doif=(.NOT.(MOD(ibnd,2)==0 .AND.(lth(ibnd,ik)&
              .OR.lth(ibnd,ik+1))) &
            .OR. (MOD(ibnd,2)==0 .AND. lth(ibnd,ik) .AND.     &
            ( .NOT.lth(ibnd,ik+1) .AND. .NOT.lth(ibnd,ik-1) )))  
         ELSEIF (ik==1) THEN
            doif=(.NOT.(MOD(ibnd,2)==0 .AND.(lth(ibnd,ik)&
              .OR.lth(ibnd,ik+1))) &
            .OR. (MOD(ibnd,2)==0 .AND. lth(ibnd,ik) .AND.     &
            ( .NOT.lth(ibnd,ik+1) ))  &
            .OR. (MOD(ibnd,2)==0 .AND. lth(ibnd,ik+1).AND.   &
            (.NOT.lth(ibnd,ik).AND..NOT.lth(ibnd,ik+2)) ) ) 
         ELSE
            doif=(.NOT.(MOD(ibnd,2)==0 .AND.(lth(ibnd,ik)&
              .OR.lth(ibnd,ik+1))) &
            .OR. (MOD(ibnd,2)==0 .AND. lth(ibnd,ik) .AND.     &
            ( .NOT.lth(ibnd,ik+1) .AND. .NOT.lth(ibnd,ik-1) ))  &
            .OR. (MOD(ibnd,2)==0 .AND. lth(ibnd,ik+1).AND.   &
            (.NOT.lth(ibnd,ik).AND..NOT.lth(ibnd,ik+2)) ) ) 
         ENDIF
         IF (doif) THEN
            x(1)=kx(ik)
            x(2)=kx(ik+1)
            x(3)=kx(ik+1)
            x(4)=kx(ik)
            y(1)=eth(ibnd,ik)
            y(2)=eth(ibnd,ik+1)
            y(3)=eth(ibnd+1,ik+1)
            y(4)=eth(ibnd+1,ik)
            plot=((eth(ibnd,ik) < emax) .OR. (eth(ibnd,ik+1) < emax) .OR. &
                 (eth(ibnd+1,ik+1) < emax) .OR. (eth(ibnd+1,ik) < emax) ) 
            plot=plot.AND.((eth(ibnd,ik)>emin) .OR. (eth(ibnd,ik+1)>emin) &
                     .OR. (eth(ibnd+1,ik+1)>emin) .OR. (eth(ibnd+1,ik)>emin) ) 
            IF (plot) CALL gnuplot_polygon(4,x,y,'pbs_opacity','pbs_color')
          ENDIF
!         ELSE IF (.NOT.(lth(ibnd,ik)).AND.lth(ibnd,ik+1).AND.&
!                           MOD(ibnd,2)==0.AND.ik/=(last_point(ilines)-1)) THEN
!            it=1
!            m1=(eth(ibnd,ik+2)-eth(ibnd,ik+1))/(kx(ik+2)-kx(ik+1))
!            m2=(eth(ibnd+it,ik+2)-eth(ibnd+it,ik+1))/(kx(ik+2)-kx(ik+1))
!            q1=(kx(ik+2)*eth(ibnd,ik+1)-kx(ik+1)*&
!                               eth(ibnd,ik+2))/(kx(ik+2)-kx(ik+1))
!            q2=(kx(ik+2)*eth(ibnd+it,ik+1)-kx(ik+1)*&
!                               eth(ibnd+it,ik+2)) /(kx(ik+2)-kx(ik+1))
!            xc=(q2-q1)/(m1-m2)
!            yc=m1*xc+q1
!            IF (xc > kx(ik+1).OR.xc < kx(ik)) THEN
!               x(3)=kx(ik)+.001_DP
!               y(3)=.5_DP*(kx(ik)*(m1+m2)+q1+q2)
!            ELSE
!               x(3)=xc
!               y(3)=yc
!            END IF
!            x(1)=kx(ik)
!            x(2)=kx(ik+1)
!            x(4)=kx(ik+1)
!            x(5)=kx(ik)
!            y(1)=eth(ibnd,ik)
!            y(2)=eth(ibnd,ik+1)
!            y(4)=eth(ibnd+1,ik+1)
!            y(5)=eth(ibnd+1,ik)
!            CALL gnuplot_polygon(5,x,y,'pbs_opacity','pbs_color')
!         ELSE IF (lth(ibnd,ik).AND..NOT.lth(ibnd,ik+1).AND. &
!                     MOD(ibnd,2)==0) THEN
!            it=1
!            m1=(eth(ibnd,ik)-eth(ibnd,ik-1))/(kx(ik)-kx(ik-1))             
!            m2=(eth(ibnd+it,ik)-eth(ibnd+it,ik-1))/(kx(ik)-kx(ik-1))
!            q1=(kx(ik)*eth(ibnd,ik-1)-kx(ik-1)*&
!                             eth(ibnd,ik))/(kx(ik)-kx(ik-1))
!            q2=(kx(ik)*eth(ibnd+it,ik-1)-kx(ik-1)*&
!                              eth(ibnd+it,ik))/(kx(ik)-kx(ik-1))
!            xc=(q2-q1)/(m1-m2)
!            yc=m1*xc+q1
!            IF (xc < kx(ik).OR.xc > kx(ik+1)) THEN
!               x(5)=kx(ik+1)-.001_DP
!               y(5)=.5_DP*(kx(ik+1)*(m1+m2)+q1+q2)
!            ELSE
!               x(5)=xc
!               y(5)=yc
!            END IF
!            x(1)=kx(ik)
!            x(2)=kx(ik+1)
!            x(3)=kx(ik+1)
!            x(4)=kx(ik)
!            y(1)=eth(ibnd,ik)
!            y(2)=eth(ibnd,ik+1)
!            y(3)=eth(ibnd+1,ik+1)
!            y(4)=eth(ibnd+1,ik)
!            CALL gnuplot_polygon(5,x,y,'pbs_opacity','pbs_color')
!         END IF
      END DO
   END DO
END DO

!
! the perimeter of PBS is plotted
!
DO ilines=1,nlines
   DO ik=start_point(ilines)+1,last_point(ilines)-1
      DO ibnd=1,2*nbnd_plot(ik)
         IF (lth(ibnd,ik).AND.lth(ibnd,ik+1)) THEN
            x(1)=kx(ik)
            x(2)=kx(ik+1)
            y(1)=eth(ibnd,ik)
            y(2)=eth(ibnd,ik+1)
            CALL adjust_beyond_limits(x,y,emax,emin,plot)
            IF (plot) CALL gnuplot_line(x,y,'pbs_border_width','front',&
                                                     'pbs_border_color')

            IF (.NOT.lth(ibnd,ik-1)) THEN
               x(2)=kx(ik)
               y(2)=eth(ibnd,ik)
               IF (MOD(ibnd,2)==0) THEN
                  it=1
               ELSE
                  it=-1
               END IF
               m1=(eth(ibnd,ik-1)-eth(ibnd,ik))/(kx(ik-1)-kx(ik))
               m2=(eth(ibnd+it,ik-1)-eth(ibnd+it,ik))/(kx(ik-1)-kx(ik))
               q1=(kx(ik-1)*eth(ibnd,ik)-kx(ik)*eth(ibnd,ik-1))&
                                                     /(kx(ik-1)-kx(ik))
               q2=(kx(ik-1)*eth(ibnd+it,ik)-kx(ik)* &
                            eth(ibnd+it,ik-1))/(kx(ik-1)-kx(ik))
               xc=(q2-q1)/(m1-m2)
               yc=m1*xc+q1
               IF (xc < kx(ik-1)) THEN
                  x(1)=kx(ik-1)+.001_DP
                  y(1)=.5_DP*(kx(ik-1)*(m1+m2)+q1+q2)
               ELSE
                  x(1)=xc
                  y(1)=yc
               END IF
               CALL adjust_beyond_limits(x,y,emax,emin,plot)
               IF (plot) CALL gnuplot_line(x,y,'pbs_border_width','front', &
                                                        'pbs_border_color')
            END IF
         ELSE IF (lth(ibnd,ik).AND.lth(ibnd,ik-1)) THEN
            x(1)=kx(ik)
            y(1)=eth(ibnd,ik)
            IF (MOD(ibnd,2)==0) THEN
               it=1
            ELSE
               it=-1
            END IF
            m1=(eth(ibnd,ik)-eth(ibnd,ik+1)) / (kx(ik)-kx(ik+1))              
            m2=(eth(ibnd+it,ik)-eth(ibnd+it,ik+1)) / (kx(ik)-kx(ik+1)) 
            q1=(kx(ik)*eth(ibnd,ik+1)-kx(ik+1)*eth(ibnd,ik))/(kx(ik)-kx(ik+1))
            q2=(kx(ik)*eth(ibnd+it,ik+1)-kx(ik+1)*&
                                    eth(ibnd+it,ik))/(kx(ik)-kx(ik+1))
            xc=(q2-q1)/(m1-m2)
            yc=m1*xc+q1
            IF (xc > kx(ik+1)) THEN
               x(2)=kx(ik+1)-.001_DP
               y(2)=.5_DP*(kx(ik+1)*(m1+m2)+q1+q2)
            ELSE
               x(2)=xc
               y(2)=yc
            END IF
            CALL adjust_beyond_limits(x,y,emax,emin,plot)
            IF (plot) CALL gnuplot_line(x,y,'pbs_border_width','front', &
                                  'pbs_border_color')
         END IF
      END DO   
   END DO
!
!  first and last point of a line must be treated in a different way
!  let us start with the first point
!
   spe=start_point(ilines)
   DO ibnd=1,2*nbnd_plot(spe)
      IF (lth(ibnd,spe).AND.lth(ibnd,spe+1)) THEN
         x(1)=kx(spe)
         x(2)=kx(spe+1)
         y(1)=eth(ibnd,spe)
         y(2)=eth(ibnd,spe+1)
         CALL adjust_beyond_limits(x,y,emax,emin,plot)

         IF (plot) CALL gnuplot_line(x,y,'pbs_border_width','front', &
                                         'pbs_border_color')

     ELSEIF (lth(ibnd,spe)) THEN
        ik=spe
        x(1)=kx(ik)
        y(1)=eth(ibnd,ik)
        IF (MOD(ibnd,2)==0) THEN
            it=1
        ELSE
            it=-1
        END IF
        m1=(eth(ibnd,ik)-eth(ibnd,ik+1)) / (kx(ik)-kx(ik+1))
        m2=(eth(ibnd+it,ik)-eth(ibnd+it,ik+1)) / (kx(ik)-kx(ik+1))
        q1=(kx(ik)*eth(ibnd,ik+1)-kx(ik+1)*eth(ibnd,ik))/(kx(ik)-kx(ik+1))
        q2=(kx(ik)*eth(ibnd+it,ik+1)-kx(ik+1)*&
                                   eth(ibnd+it,ik))/(kx(ik)-kx(ik+1))
        xc=(q2-q1)/(m1-m2)
        yc=m1*xc+q1
        IF (xc > kx(ik+1)) THEN
           x(2)=kx(ik+1)-.001_DP
           y(2)=.5_DP*(kx(ik+1)*(m1+m2)+q1+q2)
        ELSE
           x(2)=xc
           y(2)=yc
        END IF
        CALL adjust_beyond_limits(x,y,emax,emin,plot)
        IF (plot) CALL gnuplot_line(x,y,'pbs_border_width','front', &
                                        'pbs_border_color')
     ENDIF
  ENDDO
  lpe=last_point(ilines)
  DO ibnd=1,2*nbnd_plot(lpe)
     IF (lth(ibnd,lpe).AND..NOT.lth(ibnd,lpe-1)) THEN
        ik=lpe
        x(1)=kx(ik)
        y(1)=eth(ibnd,ik)
        IF (MOD(ibnd,2)==0) THEN
           it=1
        ELSE
           it=-1
        END IF
        m1=(eth(ibnd,ik)-eth(ibnd,ik-1)) / (kx(ik)-kx(ik-1))
        m2=(eth(ibnd+it,ik)-eth(ibnd+it,ik-1)) / (kx(ik)-kx(ik-1))
        q1=(kx(ik)*eth(ibnd,ik-1)-kx(ik-1)*eth(ibnd,ik))/(kx(ik)-kx(ik-1))
        q2=(kx(ik)*eth(ibnd+it,ik-1)-kx(ik-1)*&
                                    eth(ibnd+it,ik))/(kx(ik)-kx(ik-1))
        xc=(q2-q1)/(m1-m2)
        yc=m1*xc+q1
        IF (xc < kx(ik-1)) THEN
           x(2)=kx(ik-1)+.001_DP
           y(2)=.5_DP*(kx(ik-1)*(m1+m2)+q1+q2)
        ELSE
           x(2)=xc
           y(2)=yc
        END IF
        CALL adjust_beyond_limits(x,y,emax,emin,plot)
        IF (plot) CALL gnuplot_line(x,y,'pbs_border_width','front', &
                                        'pbs_border_color')
     ENDIF
   ENDDO
END DO

DEALLOCATE(et1)
DEALLOCATE(et2)
DEALLOCATE(start_point)
DEALLOCATE(last_point)
DEALLOCATE(nbnd_plot)
DEALLOCATE(eth)
DEALLOCATE(lth)
DEALLOCATE(x)
DEALLOCATE(y)

RETURN
END SUBROUTINE proj_band_structure

SUBROUTINE adjust_beyond_limits(x,y,emax,emin,plot)

USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), INTENT(IN) :: emax, emin
REAL(DP), INTENT(INOUT) :: x(2), y(2)
LOGICAL, INTENT(OUT) :: plot
REAL(DP) :: delta, ys(2)

plot=.TRUE.
ys=y
y(1)=MIN(y(1),emax)
y(2)=MIN(y(2),emax)
y(1)=MAX(y(1),emin)
y(2)=MAX(y(2),emin)
IF (.NOT. (( ABS(y(1)-emax)<1.d-8 .AND. ABS(y(2)-emax)<1.d-8 ) .OR. &
          ( ABS(y(1)-emin)<1.d-8 .AND. ABS(y(2)-emin)<1.d-8 )) ) THEN
   IF (ABS(ys(2) - ys(1))>1.d-5) THEN
      delta=(x(2) - x(1) ) / (ys(2) - ys(1))
      IF (y(1)==emax) x(1)=delta * (emax-ys(1)) + x(1)
      IF (y(2)==emax) x(2)=delta * (emax-ys(1)) + x(1)
      IF (y(1)==emin) x(1)=delta * (emin-ys(1)) + x(1)
      IF (y(2)==emin) x(2)=delta * (emin-ys(1)) + x(1)
   ENDIF
ELSE
   plot=.FALSE.
ENDIF

RETURN
END SUBROUTINE adjust_beyond_limits
