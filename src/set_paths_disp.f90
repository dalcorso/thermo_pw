!
! Copyright (C) 2013-2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE set_paths_disp()
!---------------------------------------------------------------------
!
!  This subroutine computes the paths for (band or phonon) dispersions 
!  from the label letters or the auxiliary (k or q) points given in
!  input. Note that this routine should be called after the calculation
!  of the geometry. If the latter is changing the input k or q points
!  should be given as letters, or in crystal coordinates, to be computed 
!  here with the correct geometry.
!
  USE kinds,         ONLY : DP
  USE control_paths, ONLY : nqaux, xqaux, wqaux, wqauxr, npk_label, letter, &
                            label_list, nqaux, q_in_band_form, &
                            q_in_cryst_coord, q2d, point_label_type, &
                            label_disp_q, disp_nqs, disp_q, disp_wq, &
                            nrap_plot, rap_plot, nrap_plot_in, rap_plot_in, &
                            dkmod_save, high_sym_path
  USE control_2d_bands, ONLY : nkz, sym_divide
  USE thermo_mod,    ONLY : what
  USE cell_base,     ONLY : ibrav, celldm, bg
  USE bz_form,       ONLY : transform_label_coord
  USE bz_2d_form,    ONLY : transform_2d_label_coord

  IMPLICIT NONE
  CHARACTER(len=80) :: k_points = 'tpiba'
  INTEGER :: n, i, j, ik
  REAL(DP) :: q1(3), q2(3), modq1, modq2, ps
 
  IF (nqaux==0) RETURN

  IF ( ALLOCATED(disp_q) )    DEALLOCATE (disp_q)
  IF ( ALLOCATED(disp_wq) )   DEALLOCATE (disp_wq)
  IF ( ALLOCATED(high_sym_path) )   DEALLOCATE (high_sym_path)
  IF ( ALLOCATED(nrap_plot) ) DEALLOCATE (nrap_plot)
  IF ( ALLOCATED(rap_plot) )  DEALLOCATE (rap_plot)

  IF (q_in_cryst_coord) k_points='crystal'

  IF ( npk_label > 0 ) THEN
     IF (what=='scf_2d_bands') THEN
        CALL transform_2d_label_coord(ibrav, celldm, xqaux, letter, &
                                      label_list, npk_label, nqaux, k_points )
     ELSE
        CALL transform_label_coord(ibrav, celldm, xqaux, letter, &
        label_list, npk_label, nqaux, k_points, point_label_type )
     END IF
  END IF
  IF (q_in_cryst_coord)  CALL cryst_to_cart(nqaux,xqaux,bg,+1)

  IF (q2d) THEN
     disp_nqs=wqaux(2)*wqaux(3)
     ALLOCATE(disp_q(3,disp_nqs))
     ALLOCATE(disp_wq(disp_nqs))
     CALL generate_k_in_plane(nqaux, xqaux, wqaux, disp_q, disp_wq, disp_nqs)
  ELSEIF (q_in_band_form) THEN
      CALL clean_path(nqaux, xqaux, wqaux)
      disp_nqs=SUM(wqaux(1:nqaux-1))+1
      DO i=1,nqaux-1
         IF (wqaux(i)==0) disp_nqs=disp_nqs+1
      ENDDO
      disp_nqs=disp_nqs * nkz
      ALLOCATE(disp_q(3,disp_nqs))
      ALLOCATE(disp_wq(disp_nqs))
      ALLOCATE(nrap_plot(disp_nqs))
      ALLOCATE(rap_plot(12,disp_nqs))
      nrap_plot=0
      rap_plot=0
      IF (nkz > 1.OR.sym_divide) THEN
         CALL generate_k_along_lines_kz(nqaux, xqaux, wqaux, disp_q, & 
                                  disp_wq, disp_nqs, nkz)
      ELSE
         CALL generate_k_along_lines(nqaux, xqaux, wqaux, disp_q, disp_wq, &
                               disp_nqs)
      ENDIF
      label_disp_q=0
      label_disp_q(1)=1
      DO i=2,nqaux
         label_disp_q(i)=label_disp_q(i-1)+wqaux(i-1)
         IF (wqaux(i-1)==0) label_disp_q(i)=label_disp_q(i)+1
      ENDDO
   ELSE
      disp_nqs=nqaux
      ALLOCATE(disp_q(3,disp_nqs))
      ALLOCATE(disp_wq(disp_nqs))
      ALLOCATE(nrap_plot(disp_nqs))
      ALLOCATE(rap_plot(12,disp_nqs))
      nrap_plot=0
      rap_plot=0
      disp_wq(:)=wqauxr(:)
      disp_q(:,1:disp_nqs)=xqaux(:,1:disp_nqs)
      
      DO i=1,nqaux
         label_disp_q(i)=i
      ENDDO
   ENDIF
!
!   For 2d plots we can return here.
!
   IF (q2d) RETURN
!
!  detect if some points are high symmetry points from the change
!  of direction of the path
!
  ALLOCATE(high_sym_path(disp_nqs))
  high_sym_path=.FALSE.
  dkmod_save=0.0_DP
  DO n=1,disp_nqs
     IF (n==1.OR.n==disp_nqs) THEN
!
!     first and last points are assumed to be high symmetry points
!
        high_sym_path(n) = .TRUE.
     ELSE
        q1(:) = disp_q(:,n) - disp_q(:,n-1)
        q2(:) = disp_q(:,n+1) - disp_q(:,n)
        modq1=sqrt( q1(1)*q1(1) + q1(2)*q1(2) + q1(3)*q1(3) )
        modq2=sqrt( q2(1)*q2(1) + q2(2)*q2(2) + q2(3)*q2(3) )
        dkmod_save=dkmod_save+modq2
        IF (modq1 >1.d-6 .AND. modq2 > 1.d-6) THEN
           ps = ( q1(1)*q2(1) + q1(2)*q2(2) + q1(3)*q2(3) ) / &
                   modq1 / modq2
           high_sym_path(n) = (ABS(ps-1.d0) >1.0d-4)
        ENDIF
!
!  The gamma point is a high symmetry point
!
        IF (SQRT(disp_q(1,n)**2+disp_q(2,n)**2+disp_q(3,n)**2) < 1.0d-9) &
                                              high_sym_path(n)=.TRUE.
     END IF
  END DO
  dkmod_save=dkmod_save/(disp_nqs-1)

  RETURN
END SUBROUTINE set_paths_disp
