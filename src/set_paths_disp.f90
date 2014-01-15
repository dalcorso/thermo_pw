!
! Copyright (C) 2013-2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_paths_disp()
!
!  This subroutine computes the paths for (band or phonon) dispersions 
!  from the label letters or the auxiliary (k or q) points given in
!  input. Note that this routine should be called after the calculation
!  of the geometry. If the latter is changing the input k or q points
!  should be given as letters to be computed here with the correct
!  geometry.
!
  USE control_paths, ONLY : nqaux, xqaux, wqaux, npk_label, letter, &
                            label_list, nqaux, q_in_band_form, &
                            q_in_cryst_coord, q2d, point_label_type, &
                            label_disp_q
  USE cell_base,     ONLY : ibrav, celldm, bg
  USE ifc,           ONLY : disp_nqs, disp_q, disp_wq
  USE bz_form,       ONLY : transform_label_coord

  IMPLICIT NONE
  CHARACTER(len=80) :: k_points = 'tpiba'
  INTEGER :: i, j
 
  IF (nqaux==0) CALL errore('set_paths_disp','path_not_set',1)

  IF ( ALLOCATED(disp_q) ) DEALLOCATE (disp_q)
  IF ( ALLOCATED(disp_wq) ) DEALLOCATE (disp_wq)

  IF (q_in_cryst_coord) k_points='crystal'

  IF ( npk_label > 0 ) &
     CALL transform_label_coord(ibrav, celldm, xqaux, letter, &
     label_list, npk_label, nqaux, k_points, point_label_type )
   IF (q_in_cryst_coord)  CALL cryst_to_cart(nqaux,xqaux,bg,+1)

   IF (q2d) THEN
      disp_nqs=wqaux(2)*wqaux(3)
      ALLOCATE(disp_q(3,disp_nqs))
      ALLOCATE(disp_wq(disp_nqs))
      CALL generate_k_in_plane(nqaux, xqaux, wqaux, disp_q, &
                                 disp_wq, disp_nqs)
   ELSEIF (q_in_band_form) THEN
      disp_nqs=SUM(wqaux(1:nqaux-1))+1
      DO i=1,nqaux-1
         IF (wqaux(i)==0) disp_nqs=disp_nqs+1
      ENDDO
      ALLOCATE(disp_q(3,disp_nqs))
      ALLOCATE(disp_wq(disp_nqs))
      CALL generate_k_along_lines(nqaux, xqaux, wqaux, disp_q, disp_wq, &
                                  disp_nqs)
      label_disp_q=0
      label_disp_q(1)=1
      DO i=2,nqaux
         label_disp_q(i)=label_disp_q(i-1)+wqaux(i-1)
      ENDDO
   ELSE
      disp_nqs=nqaux
      ALLOCATE(disp_q(3,disp_nqs))
      ALLOCATE(disp_wq(disp_nqs))
      disp_wq(:)=wqaux(:)
      disp_q(:,1:disp_nqs)=xqaux(:,1:disp_nqs)
      DO i=1,nqaux
         label_disp_q(i)=i
      ENDDO
   ENDIF

   RETURN
END SUBROUTINE set_paths_disp
