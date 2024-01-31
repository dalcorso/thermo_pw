!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!--------------------------------------------------------------------------
SUBROUTINE plot_2d_bz(ibrav, celldm, at, bg, xk, wk, npkt, letter, &
           letter_path, npk_label, label_list, asy_filename)
!--------------------------------------------------------------------------
!
!  This subroutine plots the BZ of a given solid. It
!  receives as input:
!  ibrav  : a Bravais lattice code specified by the routine latgen
!  celldm(6) : the six parameters that specify the Bravais lattice size
!  at,bg,    : the direct and reciprocal lattice vectors
!  xk, wk    : the coordinates of the k points of a path
!  npkt,     : the number of k points that define a path
!  letter,   : a list of letters to identify the points on the path
!  letter_path,  : the same letter list, with the letter written on the plot 
!  npk_label : the number of letters
!  label_list : an array that says which k point is given as a letter
!  asy_filename : the name of a file where an asymptote script that allows
!                 to show the BZ is written
!  It produces as output an asymptote script that can be used to generate
!  a pdf file of the BZ and a plot of the path in the BZ.
!
!
USE kinds, ONLY : DP
USE bz_2d_form, ONLY : bz_2d, init_2d_bz, allocate_2d_bz, find_2d_bz_type, &
                       deallocate_2d_bz
USE bz_asy_mod, ONLY : bz_asy, allocate_2d_bz_asy, init_2d_bz_asy, &
                       deallocate_bz_asy
USE control_paths, ONLY : npx
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: npkt, npk_label
REAL(DP), INTENT(IN) :: celldm(6), at(3,3), bg(3,3) 
REAL(DP), INTENT(INOUT) :: xk(3,npkt)
INTEGER, INTENT(IN) :: label_list(npk_label), wk(npkt) 
CHARACTER(LEN=3),INTENT(IN) :: letter(npk_label), letter_path(npk_label)
INTEGER, INTENT(IN) :: ibrav
CHARACTER(LEN=*), INTENT(IN) :: asy_filename
INTEGER :: ik
REAL(DP) :: celldm_2d(3)
TYPE(bz_2d) :: bz_2d_struc
TYPE(bz_asy) :: bz_asy_struc
INTEGER :: ibz
!
! load the information on this Brillouin zone
!
CALL allocate_2d_bz(ibz, bz_2d_struc, celldm, at, bg, celldm_2d)
IF (npx > bz_2d_struc%npx) bz_2d_struc%npx=npx
CALL init_2d_bz(bz_2d_struc, celldm_2d)

WRITE(stdout,'(/,5x,"2D Brillouin zone type",i5,/)') ibz
!
!   allocate the variables to generate the asymptote script
!
CALL allocate_2d_bz_asy(bz_2d_struc, bz_asy_struc )
CALL init_2d_bz_asy(bz_2d_struc, bz_asy_struc, ibz, celldm_2d)
CALL generate_2d_asy_figure(xk, wk, npkt, letter, letter_path, npk_label, &
                        label_list, bz_2d_struc, bz_asy_struc, asy_filename)

CALL deallocate_bz_asy(bz_asy_struc)
CALL deallocate_2d_bz(bz_2d_struc)
RETURN
END SUBROUTINE plot_2d_bz

!--------------------------------------------------------------------------
SUBROUTINE generate_2d_asy_figure(xk, wk, npkt, letter, letter_path, &
                      npk_label, label_list, bz_2d_struc, bz_asy_struc, &
                      asy_filename)
!--------------------------------------------------------------------------
!
!  This subroutine generates an input script for asymptote that
!  plots the BZ
!
USE kinds, ONLY : DP
USE asy, ONLY : asy_open2dplot, asy_closeplot, asy_write_2d_point, &
                asy_writesurface, asy_2d_join, asy_put2dlabel, &
                asy_write_2d_surface, asy_2d_plotaxis
USE bz_2d_form, ONLY : bz_2d, find_2d_letter_coordinate
USE bz_asy_mod, ONLY : bz_asy, find_2d_letter_position
USE control_paths, ONLY : q_in_band_form 

IMPLICIT NONE
INTEGER, INTENT(IN) :: npkt, npk_label
CHARACTER(LEN=3), INTENT(IN) :: letter(npk_label), letter_path(npk_label)
INTEGER, INTENT(IN) :: label_list(npk_label), wk(npkt)
REAL(DP), INTENT(INOUT) :: xk(3,npkt)
TYPE(bz_2d), INTENT(IN) :: bz_2d_struc
TYPE(bz_asy), INTENT(IN) :: bz_asy_struc
CHARACTER(LEN=*), INTENT(IN) :: asy_filename

INTEGER :: i, ik
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=7) :: label
CHARACTER(LEN=3) :: let_pos
REAL(DP) :: ak(3)
REAL(DP) :: letter_coordinates(3)

CALL asy_open2dplot(asy_filename)

DO i=1,bz_2d_struc%nvertices
   label="G"//TRIM(int_to_char(i))
   CALL asy_write_2d_point(bz_2d_struc%vertex_coord(:,i), label)
ENDDO

CALL asy_write_2d_surface(bz_2d_struc%nvertices)

DO i = 1, npk_label
   CALL find_2d_letter_coordinate(bz_2d_struc, letter(i), letter_coordinates )
   xk(:,label_list(i))=letter_coordinates
   CALL find_2d_letter_position(bz_2d_struc, bz_asy_struc, letter(i), let_pos)
   CALL asy_put2dlabel(letter_path(label_list(i)), letter_coordinates, let_pos)
END DO

DO ik=2, npkt
   IF (wk(ik-1)/=0) CALL asy_2d_join(xk(:,ik), xk(:,ik-1), 'blue')
ENDDO
!
IF (.NOT.q_in_band_form) THEN
   DO ik=1, npkt
      IF (letter_path(ik) /= '') THEN
         CALL find_2d_letter_position(bz_2d_struc, bz_asy_struc, &
         letter_path(ik), let_pos)
         CALL asy_put2dlabel(letter_path(ik), xk(:,ik), let_pos)
      END IF
   END DO
END IF
!
!  find where the coordinate axis intercept the BZ
!
ak(1)=bz_2d_struc%xi(1)
ak(2)=bz_2d_struc%yi(2)
ak(3)=0.0_DP

CALL asy_2d_plotaxis(ak)
 
CALL asy_closeplot()

RETURN
END SUBROUTINE generate_2d_asy_figure

