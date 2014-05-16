!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plot_bz(ibrav, celldm, at, bg, point_type, &
                   xk, wk, npkt, letter, letter_path, npk_label, &
                   label_list, asy_filename)
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
USE bz_form, ONLY : bz, init_bz, allocate_bz, find_letter_coordinate, &
                    find_bz_type, set_label_type
USE bz_asy_mod, ONLY : bz_asy, allocate_bz_asy, init_bz_asy
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: npkt, npk_label
REAL(DP), INTENT(IN) :: celldm(6), at(3,3), bg(3,3) 
REAL(DP), INTENT(INOUT) :: xk(3,npkt)
INTEGER, INTENT(IN) :: label_list(npk_label), wk(npkt) 
CHARACTER(LEN=3),INTENT(IN) :: letter(npk_label), letter_path(npk_label)
INTEGER, INTENT(IN) :: ibrav
CHARACTER(LEN=*), INTENT(IN) :: asy_filename
CHARACTER(LEN=10), INTENT(IN) :: point_type 
INTEGER :: ik
TYPE(bz) :: bz_struc
TYPE(bz_asy) :: bz_asy_struc
INTEGER :: bzt
!
!  find which BZ is appropriate for this Bravais lattice
!
CALL find_bz_type(ibrav, celldm, bzt)

WRITE(stdout,'(/,5x,"Brillouin zone type",i5,/)') bzt
!
! load the information on this Brillouin zone
!
CALL set_label_type(bz_struc, point_type)
CALL allocate_bz(ibrav, bzt, bz_struc, celldm, at, bg)
CALL init_bz(bz_struc)
!
!   allocate the variables to generate the asymptote script
!
CALL allocate_bz_asy(bz_struc, bz_asy_struc )
CALL init_bz_asy(bz_struc, bz_asy_struc)
CALL generate_asy_figure(xk, wk, npkt, letter, letter_path, npk_label, &
                        label_list, bz_struc, bz_asy_struc, asy_filename)

RETURN
END SUBROUTINE plot_bz

SUBROUTINE generate_asy_figure(xk, wk, npkt, letter, letter_path, &
                      npk_label, label_list, bz_struc, bz_asy_struc, &
                      asy_filename)
!
!  This subroutine generates an input script for asymptote that
!  plots the BZ
!
USE kinds, ONLY : DP
USE asy, ONLY : asy_openplot, asy_closeplot, asy_writepoint, asy_writesurface, &
                asy_plotaxis, asy_join, asy_putlabel
USE bz_form, ONLY : bz, find_letter_coordinate
USE bz_asy_mod, ONLY : bz_asy, find_letter_position
USE control_paths, ONLY : q_in_band_form 

IMPLICIT NONE
INTEGER, INTENT(IN) :: npkt, npk_label
CHARACTER(LEN=3), INTENT(IN) :: letter(npk_label), letter_path(npk_label)
INTEGER, INTENT(IN) :: label_list(npk_label), wk(npkt)
REAL(DP), INTENT(INOUT) :: xk(3,npkt)
TYPE(bz), INTENT(IN) :: bz_struc
TYPE(bz_asy), INTENT(IN) :: bz_asy_struc
CHARACTER(LEN=*), INTENT(IN) :: asy_filename

INTEGER :: i, ik
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=7) :: label
CHARACTER(LEN=3) :: let_pos
REAL(DP) :: x0(3), vect(3), ak(3), xk1(3), xmod, xmod1
REAL(DP) :: letter_coordinates(3)

CALL asy_openplot(asy_filename)

DO i=1,bz_struc%nvertices
   label="F"//TRIM(int_to_char(i))
   CALL asy_writepoint(bz_struc%vertex_coord(:,i), label)
ENDDO

DO i=1,bz_struc%nfaces
   IF (bz_asy_struc%visible(i)) CALL asy_writesurface(bz_struc%indsur(:,i))
ENDDO

DO i = 1, npk_label
   CALL find_letter_coordinate(bz_struc, letter(i), &
                                          letter_coordinates )
   xk(:,label_list(i))=letter_coordinates
   CALL find_letter_position(bz_struc, bz_asy_struc, letter(i), let_pos)
   CALL asy_putlabel(letter_path(label_list(i)), letter_coordinates, let_pos)
END DO

DO ik=2, npkt
   xmod = SQRT(xk(1,ik)**2 + xk(2,ik)**2 + xk(3,ik)**2)
   xmod1 = SQRT(xk(1,ik-1)**2 + xk(2,ik-1)**2 + xk(3,ik-1)**2)
   IF (xmod < 1.d-7 .OR. xmod1 < 1.d-7) THEN
      IF (wk(ik-1)/=0) CALL asy_join(xk(:,ik), xk(:,ik-1), 'red+dashed')
   ELSE
      IF (wk(ik-1)/=0) CALL asy_join(xk(:,ik), xk(:,ik-1), 'red')
   ENDIF
ENDDO

IF (.NOT.q_in_band_form) THEN
   DO ik=1, npkt
      IF (letter_path(ik) /= '') THEN
         CALL find_letter_position(bz_struc, bz_asy_struc, letter_path(ik), &
                                   let_pos)
         CALL asy_putlabel(letter_path(ik), xk(:,ik), let_pos)
      END IF
   END DO
END IF
!
!  find where the coordinate axis intercept the BZ
!
ak(1)=bz_struc%xi(1)
ak(2)=bz_struc%yi(2)
ak(3)=bz_struc%zi(3)

CALL asy_plotaxis(ak)
CALL asy_closeplot()

RETURN
END SUBROUTINE generate_asy_figure
