!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE color_mod
!
!  This module provides some support for using colors.
!  Presently it offers only the 
!
    USE kinds, ONLY : DP

    SAVE
    PRIVATE

    INTEGER, PARAMETER :: eight_colors=8
    CHARACTER(LEN=12) :: color(eight_colors)


    PUBLIC  color, set_colors 

CONTAINS
!--------------------------------------------------------------------
   SUBROUTINE set_colors()
!--------------------------------------------------------------------

   IMPLICIT NONE

   color(1)='color_red'
   color(2)='color_green'
   color(3)='color_blue'
   color(4)='color_yellow'
   color(5)='color_pink'
   color(6)='color_cyan'
   color(7)='color_orange'
   color(8)='color_black'

END SUBROUTINE set_colors

END MODULE color_mod
