!
! Copyright (C) 2015 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE gnuplot_color
!
!  This module provides some support for using colors in gnuplot.
!  gnuplot supports 112 colors. This modules allows to define their
!  names in the gnuplot scripts and use these names in the commands.
!  Moreover it provides routines to print these colors to test them
!  on the screen and on the 
!
    USE kinds, ONLY : DP
    USE io_global, ONLY : ionode
    USE gnuplot, ONLY : gnuplot_write_command

    SAVE
    PRIVATE

    INTEGER, PARAMETER :: nbase_colors=14 
    CHARACTER(LEN=20) :: base_colors(nbase_colors)

    INTEGER, PARAMETER :: ndark_light_colors=13
    CHARACTER(LEN=20) :: dark_light_colors(ndark_light_colors)

    INTEGER, PARAMETER :: nred_colors=9
    CHARACTER(LEN=20) :: red_colors(nred_colors)

    INTEGER, PARAMETER :: norange_colors=12
    CHARACTER(LEN=20) :: orange_colors(norange_colors)

    INTEGER, PARAMETER :: nyellow_colors=5
    CHARACTER(LEN=20) :: yellow_colors(nyellow_colors)

    INTEGER, PARAMETER :: ngreen_colors=11
    CHARACTER(LEN=20) :: green_colors(ngreen_colors)

    INTEGER, PARAMETER :: nblue_colors=13
    CHARACTER(LEN=20) :: blue_colors(nblue_colors)

    INTEGER, PARAMETER :: nviolet_colors=5
    CHARACTER(LEN=20) :: violet_colors(nviolet_colors)

    INTEGER, PARAMETER :: ngray_colors=28
    CHARACTER(LEN=20) :: gray_colors(ngray_colors)

    INTEGER, PARAMETER :: nall_colors=82
    CHARACTER(LEN=20) :: all_colors(nall_colors)

    PUBLIC  gnuplot_set_base_colors, gnuplot_set_dark_light_colors,     &
            gnuplot_set_reds, gnuplot_set_oranges,                      &
            gnuplot_set_yellows, gnuplot_set_greens, gnuplot_set_blues, &
            gnuplot_set_violets, gnuplot_set_grays,                     &
            gnuplot_set_all_colors, nbase_colors, base_colors,          &
            ndark_light_colors, dark_light_colors,                      &
            nred_colors, red_colors, norange_colors, orange_colors,     &
            nyellow_colors, yellow_colors, ngreen_colors, green_colors, &
            nblue_colors, blue_colors, gnuplot_set_blue,                &
            nviolet_colors, violet_colors, ngray_colors, gray_colors,   &
            nall_colors, all_colors, convert_color_name

CONTAINS

   SUBROUTINE gnuplot_set_base_colors()

   IMPLICIT NONE

   base_colors(1)='white'
   base_colors(2)='red'
   base_colors(3)='orange_red'
   base_colors(4)='orange'
   base_colors(5)='gold'
   base_colors(6)='yellow'
   base_colors(7)='greenyellow'
   base_colors(8)='green'
   base_colors(9)='cyan'
   base_colors(10)='blue'
   base_colors(11)='slateblue1'
   base_colors(12)='violet'
   base_colors(13)='magenta'
   base_colors(14)='black'

   CALL gnuplot_write_command('color_white="white"',.FALSE.)
   CALL gnuplot_write_command('color_red="red"',.FALSE.)
   CALL gnuplot_write_command('color_orange_red="orange-red"',.FALSE.)
   CALL gnuplot_write_command('color_orange="orange"',.FALSE.)
   CALL gnuplot_write_command('color_gold="gold"',.FALSE.)
   CALL gnuplot_write_command('color_yellow="yellow"',.FALSE.)
   CALL gnuplot_write_command('color_greenyellow="greenyellow"',.FALSE.)
   CALL gnuplot_write_command('color_green="green"',.FALSE.)
   CALL gnuplot_write_command('color_cyan="cyan"',.FALSE.)
   CALL gnuplot_write_command('color_blue="blue"',.FALSE.)
   CALL gnuplot_write_command('color_slateblue1="slateblue1"',.FALSE.)
   CALL gnuplot_write_command('color_violet="violet"',.FALSE.)
   CALL gnuplot_write_command('color_magenta="magenta"',.FALSE.)
   CALL gnuplot_write_command('color_black="black"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_base_colors

   SUBROUTINE gnuplot_set_dark_light_colors()

   IMPLICIT NONE

   dark_light_colors(1)='light_red'
   dark_light_colors(2)='dark_red'
   dark_light_colors(3)='dark_orange'
   dark_light_colors(4)='dark_yellow'
   dark_light_colors(5)='light_green'
   dark_light_colors(6)='dark_green'
   dark_light_colors(7)='light_blue'
   dark_light_colors(8)='dark_blue'
   dark_light_colors(9)='light_cyan'
   dark_light_colors(10)='dark_cyan'
   dark_light_colors(11)='dark_violet'
   dark_light_colors(12)='light_magenta'
   dark_light_colors(13)='dark_magenta'

   CALL gnuplot_write_command('color_light_red="light-red"',.FALSE.)
   CALL gnuplot_write_command('color_dark_red="dark-red"',.FALSE.)
   CALL gnuplot_write_command('color_dark_orange="dark-orange"',.FALSE.)
   CALL gnuplot_write_command('color_dark_yellow="dark-yellow"',.FALSE.)
   CALL gnuplot_write_command('color_light_green="light-green"',.FALSE.)
   CALL gnuplot_write_command('color_dark_green="dark-green"',.FALSE.)
   CALL gnuplot_write_command('color_light_blue="light-blue"',.FALSE.)
   CALL gnuplot_write_command('color_dark_blue="dark-blue"',.FALSE.)
   CALL gnuplot_write_command('color_light_cyan="light-cyan"',.FALSE.)
   CALL gnuplot_write_command('color_dark_cyan="dark-cyan"',.FALSE.)
   CALL gnuplot_write_command('color_dark_violet="dark-violet"',.FALSE.)
   CALL gnuplot_write_command('color_light_magenta="light-magenta"',.FALSE.)
   CALL gnuplot_write_command('color_dark_magenta="dark-magenta"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_dark_light_colors

   SUBROUTINE gnuplot_set_reds()

   IMPLICIT NONE

   red_colors(1)='coral'
   red_colors(2)='light_coral'
   red_colors(3)='light_pink'
   red_colors(4)='pink'
   red_colors(5)='dark_pink'
   red_colors(6)='light_salmon'
   red_colors(7)='salmon'
   red_colors(8)='dark_salmon'
   red_colors(9)='orangered4'
   
   CALL gnuplot_write_command('color_coral="coral"',.FALSE.)
   CALL gnuplot_write_command('color_light_coral="light-coral"',.FALSE.)
   CALL gnuplot_write_command('color_light_pink="light-pink"',.FALSE.)
   CALL gnuplot_write_command('color_pink="pink"',.FALSE.)
   CALL gnuplot_write_command('color_dark_pink="dark-pink"',.FALSE.)
   CALL gnuplot_write_command('color_light_salmon="light-salmon"',.FALSE.)
   CALL gnuplot_write_command('color_salmon="salmon"',.FALSE.)
   CALL gnuplot_write_command('color_dark_salmon="dark-salmon"',.FALSE.)
   CALL gnuplot_write_command('color_orangered4="orangered4"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_reds

   SUBROUTINE gnuplot_set_oranges()

   IMPLICIT NONE

   orange_colors(1)='brown'
   orange_colors(2)='brown4'
   orange_colors(3)='sienna1'
   orange_colors(4)='sienna4'
   orange_colors(5)='tan1'
   orange_colors(6)='sandybrown'
   orange_colors(7)='bisque'
   orange_colors(8)='antiquewhite'
   orange_colors(9)='beige'
   orange_colors(10)='khaki'
   orange_colors(11)='dark_khaki'
   orange_colors(12)='khaki1'
   
   CALL gnuplot_write_command('color_brown="brown"',.FALSE.)
   CALL gnuplot_write_command('color_brown4="brown4"',.FALSE.)
   CALL gnuplot_write_command('color_sienna1="sienna1"',.FALSE.)
   CALL gnuplot_write_command('color_sienna4="sienna4"',.FALSE.)
   CALL gnuplot_write_command('color_tan1="tan1"',.FALSE.)
   CALL gnuplot_write_command('color_sandybrown="sandybrown"',.FALSE.)
   CALL gnuplot_write_command('color_bisque="bisque"',.FALSE.)
   CALL gnuplot_write_command('color_antiquewhite="antiquewhite"',.FALSE.)
   CALL gnuplot_write_command('color_beige="beige"',.FALSE.)
   CALL gnuplot_write_command('color_khaki="khaki"',.FALSE.)
   CALL gnuplot_write_command('color_dark_khaki="dark-khaki"',.FALSE.)
   CALL gnuplot_write_command('color_khaki1="khaki1"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_oranges

   SUBROUTINE gnuplot_set_yellows()

   IMPLICIT NONE

   yellow_colors(1)='lemonchiffon'
   yellow_colors(2)='goldenrod'
   yellow_colors(3)='light_goldenrod'
   yellow_colors(4)='dark_goldenrod'
   yellow_colors(5)='yellow4'

   CALL gnuplot_write_command('color_lemonchiffon="lemonchiffon"',.FALSE.)
   CALL gnuplot_write_command('color_goldenrod="goldenrod"',.FALSE.)
   CALL gnuplot_write_command('color_light_goldenrod="light-goldenrod"',.FALSE.)
   CALL gnuplot_write_command('color_dark_goldenrod="dark-goldenrod"',.FALSE.)
   CALL gnuplot_write_command('color_yellow4="yellow4"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_yellows

   SUBROUTINE gnuplot_set_greens()

   IMPLICIT NONE

   green_colors(1)='olive'
   green_colors(2)='dark_olivegreen'
   green_colors(3)='chartreuse'
   green_colors(4)='dark_chartreuse'
   green_colors(5)='web_green'
   green_colors(6)='spring_green'
   green_colors(7)='dark_spring_green'
   green_colors(8)='forest_green'
   green_colors(9)='sea_green'
   green_colors(10)='seagreen'
   green_colors(11)='honeydew'

   CALL gnuplot_write_command('color_olive="olive"',.FALSE.)
   CALL gnuplot_write_command('color_dark_olivegreen="dark-olivegreen"',.FALSE.)
   CALL gnuplot_write_command('color_chartreuse="chartreuse"',.FALSE.)
   CALL gnuplot_write_command('color_dark_chartreuse="dark-chartreuse"',.FALSE.)
   CALL gnuplot_write_command('color_web_green="web-green"',.FALSE.)
   CALL gnuplot_write_command('color_spring_green="spring-green"',.FALSE.)
   CALL gnuplot_write_command('color_dark_spring_green="dark-spring-green"',&
                                                                  .FALSE.)
   CALL gnuplot_write_command('color_forest_green="forest-green"',.FALSE.)
   CALL gnuplot_write_command('color_sea_green="sea-green"',.FALSE.)
   CALL gnuplot_write_command('color_seagreen="seagreen"',.FALSE.)
   CALL gnuplot_write_command('color_honeydew="honeydew"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_greens

   SUBROUTINE gnuplot_set_blues()

   IMPLICIT NONE

   blue_colors(1)='aquamarine'
   blue_colors(2)='light_turquoise'
   blue_colors(3)='turquoise'
   blue_colors(4)='dark_turquoise'
   blue_colors(5)='skyblue'
   blue_colors(6)='slategray'
   blue_colors(7)='slategrey'
   blue_colors(8)='steelblue'
   blue_colors(9)='midnight_blue'
   blue_colors(10)='navy'
   blue_colors(11)='medium_blue'
   blue_colors(12)='web_blue'
   blue_colors(13)='royalblue'

   CALL gnuplot_write_command('color_aquamarine="aquamarine"',.FALSE.)
   CALL gnuplot_write_command('color_light_turquoise="light-turquoise"',.FALSE.)
   CALL gnuplot_write_command('color_turquoise="turquoise"',.FALSE.)
   CALL gnuplot_write_command('color_dark_turquoise="dark-turquoise"',.FALSE.)
   CALL gnuplot_write_command('color_skyblue="skyblue"',.FALSE.)
   CALL gnuplot_write_command('color_slategray="slategray"',.FALSE.)
   CALL gnuplot_write_command('color_slategrey="slategrey"',.FALSE.)
   CALL gnuplot_write_command('color_steelblue="steelblue"',.FALSE.)
   CALL gnuplot_write_command('color_midnight_blue="midnight-blue"',.FALSE.)
   CALL gnuplot_write_command('color_medium_blue="medium-blue"',.FALSE.)
   CALL gnuplot_write_command('color_navy="navy"',.FALSE.)
   CALL gnuplot_write_command('color_web_blue="web-blue"',.FALSE.)
   CALL gnuplot_write_command('color_royalblue="royalblue"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_blues

   SUBROUTINE gnuplot_set_violets()

   IMPLICIT NONE

   violet_colors(1)='orchid'
   violet_colors(2)='orchid4'
   violet_colors(3)='plum'
   violet_colors(4)='dark_plum'
   violet_colors(5)='mediumpurple3'

   CALL gnuplot_write_command('color_orchid="orchid"',.FALSE.)
   CALL gnuplot_write_command('color_orchid4="orchid4"',.FALSE.)
   CALL gnuplot_write_command('color_plum="plum"',.FALSE.)
   CALL gnuplot_write_command('color_dark_plum="dark-plum"',.FALSE.)
   CALL gnuplot_write_command('color_mediumpurple3="mediumpurple3"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_violets

   SUBROUTINE gnuplot_set_grays()

   IMPLICIT NONE

   gray_colors(1)='gray0'
   gray_colors(2)='grey0'
   gray_colors(3)='gray10'
   gray_colors(4)='grey10'
   gray_colors(5)='gray20'
   gray_colors(6)='grey20'
   gray_colors(7)='gray30'
   gray_colors(8)='grey30'
   gray_colors(9)='gray40'
   gray_colors(10)='grey40'
   gray_colors(11)='gray50'
   gray_colors(12)='grey50'
   gray_colors(13)='gray60'
   gray_colors(14)='grey60'
   gray_colors(15)='dark_gray'
   gray_colors(16)='dark_grey'
   gray_colors(17)='gray70'
   gray_colors(18)='grey70'
   gray_colors(19)='gray'
   gray_colors(20)='grey'
   gray_colors(21)='gray80'
   gray_colors(22)='grey80'
   gray_colors(23)='light_gray'
   gray_colors(24)='light_grey'
   gray_colors(25)='gray90'
   gray_colors(26)='grey90'
   gray_colors(27)='gray100'
   gray_colors(28)='grey100'

   CALL gnuplot_write_command('color_gray0="gray0"',.FALSE.)
   CALL gnuplot_write_command('color_grey0="grey0"',.FALSE.)
   CALL gnuplot_write_command('color_gray10="gray10"',.FALSE.)
   CALL gnuplot_write_command('color_grey10="grey10"',.FALSE.)
   CALL gnuplot_write_command('color_gray20="gray20"',.FALSE.)
   CALL gnuplot_write_command('color_grey20="grey20"',.FALSE.)
   CALL gnuplot_write_command('color_gray30="gray30"',.FALSE.)
   CALL gnuplot_write_command('color_grey30="grey30"',.FALSE.)
   CALL gnuplot_write_command('color_gray40="gray40"',.FALSE.)
   CALL gnuplot_write_command('color_grey40="grey40"',.FALSE.)
   CALL gnuplot_write_command('color_gray50="gray50"',.FALSE.)
   CALL gnuplot_write_command('color_grey50="grey50"',.FALSE.)
   CALL gnuplot_write_command('color_gray60="gray60"',.FALSE.)
   CALL gnuplot_write_command('color_grey60="grey60"',.FALSE.)
   CALL gnuplot_write_command('color_dark_gray="dark-gray"',.FALSE.)
   CALL gnuplot_write_command('color_dark_grey="dark-grey"',.FALSE.)
   CALL gnuplot_write_command('color_gray70="gray70"',.FALSE.)
   CALL gnuplot_write_command('color_grey70="grey70"',.FALSE.)
   CALL gnuplot_write_command('color_gray="gray"',.FALSE.)
   CALL gnuplot_write_command('color_grey="grey"',.FALSE.)
   CALL gnuplot_write_command('color_gray80="gray80"',.FALSE.)
   CALL gnuplot_write_command('color_grey80="grey80"',.FALSE.)
   CALL gnuplot_write_command('color_light_gray="light-gray"',.FALSE.)
   CALL gnuplot_write_command('color_light_grey="light-grey"',.FALSE.)
   CALL gnuplot_write_command('color_gray90="gray90"',.FALSE.)
   CALL gnuplot_write_command('color_grey90="grey90"',.FALSE.)
   CALL gnuplot_write_command('color_gray100="gray100"',.FALSE.)
   CALL gnuplot_write_command('color_grey100="grey100"',.FALSE.)

   RETURN
   END SUBROUTINE gnuplot_set_grays

   SUBROUTINE gnuplot_set_all_colors()

   IMPLICIT NONE
   INTEGER :: counter

   CALL gnuplot_set_base_colors()
   CALL gnuplot_set_dark_light_colors()
   CALL gnuplot_set_reds()
   CALL gnuplot_set_oranges()
   CALL gnuplot_set_yellows()
   CALL gnuplot_set_greens()
   CALL gnuplot_set_blues()
   CALL gnuplot_set_violets()

   all_colors(1:nbase_colors)=base_colors(:)
   counter=nbase_colors
   all_colors(counter+1:counter+ndark_light_colors)=dark_light_colors(:)
   counter=counter+ndark_light_colors
   all_colors(counter+1:counter+nred_colors)=red_colors(:)
   counter=counter+nred_colors
   all_colors(counter+1:counter+norange_colors)=orange_colors(:)
   counter=counter+norange_colors
   all_colors(counter+1:counter+nyellow_colors)=yellow_colors(:)
   counter=counter+nyellow_colors
   all_colors(counter+1:counter+ngreen_colors)=green_colors(:)
   counter=counter+ngreen_colors
   all_colors(counter+1:counter+nblue_colors)=blue_colors(:)
   counter=counter+nblue_colors
   all_colors(counter+1:counter+nviolet_colors)=violet_colors(:)

   RETURN
   END SUBROUTINE gnuplot_set_all_colors

   SUBROUTINE convert_color_name(inname, outname)
!
!   this routine substitutes the _ with a - in the color name, in order
!   to print it.
!
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN) :: inname
   CHARACTER(LEN=*), INTENT(INOUT) :: outname

   CHARACTER(LEN=1) :: c

   INTEGER :: ilen, ic

   outname=''
   DO ilen=1,LEN(inname)
      c=inname(ilen:ilen)
      ic= ICHAR(c)
      IF (ic==95) c='-'
      outname=TRIM(outname)//c
   ENDDO

   RETURN
   END SUBROUTINE convert_color_name



END MODULE gnuplot_color
