!
! Copyright (C) 2017 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE magnetic_point_group
!
!  This module contains variables and routines to deal with the 
!  crystallographic magnetic point group symmetry. 
!
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC find_mag_group_code, mag_group_name, magnetic_type, is_mag_group, &
         find_mag_group_code_ext, set_mag_group_subgroup, &
         set_mag_group_subgroup_ext, mag_group_index_from_ext, &
         find_group_subgroup_ext

CONTAINS

SUBROUTINE find_mag_group_code(group_code,subgroup_code,mag_code)
!
! This routine receives the codes of a group and of a subgroup and
! gives the code of the magnetic group that corresponds to them.
! The magnetic codes are the following:
!
!
!    1               C_1 (11')  DP        2              C_i (-11')  DP
!    3               C_s (m1')  DP        4               C_2 (21')  DP
!    5               C_3 (31')  DP        6               C_4 (41')  DP
!    7               C_6 (61')  DP        8             D_2 (2221')  DP
!    9              D_3 (321')  DP       10             D_4 (4221')  DP
!   11             D_6 (6221')  DP       12            C_2v (mm21')  DP
!   13             C_3v (3m1')  DP       14            C_4v (4mm1')  DP
!   15            C_6v (6mm1')  DP       16            C_2h (2/m1')  DP
!   17             C_3h (-61')  DP       18            C_4h (4/m1')  DP
!   19            C_6h (6/m1')  DP       20            D_2h (mmm1')  DP
!   21           D_3h (-62m1')  DP       22          D_4h (4/mmm1')  DP
!   23          D_6h (6/mmm1')  DP       24           D_2d (-42m1')  DP
!   25            D_3d (-3m1')  DP       26              S_4 (-41')  DP
!   27              S_6 (-31')  DP       28                  T (23)  DP
!   29             T_h (m-31')  DP       30            T_d (-43m1')  DP
!   31               O (4321')  DP       32            O_h (m-3m1')  DP
!   33            C_1(C_1) (1)  F        34           C_i(C_i) (-1)  F 
!   35            C_s(C_s) (m)  F        36            C_2(C_2) (2)  F 
!   37            C_3(C_3) (3)  F        38            C_4(C_4) (4)  F 
!   39            C_6(C_6) (6)  F        40          D_2(D_2) (222)  AF
!   41           D_3(D_3) (32)  AF       42          D_4(D_4) (422)  AF
!   43          D_6(D_6) (622)  AF       44        C_2v(C_2v) (mm2)  AF
!   45         C_3v(C_3v) (3m)  AF       46        C_4v(C_4v) (4mm)  AF
!   47        C_6v(C_6v) (6mm)  AF       48        C_2h(C_2h) (2/m)  F 
!   49         C_3h(C_3h) (-6)  F        50        C_4h(C_4h) (4/m)  F 
!   51        C_6h(C_3h) (6/m)  F        52        D_2h(D_2h) (mmm)  AF
!   53       D_3h(D_3h) (-62m)  AF       54      D_4h(D_4h) (4/mmm)  AF
!   55      D_6h(D_6h) (6/mmm)  AF       56       D_2d(D_2d) (-42m)  AF
!   57        D_3d(D_3d) (-3m)  AF       58           S_4(S_4) (-4)  F 
!   59           S_6(S_6) (-3)  F        60               T(T) (23)  AF
!   61          T_h(T_h) (m-3)  AF       62         T_d(T_d) (-43m)  AF
!   63              O(O) (432)  AF       64         O_h(O_h) (m-3m)  AF
!   65          C_i(C_1) (-1')  AF       66           C_2(C_1) (2')  F 
!   67           C_s(C_1) (m')  F        68        C_2h(C_2) (2/m')  AF
!   69       C_2h(C_i) (2'/m')  F        70        C_2h(C_s) (2'/m)  AF
!   71        D_2(C_2) (22'2')  F        72       C_2v(C_2) (m'm'2)  F 
!   73       C_2v(C_s) (m'm2')  F        74       D_2h(C_2v) (m'mm)  AF
!   75      D_2h(C_2h) (m'm'm)  F        76      D_2h(D_2) (m'm'm')  AF
!   77           C_4(C_2) (4')  AF       78          S_4(C_2) (-4')  AF
!   79       C_4h(C_2h) (4'/m)  AF       80        C_4h(C_4) (4/m')  AF
!   81       C_4h(S_4) (4'/m')  AF       82        D_4(D_2) (4'2'2)  AF
!   83        D_4(C_4) (42'2')  F        84      C_4v(C_2v) (4'm'm)  AF
!   85       C_4v(C_4) (4m'm')  F        86     D_2d(C_2v) (-4'2'm)  AF
!   87      D_2d(D_2) (-4'2m')  AF       88      D_2d(S_4) (-42'm')  F 
!   89     D_4h(C_4v) (4/m'mm)  AF       90    D_4h(D_2h) (4'/mm'm)  AF
!   91   D_4h(D_2d) (4'/m'm'm)  AF       92    D_4h(C_4h) (4/mm'm')  F 
!   93    D_4h(D_4) (4/m'm'm')  AF       94          S_6(C_3) (-3')  AF
!   95          D_3(C_3) (32')  F        96         C_3v(C_3) (3m')  F 
!   97       D_3d(C_3v) (-3'm)  AF       98       D_3d(D_3) (-3'm')  AF
!   99        D_3d(S_6) (-3m')  F       100           C_6(C_3) (6')  AF
!  101         C_3h(C_3) (-6')  AF      102        C_6h(S_6) (6'/m)  AF
!  103        C_6h(C_6) (6/m')  AF      104      C_6h(C_3h) (6'/m')  AF
!  105        D_6(D_3) (6'22')  AF      106        D_6(C_6) (62'2')  F 
!  107      C_6v(C_3v) (6'm'm)  AF      108       C_6v(C_6) (6m'm')  F 
!  109      D_3h(D_3) (-6'm'2)  F       110     D_3h(C_3v) (-6'm2')  AF
!  111     D_3h(C_3h) (-6m'2')  AF      112     D_6h(C_6v) (6/m'mm)  AF
!  113     D_6h(D_3h) (6/mm'm)  AF      114   D_6h(D_3d) (6'/m'm'm)  AF
!  115    D_6h(C_6h) (6/mm'm')  F       116    D_6h(D_6) (6/m'm'm')  AF
!  117            T_h(T) (m'3)  AF      118            O(T) (4'32')  AF
!  119        T_d(T) (-4'-3m')  AF      120        O_h(T_d) (m'-3m)  AF
!  121        O_h(T_h) (m-3m')  AF      122        O_h(O) (m'-3'm')  AF
!
! If the group and subgroup codes do not correspond to any magnetic
! point group the output mag_code is zero.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: group_code, subgroup_code
INTEGER, INTENT(OUT) :: mag_code
INTEGER :: mag_code_tab(32,32)

INTEGER :: group(58), subgroup(58), icode

IF (group_code<1.OR.group_code>32.OR.subgroup_code<0.OR.&
    subgroup_code>32) CALL errore('find_mag_group_code',&
          'group or subgroup codes out of range',1)

CALL set_mag_group_subgroup(group,subgroup)
mag_code_tab=0
IF (subgroup_code==group_code) THEN
!
!   This is the case of the gray groups, time reversal is an element of
!   the group
   mag_code=group_code 
ELSEIF (subgroup_code==0) THEN
!
!  These are the groups without any time reversal operation
!
   mag_code=group_code+32
ELSE
!
!  These are the magnetic groups.
!
   DO icode=1,58
      mag_code_tab(group(icode),subgroup(icode))=64+icode
   ENDDO
   mag_code=mag_code_tab(group_code, subgroup_code)
ENDIF

RETURN
END SUBROUTINE find_mag_group_code


FUNCTION mag_group_name(code)

IMPLICIT NONE
INTEGER, INTENT(IN) :: code
CHARACTER(LEN=21) :: mag_group_name

CHARACTER(LEN=21) :: gname(122)

data gname /  "C_1 (11')            ", "C_i (-11')           ", &
              "C_s (m1')            ", "C_2 (21')            ", &
              "C_3 (31')            ", "C_4 (41')            ", &
              "C_6 (61')            ", "D_2 (2221')          ", &
              "D_3 (321')           ", "D_4 (4221')          ", &
              "D_6 (6221')          ", "C_2v (mm21')         ", &
              "C_3v (3m1')          ", "C_4v (4mm1')         ", &
              "C_6v (6mm1')         ", "C_2h (2/m1')         ", &
              "C_3h (-61')          ", "C_4h (4/m1')         ", &
              "C_6h (6/m1')         ", "D_2h (mmm1')         ", &
              "D_3h (-62m1')        ", "D_4h (4/mmm1')       ", &
              "D_6h (6/mmm1')       ", "D_2d (-42m1')        ", &
              "D_3d (-3m1')         ", "S_4 (-41')           ", &
              "S_6 (-31')           ", "T (23)               ", &
              "T_h (m-31')          ", "T_d (-43m1')         ", &
              "O (4321')            ", "O_h (m-3m1')         ", &
              "C_1(C_1) (1)         ", "C_i(C_i) (-1)        ", &
              "C_s(C_s) (m)         ", "C_2(C_2) (2)         ", &
              "C_3(C_3) (3)         ", "C_4(C_4) (4)         ", &
              "C_6(C_6) (6)         ", "D_2(D_2) (222)       ", &
              "D_3(D_3) (32)        ", "D_4(D_4) (422)       ", &
              "D_6(D_6) (622)       ", "C_2v(C_2v) (mm2)     ", &
              "C_3v(C_3v) (3m)      ", "C_4v(C_4v) (4mm)     ", &
              "C_6v(C_6v) (6mm)     ", "C_2h(C_2h) (2/m)     ", &
              "C_3h(C_3h) (-6)      ", "C_4h(C_4h) (4/m)     ", &
              "C_6h(C_3h) (6/m)     ", "D_2h(D_2h) (mmm)     ", &
              "D_3h(D_3h) (-62m)    ", "D_4h(D_4h) (4/mmm)   ", &
              "D_6h(D_6h) (6/mmm)   ", "D_2d(D_2d) (-42m)    ", &
              "D_3d(D_3d) (-3m)     ", "S_4(S_4) (-4)        ", &
              "S_6(S_6) (-3)        ", "T(T) (23)            ", &
              "T_h(T_h) (m-3)       ", "T_d(T_d) (-43m)      ", &
              "O(O) (432)           ", "O_h(O_h) (m-3m)      ", &
              "C_i(C_1) (-1')       ", "C_2(C_1) (2')        ", &
              "C_s(C_1) (m')        ", "C_2h(C_2) (2/m')     ", &
              "C_2h(C_i) (2'/m')    ", "C_2h(C_s) (2'/m)     ", &
              "D_2(C_2) (22'2')     ", "C_2v(C_2) (m'm'2)    ", &
              "C_2v(C_s) (m'm2')    ", "D_2h(C_2v) (m'mm)    ", &
              "D_2h(C_2h) (m'm'm)   ", "D_2h(D_2) (m'm'm')   ", &
              "C_4(C_2) (4')        ", "S_4(C_2) (-4')       ", &
              "C_4h(C_2h) (4'/m)    ", "C_4h(C_4) (4/m')     ", &
              "C_4h(S_4) (4'/m')    ", "D_4(D_2) (4'2'2)     ", &
              "D_4(C_4) (42'2')     ", "C_4v(C_2v) (4'm'm)   ", &
              "C_4v(C_4) (4m'm')    ", "D_2d(C_2v) (-4'2'm)  ", &
              "D_2d(D_2) (-4'2m')   ", "D_2d(S_4) (-42'm')   ", &
              "D_4h(C_4v) (4/m'mm)  ", "D_4h(D_2h) (4'/mm'm) ", &
              "D_4h(D_2d) (4'/m'm'm)", "D_4h(C_4h) (4/mm'm') ", &
              "D_4h(D_4) (4/m'm'm') ", "S_6(C_3) (-3')       ", &
              "D_3(C_3) (32')       ", "C_3v(C_3) (3m')      ", &
              "D_3d(C_3v) (-3'm)    ", "D_3d(D_3) (-3'm')    ", &
              "D_3d(S_6) (-3m')     ", "C_6(C_3) (6')        ", &
              "C_3h(C_3) (-6')      ", "C_6h(S_6) (6'/m)     ", &
              "C_6h(C_6) (6/m')     ", "C_6h(C_3h) (6'/m')   ", &
              "D_6(D_3) (6'22')     ", "D_6(C_6) (62'2')     ", &
              "C_6v(C_3v) (6'm'm)   ", "C_6v(C_6) (6m'm')    ", &
              "D_3h(D_3) (-6'm'2)   ", "D_3h(C_3v) (-6'm2')  ", &
              "D_3h(C_3h) (-6m'2')  ", "D_6h(C_6v) (6/m'mm)  ", &
              "D_6h(D_3h) (6/mm'm)  ", "D_6h(D_3d) (6'/m'm'm)", &
              "D_6h(C_6h) (6/mm'm') ", "D_6h(D_6) (6/m'm'm') ", &
              "T_h(T) (m'3)         ", "O(T) (4'32')         ", &
              "T_d(T) (-4'-3m')     ", "O_h(T_d) (m'-3m)     ", &
              "O_h(T_h) (m-3m')     ", "O_h(O) (m'-3'm')     " /

mag_group_name=gname(code)

RETURN
END FUNCTION mag_group_name

FUNCTION magnetic_type(code)

IMPLICIT NONE
INTEGER, INTENT(IN) :: code
CHARACTER(LEN=2) :: magnetic_type

CHARACTER(LEN=2) :: mtype(122)

data mtype  / "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", &
              "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", &
              "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", "DP", &
              "DP", "DP", "F ", "F ", "F ", "F ", "F ", "F ", "F ", "AF", &
              "AF", "AF", "AF", "AF", "AF", "AF", "AF", "F ", "F ", "F ", &
              "F ", "AF", "AF", "AF", "AF", "AF", "AF", "F ", "F ", "AF", &
              "AF", "AF", "AF", "AF", "AF", "F ", "F ", "AF", "F ", "AF", &
              "F ", "F ", "F ", "AF", "F ", "AF", "AF", "AF", "AF", "AF", &
              "AF", "AF", "F ", "AF", "F ", "AF", "AF", "F ", "AF", "AF", &
              "AF", "F ", "AF", "AF", "F ", "F ", "AF", "AF", "F ", "AF", &
              "AF", "AF", "AF", "AF", "AF", "F ", "AF", "F ", "F ", "AF", &
              "AF", "AF", "AF", "AF", "F ", "AF", "AF", "AF", "AF", "AF", &
              "AF", "AF" /

magnetic_type=mtype(code)

RETURN
END FUNCTION magnetic_type

SUBROUTINE set_mag_group_subgroup(group,subgroup)

IMPLICIT NONE
INTEGER, INTENT(OUT) :: group(58), subgroup(58)
INTEGER :: group_(58), subgroup_(58)

DATA group_    / 2,  4,  3, 16, 16, 16,  8, 12, 12, 20, 20, 20,  6, 26, 18, &
                18, 18, 10, 10, 14, 14, 24, 24, 24, 22, 22, 22, 22, 22, 27, &
                 9, 13, 25, 25, 25,  7, 17, 19, 19, 19, 11, 11, 15, 15, 21, &
                21, 21, 23, 23, 23, 23, 23, 29, 31, 30, 32, 32, 32  / 
DATA subgroup_ / 1,  1,  1,  4,  2,  3,  4,  4,  3, 12, 16,  8,  4,  4, 16, &
                 6, 26,  8,  6, 12,  6, 12,  8, 26, 14, 20, 24, 18, 10,  5, &
                 5,  5, 13,  9, 27,  5,  5, 27,  7, 17,  9,  7, 13,  7,  9, &
                13, 17, 15, 21, 25, 19, 11, 28, 28, 28, 30, 29, 31  / 

group=group_
subgroup=subgroup_

RETURN
END SUBROUTINE set_mag_group_subgroup

LOGICAL FUNCTION is_mag_group(code_group, code_subgroup)
!
!  This function receives the codes of group and of a subgroup and
!  returns .TRUE. if the subgroup is invariant and of index=2
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group, code_subgroup

INTEGER :: group(58), subgroup(58), icode

CALL set_mag_group_subgroup(group,subgroup)

is_mag_group=.FALSE.
DO icode=1,58
   is_mag_group = (group(icode)==code_group).AND.&
                          (subgroup(icode)==code_subgroup) 
   IF (is_mag_group) EXIT
ENDDO

RETURN
END FUNCTION is_mag_group

SUBROUTINE find_group_subgroup(group_code,subgroup_code, mag_code)
!
!  This routine recieves the code of a magnetic group (between 1 and 122)
!  and sets the code of the group and of the invariant subgroup that
!  define the magnetic group. For gray groups the code of the group and
!  of the subgroups are equal, for groups that do not contain the 
!  operations with time reversal the code of the subgroup is zero.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: mag_code
INTEGER, INTENT(OUT) :: group_code, subgroup_code

INTEGER :: group_(58), subgroup_(58), aux_mag_code

CALL set_mag_group_subgroup(group_,subgroup_)
IF (mag_code<=32) THEN
   group_code=mag_code
   subgroup_code=mag_code
ELSEIF (mag_code<=72) THEN
   group_code=mag_code-32
   subgroup_code=0
ELSEIF (mag_code<=122) THEN
   aux_mag_code=mag_code-64
   group_code=group_(aux_mag_code)
   subgroup_code=subgroup_(aux_mag_code)
ELSE
   CALL errore('find_group_subgroup','mag_code not in range',1)
ENDIF

RETURN
END SUBROUTINE find_group_subgroup

SUBROUTINE set_mag_group_subgroup_ext(group_ext,subgroup_ext)
!
!  This routine sets the codes of all the possible combinations of
!  extended group - subgroups codes that define a magnetic point group.
!  Only black and white groups are considered here.
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: group_ext(317), subgroup_ext(317)
INTEGER :: group_ext_(317), subgroup_ext_(317)

DATA group_ext_ / 2,   3,   4,   5,   6,   7,   8,   9,  10,  11,       &
                 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,       &
                 22,  23,  24,  25,  26,  27,  28,  34,  35,  36,       &
                 37,  38,  38,  38,  39,  39,  39,  40,  40,  40,       &
                 41,  41,  41,  42,  42,  42,  43,  43,  43,  44,       &
                 45,  46,  47,  48,  49,  50,  50,  50,  51,  51,       &
                 51,  52,  52,  52,  53,  53,  53,  54,  54,  54,       &
                 55,  55,  55,  56,  56,  56,  57,  57,  57,  58,       &
                 58,  58,  59,  59,  59,  60,  60,  60,  61,  61,       &
                 61,  62,  62,  62,  63,  63,  63,  64,  64,  64,       &
                 65,  65,  65,  66,  66,  66,  67,  67,  67,  68,       &
                 68,  68,  69,  69,  69,  70,  70,  70,  71,  71,       &
                 71,  72,  73,  74,  75,  76,  77,  78,  78,  78,       &
                 79,  79,  79,  80,  80,  80,  81,  81,  81,  82,       &
                 82,  82,  83,  83,  83,  84,  84,  84,  85,  85,       &
                 85,  86,  86,  86,  87,  87,  87,  88,  88,  88,       &
                 89,  89,  89,  90,  90,  90,  91,  91,  91,  92,       &
                 92,  92,  93,  93,  93,  94,  94,  94,  95,  96,       &
                 96,  96,  97,  97,  97,  98,  98,  98,  99,  99,       &
                 99, 100, 100, 100, 100, 100, 100, 100, 101, 101,       &
                101, 101, 101, 101, 101, 102, 102, 102, 102, 102,       &
                102, 102, 103, 103, 103, 103, 103, 103, 103, 104,       &
                104, 104, 104, 104, 104, 104, 105, 105, 105, 105,       &
                105, 105, 105, 106, 106, 106, 107, 107, 107, 108,       &
                108, 108, 108, 108, 108, 108, 109, 109, 109, 109,       &
                109, 109, 109, 110, 110, 110, 110, 110, 110, 110,       &
                111, 111, 111, 111, 111, 111, 111, 112, 112, 112,       &
                113, 113, 113, 114, 114, 114, 115, 115, 115, 116,       &
                116, 116, 117, 117, 117, 118, 118, 118, 119, 119,       &
                119, 120, 120, 120, 121, 121, 121, 122, 122, 122,       &
                123, 123, 123, 124, 125, 126, 127, 128, 129, 130,       &
                131, 133, 134, 135, 136, 136, 136 /

DATA subgroup_ext_ /  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,  &
                      1,    1,   1,   1,   1,   1,   1,   1,   1,   1,  &
                      1,    1,   1,   1,   1,   1,   1,   2,   3,   4,  &
                      33,   2,   3,   4,   2,   5,   6,   3,   7,   8,  &
                       4,   9,  10,   2,  11,  14,   2,  12,  13,  33,  &
                      33,  29,  32,  30,  31,  34,  38,  39,  35,  38,  &
                      40,  36,  38,  41,  37,  44,  45,   4,  15,  16,  &
                       4,  22,  23,   3,  15,  17,   3,  20,  21,   2,  &
                      16,  17,   2,  18,  19,   5,  15,  19,   6,  15,  &
                      18,   7,  16,  21,   8,  16,  20,   9,  17,  23,  &
                      10,  17,  22,   2,  25,  26,   2,  24,  27,  12,  &
                      15,  26,  11,  15,  27,  14,  15,  24,  13,  15,  &
                      25,  33,  33,  29,  32,  30,  31,  34,  58,  59,  &
                      35,  56,  57,  36,  54,  55,  37,  72,  73,   2,  &
                      15,  28,   3,  16,  28,   4,  17,  28,   5,  18,  &
                      28,   6,  19,  28,   7,  20,  28,   8,  21,  28,  &
                       9,  22,  28,  10,  23,  28,  14,  27,  28,  13,  &
                      26,  28,  12,  25,  28,  11,  24,  28,  33,  34,  &
                      82, 124,  35,  83, 125,  36,  84, 126,  37,  95,  &
                     127,  38,  54,  56,  58,  82,  83,  84,  41,  55,  &
                      64,  65,  84,  89,  90,  40,  57,  62,  63,  83,  &
                      87,  88,  39,  59,  60,  61,  82,  85,  86,  42,  &
                      67,  69,  70,  82,  91,  94,  43,  66,  68,  71,  &
                      82,  92,  93,  44,  73,  95,  45,  72,  95,  50,  &
                      78,  96, 100, 103, 112, 113,  51,  79,  97, 100,  &
                     102, 114, 115,  52,  80,  98, 100, 101, 116, 117,  &
                      53,  81,  99, 106, 107, 118, 119,  38,  59, 124,  &
                      39,  58, 124,  38,  57, 125,  40,  56, 125,  38,  &
                      55, 126,  41,  54, 126,  44,  72, 127,  45,  73,  &
                     127,  46,  74, 128,  47,  75, 131,  48,  76, 129,  &
                      49,  77, 130,   2,   3,   4,  33,  29,  30,  31,  &
                      32, 132, 132, 132, 133, 134, 135 /  

group_ext=group_ext_
subgroup_ext=subgroup_ext_

RETURN
END SUBROUTINE set_mag_group_subgroup_ext

SUBROUTINE find_mag_group_code_ext(group_code_ext,subgroup_code_ext, &
                                    mag_code_ext)
!
!   This routine receives the extended codes of a group and of its invariant
!   subgroup and gives the code of the magnetic point group
!   that corresponds to them. The order of the extended magnetic point groups
!   is the following:
!
!       mag_code  mag group name            group             subgroup
!     1   273           C_2(C_1) (2')  F     2   C_2  (2)       1   C_1 (1)    
!     2   274           C_2(C_1) (2')  F     3   C_2  (2)       1   C_1 (1)    
!     3   275           C_2(C_1) (2')  F     4   C_2  (2)       1   C_1 (1)    
!     4   276           C_2(C_1) (2')  F     5   C_2  (2)       1   C_1 (1)    
!     5   277           C_2(C_1) (2')  F     6   C_2  (2)       1   C_1 (1)    
!     6   278           C_2(C_1) (2')  F     7   C_2  (2)       1   C_1 (1)    
!     7   279           C_2(C_1) (2')  F     8   C_2  (2)       1   C_1 (1)    
!     8   280           C_2(C_1) (2')  F     9   C_2  (2)       1   C_1 (1)    
!     9   281           C_2(C_1) (2')  F    10   C_2  (2)       1   C_1 (1)    
!    10   282           C_2(C_1) (2')  F    11   C_2  (2)       1   C_1 (1)    
!    11   283           C_2(C_1) (2')  F    12   C_2  (2)       1   C_1 (1)    
!    12   284           C_2(C_1) (2')  F    13   C_2  (2)       1   C_1 (1)    
!    13   285           C_2(C_1) (2')  F    14   C_2  (2)       1   C_1 (1)    
!    14   286           C_s(C_1) (m')  F    15   C_s (m)        1   C_1 (1)    
!    15   287           C_s(C_1) (m')  F    16   C_s (m)        1   C_1 (1)    
!    16   288           C_s(C_1) (m')  F    17   C_s (m)        1   C_1 (1)    
!    17   289           C_s(C_1) (m')  F    18   C_s (m)        1   C_1 (1)    
!    18   290           C_s(C_1) (m')  F    19   C_s (m)        1   C_1 (1)    
!    19   291           C_s(C_1) (m')  F    20   C_s (m)        1   C_1 (1)    
!    20   292           C_s(C_1) (m')  F    21   C_s (m)        1   C_1 (1)    
!    21   293           C_s(C_1) (m')  F    22   C_s (m)        1   C_1 (1)    
!    22   294           C_s(C_1) (m')  F    23   C_s (m)        1   C_1 (1)    
!    23   295           C_s(C_1) (m')  F    24   C_s (m)        1   C_1 (1)    
!    24   296           C_s(C_1) (m')  F    25   C_s (m)        1   C_1 (1)    
!    25   297           C_s(C_1) (m')  F    26   C_s (m)        1   C_1 (1)    
!    26   298           C_s(C_1) (m')  F    27   C_s (m)        1   C_1 (1)    
!    27   299          C_i(C_1) (-1')  AF   28   C_i (-1)       1   C_1 (1)    
!    28   300           C_4(C_2) (4')  AF   34   C_4 (4)        2   C_2  (2)   
!    29   301           C_4(C_2) (4')  AF   35   C_4 (4)        3   C_2  (2)   
!    30   302           C_4(C_2) (4')  AF   36   C_4 (4)        4   C_2  (2)   
!    31   303           C_6(C_3) (6')  AF   37   C_6 (6)       33   C_3 (3)    
!    32   304        D_2(C_2) (22'2')  F    38   D_2  (222)     2   C_2  (2)   
!    33   305        D_2(C_2) (22'2')  F    38   D_2  (222)     3   C_2  (2)   
!    34   306        D_2(C_2) (22'2')  F    38   D_2  (222)     4   C_2  (2)   
!    35   307        D_2(C_2) (22'2')  F    39   D_2  (222)     2   C_2  (2)   
!    36   308        D_2(C_2) (22'2')  F    39   D_2  (222)     5   C_2  (2)   
!    37   309        D_2(C_2) (22'2')  F    39   D_2  (222)     6   C_2  (2)   
!    38   310        D_2(C_2) (22'2')  F    40   D_2  (222)     3   C_2  (2)   
!    39   311        D_2(C_2) (22'2')  F    40   D_2  (222)     7   C_2  (2)   
!    40   312        D_2(C_2) (22'2')  F    40   D_2  (222)     8   C_2  (2)   
!    41   313        D_2(C_2) (22'2')  F    41   D_2  (222)     4   C_2  (2)   
!    42   314        D_2(C_2) (22'2')  F    41   D_2  (222)     9   C_2  (2)   
!    43   315        D_2(C_2) (22'2')  F    41   D_2  (222)    10   C_2  (2)   
!    44   316        D_2(C_2) (22'2')  F    42   D_2  (222)     2   C_2  (2)   
!    45   317        D_2(C_2) (22'2')  F    42   D_2  (222)    11   C_2  (2)   
!    46   318        D_2(C_2) (22'2')  F    42   D_2  (222)    14   C_2  (2)   
!    47   319        D_2(C_2) (22'2')  F    43   D_2  (222)     2   C_2  (2)   
!    48   320        D_2(C_2) (22'2')  F    43   D_2  (222)    12   C_2  (2)   
!    49   321        D_2(C_2) (22'2')  F    43   D_2  (222)    13   C_2  (2)   
!    50   322          D_3(C_3) (32')  F    44   D_3 (32)      33   C_3 (3)    
!    51   323          D_3(C_3) (32')  F    45   D_3 (32)      33   C_3 (3)    
!    52   324          D_3(C_3) (32')  F    46   D_3 (32)      29   C_3 (3)    
!    53   325          D_3(C_3) (32')  F    47   D_3 (32)      32   C_3 (3)    
!    54   326          D_3(C_3) (32')  F    48   D_3 (32)      30   C_3 (3)    
!    55   327          D_3(C_3) (32')  F    49   D_3 (32)      31   C_3 (3)    
!    56   328        D_4(C_4) (42'2')  F    50   D_4 (422)     34   C_4 (4)    
!    57   329        D_4(D_2) (4'2'2)  AF   50   D_4 (422)     38   D_2  (222) 
!    58   330        D_4(D_2) (4'2'2)  AF   50   D_4 (422)     39   D_2  (222) 
!    59   331        D_4(C_4) (42'2')  F    51   D_4 (422)     35   C_4 (4)    
!    60   332        D_4(D_2) (4'2'2)  AF   51   D_4 (422)     38   D_2  (222) 
!    61   333        D_4(D_2) (4'2'2)  AF   51   D_4 (422)     40   D_2  (222) 
!    62   334        D_4(C_4) (42'2')  F    52   D_4 (422)     36   C_4 (4)    
!    63   335        D_4(D_2) (4'2'2)  AF   52   D_4 (422)     38   D_2  (222) 
!    64   336        D_4(D_2) (4'2'2)  AF   52   D_4 (422)     41   D_2  (222) 
!    65   337        D_6(C_6) (62'2')  F    53   D_6 (622)     37   C_6 (6)    
!    66   338        D_6(D_3) (6'22')  AF   53   D_6 (622)     44   D_3 (32)   
!    67   339        D_6(D_3) (6'22')  AF   53   D_6 (622)     45   D_3 (32)   
!    68   340       C_2v(C_2) (m'm'2)  F    54   C_2v (mm2)     4   C_2  (2)   
!    69   341       C_2v(C_s) (m'm2')  F    54   C_2v (mm2)    15   C_s (m)    
!    70   342       C_2v(C_s) (m'm2')  F    54   C_2v (mm2)    16   C_s (m)    
!    71   343       C_2v(C_2) (m'm'2)  F    55   C_2v (mm2)     4   C_2  (2)   
!    72   344       C_2v(C_s) (m'm2')  F    55   C_2v (mm2)    22   C_s (m)    
!    73   345       C_2v(C_s) (m'm2')  F    55   C_2v (mm2)    23   C_s (m)    
!    74   346       C_2v(C_2) (m'm'2)  F    56   C_2v (mm2)     3   C_2  (2)   
!    75   347       C_2v(C_s) (m'm2')  F    56   C_2v (mm2)    15   C_s (m)    
!    76   348       C_2v(C_s) (m'm2')  F    56   C_2v (mm2)    17   C_s (m)    
!    77   349       C_2v(C_2) (m'm'2)  F    57   C_2v (mm2)     3   C_2  (2)   
!    78   350       C_2v(C_s) (m'm2')  F    57   C_2v (mm2)    20   C_s (m)    
!    79   351       C_2v(C_s) (m'm2')  F    57   C_2v (mm2)    21   C_s (m)    
!    80   352       C_2v(C_2) (m'm'2)  F    58   C_2v (mm2)     2   C_2  (2)   
!    81   353       C_2v(C_s) (m'm2')  F    58   C_2v (mm2)    16   C_s (m)    
!    82   354       C_2v(C_s) (m'm2')  F    58   C_2v (mm2)    17   C_s (m)    
!    83   355       C_2v(C_2) (m'm'2)  F    59   C_2v (mm2)     2   C_2  (2)   
!    84   356       C_2v(C_s) (m'm2')  F    59   C_2v (mm2)    18   C_s (m)    
!    85   357       C_2v(C_s) (m'm2')  F    59   C_2v (mm2)    19   C_s (m)    
!    86   358       C_2v(C_2) (m'm'2)  F    60   C_2v (mm2)     5   C_2  (2)   
!    87   359       C_2v(C_s) (m'm2')  F    60   C_2v (mm2)    15   C_s (m)    
!    88   360       C_2v(C_s) (m'm2')  F    60   C_2v (mm2)    19   C_s (m)    
!    89   361       C_2v(C_2) (m'm'2)  F    61   C_2v (mm2)     6   C_2  (2)   
!    90   362       C_2v(C_s) (m'm2')  F    61   C_2v (mm2)    15   C_s (m)    
!    91   363       C_2v(C_s) (m'm2')  F    61   C_2v (mm2)    18   C_s (m)    
!    92   364       C_2v(C_2) (m'm'2)  F    62   C_2v (mm2)     7   C_2  (2)   
!    93   365       C_2v(C_s) (m'm2')  F    62   C_2v (mm2)    16   C_s (m)    
!    94   366       C_2v(C_s) (m'm2')  F    62   C_2v (mm2)    21   C_s (m)    
!    95   367       C_2v(C_2) (m'm'2)  F    63   C_2v (mm2)     8   C_2  (2)   
!    96   368       C_2v(C_s) (m'm2')  F    63   C_2v (mm2)    16   C_s (m)    
!    97   369       C_2v(C_s) (m'm2')  F    63   C_2v (mm2)    20   C_s (m)    
!    98   370       C_2v(C_2) (m'm'2)  F    64   C_2v (mm2)     9   C_2  (2)   
!    99   371       C_2v(C_s) (m'm2')  F    64   C_2v (mm2)    17   C_s (m)    
!   100   372       C_2v(C_s) (m'm2')  F    64   C_2v (mm2)    23   C_s (m)    
!   101   373       C_2v(C_2) (m'm'2)  F    65   C_2v (mm2)    10   C_2  (2)   
!   102   374       C_2v(C_s) (m'm2')  F    65   C_2v (mm2)    17   C_s (m)    
!   103   375       C_2v(C_s) (m'm2')  F    65   C_2v (mm2)    22   C_s (m)    
!   104   376       C_2v(C_2) (m'm'2)  F    66   C_2v (mm2)     2   C_2  (2)   
!   105   377       C_2v(C_s) (m'm2')  F    66   C_2v (mm2)    25   C_s (m)    
!   106   378       C_2v(C_s) (m'm2')  F    66   C_2v (mm2)    26   C_s (m)    
!   107   379       C_2v(C_2) (m'm'2)  F    67   C_2v (mm2)     2   C_2  (2)   
!   108   380       C_2v(C_s) (m'm2')  F    67   C_2v (mm2)    24   C_s (m)    
!   109   381       C_2v(C_s) (m'm2')  F    67   C_2v (mm2)    27   C_s (m)    
!   110   382       C_2v(C_2) (m'm'2)  F    68   C_2v (mm2)    12   C_2  (2)   
!   111   383       C_2v(C_s) (m'm2')  F    68   C_2v (mm2)    15   C_s (m)    
!   112   384       C_2v(C_s) (m'm2')  F    68   C_2v (mm2)    26   C_s (m)    
!   113   385       C_2v(C_2) (m'm'2)  F    69   C_2v (mm2)    11   C_2  (2)   
!   114   386       C_2v(C_s) (m'm2')  F    69   C_2v (mm2)    15   C_s (m)    
!   115   387       C_2v(C_s) (m'm2')  F    69   C_2v (mm2)    27   C_s (m)    
!   116   388       C_2v(C_2) (m'm'2)  F    70   C_2v (mm2)    14   C_2  (2)   
!   117   389       C_2v(C_s) (m'm2')  F    70   C_2v (mm2)    15   C_s (m)    
!   118   390       C_2v(C_s) (m'm2')  F    70   C_2v (mm2)    24   C_s (m)    
!   119   391       C_2v(C_2) (m'm'2)  F    71   C_2v (mm2)    13   C_2  (2)   
!   120   392       C_2v(C_s) (m'm2')  F    71   C_2v (mm2)    15   C_s (m)    
!   121   393       C_2v(C_s) (m'm2')  F    71   C_2v (mm2)    25   C_s (m)    
!   122   394         C_3v(C_3) (3m')  F    72   C_3v (3m)     33   C_3 (3)    
!   123   395         C_3v(C_3) (3m')  F    73   C_3v (3m)     33   C_3 (3)    
!   124   396         C_3v(C_3) (3m')  F    74   C_3v (3m)     29   C_3 (3)    
!   125   397         C_3v(C_3) (3m')  F    75   C_3v (3m)     32   C_3 (3)    
!   126   398         C_3v(C_3) (3m')  F    76   C_3v (3m)     30   C_3 (3)    
!   127   399         C_3v(C_3) (3m')  F    77   C_3v (3m)     31   C_3 (3)    
!   128   400       C_4v(C_4) (4m'm')  F    78   C_4v (4mm)    34   C_4 (4)    
!   129   401      C_4v(C_2v) (4'm'm)  AF   78   C_4v (4mm)    58   C_2v (mm2) 
!   130   402      C_4v(C_2v) (4'm'm)  AF   78   C_4v (4mm)    59   C_2v (mm2) 
!   131   403       C_4v(C_4) (4m'm')  F    79   C_4v (4mm)    35   C_4 (4)    
!   132   404      C_4v(C_2v) (4'm'm)  AF   79   C_4v (4mm)    56   C_2v (mm2) 
!   133   405      C_4v(C_2v) (4'm'm)  AF   79   C_4v (4mm)    57   C_2v (mm2) 
!   134   406       C_4v(C_4) (4m'm')  F    80   C_4v (4mm)    36   C_4 (4)    
!   135   407      C_4v(C_2v) (4'm'm)  AF   80   C_4v (4mm)    54   C_2v (mm2) 
!   136   408      C_4v(C_2v) (4'm'm)  AF   80   C_4v (4mm)    55   C_2v (mm2) 
!   137   409       C_6v(C_6) (6m'm')  F    81   C_6v (6mm)    37   C_6 (6)    
!   138   410      C_6v(C_3v) (6'm'm)  AF   81   C_6v (6mm)    72   C_3v (3m)  
!   139   411      C_6v(C_3v) (6'm'm)  AF   81   C_6v (6mm)    73   C_3v (3m)  
!   140   412        C_2h(C_2) (2/m')  AF   82   C_2h (2/m)     2   C_2  (2)   
!   141   413        C_2h(C_s) (2'/m)  AF   82   C_2h (2/m)    15   C_s (m)    
!   142   414       C_2h(C_i) (2'/m')  F    82   C_2h (2/m)    28   C_i (-1)   
!   143   415        C_2h(C_2) (2/m')  AF   83   C_2h (2/m)     3   C_2  (2)   
!   144   416        C_2h(C_s) (2'/m)  AF   83   C_2h (2/m)    16   C_s (m)    
!   145   417       C_2h(C_i) (2'/m')  F    83   C_2h (2/m)    28   C_i (-1)   
!   146   418        C_2h(C_2) (2/m')  AF   84   C_2h (2/m)     4   C_2  (2)   
!   147   419        C_2h(C_s) (2'/m)  AF   84   C_2h (2/m)    17   C_s (m)    
!   148   420       C_2h(C_i) (2'/m')  F    84   C_2h (2/m)    28   C_i (-1)   
!   149   421        C_2h(C_2) (2/m')  AF   85   C_2h (2/m)     5   C_2  (2)   
!   150   422        C_2h(C_s) (2'/m)  AF   85   C_2h (2/m)    18   C_s (m)    
!   151   423       C_2h(C_i) (2'/m')  F    85   C_2h (2/m)    28   C_i (-1)   
!   152   424        C_2h(C_2) (2/m')  AF   86   C_2h (2/m)     6   C_2  (2)   
!   153   425        C_2h(C_s) (2'/m)  AF   86   C_2h (2/m)    19   C_s (m)    
!   154   426       C_2h(C_i) (2'/m')  F    86   C_2h (2/m)    28   C_i (-1)   
!   155   427        C_2h(C_2) (2/m')  AF   87   C_2h (2/m)     7   C_2  (2)   
!   156   428        C_2h(C_s) (2'/m)  AF   87   C_2h (2/m)    20   C_s (m)    
!   157   429       C_2h(C_i) (2'/m')  F    87   C_2h (2/m)    28   C_i (-1)   
!   158   430        C_2h(C_2) (2/m')  AF   88   C_2h (2/m)     8   C_2  (2)   
!   159   431        C_2h(C_s) (2'/m)  AF   88   C_2h (2/m)    21   C_s (m)    
!   160   432       C_2h(C_i) (2'/m')  F    88   C_2h (2/m)    28   C_i (-1)   
!   161   433        C_2h(C_2) (2/m')  AF   89   C_2h (2/m)     9   C_2  (2)   
!   162   434        C_2h(C_s) (2'/m)  AF   89   C_2h (2/m)    22   C_s (m)    
!   163   435       C_2h(C_i) (2'/m')  F    89   C_2h (2/m)    28   C_i (-1)   
!   164   436        C_2h(C_2) (2/m')  AF   90   C_2h (2/m)    10   C_2  (2)   
!   165   437        C_2h(C_s) (2'/m)  AF   90   C_2h (2/m)    23   C_s (m)    
!   166   438       C_2h(C_i) (2'/m')  F    90   C_2h (2/m)    28   C_i (-1)   
!   167   439        C_2h(C_2) (2/m')  AF   91   C_2h (2/m)    14   C_2  (2)   
!   168   440        C_2h(C_s) (2'/m)  AF   91   C_2h (2/m)    27   C_s (m)    
!   169   441       C_2h(C_i) (2'/m')  F    91   C_2h (2/m)    28   C_i (-1)   
!   170   442        C_2h(C_2) (2/m')  AF   92   C_2h (2/m)    13   C_2  (2)   
!   171   443        C_2h(C_s) (2'/m)  AF   92   C_2h (2/m)    26   C_s (m)    
!   172   444       C_2h(C_i) (2'/m')  F    92   C_2h (2/m)    28   C_i (-1)   
!   173   445        C_2h(C_2) (2/m')  AF   93   C_2h (2/m)    12   C_2  (2)   
!   174   446        C_2h(C_s) (2'/m)  AF   93   C_2h (2/m)    25   C_s (m)    
!   175   447       C_2h(C_i) (2'/m')  F    93   C_2h (2/m)    28   C_i (-1)   
!   176   448        C_2h(C_2) (2/m')  AF   94   C_2h (2/m)    11   C_2  (2)   
!   177   449        C_2h(C_s) (2'/m)  AF   94   C_2h (2/m)    24   C_s (m)    
!   178   450       C_2h(C_i) (2'/m')  F    94   C_2h (2/m)    28   C_i (-1)   
!   179   451         C_3h(C_3) (-6')  AF   95   C_3h (-6)     33   C_3 (3)    
!   180   452        C_4h(C_4) (4/m')  AF   96   C_4h (4/m)    34   C_4 (4)    
!   181   453       C_4h(C_2h) (4'/m)  AF   96   C_4h (4/m)    82   C_2h (2/m) 
!   182   454       C_4h(S_4) (4'/m')  AF   96   C_4h (4/m)   124   S_4 (-4)   
!   183   455        C_4h(C_4) (4/m')  AF   97   C_4h (4/m)    35   C_4 (4)    
!   184   456       C_4h(C_2h) (4'/m)  AF   97   C_4h (4/m)    83   C_2h (2/m) 
!   185   457       C_4h(S_4) (4'/m')  AF   97   C_4h (4/m)   125   S_4 (-4)   
!   186   458        C_4h(C_4) (4/m')  AF   98   C_4h (4/m)    36   C_4 (4)    
!   187   459       C_4h(C_2h) (4'/m)  AF   98   C_4h (4/m)    84   C_2h (2/m) 
!   188   460       C_4h(S_4) (4'/m')  AF   98   C_4h (4/m)   126   S_4 (-4)   
!   189   461        C_6h(C_6) (6/m')  AF   99   C_6h (6/m)    37   C_6 (6)    
!   190   462      C_6h(C_3h) (6'/m')  AF   99   C_6h (6/m)    95   C_3h (-6)  
!   191   463        C_6h(S_6) (6'/m)  AF   99   C_6h (6/m)   127   S_6 (-3)   
!   192   464      D_2h(D_2) (m'm'm')  AF  100   D_2h (mmm)    38   D_2  (222) 
!   193   465       D_2h(C_2v) (m'mm)  AF  100   D_2h (mmm)    54   C_2v (mm2) 
!   194   466       D_2h(C_2v) (m'mm)  AF  100   D_2h (mmm)    56   C_2v (mm2) 
!   195   467       D_2h(C_2v) (m'mm)  AF  100   D_2h (mmm)    58   C_2v (mm2) 
!   196   468      D_2h(C_2h) (m'm'm)  F   100   D_2h (mmm)    82   C_2h (2/m) 
!   197   469      D_2h(C_2h) (m'm'm)  F   100   D_2h (mmm)    83   C_2h (2/m) 
!   198   470      D_2h(C_2h) (m'm'm)  F   100   D_2h (mmm)    84   C_2h (2/m) 
!   199   471      D_2h(D_2) (m'm'm')  AF  101   D_2h (mmm)    41   D_2  (222) 
!   200   472       D_2h(C_2v) (m'mm)  AF  101   D_2h (mmm)    55   C_2v (mm2) 
!   201   473       D_2h(C_2v) (m'mm)  AF  101   D_2h (mmm)    64   C_2v (mm2) 
!   202   474       D_2h(C_2v) (m'mm)  AF  101   D_2h (mmm)    65   C_2v (mm2) 
!   203   475      D_2h(C_2h) (m'm'm)  F   101   D_2h (mmm)    84   C_2h (2/m) 
!   204   476      D_2h(C_2h) (m'm'm)  F   101   D_2h (mmm)    89   C_2h (2/m) 
!   205   477      D_2h(C_2h) (m'm'm)  F   101   D_2h (mmm)    90   C_2h (2/m) 
!   206   478      D_2h(D_2) (m'm'm')  AF  102   D_2h (mmm)    40   D_2  (222) 
!   207   479       D_2h(C_2v) (m'mm)  AF  102   D_2h (mmm)    57   C_2v (mm2) 
!   208   480       D_2h(C_2v) (m'mm)  AF  102   D_2h (mmm)    62   C_2v (mm2) 
!   209   481       D_2h(C_2v) (m'mm)  AF  102   D_2h (mmm)    63   C_2v (mm2) 
!   210   482      D_2h(C_2h) (m'm'm)  F   102   D_2h (mmm)    83   C_2h (2/m) 
!   211   483      D_2h(C_2h) (m'm'm)  F   102   D_2h (mmm)    87   C_2h (2/m) 
!   212   484      D_2h(C_2h) (m'm'm)  F   102   D_2h (mmm)    88   C_2h (2/m) 
!   213   485      D_2h(D_2) (m'm'm')  AF  103   D_2h (mmm)    39   D_2  (222) 
!   214   486       D_2h(C_2v) (m'mm)  AF  103   D_2h (mmm)    59   C_2v (mm2) 
!   215   487       D_2h(C_2v) (m'mm)  AF  103   D_2h (mmm)    60   C_2v (mm2) 
!   216   488       D_2h(C_2v) (m'mm)  AF  103   D_2h (mmm)    61   C_2v (mm2) 
!   217   489      D_2h(C_2h) (m'm'm)  F   103   D_2h (mmm)    82   C_2h (2/m) 
!   218   490      D_2h(C_2h) (m'm'm)  F   103   D_2h (mmm)    85   C_2h (2/m) 
!   219   491      D_2h(C_2h) (m'm'm)  F   103   D_2h (mmm)    86   C_2h (2/m) 
!   220   492      D_2h(D_2) (m'm'm')  AF  104   D_2h (mmm)    42   D_2  (222) 
!   221   493       D_2h(C_2v) (m'mm)  AF  104   D_2h (mmm)    67   C_2v (mm2) 
!   222   494       D_2h(C_2v) (m'mm)  AF  104   D_2h (mmm)    69   C_2v (mm2) 
!   223   495       D_2h(C_2v) (m'mm)  AF  104   D_2h (mmm)    70   C_2v (mm2) 
!   224   496      D_2h(C_2h) (m'm'm)  F   104   D_2h (mmm)    82   C_2h (2/m) 
!   225   497      D_2h(C_2h) (m'm'm)  F   104   D_2h (mmm)    91   C_2h (2/m) 
!   226   498      D_2h(C_2h) (m'm'm)  F   104   D_2h (mmm)    94   C_2h (2/m) 
!   227   499      D_2h(D_2) (m'm'm')  AF  105   D_2h (mmm)    43   D_2  (222) 
!   228   500       D_2h(C_2v) (m'mm)  AF  105   D_2h (mmm)    66   C_2v (mm2) 
!   229   501       D_2h(C_2v) (m'mm)  AF  105   D_2h (mmm)    68   C_2v (mm2) 
!   230   502       D_2h(C_2v) (m'mm)  AF  105   D_2h (mmm)    71   C_2v (mm2) 
!   231   503      D_2h(C_2h) (m'm'm)  F   105   D_2h (mmm)    82   C_2h (2/m) 
!   232   504      D_2h(C_2h) (m'm'm)  F   105   D_2h (mmm)    92   C_2h (2/m) 
!   233   505      D_2h(C_2h) (m'm'm)  F   105   D_2h (mmm)    93   C_2h (2/m) 
!   234   506      D_3h(D_3) (-6'm'2)  F   106   D_3h (-62m)   44   D_3 (32)   
!   235   507     D_3h(C_3v) (-6'm2')  AF  106   D_3h (-62m)   73   C_3v (3m)  
!   236   508     D_3h(C_3h) (-6m'2')  AF  106   D_3h (-62m)   95   C_3h (-6)  
!   237   509      D_3h(D_3) (-6'm'2)  F   107   D_3h (-62m)   45   D_3 (32)   
!   238   510     D_3h(C_3v) (-6'm2')  AF  107   D_3h (-62m)   72   C_3v (3m)  
!   239   511     D_3h(C_3h) (-6m'2')  AF  107   D_3h (-62m)   95   C_3h (-6)  
!   240   512    D_4h(D_4) (4/m'm'm')  AF  108   D_4h(4/mmm)   50   D_4 (422)  
!   241   513     D_4h(C_4v) (4/m'mm)  AF  108   D_4h(4/mmm)   78   C_4v (4mm) 
!   242   514    D_4h(C_4h) (4/mm'm')  F   108   D_4h(4/mmm)   96   C_4h (4/m) 
!   243   515    D_4h(D_2h) (4'/mm'm)  AF  108   D_4h(4/mmm)  100   D_2h (mmm) 
!   244   516    D_4h(D_2h) (4'/mm'm)  AF  108   D_4h(4/mmm)  103   D_2h (mmm) 
!   245   517   D_4h(D_2d) (4'/m'm'm)  AF  108   D_4h(4/mmm)  112   D_2d (-42m)
!   246   518   D_4h(D_2d) (4'/m'm'm)  AF  108   D_4h(4/mmm)  113   D_2d (-42m)
!   247   519    D_4h(D_4) (4/m'm'm')  AF  109   D_4h(4/mmm)   51   D_4 (422)  
!   248   520     D_4h(C_4v) (4/m'mm)  AF  109   D_4h(4/mmm)   79   C_4v (4mm) 
!   249   521    D_4h(C_4h) (4/mm'm')  F   109   D_4h(4/mmm)   97   C_4h (4/m) 
!   250   522    D_4h(D_2h) (4'/mm'm)  AF  109   D_4h(4/mmm)  100   D_2h (mmm) 
!   251   523    D_4h(D_2h) (4'/mm'm)  AF  109   D_4h(4/mmm)  102   D_2h (mmm) 
!   252   524   D_4h(D_2d) (4'/m'm'm)  AF  109   D_4h(4/mmm)  114   D_2d (-42m)
!   253   525   D_4h(D_2d) (4'/m'm'm)  AF  109   D_4h(4/mmm)  115   D_2d (-42m)
!   254   526    D_4h(D_4) (4/m'm'm')  AF  110   D_4h(4/mmm)   52   D_4 (422)  
!   255   527     D_4h(C_4v) (4/m'mm)  AF  110   D_4h(4/mmm)   80   C_4v (4mm) 
!   256   528    D_4h(C_4h) (4/mm'm')  F   110   D_4h(4/mmm)   98   C_4h (4/m) 
!   257   529    D_4h(D_2h) (4'/mm'm)  AF  110   D_4h(4/mmm)  100   D_2h (mmm) 
!   258   530    D_4h(D_2h) (4'/mm'm)  AF  110   D_4h(4/mmm)  101   D_2h (mmm) 
!   259   531   D_4h(D_2d) (4'/m'm'm)  AF  110   D_4h(4/mmm)  116   D_2d (-42m)
!   260   532   D_4h(D_2d) (4'/m'm'm)  AF  110   D_4h(4/mmm)  117   D_2d (-42m)
!   261   533    D_6h(D_6) (6/m'm'm')  AF  111   D_6h(6/mmm)   53   D_6 (622)  
!   262   534     D_6h(C_6v) (6/m'mm)  AF  111   D_6h(6/mmm)   81   C_6v (6mm) 
!   263   535    D_6h(C_6h) (6/mm'm')  F   111   D_6h(6/mmm)   99   C_6h (6/m) 
!   264   536     D_6h(D_3h) (6/mm'm)  AF  111   D_6h(6/mmm)  106   D_3h (-62m)
!   265   537     D_6h(D_3h) (6/mm'm)  AF  111   D_6h(6/mmm)  107   D_3h (-62m)
!   266   538   D_6h(D_3d) (6'/m'm'm)  AF  111   D_6h(6/mmm)  118   D_3d (-3m) 
!   267   539   D_6h(D_3d) (6'/m'm'm)  AF  111   D_6h(6/mmm)  119   D_3d (-3m) 
!   268   540      D_2d(D_2) (-4'2m')  AF  112   D_2d (-42m)   38   D_2  (222) 
!   269   541     D_2d(C_2v) (-4'2'm)  AF  112   D_2d (-42m)   59   C_2v (mm2) 
!   270   542      D_2d(S_4) (-42'm')  F   112   D_2d (-42m)  124   S_4 (-4)   
!   271   543      D_2d(D_2) (-4'2m')  AF  113   D_2d (-42m)   39   D_2  (222) 
!   272   544     D_2d(C_2v) (-4'2'm)  AF  113   D_2d (-42m)   58   C_2v (mm2) 
!   273   545      D_2d(S_4) (-42'm')  F   113   D_2d (-42m)  124   S_4 (-4)   
!   274   546      D_2d(D_2) (-4'2m')  AF  114   D_2d (-42m)   38   D_2  (222) 
!   275   547     D_2d(C_2v) (-4'2'm)  AF  114   D_2d (-42m)   57   C_2v (mm2) 
!   276   548      D_2d(S_4) (-42'm')  F   114   D_2d (-42m)  125   S_4 (-4)   
!   277   549      D_2d(D_2) (-4'2m')  AF  115   D_2d (-42m)   40   D_2  (222) 
!   278   550     D_2d(C_2v) (-4'2'm)  AF  115   D_2d (-42m)   56   C_2v (mm2) 
!   279   551      D_2d(S_4) (-42'm')  F   115   D_2d (-42m)  125   S_4 (-4)   
!   280   552      D_2d(D_2) (-4'2m')  AF  116   D_2d (-42m)   38   D_2  (222) 
!   281   553     D_2d(C_2v) (-4'2'm)  AF  116   D_2d (-42m)   55   C_2v (mm2) 
!   282   554      D_2d(S_4) (-42'm')  F   116   D_2d (-42m)  126   S_4 (-4)   
!   283   555      D_2d(D_2) (-4'2m')  AF  117   D_2d (-42m)   41   D_2  (222) 
!   284   556     D_2d(C_2v) (-4'2'm)  AF  117   D_2d (-42m)   54   C_2v (mm2) 
!   285   557      D_2d(S_4) (-42'm')  F   117   D_2d (-42m)  126   S_4 (-4)   
!   286   558       D_3d(D_3) (-3'm')  AF  118   D_3d (-3m)    44   D_3 (32)   
!   287   559       D_3d(C_3v) (-3'm)  AF  118   D_3d (-3m)    72   C_3v (3m)  
!   288   560        D_3d(S_6) (-3m')  F   118   D_3d (-3m)   127   S_6 (-3)   
!   289   561       D_3d(D_3) (-3'm')  AF  119   D_3d (-3m)    45   D_3 (32)   
!   290   562       D_3d(C_3v) (-3'm)  AF  119   D_3d (-3m)    73   C_3v (3m)  
!   291   563        D_3d(S_6) (-3m')  F   119   D_3d (-3m)   127   S_6 (-3)   
!   292   564       D_3d(D_3) (-3'm')  AF  120   D_3d (-3m)    46   D_3 (32)   
!   293   565       D_3d(C_3v) (-3'm)  AF  120   D_3d (-3m)    74   C_3v (3m)  
!   294   566        D_3d(S_6) (-3m')  F   120   D_3d (-3m)   128   S_6 (-3)   
!   295   567       D_3d(D_3) (-3'm')  AF  121   D_3d (-3m)    47   D_3 (32)   
!   296   568       D_3d(C_3v) (-3'm)  AF  121   D_3d (-3m)    75   C_3v (3m)  
!   297   569        D_3d(S_6) (-3m')  F   121   D_3d (-3m)   131   S_6 (-3)   
!   298   570       D_3d(D_3) (-3'm')  AF  122   D_3d (-3m)    48   D_3 (32)   
!   299   571       D_3d(C_3v) (-3'm)  AF  122   D_3d (-3m)    76   C_3v (3m)  
!   300   572        D_3d(S_6) (-3m')  F   122   D_3d (-3m)   129   S_6 (-3)   
!   301   573       D_3d(D_3) (-3'm')  AF  123   D_3d (-3m)    49   D_3 (32)   
!   302   574       D_3d(C_3v) (-3'm)  AF  123   D_3d (-3m)    77   C_3v (3m)  
!   303   575        D_3d(S_6) (-3m')  F   123   D_3d (-3m)   130   S_6 (-3)   
!   304   576          S_4(C_2) (-4')  AF  124   S_4 (-4)       2   C_2  (2)   
!   305   577          S_4(C_2) (-4')  AF  125   S_4 (-4)       3   C_2  (2)   
!   306   578          S_4(C_2) (-4')  AF  126   S_4 (-4)       4   C_2  (2)   
!   307   579          S_6(C_3) (-3')  AF  127   S_6 (-3)      33   C_3 (3)    
!   308   580          S_6(C_3) (-3')  AF  128   S_6 (-3)      29   C_3 (3)    
!   309   581          S_6(C_3) (-3')  AF  129   S_6 (-3)      30   C_3 (3)    
!   310   582          S_6(C_3) (-3')  AF  130   S_6 (-3)      31   C_3 (3)    
!   311   583          S_6(C_3) (-3')  AF  131   S_6 (-3)      32   C_3 (3)    
!   312   584            T_h(T) (m'3)  AF  133   T_h (m-3)    132   T    (23)  
!   313   585        T_d(T) (-4'-3m')  AF  134   T_d (-43m)   132   T    (23)  
!   314   586            O(T) (4'32')  AF  135   O   (432)    132   T    (23)  
!   315   587        O_h(T_h) (m-3m')  AF  136   O_h (m-3m)   133   T_h (m-3)  
!   316   588        O_h(T_d) (m'-3m)  AF  136   O_h (m-3m)   134   T_d (-43m) 
!   317   589        O_h(O) (m'-3'm')  AF  136   O_h (m-3m)   135   O   (432)  
!
!   magnetic codes from 1 to 136 are reserved to the gray groups, while
!   the codes from 137 to 272 are reserved to groups with no time reversal
!   operation. In the first case the magnetic code is the one of the group,
!   while in the second case it is 136 + the code of the group.
!
!   F means ferro or ferri-magnetic, AF means antiferromagnetic point
!   groups.
!
!   If the group-subgroup extended codes do not correspond to a magnetic
!   group the returned mag_code_ext is zero.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: group_code_ext, subgroup_code_ext
INTEGER, INTENT(OUT) :: mag_code_ext
INTEGER :: mag_code_tab_ext(136,136)

INTEGER :: group_ext(317), subgroup_ext(317), icode

IF (group_code_ext<1.OR.group_code_ext>136.OR.subgroup_code_ext<0.OR.&
    subgroup_code_ext>136) CALL errore('find_mag_group_code_ext',&
          'group or subgroup codes out of range',1)

CALL set_mag_group_subgroup_ext(group_ext,subgroup_ext)
mag_code_tab_ext=0
mag_code_ext=0
IF (subgroup_code_ext==group_code_ext) THEN
!
!   This is the case of the gray groups, time reversal is an element of
!   the group

   mag_code_ext=group_code_ext
ELSEIF (subgroup_code_ext==0) THEN
!
!  These are the groups without any time reversal operation
!
   mag_code_ext=group_code_ext+136
ELSE
!
!  These are the magnetic groups.
!
   DO icode=1,317
      mag_code_tab_ext(group_ext(icode),subgroup_ext(icode))=272+icode
   ENDDO
   mag_code_ext=mag_code_tab_ext(group_code_ext, subgroup_code_ext)
ENDIF


RETURN
END SUBROUTINE find_mag_group_code_ext

SUBROUTINE find_group_subgroup_ext(group_code_ext,subgroup_code_ext, &
                                    mag_code_ext)
!
!  This routine recieve the extended code of a magnetic group 
!  (between 1 and 589) and sets the extended code of the group and of the 
!  extended code of the invariant subgroup that define the magnetic group. 
!  For gray groups the code of the group and of the subgroups are equal, 
!  for groups that do not contain any operation combined with time reversal 
!  the code of the subgroup is conventionally set to zero.
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: group_code_ext, subgroup_code_ext
INTEGER, INTENT(IN) :: mag_code_ext

INTEGER :: aux_mag_code
INTEGER :: group_ext_(317), subgroup_ext_(317)

CALL set_mag_group_subgroup_ext(group_ext_,subgroup_ext_)
IF (mag_code_ext<=136) THEN
   group_code_ext=mag_code_ext
   subgroup_code_ext=mag_code_ext
ELSEIF (mag_code_ext<=272) THEN
   group_code_ext=mag_code_ext-136
   subgroup_code_ext=0
ELSEIF (mag_code_ext<=589) THEN
   aux_mag_code=mag_code_ext-272
   group_code_ext=group_ext_(aux_mag_code)
   subgroup_code_ext=subgroup_ext_(aux_mag_code)
ELSE
   CALL errore('find_group_subgroup_ext','mag_code_ext not in range',1)
ENDIF

RETURN
END SUBROUTINE find_group_subgroup_ext


FUNCTION mag_group_index_from_ext(mag_code_group_ext)
!
!  This routine receives the extended magnetic group code of a magnetic
!  point group and gives the magnetic_group_code
!
USE point_group, ONLY : group_index_from_ext
IMPLICIT NONE
INTEGER :: mag_group_index_from_ext
INTEGER, INTENT(IN) :: mag_code_group_ext

INTEGER :: group_ext, subgroup_ext, group_code, subgroup_code, mag_code

CALL find_group_subgroup_ext(group_ext, subgroup_ext, mag_code_group_ext)

group_code=group_index_from_ext(group_ext) 
IF (subgroup_ext>0) THEN
   subgroup_code=group_index_from_ext(subgroup_ext) 
ELSE
   subgroup_code=0
ENDIF

CALL find_mag_group_code(group_code,subgroup_code,mag_code)
mag_group_index_from_ext=mag_code

RETURN
END FUNCTION mag_group_index_from_ext


END MODULE magnetic_point_group

