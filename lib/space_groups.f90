!
! Copyright (C) 2014 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE space_groups
!
!  This module contains variables and routines to deal with the space group
!  symmetry.
!  Presently it performs the following tasks:
!  Given a bravais lattice and a set of symmetries with corresponding
!  fractional translations it finds the space group number.
!
!  Given the space group number it provides the factors that must be
!  contained in the fft grid to have a mesh compatible with the
!  fractional translations.
!
!  Given the space group name, it returns the space group number
!  Given the space group number, it returns all the space groups names
!  compatible with that number.
!  In the case of orthorhombic groups it gives the rotation needed to
!  transform the group to the one found in the ITA tables.
!
!  The point group and the Bravais lattice must be compatible, otherwise
!  this routine will not work.
!               
!
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER               :: spg_code            ! space group code
  CHARACTER(LEN=14)     :: spg_name            ! name of the space group
  REAL(DP), ALLOCATABLE :: equivalent_tau(:,:) !

  INTEGER, PARAMETER :: nsg = 1068             ! number of space groups
  INTEGER, PARAMETER :: nsg_2d = 20            ! number of 2d space groups
  CHARACTER(LEN=13) :: space_group_names_2d(nsg_2d)
  INTEGER ::  space_group_numbers_2d(nsg_2d)

  DATA space_group_names_2d  &
       / "p1           ", "p2           ", "pm          s", "pg          s", &
         "cm          s", "p2mm         ", "p2mg         ", "p2gg         ", &
         "c2mm         ", "p4           ", "p4mm         ", "p4gm         ", &
         "p3           ", "p3m1         ", "p31m         ", "p6           ", &
         "p6mm         ", "p1m1         ", "p1g1         ", "c1m1         "  /

  DATA space_group_numbers_2d &
      / 1,               2,             3,             4,            &
        5,               6,             7,             8,            &
        9,              10,            11,            12,            &
       13,              14,            15,            16,            &
       17,               3,             4,             5  /    


  CHARACTER(LEN=8) :: sym_label_sg(173)
  DATA sym_label_sg / 'E', '2z', '2y',  '2x',   '2xy', '2x-y', '4-z', '4z',    &
              '2xz',   '2x-z', '4y', '4-y',   '2yz', '2y-z', '4-x', '4x',     &
              '3-x-y-z', '3-xyz',  '3xy-z', '3x-yz', '3xyz', '3-xy-z',        &
              '3x-y-z', '3-x-yz',  '6z', '6-z', '3z', '3-z', '21-10', '2210', &
              '2010', '2110', &
              'i',   'i2z',   'i2y', 'i2x',  'i2xy', 'i2x-y', 'i4-z', 'i4z',  &
              'i2xz', 'i2x-z', 'i4y', 'i4-y', 'i2yz', 'i2y-z', 'i4-x', 'i4x', &
              'i3-x-y-z', 'i3-xyz', 'i3xy-z', 'i3x-yz', 'i3xyz', 'i3-xy-z',   &
              'i3x-y-z', 'i3-x-yz', 'i6z', 'i6-z', 'i3z', 'i3-z', 'i21-10',   &
              'i2210', 'i2010', 'i2110', '2z1', '2y1', '2x1', '2xy1', '2x-y1', &
              '4-z1', '4-z2', '4-z3', '4z1', '4z2', '4z3', '2xz1','2x-z1', &
              '4y1', '4y2', '4y3', '4-y1', '4-y2', '4-y3', '2yz1', '2y-z1', &
              '4-x1', '4-x2', '4-x3', '4x1', '4x2', '4x3', '3-x-y-z1', &
              '3-x-y-z2', '3-xyz1', '3-xyz2', '3xy-z1', '3xy-z2', '3x-yz1', &
              '3x-yz2', '3xyz1', '3xyz2', '3-xy-z1', '3-xy-z2', '3x-y-z1', &
              '3x-y-z2', '3-x-yz1', '3-x-yz2', '6z1', '6z2', '6z3', '6z4', &
              '6z5', '6-z1', '6-z2', '6-z3', '6-z4',  '6-z5', '3z1', '3z2', &
              '3-z1', '3-z2', '21-101', '22101', '20101', '21101', 'i2za', &
              'i2zb', 'i2zn', 'i2zd', 'i2ya', 'i2yc', 'i2yn', 'i2yd', &
              'i2xb', 'i2xc', 'i2xn', 'i2xd', 'i2xya', 'i2xyc', 'i2xyn', &
              'i2xyd', 'i2x-ya', 'i2x-yc', 'i2x-yn', 'i2x-yd', &
              'i2xza', 'i2xzb', 'i2xzn', 'i2xzd', 'i2x-za', 'i2x-zb', &
              'i2x-zn', 'i2x-zd', 'i2yza', 'i2yzb', 'i2yzn', 'i2yzd', &
              'i2y-za', 'i2y-zb', 'i2y-zn', 'i2y-zd', 'i21-10a', &
              'i21-10c', 'i21-10n', 'i2210a', 'i2210c', &
              'i2210n', 'i2010a', 'i2010c', 'i2010n', &
              'i2110a', 'i2110c', 'i2110n'  /

  INTEGER, PARAMETER :: mnsg=21           ! maximum number of sg names
  CHARACTER(LEN=16) :: sg_names(mnsg,230)
  CHARACTER(LEN=60) :: add_info(mnsg,230)=''

  DATA sg_names(:,1) /                                                   &
       "P1              ", "A1              ", "B1              ",       &
       "C1              ", "F1              ", "I1              ",       &
       "C_1^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,2) /                                                   &
       "P-1             ", "A-1             ", "B-1             ",       &
       "C-1             ", "F-1             ", "I-1             ",       &
       "C_i^1          S", "S_2^1          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,3) /                                                   &
       "P121            ", "P112            ", "P211            ",       &
       "A211            ", "B121            ", "C112            ",       &
       "P2             s", "C_2^1          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,4) /                                                   &
       "P12_11          ", "P112_1          ", "P2_111          ",       &
       "A2_111          ", "B12_11          ", "C112_1          ",       &
       "P2_1           s", "C_2^2          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,5) /                                                   &
       "C121            ", "A121            ", "I121            ",       &
       "A112            ", "B112            ", "I112            ",       &
       "B211            ", "C211            ", "I211            ",       &
       "C2             s", "C_2^3          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,6) /                                                   &
       "P1m1            ", "P11m            ", "Pm11            ",       &
       "Am11            ", "B1m1            ", "C11m            ",       &
       "Pm             s", "C_s^1          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,7) /                                                   &
       "P1c1            ", "P1a1            ", "P1n1            ",       &
       "P11a            ", "P11b            ", "P11m            ",       &
       "Pb11            ", "Pc11            ", "Pn11            ",       &
       "PC             s", "C_s^2          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,8) /                                                   &
       "C1m1            ", "A1m1            ", "I1m1            ",       &
       "A11m            ", "B11m            ", "I11m            ",       &
       "Bm11            ", "Cm11            ", "Im11            ",       &
       "Cm             s", "C_s^3          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,9) /                                                   &
       "C1c1            ", "A1a1            ", "I1a1            ",       &
       "A11a            ", "B11b            ", "I11a            ",       &
       "Bb11            ", "Cc11            ", "Ib11            ",       &
       "Cc             s", "C_s^4          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,10) /                                                  &
       "P12/m1          ", "P112/m          ", "P2/m11          ",       &
       "A2/m11          ", "B12/m1          ", "C112/m          ",       &
       "P2/m           s", "C_2h^1         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,11) /                                                  &
       "P12_1/m1        ", "P112_1/m        ", "P2_1/m11        ",       & 
       "A2_1/m11        ", "B12_1/m1        ", "C112_1/m        ",       &
       "P2_1/m         s", "C_2h^2         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,12) /                                                  &
       "C12/m1          ", "A12/m1          ", "I12/m1          ",       &
       "A112/m          ", "B112/m          ", "I112/m          ",       &
       "B2/m11          ", "C2/m11          ", "I2/m11          ",       &
       "C2/m           s", "C_2h^3         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,13) /                                                  &
       "P12/c1          ", "P12/a1          ", "P12/n1          ",       &
       "P112/a          ", "P112/b          ", "P112/n          ",       &
       "P2/b11          ", "P2/c11          ", "P2/n11          ",       &
       "P2/c           s", "C_2h^4         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,14) /                                                  &
       "P12_1/c1        ", "P12_1/a1        ", "P12_1/n1        ",       &
       "P112_1/a        ", "P112_1/b        ", "P112_1/n        ",       &
       "P2_1/b11        ", "P2_1/c11        ", "P2_1/n11        ",       &
       "P2_1/c         s", "C_2h^5         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,15) /                                                  &
       "C12/c1          ", "A12/n1          ", "I12/a1          ",       &
       "A12/a1          ", "C12/n1          ", "I12/c1          ",       &
       "A112/a          ", "B112/n          ", "I112/b          ",       &
       "B112/b          ", "A112/n          ", "I112/a          ",       &
       "B2/b11          ", "C2/n11          ", "I2/c11          ",       &
       "C2/c11          ", "B2/n11          ", "I2/b11          ",       &
       "C2/c           s", "C_2h^6         S", "                " /

  DATA sg_names(:,16) /                                                  &
       "P222            ", "D_2^1          S", "V^1            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,17) /                                                  &
       "P222_1          ", "P2_122          ", "P22_12          ",       &
       "D_2^2          S", "V^2            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,18) /                                                  &
       "P2_12_12        ", "P2_122_1        ", "P22_12_1        ",       &
       "D_2^3          S", "V^3            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,19) /                                                  &
       "P2_12_12_1      ", "D_2^4          S", "V^4            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,20) /                                                  &
       "C222_1          ", "A2_122          ", "B22_12          ",       &
       "D_2^5          S", "V^5            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,21) /                                                  &
       "C222            ", "A222            ", "B222            ",       &
       "D_2^6          S", "V^6            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,22) /                                                  &
       "F222            ", "D_2^7          S", "V^7            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,23) /                                                  &
       "I222            ", "D_2^8          S", "V^8            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,24) /                                                  &
       "I2_12_12_1      ", "D_2^9          S", "V^9            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,25) /                                                  &
       "Pmm2            ", "P2mm            ", "Pm2m            ",       &
       "Pmm           *s", "C_2v^1         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,26) /                                                  &
       "Pmc2_1          ", "Pcm2_1          ", "P2_1ma          ",       &
       "P2_1am          ", "Pb2_1m          ", "Pm2_1b          ",       &
       "Pmc           *s", "C_2v^2         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,27) /                                                  &
       "Pcc2            ", "P2aa            ", "Pb2b            ",       &
       "Pcc           *s", "C_2v^3         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,28) /                                                  &
       "Pma2            ", "Pbm2            ", "P2mb            ",       &
       "P2cm            ", "Pc2m            ", "Pm2a            ",       &
       "Pma           *s", "C_2v^4         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,29) /                                                  &
       "Pca2_1          ", "Pbc2_1          ", "P2_1ab          ",       &
       "P2_1ca          ", "Pc2_1b          ", "Pb2_1a          ",       &
       "Pca           *s", "C_2v^5         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,30) /                                                  &
       "Pnc2            ", "Pcn2            ", "P2na            ",       &
       "P2an            ", "Pb2n            ", "Pn2b            ",       &
       "Pnc           *s", "C_2v^6         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,31) /                                                  &
       "Pmn2_1          ", "Pnm2_1          ", "P2_1mn          ",       &
       "P2_1nm          ", "Pn2_1m          ", "Pm2_1n          ",       &
       "Pmn           *s", "C_2v^7         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,32) /                                                  &
       "Pba2            ", "Pc2a            ", "P2cb            ",       &
       "Pba           *s", "C_2v^8         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,33) /                                                  &
       "Pna2_1          ", "Pbn2_1          ", "P2_1nb          ",       &
       "P2_1cn          ", "Pc2_1n          ", "Pn2_1a          ",       &
       "Pna           *s", "C_2v^9         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /
   
  DATA sg_names(:,34) /                                                  &
       "Pnn2            ", "P2nn            ", "Pn2n            ",       &
       "Pnn           *s", "C_2v^10        S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,35) /                                                  &
       "Cmm2            ", "A2mm            ", "Bm2m            ",       &
       "Cmm           *s", "C_2v^11        S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,36) /                                                  &
       "Cmc2_1          ", "Ccm2_1          ", "A2_1ma          ",       &
       "A2_1am          ", "Bb2_1m          ", "Bm2_1b          ",       &
       "Cmc           *s", "C_2v^12        S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,37) /                                                  &
       "Ccc2            ", "A2aa            ", "Bb2b            ",       &
       "Ccc           *s", "C_2v^13        S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,38) /                                                  &
       "Amm2            ", "Bmm2            ", "B2mm            ",       &
       "C2mm            ", "Cm2m            ", "Am2m            ",       &
       "Amm           *s", "C_2v^14        S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,39) /                                                  &
       "Aem2            ", "Bme2            ", "B2em            ",       & 
       "C2me            ", "Cm2e            ", "Ae2m            ",       & 
       "Abm2          * ", "Bma2          * ", "B2am          * ",       & 
       "C2ma          * ", "Cm2a          * ", "Ab2m          * ",       & 
       "Acm2          * ", "Bmc2          * ", "B2cm          * ",       & 
       "C2mb          * ", "Cm2b          * ", "Ac2m          * ",       & 
       "Abm           *s", "C_2v^15        S", "                " /

  DATA sg_names(:,40) /                                                  &
       "Ama2            ", "Bbm2            ", "B2mb            ",       & 
       "C2cm            ", "Cc2m            ", "Am2a            ",       & 
       "Ama           *s", "C_2v^16        S", "                ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,41) /                                                  &
       "Aea2            ", "Bbe2            ", "B2eb            ",       & 
       "C2ce            ", "Cc2e            ", "Ae2a            ",       & 
       "Aba2          * ", "Bba2          * ", "B2ab          * ",       & 
       "C2ca          * ", "Cc2a          * ", "Ab2a          * ",       & 
       "Aca2          * ", "Bbc2          * ", "B2cb          * ",       & 
       "C2cb          * ", "Cc2b          * ", "Ac2a          * ",       & 
       "Aba           *s", "C_2v^17        S", "                " /

  DATA sg_names(:,42) /                                                  &
       "Fmm2            ", "F2mm            ", "Fm2m            ",       & 
       "Fmm           *s", "C_2v^18        S", "                ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,43) /                                                  &
       "Fdd2            ", "F2dd            ", "Fd2d            ",       & 
       "Fdd           *s", "C_2v^19        S", "                ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,44) /                                                  &
       "Imm2            ", "I2mm            ", "Im2m            ",       & 
       "Imm           *s", "C_2v^20        S", "                ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,45) /                                                  &
       "Iba2            ", "I2cb            ", "Ic2a            ",       & 
       "Iba           *s", "C_2v^21        S", "                ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,46) /                                                  &
       "Ima2            ", "Ibm2            ", "I2mb            ",       & 
       "I2cm            ", "Ic2m            ", "Im2a            ",       & 
       "Ima           *s", "C_2v^22        S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,47) /                                                  &
       "Pmmm            ", "P2/m2/m2/m      ", "D_2h^1         S",       & 
       "V_h^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,48) /                                                  &
       "Pnnn            ", "P2/n2/n2/n      ", "D_2h^2         S",       & 
       "V_h^2          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,49) /                                                  &
       "Pccm            ", "Pmaa            ", "Pbmb            ",       & 
       "P2/c2/c2/m      ", "D_2h^3         S", "V_h^3          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,50) /                                                  &
       "Pban            ", "Pncb            ", "Pcna            ",       & 
       "P2/b2/a2/n      ", "D_2h^4         S", "V_h^4          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,51) /                                                  &
       "Pmma            ", "Pmmb            ", "Pbmm            ",       & 
       "Pcmm            ", "Pmcm            ", "Pmam            ",       & 
       "P2_1/m2/m2/a    ", "D_2h^5         S", "V_h^5          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,52) /                                                  &
       "Pnna            ", "Pnnb            ", "Pbnn            ",       & 
       "Pcnn            ", "Pncn            ", "Pnan            ",       & 
       "P2/n2_1/n2/a    ", "D_2h^6         S", "V_h^6          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,53) /                                                  &
       "Pmna            ", "Pnmb            ", "Pbmn            ",       & 
       "Pcnm            ", "Pncm            ", "Pman            ",       & 
       "2/m2/n2_1/a     ", "D_2h^7         S", "V_h^7          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,54) /                                                  &
       "Pcca            ", "Pccb            ", "Pbaa            ",       & 
       "Pcaa            ", "Pbcb            ", "Pbab            ",       & 
       "P2_1/c2/c2/a    ", "D_2h^8         S", "V_h^8          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,55) /                                                  &
       "Pbam            ", "Pmcb            ", "Pcma            ",       & 
       "P2_1/b2_1/a2/m  ", "D_2h^9         S", "V_h^9          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,56) /                                                  &
       "Pccn            ", "Pnaa            ", "Pbnb            ",       & 
       "P2_1/c2_1/c2/n  ", "D_2h^10        S", "V_h^10         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,57) /                                                  &
       "Pbcm            ", "Pcam            ", "Pmca            ",       & 
       "Pmab            ", "Pbma            ", "Pcmb            ",       & 
       "P2/b2_1/c2_1/m  ", "D_2h^11        S", "V_h^11         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,58) /                                                  &
       "Pnnm            ", "Pmnn            ", "Pnmn            ",       & 
       "P2_1/n2_1/n2/m  ", "D_2h^12        S", "V_h^12         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,59) /                                                  &
       "Pmmn            ", "Pnmm            ", "Pmnm            ",       & 
       "P2_1/m2_1/m2/n  ", "D_2h^13        S", "V_h^13         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,60) /                                                  &
       "Pbcn            ", "Pcan            ", "Pnca            ",       & 
       "Pnab            ", "Pbna            ", "Pcnb            ",       & 
       "P2_1/b2/c2_1/n  ", "D_2h^14        S", "V_h^14         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,61) /                                                  &
       "Pbca            ", "Pcab            ", "Pbca            ",       & 
       "Pcab            ", "Pbca            ", "Pcab            ",       & 
       "P2_1/b2_1/c2_1/a", "D_2h^14        S", "V_h^14         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,62) /                                                  &
       "Pnma            ", "Pmnb            ", "Pbnm            ",       & 
       "Pcmn            ", "Pmcn            ", "Pnam            ",       & 
       "P2_1/n2_1/m2_1/a", "D_2h^15        S", "V_h^15         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,63) /                                                  &
       "Cmcm            ", "Ccmm            ", "Amma            ",       & 
       "Amam            ", "Bbmm            ", "Bmmb            ",       & 
       "C2/m2/c2_1/m    ", "D_2h^16        S", "V_h^16         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,64) /                                                  &
       "Cmce            ", "Ccme            ", "Aema            ",       & 
       "Aeam            ", "Bbem            ", "Bmeb            ",       & 
       "Cmca          * ", "Ccma          * ", "Abma          * ",       & 
       "Abam          * ", "Bbam          * ", "Bmab          * ",       & 
       "Cmcb          * ", "Ccmb          * ", "Acma          * ",       & 
       "Acam          * ", "Bbcm          * ", "Bmcb          * ",       & 
       "C2/m2/c2_1/e    ", "D_2h^18        S", "V_h^18         S"  /

  DATA sg_names(:,65) /                                                  &
       "Cmmm            ", "Ammm            ", "Bmmm            ",       & 
       "C2/m2/m2/m      ", "D_2h^19        S", "V_h^19         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,66) /                                                  &
       "Cccm            ", "Amaa            ", "Bbmb            ",       & 
       "C2/c2/c2/m      ", "D_2h^20        S", "V_h^20         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,67) /                                                  &
       "Cmme            ", "Aemm            ", "Bmem            ",       & 
       "Cmma          * ", "Abmm          * ", "Bmam          * ",       & 
       "Cmmb          * ", "Acmm          * ", "Bmcm          * ",       & 
       "C2/m2/m2/e      ", "D_2h^21        S", "V_h^21         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,68) /                                                  &
       "Ccce            ", "Aeaa            ", "Bbeb            ",       & 
       "Ccca          * ", "Abaa          * ", "Bbab          * ",       & 
       "Cccb          * ", "Acaa          * ", "Bbcb          * ",       & 
       "C2/c2/c2/e      ", "D_2h^22        S", "V_h^22         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,69) /                                                  &
       "Fmmm            ", "F2/m2/m2/m      ", "D_2h^23        S",       & 
       "V_h^23         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,70) /                                                  &
       "Fddd            ", "F2/d2/d2/d      ", "D_2h^24        S",       & 
       "V_h^24         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,71) /                                                  &
       "Immm            ", "I2/m2/m2/m      ", "D_2h^25        S",       &
       "V_h^25         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,72) /                                                  &
       "Ibam            ", "Imcb            ", "Icma            ",       &
       "I2/b2/a2/m      ", "D_2h^26        S", "V_h^26         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,73) /                                                  &
       "Ibca            ", "Icab            ", "I2_1/b2_1/c2_1/a",       & 
       "I2/b2/c2/a    # ", "D_2h^27         ", "V_h^27          ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,74) /                                                  &
       "Imma            ", "Immb            ", "Ibmm            ",       & 
       "Icmm            ", "Imcm            ", "Imam            ",       & 
       "I2_1/m2_1/m2_1/a", "I2/m2/m2/a    # ", "D_2h^28        S",       & 
       "V_h^28         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,75) /                                                  &
       "P4              ", "C4              ", "C_4^1          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,76) /                                                  &
       "P4_1            ", "C4_1            ", "C_4^2          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,77) /                                                  &
       "P4_2            ", "C4_2            ", "C_4^3          S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,78) /                                                  &
       "P4_3            ", "C4_3            ", "C_4^4          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,79) /                                                  &
       "I4              ", "F4              ", "C_4^5          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,80) /                                                  &
       "I4_1            ", "F4_1            ", "C_4^6          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,81) /                                                  &
       "P-4             ", "C-4             ", "S_4^1          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,82) /                                                  &
       "I-4             ", "F-4             ", "S_4^2          S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,83) /                                                  &
       "P4/m            ", "C4/m            ", "C_4h^1         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,84) /                                                  &
       "P4_2/m          ", "C4_2/m          ", "C_4h^2         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,85) /                                                  &
       "P4/n            ", "C4/a            ", "C_4h^3         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,86) /                                                  &
       "P4_2/n          ", "C4_2/a          ", "C_4h^4         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,87) /                                                  &
       "I4/m            ", "F4/m            ", "C_4h^5         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,88) /                                                  &
       "I4_1/a          ", "F4_1/d          ", "C_4h^6         S",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,89) /                                                  &
       "P422            ", "C422            ", "P42           *s",       & 
       "D_4^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,90) /                                                  &
       "P42_12          ", "C42_12          ", "P42_1         *s",       & 
       "D_4^2          S", "                ", "                ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,91) /                                                  &
       "P4_122          ", "C4_122          ", "P4_12         *s",       & 
       "D_4^3          S", "                ", "                ",       & 
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,92) /                                                  &
       "P4_12_12        ", "C4_12_12        ", "P4_12_1       *s",       & 
       "D_4^4          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,93) /                                                  &
       "P4_222          ", "C4_222          ", "P4_22         *s",       &
       "D_4^5          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,94) /                                                  &
       "P4_22_12        ", "C4_22_12        ", "P4_22_1       *s",       & 
       "D_4^6          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,95) /                                                  &
       "P4_322          ", "C4_322          ", "P4_32         *s",       & 
       "D_4^7          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,96) /                                                  &
       "P4_32_12        ", "C4_32_12        ", "P4_32_1       *s",       &
       "D_4^8          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,97) /                                                  &
       "I422            ", "F422            ", "I42           *s",       & 
       "D_4^9          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  DATA sg_names(:,98) /                                                  &
       "I4_122          ", "F4_122          ", "I4_12         *s",       &
       "D_4^10         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,99) /                                                  &
       "P4mm            ", "C4mm            ", "C_4v^1         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,100) /                                                 &
       "P4bm            ", "C4mg1           ", "C_4v^2         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,101) /                                                 &
       "P4_2cm          ", "C4_2mc          ", "P4cm          *s",       &
       "C_4v^3         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,102) /                                                 &
       "P4_2nm          ", "C4_2mg2         ", "P4nm          *s",       &
       "C_4v^4         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,103) /                                                 &
       "P4cc            ", "C4cc            ", "C_4v^5         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,104) /                                                 &
       "P4nc            ", "C4cg2           ", "C_4v^6         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,105) /                                                 &
       "P4_2mc          ", "C4_2cm          ", "C_4v^7         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,106) /                                                 &
       "P4_2bc          ", "C4_2cg1         ", "C_4v^8         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,107) /                                                 &
       "I4mm            ", "F4mm            ", "C_4v^9         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,108) /                                                 &
       "I4cm            ", "F4mc            ", "C_4v^10        S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,109) /                                                 &
       "I4_1md          ", "F4_1md          ", "I4md          *s",       &
       "C_4v^11        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,110) /                                                 &
       "I4_1cd          ", "F4_1dc          ", "I4cd          *s",       &
       "C_4v^12        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,111) /                                                 &
       "P-42m           ", "C-4m2           ", "D_2d^1         S",       &
       "V_d^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,112) /                                                 &
       "P-42c           ", "C-4c2           ", "D_2d^2         S",       &
       "V_d^2          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,113) /                                                 &
       "P-42_1m         ", "C-4m2_1         ", "D_2d^3         S",       &
       "V_d^3          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,114) /                                                 &
       "P-42_1c         ", "C-4c2_1         ", "D_2d^4         S",       &
       "V_d^4          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,115) /                                                 &
       "P-4m2           ", "C-42m           ", "D_2d^5         S",       &
       "V_d^5          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,116) /                                                 &
       "P-4c2           ", "C-42c           ", "D_2d^6         S",       &
       "V_d^6          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,117) /                                                 &
       "P-4b2           ", "C-42g1          ", "D_2d^7         S",       &
       "V_d^7          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,118) /                                                 &
       "P-4n2           ", "C-42g2          ", "D_2d^8         S",       &
       "V_d^8          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,119) /                                                 &
       "I-4m2           ", "F-42m           ", "D_2d^9         S",       &
       "V_d^9          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,120) /                                                 &
       "I-4c2           ", "F-42c           ", "D_2d^10        S",       &
       "V_d^10         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,121) /                                                 &
       "I-42m           ", "F-4m2           ", "D_2d^11        S",       &
       "V_d^11         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,122) /                                                 &
       "I-42d           ", "F-4d2           ", "D_2d^12        S",       &
       "V_d^12         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,123) /                                                 &
       "P4/mmm          ", "C4/mmm          ", "D_4h^1         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,124) /                                                 &
       "P4/mcc          ", "C4/mcc          ", "D_4h^2         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,125) /                                                 &
       "P4/nbm          ", "C4/amg1         ", "D_4h^3         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,126) /                                                 &
       "P4/nnc          ", "C4/acg2         ", "D_4h^4         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,127) /                                                 &
       "P4/mbm          ", "C4/mmg1         ", "D_4h^5         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,128) /                                                 &
       "P4/mnc          ", "C4/mcg2         ", "D_4h^6         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,129) /                                                 &
       "P4/nmm          ", "C4/amm          ", "D_4h^7         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,130) /                                                 &
       "P4/ncc          ", "C4/acc          ", "D_4h^8         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,131) /                                                 &
       "P4_2/mmc        ", "C4_2/mcm        ", "P4/mmc        *s",       &
       "D_4h^9         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,132) /                                                 &
       "P4_2/mcm        ", "C4_2/mmc        ", "P4/mcm        *s",       &
       "D_4h^10        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,133) /                                                 &
       "P4_2/nbc        ", "C4_2/acg1       ", "P4/nbc        *s",       &
       "D_4h^11        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,134) /                                                 &
       "P4_2/nnm        ", "C4_2/amg2       ", "P4/nbc        *s",       &
       "D_4h^12        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,135) /                                                 &
       "P4_2/mbc        ", "C4_2/mcg1       ", "P4/nbc        *s",       &
       "D_4h^13        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,136) /                                                 &
       "P4_2/mnm        ", "C4_2/mmg2       ", "P4/mnm        *s",       &
       "D_4h^14        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,137) /                                                 &
       "P4_2/nmc        ", "C4_2/acm        ", "P4/mnm        *s",       &
       "D_4h^15        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,138) /                                                 &
       "P4_2/ncm        ", "C4_2/amc        ", "P4/ncm        *s",       &
       "D_4h^16        S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,139) /                                                 &
       "I4/mmm          ", "F4/mmm          ", "D_4h^17        S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,140) /                                                 &
       "I4/mcm          ", "F4/mmc          ", "D_4h^18        S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,141) /                                                 &
       "I4_1/amd        ", "F4_1/ddm        ", "D_4h^19        S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,142) /                                                 &
       "I4_1/acd        ", "F4_1/ddc        ", "D_4h^20        S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,143) /                                                 &
       "P3              ", "H3            * ", "C3            * ",       &
       "C_3^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,144) /                                                 &
       "P3_1            ", "H3_1          * ", "C3_1          * ",       &
       "C_3^2          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,145) /                                                 &
       "P3_2            ", "H3_2          * ", "C3_2          * ",       &
       "C_3^3          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,146) /                                                 &
       "R3              ", "C_3^4          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,147) /                                                 &
       "P-3             ", "H-3           * ", "C-3           * ",       &
       "S_6^1           ", "C_3i^1         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,148) /                                                 &
       "R-3             ", "S_6^2           ", "C_3i^2         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,149) /                                                 &
       "P312            ", "H321          * ", "H32           *s",       &
       "D_3^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,150) /                                                 &
       "P321            ", "H312          * ", "C32           *s",       &
       "D_3^2          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,151) /                                                 &
       "P3_112          ", "H3_121        * ", "H3_12         *s",       &
       "D_3^3          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,152) /                                                 &
       "P3_121          ", "H3_112        * ", "C3_12         *s",       &
       "D_3^4          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,153) /                                                 &
       "P3_212          ", "H3_221        * ", "H3_22         *s",       &
       "D_3^5          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,154) /                                                 &
       "P3_221          ", "H3_212        * ", "C3_22         *s",       &
       "D_3^6          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,155) /                                                 &
       "R32             ", "D_3^7          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,156) /                                                 &
       "P3m1            ", "H31m            ", "C3m           *s",       &
       "C3m1          * ", "C_3v^1         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,157) /                                                 &
       "P31m            ", "H3m1            ", "H3m           *s",       &
       "C_3v^2         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,158) /                                                 &
       "P3c1            ", "H31c            ", "C3c           *s",       &
       "C3c1          * ", "C_3v^3         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,159) /                                                 &
       "P31c            ", "H3c1            ", "H3c           *s",       &
       "C_3v^4         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,160) /                                                 &
       "R3m             ", "C_3v^5         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,161) /                                                 &
       "R3c             ", "C_3v^6         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,162) /                                                 &
       "P-31m           ", "H-3m1           ", "H-3m          *s",       &
       "H-32/m1       * ", "D_3d^1         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,163) /                                                 &
       "P-31c           ", "H-3c1           ", "H-3c          *s",       &
       "H-32/c1       * ", "D_3d^2         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,164) /                                                 &
       "P-3m1           ", "H-31m         * ", "C-3m          *s",       &
       "C-32/m1       * ", "D_3d^3         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,165) /                                                 &
       "P-3c1           ", "H-31c         * ", "C-3c          *s",       &
       "C-32/c1       * ", "D_3d^4         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,166) /                                                 &
       "R-3m            ", "R-32/m          ", "D_3d^5         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,167) /                                                 &
       "R-3c            ", "R-32/c          ", "D_3d^6         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,168) /                                                 &
       "P6              ", "H6              ", "C6            * ",       &
       "C_6^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,169) /                                                 &
       "P6_1            ", "H6_1            ", "C6_1          * ",       &
       "C_6^2          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,170) /                                                 &
       "P6_5            ", "H6_5            ", "C6_5          * ",       &
       "C_6^3          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,171) /                                                 &
       "P6_2            ", "H6_2            ", "C6_2          * ",       &
       "C_6^4          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,172) /                                                 &
       "P6_4            ", "H6_4            ", "C6_4          * ",       &
       "C_6^5          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,173) /                                                 &
       "P6_3            ", "H6_3            ", "C6_3          * ",       &
       "C_6^6          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,174) /                                                 &
       "P-6             ", "H-6             ", "C-6           * ",       &
       "C_3h^1         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,175) /                                                 &
       "P6/m            ", "H6/m            ", "C6/m          * ",       &
       "C_6h^1         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,176) /                                                 &
       "P6_3/m          ", "H6_3/m          ", "C6_3/m        * ",       &
       "C_6h^2         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,177) /                                                 &
       "P622            ", "H622            ", "C622          * ",       &
       "C62           *s", "D_6^1          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,178) /                                                 &
       "P6_122          ", "H6_122          ", "C6_122        * ",       &
       "C6_12         *s", "D_6^2          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,179) /                                                 &
       "P6_522          ", "H6_522          ", "C6_522        * ",       &
       "C6_52         *s", "D_6^3          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,180) /                                                 &
       "P6_222          ", "H6_222          ", "C6_222        * ",       &
       "C6_22         *s", "D_6^4          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,181) /                                                 &
       "P6_422          ", "H6_422          ", "C6_422        * ",       &
       "C6_42         *s", "D_6^5          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,182) /                                                 &
       "P6_322          ", "H6_322          ", "C6_322        * ",       &
       "C6_32         *s", "D_6^6          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,183) /                                                 &
       "P6mm            ", "H6mm            ", "C6mm          * ",       &
       "C_6v^1         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,184) /                                                 &
       "P6cc            ", "H6cc            ", "C6cc          * ",       &
       "C_6v^2         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,185) /                                                 &
       "P6_3cm          ", "H6_3mc          ", "C6cm          *s",       &
       "C6_3cm        * ", "C_6v^3         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,186) /                                                 &
       "P6_3mc          ", "H6_3cm          ", "C6mc          *s",       &
       "C6_3mc        * ", "C_6v^4         S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,187) /                                                 &
       "P-6m2           ", "H-62m           ", "C-6m2         * ",       &
       "D_3h^1         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,188) /                                                 &
       "P-6c2           ", "H-62c           ", "C-6c2         * ",       &
       "D_3h^2         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,189) /                                                 &
       "P-62m           ", "H-6m2           ", "D_3h^3         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,190) /                                                 &
       "P-62c           ", "H-6c2           ", "D_3h^4         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,191) /                                                 &
       "P6/mmm          ", "H6/mmm          ", "P6/m2/m2/m      ",       &
       "C6/mmm        *s", "C6/m2/m2/m    * ", "D_6h^1         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,192) /                                                 &
       "P6/mcc          ", "H6/mcc          ", "P6/m2/c2/c      ",       &
       "C6/mcc        *s", "C6/m2/c2/c    * ", "D_6h^2         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,193) /                                                 &
       "P6_3/mcm        ", "H6_3/mmc        ", "P6_3/m2/c2/m    ",       &
       "C6_3/mcm      *s", "C6_3/m2/c2/m  * ", "D_6h^3         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,194) /                                                 &
       "P6_3/mmc        ", "H6_3/mcm        ", "P6_3/m2/m2/c    ",       &
       "C6_3/mmc      *s", "C6_3/m2/m2/c  * ", "D_6h^4         S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,195) /                                                 &
       "P23             ", "T^1            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,196) /                                                 &
       "F23             ", "T^2            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,197) /                                                 &
       "I23             ", "T^3            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,198) /                                                 &
       "P2_13           ", "T^4            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,199) /                                                 &
       "I2_13           ", "T^5            S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,200) /                                                 &
       "Pm-3            ", "P2/m-3          ", "Pm3           *s",       &
       "T_h^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,201) /                                                 &
       "Pn-3            ", "P2/n-3          ", "Pn3           *s",       &
       "T_h^2          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,202) /                                                 &
       "Fm-3            ", "F2/m-3          ", "Fm3           *s",       &
       "T_h^3          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,203) /                                                 &
       "Fd-3            ", "F2/d-3          ", "Fd3           *s",       &
       "T_h^4          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,204) /                                                 &
       "Im-3            ", "I2/m-3          ", "Im3           *s",       &
       "T_h^5          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,205) /                                                 &
       "Pa-3            ", "P2_1/a-3        ", "Pa3           *s",       &
       "T_h^6          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,206) /                                                 &
       "Ia-3            ", "I2_1/a-3        ", "Ia3           *s",       &
       "T_h^7          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,207) /                                                 &
       "P432            ", "P43           *s", "O^1            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,208) /                                                 &
       "P4_232          ", "P4_23         *s", "O^2            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,209) /                                                 &
       "F432            ", "F43           *s", "O^3            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,210) /                                                 &
       "F4_132          ", "F4_13         *s", "O^4            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,211) /                                                 &
       "I432            ", "I43           *s", "O^5            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,212) /                                                 &
       "I4_332          ", "P4_33         *s", "O^6            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,213) /                                                 &
       "P4_132          ", "P4_13         *s", "O^7            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,214) /                                                 &
       "I4_132          ", "I4_13         *s", "O^8            S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,215) /                                                 &
       "P-43m           ", "T_d^1          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,216) /                                                 &
       "F-43m           ", "T_d^2          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,217) /                                                 &
       "I-43m           ", "T_d^3          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,218) /                                                 &
       "P-43n           ", "T_d^4          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,219) /                                                 &
       "F-43c           ", "F-43c           ", "T_d^5          S",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,220) /                                                 &
       "I-43d           ", "T_d^6          S", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,221) /                                                 &
       "Pm-3m           ", "P4/m-32/m       ", "Pm3m          *s",       &
       "O_h^1          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,222) /                                                 &
       "Pn-3n           ", "P4/n-32/n       ", "Pn3n          *s",       &
       "O_h^2          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,223) /                                                 &
       "Pm-3n           ", "P4_2/m-32/n     ", "Pm3n          *s",       &
       "O_h^3          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,224) /                                                 &
       "Pn-3m           ", "P4_2/n-32/m     ", "Pn3m          *s",       &
       "O_h^4          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,225) /                                                 &
       "Fm-3m           ", "F4/m-32/m       ", "Fm3m          *s",       &
       "O_h^5          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,226) /                                                 &
       "Fm-3c           ", "F4/m-32/c       ", "Fm3c          *s",       &
       "O_h^6          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,227) /                                                 &
       "Fd-3m           ", "F4_1/d-32/m     ", "Fd3m          *s",       &
       "O_h^7          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,228) /                                                 &
       "Fd-3c           ", "F4_1/d-32/c     ", "Fd3c          *s",       &
       "O_h^8          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,229) /                                                 &
       "Im-3m           ", "I4/m-32/m       ", "Im3m          *s",       &
       "O_h^9          S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /

  DATA sg_names(:,230) /                                                 &
       "Ia-3d           ", "I4_1/a-32/d     ", "Ia3d          *s",       &
       "O_h^10         S", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                ",       &
       "                ", "                ", "                " /


  PUBLIC   spg_code, spg_name, find_space_group, sg_name, equivalent_tau, &
           sg_names, find_space_group_number, find_space_group_names, add_info,&
           set_add_info, set_fft_fact, project_frac_tran,  &
           shift_frac_tran, set_standard_sg, set_point_group_code, &
           check_code_group_ext, symmorphic_sg, sg_origin

  CONTAINS

!-----------------------------------------------------------------------
  SUBROUTINE set_add_info()
!-----------------------------------------------------------------------
  IMPLICIT NONE

!
! 17  P222_1
!
  add_info(1,17)='xyz     xyz'
  add_info(2,17)='zxy     yzx'
  add_info(3,17)='yzx     zxy'
!
! 18  P2_12_12
!
  add_info(1,18)='xyz     xyz'
  add_info(2,18)='zxy     yzx'
  add_info(3,18)='yzx     zxy'
!
! 20 C222_1
!
  add_info(1,20)='xyz     xyz'
  add_info(2,20)='zxy     yzx'
  add_info(3,20)='yzx     zxy'
!
! 21 C222
!
  add_info(1,21)='xyz     xyz'
  add_info(2,21)='zxy     yzx'
  add_info(3,21)='yzx     zxy'
!
! 25 Pmm2
!
  add_info(1,25)='xyz     xyz'
  add_info(2,25)='zxy     yzx'
  add_info(3,25)='yzx     zxy'
!
! 26 Pmc2_1
!
  add_info(1,26)='xyz     xyz'
  add_info(2,26)='yx-z    yx-z'
  add_info(3,26)='zxy     yzx'
  add_info(4,26)='-zyx    zy-x'
  add_info(5,26)='yzx     zxy'
  add_info(6,26)='x-zy    xz-y'
!
! 27 Pcc2
!
  add_info(1,27)='xyz     xyz'
  add_info(2,27)='zxy     yzx'
  add_info(3,27)='yzx     zxy'
!
! 28 Pma2
!
  add_info(1,28)='xyz     xyz'
  add_info(2,28)='yx-z    yx-z'
  add_info(3,28)='zxy     yzx'
  add_info(4,28)='-zyx    zy-x'
  add_info(5,28)='yzx     zxy'
  add_info(6,28)='x-zy    xz-y'
!
! 29 Pca2_1
!
  add_info(1,29)='xyz     xyz'
  add_info(2,29)='yx-z    yx-z'
  add_info(3,29)='zxy     yzx'
  add_info(4,29)='-zyx    zy-x'
  add_info(5,29)='yzx     zxy'
  add_info(6,29)='x-zy    xz-y'
!
! 30 Pnc2
!
  add_info(1,30)='xyz     xyz'
  add_info(2,30)='yx-z    yx-z'
  add_info(3,30)='zxy     yzx'
  add_info(4,30)='-zyx    zy-x'
  add_info(5,30)='yzx     zxy'
  add_info(6,30)='x-zy    xz-y'
!
! 31 Pmn2_1
!
  add_info(1,31)='xyz     xyz'
  add_info(2,31)='yx-z    yx-z'
  add_info(3,31)='zxy     yzx'
  add_info(4,31)='-zyx    zy-x'
  add_info(5,31)='yzx     zxy'
  add_info(6,31)='x-zy    xz-y'
!
! 32 Pba2
!
  add_info(1,32)='xyz     xyz'
  add_info(2,32)='zxy     yzx'
  add_info(3,32)='yzx     zxy'
!
! 33 Pna2_1
!
  add_info(1,33)='xyz     xyz'
  add_info(2,33)='yx-z    yx-z'
  add_info(3,33)='zxy     yzx'
  add_info(4,33)='-zyx    zy-x'
  add_info(5,33)='yzx     zxy'
  add_info(6,33)='x-zy    xz-y'
!
! 34 Pna2_1
!
  add_info(1,34)='xyz     xyz'
  add_info(2,34)='zxy     yzx'
  add_info(3,34)='yzx     zxy'
!
! 35 Cmm2
!
  add_info(1,35)='xyz     xyz'
  add_info(2,35)='zxy     yzx'
  add_info(3,35)='yzx     zxy'
!
! 36 Cmc2_1
!
  add_info(1,36)='xyz     xyz'
  add_info(2,36)='yx-z    yx-z'
  add_info(3,36)='zxy     yzx'
  add_info(4,36)='-zyx    zy-x'
  add_info(5,36)='yzx     zxy'
  add_info(6,36)='x-zy    xz-y'
!
! 37 Ccc2
!
  add_info(1,37)='xyz     xyz'
  add_info(2,37)='zxy     yzx'
  add_info(3,37)='yzx     zxy'
!
! 38 Amm2
!
  add_info(1,38)='xyz     xyz'
  add_info(2,38)='yx-z    yx-z'
  add_info(3,38)='zxy     yzx'
  add_info(4,38)='-zyx    zy-x'
  add_info(5,38)='yzx     zxy'
  add_info(6,38)='x-zy    xz-y'
!
! 39 Aem2
!
  add_info(1,39)='xyz     xyz'
  add_info(2,39)='yx-z    yx-z'
  add_info(3,39)='zxy     yzx'
  add_info(4,39)='-zyx    zy-x'
  add_info(5,39)='yzx     zxy'
  add_info(6,39)='x-zy    xz-y'
  add_info(7,39)='xyz     xyz'
  add_info(8,39)='yx-z    yx-z'
  add_info(9,39)='zxy     yzx'
  add_info(10,39)='-zyx   zy-x'
  add_info(11,39)='yzx    zxy'
  add_info(12,39)='x-zy   xz-y'
  add_info(13,39)='xyz    xyz'
  add_info(14,39)='yx-z   yx-z'
  add_info(15,39)='zxy    yzx'
  add_info(16,39)='-zyx   zy-x'
  add_info(17,39)='yzx    zxy'
  add_info(18,39)='x-zy   xz-y'
!
! 40 Ama2
!
  add_info(1,40)='xyz     xyz'
  add_info(2,40)='yx-z    yx-z'
  add_info(3,40)='zxy     yzx'
  add_info(4,40)='-zyx    zy-x'
  add_info(5,40)='yzx     zxy'
  add_info(6,40)='x-zy    xz-y'
!
! 41 Aea2
!
  add_info(1,41)='xyz     xyz'
  add_info(2,41)='yx-z    yx-z'
  add_info(3,41)='zxy     yzx'
  add_info(4,41)='-zyx    zy-x'
  add_info(5,41)='yzx     zxy'
  add_info(6,41)='x-zy    xz-y'
  add_info(7,41)='xyz     xyz'
  add_info(8,41)='yx-z    yx-z'
  add_info(9,41)='zxy     yzx'
  add_info(10,41)='-zyx   zy-x'
  add_info(11,41)='yzx    zxy'
  add_info(12,41)='x-zy   xz-y'
  add_info(13,41)='xyz    xyz'
  add_info(14,41)='yx-z   yx-z'
  add_info(15,41)='zxy    yzx'
  add_info(16,41)='-zyx   zy-x'
  add_info(17,41)='yzx    zxy'
  add_info(18,41)='x-zy   xz-y'
!
! 42 Fmm2
!
  add_info(1,42)='xyz     xyz'
  add_info(2,42)='zxy     yzx'
  add_info(3,42)='yzx     zxy'
!
! 43 Fdd2
!
  add_info(1,43)='xyz     xyz'
  add_info(2,43)='zxy     yzx'
  add_info(3,43)='yzx     zxy'
!
! 44 Imm2
!
  add_info(1,44)='xyz     xyz'
  add_info(2,44)='zxy     yzx'
  add_info(3,44)='yzx     zxy'
!
! 45 Iba2
!
  add_info(1,45)='xyz     xyz'
  add_info(2,45)='zxy     yzx'
  add_info(3,45)='yzx     zxy'
!
! 46 Ima2
!
  add_info(1,46)='xyz     xyz'
  add_info(2,46)='yx-z    yx-z'
  add_info(3,46)='zxy     yzx'
  add_info(4,46)='-zyx    zy-x'
  add_info(5,46)='yzx     zxy'
  add_info(6,46)='x-zy    xz-y'
!
! 49 Pccm
!
  add_info(1,49)='xyz     xyz'
  add_info(2,49)='zxy     yzx'
  add_info(3,49)='yzx     zxy'
!
! 50 Pban
!
  add_info(1,50)='xyz     xyz'
  add_info(2,50)='zxy     yzx'
  add_info(3,50)='yzx     zxy'
!
! 51 Pmma
!
  add_info(1,51)='xyz     xyz'
  add_info(2,51)='yx-z    yx-z'
  add_info(3,51)='zxy     yzx'
  add_info(4,51)='-zyx    zy-x'
  add_info(5,51)='yzx     zxy'
  add_info(6,51)='x-zy    xz-y'
!
! 52 Pmma
!
  add_info(1,52)='xyz     xyz'
  add_info(2,52)='yx-z    yx-z'
  add_info(3,52)='zxy     yzx'
  add_info(4,52)='-zyx    zy-x'
  add_info(5,52)='yzx     zxy'
  add_info(6,52)='x-zy    xz-y'
!
! 53 Pmna
!
  add_info(1,53)='xyz     xyz'
  add_info(2,53)='yx-z    yx-z'
  add_info(3,53)='zxy     yzx'
  add_info(4,53)='-zyx    zy-x'
  add_info(5,53)='yzx     zxy'
  add_info(6,53)='x-zy    xz-y'
!
! 54 Pcca
!
  add_info(1,54)='xyz     xyz'
  add_info(2,54)='yx-z    yx-z'
  add_info(3,54)='zxy     yzx'
  add_info(4,54)='-zyx    zy-x'
  add_info(5,54)='yzx     zxy'
  add_info(6,54)='x-zy    xz-y'
!
! 55 Pbam
!
  add_info(1,55)='xyz     xyz'
  add_info(2,55)='zxy     yzx'
  add_info(3,55)='yzx     zxy'
!
! 56 Pccn
!
  add_info(1,56)='xyz     xyz'
  add_info(2,56)='zxy     yzx'
  add_info(3,56)='yzx     zxy'
!
! 57 Pbcm
!
  add_info(1,57)='xyz     xyz'
  add_info(2,57)='yx-z    yx-z'
  add_info(3,57)='zxy     yzx'
  add_info(4,57)='-zyx    zy-x'
  add_info(5,57)='yzx     zxy'
  add_info(6,57)='x-zy    xz-y'
!
! 58 Pnnm
!
  add_info(1,58)='xyz     xyz'
  add_info(2,58)='zxy     yzx'
  add_info(3,58)='yzx     zxy'
!
! 59 Pmnm
!
  add_info(1,59)='xyz     xyz'
  add_info(2,59)='zxy     yzx'
  add_info(3,59)='yzx     zxy'
!
! 60 Pmnm
!
  add_info(1,60)='xyz     xyz'
  add_info(2,60)='yx-z    yx-z'
  add_info(3,60)='zxy     yzx'
  add_info(4,60)='-zyx    zy-x'
  add_info(5,60)='yzx     zxy'
  add_info(6,60)='x-zy    xz-y'
!
! 61 Pbca
!
  add_info(1,61)='xyz     xyz'
  add_info(2,61)='yx-z    yx-z'
  add_info(3,61)='zxy     yzx'
  add_info(4,61)='-zyx    zy-x'
  add_info(5,61)='yzx     zxy'
  add_info(6,61)='x-zy    xz-y'
!
! 62 Pnma
!
  add_info(1,62)='xyz     xyz'
  add_info(2,62)='yx-z    yx-z'
  add_info(3,62)='zxy     yzx'
  add_info(4,62)='-zyx    zy-x'
  add_info(5,62)='yzx     zxy'
  add_info(6,62)='x-zy    xz-y'
!
! 63 Cmcm
!
  add_info(1,63)='xyz     xyz'
  add_info(2,63)='yx-z    yx-z'
  add_info(3,63)='zxy     yzx'
  add_info(4,63)='-zyx    zy-x'
  add_info(5,63)='yzx     zxy'
  add_info(6,63)='x-zy    xz-y'
!
!  64 Cmce
!
  add_info(1,64)='xyz     xyz'
  add_info(2,64)='yx-z    yx-z'
  add_info(3,64)='zxy     yzx'
  add_info(4,64)='-zyx    zy-x'
  add_info(5,64)='yzx     zxy'
  add_info(6,64)='x-zy    xz-y'
  add_info(7,64)='xyz     xyz'
  add_info(8,64)='yx-z    yx-z'
  add_info(9,64)='zxy     yzx'
  add_info(10,64)='-zyx   zy-x'
  add_info(11,64)='yzx    zxy'
  add_info(12,64)='x-zy   xz-y'
  add_info(13,64)='xyz    xyz'
  add_info(14,64)='yx-z   yx-z'
  add_info(15,64)='zxy    yzx'
  add_info(16,64)='-zyx   zy-x'
  add_info(17,64)='yzx    zxy'
  add_info(18,64)='x-zy   xz-y'
!
!  65 Cmmm
!
  add_info(1,65)='xyz     xyz'
  add_info(2,65)='zxy     yzx'
  add_info(3,65)='yzx     zxy'
!
!  66 Cccm
!
  add_info(1,66)='xyz     xyz'
  add_info(2,66)='zxy     yzx'
  add_info(3,66)='yzx     zxy'
!
!  67 Cmme
!
  add_info(1,67)='xyz     xyz'
  add_info(2,67)='zxy     yzx'
  add_info(3,67)='yzx     zxy'
  add_info(4,67)='xyz     xyz'
  add_info(5,67)='zxy     yzx'
  add_info(6,67)='yzx     zxy'
  add_info(7,67)='xyz     xyz'
  add_info(8,67)='zxy     yzx'
  add_info(9,67)='yzx     zxy'
!
!  68 Ccce
!
  add_info(1,68)='xyz     xyz'
  add_info(2,68)='zxy     yzx'
  add_info(3,68)='yzx     zxy'
  add_info(4,68)='xyz     xyz'
  add_info(5,68)='zxy     yzx'
  add_info(6,68)='yzx     zxy'
  add_info(7,68)='xyz     xyz'
  add_info(8,68)='zxy     yzx'
  add_info(9,68)='yzx     zxy'
!
! 72 Ibam
!
  add_info(1,72)='xyz     xyz'
  add_info(2,72)='zxy     yzx'
  add_info(3,72)='yzx     zxy'
!
! 73 Ibca
!
  add_info(1,73)='xyz     xyz'
  add_info(2,73)='yx-z    yx-z'
!
! 74 Imma
!
  add_info(1,74)='xyz     xyz'
  add_info(2,74)='yx-z    yx-z'
  add_info(3,74)='zxy     yzx'
  add_info(4,74)='-zyx    zy-x'
  add_info(5,74)='yzx     zxy'
  add_info(6,74)='x-zy    xz-y'

  RETURN
  END SUBROUTINE set_add_info

!-----------------------------------------------------------------------
  SUBROUTINE set_fft_fact(code, unique, trig, fft_fact)
!-----------------------------------------------------------------------
!
! This routine recives as input the code of a space group and gives as
! output the factor that should be contained in the dimensions of the
! fft to accomodate the fractional translations of that space group.
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: code, unique, trig
  INTEGER, INTENT(OUT) :: fft_fact(3)
  
  fft_fact=1
  SELECT CASE (code)
     CASE (4)
        IF (unique==1) THEN
!
!  b_unique
!
          fft_fact(2)=2
        ELSE
          fft_fact(3)=2
        ENDIF
     CASE(7)
        IF (unique==1) THEN
          fft_fact(3)=2
        ELSE
          fft_fact(1)=2
        ENDIF
     CASE(9)
        IF (unique==1) THEN
          fft_fact(3)=2
        ELSE
          fft_fact(1)=2
        ENDIF
     CASE(11)
        IF (unique==1) THEN
          fft_fact(2)=2
        ELSE
          fft_fact(3)=2
        ENDIF
     CASE(12)
         fft_fact(1)=2
         fft_fact(2)=2
     CASE(13)
        IF (unique==1) THEN
          fft_fact(3)=2
        ELSE
          fft_fact(1)=2
        ENDIF
     CASE(14)
        IF (unique==1) THEN
          fft_fact(2)=2
          fft_fact(3)=2
        ELSE
          fft_fact(1)=2
          fft_fact(3)=2
        ENDIF
     CASE(15)
        IF (unique==1) THEN
          fft_fact(3)=2
        ELSE
          fft_fact(1)=2
        ENDIF
     CASE(17)
        fft_fact(3)=2
     CASE(18)
        fft_fact(1)=2
        fft_fact(2)=2
     CASE(19,24,33,34,45,46,48,52,56,58,60,61,62,68,70,72,73,74,86,94,102,&
          104,106,114,118,120,126,128,130,133,134,135,136,137,138,&
          140,198,199,201,205,206,208,218,219,222,223,224,226)
        fft_fact=2
     CASE(20,26,27,36,37,49,63,66,67,77,84,93,101,103,105,108,112,116,&
          124,131,132,158,159,163,165,173,176,182,184,185,186,188,190,&
          192,193)
        fft_fact(3)=2
     CASE(28,40,51)
        fft_fact(1)=2
     CASE(29,31,53,54)
        fft_fact(1)=2
        fft_fact(3)=2
     CASE(30,57,64)
        fft_fact(2)=2
        fft_fact(3)=2
     CASE(32,41,50,55,59,85,90,100,113,117,125,127,129)
        fft_fact(1)=2
        fft_fact(2)=2
     CASE(39)
        fft_fact(2)=2
     CASE(76,78,91,95)
        fft_fact(3)=4
     CASE(80,92,96)
        fft_fact(1)=2
        fft_fact(2)=2
        fft_fact(3)=4
     CASE(144,145,151,152,153,154,171,172,180,181)
        fft_fact(3)=3
     CASE(161,167)
        IF (trig==1) THEN
           fft_fact=2
        ELSE
           fft_fact(3)=2
        ENDIF
     CASE(169,170,178,179)
        fft_fact(3)=6
     CASE(43,88,98,109,110,122,141,142,203,210,212,213,214,220,227,228,230)
        fft_fact=4
     CASE(194)
        fft_fact(1)=3
        fft_fact(2)=3
        fft_fact(3)=2
  END SELECT 

  RETURN
  END SUBROUTINE set_fft_fact

!-----------------------------------------------------------------------
  SUBROUTINE find_space_group(ibrav, nsym, sr_in, ft_in, at, bg, &
                                     sg_number, aux_sg, s01, s02, verbosity)
!-----------------------------------------------------------------------
  !
  !  This routine receives as input: 
  !  the bravais lattice, the number of symmetries, the rotation matrices 
  !  in cartesian coordinates, the fractional translations in crystal 
  !  coordinates (true fractional translation with a - sign with respect to
  !  those used in QE), the direct and reciprocal lattice vectors.  
  !  if verbosity .true. the routine writes informations on the symmetries
  !  ibrav=0 is not supported
  !  It gives as output the space group number and the two vectors s01 and 
  !  s02 which give the position of the origin with origin choice 1 and 
  !  origin choice 2 according to ITA. 
  !  If the ITA give only one origin s02 is set to 1000.0,1000.0,1000.0
  !
  !  NB: it assumes that the space group orientation is the one of the ITA
  !      tables.
  !  aux_sg is an additional output code, used only for selected space groups.
  !      In the monoclinic and orthorhombic case it is a number from 1 to 6 
  !      that gives the orientation. It uses the same convention of table
  !      4.3.2.1 of the ITA tables, and indicates the column if all columns
  !      are different. If some columns coincide it is a pointer in the
  !      list of group names.
  !
  USE io_global, ONLY : stdout
  USE point_group, ONLY : find_group_info_ext, sym_label
  USE lattices,    ONLY : is_centered, compute_conventional, lattice_name
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: ibrav, nsym
  REAL(DP), INTENT(IN)  :: sr_in(3,3,nsym), ft_in(3,nsym), at(3,3), bg(3,3)
  LOGICAL, INTENT(IN)   :: verbosity
  INTEGER, INTENT(OUT)  :: sg_number, aux_sg
  REAL(DP), INTENT(OUT) :: s01(3), s02(3)

  LOGICAL :: is_symmorphic
  INTEGER ::                   &
             code_group,       &  ! the code of the point group
             code_group_ext,   &  ! the extended code of the point group
             group_desc(48),   &  ! the number of the symmetry of each operation
             group_desc_sg(48), & ! the number of the space group symmetry
             which_elem(48)      ! the link between standard and input order

  INTEGER :: iaxis, imirror, irot90, irot60, irot120, imirror1, imirror2, &
             idir, irot180, nsaz, isym, jsym, ipol, nft, type_group, &
             axis_dir(3), mirr_sym(3), idir1, idir2,                 &
             nsa, nmir, ncomp, nmp, ipa, gd_sg ! auxiliary

  REAL(DP) ::  &
              sr(3,3,48),                 &! symmetry matrices
              ft(3,48), ftc(3,48),        &! fractional translation crys., cart.
              ftpar(3,48), ftparc(3,48),  &! fract. transl. parallel symmetry
              ftperp(3,48), ftperpc(3,48),&! fract. transl. perp. symmetry
              s0(3,48), s0c(3,48),        &! origin of each symmetry
              ax(3), bx(3), angle, fact,  &! auxiliary
              ftrs1(3), ftrs2(3), ftrs3(3),& ! auxiliary
              atp(3,3), bgp(3,3),          & ! conventional bravais lattice
              srm(3,3), ftm(3), ftparm(3), & ! auxiliary quantities
              ftperpm(3), s0m(3),          & ! ausiliary quantities
              sp, sp1, sp2, nftd             ! auxiliary

  REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP, eps1=1.D-8
  LOGICAL :: is_axis, ft_zero, is_proper
  INTEGER :: tipo_sym, ts
  CHARACTER(LEN=11) :: group_name
  CHARACTER(LEN=40) :: latt_name
!
!  initialization
!  
  s01=1000.0_DP
  s02=1000.0_DP
  aux_sg=1
  IF (ibrav==0) THEN
     sg_number=0 
     RETURN
  ENDIF
!
!   Find the point group
!
  CALL find_group_info_ext(nsym, sr_in, code_group, code_group_ext, &
                                                  which_elem, group_desc) 

  IF (.NOT.check_code_group_ext(code_group_ext)) &
     CALL errore('find_space_group','point group orientation incorrect',1)

  IF (verbosity) THEN
     WRITE(stdout,'(/,5x,"Space group identification, ",i3," symmetries:")') &
                                                              nsym
     CALL lattice_name(ibrav, latt_name)

     WRITE(stdout,'(/,5x,"Bravais lattice ", i3, 2x, a)') ibrav, latt_name
     WRITE(stdout,'(5x,"Point group number ", i3," / ",i3,2x,a)') &
                    code_group, code_group_ext, TRIM(group_name(code_group)) 
                                                              
  ENDIF
!
!  and order the elements in the standard order for each point group
!  project each fractionary translation along or perpendicularly to
!  the symmetry element and find the tag of each space group operation.
!  Find also the origin where the operation is applied.
!
  DO isym = 1, nsym
     DO jsym = 1, nsym
        IF (which_elem(jsym) == isym) THEN
           ft(:,isym) = ft_in(:,jsym)
           sr(:,:,isym) = sr_in(:,:,jsym)
           CALL project_frac_tran(sr(1,1,isym),ft(1,isym),at,bg,&
                           ftperp(1,isym),ftpar(1,isym))
           CALL find_sg_tags(sr(1,1,isym),ft(1,isym),ibrav,at,bg, &
                             ftperp(1,isym),ftpar(1,isym),group_desc(isym),&
                             s0(1,isym),group_desc_sg(isym))
        END IF
     END DO
  END DO
!
!  if any operation is a glide plane or a screw axis, set to false the
!  is_symmorphic flag
!
  is_symmorphic=.TRUE.
  DO isym=1,nsym
     is_symmorphic = (is_symmorphic .AND. group_desc_sg(isym) < 65)
     IF (.NOT.is_symmorphic) EXIT
  END DO
!
!  if all fractional translations vanish set to .TRUE. ft_zero flag
!
  ft_zero=.NOT.(ANY(ABS(ft(:,1:nsym)) > eps1))
!
!  fractional translations in cartesian coordinates, 
!

  ftc(:,1:nsym)=ft(:,1:nsym)
  ftperpc(:,1:nsym)=ftperp(:,1:nsym)
  ftparc(:,1:nsym)=ftpar(:,1:nsym)
  s0c(:,1:nsym)=s0(:,1:nsym)

  CALL cryst_to_cart(nsym,ftc,at,1)
  CALL cryst_to_cart(nsym,ftperpc,at,1)
  CALL cryst_to_cart(nsym,ftparc,at,1)
  CALL cryst_to_cart(nsym,s0c,at,1)

  IF (is_centered(ibrav)) THEN
!
!  for centered lattices find the translations
!  in crystal coordinates of the conventional cell
!
     CALL compute_conventional(at, atp, ibrav)
     CALL recips(atp(1,1), atp(1,2), atp(1,3), bgp(1,1), bgp(1,2), bgp(1,3))
     CALL cryst_to_cart(nsym,ftc,bgp,-1)
     CALL cryst_to_cart(nsym,ftperpc,bgp,-1)
     CALL cryst_to_cart(nsym,ftparc,bgp,-1)
     CALL cryst_to_cart(nsym,s0c,bgp,-1)
  ENDIF
!
!   write a summary of the symmetry operation on output
!
  IF (verbosity) THEN
     IF (ft_zero) THEN
        WRITE(stdout,'(/,5x,"Nonsymmorphic operations not found: All &
                        &fractional translations vanish")')
        WRITE(stdout,'(5x,"Symmetries of the point group in standard order",/)')

        DO isym=1,nsym
           WRITE(stdout,'(5x,i4,a8,i4)') isym, &
              TRIM(sym_label(group_desc(isym))), group_desc(isym)
        END DO
        WRITE(stdout,*)

     ELSE
        WRITE(stdout,'(/,5x,"The coset representatives of the space group &
                        &with")') 
        WRITE(stdout,'(5x,"the standard order of operations are:")')
        WRITE(stdout,'(/,5x,"Fractional translations in crystal coordinates")')
        WRITE(stdout,'(/,6x, "PGS",4x,"Fract. transl. (all)         &
                               & Fract. transl. (no shift)",3x,"SGS",/)')
        DO isym=1,nsym
           WRITE(stdout,'(1x,a8,i4,3f8.4,4x,3f8.4,a10,i4)') &
              TRIM(sym_label(group_desc(isym))), group_desc(isym), &
              ft(1,isym), ft(2,isym), ft(3,isym), &
              ftpar(1,isym), ftpar(2,isym), ftpar(3,isym), &
              TRIM(sym_label_sg(group_desc_sg(isym))), group_desc_sg(isym) 
        END DO

        IF (is_centered(ibrav)) THEN
           WRITE(stdout,'(/,5x,"Fractional translations in conventional &
                             &crystal coordinates")')
           WRITE(stdout,'(/,6x, "PGS",6x,"Fract. transl. (all)       &
                               & Fract. transl. (no shift)",3x,"SGS",/)')
           DO isym=1,nsym
              WRITE(stdout,'(1x,a8,i4,3f8.4,4x,3f8.4,a10,i4)') &
                 TRIM(sym_label(group_desc(isym))), group_desc(isym), &
                 ftc(1,isym), ftc(2,isym), ftc(3,isym), &
                 ftparc(1,isym), ftparc(2,isym), ftparc(3,isym), &
                 TRIM(sym_label_sg(group_desc_sg(isym))), group_desc_sg(isym) 
           END DO
        END IF

        WRITE(stdout,'(/,5x,"PGS = Point group symmetry, &
                                        &SGS = space group symmetry ")')
        WRITE(stdout,'(5x,"(no shift) means the part of the fractional &
                        &translation")')
        WRITE(stdout,'(5x,"that cannot be removed by an origin shift")')

        WRITE(stdout,'(/,5x,"Fractional translation shift and origin &
                       &shift in crystal coordinates")')
        WRITE(stdout,'(/,6x, "PGS",6x,"FT shift &
          &                   &Origin shift",17x,"SGS",/)')
        DO isym=1,nsym
           WRITE(stdout,'(1x,a8,i4,3f8.4,4x,3f8.4,a10,i4)') &
              TRIM(sym_label(group_desc(isym))), group_desc(isym), &
              ftperp(1,isym), ftperp(2,isym), ftperp(3,isym), &
              s0(1,isym), s0(2,isym), s0(3,isym), &
              TRIM(sym_label_sg(group_desc_sg(isym))), group_desc_sg(isym) 
        END DO

        IF (is_centered(ibrav)) THEN
           WRITE(stdout,'(/,5x,"Fractional translation shift and origin shift &
                         &in conventional")')
           WRITE(stdout,'(5x,"crystal coordinates")')
           WRITE(stdout,'(/,6x, "PGS",6x,"FT shift &
              &                   Origin shift",17x,"SGS",/)')
           DO isym=1,nsym
              WRITE(stdout,'(1x,a8,i4,3f8.4,4x,3f8.4,a10,i4)') &
                 TRIM(sym_label(group_desc(isym))), group_desc(isym), &
                 ftperpc(1,isym), ftperpc(2,isym), ftperpc(3,isym), &
                 s0c(1,isym), s0c(2,isym), s0c(3,isym), &
                 TRIM(sym_label_sg(group_desc_sg(isym))), group_desc_sg(isym) 
           END DO
        ENDIF
        WRITE(stdout,*)
        WRITE(stdout,'(/,5x,"The origin shift gives the point of application&
                         & of the symmetry")')
        WRITE(stdout,'(5x,"Subtract this vector from all atomic coordinates")')
        WRITE(stdout,'(5x,"to have this symmetry centered at the origin")')
        WRITE(stdout,*)
     ENDIF
  ENDIF

  sg_number=0
  SELECT CASE (code_group)
      CASE(1)
!
!  1
!
         s01(:)=0.0_DP
         sg_number=1
      CASE(2)
!
!  -1
!
         s01(:)=s0(:,2)
         sg_number=2
      CASE(3)
!
!  C_s
!
          s01(:)=s0(:,2)
          IF (ibrav==13) THEN
             IF (group_desc_sg(2)==34.OR.group_desc_sg(2)==126) THEN
!
!   B11m  unique c
!      a
!
                sg_number=8
             ELSE 
                sg_number=9
             ENDIF
          ELSEIF (ibrav==-13) THEN
             IF (group_desc_sg(2)==35.OR.group_desc_sg(2)==130) THEN
!
!   C1m1  unique b
!     a
!
                sg_number=8
             ELSE 
                sg_number=9
             ENDIF
          ELSEIF(ABS(ibrav)==12) THEN
             IF (is_symmorphic) THEN
                sg_number=6
             ELSE 
                sg_number=7
             ENDIF
          ENDIF
      CASE(4)
!
!   C_2
!
         s01(:)=s0(:,2)
         IF (ABS(ibrav)==13) THEN
            sg_number=5
         ELSE
            IF (is_symmorphic) THEN
               sg_number=3
            ELSE
               sg_number=4
            ENDIF
         ENDIF
      CASE(5)
!
!  C_3
!
         s01(:)=s0(:,2)
         IF (ibrav==4) THEN
            IF (is_symmorphic) THEN
               sg_number=143
            ELSE
               IF (MOD(group_desc_sg(2),2)==0) THEN
!
!  All 3_1 operations have odd number
!
                  sg_number=144
               ELSE
                  sg_number=145
               ENDIF
            ENDIF
         ELSEIF (ibrav==5) THEN
            sg_number=146
         ENDIF
      CASE(6)
!
!   C_4
!
         s01(:)=s0(:,2)
         s01(3)=0.0_DP
         IF (ibrav==6) THEN
            IF (is_symmorphic) THEN
               sg_number=75
            ELSE
               IF (group_desc_sg(2)==73.OR.group_desc_sg(2)==78.OR. &
                   group_desc_sg(2)==89) THEN
                   sg_number=76
               ELSEIF (group_desc_sg(2)==74.OR.group_desc_sg(2)==79.OR. &
                   group_desc_sg(2)==90) THEN
                   sg_number=77
               ELSEIF (group_desc_sg(2)==75.OR.group_desc_sg(2)==80.OR. &
                   group_desc_sg(2)==91) THEN
                   sg_number=78
               ELSE
                   CALL errore('find_space_group',&
                               'symmetry not recognized C_4',1)
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (group_desc_sg(2)==8.OR.group_desc_sg(2)==74) THEN
               sg_number=79
            ELSE
               sg_number=80
               s01(:)=s0(:,3)
            ENDIF
         ENDIF
      CASE(7)
!
!  C_6
!
!
!   origin on the C_6 axis
!
         s01(:)=s0(:,2)
         IF (is_symmorphic) THEN
            sg_number=168
         ELSE
            IF (group_desc_sg(2)==108) THEN
!
!  6_1
!
               sg_number=169
            ELSEIF (group_desc_sg(2)==112) THEN
!
!  6_5
!
               sg_number=170
            ELSEIF (group_desc_sg(2)==109) THEN
!
!  6_2
!
               sg_number=171
            ELSEIF (group_desc_sg(2)==111) THEN
!
!  6_4
!
               sg_number=172
            ELSEIF (group_desc_sg(2)==110) THEN
!
!  6_3
!
               sg_number=173
            ENDIF
         ENDIF
      CASE(8)
!
!   D_2
!
         s01(:)=s0(:,2)
         s01(3)=s0(3,3)
         IF (ibrav==8) THEN
            IF (is_symmorphic) THEN
               sg_number=16
            ELSE
               nft=0
               IF (group_desc_sg(2)/=2) nft=nft+1
               IF (group_desc_sg(3)/=4) nft=nft+1
               IF (group_desc_sg(4)/=3) nft=nft+1
               IF ( nft == 0 ) CALL errore('find_space_group',&
                              'D_2: unexpected number of screw axis',1)
               IF (nft==1) THEN
                  sg_number=17
                  IF (group_desc_sg(2)/=2) aux_sg=1
                  IF (group_desc_sg(3)/=4) aux_sg=3
                  IF (group_desc_sg(4)/=3) aux_sg=5
               ENDIF
               IF (nft==2) THEN
                  sg_number=18
                  IF (group_desc_sg(2)==2) aux_sg=1
                  IF (group_desc_sg(3)==4) aux_sg=3
                  IF (group_desc_sg(4)==3) aux_sg=5
               ENDIF
               IF (nft==3) THEN
                  sg_number=19
!
!   ITA origin not programmed
!
                  s01=1000.0_DP
               ENDIF
            ENDIF
         ELSEIF (ibrav==9) THEN
            nft=0
            IF (group_desc_sg(2)/=2) nft=nft+1
            IF (group_desc_sg(3)/=4) nft=nft+1
            IF (group_desc_sg(4)/=3) nft=nft+1
            IF (nft==0.OR.nft==2) THEN
               sg_number=21
               IF (nft==0) THEN
                  aux_sg=1
               ELSE
                  IF (group_desc_sg(2)==2) aux_sg=1
                  IF (group_desc_sg(3)==4) aux_sg=3
                  IF (group_desc_sg(4)==3) aux_sg=5
               ENDIF
            ELSE
               sg_number=20
            ENDIF
         ELSEIF (ibrav==10) THEN
            sg_number=22
         ELSEIF (ibrav==11) THEN
!
!  If the three axis intesect in a point, or can be made to intersect 
!  by a direct lattice shift, the group is 23, otherwise it is 24
!
            IF (check_intersection(ftperpc(1,2), ftperpc(1,3), &
                                                 ftperpc(1,4))) THEN
               sg_number=23
            ELSE
               sg_number=24
!
!  ITA origin not programmed
!
               s01=1000.0_DP
            ENDIF
         ENDIF
      CASE(9)
!
!  D_3
!
!  origin at the intersection of the axis of order 3 and the plane that
!  contains the axis of order 2. Note that the axis 3 must be parallel to
!  z, the other cases are not programmed yet.
!
         s01(:)=s0(:,2)
         s01(3)=s0(3,4)
         type_group=100
         
         IF (code_group_ext==44) type_group=1
         IF (code_group_ext==45) type_group=0
         IF (type_group==100) CALL errore('find_space_group', &
                                            'group not recognized',1)
         IF (ibrav==4) THEN
!
!  groups of type 1 x axis is of order 2
!
            IF (type_group==1) THEN
               IF (group_desc_sg(2)==27) THEN
!
! axis 3 is proper
!
                  sg_number=150
               ELSEIF (MOD(group_desc_sg(2),2)==0) THEN
!
!   axis 3_1
!
                  sg_number=152
!
!   origin on the plane that contains the 2110 operation
!
                  s01(3)=s0(3,5)
               ELSE 
!
!  axis 3_2
!
                  sg_number=154
               ENDIF
            ELSEIF (type_group==0) THEN
!
!  groups of type 0 y axis is of order 2
!
               IF (group_desc_sg(2)==27) THEN
                  sg_number=149
               ELSEIF (MOD(group_desc_sg(2),2)==0) THEN
                  sg_number=151
!
!   origin on the plane that contains the 2210 operation (number 4)
!
               ELSE 
                  sg_number=153
               ENDIF
            ENDIF
         ELSEIF (ibrav==5) THEN
            sg_number=155
         ENDIF     
      CASE(10)
!
!  D_4
!
!
!  Origin at the intersection of axis 4 and 2x
!
         s01(:)=s0(:,2)
         s01(3)=s0(3,5)
         IF (ibrav==6) THEN
!
!   Find the axis of order 4, and the axis of order 2 parallel to x
!
            irot90=2
            irot180=5
            IF (group_desc_sg(irot90)==8) THEN
               IF (group_desc_sg(irot180)==4) THEN
                  sg_number=89
               ELSE
!
!  axis 4
!
                  sg_number=90
                  s01(:)=s0(:,3)
                  s01(3)=s0(3,5)
               ENDIF
            ELSEIF (group_desc_sg(irot90) == 73) THEN
!
!  axis 4_1
!
               IF (group_desc_sg(irot180)==4) THEN
                  sg_number=91
                  s01(:)=s0(:,2)
                  s01(3)=s0(3,7)
               ELSE
                  sg_number=92
                  s01(:)=s0(:,3)
                  s01(3)=s0(3,5)
               ENDIF
            ELSEIF (group_desc_sg(irot90) == 74) THEN
!
!  axis 4_2
!
               IF (group_desc_sg(irot180)==4) THEN
                  sg_number=93
                  s01(:)=s0(:,3)
                  s01(3)=s0(3,5)
               ELSE
                  sg_number=94
                  s01(:)=s0(:,3)
                  s01(3)=s0(3,5)
               ENDIF
            ELSEIF (group_desc_sg(irot90) == 75) THEN
!
!  axis 4_3
!
               IF (group_desc_sg(irot180)==4) THEN
                  sg_number=95
                  s01(:)=s0(:,2)
                  s01(3)=s0(3,7)
               ELSE
                  sg_number=96
                  s01(:)=s0(:,3)
                  s01(3)=s0(3,6)
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (group_desc_sg(2)==8.OR.group_desc_sg(2)==74) THEN
               sg_number=97
            ELSE 
               sg_number=98
               s01(:)=s0(:,3)
               s01(3)=s0(3,5)
            ENDIF
         ENDIF
      CASE(11)
!
!  D_6
!
         s01(:)=s0(:,2)
         s01(3)=s0(3,7)
         irot60=2
         IF (group_desc_sg(irot60) == 25) sg_number=177   ! 6
         IF (group_desc_sg(irot60) == 108) sg_number=178   ! 6_1
         IF (group_desc_sg(irot60) == 112) sg_number=179   ! 6_5
         IF (group_desc_sg(irot60) == 109) sg_number=180   ! 6_2
         IF (group_desc_sg(irot60) == 111) sg_number=181   ! 6_4
         IF (group_desc_sg(irot60) == 110) sg_number=182   ! 6_3
      CASE(12)
!
!   C_2v
!
!  identify the axis of order 2
!
         irot180=2
         imirror1=3
         imirror2=4
         ftrs1=ftc(:,3)
         ftrs2=ftc(:,4)
         IF (ibrav==8) THEN
            IF (code_group_ext==54) THEN
               idir=1
               is_proper=(group_desc_sg(2)==4) 
               ftrs1(2)=0.0_DP
               ftrs2(3)=0.0_DP
!
!   the origin on the C_2 axis
!
               s01(:)=s0(:,2)
               s01(1)=0.0_DP
               aux_sg=3
               IF (is_proper) THEN
                  IF (group_desc_sg(3)==35.AND.group_desc_sg(4)==34) THEN 
                     sg_number=25
                  ELSEIF (group_desc_sg(3)==130.AND.group_desc_sg(4)==126) THEN
                     sg_number=27
                  ELSEIF (group_desc_sg(3)==35.AND.group_desc_sg(4)==127) THEN
                     sg_number=28
                  ELSEIF (group_desc_sg(3)==131.AND.group_desc_sg(4)==36) THEN 
                     aux_sg=4
                     sg_number=28
                  ELSEIF (group_desc_sg(3)==132.AND.group_desc_sg(4)==126) THEN 
                     sg_number=30
                  ELSEIF(group_desc_sg(3)==130.AND.group_desc_sg(4)==128) THEN 
                     aux_sg=4
                     sg_number=30
                  ELSEIF (group_desc_sg(3)==131.AND.group_desc_sg(4)==127) THEN 
                     sg_number=32
                  ELSEIF (group_desc_sg(3)==132.AND.group_desc_sg(4)==138) THEN 
                     sg_number=34
                  ENDIF 
               ELSE 
                  IF (group_desc_sg(3)==35.AND.group_desc_sg(4)==126) THEN
                     sg_number=26
                  ELSEIF (group_desc_sg(3)==130.AND.group_desc_sg(4)==34) THEN 
                     aux_sg=4
                     sg_number=26
                  ELSEIF (group_desc_sg(3)==130.AND.group_desc_sg(4)==127) THEN
                     sg_number=29
                  ELSEIF (group_desc_sg(3)==131.AND.group_desc_sg(4)==126) THEN 
                     aux_sg=4
                     sg_number=29
                  ELSEIF (group_desc_sg(3)==35.AND.group_desc_sg(4)==128) THEN
                     s01(1)=s0(1,3)
                     s01(2)=s0(2,4)
                     s01(3)=0.0_DP
                     sg_number=31
                  ELSEIF (group_desc_sg(3)==132.AND.group_desc_sg(4)==34) THEN 
                     s01(1)=s0(1,3)
                     s01(2)=s0(2,4)
                     s01(3)=0.0_DP
                     aux_sg=4
                     sg_number=31
                  ELSEIF (group_desc_sg(3)==132.AND.group_desc_sg(4)==127) THEN
                     sg_number=33
                  ELSEIF (group_desc_sg(3)==131.AND.group_desc_sg(4)==128) THEN 
                     aux_sg=4
                     sg_number=33
                  ENDIF 
               ENDIF                 
            ELSEIF (code_group_ext==56) THEN
               idir=2
               is_proper=(group_desc_sg(2)==3) 
               ftrs1(3)=0.0_DP
               ftrs2(1)=0.0_DP
!
!   the origin on the C_2 axis
!
               s01(:)=s0(:,2)
               s01(2)=0.0_DP
               aux_sg=5
               IF (is_proper) THEN
                  IF (group_desc_sg(3)==34.AND.group_desc_sg(4)==36) THEN 
                     sg_number=25
                  ELSEIF (group_desc_sg(3)==127.AND.group_desc_sg(4)==134) THEN
                     sg_number=27
                  ELSEIF (group_desc_sg(3)==34.AND.group_desc_sg(4)==135) THEN
                     sg_number=28
                  ELSEIF (group_desc_sg(3)==126.AND.group_desc_sg(4)==36) THEN 
                     aux_sg=6
                     sg_number=28
                  ELSEIF (group_desc_sg(3)==128.AND.group_desc_sg(4)==134) THEN 
                     sg_number=30
                  ELSEIF(group_desc_sg(3)==127.AND.group_desc_sg(4)==136) THEN 
                     aux_sg=6
                     sg_number=30
                  ELSEIF (group_desc_sg(3)==128.AND.group_desc_sg(4)==136) THEN 
                     sg_number=34
                  ELSEIF (group_desc_sg(3)==126.AND.group_desc_sg(4)==135) THEN 
                     sg_number=32
                  ENDIF 
               ELSE 
                  IF (group_desc_sg(3)==34.AND.group_desc_sg(4)==134) THEN
                     sg_number=26
                  ELSEIF (group_desc_sg(3)==127.AND.group_desc_sg(4)==36) THEN 
                     aux_sg=6
                     sg_number=26
                  ELSEIF (group_desc_sg(3)==127.AND.group_desc_sg(4)==135) THEN
                     sg_number=29
                  ELSEIF (group_desc_sg(3)==126.AND.group_desc_sg(4)==134) THEN 
                     aux_sg=6
                     sg_number=29
                  ELSEIF (group_desc_sg(3)==34.AND.group_desc_sg(4)==136) THEN
                     s01(1)=s0(1,3)
                     s01(2)=s0(2,4)
                     s01(3)=0.0_DP
                     sg_number=31
                  ELSEIF (group_desc_sg(3)==128.AND.group_desc_sg(4)==36) THEN 
                     s01(1)=s0(1,3)
                     s01(2)=s0(2,4)
                     s01(3)=0.0_DP
                     aux_sg=6
                     sg_number=31
                  ELSEIF (group_desc_sg(3)==128.AND.group_desc_sg(4)==135) THEN
                     sg_number=33
                  ELSEIF (group_desc_sg(3)==126.AND.group_desc_sg(4)==136) THEN 
                     aux_sg=6
                     sg_number=33
                  ENDIF 
               ENDIF                 
            ELSEIF (code_group_ext==58) THEN
               idir=3
               is_proper=(group_desc_sg(2)==2) 
               ftrs1(1)=0.0_DP
               ftrs2(2)=0.0_DP
!
!   the origin on the C_2 axis
!
               s01(:)=s0(:,2)
               s01(3)=0.0_DP
               aux_sg=1
               IF (is_proper) THEN
                  IF (group_desc_sg(3)==36.AND.group_desc_sg(4)==35) THEN 
                     sg_number=25
                  ELSEIF (group_desc_sg(3)==135.AND.group_desc_sg(4)==131) THEN
                     sg_number=27
                  ELSEIF (group_desc_sg(3)==36.AND.group_desc_sg(4)==130) THEN
                     sg_number=28
                  ELSEIF (group_desc_sg(3)==134.AND.group_desc_sg(4)==35) THEN 
                     aux_sg=2
                     sg_number=28
                  ELSEIF (group_desc_sg(3)==136.AND.group_desc_sg(4)==131) THEN 
                     sg_number=30
                  ELSEIF(group_desc_sg(3)==135.AND.group_desc_sg(4)==132) THEN 
                     aux_sg=2
                     sg_number=30
                  ELSEIF (group_desc_sg(3)==136.AND.group_desc_sg(4)==132) THEN 
                     sg_number=34
                  ELSEIF (group_desc_sg(3)==134.AND.group_desc_sg(4)==130) THEN 
                     sg_number=32
                  ENDIF 
               ELSE 
                  IF (group_desc_sg(3)==36.AND.group_desc_sg(4)==131) THEN
                     sg_number=26
                  ELSEIF (group_desc_sg(3)==135.AND.group_desc_sg(4)==35) THEN 
                     aux_sg=2
                     sg_number=26
                  ELSEIF (group_desc_sg(3)==135.AND.group_desc_sg(4)==130) THEN
                     sg_number=29
                  ELSEIF (group_desc_sg(3)==134.AND.group_desc_sg(4)==131) THEN 
                     aux_sg=2
                     sg_number=29
                  ELSEIF (group_desc_sg(3)==36.AND.group_desc_sg(4)==132) THEN
                     s01(1)=s0(1,3)
                     s01(2)=s0(2,4)
                     s01(3)=0.0_DP
                     sg_number=31
                  ELSEIF (group_desc_sg(3)==132.AND.group_desc_sg(4)==35) THEN 
                     s01(1)=s0(1,3)
                     s01(2)=s0(2,4)
                     s01(3)=0.0_DP
                     aux_sg=2
                     sg_number=31
                  ELSEIF (group_desc_sg(3)==136.AND.group_desc_sg(4)==130) THEN
                     sg_number=33
                  ELSEIF (group_desc_sg(3)==134.AND.group_desc_sg(4)==132) THEN 
                     aux_sg=2
                     sg_number=33
                  ENDIF 
               ENDIF                 
            ENDIF 
         ELSEIF (ibrav==9) THEN
            IF (code_group_ext==54) THEN
               aux_sg=4
               IF (group_desc_sg(3)==35.OR.group_desc_sg(3)==132) THEN 
                  IF (group_desc_sg(4)==34.OR.group_desc_sg(4)==128) THEN
                     sg_number=38
                  ELSE
                     sg_number=39
                  ENDIF
               ELSE
                  IF (group_desc_sg(4)==34.OR.group_desc_sg(4)==128) THEN
                     sg_number=40
                  ELSE 
                     sg_number=41
                  ENDIF
               ENDIF
            ELSEIF (code_group_ext==56) THEN
               aux_sg=5
               IF (group_desc_sg(3)==34.OR.group_desc_sg(3)==128) THEN 
                  IF (group_desc_sg(4)==36.OR.group_desc_sg(4)==134) THEN
                     sg_number=38
                  ELSE
                     sg_number=40
                  ENDIF
               ELSE
                  IF (group_desc_sg(4)==36.OR.group_desc_sg(4)==134) THEN
                     sg_number=39
                  ELSE 
                     sg_number=41
                  ENDIF
               ENDIF
            ELSEIF (code_group_ext==58) THEN
               IF (group_desc_sg(2) == 65) THEN
                  sg_number=36
               ELSEIF (ABS(ftc(3,3))<eps1.AND.ABS(ftc(3,4))<eps1) THEN
                  sg_number=35
               ELSE 
                  sg_number=37
               ENDIF
            ENDIF
         ELSEIF (ibrav==91) THEN
!
!   find the mirror perpendicular to x mirror1 and to y mirror2
!
            IF (code_group_ext==54) THEN
               IF (group_desc_sg(2)/=4) THEN
                  sg_number=36
                  IF (group_desc_sg(3)==35 .OR. group_desc_sg(3)==131) THEN
                     aux_sg=3
                  ELSE
                     aux_sg=4
                  END IF
               ELSE
                  aux_sg=2
                  IF (group_desc_sg(3)==35 .OR. group_desc_sg(3)==131) THEN
                     sg_number=35
                  ELSE
                     sg_number=36
                  ENDIF 
               END IF
            ELSEIF (code_group_ext==56) THEN
               aux_sg=6
               IF (group_desc_sg(3)==34.OR.group_desc_sg(3)==127) THEN
                  IF (group_desc_sg(4)==36.OR.group_desc_sg(4)==136) THEN
                     sg_number=38
                  ELSE
                     sg_number=39
                  END IF
               ELSE
                  IF (group_desc_sg(4)==36.OR.group_desc_sg(4)==136) THEN
                     sg_number=40
                  ELSE
                     sg_number=41
                  END IF
               END IF
            ELSEIF (code_group_ext==58) THEN
               imirror1=3
               imirror2=4
               IF (group_desc_sg(imirror1)==36.OR.&
                                        group_desc_sg(imirror1)==136) THEN
                  IF (group_desc_sg(imirror2)==35.OR.&
                                        group_desc_sg(imirror2)==131) THEN
                       sg_number=38
                  ELSE
!
!    Ann2_1
!
                       sg_number=40
                  ENDIF
               ELSE
                  IF (group_desc_sg(imirror2)==35.OR.&
                                        group_desc_sg(imirror2)==131) THEN
!
!    Aec2_1
!
                     sg_number=39
                  ELSE
!
!    Aen2_1
!
                     sg_number=41
                  END IF
               END IF
            END IF
         ELSEIF (ibrav==10) THEN
            IF (code_group_ext==54) aux_sg=2
            IF (code_group_ext==56) aux_sg=3
            IF (code_group_ext==54) aux_sg=1
            IF (group_desc_sg(3)==137.AND.group_desc_sg(4)==133) THEN
               sg_number=43
            ELSE 
               sg_number=42
            ENDIF
         ELSEIF (ibrav==11) THEN
!
!   mirror1 is the mirror perpendicular to x, 
!   mirror2 is the mirror perpendicular to y, 
!
            IF (code_group_ext==54) aux_sg=2
            IF (code_group_ext==56) aux_sg=3
            IF (code_group_ext==54) aux_sg=1
            imirror1=3
            imirror2=4
            IF (group_desc_sg(imirror2)==35.OR.&
                             group_desc_sg(imirror2)==132) THEN
!
!   mirror2 is m or n
!
               sg_number=44
            ELSEIF (group_desc_sg(imirror1)==36.OR.&
                                     group_desc_sg(imirror1)==136) THEN
!
!  mirror1 is m or n
!
               sg_number=46
            ELSE
               sg_number=45
            ENDIF
         ENDIF
      CASE(13)
!
!  C_3v
!
!  the origin on axis 3
!
         s01(:)=s0(:,2)
         type_group=0
         IF (code_group_ext==72) type_group=1

         IF (ibrav==4) THEN
            IF (type_group==1) THEN
               IF (group_desc_sg(4)==36.OR.group_desc_sg(4)==134) THEN
                  sg_number=156
               ELSE
                  sg_number=158
               END IF
            ELSE     
               IF (group_desc_sg(5)==35.OR.group_desc_sg(5)==130) THEN
                  sg_number=157
               ELSE
                  sg_number=159
               ENDIF
            ENDIF
         ELSEIF (ibrav==5) THEN
            IF (group_desc_sg(4)==36.OR.group_desc_sg(4)==134) THEN
               sg_number=160
            ELSE 
               sg_number=161
            ENDIF
         ENDIF     
      CASE(14)
!
!  C_4v
!
         s01(:)=s0(:,2)
         irot90=2
         imirror=5

         IF (ibrav==6) THEN
            IF (group_desc_sg(irot90)==8) THEN
               IF (group_desc_sg(imirror)==36 ) THEN
                  sg_number=99
               ELSEIF (group_desc_sg(imirror)==134 ) THEN
                  sg_number=100
               ELSEIF (group_desc_sg(imirror)==135 ) THEN
                  sg_number=103
               ELSEIF (group_desc_sg(imirror)==136 ) THEN
                  sg_number=104
               ENDIF
            ELSE
               IF (group_desc_sg(imirror)==36 ) THEN
                  sg_number=105
               ELSEIF (group_desc_sg(imirror)==134 ) THEN
                  sg_number=106
               ELSEIF (group_desc_sg(imirror)==135 ) THEN
                  sg_number=101
               ELSEIF (group_desc_sg(imirror)==136 ) THEN
                  sg_number=102
               ENDIF
            END IF
         ELSEIF (ibrav==7) THEN
            IF (group_desc_sg(irot90)==8.OR.group_desc_sg(irot90)==74) THEN
               IF (group_desc_sg(imirror)==36.OR.&
                                      group_desc_sg(imirror)==136) THEN
                  sg_number=107
               ELSE
                  sg_number=108
               END IF
            ELSE
               IF (group_desc_sg(imirror)==36.OR.&
                                      group_desc_sg(imirror)==136) THEN
                  sg_number=109
                  s01(:) = s0(:,3)
               ELSE
                  sg_number=110
                  s01(:) = s0(:,3)
               ENDIF
            ENDIF
         ENDIF
      CASE(15)
!
!  C_6v
!
         s01(:)=s0(:,2)
         irot60=2
         imirror=7
         IF (group_desc_sg(irot60)==25) THEN
            IF (group_desc_sg(imirror)==36.OR.group_desc_sg(imirror)==134) THEN
               sg_number=183
            ELSE
               sg_number=184
            ENDIF
         ELSE
            IF (group_desc_sg(imirror)==36.OR.group_desc_sg(imirror)==134) THEN
               sg_number=186
            ELSE
               sg_number=185
            ENDIF
         ENDIF
      CASE(16)
!
!  C_2h
!
!
!   origin on the inversion center
!
          s01(:)=s0(:,3)
          IF (ibrav==-13) THEN
!
!  B112/m     unique c
!     2_1/a
!
             aux_sg=3
             IF (group_desc_sg(4)==35.OR.group_desc_sg(4)==130) THEN
                sg_number=12
             ELSE 
                sg_number=15
             ENDIF
          ELSEIF (ibrav==13) THEN
!
!  C12/m1   unique b
!    2_1/a
!
             aux_sg=1
             IF (group_desc_sg(4)==34.OR.group_desc_sg(4)==126) THEN
                sg_number=12
             ELSE 
                sg_number=15
             ENDIF
          ELSEIF(ABS(ibrav)==12) THEN
             IF (is_symmorphic) THEN
                IF (ibrav>0) THEN
                   aux_sg=3
                ELSE
                   aux_sg=1
                END IF
                sg_number=10
             ELSE

                irot180=2
                imirror=4

                IF (ibrav==-12) THEN
!
!   b unique
!
                   aux_sg=1
                   IF (group_desc_sg(irot180)==3) THEN
                      sg_number=13
                   ELSE
                      IF (group_desc_sg(imirror)==35) THEN
                         sg_number=11
                      ELSE
                         sg_number=14
                      ENDIF
                   ENDIF
                ELSE
                   aux_sg=3
                   IF (group_desc_sg(irot180)==2) THEN
                      sg_number=13
                   ELSE
                      IF (group_desc_sg(imirror)==34) THEN
                         sg_number=11
                      ELSE
                         sg_number=14
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
      CASE(17)
!
!  C_3h
!
!  origin on the inversion point of -6.
!
         s01(:)=s0(:,2)
         sg_number=174
      CASE(18)
!
!  C_4h
!
!
!  origin on inversion center
!
         s01(:)=s0(:,5)
         IF (ibrav==6) THEN
            IF (is_symmorphic) THEN
               sg_number=83
            ELSE
               irot90=2
               imirror=7
               IF (ABS(ft(3,irot90))<eps1) THEN
                  sg_number=85
!
!   origin choice 1 on the inversion center of -4
!
                  s01(:)=s0(:,6)
!
!   origin choice 2 on inversion
!
                  s02(:)=s0(:,5)
               ELSEIF (group_desc_sg(irot90) == 74) THEN
                  IF (ABS(ft(1,imirror))<eps1.AND.ABS(ft(2,imirror))<eps1) THEN
                     sg_number=84
                  ELSE
                     sg_number=86
!
!   origin choice 1 on the inversion center of -4
!
                     s01(:)=s0(:,6)
!
!   origin choice 2 on inversion
!
                     s02(:)=s0(:,5)
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF(ibrav==7) THEN
            IF (group_desc_sg(2)==8.OR.group_desc_sg(2)==74) THEN
               sg_number=87
            ELSE
               sg_number=88
!
!   origin choice 1 on the inversion center of -4
!
               s01(:)=s0(:,6)
!
!   origin choice 2 on inversion
!
               s02(:)=s0(:,5)
            ENDIF
         ENDIF
      CASE(19)
!
!  C_6h
!
!  origin on the inversion center
!
         s01(:)=s0(:,7)
         IF (is_symmorphic) THEN
            sg_number=175
         ELSE
            sg_number=176
         ENDIF
      CASE(20)
!
!  D_2h
!
!
!   by default the origin is on the inversion point
!
         s01(:)=s0(:,5)
!
!   count the number of improper axis
!
         nsa=0
         IF (group_desc_sg(2) /= 2) nsa=nsa+1
         IF (group_desc_sg(3) /= 4) nsa=nsa+1
         IF (group_desc_sg(4) /= 3) nsa=nsa+1

         IF (ibrav==8) THEN
            IF (nsa==0) THEN
               IF (is_symmorphic) THEN
                  sg_number=47
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==132) THEN
!
!  origin choice 1 at the intesection of the 2 axes
!
                        s01(:)=s0(:,2)
                        s01(3)=s0(3,3)
!
!   origin choice 2 on inversion
!
                        s02(:)=s0(:,5)
                  sg_number=48
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==131) THEN
                  sg_number=49
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=2
                  sg_number=49
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=3
                  sg_number=49
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==130) THEN
                  sg_number=50
!
!  origin choice 1 at the intesection of the 2z axis and its perpendicular 
!  mirror
!
                   s01(:)=s0(:,2)
                   s01(3)=s0(3,6)
!
!   origin choice 2 on inversion
!
                   s02(:)=s0(:,5)
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=2
                  sg_number=50
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=3
                  sg_number=50
               END IF   
            ELSEIF (nsa==1) THEN
               IF (group_desc_sg(6)==126.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==35) THEN
                  sg_number=51
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=2
                  sg_number=51
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=3
                  sg_number=51
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=4
                  sg_number=51
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=5
                  sg_number=51
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=6
                  sg_number=51
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==132) THEN
                  sg_number=52
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=2
                  sg_number=52
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=3
                  sg_number=52
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=4
                  sg_number=52
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=5
                  sg_number=52
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=6
                  sg_number=52
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==132) THEN
                  sg_number=53
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=2
                  sg_number=53
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=3
                  sg_number=53
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=4
                  sg_number=53
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=5
                  sg_number=53
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=6
                  sg_number=53
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==131) THEN
                  sg_number=54
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=2
                  sg_number=54
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=3
                  sg_number=54
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=4
                  sg_number=54
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=5
                  sg_number=54
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=6
                  sg_number=54
               END IF   

            ELSEIF (nsa==2) THEN
               IF (group_desc_sg(6)==34.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==130) THEN
                  sg_number=55
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=2
                  sg_number=55
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=3
                  sg_number=55
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==131) THEN
                  sg_number=56
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=2
                  sg_number=56
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=3
                  sg_number=56
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==131) THEN
                  sg_number=57
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=2
                  sg_number=57
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=3
                  sg_number=57
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=4
                  sg_number=57
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=5
                  sg_number=57
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=6
                  sg_number=57
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==132) THEN
                  sg_number=58
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=2
                  sg_number=58
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=3
                  sg_number=58
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==35) THEN
                  sg_number=59
!
!  origin choice 1 at the intesection of the 2z axis and its perpendicular 
!  mirror
!
                  s01(:)=s0(:,2)
                  s01(3)=s0(3,6)
!
!   origin choice 2 on inversion
!
                  s02(:)=s0(:,5)
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=2
                  sg_number=59
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=3
                  sg_number=59
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==131) THEN
                  sg_number=60
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=2
                  sg_number=60
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=3
                  sg_number=60
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=4
                  sg_number=60
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=5
                  sg_number=60
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=6
                  sg_number=60
               ENDIF
            ELSE
               IF (group_desc_sg(6)==126.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==131) THEN
                  sg_number=61
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=2
                  sg_number=61
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=3
                  sg_number=61
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=4
                  sg_number=61
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=5
                  sg_number=61
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=6
                  sg_number=61
               ELSEIF (group_desc_sg(6)==126.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==35) THEN
                  sg_number=62
               ELSEIF (group_desc_sg(6)==127.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=2
                  sg_number=62
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==134.AND. &
                       group_desc_sg(8)==132) THEN
                  aux_sg=3
                  sg_number=62
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==135.AND. &
                       group_desc_sg(8)==35) THEN
                  aux_sg=4
                  sg_number=62
               ELSEIF (group_desc_sg(6)==128.AND.group_desc_sg(7)==36.AND. &
                       group_desc_sg(8)==131) THEN
                  aux_sg=5
                  sg_number=62
               ELSEIF (group_desc_sg(6)==34.AND.group_desc_sg(7)==136.AND. &
                       group_desc_sg(8)==130) THEN
                  aux_sg=6
                  sg_number=62
               ENDIF
            ENDIF
         ELSEIF (ibrav==9) THEN
            IF (nsa>0) THEN
               IF (group_desc_sg(6)==34.OR.group_desc_sg(6)==128) THEN
                  sg_number=63 
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==134) THEN
                     aux_sg=1
                  ELSE
                     aux_sg=2
                  ENDIF
               ELSE
                  sg_number=64
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==134) THEN
                     aux_sg=1
                  ELSE
                     aux_sg=2
                  ENDIF
               ENDIF
            ELSE
               IF (group_desc_sg(6)==34.OR.group_desc_sg(6)==128) THEN
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==134) THEN
                     sg_number=65
                  ELSE
                     sg_number=66
                  ENDIF
               ELSE
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==134) THEN
                     sg_number=67
                  ELSE
                     sg_number=68
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==91) THEN
            IF (nsa>0) THEN
               IF (group_desc_sg(6)==126.OR.group_desc_sg(6)==128) THEN
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==136) THEN
                     aux_sg=3
                     sg_number=63 
                  ELSE
                     aux_sg=3
                     sg_number=64 
                  ENDIF
               ELSEIF (group_desc_sg(6)==34.OR.group_desc_sg(6)==127) THEN
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==136) THEN
                     aux_sg=4
                     sg_number=63 
                  ELSE
                     aux_sg=4
                     sg_number=64 
                  ENDIF
               ENDIF
            ELSE
               IF (group_desc_sg(6)==34.OR.group_desc_sg(6)==127) THEN
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==136) THEN
                     aux_sg=2
                     sg_number=65 
                  ELSE
                     aux_sg=2
                     sg_number=67 
                  ENDIF
               ELSE
                  IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==136) THEN
                     aux_sg=2
                     sg_number=66 
                  ELSE
                     aux_sg=2
                     sg_number=68 
                  ENDIF
               ENDIF
            END IF
         ELSEIF (ibrav==10) THEN
            IF (group_desc_sg(7)==137.AND.group_desc_sg(8)==133) THEN
               s01(:)=(s0c(:,2)+s0c(:,3)+s0c(:,4))/3.0_DP
               CALL cryst_to_cart(1,s01,bg,-1)
               s02(:)=s0(:,5)
               sg_number=70
            ELSE
               sg_number=69
            ENDIF
         ELSEIF (ibrav==11) THEN
!
!   count how many mirrors have m or n
!
            nsa=0
            IF (group_desc_sg(6)==34.OR.group_desc_sg(6)==128) nsa=nsa+1
            IF (group_desc_sg(7)==36.OR.group_desc_sg(7)==136) nsa=nsa+1
            IF (group_desc_sg(8)==35.OR.group_desc_sg(8)==132) nsa=nsa+1
            IF (nsa==0) THEN
               sg_number=73
            ELSEIF (nsa==1) THEN
               sg_number=72      
               IF (group_desc_sg(6)==34.OR.group_desc_sg(6)==128) THEN
                  aux_sg=1
               ELSEIF (group_desc_sg(7)==36.OR.group_desc_sg(7)==136) THEN
                  aux_sg=2
               ELSE
                  aux_sg=3
               ENDIF
            ELSEIF (nsa==2) THEN
               sg_number=74      
               IF (group_desc_sg(6)==126.OR.group_desc_sg(6)==127) THEN
                  aux_sg=1
               ELSEIF (group_desc_sg(7)==134.OR.group_desc_sg(7)==135) THEN
                  aux_sg=2
               ELSE
                  aux_sg=3
               ENDIF
            ELSEIF (nsa==3) THEN
               sg_number=71      
            ENDIF
         ENDIF
      CASE(21)
!
!  D_3h
!
!  origin on the inversion center of -6
!
         s01(:)=s0(:,2)

         IF (code_group_ext==107) THEN
!
!  x is perpendicular to a mirror
!
            type_group=1
            imirror=7
         ELSE
            type_group=2
            imirror=10
         ENDIF
         IF (type_group==1) THEN
            IF (group_desc_sg(imirror)==36.OR.&
                               group_desc_sg(imirror)==134) THEN
               sg_number=187
            ELSE
!
!  origin on 3c2
!
               s01(:)=s0(:,3)
               s01(3)=s0(3,10)
               sg_number=188
            ENDIF
         ELSE
            IF (group_desc_sg(imirror)==35.OR.&
                               group_desc_sg(imirror)==130) THEN
               sg_number=189
            ELSE
!
!  origin on 32c
!
               s01(:)=s0(:,3)
               s01(3)=s0(3,7)
               sg_number=190
            ENDIF
         ENDIF
      CASE(22)
!
!  D_4h
!
         s01(:)=s0(:,9)
         irot90=2
         irot180=5
         imirror=13
         IF (ibrav==6) THEN
            IF (group_desc_sg(irot90)==8) THEN
               IF (group_desc_sg(irot180)==4) THEN
                  IF (group_desc_sg(imirror)==36) THEN
                     sg_number=123
                  ELSEIF (group_desc_sg(imirror)==135) THEN
                     sg_number=124
                  ELSEIF (group_desc_sg(imirror)==134) THEN
                     sg_number=125
!
!  two origins origin 1 on the inversion point of -4, origin 2 on the inversion
!  center
!
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  ELSEIF (group_desc_sg(imirror)==136) THEN
                     sg_number=126
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  ENDIF
               ELSE
                  IF (group_desc_sg(imirror)==36) THEN
                     sg_number=129
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  ELSEIF (group_desc_sg(imirror)==135) THEN
                     sg_number=130
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  ELSEIF (group_desc_sg(imirror)==134) THEN
                     sg_number=127
                  ELSEIF (group_desc_sg(imirror)==136) THEN
                     sg_number=128
                  ENDIF
               ENDIF
            ELSE
               IF (group_desc_sg(irot180)==4) THEN
                  IF (group_desc_sg(imirror)==36) THEN
                     sg_number=131
                  ELSEIF (group_desc_sg(imirror)==135) THEN
                     sg_number=132
                  ELSEIF (group_desc_sg(imirror)==134) THEN
                     sg_number=133
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  ELSEIF (group_desc_sg(imirror)==136) THEN
                     sg_number=134
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  END IF
               ELSE
                  IF (group_desc_sg(imirror)==36) THEN
                     sg_number=137
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  ELSEIF (group_desc_sg(imirror)==135) THEN
                     sg_number=138
                     s01(:)=s0(:,10)
                     s02(:)=s0(:,9)
                  ELSEIF (group_desc_sg(imirror)==134) THEN
                     sg_number=135
                  ELSEIF (group_desc_sg(imirror)==136) THEN
                     sg_number=136
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (group_desc_sg(irot90)==8.OR.group_desc_sg(irot90)==74) THEN
               IF (group_desc_sg(imirror)==36.OR.&
                                    group_desc_sg(imirror)==136) THEN
                  sg_number=139
               ELSE
                  sg_number=140
               ENDIF
            ELSE
               IF (group_desc_sg(imirror)==36.OR.&
                                    group_desc_sg(imirror)==136) THEN
                  sg_number=141
                  s01(:)=s0(:,10)
                  s02(:)=s0(:,9)
               ELSE
                  sg_number=142
                  s01(:)=s0(:,10)
                  s02(:)=s0(:,9)
               ENDIF
            ENDIF
         ENDIF 
      CASE(23)
!
!  D_6h
!
!  origin on the inversion center of the -3_z operation
!
         s01(:)=s0(:,18)
         irot60=2
         imirror=19 ! perpendicular to x
         IF (group_desc_sg(irot60)==25) THEN
            IF (group_desc_sg(imirror)==36.OR.group_desc_sg(imirror)==134) THEN
               sg_number=191
            ELSE
               sg_number=192
            ENDIF
         ELSE
            IF (group_desc_sg(imirror)==36.OR.group_desc_sg(imirror)==134) THEN
               sg_number=194
            ELSE
               sg_number=193
            ENDIF
         ENDIF
      CASE(24)
!
!  D_2d
!
!
!  origin on the inversion point of -4
!
         s01(:)=s0(:,2)
         type_group=0
!
!  if the axis of order 2 make angles of 0,90,180,270 with x (type_group=1) 
!  if it makes angles 45,135,225,315 (type_group=0)
!
         IF (code_group_ext==112) type_group=1

         IF (type_group==1) THEN
!
!  the x axis is a two-fold rotation axis and there is a mirror 
!  perpendicular to (1,1,0). Search it
!
            irot180=5   ! x axis
            imirror1=6
         ELSE
!
!  There is a mirror perpendicular to the x axis. 
!
            imirror=5
         ENDIF

         IF (ibrav==6) THEN
            IF (type_group==1) THEN
               IF (group_desc_sg(irot180)==4) THEN
                  IF (group_desc_sg(imirror1)==139.OR.&
                                    group_desc_sg(imirror1)==140) THEN
                     sg_number=112
                  ELSE
                     sg_number=111
                  ENDIF
               ELSE
                  IF (group_desc_sg(imirror1)==139.OR.&
                                    group_desc_sg(imirror1)==140) THEN
                     sg_number=114
                  ELSE
                     sg_number=113
                  ENDIF
               ENDIF
            ELSE
               IF (group_desc_sg(imirror)==36) THEN
                  sg_number=115
               ELSEIF (group_desc_sg(imirror)==135) THEN
                  sg_number=116
               ELSEIF (group_desc_sg(imirror)==134) THEN
                  sg_number=117
               ELSEIF (group_desc_sg(imirror)==136) THEN
                  sg_number=118
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (type_group==0) THEN
               IF (group_desc_sg(imirror)==36.OR.&
                                    group_desc_sg(imirror)==136) THEN
                  sg_number=119
               ELSE
                  sg_number=120
               ENDIF
            ELSE
               IF (group_desc_sg(imirror1)==141) THEN
                  sg_number=122
               ELSE
                  sg_number=121
               ENDIF
            ENDIF
         ENDIF
      CASE(25)
!
!  D_3d
!
!  origin on the inversion point of the i3z operation
!
         s01(:)=s0(:,8)
         IF (code_group_ext==118) THEN
            type_group=1
         ELSEIF (code_group_ext==119) THEN
            type_group=2
         ELSE
            CALL errore('find_space_group','D_3d space group not here',1)
         ENDIF

         IF (ibrav==4) THEN
            IF (type_group==1) THEN
!
!   if the mirror perpendicular to x is proper or b it is 164
!
               IF (group_desc_sg(10)==36.OR.group_desc_sg(10)==134) THEN
                  sg_number=164
               ELSE
                  sg_number=165
               ENDIF
            ELSE
!
!   if all the mirror perpendicular to y is proper or a it is 162
!
               IF (group_desc_sg(11)==35.OR.group_desc_sg(10)==130) THEN
                  sg_number=162
               ELSE
                  sg_number=163
               ENDIF
            ENDIF
         ELSEIF (ibrav==5) THEN
             IF (group_desc_sg(10)==36.OR.group_desc_sg(10)==134) THEN
                sg_number=166
             ELSE
                sg_number=167
             ENDIF
         ENDIF
      CASE(26)
!
! S_4
!
!
!   origin on -4
!
         s01(:)=s0(:,2)
         IF (ibrav==6) THEN
            sg_number=81
         ELSEIF (ibrav==7) THEN
            sg_number=82
         ENDIF
      CASE(27)
!
!  S_6 (-3)
!
!  origin on the inversion center of -3
!
         s01(:)=s0(:,5)
         IF (ibrav==4) THEN
            sg_number=147
         ELSEIF (ibrav==5) THEN
            sg_number=148
         ENDIF
      CASE(28)
!
!   T
!
         s01(:)=s0(:,2)
         s01(3)=s0(3,3)
         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=195
             ELSE
!
!   origin not programmed
!
                s01=1000._DP
                sg_number=198
             ENDIF
         ELSEIF (ibrav==2) THEN
             sg_number=196
         ELSEIF (ibrav==3) THEN
!
!   If the three 2-fold axis intersect in a point or can be made to intersect by
!   a direct lattice translation, the group is 197 otherwise it is 199
!
            IF (check_intersection(ftperpc(1,2), ftperpc(1,3), &
                                                 ftperpc(1,4))) THEN
                sg_number=197
             ELSE
!
!   origin not programmed
!
                s01=1000._DP
                sg_number=199
             ENDIF
         ENDIF
      CASE(29)
!
!  T_h
!

         s01(:)=s0(:,17)
         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=200
             ELSE
                IF (group_desc_sg(14)==128) THEN
!
!   the mirror i2z is n
!
                   sg_number=201
!
!   origin choice 1 at 23
!
                   s01(:)=s0(:,2)
                   s01(3)=s0(3,5)
!
!  origin choice 2 at the center of -3
!
                   s02(:)=s0(:,17)
                ELSEIF (group_desc_sg(14)==126.OR.group_desc_sg(14)==127) THEN
                   sg_number=205
                ELSE
                   CALL errore('find_space_group',&
                               'T_h operation not recognized',1)
                ENDIF
             ENDIF
         ELSEIF (ibrav==2) THEN
             IF (group_desc_sg(14)==129) THEN
                sg_number=203
!
!   origin choice 1 at 23
!
                   s01(:)=s0(:,2)
                   s01(3)=s0(3,5)
!
!  origin choice 2 at the center of -3
!
                   s02(:)=s0(:,17)
             ELSE
                sg_number=202
             ENDIF
         ELSEIF (ibrav==3) THEN
             IF (group_desc_sg(14)==34.OR.group_desc_sg(14)==128) THEN
                sg_number=204
             ELSE
                sg_number=206
             ENDIF
         ENDIF
      CASE(30)
!
!  T_d
!
         s01(:)=s0(:,18)
         IF (ibrav==1) THEN
             IF (group_desc_sg(24)==139.OR.group_desc_sg(24)==140) THEN
                s01(:)=s0(:,4)
                s01(3)=s0(3,2)
                sg_number=218
             ELSE
                sg_number=215
             ENDIF
         ELSEIF (ibrav==2) THEN
!
!  one of the two operations i2xy or i2x-y (or both) is a real mirror in 216
!  they are all glide planes in 219. (Is this sufficient or should we
!  check the six diagonal mirrors?)
!  Same problem to distinguish 225 - 226 and 227 - 228
!
             IF (group_desc_sg(24)==37.OR.group_desc_sg(23)==38) THEN
                sg_number=216
             ELSE
                s01(:)=s0(:,4)
                s01(3)=s0(3,2)
                sg_number=219
             ENDIF
         ELSEIF (ibrav==3) THEN
             IF (group_desc_sg(24)==141) THEN
!
!    origin to be programmed
!
                s01=1000.0_DP
                sg_number=220
             ELSE
                sg_number=217
             ENDIF
         ENDIF
      CASE(31)
!
!  O
!
!   origin at the intersection of 2z and 2x
!
         s01(:)=s0(:,4)
         s01(3)=s0(3,2)
         irot90=14
         IF (ibrav==1) THEN
            IF (group_desc_sg(irot90)==16) THEN
               sg_number=207
            ELSEIF (group_desc_sg(irot90)==90) THEN
!
!    4x is 4_2 in 208, 4_1 in 213, 4_3 in 212
!
                sg_number=208
            ELSEIF (group_desc_sg(irot90)==89) THEN
                sg_number=213
!
!  origin to be programmed
!
                s01=1000.0_DP
            ELSEIF (group_desc_sg(irot90)==91) THEN
                sg_number=212
!
!  origin to be programmed
!
                s01=1000.0_DP
            ENDIF
         ELSEIF (ibrav==2) THEN
            IF (group_desc_sg(irot90)==16.OR.group_desc_sg(irot90)==90) THEN
               sg_number=209
            ELSE
               sg_number=210
            ENDIF
         ELSEIF (ibrav==3) THEN
            IF (group_desc_sg(irot90)==16.OR.group_desc_sg(irot90)==90) THEN
               sg_number=211
            ELSE
!
!  origin to be programmed
!
               s01=1000.0_DP
               sg_number=214
            ENDIF
         ENDIF
      CASE(32)
!
!  O_h
!
!   origin on the inversion point of i3xyz
!
         s01(:)=s0(:,29)
         IF (ibrav==1) THEN
            IF (group_desc_sg(18)==8) THEN
               IF (group_desc_sg(28)==34) THEN
                  sg_number=221
               ELSE
                  sg_number=222
!
!   origin choice 1 at 432
!
                  s01(:)=s0(:,18)
                  s01(3)=s0(3,4)
!
!  origin choice 2 at the inversion point of -3
!
                  s02(:)=s0(:,29)
               ENDIF
            ELSE
               IF (group_desc_sg(28)==34) THEN
                  sg_number=223
               ELSE
                  sg_number=224
!
!  origin choice 1 at -43m
!
                  s01(:)=s0(:,18)
                  s01(3)=s0(3,28)
!
!  origin choice 2 at the inversion point of -3
!
                  s02(:)=s0(:,29)
               ENDIF
            ENDIF
         ELSEIF (ibrav==2) THEN
            IF (group_desc_sg(18)==8.OR.group_desc_sg(18)==74) THEN
!
!  See above 216 - 215
!
               IF (group_desc_sg(48)==37.OR.group_desc_sg(47)==38) THEN
                  sg_number=225
               ELSE
                  sg_number=226
               ENDIF 
            ELSE
!
!  See above 216 - 215
!
               IF (group_desc_sg(48)==37.OR.group_desc_sg(47)==38) THEN
                  sg_number=227
!
!   origin choice 2 (on the inversion center)
!
                  s02(:)=s0(:,25)
!
!   origin choice 1 (on the -43m intersection point). We choose the
!   inversion point of the i4z operation
!
                  s01(:)=s0(:,42)
               ELSE
                  sg_number=228
!
!   origin choice 2 (on the inversion center)
!
                  s02(:)=s0(:,25)
!
!   origin choice 1 (on the -43m intersection point). We choose the
!   inversion point of the i4z operation
!
                  s01(:)=s0(:,42)
               ENDIF 
            ENDIF
         ELSEIF (ibrav==3) THEN
            IF (group_desc_sg(18)==8.OR.group_desc_sg(18)==74) THEN
               sg_number=229
            ELSE
               sg_number=230
            ENDIF
         ENDIF
  END SELECT

  IF (verbosity) WRITE(stdout,'(/,5x,"Space group number",i5)') sg_number

  RETURN
  END SUBROUTINE find_space_group

!-----------------------------------------------------------------------
  SUBROUTINE transform_fcc_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
!
!  QE primitive vector choice definition
!
     tau_aux(1,ia) = - aux(1,ia) - aux(2,ia) + aux(3,ia)
     tau_aux(2,ia) =   aux(1,ia) + aux(2,ia) + aux(3,ia)
     tau_aux(3,ia) = - aux(1,ia) + aux(2,ia) - aux(3,ia)
!
! More conventional primitive vector choice a/2(011) a/2(101) a/2(110) 
!
!     tau_aux(1,ia) = - aux(1,ia) + aux(2,ia) + aux(3,ia)
!     tau_aux(2,ia) =   aux(1,ia) - aux(2,ia) + aux(3,ia)
!     tau_aux(3,ia) =   aux(1,ia) + aux(2,ia) - aux(3,ia)
  ENDDO

  RETURN
  END SUBROUTINE transform_fcc_axis

!-----------------------------------------------------------------------
  SUBROUTINE transform_bcc_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia) =   aux(1,ia) + aux(3,ia)
     tau_aux(2,ia) = - aux(1,ia) + aux(2,ia) 
     tau_aux(3,ia) = - aux(2,ia) + aux(3,ia)
  ENDDO

  RETURN
  END SUBROUTINE transform_bcc_axis

!-----------------------------------------------------------------------
  SUBROUTINE transform_fco_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia) =   aux(1,ia) - aux(2,ia) + aux(3,ia)
     tau_aux(2,ia) =   aux(1,ia) + aux(2,ia) - aux(3,ia)
     tau_aux(3,ia) = - aux(1,ia) + aux(2,ia) + aux(3,ia)
  ENDDO

  RETURN
  END SUBROUTINE transform_fco_axis

!-----------------------------------------------------------------------
  SUBROUTINE transform_ct_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia) =   aux(1,ia) - aux(2,ia) 
     tau_aux(2,ia) =   aux(2,ia) + aux(3,ia) 
     tau_aux(3,ia) = - aux(1,ia) + aux(3,ia)
  ENDDO

  RETURN
  END SUBROUTINE transform_ct_axis

!-----------------------------------------------------------------------
  SUBROUTINE transform_obco_c_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia) =   aux(1,ia) + aux(2,ia) 
     tau_aux(2,ia) = - aux(1,ia) + aux(2,ia) 
  ENDDO

  RETURN
  END SUBROUTINE transform_obco_c_axis

!-----------------------------------------------------------------------
  SUBROUTINE transform_monoclinic_center_c(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia) =   aux(1,ia) - aux(3,ia) 
     tau_aux(3,ia) =   aux(1,ia) + aux(3,ia) 
  ENDDO

  RETURN
  END SUBROUTINE transform_monoclinic_center_c

!-----------------------------------------------------------------------
  SUBROUTINE transform_monoclinic_center_b(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia) =   aux(1,ia) + aux(2,ia) 
     tau_aux(2,ia) = - aux(1,ia) + aux(2,ia) 
  ENDDO

  RETURN
  END SUBROUTINE transform_monoclinic_center_b

!-----------------------------------------------------------------------
  SUBROUTINE transform_obco_a_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(2,ia) =   aux(2,ia) + aux(3,ia) 
     tau_aux(3,ia) = - aux(2,ia) + aux(3,ia) 
  ENDDO

  RETURN
  END SUBROUTINE transform_obco_a_axis

!-----------------------------------------------------------------------
  SUBROUTINE transform_bco_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia)=aux(1,ia)+aux(3,ia)
     tau_aux(2,ia)=aux(2,ia)-aux(1,ia)
     tau_aux(3,ia)=aux(3,ia)-aux(2,ia)
  ENDDO

  RETURN
  END SUBROUTINE transform_bco_axis

!-----------------------------------------------------------------------
  SUBROUTINE transform_trigonal_axis(tau_aux,naux)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: naux
  REAL(DP), INTENT(INOUT) :: tau_aux(3,naux)

  REAL(DP) :: aux(3,naux)
  INTEGER :: ia

  aux=tau_aux
  DO ia=1,naux
     tau_aux(1,ia)=   aux(1,ia) - aux(2,ia) + aux(3,ia)
     tau_aux(2,ia)=   aux(2,ia) + aux(3,ia)
     tau_aux(3,ia)= - aux(1,ia) + aux(3,ia)
  ENDDO

  RETURN
  END SUBROUTINE transform_trigonal_axis

!-----------------------------------------------------------------------
  SUBROUTINE find_space_group_number(sgroup_name, group_code)
!-----------------------------------------------------------------------
!
!  This routine returns the space group number of a given group_name if
!  the name is recognized or zero if it is not.
!
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: sgroup_name
  INTEGER, INTENT(OUT) :: group_code
  
  INTEGER :: igroup, icount

  group_code=0
  DO igroup=1,230
     DO icount=1,mnsg
        IF (TRIM(sgroup_name)==TRIM(sg_names(icount,igroup)(1:13))) THEN
           group_code=igroup
           RETURN
        ENDIF 
     ENDDO
  ENDDO
  RETURN
  END SUBROUTINE find_space_group_number

!-----------------------------------------------------------------------
  SUBROUTINE find_space_group_names(group_code)
!-----------------------------------------------------------------------
!
!  This routine writes on output all the names that correspond to a given
!  space group.
!
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_code
  
  INTEGER :: icount

  DO icount=1,mnsg
     IF (TRIM(sg_names(icount,group_code)) /='') &
        WRITE(stdout,'(5x,a,12x,a)') TRIM(sg_names(icount,group_code)), &
                                    TRIM(add_info(icount,group_code))
  ENDDO

  RETURN
  END SUBROUTINE find_space_group_names

!-----------------------------------------------------------------------
  SUBROUTINE sg_name(group_code, aux_sg, output_name)
!-----------------------------------------------------------------------
!
!  This routine receives a space group code and gives as output the 
!  space group name
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_code, aux_sg
  CHARACTER(LEN=12), INTENT(OUT) :: output_name
  
  INTEGER igroup

  IF (group_code<1.OR.group_code>230) &
     CALL errore('sg_name','group_code unknown',1)

  output_name=sg_names(aux_sg,group_code)
  RETURN
  END SUBROUTINE sg_name

!-----------------------------------------------------------------------
  SUBROUTINE project_frac_tran(sr,ft,at,bg,ftperp,ftpar)
!-----------------------------------------------------------------------
!
!  This subroutine receives as input a space group operation in the
!  form of a rotation matrix sr (in cartesian coordinates) and a 
!  fractional translation (in crystal coordinates), the direct and
!  reciprocal lattice vectors and gives as output the parts of
!  the fractional translation perpendicular or parallel to the
!  symmetry element (rotation axis or mirror plane). The output
!  fractional translations are in crystal coordinates. For improper
!  rotations or inversion ftperp is equal to ft, while ftpar is set to zero.
!
  IMPLICIT NONE 

  REAL(DP), INTENT(IN) :: sr(3,3), ft(3), at(3,3), bg(3,3)
  REAL(DP), INTENT(OUT) :: ftperp(3), ftpar(3)

  INTEGER :: ts, tipo_sym
  REAL(DP) :: ax(3), ftr(3), ftr1(3), ftr2(3), ps

  ftperp=ft
  ftpar=0.0_DP
  ts=tipo_sym(sr)
!
!  for identity, inversion or improper rotation ftperp is ft and ftpar is zero
!
  IF (ts==1 .OR. ts==2 .OR. ts==6) RETURN
     
  IF (ts==3 .OR. ts==4) THEN
!
!  rotation
!
     CALL versor(sr,ax)
  ELSE IF (ts==5) THEN
!
!  mirror
!
     CALL mirror_axis(sr, ax)
  END IF

  ftr=ft
  CALL cryst_to_cart(1, ftr, at, 1)
 
!  WRITE(6,*) 'axis', ax(1), ax(2), ax(3)
!  WRITE(6,*) 'fract', ftr(1), ftr(2), ftr(3)
  ps = ftr(1) * ax(1) + ftr(2) * ax(2) + ftr(3) * ax(3)

  ftr1(:)= ps * ax(:)
  ftr2(:)= ftr(:) - ftr1(:)
  CALL cryst_to_cart(1, ftr1, bg, -1)
  CALL cryst_to_cart(1, ftr2, bg, -1)

  IF (ts==3 .OR. ts==4) THEN
     ftpar(:)=ftr1(:)
     ftperp(:)=ftr2(:)
  ELSEIF (ts==5) THEN
     ftpar(:)=ftr2(:)
     ftperp(:)=ftr1(:)
  ENDIF

  RETURN
  END SUBROUTINE project_frac_tran

!-----------------------------------------------------------------------
  SUBROUTINE shift_frac_tran(sr,ft,at,bg,s0,ftnew)
!-----------------------------------------------------------------------
!
!  This subroutine receives as input a space group operation in the
!  form of a rotation matrix sr (in cartesian coordinates) and a 
!  fractional translation (in crystal coordinates), a shift of the
!  origin s0 (in crystal coordinates) and gives as output the fractional
!  translation (in crystal coordinates) referred to the new origin
!
  USE kinds, ONLY : DP
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: sr(3,3), ft(3), at(3,3), bg(3,3), s0(3)
  REAL(DP), INTENT(OUT) :: ftnew(3)
  
  INTEGER :: ipol, jpol
  REAL(DP) :: s0r(3), s0rr(3)

  s0r=s0
  CALL cryst_to_cart(1, s0r, at, 1)
  
  s0rr=0.0_DP
  DO ipol=1,3
     DO jpol=1,3
        s0rr(ipol)=s0rr(ipol)+sr(ipol,jpol)*s0r(jpol)
     ENDDO
  ENDDO
  
  CALL cryst_to_cart(1, s0rr, bg, -1)

  ftnew(:)=ft(:) + s0rr(:) - s0(:)

  RETURN
  END SUBROUTINE shift_frac_tran

!-----------------------------------------------------------------------
  SUBROUTINE roto_shift_sg(sr,ft,at,bg,rot,s0,srnew,ftnew)
!-----------------------------------------------------------------------
!
!  This subroutine receives as input a space group operation in the
!  form of a rotation matrix sr (in cartesian coordinates) and a 
!  fractional translation (in crystal coordinates), a shift of the
!  origin s0 (in crystal coordinates) and a rotation matrix rot
!  and gives as output the rotation matrix and the fractional
!  translation (in crystal coordinates) referred to the new origin
!  and to the rotated reference system
!
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: sr(3,3), ft(3), at(3,3), bg(3,3), s0(3), rot(3,3)
  REAL(DP), INTENT(OUT) :: srnew(3,3), ftnew(3)

  REAL(DP)  :: ft_inter(3)
  INTEGER :: ipol, jpol
!
!  first apply the origin traslation
!
  CALL shift_frac_tran(sr,ft,at,bg,s0,ft_inter)
!
!   ft_inter in cartesian coordinates
!
  CALL cryst_to_cart(1, ft_inter, at, 1)
!
!   Apply rot^-1
!
  ftnew=0.0_DP
  DO ipol=1,3
     DO jpol=1,3
        ftnew(ipol)=ftnew(ipol)+sr(jpol,ipol)*ft_inter(jpol)
     ENDDO
  ENDDO
!
!  and return to crystal coordinates
!
  CALL cryst_to_cart(1, ftnew, bg, -1)
!
!  Then rotate sr => srnew = rot^-1 sr rot
!  
  srnew=MATMUL(TRANSPOSE(rot),MATMUL(sr,rot))

  RETURN
  END SUBROUTINE roto_shift_sg


!-----------------------------------------------------------------------
  SUBROUTINE find_sg_tags(sr,ft,ibrav,at,bg,ftperp,ftpar,sym_in,s0,sg_sym_out)
!-----------------------------------------------------------------------
!
!   This routine receives as input a space group operation in the form
!   of a 3x3 matrix sr (in cartesian coordinates), 
!   a fractional translation ft, the direct and reciprocal lattice vectors, 
!   the number of the rotation, the projection of the fractional 
!   translation along the rotation axis or along the mirror plane 
!   and the projection along the perpendicular direction and 
!   gives as output the number of the space group operation and the 
!   origin shift needed to remove the part of the fractional translation 
!   perpendicular to the rotation axis or to the mirror plane. For 
!   improper rotation the new origin removes all the fractional translation.
!   ft,ftperp,ftpar and s0 are in crystal coordinates of at.
!
!   The codes of the generalized space group operations are the following
!
!   E         1 
!   2z        2   2z1      65
!   2y        3   2y1      66
!   2x        4   2x1      67
!   2xy       5   2xy1     68
!   2x-y      6   2x-y1    69
!   4-z       7   4-z1     70    4-z2     71     4-z3   72
!   4z        8   4z1      73    4z2      74     4z3    75
!   2xz       9   2xz1     76
!   2x-z     10   2x-z1    77     
!   4y       11   4y1      78    4y2      79     4y3    80
!   4-y      12   4-y1     81    4-y2     82     4-y3   83
!   2yz      13   2yz1     84
!   2y-z     14   2y-z1    85
!   4-x      15   4-x1     86    4-x2     87     4-x3   88
!   4x       16   4x1      89    4x2      90     4x3    91
!   3-x-y-z  17   3-x-y-z1 92    3-x-y-z2 93
!   3-xyz    18   3-xyz1   94    3-xyz2   95
!   3xy-z    19   3xy-z1   96    3xy-z2   97
!   3x-yz    20   3x-yz1   98    3x-yz2   99
!   3xyz     21   3xyz1   100    3xyz2   101
!   3-xy-z   22   3-xy-z1 102    3-xy-z2 103
!   3x-y-z   23   3x-y-z1 104    3x-y-z2 105
!   3-x-yz   24   3-x-yz1 106    3-x-yz2 107
!   6z       25   6z1     108    6z2     109    6z3   110   6z4   111  6z5  112
!   6-z      26   6-z1    113    6-z2    114    6-z3  115   6-z4  116  6-z5 117
!   3z       27   3z1     118    3z2     119
!   3-z      28   3-z1    120    3-z2    121
!   21-10    29   21-101  122
!   2210     30   22101   123
!   2010     31   20101   124
!   2110     32   21101   125
!   i        33  
!   i2z      34   i2za    126    i2zb    127    i2zn   128    i2zd    129 
!   i2y      35   i2ya    130    i2yc    131    i2yn   132    i2yd    133
!   i2x      36   i2xb    134    i2xc    135    i2xn   136    i2xd    137
!   i2xy     37   i2xya   138    i2xyc   139    i2xyn  140    i2xyd   141 
!   i2x-y    38   i2x-ya  142    i2x-yc  143    i2x-yn 144    i2x-yd  145
!   i4-z     39   
!   i4z      40   
!   i2xz     41   i2xza   146    i2xzb   147    i2xzn  148    i2xzd   149
!   i2x-z    42   i2x-za  150    i2x-zb  151    i2x-zn 152    i2x-zd  153
!   i4y      43
!   i4-y     44
!   i2yz     45   i2yza   154    i2yzb   155     i2yzn  156    i2yzd   157
!   i2y-z    46   i2y-za  158    i2y-zb  159     i2y-zn 160    i2y-zd  161
!   i4x      47
!   i4-x     48
!   i3-x-y-z 49  
!   i3-xyz   50 
!   i3xy-z   51  
!   i3x-yz   52 
!   i3xyz    53  
!   i3-xy-z  54 
!   i3x-y-z  55   
!   i3-x-yz  56   
!   i6z      57
!   i6-z     58
!   i3z      59
!   i3-z     60
!   i21-10   61   i21-10a 162    i21-10c  163     i21-10n  164   
!   i2210    62   i2210a  165    i2210c   166     i2210n   167   
!   i2010    63   i2010a  168    i2010c   169     i2010n   170   
!   i2110    64   i2110a  171    i2110c   172     i2110n   173   
!
!
  USE kinds, ONLY : DP
  USE linear_solvers, ONLY : linsolvx
  USE lattices, ONLY : compute_conventional
  IMPLICIT NONE
!
!   First the operations for which ft is zero
!
  INTEGER, INTENT(IN) :: sym_in, ibrav
  REAL(DP), INTENT(IN) :: sr(3,3), ft(3), at(3,3), bg(3,3), ftperp(3), &
                          ftpar(3)
  REAL(DP), INTENT(OUT) :: s0(3)
  INTEGER, INTENT(OUT) :: sg_sym_out

  REAL(DP)  :: ftmod, amat(3,3), b(3), ftpar2(3), ftpar3(3), ftpar4(3), &
               ftpar6(3), ftmod2, ftmod3, ftmod4, ftmod6, atp(3,3), bgp(3,3),&
               ftpar2m(3), ftmod2m, u(3), ps
  INTEGER :: ts, tipo_sym, iftpar2(3), iftpar3(3), iftpar4(3), iftpar6(3), &
             iftpar2m(3), imirror, ifact, ipol, jpol,idir
!
!  symmorphic operations. The space group operation coincides with the
!  point group operation
!
!  WRITE(6,*) 'find_space_group_tags'
!  DO ipol=1,3
!     WRITE(6,*) (at(jpol,ipol), jpol=1,3)
!  ENDDO
!  WRITE(6,*) 'find_bg'
!  DO ipol=1,3
!     WRITE(6,*) (bg(jpol,ipol), jpol=1,3)
!  ENDDO

  s0=0.0_DP
  ftmod = ft(1)**2 + ft(2)**2 + ft(3)**2
  IF (ftmod<1.D-6) THEN
     sg_sym_out=sym_in
     RETURN
  ENDIF
!
! symmorphic operation with displaced origin. 
! The inversion first
!
  IF (sym_in==33) THEN
     s0(:) = ft(:) * 0.5_DP
     sg_sym_out=sym_in
  END IF 

  ts=tipo_sym(sr) 
  ftmod=ftpar(1)**2 + ftpar(2)**2 + ftpar(3)**2

  IF (ftmod<1.D-6) THEN
!
! symmorphic operation with displaced origin. 
! mirror symmetry or twofold rotation axis
!
     IF (ts==4 .OR. ts==5) THEN
        s0(:)=ftperp(:) * 0.5_DP
     ENDIF
!
! improper rotation
!
     IF (ts==6) THEN
        amat(:,:)=-sr(:,:)
        DO ipol=1,3
           amat(ipol,ipol)=1.0_DP - sr(ipol,ipol)
        ENDDO
        b(:)=ftperp(:)
        CALL cryst_to_cart(1, b, at, 1)        
        CALL linsolvx(amat,3,b,s0)
        CALL cryst_to_cart(1, s0, bg, -1)        
     ENDIF
!
!  proper rotation about x, y, or z or three-fold rotation about the
!  cube diagonals
!
     IF (sym_in==7.OR.sym_in==8.OR.sym_in==25.OR.sym_in==26 &
         .OR. sym_in==27.OR.sym_in==28) THEN
        CALL find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,3)
     ELSEIF (sym_in==11.OR.sym_in==12) THEN
        CALL find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,2)
     ELSEIF (sym_in==15.OR.sym_in==16) THEN
        CALL find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,1)
     ELSEIF (sym_in>=17.AND.sym_in<=24) THEN
        CALL find_origin_3_rot(sr,ftperp,sym_in,at,bg,s0)
     ENDIF
     
     sg_sym_out=sym_in
     RETURN
  ENDIF
!
!   We first find the origin shift
!
  u=0.0_DP
  IF (ts==4 .OR. ts==5) THEN
     s0(:)=ftperp(:) * 0.5_DP
  ELSEIF (sym_in==7.OR.sym_in==8.OR.sym_in==25.OR.sym_in==26 &
         .OR. sym_in==27.OR.sym_in==28) THEN
     CALL find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,3)
     u(3)=1.0_DP
  ELSEIF (sym_in==11.OR.sym_in==12) THEN
     CALL find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,2)
     u(2)=1.0_DP
  ELSEIF (sym_in==15.OR.sym_in==16) THEN
     CALL find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,1)
     u(1)=1.0_DP
  ELSEIF (sym_in>=17.AND.sym_in<=24) THEN
     CALL find_origin_3_rot(sr,ftperp,sym_in,at,bg,s0)
     CALL versor(sr,u)
  ENDIF
!
!  and then decide which operation it is. First find if the
!  fractional translation is parallel or antiparallel to the axis
!
  b(:)=ftpar(:)
  CALL cryst_to_cart(1, b, at, 1)        
  ps=b(1)*u(1)+b(2)*u(2)+b(3)*u(3)
  idir=0
  IF (ps > 1.D-6) THEN
     idir=1
  ELSEIF (ps < -1.D-6) THEN
     idir=-1
  ENDIF   
  ftpar2(:)=ftpar(:)*2.0_DP - NINT(ftpar(:)*2.0_DP)
  ftpar3(:)=ftpar(:)*3.0_DP - NINT(ftpar(:)*3.0_DP)
  ftpar4(:)=ftpar(:)*4.0_DP - NINT(ftpar(:)*4.0_DP)
  ftpar6(:)=ftpar(:)*6.0_DP - NINT(ftpar(:)*6.0_DP)
  ftmod2 = ftpar2(1)**2 + ftpar2(2)**2 + ftpar2(3)**2
  ftmod3 = ftpar3(1)**2 + ftpar3(2)**2 + ftpar3(3)**2
  ftmod4 = ftpar4(1)**2 + ftpar4(2) **2 + ftpar4(3)**2
  ftmod6 = ftpar6(1)**2 + ftpar6(2) **2 + ftpar6(3)**2
  iftpar2(:)= NINT(ftpar(:)*2.0_DP)
  iftpar3(:)= NINT(ftpar(:)*3.0_DP)
  iftpar4(:)= NINT(ftpar(:)*4.0_DP)
  iftpar6(:)= NINT(ftpar(:)*6.0_DP)
  CALL compute_conventional(at, atp, ibrav)
  CALL recips(atp(1,1), atp(1,2), atp(1,3), bgp(1,1), bgp(1,2), bgp(1,3))
  CALL cryst_to_cart(1,b,bgp,-1)
!
! for the mirror we check the glide with respect to the vectors of the
! conventional unit cell (for centered lattices)
!
  ftpar2m(:)=b(:)*2.0_DP - NINT(b(:)*2.0_DP)
  ftmod2m = ftpar2m(1)**2 + ftpar2m(2) **2 + ftpar2m(3)**2
  iftpar2m(:)= NINT(b(:)*2.0_DP)
  SELECT CASE (sym_in)
     CASE(2)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=65
        ELSE
           GOTO 100
        ENDIF
     CASE(3)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=66
        ELSE
           GOTO 100
        ENDIF
     CASE(4)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=67
        ELSE
           GOTO 100
        ENDIF
     CASE(5)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=68
        ELSE
           GOTO 100
        ENDIF
     CASE(6)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=69
        ELSE
           GOTO 100
        ENDIF
     CASE(7)
        IF (ftmod4 < 1.D-6) THEN
           CALL find_ifact(4,iftpar4,ifact)       
           IF (ifact==1) THEN
              IF (idir==1) THEN
                 sg_sym_out=72
              ELSE
                 sg_sym_out=70
              ENDIF
           ELSEIF (ifact==2) THEN
              sg_sym_out=71
           ELSEIF (ifact==3) THEN
              IF (idir==1) THEN
                 sg_sym_out=70
              ELSE
                 sg_sym_out=72
              ENDIF
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(8)
        IF (ftmod4 < 1.D-6) THEN
           CALL find_ifact(4,iftpar4,ifact)       
           IF (ifact==1) THEN
              IF (idir==1) THEN
                 sg_sym_out=73
              ELSE
                 sg_sym_out=75
              ENDIF
           ELSEIF (ifact==2) THEN
              sg_sym_out=74
           ELSEIF (ifact==3) THEN
              IF (idir==1) THEN
                 sg_sym_out=75
              ELSE
                 sg_sym_out=73
              ENDIF
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(9)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=76
        ELSE
           GOTO 100
        ENDIF
     CASE(10)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=77
        ELSE
           GOTO 100
        ENDIF
     CASE(11)
        IF (ftmod4 < 1.D-6) THEN
           CALL find_ifact(4,iftpar4,ifact)       
           IF (ifact==1) THEN
              IF (idir==1) THEN
                 sg_sym_out=78
              ELSE
                 sg_sym_out=80
              ENDIF
           ELSEIF (ifact==2) THEN
              sg_sym_out=79
           ELSEIF (ifact==3) THEN
              IF (idir==1) THEN
                 sg_sym_out=80
              ELSE
                 sg_sym_out=78
              ENDIF
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(12)
        IF (ftmod4 < 1.D-6) THEN
           CALL find_ifact(4,iftpar4,ifact)       
           IF (ifact==1) THEN
              IF (idir==1) THEN
                 sg_sym_out=83
              ELSE
                 sg_sym_out=81
              ENDIF
           ELSEIF (ifact==2) THEN
              sg_sym_out=82
           ELSEIF (ifact==3) THEN
              IF (idir==1) THEN
                 sg_sym_out=81
              ELSE
                 sg_sym_out=83
              ENDIF
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(13)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=84
        ELSE
           GOTO 100
        ENDIF
     CASE(14)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=85
        ELSE
           GOTO 100
        ENDIF
     CASE(15)
        IF (ftmod4 < 1.D-6) THEN
           CALL find_ifact(4,iftpar4,ifact)       
           IF (ifact==1) THEN
              IF (idir==1) THEN
                 sg_sym_out=88
              ELSE
                 sg_sym_out=86
              ENDIF
           ELSEIF (ifact==2) THEN
              sg_sym_out=87
           ELSEIF (ifact==3) THEN
              IF (idir==1) THEN
                 sg_sym_out=86
              ELSE
                 sg_sym_out=88
              ENDIF
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(16)
        IF (ftmod4 < 1.D-6) THEN
           CALL find_ifact(4,iftpar4,ifact)       
           IF (ifact==1) THEN
              IF (idir==1) THEN
                 sg_sym_out=89
              ELSE
                 sg_sym_out=91
              ENDIF
           ELSEIF (ifact==2) THEN
              sg_sym_out=90
           ELSEIF (ifact==3) THEN
              IF (idir==1) THEN
                 sg_sym_out=91
              ELSE
                 sg_sym_out=89
              ENDIF
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(17)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=92
              ELSEIF (ifact==2) THEN
                 sg_sym_out=93
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=93
              ELSEIF (ifact==2) THEN
                 sg_sym_out=92
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(18)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=94
              ELSEIF (ifact==2) THEN
                 sg_sym_out=95
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=95
              ELSEIF (ifact==2) THEN
                 sg_sym_out=94
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(19)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=96
              ELSEIF (ifact==2) THEN
                 sg_sym_out=97
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=97
              ELSEIF (ifact==2) THEN
                 sg_sym_out=96
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(20)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=98
              ELSEIF (ifact==2) THEN
                 sg_sym_out=99
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=99
              ELSEIF (ifact==2) THEN
                 sg_sym_out=98
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(21)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=100
              ELSEIF (ifact==2) THEN
                 sg_sym_out=101
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=101
              ELSEIF (ifact==2) THEN
                 sg_sym_out=100
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(22)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=103
              ELSEIF (ifact==2) THEN
                 sg_sym_out=104
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=104
              ELSEIF (ifact==2) THEN
                 sg_sym_out=103
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(23)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=105
              ELSEIF (ifact==2) THEN
                 sg_sym_out=106
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=106
              ELSEIF (ifact==2) THEN
                 sg_sym_out=105
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(24)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=106
              ELSEIF (ifact==2) THEN
                 sg_sym_out=107
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=107
              ELSEIF (ifact==2) THEN
                 sg_sym_out=106
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(25)
        IF (ftmod6 < 1.D-6) THEN
           CALL find_ifact(6,iftpar6,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=108
              ELSEIF (ifact==2) THEN
                 sg_sym_out=109
              ELSEIF (ifact==3) THEN
                 sg_sym_out=110
              ELSEIF (ifact==4) THEN
                 sg_sym_out=111
              ELSEIF (ifact==5) THEN
                 sg_sym_out=112
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=112
              ELSEIF (ifact==2) THEN
                 sg_sym_out=111
              ELSEIF (ifact==3) THEN
                 sg_sym_out=110
              ELSEIF (ifact==4) THEN
                 sg_sym_out=109
              ELSEIF (ifact==5) THEN
                 sg_sym_out=108
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(26)
        IF (ftmod6 < 1.D-6) THEN
           CALL find_ifact(6,iftpar6,ifact)       
           IF (idir==-1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=113
              ELSEIF (ifact==2) THEN
                 sg_sym_out=114
              ELSEIF (ifact==3) THEN
                 sg_sym_out=115
              ELSEIF (ifact==4) THEN
                 sg_sym_out=116
              ELSEIF (ifact==5) THEN
                 sg_sym_out=117
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=117
              ELSEIF (ifact==2) THEN
                 sg_sym_out=116
              ELSEIF (ifact==3) THEN
                 sg_sym_out=115
              ELSEIF (ifact==4) THEN
                 sg_sym_out=114
              ELSEIF (ifact==5) THEN
                 sg_sym_out=113
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(27)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=118
              ELSEIF (ifact==2) THEN
                 sg_sym_out=119
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=119
              ELSEIF (ifact==2) THEN
                 sg_sym_out=118
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(28)
        IF (ftmod3 < 1.D-6) THEN
           CALL find_ifact(3,iftpar3,ifact)       
           IF (idir==1) THEN
              IF (ifact==1) THEN
                 sg_sym_out=121
              ELSEIF (ifact==2) THEN
                 sg_sym_out=120
              ENDIF   
           ELSE
              IF (ifact==1) THEN
                 sg_sym_out=120
              ELSEIF (ifact==2) THEN
                 sg_sym_out=121
              ENDIF   
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(29)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=122
        ELSE
           GOTO 100
        ENDIF
     CASE(30)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=123
        ELSE
           GOTO 100
        ENDIF
     CASE(31)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=124
        ELSE
           GOTO 100
        ENDIF
     CASE(32)
        IF (ftmod2 < 1.D-6) THEN
           sg_sym_out=125
        ELSE
           GOTO 100
        ENDIF
     CASE(34)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==1) THEN
              sg_sym_out=126
           ELSEIF (imirror==2) THEN
              sg_sym_out=127
           ELSEIF (imirror==4) THEN
              sg_sym_out=128
           ELSE
              GOTO 100
           ENDIF
        ELSE
!
!   check for d
!
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=129
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(35)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==1) THEN
              sg_sym_out=130
           ELSEIF (imirror==3) THEN
              sg_sym_out=131
           ELSEIF (imirror==4) THEN
              sg_sym_out=132
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=133
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(36)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==2) THEN
              sg_sym_out=134
           ELSEIF (imirror==3) THEN
              sg_sym_out=135
           ELSEIF (imirror==4) THEN
              sg_sym_out=136
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=137
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(37)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==3) THEN
              sg_sym_out=139
           ELSEIF (imirror==4) THEN
              sg_sym_out=140
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=141
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(38)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==3) THEN
              sg_sym_out=143
           ELSEIF (imirror==4) THEN
              sg_sym_out=144
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=145
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(41)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==2) THEN
              sg_sym_out=147
           ELSEIF (imirror==4) THEN
              sg_sym_out=148
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=149
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(42)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==2) THEN
              sg_sym_out=151
           ELSEIF (imirror==4) THEN
              sg_sym_out=152
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=153
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(45)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==1) THEN
              sg_sym_out=154
           ELSEIF (imirror==4) THEN
              sg_sym_out=156
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=157
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(46)
        IF (ftmod2m < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==1) THEN
              sg_sym_out=158
           ELSEIF (imirror==4) THEN
              sg_sym_out=160
           ELSE
              GOTO 100
           ENDIF
        ELSE
           IF (ftmod2 < 1.D-6) THEN
              sg_sym_out=157
           ELSE
              GOTO 100
           ENDIF
        ENDIF
     CASE(61)
        IF (ftmod2 < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==1.OR.imirror==2) THEN
              sg_sym_out=162
           ELSEIF (imirror==3) THEN
              sg_sym_out=163
           ELSEIF (imirror==4) THEN
              sg_sym_out=164
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(62)
        IF (ftmod2 < 1.D-6) THEN
           CALL find_mirror_type(iftpar2m,imirror)
           IF (imirror==1.OR.imirror==2) THEN
              sg_sym_out=165
           ELSEIF (imirror==3) THEN
              sg_sym_out=166
           ELSEIF (imirror==4) THEN
              sg_sym_out=167
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(63)
        IF (ftmod2 < 1.D-6) THEN
           IF (imirror==1.OR.imirror==2) THEN
              sg_sym_out=168
           ELSEIF (imirror==3) THEN
              sg_sym_out=169
           ELSEIF (imirror==4) THEN
              sg_sym_out=170
           ENDIF
        ELSE
           GOTO 100
        ENDIF
     CASE(64)
        IF (ftmod2 < 1.D-6) THEN
           IF (imirror==1.OR.imirror==2) THEN
              sg_sym_out=171
           ELSEIF (imirror==3) THEN
              sg_sym_out=172
           ELSEIF (imirror==4) THEN
              sg_sym_out=173
           ENDIF
        ELSE
           GOTO 100
        ENDIF
  CASE DEFAULT
       GOTO 100
  END SELECT
  RETURN
100 CALL errore('find_sg_tags','Some problems. The routine should not &
                                                            &arrive here',1)
  RETURN
  END SUBROUTINE find_sg_tags

!-----------------------------------------------------------------------
SUBROUTINE find_mirror_type(ifrac,imirror)
!-----------------------------------------------------------------------
!
!  imirror codes
!  1  a
!  2  b
!  3  c
!  4  n
!  All codes are appropriate for mirror d if ifrac is calculated with
!  the centered lattice vectors. Note however that a,b,c,n must be 
!  checked before d.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ifrac(3)
INTEGER, INTENT(OUT) :: imirror

IF (ifrac(1)/=0.AND.ifrac(2)==0.AND.ifrac(3)==0) imirror=1 
IF (ifrac(1)==0.AND.ifrac(2)/=0.AND.ifrac(3)==0) imirror=2 
IF (ifrac(1)==0.AND.ifrac(2)==0.AND.ifrac(3)/=0) imirror=3 
IF (ifrac(1)/=0.AND.ifrac(2)/=0.AND.ifrac(3)==0) imirror=4    
IF (ifrac(1)==0.AND.ifrac(2)/=0.AND.ifrac(3)/=0) imirror=4 
IF (ifrac(1)/=0.AND.ifrac(2)==0.AND.ifrac(3)/=0) imirror=4 
IF (ifrac(1)/=0.AND.ifrac(2)/=0.AND.ifrac(3)/=0) imirror=4

RETURN
END SUBROUTINE find_mirror_type

!-----------------------------------------------------------------------
SUBROUTINE find_ifact(axis_order,ifrac,ifact)
!-----------------------------------------------------------------------
!
!   This routine receives a vector of three integers ifrac.
!   It removes the zeros and for the remaining integers it checks
!   if they have a common factor k. k, modulus the input value
!   axis_order, is set into ifact. 
!   The use is the following. A fractional translation parallel to a
!   rotation axis must be parallel to a Bravais lattice vector
!   n_1 a_1 + n_2 a_2 + n_3 a_3. If n_1, n_2 and n_3 have no common factor
!   this is the shortest Bravais lattice vector parallel to the axis.
!   Multiplying the fractional translation written in crystal coordinates 
!   by axis_order (the order of the axis) we will obtain:
!   k n_1 a_1 + k n_2 a_2 + k n_3 a_3
!   this routine finds k (mod axis_order) receiving in input k n_1, k n_2 and
!   k n_3.
!   
!
IMPLICIT NONE
INTEGER, INTENT(IN) ::  axis_order, ifrac(3)
INTEGER, INTENT(OUT) :: ifact

INTEGER :: ifa(3), nterms, ipol, iorder

nterms=0

DO ipol=1,3
   IF (ifrac(ipol) /= 0) THEN
      nterms=nterms+1
      ifa(nterms)=ifrac(ipol)
   ENDIF
ENDDO

IF (nterms==0) CALL errore('find_ifact','called in the wrong case',1)

DO iorder=axis_order-1, 1, -1
   IF (nterms==1) THEN
      IF (MOD(ifa(1),iorder)==0) THEN
         ifact=iorder
         RETURN
      ENDIF
   ELSEIF (nterms==2) THEN
      IF ( (MOD(ifa(1),iorder)==0).AND.(MOD(ifa(2),iorder)==0)) THEN
         ifact=iorder
         RETURN
      END IF
   ELSEIF (nterms==3) THEN
      IF ((MOD(ifa(1),iorder)==0).AND. &
          (MOD(ifa(2),iorder)==0).AND. &
          (MOD(ifa(3),iorder)==0)) THEN
         ifact=iorder
         RETURN
      END IF
   ENDIF
END DO

RETURN
END SUBROUTINE find_ifact

!-----------------------------------------------------------------------
SUBROUTINE find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,ipol)
!-----------------------------------------------------------------------
!
! This routine finds the origin shift needed to bring the fractional
! translation associated to a rotation whose axis is parallel to x, y, or
! z to be parallel to the rotation axis
!
USE kinds, ONLY : DP
USE linear_solvers, ONLY : linsolvx

IMPLICIT NONE
INTEGER, INTENT(IN) :: sym_in, ipol
REAL(DP), INTENT(IN) :: sr(3,3), ftperp(3), at(3,3), bg(3,3)
REAL(DP), INTENT(OUT) :: s0(3)


REAL(DP) ::  &  
         b(3),      & ! the fractional translation in cartesian coordinates
         amat(3,3), & ! the matrix E-R
         cmat(2,2), & ! the matrix E-R reduced to 2x2
         c(2),      & ! auxiliary vector to contain a part of b
         s(2)         ! auxiliary vector solution of the linear system

INTEGER :: jpol
!
!  bring the fractional translation in cartesian coordinates
!
b(:)=ftperp(:)
CALL cryst_to_cart(1,b,at,1)
!
!  write the matrix E-R
!
amat(:,:)=-sr(:,:)
DO jpol=1,3
   amat(jpol,jpol)=1.0_DP-sr(jpol,jpol)
ENDDO
!
!  construct the linear system to invert. According to ipol choose a
!  part of the rotation matrix
!
IF (ipol==1) THEN
   cmat(1,1) = amat(2,2)
   cmat(1,2) = amat(2,3)
   cmat(2,1) = amat(3,2)
   cmat(2,2) = amat(3,3)
   c(1)=b(2)
   c(2)=b(3)
ELSEIF (ipol==2) THEN
   cmat(1,1) = amat(1,1)
   cmat(1,2) = amat(1,3)
   cmat(2,1) = amat(3,1)
   cmat(2,2) = amat(3,3)
   c(1)=b(1)
   c(2)=b(3)
ELSEIF (ipol==3) THEN
   cmat(1,1) = amat(1,1)
   cmat(1,2) = amat(1,2)
   cmat(2,1) = amat(2,1)
   cmat(2,2) = amat(2,2)
   c(1)=b(1)
   c(2)=b(2)
ENDIF
!
!  solve the linear system that gives s(1) and s(2)
!
CALL linsolvx(cmat,2,c,s)
!
!  uncomment to check the linear system
!
!WRITE(6,*) 'linear system'
!WRITE(6,*) cmat(1,1)*s(1)+ cmat(1,2)*s(2) - c(1)
!WRITE(6,*) cmat(2,1)*s(1)+ cmat(2,2)*s(2) - c(2)
!
!  compute s0
!
IF (ipol==1) THEN
   s0(1) = 0.0_DP
   s0(2) = s(1)
   s0(3) = s(2)
ELSEIF (ipol==2) THEN
   s0(1) = s(1)
   s0(2) = 0.0_DP
   s0(3) = s(2)
ELSEIF (ipol==3) THEN
   s0(1) = s(1)
   s0(2) = s(2)
   s0(3) = 0.0_DP
ENDIF
!
!  bring s0 in crystal coordinates
!
CALL cryst_to_cart(1,s0,bg,-1)

RETURN
END SUBROUTINE find_origin_xyz_rot

!-----------------------------------------------------------------------
SUBROUTINE find_origin_3_rot(sr,ftperp,sym_in,at,bg,s0)
!-----------------------------------------------------------------------
!
! This routine finds the origin shift needed to bring the fractional
! translation associated to a rotation of +-120 about a diagonal of the cube
! to be parallel to the rotation axis
!
USE kinds, ONLY : DP
USE constants, ONLY : pi, tpi
USE linear_solvers, ONLY : linsolvx

IMPLICIT NONE
INTEGER, INTENT(IN) :: sym_in
REAL(DP), INTENT(IN) :: sr(3,3), ftperp(3), at(3,3), bg(3,3)
REAL(DP), INTENT(OUT) :: s0(3)

REAL(DP), PARAMETER :: sqrt3=SQRT(3.0_DP)
REAL(DP) :: u(3),      &  ! rotation axis
            b(3),      &  ! the fractional translation
            amat(3,3), &  ! the matrix E-R
            cmat(2,2), &  ! the matrix E-R reduced to a 2x2
            alpha, beta   ! auxiliary quantities

INTEGER :: ipol
!
!  set the rotation axis
!
CALL versor(sr,u)
!
!  bring the fractional translation in cartesian coordinates
!
b(:)=ftperp(:)
CALL cryst_to_cart(1,b,at,1)
!
!  write the matrix E-R
!
amat(:,:)=-sr(:,:)
DO ipol=1,3
   amat(ipol,ipol)=1.0_DP-sr(ipol,ipol)
ENDDO
!
!  find alpha and beta  s0(3)=alpha s0(1) + beta s0(2)
!
alpha=-u(1)/u(3)
beta=-u(2)/u(3)
!
!  construct the linear system to invert
!
cmat(1,1) = amat(1,1)+alpha*amat(1,3)
cmat(1,2) = amat(1,2)+beta*amat(1,3)
cmat(2,1) = amat(2,1)+alpha*amat(2,3)
cmat(2,2) = amat(2,2)+beta*amat(2,3)
!
!  solve the linear system that gives s0(1) and s0(2)
!
CALL linsolvx(cmat,2,b,s0)
!
!  uncomment to check the linear system
!
!WRITE(6,*) 'linear system'
!WRITE(6,*) cmat(1,1)*s0(1)+ cmat(1,2)*s0(2) - b(1)
!WRITE(6,*) cmat(2,1)*s0(1)+ cmat(2,2)*s0(2) - b(2)
!
!  compute s0(3)
!
s0(3) = alpha * s0(1) + beta * s0(2)
!
!  bring s0 in crystal coordinates
!
CALL cryst_to_cart(1,s0,bg,-1)

RETURN
END SUBROUTINE find_origin_3_rot

!-----------------------------------------------------------------------
SUBROUTINE set_sg_ibrav(sgc, ibrav)
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: sgc
INTEGER, INTENT(OUT) :: ibrav

INTEGER :: sg_ibrav(230)

DATA sg_ibrav  /  14,  14,  12,  12,  13,  12,  12,  13,  13,  12,  &
                  12,  13,  12,  12,  13,   8,   8,   8,   8,   9,  &
                   9,  10,  11,  11,   8,   8,   8,   8,   8,   8,  &
                   8,   8,   8,   8,   9,   9,   9,  91,  91,  91,  &
                  91,  10,  10,  11,  11,  11,   8,   8,   8,   8,  &
                   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,  &
                   8,   8,   9,   9,   9,   9,   9,   9,  10,  10,  &
                  11,  11,  11,  11,   6,   6,   6,   6,   7,   7,  &
                   6,   7,   6,   6,   6,   6,   7,   7,   6,   6,  &
                   6,   6,   6,   6,   6,   6,   7,   7,   6,   6,  &
                   6,   6,   6,   6,   6,   6,   7,   7,   7,   7,  &
                   6,   6,   6,   6,   6,   6,   6,   6,   7,   7,  &
                   7,   7,   6,   6,   6,   6,   6,   6,   6,   6,  &
                   6,   6,   6,   6,   6,   6,   6,   6,   7,   7,  &
                   7,   7,   4,   4,   4,   5,   4,   5,   4,   4,  &
                   4,   4,   4,   4,   5,   4,   4,   4,   4,   5,  &
                   5,   4,   4,   4,   4,   5,   5,   4,   4,   4,  &
                   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,  &
                   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,  &
                   4,   4,   4,   4,   1,   2,   3,   1,   3,   1,  &
                   1,   2,   2,   3,   1,   3,   1,   1,   2,   2,  &
                   3,   1,   1,   3,   1,   2,   3,   1,   2,   3,  &
                   1,   1,   1,   1,   2,   2,   2,   2,   3,   3  /            

ibrav=sg_ibrav(sgc)

RETURN
END SUBROUTINE set_sg_ibrav

!-----------------------------------------------------------------------
SUBROUTINE set_point_group_code(sgc, group_code, group_code_ext)
!-----------------------------------------------------------------------
!
!   This routine receives the group code of a space group (1-230)
!   and gives the code of the point group, both in condensed and
!   extended form
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: sgc
INTEGER, INTENT(OUT) :: group_code, group_code_ext

INTEGER ::  group_code_list(230)

DATA group_code_list /   1,   2,   4,   4,   4,   3,   3,   3,   3,  16, &
                        16,  16,  16,  16,  16,   8,   8,   8,   8,   8, &
                         8,   8,   8,   8,  12,  12,  12,  12,  12,  12, &
                        12,  12,  12,  12,  12,  12,  12,  12,  12,  12, &
                        12,  12,  12,  12,  12,  12,  20,  20,  20,  20, &
                        20,  20,  20,  20,  20,  20,  20,  20,  20,  20, &
                        20,  20,  20,  20,  20,  20,  20,  20,  20,  20, &
                        20,  20,  20,  20,   6,   6,   6,   6,   6,   6, &
                        26,  26,  18,  18,  18,  18,  18,  18,  10,  10, &
                        10,  10,  10,  10,  10,  10,  10,  10,  14,  14, &
                        14,  14,  14,  14,  14,  14,  14,  14,  14,  14, &
                        24,  24,  24,  24,  24,  24,  24,  24,  24,  24, &
                        24,  24,  22,  22,  22,  22,  22,  22,  22,  22, &
                        22,  22,  22,  22,  22,  22,  22,  22,  22,  22, &
                        22,  22,   5,   5,   5,   5,  27,  27,   9,   9, & 
                         9,   9,   9,   9,   9,  13,  13,  13,  13,  13, & 
                        13,  25,  25,  25,  25,  25,  25,   7,   7,   7, &
                         7,   7,   7,  17,  19,  19,  11,  11,  11,  11, &
                        11,  11,  15,  15,  15,  15,  21,  21,  21,  21, &
                        23,  23,  23,  23,  28,  28,  28,  28,  28,  29, &
                        29,  29,  29,  29,  29,  29,  31,  31,  31,  31, &
                        31,  31,  31,  31,  30,  30,  30,  30,  30,  30, &
                        32,  32,  32,  32,  32,  32,  32,  32,  32,  32  /

INTEGER ::  group_code_list_ext(230)

DATA group_code_list_ext /   1,  28,   2,   2,   2,  15,  15,  15,  15,  82, &
                            82,  82,  82,  82,  82,  38,  38,  38,  38,  38, &
                            38,  38,  38,  38,  58,  58,  58,  58,  58,  58, &
                            58,  58,  58,  58,  58,  58,  58,  58,  58,  58, &
                            58,  58,  58,  58,  58,  58, 100, 100, 100, 100, &
                           100, 100, 100, 100, 100, 100, 100, 100, 100, 100, &
                           100, 100, 100, 100, 100, 100, 100, 100, 100, 100, &
                           100, 100, 100, 100,  34,  34,  34,  34,  34,  34, &
                           124, 124,  96,  96,  96,  96,  96,  96,  50,  50, &
                            50,  50,  50,  50,  50,  50,  50,  50,  78,  78, &
                            78,  78,  78,  78,  78,  78,  78,  78,  78,  78, &
                           112, 112, 112, 112, 113, 113, 113, 113, 113, 113, &
                           112, 112, 108, 108, 108, 108, 108, 108, 108, 108, &
                           108, 108, 108, 108, 108, 108, 108, 108, 108, 108, &
                           108, 108,  33,  33,  33,  33, 127, 127,  45,  44, &
                            45,  44,  45,  44,  44,  72,  73,  72,  73,  72, &
                            73, 119, 119, 118, 118, 118, 118,  37,  37,  37, &
                            37,  37,  37,  95,  99,  99,  53,  53,  53,  53, &
                            53,  53,  81,  81,  81,  81, 107, 107, 106, 106, &
                           111, 111, 111, 111, 132, 132, 132, 132, 132, 133, &
                           133, 133, 133, 133, 133, 133, 135, 135, 135, 135, &
                           135, 135, 135, 135, 134, 134, 134, 134, 134, 134, &
                           136, 136, 136, 136, 136, 136, 136, 136, 136, 136  / 


group_code=group_code_list(sgc)
group_code_ext=group_code_list_ext(sgc)

RETURN
END SUBROUTINE set_point_group_code

!-----------------------------------------------------------------------
SUBROUTINE set_ft_generators(sgc,ft)
!-----------------------------------------------------------------------
!
! This routine receives in input the space group code and sets the
! fractionary translations on the generators. The elements of ft are
! supposed to be in the standard order and only the fractional translations
! of the generators are set by this routine. All the others can be 
! generate by the routine. The generators are those of Table 13.3 of 
! S.K. Kim, Group theoretical methods and applications to molecules and
! crystals.
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: sgc
REAL(DP), INTENT(OUT) :: ft(3,48)

ft=0.0_DP
SELECT CASE (sgc)
   CASE(1,2,3,5,6,8,10,12,16,21,22,23,25,35,38,42,44,47,65,69,71,75, &
        79,81,82,83,87,89,97,99,107,111,115,119,121,123,139,143,146, &
        147,148,149,150,155,156,157,160,162,164,166,168,174,175,177,183, &
        187,189,191,195,196,197,200,202,204,207,209,211,215,216,217,221, &
        225,229)
!
!  These are the 73 symmorphic groups. All fractional translations vanish.
!
   CASE(4)
     ft(3,2)=0.5_DP
   CASE(7,9)
     ft(2,2)=0.5_DP
   CASE(11)
     ft(3,2)=0.5_DP
   CASE(13,15)
     ft(2,2)=0.5_DP
   CASE(14)
     ft(2,2)=0.5_DP
     ft(3,2)=0.5_DP
   CASE(17,20)
     ft(3,2)=0.5_DP
   CASE(18)
     ft(1,2)=0.5_DP
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
   CASE(19,24)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
   CASE(26,36,46)
     ft(3,2)=0.5_DP
     ft(3,4)=0.5_DP
   CASE(27,37,45)
     ft(3,3)=0.5_DP
     ft(3,4)=0.5_DP
   CASE(28)
     ft(1,2)=0.5_DP
     ft(1,4)=0.5_DP
   CASE(29)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(3,3)=0.5_DP
     ft(1,4)=0.5_DP
   CASE(30)
     ft(1,2)=0.5_DP
     ft(3,3)=0.5_DP
     ft(1,4)=0.5_DP
     ft(3,4)=0.5_DP
   CASE(31)
     ft(1,2)=0.5_DP
     ft(1,4)=0.5_DP
     ft(3,2)=0.5_DP
     ft(3,4)=0.5_DP
   CASE(32)
     ft(1,2)=0.5_DP
     ft(1,4)=0.5_DP
     ft(2,2)=0.5_DP
     ft(2,3)=0.5_DP
   CASE(33)
     ft(1,2)=0.5_DP
     ft(2,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,4)=0.5_DP
     ft(2,3)=0.5_DP
   CASE(34)
     ft(1,2)=0.5_DP
     ft(1,4)=0.5_DP
     ft(2,2)=0.5_DP
     ft(2,3)=0.5_DP
     ft(3,3)=0.5_DP
     ft(3,4)=0.5_DP
   CASE(39)
     ft(3,2)=0.5_DP
     ft(3,3)=0.5_DP
   CASE(40)
     ft(1,2)=0.5_DP
   CASE(41)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(3,3)=0.5_DP
   CASE(43)
     ft(1,2)=0.25_DP
     ft(2,2)=0.25_DP
     ft(3,2)=0.5_DP
     ft(2,3)=0.25_DP
     ft(3,3)=0.25_DP
   CASE(48)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(49)
     ft(3,5)=0.5_DP
   CASE(50)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(51)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
   CASE(52)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(53)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
   CASE(54)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(55)
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
   CASE(56)
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(57)
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(58)
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(59)
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(60)
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(61)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
   CASE(62)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(63)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(64)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
   CASE(66)
     ft(3,5)=0.5_DP
   CASE(67)
     ft(2,5)=0.5_DP
   CASE(68)
     ft(2,5)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(70)
     ft(1,5)=0.25_DP
     ft(2,5)=0.25_DP
     ft(3,5)=0.25_DP
   CASE(72)
     ft(3,5)=0.5_DP
   CASE(73)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
   CASE(74)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
     ft(1,3)=0.5_DP
     ft(2,3)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(76,80)
     ft(3,2)=0.25_DP
   CASE(77,84)
     ft(3,2)=0.5_DP
   CASE(78)
     ft(3,2)=0.75_DP
   CASE(85)
     ft(1,2)=0.5_DP
   CASE(86)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
   CASE(88)
     ft(1,2)=0.25_DP
     ft(2,2)=0.25_DP
     ft(3,2)=0.25_DP
   CASE(90)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(91,98)
     ft(3,2)=0.25_DP
   CASE(92)
     ft(3,2)=0.25_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(93)
     ft(3,2)=0.5_DP
   CASE(94)
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(95)
     ft(3,2)=0.75_DP
   CASE(96)
     ft(3,2)=0.75_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(100)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(101)
     ft(3,2)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(102)
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(103)
     ft(3,5)=0.5_DP
   CASE(104)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(105)
     ft(3,2)=0.5_DP
   CASE(106)
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(108)
     ft(3,5)=0.5_DP
   CASE(109)
     ft(3,2)=0.25_DP
     ft(1,5)=0.5_DP
   CASE(110)
     ft(3,2)=0.25_DP
     ft(2,5)=0.5_DP
   CASE(112,116,120)
     ft(3,5)=0.5_DP
   CASE(113,117)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(114,118)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,5)=0.5_DP
   CASE(122)
     ft(2,5)=0.5_DP
     ft(3,5)=0.25_DP
   CASE(124)
     ft(3,9)=0.5_DP
   CASE(125)
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
   CASE(126)
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(127)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(128)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(129)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
   CASE(130)
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(131)
     ft(3,2)=0.5_DP
   CASE(132)
     ft(3,2)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(133)
     ft(3,2)=0.5_DP
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
   CASE(134)
     ft(3,2)=0.5_DP
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(135)
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
   CASE(136)
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(137)
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
   CASE(138)
     ft(3,2)=0.5_DP
     ft(1,5)=0.5_DP
     ft(2,5)=0.5_DP
     ft(1,9)=0.5_DP
     ft(2,9)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(140)
     ft(3,9)=0.5_DP
   CASE(141)
     ft(3,2)=0.25_DP
     ft(2,9)=0.5_DP
     ft(3,9)=0.5_DP
   CASE(142)
     ft(3,2)=0.25_DP
     ft(2,9)=0.5_DP
   CASE(144)
     ft(3,2)=1.0_DP/3.0_DP
   CASE(145)
     ft(3,2)=2.0_DP/3.0_DP
   CASE(151,152)
     ft(3,2)=1.0_DP/3.0_DP
   CASE(153,154)
     ft(3,2)=2.0_DP/3.0_DP
   CASE(158,159,161)
     ft(3,4)=0.5_DP
   CASE(163,165,167)
     ft(3,7)=0.5_DP
   CASE(169,178)
     ft(3,2)=1.0_DP/6.0_DP
   CASE(170,179)
     ft(3,2)=5.0_DP/6.0_DP
   CASE(171,180)
     ft(3,2)=1.0_DP/3.0_DP
   CASE(172,181)
     ft(3,2)=2.0_DP/3.0_DP
   CASE(173,176,182)
     ft(3,2)=0.5_DP
   CASE(184)
     ft(3,7)=0.5_DP
   CASE(185)
     ft(3,2)=0.5_DP
     ft(3,7)=0.5_DP
   CASE(186)
     ft(3,2)=0.5_DP
   CASE(188,190)
     ft(3,7)=0.5_DP
   CASE(192)
     ft(3,13)=0.5_DP
   CASE(193)
     ft(3,2)=0.5_DP
     ft(3,13)=0.5_DP
   CASE(194)
     ft(3,2)=0.5_DP
   CASE(198,199)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
   CASE(201)
     ft(1,2)=0.5_DP
     ft(2,2)=0.5_DP
   CASE(203)
     ft(1,2)=0.25_DP
     ft(2,2)=0.25_DP
   CASE(205,206)
     ft(1,2)=0.5_DP
     ft(3,2)=0.5_DP
   CASE(208)
     ft(2,18)=-0.5_DP
     ft(3,18)=0.5_DP
   CASE(210,213,214)
     ft(2,18)=-0.25_DP
     ft(3,18)=0.25_DP
   CASE(212)
     ft(2,18)=-0.75_DP
     ft(3,18)=0.75_DP
   CASE(218,219)
     ft(1,18)=0.5_DP
     ft(2,18)=-0.5_DP
     ft(3,18)=0.5_DP
   CASE(220)
     ft(1,18)=0.25_DP
     ft(2,18)=-0.25_DP
     ft(3,18)=0.25_DP
   CASE(222,226)
     ft(1,18)=-0.5_DP
   CASE(223)
     ft(1,18)=-0.5_DP
     ft(2,18)=-0.5_DP
     ft(3,18)=0.5_DP
   CASE(224)
     ft(2,18)=-0.5_DP
     ft(3,18)=0.5_DP
   CASE(227)
     ft(1,18)=-0.5_DP
     ft(2,18)=-0.25_DP
     ft(3,18)=0.25_DP
   CASE(228)
     ft(2,18)=-0.25_DP
     ft(3,18)=0.25_DP
   CASE(230)
     ft(1,18)=-0.25_DP
     ft(2,18)=-0.25_DP
     ft(3,18)=0.25_DP
CASE DEFAULT
   CALL errore('set_ft_generators','unknown space group number',1)
END SELECT

RETURN
END SUBROUTINE set_ft_generators

!-----------------------------------------------------------------------
LOGICAL FUNCTION check_space_group_ft(ft1,ftref,nsym)
!-----------------------------------------------------------------------
!
!  This routine checks if all the fractional translations differ
!  only by integer numbers with respect to those of the standard space group.
!  Returns .true. only if all the fractionary translations match
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER :: nsym
REAL(DP), INTENT(IN) :: ft1(3,48), ftref(3,48)

INTEGER :: isym
REAL(DP) :: eps1 = 1.D-9

check_space_group_ft=.TRUE.

DO isym=1,nsym
   check_space_group_ft=check_space_group_ft.AND. &
     ABS(ft1(1,isym)-ftref(1,isym)-NINT(ft1(1,isym)-ftref(1,isym)))<eps1.AND. &
     ABS(ft1(2,isym)-ftref(2,isym)-NINT(ft1(2,isym)-ftref(2,isym)))<eps1.AND. &
     ABS(ft1(3,isym)-ftref(3,isym)-NINT(ft1(3,isym)-ftref(3,isym)))<eps1
   IF (.NOT. check_space_group_ft) EXIT
END DO

RETURN
END FUNCTION check_space_group_ft

!-----------------------------------------------------------------------
SUBROUTINE set_standard_sg(sgc,ibrav,celldm,nsym,sr,ft)
!-----------------------------------------------------------------------
!
!   This routine receives as input the sgc (1-230), the bravais lattice index
!   and the dimension of the unit cell and sets the symmetry 
!   matrices of the standard space group.  
!   sr are the rotation matrices in cartesian coordinates
!   ft are the fractional translations in crystal coordinates (of the
!      real lattice, not of the conventional one)
!
USE kinds, ONLY : DP
USE point_group, ONLY : set_group_desc, set_sym_o3
USE lattices,  ONLY : compute_conventional
IMPLICIT NONE
INTEGER, INTENT(IN) :: sgc
REAL(DP), INTENT(IN) :: celldm(6)
INTEGER, INTENT(OUT) :: nsym, ibrav
REAL(DP), INTENT(OUT) :: sr(3,3,48), ft(3,48)

REAL(DP), PARAMETER :: eps1=1.D-9
INTEGER :: group_code, group_code_ext, group_desc(48), isym, ipol
REAL(DP) :: omega, at(3,3), bg(3,3), comega, cat(3,3), cbg(3,3), &
            ftc(3,48)
!
!  find the point group of this space group and the extended point group
!  code.
!
CALL set_point_group_code(sgc, group_code, group_code_ext)
!
!  set the rotation matrices
!
CALL set_group_desc(group_desc, nsym, group_code_ext)
DO isym=1,nsym
   CALL set_sym_o3(sr(1,1,isym), group_desc(isym))
END DO
!
CALL set_sg_ibrav(sgc, ibrav)
!
!  set the fractionary translations for the generators
!
ft=0.0_DP
CALL set_ft_generators(sgc,ft)
!
!  now calculate the primitive lattice vectors and the reciprocal
!  lattice vectors. Computes also the conventional lattice vector and
!  the reciprocal conventional lattice vectors
!
CALL latgen(ibrav,celldm,at(1,1), at(1,2), at(1,3), omega)
CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
CALL compute_conventional(at, cat, ibrav)
CALL recips(cat(1,1), cat(1,2), cat(1,3), cbg(1,1), cbg(1,2), cbg(1,3))
!
!  The fractional translations are in the basis of the conventional
!  lattice vectors and are written in cartesian coordinates
!
ftc=ft
CALL cryst_to_cart(nsym,ftc,cat,1)
!
!  here we build the matrices of the space group and the fractional
!  translations from the generators
!
SELECT CASE (group_code)
   CASE(1,2,3,4)

   CASE(5)
!
!  C_3 
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))

   CASE(6,26)
!
!  C_4, S_4 
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                  sr(1,1,4), ftc(1,4))
   CASE(7,17)
!
!  C_6, C_3h
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                  sr(1,1,4), ftc(1,4))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,4), ftc(1,4), &
                  sr(1,1,5), ftc(1,5))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,5), ftc(1,5), &
                  sr(1,1,6), ftc(1,6))
   CASE(8,12,16)
!
!   D_2, C_2v, C_2h
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                   sr(1,1,4), ftc(1,4))
   CASE(9,13)
!
!   D_3, C_3v
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,4), ftc(1,4), &
                  sr(1,1,5), ftc(1,5))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,4), ftc(1,4), &
                  sr(1,1,6), ftc(1,6))

   CASE(10,14,18,24)
!
!   D_4, C_4v, C_4h, D_2d
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                  sr(1,1,4), ftc(1,4))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,5), ftc(1,5), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,5), ftc(1,5), &
                  sr(1,1,7), ftc(1,7))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,5), ftc(1,5), &
                  sr(1,1,8), ftc(1,8))

   CASE(11,15,19,21)
!
!   D_6, C_6v, C_6h, D_3h
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                  sr(1,1,4), ftc(1,4))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,4), ftc(1,4), &
                  sr(1,1,5), ftc(1,5))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,5), ftc(1,5), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,7), ftc(1,7), &
                  sr(1,1,8), ftc(1,8))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,7), ftc(1,7), &
                  sr(1,1,9), ftc(1,9))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,7), ftc(1,7), &
                  sr(1,1,10), ftc(1,10))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,7), ftc(1,7), &
                  sr(1,1,11), ftc(1,11))
      CALL space_group_multiply(sr(1,1,6), ftc(1,6), sr(1,1,7), ftc(1,7), &
                  sr(1,1,12), ftc(1,12))

   CASE(20)
!
!  D_2h
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                  sr(1,1,4), ftc(1,4))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,5), ftc(1,5), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,5), ftc(1,5), &
                  sr(1,1,7), ftc(1,7))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,5), ftc(1,5), &
                  sr(1,1,8), ftc(1,8))
   CASE(22)
!
!  D_4h
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                  sr(1,1,4), ftc(1,4))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,5), ftc(1,5), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,5), ftc(1,5), &
                  sr(1,1,7), ftc(1,7))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,5), ftc(1,5), &
                  sr(1,1,8), ftc(1,8))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,9), ftc(1,9), &
                  sr(1,1,10), ftc(1,10))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,9), ftc(1,9), &
                  sr(1,1,11), ftc(1,11))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,9), ftc(1,9), &
                  sr(1,1,12), ftc(1,12))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,9), ftc(1,9), &
                  sr(1,1,13), ftc(1,13))
      CALL space_group_multiply(sr(1,1,6), ftc(1,6), sr(1,1,9), ftc(1,9), &
                  sr(1,1,14), ftc(1,14))
      CALL space_group_multiply(sr(1,1,7), ftc(1,7), sr(1,1,9), ftc(1,9), &
                  sr(1,1,15), ftc(1,15))
      CALL space_group_multiply(sr(1,1,8), ftc(1,8), sr(1,1,9), ftc(1,9), &
                  sr(1,1,16), ftc(1,16))

   CASE(23)
!
!  D_6h
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,3), ftc(1,3), &
                  sr(1,1,4), ftc(1,4))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,4), ftc(1,4), &
                  sr(1,1,5), ftc(1,5))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,5), ftc(1,5), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,7), ftc(1,7), &
                  sr(1,1,8), ftc(1,8))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,7), ftc(1,7), &
                  sr(1,1,9), ftc(1,9))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,7), ftc(1,7), &
                  sr(1,1,10), ftc(1,10))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,7), ftc(1,7), &
                  sr(1,1,11), ftc(1,11))
      CALL space_group_multiply(sr(1,1,6), ftc(1,6), sr(1,1,7), ftc(1,7), &
                  sr(1,1,12), ftc(1,12))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,13), ftc(1,13), &
                  sr(1,1,14), ftc(1,14))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,13), ftc(1,13), &
                  sr(1,1,15), ftc(1,15))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,13), ftc(1,13), &
                  sr(1,1,16), ftc(1,16))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,13), ftc(1,13), &
                  sr(1,1,17), ftc(1,17))
      CALL space_group_multiply(sr(1,1,6), ftc(1,6), sr(1,1,13), ftc(1,13), &
                  sr(1,1,18), ftc(1,18))
      CALL space_group_multiply(sr(1,1,7), ftc(1,7), sr(1,1,13), ftc(1,13), &
                  sr(1,1,19), ftc(1,19))
      CALL space_group_multiply(sr(1,1,8), ftc(1,8), sr(1,1,13), ftc(1,13), &
                  sr(1,1,20), ftc(1,20))
      CALL space_group_multiply(sr(1,1,9), ftc(1,9), sr(1,1,13), ftc(1,13), &
                  sr(1,1,21), ftc(1,21))
      CALL space_group_multiply(sr(1,1,10), ftc(1,10), sr(1,1,13), ftc(1,13), &
                  sr(1,1,22), ftc(1,22))
      CALL space_group_multiply(sr(1,1,11), ftc(1,11), sr(1,1,13), ftc(1,13), &
                  sr(1,1,23), ftc(1,23))
      CALL space_group_multiply(sr(1,1,12), ftc(1,12), sr(1,1,13), ftc(1,13), &
                  sr(1,1,24), ftc(1,24))
   CASE(25)
!
!  D_3d
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,4), ftc(1,4), &
                  sr(1,1,5), ftc(1,5))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,4), ftc(1,4), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,7), ftc(1,7), &
                  sr(1,1,8), ftc(1,8))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,7), ftc(1,7), &
                  sr(1,1,9), ftc(1,9))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,7), ftc(1,7), &
                  sr(1,1,10), ftc(1,10))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,7), ftc(1,7), &
                  sr(1,1,11), ftc(1,11))
      CALL space_group_multiply(sr(1,1,6), ftc(1,6), sr(1,1,7), ftc(1,7), &
                  sr(1,1,12), ftc(1,12))

   CASE(27)
!
!  S_6
!
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,2), ftc(1,2), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,4), ftc(1,4), &
                  sr(1,1,5), ftc(1,5))
      CALL space_group_multiply(sr(1,1,3), ftc(1,3), sr(1,1,4), ftc(1,4), &
                  sr(1,1,6), ftc(1,6))
   CASE(28,29)
!
!  T, T_h
!
!  Note that here the generators are 2_z (2) and 3xyz (5)  
!
!  1  2   3      4     5  6    7    8   9    10   11   12
!  E  A  B^2AB  BAB^2  B  AB  ABA  BA  B^2  B^2A  BAB  AB^2
!

      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,5), ftc(1,5), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,2), ftc(1,2), &
                  sr(1,1,8), ftc(1,8))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,5), ftc(1,5), &
                  sr(1,1,9), ftc(1,9))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,9), ftc(1,9), &
                  sr(1,1,12), ftc(1,12))
      CALL space_group_multiply(sr(1,1,9), ftc(1,9), sr(1,1,2), ftc(1,2), &
                  sr(1,1,10), ftc(1,10))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,8), ftc(1,8), &
                  sr(1,1,7), ftc(1,7))
      CALL space_group_multiply(sr(1,1,8), ftc(1,8), sr(1,1,5), ftc(1,5), &
                  sr(1,1,11), ftc(1,11))
      CALL space_group_multiply(sr(1,1,9), ftc(1,9), sr(1,1,6), ftc(1,6), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,8), ftc(1,8), sr(1,1,9), ftc(1,9), &
                  sr(1,1,4), ftc(1,4))

      IF (group_code==29) THEN
         DO isym=2,12
            CALL space_group_multiply(sr(1,1,isym), ftc(1,isym), &
                      sr(1,1,13), ftc(1,13), sr(1,1,12+isym), ftc(1,12+isym))
         END DO 
      ENDIF
   CASE(30,31,32)
!
!  T_d, O, O_h
!
!    1   2        3     4  5   6       7        8    9   10     11     12
!    E BA^2B-1  ABA-1B A^2 B A-1A-1B B-1ABA-1 AB-1A B-1 B-1A^2 A^2B-1 A-1BA-1
!
!     13    14    15   16   17 18  19     20  21       22   23       24
!    B-1A  A-1B  AB-1 BA-1 A-1 A A^2B-1A  AB A^2B-1AB  BA BA^2B-1A  BA-1B

      CALL space_group_multiply(sr(1,1,18), ftc(1,18), sr(1,1,5), ftc(1,5), &
                  sr(1,1,20), ftc(1,20))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,18), ftc(1,18), &
                  sr(1,1,22), ftc(1,22))
      CALL space_group_multiply(sr(1,1,18), ftc(1,18), sr(1,1,18), ftc(1,18), &
                  sr(1,1,4), ftc(1,4))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,5), ftc(1,5), &
                  sr(1,1,9), ftc(1,9))
      CALL space_group_multiply(sr(1,1,18), ftc(1,18), sr(1,1,4), ftc(1,4), &
                  sr(1,1,17), ftc(1,17))
      CALL space_group_multiply(sr(1,1,17), ftc(1,17), sr(1,1,5), ftc(1,5), &
                  sr(1,1,14), ftc(1,14))
      CALL space_group_multiply(sr(1,1,18), ftc(1,18), sr(1,1,9), ftc(1,9), &
                  sr(1,1,15), ftc(1,15))
      CALL space_group_multiply(sr(1,1,5), ftc(1,5), sr(1,1,17), ftc(1,17), &
                  sr(1,1,16), ftc(1,16))
      CALL space_group_multiply(sr(1,1,9), ftc(1,9), sr(1,1,18), ftc(1,18), &
                  sr(1,1,13), ftc(1,13))
      CALL space_group_multiply(sr(1,1,14), ftc(1,14), sr(1,1,17), ftc(1,17), &
                  sr(1,1,12), ftc(1,12))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,9), ftc(1,9), &
                  sr(1,1,11), ftc(1,11))
      CALL space_group_multiply(sr(1,1,9), ftc(1,9), sr(1,1,4), ftc(1,4), &
                  sr(1,1,10), ftc(1,10))
      CALL space_group_multiply(sr(1,1,16), ftc(1,16), sr(1,1,5), ftc(1,5), &
                  sr(1,1,24), ftc(1,24))
      CALL space_group_multiply(sr(1,1,15), ftc(1,15), sr(1,1,18), ftc(1,18), &
                  sr(1,1,8), ftc(1,8))
      CALL space_group_multiply(sr(1,1,13), ftc(1,13), sr(1,1,16), ftc(1,16), &
                  sr(1,1,7), ftc(1,7))
      CALL space_group_multiply(sr(1,1,17), ftc(1,17), sr(1,1,14), ftc(1,14), &
                  sr(1,1,6), ftc(1,6))
      CALL space_group_multiply(sr(1,1,20), ftc(1,20), sr(1,1,14), ftc(1,14), &
                  sr(1,1,3), ftc(1,3))
      CALL space_group_multiply(sr(1,1,22), ftc(1,22), sr(1,1,15), ftc(1,15), &
                  sr(1,1,2), ftc(1,2))
      CALL space_group_multiply(sr(1,1,4), ftc(1,4), sr(1,1,13), ftc(1,13), &
                  sr(1,1,19), ftc(1,19))
      CALL space_group_multiply(sr(1,1,11), ftc(1,11), sr(1,1,20), ftc(1,20), &
                  sr(1,1,21), ftc(1,21))
      CALL space_group_multiply(sr(1,1,2), ftc(1,2), sr(1,1,18), ftc(1,18), &
                  sr(1,1,23), ftc(1,23))
      IF (group_code==32) THEN
         DO isym=2,24
            CALL space_group_multiply(sr(1,1,isym), ftc(1,isym), &
                      sr(1,1,25), ftc(1,25), sr(1,1,24+isym), ftc(1,24+isym))
         END DO 
      ENDIF

END SELECT
!
!   Bring all fractional translations in the crystal basis of the
!   real lattice
!
ft=ftc
CALL cryst_to_cart(nsym,ft,bg,-1)
!
!  remove pure translations and set translations between 0 < x <= 1
!
DO isym=1,nsym
   ft(:,isym)=ft(:,isym)-NINT(ft(:,isym))
   DO ipol=1,3
      IF (ft(ipol,isym) < -eps1) ft(ipol,isym)=ft(ipol,isym)+1.0_DP
   ENDDO
ENDDO

RETURN
END SUBROUTINE set_standard_sg

!-----------------------------------------------------------------------
SUBROUTINE space_group_multiply(sr1, ft1, sr2, ft2, sr3, ft3)
!-----------------------------------------------------------------------
!
!   This routine multiplies two space group operations
!   sr1, sr2 are two 3x3 rotation matrices in cartesian coordinates
!   ft1, ft2, are two fractional translations in cartesian coordinates
!
IMPLICIT NONE
REAL(DP), INTENT(IN)  :: sr1(3,3), sr2(3,3)
REAL(DP), INTENT(IN)  :: ft1(3), ft2(3)
REAL(DP), INTENT(OUT) :: sr3(3,3)
REAL(DP), INTENT(OUT) :: ft3(3)

REAL(DP) :: ftaux(3)
INTEGER  :: ipol

sr3 = MATMUL(sr1, sr2)

ftaux=0.0_DP
DO ipol=1,3
   ftaux(:) = ftaux(:) + sr1(:,ipol) * ft2(ipol)
END DO
ft3 = ft1 + ftaux

RETURN
END SUBROUTINE space_group_multiply

!-----------------------------------------------------------------------
SUBROUTINE find_reference_pg(code_group_ext, ref_code_group_ext)
!-----------------------------------------------------------------------
!
!  For each point group (in extended form), gives the point group of
!  the reference space group for that point group
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: code_group_ext
INTEGER, INTENT(OUT) :: ref_code_group_ext

INTEGER :: rfc(136)

DATA  rfc / 1,    2,    2,    2,    2,    2,    2,    2,    2,    2,  &
            2,    2,    2,    2,   15,   15,   15,   15,   15,   15,  &
           15,   15,   15,   15,   15,   15,   15,   28,   33,   33,  &
           33,   33,   33,   34,   34,   34,   37,   38,   38,   38,  &
           38,   38,   38,   44,   45,   44,   44,   44,   44,   50,  &
           50,   50,   53,   58,   58,   58,   58,   58,   58,   58,  &
           58,   58,   58,   58,   58,   58,   58,   58,   58,   58,  &
           58,   72,   73,   72,   72,   72,   72,   78,   78,   78,  &
           81,   82,   82,   82,   82,   82,   82,   82,   82,   82,  &
           82,   82,   82,   82,   95,   96,   96,   96,   99,  100,  &
          100,  100,  100,  100,  100,  106,  107,  108,  108,  108,  &
          111,  112,  113,  112,  112,  112,  112,  118,  119,  118,  &
          118,  118,  118,  124,  124,  124,  127,  127,  127,  127,  &
          127,  132,  133,  134,  135,  136  /


ref_code_group_ext=rfc(code_group_ext)

RETURN
END SUBROUTINE find_reference_pg

!-----------------------------------------------------------------------
LOGICAL FUNCTION check_code_group_ext(cge)
!-----------------------------------------------------------------------
!
!  this function gives .TRUE. if the point group is one of those compatible
!  with the ITA tables
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: cge
!
!  c unique, b unique, and all settings for the orthorhombic groups
!
check_code_group_ext=(cge==1.OR.cge==2.OR.cge==3.OR.cge==15.OR.cge==16      &
               .OR.cge==28.OR.cge==33.OR.cge==34.OR.cge==37.OR.cge==38      &
               .OR.cge==44.OR.cge==45.OR.cge==50.OR.cge==53.OR.cge==54      &
               .OR.cge==56.OR.cge==58.OR.cge==72.OR.cge==73.OR.cge==78      &
               .OR.cge==81.OR.cge==82.OR.cge==83.OR.cge==95.OR.cge==96      &
               .OR.cge==99.OR.cge==100.OR.cge==106.OR.cge==107.OR.cge==108  &
               .OR.cge==111.OR.cge==112.OR.cge==113.OR.cge==118.OR.cge==119 &
               .OR.cge==124.OR.cge==127.OR.cge==132.OR.cge==133.OR.cge==134 &
               .OR.cge==135.OR.cge==136)
RETURN
END FUNCTION check_code_group_ext

!-----------------------------------------------------------------------
LOGICAL FUNCTION symmorphic_sg(sgc,ssgc)
!-----------------------------------------------------------------------
!
!   This function receives a space group number (1-230) and returns .TRUE.
!   if the space group is symmorphic. It returns also the group number
!   of the symmorphic space group obtained by removing the fractional
!   translations.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: sgc
INTEGER, INTENT(OUT) :: ssgc

INTEGER :: ssgc_data(230)

DATA ssgc_data  / 1,    2,    3,     3,     5,    6,    6,    8,    8,   10, &
                 10,   12,   12,    10,    12,   16,   16,   16,   16,   21, &
                 21,   22,   23,    23,    25,   25,   25,   25,   25,   25, &
                 25,   25,   25,    25,    35,   35,   35,   38,   38,   38, &
                 38,   42,   42,    44,    44,   44,   47,   47,   47,   47, &
                 47,   47,   47,    47,    47,   47,   47,   47,   47,   47, &
                 47,   47,   65,    65,    65,   65,   65,   65,   69,   69, &
                 71,   71,   71,    71,    75,   75,   75,   75,   79,   79, &
                 81,   82,   83,    83,    83,   83,   87,   87,   89,   89, &
                 89,   89,   89,    89,    89,   89,   97,   97,   99,   99, &
                 99,   99,   99,    99,    99,   99,   107, 107,  107,  107, &
                111,  111,  111,   111,   115,  115,   115, 115,  119,  119, &
                121,  121,  123,   123,   123,  123,   123, 123,  123,  123, &
                123,  123,  123,   123,   123,  123,   123, 123,  139,  139, &
                139,  139,  143,   143,   143,  146,   147, 148,  149,  150, &
                149,  150,  149,   150,   155,  156,   157, 156,  157,  160, &
                160,  162,  162,   164,   164,  166,   166, 168,  168,  168, &
                168,  168,  168,   174,   175,  175,   177, 177,  177,  177, &
                177,  177,  183,   183,   183,  183,   187, 187,  189,  189, &
                191,  191,  191,   191,   195,  196,   197, 195,  197,  200, &
                200,  202,  202,   204,   200,  204,   207, 207,  209,  209, &
                211,  207,  207,   211,   215,  216,   217, 215,  216,  217, &
                221,  221,  221,   221,   225,  225,   225, 225,  229,  229  /
                 
ssgc=ssgc_data(sgc)
symmorphic_sg=(sgc==ssgc)

RETURN
END FUNCTION symmorphic_sg

!-----------------------------------------------------------------------
LOGICAL FUNCTION check_intersection(ft1,ft2,ft3)
!-----------------------------------------------------------------------

IMPLICIT NONE
REAL(DP), INTENT(IN) :: ft1(3), ft2(3), ft3(3)

LOGICAL :: type1, type2, type3
INTEGER :: ipol
REAL(DP) :: ft1_(3), ft2_(3), ft3_(3)
REAL(DP), PARAMETER :: eps1=1.D-8
!
!   first bring the fractional translations to be 0.0 or 0.5
!
DO ipol=1,3
   IF (ABS(ft1(ipol)-NINT(ft1(ipol)))<eps1) THEN
      ft1_(ipol)=0.0_DP
   ELSE
      ft1_(ipol)=0.5_DP
   ENDIF
   IF (ABS(ft2(ipol)-NINT(ft2(ipol)))<eps1) THEN
      ft2_(ipol)=0.0_DP
   ELSE
      ft2_(ipol)=0.5_DP
   ENDIF
   IF (ABS(ft3(ipol)-NINT(ft3(ipol)))<eps1) THEN
      ft3_(ipol)=0.0_DP
   ELSE
      ft3_(ipol)=0.5_DP
   ENDIF
END DO
!
!  If they are all 0.0,0.0 or 0.5,0.5  check_intersection is .TRUE.
!
type1=(ft1_(1)==ft1_(2))
type2=(ft2_(2)==ft2_(3))
type3=(ft3_(1)==ft3_(3))

check_intersection=type1.AND.type2.AND.type3

RETURN
END FUNCTION check_intersection

!-----------------------------------------------------------------------
SUBROUTINE sg_origin(sg_number, spaceg_name, at, s01, s02)
!-----------------------------------------------------------------------

USE kinds, ONLY : DP
USE io_global, ONLY : stdout

IMPLICIT NONE
INTEGER, INTENT(IN) :: sg_number
REAL(DP), INTENT(IN) :: at(3,3), s01(3), s02(3)
CHARACTER(LEN=12) :: spaceg_name

REAL(DP) :: s01c(3), s02c(3), s01mod, s02mod

IF (sg_number > 0) THEN
   WRITE(stdout,'(/,5x,"Space group ",a,"   (group number",i4, ").")') &
                         TRIM(spaceg_name), sg_number
   s01mod=s01(1)**2 + s01(2)**2 + s01(3)**2
   s02mod=s02(1)**2 + s02(2)**2 + s02(3)**2
   s01c=s01
   s02c=s02
   CALL cryst_to_cart(1,s01c,at,1)
   CALL cryst_to_cart(1,s02c,at,1)
   IF (s01mod < 1.D-4.AND.s02mod>1.D3) THEN
      WRITE(stdout,'(5x,"The origin coincides with the ITA tables.")')
   ELSEIF (s01mod<1.D-4.AND.s02mod<1.D3) THEN
      WRITE(stdout,'(5x,"The origin coincides with &
                   &the ITA tables (origin choice 1).")')
      WRITE(stdout,'(/,5x,"To shift to origin choice 2 &
                         &subtract to all coordinates:")')
      WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
      WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s02(:), s02c(:)
   ELSEIF (s01mod > 1.D-4 .AND. s02mod<1.D-4 ) THEN
      WRITE(stdout,'(5x,"The origin coincides with &
                  & the ITA tables (origin choice 2).")')
      WRITE(stdout,'(/,5x,"To shift to origin choice 1 &
                         &subtract to all coordinates:")')
      WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
      WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s01(:), s01c(:)
   ELSEIF (s01mod > 1.D-4 .AND. s01mod<1.d3 .AND. s02mod > 1.d3 ) THEN
      WRITE(stdout,'(5x,"The origin does not coincide with &
                  &the ITA tables.")')
      WRITE(stdout,'(/,5x,"To have the origin as in the ITA &
                     &tables subtract to all coordinates:")')
      WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
      WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s01(:), s01c(:)
   ELSEIF (s01mod > 1.D-4 .AND. s02mod < 1.d3 ) THEN
      WRITE(stdout,'(5x,"The origin does not coincide with the ITA tables.")')
      WRITE(stdout,'(/,5x,"To shift to origin choice 1 &
                           &subtract to all coordinates:")')
      WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
      WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s01(:), s01c(:)
      WRITE(stdout,'(/,5x,"to shift to origin choice 2 &
                             &subtract to all coordinates: ")')
      WRITE(stdout,'(17x,"(cryst. coord.)",18x,"(cart. coord.)")')
      WRITE(stdout,'(5x,3f12.7,2x,3f12.7)') s02(:), s02c(:)
   ENDIF
ELSE
   WRITE(stdout,'(/,5x,"Unknown space group.")') 
ENDIF

RETURN
END SUBROUTINE sg_origin

END MODULE space_groups
