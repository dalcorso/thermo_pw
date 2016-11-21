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
!  In the case of orthorombic groups it gives the rotation needed to
!  transform the group to the one found in the ITA tables.
!
!  The point group and the Bravais lattice must be compatible, otherwise
!  this routine will not work.
!
!  Limitation : there are still problems with origin shifts. The routine 
!               should recognize the space group if the origin coincides 
!               with that used on the ITA tables, but a shifted origin 
!               might still confuse the routine in particular space groups. 
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

  INTEGER, PARAMETER :: nsg = 1080             ! number of space groups
  INTEGER, PARAMETER :: nsg_2d = 20            ! number of 2d space groups
  CHARACTER(LEN=13) :: space_group_names(nsg)
  CHARACTER(LEN=13) :: space_group_names_2d(nsg_2d)
  CHARACTER(LEN=60) :: add_info(nsg)=""
  INTEGER :: space_group_numbers(nsg), space_group_numbers_2d(nsg_2d)

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

  DATA space_group_names  &
       / "P1           ", "P-1          ", "P121         ", "P12_11       ", &
         "C121         ", "P1m1         ", "P1c1         ", "C1m1         ", &
         "C1c1         ", "P12/m1       ", "P12_1/m1     ", "C12/m1       ", &
         "P12/c1       ", "P12_1/c1     ", "C12/c1       ", "B121         ", &
         "P112         ", "C112         ", "P211         ", "A211         ", &
         "B12_11       ", "P112_1       ", "C112         ", "P2_111       ", &
         "A2_111       ", "A121         ", "I121         ", "A112         ", &
         "B112         ", "I112         ", "B211         ", "C211         ", &
         "I211         ", "B1m1         ", "P11m         ", "C11m         ", &
         "Pm11         ", "Am11         ", "P1a1         ", "P1n1         ", &
         "P11a         ", "P11b         ", "P11n         ", "Pb11         ", &
         "Pc11         ", "Pn11         ", "A1m1         ", "I1m1         ", &
         "A11m         ", "B11m         ", "I11m         ", "Bm11         ", &
         "Cm11         ", "Im11         ", "A1a1         ", "I1a1         ", &
         "A11a         ", "B11b         ", "I11a         ", "Bb11         ", &
         "Cc11         ", "Ib11         ", "B12/m1       ", "P112/m       ", &
         "C112/m       ", "P2/m11       ", "A2/m11       ", "B12_1/m1     ", &
         "P112_1/m     ", "C112_1/m     ", "P2_1/m11     ", "A2_1/m11     ", &
         "A12/m1       ", "I12/m1       ", "A112/m       ", "B112/m       ", &
         "I112/m       ", "B2/m11       ", "C2/m11       ", "I2/m11       ", &
         "P12/a1       ", "P12/n1       ", "P112/a       ", "P112/b       ", &
         "P112/n       ", "P2/b11       ", "P2/c11       ", "I2/m11       ", &
         "P12_1/a1     ", "P12_1/n1     ", "             ", "             ", &
         "P112_1/a     ", "P112_1/b     ", "P112_1/n     ", "P2_1/b11     ", &
         "P2_1/c11     ", "P2_1/n11     ", "A12/a1       ", "I12/a1       ", &
         "A112/a       ", "B112/b       ", "I112/a       ", "B2/b11       ", &
         "C2/c11       ", "I2/b11       ", "P222         ", "P222_1       ", &
         "P2_12_12     ", "P2_12_12_1   ", "C222_1       ", "C222         ", &
         "F222         ", "I222         ", "I2_12_12_1   ", "Pmm2         ", &
         "Pmc2_1       ", "Pcc2         ", "Pma2         ", "Pca2_1       ", &
         "Pnc2         ", "Pmn2_1       ", "Pba2         ", "Pna2_1       ", &
         "Pnn2         ", "Cmm2         ", "Cmc2_1       ", "Ccc2         ", &
         "Amm2         ", "Aem2         ", "Ama2         ", "Aea2         ", &
         "Fmm2         ", "Fdd2         ", "Imm2         ", "Iba2         ", &
         "Ima2         ", "Pmmm         ", "Pnnn         ", "Pccm         ", &
         "Pban         ", "Pmma         ", "Pnna         ", "Pmna         ", &
         "Pcca         ", "Pbam         ", "Pccn         ", "Pbcm         ", &
         "Pnnm         ", "Pmmn         ", "Pbcn         ", "Pbca         ", &
         "Pnma         ", "Cmcm         ", "Cmce         ", "Cmmm         ", &
         "Cccm         ", "Cmme         ", "Ccce         ", "Fmmm         ", &
         "Fddd         ", "Immm         ", "Ibam         ", "Ibca         ", &
         "Imma         ", "P4           ", "P4_1         ", "P4_2         ", &
         "P4_3         ", "I4           ", "I4_1         ", "P-4          ", &
         "I-4          ", "P4/m         ", "P4_2/m       ", "P4/n         ", &
         "P4_2/n       ", "I4/m         ", "I4_1/a       ", "P422         ", &
         "P42_12       ", "P4_122       ", "P4_12_12     ", "P4_222       ", &
         "P4_22_12     ", "P4_322       ", "P4_32_12     ", "I422         ", &
         "I4_122       ", "P4mm         ", "P4bm         ", "P4_2cm       ", &
         "P4_2nm       ", "P4cc         ", "P4nc         ", "P4_2mc       ", &
         "P4_2bc       ", "I4mm         ", "I4cm         ", "I4_1md       ", &
         "I4_1cd       ", "P-42m        ", "P-42c        ", "P-42_1m      ", &
         "P-42_1c      ", "P-4m2        ", "P-4c2        ", "P-4b2        ", &
         "P-4n2        ", "I-4m2        ", "I-4c2        ", "I-42m        ", &
         "I-42d        ", "P4/mmm       ", "P4/mcc       ", "P4/nbm       ", &
         "P4/nnc       ", "P4/mbm       ", "P4/mnc       ", "P4/nmm       ", &
         "P4/ncc       ", "P4_2/mmc     ", "P4_2/mcm     ", "P4_2/nbc     ", &
         "P4_2/nnm     ", "P4_2/mbc     ", "P4_2/mnm     ", "P4_2/nmc     ", &
         "P4_2/ncm     ", "I4/mmm       ", "I4/mcm       ", "I4_1/amd     ", &
         "I4_1/acd     ", "P3           ", "P3_1         ", "P3_2         ", &
         "R3           ", "P-3          ", "R-3          ", "P312         ", &
         "P321         ", "P3_112       ", "P3_121       ", "P3_212       ", &
         "P3_221       ", "R32          ", "P3m1         ", "P31m         ", &
         "P3c1         ", "P31c         ", "R3m          ", "R3c          ", &
         "P-31m        ", "P-31c        ", "P-3m1        ", "P-3c1        ", &
         "R-3m         ", "R-3c         ", "P6           ", "P6_1         ", &
         "P6_5         ", "P6_2         ", "P6_4         ", "P6_3         ", &
         "P-6          ", "P6/m         ", "P6_3/m       ", "P622         ", &
         "P6_122       ", "P6_522       ", "P6_222       ", "P6_422       ", &
         "P6_322       ", "P6mm         ", "P6cc         ", "P6_3cm       ", &
         "P6_3mc       ", "P-6m2        ", "P-6c2        ", "P-62m        ", &
         "P-62c        ", "P6/mmm       ", "P6/mcc       ", "P6_3/mcm     ", &
         "P6_3/mmc     ", "P23          ", "F23          ", "I23          ", &
         "P2_13        ", "I2_13        ", "Pm-3         ", "Pn-3         ", &
         "Fm-3         ", "Fd-3         ", "Im-3         ", "Pa-3         ", &
         "Ia-3         ", "P432         ", "P4_232       ", "F432         ", &
         "F4_132       ", "I432         ", "P4_332       ", "P4_132       ", &
         "I4_132       ", "P-43m        ", "F-43m        ", "I-43m        ", &
         "P-43n        ", "F-43c        ", "I-43d        ", "Pm-3m        ", &
         "Pn-3n        ", "Pm-3n        ", "Pn-3m        ", "Fm-3m        ", &
         "Fm-3c        ", "Fd-3m        ", "Fd-3c        ", "Im-3m        ", &
         "Ia-3d        ", "P2_122       ", "P22_12       ", "P2_122_1     ", &
         "P22_12_1     ", "B22_12       ", "A2_122       ", "B222         ", &
         "A222         ", "Pm2m         ", "P2mm         ", "Pcm2_1       ", &
         "Pb2_1m       ", "Pm2_1b       ", "P2_1ma       ", "P2_1am       ", &
         "Pb2b         ", "P2aa         ", "Pbm2         ", "Pc2m         ", &
         "Pm2a         ", "P2mb         ", "P2cm         ", "Pbc2_1       ", &
         "Pc2_1b       ", "Pb2_1a       ", "P2_1ab       ", "P2_1ca       ", &
         "Pcn2         ", "Pb2n         ", "Pn2b         ", "P2na         ", &
         "P2an         ", "Pnm2_1       ", "Pn2_1m       ", "Pm2_1n       ", &
         "P2_1mn       ", "P2_1nm       ", "Pc2a         ", "P2cb         ", &
         "Pbn2_1       ", "Pc2_1n       ", "Pn2_1a       ", "P2_1nb       ", &
         "P2_1cn       ", "Pn2n         ", "P2nn         ", "Bm2m         ", &
         "A2mm         ", "Ccm2_1       ", "Bb2_1m       ", "Bm2_1b       ", &
         "A2_1ma       ", "A2_1am       ", "Bb2b         ", "A2aa         ", &
         "Bmm2         ", "Cm2m         ", "Am2m         ", "B2mm         ", &
         "C2mm         ", "Bme2         ", "Cm2e         ", "Ae2m         ", &
         "B2em         ", "C2me         ", "Bbm2         ", "Cc2m         ", &
         "Am2a         ", "B2mb         ", "C2cm         ", "Bbe2         ", &
         "Cc2e         ", "Ae2a         ", "B2eb         ", "C2ce         ", &
         "Fd2d         ", "F2dd         ", "Ic2a         ", "I2cb         ", &
         "Ibm2         ", "Ic2m         ", "Im2a         ", "I2mb         ", &
         "I2cm         ", "Pbmb         ", "Pmaa         ", "Pcna         ", &
         "             ", "Pncb         ", "             ", "Pmmb         ", &
         "Pmcm         ", "Pmam         ", "Pbmm         ", "Pcmm         ", &
         "Pnnb         ", "Pncn         ", "Pnan         ", "Pbnn         ", &
         "Pcnn         ", "Pnmb         ", "Pncm         ", "Pman         ", &
         "Pbmn         ", "Pcnm         ", "Pccb         ", "Pbcb         ", &
         "Pbab         ", "Pbaa         ", "Pcaa         ", "Pcma         ", &
         "Pmcb         ", "Pbnb         ", "Pnaa         ", "Pcam         ", &
         "Pbma         ", "Pcmb         ", "Pmca         ", "Pmab         ", &
         "Pnmn         ", "Pmnn         ", "Pmnm         ", "Pnmm         ", &
         "Pcan         ", "Pbna         ", "Pcnb         ", "Pnca         ", &
         "Pnab         ", "Pcab         ", "Pmnb         ", "Pmcn         ", &
         "Pnam         ", "Pbnm         ", "Pcmn         ", "Ccmm         ", &
         "Bbmm         ", "Bmmb         ", "Amma         ", "Amam         ", &
         "Ccme         ", "Bbem         ", "Bmeb         ", "Aema         ", &
         "Aeam         ", "Bbmb         ", "Amaa         ", "Bmam         ", &
         "Abmm         ", "Bbeb         ", "Aeaa         ", "Icma         ", &
         "Imcb         ", "Immb         ", "Imcm         ", "A1           ", &
         "B1           ", "C1           ", "F1           ", "I1           ", &
         "A-1          ", "B-1          ", "C-1          ", "F-1          ", &
         "I-1          ", "C4           ", "C4_1         ", "C4_2         ", &
         "C4_3         ", "F4           ", "F4_1         ", "C-4          ", &
         "F-4          ", "C4/m         ", "C4_2/m       ", "C4/a         ", &
         "C4_2/a       ", "F4/m         ", "F4_1/d       ", "C422         ", &
         "C422_1       ", "C4_122       ", "C4_122_1     ", "C4_222       ", &
         "C4_222_1     ", "C4_322       ", "C4_322_1     ", "F422         ", &
         "F4_122       ", "C4mm         ", "C4mb         ", "C4_2mc       ", &
         "C4_2mn       ", "C4cc         ", "C4cn         ", "C4_2cm       ", &
         "C4_2cb       ", "F4mm         ", "F4mc         ", "F4_1dm       ", &
         "F4_1dc       ", "C-4m2        ", "C-4c2        ", "C-4m2_1      ", &
         "C-4c2_1      ", "C-42m        ", "C-42c        ", "C-42b        ", &
         "C-42n        ", "F-42m        ", "F-42c        ", "F-4m2        ", &
         "F-4d2        ", "C4/mmm       ", "C4/mcc       ", "C4/amb       ", &
         "C4/acn       ", "C4/mmb       ", "C4/mcn       ", "C4/amm       ", &
         "C4/acc       ", "C4_2/mcm     ", "C4_2/mmc     ", "C4_2/acb     ", &
         "C4_2/anm     ", "C4_2/mcb     ", "C4_2/mmn     ", "C4_2/acm     ", &
         "C4_2/amc     ", "F4/mmm       ", "F4/mmc       ", "F4_1/ddm     ", &
         "F4_1/ddc     ", "H3         * ", "H3_1       * ", "H3_2       * ", &
         "H-3        * ", "H321       * ", "H312       * ", "H3_121     * ", &
         "H3_112     * ", "H3_221     * ", "H3_212     * ", "H31m       * ", &
         "H3m1       * ", "H31c       * ", "H3c1       * ", "H-3m1      * ", &
         "H-3c1      * ", "H-31m      * ", "H-31c      * ", "P2          s", &
         "P2_1        s", "C2          s", "Pm          s", "PC          s", &
         "Cm          s", "Cc          s", "P2/m        s", "P2_1/m      s", &
         "C2/m        s", "P2/c        s", "P2_1/c      s", "C2/c        s", &
         "Pmm        *s", "Pmc        *s", "Pcc        *s", "Pma        *s", &
         "Pca        *s", "Pnc        *s", "Pmn        *s", "Pba        *s", &
         "Pna        *s", "Pnn        *s", "Cmm        *s", "Cmc        *s", &
         "Ccc        *s", "Amm        *s", "Abm        *s", "Ama        *s", &
         "Aba        *s", "Fm2m         ", "Fdd        *s", "Im2m         ", &
         "Iba        *s", "Ima        *s", "Abm2       * ", "Aba2       * ", &
         "Cmca       * ", "Cmma       * ", "Ccca       * ", "P42        *s", &
         "P42_1      *s", "P4_12      *s", "P4_12_1    *s", "P4_22      *s", &
         "P4_22_1    *s", "P4_32      *s", "P4_32_1    *s", "I42        *s", &
         "I4_12      *s", "P4cm       *s", "P4nm       *s", "P4nc       *s", &
         "P4_bc      *s", "I4md       *s", "I4cd       *s", "C-42m      *s", &
         "C-42c      *s", "C-42b      *s", "C-42n      *s", "F-42m      *s", &
         "F-42c      *s", "C3         * ", "C3_1       * ", "C3_2       * ", &
         "C-3        * ", "H32        *s", "C32        *s", "H3_12      *s", &
         "C3_12      *s", "C3_121     * ", "H3_22      *s", "C3_22      * ", &
         "C3_221     *s", "C3m        *s", "C3m1       * ", "H3m        *s", &
         "C3c        *s", "C3c1       * ", "H3c        *s", "H-3m       *s", &
         "H-3c       *s", "C-3m       *s", "C-3m1      * ", "C-3c       *s", &
         "C-3c1      * ", "C6         * ", "C6_1       * ", "C6_5       * ", &
         "C6_2,      * ", "C6_4       * ", "C6_3       * ", "C-6        * ", &
         "C6/m       * ", "C6_3/m     * ", "C62        *s", "C6_12      *s", &
         "C6_52      *s", "C6_22      *s", "C6_42      *s", "C6_32      *s", &
         "C622       * ", "C6_122     * ", "C6_522     * ", "C6_222     * ", &
         "C6_422     * ", "C6_322     * ", "C6mm       * ", "C6cc       * ", &
         "C6cm       * ", "C6mc       * ", "C-6m2      * ", "C-6c2      * ", &
         "H-6m2      * ", "H-6c2      * ", "C6/mmm     * ", "C6/mcc     * ", &
         "C6/mcm     * ", "C6/mmc     * ", "Pm3        * ", "Pn3        * ", &
         "Fm3        * ", "Fd3        * ", "Im3        * ", "Pa3        * ", &
         "Ia3        * ", "P43        *s", "P4_23      *s", "F43        *s", &
         "F4_13      *s", "I43        *s", "P4_33      *s", "P4_13      *s", &
         "I4_13      *s", "Pm3m       *s", "Pn3n       * ", "Pm3n       * ", &
         "Pn3m       *s", "Fm3m       * ", "Fm3c       * ", "Fd3m       * ", &
         "Fd3c       *s", "Im3m       * ", "Ia3d       * ", "Bma2       * ", &
         "Cm2a       * ", "Ac2m       * ", "B2cm       * ", "C2mb       * ", &
         "Acm2       * ", "Bmc2       * ", "Cm2b       * ", "Ab2m       * ", &
         "B2am       * ", "C2ma       * ", "Bba2       * ", "Cc2a       * ", &
         "Ab2a       * ", "B2ab       * ", "C2ca       * ", "Aca2       * ", &
         "Bbc2       * ", "Cc2b       * ", "Ac2a       * ", "B2cb       * ", &
         "C2cb       * ", "F2mm         ", "Fmm        *s", "I2mm         ", &
         "Imm        *s", "Ccmb       * ", "Bbcm       * ", "Bmab       * ", &
         "Abma       * ", "Acam       * ", "Cmcb       * ", "Ccma       * ", &
         "Bbam       * ", "Bmcb       * ", "Acma       * ", "Abam       * ", &
         "Bmcm       * ", "Abmm       * ", "Cmmb       * ", "Bmam       * ", &
         "Acmm       * ", "Bbcb       * ", "Abaa       * ", "Cccb       * ", &
         "Bbab       * ", "Acaa       * ", "Icab         ", "Imam         ", &
         "Ibmm         ", "Icmm         ", "P4/mmc     *s", "P4/mcm     *s", &
         "P4/nbc     *s", "P4/nnm     *s", "P4/mbc     *s", "P4/mnm     *s", &
         "P4/nmc     *s", "P4/ncm     *s", "I4/amd     *s", "I4/acd     *s", &
         "C_1^1       S", "C_i^1       S", "C_2^1       S", "C_2^2       S", &
         "C_2^3       S", "C_s^1       S", "C_s^2       S", "C_s^3       S", &
         "C_s^4       S", "C_2h^1      S", "C_2h^2      S", "C_2h^3      S", &
         "C_2h^4      S", "C_2h^5      S", "C_2h^6      S", "D_2^1       S", &
         "D_2^2       S", "D_2^3       S", "D_2^4       S", "D_2^5       S", &
         "D_2^6       S", "D_2^7       S", "D_2^8       S", "D_2^9       S", &
         "C_2v^1      S", "C_2v^2      S", "C_2v^3      S", "C_2v^4      S", &
         "C_2v^5      S", "C_2v^6      S", "C_2v^7      S", "C_2v^8      S", &
         "C_2v^9      S", "C_2v^10     S", "C_2v^11     S", "C_2v^12     S", &
         "C_2v^13     S", "C_2v^14     S", "C_2v^15     S", "C_2v^16     S", &
         "C_2v^17     S", "C_2v^18     S", "C_2v^19     S", "C_2v^20     S", &
         "C_2v^21     S", "C_2v^22     S", "D_2h^1      S", "D_2h^2      S", &
         "D_2h^3      S", "D_2h^4      S", "D_2h^5      S", "D_2h^6      S", &
         "D_2h^7      S", "D_2h^8      S", "D_2h^9      S", "D_2h^10     S", &
         "D_2h^11     S", "D_2h^12     S", "D_2h^13     S", "D_2h^14     S", &
         "D_2h^15     S", "D_2h^16     S", "D_2h^17     S", "D_2h^18     S", &
         "D_2h^19     S", "D_2h^20     S", "D_2h^21     S", "D_2h^22     S", &
         "D_2h^23     S", "D_2h^24     S", "D_2h^25     S", "D_2h^26     S", &
         "D_2h^27     S", "D_2h^28     S", "C_4^1       S", "C_4^2       S", &
         "C_4^3       S", "C_4^4       S", "C_4^5       S", "C_4^6       S", &
         "S_4^1       S", "S_4^2       S", "C_4h^1      S", "C_4h^2      S", &
         "C_4h^3      S", "C_4h^4      S", "C_4h^5      S", "C_4h^6      S", &
         "D_4^1       S", "D_4^2       S", "D_4^3       S", "D_4^4       S", &
         "D_4^5       S", "D_4^6       S", "D_4^7       S", "D_4^8       S", &
         "D_4^9       S", "D_4^10      S", "C_4v^1      S", "C_4v^2      S", &
         "C_4v^3      S", "C_4v^4      S", "C_4v^5      S", "C_4v^6      S", &
         "C_4v^7      S", "C_4v^8      S", "C_4v^9      S", "C_4v^10     S", &
         "C_4v^11     S", "C_4v^12     S", "D_2d^1      S", "D_2d^2      S", &
         "D_2d^3      S", "D_2d^4      S", "D_2d^5      S", "D_2d^6      S", &
         "D_2d^7      S", "D_2d^8      S", "D_2d^9      S", "D_2d^10     S", &
         "D_2d^11     S", "D_2d^12     S", "D_4h^1      S", "D_4h^2      S", &
         "D_4h^3      S", "D_4h^4      S", "D_4h^5      S", "D_4h^6      S", &
         "D_4h^7      S", "D_4h^8      S", "D_4h^9      S", "D_4h^10     S", &
         "D_4h^11     S", "D_4h^12     S", "D_4h^13     S", "D_4h^14     S", &
         "D_4h^15     S", "D_4h^16     S", "D_4h^17     S", "D_4h^18     S", &
         "D_4h^19     S", "D_4h^20     S", "C_3^1       S", "C_3^2       S", &
         "C_3^3       S", "C_3^4       S", "C_3i^1      S", "C_3i^2      S", &
         "D_3^1       S", "D_3^2       S", "D_3^3       S", "D_3^4       S", &
         "D_3^5       S", "D_3^6       S", "D_3^7       S", "C_3v^1      S", &
         "C_3v^2      S", "C_3v^3      S", "C_3v^4      S", "C_3v^5      S", &
         "C_3v^6      S", "D_3d^1      S", "D_3d^2      S", "D_3d^3      S", &
         "D_3d^4      S", "D_3d^5      S", "D_3d^6      S", "C_6^1       S", &
         "C_6^2       S", "C_6^3       S", "C_6^4       S", "C_6^5       S", &
         "C_6^6       S", "C_3h^1      S", "C_6h^1      S", "C_6h^2      S", &
         "D_6^1       S", "D_6^2       S", "D_6^3       S", "D_6^4       S", &
         "D_6^5       S", "D_6^6       S", "C_6v^1      S", "C_6v^2      S", &
         "C_6v^3      S", "C_6v^4      S", "D_3h^1      S", "D_3h^2      S", &
         "D_3h^3      S", "D_3h^4      S", "D_6h^1      S", "D_6h^2      S", &
         "D_6h^3      S", "D_6h^4      S", "T^1         S", "T^2         S", &
         "T^3         S", "T^4         S", "T^5         S", "T_h^1       S", &
         "T_h^2       S", "T_h^3       S", "T_h^4       S", "T_h^5       S", &
         "T_h^6       S", "T_h^7       S", "O^1         S", "O^2         S", &
         "O^3         S", "O^4         S", "O^5         S", "O^6         S", &
         "O^7         S", "O^8         S", "T_d^1       S", "T_d^2       S", &
         "T_d^3       S", "T_d^4       S", "T_d^5       S", "T_d^6       S", &
         "O_h^1       S", "O_h^2       S", "O_h^3       S", "O_h^4       S", &
         "O_h^5       S", "O_h^6       S", "O_h^7       S", "O_h^8       S", &
         "O_h^9       S", "O_h^10      S", "S_2^1       S", "V^1         S", &
         "V^2         S", "V^3         S", "V^4         S", "V^5         S", & 
         "V^6         S", "V^7         S", "V^8         S", "V^9         S", &
         "V_h^1       S", "V_h^2       S", "V_h^3       S", "V_h^4       S", &
         "V_h^5       S", "V_h^6       S", "V_h^7       S", "V_h^8       S", &
         "V_h^9       S", "V_h^10      S", "V_h^11      S", "V_h^12      S", &
         "V_h^13      S", "V_h^14      S", "V_h^15      S", "V_h^16      S", &
         "V_h^17      S", "V_h^18      S", "V_h^19      S", "V_h^20      S", &
         "V_h^21      S", "V_h^22      S", "V_h^23      S", "V_h^24      S", &
         "V_h^25      S", "V_h^26      S", "V_h^27      S", "V_h^28      S", &
         "V_d^1       S", "V_d^2       S", "V_d^3       S", "V_d^4       S", &
         "V_d^5       S", "V_d^6       S", "V_d^7       S", "V_d^8       S", &
         "V_d^9       S", "V_d^10      S", "V_d^11      S", "V_d^12      S", &
         "S_6^1        ", "S_6^2        ", "A12_11       ", "A112_1       ", &
         "B2_111       ", "B112_1       ", "C2_111       ", "C12_11       ", &
         "I2_111       ", "I12_11       ", "I112_1       ", "I12/c1       ", &
         "I112/b       ", "I2/c11       ", "A2_12_12_1   ", "B2_12_12_1   ", &
         "C2_12_12_1   ", "F22_12_1     ", "F2_122_1     ", "F2_12_12     ", &
         "I4_2         ", "F4_2         ", "I4_3         ", "F4_3         "/


  DATA space_group_numbers &
      / 1,               2,             3,             4,            &
        5,               6,             7,             8,            &
        9,              10,            11,            12,            &
       13,              14,            15,             3,            &
        3,               3,             3,             3,            &
        4,               4,             4,             4,            &
        4,               5,             5,             5,            &
        5,               5,             5,             5,            &
        5,               6,             6,             6,            &
        6,               6,             7,             7,            &
        7,               7,             7,             7,            &
        7,               7,             8,             8,            &
        8,               8,             8,             8,            &
        8,               8,             9,             9,            &
        9,               9,             9,             9,            &
        9,               9,            10,            10,            &
       10,              10,            10,            11,            &
       11,              11,            11,            11,            &
       12,              12,            12,            12,            &
       12,              12,            12,            12,            &
       13,              13,            13,            13,            &
       13,              13,            13,            13,            &
       14,              14,             0,             0,            &
       14,              14,            14,            14,            &
       14,              14,            15,            15,            &
       15,              15,            15,            15,            &
       15,              15,            16,            17,            &
       18,              19,            20,            21,            &
       22,              23,            24,            25,            &
       26,              27,            28,            29,            &
       30,              31,            32,            33,            &
       34,              35,            36,            37,            &
       38,              39,            40,            41,            &
       42,              43,            44,            45,            &
       46,              47,            48,            49,            &
       50,              51,            52,            53,            &
       54,              55,            56,            57,            &
       58,              59,            60,            61,            &
       62,              63,            64,            65,            &
       66,              67,            68,            69,            &
       70,              71,            72,            73,            &
       74,              75,            76,            77,            &
       78,              79,            80,            81,            &
       82,              83,            84,            85,            &
       86,              87,            88,            89,            &
       90,              91,            92,            93,            &
       94,              95,            96,            97,            &
       98,              99,           100,           101,            &
      102,             103,           104,           105,            &
      106,             107,           108,           109,            &
      110,             111,           112,           113,            &
      114,             115,           116,           117,            &
      118,             119,           120,           121,            &
      122,             123,           124,           125,            &
      126,             127,           128,           129,            &
      130,             131,           132,           133,            &
      134,             135,           136,           137,            &
      138,             139,           140,           141,            &
      142,             143,           144,           145,            &
      146,             147,           148,           149,            &
      150,             151,           152,           153,            &
      154,             155,           156,           157,            &
      158,             159,           160,           161,            &
      162,             163,           164,           165,            &
      166,             167,           168,           169,            &
      170,             171,           172,           173,            &
      174,             175,           176,           177,            &
      178,             179,           180,           181,            &
      182,             183,           184,           185,            &
      186,             187,           188,           189,            &
      190,             191,           192,           193,            &
      194,             195,           196,           197,            &
      198,             199,           200,           201,            &
      202,             203,           204,           205,            &
      206,             207,           208,           209,            &
      210,             211,           212,           213,            &
      214,             215,           216,           217,            &
      218,             219,           220,           221,            &
      222,             223,           224,           225,            &
      226,             227,           228,           229,            &
      230,              17,            17,            18,            &
       18,              20,            20,            21,            &
       21,              25,            25,            26,            &
       26,              26,            26,            26,            &
       27,              27,            28,            28,            &
       28,              28,            28,            29,            &
       29,              29,            29,            29,            & 
       30,              30,            30,            30,            &
       30,              31,            31,            31,            &
       31,              31,            32,            32,            &
       33,              33,            33,            33,            &
       33,              34,            34,            35,            &
       35,              36,            36,            36,            &
       36,              36,            37,            37,            &
       38,              38,            38,            38,            &
       38,              39,            39,            39,            &
       39,              39,            40,            40,            &
       40,              40,            40,            41,            &
       41,              41,            41,            41,            &
       43,              43,            45,            45,            &
       46,              46,            46,            46,            &
       46,              49,            49,            50,            &
        0,              50,             0,            51,            &
       51,              51,            51,            51,            &
       52,              52,            52,            52,            &
       52,              53,            53,            53,            &
       53,              53,            54,            54,            &
       54,              54,            54,            55,            &
       55,              56,            56,            57,            &
       57,              57,            57,            57,            &
       58,              58,            59,            59,            &
       60,              60,            60,            60,            &
       60,              61,            62,            62,            &
       62,              62,            62,            63,            &
       63,              63,            63,            63,            &
       64,              64,            64,            64,            &
       64,              66,            66,            67,            &
       67,              68,            68,            72,            &
       72,              74,            74,             1,            &
        1,               1,             1,             1,            &
        2,               2,             2,             2,            &
        2,              75,            76,            77,            &
       78,              79,            80,            81,            &
       82,              83,            84,            85,            &
       86,              87,            88,            89,            &
       90,              91,            92,            93,            &
       94,              95,            96,            97,            &
       98,              99,           100,           101,            &
      102,             103,           104,           105,            &
      106,             107,           108,           109,            &
      110,             111,           112,           113,            &
      114,             115,           116,           117,            &
      118,             119,           120,           121,            &
      122,             123,           124,           125,            &
      126,             127,           128,           129,            &
      130,             131,           132,           133,            &
      134,             135,           136,           137,            &
      138,             139,           140,           141,            &
      142,             143,           144,           145,            &
      147,             149,           150,           151,            &
      152,             153,           154,           156,            &
      157,             158,           159,           162,            &
      163,             164,           165,             3,            & 
        4,               5,             6,             7,            &
        8,               9,            10,            11,            &
       12,              13,            14,            15,            &
       25,              26,            27,            28,            &
       29,              30,            31,            32,            &
       33,              34,            35,            36,            &
       37,              38,            39,            40,            &
       41,              42,            43,            44,            &
       45,              46,            39,            41,            &
       64,              67,            68,            89,            &
       90,              91,            92,            93,            &
       94,              95,            96,            97,            &
       98,             101,           102,           105,            &
      106,             109,           110,           115,            &
      116,             117,           118,           119,            &
      120,             143,           144,           145,            &
      147,             149,           150,           151,            &
      152,             152,           153,           153,            &
      154,             156,           156,           157,            &
      158,             158,           159,           162,            &
      163,             164,           164,           165,            &
      165,             168,           169,           170,            &
      171,             172,           173,           174,            &
      175,             176,           177,           178,            & 
      179,             180,           181,           182,            &
      177,             178,           179,           180,            &
      181,             182,           183,           184,            &
      185,             186,           187,           188,            &
      189,             190,           191,           192,            &
      193,             194,           200,           201,            &  
      202,             203,           204,           205,            &
      206,             207,           208,           209,            &
      210,             211,           212,           213,            &
      214,             221,           222,           223,            &
      224,             225,           226,           227,            &
      228,             229,           230,            39,            &
       39,              39,            39,            39,            &
       39,              39,            39,            39,            &
       39,              39,            41,            41,            &
       41,              41,            41,            41,            &          
       41,              41,            41,            41,            &         
       41,              42,            42,            44,            &
       44,              64,            64,            64,            &
       64,              64,            64,            64,            &          
       64,              64,            64,            64,            &
       67,              67,            67,            67,            &
       67,              68,            68,            68,            &
       68,              68,            73,            74,            &
       74,              74,           131,           132,            &
      133,             134,           135,           136,            &
      137,             138,           141,           142,            &
        1,               2,             3,             4,            &
        5,               6,             7,             8,            &
        9,              10,            11,            12,            &
       13,              14,            15,            16,            &
       17,              18,            19,            20,            &
       21,              22,            23,            24,            &
       25,              26,            27,            28,            &
       29,              30,            31,            32,            &
       33,              34,            35,            36,            &
       37,              38,            39,            40,            &
       41,              42,            43,            44,            &
       45,              46,            47,            48,            &
       49,              50,            51,            52,            &
       53,              54,            55,            56,            &
       57,              58,            59,            60,            &
       61,              62,            63,            64,            &
       65,              66,            67,            68,            &
       69,              70,            71,            72,            &
       73,              74,            75,            76,            &
       77,              78,            79,            80,            &
       81,              82,            83,            84,            &
       85,              86,            87,            88,            &
       89,              90,            91,            92,            &
       93,              94,            95,            96,            &
       97,              98,            99,           100,            &
      101,             102,           103,           104,            &
      105,             106,           107,           108,            &
      109,             110,           111,           112,            &
      113,             114,           115,           116,            &
      117,             118,           119,           120,            &
      121,             122,           123,           124,            &
      125,             126,           127,           128,            &
      129,             130,           131,           132,            &
      133,             134,           135,           136,            &
      137,             138,           139,           140,            &
      141,             142,           143,           144,            &
      145,             146,           147,           148,            &
      149,             150,           151,           152,            & 
      153,             154,           155,           156,            &
      157,             158,           159,           160,            &
      161,             162,           163,           164,            &
      165,             166,           167,           168,            &
      169,             170,           171,           172,            &
      173,             174,           175,           176,            &
      177,             178,           179,           180,            &
      181,             182,           183,           184,            &
      185,             186,           187,           188,            &
      189,             190,           191,           192,            &
      193,             194,           195,           196,            &
      197,             198,           199,           200,            &
      201,             202,           203,           204,            &
      205,             206,           207,           208,            & 
      209,             210,           211,           212,            &
      213,             214,           215,           216,            &
      217,             218,           219,           220,            &
      221,             222,           223,           224,            &
      225,             226,           227,           228,            &
      229,             230,             2,            16,            &
       17,              18,            19,            20,            &
       21,              22,            23,            24,            &
       47,              48,            49,            50,            &
       51,              52,            53,            54,            &
       55,              56,            57,            58,            &
       59,              60,            61,            62,            &
       63,              64,            65,            66,            &
       67,              68,            69,            70,            &
       71,              72,            73,            74,            &
      111,             112,           113,           114,            &
      115,             116,           117,           118,            &
      119,             120,           121,           122,            &
      147,             148,             5,             5,            &
        5,               5,             5,             5,            &
        5,               5,             5,            15,            &
       15,              15,            20,            20,            &
       20,              22,            22,            22,            &
       79,              79,            80,            80   /


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

 
  PUBLIC   spg_code, spg_name, set_spg_code, find_space_group, sg_name, &
           equivalent_tau, space_group_names, &
           find_space_group_number, find_space_group_names, add_info, &
           set_add_info, set_fft_fact, project_frac_tran,  &
           shift_frac_tran, set_standard_sg, set_point_group_code

  CONTAINS

  SUBROUTINE set_add_info()
  IMPLICIT NONE
!
! 17  P222_1
!
  add_info(108)='xyz   xyz'
  add_info(322)='zxy   yzx'
  add_info(323)='yzx   zxy'
!
! 18  P2_12_12
!
  add_info(109)='xyz   xyz'
  add_info(324)='yzx   zxy'
  add_info(325)='zxy   yzx'
!
! 19  P2_12_12
!
  add_info(109)='xyz   xyz'
  add_info(324)='yzx   zxy'
  add_info(325)='zxy   yzx'
!
! 20 C222_1
!
  add_info(111)='xyz   xyz'
  add_info(326)='yzx   zxy'
  add_info(327)='zxy   yzx'
!
! 21 C222
!
  add_info(112)='xyz   xyz'
  add_info(328)='yzx   zxy'
  add_info(329)='zxy   yzx'
!
! 25 Pmm2
!
  add_info(116)='xyz   xyz'
  add_info(330)='yzx   zxy'
  add_info(331)='zxy   yzx'
!
! 26 Pmc2_1
!
  add_info(117)='xyz   xyz'
  add_info(332)='yx-z  yx-z'
  add_info(333)='yzx   zxy'
  add_info(334)='x-zy  xz-y'
  add_info(335)='zxy   yzx'
  add_info(336)='-zyx  zy-x'
!
! 27 Pcc2
!
  add_info(118)='xyz   xyz'
  add_info(337)='yzx   zxy'
  add_info(338)='zxy   yzx'
!
! 28 Pma2
!
  add_info(119)='xyz   xyz'
  add_info(339)='yx-z  yx-z'
  add_info(340)='yzx   zxy'
  add_info(341)='x-zy  xz-y'
  add_info(342)='zxy   yzx'
  add_info(343)='-zyx  zy-x'
!
! 29 Pca2_1
!
  add_info(120)='xyz   xyz'
  add_info(344)='yx-z  yx-z'
  add_info(345)='yzx   zxy'
  add_info(346)='x-zy  xz-y'
  add_info(347)='zxy   yzx'
  add_info(348)='-zyx  zy-x'
!
! 30 Pnc2
!
  add_info(121)='xyz   xyz'
  add_info(349)='yx-z  yx-z'
  add_info(350)='yzx   zxy'
  add_info(351)='x-zy  xz-y'
  add_info(352)='zxy   yzx'
  add_info(353)='-zyx  zy-x'
!
! 31 Pmn2_1
!
  add_info(122)='xyz   xyz'
  add_info(354)='yx-z  yx-z'
  add_info(355)='yzx   zxy'
  add_info(356)='x-zy  xz-y'
  add_info(357)='zxy   yzx'
  add_info(358)='-zyx  zy-x'
!
! 32 Pba2
!
  add_info(123)='xyz   xyz'
  add_info(359)='yzx   zxy'
  add_info(360)='zxy   yzx'
!
! 33 Pna2_1
!
  add_info(124)='xyz   xyz'
  add_info(361)='yx-z  yx-z'
  add_info(362)='yzx   zxy'
  add_info(363)='x-zy  xz-y'
  add_info(364)='zxy   yzx'
  add_info(365)='-zyx  zy-x'
!
! 34 Pna2_1
!
  add_info(125)='xyz   xyz'
  add_info(366)='yzx   zxy'
  add_info(367)='zxy   yzx'
!
! 35 Cmm2
!
  add_info(126)='xyz   xyz'
  add_info(368)='yzx   zxy'
  add_info(369)='zxy   yzx'
!
! 36 Cmc2_1
!
  add_info(127)='xyz   xyz'
  add_info(370)='yx-z  yx-z'
  add_info(371)='yzx   zxy'
  add_info(372)='x-zy  xz-y'
  add_info(373)='zxy   yzx'
  add_info(374)='-zyx  zy-x'
!
! 37 Ccc2
!
  add_info(128)='xyz   xyz'
  add_info(375)='yzx   zxy'
  add_info(376)='zxy   yzx'
!
! 38 Amm2
!
  add_info(129)='xyz   xyz'
  add_info(377)='yx-z  yx-z'
  add_info(378)='yzx   zxy'
  add_info(379)='x-zy  xz-y'
  add_info(380)='zxy   yzx'
  add_info(381)='-zyx  zy-x'
!
! 39 Aem2
!
  add_info(130)='xyz   xyz'
  add_info(382)='yx-z  yx-z'
  add_info(383)='yzx   zxy'
  add_info(384)='x-zy  xz-y'
  add_info(385)='zxy   yzx'
  add_info(386)='-zyx  zy-x'
  add_info(608)='xyz   xyz'
  add_info(716)='yx-z  yx-z'
  add_info(717)='yzx   zxy'
  add_info(718)='x-zy  xz-y'
  add_info(719)='zxy   yzx'
  add_info(720)='-zyx  zy-x'
  add_info(721)='xyz   xyz'
  add_info(722)='yx-z  yx-z'
  add_info(723)='yzx   zxy'
  add_info(724)='x-zy  xz-y'
  add_info(725)='zxy   yzx'
  add_info(726)='-zyx  zy-x'
!
! 40 Ama2
!
  add_info(131)='xyz   xyz'
  add_info(387)='yx-z  yx-z'
  add_info(388)='yzx   zxy'
  add_info(389)='x-zy  xz-y'
  add_info(390)='zxy   yzx'
  add_info(391)='-zyx  zy-x'
!
! 41 Aea2
!
  add_info(132)='xyz   xyz'
  add_info(392)='yx-z  yx-z'
  add_info(393)='yzx   zxy'
  add_info(394)='x-zy  xz-y'
  add_info(395)='zxy   yzx'
  add_info(396)='-zyx  zy-x'
  add_info(608)='xyz   xyz'
  add_info(727)='yx-z  yx-z'
  add_info(728)='yzx   zxy'
  add_info(729)='x-zy  xz-y'
  add_info(730)='zxy   yzx'
  add_info(731)='-zyx  zy-x'
  add_info(732)='xyz   xyz'
  add_info(733)='yx-z  yx-z'
  add_info(734)='yzx   zxy'
  add_info(735)='x-zy  xz-y'
  add_info(736)='zxy   yzx'
  add_info(737)='-zyx  zy-x'
!
! 42 Fmm2
!
  add_info(133)='xyz   xyz'
  add_info(602)='yzx   zxy'
  add_info(738)='zxy   yzx'
!
! 43 Fdd2
!
  add_info(134)='xyz   xyz'
  add_info(397)='yzx   zxy'
  add_info(398)='zxy   yzx'
!
! 44 Imm2
!
  add_info(135)='xyz   xyz'
  add_info(604)='yzx   zxy'
  add_info(740)='zxy   yzx'
!
! 45 Iba2
!
  add_info(136)='xyz   xyz'
  add_info(399)='yzx   zxy'
  add_info(400)='zxy   yzx'
!
! 46 Ima2
!
  add_info(137)='xyz   xyz'
  add_info(401)='yx-z  yx-z'
  add_info(402)='yzx   zxy'
  add_info(403)='x-zy  xz-y'
  add_info(404)='zxy   yzx'
  add_info(405)='-zyx  zy-x'
!
! 49 Pccm
!
  add_info(140)='xyz   xyz'
  add_info(406)='yzx   zxy'
  add_info(407)='zxy   yzx'
!
! 50 Pban
!
  add_info(141)='xyz   xyz'
  add_info(408)='yzx   zxy'
  add_info(410)='zxy   yzx'
!
! 51 Pmma
!
  add_info(142)='xyz   xyz'
  add_info(412)='yx-z  yx-z'
  add_info(413)='yzx   zxy'
  add_info(414)='x-zy  xz-y'
  add_info(415)='zxy   yzx'
  add_info(416)='-zyx  zy-x'
!
! 52 Pmma
!
  add_info(143)='xyz   xyz'
  add_info(417)='yx-z  yx-z'
  add_info(418)='yzx   zxy'
  add_info(419)='x-zy  xz-y'
  add_info(420)='zxy   yzx'
  add_info(421)='-zyx  zy-x'
!
! 53 Pmna
!
  add_info(144)='xyz   xyz'
  add_info(422)='yx-z  yx-z'
  add_info(423)='yzx   zxy'
  add_info(424)='x-zy  xz-y'
  add_info(425)='zxy   yzx'
  add_info(426)='-zyx  zy-x'
!
! 54 Pcca
!
  add_info(145)='xyz   xyz'
  add_info(427)='yx-z  yx-z'
  add_info(428)='yzx   zxy'
  add_info(439)='x-zy  xz-y'
  add_info(430)='zxy   yzx'
  add_info(431)='-zyx  zy-x'
!
! 55 Pbam
!
  add_info(146)='xyz   xyz'
  add_info(432)='yzx   zxy'
  add_info(433)='zxy   yzx'
!
! 56 Pccn
!
  add_info(147)='xyz   xyz'
  add_info(434)='yzx   zxy'
  add_info(435)='zxy   yzx'
!
! 57 Pbcm
!
  add_info(148)='xyz   xyz'
  add_info(436)='yx-z  yx-z'
  add_info(437)='yzx   zxy'
  add_info(438)='x-zy  xz-y'
  add_info(439)='zxy   yzx'
  add_info(440)='-zyx  zy-x'
!
! 58 Pnnm
!
  add_info(149)='xyz   xyz'
  add_info(441)='yzx   zxy'
  add_info(442)='zxy   yzx'
!
! 59 Pmnm
!
  add_info(150)='xyz   xyz'
  add_info(443)='yzx   zxy'
  add_info(444)='zxy   yzx'
!
! 60 Pmnm
!
  add_info(151)='xyz   xyz'
  add_info(445)='yx-z  yx-z'
  add_info(446)='yzx   zxy'
  add_info(447)='x-zy  xz-y'
  add_info(448)='zxy   yzx'
  add_info(449)='-zyx  zy-x'
!
! 61 Pbca
!
  add_info(152)='xyz   xyz'
  add_info(450)='yx-z  yx-z'
!
! 62 Pnma
!
  add_info(153)='xyz   xyz'
  add_info(451)='yx-z  yx-z'
  add_info(452)='yzx   zxy'
  add_info(453)='x-zy  xz-y'
  add_info(454)='zxy   yzx'
  add_info(455)='-zyx  zy-x'
!
! 63 Cmcm
!
  add_info(154)='xyz   xyz'
  add_info(456)='yx-z  yx-z'
  add_info(457)='yzx   zxy'
  add_info(458)='x-zy  xz-y'
  add_info(459)='zxy   yzx'
  add_info(460)='-zyx  zy-x'
!
!  64 Cmce
!
  add_info(155)='xyz   xyz'
  add_info(461)='yx-z  yx-z'
  add_info(462)='yzx   zxy'
  add_info(463)='x-zy  xz-y'
  add_info(464)='zxy   yzx'
  add_info(465)='-zyx  zy-x'
  add_info(609)='xyz   xyz'
  add_info(742)='yx-z  yx-z'
  add_info(743)='yzx   zxy'
  add_info(744)='x-zy  xz-y'
  add_info(745)='zxy   yzx'
  add_info(746)='-zyx  zy-x'
  add_info(747)='xyz   xyz'
  add_info(748)='yx-z  yx-z'
  add_info(749)='yzx   zxy'
  add_info(750)='x-zy  xz-y'
  add_info(751)='zxy   yzx'
  add_info(752)='-zyx  zy-x'
!
!  66 Cccm
!
  add_info(157)='xyz   xyz'
  add_info(466)='yzx   zxy'
  add_info(467)='zxy   yzx'
!
!  67 Cmme
!
  add_info(158)='xyz   xyz'
  add_info(468)='yzx   zxy'
  add_info(469)='zxy   yzx'
  add_info(610)='xyz   xyz'
  add_info(753)='yzx   zxy'
  add_info(754)='zxy   xyz'
  add_info(755)='xyz   xyz'
  add_info(756)='yzx   zxy'
  add_info(757)='zxy   yzx'
!
!  68 Ccce
!
  add_info(159)='xyz   xyz'
  add_info(470)='yzx   zxy'
  add_info(471)='zxy   yzx'
  add_info(611)='xyz   xyz'
  add_info(758)='yzx   zxy'
  add_info(759)='zxy   yzx'
  add_info(760)='xyz   xyz'
  add_info(761)='yzx   zxy'
  add_info(762)='zxy   yzx'
!
! 72 Ibam
!
  add_info(163)='xyz   xyz'
  add_info(472)='yzx   zxy'
  add_info(473)='zxy   yzx'
!
! 73 Ibca
!
  add_info(164)='xyz   xyz'
  add_info(763)='yx-z  yx-z'
!
! 74 Imma
!
  add_info(165)='xyz   xyz'
  add_info(474)='yx-z  yx-z'
  add_info(475)='yzx   zxy'
  add_info(764)='x-zy  xz-y'
  add_info(765)='zxy   yzx'
  add_info(766)='-zyx  zy-x'

  RETURN
  END SUBROUTINE set_add_info

  SUBROUTINE set_fft_fact(code, unique, trig, fft_fact)
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

  SUBROUTINE find_space_group(ibrav, nsym, sr_in, ft_in, at, bg, &
                                     sg_number, s01, s02, verbosity)
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
  USE io_global, ONLY : stdout
  USE point_group, ONLY : find_group_info_ext, sym_label
  USE lattices,    ONLY : is_centered, compute_conventional, lattice_name
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: ibrav, nsym
  REAL(DP), INTENT(IN)  :: sr_in(3,3,nsym), ft_in(3,nsym), at(3,3), bg(3,3)
  LOGICAL, INTENT(IN)   :: verbosity
  INTEGER, INTENT(OUT)  :: sg_number
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
             nsa, nmir, ncomp, nmp, ipa ! auxiliary

  REAL(DP) ::  &
              sr(3,3,48),                 &! symmetry matrices
              ft(3,48), ftc(3,48),        &! fractional translation crys., cart.
              ftpar(3,48), ftparc(3,48),  &! fract. transl. parallel symmetry
              ftperp(3,48), ftperpc(3,48),&! fract. transl. perp. symmetry
              s0(3,48), s0c(3,48),        &! origin of each symmetry
              ax(3), bx(3), angle, fact,  &! auxiliary
              ftrs1(3), ftrs2(3), ftrs3(3),& ! auxiliary
              atp(3,3), bgp(3,3),          & ! conventional bravais lattice
              sp, sp1, sp2, nftd             ! auxiliary

  REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP, eps1=1.D-8
  LOGICAL :: is_axis, ft_zero
  INTEGER :: tipo_sym, ts
  CHARACTER(LEN=11) :: group_name
  CHARACTER(LEN=40) :: latt_name
  REAL(DP) :: angle_rot
!
!  initialization
!  
  s01=1000.0_DP
  s02=1000.0_DP
  IF (ibrav==0) THEN
     sg_number=0 
     RETURN
  ENDIF

!
!   Find the point group
!
  CALL find_group_info_ext(nsym, sr_in, code_group, code_group_ext, &
                                                  which_elem, group_desc) 

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
                                        SGS = space group symmetry ")')
        WRITE(stdout,'(5x,"(no shift) means the part of the fractional &
                        translation")')
        WRITE(stdout,'(5x,"that cannot be removed by an origin shift")')

        WRITE(stdout,'(/,5x,"Fractional translation shift and origin &
                       &shift in crystal coordinates")')
        WRITE(stdout,'(/,6x, "PGS",6x,"FT shift &
          &                   Origin shift",17x,"SGS",/)')
        DO isym=1,nsym
           WRITE(stdout,'(1x,a8,i4,3f8.4,4x,3f8.4,a10,i4)') &
              TRIM(sym_label(group_desc(isym))), group_desc(isym), &
              ftperp(1,isym), ftperp(2,isym), ftperp(3,isym), &
              s0(1,isym), s0(2,isym), s0(3,isym), &
              TRIM(sym_label_sg(group_desc_sg(isym))), group_desc_sg(isym) 
        END DO

        IF (is_centered(ibrav)) THEN
           WRITE(stdout,'(/,5x,"Fractional translation shift and origin shift &
                         &in conventional crystal coordinates")')
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
          IF (ABS(ibrav)==13) THEN
             IF (is_symmorphic) THEN
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
            IF (is_symmorphic) THEN
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
               IF (nft==1) sg_number=17
               IF (nft==2) sg_number=18
               IF (nft==3) THEN
                  sg_number=19
!
!   ITA origin not programmed
!
                  s01=1000.0_DP
               ENDIF
            ENDIF
         ELSEIF (ibrav==9) THEN
            IF (is_symmorphic) THEN
               sg_number=21
            ELSE
               sg_number=20
            ENDIF
         ELSEIF (ibrav==10) THEN
            sg_number=22
         ELSEIF (ibrav==11) THEN
            IF (is_symmorphic) THEN
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
            IF (is_symmorphic) THEN
               sg_number=89
            ELSE 
!
!   Find the axis of order 4, and the axis of order 2 parallel to x
!
               irot90=2
               irot180=5

               IF (ABS(ft(3,irot90)) < eps1 ) THEN
!
!  axis 4
!
                  sg_number=90
                  s01(:)=s0(:,3)
                  s01(3)=s0(3,5)
               ELSEIF (group_desc_sg(irot90) == 73) THEN
!
!  axis 4_1
!
                  IF (ABS(ft(1,irot180))<eps1) THEN
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
                  IF (ABS(ft(1,irot180))<eps1) THEN
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
                  IF (ABS(ft(1,irot180)) < eps1 ) THEN
                     sg_number=95
                     s01(:)=s0(:,2)
                     s01(3)=s0(3,7)
                  ELSE
                     sg_number=96
                     s01(:)=s0(:,3)
                     s01(3)=s0(3,6)
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (is_symmorphic) THEN
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
         IF (is_symmorphic) THEN
            sg_number=177
         ELSE
            irot60=2
            IF (group_desc_sg(irot60) == 108) sg_number=178   ! 6_1
            IF (group_desc_sg(irot60) == 112) sg_number=179   ! 6_5
            IF (group_desc_sg(irot60) == 109) sg_number=180   ! 6_2
            IF (group_desc_sg(irot60) == 111) sg_number=181   ! 6_4
            IF (group_desc_sg(irot60) == 110) sg_number=182   ! 6_3
         ENDIF
      CASE(12)
!
!   C_2v
!
!   the origin in the C_2 axis
!
         s01(:)=s0(:,2)
         s01(3)=0.0_DP
!
!  identify the axis of order 2
!
         irot180=2
         imirror1=3
         imirror2=4
         ftrs1=ftc(:,3)
         ftrs2=ftc(:,4)
         IF (code_group_ext==54) THEN
            idir=1
            ftrs1(2)=0.0_DP
            ftrs2(3)=0.0_DP
         ELSEIF (code_group_ext==56) THEN
            idir=2
            ftrs1(3)=0.0_DP
            ftrs2(1)=0.0_DP
         ELSEIF (code_group_ext==58) THEN
            idir=3
            ftrs1(1)=0.0_DP
            ftrs2(2)=0.0_DP
         ENDIF 

         IF (ibrav==8) THEN
            IF (is_symmorphic) THEN
               sg_number=25
            ELSE
!
!   find if there is a component of the fractional translation perpendicular
!   to the two fold axis
!
               sp1=0
               sp2=0
               DO ipol=1,3
                  IF (ipol/=idir) THEN
                     sp1=ftrs1(ipol)**2+sp1
                     sp2=ftrs2(ipol)**2+sp2
                  ENDIF
               ENDDO
!
!  count the true mirrors
!
               nft=0
               IF (ABS(ftrs1(1))<eps1.AND.ABS(ftrs1(2))<eps1.AND.&
                                          &ABS(ftrs1(3))<eps1) nft=nft+1
               IF (ABS(ftrs2(1))<eps1.AND.ABS(ftrs2(2))<eps1.AND.&
                                          &ABS(ftrs2(3))<eps1) nft=nft+1
!
!    Now divide the groups with a proper or improper twofold axis
!
               IF (ABS(ftc(idir,2))<eps1) THEN
!
!   proper axis 
!
                  IF (nft==1) THEN
!
!  and a true mirror
!
                     sg_number=28
                  ELSE
                     IF (ABS(ftrs1(idir))<eps1.AND.ABS(ftrs2(idir))<eps1) THEN
!
!   two mirrors with glide perpendicular to the axis
!
                        sg_number=32
                     ELSEIF (ABS(sp1)<eps1.AND.ABS(sp2)<eps1) THEN
!
!   two mirrors with glide parallel to the axis
!
                        sg_number=27
                     ELSEIF (ABS(sp1)<eps1.OR.ABS(sp2)<eps1) THEN 
!
!   one mirrors with glide perpendicular to the axis
!
                        sg_number=30
                     ELSE
!
!  something else
!
                        sg_number=34
                     ENDIF
                  ENDIF
               ELSE
!
!   In this case the two fold axis is improper
!
                  IF (nft == 1) THEN
!
!    one proper mirror
!
                     IF ( ABS(sp1) > eps1 .OR. ABS(sp2) > eps1 ) THEN
!
!  one glide parallel to the axis
!
                        s01(1)=s0(1,3)
                        s01(2)=s0(2,4)
                        s01(3)=0.0_DP
                        sg_number=31
                     ELSE
                        sg_number=26
                     ENDIF
                  ELSE
!
!  two glide planes
!
                     IF ( sp1 > eps1 .AND. sp2 > eps1 ) THEN
!
!  no glide parallel to the twofold axis
!
                        sg_number=33
                     ELSE
                        sg_number=29
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==9) THEN
            IF (is_symmorphic) THEN
               sg_number=35
            ELSE 
               IF (group_desc_sg(2) == 65) THEN
                  sg_number=36
               ELSE
                  sg_number=37
               ENDIF
            ENDIF
         ELSEIF (ibrav==91) THEN
!
!   find the mirror perpendicular to y
!
            imirror1=4
            IF (is_symmorphic) THEN
               sg_number=38
            ELSE 
               IF (ABS(ft(1,imirror1))<eps1) THEN
                  sg_number=39
               ELSE
!
!  Find the mirror perpendicular to x
!

                  imirror2=3
!
!  bring the fractional translation in the orthorombic axis
!
                  ftrs1(:) = ftc(:,imirror2)
!
!    check if it is a proper mirrror or a glide plane (in the latter case only
!    one component is different from zero)
!
                  nft=0
                  IF ((ABS(ftrs1(2)) < 1.d-6 .AND. ABS(ftrs1(3)) > 1.d-6) .OR.&
                (ABS(ftrs1(3)) < 1.d-6 .AND. ABS(ftrs1(2)) > 1.d-6) ) nft=nft+1

                  IF (nft==0) THEN
!
!   the mirror perpendicular to x is a proper mirror
!
                     sg_number=40
                  ELSE
!
!    the mirror perpendicular to x is a glide plane
!
                     sg_number=41
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==10) THEN
            IF (is_symmorphic) THEN
               sg_number=42
            ELSE 
               sg_number=43
            ENDIF
         ELSEIF (ibrav==11) THEN
            IF (is_symmorphic) THEN
               sg_number=44
            ELSE 
!
!  Find the two mirrors perpendicular to the proper axis
!

               imirror1=0
               imirror2=0
               DO isym=2, nsym
                  IF (tipo_sym(sr(1,1,isym))==5) THEN
                     IF (imirror1==0) THEN
                        imirror1=isym
                        CALL mirror_axis(sr(1,1,isym),ax)
                        DO ipol=1,3 
                           IF (is_axis(ax,ipol)) idir1=ipol
                        ENDDO
                     ELSE
                        imirror2=isym
                        CALL mirror_axis(sr(1,1,isym),ax)
                        DO ipol=1,3 
                           IF (is_axis(ax,ipol)) idir2=ipol
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
               IF ( imirror1 == 0 ) CALL errore('find_space_group',&
                                         'C_2v: no mirror found',1)
               IF ( imirror2 == 0 ) CALL errore('find_space_group',&
                                         'C_2v: second mirror not found',1)

!
!  bring the fractional translation in the orthorombic axis
!
               ftrs1(:) = ftc(:,imirror1)  
               ftrs2(:) = ftc(:,imirror2)
!
!    check if mirrror1 is a glide plane 
!
               nft=0
               IF (idir1==1) THEN
                  IF ((ABS(ftrs1(2)) < 1.d-6 .AND. ABS(ftrs1(3)) > 1.d-6) .OR. &
                 (ABS(ftrs1(3)) < 1.d-6 .AND. ABS(ftrs1(2)) > 1.d-6) ) nft=nft+1
               ELSEIF(idir1==2) THEN
                  IF ((ABS(ftrs1(1)) < 1.d-6 .AND. ABS(ftrs1(3)) > 1.d-6) .OR.&
                 (ABS(ftrs1(3)) < 1.d-6 .AND. ABS(ftrs1(1)) > 1.d-6) ) nft=nft+1
               ELSEIF(idir1==3) THEN
                  IF ((ABS(ftrs1(1)) < 1.d-6 .AND. ABS(ftrs1(2)) > 1.d-6) .OR. &
                 (ABS(ftrs1(1)) < 1.d-6 .AND. ABS(ftrs1(2)) > 1.d-6) ) nft=nft+1
               ENDIF
               IF (idir2==1) THEN
                  IF ((ABS(ftrs2(2)) < 1.d-6 .AND. ABS(ftrs2(3)) > 1.d-6) .OR. &
                 (ABS(ftrs2(3)) < 1.d-6 .AND. ABS(ftrs2(2)) > 1.d-6) ) nft=nft+1
               ELSEIF(idir2==2) THEN
                  IF ((ABS(ftrs2(1)) < 1.d-6 .AND. ABS(ftrs2(3)) > 1.d-6) .OR. &
                 (ABS(ftrs2(3)) < 1.d-6 .AND. ABS(ftrs2(1)) > 1.d-6) ) nft=nft+1
               ELSEIF(idir2==3) THEN
                  IF ((ABS(ftrs2(1)) < 1.d-6 .AND. ABS(ftrs2(2)) > 1.d-6) .OR. &
                 (ABS(ftrs2(1)) < 1.d-6 .AND. ABS(ftrs2(2)) > 1.d-6) ) nft=nft+1
               ENDIF

               IF (nft==1) THEN
!
!   one mirror is proper
!
                  sg_number=46
               ELSEIF (nft==2) THEN
!
!   both mirrors are glide planes
!
                  sg_number=45
               ENDIF
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
            IF (is_symmorphic) THEN
               IF (type_group==1) THEN
                  sg_number=156
               ELSE
                  sg_number=157
               ENDIF
            ELSE 
               IF (type_group==1) THEN
                  sg_number=158
               ELSE
                  sg_number=159
               ENDIF
            ENDIF
         ELSEIF (ibrav==5) THEN
            IF (is_symmorphic) THEN
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
            IF (is_symmorphic) THEN
               sg_number=99
            ELSE 
               IF (ABS(ft(3,irot90)) < eps1) THEN
                  IF (ABS(ft(3,imirror))<eps1) THEN
                     sg_number=100
                  ELSEIF (ABS(ft(2,imirror))<eps1) THEN
                     sg_number=103
                  ELSE
                     sg_number=104
                  ENDIF
               ELSEIF (group_desc_sg(irot90) == 74) THEN
!
!  look at the mirror perpendicular to the x axis.
!
                  IF (ABS(ft(3,imirror))<eps1.AND.&
                                        ABS(ft(2,imirror))<eps1) THEN
!
!   proper mirror
!
                     sg_number=105
                  ELSEIF (ABS(ft(3,imirror))<eps1) THEN

!   glide perpendicular to the z axis
!
                     sg_number=106
                  ELSEIF (ABS(ft(2,imirror))<eps1) THEN
!
!   glide parallel to the z axis
!
                     sg_number=101
                  ELSE
!
!   glide neither parallel nor perpendicular to the z axis
!
                     sg_number=102
                  END IF
               END IF
            END IF
         ELSEIF (ibrav==7) THEN
            IF (is_symmorphic) THEN
               sg_number=107
            ELSE 
!
!  bring the fractional translation associated to the 90 degree rotation
!  in the cartesian axis
!
               ftrs1(:) = ftc(:,irot90) 

               IF ( ABS((ABS(ftrs1(3) / at(3,1))-1.0_DP))<1.d-6 .OR. &
                     ABS(ftrs1(3)) < 1.d-6 ) THEN
!
!  the axis of order four is proper
!
                  sg_number=108
               ELSE 
!
!  the axis of order four is improprer
!
                  ftrs2(:) = ftc(:,imirror) 
                  nft=0
                  IF ((ABS(ftrs2(2)) < 1.d-6 .AND. ABS(ftrs2(3)) > 1.d-6).OR.&
                      (ABS(ftrs2(3)) < 1.d-6 .AND. ABS(ftrs2(2)) > 1.d-6) ) &
                      nft=nft+1

                  IF (nft==0) THEN
!
!   the mirror is proper. Origin on the 2 axis
!
                     sg_number=109
                     s01(:) = s0(:,3)
                  ELSE
!
!  origin on the 2 axis.
!
                     sg_number=110
                     s01(:) = s0(:,3)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      CASE(15)
!
!  C_6v
!
         s01(:)=s0(:,2)
         IF (is_symmorphic) THEN
            sg_number=183
         ELSE 
            irot60=2
            imirror=7

            IF (ABS(ft(3,irot60)) <eps1 ) THEN
               sg_number=184
            ELSEIF (group_desc_sg(irot60) == 110) THEN
               IF (ABS(ft(2,imirror))<eps1.AND.ABS(ft(3,imirror))<eps1) THEN
                  sg_number=186
               ELSE
                  sg_number=185
               ENDIF
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
          IF (ABS(ibrav)==13) THEN
             IF (is_symmorphic) THEN
                sg_number=12
             ELSE 
                sg_number=15
             ENDIF
          ELSEIF(ABS(ibrav)==12) THEN
             IF (is_symmorphic) THEN
                sg_number=10
             ELSE

                irot180=2
                imirror=4

                IF (ibrav==-12) THEN
!
!   b unique
!
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
            IF (is_symmorphic) THEN
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
         nsa=0
         nsaz=0
         iaxis=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==4) THEN
               iaxis=iaxis+1
               CALL versor(sr(1,1,isym),ax)
               DO ipol=1,3
                  IF (is_axis(ax,ipol)) THEN
                     IF (ABS(ft(ipol,isym)) >eps1) THEN
                        nsa=nsa+1
                        IF (ipol==3) nsaz=nsaz+1
                     ELSE
                        ipa=iaxis
                     ENDIF
                     axis_dir(iaxis)=ipol
                  ENDIF
               ENDDO
               DO jsym=2,nsym
                  IF (tipo_sym(sr(1,1,jsym))==5) THEN
                     CALL mirror_axis(sr(1,1,jsym),bx)
                     sp=ax(1)*bx(1) + ax(2)*bx(2) + ax(3)*bx(3)
                     IF (ABS(sp)> 1.d-6) mirr_sym(iaxis)=jsym
                  ENDIF
               ENDDO 
            ENDIF
         ENDDO

         IF (ibrav==8) THEN
!
!   Count the number of real mirrors nmir, the number of glide planes with
!   glide vector parallel to the axis nmp. type_group is used only when there
!   is only one proper axis and type_group=1 if the mirror perpendicular 
!   to the proper axis is a glide plane with glide vector parallel to the 
!   x, y or z axis.
!
            nmir=0
            nmp=0
            type_group=0
            DO iaxis=1,3
               ftrs1(:)=ft(:,mirr_sym(iaxis))
               ftrs1(axis_dir(iaxis))=0.0_DP
               nftd=0.0_DP
               ncomp=0
               DO ipol=1,3
                  nftd=nftd+ftrs1(ipol)**2
                  IF (ftrs1(ipol) /= 0) ncomp=ncomp+1
               END DO  
               IF (nftd<eps1) nmir=nmir+1
               IF (ncomp==1) nmp=nmp+1
               IF (iaxis==ipa.AND.ncomp==1) type_group=1 
            END DO

            IF (is_symmorphic) THEN
               sg_number=47
            ELSE
               IF (nsa==0) THEN
                  IF (nmir==1) THEN
                     sg_number=49
                  ELSE
                     IF (nmp==0) THEN
                        sg_number=48
!
!  origin choice 1 at the intesection of the 2 axes
!
                        s01(:)=s0(:,2)
                        s01(3)=s0(3,3)
!
!   origin choice 2 on inversion
!
                        s02(:)=s0(:,5)
                     ELSEIF (nmp==2) THEN
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

                     ENDIF
                  ENDIF
               ELSEIF (nsa==1) THEN
                 IF (nmir==1) THEN
                    sg_number=53
                 ELSEIF (nmir==2) THEN
                    sg_number=51
                 ELSEIF (nmir==0) THEN
                     IF (nmp==1) THEN
                        sg_number=52
                     ELSEIF (nmp==3) THEN
                        sg_number=54
                     ENDIF
                 ENDIF
               ELSEIF (nsa==2) THEN
                  IF (nmir==0) THEN
                     IF (type_group==1) THEN
                        sg_number=60
                     ELSE
                        sg_number=56
                     ENDIF
                  ELSEIF (nmir==1) THEN
                     IF (nmp==0) THEN
                        sg_number=58
                     ELSEIF (nmp==2) THEN
                        IF (type_group==1) THEN
                           sg_number=57
                        ELSE
                           sg_number=55
                        ENDIF
                     ENDIF
                  ELSEIF (nmir==2) THEN
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
                  ENDIF
               ELSEIF (nsa==3) THEN
                 IF (nmir==1) THEN
                    sg_number=62
                 ELSEIF (nmir==0) THEN
                    sg_number=61
                 ENDIF  
               ENDIF 
            ENDIF
         ELSEIF (ibrav==9) THEN
            IF (is_symmorphic) THEN
               sg_number=65
            ELSE
               IF (nsaz==0) THEN

                  DO iaxis=1,3
                     IF (axis_dir(iaxis)==3) ftrs1(:) = ftc(:,mirr_sym(iaxis))
                  ENDDO

                  nft=0
                  IF ((ABS(ftrs1(1)) < 1.d-6 .AND. ABS(ftrs1(2)) > 1.d-6) .OR. &
                     (ABS(ftrs1(2)) < 1.d-6 .AND. ABS(ftrs1(1)) > 1.d-6) ) &
                                                               nft=nft+1
!
!   check only the mirror perpendicular to x, the other is equal
!

                  nmir=0
                  DO iaxis=1,3
                     IF (axis_dir(iaxis)==1) THEN
                        IF (ABS(ft(3,mirr_sym(iaxis)))<eps1) nmir=nmir+1 
                     ENDIF
                  ENDDO
                  nmir=nmir*2

                  IF (nft==0) THEN
                     sg_number=66
                  ELSEIF (nft==1.AND.nmir==0) THEN
                     sg_number=68
!
!  origin choice 1 at the intesection of the 222 axis  
!
                     s01(:)=s0(:,2)
                     s01(3)=s0(3,3)
!
!   origin choice 2 on inversion
!
                     s02(:)=s0(:,5)
                  ELSEIF (nft==1.AND.nmir==2) THEN
                     sg_number=67
                  ENDIF
               ELSEIF (nsaz==1) THEN
!
!   check only the mirror parallel to the centered C face
!
                  DO iaxis=1,3
                     IF (axis_dir(iaxis)==3) ftrs1(:) = ftc(:,mirr_sym(iaxis))
                  ENDDO

                  nft=0
                  IF ((ABS(ftrs1(1)) < 1.d-6 .AND. ABS(ftrs1(2)) > 1.d-6) .OR. &
                 (ABS(ftrs1(2)) < 1.d-6 .AND. ABS(ftrs1(1)) > 1.d-6) ) nft=nft+1

                  IF (nft==0) THEN
                     sg_number=63
                  ELSEIF (nft==1) THEN
                     sg_number=64
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==10) THEN
            IF (is_symmorphic) THEN
               sg_number=69
            ELSE
               s01(:)=(s0c(:,2)+s0c(:,3)+s0c(:,4))/3.0_DP
               CALL cryst_to_cart(1,s01,bg,-1)
               s02(:)=s0(:,5)
               sg_number=70
            ENDIF
         ELSEIF (ibrav==11) THEN
            IF (is_symmorphic) THEN
               sg_number=71
            ELSE
!
!   check the three mirrors 
!
               DO iaxis=1,3
                  IF (axis_dir(iaxis)==1) ftrs1(:) = ftc(:,mirr_sym(iaxis))
                  IF (axis_dir(iaxis)==2) ftrs2(:) = ftc(:,mirr_sym(iaxis))
                  IF (axis_dir(iaxis)==3) ftrs3(:) = ftc(:,mirr_sym(iaxis))
               ENDDO

               nft=0
               IF ((ABS(ftrs1(2)) < 1.d-6 .AND. ABS(ftrs1(3)) > 1.d-6) .OR. &
               (ABS(ftrs1(3)) < 1.d-6 .AND. ABS(ftrs1(2)) > 1.d-6) ) nft=nft+1
               IF ((ABS(ftrs2(1)) < 1.d-6 .AND. ABS(ftrs2(3)) > 1.d-6) .OR. &
               (ABS(ftrs2(3)) < 1.d-6 .AND. ABS(ftrs2(1)) > 1.d-6) ) nft=nft+1
               IF ((ABS(ftrs3(2)) < 1.d-6 .AND. ABS(ftrs3(1)) > 1.d-6) .OR. &
               (ABS(ftrs3(1)) < 1.d-6 .AND. ABS(ftrs3(2)) > 1.d-6) ) nft=nft+1

               IF (nft==1) THEN
                  sg_number=74
               ELSEIF (nft==2) THEN
                  sg_number=72
               ELSEIF (nft==3) THEN
                  sg_number=73
               ENDIF
            ENDIF
         ENDIF
      CASE(21)
!
!  D_3h
!
!  origin on the inversion center of -6
!
         s01(:)=s0(:,2)

         IF (is_symmorphic) THEN
            IF (group_desc(7)==4) THEN
               sg_number=189
            ELSE
               sg_number=187
            ENDIF
         ELSE
            IF (group_desc(7)==4) THEN
!
!  origin on 32c
!
               s01(:)=s0(:,3)
               s01(3)=s0(3,7)
               sg_number=190
            ELSE
!
!  origin on 3c2
!
               s01(:)=s0(:,3)
               s01(3)=s0(3,10)
               sg_number=188
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
            IF (is_symmorphic) THEN
               sg_number=123
            ELSE
               IF (ABS(ft(3,irot90))<eps1) THEN
                  IF (ABS(ft(1,irot180))<eps1) THEN
                     IF (ABS(ft(2,imirror))<eps1) THEN
!
!  glide plane c
!
                        sg_number=124
                     ELSEIF(ABS(ft(3,imirror))<eps1) THEN
!
!  glide plane b
!
                        sg_number=125
!
!  two origins origin 1 on the inversion point of -4, origin 2 on the inversion
!  center
!
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ELSE
!
!  glide plane n
!
                        sg_number=126
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ENDIF
                  ELSEIF (group_desc_sg(irot180) == 67) THEN
                     IF (ABS(ft(2,imirror))<eps1.AND.&
                         ABS(ft(3,imirror))<eps1) THEN
!
!  mirror
!
                        sg_number=129
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ELSEIF (ABS(ft(2,imirror))<eps1) THEN
!
!  glide plane c
!
                        sg_number=130
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ELSEIF(ABS(ft(3,imirror))<eps1) THEN
!
!  glide plane b
!
                        sg_number=127
                     ELSE
!
!  glide plane n
!
                        sg_number=128
                     ENDIF
                  ENDIF
               ELSEIF (group_desc_sg(irot90) == 74) THEN
                  IF (ABS(ft(1,irot180))<eps1) THEN
                     IF (ABS(ft(2,imirror))<eps1.AND.&
                         ABS(ft(3,imirror))<eps1) THEN
!
!  mirror
!
                        sg_number=131
                     ELSEIF (ABS(ft(2,imirror))<eps1) THEN
!
!  glide plane c
!
                        sg_number=132
                     ELSEIF(ABS(ft(3,imirror))<eps1) THEN
!
!  glide plane b
!
                        sg_number=133
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ELSE
!
!  glide plane n
!
                        sg_number=134
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ENDIF
                  ELSEIF (group_desc_sg(irot180) == 67) THEN
                     IF (ABS(ft(2,imirror))<eps1.AND.&
                         ABS(ft(3,imirror))<eps1) THEN
!
!  mirror
!
                        sg_number=137
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ELSEIF (ABS(ft(2,imirror))<eps1) THEN
!
!  glide plane c
!
                        sg_number=138
                        s01(:)=s0(:,10)
                        s02(:)=s0(:,9)
                     ELSEIF(ABS(ft(3,imirror))<eps1) THEN
!
!  glide plane b
!
                        sg_number=135
                     ELSE
!
!  glide plane n
!
                        sg_number=136
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (is_symmorphic) THEN
               sg_number=139
            ELSE
!
!    Bring the fractional translation associated with the fourfold rotation
!    to cartesian axis
!
               ftrs1(:) = ftc(:,irot90)

               IF ( group_desc_sg(irot90) == 8 ) THEN
!
!  the axis of order four is proper
!
                  sg_number=140
               ELSE
!
!  the axis of order four is a screw axis
!
                  ftrs2(:) = ftc(:,imirror) 
                  nft=0
                  IF ((ABS(ftrs2(2)) < 1.d-6 .AND. ABS(ftrs2(3)) > 1.d-6).OR.&
                      (ABS(ftrs2(3)) < 1.d-6 .AND. ABS(ftrs2(2)) > 1.d-6) ) &
                      nft=nft+1

                  IF (nft==0) THEN
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
         ENDIF 
      CASE(23)
!
!  D_6h
!
!  origin on the inversion center of the -3_z operation
!

         s01(:)=s0(:,18)
         IF (is_symmorphic) THEN
            sg_number=191
         ELSE
            irot60=2
            imirror=19
            IF (ABS(ft(3,irot60))<eps1) THEN
               sg_number=192
            ELSEIF (group_desc_sg(irot60)==110) THEN
               IF (ABS(ft(2,imirror))<eps1.AND.ABS(ft(3,imirror))<eps1) THEN
                  sg_number=194
               ELSE
                  sg_number=193
               ENDIF
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
            IF (is_symmorphic) THEN
               IF (type_group==0) THEN
                  sg_number=115
               ELSE
                  sg_number=111
               ENDIF
            ELSE
               IF (type_group==1) THEN
                  IF (ABS(ft(1,irot180))<eps1) THEN
                     sg_number=112
                  ELSE
                    IF (ABS(ft(3,imirror1))<eps1) THEN
                       sg_number=113
                    ELSE
                       sg_number=114
                    ENDIF
                  ENDIF
               ELSE
                  IF (ABS(ft(3,imirror))<eps1) THEN
                     sg_number=117
                  ELSEIF (ABS(ft(2,imirror))<eps1) THEN
                     sg_number=116
                  ELSE 
                     sg_number=118
                  ENDIF 
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (is_symmorphic) THEN
               IF (type_group==1) THEN
                  sg_number=121
               ELSE
                  sg_number=119
               ENDIF
            ELSE
               IF (type_group==1) THEN
                  sg_number=122
               ELSE
                  sg_number=120
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
!   if all the mirrors are glide planes it is 165, otherwise it is 164
!
               IF (group_desc_sg(10)>64.AND.group_desc_sg(11)>64&
                   .AND.group_desc_sg(12)>64) THEN
                  sg_number=165
               ELSE
                  sg_number=164
               ENDIF
            ELSE
!
!   if all the mirrors are glide planes it is 163, otherwise it is 164
!
               IF (group_desc_sg(10)>64.AND.group_desc_sg(11)>64&
                   .AND.group_desc_sg(12)>64) THEN
                  sg_number=163
               ELSE
                  sg_number=162
               ENDIF
            ENDIF
         ELSEIF (ibrav==5) THEN
             IF (is_symmorphic) THEN
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
             IF (is_symmorphic) THEN
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
             IF (is_symmorphic) THEN
                sg_number=202
             ELSE
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
             ENDIF
         ELSEIF (ibrav==3) THEN
             IF (is_symmorphic) THEN
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
             IF (is_symmorphic) THEN
                sg_number=215
             ELSE
                s01(:)=s0(:,4)
                s01(3)=s0(3,2)
                sg_number=218
             ENDIF
         ELSEIF (ibrav==2) THEN
             IF (is_symmorphic) THEN
                sg_number=216
             ELSE
                s01(:)=s0(:,4)
                s01(3)=s0(3,2)
                sg_number=219
             ENDIF
         ELSEIF (ibrav==3) THEN
             IF (is_symmorphic) THEN
                sg_number=217
             ELSE
!
!    origin to be programmed
!
                s01=1000.0_DP
                sg_number=220
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
         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=207
             ELSE
!
!          find the 90 degree rotation about the x axis and check its
!          fractional translation along the axis.
!
                irot90=14
!
!    4x is 4_2 in 208, 4_1 in 213, 4_3 in 212
!
                IF (group_desc_sg(14)==90) sg_number=208
                IF (group_desc_sg(14)==89) THEN
                   sg_number=213
!
!  origin to be programmed
!
                   s01=1000.0_DP
                ENDIF
                IF (group_desc_sg(14)==91) THEN
                   sg_number=212
!
!  origin to be programmed
!
                   s01=1000.0_DP
                ENDIF
             ENDIF
         ELSEIF (ibrav==2) THEN
             IF (is_symmorphic) THEN
                sg_number=209
             ELSE
                sg_number=210
             ENDIF
         ELSEIF (ibrav==3) THEN
             IF (is_symmorphic) THEN
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
            IF (is_symmorphic) THEN
               sg_number=221
            ELSE
               IF (group_desc_sg(18)==8) THEN
!
!   the 4z axis has no fractional translation in space group 222
!
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
               ELSE
                  IF (group_desc_sg(28)==34) THEN
!
!  in space group 223 the mirror i2z is a mirror, in 224 a glide plane
!
                     sg_number=223
                  ELSE
                     sg_number=224
!
!   origin choice 1 at -43m
!
                     s01(:)=s0(:,18)
                     s01(3)=s0(3,28)
!
!  origin choice 2 at the inversion point of -3
!
                     s02(:)=s0(:,29)
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==2) THEN
            IF (is_symmorphic) THEN
               sg_number=225
            ELSE
               IF (group_desc_sg(28)==34) THEN
!
!  in this case i2z is a pure mirror
!
                  sg_number=226
               ELSE
                  IF (group_desc_sg(47)==38.OR.group_desc_sg(48)==37) THEN
!
!  one of i2xy or i2x-y is a pure mirror in 227 and a glide plane in 228
!
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
            ENDIF
         ELSEIF (ibrav==3) THEN
            IF (is_symmorphic) THEN
               sg_number=229
            ELSE
               sg_number=230
            ENDIF
         ENDIF
  END SELECT

  RETURN
  END SUBROUTINE find_space_group


  SUBROUTINE set_spg_code(code)
  INTEGER, INTENT(IN) :: code
  
  spg_code=code

  RETURN
  END SUBROUTINE set_spg_code


  SUBROUTINE transform_fcc_axis(tau_aux,naux)
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

  SUBROUTINE transform_bcc_axis(tau_aux,naux)
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

  SUBROUTINE transform_fco_axis(tau_aux,naux)
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

  SUBROUTINE transform_ct_axis(tau_aux,naux)
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

  SUBROUTINE transform_obco_c_axis(tau_aux,naux)
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

  SUBROUTINE transform_monoclinic_center_c(tau_aux,naux)
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

  SUBROUTINE transform_monoclinic_center_b(tau_aux,naux)
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

  SUBROUTINE transform_obco_a_axis(tau_aux,naux)
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

  SUBROUTINE transform_bco_axis(tau_aux,naux)
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

  SUBROUTINE transform_trigonal_axis(tau_aux,naux)
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

  SUBROUTINE find_space_group_number(sgroup_name, group_code)
!
!  This routine returns the space group number of a given group_name if
!  the name is recognized or zero if it is not.
!
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: sgroup_name
  INTEGER, INTENT(OUT) :: group_code
  
  INTEGER igroup

  group_code=0
  DO igroup=1,nsg
     IF (TRIM(sgroup_name)==TRIM(space_group_names(igroup)(1:11))) THEN
        group_code=space_group_numbers(igroup)
        RETURN
     ENDIF 
  ENDDO
  RETURN
  END SUBROUTINE find_space_group_number

  SUBROUTINE find_space_group_names(group_code)
!
!  This routine writes on output all the names that correspond to a given
!  space group.
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_code
  
  INTEGER igroup

  DO igroup=1,nsg
     IF (space_group_numbers(igroup)==group_code) &
!     WRITE(6,'(i5,2x,a,3x,a)') igroup, TRIM(space_group_names(igroup)), &
!                            TRIM(add_info(igroup))
     WRITE(6,'(5x,a,3x,a)') TRIM(space_group_names(igroup)), &
                            TRIM(add_info(igroup))
  ENDDO
  RETURN
  END SUBROUTINE find_space_group_names

  SUBROUTINE sg_name(group_code, output_name)
!
!  This routine receives a space group code and gives as output the 
!  space group name
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: group_code
  CHARACTER(LEN=12), INTENT(OUT) :: output_name
  
  INTEGER igroup

  DO igroup=1,nsg
     IF (space_group_numbers(igroup)==group_code) THEN
         output_name=space_group_names(igroup)
         RETURN
     ENDIF
  ENDDO
  CALL errore('sg_name','group_code unknown',1)
  RETURN
  END SUBROUTINE sg_name

  SUBROUTINE project_frac_tran(sr,ft,at,bg,ftperp,ftpar)
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

  SUBROUTINE shift_frac_tran(sr,ft,at,bg,s0,ftnew)
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

  SUBROUTINE find_sg_tags(sr,ft,ibrav,at,bg,ftperp,ftpar,sym_in,s0,sg_sym_out)
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

SUBROUTINE find_mirror_type(ifrac,imirror)
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

SUBROUTINE find_ifact(axis_order,ifrac,ifact)
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

SUBROUTINE find_origin_xyz_rot(sr,ftperp,sym_in,at,bg,s0,ipol)
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

SUBROUTINE find_origin_3_rot(sr,ftperp,sym_in,at,bg,s0)
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

SUBROUTINE set_sg_ibrav(sgc, ibrav)

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

SUBROUTINE set_point_group_code(sgc, group_code, group_code_ext)
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

SUBROUTINE set_ft_generators(sgc,ft)
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

LOGICAL FUNCTION check_space_group_ft(ft1,ftref,nsym)
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

SUBROUTINE set_standard_sg(sgc,ibrav,celldm,nsym,sr,ft)
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

SUBROUTINE space_group_multiply(sr1, ft1, sr2, ft2, sr3, ft3)
!
!   This routine multiplies two space group operations
!   sr1, sr2 are two 3x3 rotation matrices in real space
!   ft1, ft2, are two fractionary translations in real space
!
IMPLICIT NONE
REAL(DP), INTENT(IN)  :: sr1(3,3), sr2(3,3)
REAL(DP), INTENT(IN)  :: ft1(3), ft2(3)
REAL(DP), INTENT(OUT) :: sr3(3,3)
REAL(DP), INTENT(OUT) :: ft3(3)

REAL(DP) :: ftaux(3)
INTEGER  :: jpol

sr3 = MATMUL(sr1, sr2)

ftaux=0.0_DP
DO jpol=1,3
   ftaux(:) = ftaux(:) + sr1(:,jpol) * ft2(jpol)
END DO
ft3 = ft1 + ftaux

RETURN
END SUBROUTINE space_group_multiply

END MODULE space_groups
