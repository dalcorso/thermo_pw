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

  INTEGER, PARAMETER :: nsg = 1056             ! number of space groups
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
         "P12_1/a1     ", "P12_1/n1     ", "B12_1/a1     ", "B12_1/d1     ", &
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
         "D_2^2       S", "D_2^3       S", "D^2_4       S", "D_2^5       S", &
         "D_2^6       S", "D_2^7       S", "D_2^8       S", "D_2^9       S", &
         "C_2v^1      S", "C_2v^2      S", "C_2v^3      S", "C_2v^4      S", &
         "C_2v^5      S", "C_2v^6      S", "C_2v^7      S", "C_2v^8      S", &
         "C_2v^9      S", "C_2v^10     S", "C_2v^11     S", "C_2v^12     S", &
         "C_2v^13     S", "C_2v^14     S", "C_2v^15     S", "C_2v_16     S", &
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
         "S_4^1       S", "S^4^2       S", "C_4h^1      S", "C_4h^2      S", &
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
         "C_6v^3      S", "C_6v^1      S", "D_3h^1      S", "D_3h^2      S", &
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
         "V_d^9       S", "V_d^10      S", "V_d^11      S", "V_d^12      S" /

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
       14,              14,            14,            14,            &
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
      119,             120,           121,           122              /
 
  PUBLIC   spg_code, spg_name, set_spg_code, find_space_group, sg_name, &
           equivalent_tau, space_group_names, &
           find_space_group_number, find_space_group_names, add_info, &
           set_add_info, set_fft_fact

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

  SUBROUTINE find_space_group(sg_number, ibrav, code_group, nsym, s, sr, ftau, &
                              at, nr1, nr2, nr3, verbosity)
  !
  !  This routine receives as input: 
  !  the bravais lattice, the point group, the number of symmetries, the
  !  rotation matrices in crystal and in carthesian coordinates, the
  !  fractional translations in number of fft points, the fft mesh nr1, nr2, nr3
  !  and gives as output the space group number.
  !  if verbosity .true. the routine writes informations on the symmetries
  !  ibrav=0 is not supported
  !
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ibrav, code_group, nsym, nr1, nr2, nr3
  INTEGER, INTENT(IN) :: s(3,3,nsym), ftau(3,nsym)
  LOGICAL, INTENT(IN) :: verbosity
  REAL(DP) :: sr(3,3,nsym), at(3,3) ! rotation matrices in cartesian coordinates
  INTEGER, INTENT(OUT) :: sg_number
  LOGICAL :: is_symmorphic
  INTEGER :: isym, nft, iaxis, imirror, irot90, irot60, irot120, ipol, &
             imirror1, imirror2, nsaz
  INTEGER :: idir, irot180, type_group, nsp, nsa, nmir, ncomp, jsym, nmp, ipa
  INTEGER :: ft1(3), ft2(3), sp1, sp2, axis_dir(3), mirr_sym(3), idir1, idir2
  REAL(DP) :: ax(3), bx(3), angle, fact, ft(3), sp, ftrs1(3), ftrs2(3), &
              ftrs3(3), ashift(3), acryshift(3), ashift_mod
  REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP, eps=1.D-8
  LOGICAL :: is_axis
  INTEGER :: tipo_sym, ts
  REAL(DP) :: angle_rot
  
  IF (ibrav==0) THEN
     sg_number=0 
     RETURN
  ENDIF

  IF (verbosity) THEN
     WRITE(stdout,'(/,5x,"Space group identification, ",i3," symmetries:")') nsym
     WRITE(stdout,'(5x,"(Optional origin shifts are indicated for C2, mirrors,")')
     WRITE(stdout,'(5x,"and inversion symmetries)",/)')
  ENDIF
  is_symmorphic=.TRUE.
  DO isym=1, nsym
     ts=tipo_sym_sg(sr(1,1,isym), ftau(1,isym), at, nr1, nr2, nr3, ashift, acryshift)
     ashift_mod=SQRT(ashift(1)**2 + ashift(2)**2 + ashift(3)**2 )
     IF (verbosity) THEN
        WRITE(stdout,'(5x,i5,"  ",a)') isym, TRIM(label_sg(ts))
        IF (ts < 17 .AND. ashift_mod > eps) THEN
           WRITE(stdout,'(11x,"The origin is not a fixed point of this operation")')
           WRITE(6,'(11x,"To center it on the origin add to all atomic positions")')
           WRITE(6,'(11x,"(alat units)")')
           WRITE(6,'(11x,"(",3f15.10,")",/)') ashift(:)
           WRITE(6,'(11x,"(cryst. coord.)")')
           WRITE(6,'(11x,"(",3f15.10,")",/)') acryshift(:)
        ENDIF
     ENDIF
     is_symmorphic = (is_symmorphic .AND. ts < 17 )
  ENDDO

!  is_symmorphic=.NOT.(ANY(ftau(:,1:nsym) /= 0))
  sg_number=0
  SELECT CASE (code_group)
      CASE(1)
!
!  1
!
         sg_number=1
      CASE(2)
!
!  -1
!
         sg_number=2
      CASE(3)
!
!  C_s
!
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
         IF (ibrav==4) THEN
            IF (is_symmorphic) THEN
               sg_number=143
            ELSE
              irot120=0
              DO isym=2,nsym
                 IF (tipo_sym(sr(1,1,isym))==3) THEN
                    IF (ABS(angle_rot(sr(1,1,isym))-120.0_DP)< 1.D-6) &
                        irot120=isym
                 ENDIF
              ENDDO
              IF ( irot120 == 0 ) CALL errore('find_space_group',&
                                         'C_3: no 120 degrees rotation',1)
              IF (ftau(3,irot120) * 3 / nr3 == -1 .OR. &
                  ftau(3,irot120) * 3 / nr3 == 2 ) sg_number=144
              IF (ftau(3,irot120) * 3 / nr3 == 1 .OR. &
                  ftau(3,irot120) * 3 / nr3 == -2 ) sg_number=145
           ENDIF
         ELSEIF (ibrav==5) THEN
            sg_number=146
         ENDIF
      CASE(6)
!
!   C_4
!
         IF (ibrav==6) THEN
            IF (is_symmorphic) THEN
               sg_number=75
            ELSE
               irot90=0
               DO isym=2,nsym
                  IF (tipo_sym(sr(1,1,isym))==3) THEN
                     IF (ABS(angle_rot(sr(1,1,isym))-90.0_DP)< 1.D-6) &
                        irot90=isym
                  ENDIF
               ENDDO
               IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                          'C_4: no 90 degrees rotation',1)
               IF (ftau(3,irot90) * 4 / nr3 == -1 .OR. &
                                ftau(3,irot90) * 4 / nr3 == 3) sg_number=76
               IF (ABS(ftau(3,irot90) * 4 / nr3) == 2) sg_number=77
               IF (ftau(3,irot90) * 4 / nr3 == -3 .OR. &
                                ftau(3,irot90) * 4 / nr3 == 1) sg_number=78
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (is_symmorphic) THEN
               sg_number=79
            ELSE
               sg_number=80
            ENDIF
         ENDIF
      CASE(7)
!
!  C_6
!
         IF (is_symmorphic) THEN
            sg_number=168
         ELSE
            irot60=0
            DO isym=2,nsym
               IF (tipo_sym(sr(1,1,isym))==3) THEN
                  IF (ABS(angle_rot(sr(1,1,isym))-60.0_DP)< 1.D-6) &
                         irot60=isym
               ENDIF
            ENDDO
            IF ( irot60 == 0 ) CALL errore('find_space_group',&
                                          'C_6: no 60 degrees rotation',1)
            IF (ftau(3,irot60) * 6 / nr3 == -1 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == 5) sg_number=169
            IF (ftau(3,irot60) * 6 / nr3 == -2 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == 4) sg_number=171
            IF (ABS(ftau(3,irot60) * 6 / nr3) == 3) sg_number=173
            IF (ftau(3,irot60) * 6 / nr3 == 2 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == -4) sg_number=172
            IF (ftau(3,irot60) * 6 / nr3 == 1 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == -5) sg_number=170
         ENDIF
      CASE(8)
!
!   D_2
!
         IF (ibrav==8) THEN
            IF (is_symmorphic) THEN
               sg_number=16
            ELSE
               nft=0
               DO isym=2,nsym
                  CALL versor(sr(1,1,isym), ax) 
                  DO ipol=1,3
                     IF (is_axis(ax,ipol).AND.ftau(ipol,isym) /=0) nft=nft+1
                  ENDDO
               ENDDO
               IF ( nft == 0 ) CALL errore('find_space_group',&
                              'D_2: unexpected number of screw axis',1)
               IF (nft==1) sg_number=17
               IF (nft==2) sg_number=18
               IF (nft==3) sg_number=19
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
            ENDIF
         ENDIF
      CASE(9)
!
!  D_3
!
         IF (ibrav==4) THEN
            IF (is_symmorphic) THEN
               DO isym=2,nsym
                  IF (tipo_sym(sr(1,1,isym))==4) THEN
!
!  search the axis of order 2 and compute the angle with the x axis
!
                     CALL versor(sr(1,1,isym), ax)
                     angle=ACOS(ax(1))*180.0_DP / pi
                     IF (MOD(NINT(angle), 60)==0) THEN
!
!   axis of order 2 on the hexagon edges
!
                        sg_number=150
                     ELSE
!
!   axis of order 2 on the hexagon diagonals
!
                        sg_number=149
                     ENDIF
                     EXIT
                  ENDIF
               ENDDO
            ELSE
!
!           find the rotation of 120
!
               irot120=0
               DO isym=2,nsym
                  IF (tipo_sym(sr(1,1,isym))==3) THEN
                     IF (ABS(angle_rot(sr(1,1,isym))-120.0_DP)< 1.D-6) &
                         irot120=isym
                  ENDIF
               ENDDO
               IF ( irot120 == 0 ) CALL errore('find_space_group',&
                                          'D_3: no 120 degrees rotation',1)
               IF (ftau(3,irot120) * 3 / nr3 == -1 .OR. &
                   ftau(3,irot120) * 3 / nr3 == 2 ) THEN
                  DO isym=2,nsym
                     IF (tipo_sym(sr(1,1,isym))==4) THEN
!
!  search the axis of order 2 and compute the angle with the x axis
!
                        CALL versor(sr(1,1,isym), ax)
                        angle=ACOS(ax(1))*180.0_DP / pi
                        IF (MOD(NINT(angle), 60)==0) THEN
!
!   axis of order 2 on the hexagonal edges
!
                           sg_number=152
                        ELSE
                           sg_number=151
                        ENDIF
                        EXIT
                     ENDIF
                  ENDDO
               ELSEIF (ftau(3,irot120) * 3 / nr3 == 1 .OR.  &
                       ftau(3,irot120) * 3 / nr3 == -2 ) THEN
                  DO isym=2,nsym
                     IF (tipo_sym(sr(1,1,isym))==4) THEN
!
!  search the axis of order 2 and compute the angle with the x axis
!
                        CALL versor(sr(1,1,isym), ax)
                        angle=ACOS(ax(1))*180.0_DP / pi
                        IF (MOD(NINT(angle), 60)==0) THEN
!
!   axis of order 2 on the hexagonal edges
!
                           sg_number=154
                        ELSE
                           sg_number=153
                        ENDIF
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
         ELSEIF (ibrav==5) THEN
            sg_number=155
         ENDIF     
      CASE(10)
!
!  D_4
!
         IF (ibrav==6) THEN
            IF (is_symmorphic) THEN
               sg_number=89
            ELSE 
!
!   Find the axis of order 4, and the axis of order 2 parallel to x
!
               irot90=0
               irot180=0
               DO isym=2,nsym
                  IF (tipo_sym(sr(1,1,isym))==3) THEN
                     IF (ABS(angle_rot(sr(1,1,isym))-90.0_DP)< 1.D-6) &
                                                  irot90=isym
                  ENDIF
                  IF (tipo_sym(sr(1,1,isym))==4) THEN
                     CALL versor(sr(1,1,isym), ax)
                     IF (is_axis(ax,1)) irot180=isym
                  ENDIF
               ENDDO
               IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                          'D_4: no 90 degrees rotation',1)
               IF ( irot180 == 0 ) CALL errore('find_space_group',&
                                          'D_4: no 180 degrees rotation',1)

               IF (ABS(ftau(3,irot90)) == 0) THEN
!
!  axis 4
!
                  sg_number=90
               ELSEIF (ftau(3,irot90) * 4 / nr3 == -1 .OR. &
                       ftau(3,irot90) * 4 / nr3 == 3 ) THEN
!
!  axis 4_1
!
                  IF (ftau(1,irot180) == 0) THEN
                     sg_number=91
                  ELSE
                     sg_number=92
                  ENDIF
               ELSEIF (ABS(ftau(3,irot90) * 4 / nr3) == 2) THEN
!
!  axis 4_2
!
                  IF (ftau(1,irot180) == 0) THEN
                     sg_number=93
                  ELSE
                     sg_number=94
                  ENDIF
               ELSEIF (ftau(3,irot90) * 4 / nr3 == 1 .OR. &
                       ftau(3,irot90) * 4 / nr3 == -3 ) THEN
!
!  axis 4_3
!
                  IF (ftau(1,irot180) == 0) THEN
                     sg_number=95
                  ELSE
                     sg_number=96
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==7) THEN
            IF (is_symmorphic) THEN
               sg_number=97
            ELSE 
               sg_number=98
            ENDIF
         ENDIF
      CASE(11)
!
!  D_6
!
         IF (is_symmorphic) THEN
            sg_number=177
         ELSE
            irot60=0
            DO isym=2,nsym
               IF (tipo_sym(sr(1,1,isym))==3) THEN
                  IF (ABS(angle_rot(sr(1,1,isym))-60.0_DP)< 1.D-6) &
                         irot60=isym
               ENDIF
            ENDDO
            IF ( irot60 == 0 ) CALL errore('find_space_group',&
                                           'D_6: no 60 degrees rotation',1)
            IF (ftau(3,irot60) * 6 / nr3 == -1 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == 5) sg_number=178
            IF (ftau(3,irot60) * 6 / nr3 == -2 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == 4) sg_number=180
            IF (ABS(ftau(3,irot60) * 6 / nr3) == 3) sg_number=182
            IF (ftau(3,irot60) * 6 / nr3 == 2 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == -4) sg_number=181
            IF (ftau(3,irot60) * 6 / nr3 == 1 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == -5) sg_number=179
         ENDIF
      CASE(12)
!
!   C_2v
!
!
!  identify the axis of order 2
!
         irot180=0
         idir=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==4) THEN
               irot180=isym
               CALL versor(sr(1,1,isym),ax)
               DO ipol=1,3
                  IF (is_axis(ax,ipol)) idir=ipol 
               ENDDO
            ENDIF
         ENDDO

!   find if there is a component of the fractional translation perpendicular
!   to the two fold axis
!
         IF ( irot180 == 0 ) CALL errore('find_space_group',&
                                         'C_2v: no 180 degrees rotation',1)
         IF ( idir == 0 ) CALL errore('find_space_group',&
                                         'C_2v: unknown axis direction',1)

 
         IF (ibrav==8) THEN
            IF (is_symmorphic) THEN
               sg_number=25
            ELSE
!
!   identify the two mirrors and save the projection of the fractional
!   translation in the plane of the mirror
!
               imirror1=0
               DO isym=2, nsym
                  IF (tipo_sym(sr(1,1,isym))==5) THEN
                     IF (imirror1 == 0) THEN
                        imirror1=isym
                        CALL mirror_axis(sr(1,1,isym),ax)
!
!    find only the fractional translation component parallel to the mirror
!
                        ft1(:)=ftau(:,isym)
                        DO ipol=1,3
                           IF (is_axis(ax,ipol)) ft1(ipol)=0
                        ENDDO
                     ELSE

!    find only the fractional translation component parallel to the mirror
!
                        imirror2=isym
                        CALL mirror_axis(sr(1,1,isym),ax)
                        ft2(:)=ftau(:,isym)
                        DO ipol=1,3
                           IF (is_axis(ax,ipol)) ft2(ipol)=0 
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
               IF ( imirror1 == 0 ) CALL errore('find_space_group',&
                                                 'C_2v: no mirror found',1)
               IF ( imirror2 == 0 ) CALL errore('find_space_group',&
                                                 'C_2v: one missing mirror',1)
!
!   find if there is a component of the fractional translation perpendicular
!   to the two fold axis
!
               sp1=0
               sp2=0
               DO ipol=1,3
                  IF (ipol/=idir) THEN
                     sp1=ft1(ipol)**2+sp1
                     sp2=ft2(ipol)**2+sp2
                  ENDIF
               ENDDO
!
!  count the true mirrors
!
               nft=0
               IF (ft1(1)==0.AND.ft1(2)==0.AND.ft1(3)==0) nft=nft+1
               IF (ft2(1)==0.AND.ft2(2)==0.AND.ft2(3)==0) nft=nft+1
!
!    Now divide the groups with a proper or improper twofold axis
!
               IF (ftau(idir,irot180)==0) THEN
!
!   proper axis 
!
                  IF (nft==1) THEN
!
!  and a true mirror
!
                     sg_number=28
                  ELSE
                     IF (ft1(idir)==0.AND.ft2(idir)==0) THEN
!
!   two mirrors with glide perpendicular to the axis
!
                        sg_number=32
                     ELSEIF (sp1==0.AND.sp2==0) THEN
!
!   two mirrors with glide parallel to the axis
!
                        sg_number=27
                     ELSEIF (sp1==0.OR.sp2==0) THEN 
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
                     IF ( sp1 /= 0 .OR. sp2 /= 0 ) THEN
!
!  one glide parallel to the axis
!
                        sg_number=31
                     ELSE
                        sg_number=26
                     ENDIF
                  ELSE
!
!  two glide planes
!
                     IF ( sp1 /= 0 .AND. sp2 /= 0 ) THEN
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
               IF (ABS(ftau(3,irot180) * 2 / nr3) == 1) THEN
                  sg_number=36
               ELSE
                  sg_number=37
               ENDIF
            ENDIF
         ELSEIF (ibrav==91) THEN
!
!   find the mirror perpendicular to y
!
            imirror1=0
            DO isym=2, nsym
               IF (tipo_sym(sr(1,1,isym))==5) THEN
                  CALL mirror_axis(sr(1,1,isym),ax)
                  IF (is_axis(ax,2)) imirror1=isym
               ENDIF
            ENDDO
            IF ( imirror1 == 0 ) CALL errore('find_space_group',&
                                            'C_2v: no mirror found',1)
            IF (is_symmorphic) THEN
               sg_number=38
            ELSE 
               IF (ftau(1,imirror1)==0) THEN
                  sg_number=39
               ELSE
!
!  Find the mirror perpendicular to x
!

                  imirror2=0
                  DO isym=2, nsym
                     IF (tipo_sym(sr(1,1,isym))==5) THEN
                        CALL mirror_axis(sr(1,1,isym),ax)
                        IF (is_axis(ax,1)) imirror2=isym
                     ENDIF
                  ENDDO
                  IF ( imirror2 == 0 ) CALL errore('find_space_group',&
                                            'C_2v: no mirror found',1)

!
!  bring the fractional translation in the orthorombic axis
!
                  ftrs1(:) = ftau(1,imirror2)*at(:,1) / nr1 +  &
                             ftau(2,imirror2)*at(:,2) / nr2 +  &
                             ftau(3,imirror2)*at(:,3) / nr3
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
!  Find the the two mirrors perpendicular to the proper axis
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
               ftrs1(:) = ftau(1,imirror1)*at(:,1) / nr1 +  &
                          ftau(2,imirror1)*at(:,2) / nr2 +  &
                          ftau(3,imirror1)*at(:,3) / nr3
               ftrs2(:) = ftau(1,imirror2)*at(:,1) / nr1 +  &
                          ftau(2,imirror2)*at(:,2) / nr2 +  &
                          ftau(3,imirror2)*at(:,3) / nr3
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
         type_group=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==5) THEN
!
!  search the mirror and compute the angle of its normal the x axis
!
               CALL mirror_axis(sr(1,1,isym), ax)
               angle=ACOS(ax(1))*180.0_DP / pi
               IF (MOD(NINT(angle), 60)==0) type_group=1
               EXIT
            ENDIF
         ENDDO
        
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
         irot90=0
         imirror=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==3) THEN
               IF ( ABS(angle_rot(sr(1,1,isym))-90.0_DP )< 1.D-6 ) irot90=isym
            ENDIF
            IF (tipo_sym(sr(1,1,isym))==5) THEN
               CALL mirror_axis(sr(1,1,isym), ax)
               IF (is_axis(ax,1) ) imirror=isym
            ENDIF
         ENDDO
         IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                        'C_4v: no 90 degrees rotation',1)
         IF ( imirror == 0 ) CALL errore('find_space_group',&
                                        'C_4v: no mirror found',1)

         IF (ibrav==6) THEN
            IF (is_symmorphic) THEN
               sg_number=99
            ELSE 
               IF (ftau(3,irot90) == 0) THEN
                  IF (ftau(3,imirror) == 0) THEN
                     sg_number=100
                  ELSEIF (ftau(2,imirror)==0) THEN
                     sg_number=103
                  ELSE
                     sg_number=104
                  ENDIF
               ELSEIF (ABS(ftau(3,irot90) * 4 / nr3) == 2) THEN
!
!  look at the mirror perpendicular to the x axis.
!
                  IF (ftau(3,imirror) == 0.AND.ftau(2,imirror)==0) THEN
!
!   proper mirror
!
                     sg_number=105
                  ELSEIF (ftau(3,imirror) == 0) THEN

!   glide perpendicular to the z axis
!
                     sg_number=106
                  ELSEIF (ftau(2,imirror) == 0) THEN
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
               ftrs1(:) = ftau(1,irot90)*at(:,1) / nr1 + &
                          ftau(2,irot90)*at(:,2) / nr2 + &
                          ftau(3,irot90)*at(:,3) / nr3

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
                  ftrs2(:) = ftau(1,imirror) * at(:,1) / nr1 + &
                             ftau(2,imirror) * at(:,2) / nr2 + &
                             ftau(3,imirror) * at(:,3) / nr3
                  nft=0
                  IF ((ABS(ftrs2(2)) < 1.d-6 .AND. ABS(ftrs2(3)) > 1.d-6).OR.&
                      (ABS(ftrs2(3)) < 1.d-6 .AND. ABS(ftrs2(2)) > 1.d-6) ) &
                      nft=nft+1

                  IF (nft==0) THEN
!
!   the mirror is proper
!
                     sg_number=109
                  ELSE
                     sg_number=110
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      CASE(15)
!
!  C_6v
!
         IF (is_symmorphic) THEN
            sg_number=183
         ELSE 
            irot60=0
            imirror=0
            DO isym=2,nsym
               IF (tipo_sym(sr(1,1,isym))==3) THEN
                  IF (ABS(angle_rot(sr(1,1,isym))-60.0_DP)< 1.D-6) irot60=isym
               ENDIF
               IF (tipo_sym(sr(1,1,isym))==5) THEN
                  CALL mirror_axis(sr(1,1,isym),ax)
                  IF (is_axis(ax,1)) imirror=isym
               ENDIF
            ENDDO
            IF ( irot60 == 0 ) CALL errore('find_space_group',&
                                           'C_6v: no 60 degrees rotation',1)
            IF ( imirror == 0 ) CALL errore('find_space_group',&
                                           'C_6v: no mirror found',1)

            IF (ftau(3,irot60) == 0 ) THEN
               sg_number=184
            ELSEIF (ABS(ftau(3,irot60) * 6 / nr3) == 3) THEN
               IF (ftau(2,imirror)==0 .AND. ftau(3,imirror)==0) THEN
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
!
!   Identify the axis of order 2  and the mirror
!
                irot180=0
                imirror=0
                DO isym=1,nsym
                   IF (tipo_sym(sr(1,1,isym))==4) irot180=isym
                   IF (tipo_sym(sr(1,1,isym))==5) imirror=isym
                ENDDO
                IF ( irot180 == 0 ) CALL errore('find_space_group',&
                                           'C_2h: no 180 degrees rotation',1)
                IF ( imirror == 0 ) CALL errore('find_space_group',&
                                           'C_2h: no mirror found',1)

                IF (ibrav==-12) THEN
!
!   b unique
!
                   IF (ftau(2,irot180) == 0) THEN
                      sg_number=13
                   ELSE
                      IF (ftau(1,imirror)==0.AND.ftau(3,imirror)==0) THEN
                         sg_number=11
                      ELSE
                         sg_number=14
                      ENDIF
                   ENDIF
                ELSE
                   IF (ftau(3,irot180) == 0) THEN
                      sg_number=13
                   ELSE
                      IF (ftau(1,imirror)==0.AND.ftau(2,imirror)==0) THEN
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
         sg_number=174
      CASE(18)
!
!  C_4h
!
         IF (ibrav==6) THEN
            IF (is_symmorphic) THEN
               sg_number=83
            ELSE
               irot90=0
               imirror=0
               DO isym=2,nsym
                  IF (tipo_sym(sr(1,1,isym))==3) THEN
                     IF (ABS(angle_rot(sr(1,1,isym))-90.0_DP)< 1.D-6) &
                                                     irot90=isym
                  ENDIF
                  IF (tipo_sym(sr(1,1,isym))==5) THEN
                     CALL mirror_axis(sr(1,1,isym),ax)
                     IF (is_axis(ax,3)) imirror=isym
                  ENDIF
               ENDDO
               IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                           'C_4h: no 90 degrees rotation',1)
               IF ( imirror == 0 ) CALL errore('find_space_group',&
                                           'C_4h: no mirror found',1)
               IF (ftau(3,irot90)==0) THEN
                  sg_number=85
               ELSEIF (ABS(ftau(3,irot90) * 4 / nr3) == 2) THEN
                  IF (ftau(1,imirror)==0 .AND. ftau(2,imirror)==0) THEN
                     sg_number=84
                  ELSE
                     sg_number=86
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF(ibrav==7) THEN
            IF (is_symmorphic) THEN
               sg_number=87
            ELSE
               sg_number=88
            ENDIF
         ENDIF
      CASE(19)
!
!  C_6h
!
         IF (is_symmorphic) THEN
            sg_number=175
         ELSE
            sg_number=176
         ENDIF
      CASE(20)
!
!  D_2h
!
         nsa=0
         nsaz=0
         iaxis=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==4) THEN
               iaxis=iaxis+1
               CALL versor(sr(1,1,isym),ax)
               DO ipol=1,3
                  IF (is_axis(ax,ipol)) THEN
                     IF (ftau(ipol,isym) /= 0) THEN
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
               ft1(:)=ftau(:,mirr_sym(iaxis))
               ft1(axis_dir(iaxis))=0
               nft=0
               ncomp=0
               DO ipol=1,3
                  nft=nft+ft1(ipol)**2
                  IF (ft1(ipol) /= 0) ncomp=ncomp+1
               END DO  
               IF (nft==0) nmir=nmir+1
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
                     ELSEIF (nmp==2) THEN
                        sg_number=50
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
                     IF (axis_dir(iaxis)==3) THEN
                        ftrs1(:) = ftau(1,mirr_sym(iaxis))*at(:,1) / nr1 + &
                                   ftau(2,mirr_sym(iaxis))*at(:,2) / nr2 + &
                                   ftau(3,mirr_sym(iaxis))*at(:,3) / nr3  
                     ENDIF
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
                        IF (ftau(3,mirr_sym(iaxis))==0) nmir=nmir+1 
                     ENDIF
                  ENDDO
                  nmir=nmir*2

                  IF (nft==0) THEN
                     sg_number=66
                  ELSEIF (nft==1.AND.nmir==0) THEN
                     sg_number=68
                  ELSEIF (nft==1.AND.nmir==2) THEN
                     sg_number=67
                  ENDIF
               ELSEIF (nsaz==1) THEN
!
!   check only the mirror parallel to the centered C face
!
                  DO iaxis=1,3
                     IF (axis_dir(iaxis)==3) THEN
                        ftrs1(:) = ftau(1,mirr_sym(iaxis))*at(:,1) / nr1 + &
                                   ftau(2,mirr_sym(iaxis))*at(:,2) / nr2 + &
                                   ftau(3,mirr_sym(iaxis))*at(:,3) / nr3  
                    
                     ENDIF
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
               sg_number=70
            ENDIF
         ELSEIF (ibrav==11) THEN
            IF (is_symmorphic) THEN
               sg_number=71
            ELSE
               sg_number=72
!
!   check the three mirrors 
!
               DO iaxis=1,3
                  IF (axis_dir(iaxis)==1) THEN
                     ftrs1(:) = ftau(1,mirr_sym(iaxis))*at(:,1) / nr1 + &
                                ftau(2,mirr_sym(iaxis))*at(:,2) / nr2 + &
                                ftau(3,mirr_sym(iaxis))*at(:,3) / nr3
                  ENDIF
                  IF (axis_dir(iaxis)==2) THEN
                     ftrs2(:) = ftau(1,mirr_sym(iaxis))*at(:,1) / nr1 + &
                                ftau(2,mirr_sym(iaxis))*at(:,2) / nr2 + &
                                ftau(3,mirr_sym(iaxis))*at(:,3) / nr3
                  ENDIF
                  IF (axis_dir(iaxis)==3) THEN
                     ftrs3(:) = ftau(1,mirr_sym(iaxis))*at(:,1) / nr1 + &
                                ftau(2,mirr_sym(iaxis))*at(:,2) / nr2 + &
                                ftau(3,mirr_sym(iaxis))*at(:,3) / nr3
                  ENDIF
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
         type_group=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==4) THEN
!
!  search the axis of order 2 and compute the angle with the x axis
!
               CALL versor(sr(1,1,isym), ax)
               angle=ACOS(ax(1))*180.0_DP / pi
               IF (MOD(NINT(angle), 60)==0) type_group=1
               EXIT
            ENDIF
         ENDDO

         IF (is_symmorphic) THEN
            IF (type_group==1) THEN
               sg_number=189
            ELSE
               sg_number=187
            ENDIF
         ELSE
            IF (type_group==1) THEN
               sg_number=190
            ELSE
               sg_number=188
            ENDIF
         ENDIF
      CASE(22)
!
!  D_4h
!
         irot90=0
         irot180=0
         imirror=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==3) THEN
               IF ( ABS(angle_rot(sr(1,1,isym))-90.0_DP )< 1.D-6 ) irot90=isym
            ENDIF
            IF (tipo_sym(sr(1,1,isym))==4) THEN
               CALL versor(sr(1,1,isym),ax)
               IF (is_axis(ax,1)) irot180=isym
            ENDIF
            IF (tipo_sym(sr(1,1,isym))==5) THEN
               CALL mirror_axis(sr(1,1,isym), ax)
               IF ( is_axis(ax,1) ) imirror=isym
            ENDIF
         ENDDO
         IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                        'D_4h: no 90 degrees rotation',1)
         IF ( irot180 == 0 ) CALL errore('find_space_group',&
                                        'D_4h: no 180 degrees rotation',1)
         IF ( imirror == 0 ) CALL errore('find_space_group',&
                                        'D_4h: no mirror found',1)
         IF (ibrav==6) THEN
            IF (is_symmorphic) THEN
               sg_number=123
            ELSE
               IF (ftau(3,irot90)==0) THEN
                  IF (ftau(1,irot180)==0) THEN
                     IF (ftau(2,imirror)==0) THEN
!
!  glide plane c
!
                        sg_number=124
                     ELSEIF(ftau(3,imirror)==0) THEN
!
!  glide plane b
!
                        sg_number=125
                     ELSE
!
!  glide plane n
!
                        sg_number=126
                     ENDIF
                  ELSEIF (ABS(ftau(1,irot180) * 2 / nr1) == 1) THEN
                     IF (ftau(2,imirror)==0.AND.ftau(3,imirror)==0) THEN
!
!  mirror
!
                        sg_number=129
                     ELSEIF (ftau(2,imirror)==0) THEN
!
!  glide plane c
!
                        sg_number=130
                     ELSEIF(ftau(3,imirror)==0) THEN
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
               ELSEIF (ABS(ftau(3,irot90) * 4 / nr3) == 2) THEN
                  IF (ftau(1,irot180)==0) THEN
                     IF (ftau(2,imirror)==0.AND.ftau(3,imirror)==0) THEN
!
!  mirror
!
                        sg_number=131
                     ELSEIF (ftau(2,imirror)==0) THEN
!
!  glide plane c
!
                        sg_number=132
                     ELSEIF(ftau(3,imirror)==0) THEN
!
!  glide plane b
!
                        sg_number=133
                     ELSE
!
!  glide plane n
!
                        sg_number=134
                     ENDIF
                  ELSEIF (ABS(ftau(1,irot180) * 2 / nr1) == 1) THEN
                     IF (ftau(2,imirror)==0.AND.ftau(3,imirror)==0) THEN
!
!  mirror
!
                        sg_number=137
                     ELSEIF (ftau(2,imirror)==0) THEN
!
!  glide plane c
!
                        sg_number=138
                     ELSEIF(ftau(3,imirror)==0) THEN
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
               ftrs1(:) = ftau(1,irot90)*at(:,1) / nr1 + &
                          ftau(2,irot90)*at(:,2) / nr2 + &
                          ftau(3,irot90)*at(:,3) / nr3

               IF ( ABS(ABS(ftrs1(3) / at(3,1))-1.0_DP)<1.d-6.OR. &
                     ABS(ftrs1(3)) < 1.d-6 ) THEN
!
!  the axis of order four is proper
!
                  sg_number=140
               ELSE
!
!  the axis of order four is improprer
!
                  ftrs2(:) = ftau(1,imirror) * at(:,1) / nr1 + &
                             ftau(2,imirror) * at(:,2) / nr2 + &
                             ftau(3,imirror) * at(:,3) / nr3
                  nft=0
                  IF ((ABS(ftrs2(2)) < 1.d-6 .AND. ABS(ftrs2(3)) > 1.d-6).OR.&
                      (ABS(ftrs2(3)) < 1.d-6 .AND. ABS(ftrs2(2)) > 1.d-6) ) &
                      nft=nft+1

                  IF (nft==0) THEN
                     sg_number=141
                  ELSE
                     sg_number=142
                  ENDIF
               ENDIF
            ENDIF
         ENDIF 
      CASE(23)
!
!  D_6h
!
         IF (is_symmorphic) THEN
            sg_number=191
         ELSE
            irot60=0
            imirror=0
            DO isym=2,nsym
               IF (tipo_sym(sr(1,1,isym))==3) THEN
                  IF ( ABS(angle_rot(sr(1,1,isym))-60.0_DP )< 1.D-6 ) &
                                                          irot60=isym
               ENDIF
               IF (tipo_sym(sr(1,1,isym))==5) THEN
                  CALL mirror_axis(sr(1,1,isym), ax)
                  IF ( is_axis(ax,1) ) imirror=isym
               ENDIF
            ENDDO
            IF ( irot60 == 0 ) CALL errore('find_space_group',&
                                           'D_6h: no 60 degrees rotation',1)
            IF ( imirror == 0 ) CALL errore('find_space_group',&
                                           'D_6h: no mirror found',1)

            IF (ftau(3,irot60)==0) THEN
               sg_number=192
            ELSEIF (ABS(ftau(3,irot60) * 6 / nr3) == 3) THEN
               IF (ftau(2,imirror)==0.AND.ftau(3,imirror)==0) THEN
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
         type_group=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==4) THEN
!
!  search the axis of order 2 and compute the angle with the x axis,
!  it can be 0,90,180,270 (type_group=1) or 45,135,225,315 (type_group=0)
!
               CALL versor(sr(1,1,isym), ax)
               IF (ABS(ax(3)) > 1.D-6) CYCLE
               angle=ACOS(ax(1))*180.0_DP / pi
               IF (MOD(NINT(angle), 90)==0) type_group=1
               EXIT
            ENDIF
         ENDDO

         IF (type_group==1) THEN
!
!  the x axis is a two-fold rotation axis and there is a mirror 
!  perpendicular to (1,1,0). Search it
!
            irot180=0
            imirror1=0
            DO isym=2,nsym
               IF (tipo_sym(sr(1,1,isym))==4) THEN
                  CALL versor(sr(1,1,isym), ax)
                  IF ( is_axis(ax,1) ) irot180=isym
               ENDIF
               IF (tipo_sym(sr(1,1,isym))==5) THEN
                  CALL mirror_axis(sr(1,1,isym), ax)
                  sp=ax(1)+ax(2)
                  ft1(1)=ax(1)-sp * 0.5_DP
                  ft1(2)=ax(2)-sp * 0.5_DP
                  ft1(3)=0.0_DP
                  IF (ft1(1)**2+ft1(2)**2<1.d-6) imirror1=isym
               ENDIF
            ENDDO
            IF ( irot180 == 0 ) CALL errore('find_space_group',&
                                           'D_2d: no 180 degrees rotation',1)
            IF ( imirror1 == 0 ) CALL errore('find_space_group',&
                                           'D_2d: no mirror found',1)
         ELSE
!
!  There is a mirror perpendicular to the x axis. 
!
            imirror=0
            DO isym=2,nsym
               IF (tipo_sym(sr(1,1,isym))==5) THEN
                  CALL mirror_axis(sr(1,1,isym), ax)
                  IF ( is_axis(ax,1) ) imirror=isym
               ENDIF
            ENDDO
            IF ( imirror == 0 ) CALL errore('find_space_group',&
                                           'D_2d: no mirror found',1)
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
                  IF (ftau(1,irot180)==0) THEN
                     sg_number=112
                  ELSE
                    IF (ftau(3,imirror1)==0) THEN
                       sg_number=113
                    ELSE
                       sg_number=114
                    ENDIF
                  ENDIF
               ELSE
                  IF (ftau(3,imirror)==0) THEN
                     sg_number=117
                  ELSEIF (ftau(2,imirror)==0) THEN
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
         type_group=0
         DO isym=2,nsym
            IF (tipo_sym(sr(1,1,isym))==4) THEN
!
!  search the axis of order 2 and compute the angle with the x axis,
!  it can be 0,60,120 (type_group=1) or 30,90,150 (type_group=0)
!
               CALL versor(sr(1,1,isym), ax)
               angle=ACOS(ax(1))*180.0_DP / pi
               IF (MOD(NINT(angle), 60)==0) type_group=1
               EXIT
            ENDIF
         ENDDO

         IF (ibrav==4) THEN
             IF (is_symmorphic) THEN
                IF (type_group==1) THEN
                   sg_number=164
                ELSE
                   sg_number=162
                ENDIF
             ELSE
                IF (type_group==1) THEN
                   sg_number=165
                ELSE
                   sg_number=163
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
         IF (ibrav==6) THEN
            sg_number=81
         ELSEIF (ibrav==7) THEN
            sg_number=82
         ENDIF
      CASE(27)
!
!  S_6 (-3)
!
         IF (ibrav==4) THEN
            sg_number=147
         ELSEIF (ibrav==5) THEN
            sg_number=148
         ENDIF
      CASE(28)
!
!   T
!
         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=195
             ELSE
                sg_number=198
             ENDIF
         ELSEIF (ibrav==2) THEN
             sg_number=196
         ELSEIF (ibrav==3) THEN
             IF (is_symmorphic) THEN
                sg_number=197
             ELSE
                sg_number=199
             ENDIF
         ENDIF
      CASE(29)
!
!  T_h
!

         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=200
             ELSE
                irot180=0
                DO isym=2,nsym
                   IF (tipo_sym(sr(1,1,isym))==4) THEN
                      CALL versor(sr(1,1,isym), ax)
                      IF (is_axis(ax,3)) irot180=isym
                   ENDIF
                ENDDO
                IF ( irot180 == 0 ) CALL errore('find_space_group',&
                                              'T_h: no 180 degrees rotation',1)
                IF (ftau(3,irot180)==0) THEN
                   sg_number=201
                ELSE
                   sg_number=205
                ENDIF
             ENDIF
         ELSEIF (ibrav==2) THEN
             IF (is_symmorphic) THEN
                sg_number=202
             ELSE
                sg_number=203
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
         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=215
             ELSE
                sg_number=218
             ENDIF
         ELSEIF (ibrav==2) THEN
             IF (is_symmorphic) THEN
                sg_number=216
             ELSE
                sg_number=219
             ENDIF
         ELSEIF (ibrav==3) THEN
             IF (is_symmorphic) THEN
                sg_number=217
             ELSE
                sg_number=220
             ENDIF
         ENDIF
      CASE(31)
!
!  O
!
         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=207
             ELSE
!
!          find the 90 degree rotation about the x axis and check its
!          fractional translation along the axis.
!
                irot90=0
                DO isym=2,nsym
                   IF (tipo_sym(sr(1,1,isym))==3) THEN
                      CALL versor(sr(1,1,isym), ax)
                      IF (ABS(angle_rot(sr(1,1,isym))-90.0_DP)< 1.D-6 .AND. &
                         is_axis(ax,1) ) irot90=isym
                   ENDIF
                ENDDO
                IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                              'O: no 90 degrees rotation',1)
                IF (ftau(1,irot90) * 4 / nr1 == -1 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == 3) sg_number=213
                IF (ABS(ftau(1,irot90) * 4 / nr1) == 2) sg_number=208
                IF (ftau(1,irot90) * 4 / nr1 == 1 .OR.  &
                           ftau(3,irot60) * 6 / nr3 == -3) sg_number=212
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
                sg_number=214
             ENDIF
         ENDIF
      CASE(32)
!
!  O_h
!
         IF (ibrav==1) THEN
             IF (is_symmorphic) THEN
                sg_number=221
             ELSE
!
!   find the rotation of 90 degree about the xaxis and the mirror 
!   perpendicular to the x axis
!
               irot90=0
               imirror=0
               DO isym=2,nsym
                   IF (tipo_sym(sr(1,1,isym))==3) THEN
                      CALL versor(sr(1,1,isym), ax)
                      IF (ABS(angle_rot(sr(1,1,isym))-90.0_DP)< 1.D-6 .AND. &
                         is_axis(ax,1) ) irot90=isym
                   ENDIF
                   IF (tipo_sym(sr(1,1,isym))==5) THEN
                      CALL mirror_axis(sr(1,1,isym), ax)
                      IF (is_axis(ax,1)) imirror=isym
                   ENDIF
               ENDDO
               IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                              'O_h: no 90 degrees rotation',1)
               IF ( imirror == 0 ) CALL errore('find_space_group',&
                                              'O_h: no mirror found',1)
!
!  First check if the axis of order 4 has fractional translation, if not
!  the group is 222, if yes we check the fractional translation of the
!  mirror perpendicular to the x axis. In 224 this mirror is a glide plane,
!  in 223 no.
!
               IF (irot90 /= 0) THEN
                  IF (ftau(1,irot90) == 0) THEN
                     sg_number=222
                  ELSEIF (ABS(ftau(1,irot90) * 4 / nr1) == 2) THEN
                     IF (ftau(2,imirror)==0 .AND. ftau(2,imirror)==0 ) THEN
                        sg_number=223
                     ELSE
                        sg_number=224
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ibrav==2) THEN
            IF (is_symmorphic) THEN
               sg_number=225
            ELSE
!
!   find the rotation of 90 degree about the xaxis and the mirror 
!   perpendicular to the (110) axis
!
               irot90=0
               imirror=0
               DO isym=2,nsym
                  IF (tipo_sym(sr(1,1,isym))==3) THEN
                     CALL versor(sr(1,1,isym), ax)
                     IF (ABS(angle_rot(sr(1,1,isym))-90.0_DP)< 1.D-6 .AND. &
                        is_axis(ax,1) ) irot90=isym
                  ENDIF
                  IF (tipo_sym(sr(1,1,isym))==5) THEN
                     CALL mirror_axis(sr(1,1,isym), ax)
                     IF (ABS(ax(1)-ax(2))<1.d-6) imirror=isym
                  ENDIF
               ENDDO
               IF ( irot90 == 0 ) CALL errore('find_space_group',&
                                              'O_h: no 90 degrees rotation',1)
               IF ( imirror == 0 ) CALL errore('find_space_group',&
                                              'O_h: no mirror found',1)

!
!  First check if the axis of order 4 has fractional translation, if not
!  the group is 226, if yes we check the fractional translation of the
!  mirror perpendicular to the axis of order 2 (x=y). 
!  In 228 this mirror is a glide plane, in 227 no.
!
               ftrs1(:) = ftau(1,irot90) * at(:,1) / nr1 &
                        + ftau(2,irot90) * at(:,2) / nr2 & 
                        + ftau(3,irot90) * at(:,3) / nr3
               ftrs2(:) = ftau(1,imirror) * at(:,1) / nr1 &
                        + ftau(2,imirror) * at(:,2) / nr2 & 
                        + ftau(3,imirror) * at(:,3) / nr3

               IF ( ABS(ABS(ftrs1(1))-0.5_DP)<1.d-6 .OR. ABS(ftrs1(1)) &
                                                          < 1.d-6 ) THEN
                  sg_number=226
               ELSE
!
!  put the fractional translation of the mirror in cartesian coordinates
!
                  fact= ftrs2(1) + ftrs2(2)  
!
!   project on the mirror plane
!

                  ftrs2(1)=ftrs2(1) - fact / 2.0_DP
                  ftrs2(2)=ftrs2(2) - fact / 2.0_DP
                  sp = ftrs2(1)**2 + ftrs2(2)**2 + ftrs2(3)**2 
!
!  if the fractional translation vanishes it is 227, ortherwise 228
!
                  IF (sp<1.d-6) THEN
                     sg_number=227
                  ELSE
                     sg_number=228
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

  INTEGER FUNCTION tipo_sym_sg(sr, ftau, at, nr1, nr2, nr3, ashift, acryshift)
!
!   This subroutine receives as input a rotation matrix in cartesian coordinates,
!   and a fractinal translation as three integer values on an fft mesh, the direct
!   and reciprocal lattice vectors and the fft mesh and gives as output a code 
!   with the symmetry type. 1-16 are the operations of symmorphic space group,
!   non symmorphic groups have at least one of the operation 17-41.
!   
!   1   identity
!   2   inversion
!   3   mirror
!   4   rotation of 60  degrees 6
!   5   rotation of 90  degrees 4
!   6   rotation of 120 degrees 3
!   7   rotation of 180 degrees 2
!   8   rotation of 240 degrees 3^2
!   9   rotation of 270 degrees 4^3
!  10   rotation of 300 degrees 6^5
!  11   rotation of 60  degrees multiplied by inversion 6
!  12   rotation of 90  degrees multiplied by inversion 4
!  13   rotation of 120 degrees multiplied by inversion 3
!  14   rotation of 240 degrees multiplied by inversion 3^2
!  15   rotation of 270 degrees multiplied by inversion 4^3
!  16   rotation of 300 degrees multiplied by inversion 6^5
!
!  Non symmorphic space groups have some operations of the following types:
!
!  17  screw rotation 180 degrees C2 (2)     (axis 2_1) 
!  18  screw rotation 120 degrees C3 (3)     (axis 3_1) 
!  19  screw rotation 240 degrees C3^2 (3^2) (axis 3_1) 
!  20  screw rotation 120 degrees C3 (3)     (axis 3_2) 
!  21  screw rotation 240 degrees C3^2 (3)   (axis 3_2) 
!  22  screw rotation  90 degrees C4 (4)     (axis 4_1) 
!  23  screw rotation 270 degrees C4^3 (4^3) (axis 4_1) 
!  24  screw rotation  90 degrees C4 (4)     (axis 4_2) 
!  25  screw rotation 270 degrees C4^3 (4^3) (axis 4_2) 
!  26  screw rotation  90 degrees C4 (4)     (axis 4_3) 
!  27  screw rotation 270 degrees C4^3 (4^3) (axis 4_3) 
!  28  screw rotation  60 degrees C6 (6)     (axis 6_1) 
!  29  screw rotation 300 degrees C6^5 (6^5) (axis 6_1) 
!  30  screw rotation  60 degrees C6 (6)     (axis 6_2) 
!  31  screw rotation 300 degrees C6^5 (6^5) (axis 6_2) 
!  32  screw rotation  60 degrees C6 (6)     (axis 6_3) 
!  33  screw rotation 300 degrees C6^5 (6^5) (axis 6_3) 
!  34  screw rotation  60 degrees C6 (6)     (axis 6_4) 
!  35  screw rotation 300 degrees C6^5 (6^5) (axis 6_4) 
!  36  screw rotation  60 degrees C6 (6)     (axis 6_5) 
!  37  screw rotation 300 degrees C6^5 (6^5) (axis 6_5) 
!  38  glide plane a
!  39  glide plane b
!  40  glide plane c
!  41  glide plane 
!
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: sr(3,3), at(3,3)
  REAL(DP), INTENT(OUT) :: ashift(3), acryshift(3)
  INTEGER, INTENT(IN) :: ftau(3), nr1, nr2, nr3

  REAL(DP) :: fcart(3), angle, ps, fmod, ax(3), fcrys(3), bg(3,3)
  REAL(DP), PARAMETER :: f43=4.0_DP / 3.0_DP, eps=1.D-7
  REAL(DP) :: angle_rot, angle_rot_s

  INTEGER :: ts, tip_sym
  INTEGER :: tipo_sym
  
  CALL recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
  ts=tipo_sym(sr)
  fcart(:) = ftau(1) * at(:,1) / nr1 +   &
             ftau(2) * at(:,2) / nr2 +   &
             ftau(3) * at(:,3) / nr3 
  ashift=0.0_DP

  SELECT CASE (ts)
     CASE(1)
!
!    identity
!
         tipo_sym_sg=1
     CASE(2)
!
!    inversion
!
         tipo_sym_sg=2
         fmod=SQRT(fcart(1)**2 + fcart(2)**2 + fcart(3)**3)
         ashift(:) = fcart(:) * 0.5_DP
     CASE(3)
!
!    proper rotation not 180
!
        angle=angle_rot(sr)
        CALL versor(sr, ax)
        ps = ax(1) * fcart(1) + ax(2) * fcart(2) + ax(3) * fcart(3) 
        IF (ABS(ps) < eps) THEN
           IF (ABS(angle-60.0_DP) < eps) THEN
              tipo_sym_sg=4
           ELSEIF (ABS(angle-90.0_DP) < eps) THEN
              tipo_sym_sg=5
           ELSEIF (ABS(angle-120.0_DP) < eps) THEN
              tipo_sym_sg=6
           ELSEIF (ABS(angle-240.0_DP) < eps) THEN
              tipo_sym_sg=8
           ELSEIF (ABS(angle-270.0_DP) < eps) THEN
              tipo_sym_sg=9
           ELSEIF (ABS(angle-300.0_DP) < eps) THEN
              tipo_sym_sg=10
           ELSE
              CALL errore('tipo_sym_sg','angle not recognized',1)
           ENDIF 
        ELSE
           IF (ABS(angle-60.0_DP) < eps) THEN
              fcart =  ps * ax
              fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
              IF (ABS(1.2_DP * fcrys(1) - NINT(1.2_DP * fcrys(1)) ) < eps .AND. &
                  ABS(1.2_DP * fcrys(2) - NINT(1.2_DP * fcrys(2)) ) < eps .AND. &
                  ABS(1.2_DP * fcrys(3) - NINT(1.2_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=36
              ELSEIF (ABS(1.5_DP * fcrys(1) - NINT(1.5_DP * fcrys(1)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(2) - NINT(1.5_DP * fcrys(2)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(3) - NINT(1.5_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=34
              ELSEIF (ABS(2.0_DP * fcrys(1) - NINT(2.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(2) - NINT(2.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(3) - NINT(2.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=32
              ELSEIF (ABS(3.0_DP * fcrys(1) - NINT(3.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(2) - NINT(3.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(3) - NINT(3.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=30
              ELSEIF (ABS(6.0_DP * fcrys(1) - NINT(6.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(6.0_DP * fcrys(2) - NINT(6.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(6.0_DP * fcrys(3) - NINT(6.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=28
              ELSE
                 CALL errore('tipo_sym_sg','fractional translation wrong',1)
              END IF
           ELSEIF (ABS(angle-90.0_DP) < eps) THEN
              fcart =  ps * ax
              fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
              IF (ABS(f43 * fcrys(1) - NINT(f43 * fcrys(1)) ) < eps .AND. &
                 ABS(f43 * fcrys(2) - NINT(f43 * fcrys(2)) ) < eps .AND. &
                 ABS(f43 * fcrys(3) - NINT(f43 * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=26
              ELSEIF (ABS(2.0_DP * fcrys(1) - NINT(2.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(2) - NINT(2.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(3) - NINT(2.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=24
              ELSEIF (ABS(4.0_DP * fcrys(1) - NINT(4.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(4.0_DP * fcrys(2) - NINT(4.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(4.0_DP * fcrys(3) - NINT(4.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=22
              ELSE
                 CALL errore('tipo_sym_sg','fractional translation wrong',2)
              END IF
           ELSEIF (ABS(angle-120.0_DP) < eps) THEN
              fcart =  ps * ax
              fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
              IF (ABS(1.5_DP * fcrys(1) - NINT(1.5_DP * fcrys(1)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(2) - NINT(1.5_DP * fcrys(2)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(3) - NINT(1.5_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=20
              ELSEIF (ABS(3.0_DP * fcrys(1) - NINT(3.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(2) - NINT(3.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(3) - NINT(3.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=18
              ELSE
                 CALL errore('tipo_sym_sg','fractional translation wrong',3)
              END IF
           ELSEIF (ABS(angle-240.0_DP) < eps) THEN
              fcart =  ps * ax
              fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
              IF (ABS(1.5_DP * fcrys(1) - NINT(1.5_DP * fcrys(1)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(2) - NINT(1.5_DP * fcrys(2)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(3) - NINT(1.5_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=21
              ELSEIF (ABS(3.0_DP * fcrys(1) - NINT(3.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(2) - NINT(3.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(3) - NINT(3.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=19
              ELSE
                 CALL errore('tipo_sym_sg','fractional translation wrong',4)
              END IF
           ELSEIF (ABS(angle-270.0_DP) < eps) THEN
              fcart =  ps * ax
              fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
              IF (ABS(f43 * fcrys(1) - NINT(f43 * fcrys(1)) ) < eps .AND. &
                 ABS(f43 * fcrys(2) - NINT(f43 * fcrys(2)) ) < eps .AND. &
                 ABS(f43 * fcrys(3) - NINT(f43 * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=27
              ELSEIF (ABS(2.0_DP * fcrys(1) - NINT(2.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(2) - NINT(2.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(3) - NINT(2.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=25
              ELSEIF (ABS(4.0_DP * fcrys(1) - NINT(4.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(4.0_DP * fcrys(2) - NINT(4.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(4.0_DP * fcrys(3) - NINT(4.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=23
              ELSE
                 CALL errore('tipo_sym_sg','fractional translation wrong',5)
              END IF
           ELSEIF (ABS(angle-300.0_DP) < eps) THEN
              fcart =  ps * ax
              fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
              IF (ABS(1.2_DP * fcrys(1) - NINT(1.2_DP * fcrys(1)) ) < eps .AND. &
                 ABS(1.2_DP * fcrys(2) - NINT(1.2_DP * fcrys(2)) ) < eps .AND. &
                 ABS(1.2_DP * fcrys(3) - NINT(1.2_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=37
              ELSEIF (ABS(1.5_DP * fcrys(1) - NINT(1.5_DP * fcrys(1)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(2) - NINT(1.5_DP * fcrys(2)) ) < eps .AND. &
                 ABS(1.5_DP * fcrys(3) - NINT(1.5_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=35
              ELSEIF (ABS(2.0_DP * fcrys(1) - NINT(2.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(2) - NINT(2.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(2.0_DP * fcrys(3) - NINT(2.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=33
              ELSEIF (ABS(3.0_DP * fcrys(1) - NINT(3.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(2) - NINT(3.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(3.0_DP * fcrys(3) - NINT(3.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=31
              ELSEIF (ABS(6.0_DP * fcrys(1) - NINT(6.0_DP * fcrys(1)) ) < eps .AND. &
                 ABS(6.0_DP * fcrys(2) - NINT(6.0_DP * fcrys(2)) ) < eps .AND. &
                 ABS(6.0_DP * fcrys(3) - NINT(6.0_DP * fcrys(3)) ) < eps ) THEN
                 tipo_sym_sg=29
              ELSE
                 CALL errore('tipo_sym_sg','fractional translation wrong',6)
              END IF
           ELSE
              CALL errore('tipo_sym_sg','rotation angle wrong',1)
           END IF
        END IF
     CASE(4)
!
!    proper rotation 180
!
         CALL versor(sr, ax)
         ps = ax(1) * fcart(1) + ax(2) * fcart(2) + ax(3) * fcart(3) 
         IF (ABS(ps) < eps) THEN
            tipo_sym_sg=7
            ashift = 0.5_DP*fcart 
         ELSE
            ashift = 0.5_DP*(fcart - ps*ax)
            fcart(:)=2.0_DP*ps*ax(:)
            fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
            IF ( ABS( fcrys(1) - NINT(fcrys(1)) ) < eps .AND. &
                 ABS( fcrys(2) - NINT(fcrys(2)) ) < eps .AND. &
                 ABS( fcrys(3) - NINT(fcrys(3)) ) < eps ) THEN
                tipo_sym_sg=17
            ELSE
               CALL errore('tipo_sym_sg','fractional translation wrong',7)
            ENDIF
         ENDIF
     CASE(5)
!
!    mirror
!
        CALL mirror_axis(sr, ax)
        ps = ax(1) * fcart(1) + ax(2) * fcart(2) + ax(3) * fcart(3) 
        ashift(:)= ps * ax(:) * 0.5_DP
!
!    fcart is projected in the direction parallel to the mirror
!
        fcart(:)=fcart(:)-ps*ax(:)
        fmod=SQRT(fcart(1)**2 + fcart(2)**2 + fcart(3)**2)
        IF (ABS(fmod) < eps) THEN
           tipo_sym_sg=3
        ELSE
           fcart(:) = 2.0_DP * fcart(:)
           fcrys(:)= fcart(1)*bg(1,:) + fcart(2)*bg(2,:) + fcart(3)*bg(3,:)
          WRITE(6,*) 'fcart', fcart(:)
          WRITE(6,*) 'fcrys', fcrys(:)
          IF  ( ABS(fcrys(1)-NINT(fcrys(1))) < eps .AND. &
                ABS(fcrys(2)) < eps .AND. ABS(fcrys(3)) < eps ) THEN
              tipo_sym_sg=38
          ELSEIF ( ABS(fcrys(2)-NINT(fcrys(2))) < eps .AND. &
                 ABS(fcrys(1)) < eps .AND. ABS(fcrys(3)) < eps ) THEN
              tipo_sym_sg=39
          ELSEIF ( ABS(fcrys(3)-NINT(fcrys(3))) < eps .AND. &
                 ABS(fcrys(1)) < eps .AND. ABS(fcrys(2)) < eps ) THEN
              tipo_sym_sg=40
          ELSEIF ((ABS(fcrys(1)-NINT(fcrys(1)))<eps.AND. &
                   ABS(fcrys(2)-NINT(fcrys(2)))<eps.AND.ABS(fcrys(3))<eps).OR. &
                  (ABS(fcrys(1)-NINT(fcrys(1)))<eps.AND. &
                   ABS(fcrys(3)-NINT(fcrys(3)))<eps.AND.ABS(fcrys(2)) < eps ).OR. &
                  (ABS(fcrys(2)-NINT(fcrys(2)))<eps.AND. &
                   ABS(fcrys(3)-NINT(fcrys(3)))<eps.AND.ABS(fcrys(1)) < eps ).OR. &
                  (ABS(fcrys(1)-NINT(fcrys(1)))<eps.AND. & 
                   ABS(fcrys(2)-NINT(fcrys(2)))<eps.AND. &
                   ABS(fcrys(3)-NINT(fcrys(3)))<eps) ) THEN
               tipo_sym_sg=41
           ELSEIF ((ABS(2.0_DP*fcrys(1)-NINT(2.0_DP*fcrys(1)))<eps.AND. &
                    ABS(2.0_DP*fcrys(2)-NINT(2.0_DP*fcrys(2)))<eps.AND. &
                    ABS(fcrys(3))<eps).OR. &
                   (ABS(2.0_DP*fcrys(1)-NINT(2.0_DP*fcrys(1)))<eps.AND. &
                    ABS(2.0_DP*fcrys(3)-NINT(2.0_DP*fcrys(3)))<eps.AND. &
                    ABS(fcrys(2)) < eps ).OR. &
                   (ABS(2.0_DP*fcrys(2)-NINT(2.0_DP*fcrys(2)))<eps.AND. &
                    ABS(2.0_DP*fcrys(3)-NINT(2.0_DP*fcrys(3)))<eps.AND. &
                    ABS(fcrys(1)) < eps ) .OR. &
                   (ABS(2.0_DP*fcrys(1)-NINT(2.0_DP*fcrys(1)))<eps.AND. &
                    ABS(2.0_DP*fcrys(2)-NINT(2.0_DP*fcrys(2)))<eps.AND. &
                    ABS(2.0_DP*fcrys(3)-NINT(2.0_DP*fcrys(3)))<eps )) THEN
              tipo_sym_sg=42
           ELSE
              CALL errore('tipo_sym_sg','fractional translation wrong',1)
           ENDIF
        ENDIF
     CASE(6)
!
!   here the fractionary translation is not relevant to identify the space
!   group, we give only the type of the rotation
!
        angle=angle_rot(sr)
        IF (ABS(angle-60.0_DP) < eps) THEN
           tipo_sym_sg=11
        ELSEIF (ABS(angle-90.0_DP) < eps) THEN
           tipo_sym_sg=12
        ELSEIF (ABS(angle-120.0_DP) < eps) THEN
           tipo_sym_sg=13
        ELSEIF (ABS(angle-240.0_DP) < eps) THEN
           tipo_sym_sg=14
        ELSEIF (ABS(angle-270.0_DP) < eps) THEN
           tipo_sym_sg=15
        ELSEIF (ABS(angle-300.0_DP) < eps) THEN
           tipo_sym_sg=16
        ELSE
           CALL errore('tipo_sym_sg','angle not recognized',1)
        ENDIF
     CASE DEFAULT 
          CALL errore('tipo_sym_sg',' symmetry type not available',1)
     END SELECT
     acryshift(:)= ashift(1)*bg(1,:) + ashift(2)*bg(2,:) + ashift(3)*bg(3,:)
  RETURN
  END FUNCTION tipo_sym_sg 

  CHARACTER(LEN=65) FUNCTION label_sg(ts)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ts
  CHARACTER(LEN=65) :: label(41)

  IF (ts < 1 .OR. ts > 41) CALL errore('label_sg','unknown symmetry type',1 )

  label(1)='identity  E'
  label(2)='inversion i'
  label(3)='mirror    s'
  label(4)='rotation of 60  degrees    C6    (6)'
  label(5)='rotation of 90  degrees    C4    (4)'
  label(6)='rotation of 120 degrees    C3    (3)'
  label(7)='rotation of 180 degrees    C2    (2)'
  label(8)='rotation of 240 degrees    C3^2 (3^2)'
  label(9)='rotation of 270 degrees    C4^3 (4^3)'
  label(10)='rotation of 300 degrees    C6^5 (6^5)'
  label(11)='rotation of 60  degrees multiplied by inversion S3^5 (-6)'
  label(12)='rotation of 90  degrees multiplied by inversion S4^3 (-4)'
  label(13)='rotation of 120 degrees multiplied by inversion S6^5 (-3)'
  label(14)='rotation of 240 degrees multiplied by inversion S6  (-3^2)'
  label(15)='rotation of 270 degrees multiplied by inversion S4  (-4^3)'
  label(16)='rotation of 300 degrees multiplied by inversion S3  (-6^5)'
  label(17)='screw rotation 180 degrees C2    (2)  (axis 2_1)'
  label(18)='screw rotation 120 degrees C3    (3)  (axis 3_1)'
  label(19)='screw rotation 240 degrees C3^2 (3^2) (axis 3_1)'
  label(20)='screw rotation 120 degrees C3    (3)  (axis 3_2)'
  label(21)='screw rotation 240 degrees C3^2  (3)  (axis 3_2)'
  label(22)='screw rotation  90 degrees C4    (4)  (axis 4_1)'
  label(23)='screw rotation 270 degrees C4^3 (4^3) (axis 4_1)'
  label(24)='screw rotation  90 degrees C4    (4)  (axis 4_2)'
  label(25)='screw rotation 270 degrees C4^3 (4^3) (axis 4_2)'
  label(26)='screw rotation  90 degrees C4    (4)  (axis 4_3)'
  label(27)='screw rotation 270 degrees C4^3 (4^3) (axis 4_3)'
  label(28)='screw rotation  60 degrees C6    (6)  (axis 6_1)'
  label(29)='screw rotation 300 degrees C6^5 (6^5) (axis 6_1)'
  label(30)='screw rotation  60 degrees C6    (6)  (axis 6_2)'
  label(31)='screw rotation 300 degrees C6^5 (6^5) (axis 6_2)'
  label(32)='screw rotation  60 degrees C6    (6)  (axis 6_3)'
  label(33)='screw rotation 300 degrees C6^5 (6^5) (axis 6_3)'
  label(34)='screw rotation  60 degrees C6    (6)  (axis 6_4)'
  label(35)='screw rotation 300 degrees C6^5 (6^5) (axis 6_4)'
  label(36)='screw rotation  60 degrees C6    (6)  (axis 6_5)'
  label(37)='screw rotation 300 degrees C6^5 (6^5) (axis 6_5)'
  label(38)='glide plane a'
  label(39)='glide plane b'
  label(40)='glide plane c'
  label(41)='glide plane'
   
  label_sg=label(ts)

  RETURN
  END FUNCTION label_sg

!
END MODULE space_groups
