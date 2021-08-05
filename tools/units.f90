!
! Copyright (C) 2019 Andrea Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM units
!
!  This program defines the units of measure (UOM) of the physical
!  quantities in atomic units in terms of the SI units. It
!  starts from a few values of the basic quantities.
!
USE kinds,            ONLY : DP
USE mp_global,        ONLY : mp_startup, mp_global_end
USE environment,      ONLY : environment_start, environment_end
USE io_global,        ONLY : stdout, ionode
IMPLICIT NONE
INTEGER :: iunout, ios
INTEGER :: find_free_unit
REAL(DP) :: pi, zero, sqrt2

REAL(DP) :: hplanck, hbarf, cspeed, e, rydberg, alphaf, amu, me, mu0,      &
            epsilon0, abohr, bohrmag, ryev, avonum, boltz, rgas

REAL(DP) :: barl, barm, bart, barnu, barv, bara, barp, baram, barf, baru,  &
            barw, barpr, bari, barc, barrho, barcur, bare, barphi, barcap, &
            bardip, barpolar, bard, barohm, barb, barav, barwb, bary,      &
            barmu, barmag, barh, barrhom

REAL(DP) :: barifc, bardmc, baralpha, baralphap

REAL(DP) :: cspeedau, amuau, kappa, kappaa, kappadiecic, toverg, cmm1hz,    &
            hzcmm1

REAL(DP) :: cmtom, gtokg, ptop, ltol, ftof, utou, prtopr, chtoch, itoi,     &
            rhotorho, curtocur, etoe, phitophi, captocap, diptodip,         &
            polartopolar, dtod, ohmtoohm, btob, avtoav, wbtowb, ytoy,       &
            mutomu, magtomag, htoh, rhomtorhom

REAL(DP) :: alphatoalpha, alphaptoalphap, zmtozm

REAL(DP) :: barlcgs, barmcgs, barvcgs, baracgs, barpcgs, baramcgs, barfcgs, &
            barucgs, barwcgs, barprcgs, baricgs, barccgs, barrhocgs,        &
            barcurcgs, barecgs, barphicgs, barcapcgs, bardipcgs,            &
            barpolarcgs, bardcgs, barohmcgs, barbcgs, baravcgs, barwbcgs,   &
            barycgs, barmucgs, barmagcgs, barhcgs, barrhomcgs

REAL(DP) :: barmry, bartry, barcry, barnury, barvry, barary, barfry,   &
            barury, barwry, barprry, bariry, barrhory, barcurry, barery, &
            barphiry, bardipry, barpolarry, bardry, barohmry, barbry,  &
            baravry, barwbry, baryry, barmury, barmagry, barhry

REAL(DP) :: barbg, baravg, barwbg, barmug, barmagg, barhg

REAL(DP) :: rydbergerr, alphaerr, amuerr, meerr, mu0err, epsilon0err,        &
            abohrerr, hartreeerr, bohrmagerr, barterr, barnuerr, barverr,    &
            baraerr, barperr, barferr, barwerr, barprerr, barierr,           &
            barrhoerr, barcurerr, bareerr, barphierr, barcaperr, bardiperr,  &
            barpolarerr, barderr, barberr, baraverr, baryerr, barmagerr,     &
            barherr, itoierr, rhotorhoerr, curtocurerr, etoeerr,             &
            phitophierr, captocaperr, diptodiperr, polartopolarerr, dtoderr, & 
            ohmtoohmerr, kappaerr, btoberr, avtoaverr, wbtowberr, ytoyerr,   &
            mutomuerr, magtomagerr, htoherr, baricgserr, barccgserr,         &
            barrhocgserr, barcurcgserr, barecgserr, barphicgserr,            &
            barcapcgserr, bardipcgserr, barpolarcgserr, bardcgserr,          &
            barohmcgserr, barbcgserr, baravcgserr, barwbcgserr, barmuerr,    &
            barycgserr, barmucgserr, barmagcgserr, barhcgserr, cspeedauerr,  &
            amuauerr, kappaaerr, barmryerr, bartryerr, barnuryerr, barvryerr,&
            bararyerr, barfryerr, baruryerr, barwryerr, barprryerr,          &
            bariryerr, barrhoryerr, barcurryerr, bareryerr, barphiryerr,     &
            bardipryerr, barpolarryerr, bardryerr, barbryerr,                &
            baravryerr, baryryerr, barmuryerr, barmagryerr, barhryerr,       &
            barbgerr, baravgerr, barwbgerr, barmugerr, barmaggerr, barhgerr, &
            barrhomerr

CHARACTER(LEN=80) :: float_to_latex

CHARACTER(LEN=9) :: code='units'

CALL mp_startup ( start_images=.true. )
CALL environment_start ( code )

IF (ionode) THEN
   iunout=find_free_unit()
   OPEN(UNIT=iunout, FILE='units_values.tex', STATUS='unknown', &
                                            FORM='formatted', IOSTAT=ios)
ENDIF

zero=0.0_DP
pi=4.0_DP * ATAN(1.0_DP)
sqrt2=SQRT(2.0_DP)
!
!  defined constants
!
hplanck=6.62607015D-34  ! J.s       exact
hbarf= hplanck / 2.0_DP / pi  ! J.s exact
cspeed=2.99792458D8 ! m/s           exact
e=1.602176634D-19 ! C               exact
avonum=6.02214076D23 ! N_A          exact
boltz=1.380649D-23 ! k_B          exact
rgas=avonum*boltz    ! R            exact

rydberg=1.0973731568160D7 ! rydberg
alphaf=7.2973525693D-3  ! fine strcuture constants
amu=1.66053906660D-27   ! kg

WRITE(stdout,'(/,5x,"Experimental quantities exact in the SI:")') 
WRITE(stdout,'(5x,"Planck constant: ",4x,es20.8,"      J.s")') hplanck
WRITE(stdout,'(5x,"Planck constant / 2 pi:",3x,es20.13," J.s")') hbarf
WRITE(stdout,'(5x,"Speed of light: ",5x,es20.8,"      m/s")') cspeed
WRITE(stdout,'(5x,"Electron charge: ",5x,es20.9,"     C")') e
WRITE(stdout,'(5x,"Avogadro number: ",4x,es20.8)') avonum
WRITE(stdout,'(5x,"Boltzmann constant: ",es19.6,"        J/K")') boltz

WRITE(stdout,'(/,5x,"Approximate quantities determined by experiment:")') 
WRITE(stdout,'(5x,"Rydberg constant: ",8x,es20.13," 1/m")') rydberg
WRITE(stdout,'(5x,"Fine structure constant: ",es18.10)') alphaf
WRITE(stdout,'(5x,"Atomic mass unit: ",6x,es20.11,"   kg")') amu

IF (ionode) THEN
   WRITE(iunout,'("\def\hplanck{",a,"}")') TRIM(float_to_latex(hplanck,8))
   WRITE(iunout,'("\def\hbarf{",a,"}")') TRIM(float_to_latex(hbarf,13))
   WRITE(iunout,'("\def\cspeed{",a,"}")') TRIM(float_to_latex(cspeed,8))
   WRITE(iunout,'("\def\e{",a,"}")') TRIM(float_to_latex(e,9))
   WRITE(iunout,'("\def\rydberg{",a,"}")') TRIM(float_to_latex(rydberg,13))
   WRITE(iunout,'("\def\alphaf{",a,"}")') TRIM(float_to_latex(alphaf,10))
   WRITE(iunout,'("\def\amu{",a,"}")') TRIM(float_to_latex(amu,11))
   WRITE(iunout,'("\def\avonum{",a,"}")') TRIM(float_to_latex(avonum,8))
   WRITE(iunout,'("\def\kb{",a,"}")') TRIM(float_to_latex(boltz,6))
   WRITE(iunout,'("\def\rgas{",a,"}")') TRIM(float_to_latex(rgas,14))
   WRITE(iunout,*)
ENDIF
!
! Derived physical constants
!
me=rydberg * 2.0_DP * hplanck / alphaf**2 / cspeed  ! electron mass kg
abohr=alphaf/ 4.0_DP / pi / rydberg    ! m
mu0=2.0_DP*alphaf*hplanck/e**2/cspeed  ! N /A^2       
epsilon0=1.0_DP/mu0/cspeed**2          ! N m^2/C^2        
baru=2.0_DP * hplanck * cspeed * rydberg
bohrmag=hbarf*e/2.0_DP/me

WRITE(stdout,'(/,5x,"Derived Physical quantities:")') 
WRITE(stdout,'(5x,"Electron mass: ",8x,es20.10,"    kg")') me
WRITE(stdout,'(5x,"mu0: ",19x,es20.11,"   N/A^2")') mu0
WRITE(stdout,'(5x,"epsilon0: ",13x,es20.10,"    C^2/N m^2")') epsilon0
WRITE(stdout,'(5x,"E_hatree: ",16x,es20.13," J")') baru
WRITE(stdout,'(5x,"Bohr radius: ",11x,es20.11,"   m")') abohr
WRITE(stdout,'(5x,"Bohr magneton: ",8x,es20.10,"    J/T")') bohrmag

IF (ionode) THEN
   WRITE(iunout,'("\def\me{",a,"}")') TRIM(float_to_latex(me,10))
   WRITE(iunout,'("\def\abohr{",a,"}")') TRIM(float_to_latex(abohr,11))
   WRITE(iunout,'("\def\muzero{",a,"}")') TRIM(float_to_latex(mu0,11))
   WRITE(iunout,'("\def\epsilonzero{",a,"}")') &
                                            TRIM(float_to_latex(epsilon0,10))
   WRITE(iunout,'("\def\ehartree{",a,"}")') TRIM(float_to_latex(baru,13))
   WRITE(iunout,'("\def\bohrmag{",a,"}")') TRIM(float_to_latex(bohrmag,10))
   WRITE(iunout,*)
ENDIF
!
! Unit of measurement
!
barl=abohr
barm=me
bart=hbarf/baru
barrhom=me/abohr**3
barnu=1.0_DP/bart
barv=alphaf * cspeed
bara=baru * alphaf * cspeed / hbarf
barp=hbarf/abohr
baram=hbarf
barf=baru/barl
barw=baru/bart
barpr=barf/barl**2
barc=e
bari=barc * baru / hbarf
barrho=barc/barl**3
barcur=bari/barl**2
bare=barf/barc
barphi=baru/barc
barcap=barc**2/baru
bardip=barc * barl
barpolar=barc/barl**2
bard=barpolar/4.0_DP/pi
barohm=hbarf/barc**2
barb=hbarf/barl**2/barc
barav=hbarf/barl/barc
barwb=hbarf/e
bary=hbarf**2/e**2/baru
barmu=hbarf * e / me
barmag=bari/barl
barh=barmag/4.0_DP/pi

WRITE(stdout,'(/,5x,"Conversion factors (Atomic units - SI):")') 
WRITE(stdout,'(5x,"Length: \l=",13x,es20.11,"   m")') barl
WRITE(stdout,'(5x,"Mass: \m=",14x,es20.10,"    kg")') barm
WRITE(stdout,'(5x,"Mass density: \rhom=",3x,es20.10,"    kg/m^3")') barrhom
WRITE(stdout,'(5x,"Time: \t=",17x,es20.13," s")') bart
WRITE(stdout,'(5x,"Frequency: \nu=",11x,es20.13," Hz")') barnu
WRITE(stdout,'(5x,"Speed: \v=",14x,es20.11,"   m/s")') barv
WRITE(stdout,'(5x,"Acceleration: \a=",6x,es20.10,"    m/s^2")') bara
WRITE(stdout,'(5x,"Momentum: \p=",11x,es20.11,"   kg m/s")') barp
WRITE(stdout,'(5x,"Angular momentum: \L=",5x,es20.13," kg m^2/s")') baram
WRITE(stdout,'(5x,"Force: \f=",13x,es20.10,"    N")') barf
WRITE(stdout,'(5x,"Energy: \U=",15x,es20.13," J")') baru
WRITE(stdout,'(5x,"Power: \W=",16x,es20.13," W")') barw
WRITE(stdout,'(5x,"Pressure: \pr=",9x,es20.10,"    Pa")') barpr
WRITE(stdout,'(5x,"Current: \I=",13x,es20.12,"  A")') bari
WRITE(stdout,'(5x,"Charge: \C=",11x,es20.9,"     C")') barc
WRITE(stdout,'(5x,"Charge density: \rho=",3x,es20.11,"   C/m^3")') barrho
WRITE(stdout,'(5x,"Current density: \J=",4x,es20.11,"   A/m^2")') barcur
WRITE(stdout,'(5x,"Electric field: \E=",5x,es20.11,"   N/C")') bare
WRITE(stdout,'(5x,"Electric potential: \V=",3x,es20.13," V")') barphi
WRITE(stdout,'(5x,"Capacitance: \F=",9x,es20.12,"  F")') barcap
WRITE(stdout,'(5x,"Dipole moment: \dip=",3x,es20.10,"    C m")') bardip
WRITE(stdout,'(5x,"Polarization: \P=",6x,es20.10,"    C/m^2")') barpolar
WRITE(stdout,'(5x,"Electric displacement: \D=",es17.10,"    C/m^2")') bard
WRITE(stdout,'(5x,"Resistance: \R=",11x,es20.13," Ohm")') barohm
WRITE(stdout,'(5x,"Magnetic induction: \B=",1x,es20.11,"   T")') barb
WRITE(stdout,'(5x,"Vector potential: \A=",3x,es20.11,"   T m")') barav
WRITE(stdout,'(5x,"Magnetic field flux: \Phi=",es20.13," Wb")') barwb
WRITE(stdout,'(5x,"Inductance: \Y= ",9x,es20.12,"  H")') bary
WRITE(stdout,'(5x,"Magnetic dipole: \mu=",3x,es20.11,"   A m^2 (J/T)")') barmu
WRITE(stdout,'(5x,"Magnetization: \M=",6x,es20.11,"   A/m")') barmag
WRITE(stdout,'(5x,"Magnetic strength: \H=",1x,es20.10,"    A/m")') barh

IF (ionode) THEN
   WRITE(iunout,'("\def\barl{",a,"}")') TRIM(float_to_latex(barl,11))
   WRITE(iunout,'("\def\barm{",a,"}")') TRIM(float_to_latex(barm,10))
   WRITE(iunout,'("\def\barrhom{",a,"}")') TRIM(float_to_latex(barrhom,10))
   WRITE(iunout,'("\def\bart{",a,"}")') TRIM(float_to_latex(bart,13))
   WRITE(iunout,'("\def\barnu{",a,"}")') TRIM(float_to_latex(barnu,13))
   WRITE(iunout,'("\def\barv{",a,"}")') TRIM(float_to_latex(barv,11))
   WRITE(iunout,'("\def\bara{",a,"}")') TRIM(float_to_latex(bara,10))
   WRITE(iunout,'("\def\barp{",a,"}")') TRIM(float_to_latex(barp,11))
   WRITE(iunout,'("\def\baram{",a,"}")') TRIM(float_to_latex(baram,13))
   WRITE(iunout,'("\def\barf{",a,"}")') TRIM(float_to_latex(barf,10))
   WRITE(iunout,'("\def\baru{",a,"}")') TRIM(float_to_latex(baru,13))
   WRITE(iunout,'("\def\barw{",a,"}")') TRIM(float_to_latex(barw,13))
   WRITE(iunout,'("\def\barpr{",a,"}")') TRIM(float_to_latex(barpr,10))
   WRITE(iunout,'("\def\bari{",a,"}")') TRIM(float_to_latex(bari,12))
   WRITE(iunout,'("\def\barc{",a,"}")') TRIM(float_to_latex(barc,9))
   WRITE(iunout,'("\def\barrho{",a,"}")') TRIM(float_to_latex(barrho,11))
   WRITE(iunout,'("\def\barcur{",a,"}")') TRIM(float_to_latex(barcur,11))
   WRITE(iunout,'("\def\bare{",a,"}")') TRIM(float_to_latex(bare,11))
   WRITE(iunout,'("\def\barphi{",a,"}")') TRIM(float_to_latex(barphi,13))
   WRITE(iunout,'("\def\barcap{",a,"}")') TRIM(float_to_latex(barcap,12))
   WRITE(iunout,'("\def\bardip{",a,"}")') TRIM(float_to_latex(bardip,10))
   WRITE(iunout,'("\def\barpolar{",a,"}")') TRIM(float_to_latex(barpolar,10))
   WRITE(iunout,'("\def\bard{",a,"}")') TRIM(float_to_latex(bard,10))
   WRITE(iunout,'("\def\barohm{",a,"}")') TRIM(float_to_latex(barohm,13))
   WRITE(iunout,'("\def\barb{",a,"}")') TRIM(float_to_latex(barb,11))
   WRITE(iunout,'("\def\barav{",a,"}")') TRIM(float_to_latex(barav,11))
   WRITE(iunout,'("\def\barwb{",a,"}")') TRIM(float_to_latex(barwb,13))
   WRITE(iunout,'("\def\bary{",a,"}")') TRIM(float_to_latex(bary,12))
   WRITE(iunout,'("\def\barmu{",a,"}")') TRIM(float_to_latex(barmu,11))
   WRITE(iunout,'("\def\barmag{",a,"}")') TRIM(float_to_latex(barmag,11))
   WRITE(iunout,'("\def\barh{",a,"}")') TRIM(float_to_latex(barh,10))
   WRITE(iunout,*)
ENDIF


cmtom=1.D-2
gtokg=1.D-3
rhomtorhom=gtokg/cmtom**3
kappa=SQRT(1.D9/4.0_DP/pi/epsilon0)
kappaa=kappa/1.D2/cspeed

ptop=gtokg * cmtom
ltol=gtokg * cmtom**2
ftof=gtokg * cmtom
utou=gtokg * cmtom**2
prtopr=gtokg/cmtom

chtoch=1.0_DP/kappa
kappadiecic=kappa/10.0_DP/cspeed
itoi=chtoch
rhotorho=chtoch / cmtom**3
curtocur=itoi / cmtom**2
etoe= ftof / chtoch
phitophi=etoe * cmtom
captocap=1.D7/kappa**2
diptodip=1.D-2/kappa
polartopolar=chtoch / cmtom**2
dtod=polartopolar/4.0_DP/pi
ohmtoohm=kappa**2*1.D-7
btob=kappaa/1.D3
toverg=1.0_DP/btob
avtoav=1.D-5 * kappaa
wbtowb=kappaa*1D-7
ytoy=1.D-7*kappa**2
mutomu=1.D-4 / kappaa
magtomag=1.D2/kappaa
htoh=magtomag/4.0_DP/pi

WRITE(stdout,'(/,5x,"Conversion factors (c.g.s.-Gaussian - SI):")') 
WRITE(stdout,'(5x,"Length: cm=",3x,es20.1,10x,"   m")') cmtom
WRITE(stdout,'(5x,"Mass: g=",6x,es20.1,10x,"   kg")') gtokg
WRITE(stdout,'(5x,"Mass density: g/cm^3=",es13.1,12x," kg/m^3")') rhomtorhom
WRITE(stdout,'(5x,"Time: s=",19x,"1.0E+00",11x,"  s")') 
WRITE(stdout,'(5x,"Frequency: Hz=",13x,"1.0E+00",11x,"  Hz")') 

WRITE(stdout,'(5x,"Speed: cm/s=",2x,es20.1,10x,"   m/s")') cmtom
WRITE(stdout,'(5x,"Acceleration: cm/s^2=",es13.1,10x,"   m/s^2")') cmtom
WRITE(stdout,'(5x,"Momentum: g cm/s=",es17.1,10x,"   kg m/s")') ptop
WRITE(stdout,'(5x,"Angular momentum: g cm^2/s=",es7.1,10x,"   kg m^2/s")')&
                                                                ltol
WRITE(stdout,'(5x,"Force: dyne=",2x,es20.1,10x,"   N")') ftof
WRITE(stdout,'(5x,"Energy: erg=",2x,es20.1,10x,"   J")') utou
WRITE(stdout,'(5x,"Power: erg/s=",1x,es20.1,10x,"   W")') utou
WRITE(stdout,'(5x,"Pressure: Ba=",1x,es20.1,10x,"   Pa")') prtopr
WRITE(stdout,'(5x,"Current: statA=",9x,es20.11,"   A")') itoi
WRITE(stdout,'(5x,"Charge: statC=",10x,es20.11,"   C")') chtoch
WRITE(stdout,'(5x,"Charge density: statC/cm^3=",es17.11,1x,&
                                                        &"  C/m^3")') rhotorho
WRITE(stdout,'(5x,"Current density: statA/cm^2=",es17.11,1x,&
                                                        &" A/m^2")') curtocur
WRITE(stdout,'(5x,"Electric field: dyne/statC=",es17.11,1x,"  N/C")') &
                                                                     etoe
WRITE(stdout,'(5x,"Electric potential: statV=",es18.11,1x,"  V")') &
                                                                     phitophi
WRITE(stdout,'(5x,"Capacitance: cm=",8x,es20.11,1x,"  F")') captocap
WRITE(stdout,'(5x,"Dipole moment: statC cm=",es20.11,&
                                            &1x,"  C m")') diptodip
WRITE(stdout,'(5x,"Electric polarization: statC/cm^2=",es18.11,&
                                            &" C/m^2")') polartopolar
WRITE(stdout,'(5x,"Electric displacement: statC/cm^2 4 pi=",es18.11,&
                                            &" C/m^2")') dtod
WRITE(stdout,'(5x,"Resistance: s/cm=",8x,es20.12,"  Ohm")') ohmtoohm
WRITE(stdout,'(5x,"Magnetic induction: G=",3x,es20.12,"  T")') btob
WRITE(stdout,'(5x,"Vector potential: G cm=",2x,es20.12,"  T m")') avtoav
WRITE(stdout,'(5x,"Magnetic field flux: Mx=",1x,es20.12,"  Wb")') wbtowb
WRITE(stdout,'(5x,"Inductance: statH=",6x,es20.11,2x," H")') ytoy
WRITE(stdout,'(5x,"Magnetic dipole: statC cm=",es18.11,2x,&
                                            &" A m^2")') mutomu
WRITE(stdout,'(5x,"Magnetization: statC/cm^2=",es18.11,2x," A/m")') &
                                                         magtomag
WRITE(stdout,'(5x,"Magnetic strength: statC/cm^2 4 pi=",&
                                     &es18.11," A/m")') htoh

WRITE(stdout,'(/,5x,"C/statC=",es18.11," (C/statC)/10c=",es18.11)') kappa,&
                                                   kappadiecic
WRITE(stdout,'(/,5x,"mu_0/4 pi 10^-7=",3x,es20.11)') mu0/4.0_DP/pi/1.D-7

WRITE(stdout,'(/,5x,"Conversion factors (SI - c.g.s.-Gaussian):")') 
WRITE(stdout,'(5x,"Length: m=",4x,es20.1,10x,"   cm")') 1.0_DP/cmtom
WRITE(stdout,'(5x,"Mass: kg=",5x,es20.1,10x,"   g")')   1.0_DP/gtokg
WRITE(stdout,'(5x,"Time: s=",19x,"1.0E+00",11x,"  s")') 
WRITE(stdout,'(5x,"Frequency: Hz=",13x,"1.0E+00",11x,"  Hz")') 
WRITE(stdout,'(5x,"Speed: m/s=",3x,es20.1,10x,"   cm/s")') 1.0_DP/cmtom
WRITE(stdout,'(5x,"Acceleration: m/s^2=",es14.1,10x,"   cm/s^2")') &
                                                              1.0_DP/cmtom
WRITE(stdout,'(5x,"Momentum: kg m/s=",es17.1,10x"   g cm/s")') 1.0_DP/ptop
WRITE(stdout,'(5x,"Angular momentum: kg m^2/s=",es7.1,10x,"   g cm^2/s")')&
                                                                1.0_DP/ltol
WRITE(stdout,'(5x,"Force: N=",5x,es20.1,10x,"   dyne")') 1.0_DP/ftof
WRITE(stdout,'(5x,"Energy: J=",4x,es20.1,10x,"   erg")') 1.0_DP/utou
WRITE(stdout,'(5x,"Power: W=",5x,es20.1,10x,"   erg/s")') 1.0_DP/utou
WRITE(stdout,'(5x,"Pressure: Pa=",1x,es20.1,10x,"   Ba")') 1.0_DP/prtopr

WRITE(stdout,'(5x,"Current: A=",13x,es20.11,"   statA")') 1.0_DP/itoi
WRITE(stdout,'(5x,"Charge: C=",14x,es20.11,"   statC")') 1.0_DP/chtoch
WRITE(stdout,'(5x,"Charge density: C/m^3=",2x,es20.11,1x,&
                                           &"  statC/cm^3")') 1.0_DP/rhotorho
WRITE(stdout,'(5x,"Current density: A/m^2=",1x,es20.11,1x,&
                                           &"  statA/cm^2")') 1.0_DP/curtocur
WRITE(stdout,'(5x,"Electric field: N/C=",4x,es20.11,1x,"  dyne/statC")') &
                                                             1.0_DP/etoe
WRITE(stdout,'(5x,"Electric potential: V=",2x,es20.11,1x,"  statV")') &
                                                          1.0_DP/phitophi
WRITE(stdout,'(5x,"Capacitance: F=",9x,es20.11,1x,"  cm")') 1.0_DP/captocap
WRITE(stdout,'(5x,"Dipole moment: C m=",5x,es20.11,1x,"  statC cm")') 1.0_DP/diptodip
WRITE(stdout,'(5x,"Electric polarization: C/m^2=",es18.11,&
                                      &" statC/cm^2")') 1.0_DP/polartopolar
WRITE(stdout,'(5x,"Electric displacement: C/m^2=",es18.11,&
                                            &" statC/cm^2 4pi")') 1.0_DP/dtod
WRITE(stdout,'(5x,"Resistance: Ohm=",9x,es20.12,"  s/cm")') &
                                                          1.0_DP/ohmtoohm
WRITE(stdout,'(5x,"Magnetic induction: T=",3x,es20.12,&
                                            &"  G")') 1.0_DP/btob
WRITE(stdout,'(5x,"Vector potential: T m=",3x,es20.12,&
                                            &"  G cm")') 1.0_DP/avtoav
WRITE(stdout,'(5x,"Magnetic field flux: Wb=",1x,es20.12,&
                                            &"  Mx")') 1.0_DP/wbtowb
WRITE(stdout,'(5x,"Inductance: H=",10x,es20.11,1x,"  statH")') 1.0_DP/ytoy
WRITE(stdout,'(5x,"Magnetic dipole: A m^2=",1x,es20.11,1x,&
                                            &"  statC cm")') 1.0_DP/mutomu
WRITE(stdout,'(5x,"Magnetization: A/m=",5x,es20.11,1x,"  statC/cm^2")') &
                                                         1.0_DP/magtomag
WRITE(stdout,'(5x,"Magnetic strength: A/m=",&
                                &1x,es20.11,"   statC/cm^2 4pi")') 1.0_DP/htoh

IF (ionode) THEN
   WRITE(iunout,'("\def\cmtom{",a,"}")') TRIM(float_to_latex(cmtom,1))
   WRITE(iunout,'("\def\gtokg{",a,"}")') TRIM(float_to_latex(gtokg,1))

   WRITE(iunout,'("\def\ptop{",a,"}")') TRIM(float_to_latex(ptop,1))
   WRITE(iunout,'("\def\ltol{",a,"}")') TRIM(float_to_latex(ltol,1))
   WRITE(iunout,'("\def\ftof{",a,"}")') TRIM(float_to_latex(ftof,1))
   WRITE(iunout,'("\def\utou{",a,"}")') TRIM(float_to_latex(utou,1))
   WRITE(iunout,'("\def\prtopr{",a,"}")') TRIM(float_to_latex(prtopr,1))
   WRITE(iunout,'("\def\chtoch{",a,"}")') TRIM(float_to_latex(chtoch,11))
   WRITE(iunout,'("\def\kappaa{",a,"}")') TRIM(float_to_latex(kappaa,12))
   WRITE(iunout,'("\def\kappa{",a,"}")') TRIM(float_to_latex(kappa,11))
   WRITE(iunout,'("\def\kappadiecic{",a,"}")') &
                                         TRIM(float_to_latex(kappadiecic,11))
   WRITE(iunout,'("\def\itoi{",a,"}")') TRIM(float_to_latex(itoi,11))
   WRITE(iunout,'("\def\rhotorho{",a,"}")') TRIM(float_to_latex(rhotorho,11))
   WRITE(iunout,'("\def\curtocur{",a,"}")') TRIM(float_to_latex(curtocur,11))
   WRITE(iunout,'("\def\etoe{",a,"}")') TRIM(float_to_latex(etoe,11))
   WRITE(iunout,'("\def\phitophi{",a,"}")') TRIM(float_to_latex(phitophi,11))
   WRITE(iunout,'("\def\captocap{",a,"}")') TRIM(float_to_latex(captocap,11))
   WRITE(iunout,'("\def\diptodip{",a,"}")') TRIM(float_to_latex(diptodip,11))
   WRITE(iunout,'("\def\polartopolar{",a,"}")') &
                                      TRIM(float_to_latex(polartopolar,11))
   WRITE(iunout,'("\def\dtod{",a,"}")') TRIM(float_to_latex(dtod,11))
   WRITE(iunout,'("\def\ohmtoohm{",a,"}")') TRIM(float_to_latex(ohmtoohm,10))
   WRITE(iunout,'("\def\btob{",a,"}")') TRIM(float_to_latex(btob,12))
   WRITE(iunout,'("\def\avtoav{",a,"}")') TRIM(float_to_latex(avtoav,12))
   WRITE(iunout,'("\def\wbtowb{",a,"}")') TRIM(float_to_latex(wbtowb,12))
   WRITE(iunout,'("\def\ytoy{",a,"}")') TRIM(float_to_latex(ytoy,10))
   WRITE(iunout,'("\def\mutomu{",a,"}")') TRIM(float_to_latex(mutomu,11))
   WRITE(iunout,'("\def\magtomag{",a,"}")') TRIM(float_to_latex(magtomag,11))
   WRITE(iunout,'("\def\htoh{",a,"}")') TRIM(float_to_latex(htoh,11))
   WRITE(iunout,'("\def\toverg{",a,"}")') TRIM(float_to_latex(toverg,11))
   WRITE(iunout,*)
ENDIF

barlcgs=barl/cmtom
barmcgs=me/gtokg
barrhomcgs=barrhom/rhomtorhom
barvcgs=barv / cmtom
baracgs=bara / cmtom
barpcgs=barp / ptop
baramcgs=baram/ltol
barfcgs=barf/ftof
barucgs=baru/utou
barwcgs=barw/utou
barprcgs=barpr/prtopr
barccgs=barc/chtoch
baricgs=bari/itoi
barrhocgs=barrho/rhotorho
barcurcgs=barcur/curtocur
barecgs=bare/etoe
barphicgs=barphi/phitophi
bardipcgs=bardip/diptodip
barpolarcgs=barpolar/polartopolar
bardcgs=bard/dtod
barbcgs=barb/btob
baravcgs=barav/avtoav
barwbcgs=barwb/wbtowb
barmucgs=barmu/mutomu
barmagcgs=barmag/magtomag
barhcgs=barmag/htoh/4.0_DP/pi
barohmcgs=barohm/ohmtoohm
barcapcgs=barcap/captocap
barycgs=bary/ytoy

WRITE(stdout,'(/,5x,"Conversion factors (Atomic units - c.g.s.-Gaussian):")') 
WRITE(stdout,'(5x,"Length:",17x,es20.11,"   cm")') barlcgs
WRITE(stdout,'(5x,"Mass:",18x,es20.10,"    g")') barmcgs
WRITE(stdout,'(5x,"Mass density:",10x,es20.10,"    g/cm^3")') barrhomcgs
WRITE(stdout,'(5x,"Time:",21x,es20.13," s")') bart
WRITE(stdout,'(5x,"Frequency:",16x,es20.13," Hz")') barnu
WRITE(stdout,'(5x,"Speed:",18x,es20.11,"   cm/s")') barvcgs
WRITE(stdout,'(5x,"Acceleration:",10x,es20.10,"    cm/s^2")') baracgs
WRITE(stdout,'(5x,"Momentum:",15x,es20.11,"   g cm/s")') barpcgs
WRITE(stdout,'(5x,"Angular momentum:",9x,es20.13," g cm^2/s")') baramcgs
WRITE(stdout,'(5x,"Force:",17x,es20.10,"    dyne")') barfcgs
WRITE(stdout,'(5x,"Energy:",19x,es20.13," erg")') barucgs
WRITE(stdout,'(5x,"Power:",20x,es20.13," erg/s")') barwcgs
WRITE(stdout,'(5x,"Pressure:",14x,es20.10,"    Ba")') barprcgs
WRITE(stdout,'(5x,"Current:",16x,es20.11,"   statA")') baricgs
WRITE(stdout,'(5x,"Charge:",17x,es20.11,"   statC")') barccgs
WRITE(stdout,'(5x,"Charge density:",8x,es20.10,3x," statC/cm^3")') barrhocgs
WRITE(stdout,'(5x,"Current density:",7x,es20.10,3x," statA/cm^2")') barcurcgs
WRITE(stdout,'(5x,"Electric field:",9x,es20.11,2x," dyne/statC")') barecgs
WRITE(stdout,'(5x,"Electric potential:",5x,es20.11,2x," statV")') barphicgs
WRITE(stdout,'(5x,"Capacitance:",12x,es20.11,2x," cm")') barcapcgs
WRITE(stdout,'(5x,"Dipole moment:",10x,es20.11,2x," statC cm")') bardipcgs
WRITE(stdout,'(5x,"Polarization:",11x,es20.11,2x," statC/cm^2")') barpolarcgs
WRITE(stdout,'(5x,"Electric displacement:",2x,es20.11,&
                                           &"   statC/cm^2 4pi")') bardcgs
WRITE(stdout,'(5x,"Resistance:",13x,es20.11,2x," s/cm")') barohmcgs
WRITE(stdout,'(5x,"Magnetic induction:",5x,es20.11,2x," G")') barbcgs
WRITE(stdout,'(5x,"Vector potential:",7x,es20.11,2x," G cm")') baravcgs
WRITE(stdout,'(5x,"Magnetic field flux:",4x,es20.11,2x," Mx")') barwbcgs
WRITE(stdout,'(5x,"Inductance:",13x,es20.11,"   statH")') barycgs
WRITE(stdout,'(5x,"Magnetic dipole:",8x,es20.11,"   statC cm")') barmucgs
WRITE(stdout,'(5x,"Magnetization:",10x,es20.11,&
                                           &"   statC/cm^2")') barmagcgs
WRITE(stdout,'(5x,"Magnetic strength:",6x,es20.11,"   statC/cm^2 4pi ")') &
                                                               barhcgs

IF (ionode) THEN
   WRITE(iunout,'("\def\barlcgs{",a,"}")') TRIM(float_to_latex(barlcgs,11))
   WRITE(iunout,'("\def\barmcgs{",a,"}")') TRIM(float_to_latex(barmcgs,10))
   WRITE(iunout,'("\def\barrhomcgs{",a,"}")') &
                                           TRIM(float_to_latex(barrhomcgs,10))
   WRITE(iunout,'("\def\barvcgs{",a,"}")') TRIM(float_to_latex(barvcgs,11))
   WRITE(iunout,'("\def\baracgs{",a,"}")') TRIM(float_to_latex(baracgs,10))
   WRITE(iunout,'("\def\barpcgs{",a,"}")') TRIM(float_to_latex(barpcgs,11))
   WRITE(iunout,'("\def\baramcgs{",a,"}")') TRIM(float_to_latex(baramcgs,13))
   WRITE(iunout,'("\def\barfcgs{",a,"}")') TRIM(float_to_latex(barfcgs,10))
   WRITE(iunout,'("\def\barucgs{",a,"}")') TRIM(float_to_latex(barucgs,13))
   WRITE(iunout,'("\def\barwcgs{",a,"}")') TRIM(float_to_latex(barwcgs,13))
   WRITE(iunout,'("\def\barprcgs{",a,"}")') TRIM(float_to_latex(barprcgs,10))
   WRITE(iunout,'("\def\baricgs{",a,"}")') TRIM(float_to_latex(baricgs,11))
   WRITE(iunout,'("\def\barccgs{",a,"}")') TRIM(float_to_latex(barccgs,11))
   WRITE(iunout,'("\def\barrhocgs{",a,"}")') TRIM(float_to_latex(barrhocgs,10))
   WRITE(iunout,'("\def\barcurcgs{",a,"}")') TRIM(float_to_latex(barcurcgs,10))
   WRITE(iunout,'("\def\barecgs{",a,"}")') TRIM(float_to_latex(barecgs,11))
   WRITE(iunout,'("\def\barphicgs{",a,"}")') TRIM(float_to_latex(barphicgs,11))
   WRITE(iunout,'("\def\barcapcgs{",a,"}")') TRIM(float_to_latex(barcapcgs,11))
   WRITE(iunout,'("\def\bardipcgs{",a,"}")') TRIM(float_to_latex(bardipcgs,11))
   WRITE(iunout,'("\def\barpolarcgs{",a,"}")') TRIM(float_to_latex(barpolarcgs,11))
   WRITE(iunout,'("\def\bardcgs{",a,"}")') TRIM(float_to_latex(bardcgs,11))
   WRITE(iunout,'("\def\barohmcgs{",a,"}")') TRIM(float_to_latex(barohmcgs,11))
   WRITE(iunout,'("\def\barbcgs{",a,"}")') TRIM(float_to_latex(barbcgs,11))
   WRITE(iunout,'("\def\baravcgs{",a,"}")') TRIM(float_to_latex(baravcgs,11))
   WRITE(iunout,'("\def\barwbcgs{",a,"}")') TRIM(float_to_latex(barwbcgs,11))
   WRITE(iunout,'("\def\barycgs{",a,"}")') TRIM(float_to_latex(barycgs,11))
   WRITE(iunout,'("\def\barmucgs{",a,"}")') TRIM(float_to_latex(barmucgs,11))
   WRITE(iunout,'("\def\barmagcgs{",a,"}")') TRIM(float_to_latex(barmagcgs,11))
   WRITE(iunout,'("\def\barhcgs{",a,"}")') TRIM(float_to_latex(barhcgs,11))
   WRITE(iunout,*)
ENDIF

barmry=2.0_DP * me
bartry=2.0_DP * bart
barcry=barc/sqrt2
barnury=0.5_DP * barnu
barvry=0.5_DP * barv
barary=0.5_DP * bara
barfry=0.5_DP * barf
barury=0.5_DP * baru
barwry=0.25_DP * barw
barprry=0.5_DP * barpr
bariry=0.5_DP * bari / sqrt2
barrhory=barrho / sqrt2
barcurry=0.5_DP * barcur / sqrt2
barery=bare / sqrt2
barphiry=barphi / sqrt2
bardipry=bardip / sqrt2
barpolarry=barpolar / sqrt2
bardry=bard / sqrt2
barohmry=2.0_DP * barohm
barbry=sqrt2 * barb
baravry=sqrt2 * barav
barwbry=sqrt2 * barwb
baryry=4.0_DP * bary
barmury=0.5_DP * barmu / sqrt2
barmagry=0.5_DP * barmag / sqrt2
barhry=0.5_DP * barh / sqrt2

WRITE(stdout,'(/,5x,"Conversion factors (Rydberg atomic units - SI):")') 
WRITE(stdout,'(5x,"Length: \l_R=",11x,es20.11,"   m")') barl
WRITE(stdout,'(5x,"Mass: \m_R=",13x,es20.11,"   kg")') barmry
WRITE(stdout,'(5x,"Time: \t_R=",15x,es20.13," s")') bartry
WRITE(stdout,'(5x,"Frequency: \nu_R=",9x,es20.13," Hz")') barnury
WRITE(stdout,'(5x,"Speed: \v_R=",12x,es20.11,"   m/s")') barvry
WRITE(stdout,'(5x,"Acceleration: \a_R=",5x,es20.11,"   m/s^2")') barary
WRITE(stdout,'(5x,"Momentum: \p_R=",9x,es20.11,"   kg m/s")') barp
WRITE(stdout,'(5x,"Angular momentum: \L_R=",3x,es20.13," kg m^2/s")') baram
WRITE(stdout,'(5x,"Force: \f_R=",12x,es20.11,"   N")') barfry
WRITE(stdout,'(5x,"Energy: \U_R=",13x,es20.13," J")') barury
WRITE(stdout,'(5x,"Power: \W_R=",13x,es20.12,"  W")') barwry
WRITE(stdout,'(5x,"Pressure: \pr_R=",8x,es20.11,"   Pa")') barprry
WRITE(stdout,'(5x,"Current: \I_R=",12x,es20.13," A")') bariry
WRITE(stdout,'(5x,"Charge: \C_R=",13x,es20.13," C")') barcry
WRITE(stdout,'(5x,"Charge density: \rho_R=",es20.10,"    C/m^3")') barrhory
WRITE(stdout,'(5x,"Current density: \J_R=",2x,es20.11,"   A/m^2")') barcurry
WRITE(stdout,'(5x,"Electric field: \E_R=",3x,es20.11,"   N/C")') barery
WRITE(stdout,'(5x,"Electric potential: \V_R=",1x,es20.13," V")') barphiry
WRITE(stdout,'(5x,"Capacitance: \F_R=",7x,es20.12,"  F")') barcap
WRITE(stdout,'(5x,"Dipole moment: \dip_R=",2x,es20.11,"   C m")') bardipry
WRITE(stdout,'(5x,"Polarization: \P_R=",4x,es20.10,"    C/m^2")') barpolarry
WRITE(stdout,'(5x,"Electric displacement: \D_R=",es17.11,"   C/m^2")') &
                                                                  bardry
WRITE(stdout,'(5x,"Resistance: \R_R=",9x,es20.13," Ohm")') barohmry
WRITE(stdout,'(5x,"Magnetic induction: \B_R=",es18.10,"    T")') barbry
WRITE(stdout,'(5x,"Vector potential: \A_R=",1x,es20.11,"   T m")') baravry
WRITE(stdout,'(5x,"Magnetic field flux: \Phi_R=",es19.13," Wb")') barwbry
WRITE(stdout,'(5x,"Inductance: \Y_R= ",8x,es20.13," H")') baryry
WRITE(stdout,'(5x,"Magnetic dipole: \mu_R=",es20.10,"    A m^2 (J/T)")') &
                                                                    barmury
WRITE(stdout,'(5x,"Magnetization: \M_R=",4x,es20.11,"   A/m")') barmagry
WRITE(stdout,'(5x,"Magnetic field: \H_R=",3x,es20.11,"   A/m")') barhry

IF (ionode) THEN
   WRITE(iunout,'("\def\barmry{",a,"}")') TRIM(float_to_latex(barmry,11))
   WRITE(iunout,'("\def\bartry{",a,"}")') TRIM(float_to_latex(bartry,13))
   WRITE(iunout,'("\def\barnury{",a,"}")') TRIM(float_to_latex(barnury,13))
   WRITE(iunout,'("\def\barvry{",a,"}")') TRIM(float_to_latex(barvry,11))
   WRITE(iunout,'("\def\barary{",a,"}")') TRIM(float_to_latex(barary,11))
   WRITE(iunout,'("\def\barfry{",a,"}")') TRIM(float_to_latex(barfry,11))
   WRITE(iunout,'("\def\barury{",a,"}")') TRIM(float_to_latex(barury,13))
   WRITE(iunout,'("\def\barwry{",a,"}")') TRIM(float_to_latex(barwry,12))
   WRITE(iunout,'("\def\barprry{",a,"}")') TRIM(float_to_latex(barprry,11))
   WRITE(iunout,'("\def\bariry{",a,"}")') TRIM(float_to_latex(bariry,13))
   WRITE(iunout,'("\def\barcry{",a,"}")') TRIM(float_to_latex(barcry,13))
   WRITE(iunout,'("\def\barrhory{",a,"}")') TRIM(float_to_latex(barrhory,10))
   WRITE(iunout,'("\def\barcurry{",a,"}")') TRIM(float_to_latex(barcurry,11))
   WRITE(iunout,'("\def\barery{",a,"}")') TRIM(float_to_latex(barery,11))
   WRITE(iunout,'("\def\barphiry{",a,"}")') TRIM(float_to_latex(barphiry,13))
   WRITE(iunout,'("\def\bardipry{",a,"}")') TRIM(float_to_latex(bardipry,11))
   WRITE(iunout,'("\def\barpolarry{",a,"}")') TRIM(float_to_latex(barpolarry,10))
   WRITE(iunout,'("\def\bardry{",a,"}")') TRIM(float_to_latex(bardry,11))
   WRITE(iunout,'("\def\barohmry{",a,"}")') TRIM(float_to_latex(barohmry,13))
   WRITE(iunout,'("\def\barbry{",a,"}")') TRIM(float_to_latex(barbry,10))
   WRITE(iunout,'("\def\baravry{",a,"}")') TRIM(float_to_latex(baravry,11))
   WRITE(iunout,'("\def\barwbry{",a,"}")') TRIM(float_to_latex(barwbry,13))
   WRITE(iunout,'("\def\baryry{",a,"}")') TRIM(float_to_latex(baryry,13))
   WRITE(iunout,'("\def\barmury{",a,"}")') TRIM(float_to_latex(barmury,10))
   WRITE(iunout,'("\def\barmagry{",a,"}")') TRIM(float_to_latex(barmagry,11))
   WRITE(iunout,'("\def\barhry{",a,"}")') TRIM(float_to_latex(barhry,11))
   WRITE(iunout,*)
ENDIF

barbg = barb * alphaf
baravg = barav * alphaf
barwbg = barwb * alphaf
barmug = barmu / alphaf
barmagg = barmag / alphaf
barhg = barh / alphaf

WRITE(stdout,'(/,5x,"Conversion factors (Gaussian atomic units - SI):")') 
WRITE(stdout,'(5x,"Magnetic induction: \B_G=",es19.11,"   T")') barbg
WRITE(stdout,'(5x,"Vector potential: \A_G=",es20.10,"    T m")') baravg
WRITE(stdout,'(5x,"Magnetic field flux: \Phi_G=",es17.11,"  Wb")') barwbg
WRITE(stdout,'(5x,"Magnetic dipole: \mu_G=",es20.10,"    A m^2 (J/T)")') &
                                                                    barmug
WRITE(stdout,'(5x,"Magnetization: \M_G=",4x,es20.11,"   A/m")') barmagg
WRITE(stdout,'(5x,"Magnetic field: \H_G=",3x,es20.11,"   A/m")') barhg


IF (ionode) THEN
   WRITE(iunout,'("\def\barbg{",a,"}")') TRIM(float_to_latex(barbg,11))
   WRITE(iunout,'("\def\baravg{",a,"}")') TRIM(float_to_latex(baravg,10))
   WRITE(iunout,'("\def\barwbg{",a,"}")') TRIM(float_to_latex(barwbg,11))
   WRITE(iunout,'("\def\barmug{",a,"}")') TRIM(float_to_latex(barmug,10))
   WRITE(iunout,'("\def\barmagg{",a,"}")') TRIM(float_to_latex(barmagg,11))
   WRITE(iunout,'("\def\barhg{",a,"}")') TRIM(float_to_latex(barhg,11))
   WRITE(iunout,*)
ENDIF

cspeedau=cspeed/barv
amuau=amu/barm
ryev=barphi*0.5_DP
hzcmm1=1.D-2 / cspeed
cmm1hz=1.D2 * cspeed

WRITE(stdout,'(/,5x,"Physical constants in Hartree atomic units:")') 
WRITE(stdout,'(5x,"Speed of light:",9x,es20.11)') cspeedau
WRITE(stdout,'(5x,"atomic mass unit:",6x,es20.10)') amuau

WRITE(stdout,'(/,5x,"Physical constants in eV:")') 
WRITE(stdout,'(5x,"Hartree in eV:",9x,es20.10)') barphi
WRITE(stdout,'(5x,"Rydberg in eV:",9x,es20.10)') ryev

WRITE(stdout,'(/,5x,"Frequency conversion:")') 
WRITE(stdout,'(5x,"Hz in cm^-1:",15x,es20.14)') hzcmm1
WRITE(stdout,'(5x,"cm^-1 in Hz:",9x,es20.8)') cmm1hz

IF (ionode) THEN
   WRITE(iunout,'("\def\cspeedau{",a,"}")') TRIM(float_to_latex(cspeedau,11))
   WRITE(iunout,'("\def\amuau{",a,"}")') TRIM(float_to_latex(amuau,10))
   WRITE(iunout,'("\def\ryev{",a,"}")') TRIM(float_to_latex(ryev,13))
   WRITE(iunout,'("\def\hzcmm1{",a,"}")') TRIM(float_to_latex(hzcmm1,14))
   WRITE(iunout,'("\def\cmm1hz{",a,"}")') TRIM(float_to_latex(cmm1hz,8))
   WRITE(iunout,*)
ENDIF

barifc=baru/barl**2
bardmc=e*baru*abohr/hbarf
baralpha=1.0_DP/barv
baralphap=4.0_DP * pi * alphaf**2 * baralpha
WRITE(stdout,'(/,5x,"Material properties (a.u. - SI):")') 
WRITE(stdout,'(5x,"Int. force const.: \ifc=",es20.11,"   J/m^2")') barifc
WRITE(stdout,'(5x,"Dyn. mag. charge: \dmc=",1x,es20.11,"   A m")') bardmc
WRITE(stdout,'(5x,"ME tensor (I): \alpha=",2x,es20.11,"   s/m")') baralpha
WRITE(stdout,'(5x,"ME tensor (II): \alpha''=",es20.11,"   s/m")') baralphap

IF (ionode) THEN
   WRITE(iunout,'("\def\barifc{",a,"}")') TRIM(float_to_latex(barifc,11))
   WRITE(iunout,'("\def\bardmc{",a,"}")') TRIM(float_to_latex(bardmc,11))
   WRITE(iunout,'("\def\baralpha{",a,"}")') TRIM(float_to_latex(baralpha,11))
   WRITE(iunout,'("\def\baralphap{",a,"}")') TRIM(float_to_latex(baralphap,11))
ENDIF

!ifctoifc=baru/barl**2
!dmctodmc=e*baru*abohr/hbarf
zmtozm=1.D-2/kappaa
alphatoalpha=mu0 * magtomag / etoe / 4.0_DP / pi
alphaptoalphap=mu0 * magtomag / etoe 
WRITE(stdout,'(/,5x,"Material properties (c.g.s. - SI):")') 
!WRITE(stdout,'(5x,"Int. force const.: \ifc=",13x,es20.11,"   J/m^2")') barifc
WRITE(stdout,'(5x,"Dyn. mag. charge: \dmc=",1x,es20.11,"   A m")') zmtozm 
WRITE(stdout,'(5x,"ME tensor (I): \alpha=",2x,es20.11,"   s/m")') alphatoalpha
WRITE(stdout,'(5x,"ME tensor (II): \alpha''=",es20.11,"   s/m")') &
                                                                 alphaptoalphap

IF (ionode) THEN
!   WRITE(iunout,'("\def\barifc{",a,"}")') TRIM(float_to_latex(barifc,11))
!   WRITE(iunout,'("\def\bardmc{",a,"}")') TRIM(float_to_latex(bardmc,11))
   WRITE(iunout,'("\def\zmtozm{",a,"}")') TRIM(float_to_latex(zmtozm,11))
   WRITE(iunout,'("\def\alphatoalpha{",a,"}")') &
                                 TRIM(float_to_latex(alphatoalpha,11))
   WRITE(iunout,'("\def\alphaptoalphap{",a,"}")') &
                                 TRIM(float_to_latex(alphaptoalphap,11))
ENDIF

WRITE(stdout,'(/,5x,70("-"))') 
WRITE(stdout,'(5x,"Errors:                     Absolute                  Relative",/)') 

rydbergerr=0.0000000000021D7
alphaerr=0.0000000011D-3
amuerr=0.00000000050D-27

WRITE(stdout,'(5x,"rydberg = ",5x,es20.2," 1/m",4x,es18.2)') rydbergerr, rydbergerr/rydberg
WRITE(stdout,'(5x,"alpha = ",7x,es20.2,2x,4x,es20.2)') alphaerr, alphaerr/alphaf
WRITE(stdout,'(5x,"amu = ",9x,es20.2,2x,4x,es20.2)') amuerr, amuerr/amu
WRITE(stdout,*)

meerr=( rydbergerr / rydberg + 2.0_DP * alphaerr/alphaf) * me
abohrerr=(rydbergerr / rydberg + alphaerr/alphaf) * abohr 
mu0err=alphaerr * mu0 / alphaf
epsilon0err= mu0err * epsilon0 / mu0
hartreeerr=(rydbergerr / rydberg ) * baru
bohrmagerr=meerr/me * bohrmag

WRITE(stdout,'(5x,"me = ",10x,es20.2," kg",4x,es19.2)') meerr, meerr/me
WRITE(stdout,'(5x,"abohr = ",7x,es20.2," m",4x,es20.2)') &
                                                    abohrerr, abohrerr/abohr
WRITE(stdout,'(5x,"mu0 = ",9x,es20.2," N/A^2",4x,es16.2)') mu0err, mu0err/mu0
WRITE(stdout,'(5x,"epsilon0 = ",4x,es20.2," C^2/Nm^2",4x,es13.2)') &
                                       epsilon0err, epsilon0err/epsilon0
WRITE(stdout,'(5x,"hartree = ",5x,es20.2," J", 4x,es20.2)') hartreeerr, &
                                                    hartreeerr/baru
WRITE(stdout,'(5x,"bohr mag = ",4x,es20.2," J/T", 4x,es18.2)') bohrmagerr, &
                                                    bohrmagerr/bohrmag
WRITE(stdout,*)

barrhomerr=(meerr/me + 3.0_DP * abohrerr/abohr)*barrhom
barterr=hartreeerr/baru * bart
barnuerr=barterr/bart * barnu
barverr=(alphaerr/alphaf)*barv
baraerr=(hartreeerr/baru + alphaerr/alphaf) * bara
barperr=abohrerr/barl * barp
barferr=(abohrerr/barl + hartreeerr/baru )* barf
barwerr=(hartreeerr/baru + barterr/bart) * barw
barprerr=(2.0_DP * abohrerr/abohr + barferr/barf ) * barpr
barierr=(hartreeerr/baru) * bari
barrhoerr=3.0_DP*abohrerr/abohr*barrho
barcurerr=(barierr/bari + 2.0_DP*abohrerr/abohr)*barcur
bareerr= barferr/ barf * bare
barphierr=hartreeerr/baru * barphi
barcaperr=hartreeerr/baru*barcap
bardiperr=abohrerr / abohr *bardip
barpolarerr=2.0_DP * abohrerr/abohr * barpolar
barderr=barpolarerr/barpolar * bard
barberr=2.0_DP * abohrerr/abohr * barb
baraverr=abohrerr/abohr * barav
baryerr=hartreeerr/baru * bary
barmuerr=meerr/me * barmu
barmagerr=(barierr/bari + abohrerr/abohr)*barmag 
barherr=barmagerr/barmag * barh

WRITE(stdout,'(5x,"Errors of conversion factors (atomic units - SI):")') 
WRITE(stdout,'(5x,"\l = ",10x,es20.2," m", 4x,es20.2)') abohrerr, abohrerr/barl
WRITE(stdout,'(5x,"\m = ",10x,es20.2," kg", 3x,es20.2)') meerr, meerr/me
WRITE(stdout,'(5x,"\rhom = ",7x,es20.2," kg/m^3", es19.2)') barrhomerr, &
                                                           barrhomerr/barrhom
WRITE(stdout,'(5x,"\t = ",10x,es20.2," s", 4x,es20.2)') barterr, barterr/bart
WRITE(stdout,'(5x,"\nu = ",9x,es20.2," s", 4x,es20.2)') barnuerr, &
                                                        barnuerr/barnu
WRITE(stdout,'(5x,"\v = ",10x,es20.2," m/s", 4x,es18.2)') barverr, barverr/barv
WRITE(stdout,'(5x,"\a = ",10x,es20.2," m/s^2", 4x,es16.2)') baraerr, &
                                                        baraerr/bara
WRITE(stdout,'(5x,"\p = ",10x,es20.2," kg m/s", 4x,es15.2)') barperr, &
                                                        barperr/barp
WRITE(stdout,'(5x,"\L = ",10x,es20.2," kg m^2/s",es17.2)') zero, zero
WRITE(stdout,'(5x,"\f = ",10x,es20.2," N", 4x,es20.2)') barferr, barferr/barf
WRITE(stdout,'(5x,"\U = ",10x,es20.2," J", 4x,es20.2)') hartreeerr, &
                                                        hartreeerr/baru
WRITE(stdout,'(5x,"\W = ",10x,es20.2," W", 4x,es20.2)') barwerr, barwerr/barw
WRITE(stdout,'(5x,"\pr = ",9x,es20.2," Pa", 4x,es19.2)') barprerr, &
                                                        barprerr/barpr
WRITE(stdout,'(5x,"\I = ",10x,es20.2," A", 4x,es20.2)') barierr, barierr/bari
WRITE(stdout,'(5x,"\C = ",10x,es20.2," C", 4x,es20.2)') zero, zero
WRITE(stdout,'(5x,"\rho = ",8x,es20.2," C/m^3", 4x,es16.2)') barrhoerr, &
                                                        barrhoerr/barrho
WRITE(stdout,'(5x,"\J = ",11x,es19.2," A/m^2", 4x,es16.2)') barcurerr, &
                                                        barcurerr/barcur
WRITE(stdout,'(5x,"\E = ",10x,es20.2," V/m", 4x,es18.2)') bareerr, &
                                                        bareerr/bare
WRITE(stdout,'(5x,"\V = ",10x,es20.2," V", 4x,es20.2)') barphierr, &
                                                        barphierr/barphi
WRITE(stdout,'(5x,"\F = ",10x,es20.2," F", 4x,es20.2)') barcaperr, &
                                                        barcaperr/barcap
WRITE(stdout,'(5x,"\dip = ",6x,es22.2," C m", 4x,es18.2)') bardiperr, &
                                                        bardiperr/bardip
WRITE(stdout,'(5x,"\P = ",10x,es20.2," C/m^2", 4x,es16.2)') barpolarerr, &
                                                        barpolarerr/barpolar
WRITE(stdout,'(5x,"\D = ",10x,es20.2," C/m^2",4x, es16.2)') barderr, &
                                                        barderr/bard
WRITE(stdout,'(5x,"\R = ",10x,es20.2," Ohm", 2x,es20.2)') zero, zero
WRITE(stdout,'(5x,"\B = ",10x,es20.2," T", 4x,es20.2)') barberr, barberr/barb
WRITE(stdout,'(5x,"\A = ",10x,es20.2," T m", 2x,es20.2)') baraverr, &
                                                          baraverr/barav
WRITE(stdout,'(5x,"\phi = ",8x,es20.2," Wb", 3x,es20.2)') zero, zero
WRITE(stdout,'(5x,"\Y = ",10x,es20.2," H", 4x,es20.2)') baryerr, baryerr/bary
WRITE(stdout,'(5x,"\mu = ",9x,es20.2," J/T", 4x,es18.2)') barmuerr, &
                                                        barmuerr/barmu
WRITE(stdout,'(5x,"\M = ",10x,es20.2," A/m", 4x,es18.2)') barmagerr, &
                                                        barmagerr/barmag
WRITE(stdout,'(5x,"\H = ",10x,es20.2," A/m", 4x,es18.2)') barherr, &
                                                        barherr/barh
WRITE(stdout,*)
cspeedauerr=barverr/barv * cspeedau
amuauerr=(amuerr/amu+meerr/me) * amuau
WRITE(stdout,'(5x,"cspeed = ",6x,es20.2," \v",5x,es18.2)') cspeedauerr, &
                                              cspeedauerr/cspeedau 
WRITE(stdout,'(5x,"amu = ",9x,es20.2," \m",5x,es18.2)') amuauerr, &
                                              amuauerr/amuau 

WRITE(stdout,*)
kappaerr=0.5_DP * epsilon0err/epsilon0 * kappa
kappaaerr=kappaerr/kappa * kappaa
itoierr=kappaerr/kappa * itoi
rhotorhoerr=kappaerr/kappa*rhotorho
curtocurerr=kappaerr/kappa * curtocur
phitophierr=kappaerr/kappa * phitophi
etoeerr=kappaerr/kappa * etoe
captocaperr=2.0_DP * kappaerr / kappa * captocap
diptodiperr=kappaerr/kappa * diptodip
polartopolarerr=kappaerr/kappa * polartopolar
dtoderr=polartopolarerr/polartopolar * dtod
ohmtoohmerr=2.0_DP * kappaerr / kappa * ohmtoohm
btoberr=kappaerr/kappa*btob
avtoaverr=kappaerr/kappa*avtoav
wbtowberr=kappaerr/kappa * wbtowb
ytoyerr=2.0_DP * kappaerr/kappa * ytoy
mutomuerr=kappaerr/kappa * mutomu
magtomagerr=kappaerr/kappa * magtomag
htoherr=magtomagerr/magtomag * htoh

WRITE(stdout,'(/,5x,"Errors of conversion factors (c.g.s.-Gaussian - SI):")') 
WRITE(stdout,'(5x,"Current = ",5x,es20.2," A", 4x,es20.2)') itoierr, &
                                                           itoierr/itoi
WRITE(stdout,'(5x,"Charge = ",6x,es20.2," C", 4x,es20.2)') itoierr, &
                                                           itoierr/itoi
WRITE(stdout,'(5x,"Charge density = ",es18.2," C/m^3", 4x,es16.2)') &
                                      rhotorhoerr, rhotorhoerr/rhotorho
WRITE(stdout,'(5x,"Current density = ",es17.2," A/m^2", 4x,es16.2)') &
                                      curtocurerr, curtocurerr/curtocur
WRITE(stdout,'(5x,"Electric field = ",es18.2," V/m", 4x,es18.2)') etoeerr, &
                                      etoeerr/etoe
WRITE(stdout,'(5x,"Electric potential = ",es14.2," V", 4x,es20.2)') &
                                      phitophierr, phitophierr/phitophi
WRITE(stdout,'(5x,"Capacitance = ",3x,es18.2," F", 4x,es20.2)') captocaperr, &
                                      captocaperr/captocap
WRITE(stdout,'(5x,"Dipole moment = ",1x,es18.2," C m", 4x,es18.2)') &
                                      diptodiperr, diptodiperr/diptodip
WRITE(stdout,'(5x,"Electric Polarization = ",es11.2," C/m^2", 4x,es16.2)') &
                            polartopolarerr, polartopolarerr/polartopolar
WRITE(stdout,'(5x,"Electric Displacement = ",es11.2," C/m^2",4x, es16.2)') &
                                                 dtoderr, dtoderr/dtod
WRITE(stdout,'(5x,"Resistance = ",4x,es18.2," Ohm", 4x,es18.2)') ohmtoohmerr, &
                                              ohmtoohmerr/ohmtoohm
WRITE(stdout,'(5x,"Magnetic induction = ",es14.2," T", 4x,es20.2)') &
                                              btoberr, btoberr/btob
WRITE(stdout,'(5x,"Vector potential = ",es16.2," T m",2x,es20.2)') &
                                              avtoaverr, avtoaverr/avtoav
WRITE(stdout,'(5x,"Magnetic field flux = ",es13.2," Wb", 4x,es19.2)') &
                                             wbtowberr, wbtowberr/wbtowb
WRITE(stdout,'(5x,"Inductance = ",2x,es20.2," H", 4x,es20.2)') ytoyerr, &
                                               ytoyerr/ytoy
WRITE(stdout,'(5x,"Magnetic dipole = ",es17.2," J/T", 4x,es18.2)') &
                                         mutomuerr, mutomuerr/mutomu
WRITE(stdout,'(5x,"Magnetization = ",es19.2," A/m", 4x,es18.2)') &
                                           magtomagerr, magtomagerr/magtomag
WRITE(stdout,'(5x,"Magnetic strength = ",es15.2," A/m", 4x,es18.2)') &
                                           htoherr, htoherr/htoh

WRITE(stdout,*)

baricgserr=(barierr/bari+kappaerr/kappa)*baricgs
barccgserr=(kappaerr/kappa)*barccgs
barrhocgserr=(barrhoerr/barrho + rhotorhoerr/rhotorho)*barrhocgs
barcurcgserr=(barcurerr/barcur + curtocurerr/curtocur)*barcurcgs
barecgserr=(bareerr/bare+etoeerr/etoe)*barecgs
barphicgserr=(barphierr/barphi+phitophierr/phitophi)*barphicgs
barcapcgserr=(barcaperr/barcap + captocaperr/captocap) * barcapcgs
bardipcgserr=(bardiperr/bardip + diptodiperr/diptodip) * bardipcgs
barpolarcgserr=(barpolarerr/barpolar+polartopolarerr/polartopolar)*barpolarcgs
bardcgserr=(barderr/bard+dtoderr/dtod) * bardcgs
barohmcgserr=ohmtoohmerr/ohmtoohm * barohmcgs
barbcgserr=(barberr/barb + btoberr/btob) * barbcgs
baravcgserr=(baraverr/barav + avtoaverr/avtoav) * baravcgs
barwbcgserr=(wbtowberr/wbtowb)*barwbcgs
barycgserr=(baryerr/bary + ytoyerr/ytoy)*barycgs
barmucgserr=(barmuerr/barmu + mutomuerr/mutomu)*barmucgs
barmagcgserr=(barmagerr/barmag + magtomagerr/magtomag)*barmagcgs
barhcgserr=(barherr/barh + htoherr/htoh)*barhcgs

WRITE(stdout,'(/,5x,"Errors of conversion factors &
                                      &(atomic units - c.g.s.-Gaussian):")') 

WRITE(stdout,'(5x,"Current = ",5x,es20.2," statA", 4x,es16.2)')  &
                                   baricgserr, baricgserr/baricgs
WRITE(stdout,'(5x,"Charge = ",6x,es20.2," statC", 4x,es16.2)') &
                                   barccgserr, barccgserr/barccgs
WRITE(stdout,'(5x,"Charge density = ",es18.2," statC/cm^3", 4x,es11.2)') &
                                   barrhocgserr, barrhocgserr/barrhocgs
WRITE(stdout,'(5x,"Current density = ",es17.2," statA/cm^2", 4x,es11.2)') &
                                   barcurcgserr, barcurcgserr/barcurcgs
WRITE(stdout,'(5x,"Electric field = ",es18.2," statV/cm",4x,es13.2)') &
                                   barecgserr, barecgserr/barecgs
WRITE(stdout,'(5x,"Electric potential = ",es14.2," statV", 4x,es16.2)') &
                                   barphicgserr, barphicgserr/barphicgs
WRITE(stdout,'(5x,"Capacitance = ",2x,es19.2," cm", 4x,es19.2)') &
                                   barcapcgserr, barcapcgserr/barcapcgs
WRITE(stdout,'(5x,"Dipole moment = ",es19.2," statC cm", 4x,es13.2)') &
                                   bardipcgserr, bardipcgserr/bardipcgs
WRITE(stdout,'(5x,"Electric Polarization = ",es11.2," statC/cm^2", &
                      &4x,es11.2)') barpolarcgserr, barpolarcgserr/barpolarcgs
WRITE(stdout,'(5x,"Electric displacement = ",es11.2," statC/cm^2 4pi", &
                      &es11.2)') bardcgserr, bardcgserr/bardcgs
WRITE(stdout,'(5x,"Resistance = ",3x,es19.2," s/cm",4x,es17.2)') &
                                 barohmcgserr, barohmcgserr/barohmcgs
WRITE(stdout,'(5x,"Magnetic induction = ",es14.2," G",4x,es20.2)') &
                                 barbcgserr, barbcgserr/barbcgs
WRITE(stdout,'(5x,"Vector potential = ",es16.2," G cm",1x,es20.2)') &
                                 baravcgserr, baravcgserr/baravcgs
WRITE(stdout,'(5x,"Magnetic field flux = ",es13.2," Mx",3x,es20.2)') &
                                 barwbcgserr, barwbcgserr/barwbcgs
WRITE(stdout,'(5x,"Inductance = ",2x,es20.2," statH",es20.2)') &
                                 barycgserr, barycgserr/barycgs
WRITE(stdout,'(5x,"Magnetic dipole = ",es17.2," statC cm",4x,es13.2)') &
                                 barmucgserr, barmucgserr/barmucgs
WRITE(stdout,'(5x,"Magnetization = ",es19.2," statC/cm^2",es15.2)') &
                                 barmagcgserr, barmagcgserr/barmagcgs
WRITE(stdout,'(5x,"Magnetic strength = ",es15.2," statC/cm^2 4pi",es11.2)') &
                                 barhcgserr, barhcgserr/barhcgs

barmryerr= meerr / me * barmry
bartryerr= barterr / bart * bartry
barnuryerr=barnuerr / barnu * barnury
barvryerr= barverr / barv * barvry
bararyerr= baraerr / bara * barary
barfryerr= barferr / barf * barfry
baruryerr= hartreeerr / baru * barury
barwryerr= barwerr / barw * barwry
barprryerr= barprerr / barpr * barprry
bariryerr= barierr / bari * bariry
barrhoryerr=barrhoerr / barrho * barrhory
barcurryerr=barcurerr/barcur * barcurry
bareryerr=bareerr / bare * barery
barphiryerr=barphierr/ barphi * barphi 
bardipryerr=bardiperr / bardip * bardipry
barpolarryerr=barpolarerr/barpolar * barpolarry 
bardryerr=barderr / bard * bardry
barbryerr= barberr/barb * barbry
baravryerr= baraverr/barav * baravry
baryryerr= baryerr/bary * baryry
barmuryerr=barmuerr/barmu * barmury
barmagryerr=barmagerr / barmag * barmagry
barhryerr=barherr / barh * barhry

WRITE(stdout,'(/,5x,"Errors of conversion factors &
                                          &(Rydberg atomic units - SI):")') 
WRITE(stdout,'(5x,"\l_R = ",8x,es20.2," m", 4x,es20.2)') abohrerr, &
                                                          abohrerr/barl
WRITE(stdout,'(5x,"\m_R = ",8x,es20.2," kg", 3x,es20.2)') barmryerr, &
                                                          barmryerr/barmry
WRITE(stdout,'(5x,"\t_R = ",8x,es20.2," s", 4x,es20.2)') bartryerr, &
                                                          bartryerr/bart
WRITE(stdout,'(5x,"\nu_R = ",7x,es20.2," s", 4x,es20.2)') barnuryerr, &
                                                          barnuryerr/barnury
WRITE(stdout,'(5x,"\v_R = ",8x,es20.2," m/s", 4x,es18.2)') barvryerr, &
                                                          barvryerr/barvry
WRITE(stdout,'(5x,"\a_R = ",8x,es20.2," m/s^2", 4x,es16.2)') bararyerr, &
                                                          bararyerr/barary
WRITE(stdout,'(5x,"\p_R = ",8x,es20.2," kg m/s", 4x,es15.2)') barperr, &
                                                          barperr/barp
WRITE(stdout,'(5x,"\L_R = ",8x,es20.2," kg m^2/s",es17.2)') zero, zero
WRITE(stdout,'(5x,"\f_R = ",8x,es20.2," N", 4x,es20.2)') barfryerr, &
                                                          barfryerr/barfry
WRITE(stdout,'(5x,"\U_R = ",8x,es20.2," J", 4x,es20.2)') baruryerr, &
                                                          baruryerr/barury
WRITE(stdout,'(5x,"\W_R = ",8x,es20.2," W", 4x,es20.2)') barwryerr, &
                                                          barwryerr/barwry
WRITE(stdout,'(5x,"\pr_R = ",7x,es20.2," Pa", 4x,es19.2)') barprryerr, &
                                                          barprryerr/barprry
WRITE(stdout,'(5x,"\I_R = ",8x,es20.2," A", 4x,es20.2)') bariryerr, &
                                                          bariryerr/bariry
WRITE(stdout,'(5x,"\C_R = ",8x,es20.2," C", 4x,es20.2)') zero, zero
WRITE(stdout,'(5x,"\rho_R = ",6x,es20.2," C/m^3", 4x,es16.2)') barrhoryerr, &
                                                          barrhoryerr/barrhory
WRITE(stdout,'(5x,"\J_R = ",9x,es19.2," A/m^2", 4x,es16.2)') barcurryerr, &
                                                          barcurryerr/barcurry
WRITE(stdout,'(5x,"\E_R = ",8x,es20.2," V/m", 4x,es18.2)') bareryerr, &
                                                          bareryerr/barery
WRITE(stdout,'(5x,"\V_R = ",8x,es20.2," V", 4x,es20.2)') barphiryerr, &
                                                          barphiryerr/barphiry
WRITE(stdout,'(5x,"\F_R = ",8x,es20.2," F", 4x,es20.2)') barcaperr, &
                                                        barcaperr/barcap
WRITE(stdout,'(5x,"\dip_R = ",4x,es22.2," C m", 4x,es18.2)') bardipryerr, &
                                                        bardipryerr/bardipry
WRITE(stdout,'(5x,"\P_R = ",8x,es20.2," C/m^2", 4x,es16.2)') barpolarryerr, &
                                                     barpolarryerr/barpolarry
WRITE(stdout,'(5x,"\D_R = ",8x,es20.2," C/m^2",4x, es16.2)') bardryerr, &
                                                        bardryerr/bardry
WRITE(stdout,'(5x,"\R_R = ",8x,es20.2," Ohm", 2x,es20.2)') zero, zero
WRITE(stdout,'(5x,"\B_R = ",8x,es20.2," T", 4x,es20.2)') barbryerr, &
                                                        barbryerr/barbry
WRITE(stdout,'(5x,"\A_R = ",8x,es20.2," T m", 2x,es20.2)') baravryerr, &
                                                          baravryerr/baravry
WRITE(stdout,'(5x,"\phi_R = ",6x,es20.2," Wb", 3x,es20.2)') zero, zero
WRITE(stdout,'(5x,"\Y_R = ",8x,es20.2," H", 4x,es20.2)') baryryerr, &
                                                          baryryerr/baryry
WRITE(stdout,'(5x,"\mu_R = ",7x,es20.2," J/T", 4x,es18.2)') barmuryerr, &
                                                        barmuryerr/barmury
WRITE(stdout,'(5x,"\M_R = ",8x,es20.2," A/m", 4x,es18.2)') barmagryerr, &
                                                        barmagryerr/barmagry
WRITE(stdout,'(5x,"\H_R = ",8x,es20.2," A/m", 4x,es18.2)') barhryerr, &
                                                        barhryerr/barhry
WRITE(stdout,*)

barbgerr=(barberr/barb+alphaerr/alphaf)*barbg
baravgerr=(baraverr/barav+alphaerr/alphaf)*baravg
barwbgerr=alphaerr/alphaf * barwbg
barmugerr=(barmuerr/barmu+alphaerr/alphaf)*barmug
barmaggerr=(barmagerr/barmag+alphaerr/alphaf)*barmagg
barhgerr=(barherr/barh+alphaerr/alphaf)*barhg


WRITE(stdout,'(/,5x,"Errors of conversion factors &
                                          &(Gaussian atomic units - SI):")') 
WRITE(stdout,'(5x,"\B_G = ",8x,es20.2," T", 4x,es20.2)') barbgerr, &
                                                        barbgerr/barbg
WRITE(stdout,'(5x,"\A_G = ",8x,es20.2," T m", 2x,es20.2)') baravgerr, &
                                                          baravgerr/baravg
WRITE(stdout,'(5x,"\phi_G = ",6x,es20.2," Wb", 3x,es20.2)') barwbgerr, &
                                                          barwbgerr/barwbg
WRITE(stdout,'(5x,"\mu_G = ",7x,es20.2," J/T", 4x,es18.2)') barmugerr, &
                                                        barmugerr/barmug
WRITE(stdout,'(5x,"\M_G = ",8x,es20.2," A/m", 4x,es18.2)') barmaggerr, &
                                                        barmaggerr/barmagg
WRITE(stdout,'(5x,"\H_G = ",8x,es20.2," A/m", 4x,es18.2)') barhgerr, &
                                                        barhgerr/barhg
WRITE(stdout,*)




IF (ionode) CLOSE(UNIT=iunout, status='KEEP')

CALL environment_end( code )
CALL mp_global_end ()
END PROGRAM units

FUNCTION float_to_latex(a,n)
USE kinds, ONLY : DP
IMPLICIT NONE
CHARACTER(LEN=80) :: float_to_latex
REAL(DP) :: a
INTEGER :: n
CHARACTER(LEN=80) :: aux, formt
CHARACTER(LEN=3) :: expon
CHARACTER(LEN=76) :: mantis
CHARACTER(LEN=6) :: int_to_char
INTEGER :: n1

formt="(es80."//TRIM(int_to_char(n))//")"
WRITE(aux,formt) a
expon=aux(78:80)
!
!  correct for the presence of a zero in the exponent
!
IF (expon(2:2)=='0') THEN
   expon(2:2)=expon(3:3)
   expon(3:3)=' '
ENDIF
!
! remove the + in the exponent
!
IF (expon(1:1)=='+') THEN
   expon(1:1)=expon(2:2)
   expon(2:2)=expon(3:3)
   expon(3:3)=' '
ENDIF

n1=76
mantis=TRIM(ADJUSTL(aux(1:n1)))
IF (TRIM(expon)/='0') THEN
   float_to_latex=TRIM(mantis)//"\times 10^{"//TRIM(expon)//"}"
ELSE
   float_to_latex=TRIM(mantis)
ENDIF

END FUNCTION float_to_latex
