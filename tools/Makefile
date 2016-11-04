# Makefile for thermo_pw 

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib  \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH $(MOD_FLAG)../lib  $(MOD_FLAG).

IFLAGS=

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PHonon/PH/libph.a ../../PP/src/libpp.a ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

LIBOBJS	= ../lib/thermo_lib.a ../qe/thermo_qe.a ../../PHonon/PH/libph.a ../../clib/clib.a ../../iotk/src/libiotk.a

TLDEPS= bindir mods libs pw

all : group_name.x group_number.x crystal_point_group.x gener_2d_slab.x gener_3d_slab.x gener_nanowire.x supercell.x plot_sur_states.x translate.x elastic.x test_colors.x epsilon_tpw.x kovalev.x hex_trig.x bravais_lattices.x

libs-ph : libph.a libphaux.a

group_name.x : group_name.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ group_name.o  \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/group_name.x . )

group_number.x : group_number.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ group_number.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/group_number.x . )

crystal_point_group.x : crystal_point_group.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ crystal_point_group.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/crystal_point_group.x . )

hex_trig.x : hex_trig.o
	$(LD) $(LDFLAGS) -o $@ hex_trig.o 
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/hex_trig.x . )

bravais_lattices.x : bravais_lattices.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ bravais_lattices.o \
                 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/bravais_lattices.x . )

kovalev.x : kovalev.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ kovalev.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/kovalev.x . )

supercell.x : supercell.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ supercell.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/supercell.x . )

translate.x : translate.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ translate.o \
	         $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/translate.x . )

elastic.x : elastic.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ elastic.o \
	         $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/elastic.x . )

plot_sur_states.x : plot_sur_states.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ plot_sur_states.o  \
	$(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/plot_sur_states.x . )


gener_2d_slab.x : gener_2d_slab.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_2d_slab.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_2d_slab.x . )

gener_3d_slab.x : gener_3d_slab.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_3d_slab.o  \
	$(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_3d_slab.x . )

gener_nanowire.x : gener_nanowire.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_nanowire.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_nanowire.x . )

test_colors.x : test_colors.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ test_colors.o \
	         $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/test_colors.x . )

epsilon_tpw.x : epsilon_tpw.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLtAGS) -o $@ epsilon_tpw.o \
	         $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LRMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/epsilon_tpw.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/thermo_pw.x 
	- /bin/rm -rf ../../bin/group_name.x 
	- /bin/rm -rf ../../bin/group_number.x 
	- /bin/rm -rf ../../bin/supercell.x 
	- /bin/rm -rf ../../bin/translate.x 
	- /bin/rm -rf ../../bin/elastic.x 
	- /bin/rm -rf ../../bin/plot_sur_states.x 
	- /bin/rm -rf ../../bin/gener_2d_slab.x 
	- /bin/rm -rf ../../bin/gener_3d_slab.x 
	- /bin/rm -rf ../../bin/gener_nanowire.x 
	- /bin/rm -rf ../../bin/test_colors.x 
	- /bin/rm -rf ../../bin/epsilon_tpw.x 

include make.depend
# DO NOT DELETE
