# Makefile for thermo_pw 

include ../../make.sys

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib $(MOD_FLAG).
IFLAGS=

PWOBJS = ../../PHonon/PH/libph.a ../../PP/src/libpp.a ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a

LIBOBJS	= ../lib/thermo_lib.a ../../PHonon/PH/libph.a ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

TLDEPS= bindir mods libs pw

all : group_name.x group_number.x gener_2d_slab.x gener_nanowire.x supercell.x

libs-ph : libph.a libphaux.a

group_name.x : group_name.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ group_name.o  \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/group_name.x . )

group_number.x : group_number.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ group_number.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/group_number.x . )

supercell.x : supercell.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ supercell.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/supercell.x . )

translate.x : translate.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ translate.o \
	         $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/translate.x . )

gener_2d_slab.x : gener_2d_slab.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_2d_slab.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_2d_slab.x . )

gener_nanowire.x : gener_nanowire.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_nanowire.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_nanowire.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/thermo_pw.x 

include make.depend
# DO NOT DELETE
