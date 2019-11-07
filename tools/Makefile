# Makefile for thermo_pw/tools

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib 

IFLAGS=

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PHonon/PH/libph.a ../../PP/src/libpp.a ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/Davidson/libdavid.a \
         ../../KS_Solvers/CG/libcg.a ../../FFTXlib/libqefft.a \
         ../../LAXlib/libqela.a ../../UtilXlib/libutil.a
LIBOBJS	= ../lib/thermo_lib.a ../qe/thermo_qe.a ../../PHonon/PH/libph.a \
          ../../PW/src/libpw.a ../../clib/clib.a ../../iotk/src/libiotk.a

all : crystal_point_group.x mag_point_group.x gener_2d_slab.x gener_3d_slab.x \
      gener_nanowire.x supercell.x plot_sur_states.x translate.x elastic.x \
      test_colors.x epsilon_tpw.x kovalev.x hex_trig.x bravais_lattices.x \
      space_groups.x debye.x units.x

space_groups.x : space_groups.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ space_groups.o \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/space_groups.x . )

crystal_point_group.x : crystal_point_group.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ crystal_point_group.o \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/crystal_point_group.x . )

mag_point_group.x : mag_point_group.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ mag_point_group.o \
	    $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/mag_point_group.x . )

units.x : units.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ units.o \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/units.x . )

hex_trig.x : hex_trig.o
	$(LD) $(LDFLAGS) -o $@ hex_trig.o 
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/hex_trig.x . )

bravais_lattices.x : bravais_lattices.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ bravais_lattices.o \
        $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/bravais_lattices.x . )

kovalev.x : kovalev.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ kovalev.o \
		 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/kovalev.x . )

supercell.x : supercell.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ supercell.o \
        $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/supercell.x . )

translate.x : translate.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ translate.o \
	         $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/translate.x . )

elastic.x : elastic.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ elastic.o \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/elastic.x . )

plot_sur_states.x : plot_sur_states.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ plot_sur_states.o  \
	$(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/plot_sur_states.x . )

gener_2d_slab.x : gener_2d_slab.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_2d_slab.o \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_2d_slab.x . )

gener_3d_slab.x : gener_3d_slab.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_3d_slab.o  \
	$(PWOBJS) $(QEMODS) $(LIBOBJS) $(QEMODS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_3d_slab.x . )

gener_nanowire.x : gener_nanowire.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ gener_nanowire.o \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_nanowire.x . )

debye.x : debye.o $(PWOBJS) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ debye.o  \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/debye.x . )

rotate_tensors.x : rotate_tensors.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ rotate_tensors.o \
	        $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/rotate_tensors.x . )

test_colors.x : test_colors.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ test_colors.o \
	 $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/test_colors.x . )

epsilon_tpw.x : epsilon_tpw.o $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ epsilon_tpw.o \
        $(LIBOBJS) $(PWOBJS) $(QEMODS) $(LRMODS) $(LIBOBJS) $(PWOBJS) $(LIBS) \
        $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/epsilon_tpw.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/space_groups.x 
	- /bin/rm -rf ../../bin/crystal_point_group.x 
	- /bin/rm -rf ../../bin/mag_point_group.x 
	- /bin/rm -rf ../../bin/hex_trig.x 
	- /bin/rm -rf ../../bin/bravais_lattices.x 
	- /bin/rm -rf ../../bin/kovalev.x 
	- /bin/rm -rf ../../bin/supercell.x 
	- /bin/rm -rf ../../bin/translate.x 
	- /bin/rm -rf ../../bin/elastic.x 
	- /bin/rm -rf ../../bin/plot_sur_states.x 
	- /bin/rm -rf ../../bin/gener_2d_slab.x 
	- /bin/rm -rf ../../bin/gener_3d_slab.x 
	- /bin/rm -rf ../../bin/gener_nanowire.x 
	- /bin/rm -rf ../../bin/rotate_tensors.x 
	- /bin/rm -rf ../../bin/test_colors.x 
	- /bin/rm -rf ../../bin/epsilon_tpw.x 

include make.depend
# DO NOT DELETE
