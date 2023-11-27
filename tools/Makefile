# Makefile for thermo_pw/tools

include ../../make.inc

# location of needed modules and included files 
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib \
          $(MOD_FLAG)../../dft-d3/ \
          $(MOD_FLAG). 

#libraries of QE
PWOBJS = ../../PHonon/PH/libph.a \
         ../../PP/src/libpp.a \
         ../../LR_Modules/liblrmod.a \
         ../../PW/src/libpw.a 

#library modules of QE
QEMODS = ../../Modules/libqemod.a \
         ../../KS_Solvers/libks_solvers.a \
         ../../dft-d3/libdftd3qe.a

#libraries of thermo_pw
LIBTHERMO= ../lapack/lapack_qe.a \
           ../fftpack5.1/fftpack5.1.a  \
           ../lib/thermo_lib.a \
           ../qe/thermo_qe.a 

MODULES = $(LIBTHERMO) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(BASEMODS)

all : space_groups.x crystal_point_group.x mag_point_group.x units.x \
      hex_trig.x bravais_lattices.x kovalev.x supercell.x translate.x \
      elastic.x plot_sur_states.x gener_2d_slab.x gener_3d_slab.x \
      gener_nanowire.x debye.x rotate_tensors.x test_colors.x epsilon_tpw.x \
      pdec.x optical.x change_e_name.x change_tel_name.x \
      change_dynmat_name.x merge_interp.x template.x test_eos.x emp_f.x \
      emp_g.x average_grun.x density.x
      

space_groups.x : space_groups.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ space_groups.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/space_groups.x . )

crystal_point_group.x : crystal_point_group.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ crystal_point_group.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/crystal_point_group.x . )


mag_point_group.x : mag_point_group.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ mag_point_group.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/mag_point_group.x . )

units.x : units.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ units.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/units.x . )

hex_trig.x : hex_trig.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ hex_trig.o 
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/hex_trig.x . )

bravais_lattices.x : bravais_lattices.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ bravais_lattices.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/bravais_lattices.x . )

kovalev.x : kovalev.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ kovalev.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/kovalev.x . )

supercell.x : supercell.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ supercell.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/supercell.x . )

translate.x : translate.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ translate.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/translate.x . )

elastic.x : elastic.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ elastic.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/elastic.x . )

optical.x : optical.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ optical.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/optical.x . )

plot_sur_states.x : plot_sur_states.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ plot_sur_states.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/plot_sur_states.x . )

gener_2d_slab.x : gener_2d_slab.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ gener_2d_slab.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_2d_slab.x . )

gener_3d_slab.x : gener_3d_slab.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ gener_3d_slab.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_3d_slab.x . )

gener_nanowire.x : gener_nanowire.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ gener_nanowire.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/gener_nanowire.x . )

debye.x : debye.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ debye.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/debye.x . )

rotate_tensors.x : rotate_tensors.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ rotate_tensors.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/rotate_tensors.x . )

test_colors.x : test_colors.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ test_colors.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/test_colors.x . )

epsilon_tpw.x : epsilon_tpw.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ epsilon_tpw.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/epsilon_tpw.x . )

pdec.x : pdec.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ pdec.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/pdec.x . )

change_dynmat_name.x : change_dynmat_name.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ change_dynmat_name.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/change_dynmat_name.x . )

change_e_name.x : change_e_name.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ change_e_name.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/change_e_name.x . )

change_tel_name.x : change_tel_name.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ change_tel_name.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/change_tel_name.x . )

merge_interp.x : merge_interp.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ merge_interp.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/merge_interp.x . )

template.x : template.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ template.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/template.x . )

test_eos.x : test_eos.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ test_eos.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/test_eos.x . )

emp_f.x : emp_f.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ emp_f.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/emp_f.x . )

emp_g.x : emp_g.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ emp_g.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/emp_g.x . )

average_grun.x : average_grun.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ average_grun.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/average_grun.x . )

density.x : density.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ density.o $(MODULES) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../thermo_pw/tools/density.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/space_groups.x 
	- /bin/rm -rf ../../bin/crystal_point_group.x 
	- /bin/rm -rf ../../bin/mag_point_group.x 
	- /bin/rm -rf ../../bin/units.x 
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
	- /bin/rm -rf ../../bin/debye.x 
	- /bin/rm -rf ../../bin/rotate_tensors.x 
	- /bin/rm -rf ../../bin/test_colors.x 
	- /bin/rm -rf ../../bin/epsilon_tpw.x 
	- /bin/rm -rf ../../bin/pdec.x 
	- /bin/rm -rf ../../bin/optical.x 
	- /bin/rm -rf ../../bin/change_dynmat_name.x 
	- /bin/rm -rf ../../bin/change_e_name.x 
	- /bin/rm -rf ../../bin/change_tel_name.x 
	- /bin/rm -rf ../../bin/test_eos.x 
	- /bin/rm -rf ../../bin/emp_f.x 
	- /bin/rm -rf ../../bin/emp_g.x 
	- /bin/rm -rf ../../bin/average_grun.x 
	- /bin/rm -rf ../../bin/merge_interp.x 
	- /bin/rm -rf ../../bin/density.x 
	- /bin/rm -rf ../../bin/template.x 

include make.depend
# DO NOT DELETE
