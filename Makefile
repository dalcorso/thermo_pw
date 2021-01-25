# Makefile for thermo_pw
sinclude ../make.sys

default: all

all: thermo_pw thermo_tools 

thermo_tools: thermo_lib
	( cd tools ; $(MAKE) all || exit 1 )

thermo_pw: thermo_lib thermo_qe
	( cd src ; $(MAKE) all || exit 1 )

thermo_lib: 
	( cd lib ; $(MAKE) all || exit 1 )

thermo_qe: 
	( cd qe ; $(MAKE) all || exit 1 )

join_qe:
	if test -f main_Makefile ; then mv ../Makefile Makefile_qe ; \
           mv main_Makefile ../Makefile ; fi
	if test -f plugins_makefile ; then mv ../install/plugins_makefile \
           plugins_makefile_qe ; mv plugins_makefile ../install ; fi
	if test -f makedeps.sh ; then mv ../install/makedeps.sh \
           makedeps.sh_qe ; mv makedeps.sh ../install/ ; fi
	if test -f ph_restart.f90 ; then mv ../PHonon/PH/ph_restart.f90 \
           ph_restart.f90_qe ; mv ph_restart.f90 ../PHonon/PH ; fi

leave_qe:
	if test -f Makefile_qe ; then mv ../Makefile main_Makefile ; \
           mv Makefile_qe ../Makefile ; fi
	if test -f plugins_makefile_qe ; then \
           mv ../install/plugins_makefile . ; \
           mv plugins_makefile_qe ../install/plugins_makefile ; fi
	if test -f makedeps.sh_qe ; then mv ../install/makedeps.sh . ; \
           mv makedeps.sh_qe ../install/makedeps.sh ; fi
	if test -f ph_restart.f90_qe ; then mv ../PHonon/PH/ph_restart.f90 \
           ph_restart.f90 ; mv ph_restart.f90_qe ../PHonon/PH/ph_restart.f90 ; fi

tpw_gpu:
	if test -f cpu_gpu/gpu_setup.f90 ; then mv qe/setup.f90 \
           cpu_gpu/cpu_setup.f90 ; \
           mv cpu_gpu/gpu_setup.f90 qe/setup.f90 ; fi
	if test -f cpu_gpu/gpu_electrons.f90 ; then mv qe/electrons.f90 \
           cpu_gpu/cpu_electrons.f90 ; \
           mv cpu_gpu/gpu_electrons.f90 qe/electrons.f90 ; fi
	if test -f cpu_gpu/gpu_non_scf.f90 ; then mv qe/non_scf.f90 \
           cpu_gpu/cpu_non_scf.f90 ; \
           mv cpu_gpu/gpu_non_scf.f90 qe/non_scf.f90 ; fi
	if test -f cpu_gpu/gpu_c_bands.f90 ; then mv qe/c_bands.f90 \
           cpu_gpu/cpu_c_bands.f90 ; \
           mv cpu_gpu/gpu_c_bands.f90 qe/c_bands.f90 ; fi
	if test -f cpu_gpu/gpu_stress.f90 ; then mv qe/stress.f90 \
           cpu_gpu/cpu_stress.f90 ; \
           mv cpu_gpu/gpu_stress.f90 qe/stress.f90 ; fi

tpw_cpu:
	if test -f cpu_gpu/cpu_setup.f90 ; then mv qe/setup.f90 \
           cpu_gpu/gpu_setup.f90 ; \
           mv cpu_gpu/cpu_setup.f90 qe/setup.f90 ; fi
	if test -f cpu_gpu/cpu_electrons.f90 ; then mv qe/electrons.f90 \
           cpu_gpu/gpu_electrons.f90 ; \
           mv cpu_gpu/cpu_electrons.f90 qe/electrons.f90 ; fi
	if test -f cpu_gpu/cpu_non_scf.f90 ; then mv qe/non_scf.f90 \
           cpu_gpu/gpu_non_scf.f90 ; \
           mv cpu_gpu/cpu_non_scf.f90 qe/non_scf.f90 ; fi
	if test -f cpu_gpu/cpu_c_bands.f90 ; then mv qe/c_bands.f90 \
           cpu_gpu/gpu_c_bands.f90 ; \
           mv cpu_gpu/cpu_c_bands.f90 qe/c_bands.f90 ; fi
	if test -f cpu_gpu/cpu_stress.f90 ; then mv qe/stress.f90 \
           cpu_gpu/gpu_stress.f90 ; \
           mv cpu_gpu/cpu_stress.f90 qe/stress.f90 ; fi


clean: thermo_tools_clean thermo_pw_clean thermo_lib_clean thermo_qe_clean examples_clean examples_qe_clean doc_clean

thermo_pw_clean:
	( cd src ; $(MAKE) clean )

thermo_tools_clean:
	( cd tools ; $(MAKE) clean )

thermo_lib_clean:
	( cd lib ; $(MAKE) clean )

thermo_qe_clean:
	( cd qe ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi

examples_qe_clean:
	if test -d examples_qe ; then \
	( cd examples_qe ; ./clean_all ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi
doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi


distclean: clean 
