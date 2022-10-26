# Makefile for thermo_pw
sinclude ../make.inc

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
	if test -f ./main_Makefile ; then mv ../Makefile Makefile_qe ; \
           mv ./main_Makefile ../Makefile ; fi
	if test -f ./plugins_makefile ; then mv ../install/plugins_makefile \
           ./plugins_makefile_qe ; mv ./plugins_makefile ../install ; fi
	if test -f ./makedeps.sh ; then mv ../install/makedeps.sh \
           ./makedeps.sh_qe ; mv ./makedeps.sh ../install/ ; fi
	if test -f ./main_CMakeLists.txt ; then mv ../CMakeLists.txt \
          ./CMakeLists.txt_qe ; mv ./main_CMakeLists.txt ../CMakeLists.txt ; fi
	if test -f ./QESUB/lmdif.f90 ; then mv ../Modules/lmdif.f90 \
           ./QESUB/lmdif.f90_qe ; mv ./QESUB/lmdif.f90 ../Modules ; fi
	if test -f ./QESUB/v_of_rho.f90 ; then mv ../PW/src/v_of_rho.f90 \
           ./QESUB/v_of_rho.f90_qe ; \
           mv ./QESUB/v_of_rho.f90 ../PW/src/v_of_rho.f90 ; fi
	if test -f ./QESUB/write_upf_new.f90 ; then \
           mv ../upflib/write_upf_new.f90 ./QESUB/write_upf_new.f90_qe ; \
           mv ./QESUB/write_upf_new.f90 ../upflib/write_upf_new.f90 ; fi
	if test -f ./QESUB/qexsd_init.f90 ; then \
           mv ../Modules/qexsd_init.f90 ./QESUB/qexsd_init.f90_qe ; \
           mv ./QESUB/qexsd_init.f90 ../Modules/qexsd_init.f90 ; fi
leave_qe:
	if test -f ./Makefile_qe ; then mv ../Makefile ./main_Makefile ; \
           mv ./Makefile_qe ../Makefile ; fi
	if test -f ./plugins_makefile_qe ; then \
           mv ../install/plugins_makefile . ; \
           mv ./plugins_makefile_qe ../install/plugins_makefile ; fi
	if test -f ./makedeps.sh_qe ; then mv ../install/makedeps.sh . ; \
           mv ./makedeps.sh_qe ../install/makedeps.sh ; fi
	if test -f ./CMakeLists.txt_qe ; then mv ../CMakeLists.txt \
           ./main_CMakeLists.txt ; \
           mv ./CMakeLists.txt_qe ../CMakeLists.txt ; fi
	if test -f ./QESUB/lmdif.f90_qe ; then mv ../Modules/lmdif.f90 \
           ./QESUB/lmdif.f90 ; \
           mv ./QESUB/lmdif.f90_qe ../Modules/lmdif.f90 ; fi
	if test -f ./QESUB/v_of_rho.f90_qe ; then mv ../PW/src/v_of_rho.f90 \
           ./QESUB/v_of_rho.f90 ; \
           mv ./QESUB/v_of_rho.f90_qe ../PW/src/v_of_rho.f90 ; fi
	if test -f ./QESUB/write_upf_new.f90_qe ; then \
           mv ../upflib/write_upf_new.f90 ./QESUB/write_upf_new.f90 ; \
           mv ./QESUB/write_upf_new.f90_qe ../upflib/write_upf_new.f90 ; fi
	if test -f ./QESUB/qexsd_init.f90_qe ; then \
           mv ../Modules/qexsd_init.f90 ./QESUB/qexsd_init.f90 ; \
           mv ./QESUB/qexsd_init.f90_qe ../Modules/qexsd_init.f90 ; fi

clean: thermo_tools_clean thermo_pw_clean thermo_lib_clean thermo_qe_clean examples_clean examples_qe_clean space_groups_clean doc_clean

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

space_groups_clean:
	if test -d space_groups ; then \
	( cd space_groups/examples ; ./clean_all ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi
doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi


distclean: clean 
