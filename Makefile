# Makefile for PH
sinclude ../make.sys

default: all

#all: phonon phgamma_only third_order third_order_q
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
	mv ../Makefile Makefile_qe
	mv main_Makefile ../Makefile
	mv ../install/plugins_makefile plugins_makefile_qe
	mv plugins_makefile ../install
	mv ../install/makedeps.sh makedeps.sh_qe
	mv makedeps.sh ../install/

leave_qe:
	mv ../Makefile main_Makefile
	mv Makefile_qe ../Makefile
	mv ../install/plugins_makefile .
	mv plugins_makefile_qe ../install/plugins_makefile
	mv ../install/makedeps.sh .
	mv makedeps.sh_qe ../install/makedeps.sh

clean: thermo_tools_clean thermo_pw_clean thermo_lib_clean thermo_qe_clean examples_clean doc_clean

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

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi
doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi


distclean: clean doc_clean
