# Makefile for PH
sinclude ../make.sys

default: all

#all: phonon phgamma_only third_order third_order_q
all: thermo_pw thermo_red

thermo_red: thermo_lib
	( cd src ; $(MAKE) all || exit 1 )

thermo_pw: thermo_lib
	( cd src ; $(MAKE) all || exit 1 )

thermo_lib: 
	( cd lib ; $(MAKE) all || exit 1 )

join_qe:
	mv ../Makefile Makefile_qe
	mv main_Makefile ../Makefile
	mv ../install/plugins_makefile plugins_makefile_qe
	mv plugins_makefile ../install

leave_qe:
	mv ../Makefile main_Makefile
	mv Makefile_qe ../Makefile
	mv ../install/plugins_makefile .
	mv plugins_makefile_qe ../install/plugins_makefile

clean: thermo_pw_clean

thermo_pw_clean:
	( cd src ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
