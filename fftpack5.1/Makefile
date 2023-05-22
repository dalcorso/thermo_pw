# Makefile for thermo_pw/qe

include ../../make.inc

MODULES=      \
   c1f2kb_dev.o  \
   c1f2kf_dev.o  \
   c1f3kb_dev.o  \
   c1f3kf_dev.o  \
   c1f4kb_dev.o  \
   c1f4kf_dev.o  \
   c1f5kb_dev.o  \
   c1f5kf_dev.o  \
   c1fgkb_dev.o  \
   c1fgkf_dev.o  \
   c1fm1b_dev.o  \
   c1fm1f_dev.o  \
   cfft1b_dev.o  \
   cfft1f_dev.o  \
   cfft1i_dev.o  \
   r8_factor_dev.o  \
   r8_mcfti1_dev.o  \
   r8_tables_dev.o  \
   xerfft_dev.o  

all : fftpack5.1.a

fftpack5.1.a: $(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

