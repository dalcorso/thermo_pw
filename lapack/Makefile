# Makefile for thermo_pw/qe

include ../../make.inc

ifneq ($(strip $(CUDA_F90FLAGS)),)
ifneq ($(findstring Mcuda, $(CUDA_F90FLAGS)),)
.f90.o:
	$(MPIF90) -Mcuda=charstring $(F90FLAGS) -c $<
else
ifneq ($(findstring cuda, $(CUDA_F90FLAGS)),)
.f90.o:
	$(MPIF90) -cuda=charstring $(F90FLAGS) -c $<
endif
endif
endif

MODULES=      \
la_constants_dev.o \
dcabs1_dev.o \
dcopy_dev.o \
disnan_dev.o \
dladiv1_dev.o \
dladiv2_dev.o \
dladiv_dev.o \
dlae2_dev.o \
dlaebz_dev.o \
dlaev2_dev.o \
dlagtf_dev.o \
dlagts_dev.o \
dlaisnan_dev.o \
dlamc3_dev.o \
dlamch_dev.o \
dlanst_dev.o \
dlapy2_dev.o \
dlapy3_dev.o \
dlarnv_dev.o \
dlartg_dev.o \
dlaruv_dev.o \
dlascl_dev.o \
dlasrt_dev.o \
dlassq_dev.o \
dnrm2_dev.o \
dscal_dev.o \
dsterf_dev.o \
dstebz_dev.o \
dznrm2_dev.o \
idamax_dev.o \
ieeeck_dev.o \
ilaenv_dev.o \
ilazlc_dev.o \
ilazlr_dev.o \
iparmq_dev.o \
la_xisnan_dev.o \
lsame_dev.o \
xerbla_dev.o \
zaxpy_dev.o \
zcopy_dev.o \
zdotc_dev.o \
zdscal_dev.o \
zgemm_dev.o \
zgemv_dev.o \
zgerc_dev.o \
zheev_dev.o \
zhegs2_dev.o \
zhegst_dev.o \
zhegv_dev.o \
zheevx_dev.o \
zhegvx_dev.o \
zhemm_dev.o \
zhemv_dev.o \
zher2_dev.o \
zher2k_dev.o \
zherk_dev.o \
zhetd2_dev.o \
zhetrd_dev.o \
zlacgv_dev.o \
zlacpy_dev.o \
zladiv_dev.o \
zlanhe_dev.o \
zlarf_dev.o \
zlarfb_dev.o \
zlarfg_dev.o \
zlarft_dev.o \
zlartg_dev.o \
zlascl_dev.o \
zlaset_dev.o \
zlasr_dev.o \
zlassq_dev.o \
zlatrd_dev.o \
zpotrf2_x_dev.o \
zpotrf_dev.o \
zscal_dev.o \
zstein_dev.o \
zsteqr_dev.o \
zswap_dev.o \
ztrmm_dev.o \
ztrmv_dev.o \
ztrsm_dev.o \
ztrsv_dev.o \
zung2l_dev.o \
zung2r_dev.o \
zungql_dev.o \
zungqr_dev.o \
zungtr_dev.o \
zunm2l_dev.o \
zunm2r_dev.o \
zunmql_dev.o \
zunmqr_dev.o \
zunmtr_dev.o 


all : lapack_qe.a

lapack_qe.a: $(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

include make.depend

