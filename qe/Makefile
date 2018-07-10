# Makefile for thermo_pw/qe

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS)  \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../lib  \
          $(MOD_FLAG).

MODULES = \
add_dkmds.o \
addnlcc_zstar_eu_us.o \
adddvepsi_us.o \
adddvscf_tran.o \
addusddenseq.o \
add_zstar_ue.o \
add_zstar_us.o \
allocate_cg.o \
allocate_phq.o \
allocate_lanczos.o \
apply_ac.o \
bcast_qes_module.o \
c_bands.o \
calc_chi.o \
check_bands.o \
check_initial_status.o \
check_initial_geometry.o \
clean_all_pw.o \
collect_grid_files.o \
compute_intq.o \
compute_drhous.o \
compute_drhous_nc.o \
control_qe.o \
ccg_psi.o \
deallocate_ph_tpw.o \
divide_et_impera.o \
do_cg_e.o \
do_cg_ph.o \
do_lanczos.o \
do_lanczos_psh.o \
drho.o \
dveqpsi_us.o \
dveqpsi_us_only.o \
dv_of_drho_clf.o \
dv_of_drho_tran.o \
electrons.o \
ev_xml.o \
find_band_sym.o \
find_mode_sym.o \
h_pcg_step.o \
initialize_ph.o \
kpoint_grid.o \
irrek.o \
lanczos_write_restart.o \
lr_dot.o \
lr_restart.o \
lr_sm1_psi.o \
lr_us.o  \
mix_pot.o \
nh_lanczos_step.o \
non_scf.o \
openfilq.o \
orthogonalize_omega.o \
paw_add_onecenter.o \
paw_add_symmetry.o \
phescf.o \
phqscf.o \
pheqscf.o \
phq_init.o \
phq_readin.o \
phq_summary.o \
polarizc.o \
prepare_q.o \
print_symmetries.o \
psh_lanczos_step.o \
psyme_fpol.o \
psym_dmageq.o \
psymeq.o \
psym_dmage_fpol.o \
q_points.o \
run_nscf.o \
read_dyn_from_file.o \
read_lanczos_chain.o \
realspace_grid_init.o \
rotate_all_psi.o \
ruotaijk.o \
scal_prod.o \
set_defaults_pw.o \
set_fxc_tran.o \
set_intq_nc.o \
set_kplusq.o \
set_kup_and_kdw.o \
setup.o \
solve_e.o \
solve_e_fpolc.o \
solve_eq.o \
solve_eq_tran.o \
solve_linter.o  \
syme_fpol.o     \
symeq.o \
sym_dmage_fpol.o \
sym_dmageq.o \
transform_intq_nc.o \
transform_intq_so.o  \
write_group_info_ph.o \
zstar_eu.o          \
zstar_eu_loc.o     \
zstar_eu_us.o    

all : thermo_qe.a

thermo_qe.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

include make.depend
