# THERMO_PW QUICK HELP

<p align="justify"> In order to use <code>thermo_pw</code> you need a
working version of the <a href="http://www.quantum-espresso.org/">Quantum 
ESPRESSO (QE)</a> package. <code>Thermo_pw</code> can be downloaded from its 
<a href="http://dalcorso.github.io/thermo_pw/">main page</a> as 
a <code>.tar.gz</code> file. The current production version is 
<code>1.7.1</code> compatible with 
<code>QE-7.1</code>. The <code>thermo_pw</code> file should be copied
in the main (QE) directory and unpacked with the command:</p>
<p align="center">
<code>tar -xzvf thermo_pw.1.7.1.tar.gz</code>
</p>
<p align="justify">This command produces a directory called 
<code>thermo_pw</code>. To compile the code you need a Fortran compiler, for 
instance the <code>gcc</code> package and <code>gfortran</code> (or 
<code>gcc-fortran</code> in some distributions),
and the same libraries required by <code>QE</code>. 
After getting the directory <code>thermo_pw</code> in the main <code>QE</code> 
directory, cd to the directory <code>thermo_pw</code> and give the command 
<code>make join_qe</code>. Then cd to the main <code>QE</code> 
directory and compile <code>thermo_pw</code> with the command:</p>
<p align="center"><code>
make thermo_pw
</code></p>
<p align="justify"> Instead, to use <code>cmake</code>, you 
create the directory <code>build</code> enter there and give the command:</p>
<p align="center"><code>
cmake -DCMAKE_C_COMPILER=c_compiler -DCMAKE_Fortran_COMPILER=fortran_compiler ../
</code></p>
<p align="justify"> Then the command 
<code>make</code> produces also the <code>thermo_pw.x</code> executable.</p>
<p align="justify"> To run <code>thermo_pw</code> it is useful to have the 
<a href="http://www.gnuplot.info/">gnuplot</a> package, and to plot 
the Brillouin zone
you need the <a href="http://asymptote.sourceforge.net/">asymptote</a> package. 
Both are available as precompiled packages in many distributions.
For further information, please refer to the user guide
available in the <code>thermo_pw/Doc</code> directory.
Please report any problem to
<a href="mailto:dalcorso .at. sissa.it"> dalcorso .at. sissa.it</a>.</p>

<p align="justify"> The development version of <code>thermo_pw</code> is hosted at <a href="https://github.com/dalcorso/thermo_pw">https://github.com/dalcorso/thermo_pw</a>. To download it you need the <code>git</code>
package. Then you can give the command:</p>
<p align="center">
<code>git clone https://github.com/dalcorso/thermo_pw</code>
</p>
<p align="justify"> and you should get a directory called <code>thermo_pw</code> that contains the source code.
The <code>git</code> version can be used only together with the version of <code>QE</code> reported here: <code>7.1</code>. Please note that sometimes the <code>git</code> version is not working properly and in any case its use is not recommended.</p> 

<p align="justify"> Although <code>thermo_pw</code> has been 
used for several years and can be considered reasonably stable, it remains an 
experimental code given as it is.
A version of <code>QE</code> older than <code>7.1</code>
can still be used with <code>thermo_pw</code> matching carefully the
versions of <code>thermo_pw</code> and of <code>QE</code> as explained in the
main <code>thermo_pw</code> page.</p>

<p align="justify"> Before using <code>thermo_pw</code>, please apply the 
patches given below.</p>

**Patches for thermo_pw.1.7.1**:
<br>
**Patches for thermo_pw.1.7.0**:
<br>
<br>
**Patches for thermo_pw.1.6.1**:
<br>
* There is a problem with the GPU version and metals. Correct as in commit 
<code>7d344d0</code> of 18/01/2022.
* Correct as in the commit of 04/03/2022 if you have problem with
magnetic systems.
* Problems with spin-orbit. Correct the routine PW/src/v_of_rho.f90 
adding the instruction v(:,:)=0.0_DP after line 208 and after line
476 and recompile.

**Patches for thermo_pw.1.6.0**:
<br>
**Patches for thermo_pw.1.5.1**:
<br>
* tools/epsilon_tpw.f90 was not updated to the QE68 conventions.
Please change as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/cd4353f48263e6015b770ef7488337f75a3184c4">commit_cd4353f</a> of 13/08/2021.
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.8.
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.
* At line 723 of upflib/write_upf_new.f90 of QE7.0 change PP_AEWFC_rel with
PP_AEWFC_REL. See also the FAQ 34.

**Patches for thermo_pw.1.5.0**:
<br>
**Patches for thermo_pw.1.4.1**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
* At line 576 of upflib/read_upf_new.f90 of QE6.7 change PP_AEWFC_rel with
PP_AEWFC_REL.
* I usually change line 400 of Modules/read_namelists.f90 of QE6.7
restoring the old default diago_david_ndim=4.
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.7. First call radial_grid_copy.
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.

**Patches for thermo_pw.1.4.0**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
* Problems with examples 10,11, and 18. Problems with electric fields and
FR-PP. Please apply the changes as in commit:
<a href="https://github.com/dalcorso/thermo_pw/commit/743148245d3ee9ce524afe3edc323d5ff3b31a92">commit_7431482</a>.
* To reproduce ph_example07 it is necessary to change line 1049 of 
file PW/src/pw_restart_new.f90 of QE6.6 as explained for version 1.3.0.
* At line 557 of upflib/read_upf_new.f90 of QE6.6 change PP_AEWFC_rel with
PP_AEWFC_REL.
* I usually change line 400 of Modules/read_namelists.f90 of QE6.6 
restoring the old default diago_david_ndim=4.
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.6. First call radial_grid_copy.
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.

**Patches for thermo_pw.1.3.1 and thermo_pw.1.3.2**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
* The phonon calculation with US-PP and PAW-PP is unstable. Correct as in
commit: 
<a href="https://github.com/dalcorso/thermo_pw/commit/51b600a25dd1d46c2ea004a509beb830d82a5811">commit_51b600a</a>.
* To reproduce ph_example07 it is necessary to change line 1049 of 
file PW/src/pw_restart_new.f90 of QE6.6 as explained for version 1.3.0.
* QE6.6 does not stop any longer if some pools have no k point. thermo_pw
is not working in this case. See in the FAQ 24 to solve this problem.
* tools/pdec.f90 does not compile with some compilers. Take the git version
of this file and recompile.
* I usually change line 400 of Modules/read_namelists.f90 of QE6.6 restoring 
the old default diago_david_ndim=4.
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.6. First call radial_grid_copy.
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.

**Patches for thermo_pw.1.3.0**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
* To reproduce ph_example07 it is necessary to change line 999 of 
file PW/src/pw_restart_new.f90 of QE6.5. Substitute angle1, angle2,
starting_magnetization with starting_magnetization, angle1, angle2.
* At line 170 of atomic/src/import_upf.f90 exchange the calls to 
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.5. First call radial_grid_copy.

**Patches for thermo_pw.1.2.1**:
<br>
* The phonon dispersions have a wrong scale when the dynamical matrices 
are written in the old format. Change as described in: 
<a href="https://github.com/dalcorso/thermo_pw/commit/3bbabc9b903a1255012f6cc0927cac705645f05b">commit 3bbabc9</a>, 
in 
<a href="https://github.com/dalcorso/thermo_pw/commit/42d4b2d4c0247a4a332f12f9d3abbaf2e77b7f38"> commit 42d4b2d</a> and in
<a href="https://github.com/dalcorso/thermo_pw/commit/35f6defc268508f8ff6dc422581ab3caa65f6e25">commit 35fedef</a> (Present also in version 1.2.0).

**Patches for thermo_pw.1.2.0**:
<br>
* When what='elastic_constants_t' a bug in QE prevents the use of 
use_free_energy=.TRUE. and elastic_algorithm='energy_std'.
Add the instruction:
ibrav_ => NULL()
at line 490 of Modules/qexsd.f90 of QE version 6.4.1.

**Patches for thermo_pw.1.1.1**:
<br>
* A bug in src/initialize_thermo_work.f90 could give wrong
geometries for mur_lc_t for odd ngeo. Change line 607 of this 
file to delta=0.0_DP. (Only in versions 1.1.1 and 1.1.0).

**Patches for thermo_pw.1.0.0**:
<br>
* Grimme-d3 not implemented. 
* zeu+US+pools not working. Apply the changes described in 
<a href="https://github.com/dalcorso/thermo_pw/commit/6c70c8f68abb017b90da9d4ce4ab0fb7620a3308">commit 6c70c8f</a>. 
* The plotted Gruneisen parameters have the wrong sign when
lmurn=.TRUE.. Apply the change described in <a href="https://github.com/dalcorso/thermo_pw/commit/d78859e8719894646ee4b416a401676c40ff8eff">commit d78859e</a>.
* Ionic relaxations are working only the first time pw.x is called. Apply the change described in <a href="https://github.com/dalcorso/thermo_pw/commit/6a1a5d8464be36d9d4ce9435f2165bd8484a6acf">commit 6a1a5d8</a>.
<br>

**Known problems of thermo_pw.0.9.0**:
<br>
* Phonons + tetrahedra are not working (not implemented yet).
* Phonons + lsda are not working (use one of previous versions).
<br>

**Patches for thermo_pw.0.9.0**:
<br>
* EELS with US-PP has still a bug. At line 75 of 
<code>qe/addusddenseq.f90</code> change <code>gg</code> with <code>qmod</code>.
* There is a problem with <code>sym_for_diago=.TRUE.</code> in the 
noncollinear/spin-orbit case. Use <code>sym_for_diago=.FALSE.</code> or use 
the following file for 
<a href="http://people.sissa.it/~dalcorso/thermo_pw/c_bands.f90">qe/c_bands.f90</a>
* The code stops using the old diagonalization in the phonon with the flag
<code>sym_for_diago=.FALSE</code>. At line 144 of 
<code>qe/set_defaults_pw.f90</code> remove the <code>_tpw</code> from the 
call to <code>set_kplusq</code>.
* Some compilers could have problems to compile the routine <code>thermo_pw/qe/set_kplusq.f90</code>. Use the following <a href="http://people.sissa.it/~dalcorso/thermo_pw/set_kplusq.f90">set_kplusq.f90</a>.
* The plotted Gruneisen parameters have the wrong sign when lmurn=.TRUE.. 
Apply the change described in <a href="https://github.com/dalcorso/thermo_pw/commit/d78859e8719894646ee4b416a401676c40ff8eff">commit d78859e</a>.  
<br>

**Patches for thermo_pw.0.8.0**:
<br>
* A bug might create some differences for phonons and US and PAW-PP calculated 
with images, so for these cases update to <code>thermo_pw.0.9.0</code> is 
recommended. 
It might affect also previous versions that use the new xml output.
<br>

**Patches for thermo_pw.0.7.0**:
<br>
* With pools, all bands are red in the band plot. At line 550 of 
<code>src/sym_band_sub.f90</code> substitute <code>nks</code> with 
<code>nkstot</code> and recompile.
* Some problems with Intel compiler can be solved as described in the patch
<a href="https://github.com/dalcorso/thermo_pw/commit/68ec9d73fb110f9a10e36b01bab02a01f80b4968">68ec9d7</a>
* When <code>emin_input</code> and <code>emax_input</code> are given in the 
<code>thermo_control</code> namelist,
some bands could be missing. The problem can be solved as described
in the patch 
<a href="https://github.com/dalcorso/thermo_pw/commit/c39f6f063433a4766dbd1999a51316f19adeebbd">c39f6f0</a>
<br>

**Patches for thermo_pw.0.6.0**:
<br>
* There is a problem with anharmonic properties calculated recovering the run
with <code>after_disp=.TRUE.</code> introduced in this version. Take the file
<a href="https://people.sissa.it/~dalcorso/q2r_sub.f90">thermo_pw/src/q2r_sub.f90</a>, substitute the one 
of <code>thermo_pw.0.6.0</code> and recompile.<br>
Moreover, at lines 11307 and 11336 of <code>lib/point_group.f90</code>,
change <code>1D-8</code> with <code>1D-5</code>.
* Modules/clocks.f90 : line 41 set <code>maxclock=200</code> otherwise 
<code>thermo_pw</code> might run out of clocks.
<br>

**Patches for thermo_pw.0.5.0**:
<br>
* Modules/clocks.f90 : line 41 set <code>maxclock=200</code> otherwise 
<code>thermo_pw</code> might run out of clocks.
<br>

**Patches for thermo_pw.0.4.0**:
<br>
* A problem with <code>max_geometries</code>: this is a bug. Add the instruction
<code>ph_geometries=0</code> at the line 431 of the file 
<code>src/thermo_pw.f90</code> and recompile.
* Compilation problem of <code>tools/test_colors.f90</code>: remove the RETURN command at the end of the file.
* Error from <code>find_aux_ind_two_groups</code> in a phonon plot. Please 
check commit
<a href="https://people.sissa.it/~dalcorso/matdyn_sub.f90">122688</a>
and make the same changes to <code>src/matdyn_sub.f90</code>.
<br>

**Patches for thermo_pw.0.3.0**:
<br>
* With some compilers the code crashes at the beginning.
Please change line 571 of <code>src/thermo_readin.f90</code> from 
<code>CALL clean_ngeo()</code> to
<code>CALL clean_ngeo(ngeo,ibrav)</code>. 
* Anharmonic properties can be calculated only with the dynamical matrix in
<code>.xml</code> format. Old format is not working. (See commit 110778).
* The code is not recovering correctly and gives an error 
<code>check_stop_init</code> not initialized. (Please apply commit 110838).
<br>

**Patches for thermo_pw.0.2.0**:
<br>
* Problem in anharmonic properties: update to a newer version.
* Modules/clocks.f90 : line 41 set <code>maxclock=200</code> otherwise 
<code>thermo_pw</code> might run out of clocks.
* Bug fix: In anharmonic calculations some vertical lines in phonon dispersion 
plots are double. Update to a newer version. 
<br>

**Patches for thermo_pw.0.1.0**:
<br>
* src/Makefile : line 83 change <code>THERMO_PW</code> with 
<code>thermo_pw</code>.
* outdir: must end with a <code>/</code>, the other case is not dealt with 
correctly.
<br> 

**FAQ:**
<br><br>
1. How can I learn to use <code>thermo_pw</code>?
<br>
Please learn the basic use of <code>Quantum ESPRESSO</code> first. 
Then you can read the <code>thermo_pw</code> tutorial and user's guide and run 
the examples. These <code>FAQ</code> assume that you have a basic 
understanding of <code>thermo_pw</code> and contain miscellaneous information
not available in the user's guide.
<br><br>
2. Can I study the thermal expansion of anisotropic solids using <code>thermo_pw</code>?
<br>
For some crystal systems, yes, but not all systems are supported
or tested. Read carefully the user's guide and use a version higher than 
<code>0.3.0</code>. Also use dynamical matrices in <code>.xml</code> format 
or the calculation
of thermal expansion with Gruneisen parameters will not work with all 
versions previous to <code>0.5.0</code>.
<br><br>
3. Can I calculate the temperature dependence of the band gap or
in general of the band structure using <code>thermo_pw</code>?
<br>
You can calculate the band structure at the crystal geometry that corresponds
to a given temperature. In this way you evaluate the effect of thermal
expansion on the band structure or on the gap. However an important 
temperature dependence of the band gaps and of the band structure
comes from the electron-phonon interactions that are not included in 
<code>thermo_pw</code>. 
For this purpose you should use another package.
<br><br>
4. Can I calculate the equilibrium geometry of a solid at a given temperature using <code>thermo_pw</code>?
<br>
Yes, but the calculation is as heavy as computing the anharmonic properties
and it will take a lot of time and resources. You need to learn how to use 
<code>thermo_pw</code> before starting such a complex calculation.
<br><br>
5. Which is the difference between <code>examples</code> and <code>inputs</code>?
<br>
<code>Examples</code> illustrate the features of <code>thermo_pw</code> and are fast, but are not converged. <code>inputs</code> are more realistic examples.
<br><br>
6. Sometimes the examples of <code>thermo_pw</code> run correctly, sometimes they crash. Which is the problem?
<br>
The most probable reason is that you have not removed the <code>results</code> 
directory produced by a previous run of the example script.
<br><br>
7. <code>make thermo_pw</code> is not working. Compilation stops with some missing routines error. 
<br>
Most probably you have not matched the versions of <code>QE</code> and of <code>thermo_pw</code>.
<br><br>
8. I have compiled <code>thermo_pw</code> but as I run it, it stops immediately.
I am using <code>thermo_pw.0.3.0</code>.
<br>
Most probably you have not applied the patch described above. Update to a
newer version.
<br><br>
9. After unpacking the tar file there is no <code>thermo_pw</code> directory.
<br>
The directory obtained unpacking the source files obtained from the github 
releases web page is called <code>thermo_pw-#version number</code>. Just 
change the name of this directory to <code>thermo_pw</code>.
<br><br>
10. I cannot run the examples. I have problems using images. What 
should I do?
<br>
If you want to run the examples without images 
edit the file <code>environment_variables</code> in the main <code>QE</code>
directory. Search the two variables <code>PARA_IMAGE_PREFIX</code> and
<code>PARA_IMAGE_POSTFIX</code> and set <code>-ni 1</code>. 
<br><br>
11. I have not a parallel computer. I do not know what <code>mpi</code> is. 
Can I run <code>thermo_pw</code>?
<br>
Only <code>thermo_pw.0.5.0</code> or later versions can be compiled in serial. All previous versions must be compiled together with <code>mpi</code>.
<br><br>
12. An ionic relaxation converges with <code>pw.x</code> but not with <code>thermo_pw.x</code> (version <code>0.4.0</code>).
<br>
This is a bug of version <code>0.4.0</code>. Please  
update to a newer version.
<br><br>
13. The plot of the phonon dispersions is very strange with several disjoint
parts. Moreover the modes are not classified using symmetry. Why?
<br>
The mode symmetry analysis requires dynamical matrices in <code>.xml</code> 
format. Please put the  <code>.xml</code> extension in the <code>fildyn</code>
variable in the <code>ph.x</code> input.
Symmetry matrices are needed also to recognize symmetry equivalent point
on the Brillouin zone.
<br><br>
14. The plot of the Gruneisen parameters has strange crossings in some points.
Why?
<br>
In some cases the plot of the Gruneisen parameters needs more accuracy 
on the symmetry analysis than the phonon plot. Accidentally degenerate 
frequencies might have very different Gruneisen parameters. Change the 
parameter <code>5.D-2</code> at line 148 of 
<code>PHonon/PH/find_mode_sym.f90</code> to <code>1.D-2</code> or 
less and recompile <code>thermo_pw</code>.
<br><br>
15. Thermo_pw documentation does not compile and stops with an error 
saying that
<code>html.sty</code> is missing or <code>latex2html</code> is missing.
<br>
This is not a problem of <code>thermo_pw</code>. In order to compile the
documentation <code>thermo_pw</code> needs a quite complete 
<code>latex</code> package. You can find <code>html.sty</code> on
the web and copy it inside <code>thermo_pw/Doc</code> directory and you
can install the package <code>latex2html</code>. Even if you do not solve this problem, 
<code>thermo_pw.x</code> will be available in the <code>bin</code> directory 
of QE. Only the documentation will be missing. 
<br><br>
16. The plot of the projected band structure has some problems. Some gaps
have the same color of the projected band structure. 
<br>
This is a problem of old versions of gnuplot. Update to gnuplot 5.0 or higher.
<br><br>
17. The phonon dispersion plot seems strange, some branches are missing.
<br>
Please check that you used enough digits for the atomic positions. A typical
problem appears when you write 1/3 and 2/3 in single precision. The 
<code>pw.x</code> code finds more symmetries than those that are
actually present in the final modes and the routine that identifies the 
mode symmetry gets 
confused.
<br><br>
18. The code fails to identify the space group and stops with an error
''point group orientation incorrect''.
<br>
Most probably you are simulating a noncollinear magnetic system. Magnetic
space group identification is not implemented but no check is done in versions
up to 0.9.0. Please make the same changes as commit a68e6cb of 18 January 2018.
If you find this error, you are using ibrav/=0, and your system is collinear,
please send me your input.
<br><br>
19. <code>what='scf_disp'</code> and partial phonon computations with
<code>start_q</code>, <code>last_q</code> or <code>start_irr</code> 
<code>last_irr</code> gives strange error messages.
<br>
The option <code>what='scf_disp'</code> requires all the dynamical matrices
files in the <code>dynamical_matrices</code> directory. Use 
<code>what='scf_ph'</code> until you collect all the dynamical matrices
and do a final run with <code>what='scf_disp'</code>.
<br><br>
20. I am computing a phonon dispersion but some <B>q</B> points are not 
computed.
<br>
Most probably you have not cleaned the <code>outdir</code> directory. Note that the 
<code>thermo_pw</code> always tries to use the content of the 
<code>outdir</code> directory if present.
<br><br>
21. Is it possible to increase the temperature range? 
<br>
Yes, you have to remove the <code>therm_files</code> directory
while keeping the <code>dynamical_matrices</code> and the 
<code>restart</code> directories. If you removed the <code>outdir</code> 
directory, use <code>after_disp=.TRUE.</code> and set <code>fildyn</code> 
with the name of the dynamical matrices.
<br><br>
22. Is it possible to increase the number of points used to compute the
phonon dos?
<br>
Yes, you have to remove both the <code>phdisp_files</code> and the
<code>therm_files</code> directories
while keeping the <code>dynamical_matrices</code> and the 
<code>restart</code> directories. 
<br><br>
23. I made a calculation with <code>with_eigen=.FALSE.</code>. Is it possible
to restart with <code>with_eigen=.TRUE.</code>?
<br>
Yes, but you have to remove both the <code>phdisp_files</code> and the 
<code>therm_files</code> directories, while keeping the 
<code>dynamical_matrices</code> and the <code>restart</code> directories.
<br><br>
24. I am using <code>thermo_pw 1.3.1</code> with<code>QE6.6</code> but 
the code hangs or stops in random places when computing phonon dispersions. 
<br>
Be careful with the use of pools. Since version <code>6.6 QE</code> 
does not stop any longer if some pools have no k points, but 
<code>thermo_pw</code> cannot deal with 
this case. In order to check if you are in this situation search the
string 'suboptimal parallelization: some nodes have no k-points' in your
output. The solution is to decrease the number of pools until the
message disappear.
If you want a permanent check of this problem use <code>thermo_pw.1.3.2</code>.
<br><br>
25. <code>tmp_dir</code> cannot be opened.
<br>
Check your <code>outdir</code> directory in the <code>pw.x</code> input.
Usually this error indicates a missing parent directory. You might have
an error in the path indicated in <code>outdir</code> or
you are not allowed to write or execute the parent directory.
<br><br>
26. Error in namelist.
<br>
Most probably there is a mistake in a variable of the namelist. Please
check accurately the user guide. The other possibility is that your
editor added some hidden characters in the input file. Please check
for it for instance with the command <code>cat -A input_file</code>.
Another possibility is that you are reading the user guide of a version
of <code>thermo_pw</code> different from the one you are using and the
variable is not yet available. Please match the versions of the user guide
and of <code>thermo_pw</code>.
<br><br>
27. Point group incompatible with the Bravais lattice.
<br>
This means that your point group is different from the point groups compatible
with a given Bravais lattice. The calculation is still possible but
<code>thermo_pw</code> will not be able to find the space group and 
will not use the symmetries to simplify the calculation of the physical 
properties. Please check if you can find why some symmetries are missing,
or why you have too many symmetries. Try to use one of the Bravais
lattices suggested by the code. The message might also indicate that
you have a supercell. If this is what you want, just ignore the message and
continue the calculation, otherwise simplify your cell.
<br><br>
28. Is <code>thermo_pw</code> compatible with the GPU version of QE?
<br>
Presently there is no support for GPU in <code>thermo\_pw</code> routines,
but with QE6.7 and with version 1.5.0 
you can give the command <code>make tpw_gpu</code>
to obtain a version of the code that can be compiled with 
<code>q-e-gpu.6.7</code>. Version 1.5.1 and QE6.8 or later versions
can be compiled with CUDA enabled. 
<br><br>
29. The band or phonon symmetry is not indicated. There are many
question marks instead of the names of the irreducible representations.
<br>
The question marks indicate that the algorithm that finds the symmetry is
confused and is unable to find a well defined symmetry. There are several 
possible reasons: 
* You might have a too small cut-off, or
a too large threshold for the self consistence and the symmetry is not
accurate enough. Please modify these parameters.
* Your atomic positions are quite close to a symmetry position but
not exactly there. This is a common problem when using single precision
atomic coordinates. Please correct the atomic coordinates adding more digits. 
* There might be some problem with the pseudopotential and
there is some ghost state. Please check other pseudopotentials to
see if the problem disappears. 
* If none of the above applies, there might be a problem in the algorithm 
that finds the symmetry. Please send me your input or post it to one of
the forum mailing lists.
<br><br>
30. Can I compute the temperature dependent elastic constants with 
<code>thermo_pw</code>?
* Quasi-static elastic constants are available from version <code>0.6.0</code>.Quasi-harmonic elastic constants require version <code>1.2.0</code> or
later. The electronic contribution is implemented only starting for
version <code>1.4.0</code>. 
* Note that this feature is still work in progress and there are still 
some limitations. For instance, atomic coordinates are relaxed only at 
zero temperature.
* Note also that the quasi-harmonic calculation is very time consuming.
As an order of magnitude it requires hundreds of phonon dispersion 
calculations.
<br><br>
31. Laue class not available when computing elastic constants.
* This error usually means that your system has less symmetry than 
expected from the Bravais lattice. For instance there is no Laue class 
available for a solid with a cubic Bravais lattice and point group symmetry 
different from T, T_d, T_h, O, or O_h.
* In this case the <code>thermo_pw</code> output writes that the point 
group and the Bravais lattice are not compatible and gives a set of 
Bravais lattices compatible with the symmetry. If you think that the 
symmetry of your system is correct, then you should use one of the 
Bravais lattices suggested by <code>thermo_pw</code>. Instead if some 
symmetry is missing for other reasons (see the Quantum Espresso PW 
user's guide for possible reasons), 
then correct the problem before running the elastic constant calculation.
* If you are using supercells or you have a low dimensional system in a 
supercell, then probably thermo_pw is not yet suited to compute automatically
the elastic constants of your system.
<br><br>
32. <code>Thermo_pw</code> does not compile (with Quantum ESPRESSO version 7.0 or later). There is an error no rule to make file <code>make.depend</code>.
* After writing <code>make join_qe</code> and returning to QE root directory
you need to rerun <code>./configure</code> before <code>make thermo_pw</code>.
(Thanks to H. Zhao for reporting the problem).
33. <code>Thermo_pw</code> has problems with fully relativistic PAW. 
* Before reporting such problems check the correct matching of the 
PP_AEWFC_REL tag in the UPF file and in the upf reading routine. 
Starting from version 6.5 of QE and until version 6.7 the xml tag for 
the small component of the all electron partial waves has been called 
PP_AEWFC_rel, while in previous versions of QE it was called PP_AEWFC_REL. 
As such the fully relativistic
pseudopotentials created with a version of QE older than 6.5 (as many of 
the PPs distributed in the QE site) were no more read correctly. The code
does not stop and most of the time produces only slightly uncorrect results 
expecially in the PP test. To read a PP that contains the 
PP_AEWFC_REL tag you can change the file upflib/read_upf_new.f90 and 
upflib/write_upf_new.f90 searching
for the string PP_AEWFC_rel in both files and changing it into PP_AEWFC_REL.
In QE6.8, QE7.0, and QE7.1 the routine upflib/read_upf_new.f90 
has been corrected and these versions of QE read correctly UPF PPs with the tag
PP_AEWFC_REL, but unfortunately continue to write UPF PP with the tag 
PP_AEWFC_rel. The routine upflib/write_upf_new.f90 must be changed to 
make the code consistent and to read correctly the pseudopotential generated
by the same version of QE. For version QE7.0 apply also the correction to
PW/src/v_of_rho.f90 described above. 
