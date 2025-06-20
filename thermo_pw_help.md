# THERMO_PW QUICK HELP

<p align="justify"> 
To use <code>thermo_pw</code>, a functional version of the <a href="http://www.quantum-espresso.org/">Quantum ESPRESSO (QE)</a> package is required. <code>Thermo_pw</code> can be downloaded from its <a href="http://dalcorso.github.io/thermo_pw/">main page</a> as a <code>.tar.gz</code> file. The current production version is <code>2.0.2</code>, compatible with <code>QE-7.4</code>.</p>

<p align="justify"> 
Copy the <code>thermo_pw</code> archive into the main <code>QE</code> directory and unpack it using the command:</p>
<p align="center"> 
<code>tar -xzvf thermo_pw.2.0.2.tar.gz</code>
</p>
<p align="justify"> 
This command will create a directory named <code>thermo_pw</code>.
To compile the code, you will need a Fortran compiler (e.g., the <code>gcc</code> package with <code>gfortran</code>, or <code>gcc-fortran</code> in some distributions) and the same libraries required by <code>QE</code>.</p>

<p align="justify">
Once the <code>thermo_pw</code> directory is present in the main 
<code>QE</code> directory, cd into <code>thermo_pw</code> and write:</p>
<p align="center">
<code>make join_qe</code>
</p>
<p align="justify">
Then, cd back to the main <code>QE</code> directory and compile <code>thermo_pw</code> with the commands:</p>
<p align="center">
<code>./configure</code> </p>
<p align="center">
<code>make thermo_pw</code></p>

<p align="justify">
Alternatively, to use <code>cmake</code>, create a <code>build</code> directory, enter it, and execute:</p>
<p align="center">
<code>cmake -DCMAKE_C_COMPILER=c_compiler -DCMAKE_Fortran_COMPILER=fortran_compiler ../ </code></p>
<p align="justify">
or, to compile for GPU: </p>
<p align="center">
<code>cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DQE_ENABLE_CUDA=ON -DQE_ENABLE_OPENACC=ON -DQE_ENABLE_MPI_GPU_AWARE=ON ../</code></p>

<p align="justify">
Subsequently, the <code>make</code> command will also generate the <code>thermo_pw.x</code> executable.</p>

<p align="justify">
For running <code>thermo_pw</code>, the <a href="http://www.gnuplot.info/">gnuplot</a> package is beneficial. To plot the Brillouin zone, the <a href="http://asymptote.sourceforge.net/">asymptote</a> package is necessary. Both are typically available as precompiled packages in many distributions.</p>

<p align="justify">
For further details, please consult the user guide located in the <code>thermo_pw/Doc</code> directory. Please report any issues to <a href="mailto:dalcorso .at. sissa.it">dalcorso.at.sissa.it</a>.</p>

**Development Version (Git):**
<p align="justify"> 
The development version of <code>thermo_pw</code> is hosted on <a href="https://github.com/dalcorso/thermo_pw">https://github.com/dalcorso/thermo_pw</a>. To download it, the git package is required. Then, you can run the command:</p>
<p align="center"> 
<code>git clone https://github.com/dalcorso/thermo_pw</code></p>

<p align="justify"> 
This will create a thermo_pw directory containing the source code.
Important Note: The Git version can only be used with <code>QE</code> version <code>7.4</code>. Be aware that the Git version may not always function correctly, and its use is generally not recommended.
Although thermo_pw has been in use for several years and is considered reasonably stable, it remains an experimental code provided "as is."</p>

**Compatibility with Older QE Versions:**

<p align="justify"> 
Older versions of <code>QE</code> (earlier than <code>7.4</code>) can still be used with <code>thermo_pw</code> by carefully matching the <code>thermo_pw</code> and <code>QE</code> versions, as detailed on the main <code>thermo_pw</code> page.</p>

<p align="justify"> Before using <code>thermo_pw</code>, please apply the 
patches given below.</p>

**Patches for thermo_pw.2.0.2**:
<br>

**Patches for thermo_pw.2.0.1**:
<br>
* <code>many_k=.TRUE.</code> not working with LSDA. Substitute the file
<code>thermo_pw/qe/incdrhoscf_dev.f90</code>
with this <a href="https://people.sissa.it/~dalcorso/incdrhoscf_dev.f90">file</a>.
<br>

**Patches for thermo_pw.2.0.0**:
<br>
* To compile with <code>cmake</code> copy in <code>thermo_pw/CMakeLists.txt</code> 
the file that you find <a href="https://people.sissa.it/~dalcorso/thermo_pw/CMakeLists.txt">here</a>.
<br>
* <code>many_k=.TRUE.</code> not working with LSDA. Correct as described in version <code>2.0.1</code>.
<br>

**Patches for thermo_pw.1.9.1**:
<br>
* The code hangs when using <code>start_q</code> and <code>last_q</code>
with <code>what='elastic_constants_geo'</code>. Correct
as in commit <code>48b77cc</code> of Mar. 11, 2024.
<br>
* <code>many_k=.TRUE.</code> not working with LSDA. Correct as described in version <code>2.0.1</code>.
<br>

**Patches for thermo_pw.1.9.0**:
<br>
* At line 150 of <code>qe/many_k_ph.f90</code> continuation line
'&' is missing.
<br>
* Thermo_pw was not working with scalapak.
See bug fix of Feb. 1, 2024.
<br>
* A problem with calculation of EELS spectrum with 
Sternheimer method.
See bug fix of Jan. 31, 2024
<br>

**Patches for thermo_pw.1.8.1**:
<br>
* Thermo_pw was not working with scalapak.
See bug fix of Feb. 1, 2024
<br>
* A problem with calculation of EELS spectrum with 
Sternheimer method.
See bug fix of Jan. 31, 2024
<br>
* A problem with calculation of electronic free energy 
See bug fix aeedce4 of Jul. 3, 2023.
<br>

**Patches for thermo_pw.1.8.0**:
<br>
* To compile using cmake copy the thermo_pw/CMakeLists.txt of thermo_pw.1.8.1
in the thermo_pw directory.
<br>
* A problem with calculation of electronic free energy 
See bug fix aeedce4 of Jul. 3, 2023.
<br>

**Patches for thermo_pw.1.7.1**:
<br>

**Patches for thermo_pw.1.7.0**:
<br>

**Patches for thermo_pw.1.6.1**:
<br>
* There is a problem with the GPU version and metals. Correct as in commit 
<code>7d344d0</code> of 18/01/2022.
<br>
* Correct as in the commit of 04/03/2022 if you have problem with
magnetic systems.
<br>
* Problems with spin-orbit. Correct the routine PW/src/v_of_rho.f90 
adding the instruction v(:,:)=0.0_DP after line 208 and after line
476 and recompile.
<br>

**Patches for thermo_pw.1.6.0**:
<br>

**Patches for thermo_pw.1.5.1**:
<br>
* tools/epsilon_tpw.f90 was not updated to the QE68 conventions.
Please change as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/cd4353f48263e6015b770ef7488337f75a3184c4">commit_cd4353f</a> of 13/08/2021.
<br>
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.8.
<br>
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.
<br>
* At line 723 of upflib/write_upf_new.f90 of QE7.0 change PP_AEWFC_rel with
PP_AEWFC_REL. See also the FAQ 34.
<br>

**Patches for thermo_pw.1.5.0**:
<br>

**Patches for thermo_pw.1.4.1**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
<br>
* At line 576 of upflib/read_upf_new.f90 of QE6.7 change PP_AEWFC_rel with
PP_AEWFC_REL.
<br>
* I usually change line 400 of Modules/read_namelists.f90 of QE6.7
restoring the old default diago_david_ndim=4.
<br>
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.7. First call radial_grid_copy.
<br>
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.
<br>

**Patches for thermo_pw.1.4.0**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
<br>
* Problems with examples 10,11, and 18. Problems with electric fields and
FR-PP. Please apply the changes as in commit:
<a href="https://github.com/dalcorso/thermo_pw/commit/743148245d3ee9ce524afe3edc323d5ff3b31a92">commit_7431482</a>.
<br>
* To reproduce ph_example07 it is necessary to change line 1049 of 
file PW/src/pw_restart_new.f90 of QE6.6 as explained for version 1.3.0.
<br>
* At line 557 of upflib/read_upf_new.f90 of QE6.6 change PP_AEWFC_rel with
PP_AEWFC_REL.
<br>
* I usually change line 400 of Modules/read_namelists.f90 of QE6.6 
restoring the old default diago_david_ndim=4.
<br>
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.6. First call radial_grid_copy.
<br>
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.
<br>

**Patches for thermo_pw.1.3.1 and thermo_pw.1.3.2**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
<br>
* The phonon calculation with US-PP and PAW-PP is unstable. Correct as in
commit: 
<a href="https://github.com/dalcorso/thermo_pw/commit/51b600a25dd1d46c2ea004a509beb830d82a5811">commit_51b600a</a>.
<br>
* To reproduce ph_example07 it is necessary to change line 1049 of 
file PW/src/pw_restart_new.f90 of QE6.6 as explained for version 1.3.0.
<br>
* QE6.6 does not stop any longer if some pools have no k point. thermo_pw
is not working in this case. See in the FAQ 24 to solve this problem.
<br>
* tools/pdec.f90 does not compile with some compilers. Take the git version
of this file and recompile.
<br>
* I usually change line 400 of Modules/read_namelists.f90 of QE6.6 restoring 
the old default diago_david_ndim=4.
<br>
* At line 170 of atomic/src/import_upf.f90 exchange the calls to
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.6. First call radial_grid_copy.
<br>
* At line 131 of thermo_pw/qe/pheqscf.f90 remove tpiba2 to have example21
working again.

**Patches for thermo_pw.1.3.0**:
<br>
* Still a missing transformation of ftau into ft. Please change
as in commit <a href="https://github.com/dalcorso/thermo_pw/commit/3e3953152e2c81d301eb6a596de97eba07f6841d">commit_3e39531</a> of 7/5/2021.
<br>
* To reproduce ph_example07 it is necessary to change line 999 of 
file PW/src/pw_restart_new.f90 of QE6.5. Substitute angle1, angle2,
starting_magnetization with starting_magnetization, angle1, angle2.
<br>
* At line 170 of atomic/src/import_upf.f90 exchange the calls to 
set_pawsetup and radial_grid_copy to have the atomic
paw tests working again in QE6.5. First call radial_grid_copy.
<br>

**Patches for thermo_pw.1.2.1**:
<br>
* The phonon dispersions have a wrong scale when the dynamical matrices 
are written in the old format. Change as described in: 
<a href="https://github.com/dalcorso/thermo_pw/commit/3bbabc9b903a1255012f6cc0927cac705645f05b">commit 3bbabc9</a>, 
in 
<a href="https://github.com/dalcorso/thermo_pw/commit/42d4b2d4c0247a4a332f12f9d3abbaf2e77b7f38"> commit 42d4b2d</a> and in
<a href="https://github.com/dalcorso/thermo_pw/commit/35f6defc268508f8ff6dc422581ab3caa65f6e25">commit 35fedef</a> (Present also in version 1.2.0).
<br>

**Patches for thermo_pw.1.2.0**:
<br>
* When what='elastic_constants_t' a bug in QE prevents the use of 
use_free_energy=.TRUE. and elastic_algorithm='energy_std'.
Add the instruction:
ibrav_ => NULL()
at line 490 of Modules/qexsd.f90 of QE version 6.4.1.
<br>

**Patches for thermo_pw.1.1.1**:
<br>
* A bug in src/initialize_thermo_work.f90 could give wrong
geometries for mur_lc_t for odd ngeo. Change line 607 of this 
file to delta=0.0_DP. (Only in versions 1.1.1 and 1.1.0).
<br>

**Patches for thermo_pw.1.0.0**:
<br>
* Grimme-d3 not implemented. 
<br>
* zeu+US+pools not working. Apply the changes described in 
<a href="https://github.com/dalcorso/thermo_pw/commit/6c70c8f68abb017b90da9d4ce4ab0fb7620a3308">commit 6c70c8f</a>. 
<br>
* The plotted Grüneisen parameters have the wrong sign when
lmurn=.TRUE.. Apply the change described in <a href="https://github.com/dalcorso/thermo_pw/commit/d78859e8719894646ee4b416a401676c40ff8eff">commit d78859e</a>.
<br>
* Ionic relaxations are working only the first time pw.x is called. Apply the change described in <a href="https://github.com/dalcorso/thermo_pw/commit/6a1a5d8464be36d9d4ce9435f2165bd8484a6acf">commit 6a1a5d8</a>.
<br>

**Known problems of thermo_pw.0.9.0**:
<br>
* Phonons + tetrahedra are not working (not implemented yet).
<br>
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
<br>
* Some compilers could have problems to compile the routine <code>thermo_pw/qe/set_kplusq.f90</code>. Use the following <a href="http://people.sissa.it/~dalcorso/thermo_pw/set_kplusq.f90">set_kplusq.f90</a>.
<br>
* The plotted Grüneisen parameters have the wrong sign when lmurn=.TRUE.. 
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
<br>
* Some problems with Intel compiler can be solved as described in the patch
<a href="https://github.com/dalcorso/thermo_pw/commit/68ec9d73fb110f9a10e36b01bab02a01f80b4968">68ec9d7</a>
<br>
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
<br>
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
<br>
* Compilation problem of <code>tools/test_colors.f90</code>: remove the RETURN command at the end of the file.
<br>
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
<br>
* Anharmonic properties can be calculated only with the dynamical matrix in
<code>.xml</code> format. Old format is not working. (See commit 110778).
<br>
* The code is not recovering correctly and gives an error 
<code>check_stop_init</code> not initialized. (Please apply commit 110838).
<br>

**Patches for thermo_pw.0.2.0**:
<br>
* Problem in anharmonic properties: update to a newer version.
<br>
* Modules/clocks.f90 : line 41 set <code>maxclock=200</code> otherwise 
<code>thermo_pw</code> might run out of clocks.
<br>
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
<br>
1. How can I learn to use <code>thermo_pw</code>?
<br>
To begin, please familiarize yourself with the basic use of <code>Quantum ESPRESSO</code>. Afterward, you can read the <code>thermo_pw</code> tutorial and user's guide and run the examples. These FAQs assume a basic understanding of <code>thermo_pw</code> and contain miscellaneous information not available in the user's guide.
<br><br>
2. Can I study the thermal expansion of anisotropic solids using <code>thermo_pw</code>?
<br>
Yes, for certain crystal systems, but not all are currently supported or extensively tested. Please read the user's guide carefully and ensure you are using a version higher than <code>0.3.0</code>. Additionally, use dynamical matrices in <code>.xml</code> format, as calculations of thermal expansion with Grüneisen parameters will not work with all versions prior to <code>0.5.0</code> otherwise.
<br><br>
3. Can I calculate the temperature dependence of the band gap or
in general of the band structure using <code>thermo_pw</code>?
<br>
You can calculate the band structure at the crystal geometry corresponding to a specific temperature. This allows you to evaluate the effect of thermal expansion on the band structure or band gap. However, a significant component of the temperature dependence of band gaps and band structures arises from electron-phonon interactions, which are not included in <code>thermo_pw</code>. For such calculations, you will need to use a different package.
<br><br>
4. Can I calculate the equilibrium geometry of a solid at a given temperature using <code>thermo_pw</code>?
<br>
Yes, this is possible. However, the calculation is as computationally intensive as determining quasi-harmonic properties and will require significant time and resources. We strongly recommend that users gain a thorough understanding of <code>thermo_pw</code>'s usage before attempting such a complex calculation.
<br><br>
5. Which is the difference between <code>examples</code> and <code>inputs</code>?
<br>
<code>Examples</code> are designed to quickly illustrate <code>thermo_pw</code>'s features; however, they are not fully converged for accuracy. In contrast, <code>inputs</code> provide more realistic, production-ready examples.
<br><br>
6. Why do <code>thermo_pw</code> examples sometimes run correctly and sometimes crash?
<br>
The most probable reason is that you have not removed the results directory generated by a previous run of the example script.
<br><br>
7. <code>make thermo_pw</code> is not working. The compilation stops with a missing routines error. What's wrong?
<br>
The most probable reason is an incompatibility between the versions of <code>QE</code> and <code>thermo_pw</code>. Please ensure they are compatible.
<br><br>
8. After unpacking the tar file, I don't see a <code>thermo_pw</code> directory. Where is it?
<br>
When you unpack the source files downloaded from the GitHub releases page, the directory is named <code>thermo_pw-#version_number</code> (e.g., <code>thermo_pw-1.0.0</code>). Simply rename this directory to <code>thermo_pw</code> to proceed.
<br><br>
9. I'm having trouble running the examples, specifically with issues 
related to images. What should I do?
<br>
If you wish to run the examples without using images, you need to modify the <code>environment_variables</code> file located in your main <code>QE</code> directory. Locate the variable <code>PARA_IMAGE_POSTFIX</code> and set <code>-ni 1</code>.
<br><br>
10. I don't have a parallel computer and am unfamiliar with <code>MPI</code>. Can I still run <code>thermo_pw</code>? 
<br>
Yes, you can, but only <code>thermo_pw</code> version 0.5.0 and later can be compiled and run in serial (single-processor) mode. Earlier versions require <code>MPI</code> for compilation and execution.
<br><br>
11. The phonon dispersion plot appears strange, with disjoint parts and modes not classified by symmetry. What is the cause?
<br>
This issue arises because the mode symmetry analysis requires dynamical matrices in <code>.xml</code> format. To resolve this, please ensure you include the <code>.xml</code> extension in the <code>fildyn</code> variable within your <code>ph_control</code> input file. Symmetry matrices are also crucial for recognizing symmetry-equivalent points in the Brillouin zone.
<br><br>
12. Why does the plot of Grüneisen parameters show strange crossings at certain points?
<br>
In some cases, plotting Grüneisen parameters requires higher accuracy in symmetry analysis than the phonon plot. Accidentally degenerate frequencies might have significantly different Grüneisen parameters. To address this, open <code>PHonon/PH/find_mode_sym.f90</code>, change the parameter <code>5.D-2</code> at line 148 to <code>1.D-2</code> (or a smaller value), and then recompile <code>thermo_pw</code>.
<br><br>
13. The <code>thermo_pw</code> documentation fails to compile with errors indicating <code>html.sty</code> or <code>latex2html</code> are missing. How can I fix this?
<br>
This is not an issue with <code>thermo_pw</code> itself. Compiling the <code>thermo_pw</code> documentation requires a comprehensive <code>LaTeX</code> distribution. To resolve the errors, you can:
<br>
 * Download <code>html.sty</code> from the web and copy it into the <code>thermo_pw/Doc</code> directory.
<br>
 * Install the <code>latex2html</code> package.
<br>
Even if you do not resolve these documentation compilation issues, the <code>thermo_pw.x</code> executable will still be available in the <code>bin</code> directory of your <code>QE</code> installation. Only the documentation will be inaccessible.
<br><br>
14.  I'm experiencing issues with the projected band structure plot; some gaps have the same color as the projected band structure. What's the problem?
<br>
This is a known issue with older versions of <code>Gnuplot</code>. Please update to <code>Gnuplot 5.0</code> or higher to resolve it.
<br><br>
15.  My phonon dispersion plot looks strange, with some branches missing. What could be the issue?
<br>
This often occurs if you haven't used enough digits for the atomic positions in your input. A common problem arises when positions like 1/3 and 2/3 are specified with insufficient precision. In such cases, the <code>pw.x</code> code may incorrectly identify more symmetries than are actually present in the final modes, confusing the routine that identifies mode symmetry.
<br><br>
16. The code fails to identify the space group and stops with an error: "point group orientation incorrect." Why?
<br>
This error most likely indicates you are simulating a noncollinear magnetic system. Keep in mind that magnetic space group identification isn't implemented in <code>thermo_pw</code>. Versions up to <code>0.9.0</code> don't perform a check for this, which can lead to the error.
To resolve this, please apply the changes introduced in commit <code>a68e6cb</code> (January 18, 2018).
If you encounter this error while using <code>ibrav/=0</code> with a collinear system, please send your input file for further investigation.
<br><br>
17. When using <code>what='scf_disp'</code> with partial phonon computations (e.g., <code>start_q</code>, <code>last_q</code>, <code>start_irr</code>, <code>last_irr</code>), I receive strange error messages. Why?
<br>
The <code>what='scf_disp'</code> option requires all the dynamical matrix files to be present in the <code>dynamical_matrices</code> directory. If you are performing partial phonon calculations, use <code>what='scf_ph'</code> until all necessary dynamical matrices have been generated and collected. Once all files are present, you can then run with <code>what='scf_disp'</code> for the final computation.
<br><br>
18. I'm computing a phonon dispersion, but some <B>q</B> points are not being computed. What's the reason? 
<br>
The most probable reason is that you have not cleaned your <code>outdir</code> directory. Note that <code>thermo_pw</code> always attempts to reuse content from the <code>outdir</code> directory if it finds any.
<br><br>
19. Is it possible to increase the temperature range of my calculation?
<br>
Yes, you can. To do so, you must remove the <code>therm_files</code> directory while preserving the <code>dynamical_matrices</code> and <code>restart</code> directories. If you have also removed the <code>outdir</code> directory, set <code>after_disp=.TRUE.</code> in your input and specify the name of the dynamical matrices using the <code>fildyn</code> variable.
<br><br>
20.  Is it possible to increase the number of points used to compute the phonon DOS?
<br>
Yes, you can. To do so, you must remove both the <code>phdisp_files</code> and <code>therm_files</code> directories while preserving the <code>dynamical_matrices</code> and <code>restart</code> directories.
<br><br>
21. I performed calculations with <code>with_eigen=.FALSE.</code>. Is it 
possible to restart with <code>with_eigen=.TRUE.</code>?
<br>
Yes, you can. To do so, you must remove both the <code>phdisp_files</code> and <code>therm_files</code> directories while preserving the <code>dynamical_matrices</code> and <code>restart</code> directories.
<br><br>
22. I'm using <code>thermo_pw 1.3.1</code> with <code>QE 6.6</code>, and the code hangs or stops randomly during phonon dispersion calculations. What's the issue? 
<br>
Be cautious when using parallel pools. Since version <code>6.6</code>, <code>QE</code> no longer halts if some pools lack <B>k</B>-points, but <code>thermo_pw</code> (version <code>1.3.1</code>) cannot handle this scenario, leading to hangs or crashes.
To check if this is your situation, search your output file for the string: <code>'suboptimal parallelization: some nodes have no k-points'</code>.
Solution: Decrease the number of pools until this message no longer appears.
For a permanent check and resolution of this problem, consider upgrading to <code>thermo_pw 1.3.2</code> or higher.
<br><br>
23.  I'm receiving an error message that <code>tmp_dir</code> cannot be opened. What should I do?
<br>
This error typically indicates a problem with the <code>outdir</code> directory specified in your <code>pw.x</code> input file. Most probable causes include:
<br>
  * A missing parent directory in the path specified for <code>outdir</code>.
<br>
  * Insufficient permissions to write to or execute in the <code>outdir</code>'s parent directory.
<br>
Please check the <code>outdir</code> path and your directory permissions.
<br><br>
24.  I'm receiving an "Error in namelist." What should I do?
<br>
This error most likely indicates a mistake in a variable within your namelist. Please consult the user guide carefully to verify all variable spellings and expected values.
<br>
Other common causes include:
<br>
* Hidden characters: Your text editor might have introduced hidden characters into the input file. You can check for these using a command like <code>cat -A input_file</code>.
<br>
* Version mismatch: You might be using a user guide from a different <code>thermo_pw</code> version than the one you are running, leading to unrecognized variables. Always ensure that the version of your user guide matches your <code>thermo_pw</code> installation.
<br><br>
25. I'm receiving the error "Point group incompatible with the Bravais lattice." What does this mean?
<br>
This message indicates that the point group identified for your system is not compatible with the specified Bravais lattice. While the calculation can still proceed, <code>thermo_pw</code> will not be able to automatically determine the space group, nor will it utilize symmetries to simplify the calculation of physical properties.
<br>
Possible causes and solutions:
<br>
Incorrect symmetry:
  * Missing symmetries: Investigate why some symmetries might be absent from your input (e.g., imprecise atomic coordinates).
<br>
  * Excess symmetries: Conversely, you might have defined too many symmetries that are not truly present in the actual structure.
Action: Try using one of the Bravais lattices suggested by the code, or review your structure for subtle deviations from the intended symmetry.
<br>
Supercell:
<br>
  * The message might also appear if you are intentionally simulating a supercell. Action: If this is your intent, you can safely ignore this message and continue the calculation. Otherwise, simplify your unit cell to match the primitive cell.
<br><br>
26. Is <code>thermo_pw</code> compatible with the GPU version of QE?
<br>
Partial compatibility exists.
    With <code>QE 6.7</code> and <code>thermo_pw</code> version <code>1.5.0</code>, you can compile a GPU-compatible version of <code>thermo_pw</code> by running <code>make tpw_gpu</code>. This is designed to work with <code>q-e-gpu.6.7</code>.
For <code>thermo_pw</code> version <code>1.5.1</code> and <code>QE 6.8</code> (or later versions), <code>thermo_pw</code> can be compiled with CUDA support using the same commands you would use to enable CUDA in standard QE. This typically involves configuring QE with specific CUDA options.
<code>Thermo_pw</code> versions from <code>1.9.0</code> to <code>2.0.2</code> have been specifically tested on the Leonardo supercomputer at CINECA using the NVIDIA Fortran compiler.
These later versions also include custom GPU routines, exclusive to <code>thermo_pw</code>, that enable the simultaneous calculation of numerous <B>k</B>-points on the GPU. These routines are optimized for metallic systems with small cells.
<br><br>
27. The band or phonon symmetry is not indicated; instead, many question marks appear in place of irreducible representation names. Why?
<br>
The presence of question marks signifies that the symmetry-finding algorithm is unable to definitively determine the symmetry. This can be due to several reasons:
<br>
* Insufficient accuracy: Your calculation might have a cut-off that is too small, or a self-consistency threshold that is too large, resulting in insufficient accuracy for proper symmetry identification. Action: Please adjust these parameters.
<br>
* Imprecise atomic positions: Your atomic positions may be very close to, but not exactly on, a true symmetry position. This is a common issue when using single-precision atomic coordinates. Action: Correct your atomic coordinates by adding more digits of precision.
<br>
* Pseudopotential issues: There might be a problem with the pseudopotential, possibly indicating a "ghost state". Action: Try using different pseudopotentials to see if the issue resolves.
<br>
* Algorithm problem: If none of the above solutions apply, there might be an issue with the symmetry-finding algorithm itself. Action: In this case, please send me your input file or post it to one of the forum mailing lists.
<br><br>
28. Can I compute the temperature dependent elastic constants with 
<code>thermo_pw</code>?
<br>
Yes, with version-specific capabilities:
* Quasi-static elastic constants are available starting from version <code>0.6.0</code>.
<br>
* Quasi-harmonic elastic constants require version <code>1.2.0</code> or later.
<br>
* The electronic contribution to elastic constants is implemented only from version <code>1.4.0</code> onwards.
<br>
Important notes:
<br>
* This feature is still under development and has some limitations. For instance, atomic coordinates are currently relaxed only at zero temperature (within ZSISA). Full free energy minimization is possible only when there is only one internal degree of freedom.
<br>
* Be aware that quasi-harmonic calculations are highly time-consuming, often requiring, as an order of magnitude, hundreds of phonon dispersion computations.
<br><br>
29. I'm receiving a "Laue class not available" error when computing elastic constants. What does this mean?
<br>
* This error typically indicates that your system possesses fewer symmetries than expected for its specified Bravais lattice. For example, for a solid with a cubic Bravais lattice, a Laue class cannot be assigned if its point group symmetry is different from T, T_d, T_h, O, or O_h.
<br>
* In such cases, <code>thermo_pw</code>'s output will state that the point group and Bravais lattice are incompatible and will suggest alternative compatible Bravais lattices.
<br>
Recommended actions:
<br>
* If you are confident in your system's symmetry: Consider using one of the Bravais lattices suggested by <code>thermo_pw</code>.
<br>
* If symmetries are missing due to other reasons: (e.g., imprecise atomic coordinates, as discussed in the <code>QE</code> <code>PW</code> user's guide), you must correct these issues before proceeding with the elastic constant calculation.
<br>
* For supercells or low-dimensional systems: If you are using supercells, or have a low-dimensional system within a supercell, <code>thermo_pw</code> might not yet be automatically suited for computing the elastic constants of your system.
<br><br>
30. <code>thermo_pw</code> fails to compile with <code>QE</code> version <code>7.0</code> or later, showing a "no rule to make file make.depend" error. How can I fix this?
<br>
After running <code>make join_qe</code> and returning to the <code>QE</code> root directory, you must re-run <code>./configure</code> before attempting <code>make thermo_pw</code>. This step is crucial for regenerating the necessary build dependencies.
(Thanks to H. Zhao for reporting the problem).
<br><br>
31. <code>Thermo_pw</code> is having problems with fully relativistic PAW pseudopotentials. What should I do?
<br>
Before reporting any issues, please check for a mismatch in the <code>PP_AEWFC_REL</code> tag between your <code>UPF</code> file and <code>thermo_pw</code>'s <code>UPF</code> reading routines. This is a common source of problems.
<br>
Understanding the tag mismatch:
<br>
* From <code>QE</code> versions <code>6.5</code> to <code>6.7</code>, the XML tag for the small component of all-electron partial waves was named <code>PP_AEWFC_rel</code>. In previous <code>QE</code> versions, it was <code>PP_AEWFC_REL</code>.
This change means fully relativistic pseudopotentials created with <code>QE</code> versions older than <code>6.5</code> (including many distributed on the <code>QE</code> site) might no longer be read correctly by <code>thermo_pw</code> versions relying on the <code>PP_AEWFC_rel</code> tag. The code often doesn't stop but might produce subtly incorrect results, especially during pseudopotential tests.
<br>
Solutions for tag mismatch:
<br>
* For PPs with <code>PP_AEWFC_REL</code> (older <code>QE</code>): If your pseudopotential contains the <code>PP_AEWFC_REL</code> tag, manually edit <code>upflib/read_upf_new.f90</code> and <code>upflib/write_upf_new.f90</code>. In both files, search for the string <code>PP_AEWFC_rel</code> and change it to <code>PP_AEWFC_REL</code>.
<br>
Consistency with newer <code>QE</code> versions:
</br>
* <code>QE</code> 6.8, 7.0, and 7.1: The <code>upflib/read_upf_new.f90</code> routine in these QE versions correctly reads <code>UPF</code> PPs with the <code>PP_AEWFC_REL</code> tag. However, they continue to write <code>UPF</code> PPs with the <code>PP_AEWFC_rel</code> tag, creating an inconsistency. To ensure <code>thermo_pw</code> correctly reads pseudopotentials generated by the same QE version, you must modify <code>upflib/write_upf_new.f90</code> (changing <code>PP_AEWFC_rel</code> to <code>PP_AEWFC_REL</code>).
For <code>QE</code> 7.0 only: Additionally, apply the correction described above to <code>PW/src/v_of_rho.f90</code> for <code>QE</code> 7.0.
