# THERMO_PW QUICK HELP

<p align="justify"> In order to use <code>thermo_pw</code> you need a
working version of the <a href="http://www.quantum-espresso.org/">Quantum 
ESPRESSO (QE)</a> package. <code>Thermo_pw</code> can be downloaded from its 
<a href="http://dalcorso.github.io/thermo_pw/">main page</a> as 
a <code>.tar.gz</code> file. The current production version is 
<code>0.9.0</code> compatible with 
<code>QE-6.2.1</code>. The <code>thermo_pw</code> file should be copied
in the main (QE) directory and unpacked with the command:</p>
<p align="center">
<code>tar -xzvf thermo_pw.0.9.0.tar.gz</code>
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
<p align="justify"> To run <code>thermo_pw</code> it is useful to have the 
<a href="http://www.gnuplot.info/">gnuplot</a> package, and to plot 
the Brillouin zone
you need the <a href="http://asymptote.sourceforge.net/">asymptote</a> package. 
Both are available as precompiled packages in many distributions.
For further information, please refer to the user's guide
available in the <code>thermo_pw/Doc</code> directory.
Please report any problem to
<a href="mailto:dalcorso .at. sissa.it"> dalcorso .at. sissa.it</a>.</p>

<p align="justify"> The development version of <code>thermo_pw</code> is hosted at <a href="https://github.com/dalcorso/thermo_pw">https://github.com/dalcorso/thermo_pw</a>. To download it you need the <code>git</code>
package. Then you can give the command:</p>
<p align="center">
<code>git clone https://github.com/dalcorso/thermo_pw</code>
</p>
<p align="justify"> and you should get a directory called <code>thermo_pw</code> that contains the source code.
The <code>git</code> version can be used only together with the version of <code>QE</code> reported here: <code>6.2.1</code>. Please note that sometimes the <code>git</code> version is not working properly and in any case its use is not recommended.</p> 

<p align="justify"> Although <code>thermo_pw</code> has been 
used for several years and can be considered reasonably stable, it remains an 
experimental code given as it is.
If you are running a version of <code>QE</code> older than <code>6.2.1</code>
you can still use <code>thermo_pw</code> but you should carefully match the
versions of <code>thermo_pw</code> and of <code>QE</code> as explained in the
main <code>thermo_pw</code> page.</p>

<p align="justify"> Before using <code>thermo_pw</code>, please apply the 
patches given below.</p>

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
<a href="./q2r_sub.f90">thermo_pw/src/q2r_sub.f90</a>, substitute the one 
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
<a href="https:/people.sissa.it/~dalcorso/matdyn_sub.f90">122688</a>
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
* Problem in anharmonic properties: use <code>thermo_pw.0.3.0</code> or higher, 
or make the same changes as in commit 
<a href="http://www.qe-forge.org/gf/project/thermo_pw/scmcvs/?action=ScmCommitDetail&scm_commit_id=19508">
19508.</a>
* Modules/clocks.f90 : line 91 set <code>maxclock=200</code> otherwise 
<code>thermo_pw</code> might run out of clocks.
* Bug fix: In anharmonic calculations some vertical lines in phonon dispersion 
plots are double. Change as in commit 
<a href="http://www.qe-forge.org/gf/project/thermo_pw/scmcvs/?action=ScmCommitDetail&scm_commit_id=19386">
19386</a>. 
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
0. How can I learn to use <code>thermo_pw</code>?
<br>
Please learn the basic use of <code>Quantum ESPRESSO</code> first. 
Then you can read the <code>thermo_pw</code> user's guide and run 
the examples. These <code>FAQ</code> assume that you have a basic 
understanding of <code>thermo_pw</code> and contain miscellaneous information
not available in the user's guide.
<br><br>
1. Can I study the thermal expansion of anisotropic solids using <code>thermo_pw</code>?
<br>
For some crystal systems the answer is yes, but not all systems are supported
or tested. Read carefully the user's guide and use version <code>0.3.0</code> 
or higher. Also use dynamical matrices in .xml format or the calculation
of thermal expansion with Gruneisen parameters will not work with all 
versions previous to <code>0.5.0</code>.
<br><br>
2. Can I introduce a finite pressure?
<br>
You can calculate the equilibrium structure at a given pressure.
For other physical properties some experiments started with version
<code>0.4.0</code> but in general this part of the code is still quite 
experimental.
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
9. I cannot run the examples. I have problems using images. What 
should I do?
<br>
If you want to run the examples without images 
edit the file <code>environment_variables</code> in the main <code>QE</code>
directory. Search the two variables <code>PARA_IMAGE_PREFIX</code> and
<code>PARA_IMAGE_POSTFIX</code> and set <code>-ni 1</code>. 
<br><br>
10. I have not a parallel computer. I do not know what <code>mpi</code> is. 
Can I run <code>thermo_pw</code>?
<br>
Only <code>thermo_pw.0.5.0</code> or later versions can be compiled in serial. All previous versions must be compiled together with <code>mpi</code>.
<br><br>
11. An ionic relaxation converges with <code>pw.x</code> but not with <code>thermo_pw.x</code> (version <code>0.4.0</code>).
<br>
This is a bug of version <code>0.4.0</code>. Please apply the same changes 
as in the
commit <a href="http://www.qe-forge.org/gf/project/thermo_pw/scmcvs/?action=ScmCommitDetail&scm_commit_id=197343">197343</a> or update to a newer version.
<br><br>
12. The plot of the phonon dispersions is very strange with several disjoint
parts. Moreover the modes are not classified using symmetry. Why?
<br>
The mode symmetry analysis requires dynamical matrices in <code>.xml</code> 
format. Please put the  <code>.xml</code> extension in <code>fildyn</code>.
Symmetry matrices are needed also to recognize symmetry equivalent point
on the Brillouin zone.
<br><br>
13. The plot of the Gruneisen parameters has strange crossings in some points.
Why?
<br>
In some cases the plot of the Gruneisen parameters needs more accuracy 
on the symmetry analysis than the phonon plot. Accidentally degenerate 
frequencies might have very different Gruneisen parameters. Change the 
parameter <code>5.D-2</code> at line 148 of 
<code>PHonon/PH/find_mode_sym.f90</code> to <code>1.D-2</code> or 
less and recompile <code>thermo_pw</code>.
<br><br>
14. Thermo_pw does not compile and stops with an error saying that
<code>html.sty</code> is missing or <code>latex2html</code> is missing.
<br>
This is not a problem of <code>thermo_pw</code>. In order to compile the
documentation <code>thermo_pw</code> needs a quite complete 
<code>latex</code> package. You can find <code>html.sty</code> on
the web and copy it inside <code>thermo_pw/Doc</code> directory and you
can install <code>latex2html</code>. Even if you do not solve this problem, 
<code>thermo_pw.x</code> will be available in the <code>bin</code> directory 
of QE. Only the documentation will be missing. 
You can also remove the error editing <code>thermo_pw/Makefile</code> 
and removing the string <code>doc</code> at line 7. 
<br><br>
15. The plot of the projected band structure has some problems. Some gaps
have the same color of the projected band structure. 
<br>
This is a problem of old versions of gnuplot. Update to gnuplot 5.0 or higher.
<br><br>
16. The phonon dispersion plot seems strange, some branches are missing.
<br>
Please check that you used enough digits for the atomic positions. A typical
problem appear when you write 1/3 and 2/3 in single precision. The 
<code>pw.x</code> code finds more symmetry that it is actually present in 
the final modes and the routine that identifies the mode symmetry gets 
confused.
