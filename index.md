<p align="justify"> Thermo_pw is a set of Fortran drivers for the parallel and/or automatic 
computation of materials properties using Quantum ESPRESSO (QE) routines 
as the underlying engine. It provides an alternative organization of the 
QE work-flow for the most common tasks exploiting, when possible, 
an asynchronous image parallelization. Moreover, the code has a set of 
pre-processing tools to reduce the input information given by the user 
and a set of 
post-processing tools to produce plots directly comparable with experiment.</p>
<p align="justify"> A quick introduction to the <code>thermo_pw</code> code can be found 
<a href="https://dalcorso.github.io/thermo_pw/thermo_pw_help.html">here</a>,
a brief tutorial is available <a href="https://people.sissa.it/dalcorso/thermo_pw/tutorial/tutorial.html">here</a>,
while the user's guide of <code>thermo_pw</code> version <code>1.6.1</code> 
can be found <a href="https://people.sissa.it/dalcorso/thermo_pw/user_guide/index.html">here</a>.</p>
<p align="justify"> Presently there is no reference work for citing <code>thermo_pw</code>. If you want to mention it in your work, you can put a reference to this web page.</p>
<p align="justify">The following papers describe new
features implemented in <code>thermo_pw</code>:</p>

12. C. Malica and A. Dal Corso,
Quasi-harmonic thermoelasticity of palladium, platinum, copper and gold
from first principles,
<a href="https://iopscience.iop.org/article/10.1088/1361-648X/ac2041">
J. of Phys.: Condens. Matter <B>33</B>, 475901 (2021)</a>

11. O. Motornyi, N. Vast, I. Timrov, O. Baseggio, S. Baroni, and A. Dal Corso,
Electron energy loss spectroscopy of bulk gold with ultrasoft pseudopotentials and the Liouville-Lanczos method,
<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.035156">
Phys. Rev. B <B>102</B>, 035156 (2020).</a>

10. C. Malica and A. Dal Corso,
Quasi-harmonic temperature dependent elastic constants: Applications 
to Silicon, Aluminum, and Silver,
<a href="https://iopscience.iop.org/article/10.1088/1361-648X/ab8426/meta">
J. of Phys.: Condens. Matter <B>32</B>, 315902 (2020).</a>
<br>
<br>
9. A. Urru and A. Dal Corso,
Density functional perturbation theory for lattice dynamics with fully
relativistic ultrasoft pseudopotentials: The magnetic case,
<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.045115">Phys. Rev. B <B>100</B>, 045115 (2019).</a>
<br>
<br>
8. C. Malica and A. Dal Corso,
Temperature dependent atomic B-factor: an ab-initio calculation,
<a href="https://doi.org/10.1107/S205327331900514X">Acta Cryst. A <B>75</B>, 624 (2019).</a>
<br>
<br>
7. A. Urru and A. Dal Corso,
Spin-polarized electronic surface states of Re(0001): an ab-initio investigation,
<a href="https://doi.org/10.1016/j.susc.2019.03.008">Surf. Sci. <B>686</B>, 22 (2019).</a>
<br>
<br>
6. O. Motornyi, M. Raynaud, A. Dal Corso, and N. Vast,
Simulation of electron energy loss spectra with the turboEELS and
thermo&#95;pw codes,
<a href="https://iopscience.iop.org/article/10.1088/1742-6596/1136/1/012008/meta">J. Phys.: Conf. Ser. <B>1136</B>, 012008 (2018).</a>
<br>
<br>
5. A. Urru and A. Dal Corso,
Clean Os(0001) electronic surface states: a first-principle fully relativistic investigation,
<a href="https://www.sciencedirect.com/science/article/pii/S0039602817309469">Surf. Sci. <B> 671</B>, 17 (2018).</a>
<br>
<br>
4. M. Palumbo and A. Dal Corso,
Lattice dynamics and thermophysical properties of h.c.p. Os and Ru from
the quasi-harmonic approximation,
<a href="http://iopscience.iop.org/article/10.1088/1361-648X/aa7dca">
J. of Phys.: Condens. Matter <B>29</B>, 395401 (2017).
</a>
<br>
<br>
3. M. Palumbo and A. Dal Corso,
Lattice dynamics and thermophysical properties of h.c.p. Re and Tc from
the quasi-harmonic approximation,
<a href="http://dx.doi.org/10.1002/pssb.201700101">Physica Status Solidi B:
Basic Solid State Physics <B>254</B>, 1700101 (2017).
</a>
<br>
<br>
2. A. Dal Corso,
Elastic constants of Beryllium: a first-principles investigation,
<a href="http://dx.doi.org/10.1088/0953-8984/28/7/075401"> J. Phys.: Condens. Matter <B>28</B>, 075401 (2016) </a>.
<br>
<br>
1. A. Dal Corso,
Clean Ir(111) and Pt(111) electronic surface states: a first-principle fully relativistic investigation,
<a href="http://www.sciencedirect.com/science/article/pii/S0039602815000734"> Surf. Sci. <B>637-638</B>, 106 (2015)</a>.
<br>
<br>
<p>Refs.1,3 describe how <code>thermo_pw</code> computes temperature
dependent elastic constants, Ref.11 elastic constants at T=0 K,
Ref.12,8 surface band structures, Refs.2,7 eels and optical properties,
Ref.4 DFPT with fully relativistic pseudopotentials, and Ref.9,10 
thermodynamic anharmonic properties.</p>

<p align="justify">The following works contain some calculations made by <code>thermo_pw</code>:</p>
<br>
5. C. Malica,
From ab-initio thermodynamics to quasi-harmonic thermoelastic properties 
of crystals: A new workflow and selected applications,
<a href="https://iris.sissa.it/retrieve/handle/20.500.11767/125489/148394/Cristiano_Malica_Thesis.pdf">PhD thesis.</a>
<br>
<br>
4. A. Urru,
Lattice dynamics with Fully Relativistic Pseudopotentials for magnetic systems,
with selected applications,
<a href="https://iris.sissa.it/handle/20.500.11767/115671#.X8uStKYo-cM">PhD thesis.</a>
<br>
<br>
3. A. Urru and A. Dal Corso,
Lattice dynamics effects on the magnetocrystalline anisotropy
energy: application to MnBi,
<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.115126">
Phys. Rev. B <B>102</B>, 115126 (2020).</a>
<br>
<br>
2. C. Malica and A. Dal Corso,
Temperature dependent elastic constants and thermodynamic properties 
of BAs: an ab-initio investigation, 
<a href="https://aip.scitation.org/doi/10.1063/5.0011111">
Jour. of Applied Phys. <B>127</B>, 245103 (2020).</a>
<br>
<br>
1. S. Poncé, D. Jena, and F. Giustino,
Hole mobility of strained GaN from first principles,
<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.085204">
Phys. Rev. B <B>100</B>, 085204 (2019)</a>.

See also the presentation given at the Quantum-ESPRESSO developers meeting 2017:
<br>
<br>
<a href="https://people.sissa.it/~dalcorso/thermo_pw_2017.pdf">Thermo_pw_2017.pdf</a>
<br>
<br>
and at the Quantum-ESPRESSO developers meeting 2018:
<br>
<br>
<a href="https://people.sissa.it/~dalcorso/thermo_pw_2018.pdf">Thermo_pw_2018.pdf</a>

<p align="justify">
The latest developments of the <code>thermo_pw</code> software can be
followed <a href="https://github.com/dalcorso/thermo_pw/commits/master">here</a>.</p>

<p align="justify">For problems to compile or run <code>thermo_pw</code> or if you think
that you have found a bug, please check the quick-help page mentioned above, 
apply all the patches and if your problem is not solved, post it to the
<a href="mailto:thermo_pw-forum@lists.quantum-espresso.org">thermo_pw-forum mailing list</a> or e-mail me: <a href="mailto:dalcorso .at. sissa.it">dalcorso .at. sissa.it</a>. To subscribe to the <code>thermo_pw-forum</code> mailing list
click <a href="https://lists.quantum-espresso.org/mailman/listinfo/thermo_pw-forum">here</a>.
<br>
Please do not send me questions about the input of <code>pw.x</code>.
<br>
Errors such as: 
<br>
1. tmp_dir cannot be opened
<br>
2. error in namelist 
<br>
are exactly what the error message says. Some possible solutions
are indicated in the FAQ (25, 26) 
<a href="https://dalcorso.github.io/thermo_pw/thermo_pw_help.html">here.</a>
If these do not solve your problem I have no better answer. 
You can also find help in the QE user guide, at the
<code>users@lists.quantum-espresso.org</code> mailing list or search
in the examples directories.</p>

**Thermo_pw downloads**:
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.6.1.tar.gz">
thermo_pw.1.6.1.tar.gz</a>  (released 10-01-2022) compatible with QE-7.0.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.6.0.tar.gz">
thermo_pw.1.6.0.tar.gz</a>  (released 27-12-2021) compatible with QE-6.8.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.5.1.tar.gz">
thermo_pw.1.5.1.tar.gz</a>  (released 22-07-2021) compatible with QE-6.7.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.5.0.tar.gz">
thermo_pw.1.5.0.tar.gz</a>  (released 19-07-2021) compatible with QE-6.7.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.4.1.tar.gz">
thermo_pw.1.4.1.tar.gz</a>  (released 29-12-2020) compatible with QE-6.7.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.4.0.tar.gz">
thermo_pw.1.4.0.tar.gz</a>  (released 22-12-2020) compatible with QE-6.6.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.3.2.tar.gz">
thermo_pw.1.3.2.tar.gz</a>  (released 17-8-2020) compatible with QE-6.6.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.3.1.tar.gz">
thermo_pw.1.3.1.tar.gz</a>  (released 13-8-2020) compatible with QE-6.6.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.3.0.tar.gz">
thermo_pw.1.3.0.tar.gz</a>  (released 12-8-2020) compatible with QE-6.5.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.2.1.tar.gz">
thermo_pw.1.2.1.tar.gz</a>  (released 23-1-2020) compatible with QE-6.5.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.2.0.tar.gz">
thermo_pw.1.2.0.tar.gz</a>  (released 28-12-2019) compatible with QE-6.4.1.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.1.1.tar.gz">
thermo_pw.1.1.1.tar.gz</a>  (released 16-04-2019) compatible with QE-6.4.1.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.1.0.tar.gz">
thermo_pw.1.1.0.tar.gz</a>  (released 16-04-2019) compatible with QE-6.4.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.0.9.tar.gz">
thermo_pw.1.0.9.tar.gz</a>  (released 06-03-2019) compatible with QE-6.3.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.1.0.0.tar.gz">
thermo_pw.1.0.0.tar.gz</a>  (released 17-07-2018) compatible with QE-6.3.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw-0.9.9.tar.gz">
thermo_pw-0.9.9.tar.gz</a>  (released 05-07-2018) compatible with QE-6.2.1.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw-0.9.0.tar.gz">
thermo_pw-0.9.0.tar.gz</a>  (released 20-12-2017) compatible with QE-6.2.1.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw-0.8.0.tgz">
thermo_pw-0.8.0.tar.gz</a>  (released 24-10-2017) compatible with QE-6.2.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw-0.8.0-beta.tgz">
thermo_pw-0.8.0-beta.tar.gz</a>  (released 31-08-2017) compatible with QE-6.2-beta.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.7.9.tgz">
thermo_pw.0.7.9.tgz</a>  (released 06-07-2017) compatible with QE-6.1.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.7.0.tgz">
thermo_pw.0.7.0.tgz</a>  (released 18-03-2017) compatible with QE-6.1.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.6.0.tgz">
thermo_pw.0.6.0.tgz</a>  (released 05-10-2016) compatible with QE-6.0.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.5.0.tar.gz">
thermo_pw.0.5.0.tar.gz</a>  (released 26-04-2016) compatible with QE-5.4.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.4.0.tar.gz">
thermo_pw.0.4.0.tar.gz</a>  (released 23-01-2016) compatible with QE-5.3.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.3.0.tar.gz">
thermo_pw.0.3.0.tar.gz</a>  (released 23-06-2015) compatible with QE-5.2.0 and QE-5.2.1.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.2.0.tar.gz">
thermo_pw.0.2.0.tar.gz</a>   (released 13-03-2015) compatible with QE-5.1.2.
<br>
<br>
- <a href="http://people.sissa.it/%7Edalcorso/thermo_pw/thermo_pw.0.1.0.tar.gz">
thermo_pw.0.1.0.tar.gz</a>   (released 28-11-2014) compatible with QE-5.1.1.
<br>
<br>
<p align="justify">Please note that the versions of <code>thermo_pw</code> and 
of <code>QE</code> must be carefully matched as written above. Mixing two 
unmatched versions leads to a compilation error.</p>
<br>
