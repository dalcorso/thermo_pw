![thermo-logo](Doc/thermo_pw.jpg)

> This is the distribution of the thermo\_pw package 
> (see https://dalcorso.github.io/thermo_pw). This program is
> a driver to compute the physical properties of materials, using
> Quantum ESPRESSO (QE) (see http://www.quantum-espresso.org)
> as the underlying engine.

Thermo\_pw reads the same input as the pw.x code of QE and produces postscript
figures of some material properties. While less flexible than QE,
for properties such as the electronic band structures, the phonon
dispersions, or the harmonic and anharmonic thermodynamic quantities,
it is simpler to use and faster to learn. Moreover it can run in parallel
creating several images of itself, carrying out asynchronous tasks.
See [Doc/tutorial.pdf](https://people.sissa.it/~dalcorso/thermo_pw/tutorial.pdf) (or [here](https://people.sissa.it/~dalcorso/thermo_pw/tutorial/tutorial.html)) for an overview of the code
and the file [Doc/user\_guide.pdf](https://people.sissa.it/~dalcorso/thermo_pw/user_guide.pdf) (or [here](https://people.sissa.it/~dalcorso/thermo_pw/user_guide/user_guide.html)) for a detailed description of its 
options.

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## USAGE

Put this directory inside the main Quantum ESPRESSO (QE) directory,
cd here and type

```
make join_qe
```

cd to the main QE directory and type

```
./configure
make thermo_pw
```

Now you should have the file thermo\_pw.x in the QE/bin directory and
you can run the examples.

Please check the [quick-help page](https://dalcorso.github.io/thermo_pw/thermo_pw_help.html) before running, for possible patches to this version.

Uninstal:
cd here and type
```
make leave_qe
```
Then remove this directory.

NB: This code substitutes the main Makefile, and the files install/makedeps.sh
and install/plugins\_makefile of the QE distribution. Only after typing
make leave\_qe the original files are copied in the QE distribution. If
you just remove the thermo\_pw directory the files of the QE package are not
restored and you could have problems to reinstall thermo\_pw.

## LICENSE

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
