![thermo-logo](Doc/thermo_pw.jpg)

> This is the distribution of the thermo\_pw package 
> (see https://dalcorso.github.io/thermo_pw). This program is
> a driver to compute the physical properties of materials, using
> Quantum ESPRESSO (QE) (see http://www.quantum-espresso.org)
> as the underlying engine.

Thermo\_pw uses the same input format as Quantum ESPRESSO's <code>pw.x</code> code and generates PostScript figures of various material properties. While less flexible than QE for certain aspects, it offers a simpler and faster learning curve for properties such as electronic band structures, phonon dispersions, and both harmonic and anharmonic thermodynamic quantities. Additionally, it supports parallel execution, allowing multiple instances to run asynchronously for image generation. For an overview, refer to [Doc/tutorial.pdf](https://people.sissa.it/~dalcorso/thermo_pw/tutorial.pdf) (or [here](https://people.sissa.it/~dalcorso/thermo_pw/tutorial/tutorial.html)), and for a detailed description of its options, consult [Doc/user\_guide.pdf](https://people.sissa.it/~dalcorso/thermo_pw/user_guide.pdf) (or [here](https://people.sissa.it/~dalcorso/thermo_pw/user_guide/user_guide.html)).

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

NB: This code replaces the main <code>Makefile</code>, <code>CMakeLists.txt</code>, the <code>install/makedeps.sh</code> and <code>install/plugins_makefile</code> files within your Quantum ESPRESSO (QE) distribution. The original QE files are only restored when you type <code>make leave_qe</code>. If you simply delete the <code>thermo_pw</code> directory, the QE package files will not be restored, which could cause issues if you try to reinstall <code>thermo_pw</code>. Additionally, all files in the <code>QESUB</code> directory are replaced. Please refer to the <code>AAAREADME</code> file in that directory for details on the modifications made to QE.

## PEOPLE

The THERMO\_PW code is primarily designed, written, and maintained by Andrea
Dal Corso
(SISSA - Trieste).

Some routines have been contributed by SISSA PhD students and post-docs.
Among them I mention M. Palumbo, O. Motornyi, A. Urru, C. Malica, X. Gong,
B. Thakur, and A. Ahmed.

I would like also to thank all the people that contributed with comments,
requests for improvements, and bug reports.

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
