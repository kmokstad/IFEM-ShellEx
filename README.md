# IFEM-ShellEx

A shell FE application using external element routines.
This application is based on unstructured Lagrange interpolation
(the [ASMu2DLag](https://github.com/OPM/IFEM/blob/master/src/ASM/ASMu2DLag.h) class),
and can read FE models from Nastran bulk data files.

The Fortran source code of ANDES shell element used here are from the
[fedem-solvers](https://github.com/openfedem/fedem-solvers) project,
whereas the code in the FFlLib sub-folder are from the related
[fedem-foundation](https://github.com/openfedem/fedem-foundation) project.
