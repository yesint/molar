# MolAR

**MolAR** is a **Mol**ecular **A**nalysis and modeling library for **R**ust.

## What is molar?

Molar is a rust library for molecular analysis and modeling. It is started as a successor of [Pteros](https://github.com/yesint/pteros) molecular modeling library, which is written in C++ and become hard to develop and maintain due to all C++ idiosyncrasies. Eventually molar may become a "Pteros 3.0".

Molar is designed to simplify the analysis of molecular dynamics trajectories and to implement new analysis algorithms. Molar is intended to provide facilities, which are routinely used in all molecular analysis programs, namely input/output of popular file formats, powerful and flexible atom selections, geometry transformations, RMSD fitting and alignment, etc.

## Features
* Reading and writing PDB, GRO, XTC files.
* Reading Gromacs TPR files.
* Selections using the syntaxis similar to VMD and Pteros.
* Subselections and splitting selections.
* SASA calculations with the fastest PowerSasa method.
* RMSD fitting and alignment
* Basic algorithm (center of mass, center of geometry, etc.)
* PBC unwrapping
* Automatic seamless compiling and linking with VMD molfile and xdrfile libraries and the git version of Gromacs.

## Current status
Molar is currently close to be usable in useful projects. Documentation is still missing.

## Dependencies
Molar depends on C/C++ compiler and CMake for compiling third-party libraries.

## Change log
### v0.3.0
* Selection kinds are introduced allowing for chosing between overlapping single-threaded mutable, overlapping multi-threaded immmutable and non-overlapping multi-threaded mutable kinds.
* Distance seacrh is parallelized automatically.
