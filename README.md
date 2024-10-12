**MolAR** is a **Mol**ecular **A**nalysis and modeling library for **R**ust.

# Table of contents
- [What is MolAR?](#what-is-molar)
- [Features](#features)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [Current status](#current-status)
- [Change log](CHANGELOG.md)
- [Design decisions](#design-decisions)

# What is molar?
Molar is a rust library for molecular analysis and modeling. It is started as a successor of [Pteros](https://github.com/yesint/pteros) molecular modeling library, which is written in C++ and become hard to develop and maintain due to all C++ idiosyncrasies. Eventually molar will become a "Pteros 3.0".

Molar is designed to simplify the analysis of molecular dynamics trajectories and to implement new analysis algorithms. Molar is intended to provide facilities, which are routinely used in all molecular analysis programs, namely input/output of popular file formats, powerful and flexible atom selections, geometry transformations, RMSD fitting and alignment, etc.

# Features
* Reading and writing PDB, GRO, XYZ, XTC files
    * Recognizes any VMD molfile plugins. 
    * Reading and writing Gromacs XTC format with random access.
    * Reading Gromacs TPR files in Gromacs is installed.
* Selections using the syntax similar to VMD and Pteros.
    * Memory-safe selections for serial and parallel analysis tasks.
    * Powerful subselections and selection splitting.
* SASA calculations with the fastest PowerSasa method.
* RMSD fitting and alignment.
* Basic algorithm (center of mass, center of geometry, etc.).
* Seamless PBC treatment.

# Current status
Molar is close to be feature complete and usable in useful projects. Documentation is still rudimentary.

# Installation
Molar requires Rust 1.80 or above and a C/C++ compiler for compiling third-party libraries. Any sufficiently modern gcc or clang compiler should work.

## Linking to Gromacs
In order to be able to read Gromacs TPR files MolAR should link to locally installed Gromacs. Unfortunately, modern versions of Gromacs do not expose all needed functionality in the public API, so MolAR has to hack into the internals and thus requires an access to the whole Gromacs source and build directories. This means that you have to _compile_ Gromacs on your local machine from source.

In order to link with Gromacs create a `.cargo/config.toml` file in the root directory of your project with the following content:
```toml
[env]
# Location of Gromacs source tree
GROMACS_SOURCE_DIR = "<path-to-gromacs-source>/gromacs-2023"
# Location of Gromacs *build* directory (for generated headers)
GROMACS_BINARY_DIR = "<path-to-gromacs-source>/gromacs-2023>/build"
# Location of installed gromacs libraries (where libgromacs.so is located)
GROMACS_LIB_DIR = "<path-to-installed-gromacs>/lib64"
```
You may use a template: `mv config.toml.template config.toml`.

When configuring you project set the `gromacs` feature in Cargo.toml:
```toml
[dependencies]
molar = {version="0.5", features=["gromacs"]}
```

# Tutorial
TODO

# Design decisions
Molecular analysis typically involves multiple views of the arrays of atoms and
coordinates (aka Selections). Selections could overlap arbitraryly - for example
one selection may represent the whole protein and the other only the resiudes of its
active site. Selections are mutable in the sense that their atoms could be 
manipulated - translated, rotated, aligned, renamed, etc. Such cnahges are expected 
to be immediately picked up by all other selections involving affected atoms. 

This concept doesn't play well with Rust ownership rules, where either a single
exclusive reference (`mut&`) _or_ multiple immutable references (`&`) could exists
at the same time. If one selection holds a `mut&` of the underlying array of atoms 
then no other selection is allowed to access it neither for reading nor for writing.

If we want to obey the single ownership rule we need to create and drop selections 
for every single operation to make sure that the atoms array is never aliased mutably
by several selections. This appears to be extremely inconveniet in practice and
contrasts with _all_ existing molecular analysis software.

Thus Molar uses unsafe rust and interior mutability pattern 
internally to sidestep aliasing restriction while still having safe user 
accessible API.
In Molar one can have multiple immutable references `&Sel` to selections
while still being able to invoke selection methods that mutate underlying arrays
of atoms and coordinates, such as `sel.translate()` or `sel.rotate()`. Changes
made by one selection are immediately visible to all other selections that
point to the same atoms.
