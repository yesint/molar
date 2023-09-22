# MolAR

**MolAR** is a **Mol**ecular **A**nalysis and modeling library for **R**ust.

### What is molar?

Molar is a rust library for molecular analysis and modeling. It is started as a successor of [Pteros](https://github.com/yesint/pteros) molecular modeling library, which is written in C++ and become hard to develop and maintain due to all C++ idiosyncrasies. Eventually molar may become a "Pteros 3.0".

Molar is designed to simplify the analysis of molecular dynamics trajectories and to implement new analysis algorithms. Molar is intended to provide facilities, which are routinely used in all molecular analysis programs, namely input/output of popular file formats, powerful and flexible atom selections, geometry transformations, RMSD fitting and alignment, etc.

### Current status

Molar is currently in early alpha stage. 

#### Currently implemented features
* Reading and writing PDB and XTC files.
* Reading Gromacs TPR files.
* Selections using the syntaxys similar to VMD and Pteros.
* Automatic seamless compiling and linking with VMD molfile and xdrfile libraries and the git version of Gromacs.

