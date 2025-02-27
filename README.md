**MolAR** is a **Mol**ecular **A**nalysis and modeling library for **R**ust.

# Table of contents
- [What is MolAR?](#what-is-molar)
- [Features](#features)
- [Current status](#current-status)
- [Design and performance](#design-and-performance)
- [Installation](#installation)
- [Tutorial](#tutorial)

# What is molar?
MolAR is a library for molecular modeling and analysis written in Rust with an emphasis on memory safety and performance. 

Molar is designed to simplify the analysis of molecular dynamics trajectories and to implement new analysis algorithms. Molar is intended to provide facilities, which are routinely used in all molecular analysis programs, namely input/output of popular file formats, powerful and flexible atom selections, geometry transformations, RMSD fitting and alignment, etc.

MolAR is a logical successor of [Pteros](https://github.com/yesint/pteros) molecular modeling library, which is written in C++ and become hard to develop and maintain due to all C++ idiosyncrasies.

# Features
* Reading and writing PDB, GRO, XYZ, XTC, TPR files
    * Recognizes any VMD molfile plugins. 
    * Reading and writing Gromacs XTC format with random access.
    * Reading Gromacs TPR files if Gromacs is installed.
* Selections using the syntax similar to VMD and Pteros.
    * Memory-safe selections for serial and parallel analysis tasks.
    * Powerful subselections and selection splitting.
* SASA calculations with the fastest PowerSasa method.
* RMSD fitting and alignment.
* Basic algorithm (center of mass, center of geometry, etc.).
* Seamless PBC treatment.

# Design and Performance
Please refer to the [MolAR paper](https://onlinelibrary.wiley.com/doi/10.1002/jcc.27536).

# Current status
Molar is close to be feature complete and usable in useful projects. Documentation is still rudimentary.

# Installation
Molar requires Rust 1.80 or above and a C/C++ compiler for compiling third-party libraries. Any sufficiently modern gcc or clang compiler should work.

To add MolAR to your Rust project just use `cargo add molar`.

## Linking to Gromacs
In order to be able to read Gromacs TPR files MolAR should link to locally installed Gromacs. Unfortunately, modern versions of Gromacs do not expose all needed functionality in the public API, so MolAR has to hack into the internals and thus requires an access to the whole Gromacs source and build directories. This means that you have to _compile_ Gromacs on your local machine from source.

In order to link with Gromacs create a `.cargo/config.toml` file in the root directory of your project with the following content:
```toml
[env]
# Location of Gromacs source tree
GROMACS_SOURCE_DIR = "<path-to-gromacs-source>/gromacs-2023"
# Location of Gromacs *build* directory (for generated headers)
GROMACS_BUILD_DIR = "<path-to-gromacs-source>/gromacs-2023>/build"
# Location of installed gromacs libraries (where libgromacs.so is located)
GROMACS_LIB_DIR = "<path-to-installed-gromacs>/lib64"
```
You may use a template: `mv config.toml.template config.toml`.

When configuring you project enable the `gromacs` feature for MolAR in Cargo.toml:
```toml
[dependencies]
molar = {version="*", features=["gromacs"]}
```

# Tutorial
We will write an example program that reads a file of some molecular system containing TIP3P water molecules, convert all water to TIP4P and saves this as a new file. TIP3P water has 3 particles (oxygen and two hydrogens), while TIP4P [has 4]((http://www.sklogwiki.org/SklogWiki/index.php/TIP4P/2005_model_of_water)) (oxygen, two hydrogens and a dummy particle). Our goal is to add these dummy particles to each water molecule.

## Preparing the stage
First, let's create a new Rust project called `tip3to4`: 
```shell
cargo new tip3to4
```

Then, let's add dependencies: the MolAR itself and anyhow crate for easy error reporting:
```shell
cargo add molar anyhow
```

In `src/main.rs` add needed boilerplate:
```rust,ignore
// For processing command line arguments
use std::env;
// Import all baic things from molar
use molar::prelude::*;
// For error handling
use anyhow::Result;

fn main() -> Result<()> {
    // Get the command line arguments
    let args: Vec<String> = env::args().collect();

    // Here our program is going to be written

    // Report successful completion of the program
    Ok(())
}
```

Now we can start writing our program.

## Reading an input file
The simples way of loading the molecular system in MolAR is to use a `Source` - an object that holds a `Topology` and `State` of the system and is used to create atom selections for manipulating this data:

```rust,ignore
// Load the source file from the first command line argument
let src = Source::serial_from_file(&args[0])?;
```

Unlike other molecular analysis libraries, there are four kinds of sources and atom selections in MolAR: `serial`, `serial builder`, `parallel mutable` and `parallel immutable`. This is required to enforce memory safety and to guarantee the absense of data races in paralell programs. For now we will just work with the simplest `serial` kind of sources and selections, which behave in the most intuitive way similar to what you see in other analysis libraries. `Source::serial_from_file()` creates such serial `Source` by reading a file specified in the first command line argument.

The file type (PDB, GRO, etc) is automatically recognized by its extention.

If reading the file fails for whatever reason the `?` operator will return an error, which will be nicely printed by `anyhow` crate.

## Making selections
Now we need to select all waters that are going to be converted to TIP4. We also need to select all non-water part of the system to keep it as is.

```rust,ignore
let water = src.select_str("resname TIP3")?;
let non_water = src.select_str("not resname TIP3")?;
```

Selections are created with the syntax that is very similar to one used in VMD Pteros and Gromacs. Here we select water and non-water by residue name.

In MolAR empty selections are not permitted, so if no atoms are selected (or if anything else goes wrong) the error will be reported.

## Looping over indiviudual water molecules
We selected all water molecules as a single selection but we need to loop over individual water molecules to add an additional dummy particle to each of them. In order to do this we are splitting a selection to fragments by the residue index:

```rust,ignore
// Go over water molecules one by one                   
for mol in water.into_iter_contig_resindex() {
    // Do something with mol
}
```

The method `split_resindex_into_iter()` returns a Rust iterator, which produces contigous selections containig distinct residue index each. There are many other ways of splitting selections into parts using arbitrary logic in MolAR, but this simplest one is what we need now. 

## Working with coordinates
Now we need to get the coordinates of atoms for current water molecules and compute a position of the dummy atom.

```rust,ignore
// Go over water molecules one by one                   
for mol in water.split_resindex_into_iter() {
    // TIP3 is arranged as O->H->H
    // so atom 0 is O, atoms 1 and 2 are H
    // Get cooridnates
    let o_pos = mol.nth_pos(0).unwrap();
    let h1_pos = mol.nth_pos(1).unwrap();
    let h2_pos = mol.nth_pos(2).unwrap();
    // Get center of masses of H
    let hc = 0.5*(h1_pos.coords + h2_pos.coords);
    // Unit vector from o to hc
    let v = (hc-o_pos.coords).normalize();	
    // Position of the M dummy particle in TIP4
    let m_pos = o_pos + v*0.01546;
    // Dummy atom M
    let m_at = Atom {   
        resname: "TIP4".into(),
        name: "M".into(),
        ..mol.first_particle().atom.clone()
    };
    println!("{:?} {:?}",m_at,m_pos);
}
```

First, we are getting the coordinates of oxigen and two hydrogens. `nth_pos(n)` returns the position of n-th atom in selection. Since n may potentially be out of range, it returns an `Option<&Pos>`. We are sure that there are just 3 atoms in water molecule, so we just unwrapping an option.

Then we are computing the position of the dummy atom, which is on the bissection of H-O-H angle at the distance of 0.01546 from the oxygen.

Finally, we are constructing a new `Atom` with name 'M', residue name 'TIP4' and all other fields (resid,resindex, etc) the same as in our water molecule.

We are printing our new dummy atom and its position just to be sure that everything works as intended.

## Constructing output system
All this is fine, but we still have no system to write our converted water molecules to. Let's fix this and modify the beginning of our main funtion like this:

```rust,ignore
// Load the source file from the first command line argument
let src = Source::serial_from_file(&args[0])?;

// Make empty output system
let out = Source::empty_builder();
```

Here we are creating new empty `Source` of kind `builder`. This means that we will be able to add and delete the atoms to this source. Conventional `serial` source can access and alter existing atoms, but can't add or delete them. Such a distinction is dictated by performance and memory safety reasons - `builder` sources and selections require additional range checks, which make them a bit slower, so it only makes sense to use them when you actually need to add or delete the atoms.

The first thing that we add to out new empty system is all non-water atoms:
```rust,ignore
// Add non-water atoms to the output
out.append(&non_water);
```

Now, at the end of our loop over water molecules, we can add new dummy atoms properly to the new system:
```rust,ignore
// Add new converted water molecule
// We assume that the dummy is the last atom.
out.append_atoms(
    mol.iter_atoms().cloned().chain(std::iter::once(m_at)),
    mol.iter_pos().cloned().chain(std::iter::once(m_pos)),
);
```

This code snippet may look a bit puzzling for non-rustaceans, so let's go through it.
- `append_atoms()` method accepts two iterators: the first yielding atoms and the second yielding their corresponding coordinates. 
- Our selected water molecule `mol` has methods `iter_atoms()` and `iter_pos()` for getting these iterators. 
- `cloned()` adaptor is used to get copies of existing atoms and coordinates instead of references to them. 
- We add our new dummy atom at the end of water molecule by "chaining" another iterator at the end of the current one. `std::iter::once(value)` returns an iterator yielding a single value and allows us to add newly constructed `m_at` and `m_pos` to the corrsponding iterators.

## Writing the output file
Out output system is now fully constructed but it still lacks an important element - the periodic box description. Most molecular systems originating from MD are periodic and the information about the periodic box has to be copied to our newly constructed system:

```rust,ignore
// Transfer the box from original file
out.set_box_from(&src);
```

Here we provide a reference to the input system, so the box is cloned from it to the output system.

Finally we are ready to write the output file:

```rust,ignore
// Write out new system
out.save(&args[1])?;
```

Again, file format will be determined by extension. The file name is provided by the second command line argument.

## The final result
The complete program looks like this:
```rust,no_run
// For processing command line arguments
use std::env;
// Import all baic things from molar
use molar::prelude::*;
// For error handling
use anyhow::Result;

fn main() -> Result<()> {
    // Get the command line arguments
    let args: Vec<String> = env::args().collect();

    // Load the source file from the first command line argument
    let src = Source::serial_from_file(&args[0])?;

    // Make empty output system
    let out = Source::empty_builder();

    let water = src.select_str("resname TIP3")?;
    let non_water = src.select_str("not resname TIP3")?;

    // Add non-water atoms to the output
    out.append(&non_water);

    // Go over water molecules one by one                   
    for mol in water.split_resindex_into_iter() {
        // TIP3 is arranged as O->H->H
        // so atom 0 is O, atoms 1 and 2 are H
	    // Get cooridnates
        let o_pos = mol.nth_pos(0).unwrap();
        let h1_pos = mol.nth_pos(1).unwrap();
        let h2_pos = mol.nth_pos(2).unwrap();
	    // Get center of masses of H
	    let hc = 0.5*(h1_pos.coords + h2_pos.coords);
	    // Unit vector from o to hc
	    let v = (hc-o_pos.coords).normalize();	
	    // Position of the M dummy particle in TIP4
	    let m_pos = o_pos + v*0.01546;
        // Dummy atom M
        let m_at = Atom {   
            resname: "TIP4".into(),
            name: "M".into(),
            ..mol.first_particle().atom.clone()
        };

        // Add new converted water molecule
        // We assume that the dummy is the last atom.
        out.append_atoms(
            mol.iter_atoms().cloned().chain(std::iter::once(m_at)),
            mol.iter_pos().cloned().chain(std::iter::once(m_pos)),
        );

    }

    // Transfer the box
    out.set_box_from(&src);

    // Write out new system
    out.save(&args[1])?;

    // Report successful completion of the program
    Ok(())
}
```