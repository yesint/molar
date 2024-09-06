import mdtraj as md
import numpy as np
import time

t = time.process_time()

# Load the reference structure (PDB file)
reference = md.load('protein.pdb')  # Replace with your PDB file path

# Load the trajectory file (XTC) along with the reference structure topology
trajectory = md.load('protein.xtc', top='protein.pdb')  # Replace with your XTC and PDB file paths

# Select the protein atoms for alignment and RMSD calculation (e.g., CA atoms)
atom_indices = reference.topology.select('all')

# Align the trajectory to the reference structure using the selected atoms
trajectory.superpose(reference, atom_indices=atom_indices)

# Calculate RMSD of each frame to the reference structure
rmsd = md.rmsd(trajectory, reference, atom_indices=atom_indices)

print(rmsd)
print(f"Elapsed: {time.process_time() - t}")