import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import numpy as np
import time

t = time.process_time()

# Load the reference structure (PDB file)
ref = mda.Universe('protein.pdb')  # Replace with your PDB file path

# Load the trajectory file (XTC) using the same topology
u = mda.Universe('protein.pdb', 'protein.xtc')  # Replace with your files

# Select the protein atoms (common selection can be changed if needed)
ref_atoms = ref.select_atoms('protein')
u_atoms = u.select_atoms('protein')

# Align each frame of the trajectory to the reference structure
# The alignment is done using the RMSD alignment procedure from the align module

# Prepare RMSD calculation
rmsd = []

# Iterate through each frame of the trajectory
for ts in u.trajectory:
    # Align the mobile atoms to the reference
    align.alignto(u_atoms, ref_atoms)
    
    # Calculate RMSD for the aligned structure against the reference
    rmsd_value = rms.rmsd(u_atoms.positions, ref_atoms.positions)
    rmsd.append(rmsd_value)

print(rmsd)
print(f"Elapsed: {time.process_time() - t}")