import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import time

t = time.process_time()
    
# Load the reference structure (PDB) and trajectory (XTC)
u = mda.Universe('protein.pdb', 'protein.xtc')  # Replace with your file paths

# Select the residue of interest (residue 1 in this case)
residue_1 = u.select_atoms('resid 560')  # Selecting residue 1

cm = []

# Iterate through each frame of the trajectory
nearby_atoms = u.select_atoms("around 10.0 resid 560", updating=True)
for frame_index, ts in enumerate(u.trajectory):
    # Compute the center of mass for the selected atoms
    cm.append( nearby_atoms.center_of_mass() )

print(cm)
print(f"Elapsed: {time.process_time() - t}")