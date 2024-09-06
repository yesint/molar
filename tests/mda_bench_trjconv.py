import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import time

t = time.process_time()

# Load the reference structure (PDB) and trajectory (XTC)
u = mda.Universe('protein.pdb', 'protein.xtc')  # Replace with your file paths

selection = u.select_atoms("resid 560")

# Create an XTCWriter to write the new trajectory
with XTCWriter('.extracted.xtc', n_atoms=selection.n_atoms) as writer:
    # Iterate over each frame of the original trajectory
    for ts in u.trajectory:
        # Write only the selected atoms to the new trajectory
        writer.write(selection)

print(f"Elapsed: {time.process_time() - t}")