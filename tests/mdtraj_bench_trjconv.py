import mdtraj as md
import time

t = time.process_time()

# Load the reference structure and trajectory
trajectory = md.load('protein.xtc', top='protein.pdb')  # Replace with your file paths

atom_indices = trajectory.topology.select("residue 560")

# Slice the trajectory to include only the selected atoms
selected_trajectory = trajectory.atom_slice(atom_indices)

# Save the selected atoms' trajectory to a new XTC file
selected_trajectory.save_xtc('.extracted.xtc')

print(f"Elapsed: {time.process_time() - t}")