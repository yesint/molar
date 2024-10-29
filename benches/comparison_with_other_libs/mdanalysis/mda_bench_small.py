import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis import align, rms
import time
import numpy as np



def trjconv_bench():
    u = mda.Universe('tests/protein.pdb', 'tests/protein.xtc')
    selection = u.select_atoms("resid 560")
    with XTCWriter('tests/.extracted.xtc', n_atoms=selection.n_atoms) as writer:
        for ts in u.trajectory:
            writer.write(selection)


def align_bench():
    ref = mda.Universe('tests/protein.pdb')
    u = mda.Universe('tests/protein.pdb', 'tests/protein.xtc')
    ref_atoms = ref.select_atoms('protein')
    u_atoms = u.select_atoms('protein')
    # Prepare RMSD calculation
    rmsd = []
    # Iterate through each frame of the trajectory
    for ts in u.trajectory:
        # Align the mobile atoms to the reference
        align.alignto(u_atoms, ref_atoms)
        # Calculate RMSD for the aligned structure against the reference
        rmsd_value = rms.rmsd(u_atoms.positions, ref_atoms.positions)
        rmsd.append(rmsd_value)
    #print(rmsd)


def within_bench():
    u = mda.Universe('tests/protein.pdb', 'tests/protein.xtc')  # Replace with your file paths
    cm = []
    # Iterate through each frame of the trajectory
    nearby_atoms = u.select_atoms("around 10.0 resid 560", updating=True)
    for _,_ in enumerate(u.trajectory):
        cm.append( nearby_atoms.center_of_mass() )
    #print(cm)


def run(name,func,Nwarm,N):
    times = np.zeros(N)
    for i in range(Nwarm):
        func()
    for i in range(N):
        t = time.process_time()
        func()
        times[i] = time.process_time() - t
    print(f"{name}: {times.mean()}±{times.std()}")


run('align',align_bench,5,10)
run('within',within_bench,5,10)
run('trjconv',trjconv_bench,5,10)