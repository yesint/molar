import mdtraj as md
import time
import numpy as np


def trjconv_bench():
    trajectory = md.load('tests/protein.xtc', top='tests/protein.pdb')
    atom_indices = trajectory.topology.select("residue 560")
    selected_trajectory = trajectory.atom_slice(atom_indices)
    selected_trajectory.save_dcd('target/.extracted.dcd')


def align_bench():
    reference = md.load('tests/protein.pdb')
    trajectory = md.load('tests/protein.xtc', top='tests/protein.pdb')
    atom_indices = reference.topology.select('all')
    trajectory.superpose(reference, atom_indices=atom_indices)
    rmsd = md.rmsd(trajectory, reference, atom_indices=atom_indices)


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
run('trjconv',trjconv_bench,5,10)