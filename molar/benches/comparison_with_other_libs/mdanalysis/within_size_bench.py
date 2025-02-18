import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis import align, rms
import time
import numpy as np


def within_bench(sys,n_res,d):
    for start_res in range(0,100,10):
        sel = sys.select_atoms(f"around {d} resid {start_res} to {start_res+n_res}")


def within_pbc_bench(sys,n_res,d):
    for start_res in range(0,100,10):
        sel = sys.select_atoms(f"around {d} resid {start_res} to {start_res+n_res}",periodic=True)


sys = mda.Universe('tests/protein.pdb')
N = 1

def run(bench,sys,n_res,d,out):
    times = np.zeros(N)
    for i in range(N):
        t = time.process_time()
        bench(sys,n_res,d)
        times[i] = time.process_time() - t
    print(f"{n_res} {d} {times.mean()}±{times.std()}")
    out.write(f"{d} {times.mean()}\n")


for n_res in [1,20,40,60]:
    with open(f"target/mda_{n_res}.dat","w") as out:
        for d in [0.3+0.1*val for val in range(40)]:
            run(within_bench,sys,n_res,10.0*d,out)
    
    with open(f"target/mda_pbc_{n_res}.dat","w") as out:
        for d in [0.3+0.1*val for val in range(40)]:
            run(within_pbc_bench,sys,n_res,10.0*d,out)
            