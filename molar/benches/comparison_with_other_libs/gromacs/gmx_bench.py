from subprocess import call
import time
import numpy as np


def trjconv_bench():
    call("echo 3 | gmx trjconv -f tests/protein.xtc -s tests/protein.pdb -o target/extracted.xtc",shell=True)

def align_bench():
    call("echo 0 0 | gmx rms -f tests/protein.xtc -s tests/protein.pdb",shell=True)


def run(name,func,Nwarm,N):
    times = np.zeros(N)
    for i in range(Nwarm):
        func()
    for i in range(N):
        t = time.time_ns()
        func()
        times[i] = time.time_ns() - t
    times /= 1.0e9;
    print(f"{name}: {times.mean():.3}±{times.std():.3}")


run('align',align_bench,5,10)
#run('trjconv',trjconv_bench,5,10)