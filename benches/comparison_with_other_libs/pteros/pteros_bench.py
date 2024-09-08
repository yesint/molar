from pteros import *
import time
import numpy as np


def trjconv_bench():
    s=System("tests/protein.pdb")
    s.load("tests/protein.xtc")
    sel = s("resid 560")
    sel.write("tests/.extracted.dcd")


def align_bench():
    s=System("tests/protein.pdb")
    ref_sel = s();
    cur_sel = s();
    rmsds = []
    def do_align(sys,fr):
        cur_sel.set_frame(fr)
        tr = fit_transform(cur_sel,ref_sel)
        cur_sel.apply_transform(tr)
        rmsds.append( rmsd(cur_sel,ref_sel) )
        return True
    s.load("tests/protein.xtc", on_frame=do_align)


def within_bench():
    s=System("tests/protein.pdb")
    sel = s("within 1.0 of resid 560")
    cm = []
    def do_cm(sys,fr):
        sel.set_frame(fr)
        cm.append( sel.center(True) )
        return True
    s.load("tests/protein.xtc", on_frame=do_cm)


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