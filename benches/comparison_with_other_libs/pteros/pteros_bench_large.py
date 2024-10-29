from pteros import *
import time
import numpy as np

PDB = "inp.pdb"
XTC = "traj_comp.xtc"


def trjconv_bench():
    s=System(PDB)
    s.load(XTC,e=500)
    sel = s("protein")
    sel.write("tests/.extracted.dcd")


def align_bench():
    s=System(PDB)
    ref_sel = s('protein');
    cur_sel = s('protein');
    rmsds = []
    def do_align(sys,fr):
        cur_sel.set_frame(fr)
        tr = fit_transform(cur_sel,ref_sel)
        cur_sel.apply_transform(tr)
        rmsds.append( rmsd(cur_sel,ref_sel) )
        return True
    s.load(XTC, on_frame=do_align, e=500)


def within_bench():
    s=System(PDB)
    sel = s("within 1.0 of protein")
    cm = []
    def do_cm(sys,fr):
        sel.set_frame(fr)
        cm.append( sel.center(True) )
        return True
    s.load(XTC, on_frame=do_cm, e=500)


def run(name,func,Nwarm,N):
    times = np.zeros(N)
    for i in range(Nwarm):
        func()
    for i in range(N):
        t = time.process_time()
        func()
        times[i] = time.process_time() - t
    print(f"{name}: {times.mean()}±{times.std()}")


#run('align',align_bench,2,10)
run('within',within_bench,1,5)
#run('trjconv',trjconv_bench,1,5)