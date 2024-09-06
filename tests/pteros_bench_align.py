from pteros import *
import time

t = time.process_time()

s=System("protein.pdb")

ref_sel = s();
cur_sel = s();

rmsds = []

def do_align(sys,fr):
    cur_sel.set_frame(fr)
    tr = fit_transform(cur_sel,ref_sel)
    cur_sel.apply_transform(tr)
    rmsds.append( rmsd(cur_sel,ref_sel) )
    return True

s.load("protein.xtc", on_frame=do_align)

print(f"Elapsed: {time.process_time() - t}")
#print(f"{rmsds}")