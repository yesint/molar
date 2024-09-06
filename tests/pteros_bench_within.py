from pteros import *
import time

t = time.process_time()

s=System("protein.pdb")

sel = s("within 1.0 of resid 560");

cm = []

def do_cm(sys,fr):
    sel.set_frame(fr)
    cm.append( sel.center(True) )
    return True

s.load("protein.xtc", on_frame=do_cm)

print(f"{cm}")
print(f"Elapsed: {time.process_time() - t}")