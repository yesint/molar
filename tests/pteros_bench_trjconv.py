from pteros import *
import time

t = time.process_time()

s=System("protein.pdb")
s.load("protein.xtc")
sel = s("resid 560")

sel.write(".extracted.xtc")

print(f"Elapsed: {time.process_time() - t}")