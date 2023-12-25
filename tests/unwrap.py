from pteros import *

s=System('triclinic.pdb')
prot = s('not resname TIP3 POT CLA')
l = prot.unwrap_bonds(0.2)
print(l)
