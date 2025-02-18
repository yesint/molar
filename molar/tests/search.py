from pteros import *

s=System('no_ATP.pdb')
prot = s()
pairs,dist = search_contacts(0.3,prot,pbc=[1,1,1])
print(len(pairs))
