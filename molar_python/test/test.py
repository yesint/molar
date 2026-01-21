from pymolar import *
from sys import getrefcount
import numpy as np

#sel = System('../../molar/tests/protein.pdb').select("resid 5:600")

def test_set_state():
    top,st1 = FileHandler('../../molar/tests/protein.pdb','r').read()
    st2 = FileHandler('../../molar/tests/protein.pdb','r').read_state()
    st2.time = 100
    print(f'st1: {st1.time}, st2: {st2.time}')

    sys = System(top,st1)
    sel1 = sys('name CA')
    sel2 = sys('name CB')

    print(f'Before: sys: {sys.time}, sel1: {sel1.time}, sel2: {sel2.time}')

    old = sel1.set_state(st2)

    print(f'st1: {st1.time}, st2: {st2.time}')
    print(f'After: sys: {sys.time}, sel1: {sel1.time}, sel2: {sel2.time}')
    print(f'old: {old.time}')
    

def test1():
    sel = System('../../molar/tests/protein.pdb')("resid 5:600")
    trj = FileHandler('../../molar/tests/protein.xtc','r')
    for st in trj:
        sel.set_state(st)
        print(st.time,sel.com())

def test2():
    src = System('../../molar/tests/protein.pdb')
    sel = src("resid 5:600")
    pos0 = sel[0].pos
    print(f"ref: {getrefcount(pos0)} pos0: {pos0}")
    pos0[1]+=1
    print(f"ref: {getrefcount(pos0)} pos0: {sel[0].pos}")


def test3():
    src = System('../../molar/tests/protein.pdb')
    sel = src("resid 5:600")
    print("[0] before:",sel[0].pos)
    sel[0].pos = [100,100,3]
    print("[0] after [0]=[1]:",sel[0].pos)

    print("[0].name before:",sel[0].name,sel[1].name)
    sel[0].name = "AAA"
    print("[0].name after [0]=[1]:",sel[0].name, sel[1].name)

    print("[0].resid before:",sel[0].resid,sel[100].resid)
    sel[0].resid = sel[100].resid
    print("[0].resid after [0]=[1]:",sel[0].resid, sel[100].resid)

    print("[0].resid before:",sel[0].resname,sel[50].resname)
    sel[0].resname = sel[50].resname
    print("[0].resid after [0]=[1]:",sel[0].resname, sel[50].resname)

def test4():
    print(f'Size: {len(sel)}')
    for p in sel:
        print(p.id, p.pos)


def test5():
    src = System('../../molar/tests/protein.pdb')
    sel = src("resid 5:600")
    print(sel[100].pos, sel[0].name)
    print(sel[-100].atom.name)
    print("x=",sel[-100].x)
    sel[-100].x = 42
    print("x=",sel[-100].x)
    

def test6():
    src = System('../../molar/tests/protein.pdb')
    sel = src("resid 5:600")
    subsel = sel("name CA")
    print(len(sel), sel[0].name)
    print(len(subsel), subsel[0].name)

def test7():
    src = System('../../molar/tests/protein.pdb')
    sel = src("resid 5:600")
    crd = sel.get_coord()
    print("before:",sel[0].pos, crd[:,0])
    crd[0,0] = 42
    sel.set_coord(crd)
    print("after:",sel[0].pos, crd[:,0])

    arr = np.zeros((3, len(sel)), dtype=np.float32)
    sel.set_coord(arr)
    print(sel[0].pos, crd[:,0])

    print(sel[5].resname)

def test8():
    src = System('../../molar/tests/protein.pdb')
    sel = []
    sel.append( src("resid 500:600") )
    sel.append( src(None) )
    sel.append( src() )
    sel.append( src((0,199)) )

    sel.append( src([1,3,4,5,6,7]) )
    for s in sel:
        print(len(s))

def test9():
    #b = PeriodicBox([[1,0,0],[0,1,0],[0,0,1]])
    b = PeriodicBox([1,2,3],[90,90,90])
    print(b.to_vectors_angles())
    print(b.shortest_vector([0.9,0.5,0.6]))


def test_distance_search():
    s = System('../../molar/tests/topol.tpr')
    sel1 = s("resid 5:100")
    sel2 = s("resid 101:200")
    pairs,dist = distance_search('vdw',sel1,sel2)
    print(len(pairs),len(dist))
    print(pairs,dist)

def test_append():
    sys = System('../../molar/tests/protein.pdb')
    sel = sys("resid 550")
    sys.append(sel)


test7()
#test7()
#test2()
#test_distance_search()
#test_distance_search()
#test_set_state()
#test_append()