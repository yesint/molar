from python import *
from sys import getrefcount
import numpy as np

sel = Source('../../tests/protein.pdb').select_str("resid 5:600")

def test1():
    sel = Source('../../tests/protein.pdb')("resid 5:600")
    trj = FileHandler('../../tests/protein.xtc')
    for st in trj:
        sel.set_state(st)
        print(st.time,sel.com())

def test2():
    fh = FileHandler('../../tests/protein.pdb')
    top = fh.read_topology()
    st = fh.read_state()
    src = Source(top,st)
    sel = src("resid 5:600")
    del top,st,src,fh
    pos0 = sel[0].pos
    pos1 = sel[1].pos
    pos2 = sel[2].pos
    print(f"ref: {getrefcount(sel)-1} {getrefcount(pos0)-1}")
    del pos2
    print(f"ref: {getrefcount(sel)-1} {getrefcount(pos0)-1}")
    del pos1
    print(f"ref: {getrefcount(sel)-1} {getrefcount(pos0)-1}")
    del sel
    print(f"ref: {getrefcount(pos0)-1} pos0: {pos0}")


def test3():
    print("[0] before:",sel[0].pos)
    print("[1] before]:",sel[1].pos)
    sel[0].pos = sel[1].pos
    print("[0] after [0]=[1]:",sel[0].pos)

    print("[0].name before:",sel[0].name)
    print("[1].name before:",sel[1].name)
    print("[0].name after [0]=[1]:",sel[0].name)


def test4():
    print(f'Size: {len(sel)}')
    for p in sel:
        print(p.id, p.pos)


def test5():
    print(sel[100].pos, sel[0].name)
    print(sel[-100].atom.name)
    print("x=",sel[-100].x)
    sel[-100].x = 42
    print("x=",sel[-100].x)
    

def test6():
    subsel = sel("name CA")
    print(len(sel), sel[0].name)
    print(len(subsel), subsel[0].name)

def test7():
    crd = sel.get_coord()
    print(sel[0].pos, crd[:,0])
    crd[0,0] = 42
    sel.set_coord(crd)
    print(sel[0].pos, crd[:,0])

    arr = np.zeros((3, len(sel)), dtype=np.float32)
    sel.set_coord(arr)
    print(sel[0].pos, crd[:,0])

    print(sel[5].atom.resname)

def test8():
    src = Source.from_file('../../tests/protein.pdb')
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
    print(b.shortest_vector([0.5,0.5,0.6]))

#test3()
#test2()
test1()
