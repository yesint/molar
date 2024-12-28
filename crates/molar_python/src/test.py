from molar_python import *
from sys import getrefcount
import numpy as np

sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")

def test1():
    topst = f.read()
    src = Source(*topst)

    sel1 = src.select_all()
    sel2 = src.select_str("resid 5:600")

    print("com1:",sel1.com())
    print("com2:",sel2.com())

    print("py1:",sel2[0].pos)
    sel2[0].pos[0] = 42
    print("py2:",sel2[0].pos)


def test2():
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

#test3()
#test2()
test7()
