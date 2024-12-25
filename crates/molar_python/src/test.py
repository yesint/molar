from molar_python import *
from sys import getrefcount

def test1():
    f = FileHandler.open('/home/semen/work/Projects/molar/tests/protein.pdb')
    topst = f.read()
    src = Source(*topst)

    sel1 = src.select_all()
    sel2 = src.select_str("resid 5:600")

    print("com1:",sel1.com())
    print("com2:",sel2.com())

    print("py1:",sel2.nth_pos(0))
    sel2[0].pos[0] = 42
    print("py2:",sel2.nth_pos(0))


def test2():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    pos0 = sel.nth_pos(0)
    pos1 = sel.nth_pos(0)
    pos2 = sel.nth_pos(0)
    print(f"ref: {getrefcount(sel)-1} {getrefcount(pos0)-1}")
    del pos2
    print(f"ref: {getrefcount(sel)-1} {getrefcount(pos0)-1}")
    del pos1
    print(f"ref: {getrefcount(sel)-1} {getrefcount(pos0)-1}")
    del pos0
    print(f"ref: {getrefcount(sel)-1} {getrefcount(pos0)-1}")

def test3():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    print("[0] before:",sel[0].pos)
    print("[1] before]:",sel[1].pos)
    sel[0].pos = sel[1].pos
    print("[0] after [0]=[1]:",sel[0].pos)

    print("[0].name before:",sel[0].name)
    print("[1].name before:",sel[1].name)
    print("[0].name after [0]=[1]:",sel[0].name)


def test4():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    print(f'Size: {sel.len()}')
    for p in sel:
        print(p.id, p.pos)


def test5():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    print(sel[100].pos, sel[0].name)
    print(sel[-1000].atom.name)
    

def test6():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    subsel = sel("name CA")
    print(len(sel), sel[0].name)
    print(len(subsel), subsel[0].name)

#test3()
#test2()
test2()
