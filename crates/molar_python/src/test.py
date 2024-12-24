from molar_python import *

def test1():
    f = FileHandler.open('/home/semen/work/Projects/molar/tests/protein.pdb')
    topst = f.read()
    src = Source(*topst)

    sel1 = src.select_all()
    sel2 = src.select_str("resid 5:600")

    print("com1:",sel1.com())
    print("com2:",sel2.com())

    print("py1:",sel2.nth_pos(0))
    sel2.nth_pos(0)[0] = 42.0
    print("py2:",sel2.nth_pos(0))


def test2():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    n = sel.len()
    print(f"len: {n} {type(n)}")
    com = sel.com()
    print(f"com: {com} {type(com)}")
    pos = sel.nth_pos(0)
    print(f"py1: {pos} {type(pos)}")
    pos[0] = 42.0
    print("py2:",pos)
    del sel
    print("py3:",pos)


def test3():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    print("[0] before:",sel.nth_particle(0).pos)
    print("[1] before]:",sel.nth_particle(1).pos)
    sel.nth_particle(0).pos = sel.nth_particle(1).pos
    print("[0] after [0]=[1]:",sel.nth_particle(0).pos)

    print("[0].name before:",sel.nth_particle(0).name)
    print("[1].name before:",sel.nth_particle(1).name)
    print("[0].name after [0]=[1]:",sel.nth_particle(0).name)

def test4():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    print(f'Size: {sel.len()}')
    for p in sel:
        print(p.id, p.pos)

def test5():
    sel = Source.from_file('../../tests/protein.pdb').select_str("resid 5:600")
    print(sel[:].pos)
    

#test1()
#test2()
test5()
