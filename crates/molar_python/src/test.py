from molar_python import *

def test1():
    f = FileHandler.open('../molar/tests/protein.pdb')
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
    sel = Source.from_file('../molar/tests/protein.pdb').select_str("resid 5:600")
    n = sel.len()
    print(f"len: {n} {type(n)}")
    com = sel.com()
    print(f"com: {com} {type(com)}")
    pos = sel.nth_pos(0)
    print(f"py1: {pos} {type(pos)}")
    pos[0] = 42.0
    print("py2:",pos)


def test3():
    sel = Source.from_file('../molar/tests/protein.pdb').select_str("resid 5:600")
    print("py1:",sel.nth_pos(100))

#test1()
test2()
#test3()
