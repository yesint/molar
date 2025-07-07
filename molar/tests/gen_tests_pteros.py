from pteros import *
sys = System('tests/albumin.pdb')

# Format: "pteros selection"  "molar selection"
selections = [
    "name CA",   "name CA",
    "resid 555",  "resid 555",
    "by residue (name CA and resid 555)",    "same residue as (name CA and resid 555)",
    "within 0.5 of resid 555",  "within 0.5 of resid 555",
    "within 0.5 pbc of resid 555",  "within 0.5 pbc yyy of resid 555",
    "x<4", "x<4",
    "dist from 4.1 3.9 5.1 < 0.4", "dist point 4.1 3.9 5.1<0.4", 
]

rust_code = ""
n = 0
#for {set i 0} {$i < [llength $selections]} {incr i 2} {
for i in range(0,len(selections),2):
    print(selections[i])
    res = sys(selections[i]).get_index()
    s = selections[i+1]
    code="""
    #[test]
    fn selection_pteros_test_{n} () {{
        let answer: Vec<usize> = vec!{res};
        assert_eq!(get_selection_index2(\"{s}\"),answer);
    }}
    """.format(n=n,res=res,s=s)
    rust_code += code
    n+=1


with open("tests/generated_pteros_tests.in","w") as f:
    f.write(rust_code)
