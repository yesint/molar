mol new tests/albumin.pdb type pdb

# Format: "vmd selection"  "molar selection"
set selections {
    "name CA"   "name CA"
    "resid 10"  "resid 10"
    "same residue as (name CA and resid 10)"    "same residue as (name CA and resid 10)"
    "within 5 of resid 10"  "within 0.5 of resid 10"
    "within 3 of resid 20"  "within 0.3 of resid 20"
} 
#"within 10 of resid 1 to 200"  "within 1.0 of resid 1:200"
#"within 30 of resid 1 to 200"  "within 3.0 of resid 1:200"

set rust_code ""
set n 0
for {set i 0} {$i < [llength $selections]} {incr i 2} {
    set res [[atomselect top [lindex $selections $i]] get index]

    set s [lindex $selections [expr $i+1]]
    append rust_code "
    \#\[test\]
    fn selection_test_$n () {
        let answer: Vec<usize> =
        \"$res\"
        .split(\" \").map(|x| x.parse::<usize>().unwrap()).collect();
        assert_eq!(get_selection_index(\"$s\"),answer);
    }
    "
    incr n
}

set f [open "tests/generated_selection_tests.in" w]
puts $f $rust_code
close $f

exit