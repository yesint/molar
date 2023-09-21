mol new triclinic.pdb type pdb

set selections {
    "name CA"
    "resid 10"
    "same residue as (name CA and resid 10)"
} 

set rust_code ""
set n 0
foreach s $selections {
    set res [[atomselect top $s] get index]

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

set f [open "generated_selection_tests.in" w]
puts $f $rust_code
close $f

exit