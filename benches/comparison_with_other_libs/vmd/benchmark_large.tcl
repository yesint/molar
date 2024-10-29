source ./bigdcd.tcl

# Function to compute mean and standard deviation
proc compute_stats {times} {
    set n [llength $times]
    set mean [expr {[tcl::mathop::+ {*}$times] / double($n)}]

    set variance 0
    foreach time $times {
        set variance [expr {$variance + ($time - $mean) * ($time - $mean)}]
    }
    set std [expr {sqrt($variance / double($n))}]

    return [list $mean $std]
}

# Run given benchmark
proc runner {name func n_warming n_bench} {
    puts "Benchmarking $name"
    puts "Warming up..."
    # Warmng up
    for {set i 0} {$i < $n_warming} {incr i} {
        $func
    }
    puts "Collecting statistics..."
    # Benchmarking
    set times {}
    for {set i 0} {$i < $n_bench} {incr i} {
        set start [clock seconds]
        $func
        set end [clock seconds]
        lappend times [expr {$end - $start}]
    }

    # Compute statistics
    set stats [compute_stats $times]
    puts "$name: [format %.6f [lindex $stats 0]]±[format %.6f [lindex $stats 1]]"
}

#--------------------------------------
# Align
#--------------------------------------

proc bench_align_on_frame {frame} {
    global cur_sel ref_sel rmsd_list
    # Select the protein atoms in the current loaded frame
    $cur_sel frame $frame
    # Align the selection to the reference
    $cur_sel move [measure fit $cur_sel $ref_sel]
    # Calculate RMSD after alignment
    set rmsd_val [measure rmsd $cur_sel $ref_sel]
    # Append RMSD value to the list
    lappend rmsd_list $rmsd_val
}

proc bench_align {} {
    set PDB "inp.pdb"
    set XTC "traj_comp.xtc"

    # Load strutcure
    mol new $PDB waitfor all
    # Select the protein atoms from the reference structure
    set ::ref_sel [atomselect top "protein" frame 0]
    set ::cur_sel [atomselect top "protein" frame 0]
    set ::rmsd_list {}
    # Read trajectory frame bu frame and call benchmark function
    bigdcd bench_align_on_frame xtc $XTC
    bigdcd_wait
    #puts $::rmsd_list
}

#-------------------------------------------
# Within
#--------------------------------------

proc bench_within_on_frame {frame} {
    global within_sel cm_list
    $within_sel frame $frame
    lappend cm_list [measure center $within_sel weight mass]
}

proc bench_within {} {
    set PDB "inp.pdb"
    set XTC "traj_comp.xtc"

    # Load strutcure
    mol new $PDB waitfor all
    # Select the protein atoms from the reference structure
    set ::within_sel [atomselect top "within 10.0 of protein" frame 0]
    set ::cm_list {}
    # Read trajectory frame by frame and call benchmark function
    bigdcd bench_within_on_frame xtc $XTC
    bigdcd_wait
    unset ::within_sel
    #puts $::rmsd_list
}

#-------------------------------------------
# trjconv
#--------------------------------------
proc bench_trjconv {} {
    set PDB "inp.pdb"
    set XTC "traj_comp.xtc"

    mol new $PDB
    mol addfile $XTC waitfor all
    set sel [atomselect top "protein"]
    $sel writedcd ../../../target/extracted.dcd
}

#-------------------------------------------
# Run benchmarks
# Warm up for 5 iterations
#-------------------------------------------
runner "align" bench_align 2 10
#runner "within" bench_within 5 10
#runner "trjconv" bench_trjconv 0 2

exit
