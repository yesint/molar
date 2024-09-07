# Callback function to handle each frame
proc process_frame_align {} {
    global ref_sel rmsd_list current_frame

    # Select the protein atoms in the current loaded frame
    set sel [atomselect top "protein" frame [molinfo top get frame]]
    # Align the selection to the reference
    $sel move [measure fit $sel $ref_sel]
    # Calculate RMSD after alignment
    set rmsd_val [measure rmsd $sel $ref_sel]
    # Append RMSD value to the list
    lappend rmsd_list $rmsd_val
    # Increment frame counter
    incr current_frame
}

proc bench_align {} {
    # Load the reference structure (PDB file)
    mol new tests/protein.pdb type pdb
    # Select the protein atoms from the reference structure
    set ref_sel [atomselect top "protein" frame 0]
    # Create a list to store RMSD values
    set rmsd_list {}
    # Variable to track the frame index
    set current_frame 0
    # Set up the callback to call `process_frame` whenever a new frame is loaded
    trace variable ::vmd_frame_change w process_frame_align
    # Load the trajectory file, triggering the callback on each frame
    mol addfile tests/protein.xtc type xtc waitfor all
    # Disable the callback after the trajectory is processed
    trace vdelete ::vmd_frame_change w process_frame
}

#----------------------------------------------------

# Callback function to handle each frame
proc process_frame_within {} {
    global cm_list nearby_sel
    # Update the selection to reflect the current frame
    $nearby_sel update
    # Compute the center of mass of the selected atoms
    if {[$nearby_sel num] > 0} {
        set com [$nearby_sel measure center weight mass]
        lappend cm_list $com
    }
}

proc bench_within {} {
    mol new tests/protein.pdb type pdb
    set nearby_sel [atomselect top "within 10.0 of resid 560"]
    # Initialize a list to store centers of mass
    set cm_list {}
    # Set up the callback to call `process_frame` whenever a new frame is loaded
    trace variable ::vmd_frame_change w process_frame
    # Load the trajectory file, triggering the callback on each frame
    mol addfile tests/protein.xtc type xtc waitfor all
    # Disable the callback after the trajectory is processed
    trace vdelete ::vmd_frame_change w process_frame_within
}
#----------------------------------------------------

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


proc run {name func} {
    set times {}
    for {set i 0} {$i < 3} {incr i} {
        set start [clock seconds]
        $func
        set end [clock seconds]
        lappend times [expr {$end - $start}]
    }

    # Compute and display the statistics
    set stats [compute_stats $times]
    puts "$name: [format %.6f [lindex $stats 0]]±[format %.6f [lindex $stats 1]]"
}

run "align" bench_align
#run "within" bench_within

exit