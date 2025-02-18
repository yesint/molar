# Function to perform a selection within a distance `d` of residues
proc within_bench {n_res d} {
    for {set start_res 0} {$start_res < 100} {incr start_res 10} {
        puts "within $d of resid $start_res to [expr $start_res + $n_res]"
        set sel [atomselect top "within $d of resid $start_res to [expr $start_res + $n_res]"]
        $sel delete
    }
}

# Function to run the benchmark and measure execution time
proc run {n_res d out_file} {
    set start_time [clock milliseconds]
    within_bench $n_res $d
    set end_time [clock milliseconds]
    set t [expr ($end_time - $start_time)/1000.0]
    puts "$n_res $d $t"
    puts $out_file "[expr $d/10.0] $t"
}

# Main script execution
mol new tests/albumin.pdb

# List of residue counts and distances
set res_counts {1 20 40 60}

foreach n_res $res_counts {
    set out_file [open "target/vmd_${n_res}.dat" "w"]
    for {set d 3.0} {$d <= 43.0} {set d [expr {$d + 1.0}]} {
        run $n_res $d $out_file
    }
    close $out_file
}


# for {set d 3.0} {$d <= 50.0} {set d [expr {$d + 1.0}]} {
#     set start_time [clock milliseconds]

#     set sel [atomselect top "within $d of residue 0 to 200"]
    
#     set end_time [clock milliseconds]
    
#     set t [expr ($end_time - $start_time)/1000.0]
#     puts "$d $t [$sel num]"
# }

# Quit VMD after completion
quit
