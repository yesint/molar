# VMD Tcl Script to Calculate Center of Mass for Nearby Atoms Using Pre-Created Selection and Callbacks
# Start high-precision timing in milliseconds
set t_start [clock milliseconds]

# Load the reference structure (PDB file)
mol new protein.pdb type pdb

# Pre-create the selection of atoms within 10 Å of residue 560
# The 'update' keyword allows this selection to refresh for each frame
set nearby_sel [atomselect top "within 10.0 of resid 560"]

# Initialize a list to store centers of mass
set cm_list {}


# Callback function to process each frame
proc process_frame {} {
    global cm_list nearby_sel

    # Update the selection to reflect the current frame
    $nearby_sel update

    # Compute the center of mass of the selected atoms
    if {[$nearby_sel num] > 0} {
        set com [$nearby_sel measure center weight mass]
        lappend cm_list $com
    }
}

# Set up the callback to call `process_frame` whenever a new frame is loaded
trace variable ::vmd_frame_change w process_frame

# Load the trajectory file, triggering the callback on each frame
mol addfile protein.xtc type xtc waitfor all

# Disable the callback after the trajectory is processed
trace vdelete ::vmd_frame_change w process_frame

# End high-precision timing in milliseconds
set t_end [clock milliseconds]

# Calculate elapsed time in seconds with fractional precision
set elapsed_time [expr {($t_end - $t_start) / 1000.0}]

# Output the list of centers of mass and elapsed time
puts "Centers of Mass: $cm_list"
puts "Elapsed: $elapsed_time seconds"
exit