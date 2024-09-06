# VMD Tcl Script to Align and Calculate RMSD using Callbacks for Trajectory Frames
set t_start [clock milliseconds]
# Load the reference structure (PDB file)
mol new protein.pdb type pdb

# Select the protein atoms from the reference structure
set ref_sel [atomselect top "protein" frame 0]
set ref_coords [$ref_sel get {x y z}]

# Create a list to store RMSD values
set rmsd_list {}

# Variable to track the frame index
set current_frame 0

# Callback function to handle each frame
proc process_frame {} {
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

# Set up the callback to call `process_frame` whenever a new frame is loaded
trace variable ::vmd_frame_change w process_frame

# Load the trajectory file, triggering the callback on each frame
mol addfile protein.xtc type xtc waitfor all

# Disable the callback after the trajectory is processed
trace vdelete ::vmd_frame_change w process_frame

# Output RMSD values
puts "RMSD values: $rmsd_list"

# Calculate elapsed time
set elapsed_time [expr {([clock milliseconds] - $t_start) / 1000.0}]
puts "Elapsed: $elapsed_time seconds"
exit