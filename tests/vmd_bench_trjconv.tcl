# Start high-precision timing in milliseconds
set t_start [clock milliseconds]

# Load the reference structure (PDB file)
mol new protein.pdb type pdb

# Load the trajectory file (entirely) using the same topology
mol addfile protein.xtc type xtc waitfor all

# Select the protein atoms from the loaded trajectory based on indices
set sel [atomselect top "resid 560"]

# Write the selected atoms to a new trajectory
animate write xtc extracted.xtc sel $sel

# End high-precision timing in milliseconds
set t_end [clock milliseconds]

# Calculate elapsed time in seconds with fractional precision
set elapsed_time [expr {($t_end - $t_start) / 1000.0}]

# Output the elapsed time
puts "Elapsed: $elapsed_time seconds"
exit