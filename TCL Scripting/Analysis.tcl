# Load trajectory

mol new PR_protein_bgc_41_center.gro
mol addfile PR_protein_bgc_41_center.xtc first 0 last -1 step 10 waitfor all


set reference [atomselect top "protein and backbone" frame 0]

# the frame to be  compared
set compare [atomselect top "protein and backbone"]

#get the number of frames of the trajectory 

set num_steps [molinfo top get numframes]

# get the correct num as the trajectory start from zero

##set outfile [open rmsd2o9t_500new.dat w]
set outfile [open PR_protein_bgc_41_rmsd.dat w]

for {set frame 0} {$frame < $num_steps} {incr frame} {
                # get the correct frame
                $compare frame $frame
                # compute the 4*4 matrix transformation that takes one set of coordinates onto the other 
                set trans_mat [measure fit $compare $reference]

                # do the alignment
                $compare move $trans_mat
                # compute the RMSD
                set rmsd [measure rmsd $compare $reference ]
                # print the RMSD
   puts $outfile "$frame    $rmsd"
}
close $outfile

#RMSF Calculation 
set num [expr {$num_steps - 1}]
# RMSF calculation
##set outfile [open .dat w]
set outfile [open PR_protein_bgc_41_rmsf.dat w]

set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
for {set i 0} {$i < [$sel num]} {incr i} {
  puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
}
close $outfile

#-detailout <details output file> (default: stdout)

## SASA 
# selection
set sel [atomselect top "backbone"]
set n [molinfo top get numframes]
##set output [open "SASA2o9t_500new.dat" w]
set output [open "PR_protien_bgc_41_sasa.dat" w]
# sasa calculation loop
for {set i 0} {$i < $n} {incr i} {
        molinfo top set frame $i
        set sasa [measure sasa 1.4 $sel -restrict $sel]
        puts "\t \t progress: $i/$n"
        puts $output "$i $sasa"
}
puts "\t \t progress: $n/$n"
puts "Done."
puts "output file: SASA_O08_95mM.dat"
close $output

source rog_loop_1us.tcl
