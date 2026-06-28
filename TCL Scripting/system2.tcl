# Open output file once
set out [open "PR_434_ssp.dat" w]

# Set scale for phi/psi normalization
set scale 180.0

# Load reference and trajectory system
set sysref [mol new PR_Protein_bgc_434_center.gro]
set sys [mol new PR_Protein_bgc_434_center.gro]
mol addfile PR_Protein_bgc_434_center.xtc first 0 last -1 step 1 waitfor all

# Get total number of frames
set nf [molinfo $sys get numframes]
puts "Total frames: $nf"

# Divide into 3 segments
set nf1 [expr {int($nf / 3)}]
set nf2 [expr {int(($nf * 2) / 3)}]
set nf3 $nf

# Frame segments
set segments [list [list 0 $nf1] [list [expr {$nf1 + 1}] $nf2] [list [expr {$nf2 + 1}] $nf3]]

# Loop over segments
foreach seg $segments {
    set start [lindex $seg 0]
    set end [lindex $seg 1]
    
    puts "Processing segment from frame $start to $end"

    for {set f $start} {$f <= $end} {incr f 10} {
        molinfo $sys set frame $f
        set t [expr {($f*1.0)/1.0}]
        puts "Frame: $t"

        set sum 0.0
        set fres_b1 1
        set lres_b1 418
        set nres_b1 [expr {($lres_b1 - $fres_b1) + 1}]

        for {set d1 1} {$d1 <= 418} {incr d1} {
            set sel_ref1 [atomselect $sysref "resid $d1 and name CA"]
            set n [$sel_ref1 get resname]

            if {$n != "GLY"} {
                set ref_phi1 [$sel_ref1 get phi]
                set ref_psi1 [$sel_ref1 get psi]
                $sel_ref1 delete

                set sel_b1 [atomselect $sys "resid $d1 and name CA"]
                set curr_phi1 [$sel_b1 get phi]
                set curr_psi1 [$sel_b1 get psi]

                set delphi1 [expr {$curr_phi1 - $ref_phi1}]
                set delpsi1 [expr {$curr_psi1 - $ref_psi1}]
                set mphi1 [expr {$delphi1 / $scale}]
                set mpsi1 [expr {$delpsi1 / $scale}]
                set np1 [expr {abs($mphi1)}]
                set ns1 [expr {abs($mpsi1)}]

                set sum [expr {$sum + (exp(-1*$np1) * exp(-1*$ns1))}]

                $sel_b1 delete

                # Clean up
                unset ref_phi1 ref_psi1 curr_phi1 curr_psi1 delphi1 delpsi1 mphi1 mpsi1 np1 ns1
            } else {
                $sel_ref1 delete
            }
        }

        # Average score for this frame
        set pval [expr {$sum / (1.0 * $nres_b1)}]
        puts $out "$t $pval"
    }
}

close $out
exit

