
source gyr_radius.tcl
source center_of_mass.tcl

set outfile [open PR_protein_bgc_41_rog.dat w]
puts $outfile "i rad_of_gyr"
set nf [molinfo top get numframes]
set i 0

set prot [atomselect top "backbone"]
while {$i < $nf} {

    $prot frame $i
    $prot update

    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]

    puts $outfile "$i $rog"

}

close $outfile
#exit
