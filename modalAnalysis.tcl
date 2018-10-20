# --------------------------------------
# Eigen value analysis
# --------------------------------------
# Copyright by Gerard J. O'Reilly, 2017
# Written: Gerard J. O'Reilly
# Date: February 2015
# Last Updated: March 2017
# --------------------------------------
proc modalAnalysis {numModes {pflag 0} {tag ""} {modaldir ""} {rdamp {0}} {xi 0.05}} {
# numModes 		Number of Modes (e.g. 2)
# pflag 		Print flag - 1 to show the modal info onscreen or 0 for nothing (default: 0)
# tag 		Tag to append to period and eigenvector filenames (e.g. "_Frame") (default "")
# modaldir 		Directory to put output files (e.g. "outs/") (default "")
# rdamp		Modes to which xi ratio of Rayleigh is applied are applied, specified as {a b} (e.g. {1 3} for the 1st and 3rd mode)
# xi			Ratio of Critical Damping to be applied to the previously listed modes (e.g. 0.05 for 5%)

# Please note: This only computes and spits out what the damping at each mode is, it doesnt actually apply it.
# This is just to see this information quickly while are computing modal properties here.

# If the rdamp list or the xi value are not input, then this part wont be conducted

# This makes the Prd and omega variables available outside the proc
upvar 1 Prd Prd
upvar 1 omega omega

# Solve for lambda
set lambda [eigen -genBandArpack $numModes]; # Default

# If this solver doesn't work, try another
if {[llength $lambda]==0} {set lambda [eigen -fullGenLapack $numModes]};
if {[llength $lambda]==0} {set lambda [eigen -symmBandLapack $numModes]};

# Record the eigenvectors
record

# Extract the eigenvalues to the appropriate arrays
set omega {}
set freq {}
set Prd {}
foreach lam $lambda {
	lappend omega [expr sqrt($lam)]
	lappend freq [expr sqrt($lam)/(2*3.14159)]
	lappend Prd [expr (2*3.14159)/sqrt($lam)]
}

# Print the Periods to a text file
file delete "${modaldir}Periods${tag}.txt"
set period "${modaldir}Periods${tag}.txt"
set Periods [open $period "w"];
for {set i 1} {$i<=$numModes} {incr i 1} {
	puts $Periods [format "%.3f" [lindex $Prd $i-1]]
}
close $Periods

# Print the eigenvectors to a text file
file delete "${modaldir}eigenVectors${tag}.txt"
print "${modaldir}eigenVectors${tag}.txt" -node


# Compute the Rayleigh damping % if requested
set n [llength $rdamp]
if {$n>1} {
	set wi [lindex $omega [lindex $rdamp 0]-1];
	set wj [lindex $omega [lindex $rdamp 1]-1];
	set a0 [expr $xi*2.0*$wi*$wj/($wi+$wj)];
	set a1 [expr 2.0*$xi/($wi+$wj)];
	set xi_modes {}
	for {set i 1} {$i<=$numModes} {incr i 1} {
		set wn [lindex $omega $i-1];
		lappend xi_modes [expr 0.5*($a0/$wn+$a1*$wn)];
	}
} elseif {$n==1} {
	lappend xi_modes 0.0
}

# Create the on screen output
if {$pflag==1} {
	for {set i 1} {$i<=$numModes} {incr i 1} {
		if {$n>1} {
			puts [format "Mode %d - T=%.2fs   f=%.2fHz  omega:%.1frad/s    xi=%.2f%%" $i [lindex $Prd $i-1] [lindex $freq $i-1] [lindex $omega $i-1] [expr [lindex $xi_modes $i-1]*100.0]]
		} else {
			puts [format "Mode %d - T=%.2fs   f=%.2fHz  omega:%.1frad/s " $i [lindex $Prd $i-1] [lindex $freq $i-1] [lindex $omega $i-1]]
		}
	}
}


set ::omega $omega
set ::freq $freq
set ::Prd $Prd

}
