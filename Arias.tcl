# This is a function to compute the Arias Intensity of a ground motion
proc Arias {EQ dt I1 I2} {
	upvar 1 I I
	upvar 1 t12 t12

	# Import the ground motion
	set imp [open $EQ "r"];
	set accg [read $imp];
	close $imp

	# Create time vector
	set t {0.0};
	for {set i 1} {$i < [llength $accg]} {incr i 1} {
		lappend t [expr [lindex $t $i-1]+$dt];
	}

	# Compute the Arias Intensity
	set Ia {0.0};
	for {set i 1} {$i < [llength $accg]} {incr i 1} {
	    lappend Ia [expr [lindex $Ia $i-1]+$dt*0.5*(pow([lindex $accg $i],2.0)+pow([lindex $accg $i-1],2.0))*3.14/2.0/9.81];
    	}

	# Get the Arias intensit of the record
	set I [lindex $Ia [llength $accg]-1];

	# Find the time span between the specified percentages of I
	for {set i 1} {$i <= [llength $Ia]} {incr i 1} {
		if {[expr $I1*$I]>[lindex $Ia $i-1]} {set c1 $i};
		if {[expr $I2*$I]>[lindex $Ia $i-1]} {set c2 $i};
	}
	set t1 [lindex $t $c1-1];
	set t2 [lindex $t $c2-1];
	set t12 [expr $t2-$t1];

}
