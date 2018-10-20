# ------------------------------------
# -- Script to Compute Sa(T,xi%) -----
# ------------------------------------
# Gerard O'Reilly
# EUCENTRE/IUSSPavia
# Created: November 2014
# Last Updated: February 2015
#
# This is a script that will return the Sa(T1) of a given record,
# for a specified value of period T using the Newmark Average
# Acceleration. This is based on the Newmark.m matlab function.

# -----------------------------------
# Inputs:
# -----------------------------------
# EQ:   Filename which is a single column file in units g (e.g "Eq.txt")
# dt:   Time step in seconds (e.g 0.01)
# xi:   Elastic damping (e.g 0.05)
# T:    Period in seconds (e.g 1.0)

# -----------------------------------
# Outputs:
#------------------------------------
# Sa:	Sa(T,%) - Pseudo-Spectral Acceleration in g
# Sv: 	Sv(T,%) - Pseudo-Spectral Velocity in m/s
# Sd: 	Sd(T,%) - Spectral Displacement in m

proc getSaT {EQ dt T xi} {
	# -----------------------------------
	# Inputs:
	# -----------------------------------
	# EQ:   Filename which is a single column file in units g (e.g "Eq.txt")
	# dt:   Time step in seconds (e.g 0.01)
	# xi:   Elastic damping (e.g 0.05)
	# T:    Period in seconds (e.g 1.0)

	# -----------------------------------
	# Outputs:
	#------------------------------------
	# Sa:	Sa(T,%) - Pseudo-Spectral Acceleration in g
	# Sv: 	Sv(T,%) - Pseudo-Spectral Velocity in m/s
	# Sd: 	Sd(T,%) - Spectral Displacement in m

	# Make the outputs available outside the procedure
	upvar 1 Sa Sa
	upvar 1 Sv Sv
	upvar 1 Sd Sd
	upvar 1 pga pga

	# Import the ground motion
	set imp [open $EQ "r"];
	set accg [read $imp];
	close $imp

	# Set up some basics
	set g 	9.81;		# Define gravity in metric, because this is not the stone age anymore....
	set gamma 	0.5;		# Newmark terms (Set for Average Acceleration method)
	set beta 	0.25;		# Newmark terms (Set for Average Acceleration method)
	set ms	1.0; 		# Set the mass to 1kg
	set acc 	{};
	set p 	{};
	for {set i 1} {$i <= [llength $accg]} {incr i 1} {
		lappend acc [expr [lindex $accg $i-1]*$g]; 	# Change the units of the record to m/s2
		lappend p	[expr -$ms*[lindex $acc $i-1]];	# Create the force in N
	}

	# Create time vector
	set t {0.0};
	for {set i 1} {$i < [llength $acc]} {incr i 1} {
		lappend t [expr [lindex $t $i-1]+$dt];
	}

	# Calculate the initial values
	set k 		[expr $ms*pow((2.0*3.141592654/$T),2)];		# Stiffness in N/m (which will give T following assumption of mass)
	set w 		[expr pow(($k/$ms),0.5)];				# Circular frequency
	set c 		[expr 2*$xi*$ms*$w];					# Damping coefficient
	set a0 		[expr [lindex $p 0]/$ms];				# Initial acceleration in m/s2
	set k_bar	[expr $k+$gamma*$c/$beta/$dt+$ms/$beta/$dt/$dt];	# Stiffness term (see Chopra book)
	set A		[expr $ms/$beta/$dt+$gamma*$c/$beta];			# Constant term (see Chopra book)
	set B		[expr $ms/2/$beta+$dt*$c*($gamma/2/$beta-1)];		# Constant term (see Chopra book)

	# Initialise some vectors
	set u {0.0};
	set v {0.0};
	set a {$a0};
	set du {};
	set dv {};
	set da {};
	set dp {};
	set dp_bar {};

	# Loop the terms (This loop is to length minus 1)
	for {set i 1} {$i < [llength $t]} {incr i 1} {
		lappend dp 		[expr [lindex $p $i]-[lindex $p $i-1]];
		lappend dp_bar 	[expr [lindex $dp $i-1]+$A*[lindex $v $i-1]+$B*[lindex $a $i-1]];
		lappend du		[expr [lindex $dp_bar $i-1]/$k_bar];
		lappend dv		[expr $gamma*[lindex $du $i-1]/$beta/$dt-$gamma*[lindex $v $i-1]/$beta+$dt*(1-$gamma/2/$beta)*[lindex $a $i-1]];
		lappend da		[expr [lindex $du $i-1]/$beta/$dt/$dt-[lindex $v $i-1]/$beta/$dt-[lindex $a $i-1]/2/$beta];
		lappend u		[expr [lindex $u $i-1]+[lindex $du $i-1]];
		lappend v		[expr [lindex $v $i-1]+[lindex $dv $i-1]];
		lappend a		[expr [lindex $a $i-1]+[lindex $da $i-1]];
	}

	# I'm not coding the other stuff for getting total acceleration etc. like the Matlab function, dont need it here
	set umax 0.0;
	set pga 0.0;
	# Loop the terms (This loop is to length)
	for {set i 1} {$i <= [llength $t]} {incr i 1} {
		set temp1 [expr abs([lindex $u $i-1])]; 		# Get the absolute current displacement
		set temp2 [expr abs([lindex $accg $i-1])]; 	# Get the absolute current ground acceleration
		if {$temp1>$umax} {set umax $temp1};		# Find the peak displacement
		if {$temp2>$pga} {set pga $temp2};			# Find the pga

	}

	# Calculate Spectral Values
	set Sd 	$umax;			# Spectral acceleration in m
	set Sv	[expr $Sd*$w];		# Pseudo spectral velocity in m/s
	set Sa	[expr $Sd*$w*$w/$g];	# Pseudo spectral acceleration in g

}
