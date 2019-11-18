# ------------------------------------
# -- Script to Compute Sa(T,xi%) -----
# ------------------------------------
# Gerard O'Reilly
# EUCENTRE/IUSSPavia
# Created: November 2014
# Last Updated: November 2019
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
	# pga: 	Peak ground acceelration in g
	# pgv:	Peak ground velocity in m/s
	# pgd:	Peak ground displacement in m

	# Make the outputs available outside the procedure
	upvar 1 Sa Sa
	upvar 1 Sv Sv
	upvar 1 Sd Sd
	upvar 1 pga pga
	upvar 1 pgv pgv
	upvar 1 pgd pgd

	# Import the ground motion
	set imp [open $EQ "r"];
	set accg [read $imp];
	close $imp

	# Set up some basics
	set g 			9.81;		# Define gravity in metric, because this is not the stone age anymore....
	set gamma 	0.5;		# Newmark terms (Set for Average Acceleration method)
	set beta 		0.25;		# Newmark terms (Set for Average Acceleration method)
	set ms			1.0; 		# Set the mass to 1kg
	set acc 		{};
	set p 			{};
	set vg 			{0.0}
	set dg 			{0.0}

	for {set i 0} {$i < [llength $accg]} {incr i 1} {
		lappend acc [expr [lindex $accg $i]*$g]; 	# Change the units of the record to m/s2
		lappend p	[expr -$ms*[lindex $acc $i]];	# Create the force in N
		lappend vg [expr [lindex $vg $i]+0.5*$dt*$g*([lindex $accg $i+1]+[lindex $accg $i])]
	  lappend dg [expr [lindex $dg $i]+0.5*$dt*([lindex $vg $i+1]+[lindex $vg $i])]
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
	for {set i 0} {$i < [llength $accg]-1} {incr i 1} {
		lappend dp 		[expr [lindex $p $i+1]-[lindex $p $i]];
		lappend dp_bar 	[expr [lindex $dp $i]+$A*[lindex $v $i]+$B*[lindex $a $i]];
		lappend du		[expr [lindex $dp_bar $i]/$k_bar];
		lappend dv		[expr $gamma*[lindex $du $i]/$beta/$dt-$gamma*[lindex $v $i]/$beta+$dt*(1-$gamma/2/$beta)*[lindex $a $i]];
		lappend da		[expr [lindex $du $i]/$beta/$dt/$dt-[lindex $v $i]/$beta/$dt-[lindex $a $i]/2/$beta];
		lappend u		[expr [lindex $u $i]+[lindex $du $i]];
		lappend v		[expr [lindex $v $i]+[lindex $dv $i]];
		lappend a		[expr [lindex $a $i]+[lindex $da $i]];
	}

	# I'm not coding the other stuff for getting total acceleration etc. like the Matlab function, dont need it here
	set umax 0.0;
	set pga 0.0;
	set pgv 0.0;
	set pgd 0.0;

	# Loop the terms (This loop is t length)
	for {set i 0} {$i < [llength $accg]} {incr i 1} {
		set temp1 [expr abs([lindex $u $i])]; 		# Get the absolute current displacement
		set temp2 [expr abs([lindex $accg $i])]; 	# Get the absolute current ground acceleration
		set temp3 [expr abs([lindex $vg $i])]; 		# Get the absolute current ground velocity
	  set temp4 [expr abs([lindex $dg $i])]; 		# Get the absolute current ground displacement

		if {$temp1>$umax} {set umax $temp1};		# Find the peak displacement
		if {$temp2>$pga} {set pga $temp2};			# Find the pga
		if {$temp3>$pgv} {set pgv $temp3};			# Find the pgv
	  if {$temp4>$pgd} {set pgd $temp4};			# Find the pgd

	}

	# Calculate Spectral Values
	set Sd 	$umax;			# Spectral acceleration in m
	set Sv	[expr $Sd*$w];		# Pseudo spectral velocity in m/s
	set Sa	[expr $Sd*$w*$w/$g];	# Pseudo spectral acceleration in g

}
