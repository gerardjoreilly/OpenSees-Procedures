# ----------------------------------------------------------------
# -- Script to Conduct 3D Non-linear Response History Analysis ---
# ----------------------------------------------------------------
# Copyright by Gerard J. O'Reilly, 2017
# EUCENTRE and IUSS Pavia, Italy
# Date Created: April 2012
# Last Updated: May 2017

# This procedure is a simple script that executes the NRHA of a 3D model. It
# requires that the model has the dynamic analysis objects defined and just the
# 'analyze' of a regular OpenSees model needs to be executed. Therefore, the model
# needs to be built to the point where a modal analysis can be conducted. The
# ground motion timeSeries and pattern need to be setup and the constraints,
# numberer and system analysis objects also need to be predefined.

# When conducting the NRHA, this proc will try different options to achieve
# convergence and finish the ground motion. This allows for a relatively robust
# analysis algorithm to be implemented with a single command.

# In addition, the analysis requires that a deformation capacity be specified
# to define a limit that upon exceedance, will stop the NRHA and mark the
# analysis as a collapse case. It monitors the current deformation of a number
# of specified nodes and flags collapse based on their deforormation. This
# prevents the models getting 'stuck' trying to converge a model that has
# essentially collapsed, or would be deemed a collapse case in post processing.
# These are defined in terms of both the local storey drifts and also the roof
# drift in either direction. For 3D analysis, the absolute maximum drift in
# either direction is used. Other definitions are possible but not yet implemented.

# Lastly, a log file identifier is also required in the input. This will log
# all of the essential information about the maximum storey drifts. This script
# was developed for analysing buildings so the deformation capacity typically
# corresponds to a drift capacity and the top and bottom nodes would typically
# correspond to the centreline nodes of the floorslabs.

# WARNING: The acceleration values that are ouput from this are the relative
# accelerations. For some reason, you cant get the current absolute values with
# inline commands.

proc runNRHA3D {Dt Tmax Dc tNode bNode log {pflag 0}} {
	if {$pflag>0} {puts "Starting runNRHA"};
	# --------------------------------------------------
	# Description of Parameters
	# --------------------------------------------------
	# Dt:		Analysis time step
	# Tmax:	Length of the record (including padding of 0's)
	# Dc:		Drift capacity for both storey and roof drift (%)
	# tNode:	List of top nodes (e.g. {2 3 4 5})
	# bNode:	List of bottom node (e.g. {1 2 3 4})
	# log:	File handle of the logfile
	# --------------------------------------------------

	# Make the control index available outside the proc
	upvar 1 cIndex cIndex

	# Define the initial Analysis parameters
	set 	testType	NormDispIncr;	# Set the test type
	set 	tolInit 	1.0e-7;		# Set the initial tolerance, so it can be referred back to
	set 	iterInit 	50;			# Set the initial max number of iterations
	set 	algorithmType KrylovNewton;	# Set the algorithm type
	test 		$testType $tolInit $iterInit
	algorithm 	$algorithmType
	integrator 	Newmark 0.5 0.25
	analysis 	Transient

	# Set up analysis parameters
	set Nsteps [expr int($Tmax/$Dt)]; 	# integer of the number of steps needed
	set cIndex 		0;			# Initially define the control index (-1 for non-converged, 0 for stable, 1 for global collapse, 2 for local collapse)
	set controlTime 	0.0;			# Start the controlTime
	set ok		0;			# Set the convergence to 0 (initially converged)
	set mdrft 		0.0;			# Set the initial storey drift
	set mflr		0;			# Set the initial storey collapse location

	# Set up the storey drift and acceleration values
	set h 		{};
	set mdrftX 		{};
	set mdrftY 		{};
	set maccelX 	{0.0};
	set maccelY 	{0.0};
	set bdg_h 		0.0;
	for {set i 1} {$i<=[llength $tNode]} {incr i 1} {
		# Find the coordinates of the nodes in Global Z (3)
		set top2	[nodeCoord [lindex $tNode $i-1] 3];
		set bot2	[nodeCoord [lindex $bNode $i-1] 3];
		set dist [expr $top2-$bot2];

		# Calculate the building height as the running sum of the distances
		# between nodes specified
		set bdg_h [expr $bdg_h+$dist];

		# This means we take the distance in Z (3) in my coordinates systems at least. This is X-Y/Z| so X=1 Y=2 Z=3. (gli altri vacca gare)
		lappend h $dist;
		lappend mdrftX 0.0; # We will populate the lists with zeros initially
		lappend mdrftY 0.0;
		lappend maccelX 0.0;
		lappend maccelY 0.0;
		if {$dist==0} {puts "WARNING: Zerolength found in drift check"};


	}

	# Run the actual analysis now
	while {$cIndex==0 && $controlTime <= $Tmax && $ok==0} {
		# Runs while the building is stable, time is less
		# than that of the length of the record (plus buffering)
		# and the analysis is still converging

		# Do the analysis
		set ok [analyze 1 $Dt];		# Run a step of the analysis
		set controlTime [getTime];	# Update the control time
		if {$pflag>1} {puts "Completed  $controlTime of $Tmax seconds"}

		# If the analysis fails, try the following changes to achieve convergence
		# Analysis will be slower in here though...

		# First changes are to change algorithm to achieve convergence...
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Reduced timestep by half......"
			set Dtt [expr 0.5*$Dt]
			set ok [analyze 1 $Dtt]
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Reduced timestep by quarter......"
			set Dtt [expr 0.25*$Dt]
			set ok [analyze 1 $Dtt]
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying Broyden......"
			algorithm Broyden 8
			set ok [analyze 1 $Dt]
			algorithm 	$algorithmType
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying Newton with Initial Tangent ......"
			algorithm Newton -initial
			set ok [analyze 1 $Dt]
			algorithm 	$algorithmType
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying NewtonWithLineSearch......"
			algorithm NewtonLineSearch .8
			set ok [analyze 1 $Dt]
			algorithm 	$algorithmType
		}
		# Next change both algorithm and tolerance to achieve convergence....
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying Newton with Initial Tangent & relaxed convergence......"
			test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
			algorithm Newton -initial
			set ok [analyze 1 $Dt]
			test 		$testType $tolInit $iterInit
			algorithm 	$algorithmType
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying Newton with Initial Tangent & relaxed convergence......"
			test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
			algorithm Newton -initial
			set ok [analyze 1 $Dt]
			test 		$testType $tolInit $iterInit
			algorithm 	$algorithmType
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying NewtonWithLineSearch & relaxed convergence......"
			test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
			algorithm NewtonLineSearch .8
			set ok [analyze 1 $Dt]
			test 		$testType $tolInit $iterInit
			algorithm 	$algorithmType
		}
		# Next half the timestep with both algorithm and tolerance reduction, if this doesn't work - in bocca al lupo
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying Newton with Initial Tangent, reduced timestep & relaxed convergence......"
			test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
			algorithm Newton -initial
			set Dtt [expr 0.5*$Dt]
			set ok [analyze 1 $Dtt]
			test 		$testType $tolInit $iterInit
			algorithm 	$algorithmType
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying Newton with Initial Tangent, reduced timestep & relaxed convergence......"
			test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
			algorithm Newton -initial
			set Dtt [expr 0.5*$Dt]
			set ok [analyze 1 $Dtt]
			test 		$testType $tolInit $iterInit
			algorithm 	$algorithmType
		}
		if {$ok != 0} {
			puts " ~~~ Failed at $controlTime - Trying NewtonWithLineSearch, reduced timestep & relaxed convergence......"
			test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
			algorithm NewtonLineSearch .8
			set Dtt [expr 0.5*$Dt]
			set ok [analyze 1 $Dtt]
			test 		$testType $tolInit $iterInit
			algorithm 	$algorithmType
		}
		# Game over......
		if {$ok !=0} {
			puts " ~~~ Failed at $controlTime - exit analysis......"
			# Failed to converge, exit analysis
			# wipe;
			set cIndex -1;
		}

		# Check the actual state of the model with respect to the limits provided
		# Need to get the PGA (this is actually zero since nodeAcce returns relative not absolute values)
		set base_accelX 	[expr [nodeAccel [lindex $bNode 0] 1]/9.81]; 	# Current base node accel in X in g
		set base_accelY 	[expr [nodeAccel [lindex $bNode 0] 2]/9.81]; 	# Current base node accel in Y in g

		if {$base_accelX>[lindex $maccelX 0]} {set maccelX [lreplace $maccelX 0 0 $base_accelX]};
		if {$base_accelY>[lindex $maccelY 0]} {set maccelY [lreplace $maccelY 0 0 $base_accelY]};

		# Check the storey drifts and accelerations
		for {set i 1} {$i<=[llength $tNode]} {incr i 1} {
			set tNode_dispX 	[nodeDisp [lindex $tNode $i-1] 1]; 	# Current top node disp in X
			set tNode_dispY 	[nodeDisp [lindex $tNode $i-1] 2]; 	# Current top node disp in Y
			set bNode_dispX 	[nodeDisp [lindex $bNode $i-1] 1]; 	# Current bottom node disp in X
			set bNode_dispY 	[nodeDisp [lindex $bNode $i-1] 2]; 	# Current bottom node disp in Y
			set cHt				[lindex $h $i-1];					# Current storey height
			set cdrftX			[expr 100.0*abs($tNode_dispX-$bNode_dispX)/$cHt];	# Current storey drift in X at the current floor in %
			set cdrftY			[expr 100.0*abs($tNode_dispY-$bNode_dispY)/$cHt];	# Current storey drift in X at the current floor in %
			if {$cdrftX>=[lindex $mdrftX $i-1]} {set mdrftX [lreplace $mdrftX $i-1 $i-1 $cdrftX]};
			if {$cdrftY>=[lindex $mdrftY $i-1]} {set mdrftY [lreplace $mdrftY $i-1 $i-1 $cdrftY]};

	# 		set cdrft	[expr sqrt(pow($cdrftX,2)+pow($cdrftY,2))]; # square root sum of squares
			if {$cdrftX>=$cdrftY} {set cdrft $cdrftX} else {set cdrft $cdrftY}; # maximum of the two components
			if {$cdrft>$mdrft} {set mdrft $cdrft; set mflr $i}; # Update the current maximum storey drift and where it is

			# Now get the accelerations
			set Node_accelX 	[expr [nodeAccel [lindex $tNode $i-1] 1]/9.81]; 	# Current top node accel in X in g
			set Node_accelY 	[expr [nodeAccel [lindex $tNode $i-1] 2]/9.81]; 	# Current top node accel in Y in g
			set caccelX			[expr 1.0*abs($Node_accelX)];
			set caccelY			[expr 1.0*abs($Node_accelY)];
			if {$caccelX>=[lindex $maccelX $i]} {set maccelX [lreplace $maccelX $i $i $caccelX]};
			if {$caccelY>=[lindex $maccelY $i]} {set maccelY [lreplace $maccelY $i $i $caccelY]};
		}

		if {$mdrft>=$Dc} {set cIndex 1; set mdrft $Dc; wipe}; 			# Set the state of the model to local collapse (=1)

	}

	# Create some output
	puts $log [format "FinalState:%d at %.3f of %.3f seconds" $cIndex $controlTime $Tmax]; # Print to the logfile the final state
	puts $log [format "PeakStoreyDrift:%.2f%% at %d" $mdrft $mflr]; 	# Print to the max interstorey drift and where it is


	# Print to the max interstorey drifts
	puts -nonewline $log "PeakStoreyDriftX: ";
	for {set ii 1} {$ii <= [llength $mdrftX]} {incr ii 1} {
		puts -nonewline $log [format "%.2f " [lindex $mdrftX $ii-1]]
	}
	puts $log "%"; 	# Print to the max storey drifts
	puts -nonewline $log "PeakStoreyDriftY: ";
	for {set ii 1} {$ii <= [llength $mdrftY]} {incr ii 1} {
		puts -nonewline $log [format "%.2f " [lindex $mdrftY $ii-1]]
	}
	puts $log "%"; 	# Print to the max storey drifts


	# Print to the max max floor accelerations
	puts -nonewline $log "PeakFloorAccelerationX: ";
	for {set ii 1} {$ii <= [llength $maccelX]} {incr ii 1} {
		puts -nonewline $log [format "%.2f " [lindex $maccelX $ii-1]]
	}
	puts $log "g"; 	# Print to the max floor accelerations
	puts -nonewline $log "PeakFloorAccelerationY: ";
	for {set ii 1} {$ii <= [llength $maccelY]} {incr ii 1} {
		puts -nonewline $log [format "%.2f " [lindex $maccelY $ii-1]]
	}
	puts $log "g"; 	# Print to the max floor accelerations


	if {$pflag>0} {
		puts [format "PeakStoreyDrift:%.2f %% at %d" $mdrft $mflr]; 	# Print to the max roof drift
	}
	if {$cIndex == -1} {puts ":::::: ANALYSIS FAILED TO CONVERGE at $controlTime of $Tmax :::::"}
	if {$cIndex == 0} {puts  "######## ANALYSIS COMPLETED SUCCESSFULLY #####"}
	if {$cIndex == 1} { puts "========== LOCAL STRUCTURE COLLAPSE =========="}

}
