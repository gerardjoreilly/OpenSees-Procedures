# Procedure to carry out a dynamic analysis of a model for 1D models
# Gerard O'Reilly
# EUCENTRE/IUSSPavia
# Date Created: April 2012
# Last Updated: December 2016
# The coding is still pretty slap dash, but it does the job

proc runNRHA1D {Dt Tmax Dc aNode bNode log {pflag 0}} {
# --------------------------------------------------
# Description of Parameters
# --------------------------------------------------
# Dt:		Timestep of the analysis
# Tmax:		Length of the Record (Including Padding of 0's)
# Dc:		Displacement Capacity for both Interstorey and Roof Drift(rad)
# aNode:	Free Node
# bNode:	Fixed node
# log:		File handle of the logfile
# --------------------------------------------------
upvar 1 cIndex cIndex
if {$pflag>0} {
	puts "Starting runNRHA"
}
# Define the Initial Analysis Parameters
set 	testType	NormDispIncr;	# Set the test type
set 	tolInit 	1.0e-7;			# Set the initial Tolerance, so it can be referred back to
set 	iterInit 	50;				# Set the initial Max Number of Iterations
set 	algorithmType KrylovNewton;	# Set the algorithm type

test 		$testType $tolInit $iterInit
algorithm 	$algorithmType
integrator 	Newmark 0.5 0.25
analysis 	Transient

# Set up analysis parameters
set Nsteps [expr int($Tmax/$Dt)]; 	# integer of the number of steps needed
set cIndex 		0;					# Initially define the control index (-1 for non-converged, 0 for stable, 1 for global collapse, 2 for local collapse)
set controlTime 0.0;				# Start the controlTime
set ok			0;					# Set the convergence to 0 (converged)
set Dtt $Dt
# Set up initials
set mdef 0.0;

# Run the actual analysis now
while {$cIndex==0 && $controlTime <= $Tmax && $ok==0} {
	# Runs while the building is stable, time is less
	# than that of the length of the record (plus buffering)
	# and the analysis is still converging

	# Do the analysis
	set ok [analyze 1 $Dtt];		# Run a step of the analysis
	set controlTime [getTime];	# Update the control time
	if {$pflag>1} {
		puts "Completed  $controlTime of $Tmax seconds"
	}

	# Check the actual state of the model with respect to the limits provided

	# Check the drifts
	set aNode_def 	[nodeDisp $aNode 1]; 	# Current top node disp
	set bNode_def 	[nodeDisp $bNode 1]; 	# Current bottom node disp
	set cdef 		[expr abs($aNode_def-$bNode_def)];	# Current disp
	if {$cdef>=$mdef} {set mdef $cdef};
	if {$cdef>$mdef} {set mdef $cdef; set mflr $i};			# Update the current maximum

	if {$mdef>=$Dc} {set cIndex 1; set mdef $Dc; wipe}; 			# Set the state of the model to local collapse (=2)

	# If the analysis fails, try the following changes to achieve convergence
	# Analysis will be slower in here though...

	# First changes are to change algorithm to achieve convergence...
	# Next half the timestep with both algorithm and tolerance reduction, if this doesn't work - in bocca al lupo
	if {$ok != 0} {
		puts "Failed at $controlTime - Reduced timestep by half......"
		set Dtt [expr 0.5*$Dt]
		set ok [analyze 1 $Dtt]
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Reduced timestep by quarter......"
		set Dtt [expr 0.25*$Dtt]
		set ok [analyze 1 $Dtt]
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying Broyden......"
		algorithm Broyden 8
		set ok [analyze 1 $Dtt]
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying Newton with Initial Tangent ......"
		algorithm Newton -initial
		set ok [analyze 1 $Dtt]
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying NewtonWithLineSearch......"
		algorithm NewtonLineSearch .8
		set ok [analyze 1 $Dtt]
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	# Next change both algorithm and tolerance to achieve convergence....
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying Newton with Initial Tangent & relaxed convergence......"
		test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
		algorithm Newton -initial
		set ok [analyze 1 $Dtt]
		test 		$testType $tolInit $iterInit
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying Newton with Initial Tangent & relaxed convergence......"
		test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
		algorithm Newton -initial
		set ok [analyze 1 $Dtt]
		test 		$testType $tolInit $iterInit
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying NewtonWithLineSearch & relaxed convergence......"
		test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
		algorithm NewtonLineSearch .8
		set ok [analyze 1 $Dtt]
		test 		$testType $tolInit $iterInit
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	# Next half the timestep with both algorithm and tolerance reduction, if this doesn't work - in bocca al lupo
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying Newton with Initial Tangent, reduced timestep & relaxed convergence......"
		test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
		algorithm Newton -initial
		set Dtt [expr 0.5*$Dt]
		set ok [analyze 1 $Dtt]
		test 		$testType $tolInit $iterInit
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying Newton with Initial Tangent, reduced timestep & relaxed convergence......"
		test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
		algorithm Newton -initial
		set Dtt [expr 0.5*$Dt]
		set ok [analyze 1 $Dtt]
		test 		$testType $tolInit $iterInit
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	if {$ok != 0} {
		puts "Failed at $controlTime - Trying NewtonWithLineSearch, reduced timestep & relaxed convergence......"
		test NormDispIncr [expr $tolInit*0.1] [expr $iterInit*50]
		algorithm NewtonLineSearch .8
		set Dtt [expr 0.5*$Dt]
		set ok [analyze 1 $Dtt]
		test 		$testType $tolInit $iterInit
		algorithm 	$algorithmType
		set Dtt [expr $Dt]
	}
	# One last hail mary could be to try and introduce some numerical damping with the HHT
	if {$ok != 0} {
		puts "Failed at $controlTime - Introduce some numerical damping with HHT integrator......"
		integrator 	HHT 0.8
		set ok [analyze 1 $Dtt]
		integrator 	Newmark 0.5 0.25
		set Dtt [expr $Dt]
	}

	if {$ok !=0} {
		# Failed to converge, exit analysis
		wipe;
		set cIndex -1;
	}
}

# Create some output
puts $log [format "FinalState:%d at %.3f of %.3f seconds" $cIndex $controlTime $Tmax]; # Print to the logfile the final state
puts $log [format "Peak Deformation:%.4f rad" $mdef]; 	# Print to the max

if {$pflag>0} {
	puts [format "Peak Deformation: %.4f" $mdef]; 	# Print to the max
}
if {$cIndex == -1} {puts ":::::: ANALYSIS FAILED TO CONVERGE at $controlTime of $Tmax :::::"}
if {$cIndex == 0} {puts  "######## ANALYSIS COMPLETED SUCCESSFULLY #####"}
if {$cIndex == 1} { puts "========== LOCAL STRUCTURE COLLAPSE =========="}


}
