# Procedure to carry out a pushover of a model
# Copyright by Gerard J. O'Reilly, 2017

# EUCENTRE/IUSSPavia
# Date Created: April 2012
# Last Updated: November 2015


proc singlePush {dref mu ctrlNode dispDir nSteps {IOflag 1} {PrintFlag 0}} {
# --------------------------------------------------
# Description of Parameters
# --------------------------------------------------
# dref:			Reference displacement to which cycles are run. Corresponds to yield or equivalent other, such as 1mm
# mu:			Multiple of dref to which the push is run. So pushover can be run to a specified ductility or displacement
# ctrlNode:		Node to control with the displacement integrator.
# dispDir:		DOF the loading is applied.
# nSteps:		Number of steps.
# IOflag:		Option to print details on screen. 2 for print of each step, 1 for basic info (default), 0 for off
# ---------------------------------------------------

# Set up the initial analysis
# Define the Initial Analysis Parameters
# set testType	NormUnbalance;				# Dont use with Penalty constraints
set testType	NormDispIncr
# set testType	EnergyIncr;					# Dont use with Penalty constraints
# set testType	RelativeNormUnbalance;		# Dont use with Penalty constraints
# set testType	RelativeNormDispIncr;		# Dont use with Lagrange constraints
# set testType	RelativeTotalNormDispIncr;	# Dont use with Lagrange constraints
# set testType	RelativeEnergyIncr;			# Dont use with Penalty constraints

set 	tolInit 	1.0e-7;			# Set the initial Tolerance, so it can be referred back to
set 	iterInit 	50;				# Set the initial Max Number of Iterations

set 	algorithmType KrylovNewton;		# Set the algorithm type
# set 	algorithmType Newton;		# Set the algorithm type
# set 	algorithmType Newton;		# Set the algorithm type

test 		$testType $tolInit $iterInit

algorithm 	$algorithmType
set 		disp 				[expr $dref*$mu];
set 		dU 					[expr $disp/(1.0*$nSteps)];
integrator 	DisplacementControl $ctrlNode $dispDir $dU
analysis 	Static

# Print values
if {$IOflag >= 1} {
	puts "singlePush: Push $ctrlNode to $mu"
}

# Set the initial values to start the while loop
set ok 		0;
set step 	1;
set loadf 	1.0;
# This feature of disabling the possibility of having a negative loading has been included.
# This has been adapted from a similar script by Prof. Garbaggio

while {$step<=$nSteps && $ok==0 && $loadf>0} {
	set ok 		[analyze 1];
	set loadf 	[getTime];
	set temp 	[nodeDisp $ctrlNode $dispDir];

	# Print the current displacement
	if {$IOflag >=2} {
		puts "Pushed $ctrlNode in $dispDir to $temp with $loadf"
	}

	# If the analysis fails, try the following changes to achieve convergence
	# Analysis will be slower in here though...
	if {$ok != 0} {
		puts "Trying relaxed convergence.."
		test 		$testType [expr $tolInit*0.01] [expr $iterInit*50]
		set ok [analyze 1]
		test 		$testType $tolInit $iterInit
	}
	if {$ok != 0} {
		puts "Trying Newton with initial then current .."
		test 		$testType [expr $tolInit*0.01] [expr $iterInit*50]
		algorithm Newton -initialThenCurrent
		set ok [analyze 1]
		algorithm 	$algorithmType
		test 		$testType $tolInit $iterInit
	}
	if {$ok != 0} {
		puts "Trying ModifiedNewton with initial .."
		test 		$testType [expr $tolInit*0.01] [expr $iterInit*50]
		algorithm ModifiedNewton -initial
		set ok [analyze 1]
		algorithm 	$algorithmType
		test 		$testType $tolInit $iterInit
	}
	if {$ok != 0} {
		puts "Trying KrylovNewton .."
		test 		$testType [expr $tolInit*0.01] [expr $iterInit*50]
		algorithm KrylovNewton
		set ok [analyze 1]
		algorithm 	$algorithmType
		test 		$testType $tolInit $iterInit
	}
	if {$ok != 0} {
		puts "Perform a Hail Mary ...."
		test 	FixedNumIter $iterInit
		set ok [analyze 1]
	}


	set temp 	[nodeDisp $ctrlNode $dispDir];
	set loadf 	[getTime];
	incr step 1;


}; # Close the while loop

if {$ok != 0} {
	puts "DispControl Analysis FAILED"
    #puts "Do you wish to continue y/n ?"; # include if want to pause at analysis failure
    #gets stdin ans; # not recommended in parameter study
    #if {$ans == "n"} done; # as it interrupts batch file
} else {
    puts "DispControl Analysis SUCCESSFUL"
}
if {$loadf<=0} {
    puts "Stopped because of Load factor below zero: $loadf"
}
if {$PrintFlag} {
	file delete singlePush.txt
	print singlePush.txt
}

}
