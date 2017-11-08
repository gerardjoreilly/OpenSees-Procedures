proc cyclicPush {dref mu numCycles ctrlNode dispDir dispIncr IOflag {PrintFlag 0}} {
# Procedure to carry out a cyclic pushover of a model
# Copyright by Gerard J. O'Reilly, 2017
# Created: March 2014
# --------------------------------------------------
# Description of Parameters
# --------------------------------------------------
# Command:		cyclicPush
# dref:		Reference displacement to which cycles are run. Corresponds to yield or equivalent other.
# mu:			Multiple of dref to which the cycles is run.
# numCycles:	No. of cycles. Valid options either 1,2,3,4,5,6
# ctrlNode:		Node to control with the displacement integrator
# dispDir:		Direction the loading is applied.
# dispIncr:		Number of displacement increments.
# IOflag:		Option to print cycle details on screen. 1 for on, 0 for off
# PrintFlag:	Optional flag for printing nodes/elements at max cycle
# ---------------------------------------------------
# cyclicPush dref mu numCycles ctrlNode dir dispIncr IOflag

# Set up the initial analysis
# Define the Initial Analysis Parameters
# set testType	NormUnbalance;			# Dont use with Penalty constraints
set testType	NormDispIncr
# set testType	EnergyIncr;				# Dont use with Penalty constraints
# set testType	RelativeNormUnbalance;		# Dont use with Penalty constraints
# set testType	RelativeNormDispIncr;		# Dont use with Lagrange constraints
# set testType	RelativeTotalNormDispIncr;	# Dont use with Lagrange constraints
# set testType	RelativeEnergyIncr;		# Dont use with Penalty constraints

set 	tolInit 	1.0e-7;			# Set the initial Tolerance, so it can be referred back to
set 	iterInit 	500;				# Set the initial Max Number of Iterations

set 	algorithmType KrylovNewton;		# Set the algorithm type
# set 	algorithmType Newton;		# Set the algorithm type
# set 	algorithmType Newton;		# Set the algorithm type

test 		$testType $tolInit $iterInit

algorithm 	$algorithmType


# Create the list of displacements
if {$numCycles == 1} {
	set dispList [list [expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]]
	set dispNoMax 3
} elseif {$numCycles == 2} {
	set dispList [list  [expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]]
	set dispNoMax 6
} elseif {$numCycles == 3} {
	set dispList [list  [expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]]
	set dispNoMax 9
} elseif {$numCycles == 4} {
	set dispList [list  [expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]]
	set dispNoMax 12
} elseif {$numCycles == 5} {
	set dispList [list  [expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]]
	set dispNoMax 15
} elseif {$numCycles == 6} {
	set dispList [list  [expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]\
						[expr $dref*$mu] [expr -2*$dref*$mu] [expr $dref*$mu]]
	set dispNoMax 18
} else {puts "ERROR: Value for numCycles not a valid choice. Choose between 1 and 6"}

# Print values
if {$IOflag == 1} {
	puts "cyclicPush: $numCycles cycles to $mu at $ctrlNode"
}

# Carry out loading
for {set d 1} {$d <= $dispNoMax} {incr d} {
	set numIncr $dispIncr
	set dU [expr [lindex $dispList [expr $d-1]]/(1.0*$numIncr)]
	integrator DisplacementControl $ctrlNode $dispDir $dU
	analysis Static
	for {set l 0} {$l < $numIncr} {incr l} {
		# puts "Analysis step: $l"
		set ok [analyze 1]
		if {$ok !=0} {
			puts "@@@@@@@@@@@@@@@@@@@@@@@@@"
			puts "@@@@ ANALYSIS FAILED @@@@"
			puts "@@@@@@@@@@@@@@@@@@@@@@@@@"
			puts "Analysis failed at cycle: $d and dispIncr: $l"
			exit;
		}
	}
	if {$PrintFlag ==1 && $d==1 } {
		print nodes_$mu.txt -node
		print ele_$mu.txt -ele
	}
}
}
