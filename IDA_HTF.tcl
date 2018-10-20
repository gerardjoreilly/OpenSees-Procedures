# --------------------------------------------------------------------------------------------------
# -- Script to Conduct Incremental Dynamic Analysis Using Hunt, Trace and Fill Algorithm -----------
# --------------------------------------------------------------------------------------------------
# Gerard O'Reilly
# EUCENTRE/IUSSPavia
# Created: November 2015
# Last Updated: May 2017

# --------------------------------------------------------------------------------------------------
# ----------------------------------------- Overview -----------------------------------------------
# --------------------------------------------------------------------------------------------------
# This is a script that will conduct an Incremental Dynamic Analysis (IDA) of a given structure
# using hunting, tracing and filling (HTF), where the user needs to provide a list of ground motion
# records and specify the increments, steps and maximum number of runs to conduct per record.
#
# The algorithm is inspired by that described in the Vamvatsikos & Cornell [2004] paper in
# Earthquake Spectra, but the main difference here is that the script is entirely TCL-based. This
# means that the procedure does not need Matlab to work, so the back and forth between Matlab and
# OpenSees during analysis is removed. The number of inputs is also reduced where just the initial
# intensity, the increment by which to increment the record and the maximum number of runs to be
# conducted per record are specified.
#
# The algorithm works by conducting an initial "hunting" phase, where the record is scaled up
# relatively quickly to find a collapse - hence, we are hunting out the collapse. This collapsed run
# is then discarded and this run is re-run at an intensity of the last non-collapsing run plus a
# fraction of the difference between the hunted collapse's intensity and the last non-collapsing
# intensity. This fraction is currently set at 0.25, but can be modified in the code. Once we go
# back to the last non-collapsing run and start slowly incrementing, we are in the "tracing" phase
# of the algorithm. Once collapse has been found from tracing, the remainder of the available runs
# are then used to go back and fill in the biggest gaps in the hunting phase's intensity steps,
# which is known as the "filling" phase.
#
# Some of the warning outputs:
# 	"----- WARNING: Collapsed achieved on first increment, reduce increment..."
# 	This means that the first step following the initial elastic run resuted in a collapse. In
# 	short, you have well over-shot the runway with the increment step. This is the second run,
# 	but if the building is collapsing on the very first run (i.e. the user defined intensity),
# 	well it seems you have bigger problems.
#
# 	"--- WARNING: First trace for collapse resulted in collapse..."
# 	Once the collapse has been hunted out, we go back to tracing the collapse more carefully.
# 	If this first trace results in collapse, either the last non-collapsing run was quite close
# 	to the collapsing intensity or the
# 	fraction of the difference is too big. Ideally, we would like to reduce this fraction
# 	automatically, but it's a bitch to code at the minute. Flagged for follow-up.....
#
# 	"----- WARNING: Collapse not achieved, increase increment or number of runs..."
# 	This means the collapse hasn't been hunted out yet, so either the incrementing to reach
# 	collapse is increased	or the maximum nuber of runs at the current increment is increased
# 	to allow it to go further. This warning	could be used as a way of gauging when collapse
# 	occurs, such that the max runs can be specified a priori such that there are enough runs left
# 	for tracing and filling.
#
# 	"----- WARNING: No filling, algorithm still tracing for collapse (reduce increment & increase runs)..."
# 	The algorithm is still tracing for collapse, either because there are not enough runs or
# 	because the increment used during the hunting was too big with respect to the fraction of
# 	the difference being used to trace.
#
# --------------------------------------------------------------------------------------------------
# Remaining issues
# --------------------------------------------------------------------------------------------------
#	Need to convert the names file into two column file for each direction. This current two
# 	file setup is a pain in planar analysis
#
# 	Handle collapse upon first trace. Try an exit command and print an error ti the file
#
# 	The algorithm hasnt been checked against how it copes with intermediate collapsing systems.
#
# 	The occurrence of non-converging has yet to be handled.
#
# 	Could implement a method, whereby 1 more intensity after first collapse is tried just to be
# 	sure
#
#  	The use of a single IM in 3D analysis is done by taking the IMX and IMY and getting the
# 	geomtric mean of this. Although the relative scale between the two components is maintained.
#
#	The periods are required to be know a priori which is a bit of a pain. An older version I
# 	had written of this function didnt require this because it loaded the model and did a modal
# 	analysis itself to get T1. But this can be problemnatic in the future when more complicated
# 	buildings are being analysed since OpenSees doesn't tell you explicitly which are the
# 	T1X T1Y T2X T2Y etc. so you might be mixing the periods from two different directions.
#
# 	There is bug that does a maxRun+1 run when collapse is reached at exactly maxRun number of runs
#
# --------------------------------------------------------------------------------------------------
# References
# --------------------------------------------------------------------------------------------------
# Vamvatsikos, D., and Cornell, C. A. [2004] “Applied Incremental Dynamic Analysis,”
# Earthquake Spectra, Vol. 20, No.2, pp. 523–553.
#
# Kazantzi, A. K., and Vamvatsikos, D. [2015a] “Intensity measure selection for vulnerability studies
# of building classes,” Earthquake Engineering & Structural Dynamics, Vol. 44, No.15, pp. 2677–2694.
#
# Kazantzi, A. K., and Vamvatsikos, D. [2015b] “A Next Generation Scalar Intensity Measure for
# Analytical Vulnerability Studies,” COMPDYN 2015 - 5th ECCOMAS Thematic Conference on Computational
# Methods in Structural Dynamics and Earthquake Engineering, Crete Island, Greece.
# --------------------------------------------------------------------------------------------------
# Inputs:
# --------------------------------------------------------------------------------------------------
# firstInt:		This is the first intensity to run the elastic run (e.g. 0.05g)
# incrStep:		This is the increment used during the hunting phase (e.g. 0.10g)
# maxRuns:		This is the maximum number of runs to use (e.g. 20)
# IMtype:		Intensity measure with which to conduct the IDA.
# 				1: 	PGA
# 				2: 	Sa at a given period (e.g. Sa(T1)), in 3D analysis this will take
# 					the geometric mean of the two record Sa(T) as the IM
# 				3:	From Kazantzi & Vamvatsikos [2015a], the intensity measure (d) which is
# 					given as the geometric mean (Sa,gm) at four periods defined as
# 					[T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m], where T1m and T2m are
# 					the mean first and second mode periods of the building class. Again
# 					the 3D IM will be the Sa,gm of the two components.
# 				4:	From Kazantzi & Vamvatsikos [2015a], the intensity measure (e) which is
# 					given as the geometric mean (Sa,gm) at five periods defined as
# 					[T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m, 2T1m], where T1m and T2m
# 					are the mean first and second mode periods of the building class. Again
# 					the 3D IM will be the Sa,gm of the two components.
# 				5:	From Kazantzi & Vamvatsikos [2015b], the intensity measure the same as
# 					4 above but with the inclusion of the duration via time between Arias
# 					intensity of 5% and 75% raised to the power 0.2. This gives the IM as
# 					Sa,gm*(t75%-t5%)^0.2
# Tinfo:		List of period info required by specified IM (e.g {1.0 2.0})
# 				1 - 	Dont need anything, just specify an empty list {}. It will ignore
# 					any entries if present
#				2 - 	Single value of period to condition to, like {1.0}
# 				3 - 	List of the two periods T1m and T2m, like {1.0 0.5}
# 				4 - 	List of the two periods T1m and T2m, like {1.0 0.5}
# 				5 - 	List of the two periods T1m and T2m, like {1.0 0.5}
# xi:			Elastic damping, typically 0.05
# Dt:			Analysis time step
# dCap:		Drift capacity in %
# nmsfileX:		Text file with the names of the X direction records in the form "*.txt"
# nmsfileY:		Text file with the names of the Y direction records in the form "*.txt"
# dtsfile:		Text file with the time steps of the records
# dursfile:		Text file with the durations of the records
# outsdir:		Where to print the outputs
# mdlfile:		OpenSees model (in 3D), the model must have the EQ pattern and response recorders
# 			already assigned inside model
# pflag:		Print flag (optional) put greater than 1 to print output on screen during analysis

# --------------------------------------------------------------------------------------------------
# Outputs:
# --------------------------------------------------------------------------------------------------
# The procedure will then print the intensities into the file:
# $outsdir/IM_Record${record number}.txt
#
# and also print a log file with the NRHA results into:
# $outsdir/log_IDA_Record${record number}_Run${run number}.txt

# --------------------------------------------------------------------------------------------------
# REQUIREMENTS:
# --------------------------------------------------------------------------------------------------
# To run this script, the following procedures are also required:
# runNRHA3D: This is to run the actual non-linear response history analysis

proc IDA_HTF {firstInt incrStep maxRuns IMtype Tinfo xi Dt dCap nmsfileX nmsfileY dtsfile dursfile outsdir mdlfile {pflag 0}} {
	# Create the output directory
	file mkdir $outsdir

	# Open an error file that will log the IDA_HTF errors
	set error_log [open $outsdir/IDA_HTF_error_log.txt "w"];
	puts "^^^^^^^^ STARTING IDA HTF ^^^^^^^^"
	puts $error_log "^^^^^^^^ STARTING IDA HTF ^^^^^^^^"

	# Get the ground motion set information
	set eqnms_listX 	[read [open $nmsfileX "r"]];
	set eqnms_listY 	[read [open $nmsfileY "r"]];
	set dts_list 	[read [open $dtsfile "r"]];
	set durs_list 	[read [open $dursfile "r"]];
	set nrecs 		[llength $dts_list];
	set g 9.81;

	for {set i 1} {$i <= $nrecs} {incr i 1} {
		set IM_log [open $outsdir/IM_${i}.txt "w"];

		# Load the info about the record
		set EQnameX [lindex $eqnms_listX $i-1]; 	# Get the name of the record1
		set EQnameY [lindex $eqnms_listY $i-1]; 	# Get the name of the record2
		set dt	[lindex $dts_list $i-1];	# Current dt
		set dur	[lindex $durs_list $i-1];	# Current duration

		# Establish the IM
		if {$IMtype==1} {
			# IM is PGA
			# Now get the spectral ordinates
			getSaT $EQnameX [lindex $dts_list $i-1] 0.0 $xi; # Get the PGA of the record in the X direction
			set IMX $pga
			getSaT $EQnameY [lindex $dts_list $i-1] 0.0 $xi; # Get the PGA of the record in the Y direction
			set IMY $pga

			# Now we have the IMX and IMY. In IDA we will use the geomen of these to get the
			# "current" IM. This way, the same scale factor will be applied to both
			set IMgeomean [expr pow($IMX*$IMY,0.5)];
		} elseif {$IMtype==2} {
			# IM is Sa at a specified period
			# Need to get the conditioning period
			set Tcond [lindex $Tinfo 0]; # It will be the first entry in the Tinfo list

			# Now get the spectral ordinates
			getSaT $EQnameX [lindex $dts_list $i-1] $Tcond $xi; # Get the Sa(T1,5%) of the record in the X direction
			set IMX $Sa
			getSaT $EQnameY [lindex $dts_list $i-1] $Tcond $xi; # Get the Sa(T1,5%) of the record in the Y direction
			set IMY $Sa

			# Now we have the IMX and IMY. In IDA we will use the geomen of these to get the
			# "current" IM. This way, the same scale factor will be applied to both
			set IMgeomean [expr pow($IMX*$IMY,0.5)];
		} elseif {$IMtype==3} {
			# IM is Sa,gm at a specified periods
			set T1m [lindex $Tinfo 0]; # It will be the first entry in the Tinfo list
			set T2m [lindex $Tinfo 1]; # It will be the second entry in the Tinfo list

			# Period list [T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m]
			set p_list {};
			lappend p_list $T2m;
			if {[expr 0.5*($T2m+$T1m)]<[expr 1.5*$T2m]} {
				lappend p_list [expr 0.5*($T2m+$T1m)];
			} else {
				lappend p_list [expr 1.5*$T2m];
			}
			lappend p_list $T1m;
			lappend p_list [expr 1.5*$T1m];

			# Get Spectral values at each
			set Sa_listX {};
			set Sa_listY {};
			for {set pn 1} {$pn<=[llength $p_list]} {incr pn 1} {
				getSaT $EQnameX [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the X direction
				lappend Sa_listX $Sa
				getSaT $EQnameY [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the Y direction
				lappend Sa_listY $Sa
			}

			# Get the geometric mean of these
			set IMX [expr pow([lindex $Sa_listX 0]*[lindex $Sa_listX 1]*[lindex $Sa_listX 2]*[lindex $Sa_listX 3],1/4.0)];
			set IMY [expr pow([lindex $Sa_listY 0]*[lindex $Sa_listY 1]*[lindex $Sa_listY 2]*[lindex $Sa_listY 3],1/4.0)];

			# Now we have the IMX and IMY. In IDA we will use the geomen of these to get the
			# "current" IM. This way, the same scale factor will be applied to both
			set IMgeomean [expr pow($IMX*$IMY,0.5)];
		} elseif {$IMtype==4} {
			# IM is Sa,gm at a specified periods
			set T1m [lindex $Tinfo 0]; # It will be the first entry in the Tinfo list
			set T2m [lindex $Tinfo 1]; # It will be the second entry in the Tinfo list

			# Period list [T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m, 2T1m]
			set p_list {};
			lappend p_list $T2m;
			if {[expr 0.5*($T2m+$T1m)]<[expr 1.5*$T2m]} {
				lappend p_list [expr 0.5*($T2m+$T1m)];
			} else {
				lappend p_list [expr 1.5*$T2m];
			}
			lappend p_list $T1m;
			lappend p_list [expr 1.5*$T1m];
			lappend p_list [expr 2.0*$T1m];

			# Get Spectral values at each
			set Sa_listX {};
			set Sa_listY {};
			for {set pn 1} {$pn<=[llength $p_list]} {incr pn 1} {
				getSaT $EQnameX [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the X direction
				lappend Sa_listX $Sa
				getSaT $EQnameY [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the Y direction
				lappend Sa_listY $Sa
			}

			# Get the geometric mean of these
			set IMX [expr pow([lindex $Sa_listX 0]*[lindex $Sa_listX 1]*[lindex $Sa_listX 2]*[lindex $Sa_listX 3]*[lindex $Sa_listX 4],1/5.0)];
			set IMY [expr pow([lindex $Sa_listY 0]*[lindex $Sa_listY 1]*[lindex $Sa_listY 2]*[lindex $Sa_listY 3]*[lindex $Sa_listY 4],1/5.0)];

			# Now we have the IMX and IMY. In IDA we will use the geomen of these to get the
			# "current" IM. This way, the same scale factor will be applied to both
			set IMgeomean [expr pow($IMX*$IMY,0.5)];
		} elseif {$IMtype==5} {
			# IM is Sa,gm at a specified periods
			set T1m [lindex $Tinfo 0]; # It will be the first entry in the Tinfo list
			set T2m [lindex $Tinfo 1]; # It will be the second entry in the Tinfo list

			# Period list [T2m, min[(T2m+T1m)/2, 1.5T2m], T1m, 1.5T1m, 2T1m]
			set p_list {};
			lappend p_list $T2m;
			if {[expr 0.5*($T2m+$T1m)]<[expr 1.5*$T2m]} {
				lappend p_list [expr 0.5*($T2m+$T1m)];
			} else {
				lappend p_list [expr 1.5*$T2m];
			}
			lappend p_list $T1m;
			lappend p_list [expr 1.5*$T1m];
			lappend p_list [expr 2.0*$T1m];

			# Get Spectral values at each
			set Sa_listX {};
			set Sa_listY {};
			for {set pn 1} {$pn<=[llength $p_list]} {incr pn 1} {
				getSaT $EQnameX [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the X direction
				lappend Sa_listX $Sa
				getSaT $EQnameY [lindex $dts_list $i-1] [lindex $p_list $pn-1] $xi; # Get the Sa(T1,5%) of the record in the Y direction
				lappend Sa_listY $Sa
			}

			# Get the geometric mean of these
			set SagmX [expr pow([lindex $Sa_listX 0]*[lindex $Sa_listX 1]*[lindex $Sa_listX 2]*[lindex $Sa_listX 3]*[lindex $Sa_listX 4],1/5.0)];
			set SagmY [expr pow([lindex $Sa_listY 0]*[lindex $Sa_listY 1]*[lindex $Sa_listY 2]*[lindex $Sa_listY 3]*[lindex $Sa_listY 4],1/5.0)];

			# Get the time between 5 and 75% of Arias Intensity
			Arias $EQnameX [lindex $dts_list $i-1] 0.05 0.75;
			set t12X $t12;
			Arias $EQnameY [lindex $dts_list $i-1] 0.05 0.75;
			set t12Y $t12;

			# Get the IMs in the two directions
			set IMX [expr $SagmX*pow($t12X,0.2)];
			set IMY [expr $SagmY*pow($t12Y,0.2)];

			# Now we have the IMX and IMY. In IDA we will use the geomen of these to get the
			# "current" IM. This way, the same scale factor will be applied to both
			set IMgeomean [expr pow($IMX*$IMY,0.5)];
		}

		# Set up the initial indices for HTF
		set j 	1;
		set IM	{};		# Initialise the list of IM used for printing
		set IMlist 	{};		# This is just a list that will be used in filling
		set hFlag	1;		# Hunting flag (1 for when we're hunting)
		set tFlag	0;		# Tracing flag (0 at first)
		set fFlag	0;		# Filling flag (0 at first)

		while {$j<=$maxRuns} {
			# As long as the hunting flag is 1, meaning we havent reached a collapse
			if {$hFlag==1} {
				# Determine the intensity to run at during the hunting (Andiamo a cacciare!)
				if {$j==1} {
					lappend IM $firstInt;
				} else {
					lappend IM [expr [lindex $IM $j-2]+($j-1)*$incrStep]; # Ramp it up!
				}
				# Determine the scale factor that needs to be applied to the record
				set sfX 	[expr [lindex $IM $j-1]/$IMgeomean*$g];
				set sfY 	[expr [lindex $IM $j-1]/$IMgeomean*$g];
				set 	run "Record${i}_Run${j}";			# This is a tag that outputs will be labelled with
				set 	log [open $outsdir/log_IDA_${run}.txt "w"];

				# The hunting intensity has been determined, now we can analyse
				source 	$mdlfile
				if {$pflag>0} {puts [format "Record:$i  Run:$j IM:%.3f" [lindex $IM $j-1] ]}
				runNRHA3D $Dt $dur $dCap $tNode $bNode $log $pflag
				close $log
				incr j 1;

				# Check the hunted run for collapse
				if {$cIndex>0} {
					set 	hFlag 	0;  # Stop hunting
					set 	tFlag 	1;  # Start tracing
					incr 	j 		-1; # Reduce by 1 because j was increased at end of hunting and we want to redo that point
					set 	jhunt	$j; # The value of j we hunted to
					if {$jhunt==2} {puts $error_log "WARNING: ${run} - Collapsed achieved on first increment, reduce increment..."};
				} else {
					puts $IM_log [format "%.3f" [lindex $IM $j-2]]; # j-2 because we've already increased j, but need to know if collapsed
				}
				wipe;
			}; # Close hunting

			# When the first collapse is reached, we start tracing between last convergence and the first collapse
			if {$tFlag==1} {
				# The first phase is to trace the last DeltaIM to get within the resolution
				if {$j==$jhunt} {
					set firstC 	[lindex $IM $j-1];				# This is the IM of the hunting collapse
					set IM 	[lreplace $IM $j-1 $j-1];		# Remove that value of IM from the array (it's already been appended)
				}
				set diff 	[expr $firstC-[lindex $IM $j-2]];	# Determine the difference between the hunting's noncollapse and collapse IM
				set inctr	[expr 0.20*$diff]; # Take 0.2 of the difference
				if {$inctr<0.05} {set inctr 0.025}; # Place a lower threshold on the increment so it doesnt start tracing too fine
				set IMtr 	[expr [lindex $IM $j-2]+$inctr];	# Calculate new tracing IM, which is previous noncollapse plus increment
				lappend IM $IMtr
				set sfX 	[expr [lindex $IM $j-1]/$IMgeomean*$g];
				set sfY 	[expr [lindex $IM $j-1]/$IMgeomean*$g];
				puts $IM_log  [format "%.3f" $IMtr];
				set 	run "Record${i}_Run${j}";			# This is a tag that outputs will be labelled with
				set 	log [open $outsdir/log_IDA_${run}.txt "w"];

				# The trace intensity has been determined, now we can analyse
				source 	$mdlfile
				if {$pflag>0} {puts [format "Record:$i  Run:$j IM:%.3f" $IMtr ]}
				runNRHA3D $Dt $dur $dCap $tNode $bNode $log $pflag

				close $log

				if {$cIndex>0} {
					# Not sure if this is the best way, to just trace back up to collapse again
					set tFlag 0; # Stop tracing
					set fFlag 1; # Start filling
					set jtrace	$j; # The value of j we traced to
					set IMlist $IM; # Get the list of IMs
					if {$j==$jhunt} {
						# This means the first trace collapsed, should reduce the increment
						puts $error_log "WARNING: ${run} - First trace for collapse resulted in collapse..."
					}
				}
				incr j 1;
				wipe;
			}; # Close the tracing

			# When the required resolution is reached, we start filling
			if {$fFlag==1} {
				# Reorder the list so we can account for filled runs
				set IMlist [lsort -real $IMlist];

				# Determine the biggest gap in IM for the hunted runs
				set gap 0.0;
				# We go to the end of the list minus 1 because, if not we would be filling between a noncollapsing and a collapsing run,
				# for which we are not sure if that filling run would be a non collapse - In short, does away with collapsing fills
				for {set ii 1} {$ii<[expr [llength $IMlist]-1]} {incr ii 1} {
					set temp [expr [lindex $IMlist $ii]-[lindex $IMlist $ii-1]]; 	# Find the running gap of hunted runs
					if {$temp>$gap} {
						set gap $temp
						set IMfil [expr [lindex $IMlist $ii-1]+$gap/2];		# Determine new filling IM
					};	# Update to maximum gap
				}

				lappend IM $IMfil
				lappend IMlist $IMfil
				set sfX 	[expr [lindex $IM $j-1]/$IMgeomean*$g];
				set sfY 	[expr [lindex $IM $j-1]/$IMgeomean*$g];
				puts $IM_log  [format "%.3f" $IMfil];
				set 	run "Record${i}_Run${j}";			# This is a tag that outputs will be labelled with
				set 	log [open $outsdir/log_IDA_${run}.txt "w"];

				# The trace intensity has been determined, now we can analyse
				source 	$mdlfile
				if {$pflag>0} {puts [format "Record:$i  Run:$j IM:%.3f" $IMfil ]}
				runNRHA3D $Dt $dur $dCap $tNode $bNode $log $pflag

				close $log
				incr j 1;
				wipe;
			}; # Close the filling

			# Wrap it up and finish
			if {$j==$maxRuns && $hFlag==1} {
				puts $error_log "WARNING: ${run} - Collapse not achieved, increase increment or number of runs..."
			};
			if {$j==$maxRuns && $fFlag==0} {
				puts $error_log "WARNING: ${run} - No filling, algorithm still tracing for collapse (reduce increment & increase runs)..."
			};
			wipe;
		}; # Close the maxRuns while loop
		close $IM_log
	}; # Close the ground motion loop

	puts "^^^^^^^^ FINISHED IDA HTF ^^^^^^^^"
	puts $error_log "^^^^^^^^ FINISHED IDA HTF ^^^^^^^^"
	close $error_log
};

# --------------------------------------------------------------------------------------------------
# -- Script to Compute Sa(T,xi%) -----
# --------------------------------------------------------------------------------------------------
#
# This is a script that will return the Sa(T1) of a given record,
# for a specified value of period T using the Newmark Average
# Acceleration. This is based on the Newmark.m matlab function.
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
	# Sv: Sv(T,%) - Pseudo-Spectral Velocity in m/s
	# Sd: Sd(T,%) - Spectral Displacement in m
	# PGA - Peak Ground Acceleration in g

	# Import the ground motion
	set imp [open $EQ "r"];
	set accg [read $imp];
	close $imp

	if {$T==0.0} {
		# Make the outputs available outside the procedure
		upvar 1 pga pga

		set pga 0.0;
		# Loop the terms (This loop is to length)
		for {set i 1} {$i <= [llength $accg]} {incr i 1} {
			set temp2 [expr abs([lindex $accg $i-1])]; 	# Get the absolute current ground acceleration
			if {$temp2>$pga} {set pga $temp2};			# Find the pga
		}
	} else {
		# Make the outputs available outside the procedure
		upvar 1 Sa Sa
		upvar 1 Sv Sv
		upvar 1 Sd Sd
		upvar 1 pga pga

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

		set umax 0.0;
		# Loop the terms (This loop is to length)
		for {set i 1} {$i <= [llength $t]} {incr i 1} {
			set temp1 [expr abs([lindex $u $i-1])]; 		# Get the absolute current displacement
			if {$temp1>$umax} {set umax $temp1};		# Find the peak displacement
		}

		set pga 0.0;
		# Loop the terms (This loop is to length)
		for {set i 1} {$i <= [llength $accg]} {incr i 1} {
			set temp2 [expr abs([lindex $accg $i-1])]; 	# Get the absolute current ground acceleration
			if {$temp2>$pga} {set pga $temp2};			# Find the pga
		}

		# Calculate Spectral Values
		set Sd 	$umax;			# Spectral acceleration in m
		set Sv	[expr $Sd*$w];		# Pseudo spectral velocity in m/s
		set Sa	[expr $Sd*$w*$w/$g];	# Pseudo spectral acceleration in g
	}

}

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
