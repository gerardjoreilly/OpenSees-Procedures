# ---------------------------------------------------
# -- Script to Create a Lumped Hinge Element ---
# ---------------------------------------------------
# Gerard O'Reilly
# IUSS Pavia
# Created: January 2019
proc LumpedHingeElement {ET GT iNode jNode My1 My2 mu Lp fy fc b h ap app} {
# --------------------------------------
# -- Description of the Parameters -----
# --------------------------------------
# ET		# Element Tag
# GT		# Geometric Transf Tag
# iNode	# i Node
# jNode	# j Node
# My1		# Yield moment at end 1
# My2		# Yield moment at end 2
# mu		# Curvature ductility to peak
# Lp		# Plastic hinge length (m)
# fy		# Rebar yield strength (MPa)
# fc		# Concrete c strength (MPa)
# b		# Section width (m)
# h		# Section height (m)
# ap		# Peak strength ratio
# app		# Post-peak strength stiffness ratio

# --------------------------------------
# ------ Section Parameters ------------
# --------------------------------------
set epsy	[expr $fy/200e3]; # Assume E=200GPa for steel
set Ec	[expr (3320*sqrt($fc)+6900)*1.0e3];; # Assume for concrete (MPa)
set Ag	[expr $b*$h];
set Iz 	[expr $b*pow($h,3)/12];
set EIg	[expr $Ec*$Iz];
set phiy	[expr 2.1*$epsy/$h];
set EIcr1	[expr $My1/$phiy];
set EIcr2	[expr $My2/$phiy];
set cr1	[expr $EIcr1/$EIg];
set cr2	[expr $EIcr2/$EIg];
set Icr1	[expr $cr1*$Iz];
set Icr2	[expr $cr2*$Iz];

# --------------------------------------
# ------ Hysteretic Parameters ---------
# --------------------------------------
set pinchX 	0.8;
set pinchY 	0.5;
set damage1 0.0;
set damage2 0.0;
set beta 0.0;

# --------------------------------------
# -- Create the Material Model ---------
# --------------------------------------
set hingeMTag1 101${ET}; # Create an nonlinear material for the flexural hinge to be aggregated to the actual PH
set hingeMTag2 102${ET}; # Create an nonlinear material for the flexural hinge to be aggregated to the actual PH

set Mp1 [expr $ap*$My1];
set Mp2 [expr $ap*$My2];
set Mu1 [expr 0.1*$My1];
set Mu2 [expr 0.1*$My2];

set phip	[expr $phiy*$mu];
set phiu1	[expr $phip+($Mp1-$Mu1)/($app*$My1/$phiy)];
set phiu2	[expr $phip+($Mp2-$Mu2)/($app*$My2/$phiy)];

uniaxialMaterial Hysteretic $hingeMTag1 $My1 $phiy $Mp1 $phip $Mu1 $phiu1 -$My1 -$phiy -$Mp1 -$phip -$Mu1 -$phiu1 $pinchX $pinchY $damage1 $damage2 $beta
uniaxialMaterial Hysteretic $hingeMTag2 $My2 $phiy $Mp2 $phip $Mu2 $phiu2 -$My2 -$phiy -$Mp2 -$phip -$Mu2 -$phiu2 $pinchX $pinchY $damage1 $damage2 $beta
# --------------------------------------
# -- Create the Element  ---------------
# --------------------------------------
set intTag 		105${ET};	# Internal elastic section tag
set phTag1  	106${ET}; 	# Create an aggregated section with this tag for the actual element
set phTag2  	107${ET}; 	# Create an aggregated section with this tag for the actual element

# Create the internal elastic element behaviour
# section Elastic $intTag $Ec $Ag $Izz $Iyy $Gc 0.01
section Elastic $intTag $Ec $Ag $Iz $Iz [expr 0.4*$Ec] 0.01
puts "Cracked section to be fixed"

# Create the plastic hinge section
section Uniaxial $phTag1 $hingeMTag1 Mz; 							# Create the PH flexural section about ZZ
section Uniaxial $phTag2 $hingeMTag2 Mz; 							# Create the PH flexural section about ZZ

# Create the whole element
# element forceBeamColumn $eleTag $iNode $jNode $transfTag "HingeRadau $secTagI $LpI $secTagJ $LpJ $secTagInterior" <-mass $massDens> <-iter $maxIters $tol>
element forceBeamColumn $ET $iNode $jNode $GT "HingeRadau $phTag1 $Lp $phTag2 $Lp $intTag"
# element elasticBeamColumn $ET $iNode $jNode $Ag $Ec $Gc 0.01 $Iyy $Izz $GT

}
