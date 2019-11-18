# -------------------------------------------------------------------------------
# -- Script to Create a Lumped Plasticity Element for Ductile RC Beam Columns ---
# -------------------------------------------------------------------------------
# Gerard O'Reilly
# EUCENTRE/IUSSPavia
# Created: November 2014
# Last Updated: April 2016


# First create a function that performs Moment Curvature Analysis
proc MomentCurvature {index h b cv dbL dbV fc Ec P fyL Es rho1 rho2 rho3 {pflag 0}} {
# --------------------------------------------------------
# -- Set some stuff that makes My available above   ------
# --------------------------------------------------------
upvar Myp${index} Mp
upvar Myn${index} Mn
upvar cp${index} cp
upvar cn${index} cn
# --------------------------------------
# -- Compute Some General Stuff   ------
# --------------------------------------
set n_c		[expr 0.8+$fc/18];					# n term for concrete
set e_c		[expr 1.0*$fc/$Ec*($n_c/($n_c-1))]; 	# epsilon_c' for concrete
set e_s		[expr $fyL/$Es];					# Steel yield strain

# --------------------------------------
# -- Compute the Yield Curvature -------
# --------------------------------------
set phiY [expr 2.1*$fyL/$Es/$h]; # Yield Curvature (rad/m)

# --------------------------------------
# -- Compute the Yield Moment ----------
# --------------------------------------
# Do a Moment Curvature Analysis
set d1 		[expr $cv+$dbV+$dbL/2];		# Depth to Top Bars (m)
set d2 		[expr $h/2];				# Depth to Middle Bars (m)
set d3 		[expr $h-$cv-$dbV-$dbL/2];	# Depth to Bottom Bars (m)

# Positive Bending
set c 		[expr $h/2]; 				# initial trial of NA depth (m)
set count 	0;
set err		0.5;
while {$err > 0.001 && $count < 1000 } {
	if {$pflag>=2} {
		puts "Iteration $count c:$c  Error: $err kN";
	}
	# Compute the strains at each level
	set e_s1	[expr ($c-$d1)*$phiY]; 	# Strain in top steel (in strains)
	set e_s2	[expr ($d2-$c)*$phiY]; 	# Strain in middle steel
	set e_s3	[expr ($d3-$c)*$phiY]; 	# Strain in middle steel
	set e_top	[expr $c*$phiY];		# Strain in top of section

	# Compute the steel stuff
	if {$e_s1 < $e_s} {set f_s1 [expr $e_s1*$Es]} else {set f_s1 $fyL}; # 	in MPa
	if {$e_s2 < $e_s} {set f_s2 [expr $e_s2*$Es]} else {set f_s2 $fyL}; # 	in MPa
	if {$e_s3 < $e_s} {set f_s3 [expr $e_s3*$Es]} else {set f_s3 $fyL}; # 	in MPa

	set Fs1		[expr $f_s1*$rho1*$b*$d3*1000]; # (kN)
	set Fs2		[expr $f_s2*$rho2*$b*$d3*1000]; # (kN)
	set Fs3		[expr $f_s3*$rho3*$b*$d3*1000]; # (kN)

	# Compute concrete stuff
	set a1b1	[expr ($e_top/$e_c)-pow($e_top/$e_c,2)/3];		# alpha1beta1 term
	set b1		[expr (4-($e_top/$e_c))/(6-2*($e_top/$e_c))];	# beta1
	set Fc 		[expr $a1b1*$c*$fc*$b*1000];					# Concrete block force (kN)

	# Section force
	set Psec 	[expr $P+$Fs2+$Fs3-$Fc-$Fs1];					# Section Force (kN)

	# Adjust NA depth to balance section forces
	if {$Psec < 0} {
		set c0 $c
		set c [expr $c-0.001];
	} elseif {$Psec > 0} {
		set c0 $c
		set c [expr $c+0.001];
	}
	set err [expr abs($Psec)];

	if {$err < 5} {
	break
	}
	incr count 1
}
# Compute the moment
set Mp [expr $P*(0.5*$h-$c)+$Fs1*($c-$d1)+$Fs3*($d3-$c)+$Fs2*($d2-$c)+$Fc*$c*(1-$b1/2)];
set cp $c

# Negative Bending
set c 		[expr $h/2]; 				# initial trial of NA depth (m)
set count 	0;
set err		0.5;
while {$err > 0.001 && $count < 1000 } {
	if {$pflag>=2} {
		puts "Iteration $count c:$c  Error: $err kN";
	}
	# Compute the strains at each level
	set e_s1	[expr ($c-$d1)*$phiY]; 	# Strain in top steel (in strains)
	set e_s2	[expr ($d2-$c)*$phiY]; 	# Strain in middle steel
	set e_s3	[expr ($d3-$c)*$phiY]; 	# Strain in middle steel
	set e_top	[expr $c*$phiY];		# Strain in top of section

	# Compute the steel stuff
	if {$e_s1 < $e_s} {set f_s1 [expr $e_s1*$Es]} else {set f_s1 $fyL}; # 	in MPa
	if {$e_s2 < $e_s} {set f_s2 [expr $e_s2*$Es]} else {set f_s2 $fyL}; # 	in MPa
	if {$e_s3 < $e_s} {set f_s3 [expr $e_s3*$Es]} else {set f_s3 $fyL}; # 	in MPa

	set Fs1		[expr $f_s1*$rho3*$b*$d3*1000]; # (kN)
	set Fs2		[expr $f_s2*$rho2*$b*$d3*1000]; # (kN)
	set Fs3		[expr $f_s3*$rho1*$b*$d3*1000]; # (kN)

	# Compute concrete stuff
	set a1b1	[expr ($e_top/$e_c)-pow($e_top/$e_c,2)/3];		# alpha1beta1 term
	set b1		[expr (4-($e_top/$e_c))/(6-2*($e_top/$e_c))];	# beta1
	set Fc 		[expr $a1b1*$c*$fc*$b*1000];					# Concrete block force (kN)

	# Section force
	set Psec 	[expr $P+$Fs2+$Fs3-$Fc-$Fs1];					# Section Force (kN)

	# Adjust NA depth to balance section forces
	if {$Psec < 0} {
		set c0 $c
		set c [expr $c-0.001];
	} elseif {$Psec > 0} {
		set c0 $c
		set c [expr $c+0.001];
	}
	set err [expr abs($Psec)];

	if {$err < 5} {
	break
	}
	incr count 1
}
# Compute the moment
set Mn [expr $P*(0.5*$h-$c)+$Fs1*($c-$d1)+$Fs3*($d3-$c)+$Fs2*($d2-$c)+$Fc*$c*(1-$b1/2)];
set cn $c
# set Myp${index} $Mp
# set Myn${index} $Mn

if {$pflag>=1} {
	puts "Myp${index}: $Mp kNm"
	puts "Myn${index}: $Mn kNm"
}

}


# Ductile Element based on Hasleton Parameters
proc rcBC_duct {ET GT iNode jNode fyL fyV Es fc Ec b h s cv dbL dbV k asl P Ls rho_shr rho_top1zz rho_mid1zz rho_bot1zz rho_top2zz rho_mid2zz rho_bot2zz rho_top1yy rho_mid1yy rho_bot1yy rho_top2yy rho_mid2yy rho_bot2yy pfile {pflag 0}} {
global MPa
global mm
if {$pflag==1} {
	puts ""
	puts ""
	puts "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	puts "Element $ET between nodes $iNode and $jNode"
}
# --------------------------------------
# -- Description of the Parameters -----
# --------------------------------------
# ET		# Element Tag
# GT		# Geometric Transf Tag
# iNode	# i Node
# jNode	# j Node
# fyL 	# Steel yield strength long. bars (MPa)
# fyV 	# Steel yield strength shear. bars (MPa)
# Es 		# Steel elastic modulus (MPa)
# fc		# Concrete compressive strength (MPa)
# Ec		# Concrete elastic modulus (MPa)
# b		# Section width (m)
# h		# Section height (m)
# s		# Shear rft spacing
# cv		# Cover (m)
# dbL		# Long. bar diameter (m)
# dbV		# Shear bar diameter (m)
# k		# Confined strength ratio (1.0 for unconfined concrete) fcc'=k*fc
# asl		# Bond slip parameter (1 when possible, 0 when not)
# P			# Axial force (kN) (+ Compression, - Tension)
# Ls			# Shear span
# rho_shr		# Ratio of shear rft =Ash/bs
# rho_top1zz	# Ratio of top rft = Astop/bd (END 1 about zz axis)
# rho_mid1zz	# Ratio of web rft =Asmid/bd (END 1 about zz axis)
# rho_bot1zz	# Ratio of bottom rft = Asbot/bd (END 1 about zz axis)
# rho_top2zz	# Ratio of top rft = Astop/bd (END 2 about zz axis)
# rho_mid2zz	# Ratio of web rft =Asmid/bd (END 2 about zz axis)
# rho_bot2zz	# Ratio of bottom rft = Asbot/bd (END 2 about zz axis)
# rho_top1yy	# Ratio of top rft = Astop/bd (END 1 about yy axis)
# rho_mid1yy	# Ratio of web rft =Asmid/bd (END 1 about yy axis)
# rho_bot1yy	# Ratio of bottom rft = Asbot/bd (END 1 about yy axis)
# rho_top2yy	# Ratio of top rft = Astop/bd (END 2 about yy axis)
# rho_mid2yy	# Ratio of web rft =Asmid/bd (END 2 about yy axis)
# rho_bot2yy	# Ratio of bottom rft = Asbot/bd (END 2 about yy axis)
# pflag		# Print flag (optional - 1: details about hinge are printed)
# pfile		Print file -  prints the backbone properties to a specified file

# --------------------------------------
# -- Compute Some General Stuff   ------
# --------------------------------------
set fc 		[expr $fc*1.0]; 					# Change to real number incase it is integer
set nu 		[expr $P/($b*$h*$fc*$MPa)];			# Normalised Axial Load Ratio (-)
set dyy		[expr $h-$dbV-$cv-$dbL/2];			# Depth to bottom bars zz (m)
set dzz		[expr $b-$dbV-$cv-$dbL/2];			# Depth to bottom bars zz (m)
set d1		[expr $dbV-$cv-$dbL/2];	    			# Depth to top bars zz (m)
set Ag 		[expr $b*$h];					# Gross Cross Section Area (mm2)
set Izz 		[expr $b*pow($h,3)/12]; 			# I about local zz axis (mm4)
set Iyy 		[expr $h*pow($b,3)/12]; 			# I about local yy axis (mm4)
set EIzz		[expr $Ec*$MPa*$Izz];				# EI about zz (kNm2)
set EIyy		[expr $Ec*$MPa*$Iyy];				# EI about yy (kNm2)
set EA		[expr $Ec*$MPa*$Ag];				# EA (kN)
set Gc		[expr 0.4*$Ec*$MPa];				# Shear Modulus of Concrete (kN/m2)
set Kshear		[expr $Gc*$Ag];					# Shear Stiffness of Section

# Torsional Stiffness
if {$h>=$b} {set J [expr $h*pow($b,3)*(0.333-0.21*($h/$b)*(1-pow($b/$h,4)/12))];};
if {$b>$h} {set J [expr $b*pow($h,3)*(0.333-0.21*($b/$h)*(1-pow($h/$b,4)/12))];};

if {$pflag==1} {
	puts ""
	puts [format "b: %.3fm   h: %.3fm  N:%.1fkN" $b $h $P];
	puts [format "nu: %.3f" $nu];
	puts [format "Izz: %.6fm4   Iyy: %.6fm4" $Izz $Iyy];
	puts [format "EIyy: %.1fkNm2   EIzz: %.1fkNm2"  $EIyy $EIzz];
}

# --------------------------------------
# -- Compute the Yield Curvature -------
# --------------------------------------
# Assume square section and use PCK2007 expression
set phiYzz [expr 2.1*$fyL/$Es/$h]; # Yield Curvature (rad/m)
set phiYyy [expr 2.1*$fyL/$Es/$b]; # Yield Curvature (rad/m)

if {$pflag==1} {
	puts ""
	puts [format "phiYzz: %.4f rad/m      phiYyy: %.4f rad/m" $phiYzz $phiYyy];

}


# --------------------------------------
# -- Compute the Nominal Moment --------
# --------------------------------------
MomentCurvature 1zz $h $b $cv $dbL $dbV [expr $k*$fc] $Ec $P $fyL $Es $rho_top1zz $rho_mid1zz $rho_bot1zz
MomentCurvature 2zz $h $b $cv $dbL $dbV [expr $k*$fc] $Ec $P $fyL $Es $rho_top2zz $rho_mid2zz $rho_bot2zz
MomentCurvature 1yy $b $h $cv $dbL $dbV [expr $k*$fc] $Ec $P $fyL $Es $rho_top1yy $rho_mid1yy $rho_bot1yy
MomentCurvature 2yy $b $h $cv $dbL $dbV [expr $k*$fc] $Ec $P $fyL $Es $rho_top2yy $rho_mid2yy $rho_bot2yy

set Kizz		[expr $Myp1zz/$phiYzz];		# Initial cracked section stiffness about zz
set Kiyy		[expr $Myp1yy/$phiYyy];		# Initial cracked section stiffness about yy
set EIrzz 		[expr $Kizz/$EIzz];		# Ratio of gross EI to get cracked section EI
set EIryy 		[expr $Kiyy/$EIyy];		# Ratio of gross EI to get cracked section EI
set EIzze		[expr $EIrzz*$EIzz];		# Cracked EI about zz
set EIyye		[expr $EIryy*$EIyy];		# Cracked EI about yy
set Izze		[expr $EIrzz*$Izz];		# Cracked I about zz
set Iyye		[expr $EIryy*$Iyy];		# Cracked I about yy

if {$pflag==1} {
	puts ""
	puts [format "EIrzz: %.3f   EIryy: %.3f"  $EIrzz $EIryy];
	puts ""
	puts [format "Myp1zz: %.1fkNm  Myp2zz: %.1fkNm     Myp1yy: %.1fkNm  Myp2yy: %.1fkNm" $Myp1zz $Myp2zz $Myp1yy $Myp2yy];
	puts [format "Myn1zz: %.1fkNm  Myn2zz: %.1fkNm     Myn1yy: %.1fkNm  Myn2yy: %.1fkNm" $Myn1zz $Myn2zz $Myn1yy $Myn2yy];
}

# --------------------------------------
# -- Compute the Capping Moment --------
# --------------------------------------
set Mcp1zz 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myp1zz]; # p59
set Mcp2zz 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myp2zz];
set Mcp1yy 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myp1yy];
set Mcp2yy 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myp2yy];

set Mcn1zz 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myn1zz]; # p59
set Mcn2zz 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myn2zz];
set Mcn1yy 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myn1yy];
set Mcn2yy 		[expr 1.25*pow(0.89,$nu)*pow(0.91,0.01*$fc)*$Myn2yy];

if {$pflag==1} {
	puts ""
	puts [format "Mcp1zz: %.1fkNm  Mcp2zz: %.1fkNm     Mcp1yy: %.1fkNm  Mcp2yy: %.1fkNm" $Mcp1zz $Mcp2zz $Mcp1yy $Mcp2yy];
	puts [format "Mcn1zz: %.1fkNm  Mcn2zz: %.1fkNm     Mcn1yy: %.1fkNm  Mcn2yy: %.1fkNm" $Mcn1zz $Mcn2zz $Mcn1yy $Mcn2yy];
}

# --------------------------------------
# -- Compute the Ultimate Moment -------
# --------------------------------------
set Mup1zz [expr 0.1*$Myp1zz];
set Mup2zz [expr 0.1*$Myp2zz];
set Mup1yy [expr 0.1*$Myp1yy];
set Mup2yy [expr 0.1*$Myp2yy];

set Mun1zz [expr 0.1*$Myn1zz];
set Mun2zz [expr 0.1*$Myn2zz];
set Mun1yy [expr 0.1*$Myn1yy];
set Mun2yy [expr 0.1*$Myn2yy];

if {$pflag==1} {
	puts ""
	puts [format "Mup1zz: %.1fkNm  Mup2zz: %.1fkNm     Mup1yy: %.1fkNm  Mup2yy: %.1fkNm" $Mup1zz $Mup2zz $Mup1yy $Mup2yy];
	puts [format "Mun1zz: %.1fkNm  Mun2zz: %.1fkNm     Mun1yy: %.1fkNm  Mun2yy: %.1fkNm" $Mun1zz $Mun2zz $Mun1yy $Mun2yy];
}

# --------------------------------------
# -- Compute the Plastic Hinge Length --
# --------------------------------------
set Lp [expr 0.08*$Ls+0.022*$fyL*$dbL]; # Priestley + Park (1992)

if {$pflag==1} {
	puts ""
	puts [format "Lp: %.4fm" $Lp];
}

# --------------------------------------
# -- Compute the Ultimate Curvature ----
# --------------------------------------
set sn 		[expr $s/$dbL*pow(0.01*$fyL,0.5)]; # p xv

set thetaCapPl1zz [expr 0.12*(1+$asl*0.55)*pow(0.16,$nu)*pow(0.02+40*$rho_shr,0.43)*pow(0.54,0.01*$fc)*pow(0.66,0.1*$sn)*pow(2.27,10*($rho_top1zz+$rho_mid1zz+$rho_bot1zz))]; # p40
set thetaCapPl2zz [expr 0.12*(1+$asl*0.55)*pow(0.16,$nu)*pow(0.02+40*$rho_shr,0.43)*pow(0.54,0.01*$fc)*pow(0.66,0.1*$sn)*pow(2.27,10*($rho_top2zz+$rho_mid2zz+$rho_bot2zz))];
set thetaCapPl1yy [expr 0.12*(1+$asl*0.55)*pow(0.16,$nu)*pow(0.02+40*$rho_shr,0.43)*pow(0.54,0.01*$fc)*pow(0.66,0.1*$sn)*pow(2.27,10*($rho_top1yy+$rho_mid1yy+$rho_bot1yy))];
set thetaCapPl2yy [expr 0.12*(1+$asl*0.55)*pow(0.16,$nu)*pow(0.02+40*$rho_shr,0.43)*pow(0.54,0.01*$fc)*pow(0.66,0.1*$sn)*pow(2.27,10*($rho_top2yy+$rho_mid2yy+$rho_bot2yy))];

set phiC1zz 	[expr $thetaCapPl1zz/$Lp+$phiYzz];
set phiC2zz 	[expr $thetaCapPl2zz/$Lp+$phiYzz];
set phiC1yy 	[expr $thetaCapPl1yy/$Lp+$phiYyy];
set phiC2yy 	[expr $thetaCapPl2yy/$Lp+$phiYyy];

if {$pflag==1} {
	puts ""
	puts [format "phiC1zz: %.3frad/m  phiC2zz: %.3frad/m     phiC1yy: %.3frad/m  phiC2yy: %.3frad/m" $phiC1zz $phiC2zz $phiC1yy $phiC2yy];
}

# --------------------------------------
# -- Compute the Capping Curvature -----
# --------------------------------------
set thetaPc [expr 0.76*pow(0.031,$nu)*pow(0.02+40*$rho_shr,1.02)]; # p54
if {$thetaPc > 0.10 } {set thetaPc 0.1}

set phiU1zz	[expr $phiC1zz+$thetaPc/$Lp];
set phiU2zz	[expr $phiC2zz+$thetaPc/$Lp];
set phiU1yy	[expr $phiC1yy+$thetaPc/$Lp];
set phiU2yy	[expr $phiC2yy+$thetaPc/$Lp];

if {$pflag==1} {
	puts ""
	puts [format "phiU1zz: %.3frad/m  phiU2zz: %.3frad/m     phiU1yy: %.3frad/m  phiU2yy: %.3frad/m" $phiU1zz $phiU2zz $phiU1yy $phiU2yy];
}

# --------------------------------------
# -- Compute Stiffness Deterioration ---
# --------------------------------------
set Lambdazz [expr 170.7*pow(0.27,$nu)*pow(0.1,$s/$dzz)]; # p65
set Lambdayy [expr 170.7*pow(0.27,$nu)*pow(0.1,$s/$dyy)];
# also see comments on p9
if {$pflag==1} {
	puts ""
	puts [format "Lambdazz: %.3f  Lambdayy: %.3f" $Lambdazz $Lambdayy];
}

# --------------------------------------
# -- Compute Hardening Ratio  ----------
# --------------------------------------
set Ksp1zz [expr ($Mcp1zz-$Myp1zz)/($phiC1zz-$phiYzz)/$Myp1zz*$phiYzz];
set Ksp2zz [expr ($Mcp2zz-$Myp2zz)/($phiC2zz-$phiYzz)/$Myp2zz*$phiYzz];
set Ksp1yy [expr ($Mcp1yy-$Myp1yy)/($phiC1yy-$phiYyy)/$Myp1yy*$phiYyy];
set Ksp2yy [expr ($Mcp2yy-$Myp2yy)/($phiC2yy-$phiYyy)/$Myp2yy*$phiYyy];

set Ksn1zz [expr ($Mcn1zz-$Myn1zz)/($phiC1zz-$phiYzz)/$Myn1zz*$phiYzz];
set Ksn2zz [expr ($Mcn2zz-$Myn2zz)/($phiC2zz-$phiYzz)/$Myn2zz*$phiYzz];
set Ksn1yy [expr ($Mcn1yy-$Myn1yy)/($phiC1yy-$phiYyy)/$Myn1yy*$phiYyy];
set Ksn2yy [expr ($Mcn2yy-$Myn2yy)/($phiC2yy-$phiYyy)/$Myn2yy*$phiYyy];

if {$pflag==1} {
	puts ""
	puts [format "Ksp1zz: %.3f  Ksp2zz: %.3f     Ksp1yy: %.3f  Ksp2yy: %.3f" $Ksp1zz $Ksp2zz $Ksp1yy $Ksp2yy];
	puts [format "Ksn1zz: %.3f  Ksn2zz: %.3f     Ksn1yy: %.3f  Ksn2yy: %.3f" $Ksn1zz $Ksn2zz $Ksn1yy $Ksn2yy];
}

# --------------------------------------
# -- Print Some Output  ----------------
# --------------------------------------
puts $pfile [format "Element %d between nodes %d and %d. Myp1zz: %.1f kNm Myp2zz: %.1f kNm Myp1yy: %.1f kNm Myp2yy: %.1f kNm Myn1zz: %.1f kNm Myn2zz: %.1f kNm Myn1yy: %.1f kNm Myn2yy: %.1f kNm Mcp1zz: %.1f kNm Mcp2zz: %.1f kNm Mcp1yy: %.1f kNm Mcp2yy: %.1f kNm Mcn1zz: %.1f kNm Mcn2zz: %.1f kNm Mcn1yy: %.1f kNm Mcn2yy: %.1f kNm Mup1zz: %.1f kNm Mup2zz: %.1f kNm Mup1yy: %.1f kNm Mup2yy: %.1f kNm Mun1zz: %.1f kNm Mun2zz: %.1f kNm Mun1yy: %.1f kNm Mun2yy: %.1f kNm phiYzz: %.4f rad/m phiYyy: %.4f rad/m phiC1zz: %.3f rad/m phiC2zz: %.3f rad/m phiC1yy: %.3f rad/m phiC2yy: %.3f rad/m phiU1zz: %.3f rad/m phiU2zz: %.3f rad/m phiU1yy: %.3f rad/m phiU2yy: %.3f rad/m thetaCapPl1zz: %.3f rad thetaCapPl2zz: %.3f rad thetaCapPl1yy: %.3f rad thetaCapPl2yy: %.3f rad thetaPc: %.3f rad Lp: %.4f m Lambdazz: %.3f Lambdayy: %.3f nu: %.3f EIrzz: %.3f EIryy: %.3f" $ET $iNode $jNode $Myp1zz $Myp2zz $Myp1yy $Myp2yy $Myn1zz $Myn2zz $Myn1yy $Myn2yy $Mcp1zz $Mcp2zz $Mcp1yy $Mcp2yy $Mcn1zz $Mcn2zz $Mcn1yy $Mcn2yy $Mup1zz $Mup2zz $Mup1yy $Mup2yy $Mun1zz $Mun2zz $Mun1yy $Mun2yy $phiYzz $phiYyy $phiC1zz $phiC2zz $phiC1yy $phiC2yy $phiU1zz $phiU2zz $phiU1yy $phiU2yy $thetaCapPl1zz $thetaCapPl2zz $thetaCapPl1yy $thetaCapPl2yy $thetaPc $Lp $Lambdazz $Lambdayy $nu $EIrzz $EIryy];

# --------------------------------------
# -- Create the Material Model ---------
# --------------------------------------
set hingeMTag1zz 101${ET}; # Create an nonlinear material for the flexural hinge to be aggregated to the actual PH
set hingeMTag2zz 102${ET}; # Create an nonlinear material for the flexural hinge to be aggregated to the actual PH
set hingeMTag1yy 103${ET}; # Create an nonlinear material for the flexural hinge to be aggregated to the actual PH
set hingeMTag2yy 104${ET}; # Create an nonlinear material for the flexural hinge to be aggregated to the actual PH

# Material Model
set K0 			[expr $Myp1zz/$phiYzz]; 	# elastic stiffness
set as_Plus 		$Ksp1zz; 			# strain hardening ratio for positive loading direction
set as_Neg 			$Ksp1zz; 			# strain hardening ratio for negative loading direction
set My_Plus 		$Myp1zz; 			# effective yield strength for positive loading direction
set My_Neg 			-$Myn1zz; 			# effective yield strength for negative loading direction (negative value)
set Lamda_S 		$Lambdazz; 			# Cyclic deterioration parameter for strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_C 		$Lambdazz; 			# Cyclic deterioration parameter for post-capping strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_A 		0.0; 				# Cyclic deterioration parameter for acceleration reloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_K 		0.0; 				# Cyclic deterioration parameter for unloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set c_S 			1.0; 				# rate of strength deterioration. The default value is 1.0.
set c_C 			1.0; 				# rate of post-capping strength deterioration. The default value is 1.0.
set c_A 			0.0; 				# rate of accelerated reloading deterioration. The default value is 1.0.
set c_K 			0.0; 				# rate of unloading stiffness deterioration. The default value is 1.0.
set phi_p_Plus		[expr $phiC1zz-$phiYzz]; 	# pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
set phi_p_Neg		[expr $phiC1zz-$phiYzz]; 	# pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (must be defined as a positive value)
set phi_pc_Plus		[expr $phiU1zz-$phiC1zz]; 	# post-capping rotation for positive loading direction
set phi_pc_Neg		[expr $phiU1zz-$phiC1zz]; 	# post-capping rotation for negative loading direction (must be defined as a positive value)
set Res_Pos			[expr $Mup1zz/$Myp1zz]; 	# residual strength ratio for positive loading direction
set Res_Neg			[expr $Mun1zz/$Myn1zz];  	# residual strength ratio for negative loading direction (must be defined as a positive value)
set phi_u_Plus		[expr 2*$phiU1zz]; 		# ultimate rotation capacity for positive loading direction
set phi_u_Neg		[expr 2*$phiU1zz]; 		# ultimate rotation capacity for negative loading direction (must be defined as a positive value)
set D_Plus			1.0; 				# rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
set D_Neg			1.0; 				# rate of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
uniaxialMaterial ModIMKPeakOriented $hingeMTag1zz $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $phi_p_Plus $phi_p_Neg $phi_pc_Plus $phi_pc_Neg $Res_Pos $Res_Neg $phi_u_Plus $phi_u_Neg $D_Plus $D_Neg

set K0 			[expr $Myp2zz/$phiYzz]; 	# elastic stiffness
set as_Plus 		$Ksp2zz; 			# strain hardening ratio for positive loading direction
set as_Neg 			$Ksp2zz; 			# strain hardening ratio for negative loading direction
set My_Plus 		$Myp2zz; 			# effective yield strength for positive loading direction
set My_Neg 			-$Myn2zz; 			# effective yield strength for negative loading direction (negative value)
set Lamda_S 		$Lambdazz; 			# Cyclic deterioration parameter for strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_C 		$Lambdazz; 			# Cyclic deterioration parameter for post-capping strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_A 		0.0; 				# Cyclic deterioration parameter for acceleration reloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_K 		0.0; 				# Cyclic deterioration parameter for unloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set c_S 			1.0; 				# rate of strength deterioration. The default value is 1.0.
set c_C 			1.0; 				# rate of post-capping strength deterioration. The default value is 1.0.
set c_A 			0.0; 				# rate of accelerated reloading deterioration. The default value is 1.0.
set c_K 			0.0; 				# rate of unloading stiffness deterioration. The default value is 1.0.
set phi_p_Plus		[expr $phiC2zz-$phiYzz]; 	# pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
set phi_p_Neg		[expr $phiC2zz-$phiYzz]; 	# pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (must be defined as a positive value)
set phi_pc_Plus		[expr $phiU2zz-$phiC2zz]; 	# post-capping rotation for positive loading direction
set phi_pc_Neg		[expr $phiU2zz-$phiC2zz]; 	# post-capping rotation for negative loading direction (must be defined as a positive value)
set Res_Pos			[expr $Mup2zz/$Myp2zz]; 	# residual strength ratio for positive loading direction
set Res_Neg			[expr $Mun2zz/$Myn2zz];  	# residual strength ratio for negative loading direction (must be defined as a positive value)
set phi_u_Plus		[expr 2*$phiU2zz]; 		# ultimate rotation capacity for positive loading direction
set phi_u_Neg		[expr 2*$phiU2zz]; 		# ultimate rotation capacity for negative loading direction (must be defined as a positive value)
set D_Plus			1.0; 				# rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
set D_Neg			1.0; 				# rate of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
uniaxialMaterial ModIMKPeakOriented $hingeMTag2zz $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $phi_p_Plus $phi_p_Neg $phi_pc_Plus $phi_pc_Neg $Res_Pos $Res_Neg $phi_u_Plus $phi_u_Neg $D_Plus $D_Neg

set K0 			[expr $Myp1yy/$phiYyy]; 	# elastic stiffness
set as_Plus 		$Ksp1yy; 			# strain hardening ratio for positive loading direction
set as_Neg 			$Ksp1yy; 			# strain hardening ratio for negative loading direction
set My_Plus 		$Myp1yy; 			# effective yield strength for positive loading direction
set My_Neg 			-$Myn1yy; 			# effective yield strength for negative loading direction (negative value)
set Lamda_S 		$Lambdayy; 			# Cyclic deterioration parameter for strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_C 		$Lambdayy; 			# Cyclic deterioration parameter for post-capping strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_A 		0.0; 				# Cyclic deterioration parameter for acceleration reloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_K 		0.0; 				# Cyclic deterioration parameter for unloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set c_S 			1.0; 				# rate of strength deterioration. The default value is 1.0.
set c_C 			1.0; 				# rate of post-capping strength deterioration. The default value is 1.0.
set c_A 			0.0; 				# rate of accelerated reloading deterioration. The default value is 1.0.
set c_K 			0.0; 				# rate of unloading stiffness deterioration. The default value is 1.0.
set phi_p_Plus		[expr $phiC1yy-$phiYyy]; 	# pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
set phi_p_Neg		[expr $phiC1yy-$phiYyy]; 	# pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (must be defined as a positive value)
set phi_pc_Plus		[expr $phiU1yy-$phiC1yy]; 	# post-capping rotation for positive loading direction
set phi_pc_Neg		[expr $phiU1yy-$phiC1yy]; 	# post-capping rotation for negative loading direction (must be defined as a positive value)
set Res_Pos			[expr $Mup1yy/$Myp1yy]; 	# residual strength ratio for positive loading direction
set Res_Neg			[expr $Mun1yy/$Myn1yy];  	# residual strength ratio for negative loading direction (must be defined as a positive value)
set phi_u_Plus		[expr 2*$phiU1yy]; 		# ultimate rotation capacity for positive loading direction
set phi_u_Neg		[expr 2*$phiU1yy]; 		# ultimate rotation capacity for negative loading direction (must be defined as a positive value)
set D_Plus			1.0; 				# rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
set D_Neg			1.0; 				# rate of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
uniaxialMaterial ModIMKPeakOriented $hingeMTag1yy $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $phi_p_Plus $phi_p_Neg $phi_pc_Plus $phi_pc_Neg $Res_Pos $Res_Neg $phi_u_Plus $phi_u_Neg $D_Plus $D_Neg

set K0 			[expr $Myp2yy/$phiYyy]; 	# elastic stiffness
set as_Plus 		$Ksp2yy; 			# strain hardening ratio for positive loading direction
set as_Neg 			$Ksp2yy; 			# strain hardening ratio for negative loading direction
set My_Plus 		$Myp2yy; 			# effective yield strength for positive loading direction
set My_Neg 			-$Myn2yy; 			# effective yield strength for negative loading direction (negative value)
set Lamda_S 		$Lambdayy; 			# Cyclic deterioration parameter for strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_C 		$Lambdayy; 			# Cyclic deterioration parameter for post-capping strength deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_A 		0.0; 				# Cyclic deterioration parameter for acceleration reloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set Lamda_K 		0.0; 				# Cyclic deterioration parameter for unloading stiffness deterioration [see definitions in Lignos and Krawinkler (2011)]
set c_S 			1.0; 				# rate of strength deterioration. The default value is 1.0.
set c_C 			1.0; 				# rate of post-capping strength deterioration. The default value is 1.0.
set c_A 			0.0; 				# rate of accelerated reloading deterioration. The default value is 1.0.
set c_K 			0.0; 				# rate of unloading stiffness deterioration. The default value is 1.0.
set phi_p_Plus		[expr $phiC2yy-$phiYyy]; 	# pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
set phi_p_Neg		[expr $phiC2yy-$phiYyy]; 	# pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (must be defined as a positive value)
set phi_pc_Plus		[expr $phiU2yy-$phiC2yy]; 	# post-capping rotation for positive loading direction
set phi_pc_Neg		[expr $phiU2yy-$phiC2yy]; 	# post-capping rotation for negative loading direction (must be defined as a positive value)
set Res_Pos			[expr $Mup2yy/$Myp2yy]; 	# residual strength ratio for positive loading direction
set Res_Neg			[expr $Mun2yy/$Myn2yy];  	# residual strength ratio for negative loading direction (must be defined as a positive value)
set phi_u_Plus		[expr 2*$phiU2yy]; 		# ultimate rotation capacity for positive loading direction
set phi_u_Neg		[expr 2*$phiU2yy]; 		# ultimate rotation capacity for negative loading direction (must be defined as a positive value)
set D_Plus			1.0; 				# rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
set D_Neg			1.0; 				# rate of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
uniaxialMaterial ModIMKPeakOriented $hingeMTag2yy $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $phi_p_Plus $phi_p_Neg $phi_pc_Plus $phi_pc_Neg $Res_Pos $Res_Neg $phi_u_Plus $phi_u_Neg $D_Plus $D_Neg


# --------------------------------------
# -- Create the Element  ---------------
# --------------------------------------
set intTag 		105${ET};	# Internal elastic section tag
set fTag1zz 	106${ET}; 	# Section tag Mzz 1
set fTag2zz 	107${ET}; 	# Section tag Mzz 2
set fTag1yy 	108${ET}; 	# Section tag Myy 1
set fTag2yy 	109${ET}; 	# Section tag Myy 2
set axialMTag  	110${ET}; 	# Create an elastic material for the actual hinge to be aggregated to the actual PH

set phTag1  	111${ET}; 	# Create an aggregated section with this tag for the actual element
set phTag2  	112${ET}; 	# Create an aggregated section with this tag for the actual element
set tempTag		113${ET};	# temporary section tag to be used in creating bidirectional response

# Create the internal elastic element behaviour
section Elastic $intTag [expr $Ec*1000] $Ag $Izze $Iyye $Gc $J

# Create the plastic hinge section
section Uniaxial $fTag1zz $hingeMTag1zz Mz; 			# Create the PH flexural section about ZZ
section Uniaxial $fTag2zz $hingeMTag2zz Mz; 			# Create the PH flexural section about ZZ
uniaxialMaterial Elastic $axialMTag [expr $Ec*$Ag*1000];	# Create the PH axial material

section Aggregator $phTag1 $axialMTag P $hingeMTag1yy My -section $fTag1zz;		# Aggregate P and Myy behaviour to Mzz behaviour
section Aggregator $phTag2 $axialMTag P $hingeMTag2yy My -section $fTag2zz;		# Aggregate P and Myy behaviour to Mzz behaviour

# Create the whole element
# element forceBeamColumn $eleTag $iNode $jNode $transfTag "HingeRadau $secTagI $LpI $secTagJ $LpJ $secTagInterior" <-mass $massDens> <-iter $maxIters $tol>
element forceBeamColumn $ET $iNode $jNode $GT "HingeRadau $phTag1 $Lp $phTag2 $Lp $intTag"
# element elasticBeamColumn $ET $iNode $jNode $A [expr $Ec*1000] $Gc 0.01 [expr $Iyy*$EIratio] [expr $Izz*$EIratio] $GT

}
