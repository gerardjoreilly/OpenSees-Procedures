# Units
# --------------
# Standardised Units
# --------------
set pi 3.141592654
set g 9.81

# Length
set m 1.0;
set mm [expr $m/1000.0];
set cm [expr $m/100.0];
set in [expr 25.4*$mm];
set ft [expr 12*$in];

# Area
set m2 [expr $m*$m];
set mm2 [expr $mm*$mm];
set cm2 [expr $cm*$cm];
set in2 [expr $in*$in];

# Second Moment of Area
set m4 [expr $m*$m*$m*$m];
set cm4 [expr $cm*$cm*$cm*$cm];
set mm4 [expr $mm*$mm*$mm*$mm];
set in4 [expr $in*$in*$in*$in];

# Force
set kN 1.0;
set N [expr $kN/1000.0];
set kips [expr $kN*4.448221615];

# Moment
set kNm [expr $kN*$m];

# Mass
set kg [expr $N/$g];
set tonne [expr $kg*1000];

# Stress
set Pa [expr $N/($m*$m)];
set kPa [expr $Pa*1.0e3];
set MPa [expr $Pa*1.0e6];
set Nmm2 [expr $N/($mm*$mm)];
set kNmm2 [expr $Nmm2*1.0e3];
set GPa [expr $Pa*1.0e9];
set ksi [expr 6.8947573*$MPa];

# Angles
set degrees [expr $pi/180];
