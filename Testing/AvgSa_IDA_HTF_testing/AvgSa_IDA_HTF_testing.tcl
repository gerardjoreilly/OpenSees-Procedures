# This is a testing script for the computation of AvgSa

# Need to load the getSaT function
set procdir "../../"
source $procdir/getSaT.tcl

set IMtype 6
set Tinfo {0.27 0.46 0.65 0.83 1.02 1.21 1.40 1.58 1.77 1.96}

set gmsdir "../"
set EQnameX "eq_1g_0.01s.txt"
set EQnameY "eq_1g_0.01s.txt"
set dts_list [list 0.01]
set xi 0.05
set i 1


if {$IMtype==6} {
  # Get the length of the periods
  set nT [llength $Tinfo]

  # Get the spectral accelerations at each period
  set Sa_listX {};
  set Sa_listY {};
  for {set pn 0} {$pn < $nT} {incr pn 1} {
    getSaT $gmsdir/$EQnameX [lindex $dts_list $i-1] [lindex $Tinfo $pn] $xi; # Get the Sa(T1,5%) of the record in the X direction
    lappend Sa_listX $Sa
    getSaT $gmsdir/$EQnameY [lindex $dts_list $i-1] [lindex $Tinfo $pn] $xi; # Get the Sa(T1,5%) of the record in the Y direction
    lappend Sa_listY $Sa
  }

  # Compute the geometric mean
  set SaXsumprod [lindex $Sa_listX 0]
  set SaYsumprod [lindex $Sa_listY 0]
  for {set pn 1} {$pn < $nT} {incr pn 1} {
    set SaXsumprod [expr [lindex $Sa_listX $pn]*$SaXsumprod]
    set SaYsumprod [expr [lindex $Sa_listX $pn]*$SaYsumprod]
  }
  set SaXgm [expr pow($SaXsumprod,1/($nT*1.0))]
  set SaYgm [expr pow($SaYsumprod,1/($nT*1.0))]

  # Using the geoetreic mean of the two compoents AvgSa
  # This is the same as just takeing the AvgSa of the combined set of X and Y (i.e. at Tinfo x 2 periods)
  set IMgeomean [expr pow($SaXgm*$SaYgm,0.5)]

}

puts $Sa_listX
puts $SaXgm
puts $SaYgm

puts "Finished!"
