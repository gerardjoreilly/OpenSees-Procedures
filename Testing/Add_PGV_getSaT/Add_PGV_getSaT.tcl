# Test and verify the addition of PGV and PGD calculation to the getSaT function

# Define a ground motion
set EQ "../eq_1g_0.01s.txt"
set dt 0.01
set xi 0.05
set g 9.81

# Read the signal
set imp [open $EQ "r"];
set accg [read $imp];
close $imp

# Integrate it for ground displacement and velocity
set acc 	{};
set vg {0.0}
set dg {0.0}
for {set i 0} {$i < [llength $accg]} {incr i 1} {
  lappend vg [expr [lindex $vg $i]+0.5*$dt*$g*([lindex $accg $i+1]+[lindex $accg $i])]
  lappend dg [expr [lindex $dg $i]+0.5*$dt*([lindex $vg $i+1]+[lindex $vg $i])]
}

# Find the maximum values
set pgv 0.0;
set pgd 0.0;

file delete ground_response.txt
set gr [open ground_response.txt a]

# Loop the terms (This loop is t length)
for {set i 0} {$i < [llength $accg]} {incr i 1} {
  set temp1 [expr abs([lindex $vg $i])]; 		# Get the absolute current ground velocity
  set temp2 [expr abs([lindex $dg $i])]; 	# Get the absolute current ground displacement
  if {$temp1>$pgv} {set pgv $temp1};		# Find the pgv
  if {$temp2>$pgd} {set pgd $temp2};			# Find the pgd


  puts $gr "[lindex $accg $i] [lindex $vg $i] [lindex $dg $i]"

}
close $gr

puts $pgv
puts $pgd

source ../../getSaT.tcl
getSaT $EQ $dt 1.0 $xi
# puts $Sa
# puts $pga
puts $pgv
puts $pgd
