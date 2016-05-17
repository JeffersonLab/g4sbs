#!/usr/bin/perl

## For easier debuging, require strict perl rules
use strict;

## Specify the number of test events (default is 100)
my $numOfEvents = 100;

## Configure the kinematic generator
my $macTemplate = "GMn/template.mac";
my $macPrefix = "../gmn_";
my $kinFile = "GMn/kinematics.txt";

## Some global variables that are needed
my $kinNum = 0;


## Now open the kinematics file and generate a sample macro file for
## each kinematic (Q^2) point.
open(my $gmnFH, "<", $kinFile ) or die printError($kinFile,$!);
while(my $line = <$gmnFH>) {
  ## Strip end of line character
  chomp($line);
  ## Only accept lines that do not begin with a #
  if( !($line =~ /^\#/) ) {
    # Found a valid line, let's process it
    my $kin = sprintf("%02d",++$kinNum);

    my ($tmp,$q2,$beam_e,$theta_bb,$theta_sbs,$dist_bb,$dist_mag,$dist_hcal,$bdl,
      $voff_hcal,$th_min,$th_max,$ph_min,$ph_max) = split(" ",$line);
    $th_min=$theta_bb-10;
    $th_max=$theta_bb+10;

    ## Determine magnetic field using the standard 48 inch length of 48D48
    my $field_mag = $bdl/1.2192; ## Given in Tesla
    ## Round off to 2 decimal places
    $field_mag = sprintf("%.2f",$field_mag);

    print "Generating sample macro file for Kin$kin. (Q^2=$q2 GeV^2)\n";
    ## Define the name of the macro file for this kinematic
    my $macFileName = "${macPrefix}kin${kin}.mac";

    ## Open up the macro file for output
    open(my $outMacFH, ">",$macFileName)
      or die printError($macFileName,$!);

    ## Generate the sample macro for this kinematic
    my $tmpBuffer;

    ## Produce the general macro
    open(my $inMacFH, "<", "$macTemplate") or die printError("$macTemplate");
    while(<$inMacFH>) {
      s/%%(\$\w+)%%/$1/eeg;
      print $outMacFH $_;
    }
    close($inMacFH);

    ## Close output macro file
    close($outMacFH);
  }
}
close($gmnFH);

sub printError
{
  print "Could not open file $_[0]: $_[1]\n";
}
