#!/usr/bin/perl

## For easier debuging, require strict perl rules
use strict;

## Specify the number of test events (default is 100)
my $numOfEvents = 100;

## Configure the kinematic generator
my $macTemplate = "GEn/template.mac";
my $macPrefix = "../gen_";
my $kinFile = "GEn/kinematics.txt";

## Some global variables that are needed
my $kinNum = 0;


## Now open the kinematics file and generate a sample macro file for
## each kinematic (Q^2) point.
open(my $genFH, "<", $kinFile ) or die printError($kinFile,$!);
while(my $line = <$genFH>) {
  ## Strip end of line character
  chomp($line);
  ## Only accept lines that do not begin with a #
  if( !($line =~ /^\#/) ) {
    # Found a valid line, let's process it
    my $kin = sprintf("%02d",++$kinNum);

    my ($q2,$beam_e,$theta_bb,$theta_sbs,$th_min,$th_max) = split(" ",$line);

    print "Generating sample macro file for Kin$kin. (Q^2=$q2 GeV^2)\n";
    ## Define the name of the macro file for this kinematic
    my $macFileName = "${macPrefix}${q2}GeV2.mac";

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
close($genFH);

sub printError
{
  print "Could not open file $_[0]: $_[1]\n";
}
