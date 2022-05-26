#!/usr/bin/perl -w 

open(FP,"Results/IndParam_impact_1000_Iterations.txt") or die; 
open(WR,">Results/Impactful_parameter_list_DM.txt") or die; 
$mpc=0; 
while($fp=<FP>)
{
  chomp($fp); 
  if($fp=~/^Unfiltered/ || $fp=~/^Variable/ || $fp=~/^-+/)
  {
    next; 
  }
  @arr=split(/\s+/,$fp); 
  @qw1=split(/\±/,$arr[1]);
  @qw2=split(/\±/,$arr[2]);
 
  if($qw1[0]>=0.025 || $qw2[0]>=0.025 || ($qw1[0]>=0.025 && $qw2[0]>=0.025)) 
  {
    $mpc++;
    print WR "$arr[0]\n"; 
    ##print "$arr[0]+";
  }
 
  undef(@qw2);
  undef(@qw1);
  undef(@arr);
}
close(WR);
close(FP);

##print "\n";

print "Impactful Params no.: $mpc\n"; 

open(FP,"Results/IndParam_impact_1000_Iterations_nondup.txt") or die; 
open(WR,">Results/Impactful_parameter_list_DM_nondup.txt") or die; 
$mpc=0; 
while($fp=<FP>)
{
  chomp($fp); 
  if($fp=~/^Unfiltered/ || $fp=~/^Variable/ || $fp=~/^-+/)
  {
    next; 
  }
  @arr=split(/\s+/,$fp); 
  @qw1=split(/\±/,$arr[1]);
  @qw2=split(/\±/,$arr[2]);
 
  if($qw1[0]>=0.025 || $qw2[0]>=0.025 || ($qw1[0]>=0.025 && $qw2[0]>=0.025)) 
  {
    $mpc++;
    print WR "$arr[0]\n"; 
    ##print "$arr[0]+";
  }
 
  undef(@qw2);
  undef(@qw1);
  undef(@arr);
}
close(WR);
close(FP);

##print "\n";

print "Impactful Params no. nondup: $mpc\n"; 

