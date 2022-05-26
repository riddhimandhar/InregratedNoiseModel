#!/usr/bin/perl -w 
#

## Summarize the results of the onrate and offrate explortation 
#

system "ls Plots_sim_onrate_offrate_expl/NoiseData_*.txt > TMP1"; 

open(TM,"TMP1") or die; 

while($tm=<TM>) 
{
   chomp($tm); 
   @qw1=split(/\//,$tm); 
   @qw2=split(/\./,$qw1[1]); 
   @we=split(/\_/,$qw2[0]); 
   @we2=split(/binding/,$we[1]); 
   $id=$we2[0]."_$we[2]"."_$we[3]"; 

   $av1=0; $av2=0; $sd1=0; $sd2=0; 
   $cc=0; 

   print "$tm\n"; 
   open(FP,$tm) or die; 
   while($fp=<FP>)
   {
      chomp($fp); 
      if($fp=~/^Timepoint/) { next; } 
      @arr=split(/\s+/,$fp); 
      if($arr[7] ne "NA" && $arr[9] ne "NA") 
      {
        $av1+=$arr[7];
        $sd1+=($arr[7])**2; 
        $av2+=$arr[9];
        $sd2+=($arr[9])**2; 
	$cc++; 
      }
      undef(@arr); 
   }
   close(FP); 

   undef(@we2); 
   undef(@we); 
   undef(@qw2); 
   undef(@qw1); 

   if($cc!=0) 
   {
     $sd1=(($sd1/$cc)-($av1/$cc)**2); 
     if($sd1<0) { $sd1=0; }
     $sd1=sprintf("%0.2f",sqrt($sd1)); 
     $av1=sprintf("%0.2f",($av1/$cc)); 

     $sd2=(($sd2/$cc)-($av2/$cc)**2); 
     if($sd2<0) { $sd2=0; } 
     $sd2=sprintf("%0.4f",sqrt($sd2)); 
     $av2=sprintf("%0.4f",($av2/$cc)); 
   }
   else
   {
     $av1="NA"; $sd1="NA"; $av2="NA"; $sd2="NA"; 
   }
   $list{$id}[0]=$av1; 
   $list{$id}[1]=$sd1; 
   $list{$id}[2]=$av2; 
   $list{$id}[3]=$sd2; 

   undef $id; 
}
close(TM); 

system "rm TMP1"; 

system "ls Plots_sim_onrate_offrate_expl/*_Burst_SizeFreq*.txt > TMP2"; 

open(TM,"TMP2") or die; 

while($tm=<TM>) 
{
   chomp($tm); 
   @qw1=split(/\//,$tm); 
   @qw2=split(/\_/,$qw1[1]); 
   @we=split(/\./,$qw2[4]); 
   $id=$qw2[0]."_$qw2[3]"."_$we[0]"; 

   $av1=0; $av2=0; $sd1=0; $sd2=0; 
   $cc=0; 

   print "$tm\n"; 
   open(FP,$tm) or die; 
   while($fp=<FP>)
   {
      chomp($fp); 
      if($fp=~/^CellID/) { next; } 
      @arr=split(/\s+/,$fp); 
      if($arr[1] ne "NA" && $arr[2] ne "NA") 
      {
        $av1+=$arr[1];
        $sd1+=($arr[1])**2; 
        $av2+=$arr[2];
        $sd2+=($arr[2])**2; 
	$cc++; 
      }
      undef(@arr); 
   }
   close(FP); 

   undef(@we); 
   undef(@qw2); 
   undef(@qw1); 

   if($cc!=0) 
   {
     $sd1=(($sd1/$cc)-($av1/$cc)**2); 
     if($sd1<0) { $sd1=0; }
     $sd1=sprintf("%0.2f",sqrt($sd1)); 
     $av1=sprintf("%0.2f",($av1/$cc)); 

     $sd2=(($sd2/$cc)-($av2/$cc)**2); 
     if($sd2<0) { $sd2=0; } 
     $sd2=sprintf("%0.4f",sqrt($sd2)); 
     $av2=sprintf("%0.4f",($av2/$cc)); 
   }
   else
   {
     $av1="NA"; $sd1="NA"; $av2="NA"; $sd2="NA"; 
   }

   $list{$id}[4]=$av1; 
   $list{$id}[5]=$sd1; 
   $list{$id}[6]=$av2; 
   $list{$id}[7]=$sd2; 

   undef $id; 
}
close(TM); 

system "rm TMP2"; 

open(WR,">Plots_sim_onrate_offrate_expl/Overall_Results_Stats_Onrate_Offrate_expl.txt") or die; 
print WR "TFbinding OnRate OffRate Mean_AvgProt Sd_AvgProt Mean_CVProt Sd_CVProt Mean_BurstSize Sd_BurstSize Mean_BurstFreq Sd_BurstFreq\n"; 

foreach $key (sort keys %list) 
{
  @we=split(/\_/,$key); 
  print WR "$we[0] $we[1] $we[2] $list{$key}[0] $list{$key}[1] $list{$key}[2] $list{$key}[3] $list{$key}[4] $list{$key}[5] $list{$key}[6] $list{$key}[7]\n";  undef(@we); 
}
close(WR); 

undef %list; 
