#!/usr/bin/perl -w 

## Filter out non-correlating parameters with MED_DM_SD and IQR_DM_SD 

open(FP,"../../datasets/WolfeLab_Ohnolog_list.csv") or die; 
while($fp=<FP>)
{
   chomp($fp); 
   @arr=split(/\,/,$fp); 
   $duplist{$arr[3]}=$arr[5]; 
   $duplist{$arr[5]}=$arr[3]; 
   undef(@arr); 
}
close(FP); 


$lnc=0;
open(WR,">Filtered_nondup_MERGED_DATA_Yeast.txt") or die;
open(FP,"Filtered_MERGED_DATA_Yeast.txt") or die;
while($fp=<FP>)
{
  chomp($fp);
  @arr=split(/\s+/,$fp);
  if($fp=~/^Gene/)
  {
    print WR "$fp\n";  
    next;
  }

  if(exists $duplist{$arr[0]}) 
  {
     $id=$duplist{$arr[0]};
     if(!exists $exlist{$id}) 
     {
        print WR "$fp\n"; 
	$exlist{$arr[0]}=1; 
     }
     else
     {
	$lnc++; 
     }
     undef $id; 
  }
  else
  {
     print WR "$fp\n"; 
  }
  undef(@arr);
}
close(FP);
close(WR); 

undef %exlist; 
undef %duplist; 

print "Duplicates discarded: $lnc\n"; 
