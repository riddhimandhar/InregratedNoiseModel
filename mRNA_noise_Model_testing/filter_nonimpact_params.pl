#!/usr/bin/perl -w 

## Filter out parameters that do not show any impact on MED_DM 

open(FP,"Results/Impactful_parameter_list_DM.txt") or die; 
while($fp=<FP>)
{
  chomp($fp); 
  $implist{$fp}=1; 
}
close(FP);

open(WR,">ImpFiltered_MERGED_DATA_Yeast.txt") or die; 
open(FP,"Filtered_MERGED_DATA_Yeast.txt") or die;
while($fp=<FP>)
{
  chomp($fp);
  @arr=split(/\s+/,$fp);
  $no=@arr;
  if($fp=~/^Gene/)
  {
    for($i=0;$i<6;$i++)
    {
      $flag{$i}=1; 
      print WR "$arr[$i]\t"; 
    }
    for($i=6;$i<$no;$i++)
    {
      if(exists $implist{$arr[$i]})
      {
        $flag{$i}=1; 
        print WR "$arr[$i]\t"; 
      }
    }
    print WR "\n"; 
    next;
  }

  for($i=0;$i<$no;$i++)
  {
    if(exists $flag{$i})
    {
      print WR "$arr[$i]\t"; 
    }
  }
  print WR "\n"; 

  undef $no;
  undef(@arr);
}
close(WR);
close(FP);

undef %flag; 
undef %implist; 

## For nondup data 

open(FP,"Results/Impactful_parameter_list_DM_nondup.txt") or die; 
while($fp=<FP>)
{
  chomp($fp); 
  $implist{$fp}=1; 
}
close(FP);

open(WR,">ImpFiltered_nondup_MERGED_DATA_Yeast.txt") or die; 
open(FP,"Filtered_nondup_MERGED_DATA_Yeast.txt") or die;
while($fp=<FP>)
{
  chomp($fp);
  @arr=split(/\s+/,$fp);
  $no=@arr;
  if($fp=~/^Gene/)
  {
    for($i=0;$i<6;$i++)
    {
      $flag{$i}=1; 
      print WR "$arr[$i]\t"; 
    }
    for($i=6;$i<$no;$i++)
    {
      if(exists $implist{$arr[$i]})
      {
        $flag{$i}=1; 
        print WR "$arr[$i]\t"; 
      }
    }
    print WR "\n"; 
    next;
  }

  for($i=0;$i<$no;$i++)
  {
    if(exists $flag{$i})
    {
      print WR "$arr[$i]\t"; 
    }
  }
  print WR "\n"; 

  undef $no;
  undef(@arr);
}
close(WR);
close(FP);

undef %flag; 
undef %implist; 
