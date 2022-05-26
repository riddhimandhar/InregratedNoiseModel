#!/usr/bin/perl -w 

## Filter MERGED_DATA_Yeast.txt file for NAs 

use Statistics::R;
$R=Statistics::R->new();
$R->startR;

open(FP,"../../MERGED_DATA_Yeast.txt") or die; 
open(WR,">ColNAremoved_MERGED_DATA_Yeast.txt") or die; 
while($fp=<FP>)
{
  chomp($fp); 
  if($fp=~/^Gene/)
  {
    print WR "$fp\n"; 
    next;
  }

  @arr=split(/\s+/,$fp); 
  if($arr[3]!~/NA/)
  {
    print WR "$fp\n"; 
  }
  undef(@arr); 
}
close(WR);
close(FP); 

$oritab="ColNAremoved_MERGED_DATA_Yeast.txt"; 

$R->run(qq'tab=read.table("$oritab",header=T)'); 
##$R->run(q'df=colSums(is.na(tab))');
$R->run(q'newtab=tab[,colSums(is.na(tab))<0.8*nrow(tab)]');
$R->run(q'write.table(newtab,"ColNAremoved_MERGED_DATA_Yeast.txt",row.names=F,sep="\t",quote=F)'); 

open(FP,"ColNAremoved_MERGED_DATA_Yeast.txt") or die; 
open(WR,">Filtered_MERGED_DATA_Yeast.txt") or die; 
while($fp=<FP>)
{
  chomp($fp); 
  if($fp=~/^Gene/)
  {
    print WR "$fp\n"; 
    next;
  }

  @arr=split(/\s+/,$fp); 
  $no=@arr;
  $masterfl=0; 
  if($arr[5]!~/NA/)
  {
    for($i=6;$i<$no;$i++)
    {
      if($arr[$i]!~/NA/ && $arr[$i]!=0) { $masterfl++; }
    }
  }
  if($masterfl/$no>0.1)
  {
    print WR "$fp\n"; 
  }

  undef $masterfl; 
  undef $no; 
  undef(@arr); 
}
close(WR);
close(FP); 

$R->stopR; 
