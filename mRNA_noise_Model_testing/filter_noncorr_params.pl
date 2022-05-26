#!/usr/bin/perl -w 

## Filter out non-correlating parameters with DM_YPD  

use Statistics::R;
$R=Statistics::R->new();
$R->startR;

$lnc=0;
open(FP,"Filtered_MERGED_DATA_Yeast.txt") or die;
while($fp=<FP>)
{
  chomp($fp);
  @arr=split(/\s+/,$fp);
  $no=@arr;
  if($fp=~/^Gene/)
  {
    for($i=0;$i<$no;$i++)
    {
      $head[1][$i]=$arr[$i];
      $head[0]=$no;
    }
    next;
  }

  for($i=0;$i<$no;$i++)
  {
    $list[$lnc][$i]=$arr[$i];
  }
  $lnc++;

  undef $no;
  undef(@arr);
}
close(FP);

$ptcol=3;
for($mm=0;$mm<1;$mm++)
{
 $ofs=0; 
 if($head[1][$ptcol+$mm] eq "DM_YPD") { $ofs=1; }

 for($j=$ptcol+$mm+2+$ofs;$j<$head[0];$j++)
 {
   open(GR,">TmpFilt.txt") or die;
   print "$head[1][$j] $head[1][$ptcol+$mm]\t";
   $gt=0;
   for($k=0;$k<$lnc;$k++)
   {
     if($list[$k][$j]!~/NA/ && $list[$k][$ptcol+$mm]!~/NA/) 
     {
       print GR "$list[$k][$j] $list[$k][$ptcol+$mm]\n";
       $gt++;
     }
   }
   close(GR);

   if($gt>2)
   {
   $R->run(q'tab=read.table("TmpFilt.txt")');
   $R->run(q'rv=cor.test(tab$V1,tab$V2)');
   $R->run(q'rves=rv$estimate');
   $R->run(q'rvp=rv$p.value');
   $rves=$R->get('rves');
   $rvp=$R->get('rvp');
   $R->run(q'rm(rv)');
   $R->run(q'rm(rves)');
   $R->run(q'rm(rvp)');
   $R->run(q'rm(tab)');
   system "rm TmpFilt.txt";
 
   print "$rves $rvp\n";
  
   if($rvp!~/NA/ && $rvp<0.05) 
   {
     $flag{$j}=1;  
   }  
   undef $rves;
   undef $rvp;
   }
   else { print "NA NA\n"; }
 }

 if($head[1][$ptcol+$mm] eq "DM_YPD") 
 { 
   open(WR,">CorrFiltered_MERGED_DATA_Yeast.txt") or die; 
 }

 for($i=0;$i<=$ptcol+$mm;$i++)
 {
   print WR "$head[1][$i]\t";
 }
 for($i=$ptcol+$mm+2+$ofs;$i<$head[0];$i++)
 {
   if(exists $flag{$i})  
   {
     print WR "$head[1][$i]\t"; 
   }
 }
 print WR "\n"; 

 for($k=0;$k<$lnc;$k++)
 {
   for($i=0;$i<=$ptcol+$mm;$i++)
   {
     print WR "$list[$k][$i]\t";
   }
   for($i=$ptcol+$mm+2+$ofs;$i<$head[0];$i++)
   {
     if(exists $flag{$i})  
     {
       print WR "$list[$k][$i]\t"; 
     }
   }
   print WR "\n"; 
 }
 close(WR);
 undef %flag; 
}

undef(@head);
undef(@list);

## For nondup data 

$lnc=0;
open(FP,"Filtered_nondup_MERGED_DATA_Yeast.txt") or die;
while($fp=<FP>)
{
  chomp($fp);
  @arr=split(/\s+/,$fp);
  $no=@arr;
  if($fp=~/^Gene/)
  {
    for($i=0;$i<$no;$i++)
    {
      $head[1][$i]=$arr[$i];
      $head[0]=$no;
    }
    next;
  }

  for($i=0;$i<$no;$i++)
  {
    $list[$lnc][$i]=$arr[$i];
  }
  $lnc++;

  undef $no;
  undef(@arr);
}
close(FP);

$ptcol=3;
for($mm=0;$mm<1;$mm++)
{
 $ofs=0; 
 if($head[1][$ptcol+$mm] eq "DM_YPD") { $ofs=1; }

 for($j=$ptcol+$mm+2+$ofs;$j<$head[0];$j++)
 {
   open(GR,">TmpFilt.txt") or die;
   print "$head[1][$j] $head[1][$ptcol+$mm]\t";
   $gt=0;
   for($k=0;$k<$lnc;$k++)
   {
     if($list[$k][$j]!~/NA/ && $list[$k][$ptcol+$mm]!~/NA/) 
     {
       print GR "$list[$k][$j] $list[$k][$ptcol+$mm]\n";
       $gt++;
     }
   }
   close(GR);

   if($gt>2)
   {
   $R->run(q'tab=read.table("TmpFilt.txt")');
   $R->run(q'rv=cor.test(tab$V1,tab$V2)');
   $R->run(q'rves=rv$estimate');
   $R->run(q'rvp=rv$p.value');
   $rves=$R->get('rves');
   $rvp=$R->get('rvp');
   $R->run(q'rm(rv)');
   $R->run(q'rm(rves)');
   $R->run(q'rm(rvp)');
   $R->run(q'rm(tab)');
   system "rm TmpFilt.txt";
 
   print "$rves $rvp\n";
  
   if($rvp!~/NA/ && $rvp<0.05) 
   {
     $flag{$j}=1;  
   }  
   undef $rves;
   undef $rvp;
   }
   else { print "NA NA\n"; }
 }

 if($head[1][$ptcol+$mm] eq "DM_YPD") 
 { 
   open(WR,">CorrFiltered_nondup_MERGED_DATA_Yeast.txt") or die; 
 }

 for($i=0;$i<=$ptcol+$mm;$i++)
 {
   print WR "$head[1][$i]\t";
 }
 for($i=$ptcol+$mm+2+$ofs;$i<$head[0];$i++)
 {
   if(exists $flag{$i})  
   {
     print WR "$head[1][$i]\t"; 
   }
 }
 print WR "\n"; 

 for($k=0;$k<$lnc;$k++)
 {
   for($i=0;$i<=$ptcol+$mm;$i++)
   {
     print WR "$list[$k][$i]\t";
   }
   for($i=$ptcol+$mm+2+$ofs;$i<$head[0];$i++)
   {
     if(exists $flag{$i})  
     {
       print WR "$list[$k][$i]\t"; 
     }
   }
   print WR "\n"; 
 }
 close(WR);
 undef %flag; 
}

undef(@head);
undef(@list);
$R->stopR; 
