#!/usr/bin/perl -w 

## Summarize results of all paramter explorations - Coop, Coop2, Comp 

system "ls Cooperative_TFs/*_NoiseData.txt > ST1"; 
system "ls Cooperative_TFs/*_Burst_SizeFreq.txt > ST2"; 

print "Stats for Cooperative TF binding ...\n"; 
open(GR1,">Stats_CooperativeTFbinding_Noise_data.txt") or die;

open(ST,"ST1") or die; 
while($st=<ST>)
{
  chomp($st); 

  $av1=0; $av2=0; $av3=0; $av4=0; 
  $sd1=0; $sd2=0; $sd3=0; $sd4=0; 

  $pp=0; 
  $pp1=0; $pp2=0; $pp3=0; $pp4=0;  

  open(FP,$st) or die; 
  $ln=0; 
  while($fp=<FP>)
  {
     chomp($fp); 
     if($ln==0 || $ln==1) 
     {
	if($ln==0) 
	{ 
	  print GR1 "$fp\t"; 
          print "   $fp\t"; 
        } 
	$ln++; 
	next; 
     }
     
     $pp++; 

     @arr=split(/\s+/,$fp); 
     if($arr[7]!~/NA/)
     {
       $av1+=$arr[7]; 
       $sd1+=($arr[7])**2; 
       $pp1++;
     }
  
     if($arr[9]!~/NA/)
     {
       $av2+=$arr[9];
       $sd2+=($arr[9])**2; 
       $pp2++;
     }

     if($arr[17]!~/NA/)
     {
       $av3+=$arr[17];
       $sd3+=($arr[17])**2; 
       $pp3++;
     }

     if($arr[19]!~/NA/)
     {
       $av4+=$arr[19];
       $sd4+=($arr[19])**2; 
       $pp4++; 
     }
     undef(@arr); 
  }
  close(FP); 

  $sd1=(($sd1/$pp1)-(($av1/$pp1)**2));
  $sd2=(($sd2/$pp2)-(($av2/$pp2)**2));
  $sd3=(($sd3/$pp3)-(($av3/$pp3)**2));
  $sd4=(($sd4/$pp4)-(($av4/$pp4)**2));


  if($sd1<0) { $sd1=0; } 
  if($sd2<0) { $sd2=0; } 
  if($sd3<0) { $sd3=0; } 
  if($sd4<0) { $sd4=0; } 

  $sd1=sprintf("%0.4f",sqrt($sd1));
  $sd2=sprintf("%0.4f",sqrt($sd2));
  $sd3=sprintf("%0.4f",sqrt($sd3));
  $sd4=sprintf("%0.4f",sqrt($sd4));

  $av1=sprintf("%0.4f",$av1/$pp1);
  $av2=sprintf("%0.4f",$av2/$pp2);
  $av3=sprintf("%0.4f",$av3/$pp3);
  $av4=sprintf("%0.4f",$av4/$pp4);

  print GR1 "Si_Mean_AvgProt $av1  Si_Sd_AvgProt $sd1  Si_Mean_CVProt $av2  Si_Sd_CVProt $sd2\t";
  print GR1 "Cp_Mean_AvgProt $av3  Cp_Sd_AvgProt $sd3  Cp_Mean_CVProt $av4  Cp_Sd_CVProt $sd4\n";

  $n1=sprintf("%0.2f",($pp-$pp1)/$pp*100); 
  $n2=sprintf("%0.2f",($pp-$pp2)/$pp*100); 
  $n3=sprintf("%0.2f",($pp-$pp3)/$pp*100); 
  $n4=sprintf("%0.2f",($pp-$pp4)/$pp*100); 
  print "%NAs $n1 $n2 $n3 $n4\n"; 
}
close(ST); 
close(GR1); 

open(GR2,">Stats_CooperativeTFbinding_BurstSizeFreq.txt") or die;

open(ST,"ST2") or die; 
while($st=<ST>)
{
  chomp($st); 

  $av1=0; $av2=0; $av3=0; $av4=0; 
  $sd1=0; $sd2=0; $sd3=0; $sd4=0; 
  $numcell=0; 
  $pp1=0; $pp2=0; $pp3=0; $pp4=0; 

  open(FP,$st) or die; 
  $ln=0; 
  while($fp=<FP>)
  {
     chomp($fp); 
     if($ln==0 || $ln==1) 
     {
	if($ln==0) { print GR2 "$fp\t"; } 
	if($ln==0) { print "  $fp\t"; } 
	$ln++; 
	next; 
     }
     
     $numcell++; 

     @arr=split(/\s+/,$fp); 
     if($arr[1]!~/NA/)
     {
       $av1+=$arr[1]; 
       $sd1+=($arr[1])**2; 
       $pp1++; 
     }
  
     if($arr[2]!~/NA/)
     {
       $av2+=$arr[2];
       $sd2+=($arr[2])**2; 
       $pp2++; 
     }

     if($arr[4]!~/NA/)
     {
       $av3+=$arr[4];
       $sd3+=($arr[4])**2; 
       $pp3++; 
     }

     if($arr[5]!~/NA/)
     {
       $av4+=$arr[5];
       $sd4+=($arr[5])**2; 
       $pp4++; 
     }
     undef(@arr); 
  }
  close(FP); 
  
  $sd1=(($sd1/$pp1)-(($av1/$pp1)**2));
  $sd2=(($sd2/$pp2)-(($av2/$pp2)**2));
  $sd3=(($sd3/$pp3)-(($av3/$pp3)**2));
  $sd4=(($sd4/$pp4)-(($av4/$pp4)**2));

  if($sd1<0) { $sd1=0; } 
  if($sd2<0) { $sd2=0; } 
  if($sd3<0) { $sd3=0; } 
  if($sd4<0) { $sd4=0; } 

  $sd1=sprintf("%0.4f",sqrt($sd1));
  $sd2=sprintf("%0.4f",sqrt($sd2));
  $sd3=sprintf("%0.4f",sqrt($sd3));
  $sd4=sprintf("%0.4f",sqrt($sd4));

  $av1=sprintf("%0.4f",$av1/$pp1); 
  $av2=sprintf("%0.4f",$av2/$pp2); 
  $av3=sprintf("%0.4f",$av3/$pp3); 
  $av4=sprintf("%0.4f",$av4/$pp4); 

  print GR2 "Si_Mean_BurstSize $av1  Si_Sd_BurstSize $sd1  Si_Mean_BurstFreq $av2  Si_Sd_BurstFreq $sd2\t";
  print GR2 "Cp_Mean_BurstSize $av3  Cp_Sd_BurstSize $sd3  Cp_Mean_BurstFreq $av4  Cp_Sd_BurstFreq $sd4\n";

  $n1=sprintf("%0.2f",($numcell-$pp1)/$numcell*100); 
  $n2=sprintf("%0.2f",($numcell-$pp2)/$numcell*100); 
  $n3=sprintf("%0.2f",($numcell-$pp3)/$numcell*100); 
  $n4=sprintf("%0.2f",($numcell-$pp4)/$numcell*100); 
  print "%NAs $n1 $n2 $n3 $n4\n"; 
}
close(ST); 

close(GR2);

system "rm ST1";
system "rm ST2";

## For Cooperative2 
#
print "Stats2 for Cooperative TF binding ...\n"; 

system "ls Cooperative_TFs2/*_NoiseData.txt > ST1"; 
system "ls Cooperative_TFs2/*_Burst_SizeFreq.txt > ST2"; 

open(GR1,">Stats2_CooperativeTFbinding_Noise_data.txt") or die;

open(ST,"ST1") or die; 
while($st=<ST>)
{
  chomp($st); 

  $av1=0; $av2=0; $av3=0; $av4=0; 
  $sd1=0; $sd2=0; $sd3=0; $sd4=0; 
  $pp=0; 
  $pp1=0; $pp2=0; $pp3=0; $pp4=0; 

  open(FP,$st) or die; 
  $ln=0; 
  while($fp=<FP>)
  {
     chomp($fp); 
     if($ln==0 || $ln==1) 
     {
	if($ln==0) { print GR1 "$fp\t"; } 
	if($ln==0) { print "   $fp\t"; } 
	$ln++; 
	next; 
     }
     
     $pp++; 

     @arr=split(/\s+/,$fp); 
     if($arr[7]!~/NA/)
     {
       $av1+=$arr[7]; 
       $sd1+=($arr[7])**2; 
       $pp1++; 
     }
  
     if($arr[9]!~/NA/)
     {
       $av2+=$arr[9];
       $sd2+=($arr[9])**2; 
       $pp2++; 
     }

     if($arr[17]!~/NA/)
     {
       $av3+=$arr[17];
       $sd3+=($arr[17])**2; 
       $pp3++; 
     }

     if($arr[19]!~/NA/)
     {
       $av4+=$arr[19];
       $sd4+=($arr[19])**2; 
       $pp4++; 
     }
     undef(@arr); 
  }
  close(FP); 

  $sd1=(($sd1/$pp1)-(($av1/$pp1)**2));
  $sd2=(($sd2/$pp2)-(($av2/$pp2)**2));
  $sd3=(($sd3/$pp3)-(($av3/$pp3)**2));
  $sd4=(($sd4/$pp4)-(($av4/$pp4)**2));

  if($sd1<0) { $sd1=0; } 
  if($sd2<0) { $sd2=0; } 
  if($sd3<0) { $sd3=0; } 
  if($sd4<0) { $sd4=0; } 

  $sd1=sprintf("%0.4f",sqrt($sd1));
  $sd2=sprintf("%0.4f",sqrt($sd2));
  $sd3=sprintf("%0.4f",sqrt($sd3));
  $sd4=sprintf("%0.4f",sqrt($sd4));

  $av1=sprintf("%0.4f",$av1/$pp1);
  $av2=sprintf("%0.4f",$av2/$pp2);
  $av3=sprintf("%0.4f",$av3/$pp3);
  $av4=sprintf("%0.4f",$av4/$pp4);

  print GR1 "Si_Mean_AvgProt $av1  Si_Sd_AvgProt $sd1  Si_Mean_CVProt $av2  Si_Sd_CVProt $sd2\t";
  print GR1 "Cp_Mean_AvgProt $av3  Cp_Sd_AvgProt $sd3  Cp_Mean_CVProt $av4  Cp_Sd_CVProt $sd4\n";

  $n1=sprintf("%0.2f",($pp-$pp1)/$pp*100); 
  $n2=sprintf("%0.2f",($pp-$pp2)/$pp*100); 
  $n3=sprintf("%0.2f",($pp-$pp3)/$pp*100); 
  $n4=sprintf("%0.2f",($pp-$pp4)/$pp*100); 
  print "%NAs $n1 $n2 $n3 $n4\n"; 
}
close(ST); 
close(GR1); 


open(GR2,">Stats2_CooperativeTFbinding_BurstSizeFreq.txt") or die;

open(ST,"ST2") or die; 
while($st=<ST>)
{
  chomp($st); 

  $av1=0; $av2=0; $av3=0; $av4=0; 
  $sd1=0; $sd2=0; $sd3=0; $sd4=0; 
  $numcell=0; 
  $pp1=0; $pp2=0; $pp3=0; $pp4=0; 

  open(FP,$st) or die; 
  $ln=0; 
  while($fp=<FP>)
  {
     chomp($fp); 
     if($ln==0 || $ln==1) 
     {
	if($ln==0) { print GR2 "$fp\t"; } 
	if($ln==0) { print "   $fp\t"; } 
	$ln++; 
	next; 
     }
     
     $numcell++; 

     @arr=split(/\s+/,$fp); 
     if($arr[1]!~/NA/)
     {
       $av1+=$arr[1]; 
       $sd1+=($arr[1])**2; 
       $pp1++; 
     }
  
     if($arr[2]!~/NA/)
     {
       $av2+=$arr[2];
       $sd2+=($arr[2])**2; 
       $pp2++; 
     }

     if($arr[4]!~/NA/)
     {
       $av3+=$arr[4];
       $sd3+=($arr[4])**2; 
       $pp3++; 
     }

     if($arr[5]!~/NA/)
     {
       $av4+=$arr[5];
       $sd4+=($arr[5])**2; 
       $pp4++; 
     }
     undef(@arr); 
  }
  close(FP); 
  
  $sd1=(($sd1/$pp1)-(($av1/$pp1)**2));
  $sd2=(($sd2/$pp2)-(($av2/$pp2)**2));
  $sd3=(($sd3/$pp3)-(($av3/$pp3)**2));
  $sd4=(($sd4/$pp4)-(($av4/$pp4)**2));

  if($sd1<0) { $sd1=0; } 
  if($sd2<0) { $sd2=0; } 
  if($sd3<0) { $sd3=0; } 
  if($sd4<0) { $sd4=0; } 

  $sd1=sprintf("%0.4f",sqrt($sd1));
  $sd2=sprintf("%0.4f",sqrt($sd2));
  $sd3=sprintf("%0.4f",sqrt($sd3));
  $sd4=sprintf("%0.4f",sqrt($sd4));

  $av1=sprintf("%0.4f",$av1/$pp1); 
  $av2=sprintf("%0.4f",$av2/$pp2); 
  $av3=sprintf("%0.4f",$av3/$pp3); 
  $av4=sprintf("%0.4f",$av4/$pp4); 

  print GR2 "Si_Mean_BurstSize $av1  Si_Sd_BurstSize $sd1  Si_Mean_BurstFreq $av2  Si_Sd_BurstFreq $sd2\t";
  print GR2 "Cp_Mean_BurstSize $av3  Cp_Sd_BurstSize $sd3  Cp_Mean_BurstFreq $av4  Cp_Sd_BurstFreq $sd4\n";

  $n1=sprintf("%0.2f",($numcell-$pp1)/$numcell*100); 
  $n2=sprintf("%0.2f",($numcell-$pp2)/$numcell*100); 
  $n3=sprintf("%0.2f",($numcell-$pp3)/$numcell*100); 
  $n4=sprintf("%0.2f",($numcell-$pp4)/$numcell*100); 
  print "%NAs $n1 $n2 $n3 $n4\n"; 
}
close(ST); 

close(GR2);

system "rm ST1";
system "rm ST2";

## For Competitive TFbinding 

print "Stats for Competitive TF binding ...\n"; 

system "ls Competitive_TFs/*_NoiseData.txt > ST1"; 
system "ls Competitive_TFs/*_Burst_SizeFreq.txt > ST2"; 

open(GR1,">Stats_CompetitiveTFbinding_Noise_data.txt") or die;

open(ST,"ST1") or die; 
while($st=<ST>)
{
  chomp($st); 

  $av1=0; $av2=0; $av3=0; $av4=0; 
  $sd1=0; $sd2=0; $sd3=0; $sd4=0; 
  $pp=0; 
  $pp1=0; $pp2=0; $pp3=0; $pp4=0; 

  open(FP,$st) or die; 
  $ln=0; 
  while($fp=<FP>)
  {
     chomp($fp); 
     if($ln==0 || $ln==1) 
     {
	if($ln==0) { print GR1 "$fp\t"; } 
	if($ln==0) { print "   $fp\t"; } 
	$ln++; 
	next; 
     }
     
     $pp++; 

     @arr=split(/\s+/,$fp); 

     if($arr[7]!~/NA/)
     {
       $av1+=$arr[7]; 
       $sd1+=($arr[7])**2; 
       $pp1++; 
     }
  
     if($arr[9]!~/NA/)
     {
       $av2+=$arr[9];
       $sd2+=($arr[9])**2; 
       $pp2++; 
     }

     if($arr[17]!~/NA/)
     {
       $av3+=$arr[17];
       $sd3+=($arr[17])**2; 
       $pp3++; 
     }

     if($arr[19]!~/NA/)
     {
       $av4+=$arr[19];
       $sd4+=($arr[19])**2; 
       $pp4++; 
     }
     undef(@arr); 
  }
  close(FP); 

  $sd1=(($sd1/$pp1)-(($av1/$pp1)**2));
  $sd2=(($sd2/$pp2)-(($av2/$pp2)**2));
  $sd3=(($sd3/$pp3)-(($av3/$pp3)**2));
  $sd4=(($sd4/$pp4)-(($av4/$pp4)**2));

  if($sd1<0) { $sd1=0; } 
  if($sd2<0) { $sd2=0; } 
  if($sd3<0) { $sd3=0; } 
  if($sd4<0) { $sd4=0; } 

  $sd1=sprintf("%0.4f",sqrt($sd1));
  $sd2=sprintf("%0.4f",sqrt($sd2));
  $sd3=sprintf("%0.4f",sqrt($sd3));
  $sd4=sprintf("%0.4f",sqrt($sd4));

  $av1=sprintf("%0.4f",$av1/$pp1);
  $av2=sprintf("%0.4f",$av2/$pp2);
  $av3=sprintf("%0.4f",$av3/$pp3);
  $av4=sprintf("%0.4f",$av4/$pp4);

  print GR1 "Si_Mean_AvgProt $av1  Si_Sd_AvgProt $sd1  Si_Mean_CVProt $av2  Si_Sd_CVProt $sd2\t";
  print GR1 "Cm_Mean_AvgProt $av3  Cm_Sd_AvgProt $sd3  Cm_Mean_CVProt $av4  Cm_Sd_CVProt $sd4\n";

  $n1=sprintf("%0.2f",($pp-$pp1)/$pp*100); 
  $n2=sprintf("%0.2f",($pp-$pp2)/$pp*100); 
  $n3=sprintf("%0.2f",($pp-$pp3)/$pp*100); 
  $n4=sprintf("%0.2f",($pp-$pp4)/$pp*100); 
  print "%NAs $n1 $n2 $n3 $n4\n"; 
}
close(ST); 
close(GR1); 


open(GR2,">Stats_CompetitiveTFbinding_BurstSizeFreq.txt") or die;

open(ST,"ST2") or die; 
while($st=<ST>)
{
  chomp($st); 

  $av1=0; $av2=0; $av3=0; $av4=0; 
  $sd1=0; $sd2=0; $sd3=0; $sd4=0; 
  $numcell=0; 
  $pp1=0; $pp2=0; $pp3=0; $pp4=0; 

  open(FP,$st) or die; 
  $ln=0; 
  while($fp=<FP>)
  {
     chomp($fp); 
     if($ln==0 || $ln==1) 
     {
	if($ln==0) { print GR2 "$fp\t"; } 
	if($ln==0) { print "   $fp\t"; } 
	$ln++; 
	next; 
     }
     
     $numcell++; 

     @arr=split(/\s+/,$fp); 
     if($arr[1]!~/NA/)
     {
       $av1+=$arr[1]; 
       $sd1+=($arr[1])**2; 
       $pp1++; 
     }
  
     if($arr[2]!~/NA/)
     {
       $av2+=$arr[2];
       $sd2+=($arr[2])**2; 
       $pp2++; 
     }

     if($arr[4]!~/NA/)
     {
       $av3+=$arr[4];
       $sd3+=($arr[4])**2; 
       $pp3++; 
     }

     if($arr[5]!~/NA/)
     {
       $av4+=$arr[5];
       $sd4+=($arr[5])**2; 
       $pp4++; 
     }
     undef(@arr); 
  }
  close(FP); 
  
  $sd1=(($sd1/$pp1)-(($av1/$pp1)**2));
  $sd2=(($sd2/$pp2)-(($av2/$pp2)**2));
  $sd3=(($sd3/$pp3)-(($av3/$pp3)**2));
  $sd4=(($sd4/$pp4)-(($av4/$pp4)**2));

  if($sd1<0) { $sd1=0; } 
  if($sd2<0) { $sd2=0; } 
  if($sd3<0) { $sd3=0; } 
  if($sd4<0) { $sd4=0; } 

  $sd1=sprintf("%0.4f",sqrt($sd1));
  $sd2=sprintf("%0.4f",sqrt($sd2));
  $sd3=sprintf("%0.4f",sqrt($sd3));
  $sd4=sprintf("%0.4f",sqrt($sd4));

  $av1=sprintf("%0.4f",$av1/$pp1); 
  $av2=sprintf("%0.4f",$av2/$pp2); 
  $av3=sprintf("%0.4f",$av3/$pp3); 
  $av4=sprintf("%0.4f",$av4/$pp4); 

  print GR2 "Si_Mean_BurstSize $av1  Si_Sd_BurstSize $sd1  Si_Mean_BurstFreq $av2  Si_Sd_BurstFreq $sd2\t";
  print GR2 "Cm_Mean_BurstSize $av3  Cm_Sd_BurstSize $sd3  Cm_Mean_BurstFreq $av4  Cm_Sd_BurstFreq $sd4\n";

  $n1=sprintf("%0.2f",($numcell-$pp1)/$numcell*100); 
  $n2=sprintf("%0.2f",($numcell-$pp2)/$numcell*100); 
  $n3=sprintf("%0.2f",($numcell-$pp3)/$numcell*100); 
  $n4=sprintf("%0.2f",($numcell-$pp4)/$numcell*100); 
  print "%NAs $n1 $n2 $n3 $n4\n"; 
}
close(ST); 

close(GR2);

system "rm ST1";
system "rm ST2";
