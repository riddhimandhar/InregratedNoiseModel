#!/usr/bin/perl -w 

## Explore the impact of paramter values, number of TFs in cooperative TF binding in 
## simulating TF cooperative binding 
## Calculate impact on gene expression 

print "Calculating and plotting time-variant dynamics ....\n"; 

use Statistics::R;

$numevt=500; ## Num of events 
$bas=0; ## Basal transcription rate  
$tottime=0.0500; ## Total time of simulation 
$uot=0.0001;  ## Time interval 

$inimrna=10; 
$iniprot=1000; 

## Variables - No. of TFs, synthesis rate, onrate, offrate, degradation rates 
## --------------------------
 
$mrnarate=100; ## mRNA production rate  

$numtf=2; 

$onrate=1000; ## rate parameter exponential distribution for onrate 
$offrate=1200; ## rate parameter exponential distribution for offrate 

$protrate=100; ## protein synthesis rate per mRNA molecule per unit time 

$degm=0.1; ## mRNA degradation rate per mRNA molecule per unit time 
$degp=0.02; ## protein degradation rate per protein moleculer per unit time 

$coopfac=2.25; 

## Parameter exploration - combination of numTF and coopfac  
#

$nmt[0]=8; ## $numtf 
$nmt[1][0]=2; $nmt[1][1]=3; $nmt[1][2]=4; $nmt[1][3]=5; $nmt[1][4]=6; $nmt[1][5]=7; $nmt[1][6]=8; $nmt[1][7]=9; ## 2 to 9 
## $rates -> average rate $mrnarate 

$cpf[0]=10; ## $coopfac
$cpf[1][0]=1; $cpf[1][1]=2; $cpf[1][2]=3; $cpf[1][3]=4; $cpf[1][4]=5; $cpf[1][5]=6; $cpf[1][6]=7; $cpf[1][7]=8; $cpf[1][8]=9; $cpf[1][9]=10;  


$runfl=0; 
$modrun=1; 

$mp=0; 

while($runfl==0)
{
  ## Parameter value reset 
  
  $mrnarate=100; ## mRNA production rate  
  $numtf=2; 
  $onrate=1000; ## rate parameter exponential distribution for onrate 
  $offrate=1200; ## rate parameter exponential distribution for offrate 
  $protrate=100; ## protein synthesis rate per mRNA molecule per unit time 
  $degm=0.1; ## mRNA degradation rate per mRNA molecule per unit time 
  $degp=0.02; ## protein degradation rate per protein moleculer per unit time 

  for($li=0;$li<$nmt[0]; $li++)
  {
     $numtf=$nmt[1][$li];
     for($lj=0;$lj<$cpf[0]; $lj++)
     {
	 $coopfac=$cpf[1][$lj];      

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 $modlist[$mp][8]=$mrnarate; 
	 $modlist[$mp][9]=$coopfac; 
	 $mp++;

	 $modrun++;
     }
  }
  last; 
}

use LWP::Simple;
use Parallel::ForkManager;

$pm = Parallel::ForkManager->new(20);

LINKS:

for($i=0;$i<$mp;$i++)
{
   $pm->start and next LINKS; # do the fork

   $modrun=$modlist[$i][0]; 
   $mrnarate=$modlist[$i][1];
   $onrate=$modlist[$i][2];
   $offrate=$modlist[$i][3];
   $protrate=$modlist[$i][4];
   $degm=$modlist[$i][5];
   $degp=$modlist[$i][6];
   $numtf=$modlist[$i][7];
   $coopfac=$modlist[$i][9];

   print "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp NumTF $numtf Transcr_Rate $mrnarate DivFactor $coopfac\n";
   simulate($modrun); 

   $pm->finish; # do the exit in the child process
}

$pm->wait_all_children;

print "Done...\n"; 

undef(@modlist); 



sub simulate
{

$R=Statistics::R->new();
$R->startR;

local($modno)=@_; 
## -------------------------------------------------------------------
##  Noise calculation in protein levels from 10000 cells 
## -------------------------------------------------------------------

## Single TF 

print "Noise calculation across cells ...\n"; 
print "   SingleTF binding\n"; 

$tottime=0.1500;
$numcell=10000; 

for($cn=0;$cn<$numcell;$cn++)
{
  $R->run(qq'di=rexp($numevt,rate=$onrate)');
  $ref1=$R->get('di');
  $R->run(q'rm(di)');
  @arr1=@{$ref1};

  $R->run(qq'di=rexp($numevt,rate=$offrate)');
  $ref2=$R->get('di');
  $R->run(q'rm(di)');
  @arr2=@{$ref2};

  $no=@arr1;

  $time=0;

  for($i=0;$i<$no;$i++)
  {
    $rn=0;
    $time+=$arr1[$i];
    $time=sprintf("%0.4f",$time);
    $ttlist{$time}=$mrnarate;
    $time+=$arr2[$i];
    $time=sprintf("%0.4f",$time);
    $ttlist{$time}=$bas;
    undef $rn;
  }
  
  $prev=$bas;
  $mrna=$inimrna;
  $prot=$iniprot; 

  $cl=0; 
  for($i=0;$i<$tottime;$i+=$uot)
  {
    $i=sprintf("%0.4f",$i);
    $timelist[$cl]=$i;
    if(exists $ttlist{$i})
    {
      $trlist[$cn][$cl]=$ttlist{$i}; 

      if($prev==$ttlist{$i})
      {
        $prot+=($protrate*$mrna-$degp*$prot);
        $mrna+=($mrnarate-$degm*$mrna); 
      }
      else
      {
        $prot+=(-$degp*$prot);
        $mrna+=(-$degm*$mrna); 
      }
      if($mrna<0) { $mrna=0; }
      if($prot<0) { $prot=0; }

      $mrnalist[$cn][$cl]=$mrna;
      $protlist[$cn][$cl]=$prot;

      $prev=$ttlist{$i};
    }
    else
    {
      if($prev>$bas)
      {
        $prot+=($protrate*$mrna-$degp*$prot);
        $mrna+=($mrnarate-$degm*$mrna); 
      }
      else
      {
        $prot+=(-$degp*$prot);
        $mrna+=(-$degm*$mrna); 
      }

      if($mrna<0) { $mrna=0; }
      if($prot<0) { $prot=0; }

      $trlist[$cn][$cl]=$prev; 
      $mrnalist[$cn][$cl]=$mrna;
      $protlist[$cn][$cl]=$prot;
       
    }
    $cl++;
  }
  undef $prev;
  undef $ref1;
  undef $ref2;
  undef $no;
  undef(@arr1);
  undef(@arr2);
  
  undef %ttlist; 
}


$pp=0; 
for($i=0;$i<$cl;$i++)
{
  $avgtr=0; $sdtr=0; 
  $avgm=0; $sdm=0; 
  $avgp=0; $sdp=0; 

  for($cn=0;$cn<$numcell;$cn++)
  {
    $avgtr+=$trlist[$cn][$i];
    $sdtr+=($trlist[$cn][$i])**2;
    $avgm+=$mrnalist[$cn][$i];
    $sdm+=($mrnalist[$cn][$i])**2;
    $avgp+=$protlist[$cn][$i];
    $sdp+=($protlist[$cn][$i])**2;
  }

  $sdm=(($sdm/$numcell)-(($avgm/$numcell)**2));
  if($sdm<0) { $sdm=0; }
  $sdm=sprintf("%0.4f",sqrt($sdm));
  $avgm=sprintf("%0.4f",$avgm/$numcell);

  if($avgm!=0)
  { 
     $cvm=sprintf("%0.4f",$sdm/$avgm);
  }
  else
  {
     $cvm="NA"; 
  }

  $sdp=(($sdp/$numcell)-(($avgp/$numcell)**2));
  if($sdp<0) { $sdp=0; }
  $sdp=sprintf("%0.4f",sqrt($sdp));
  $avgp=sprintf("%0.4f",$avgp/$numcell);
  if($avgp!=0)
  { 
     $cvp=sprintf("%0.4f",$sdp/$avgp);
  }
  else
  {
     $cvp="NA"; 
  }

  $sdtr=(($sdtr/$numcell)-(($avgtr/$numcell)**2));
  if($sdtr<0) { $sdtr=0; }
  $sdtr=sprintf("%0.4f",sqrt($sdtr));
  $avgtr=sprintf("%0.4f",$avgtr/$numcell);
  if($avgtr!=0) 
  {
    $cvtr=sprintf("%0.4f",$sdtr/$avgtr);
  }
  else
  {
    $cvtr="NA"; 
  }

  if($i>150) 
  {
      ## print WR "$i\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; }
      $silist[$pp][0]=$i; 
      $silist[$pp][1]=$avgtr; 
      $silist[$pp][2]=$sdtr; 
      $silist[$pp][3]=$cvtr; 
      $silist[$pp][4]=$avgm; 
      $silist[$pp][5]=$sdm; 
      $silist[$pp][6]=$cvm; 
      $silist[$pp][7]=$avgp; 
      $silist[$pp][8]=$sdp; 
      $silist[$pp][9]=$cvp; 
      $pp++;
  }

  undef $avgtr;
  undef $sdtr;
  undef $cvtr; 
  undef $avgm;
  undef $sdm;
  undef $cvm;
  undef $avgp;
  undef $sdp;
  undef $cvp;

}


for($cn=0;$cn<$numcell;$cn++)
{
  $bsize=0; $bfreq=0; 
  for($i=0;$i<$cl;$i++)
  {
    if($trlist[$cn][$i]>$bas) 
    {
      $bsize+=$trlist[$cn][$i];
      $bfreq++; 
    }
  }
  if($bfreq!=0)
  {
     $bsize=sprintf("%0.4f",$bsize/$bfreq); 
  }
  else
  {
     $bsize="NA"; 
  }
  $bfreq=sprintf("%0.4f",$bfreq/$cl); 
  $silist2[$cn][0]=$bsize; 
  $silist2[$cn][1]=$bfreq; 
  undef $bsize; 
  undef $bfreq; 
}

undef $cl;
undef(@trlist);
undef(@timelist); 
undef(@mrnalist); 
undef(@protlist); 


## Cooperative TFs
print "   CooperativeTF binding\n"; 


for($cn=0;$cn<$numcell;$cn++)
{

   for($lk=0;$lk<$numtf;$lk++)
   {
     $R->run(qq'di=rexp($numevt,rate=$onrate)');
     $ref1[$lk]=$R->get('di');
     $R->run(q'rm(di)');

     $R->run(qq'di=rexp($numevt,rate=$offrate/$coopfac)');
     $ref2[$lk]=$R->get('di');
     $R->run(q'rm(di)');
   
     $time[$lk]=0;
   }
   @arrt=@{$ref1[0]};
   $no=@arrt;
   undef(@arrt);

   for($i=0;$i<$no;$i++)
   {
     for($lk=0;$lk<$numtf;$lk++)
     {
       @arrt1=@{$ref1[$lk]}; 
       @arrt2=@{$ref2[$lk]}; 

       $time[$lk]+=$arrt1[$i];
       $time[$lk]=sprintf("%0.4f",$time[$lk]);
       $ttlist{$time[$lk]}[$lk]=1;
       $time[$lk]+=$arrt2[$i];
       $time[$lk]=sprintf("%0.4f",$time[$lk]);
       $ttlist{$time[$lk]}[$lk]=0;

       undef(@arrt1);
       undef(@arrt2); 
     }
   }

   $cl=0;
   for($lk=0;$lk<$numtf;$lk++)
   {
     $pk[$lk]=0; 
   }
   $mrna=$inimrna;
   $prot=$iniprot;
   $prev=$bas;

   for($i=0;$i<$tottime;$i+=$uot)
   {
      $i=sprintf("%0.4f",$i);
      $timelist[$cl]=$i;

      $lfl=0;
      $ccn=0;  $ccn2=0;  
      for($lk=0;$lk<$numtf;$lk++)
      {
	 if(exists $ttlist{$i}[$lk] && exists $ttlist{$i}[$lk]==1)
	 { 
	    $ccn++; 
	 }
	 if($pk[$lk]==1)
	 { 
            $ccn2++; 
	 }
      }
      if($ccn==$numtf) { $lfl=1; } 
      if($ccn2==$numtf) { $lfl=2; } 

      if($lfl==1) 
      {
           $trlist[$cn][$cl]=$mrnarate;
           if($prev==$mrnarate)
           {
              $prot+=($protrate*$mrna-$degp*$prot);
              $mrna+=($mrnarate-$degm*$mrna);
           }
           else
           {
              $prot+=(-$degp*$prot);
              $mrna+=(-$degm*$mrna);
           }
           if($mrna<0) { $mrna=0; }
           if($prot<0) { $prot=0; }
           $mrnalist[$cn][$cl]=$mrna;
           $protlist[$cn][$cl]=$prot;

	   $prev=$mrnarate;
      }
      elsif($lfl==2)
      {
           $trlist[$cn][$cl]=$mrnarate;
           if($prev==$mrnarate)
           {
              $prot+=($protrate*$mrna-$degp*$prot);
              $mrna+=($mrnarate-$degm*$mrna);
           }
           else
           {
              $prot+=(-$degp*$prot);
              $mrna+=(-$degm*$mrna);
           }
           if($mrna<0) { $mrna=0; }
           if($prot<0) { $prot=0; }

           $mrnalist[$cn][$cl]=$mrna;
           $protlist[$cn][$cl]=$prot;
           $prev=$mrnarate;
      }
      else
      {
	   $trlist[$cn][$cl]=$bas;
           $prot+=(-$degp*$prot);
           $mrna+=(-$degm*$mrna);
           if($mrna<0) { $mrna=0; }
           if($prot<0) { $prot=0; }
           $mrnalist[$cn][$cl]=$mrna;
           $protlist[$cn][$cl]=$prot;
           $prev=$bas;
      }
 
      for($lk=0;$lk<$numtf;$lk++)
      {
	 if(exists $ttlist{$i}[$lk])
	 { 
	    $pk[$lk]=$ttlist{$i}[$lk]; 
	 }
      }
      $cl++;
   }
   undef(@pk);
   undef(@ref1);
   undef(@ref2);
   undef $no; 
   undef %ttlist;
}

$pp=0; 
for($i=0;$i<$cl;$i++)
{
  $avgtr=0; $sdtr=0;
  $avgm=0; $sdm=0;
  $avgp=0; $sdp=0;

  for($cn=0;$cn<$numcell;$cn++)
  {
    $avgtr+=$trlist[$cn][$i];
    $sdtr+=($trlist[$cn][$i])**2;
    $avgm+=$mrnalist[$cn][$i];
    $sdm+=($mrnalist[$cn][$i])**2;
    $avgp+=$protlist[$cn][$i];
    $sdp+=($protlist[$cn][$i])**2;
  }

  $sdm=(($sdm/$numcell)-(($avgm/$numcell)**2));
  if($sdm<0) { $sdm=0; }
  $sdm=sprintf("%0.4f",sqrt($sdm));
  $avgm=sprintf("%0.4f",$avgm/$numcell);

  if($avgm!=0)
  {
    $cvm=sprintf("%0.4f",$sdm/$avgm);
  }
  else
  {
    $cvm="NA"; 
  }

  $sdp=(($sdp/$numcell)-(($avgp/$numcell)**2));
  if($sdp<0) { $sdp=0; }
  $sdp=sprintf("%0.4f",sqrt($sdp));
  $avgp=sprintf("%0.4f",$avgp/$numcell);
  if($avgp!=0)
  {
    $cvp=sprintf("%0.4f",$sdp/$avgp);
  }
  else
  {
    $cvp="NA"; 
  }

  $sdtr=(($sdtr/$numcell)-(($avgtr/$numcell)**2));
  if($sdtr<0) { $sdtr=0; }
  $sdtr=sprintf("%0.4f",sqrt($sdtr));
  $avgtr=sprintf("%0.4f",$avgtr/$numcell);
  if($avgtr!=0)
  {
    $cvtr=sprintf("%0.4f",$sdtr/$avgtr);
  }
  else
  {
    $cvtr="NA"; 
  }

  if($i>150) 
  { 
     ## print WR "$i\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; 
     $cmlist[$pp][0]=$i;
     $cmlist[$pp][1]=$avgtr; 
     $cmlist[$pp][2]=$sdtr; 
     $cmlist[$pp][3]=$cvtr;
     $cmlist[$pp][4]=$avgm; 
     $cmlist[$pp][5]=$sdm; 
     $cmlist[$pp][6]=$cvm; 
     $cmlist[$pp][7]=$avgp;
     $cmlist[$pp][8]=$sdp; 
     $cmlist[$pp][9]=$cvp; 
     $pp++; 
  }

  undef $avgtr;
  undef $sdtr;
  undef $cvtr;
  undef $avgm;
  undef $sdm;
  undef $cvm;
  undef $avgp;
  undef $sdp;
  undef $cvp;
}

for($cn=0;$cn<$numcell;$cn++)
{
  $bsize=0; $bfreq=0;
  for($i=0;$i<$cl;$i++)
  {
    if($trlist[$cn][$i]>$bas)
    {
      $bsize+=$trlist[$cn][$i];
      $bfreq++;
    }
  }
  if($bfreq!=0)
  {
     $bsize=sprintf("%0.4f",$bsize/$bfreq);
  }
  else
  {
     $bsize="NA"; 
  }
  $bfreq=sprintf("%0.4f",$bfreq/$cl);
  $cmlist2[$cn][0]=$bsize; 
  $cmlist2[$cn][1]=$bfreq; 
  undef $bsize;
  undef $bfreq;
}

$ff="Cooperative_TFs2/$modno"."_NoiseData.txt";
open(WR1,">$ff") or die; 

print WR1 "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp NumTF $numtf Transcr_rate $mrnarate DivFactor $coopfac\n";

print WR1 "Si_Timepoint\tSi_Avg_TR\tSi_Sd_TR\tSi_CV_TR\tSi_Avg_mRNA\tSi_Sd_mRNA\tSi_CV_mRNA\tSi_Avg_Prot\tSi_Sd_Prot\tSi_CV_Prot\t"; 
print WR1 "Cp_Timepoint\tCp_Avg_TR\tSi_Sd_TR\tCp_CV_TR\tCp_Avg_mRNA\tCp_Sd_mRNA\tCp_CV_mRNA\tCp_Avg_Prot\tCp_Sd_Prot\tCp_CV_Prot\n"; 

for($i=0;$i<$pp;$i++)
{
  for($lo=0;$lo<10;$lo++)
  {
	  print WR1 "$silist[$i][$lo]\t"; 
  }
  for($lo=0;$lo<10;$lo++)
  {
	  print WR1 "$cmlist[$i][$lo]\t"; 
  }
  print WR1 "\n"; 
}
close(WR1);


$ff="Cooperative_TFs2/$modno"."_Burst_SizeFreq.txt";
open(WR2,">$ff") or die; 

print WR2 "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp NumTF $numtf Transcr_Rate $mrnarate DivFactor $coopfac\n";

print WR2 "Si_CellID\tSi_BurstSize\tSi_BurstFreq\t"; 
print WR2 "Cp_CellID\tCp_BurstSize\tCp_BurstFreq\n"; 

for($i=0;$i<$numcell;$i++)
{
  print WR2 "$i\t$silist2[$i][0]\t$silist2[$i][1]\t$i\t$cmlist2[$i][0]\t$cmlist2[$i][1]\n"; 
}
close(WR2); 

undef(@silist);
undef(@silist2);
undef(@cmlist);
undef(@cmlist2);

undef $cl;
undef(@trlist);
undef(@timelist);
undef(@mrnalist);
undef(@protlist);

$R->stopR;

}
