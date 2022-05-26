#!/usr/bin/perl -w 

## Explore the impact of paramter values, number of TFs in competitve TF binding in 
## simulating TF competitive binding 
## Calculate impact on gene expression 

use Statistics::R;

print "Calculating and plotting time-variant dynamics ....\n"; 

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


## Parameter exploration - one at a time 
#

$mrn[0]=10; # $mrnarate 
$mrn[1][0]=10; $mrn[1][1]=20; $mrn[1][2]=50; $mrn[1][3]=100; $mrn[1][4]=200; $mrn[1][5]=500; $mrn[1][6]=1000; $mrn[1][7]=2000; $mrn[1][8]=5000; $mrn[1][9]=10000; ## 10 to 10000; 

$rt[0]=9; ## rate +-
$rt[1][0]=0.1; $rt[1][1]=0.2; $rt[1][2]=0.3; $rt[1][3]=0.4; $rt[1][4]=0.5; $rt[1][5]=0.6; $rt[1][6]=0.7; $rt[1][7]=0.8; $rt[1][8]=0.9; 

$nmt[0]=8; ## $numtf 
$nmt[1][0]=2; $nmt[1][1]=3; $nmt[1][2]=4; $nmt[1][3]=5; $nmt[1][4]=6; $nmt[1][5]=7; $nmt[1][6]=8; $nmt[1][7]=9; ## 2 to 9 
## $rates -> average rate $mrnarate 

$ont[0]=6; ## $onrate 
$ont[1][0]=100; $ont[1][1]=200; $ont[1][2]=500; $ont[1][3]=1000; $ont[1][4]=2000; $ont[1][5]=5000;  ## 100 to 5000

$oft[0]=7; ## $offrate
$oft[1][0]=100; $oft[1][1]=200; $oft[1][2]=500; $oft[1][3]=1000; $oft[1][4]=2000; $oft[1][5]=4000; $oft[1][6]=6000; ## 100 to 6000

$prt[0]=10; ## $protrate 
$prt[1][0]=10; $prt[1][1]=20; $prt[1][2]=50; $prt[1][3]=100; $prt[1][4]=200; $prt[1][5]=500; $prt[1][6]=1000; $prt[1][7]=2000; $prt[1][8]=5000; $prt[1][9]=10000; ## 10 to 10000; 

$dm[0]=8; ## $degm
$dm[1][0]=0.01; $dm[1][1]=0.02; $dm[1][2]=0.05; $dm[1][3]=0.1; $dm[1][4]=0.2; $dm[1][5]=0.3; $dm[1][6]=0.4; $dm[1][7]=0.5;  ## 0.01 to 0.5;

$dp[0]=8; ## $degp; 
$dp[1][0]=0.01; $dp[1][1]=0.02; $dp[1][2]=0.05; $dp[1][3]=0.1; $dp[1][4]=0.2; $dp[1][5]=0.3; $dp[1][6]=0.4; $dp[1][7]=0.5; ## 0.01 to 0.5; 


$runfl=0; 
$tcl=0;
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

   if($tcl==0)
   {
      for($li=0;$li<$mrn[0]; $li++)
      {
	 $mrnarate=$mrn[1][$li]; 
         $rate[0]=1.3*$mrnarate;
         $rate[1]=0.7*$mrnarate;

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

	 $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==1) 
   {
      for($li=0;$li<$rt[0]; $li++)
      {
	 $rate[0]=$mrnarate*(1+$rt[1][$li]);
	 $rate[1]=$mrnarate*(1-$rt[1][$li]);

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

	 $modrun++;
      }
      $tcl++; 
      next; 
   }
   if($tcl==2) 
   {
      for($li=0;$li<$nmt[0]; $li++)
      {
	 if($nmt[1][$li]==2) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	 }
	 elsif($nmt[1][$li]==3) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1*$mrnarate;
	 }
	 elsif($nmt[1][$li]==4) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.5*$mrnarate;
	   $rate[3]=0.5*$mrnarate;
	 }
	 elsif($nmt[1][$li]==5) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.5*$mrnarate;
	   $rate[3]=0.5*$mrnarate;
	   $rate[4]=1*$mrnarate;
	 }
	 elsif($nmt[1][$li]==6) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.5*$mrnarate;
	   $rate[3]=0.5*$mrnarate;
	   $rate[4]=1.7*$mrnarate;
	   $rate[5]=0.3*$mrnarate;
	 }
	 elsif($nmt[1][$li]==7) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.5*$mrnarate;
	   $rate[3]=0.5*$mrnarate;
	   $rate[4]=1.7*$mrnarate;
	   $rate[5]=0.3*$mrnarate;
	   $rate[6]=1*$mrnarate;
	 }
	 elsif($nmt[1][$li]==8) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.5*$mrnarate;
	   $rate[3]=0.5*$mrnarate;
	   $rate[4]=1.7*$mrnarate;
	   $rate[5]=0.3*$mrnarate;
	   $rate[6]=1.9*$mrnarate;
	   $rate[7]=0.1*$mrnarate;
	 }
	 elsif($nmt[1][$li]==9) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.5*$mrnarate;
	   $rate[3]=0.5*$mrnarate;
	   $rate[4]=1.7*$mrnarate;
	   $rate[5]=0.3*$mrnarate;
	   $rate[6]=1.9*$mrnarate;
	   $rate[7]=0.1*$mrnarate;
	   $rate[8]=1*$mrnarate;
	 }
	 $numtf=$nmt[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

	 $modrun++;
      }
      $tcl++; 
      next; 
   }
   if($tcl==3) 
   {
      for($li=0;$li<$nmt[0]; $li++)
      {
	 if($nmt[1][$li]==2) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	 }
	 elsif($nmt[1][$li]==3) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1*$mrnarate;
	 }
	 elsif($nmt[1][$li]==4) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.3*$mrnarate;
	   $rate[3]=0.7*$mrnarate;
	 }
	 elsif($nmt[1][$li]==5) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.3*$mrnarate;
	   $rate[3]=0.7*$mrnarate;
	   $rate[4]=1*$mrnarate;
	 }
	 elsif($nmt[1][$li]==6) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.3*$mrnarate;
	   $rate[3]=0.7*$mrnarate;
	   $rate[4]=1.3*$mrnarate;
	   $rate[5]=0.7*$mrnarate;
	 }
	 elsif($nmt[1][$li]==7) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.3*$mrnarate;
	   $rate[3]=0.7*$mrnarate;
	   $rate[4]=1.3*$mrnarate;
	   $rate[5]=0.7*$mrnarate;
	   $rate[6]=1*$mrnarate;
	 }
	 elsif($nmt[1][$li]==8) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.3*$mrnarate;
	   $rate[3]=0.7*$mrnarate;
	   $rate[4]=1.3*$mrnarate;
	   $rate[5]=0.7*$mrnarate;
	   $rate[6]=1.3*$mrnarate;
	   $rate[7]=0.7*$mrnarate;
	 }
	 elsif($nmt[1][$li]==9) 
	 {
	   $rate[0]=1.3*$mrnarate;
	   $rate[1]=0.7*$mrnarate;
	   $rate[2]=1.3*$mrnarate;
	   $rate[3]=0.7*$mrnarate;
	   $rate[4]=1.3*$mrnarate;
	   $rate[5]=0.7*$mrnarate;
	   $rate[6]=1.3*$mrnarate;
	   $rate[7]=0.7*$mrnarate;
	   $rate[8]=1*$mrnarate;
	 }
	 $numtf=$nmt[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

	 $modrun++;
      }
      $tcl++; 
      next; 
   }
   if($tcl==4) 
   {
      for($li=0;$li<$ont[0]; $li++)
      {
         $onrate=$ont[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==5) 
   {
      for($li=0;$li<$oft[0]; $li++)
      {
         $offrate=$oft[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==6) 
   {
      for($li=0;$li<$prt[0]; $li++)
      {
         $protrate=$prt[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==7) 
   {
      for($li=0;$li<$dm[0]; $li++)
      {
         $degm=$dm[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==8) 
   {
      for($li=0;$li<$dp[0]; $li++)
      {
         $degp=$dp[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $modlist[$mp][7]=$numtf; 
	 for($ri=0;$ri<$numtf;$ri++)
	 {
	   $modlist[$mp][8][$ri]=$rate[$ri]; 
	 }
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      last; 
   } 
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

   print "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp NumTF $numtf \t";

   for($ri=0;$ri<$numtf;$ri++)
   {
      $rate[$ri]=$modlist[$i][8][$ri];
      print "$rate[$ri]\t";
   }
   print "\n"; 

   simulate($modrun); 
   undef(@rate); 

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
##  Noise calculation in protein levels from 1000 cells 
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


## Competitive TFs
print "   CompetitiveTF binding\n"; 

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
    $rn=int(rand($numtf));
    $time+=$arr1[$i];
    $time=sprintf("%0.4f",$time);
    $ttlist{$time}=$rate[$rn];
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
        $mrna+=($prev-$degm*$mrna);
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
      $trlist[$cn][$cl]=$prev; 

      if($prev>$bas)
      {
        $prot+=($protrate*$mrna-$degp*$prot);
        $mrna+=($prev-$degm*$mrna);
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

$ff="Competitive_TFs/$modno"."_NoiseData.txt";
open(WR1,">$ff") or die; 

print WR1 "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp NumTF $numtf\t";
for($k=0;$k<$numtf;$k++)
{
  print WR1 "$rate[$k]\t";
}
print WR1 "\n";

print WR1 "Si_Timepoint\tSi_Avg_TR\tSi_Sd_TR\tSi_CV_TR\tSi_Avg_mRNA\tSi_Sd_mRNA\tSi_CV_mRNA\tSi_Avg_Prot\tSi_Sd_Prot\tSi_CV_Prot\t"; 
print WR1 "Cm_Timepoint\tCm_Avg_TR\tSi_Sd_TR\tCm_CV_TR\tCm_Avg_mRNA\tCm_Sd_mRNA\tCm_CV_mRNA\tCm_Avg_Prot\tCm_Sd_Prot\tCm_CV_Prot\n"; 

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


$ff="Competitive_TFs/$modno"."_Burst_SizeFreq.txt";
open(WR2,">$ff") or die; 

print WR2 "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp NumTF $numtf\t";
for($k=0;$k<$numtf;$k++)
{
   print WR2 "$rate[$k]\t";
}
print WR2 "\n";

print WR2 "Si_CellID\tSi_BurstSize\tSi_BurstFreq\t"; 
print WR2 "Cm_CellID\tCm_BurstSize\tCm_BurstFreq\n"; 

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
