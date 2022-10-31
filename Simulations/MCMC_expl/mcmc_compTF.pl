#!/usr/bin/perl -w 

# Explore parameter space using MCMC sampling for single TF binding
# 
# Mean expression target 115000-120000
# 

$lowlim=1000000; 
$uplim=1500000; 

$testrun=1; 

$mmrun=10000; 

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

$mrn[0]=10; # $mrnarate step  
$mrn[1]=10; ## Minimum 
$mrn[2]=1000; ## Maximum, range - 10 to 10000; 

$nmt[0]=1; ## $numtf step 
$nmt[1]=2; ## Min 
$nmt[2]=9; ## Max 
## $rates -> average rate $mrnarate 

$rt[0]=0.1; ## rate +- Step for multi-TF regulation 
$rt[1]=0.1; ## Min 
$rt[2]=0.9; ## Max 

$ont[0]=100; ## $onrate step 
$ont[1]=100; ## Min 
$ont[2]=5000; ## Max 

$oft[0]=100; ## $offrate step 
$oft[1]=100; ## Min 
$oft[2]=6000; ## Max 

$prt[0]=10; ## $protrate step 
$prt[1]=10; ## Min 
$prt[2]=1000; ## Max  

$dm[0]=0.01; ## $degm step
$dm[1]=0.01; ## Min 
$dm[2]=0.25;  ## Max 

$dp[0]=0.01; ## $degp step 
$dp[1]=0.01; ## Min 
$dp[2]=0.25;  ## Max 

print "Noise calculation across cells ...\n"; 
print "   CompetitiveTF binding\n"; 

$runlim=50; 

use LWP::Simple;
use Parallel::ForkManager;

$pm = Parallel::ForkManager->new(20);

LINKS:

for($rnd=1;$rnd<=$mmrun;$rnd++)
{

## Initiate with a random set of parameter values 

$retexp=0; 

while($retexp<$lowlim/5 || $retexp>$uplim*5)
{
  $mrnarate=$mrn[1]+(int(rand(($mrn[2]-$mrn[1])/$mrn[0]-1)))*$mrn[0]; ## mRNA production rate  
  $numtf=$nmt[1]+(int(rand(($nmt[2]-$nmt[1])/$nmt[0]-1)))*$nmt[0]; 
  $onrate=$ont[1]+(int(rand(($ont[2]-$ont[1])/$ont[0]-1)))*$ont[0]; ## rate parameter exponential distribution for onrate 
  $offrate=$oft[1]+(int(rand(($oft[2]-$oft[1])/$oft[0]-1)))*$oft[0]; ## rate parameter exponential distribution for offrate 
  $protrate=$prt[1]+(int(rand(($prt[2]-$prt[1])/$prt[0]-1)))*$prt[0]; ## protein synthesis rate per mRNA molecule per unit time 
  $degm=$dm[1]+(int(rand(($dm[2]-$dm[1])/$dm[0]-1)))*$dm[0]; ## mRNA degradation rate per mRNA molecule per unit time 
  $degp=$dp[1]+(int(rand(($dp[2]-$dp[1])/$dp[0]-1)))*$dp[0]; ## protein degradation rate per protein moleculer per unit time 

  for($ri=0;$ri<$numtf;$ri++)
  {
    $rde=$rt[1]+rand(($rt[2]-$rt[1])); 
    if($rde<$rt[1]) { $rde=$rt[1]; }
    if($rde>$rt[2]) { $rde=$rt[2]; }
    $sig=rand(1); 
    if($sig<=0.5) { $sign=1; } 
    else { $sign=-1; } 
    $rt=(1+$sign*$rde); 
    $rate[$ri]=sprintf("%0.4f",$rt); 
    #print "$sign $rde $rate[$ri]\n"; 
    undef $rt;
    undef $sign; 
    undef $sig; 
    undef $rde; 
  }
  $retexp=simulate(-1,10,-1); 
}

$pm->start and next LINKS; # do the fork

$runfl=0;
$prevexp=0; 

while($runfl<=$runlim)
{
  $meanexp=simulate($testrun,1000,$rnd);
  print "$rnd mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp NumTF $numtf  MEANEXP: $meanexp\n";
  #print "MEANEXP: $meanexp\n"; 

  if($meanexp<$lowlim) 
  {
     if($lowlim-$meanexp<$lowlim-$prevexp && $runfl!=0) 
     {
	if($prevchng[1]==0)
	{
          $prevchng=randomchange(); ## Randomly pick another parameter and change in random direction 
	  @qw=split(/\*/,$prevchng); 
          $prevchng[0]=$qw[0]; ## Parameter  
          $prevchng[1]=$qw[1]; ## Change direction +ve or -ve
	  $prevchng[2]=$qw[2]; ## Fold-change
	  undef(@qw); 
	  undef $prevchng; 
	  #print "RANDOMCH\n";
	}
	else
	{
	  $prevchng=targetedchange($prevchng[0],$prevchng[1]); 
	  @qw=split(/\*/,$prevchng); 
          $prevchng[0]=$qw[0]; ## Parameter  
          $prevchng[1]=$qw[1]; ## Change direction +ve or -ve
	  $prevchng[2]=$qw[2]; ## Fold-change
	  undef(@qw); 
	  undef $prevchng; 
	  #print "TARGETEDCH\n";
	}
     }
     else
     {
	if($runfl!=0 && $meanexp<$lowlim/2)
	{
	   restoreprevparam($prevchng[0],$prevchng[1],$prevchng[2]); 
	   $meanexp=$prevexp;
	}
        $prevchng=randomchange(); ## Randomly pick another parameter and change in random direction 
	@qw=split(/\*/,$prevchng); 
        $prevchng[0]=$qw[0]; ## Parameter  
        $prevchng[1]=$qw[1]; ## Change direction +ve or -ve
	$prevchng[2]=$qw[2]; ## Fold-change
	undef(@qw); 
	undef $prevchng; 
	#print "RANDOMCH\n";
     } 
  }
  elsif($meanexp>$uplim) 
  {
     if($meanexp-$uplim<$prevexp-$uplim && $runfl!=0) 
     {
	  $prevchng=targetedchange($prevchng[0],$prevchng[1]); 
	  @qw=split(/\*/,$prevchng); 
          $prevchng[0]=$qw[0]; ## Parameter  
          $prevchng[1]=$qw[1]; ## Change direction +ve or -ve
	  $prevchng[2]=$qw[2]; ## Fold-change
	  undef(@qw); 
	  undef $prevchng; 
	  #print "TARGETEDCH\n";
     }
     else
     {
	if($runfl!=0 && $meanexp>$uplim*2)
	{
	   restoreprevparam($prevchng[0],$prevchng[1],$prevchng[2]); 
	   $meanexp=$prevexp;
	}
        $prevchng=randomchange(); ## Randomly pick another parameter and change in random direction 
	@qw=split(/\*/,$prevchng); 
        $prevchng[0]=$qw[0]; ## Parameter  
        $prevchng[1]=$qw[1]; ## Change direction +ve or -ve
	$prevchng[2]=$qw[2]; ## Fold-change
	undef(@qw); 
	undef $prevchng; 
	#print "RANDOMCH\n";
     }
  }
  else ## $eanexp>=$lowlim && $meanexp<=$uplim
  {
     last; 
  }
  $prevexp=$meanexp; 
  undef $meanexp; 
  $runfl++; 
}

undef(@rate); 
undef(@prevchng); 

 $pm->finish; # do the exit in the child process
}

$pm->wait_all_children;
print "Done...\n"; 


$ff="Comp/$testrun"."_NoiseData.txt";
open(WR,">$ff") or die; 
print WR "mRNArate\tOnrate\tOffrate\tProtrate\tDegmRNA\tDegProt\tNumTF\t\tSi_Timepoint\tSi_Avg_TR\tSi_Sd_TR\tSi_CV_TR\tSi_Avg_mRNA\tSi_Sd_mRNA\tSi_CV_mRNA\tSi_Avg_Prot\tSi_Sd_Prot\tSi_CV_Prot\n"; 

system "ls Comp/TMP_*_Noise.txt > TMPComp.txt"; 
open(TM,"TMPComp.txt") or die; 
while($tm=<TM>)
{
   chomp($tm); 
   open(FP,$tm) or die; 
   while($fp=<FP>)
   {
     chomp($fp); 
     print WR "$fp\n";
   }
   close(FP); 
}
close(TM); 
close(WR); 
system "rm Comp/TMP_*_Noise.txt"; 

$ff2="Comp/$testrun"."_Burst_SizeFreq.txt";
open(WR,">$ff2") or die; 
print WR "mRNArate\tOnrate\tOffrate\tProtrate\tDegmRNA\tDegProt\tNumTF\t\tSi_BurstSize\tSi_BurstFreq\n"; 
system "ls Comp/TMP_*_Burst_SizeFreq.txt > TMPComp.txt"; 
open(TM,"TMPComp.txt") or die; 
while($tm=<TM>)
{
   chomp($tm); 
   open(FP,$tm) or die; 
   while($fp=<FP>)
   {
     chomp($fp); 
     print WR "$fp\n";
   }
   close(FP); 
}
close(TM); 
close(WR); 
system "rm Comp/TMP_*_Burst_SizeFreq.txt"; 
close(WR); 
undef $ff; 
undef $ff2; 


sub randomchange
{
  $var=int(rand(8)); 
  $dr=rand(1); 

  if($dr<0.5) { $dirch=1; }
  else { $dirch=-1; } 

  $lfl=1+int(rand(9)); 

        if($var==0) 
        {
	   $mrnarate+=($dirch*$mrn[0]*$lfl); 
	   if($mrnarate<$mrn[1]) { $mrnarate=$mrn[1]; $dirch=0; }
	   if($mrnarate>$mrn[2]) { $mrnarate=$mrn[2]; $dirch=0; }
	}
	if($var==1) 
        {
	   $lprn=$numtf; 
	   $lfl=1; 
	   $numtf+=($dirch*$nmt[0]); 
	   if($numtf<$nmt[1]) { $numtf=$nmt[1]; $dirch=0; }
	   if($numtf>$nmt[2]) { $numtf=$nmt[2]; $dirch=0; }
           for($ri=$lprn;$ri<$numtf;$ri++)
           {
              $rde=$rt[1]+rand(($rt[2]-$rt[1])); 
              if($rde<$rt[1]) { $rde=$rt[1]; $dirch=0; }
              if($rde>$rt[2]) { $rde=$rt[2]; $dirch=0; }
              $sig=rand(1); 
              if($sig<=0.5) { $sign=1; } 
              else { $sign=-1; } 
              $rt=(1+$sign*$rde); 
              #print "$sign $rde $rate[$ri]\n"; 
              $rate[$ri]=sprintf("%0.4f",$rt); 
              undef $rt;
              undef $sign; 
              undef $sig; 
              undef $rde; 
           }
	   undef $lprn; 
	}
	if($var==2) 
	{
	   $tfid=int(rand($numtf-1)); 
	   $rate[$tfid]+=($dirch*$rt[0]); 
	   if($rate[$tfid]<$rt[1]) { $rate[$tfid]=$rt[1]; $dirch=0; }
	   if($rate[$tfid]>$rt[2]) { $rate[$tfid]=$rt[2]; $dirch=0; }
	}
	if($var==3)
	{
	   $onrate+=($dirch*$ont[0]*$lfl); 
	   if($onrate<$ont[1]) { $onrate=$ont[1]; $dirch=0; }
	   if($onrate>$ont[2]) { $onrate=$ont[2]; $dirch=0; }
	}
	if($var==4)
	{
	   $offrate+=($dirch*$oft[0]*$lfl); 
	   if($offrate<$oft[1]) { $offrate=$oft[1]; $dirch=0; }
	   if($offrate>$oft[2]) { $offrate=$oft[2]; $dirch=0; }
	}
	if($var==5)
	{
	   $protrate+=($dirch*$prt[0]*$lfl); 
	   if($protrate<$prt[1]) { $protrate=$prt[1]; $dirch=0; }
	   if($protrate>$prt[2]) { $protrate=$prt[2]; $dirch=0; }
	}
	if($var==6)
	{
	   $degm+=($dirch*$dm[0]*$lfl); 
	   if($degm<$dm[1]) { $degm=$dm[1]; $dirch=0; }
	   if($degm>$dm[2]) { $degm=$dm[2]; $dirch=0; }
	}
	if($var==7)
	{
	   $degp+=($dirch*$dp[0]*$lfl); 
	   if($degp<$dp[1]) { $degp=$dp[1]; $dirch=0; }
	   if($degp>$dp[2]) { $degp=$dp[2]; $dirch=0; }
	}

        $lprev0=$var; ## Parameter  
	if($var==2) 
	{
	   $lprev0.=$tfid; 
	   undef $tfid; 
	}
        $lprev1=$dirch; ## Change direction +ve or -ve
	undef $tfid; 
	undef $var; 
	undef $dr; 
	undef $dirch; 

    $retval=$lprev0."*$lprev1*$lfl"; 
    return $retval; 
}

sub targetedchange 
{
   local($lprevchng0,$ldr)=@_; 

  $lfl=1+int(rand(9)); 

          if($lprevchng0==0) 
	  {
	     $mrnarate+=($ldr*$mrn[0]*$lfl); 
	     if($mrnarate<$mrn[1]) { $mrnarate=$mrn[1]; $ldr=0; }
	     if($mrnarate>$mrn[2]) { $mrnarate=$mrn[2]; $ldr=0; }
	  }
	  if($lprevchng0==1)
	  {
	     $prn=$numtf; 
	     $lfl=1; 
	     $numtf+=($ldr*$nmt[0]); 
	     if($numtf<$nmt[1]) { $numtf=$nmt[1]; $ldr=0; }
	     if($numtf>$nmt[2]) { $numtf=$nmt[2]; $ldr=0; }
             for($ri=$prn;$ri<$numtf;$ri++)
             {
               $rde=$rt[1]+rand(($rt[2]-$rt[1])); 
               if($rde<$rt[1]) { $rde=$rt[1]; $ldr=0; }
               if($rde>$rt[2]) { $rde=$rt[2]; $ldr=0; }
               $sig=rand(1); 
               if($sig<=0.5) { $sign=1; } 
               else { $sign=-1; } 
               $rt=(1+$sign*$rde); 
               #print "$sign $rde $rate[$ri]\n"; 
               $rate[$ri]=sprintf("%0.4f",$rt); 
               undef $rt;
               undef $sign; 
               undef $sig; 
               undef $rde; 
             }
	     undef $prn; 
	  }
	  if($lprevchng0>=20)
	  {
	     $ls=substr($lprevchng0,1,1); 
	     $rate[$ls]+=($ldr*$rt[0]); 
	     if($rate[$ls]<$rt[1]) { $rate[$ls]=$rt[1]; $ldr=0; } 
	     if($rate[$ls]>$rt[2]) { $rate[$ls]=$rt[2]; $ldr=0; } 
	     undef $ls; 
          }
          if($lprevchng0==3) 
	  {
	     $onrate+=($ldr*$ont[0]*$lfl); 
	     if($onrate<$ont[1]) { $onrate=$ont[1]; $ldr=0; }
	     if($onrate>$ont[2]) { $onrate=$ont[2]; $ldr=0; }
	  }
          if($lprevchng0==4) 
	  {
	     $offrate+=($ldr*$oft[0]*$lfl); 
	     if($offrate<$oft[1]) { $offrate=$oft[1]; $ldr=0; }
	     if($offrate>$oft[2]) { $offrate=$oft[2]; $ldr=0; }
	  }
          if($lprevchng0==5) 
	  {
	     $protrate+=($ldr*$prt[0]*$lfl); 
	     if($protrate<$prt[1]) { $protrate=$prt[1]; $ldr=0; }
	     if($protrate>$prt[2]) { $protrate=$prt[2]; $ldr=0; }
	  }
          if($lprevchng0==6) 
	  {
	     $degm+=($ldr*$dm[0]*$lfl); 
	     if($degm<$dm[1]) { $degm=$dm[1]; $ldr=0; }
	     if($degm>$dm[2]) { $degm=$dm[2]; $ldr=0; }
	  }
          if($lprevchng0==7) 
	  {
	     $degp+=($ldr*$dp[0]*$lfl); 
	     if($degp<$dp[1]) { $degp=$dp[1]; $ldr=0; }
	     if($degp>$dp[2]) { $degp=$dp[2]; $ldr=0; }
	  }

	  $lretval=$lprevchng0."*$ldr*$lfl"; 
	  return $lretval;
}

sub restoreprevparam
{
   local($lprev0,$lprev1,$lfld)=@_; 

          if($lprev0==0) 
	  {
	     $mrnarate-=($lprev1*$mrn[0]*$lfld); 
	  }
	  if($lprev0==1)
	  {
	     $prn=$numtf; 
	     $numtf-=($lprev1*$nmt[0]); 
             for($ri=$prn;$ri<$numtf;$ri++)
             {
               $rate[$ri]=1; 
             }
	     undef $prn; 
	  }
	  if($lprev0>=20)
	  {
	     $ls=substr($lprev0,1,1); 
	     $rate[$ls]-=($lprev1*$rt[0]); 
	     undef $ls; 
          }
          if($lprev0==3) 
	  {
	     $onrate-=($lprev1*$ont[0]*$lfld); 
	  }
          if($lprev0==4) 
	  {
	     $offrate-=($lprev1*$oft[0]*$lfld); 
	  }
          if($lprev0==5) 
	  {
	     $protrate-=($lprev1*$prt[0]*$lfld); 
	  }
          if($lprev0==6) 
	  {
	     $degm-=($lprev1*$dm[0]*$lfld); 
	  }
          if($lprev0==7) 
	  {
	     $degp-=($lprev1*$dp[0]*$lfld); 
	  }

	  return;
}


sub simulate
{

$R=Statistics::R->new();
$R->start_sharedR;

local($modno,$lcell,$lrn)=@_; 
local($ltmeanexp); 

#
## -------------------------------------------------------------------
##  Noise calculation in protein levels from 1000 cells 
## -------------------------------------------------------------------

## Competitive TF 


$tottime=0.1000; ## 0.1500
$numcell=$lcell;  ## 10000 


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
    $ttlist{$time}=$rate[$rn]*$mrnarate;
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
for($i=0;$i<205;$i++) ## $cl
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

  if($i==204) 
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
      $ltmeanexp=$avgp; 
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

if($modno==-1)
{
  undef(@silist);
  undef $cl;
  undef(@trlist);
  undef(@timelist);
  undef(@mrnalist);
  undef(@protlist);
  return $ltmeanexp;
}

$avbsize=0; $avbfreq=0; 
$lpr=0; 
for($cn=0;$cn<$numcell;$cn++)
{
  $bsize=0; $bfreq=0; 
  for($i=0;$i<205;$i++)
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
     $avbsize+=$bsize; 
     $lpr++; 
  }
  else
  {
     $bsize="NA"; 
  }
  $bfreq=sprintf("%0.4f",$bfreq/205); 
  $avbfreq+=$bfreq; 
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


$ff="Comp/TMP"."_$lrn"."_Noise.txt";

open(WR1,">$ff") or die; 

$wrfl=0;
for($i=0;$i<$pp;$i++)
{
  if($silist[$i][7]>=$lowlim && $silist[$i][7]<=$uplim)
  {
     print WR1 "$mrnarate\t$onrate\t$offrate\t$protrate\t$degm\t$degp\t$numtf\t";
     for($k=0;$k<$numtf;$k++)
     {
       print WR1 "$rate[$k] ";
     }
     print WR1 "\t\t";

     for($lo=0;$lo<10;$lo++)
     {
	  print WR1 "$silist[$i][$lo]\t"; 
     }
     print WR1 "\n"; 
     $wrfl=1;
  }
}
close(WR1);


$ff="Comp/TMP"."_$lrn"."_Burst_SizeFreq.txt";

open(WR2,">$ff") or die;

$avbsize=sprintf("%0.4f",$avbsize/$lpr);
$avbfreq=sprintf("%0.4f",$avbfreq/$lpr);

if($wrfl==1)
{
	#for($i=0;$i<$numcell;$i++)
	#{
    print WR2 "$mrnarate\t$onrate\t$offrate\t$protrate\t$degm\t$degp\t$numtf\t";
    for($k=0;$k<$numtf;$k++)
    {
      print WR2 "$rate[$k] ";
    }
    print WR2 "\t\t$avbsize\t$avbfreq\n";
    #print WR2 "$i\t$silist2[$i][0]\t$silist2[$i][1]\n"; 
    #}
}

close(WR2); 

undef(@silist);
undef(@silist2);

$R->stopR;

return $ltmeanexp; 

}
