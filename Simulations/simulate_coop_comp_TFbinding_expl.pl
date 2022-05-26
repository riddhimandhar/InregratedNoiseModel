#!/usr/bin/perl -w 

## Simulate TF co-operativity and competitive binding and their impact on gene expression 
## Explore appropriate on-rate and off-rate ratio and plot
#

print "Calculating and plotting time-variant dynamics ....\n"; 

use LWP::Simple;
use Parallel::ForkManager;

$numevt=500; 
$bas=20; 
$mrnarate=100; 
$rate[0]=$mrnarate;
$rate[1]=1.3*$mrnarate;
$rate[2]=0.7*$mrnarate;

$onrate=1000;
$offrate=1200; 

$tottime=0.0500; 
$uot=0.0001;

$degm=0.1; ## per mRNA molecule per unit time 

$protrate=100; ## per mRNA molecule per unit time 

$degp=0.02; ## per protein moleculer per unit time 

$inimrna=10; 
$iniprot=1000; 

$coopfac1=1; 
$coopfac2=1; 

$oncc=14; 
$ont[0]=10; $ont[1]=50; $ont[2]=100; $ont[3]=200; $ont[4]=500; $ont[5]=1000; $ont[6]=1500; $ont[7]=2000; $ont[8]=2500; $ont[9]=3000; $ont[10]=3500; 
$ont[11]=4000; $ont[12]=4500; $ont[13]=5000; 

$offcc=16; 
$oft[0]=10; $oft[1]=50; $oft[2]=100; $oft[3]=200; $oft[4]=500; $oft[5]=1000; $oft[6]=1500; $oft[7]=2000; $oft[8]=2500; $oft[9]=3000; $oft[10]=3500; 
$oft[11]=4000; $oft[12]=4500; $oft[13]=5000; $oft[14]=5500; $oft[15]=6000; 

$u=0; 

for($kon=0;$kon<$oncc;$kon++) 
{
  for($koff=0; $koff<$offcc;$koff++)
  {
    $val[$u][0]=$ont[$kon]; 
    $val[$u][1]=$oft[$koff]; 
    $u++;
  }
}


$pm = Parallel::ForkManager->new(20);

LINKS:

for($k=0;$k<$u;$k++) 
{
   $pm->start and next LINKS; # do the fork

   use Statistics::R;
   $R=Statistics::R->new();
   $R->startR;

   $onrate=$val[$k][0]; 
   $offrate=$val[$k][1]; 

   print "On-rate: $onrate Off-rate: $offrate\n"; 

## Time varriations
#
$R->run(q'set.seed(1234)');

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
  $ttlist{$time}=$rate[$rn];
  $time+=$arr2[$i];
  $time=sprintf("%0.4f",$time);
  $ttlist{$time}=$bas;
  undef $rn;
}

print "   SingleTF\t$no\t$time\n";

$prev=$bas; 
$mrna=$inimrna;
$prot=$iniprot; 

$tm1="SingleTF_expl_$onrate"."_$offrate.txt";
$tm2="SingleTF_mRNA_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;
print TM "Time\tTrrate\n";
print TM2 "Time\tmRNA\tprot\n";
for($i=0;$i<$tottime;$i+=$uot)
{
  $i=sprintf("%0.4f",$i);
  if(exists $ttlist{$i})
  {
    print TM "$i\t$ttlist{$i}\n"; 
  
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

    print TM2 "$i\t$mrna\t$prot\n"; 

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

    print TM2 "$i\t$mrna\t$prot\n"; 
    print TM "$i\t$prev\n"; 
  }
}
close(TM);
close(TM2); 

undef $prev; 
undef $ref1;
undef $ref2; 
undef $no;
undef(@arr1);
undef(@arr2);

$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');


$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Trrate) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/SingleTF_Timevar_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()+ylim(0,180)');
$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm1";
undef %ttlist; 

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, mRNA) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/SingleTF_Timevar_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'df <- tab %>% select(Time, prot) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/SingleTF_Timevar_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'rm(tab)');
system "rm $tm2";

undef $tm1;
undef $tm2;

#Multi TF 
#
$R->run(q'set.seed(1234)');

$R->run(qq'di=rexp($numevt,rate=$onrate)');
$ref1=$R->get('di');
$R->run(q'rm(di)');
@arr1=@{$ref1};


$R->run(qq'di=rexp($numevt,rate=$offrate)');
$ref2=$R->get('di');
$R->run(q'rm(di)');
@arr2=@{$ref2};

$R->run(qq'di=rexp($numevt,rate=$onrate)');
$ref3=$R->get('di');
$R->run(q'rm(di)');
@arr3=@{$ref3};

$R->run(qq'di=rexp($numevt,rate=$offrate)');
$ref4=$R->get('di');
$R->run(q'rm(di)');
@arr4=@{$ref4};

$no=@arr1;

$time1=0;
$time2=0; 

for($i=0;$i<$no;$i++)
{
  $time1+=$arr1[$i];
  $time1=sprintf("%0.4f",$time1);
  $ttlist1{$time1}=$rate[1];
  $time1+=$arr2[$i];
  $time1=sprintf("%0.4f",$time1);
  $ttlist1{$time1}=$bas;

  $time2+=$arr3[$i];
  $time2=sprintf("%0.4f",$time2);
  $ttlist2{$time2}=$rate[2];
  $time2+=$arr4[$i];
  $time2=sprintf("%0.4f",$time2);
  $ttlist2{$time2}=$bas;
}

if($time1>$time2) { $time=$time1; } 
else { $time=$time2; }

print "   MultiTF\t$no\t$time\n";

$prev=$bas; 
$mrna=$inimrna;
$prot=$iniprot; 

$tm1="MultiTF_expl_$onrate"."_$offrate.txt";
$tm2="MultiTF_mRNA_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;
print TM "Time\tTrrate\n";
print TM2 "Time\tmRNA\tprot\n";
for($i=0;$i<$tottime;$i+=$uot)
{
  $i=sprintf("%0.4f",$i);
  if(exists $ttlist1{$i} && exists $ttlist2{$i})
  {
    $rt=($ttlist1{$i}+$ttlist2{$i})*0.77;
    print TM "$i\t$rt\n"; 
  
    if($prev>2*$bas)
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

    print TM2 "$i\t$mrna\t$prot\n"; 

    $prev=$rt;
  }
  elsif(exists $ttlist1{$i})
  {
    $rt=$ttlist1{$i}*0.77;
    print TM "$i\t$rt\n"; 
  
    if($prev>2*$bas)
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

    print TM2 "$i\t$mrna\t$prot\n"; 

    $prev=$rt;
  }
  elsif(exists $ttlist2{$i})
  {
    $rt=$ttlist2{$i}*0.77;
    print TM "$i\t$rt\n"; 
  
    if($prev>2*$bas)
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

    print TM2 "$i\t$mrna\t$prot\n"; 

    $prev=$rt;
  }
  else
  {
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

    print TM2 "$i\t$mrna\t$prot\n"; 
    print TM "$i\t$prev\n"; 
  }
}
close(TM);
close(TM2); 

undef $prev; 
undef $ref1;
undef $ref2; 
undef $ref3; 
undef $ref4; 
undef $no;
undef(@arr1);
undef(@arr2);
undef(@arr3);
undef(@arr4);

$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');


$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Trrate) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/MultiTF_Timevar_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()+ylim(0,180)');
$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm1";
undef %ttlist1; 
undef %ttlist2; 

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, mRNA) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/MultiTF_Timevar_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'df <- tab %>% select(Time, prot) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/MultiTF_Timevar_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'rm(tab)');
system "rm $tm2";

undef $tm1;
undef $tm2;

## Cooperative TF 
$R->run(q'set.seed(1234)');

$R->run(qq'di=rexp($numevt,rate=$onrate*$coopfac1)');
$ref1=$R->get('di');
$R->run(q'rm(di)');
@arr1=@{$ref1};

$R->run(qq'di=rexp($numevt,rate=$offrate/$coopfac2)');
$ref2=$R->get('di');
$R->run(q'rm(di)');
@arr2=@{$ref2};

$R->run(qq'di=rexp($numevt,rate=$onrate*$coopfac1)');
$ref3=$R->get('di');
$R->run(q'rm(di)');
@arr3=@{$ref3};

$R->run(qq'di=rexp($numevt,rate=$offrate/$coopfac2)');
$ref4=$R->get('di');
$R->run(q'rm(di)');
@arr4=@{$ref4};

$no=@arr1;

$time1=0;
$time2=0; 

for($i=0;$i<$no;$i++)
{
  $time1+=$arr1[$i];
  $time1=sprintf("%0.4f",$time1);
  $ttlist1{$time1}=1;
  $time1+=$arr2[$i];
  $time1=sprintf("%0.4f",$time1);
  $ttlist1{$time1}=0;

  $time2+=$arr3[$i];
  $time2=sprintf("%0.4f",$time2);
  $ttlist2{$time2}=1;
  $time2+=$arr4[$i];
  $time2=sprintf("%0.4f",$time2);
  $ttlist2{$time2}=0;
}

if($time1>$time2) { $time=$time1; } 
else { $time=$time2; }

print "   CooperativeTF\t$no\t$time\n";

$p1=0; $p2=0; 
$mrna=$inimrna;
$prot=$iniprot;
$prev=$bas;

$tm1="CooperativeTF_expl_$onrate"."_$offrate.txt";
$tm2="CooperativeTF_mRNA_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;
print TM "Time\tTrrate\n";
print TM2 "Time\tmRNA\tprot\n";

for($i=0;$i<$tottime;$i+=$uot)
{
  $i=sprintf("%0.4f",$i);
  if(exists $ttlist1{$i} && exists $ttlist2{$i})
  {
    if($ttlist1{$i}==1 && $ttlist2{$i}==1)
    {
      print TM "$i\t$rate[0]\n";

      if($prev==$rate[0])
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

      print TM2 "$i\t$mrna\t$prot\n"; 

      $prev=$rate[0];
    }
    else
    {
      print TM "$i\t$bas\n"; 

      $prot+=(-$degp*$prot);
      $mrna+=(-$degm*$mrna);
      if($mrna<0) { $mrna=0; }
      if($prot<0) { $prot=0; }

      print TM2 "$i\t$mrna\t$prot\n"; 

      $prev=$bas; 
    }
    
  }
  else
  { 
    if($p1==1 && $p2==1) 
    {
      print TM "$i\t$rate[0]\n";

      if($prev==$rate[0])
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

      print TM2 "$i\t$mrna\t$prot\n"; 

      $prev=$rate[0];
    }
    else
    {
      print TM "$i\t$bas\n";

      $prot+=(-$degp*$prot);
      $mrna+=(-$degm*$mrna);
      if($mrna<0) { $mrna=0; }
      if($prot<0) { $prot=0; }

      print TM2 "$i\t$mrna\t$prot\n"; 

      $prev=$bas;
    }
  }
  if(exists $ttlist1{$i})
  {
    $p1=$ttlist1{$i};
  }
  if(exists $ttlist2{$i})
  {
    $p2=$ttlist2{$i};
  }
}
close(TM);
close(TM2); 

undef $p1;
undef $p2;
undef $ref1;
undef $ref2;
undef $ref3;
undef $ref4;
undef $no;
undef(@arr1);
undef(@arr2);
undef(@arr3);
undef(@arr4);

$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');


$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Trrate) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CooperativeTF_Timevar_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()+ylim(0,180)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm1";
undef %ttlist1; 
undef %ttlist2; 

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, mRNA) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CooperativeTF_Timevar_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'df <- tab %>% select(Time, prot) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CooperativeTF_Timevar_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'rm(tab)');
system "rm $tm2";

undef $tm1; 
undef $tm2; 

## Competitive TFs 
$R->run(q'set.seed(1234)');

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
  $rn=1+int(rand(2));
  $time+=$arr1[$i];
  $time=sprintf("%0.4f",$time);
  $ttlist{$time}=$rate[$rn];
  $time+=$arr2[$i];
  $time=sprintf("%0.4f",$time);
  $ttlist{$time}=$bas;
  undef $rn;
}

print "   CompetitiveTF\t$no\t$time\n";

$prev=$bas;
$prev=$bas;
$mrna=$inimrna;
$prot=$iniprot;

$tm1="CompetitiveTF_expl_$onrate"."_$offrate.txt";
$tm2="CompetitiveTF_mRNA_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;

print TM "Time\tTrrate\n";
print TM2 "Time\tmRNA\tprot\n";

for($i=0;$i<$tottime;$i+=$uot)
{
  $i=sprintf("%0.4f",$i);
  if(exists $ttlist{$i})
  {
    print TM "$i\t$ttlist{$i}\n";

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

    print TM2 "$i\t$mrna\t$prot\n";

    $prev=$ttlist{$i};
  }
  else
  {
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

    print TM2 "$i\t$mrna\t$prot\n";

    print TM "$i\t$prev\n";
  }
}
close(TM);
close(TM2);

undef $prev;
undef $ref1;
undef $ref2;
undef $no;
undef(@arr1);
undef(@arr2);

$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');


$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Trrate) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CompetitiveTF_Timevar_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()+ylim(0,180)');
$R->run(q'rm(tab)');
$R->run(q'rm(df)');
$R->run(q'dev.off()');

system "rm $tm1";
undef %ttlist; 

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, mRNA) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CompetitiveTF_Timevar_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'df <- tab %>% select(Time, prot) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CompetitiveTF_Timevar_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=20,height=6,useDingbats=F)');
undef $wrf; 

$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) +
          theme_minimal()');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

$R->run(q'rm(tab)');
system "rm $tm2";

undef $tm1; 
undef $tm2; 

## ----------------------------------------------------------------------------------------

## Inter-Individual variations and time-averaged variation in burst size 

## -----------------------------------------------------------------------------------------

## Single TF 

print "Plotting Cell-to-Cell variation ...\n"; 
print "   SingleTF binding\n"; 

for($cn=0;$cn<10;$cn++)
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
  close(TM);

  undef $prev;
  undef $ref1;
  undef $ref2;
  undef $no;
  undef(@arr1);
  undef(@arr2);
  
  undef %ttlist; 
}

$tm1="TMP_expl_$onrate"."_$offrate.txt";
$tm2="TMP_mRNA_expl_$onrate"."_$offrate.txt";
$tm3="TMP_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;
open(TM3,">$tm3") or die;

print TM "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM2 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM3 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
for($i=0;$i<$cl;$i++)
{
  print TM "$timelist[$i]\t";
  print TM2 "$timelist[$i]\t";
  print TM3 "$timelist[$i]\t";
  for($cn=0;$cn<10;$cn++)
  {
    $val=$trlist[$cn][$i]+250*$cn;
    $val2=$mrnalist[$cn][$i]+2500*$cn;
    $val3=$protlist[$cn][$i]+5000000*$cn;
    print TM "$val\t";
    print TM2 "$val2\t";
    print TM3 "$val3\t";
    undef $val; 
    undef $val2;
    undef $val3; 
  }
  print TM "\n";
  print TM2 "\n";
  print TM3 "\n";
}
close(TM);
close(TM2);
close(TM3);

undef $cl;
undef(@trlist);
undef(@timelist); 
undef(@mrnalist); 
undef(@protlist); 


$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');

$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/SingleTF_CelltoCell_var_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,2500)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E")) + theme_minimal()+ylim(0,1500)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm1";

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/SingleTF_CelltoCell_var_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,25000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E"  )) + theme_minimal()+ylim(0,15000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm2";

$R->run(qq'tab=read.table("$tm3",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/SingleTF_CelltoCell_var_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,50000000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" )) + theme_minimal()+ylim(0,30000000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm3";

undef $tm1; 
undef $tm2; 
undef $tm3; 

#Multi TF 
print "   MultiTF binding\n"; 

for($cn=0;$cn<10;$cn++)
{
  $R->run(qq'di=rexp($numevt,rate=$onrate)');
  $ref1=$R->get('di');
  $R->run(q'rm(di)');
  @arr1=@{$ref1};

  $R->run(qq'di=rexp($numevt,rate=$offrate)');
  $ref2=$R->get('di');
  $R->run(q'rm(di)');
  @arr2=@{$ref2};

  $R->run(qq'di=rexp($numevt,rate=$onrate)');
  $ref3=$R->get('di');
  $R->run(q'rm(di)');
  @arr3=@{$ref3};

  $R->run(qq'di=rexp($numevt,rate=$offrate)');
  $ref4=$R->get('di');
  $R->run(q'rm(di)');
  @arr4=@{$ref4};

  $no=@arr1;

  $time1=0;
  $time2=0; 

  for($i=0;$i<$no;$i++)
  {
    $time1+=$arr1[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=$rate[1];
    $time1+=$arr2[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=$bas;

    $time2+=$arr3[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=$rate[2];
    $time2+=$arr4[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=$bas;
  }

  if($time1>$time2) { $time=$time1; } 
  else { $time=$time2; }
  
  $prev=$bas;
  $mrna=$inimrna;
  $prot=$iniprot; 

  $cl=0; 
  for($i=0;$i<$tottime;$i+=$uot)
  {
    $i=sprintf("%0.4f",$i);
    $timelist[$cl]=$i;
    if(exists $ttlist1{$i} && exists $ttlist2{$i})
    {
      $rt=($ttlist1{$i}+$ttlist2{$i})*0.77; 
      $trlist[$cn][$cl]=$rt; 

      if($prev>2*$bas)
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

      $prev=$rt;
    }
    elsif(exists $ttlist1{$i})
    {
      $rt=$ttlist1{$i}*0.77; 
      $trlist[$cn][$cl]=$rt; 

      if($prev>2*$bas)
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

      $prev=$rt;
    }
    elsif(exists $ttlist2{$i})
    {
      $rt=$ttlist2{$i}*0.77; 
      $trlist[$cn][$cl]=$rt; 

      if($prev>2*$bas)
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

      $prev=$rt;
    }
    else
    {
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

      $trlist[$cn][$cl]=$prev; 
      $mrnalist[$cn][$cl]=$mrna;
      $protlist[$cn][$cl]=$prot;
       
    }
    $cl++;
  }
  close(TM);

  undef $prev;
  undef $ref1;
  undef $ref2;
  undef $ref3;
  undef $ref4;
  undef $no;
  undef(@arr1);
  undef(@arr2);
  undef(@arr3);
  undef(@arr4);
  
  undef %ttlist1; 
  undef %ttlist2; 
}

$tm1="TMP_expl_$onrate"."_$offrate.txt";
$tm2="TMP_mRNA_expl_$onrate"."_$offrate.txt";
$tm3="TMP_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;
open(TM3,">$tm3") or die;

print TM "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM2 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM3 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
for($i=0;$i<$cl;$i++)
{
  print TM "$timelist[$i]\t";
  print TM2 "$timelist[$i]\t";
  print TM3 "$timelist[$i]\t";
  for($cn=0;$cn<10;$cn++)
  {
    $val=$trlist[$cn][$i]+250*$cn;
    $val2=$mrnalist[$cn][$i]+2500*$cn;
    $val3=$protlist[$cn][$i]+5000000*$cn;
    print TM "$val\t";
    print TM2 "$val2\t";
    print TM3 "$val3\t";
    undef $val; 
    undef $val2;
    undef $val3; 
  }
  print TM "\n";
  print TM2 "\n";
  print TM3 "\n";
}
close(TM);
close(TM2);
close(TM3);

undef $cl;
undef(@trlist);
undef(@timelist); 
undef(@mrnalist); 
undef(@protlist); 


$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');

$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/MultiTF_CelltoCell_var_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,2500)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E")) + theme_minimal()+ylim(0,1500)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm1";

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/MultiTF_CelltoCell_var_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,25000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E"  )) + theme_minimal()+ylim(0,15000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm2";

$R->run(qq'tab=read.table("$tm3",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/MultiTF_CelltoCell_var_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,50000000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" )) + theme_minimal()+ylim(0,30000000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm3";

undef $tm1; 
undef $tm2; 
undef $tm3; 
##  Cooperative TF 
print "   CooperativeTF binding\n"; 

for($cn=0;$cn<10;$cn++)
{
  $R->run(qq'di=rexp($numevt,rate=$onrate*$coopfac1)');
  $ref1=$R->get('di');
  $R->run(q'rm(di)');
  @arr1=@{$ref1};

  $R->run(qq'di=rexp($numevt,rate=$offrate/$coopfac2)');
  $ref2=$R->get('di');
  $R->run(q'rm(di)');
  @arr2=@{$ref2};

  $R->run(qq'di=rexp($numevt,rate=$onrate*$coopfac1)');
  $ref3=$R->get('di');
  $R->run(q'rm(di)');
  @arr3=@{$ref3};

  $R->run(qq'di=rexp($numevt,rate=$offrate/$coopfac2)');
  $ref4=$R->get('di');
  $R->run(q'rm(di)');
  @arr4=@{$ref4};

  $no=@arr1;

  $time1=0;
  $time2=0;

  for($i=0;$i<$no;$i++)
  {
    $time1+=$arr1[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=1;
    $time1+=$arr2[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=0;

    $time2+=$arr3[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=1;
    $time2+=$arr4[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=0;
  }

  $cl=0; 
  $p1=0; $p2=0; 
  $mrna=$inimrna;
  $prot=$iniprot;
  $prev=$bas;
  for($i=0;$i<$tottime;$i+=$uot)
  {
    $i=sprintf("%0.4f",$i);
    $timelist[$cl]=$i;
    if(exists $ttlist1{$i} && exists $ttlist2{$i})
    {
      if($ttlist1{$i}==1 && $ttlist2{$i}==1)
      {
        $trlist[$cn][$cl]=$rate[0];

        if($prev==$rate[0])
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

        $prev=$rate[0]; 
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
    }
    else
    { 
      if($p1==1 && $p2==1) 
      {
        $trlist[$cn][$cl]=$rate[0];
        if($prev==$rate[0])
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

        $prev=$rate[0];
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
    }
    if(exists $ttlist1{$i})
    {
      $p1=$ttlist1{$i};
    }
    if(exists $ttlist2{$i})
    {
      $p2=$ttlist2{$i};
    }
    $cl++; 
  }
  close(TM);

  undef $p1;
  undef $p2;
  undef $ref1;
  undef $ref2;
  undef $ref3;
  undef $ref4;
  undef $no;
  undef(@arr1);
  undef(@arr2);
  undef(@arr3);
  undef(@arr4);
  undef %ttlist1;
  undef %ttlist2;
}

$tm1="TMP_expl_$onrate"."_$offrate.txt";
$tm2="TMP_mRNA_expl_$onrate"."_$offrate.txt";
$tm3="TMP_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;
open(TM3,">$tm3") or die;

print TM "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM2 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM3 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
for($i=0;$i<$cl;$i++)
{
  print TM "$timelist[$i]\t";
  print TM2 "$timelist[$i]\t";
  print TM3 "$timelist[$i]\t";
  for($cn=0;$cn<10;$cn++)
  {
    $val=$trlist[$cn][$i]+250*$cn;
    $val2=$mrnalist[$cn][$i]+2500*$cn;
    $val3=$protlist[$cn][$i]+5000000*$cn;
    print TM "$val\t";
    print TM2 "$val2\t";
    print TM3 "$val3\t";
    undef $val;
  }
  print TM "\n";
  print TM2 "\n";
  print TM3 "\n";
}
close(TM);
close(TM2);
close(TM3);

undef $cl;
undef(@trlist);
undef(@timelist);
undef(@mrnalist);
undef(@protlist); 


$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');

$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CooperativeTF_CelltoCell_var_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,2500)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" )) + theme_minimal()+ylim(0,1500)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm1";

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CooperativeTF_CelltoCell_var_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,25000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E"  )) + theme_minimal()+ylim(0,15000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm2";

$R->run(qq'tab=read.table("$tm3",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CooperativeTF_CelltoCell_var_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,50000000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" )) + theme_minimal()+ylim(0,30000000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm3";

undef $tm1; 
undef $tm2; 
undef $tm3; 

## Competitive TFs
print "   CompetitiveTF binding\n"; 

for($cn=0;$cn<10;$cn++)
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
    $rn=1+int(rand(2));
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

$tm1="TMP_expl_$onrate"."_$offrate.txt";
$tm2="TMP_mRNA_expl_$onrate"."_$offrate.txt";
$tm3="TMP_prot_expl_$onrate"."_$offrate.txt";

open(TM,">$tm1") or die;
open(TM2,">$tm2") or die;
open(TM3,">$tm3") or die;

print TM "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM2 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
print TM3 "Time\tCell1\tCell2\tCell3\tCell4\tCell5\tCell6\tCell7\tCell8\tCell9\tCell10\n";
for($i=0;$i<$cl;$i++)
{
  print TM "$timelist[$i]\t";
  print TM2 "$timelist[$i]\t";
  print TM3 "$timelist[$i]\t";
  for($cn=0;$cn<10;$cn++)
  {
    $val=$trlist[$cn][$i]+250*$cn;
    $val2=$mrnalist[$cn][$i]+2500*$cn;
    $val3=$protlist[$cn][$i]+5000000*$cn;
    print TM "$val\t";
    print TM2 "$val2\t";
    print TM3 "$val3\t";
    undef $val; 
    undef $val2; 
    undef $val3; 
  }
  print TM "\n";
  print TM2 "\n";
  print TM3 "\n";
}
close(TM);
close(TM2);
close(TM3);

undef $cl;
undef(@trlist);
undef(@timelist); 
undef(@mrnalist); 
undef(@protlist); 


$R->run(q'library(ggplot2)');
$R->run(q'library(tidyr)');
$R->run(q'library(dplyr)');

$R->run(qq'tab=read.table("$tm1",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CompetitiveTF_CelltoCell_var_trarate"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,2500)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E")) + theme_minimal()+ylim(0,1500)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm1";

$R->run(qq'tab=read.table("$tm2",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CompetitiveTF_CelltoCell_var_mRNA"."_$onrate"."_$offrate.pdf";
$R->run(q'pdf(file="CompetitiveTF_CelltoCell_var_mRNA.pdf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,25000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E")) + theme_minimal()+ylim(0,15000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm2";

$R->run(qq'tab=read.table("$tm3",header=T)');

$R->run(q'df <- tab %>% select(Time, Cell1, Cell2, Cell3, Cell4, Cell5, Cell6) %>% gather(key = "variable", value = "value", -Time)');

$wrf="Plots_sim_onrate_offrate_expl/CompetitiveTF_CelltoCell_var_prot"."_$onrate"."_$offrate.pdf";
$R->run(qq'pdf(file="$wrf",width=12,height=6,useDingbats=F)');
undef $wrf; 

##$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E" , "#808000" , "#FA8072" , "#0000FF" , "#DC143C" )) + theme_minimal()+ylim(0,50000000)');
$R->run(q'ggplot(df, aes(x = Time, y = value)) + geom_line(aes(color = variable), size = 1) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FF4500" , "#228B22" , "#C71585" , "#D2691E")) + theme_minimal()+ylim(0,30000000)');

$R->run(q'rm(tab)');
$R->run(q'rm(df)');

$R->run(q'dev.off()');

system "rm $tm3";

undef $tm1; 
undef $tm2; 
undef $tm3; 

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

$wrf=">Plots_sim_onrate_offrate_expl/NoiseData_SingleTFbinding"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die; 
undef $wrf; 
print WR "Timepoint\tAvg_TR\tSd_TR\tCV_TR\tAvg_mRNA\tSd_mRNA\tCV_mRNA\tAvg_Prot\tSd_Prot\tCV_Prot\n"; 


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
  $sdm=sprintf("%0.2f",sqrt($sdm));
  $avgm=sprintf("%0.2f",$avgm/$numcell);
  if($avgm!=0) 
  { 
     $cvm=sprintf("%0.2f",$sdm/$avgm);
  }
  else 
  { 
     $cvm="NA"; 
  } 

  $sdp=(($sdp/$numcell)-(($avgp/$numcell)**2));
  if($sdp<0) { $sdp=0; }
  $sdp=sprintf("%0.2f",sqrt($sdp));
  $avgp=sprintf("%0.2f",$avgp/$numcell);
  if($avgp!=0) 
  { 
     $cvp=sprintf("%0.2f",$sdp/$avgp);
  }
  else
  {
     $cvp="NA"; 
  }

  $sdtr=(($sdtr/$numcell)-(($avgtr/$numcell)**2));
  if($sdtr<0) { $sdtr=0; }
  $sdtr=sprintf("%0.2f",sqrt($sdtr));
  $avgtr=sprintf("%0.2f",$avgtr/$numcell);
  if($avgtr!=0) 
  { 
     $cvtr=sprintf("%0.2f",$sdtr/$avgtr);
  }
  else
  {
     $cvtr="NA"; 
  }

  if($i>150) { print WR "$i\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; }

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
close(WR); 

$wrf=">Plots_sim_onrate_offrate_expl/SingleTF_Burst_SizeFreq"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die; 
undef $wrf; 

print WR "CellID\tBurstSize\tBurstFreq\n"; 
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
    $bsize=sprintf("%0.2f",$bsize/$bfreq); 
  }
  else
  {
    $bsize="NA"; 
  }
  $bfreq=sprintf("%0.2f",$bfreq/$cl); 
  print WR "$cn\t$bsize\t$bfreq\n"; 
  undef $bsize; 
  undef $bfreq; 
}
close(WR); 

undef $cl;
undef(@trlist);
undef(@timelist); 
undef(@mrnalist); 
undef(@protlist); 

## Multi TF
print "   MultiTF binding\n"; 

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

  $R->run(qq'di=rexp($numevt,rate=$onrate)');
  $ref3=$R->get('di');
  $R->run(q'rm(di)');
  @arr3=@{$ref3};

  $R->run(qq'di=rexp($numevt,rate=$offrate)');
  $ref4=$R->get('di');
  $R->run(q'rm(di)');
  @arr4=@{$ref4};

  $no=@arr1;

  $time1=0;
  $time2=0; 

  for($i=0;$i<$no;$i++)
  {
    $time1+=$arr1[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=$rate[1];
    $time1+=$arr2[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=$bas;

    $time2+=$arr3[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=$rate[2];
    $time2+=$arr4[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=$bas;
  }

  if($time1>$time2) { $time=$time1; } 
  else { $time=$time2; }

  print "   MultiTF\t$no\t$time\n";
  
  $prev=$bas;
  $mrna=$inimrna;
  $prot=$iniprot; 

  $cl=0; 
  for($i=0;$i<$tottime;$i+=$uot)
  {
    $i=sprintf("%0.4f",$i);
    $timelist[$cl]=$i;
    if(exists $ttlist1{$i} && exists $ttlist2{$i})
    {
      $rt=($ttlist1{$i}+$ttlist2{$i})*0.77; 
      $trlist[$cn][$cl]=$rt; 

      if($prev>2*$bas)
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

      $prev=$rt;
    }
    elsif(exists $ttlist1{$i})
    {
      $rt=$ttlist1{$i}*0.77; 
      $trlist[$cn][$cl]=$rt; 

      if($prev>2*$bas)
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

      $prev=$rt;
    }
    elsif(exists $ttlist2{$i})
    {
      $rt=$ttlist2{$i}*0.77; 
      $trlist[$cn][$cl]=$rt; 

      if($prev>2*$bas)
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

      $prev=$rt;
    }
    else
    {
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

      $trlist[$cn][$cl]=$prev; 
      $mrnalist[$cn][$cl]=$mrna;
      $protlist[$cn][$cl]=$prot;
       
    }
    $cl++;
  }
  undef $prev;
  undef $ref1;
  undef $ref2;
  undef $ref3;
  undef $ref4;
  undef $no;
  undef(@arr1);
  undef(@arr2);
  undef(@arr3);
  undef(@arr4);
  
  undef %ttlist1; 
  undef %ttlist2; 
}

$wrf=">Plots_sim_onrate_offrate_expl/NoiseData_MultiTFbinding"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die; 
undef $wrf; 
print WR "Timepoint\tAvg_TR\tSd_TR\tCV_TR\tAvg_mRNA\tSd_mRNA\tCV_mRNA\tAvg_Prot\tSd_Prot\tCV_Prot\n"; 


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
  $sdm=sprintf("%0.2f",sqrt($sdm));
  $avgm=sprintf("%0.2f",$avgm/$numcell);
  if($avgm!=0) 
  { 
     $cvm=sprintf("%0.2f",$sdm/$avgm);
  }
  else 
  { 
     $cvm="NA"; 
  } 

  $sdp=(($sdp/$numcell)-(($avgp/$numcell)**2));
  if($sdp<0) { $sdp=0; }
  $sdp=sprintf("%0.2f",sqrt($sdp));
  $avgp=sprintf("%0.2f",$avgp/$numcell);
  if($avgp!=0) 
  { 
     $cvp=sprintf("%0.2f",$sdp/$avgp);
  }
  else
  {
     $cvp="NA"; 
  }

  $sdtr=(($sdtr/$numcell)-(($avgtr/$numcell)**2));
  if($sdtr<0) { $sdtr=0; }
  $sdtr=sprintf("%0.2f",sqrt($sdtr));
  $avgtr=sprintf("%0.2f",$avgtr/$numcell);
  if($avgtr!=0) 
  { 
     $cvtr=sprintf("%0.2f",$sdtr/$avgtr);
  }
  else
  {
     $cvtr="NA"; 
  }

  if($i>150) { print WR "$i\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; }

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
close(WR); 

$wrf=">Plots_sim_onrate_offrate_expl/MultiTF_Burst_SizeFreq"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die; 
undef $wrf; 

print WR "CellID\tBurstSize\tBurstFreq\n"; 
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
    $bsize=sprintf("%0.2f",$bsize/$bfreq); 
  }
  else
  {
    $bsize="NA"; 
  }
  $bfreq=sprintf("%0.2f",$bfreq/$cl); 
  print WR "$cn\t$bsize\t$bfreq\n"; 
  undef $bsize; 
  undef $bfreq; 
}
close(WR); 

undef $cl;
undef(@trlist);
undef(@timelist); 
undef(@mrnalist); 
undef(@protlist); 

## Cooperative TF 

print "   CooperativeTF binding\n"; 

for($cn=0;$cn<$numcell;$cn++)
{
  $R->run(qq'di=rexp($numevt,rate=$onrate*$coopfac1)');
  $ref1=$R->get('di');
  $R->run(q'rm(di)');
  @arr1=@{$ref1};

  $R->run(qq'di=rexp($numevt,rate=$offrate/$coopfac2)');
  $ref2=$R->get('di');
  $R->run(q'rm(di)');
  @arr2=@{$ref2};

  $R->run(qq'di=rexp($numevt,rate=$onrate*$coopfac1)');
  $ref3=$R->get('di');
  $R->run(q'rm(di)');
  @arr3=@{$ref3};

  $R->run(qq'di=rexp($numevt,rate=$offrate/$coopfac2)');
  $ref4=$R->get('di');
  $R->run(q'rm(di)');
  @arr4=@{$ref4};

  $no=@arr1;

  $time1=0;
  $time2=0;

  for($i=0;$i<$no;$i++)
  {
    $time1+=$arr1[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=1;
    $time1+=$arr2[$i];
    $time1=sprintf("%0.4f",$time1);
    $ttlist1{$time1}=0;

    $time2+=$arr3[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=1;
    $time2+=$arr4[$i];
    $time2=sprintf("%0.4f",$time2);
    $ttlist2{$time2}=0;
  }

  $cl=0; 
  $p1=0; $p2=0; 
  $mrna=$inimrna;
  $prot=$iniprot;
  $prev=$bas;
  for($i=0;$i<$tottime;$i+=$uot)
  {
    $i=sprintf("%0.4f",$i);
    $timelist[$cl]=$i;
    if(exists $ttlist1{$i} && exists $ttlist2{$i})
    {
      if($ttlist1{$i}==1 && $ttlist2{$i}==1)
      {
        $trlist[$cn][$cl]=$rate[0];

        if($prev==$rate[0])
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

        $prev=$rate[0]; 
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
    }
    else
    { 
      if($p1==1 && $p2==1) 
      {
        $trlist[$cn][$cl]=$rate[0];
        if($prev==$rate[0])
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

        $prev=$rate[0];
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
    }
    if(exists $ttlist1{$i})
    {
      $p1=$ttlist1{$i};
    }
    if(exists $ttlist2{$i})
    {
      $p2=$ttlist2{$i};
    }
    $cl++; 
  }
  undef $p1;
  undef $p2;
  undef $ref1;
  undef $ref2;
  undef $ref3;
  undef $ref4;
  undef $no;
  undef(@arr1);
  undef(@arr2);
  undef(@arr3);
  undef(@arr4);
  undef %ttlist1;
  undef %ttlist2;
}

$wrf=">Plots_sim_onrate_offrate_expl/NoiseData_CooperativeTFbinding"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die; 
undef $wrf; 

print WR "Timepoint\tAvg_TR\tSd_TR\tCV_TR\tAvg_mRNA\tSd_mRNA\tCV_mRNA\tAvg_Prot\tSd_Prot\tCV_Prot\n"; 
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
  $sdm=sprintf("%0.2f",sqrt($sdm));
  $avgm=sprintf("%0.2f",$avgm/$numcell);
  if($avgm!=0)
  {
     $cvm=sprintf("%0.2f",$sdm/$avgm);
  }
  else
  {
     $cvm="NA"; 
  }

  $sdp=(($sdp/$numcell)-(($avgp/$numcell)**2));
  if($sdp<0) { $sdp=0; }
  $sdp=sprintf("%0.2f",sqrt($sdp));
  $avgp=sprintf("%0.2f",$avgp/$numcell);
  if($avgp!=0)
  {
     $cvp=sprintf("%0.2f",$sdp/$avgp);
  }
  else
  {
     $cvp="NA"; 
  }

  $sdtr=(($sdtr/$numcell)-(($avgtr/$numcell)**2));
  if($sdtr<0) { $sdtr=0; }
  $sdtr=sprintf("%0.2f",sqrt($sdtr));
  $avgtr=sprintf("%0.2f",$avgtr/$numcell);
  if($avgtr!=0)
  {
     $cvtr=sprintf("%0.2f",$sdtr/$avgtr);
  }
  else
  {
     $cvtr="NA"; 
  }

  if($i>150) { print WR "$i\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; }

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
close(WR); 

$wrf=">Plots_sim_onrate_offrate_expl/CooperativeTF_Burst_SizeFreq"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die;
undef $wrf; 

print WR "CellID\tBurstSize\tBurstFreq\n";
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
     $bsize=sprintf("%0.2f",$bsize/$bfreq);
  }
  else
  {
     $bsize="NA"; 
  }
  $bfreq=sprintf("%0.2f",$bfreq/$cl);
  print WR "$cn\t$bsize\t$bfreq\n";
  undef $bsize;
  undef $bfreq;
}
close(WR);

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
    $rn=1+int(rand(2));
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

$wrf=">Plots_sim_onrate_offrate_expl/NoiseData_CompetitiveTFbinding"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die; 
undef $wrf; 

print WR "Timepoint\tAvg_TR\tSd_TR\tCV_TR\tAvg_mRNA\tSd_mRNA\tCV_mRNA\tAvg_Prot\tSd_Prot\tCV_Prot\n"; 
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
  $sdm=sprintf("%0.2f",sqrt($sdm));
  $avgm=sprintf("%0.2f",$avgm/$numcell);
  if($avgm!=0)
  {
     $cvm=sprintf("%0.2f",$sdm/$avgm);
  }
  else
  {
     $cvm="NA"; 
  }

  $sdp=(($sdp/$numcell)-(($avgp/$numcell)**2));
  if($sdp<0) { $sdp=0; }
  $sdp=sprintf("%0.2f",sqrt($sdp));
  $avgp=sprintf("%0.2f",$avgp/$numcell);
  if($avgp!=0)
  {
     $cvp=sprintf("%0.2f",$sdp/$avgp);
  }
  else
  {
     $cvp="NA"; 
  }

  $sdtr=(($sdtr/$numcell)-(($avgtr/$numcell)**2));
  if($sdtr<0) { $sdtr=0; }
  $sdtr=sprintf("%0.2f",sqrt($sdtr));
  $avgtr=sprintf("%0.2f",$avgtr/$numcell);
  if($avgtr!=0)
  {
     $cvtr=sprintf("%0.2f",$sdtr/$avgtr);
  }
  else
  {
     $cvtr="NA"; 
  }

  if($i>150) { print WR "$i\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; }

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
close(WR);

$wrf=">Plots_sim_onrate_offrate_expl/CompetitiveTF_Burst_SizeFreq"."_$onrate"."_$offrate.txt";
open(WR,">$wrf") or die;
undef $wrf; 

print WR "CellID\tBurstSize\tBurstFreq\n";
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
     $bsize=sprintf("%0.2f",$bsize/$bfreq);
  }
  else
  {
     $bsize="NA"; 
  }
  $bfreq=sprintf("%0.2f",$bfreq/$cl);
  print WR "$cn\t$bsize\t$bfreq\n";
  undef $bsize;
  undef $bfreq;
}
close(WR);

undef $cl;
undef(@trlist);
undef(@timelist);
undef(@mrnalist);
undef(@protlist);

$R->stopR;
$pm->finish; # do the exit in the child process

}

$pm->wait_all_children;

undef(@val); 

undef(@ont); 
undef(@oft); 
