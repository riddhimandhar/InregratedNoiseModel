#!/usr/bin/perl -w 

#$unfilttab="Filtered_MERGED_DATA_Yeast.txt"; 
$unfilttab="Filtered_nondup_MERGED_DATA_Yeast.txt"; 
$iter=1000; 

use Statistics::R;
$R=Statistics::R->new();
$R->startR;

$R->run(qq'tab=read.table("$unfilttab",header=T)'); 
##$R->run(q'inputData=tab[,c(4:229)]'); ## Filtered MERGED data OLD
##$R->run(q'inputData=tab[,c(6:241)]'); ## Filtered MERGED data 
$R->run(q'inputData=tab[,c(4,7:288)]');
#$R->run(q'inputData=inputData[,which(colSums(abs(inputData))!=0)]');
$R->run(q'inputData=inputData[,colSums(is.na(inputData))<0.8*nrow(inputData)]');
$R->run(q'inputData=inputData[,which(apply(inputData, 2, var, na.rm=TRUE) != 0)]');
$R->run(q'inputData=as.data.frame(scale(inputData))');
$R->run(q'write.table(inputData,"UNFILTGENE_nondup")');

$cn=0; 
open(FP,"UNFILTGENE_nondup") or die;
while($fp=<FP>)
{
  chomp($fp); 
  if($cn==0)
  {
    @arr=split(/\s+/,$fp); 
    $no=@arr;
    for($i=1;$i<$no;$i++)
    {
      while($arr[$i]=~/\"/) { $arr[$i]=~s/\"//; }
      $var[$i-1]=$arr[$i];
    }
    $varc=$no-1; 
    undef $no; 
    undef(@arr);
    $cn++; 
    last; 
  }
}
close(FP);

for($k=0;$k<$varc;$k++)
{
print "$var[$k]\n"; 
$rdpval[$k]=0;  

for($i=0;$i<$iter;$i++)
{
  $R->run(q'trainingIndex <- sample(1:nrow(inputData), 0.80*nrow(inputData))'); # indices for 80% training data
  $R->run(q'trainingData <- inputData[trainingIndex,]'); # training data
  $R->run(q'testData <- inputData[-trainingIndex,]'); # test data

  $R->run(qq'lmmod <- lm(DM_YPD ~ $var[$k],data=trainingData)');  
  $R->run(q'expv=summary(lmmod)[9][[1]]'); 
  $expval=$R->get('expv');
  $R->run(q'preds <- predict(lmmod, testData, na.action=na.pass)');  # predict on test data
  $R->run(q'actual <- testData$DM_YPD');

  $R->run(q'full <- as.data.frame(cbind(actual,preds))');
  $R->run(q'actual=full$actual[which(full$preds!="NA")]');
  $R->run(q'preds=full$preds[which(full$preds!="NA")]');

  $R->run(q'rss <- sum(na.omit(preds - actual) ^ 2)/length(na.omit(preds-actual))');
  $R->run(q'tss <- sum(na.omit(actual - mean(actual)) ^ 2/length(na.omit(actual-mean(actual))))');
  $R->run(q'rsq <- 1 - rss/tss');
  $rsq=$R->get('rsq');
  $rmse=$R->get('rss');
  $rmse=sqrt($rmse);
  
  $R->run(q'cmv=as.data.frame(cbind(actual,preds))');
  $R->run(q'x=cmv$actual[which(is.na(cmv$actual)=="FALSE" & is.na(cmv$preds)=="FALSE")]');
  $R->run(q'y=cmv$preds[which(is.na(cmv$actual)=="FALSE" & is.na(cmv$preds)=="FALSE")]');

  $R->run(q'ln1=length(x)');
  $R->run(q'ln2=length(y)');
  $ln1=$R->get(q'ln1');
  $ln2=$R->get(q'ln2');
  if($ln1>=3 && $ln2>=3)
  {
    $R->run(q'rv=cor.test(actual,preds)');
    $R->run(q'rves=rv$estimate');
    $R->run(q'rvp=rv$p.value');
    $rves=$R->get('rves');
    $rvp=$R->get('rvp');
    $R->run(q'rm(rv)');
    $R->run(q'rm(rves)');
    $R->run(q'rm(rvp)');
  }
  else
  {
    $rves="NA";
    $rvp="NA";
  }

  $R->run(q'rm(full)');
  $R->run(q'rm(cmv)');
  $R->run(q'rm(x)');
  $R->run(q'rm(y)');
  $R->run(q'rm(ln1)');
  $R->run(q'rm(ln2)');
  $R->run(q'rm(lmmod)');
  $R->run(q'rm(expv)');
  $R->run(q'rm(actual)');
  $R->run(q'rm(preds)');
  $R->run(q'rm(rss)');
  $R->run(q'rm(tss)');
  $R->run(q'rm(rsq)');
  
  $expval=sprintf("%0.3f",$expval); 

  $rdexpli[$k][$i]=$expval;
  $rdrsqli[$k][$i]=$rsq;
  $rdcorli[$k][$i]=$rves;
  $rdrmse[$k][$i]=$rmse;
  if($rvp!~/NA/ && $rvp<0.05) { $rdpval[$k]++; }
  undef $ln1;
  undef $ln2;
  undef $expval;
  undef $rsq;
  undef $rves;
  undef $rvp;
  undef $rmse;

  $R->run(q'rm(testData)');
  $R->run(q'rm(trainingData)');
  $R->run(q'rm(trainingIndex)');
  $R->run(q'rm(tab)');
}
}

$R->run(q'rm(inputData)');
system "rm UNFILTGENE_nondup"; 

$wrfile="IndParam_impact_".$iter."_Iterations_nondup.txt"; 
open(WR,">Results/$wrfile") or die; 
print WR "Unfiltered data $iter iterations\nVariable\tAvg_frVarExpl\tSd_frVarExpl\tAvg_predRsq\tSd_predRsq\tPred_RMSE\tPred_Cor\tPred_Sig_Pval\n------------------------------------------------------------------------------------------------------\n";
open(GR,">Results/Impactful_parameter_list_DM_nondup.txt") or die; 

$imp=0; 
for($k=0;$k<$varc;$k++)
{
  $s1=0; $s2=0; $s3=0; $s4=0;  
  $sd1=0; $sd2=0; $sd3=0; $sd4=0; 

  $c1=0; $c2=0; $c3=0; $c4=0; 

  for($i=0;$i<$iter;$i++)
  {
    if($rdrsqli[$k][$i]!~/NA/ && $rdrsqli[$k][$i]!~/-Inf/)
    {
      $s1+=($rdrsqli[$k][$i]);  
      $sd1+=($rdrsqli[$k][$i])**2;  
      $c1++;
    }
    if($rdcorli[$k][$i]!~/NA/)
    {
      $s2+=($rdcorli[$k][$i]);  
      $sd2+=($rdcorli[$k][$i])**2;  
      $c2++;
    }
    if($rdexpli[$k][$i]!~/NA/ && $rdexpli[$k][$i]!~/-Inf/)
    {
      $s3+=($rdexpli[$k][$i]);  
      $sd3+=($rdexpli[$k][$i])**2;  
      $c3++;
    }
    if($rdrmse[$k][$i]!~/NA/ && $rdrmse[$k][$i]!~/-Inf/)
    {
      $s4+=($rdrmse[$k][$i]);  
      $sd4+=($rdrmse[$k][$i])**2;  
      $c4++;
    }
  }

  if($c1!=0)
  {
    $sd1=sprintf("%0.3f",sqrt(($sd1/$c1)-(($s1/$c1)**2)));
    $s1=sprintf("%0.3f",$s1/$c1); 
  }
  else { $sd1="NA"; $s1="NA"; }
  if($c2!=0)
  {
    $sd2=sprintf("%0.3f",sqrt(($sd2/$c2)-(($s2/$c2)**2)));
    $s2=sprintf("%0.3f",$s2/$c2); 
  }
  else { $sd2="NA"; $s2="NA"; }
  if($c3!=0)
  {
    $sd3=sprintf("%0.3f",sqrt(($sd3/$c3)-(($s3/$c3)**2)));
    $s3=sprintf("%0.3f",$s3/$c3); 
  }
  else { $sd3="NA"; $s3="NA"; } 
  if($c4!=0)
  {
    $sd4=sprintf("%0.3f",sqrt(($sd4/$c4)-(($s4/$c4)**2)));
    $s4=sprintf("%0.3f",$s4/$c4); 
  }
  else { $sd4="NA"; $s4="NA"; } 

  $rdpval[$k]=sprintf("%0.3f",$rdpval[$k]/$iter*100);

  print WR "$var[$k]\t$s3\t$sd3\t$s1\t$sd1\t$s4±$sd4\t$s2±$sd2\t$rdpval[$k]%\n"; 
  
  if($s3>=0.05 || $s1>=0.05 || ($s3>=0.05 && $s1>=0.05)) 
  { 
    print GR "$var[$k]\n"; 
    $imp++;
  }
}
print WR "------------------------------------------------------------------------------------------------------\n";
close(GR); 
close(WR);

undef(@rdexpli);
undef(@rdrsqli);
undef(@rdcorli);
undef(@rdrmse);
undef(@rdpval);

print "Impactful parameters no. : $imp\n"; 

$R->stopR;
