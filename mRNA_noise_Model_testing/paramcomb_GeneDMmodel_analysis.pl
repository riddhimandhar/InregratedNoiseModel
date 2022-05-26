#!/usr/bin/perl -w 

$iter=1000; 

$datatype[0]="Filtered_MERGED_DATA_Yeast.txt"; 
$datatype[1]="CorrFiltered_MERGED_DATA_Yeast.txt"; 
$datatype[2]="ImpFiltered_MERGED_DATA_Yeast.txt"; 

use Statistics::R;
$R=Statistics::R->new();
$R->startR;

$R->run(qq'whtab=read.table("$datatype[0]",header=T)'); 

$R->run(q'library(pls)');
$R->run(q'library(olsrr)');
$R->run(q'library(ridge)');
$R->run(q'library(ISLR)');
$R->run(q'library(glmnet)');
$R->run(q'library(dplyr)');
$R->run(q'library(tidyr)');
$R->run(q'library(randomForestSRC)');

$wrfile="SubModel_ResultsGeneDM_".$iter."_Iterations.txt";
$wrfile2="Formatted_frvarExpl_SubModelGeneDM_".$iter."_Iterations.txt";
$wrfile3="Formatted_predRsq_SubModelGeneDM_".$iter."_Iterations.txt";
open(WR,">Results/$wrfile") or die;
open(WR2,">Results/$wrfile2") or die;
open(WR3,">Results/$wrfile3") or die;

print WR2 "Model_no\tModel\tMean_frv_U_Ridge\tMean_frv_U_RidgeAIC\tMean_frv_U_Lasso\tMean_frv_U_LassoAIC\tMean_frv_U_RF_wi\tMean_frv_U_RF_wo\t";
print WR2 "Mean_frv_F_Ridge\tMean_frv_F_RidgeAIC\tMean_frv_F_Lasso\tMean_frv_F_LassoAIC\tMean_frv_F_RF_wi\tMean_frv_F_RF_wo\t";
print WR2 "Mean_frv_F2_Ridge\tMean_frv_F2_RidgeAIC\tMean_frv_F2_Lasso\tMean_frv_F2_LassoAIC\tMean_frv_F2_RF_wi\tMean_frv_F2_RF_wo\t";
print WR2 "Sd_frv_U_Ridge\tSd_frv_U_RidgeAIC\tSd_frv_U_Lasso\tSd_frv_U_LassoAIC\tSd_frv_U_RF_wi\tSd_frv_U_RF_wo\t";
print WR2 "Sd_frv_F_Ridge\tSd_frv_F_RidgeAIC\tSd_frv_F_Lasso\tSd_frv_F_LassoAIC\tSd_frv_F_RF_wi\tSd_frv_F_RF_wo\t";
print WR2 "Sd_frv_F2_Ridge\tSd_frv_F2_RidgeAIC\tSd_frv_F2_Lasso\tSd_frv_F2_LassoAIC\tSd_frv_F2_RF_wi\tSd_frv_F2_RF_wo\n";

print WR3 "Model_no\tModel\tMean_prRsq_U_Ridge\tMean_prRsq_U_RidgeAIC\tMean_prRsq_U_Lasso\tMean_prRsq_U_LassoAIC\tMean_prRsq_U_RF_wi\tMean_prRsq_U_RF_wo\t";print WR3 "Mean_prRsq_F_Ridge\tMean_prRsq_F_RidgeAIC\tMean_prRsq_F_Lasso\tMean_prRsq_F_LassoAIC\tMean_prRsq_F_RF_wi\tMean_prRsq_F_RF_wo\t";
print WR3 "Mean_prRsq_F2_Ridge\tMean_prRsq_F2_RidgeAIC\tMean_prRsq_F2_Lasso\tMean_prRsq_F2_LassoAIC\tMean_prRsq_F2_RF_wi\tMean_prRsq_F2_RF_wo\t";
print WR3 "Sd_prRsq_U_Ridge\tSd_prRsq_U_RidgeAIC\tSd_prRsq_U_Lasso\tSd_prRsq_U_LassoAIC\tSd_prRsq_U_RF_wi\tSd_prRsq_U_RF_wo\t";
print WR3 "Sd_prRsq_F_Ridge\tSd_prRsq_F_RidgeAIC\tSd_prRsq_F_Lasso\tSd_prRsq_F_LassoAIC\tSd_prRsq_F_RF_wi\tSd_prRsq_F_RF_wo\t";
print WR3 "Sd_prRsq_F2_Ridge\tSd_prRsq_F2_RidgeAIC\tSd_prRsq_F2_Lasso\tSd_prRsq_F2_LassoAIC\tSd_prRsq_F2_RF_wi\tSd_prRsq_F2_RF_wo\n";

$model[0]="Model2 DM_YPD~PromNucOcc";
$model[1]="Model3 DM_YPD~TATAbox+PromNucOcc";  
$model[2]="Model4 DM_YPD~Histone_Mods";  
$model[3]="Model5 DM_YPD~TATAbox+PromNucOcc+Histone_Mods";  
$model[4]="Model6 DM_YPD~3D_interactions";  
$model[5]="Model7 DM_YPD~TF_binding";  
$model[6]="Model8 DM_YPD~TATAbox+PromNucOcc+Histone_Mods+TF_binding";  
$model[7]="Model9 DM_YPD~TATAbox+PromNucOcc+Histone_Mods+TF_binding+mRNAprop";  
$model[8]="Model10 DM_YPD~SAGA_TFIIDprops";
$model[9]="Model11 DM_YPD~TBPprop_NSMB2009";
$model[10]="Model12 DM_YPD~Nuc_Assoc_prot";
$model[11]="Model13 DM_YPD~Hist_dynamics";
$model[12]="Model14 DM_YPD~Model11-13";
$model[13]="Model15 DM_YPD~Model9+Model11-13";
$model[14]="Model16 DM_YPD~TATAbox+PromNucOcc+Histone_Mods+TF_binding+mRNAprop+PTM";  
$model[15]="Model17 DM_YPD~Model9+Model11-13+PTM";

for($modt=0;$modt<16;$modt++)
{
print WR2 "$model[$modt]\t";
print WR3 "$model[$modt]\t";

print WR "\n$model[$modt]\n"; 

for($dtc=0;$dtc<3;$dtc++)
{
## DATA TYPE

$R->run(qq'tab=read.table("$datatype[$dtc]",header=T)'); 

if($dtc==0) 
{ 
   print "Unfiltered Data ...$model[$modt]\n";  
   $R->run(q'tab=tab[,which(colSums(is.na(tab))<2000)]');
   $R->run(q'tab_filt1 <- tab[,c(4,7:242)]');
   $R->run(q'tab_filt2 <- tab_filt1[,which(apply(tab_filt1, 2, var, na.rm=TRUE) != 0)]'); 
   if($modt==0) { $R->run(q'inputData=tab_filt2[,c(1,4:32)]'); } ## PromNucOcc 
   if($modt==1) { $R->run(q'inputData=tab_filt2[,c(1,3:32)]'); } ## TATAbox + PromNucOcc
   if($modt==2) { $R->run(q'inputData=tab_filt2[,c(1,167:218)]'); } ## Histone_mods
   if($modt==3) { $R->run(q'inputData=tab_filt2[,c(1,3:32,167:218)]'); }  ##TATAbox + PromNucOcc + Histone_mods
   if($modt==4) 
   { 
     print WR "---------------------------------------------------------\n";

     print WR2 "---\t---\t---\t---\t---\t---\t";
     print WR3 "---\t---\t---\t---\t---\t---\t";
     for($mu=0;$mu<12;$mu++)
     {
       $res[$dtc][$mu]="---";
     }
     next; 
   }
   if($modt==5) { $R->run(q'inputData=tab_filt2[,c(1,33:157)]'); }   ## TF_binding
   if($modt==6) { $R->run(q'inputData=tab_filt2[,c(1,3:157,167:218)]'); }  ## TATAbox+PromNucOcc+Histone_mods+TF_binding 
   if($modt==7) { $R->run(q'inputData=tab_filt2[,c(1,3:157,167:218,223)]'); }  ## TATAbox+PromNucOcc+Histone_mods+TF_binding+mRNAprop
   if($modt==8 || $modt==9 || $modt==10 || $modt==11 || $modt==12 || $modt==13 || $modt==15) ## 
   {
     print WR "---------------------------------------------------------\n";
     print WR2 "---\t---\t---\t---\t---\t---\t";
     print WR3 "---\t---\t---\t---\t---\t---\t";
     for($mu=0;$mu<12;$mu++)
     {
       $res[$dtc][$mu]="---";
     }
     next;
   } 
   if($modt==14) { $R->run(q'inputData=tab_filt2[,c(1,3:157,167:218,223:237)]'); }  ## TATAbox+PromNucOcc+Histone_mods+TF_binding+mRNAprop+PTM
    ##print "$name $pval\n"; 
}
if($dtc==1)
{
   print "CorrFiltered Data ...$model[$modt]\n";
   $R->run(q'tab_filt1 <- tab[,c(4:168)]');
   $R->run(q'tab_filt2 <- tab_filt1[,which(apply(tab_filt1, 2, var, na.rm=TRUE) != 0)]');
   if($modt==0) { $R->run(q'inputData=tab_filt2[,c(1,4:27)]'); } ## PromNucOcc
   if($modt==1) { $R->run(q'inputData=tab_filt2[,c(1,3:27)]'); } ## TATAbox + PromNucOcc
   if($modt==2) { $R->run(q'inputData=tab_filt2[,c(1,111:141)]'); } ## Histone_mods
   if($modt==3) { $R->run(q'inputData=tab_filt2[,c(1,3:27,111:141)]'); }  ##TATAbox + PromNucOcc + Histone_mods
   if($modt==4 || $modt==8 || $modt==9 || $modt==10 || $modt==11 || $modt==12)
   {
     print WR "---------------------------------------------------------\n";
     print WR2 "---\t---\t---\t---\t---\t---\t";
     print WR3 "---\t---\t---\t---\t---\t---\t";
     for($mu=0;$mu<12;$mu++)
     {
       $res[$dtc][$mu]="---";
     }
     next;
   }
   if($modt==5) { $R->run(q'inputData=tab_filt2[,c(1,28:101)]'); }   ## TF_binding
   if($modt==6) { $R->run(q'inputData=tab_filt2[,c(1,3:101,111:141)]'); }  ## TATAbox+PromNucOcc+Histone_mods+TF_binding
   if($modt==7) { $R->run(q'inputData=tab_filt2[,c(1,3:101,111:141,144:150)]'); }  ## TATAbox+PromNucOcc+Histone_mods+TF_binding+mRNAprop
   if($modt==13) { $R->run(q'inputData=tab_filt2[,c(1,3:101,111:141,144:160)]'); }  ## Combinations 
   if($modt==14) { $R->run(q'inputData=tab_filt2[,c(1,3:101,111:141,144:150,161:165)]'); }  ## TATAbox+PromNucOcc+Histone_mods+TF_binding+mRNAprop+PTM
   if($modt==15) { $R->run(q'inputData=tab_filt2[,c(1,3:101,111:141,144:165)]'); }  ## Combinations+PTM 
}
if($dtc==2) ###
{
   print "ImpFiltered Data ...$model[$modt]\n";
   $R->run(q'tab_filt1 <- tab[,c(4,7:33)]');
   $R->run(q'tab_filt2 <- tab_filt1[,which(apply(tab_filt1, 2, var, na.rm=TRUE) != 0)]');
   if($modt>=0 && $modt<=4)
   {
     print WR "---------------------------------------------------------\n";
     print WR2 "---\t---\t---\t---\t---\t---\t";
     print WR3 "---\t---\t---\t---\t---\t---\t";
     for($mu=0;$mu<12;$mu++)
     {
       $res[$dtc][$mu]="---";
     }
     next;
   }
   if($modt==5) { $R->run(q'inputData=tab_filt2[,c(1,3:25)]'); }   ## TF_binding
   if($modt==6) { $R->run(q'inputData=tab_filt2[,c(1:26)]'); }  ## TATAbox+PromNucOcc+Histone_mods+TF_binding
   if(($modt>=7 && $modt<=12) || ($modt==14) || ($modt==15))
   {
     print WR "---------------------------------------------------------\n";
     print WR2 "---\t---\t---\t---\t---\t---\t";
     print WR3 "---\t---\t---\t---\t---\t---\t";
     for($mu=0;$mu<12;$mu++)
     {
       $res[$dtc][$mu]="---";
     }
     next;
   }
   if($modt==13) { $R->run(q'inputData=tab_filt2[,c(1:28)]'); }  ## Combinations 
}

$R->run(q'inputData=as.data.frame(scale(inputData))');

$R->run(q'write.table(inputData,"TMP.txt")');

open(TM,"TMP.txt") or die; 
open(TR,">inputData.txt") or die; 
$cc=0; 

while($tm=<TM>)
{
  chomp($tm); 
  @arr=split(/\s+/,$tm);
  $no=@arr;
  if($tm=~/^\"DM/)
  {
    for($i=0;$i<$no;$i++)
    {
      @qw=split(/\"/,$arr[$i]);
      print TR "$qw[1]\t"; 
      undef(@qw);  
    }
    print TR "\n"; 
    next; 
  }
  
  $masterfl=1;
  for($i=2;$i<$no;$i++)
  {
    if($arr[$i]!~/NA/ && $arr[$i]!=0) { $masterfl++; }
  }
  if($masterfl/$no>0.1)
  {
    for($i=1;$i<$no;$i++)
    {
      print TR "$arr[$i]\t";
    }
    print TR "\n"; 
    $cc++; 
  }
  undef $masterfl;
  undef $no;
  undef(@arr);
}
close(TR); 
close(TM);

print "$model[$modt] $cc\n"; 

system "rm TMP.txt"; 

$R->run(q'rm(inputData)');
$R->run(q'inputData=read.table("inputData.txt",header=T)');

system "rm inputData.txt"; 

$R->run(q'lmData=inputData[,which(colSums(abs(inputData))!=0)]');
$R->run(q'lmData=lmData[,which(colSums(is.na(lmData))<0.9*nrow(lmData))]');
$R->run(q'lmData=lmData[,which(colSums(lmData==0)<0.9*nrow(lmData))]');
$R->run(q'lmData=lmData[,which(apply(lmData, 2, var, na.rm=TRUE) != 0)]');
$R->run(q'lmData=as.data.frame(scale(lmData))');

$R->run(q'inputData=inputData[,which(colSums(abs(inputData),na.rm=T)!=0)]');
$R->run(q'inputData=inputData[,which(colSums(is.na(inputData))<0.9*nrow(inputData))]');
$R->run(q'inputData=inputData[,which(colSums(inputData==0,na.rm=T)<0.9*nrow(inputData))]');
$R->run(q'inputData=inputData[,which(apply(inputData, 2, var, na.rm=TRUE) != 0)]');

##Ridge regression on unfiltered data 

$R->run(q'linRidgeMod <- linearRidge(DM_YPD ~ ., data = lmData)');
$R->run(q'sum=summary(linRidgeMod)'); 
$R->run(q'npc=sum$chosen.nPCs');
$npc=$R->get('npc');
$R->run(qq'len=length(sum\$summaries\$summary$npc\$coefficients\[,5\])');
$len=$R->get('len');
##print "$npc $len\n"; 
$mod1ln=0; 
$modvar1="";
for($kk=2;$kk<=$len;$kk++)
{
  $R->run(qq'pval=sum\$summaries\$summary$npc\$coefficients\[,5\][[$kk]]');
  $pval=$R->get('pval');
  $R->run(qq'name=names(sum\$summaries\$summary$npc\$coefficients\[,5\])[[$kk]]');
  $name=$R->get('name');
  if($pval<0.05)
  {
    ##print "$name $pval\n"; 
    $modvar1.=$name."+"; 
    $mod1ln++; 
  } 
  undef $name; 
  undef $pval; 
}
$R->run(q'rm(len)');
$R->run(q'rm(npc)');
$R->run(q'rm(sum)');
$R->run(q'rm(linRidgeMod)');
$R->run(q'rm(lmData)');

undef $len; 
undef $npc; 
  
chop($modvar1);

$mod3ln=0;
if($mod1ln==1) { $modvar3=$modvar1; }
else
{
$R->run(qq'modlm <- lm(DM_YPD ~ $modvar1 , data=inputData)');
$R->run(q'modaic=ols_step_both_aic(modlm)');
$ref1=$R->get(q'modaic$predictors');
$ref2=$R->get(q'modaic$method');
@arref1=@{$ref1};
@arref2=@{$ref2};
$len=@arref1;

$mod3ln=0;
$modvar3="";
for($kk=0;$kk<$len;$kk++)
{
  ## print "$arref1[$kk] $arref2[$kk]\n";
  if($arref2[$kk]=~/addition/)
  {
    ##print "$name $pval\n";
    if($arref1[$kk]!~/NA/ && $arref1[$kk]!~/\[/)
    {
      $varlist{$arref1[$kk]}=1;
    }
  }
  if($arref2[$kk]=~/removal/)
  {
    if($arref1[$kk]!~/NA/ && $arref1[$kk]!~/\[/)
    {
      if(exists $varlist{$arref1[$kk]})
      {
        $varlist{$arref1[$kk]}=0;
      }
    }
  }
}
foreach $key (keys %varlist)
{
  if($varlist{$key}==1)
  {
    $modvar3.=$key."+";
    $mod3ln++;
  }
}
undef %varlist;

$R->run(q'rm(modaic)');
$R->run(q'rm(modlm)');

undef $len;
undef $ref1;
undef $ref2;
undef(@arref1);
undef(@arref2);

chop($modvar3);
}

##$R->run(q'datatab=na.omit(inputData)'); 
$R->run(q'library(mice)');
$R->run(q'set.seed(1234)');
$R->run(q'miceMod <- mice(inputData[, !names(inputData) %in% "DM_YPD"], method="rf",rfPackage = "randomForest")');  # perform mice imputation, based on random forests.
$R->run(q'newdata <- complete(miceMod)');  # generate the completed data.
$R->run(q'DM_YPD = inputData$DM_YPD');
$R->run(q'datatab <- as.data.frame(cbind(DM_YPD,newdata))');  # generate the completed data.
$R->run(q'rm(miceMod)');
$R->run(q'rm(newdata)');
$R->run(q'rm(DM_YPD)');

$R->run(q'len=length(names(datatab))');
$len=$R->get('len'); 
$R->run(q'x = model.matrix(DM_YPD~., datatab)[,-1]'); 
$R->run(q'y = datatab %>% select(DM_YPD) %>% unlist() %>% as.numeric()');
$R->run(q'grid = 10^seq(10, -2, length = 100)');
$R->run(q'cv.out = cv.glmnet(x, y, alpha = 1)');
$R->run(q'bestlam = cv.out$lambda.min'); 
$R->run(q'out = glmnet(x, y, alpha = 1, lambda = grid)'); 
$R->run(qq'lasso_coef = predict(out, type = "coefficients", s = bestlam)[1:$len,]'); 
$R->run(q'impvar=lasso_coef[lasso_coef != 0]'); 
$R->run(q'implen=length(names(impvar))'); 
$implen=$R->get(q'implen');

$mod2ln=0; 
$modvar2="";
for($kk=2;$kk<=$implen;$kk++)
{
  $R->run(qq'varname=names(impvar)[[$kk]]'); 
  $varname=$R->get('varname');
  $modvar2.="$varname+";
  $mod2ln++;
  undef $varname; 
}
chop($modvar2); 

$R->run(q'rm(datatab)');
$R->run(q'rm(len)');
$R->run(q'rm(x)');
$R->run(q'rm(y)');
$R->run(q'rm(grid)');
$R->run(q'rm(cv.out)');
$R->run(q'rm(out)');
$R->run(q'rm(bestlam)');
$R->run(q'rm(lasso.coef)');
$R->run(q'rm(impvar)');
$R->run(q'rm(implen)');

$mod4ln=0;
if($mod2ln==1) { $modvar4=$modvar2; }
else
{
$R->run(qq'modlm <- lm(DM_YPD ~ $modvar2 , data=inputData)');
$R->run(q'modaic=ols_step_both_aic(modlm)');
$ref1=$R->get(q'modaic$predictors');
$ref2=$R->get(q'modaic$method');
@arref1=@{$ref1};
@arref2=@{$ref2};
$len=@arref1;

$mod4ln=0;
$modvar4="";
for($kk=0;$kk<$len;$kk++)
{
  ## print "$arref1[$kk] $arref2[$kk]\n";
  if($arref2[$kk]=~/addition/)
  {
    ##print "$name $pval\n";
    if($arref1[$kk]!~/NA/ && $arref1[$kk]!~/\[/)
    {
      $varlist{$arref1[$kk]}=1;
    }
  }
  if($arref2[$kk]=~/removal/)
  {
    if($arref1[$kk]!~/NA/ && $arref1[$kk]!~/\[/)
    {
      if(exists $varlist{$arref1[$kk]})
      {
        $varlist{$arref1[$kk]}=0;
      }
    }
  }
}
foreach $key (keys %varlist)
{
  if($varlist{$key}==1)
  {
    $modvar4.=$key."+";
    $mod4ln++;
  }
}
undef %varlist;

$R->run(q'rm(modaic)');
$R->run(q'rm(modlm)');

undef $len;
undef $ref1;
undef $ref2;
undef(@arref1);
undef(@arref2);

chop($modvar4);
}

print "Ridge: $modvar1\n";
print "Ridge_aic: $modvar3\n";
print "Lasso: $modvar2\n";
print "Lasso_aic: $modvar4\n";


$rdc=0; $rdpval=0; 
$lsc=0; $lspval=0; 
$rac=0; $rapval=0; 
$lac=0; $lapval=0;
$rfc=0; $rfc2=0;

for($i=0;$i<$iter;$i++)
{
  
  $R->run(q'trainingIndex <- sample(1:nrow(inputData), 0.80*nrow(inputData))'); # indices for 80% training data
  $R->run(q'trainingData <- inputData[trainingIndex,]'); # training data
  $R->run(q'testData <- inputData[-trainingIndex,]'); # test data

  if($mod1ln!=0)
  {
  $R->run(qq'lmmod <- lm(DM_YPD ~ $modvar1 , data=trainingData)'); 
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

  $R->run(q'cb=cbind(actual,preds)');
  $R->run(q'nav=length(which(rowSums(cb)!="NA"))');
  $nav=$R->get('nav');
  $R->run(q'rm(nav)');
  $R->run(q'rm(cb)');
  if($nav>2)
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
  else { $rves=0; $rvp=1; }
  undef $nav;

  $expval=sprintf("%0.2f",$expval*100); 
  print "Ridge $expval $rsq $rves $rvp\n";

  $R->run(q'rm(full)');
  $R->run(q'rm(lmmod)');
  $R->run(q'rm(expv)');
  $R->run(q'rm(actual)');
  $R->run(q'rm(preds)');
  $R->run(q'rm(rss)');
  $R->run(q'rm(tss)');
  $R->run(q'rm(rsq)');

  $rdexpli[$rdc]=$expval;
  $rdrsqli[$rdc]=$rsq;
  $rdcorli[$rdc]=$rves;
  $rdrmse[$rdc]=$rmse;
  $rdc++;
  if($rvp<0.05) { $rdpval++; }
  
  undef $expval; 
  undef $rsq;
  undef $rves;
  undef $rvp; 
  undef $rmse; 
  }

  if($mod2ln!=0)
  {
  $R->run(qq'lmmod <- lm(DM_YPD ~ $modvar2, data=trainingData)');  

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

  $R->run(q'cb=cbind(actual,preds)');
  $R->run(q'nav=length(which(rowSums(cb)!="NA"))');
  $nav=$R->get('nav');
  $R->run(q'rm(nav)');
  $R->run(q'rm(cb)');
  if($nav>2)
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
  else { $rves=0; $rvp=1; }
  undef $nav;

  $expval=sprintf("%0.2f",$expval*100);
  print "Lasso $expval $rsq $rves $rvp\n";

  $R->run(q'rm(full)');
  $R->run(q'rm(lmmod)');
  $R->run(q'rm(expv)');
  $R->run(q'rm(actual)');
  $R->run(q'rm(preds)');
  $R->run(q'rm(rss)');
  $R->run(q'rm(tss)');
  $R->run(q'rm(rsq)');
 
  $lsexpli[$lsc]=$expval;
  $lsrsqli[$lsc]=$rsq;
  $lscorli[$lsc]=$rves;
  $lsrmse[$lsc]=$rmse;
  $lsc++;
  if($rvp<0.05) { $lspval++; }

  undef $expval;
  undef $rsq;
  undef $rves;
  undef $rvp;
  undef $rmse;
  }
  if($mod3ln!=0)
  {
  $R->run(qq'lmmod <- lm(DM_YPD ~ $modvar3, data=trainingData)');
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

  $R->run(q'cb=cbind(actual,preds)');
  $R->run(q'nav=length(which(rowSums(cb)!="NA"))');
  $nav=$R->get('nav');
  $R->run(q'rm(nav)');
  $R->run(q'rm(cb)');
  if($nav>2)
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
  else { $rves=0; $rvp=1; }
  undef $nav;

  $expval=sprintf("%0.2f",$expval*100);
  print "Ridge_aic $expval $rsq $rves $rvp\n";

  $R->run(q'rm(full)');
  $R->run(q'rm(lmmod)');
  $R->run(q'rm(expv)');
  $R->run(q'rm(actual)');
  $R->run(q'rm(preds)');
  $R->run(q'rm(rss)');
  $R->run(q'rm(tss)');
  $R->run(q'rm(rsq)');

  $raexpli[$rac]=$expval;
  $rarsqli[$rac]=$rsq;
  $racorli[$rac]=$rves;
  $rarmse[$rac]=$rmse;
  $rac++;
  if($rvp<0.05) { $rapval++; }

  undef $expval;
  undef $rsq;
  undef $rves;
  undef $rvp;
  undef $rmse;
  }

  if($mod4ln!=0)
  {
  $R->run(qq'lmmod <- lm(DM_YPD ~ $modvar4, data=trainingData)');
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

  $R->run(q'cb=cbind(actual,preds)');
  $R->run(q'nav=length(which(rowSums(cb)!="NA"))');
  $nav=$R->get('nav');
  $R->run(q'rm(nav)');
  $R->run(q'rm(cb)');
  if($nav>2)
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
  else { $rves=0; $rvp=1; }
  undef $nav;

  $expval=sprintf("%0.2f",$expval*100);
  print "Lasso_aic $expval $rsq $rves $rvp\n";

  $R->run(q'rm(full)');
  $R->run(q'rm(lmmod)');
  $R->run(q'rm(expv)');
  $R->run(q'rm(actual)');
  $R->run(q'rm(preds)');
  $R->run(q'rm(rss)');
  $R->run(q'rm(tss)');
  $R->run(q'rm(rsq)');

  $laexpli[$lac]=$expval;
  $larsqli[$lac]=$rsq;
  $lacorli[$lac]=$rves;
  $larmse[$lac]=$rmse;
  $lac++;
  if($rvp<0.05) { $lapval++; }

  undef $expval;
  undef $rsq;
  undef $rves;
  undef $rvp;
  undef $rmse;
  }
 
  ## Random Forest 

  $R->run(q'rf=rfsrc(DM_YPD~.,data=trainingData,na.action="na.impute")');
  $R->run(q'err.rate1=rf$err.rate');
  $R->run(q'per.var1=round(100 * (1 - err.rate1[length(err.rate1)]/var(rf$yvar,na.rm=T)),2)');
  $rf_perv=$R->get('per.var1');

  $R->run(q'pred=predict(rf,newdata=testData)');
  $R->run(q'err.rate2=pred$err.rate');
  $R->run(q'per.var2=round(100 * (1 - err.rate2[length(err.rate2)]/var(pred$yvar,na.rm=T)),2)');
  $rf_pred=$R->get('per.var2');

  if($rf_perv!~/NA/ && $rf_pred!~/NA/)
  {
     $rfexp[$rfc]=$rf_perv;
     $rfpred[$rfc]=$rf_pred;
     $rfc++;
  }

  print "RandomForestSRC_wiImput\t$rf_perv\t$rf_pred\n";
  undef $rf_perv;
  undef $rf_pred;

  $R->run(q'rm(rf)');
  $R->run(q'rm(err.rate1)');
  $R->run(q'rm(per.var1)');
  $R->run(q'rm(err.rate2)');
  $R->run(q'rm(per.var2)');
  $R->run(q'rm(pred)');

  $R->run(q'rf=rfsrc(DM_YPD~.,data=trainingData)');
  $R->run(q'err.rate1=rf$err.rate');
  $R->run(q'per.var1=round(100 * (1 - err.rate1[length(err.rate1)]/var(rf$yvar,na.rm=T)),2)');
  $rf_perv=$R->get('per.var1');

  $R->run(q'pred=predict(rf,newdata=testData)');
  $R->run(q'err.rate2=pred$err.rate');
  $R->run(q'per.var2=round(100 * (1 - err.rate2[length(err.rate2)]/var(pred$yvar,na.rm=T)),2)');
  $rf_pred=$R->get('per.var2');

  if($rf_perv!~/NA/ && $rf_pred!~/NA/)
  {
     $rfexp2[$rfc2]=$rf_perv;
     $rfpred2[$rfc2]=$rf_pred;
     $rfc2++;
  }

  print "RandomForestSRC_woImput\t$rf_perv\t$rf_pred\n";
  undef $rf_perv;
  undef $rf_pred;  

  $R->run(q'rm(rf)');
  $R->run(q'rm(err.rate1)');
  $R->run(q'rm(per.var1)');
  $R->run(q'rm(err.rate2)');
  $R->run(q'rm(per.var2)');
  $R->run(q'rm(pred)');

  $R->run(q'rm(testData)');
  $R->run(q'rm(trainingData)');
  $R->run(q'rm(trainingIndex)');
}

$R->run(q'rm(inputData)');
$R->run(q'rm(tab_filt2)');
$R->run(q'rm(tab_filt1)');
$R->run(q'rm(tab)');

$p1=0; $sp1=0; $p2=0; $sp2=0;
$s1=0; $s2=0; $s3=0; $s4=0; $s5=0; $s6=0;  
$sd1=0; $sd2=0; $sd3=0; $sd4=0; $sd5=0; $sd6=0;

$ra1=0; $ra2=0; $ra3=0; $ra4=0; $ra5=0; $ra6=0; $ra7=0; $ra8=0;
$la1=0; $la2=0; $la3=0; $la4=0; $la5=0; $la6=0; $la7=0; $la8=0;


$rm1=0; $rm2=0; $rm3=0; $rm4=0;
$rd1=0; $rd2=0; $rd3=0; $rd4=0;

for($i=0;$i<$rdc;$i++)
{
  $s1+=($rdrsqli[$i]);
  $sd1+=($rdrsqli[$i])**2;
  $s2+=($rdcorli[$i]);
  $sd2+=($rdcorli[$i])**2;
  $p1+=($rdexpli[$i]);
  $sp1+=($rdexpli[$i])**2;
  $s5+=($rdrmse[$i]);
  $sd5+=($rdrmse[$i])**2;
}
for($i=0;$i<$lsc;$i++)
{
  $s3+=($lsrsqli[$i]);
  $sd3+=($lsrsqli[$i])**2;
  $s4+=($lscorli[$i]);
  $sd4+=($lscorli[$i])**2;
  $p2+=($lsexpli[$i]);
  $sp2+=($lsexpli[$i])**2;
  $s6+=($lsrmse[$i]);
  $sd6+=($lsrmse[$i])**2;
}

for($i=0;$i<$rac;$i++)
{
  $ra1+=($rarsqli[$i]);
  $ra2+=($rarsqli[$i])**2;
  $ra3+=($racorli[$i]);
  $ra4+=($racorli[$i])**2;
  $ra5+=($raexpli[$i]);
  $ra6+=($raexpli[$i])**2;
  $ra7+=($rarmse[$i]);
  $ra8+=($rarmse[$i])**2;
}
for($i=0;$i<$lac;$i++)
{
  $la1+=($larsqli[$i]);
  $la2+=($larsqli[$i])**2;
  $la3+=($lacorli[$i]);
  $la4+=($lacorli[$i])**2;
  $la5+=($laexpli[$i]);
  $la6+=($laexpli[$i])**2;
  $la7+=($larmse[$i]);
  $la8+=($larmse[$i])**2;
}

for($i=0;$i<$rfc;$i++)
{
  $rm1+=($rfexp[$i]);
  $rd1+=($rfexp[$i])**2;
  $rm2+=($rfpred[$i]);
  $rd2+=($rfpred[$i])**2;
}

for($i=0;$i<$rfc2;$i++)
{
  $rm3+=($rfexp2[$i]);
  $rd3+=($rfexp2[$i])**2;
  $rm4+=($rfpred2[$i]);
  $rd4+=($rfpred2[$i])**2;
}

if($rdc!=0)
{
$sd1=sprintf("%0.2f",sqrt(($sd1/$rdc)-(($s1/$rdc)**2)));
$sd2=sprintf("%0.2f",sqrt(($sd2/$rdc)-(($s2/$rdc)**2)));
$sp1=sprintf("%0.2f",sqrt(($sp1/$rdc)-(($p1/$rdc)**2)));
$sd5=sprintf("%0.2f",sqrt(($sd5/$rdc)-(($s5/$rdc)**2)));

$s1=sprintf("%0.2f",$s1/$rdc);
$s2=sprintf("%0.2f",$s2/$rdc);
$p1=sprintf("%0.2f",$p1/$rdc);
$s5=sprintf("%0.2f",$s5/$rdc);

$rdpval=sprintf("%0.2f",$rdpval/$rdc*100);
}

if($lsc!=0)
{
$sd3=sprintf("%0.2f",sqrt(($sd3/$lsc)-(($s3/$lsc)**2)));
$sd4=sprintf("%0.2f",sqrt(($sd4/$lsc)-(($s4/$lsc)**2)));
$sp2=sprintf("%0.2f",sqrt(($sp2/$lsc)-(($p2/$lsc)**2)));
$sd6=sprintf("%0.2f",sqrt(($sd6/$lsc)-(($s6/$lsc)**2)));

$s3=sprintf("%0.2f",$s3/$lsc);
$s4=sprintf("%0.2f",$s4/$lsc);
$p2=sprintf("%0.2f",$p2/$lsc);
$s6=sprintf("%0.2f",$s6/$lsc);

$lspval=sprintf("%0.2f",$lspval/$lsc*100);
}

if($rac!=0)
{
$ra2=sprintf("%0.2f",sqrt(($ra2/$rac)-(($ra1/$rac)**2)));
$ra4=sprintf("%0.2f",sqrt(($ra4/$rac)-(($ra3/$rac)**2)));
$ra6=sprintf("%0.2f",sqrt(($ra6/$rac)-(($ra5/$rac)**2)));
$ra8=sprintf("%0.2f",sqrt(($ra8/$rac)-(($ra7/$rac)**2)));

$ra1=sprintf("%0.2f",$ra1/$rac);
$ra3=sprintf("%0.2f",$ra3/$rac);
$ra5=sprintf("%0.2f",$ra5/$rac);
$ra7=sprintf("%0.2f",$ra7/$rac);

$rapval=sprintf("%0.2f",$rapval/$rac*100);
}

if($lac!=0)
{
$la2=sprintf("%0.2f",sqrt(($la2/$lac)-(($la1/$lac)**2)));
$la4=sprintf("%0.2f",sqrt(($la4/$lac)-(($la3/$lac)**2)));
$la6=sprintf("%0.2f",sqrt(($la6/$lac)-(($la5/$lac)**2)));
$la8=sprintf("%0.2f",sqrt(($la8/$lac)-(($la7/$lac)**2)));

$la1=sprintf("%0.2f",$la1/$lac);
$la3=sprintf("%0.2f",$la3/$lac);
$la5=sprintf("%0.2f",$la5/$lac);
$la7=sprintf("%0.2f",$la7/$lac);

$lapval=sprintf("%0.2f",$lapval/$lac*100);
}

if($rfc!=0)
{
$rd1=sprintf("%0.2f",sqrt(($rd1/$rfc)-(($rm1/$rfc)**2)));
$rd2=sprintf("%0.2f",sqrt(($rd2/$rfc)-(($rm2/$rfc)**2)));
$rm1=sprintf("%0.2f",$rm1/$rfc);
$rm2=sprintf("%0.2f",$rm2/$rfc);
}
if($rfc2!=0)
{
$rd3=sprintf("%0.2f",sqrt(($rd3/$rfc2)-(($rm3/$rfc2)**2)));
$rd4=sprintf("%0.2f",sqrt(($rd4/$rfc2)-(($rm4/$rfc2)**2)));
$rm3=sprintf("%0.2f",$rm3/$rfc2);
$rm4=sprintf("%0.2f",$rm4/$rfc2);
}

if($dtc==0) { print WR "Unfiltered data $iter iterations\n------------------------------------------------------------------------------------------------------\n"; }

if($dtc==1) { print WR "Filtered data $iter iterations\n------------------------------------------------------------------------------------------------------\n"; }
if($dtc==2) { print WR "Filtered2 data $iter iterations\n------------------------------------------------------------------------------------------------------\n"; }

print WR "Ridge    : $modvar1\n";
print WR "Ridge_aic: $modvar3\n";
print WR "Lasso    : $modvar2\n";
print WR "Lasso_aic: $modvar4\n\n";

print WR "Ridge regression: %VarExpl $p1±$sp1  Pred.Rsq $s1±$sd1 RMSE: $s5±$sd5 Cor: $s2±$sd2  Sig. Pval: $rdpval %\n";
print WR "Ridge aic       : %VarExpl $ra5±$ra6  Pred.Rsq $ra1±$ra2 RMSE: $ra7±$ra8 Cor: $ra3±$ra4  Sig. Pval: $rapval %\n";
print WR "Lasso regression: %VarExpl $p2±$sp2  Pred.Rsq $s3±$sd3 RMSE: $s6±$sd6 Cor: $s4±$sd4  Sig. Pval: $lspval %\n";
print WR "Lasso aic       : %VarExpl $la5±$la6  Pred.Rsq $la1±$la2 RMSE: $la7±$la8 Cor: $la3±$la4  Sig. Pval: $lapval %\n";
print WR "RF_wiImp        : %VarExpl $rm1±$rd1  PredVar $rm2±$rd2\n";
print WR "RF_woImput      : %VarExpl $rm3±$rd3  PredVar $rm4±$rd4\n";
print WR "------------------------------------------------------------------------------------------------------\n";

$p1=sprintf("%0.2f",$p1/100);
$p2=sprintf("%0.2f",$p2/100);
$sp1=sprintf("%0.2f",$sp1/100);
$sp2=sprintf("%0.2f",$sp2/100);

$rm1=sprintf("%0.2f",$rm1/100);
$rm2=sprintf("%0.2f",$rm2/100);
$rd1=sprintf("%0.2f",$rd1/100);
$rd2=sprintf("%0.2f",$rd2/100);

$rm3=sprintf("%0.2f",$rm3/100);
$rm4=sprintf("%0.2f",$rm4/100);
$rd3=sprintf("%0.2f",$rd3/100);
$rd4=sprintf("%0.2f",$rd4/100);

$ra5=sprintf("%0.2f",$ra5/100);
$ra6=sprintf("%0.2f",$ra6/100);
$la5=sprintf("%0.2f",$la5/100);
$la6=sprintf("%0.2f",$la6/100);

print WR2 "$p1\t$ra5\t$p2\t$la5\t$rm1\t$rm3\t";
print WR3 "$s1\t$ra1\t$s3\t$la1\t$rm2\t$rm4\t";
$res[$dtc][0]=$sp1; $res[$dtc][1]=$ra6; $res[$dtc][2]=$sp2;
$res[$dtc][3]=$la6; $res[$dtc][4]=$rd1; $res[$dtc][5]=$rd3;

$res[$dtc][6]=$sd1; $res[$dtc][7]=$ra2; $res[$dtc][8]=$sd3;
$res[$dtc][9]=$la2; $res[$dtc][10]=$rd2; $res[$dtc][11]=$rd4;

undef $mod1ln;
undef $mod2ln;
undef $mod3ln;
undef $mod4ln;
undef $modvar1;
undef $modvar2;
undef $modvar3;
undef $modvar4;

undef(@rfexp);
undef(@rfpred);
undef(@rfexp2);
undef(@rfpred2);
undef(@rdrsqli);
undef(@rdcorli);
undef(@rdexpli);
undef(@rdrmse);
undef(@rarsqli);
undef(@racorli);
undef(@raexpli);
undef(@rarmse);
undef(@lsrsqli);
undef(@lscorli);
undef(@lsexpli);
undef(@lsrmse);
undef(@larsqli);
undef(@lacorli);
undef(@laexpli);
undef(@larmse);

}
for($u=0;$u<$dtc;$u++)
{
  for($t=0;$t<6;$t++)
  {
    print WR2 "$res[$u][$t]\t"; 
  }
}
print WR2 "\n";

for($u=0;$u<$dtc;$u++)
{
  for($t=6;$t<12;$t++)
  {
    print WR3 "$res[$u][$t]\t";
  }
}
print WR3 "\n";

undef(@res);

}
close(WR3); 
close(WR2);
close(WR);

undef(@model);

$R->run('rm(whtab)');
$R->stopR; 
