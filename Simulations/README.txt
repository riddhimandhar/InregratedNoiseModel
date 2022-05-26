Codes
======

1. simulate_coop_comp_TFbinding.pl 

   - Simulate TF co-operativity and competitive binding and their impact on gene expression 
     comparison with regulation by single TF and two independent TFs
   - Output files: 
     SingleTF_Timevar_trarate.pdf
     SingleTF_Timevar_mRNA.pdf
     SingleTF_Timevar_prot.pdf
     MultiTF_Timevar_trarate.pdf
     MultiTF_Timevar_mRNA.pdf
     MultiTF_Timevar_prot.pdf
     CooperativeTF_Timevar_trarate.pdf
     CooperativeTF_Timevar_mRNA.pdf
     CooperativeTF_Timevar_prot.pdf
     CompetitiveTF_Timevar_trarate.pdf
     CompetitiveTF_Timevar_mRNA.pdf
     CompetitiveTF_Timevar_prot.pdf
     SingleTF_CelltoCell_var_trarate.pdf
     SingleTF_CelltoCell_var_mRNA.pdf
     SingleTF_CelltoCell_var_prot.pdf
     MultiTF_CelltoCell_var_trarate.pdf
     MultiTF_CelltoCell_var_mRNA.pdf
     MultiTF_CelltoCell_var_prot.pdf
     CooperativeTF_CelltoCell_var_trarate.pdf
     CooperativeTF_CelltoCell_var_mRNA.pdf
     CooperativeTF_CelltoCell_var_prot.pdf
     CompetitiveTF_CelltoCell_var_trarate.pdf
     CompetitiveTF_CelltoCell_var_mRNA.pdf
     CompetitiveTF_CelltoCell_var_prot.pdf
     NoiseData_SingleTFbinding.txt
     NoiseData_MultiTFbinding.txt
     NoiseData_CooperativeTFbinding.txt
     NoiseData_CompetitiveTFbinding.txt

2. simulate_coop_comp_TFbinding_expl.pl
   
   - Explore impact of on-rate and off-rate values on noise in SingleTF, MultiTF, CooperativeTF and CompetitiveTF binding 

3. process_results_OnOffRate_expl.pl 

   - summarize results obtained in step 2 

4. Parameter_exploration/expl_par_sim_compTFs.pl 

   - Explore the impact of variation in parameter values of the model on noise in case of competitve TF binding 

5. Parameter_exploration/expl_par_sim_coopTFs.pl 

   - Explore the impact of variation in parameter values of the model on noise in case of cooperative TF binding 

6. Parameter_exploration/expl2_par_sim_compTFs.pl 

   - Explore the impact of variation in number of TFs and off rate reduction in the model on noise in case of cooperative TF binding 

7. summarize_all_expl.pl 

   - Summarize results of all parameter exploration for all types of binding 
   - Output files: 
	Stats_CooperativeTFbinding_Noise_data.txt
	Stats_CooperativeTFbinding_BurstSizeFreq.txt
	Stats2_CooperativeTFbinding_Noise_data.txt
	Stats2_CooperativeTFbinding_BurstSizeFreq.txt
	Stats_CompetitiveTFbinding_Noise_data.txt
	Stats_CompetitiveTFbinding_BurstSizeFreq.txt
