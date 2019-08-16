python2 phylowgs/parser/create_phylowgs_inputs.py --vcf-type S1=mutect_smchet S1=Documents/RAP_ANALYSIS/LY_RAP_0003_Dia_FoT_05_mutect_patient_results_all.vcf_filtered.vcf --cnvs S1=Documents/RAP_ANALYSIS/optimalClusterSolution/tumor_sample_20_cluster1.segs.txt.parsed.txt 



#test version 


less Documents/RAP_ANALYSIS/LY_RAP_0003_Dia_FoT_05_mutect_patient_results_all.vcf_filtered.vcf | head -3000 > test.vcf
less Documents/RAP_ANALYSIS/optimalClusterSolution/tumor_sample_20_cluster1.segs.txt.parsed.txt | head -1000 > test_cnv.txt

python2 phylowgs/parser/create_phylowgs_inputs.py --vcf-type S1=mutect_smchet S1=test.vcf --cnvs S1=test_cnv.txt --output-cnvs test_cnv_data.txt --output-variants test_ssm_data.txt 


#run phylowgs
phylowgs=phylowgs/multievolve.py 
python2 $phylowgs --num-chains 4 --ssms test_ssm_data.txt --cnvs test_cnv_data.txt 

#To write JSON results, please run 
python2 phylowgs/write_results.py run_name chains/trees.zip run_name.summ.json.gz run_name.muts.json.gz run_name.mutass.zip