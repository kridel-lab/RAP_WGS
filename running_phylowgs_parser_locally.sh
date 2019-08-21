module load python2

less LY_RAP_0003_Aut_FzT_16_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input.vcf.gz | head -1000 > test.vcf
tum1=test.vcf
cnv1=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_20_cluster1.segs.txt.parsed.txt_plusone.txt
parser=/cluster/home/kisaev/phylowgs/parser/create_phylowgs_inputs.py 

less LY_RAP_0003_Aut_FzT_17_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input.vcf.gz | head -1000 > test2.vcf
tum2=test2.vcf
cnv2=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_15_cluster1.segs.txt.parsed.txt_plusone.txt

python2 $parser --vcf-type S2=mutect_smchet S2=$tum2 --cnvs S2=$cnv2
python2 $parser --vcf-type S2=mutect_smchet S2=$tum1 --cnvs S2=$cnv1
python2 $parser --vcf-type S1=mutect_smchet --vcf-type S2=mutect_smchet S1=$tum1 S2=$tum2 --cnvs S1=$cnv1 --cnvs S2=$cnv2

less LY_RAP_0003_Dia_FoT_03_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input.vcf.gz | head -250 > test_auto.vcf
less LY_RAP_0003_Aut_FzT_15_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input.vcf.gz | head -250 > test_frozen.vcf
tum1=test.vcf
python2 $parser --vcf-type S1=mutect_smchet S1=test_frozen.vcf --regions=all 




#test version 

less Documents/RAP_ANALYSIS/LY_RAP_0003_Dia_FoT_05_mutect_patient_results_all.vcf_filtered.vcf | head -3000 > test.vcf
less Documents/RAP_ANALYSIS/optimalClusterSolution/tumor_sample_20_cluster1.segs.txt.parsed.txt | head -1000 > test_cnv.txt

python2 phylowgs/parser/create_phylowgs_inputs.py --vcf-type S1=mutect_smchet S1=test.vcf --cnvs S1=test_cnv.txt --output-cnvs test_cnv_data.txt --output-variants test_ssm_data.txt 

#run phylowgs
phylowgs=phylowgs/multievolve.py 
python2 $phylowgs --num-chains 4 --ssms test_ssm_data.txt --cnvs test_cnv_data.txt 

#To write JSON results, please run 
python2 phylowgs/write_results.py run_name chains/trees.zip run_name.summ.json.gz run_name.muts.json.gz run_name.mutass.zip