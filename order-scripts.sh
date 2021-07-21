#Last updated: July 14th, 2021
#By: Karin Isaev

#READ ME
#This file contains the order of scripts run by KI for all components of the Rapid autopsy project
#if any questions please email: karin.isaev@gmail.com

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[1] prep raw data files
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#----navigate to main directory-------------------------------------------------

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#collect all sample BAM files into file
#ls */*/*recal.cram* > all_bam_files_raw.txt #with some manual clean up to only include tumour cram files

less all_bam_files_raw.txt | wc -l #this file contains the path to the 27 tumour CRAM files

#collect all control samples
less all_control_samples.txt | wc -l #this file contains the path to the 3 control sample files

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[2] run strelka and manta on WGS tumour samples with normal controls
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#----run Manta------------------------------------------------------------------

#Tool used to call translocations and large insertions and deletions

#[1] Run Manta
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_001_manta.sh

#[2] Next process and summmarize manta results
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_007_processing_manta_results.sh

#[3] Combine all SV calls from all samples into one file
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_008_processing_manta_results.R

#[4] Clean up SV calls and sample names so it's more readable
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_009_processing_manta_results.R

#----run Strelka----------------------------------------------------------------

#Tool uses to call SNVs and indels

#[1] call variants
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_002_strelka.sh

#[2] select variants that pass filters
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_003_SelectVariants_post_strelka.sh

#[3] select variants that pass filters (indels)
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/Strelka_003_SelectVariants_post_strelka_indels.sh

#[4] prepare bed files to merge with Mutect2
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/strelka_prepare_bed_files.R

#[5] prepare bed files to merge with Mutect2 (indels)
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Strelka/strelka_prepare_bed_files_indels.R

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[3] run Mutect2 on WGS tumour samples with normal controls
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. first split normal samples by chromosomes
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_001_start_running_bam_splits_normal.sh

#2. split bam files by chromosomes
#requested temp dirctory for running mutect2
#can do this by emailing Zhibin and he can give you temporary project
#space for two weeks for example
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_001_start_running_bam_splits.sh

#3. run mutect2 mutation calling on individual chr based bam files
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_001_running_mutect2.sh

#4. run getpileup on tumour samples
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_002A_running_GetPileupSummaries.sh

#5. run getpileup on normal samples
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_002A_running_GetPileupSummaries_control.sh

#6. get contamination tables
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_003A_running_CalculateContamination.sh

#7. run filter variant calls
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_004A_FilterVariants.sh

#8. run gather VCFs from across chromosomes
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_005A_running_gatherVCFs.sh

#9. select variants from combined VCF files to only include variants that passed
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_006_SelectVariants_post_mutect2.sh

#10. select variants from combined VCF files to only include variants that passed
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_006B_SelectVariants_post_mutect2_save_INDELS.sh

#11. soft filtering of mutect2 variants in R before merging with strelka
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_007_soft_filtering_mutect2_filtered_variants.R

#12. soft filtering of mutect2 variants in R before merging with strelka (indels)
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Mutect2/pipeline_007_soft_filtering_mutect2_filtered_variants_indels.R

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[4] merge variants from Strelka and Mutect2
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. get list of mutations called by both tools
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_001_evaluating_overlap_between_callers.R

#2. get list of mutations called by both tools (indels)
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_001_evaluating_overlap_between_callers_indels.R

#3. use original mutect2 VCF files and keep only the mutations found by both callers
#run annovar
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_002_filter_variants_only_merged_run_annovar.sh

#4. use original mutect2 VCF files and keep only the mutations found by both callers (indels)
#run annovar
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_002_filter_variants_only_merged_run_annovar_indels.sh

#5. remove extra AF column from VCF files otherwise won't load properly into R
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003B_remove_extra_AF_column.sh

#6. run annovar annotated VCF files through soft filtering
#vcfR functions don't work in 3.5.0
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003_processing_merged_annovar.R

#7. run annovar annotated VCF files through soft filtering (indels)
#vcfR functions don't work in 3.5.0
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_003_processing_merged_annovar_indels.R

#8. some additional soft filtering and also the script that looks at overlaps
#between our gene mutations and several gene panels
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_004_additional_soft_filters_applied.R

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[5] Run Sequenza
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. Run once (won't need to re-run to re-run Sequenza)
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/001_fastq_fc_wiggle.sh

#2. To run sequenza you actually need BAM files but we got CRAM files uploaded
#from sickkids so these scripts are run to convert CRAM -> BAM
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/001_p001_p002_get_BAM_files.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/002_p001_p002_get_BAM_files_control_samples.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/p003_get_BAM_files.sh

#3. Now run Sequenza steps
#First run the sequenza-utils bam2seqz step (once it's done I commmented out this line and
#reran the same script to run the seqz_binning step)

#sequenza-utils bam2seqz (comment out seqz_binning)
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/002_get_individual_sequenza_files.sh

#sequenza-utils seqz_binning (comment out bam2seqz)
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/002_get_individual_sequenza_files.sh

#4. Run final Sequenza command which is an R script (this produces the final Sequenza output files)
sbatch /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/003_R_script_job_submission.sh

#5. Assemble all Sequenza output together
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/prepare_Sequenza_CNA_input_Palimpsest.R

#6. intersect mutations with CNAs - the output from these files get read in
#by config-file.R to generate the final mutation file
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_prepare_bed_files_to_fillter_VCFs_by_soft_filters_using_sequenza.R
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_prepare_bed_files_to_fillter_VCFs_by_soft_filters_indels_using_sequenza.R

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[6] Run Pyclone-VI
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pyclone_folder=/cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP

#prepare input for bamreadcount
Rscript $pyclone_folder/pyclone-vi_001_prepare_subset_and_all_mutation_files.R

#run bamreadcount
sbatch $pyclone_folder/pyclone-vi_002_Pyclone_get_missing_mutations.sh

#assemble results from bamreadcount into pyclone mutation list
sbatch $pyclone_folder/pyclone-vi_003_Merged_008_get_pyclone_missing_counts.sh

#prepare input for pyclone-vi
Rscript $pyclone_folder/pyclone-vi_004_make_input_files_sequenza.R

#run pyclone-vi all muts
sbatch $pyclone_folder/pyclone-vi_005_run_main_program_all_muts.sh

#evaluate mutation signatures across pyclone clusters
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone-vi_006_pyclone_mutation_signatures_in_clusters.R

#prepare input for mapscape (edges are prepared manually seperatley and locally)
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone-vi_007_prep_for_mapscape.R

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[7] Run Pairtree
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#generate input files for pairtree using pyclone input and output files
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pairtree/001_create_ssm_input_files.R

#run pairtree with default parameters
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pairtree/002_run_Pairtree.sh

#run pairtree on specific tree ID chosen as model
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pairtree/003_run_Pairtree_choose_tree_model.sh

#post pairtree
module purge
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pairtree/004_poss_pairtree_analysis.R

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[8] Run Treeomics
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. make patient specific bed files from which to pull mutations in vcf files for treeomics
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_treeomics_prepare_bed_files_to_fillter_VCFs_by_soft_filters.R

#2. clean up VCF first for each sample to only inlcude mutations that want to include in analysis
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Treeomics/treeomics_001_A_input_VCF_files.sh

#3. change VCF file names in treeomics input to sample origin / organ
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Treeomics/treeomics_001_B_rename_sample_names_in_VCF_file.sh

#4. Run Treeomics on three patients
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Treeomics/treeomics_002_PCG_mutations_only_all.sh

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[9] Run Lymphgen Staudt Classifier
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. set up sample file
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/001_setting_up_sample_annotation_file.R

#2. set up mutation file
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/002_setting_up_mutation_file.R

#3. set up CNA file
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/003_setting_up_CNA_file.R

#4. summarize results from Staudt classifier
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Staudt_subtyping_samples/004_summarize_LymphGen_results.R

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[10] ctDNA analysis
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. Run Consensus Cruncher
#NOTE 1: the pipeline (snakemake) in the RAP directory isn't the most up to date
#pipeline that Gabrielle has worked on. Once her pipeline is fully tested and
#complete you could replace the one I used with hers and re-run the analysis that way
#the way mine is currently set up should also be appropriate to use
#NOTE 2: also note that the final data used for ctDNA includes fastq files
#from two different lanes that were combined more details are available
#in the README file on the cluster

#First combine fastq files from the first upload of data and the
#second upload of data
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/001_combine_fastq_files.sh

#Now set up config file to indicate where the fastq files are
#edit this file locally and push changes /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/config/config_all_samples_combined.json

#Now edit output destination in the snakemake itself and submit script to run
#snakemake pipeline for first steps of consensuscruncher
#edit this file if needed to change output folders:
# /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Consensus_Cruncher_Pipeline_Combined.snakefile

sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/snakemake_submit.sh

#2. Run Picard tools to get coverage across targets
#NOTE: here you can run these scripts on the uncollapsed BAM files and also the
#BAM files produced with the different corrections (SSCS vs DCS for example)

#for collapsed bam files (in the end we use the dcs_sc correction)
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs_sc/001_CollectTargetsPCRMetrics.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/sscs/001_CollectTargetsPCRMetrics.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/ssc_sc/001_CollectTargetsPCRMetrics.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs/001_CollectTargetsPCRMetrics.sh

#3. Run Mutect2 on each type of BAM files
#NOTE: here I ran the scripts related to mutation calling and annotation
#seperatley for each of the four corrections but if you run Gabrielle's latest
#pipeline it will automaitcally run Mutect2 on all 4 corrections... I just
#didn't have that much time and set it up this way but it's definitely not the
#cleanest
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs_sc/002_run_Mutect2.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/sscs/002_run_Mutect2.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/ssc_sc/002_run_Mutect2.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs/002_run_Mutect2.sh

sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs_sc/003_Mutect2_filter_variants.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/sscs/003_Mutect2_filter_variants.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/ssc_sc/003_Mutect2_filter_variants.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs/003_Mutect2_filter_variants.sh

sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs_sc/004_Mutect2_annotate_wAnnovar.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/sscs/004_Mutect2_annotate_wAnnovar.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/ssc_sc/004_Mutect2_annotate_wAnnovar.sh
sbatch /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/dcs/004_Mutect2_annotate_wAnnovar.sh

#4. Combine mutations into one file and clean up output from Annovar
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Mutect2/Collect_annovar_mutations_into_matrix.R

#5. Look at overlap of mutations across different corrections and also with RAP
#samples
module load R/4.0.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Mutect2/compare_ctDNA_to_bulkDNA.R

#6. Summarize mutations called in ctDNA and make a plot


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[11] Scripts for figures in manuscript
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Figure 1+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. prep CNAs for Figure 1 barplot (on the cluster)
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure1_plots/001_CNAs_summary_plots.R

#2. prep SNVs and indels for Figure 1 barplot (on the cluster)
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure1_plots/001_SNVs_indels_summary.R

#3. prep Structural Variants (SVs) for Figure 1 barplot (on the cluster)
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure1_plots/001_SVs_summary_plots.R

#4. compile final plot for Figure 1B summarizing all events (locally using files
#on UHN teams downloaded from previous scripts output, directory is in the scripts)
Rscript /Users/kisaev/github/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure1_plots/Main_Figure_All_Facets.R

#Supp Figure 1A+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure1_plots/001_SNVs_indels_private_vs_purity.R

#Figure 2+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#1. make circos plots (easier to run this script in Rstudio interactively in
#order to save the plot produced by the function)
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure2_plots/001_circos_plots_local.R

#2. Prepare data for Figure 2 D, E, F
#prep SNVs numbers for line plots
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure2_plots/001_SNVs_indels_lineplot.R
#prep CNAs summarize number of genes affected by CNAs in each sample
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure2_plots/001_CNAs_lineplot.R
#prep SVs summarize how many affected in each sample
Rscript /cluster/home/kisaev/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure2_plots/001_SVs_summary_plots_lineplot.R
#combine all together (run locally and download outputs from scripts above to
#make final figure 2)
Rscript /Users/kisaev/github/RAP_WGS/Part_3_Analysis_variants_preliminary_findings/Figure2_plots/Main_Figure_2.R

#DONE
