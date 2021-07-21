# Rapid Autopsy project

Input: BAM files provided by TCAG for 3 normal, 3 diagnostic samples and 24 autopsy samples

## Part 1: Somatic variant calling
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_1_Somatic_variant_calling
- Ran Mutect2 which is optimized for identifying low frequency somatic variants
- Ran Manta to identify somatic structural variants for each sample
- Ran Strelka to idenify somatic SNVs for each sample
- Merge variants from Mutect2 and Strelka and apply further soft filters
  - minimum depth of 60
  - no indels
  - no chromosome X and Y
  - no population variants (Gnomad)
  - minimum VAF of 0.1
  - remove variants overlapping blacklisted regions (ENCODE)

## Part 2: Somatic copy number calling
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_2_Somatic_copy_number_calling
- Ran Sequenza to obtain overall copy number status, purity, ploidy and major minor CN status

## Part 3: Analysis of variants, preliminary findings, basic summaries and visualizations
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_3_Analysis_variants_preliminary_findings
- basic summary of variants per sample
- scripts for figures

## Part 4: Phylogeny analysis
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_4_Phylogeny_analysis
- Treeomics
- Pyclone-VI
- Pairtree

## Part 5: ctDNA analysis
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_5_ctDNA
- Ran ConsensusCruncher to process UMIs and generate BAM files
- Ran Mutect2 in tumour/normal mode to detect SNVs in plasma ctDNA (using WGS germline samples as controls)

## Notes
- To easily load processed mutation and CNA data to work with on the cluster, use this script: https://github.com/kridel-lab/RAP_WGS/blob/master/config-file.R 
- Full order of scripts in which they were run is here: https://github.com/kridel-lab/RAP_WGS/blob/master/order-scripts.sh
