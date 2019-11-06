# Rapid Autopsy project 

Input: BAM files provided by TCAG for 1 normal, 3 diagnostic samples and 17 autopsy samples 

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
  - remove Synonymous variants 
  - remove variants overlapping blacklisted regions (ENCODE) 

## Part 2: Somatic copy number calling
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_2_Somatic_copy_number_calling 
- Ran TitanCNA snakemake files (minor adjustments were done to TitanCNA scripts so make sure to use edited scripts if need to in the future, edited scripts here: 

## Part 3: Analysis of variants, preliminary findings, basic summaries and visualizations
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_3_Analysis_variants_preliminary_findings
- basic summary of variants per sample 
- founder mutations versus sample specific variants 
- preliminary analysis of mutation signatures 

## Part 4: Phylogeny analysi 
https://github.com/kridel-lab/RAP_WGS/tree/master/Part_4_Phylogeny_analysis
- Treeomics 
- PhyloWGS 
- Palimpsest

