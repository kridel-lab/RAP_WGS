# Rapid Autopsy project 

Input: BAM files provided by TCAG for 1 normal, 3 diagnostic samples and 17 autopsy samples 

## Part 1: Somatic variant calling 
- Ran Mutect2 which is optimized for identifying low frequency somatic variants scripts here: 
- Ran Manta to identify somatic structural variants for each sample, script here: 
- Ran Strelka to idenify somatic SNVs for each sample, script here: 
- Merge variants from Mutect2 and Strelka and apply further soft filters 
  - minimum depth of 60
  - no indels 
  - no chromosome X and Y 
  - no population variants (Gnomad) 
  - minimum VAF of 0.1 
  - remove Synonymous variants 
  - remove variants overlapping blacklisted regions (ENCODE) 

## Part 2: Somatic copy number calling 
- Ran TitanCNA snakemake files (minor adjustments were done to TitanCNA scripts so make sure to use edited scripts if need to in the future, edited scripts here: 

## Part 3: Analysis of variants, preliminary findings, basic summaries and visualizations 
- basic summary of variants per sample 
- founder mutations versus sample specific variants 
- preliminary analysis of mutation signatures 

## Part 4: Phylogenetic analysis 
- Treeomics 
- PhyloWGS 
- Palimpsest

