
#use mutations evaluted by pyclone-VI

#this is what it needs to look like:
#id	name	var_reads	total_reads	var_read_prob
#s0	S_0	54,175,196	1000,1000,1000	0.5,0.5,0.5
#s1	S_1	24,10,98	1000,1000,1000	0.5,0.5,0.5
#s2	S_2	241,327,397	1000,1000,1000	0.5,0.5,0.5

#keep track of ordering of samples by commas

#var_reads: a comma-separated vector of how many whole-number-valued
#genomic reads at this mutation's locus corresponded to the variant allele
#in each tissue sample. Tissue samples may be provided in any order, so long as
#this order is consistent for all mutations. The names of the associated tissue
#samples are given in the .params.json file, detailed below.

#var_read_prob: f a copy-number aberration (CNA) duplicated the reference allele
#in the lineage bearing mutation j prior to j occurring, there will be two
#reference alleles and a single variant allele in all cells bearing j, such that \omega_{js} = 0.3333
