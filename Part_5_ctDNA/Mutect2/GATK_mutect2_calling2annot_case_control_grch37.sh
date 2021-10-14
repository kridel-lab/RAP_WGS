#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=46000M
#SBATCH -t 0-10:00 # Runtime in D-HH:MM
#SBATCH -J dcs_sc
#SBATCH --array=0-2 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

#set -o pipefail
#set -e

work_path=$SLURM_SUBMIT_DIR
##the slurm output will be written to your current directory

gatk_bundle="/cluster/projects/kridelgroup/resources/GATK_bundle/homo_sapiens_v37"
#resources to run the gatk, no slash at the end

file_path="${work_path}/bam_files"
control_path="/cluster/projects/kridelgroup/LIBERATE/RAP_Ting_test_2021/germline_control"
#path to the fq files, no slash at the end
#file=$file_path/*.fastq.gz
##use the prefix of the fq file will be used as the output prefix

module load bwa/0.7.15
module load python/3.4.3
module load gatk/4.1.8.1
module load samtools/1.9
module load annovar/20180416
module load tabix
module load vt/0.577

##output directory

if [ ! -d "${work_path}/VCFs" ]; then
        mkdir -p "${work_path}/VCFs"
fi

if [ ! -d "${work_path}/coverage" ]; then
        mkdir -p "${work_path}/coverage"
fi

if [ ! -d "${work_path}/Annovar_annot" ]; then
        mkdir -p "${work_path}/Annovar_annot"
fi

ls ${file_path}/*.bam  > ${work_path}/bam_files_list.txt
#create a list to save the bam files

##assign each sample into different ARRAYID
samples=${work_path}/bam_files_list.txt
echo $sample
prefixs=($(cat $samples))
echo $prefixs
input_bam=${prefixs[${SLURM_ARRAY_TASK_ID}]}

##get the tumor samplename to input into the MuTect2
tumor_samplename=($(samtools view -H ${input_bam} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))

#coordinates of probes (coding and non-coding) used as provided by IDT and Robert
amplicon_interval_list=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/picard_tools_amps_input.bed
#coordinates of targets (coding and non-coding) used as provided by IDT and Robert
targets_interval_list=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/picard_tools_targets_input.bed

###traverse the file_path/.fq folder, then get the file prefix from the r1 reads
###the script will not check if the r2 reads exists or not, please make sure the paired reads is placed in the same ./cleanfq folder

##get the prefix of the bam files, use the prefix to indentify the output
if [[ "$input_bam" =~ ${file_path}/(.*).sorted.*.sorted.bam ]]
 then
	echo "$input_bam"
	prefix="${BASH_REMATCH[1]}"
	echo "$prefix"

##get the file name and sample name of the tumor sample
cont_sample="${tumor_samplename##*_}"
echo $cont_sample
cd $control_path
normal_file=$(ls LY_RAP_${cont_sample}_Ctl*.bam)
echo $normal_file

##get the sample name of the normal sample
normal_samplename=($(samtools view -H ${control_path}/${normal_file} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo $normal_samplename

#somatic mutation calling
#https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136#ref4

##create interval files
gatk BedToIntervalList \
-I "$amplicon_interval_list" \
-O "${work_path}/coverage/${prefix}_amplicon.interval_list" \
-SD "$input_bam" && echo "** bed2intervallist $prefix done **"

gatk BedToIntervalList \
-I "$targets_interval_list" \
-O "${work_path}/coverage/${prefix}_targets.interval_list" \
-SD "$input_bam" && echo "** bed2intervallist targets $prefix done **"

#save new interval lists as variables which will be used as input for final picard function
amps=${work_path}/coverage/${prefix}_amplicon.interval_list
ints=${work_path}/coverage/${prefix}_targets.interval_list

gatk CollectTargetedPcrMetrics \
-I "$input_bam" \
-O "${work_path}/coverage/${prefix}.output_pcr_metrics.txt" \
-R "${gatk_bundle}/human_g1k_v37_decoy.fasta" \
--PER_TARGET_COVERAGE "${work_path}/coverage/${prefix}.per_target_coverage.txt" \
--AMPLICON_INTERVALS "$amps" \
--TARGET_INTERVALS "$ints"  && echo "** CollectTargetedPcrMetrics $prefix done **"

##calling variants
gatk Mutect2 \
-R "${gatk_bundle}/human_g1k_v37_decoy.fasta" \
-I "$input_bam"  \
-I "${control_path}/${normal_file}"  \
-tumor ${tumor_samplename} \
-normal ${normal_samplename} \
-L "${work_path}/coverage/${prefix}_targets.interval_list" \
-O "${work_path}/VCFs/${prefix}.vcf.gz" \
--germline-resource "/cluster/projects/kridelgroup/resources/GATK_bundle/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz" && echo "** $prefix Mutect2 variant calling done **"

#Retrieve contamination and filter SNVs called by MuTect2
##https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136
gatk GetPileupSummaries \
-I "$input_bam" \
-V "/cluster/projects/kridelgroup/resources/GATK_bundle/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz"  \
-L "${work_path}/coverage/${prefix}_targets.interval_list" \
-O "${work_path}/VCFs/${prefix}.pileups.table" && echo "** $prefix pile up summaries done **"

gatk CalculateContamination \
-I "${work_path}/VCFs/${prefix}.pileups.table" \
-O "${work_path}/VCFs/${prefix}.contamination.table" && echo "** $prefix Calculation Contamination done **"

gatk FilterMutectCalls \
-R "${gatk_bundle}/human_g1k_v37_decoy.fasta" \
-V "${work_path}/VCFs/${prefix}.vcf.gz" \
-O "${work_path}/VCFs/${prefix}.filtered.vcf.gz"  && echo "** $prefix filtering VCFs done **"

#Normalizing and decomposing filtered VCF files
vt normalize "${work_path}/VCFs/${prefix}.filtered.vcf.gz" -r "${gatk_bundle}/human_g1k_v37_decoy.fasta" -o "${work_path}/VCFs/${prefix}.filtered_norm.vcf.gz" && echo "** $prefix Normalizing done  **"

#split multiallelic variants to biallelic
vt decompose -s "${work_path}/VCFs/${prefix}.filtered_norm.vcf.gz" -o "${work_path}/VCFs/${prefix}.filtered_norm_decomp.vcf.gz"  && echo "** $prefix decomposed done  **"
#Annotating variants with Annovar

table_annovar.pl --buildver hg19 "${work_path}/VCFs/${prefix}.filtered_norm_decomp.vcf.gz" /cluster/tools/software/annovar/humandb \
--protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f \
--outfile "${work_path}/Annovar_annot/${prefix}.annot.norm_decomp_filtered" --vcfinput && echo "** $prefix annovar annotation done  **"

fi
