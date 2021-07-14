module load samtools

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing

samtools view -b -h LY_0003.bam "7:148495437-148515437" > LY_0003_ctDNA_EZH2_region.bam
samtools index LY_0003_ctDNA_EZH2_region.bam

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files
output=/cluster/projects/kridelgroup/RAP_ANALYSIS/IGV_bam_files

samtools view -b -h LY_RAP_0003_Dia_FoT_01.bam "7:148495437-148515437" > $output/LY_RAP_0003_Dia_FoT_01_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_07.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_07_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_12.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_12_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_05.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_05_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Dia_FoT_05.bam "7:148495437-148515437" > $output/LY_RAP_0003_Dia_FoT_05_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_06.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_06_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_02.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_02_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_16.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_16_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_04.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_04_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_18.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_18_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_15.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_15_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_09.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_09_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_17.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_17_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_03.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_03_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_10.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_10_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_14.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_14_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_13.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_13_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_11.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_11_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Aut_FzT_01.bam "7:148495437-148515437" > $output/LY_RAP_0003_Aut_FzT_01_ctDNA_EZH2_region.bam
samtools view -b -h LY_RAP_0003_Dia_FoT_03.bam "7:148495437-148515437" > $output/LY_RAP_0003_Dia_FoT_03_ctDNA_EZH2_region.bam

cd $output
for i in *.bam

do

echo "Indexing: "$i

samtools index $i

done
