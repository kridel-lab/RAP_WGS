module load vt

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_RESULTS/strelka_filtered

for filename in *.normalized.vcf.gz
   do vt rminfo $filename -t DP -o no_info_DP_$filename
done
