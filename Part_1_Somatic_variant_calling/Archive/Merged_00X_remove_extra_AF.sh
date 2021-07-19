module load vt

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_RESULTS/strelka_filtered

for filename in *.hg19_multianno.vcf
   do vt rminfo $filename -t AF -o no_info_DP_$filename
done


