module load vt

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcfs_annovar_annotated

for filename in *.hg19_multianno.vcf
   do vt rminfo $filename -t AF -o no_info_AF_$filename
done


