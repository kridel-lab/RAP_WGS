module load vt

for filename in *.vcf
   do vt rminfo $filename -t AF -o no_info_AF_$filename
done


