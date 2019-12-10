module load bedops 

#go to location where gtf file from telescope manuscript is saved 
#same file that was actually used for telescope analysis
#convert to bed file so that can look at overlap between highly expressed ERVs and 
#ChIP-seq peaks 

cd /cluster/home/kisaev/data

scp HERV_TE_annotations_S1_Table_Matthew_L_Bendall_PLOS_comp_bio_2019.no.chr.gtf telescope_erv_gtf_to_bed.gtf

#remove ## from lines of file
awk '!/##/' telescope_erv_gtf_to_bed.gtf > temp && mv temp telescope_erv_gtf_to_bed.gtf

gtf2bed < telescope_erv_gtf_to_bed.gtf > telescope_erv_gtf_to_bed.bed