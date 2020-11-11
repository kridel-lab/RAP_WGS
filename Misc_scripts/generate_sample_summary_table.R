library(data.table)

setwd("/Users/kisaev/UHN/kridel-lab - Documents/RAP_WGS/Data-Files")

#files downloaded from RAP samples master sheet
ff = fread("FF_biopsies_extracted_RAP.csv")
dna = fread("DNA_sequencing_samples.csv")
path = fread("pathology_review_RAP.csv")

#data wrangling
which(dna$Specimen_barcode_FF %in% ff$Specimen_Barcode_Frozen)
dna$Specimen_barcode_FF = as.numeric(dna$Specimen_barcode_FF)
colnames(ff)[which(colnames(ff)=="Specimen_Barcode_Frozen")] = "Specimen_barcode_FF"

merged = merge(ff, dna, by="Specimen_barcode_FF")

#need to add FFPE diagnostic samples
ffpe_samples = c("LY_RAP_0003_Dia_FoT_01" , "LY_RAP_0003_Dia_FoT_03" ,"LY_RAP_0003_Dia_FoT_05")
ffpe = as.data.table(matrix(ncol=ncol(merged), nrow=3))
colnames(ffpe) = colnames(merged)
ffpe = as.data.frame(ffpe)

#ffpe$Indiv = unique(muts$vcf_sample[which(!(muts$vcf_sample %in% dna$Indiv))])
ffpe$LY_RAP_ID = c("LY_RAP_0003_Dia_FoT_05", "LY_RAP_0003_Dia_FoT_01" ,"LY_RAP_0003_Dia_FoT_03")
ffpe$Specimen_Barcode_FFPET =c("15:S12966E", "15:S12966A", "15:S12966C")
ffpe$Tissue_Site = c("left_breast", "right_neck_LN", "left_axilla_LN")
ffpe$Specimen_Type = "FFPE"
ffpe$STUDY_PATIENT_ID = "LY_RAP_0003"
merged = rbind(merged, ffpe)

#which(merged$Specimen_Barcode_FFPET %in% path$Barcode)
#colnames(merged)[which(colnames(merged) == "Specimen_Barcode_FFPET")] = "Barcode"

#merged = merge(merged, path, by = c("Barcode", "STUDY_PATIENT_ID"))
#colnames(merged)[13]="tum_perc"
#colnames(merged)[14]="necrosis_perc"
#colnames(merged)[15]="benign_normal_perc"

merged = unique(merged[,c("STUDY_PATIENT_ID", "Tissue_Site",
"Specimen_Type", "LY_RAP_ID")])

#> unique(all_muts$Indiv)
# [1] "LY_RAP_0001_Aut_FzT_02" "LY_RAP_0001_Aut_FzT_05" "LY_RAP_0001_Aut_FzT_08"
 #[4] "LY_RAP_0002_Aut_FzT_02" "LY_RAP_0002_Aut_FzT_03" "LY_RAP_0002_Aut_FzT_14"
 #[7] "LY_RAP_0002_Aut_FzT_15"

 #"LY_RAP_0003_Aut_FzT_01" "LY_RAP_0003_Aut_FzT_02"
 #"LY_RAP_0003_Aut_FzT_03" "LY_RAP_0003_Aut_FzT_04" "LY_RAP_0003_Aut_FzT_05"
#[13] "LY_RAP_0003_Aut_FzT_06" "LY_RAP_0003_Aut_FzT_07" "LY_RAP_0003_Aut_FzT_09"
#[16] "LY_RAP_0003_Aut_FzT_10" "LY_RAP_0003_Aut_FzT_11" "LY_RAP_0003_Aut_FzT_12"
#[19] "LY_RAP_0003_Aut_FzT_13" "LY_RAP_0003_Aut_FzT_14" "LY_RAP_0003_Aut_FzT_15"
#[22] "LY_RAP_0003_Aut_FzT_16" "LY_RAP_0003_Aut_FzT_17" "LY_RAP_0003_Aut_FzT_18"

#[25] "LY_RAP_0003_Dia_FoT_01" "LY_RAP_0003_Dia_FoT_03" "LY_RAP_0003_Dia_FoT_05"

saveRDS(merged, file="copy_RAP_masterlist_samples.rds")
