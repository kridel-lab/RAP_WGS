#-------------------------------------------------------------------------------
#prep_SNVs_for_CCF_calculation.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")

#-------------------------------------------------------------------------------
#Purpose
#-------------------------------------------------------------------------------

#prepare SNV data in format so that it can be run through Hatchet's script
#to extract CCF for every mutation using copy number calls form Hatchet
#ccf=cancer cell fraction
#these CCFs can then be converted into VAF via CCF/2
#these VAF values can then be run through Sciclone to get clusters
#and then through cloneevol to get a tree

#-------------------------------------------------------------------------------
#Analysis
#-------------------------------------------------------------------------------

get_muts = function(patient){
    #subset mutation data for given patient
    #format required
    #CSV file with the following fields
    #(whose names must be specified in the first-row header)
    #chrom, position, Patient, Sample, somatic_status (Somatic),
    #tumor_var_freq, tumor_reads1 (REF count), tumor_reads2 (ALT count)

    #start
    print(patient)

    #get data
    pat_dat = as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
    pat_dat = pat_dat[,c("CHROM", "POS", "STUDY_PATIENT_ID", "Sample", "biotype",
      "gt_AF", "Ref_counts", "alt_counts")]
    colnames(pat_dat)[which(colnames(pat_dat)=="biotype")] = "somatic_status"
    pat_dat$somatic_status = "Somatic"
    #remove "chr"
    pat_dat$CHROM = sapply(pat_dat$CHROM, function(x){unlist(strsplit(x, "chr"))[2]})
    colnames(pat_dat)=c("chrom", "position", "Patient", "Sample", "somatic_status",
    "tumor_var_freq", "tumor_reads1", "tumor_reads2")

    #prep output file
    output_file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Hatchet/", patient,".csv", sep="")

    #write file
    write.csv(pat_dat, output_file, quote=F, row.names=F)
    print("done")
}

patients = c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")
llply(patients, get_muts)

print("done all samples prep")
