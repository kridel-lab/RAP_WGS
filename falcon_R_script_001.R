#falcon_R_script_001.R

library(MARATHON)

#load in merged vcf file (including all 21 samples = 3 diagnostic + 17 RAP + 1 control)

vcfFile = "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/multisample.vcf.gz"

coverageDataDemo = readVCFforFalcon(vcfFile)
saveRDS(coverageDataDemo, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/marathon_intermediate_files/coverageData_ALL_RAP.rds")

head(coverageDataDemo)
print("done generating coverage file")

