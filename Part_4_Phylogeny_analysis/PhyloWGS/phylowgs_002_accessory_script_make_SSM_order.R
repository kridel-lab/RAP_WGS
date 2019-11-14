#need to rename the mutation IDs in SSM file for phylowgs to run 

library(data.table)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS")

# Get all arguments
#sample seed and print it 
seed = sample(1:100000, 1)
print(seed)
set.seed(seed)

args<-commandArgs(TRUE)

args = commandArgs(trailingOnly = TRUE)
index = as.integer(args[1])
print(index)

#read in main SSM file
f = fread("ssm_data.txt")

#randomly sample file rows and keep header
z = sample(f$id, 1000)
random_data = f[which(f$id %in% z),]

random_data$id = paste("s", 0:(nrow(random_data)-1), sep="")
write.table(random_data, file=paste("ssm_data_ordered", "_", index, ".txt", sep=""), sep="\t", quote=F, row.names=F)

print("done")



