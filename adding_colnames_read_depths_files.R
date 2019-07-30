getwd()
#/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/tumCounts

library(data.table)

files = list.files(pattern="tumCounts.txt")
for(i in 1:length(files)){
	f=fread(files[i])
	f$V7 = NULL
	colnames(f) = c("chr","posn","ref","refOriginal","nonRef","tumDepth")
	write.table(f, file=files[i], sep="\t", quote=F, row.names=F, col.names=T)
}