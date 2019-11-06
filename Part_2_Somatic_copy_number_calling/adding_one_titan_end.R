library(data.table)
library(stringr)

files=list.files(pattern="segs.txt.parsed.txt")
for(i in 1:length(files)){
	file=files[i]
	if(length(which(str_detect(file, "plusone"))) ==0){
		file_name = paste(file,"plusone.txt", sep="_")
		f = fread(file)
		f$end = f$end+1
		write.table(f, file=file_name, quote=F, row.names=F, sep="\t")
	}
}