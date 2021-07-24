library(sequenza)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Sequenza")

print("hello")

#get arguments
args = commandArgs(trailingOnly = TRUE) #patient ID
index = args[1]
print(index)

data.file <-  index
print(data.file)

#1. sequenza.extract: process seqz data, normalization and segmentation
test <- sequenza.extract(data.file, verbose = TRUE)

#2. sequenza.fit: run grid-search approach to estimate cellularity and ploidy
#CP <- sequenza.fit(test)

#3. sequenza.results: write files and plots using suggested or selected solution
sample_id = paste(unlist(strsplit(unlist(strsplit(data.file, ".small"))[1], "_"))[4:6], collapse="_")
out_dir = unlist(strsplit(data.file, ".small"))[1]

sequenza.results(sequenza.extract = test,
    ploidy = 3.7, cellularity = 0.69, sample.id = sample_id,
    out.dir = out_dir)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#plots and outputs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Grid search maximum likelihood
#pdf(paste(out_dir, "/", "cp_plot.pdf", sep=""))
#cp.plot(CP)
#cp.plot.contours(CP, add = TRUE,
#   likThresh = c(0.999, 0.95),
#   col = c("lightsalmon", "red"), pch = 20)
#dev.off()

#Chromosome view
#Chromosome view is the visualization that displays chromosome by chromosome,
#mutations, B-allele frequency and depth-ratio.
#The visualization makes it easier to inspect the segmentation results,
#comparing to a binned profile of the raw data. It also visualize the copy
#number calling using the cellularity and ploidy solution, making useful to
#asses if the copy number calling is acurate. In addition it provides a
#visualization of the mutation frequency that can also help to corroborate the solution.

#pdf(paste(out_dir, "/", "chrom_view_plot.pdf", sep=""))
#chromosome.view(mut.tab = test$mutations[[1]], baf.windows = test$BAF[[1]],
#                   ratio.windows = test$ratio[[1]],  min.N.ratio = 1,
#                   segments = test$segments[[1]],
#                   main = test$chromosomes[1],
#                   cellularity = 0.89, ploidy = 1.9,
#                   avg.depth.ratio = 1)
#dev.off()
