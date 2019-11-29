library(biomaRt)
library(data.table)

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
all_genes <- getBM(attributes=c('ensembl_gene_id',
'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = grch37)

all_genes = as.data.table(all_genes)

saveRDS(all_genes, file="/cluster/home/kisaev/data/ensembl_biomaRt_coordinate_data.rds")