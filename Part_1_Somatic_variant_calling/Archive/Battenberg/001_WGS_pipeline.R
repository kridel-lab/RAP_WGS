library(Battenberg)
library(optparse)

option_list = list(
  make_option(c("-t", "--tumourname"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("-n", "--normalname"), type="character", default=NULL, help="Samplename of the normal", metavar="character"),
  make_option(c("--tb"), type="character", default=NULL, help="Tumour BAM file", metavar="character"),
  make_option(c("--nb"), type="character", default=NULL, help="Normal BAM file", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  make_option(c("--skip_allelecount"), type="logical", default=FALSE, action="store_true", help="Provide when alleles don't have to be counted. This expects allelecount files on disk", metavar="character"),
  make_option(c("--skip_preprocessing"), type="logical", default=FALSE, action="store_true", help="Provide when pre-processing has previously completed. This expects the files on disk", metavar="character"),
  make_option(c("--skip_phasing"), type="logical", default=FALSE, action="store_true", help="Provide when phasing has previously completed. This expects the files on disk", metavar="character"),
  make_option(c("--cpu"), type="numeric", default=8, help="The number of CPU cores to be used by the pipeline (Default: 8)", metavar="character"),
  make_option(c("--bp"), type="character", default=NULL, help="Optional two column file (chromosome and position) specifying prior breakpoints to be used during segmentation", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#TUMOURNAME = opt$tumourname
#NORMALNAME = opt$normalname
TUMOURNAME = "LY_RAP_0003_Aut_FzT_01"
NORMALNAME = "LY_RAP_0003_Ctl_FzG_01"

#NORMALBAM = opt$nb
NORMALBAM = "/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam"

#TUMOURBAM = opt$tb
TUMOURBAM = "/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Aut_FzT_01_files/gatk/LY_RAP_0003_Aut_FzT_01.sorted.dup.recal.cram.bam"

#IS.MALE = opt$sex=="male" | opt$sex=="Male"
IS.MALE = TRUE

#RUN_DIR = opt$output
RUN_DIR = "/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Battenberg"

#SKIP_ALLELECOUNTING = opt$skip_allelecount
#SKIP_PREPROCESSING = opt$skip_preprocessing
SKIP_ALLELECOUNTING = FALSE
SKIP_PREPROCESSING = FALSE

#SKIP_PHASING = opt$skip_phasing
SKIP_PHASING = FALSE

#NTHREADS = opt$cpu
NTHREADS = 8

#PRIOR_BREAKPOINTS_FILE = opt$bp
PRIOR_BREAKPOINTS_FILE = NULL

###############################################################################
# 2018-11-01
# A pure R Battenberg v2.2.9 WGS pipeline implementation.
# sd11 [at] sanger.ac.uk
###############################################################################

# General static
IMPUTEINFOFILE = "/cluster/projects/kridelgroup/RAP_ANALYSIS/Battenberg_files/battenberg_impute_v3/impute_info.txt"
G1000PREFIX = "/cluster/projects/kridelgroup/RAP_ANALYSIS/Battenberg_files/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr"
G1000PREFIX_AC = "/cluster/projects/kridelgroup/RAP_ANALYSIS/Battenberg_files/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr"
GCCORRECTPREFIX = "/cluster/projects/kridelgroup/RAP_ANALYSIS/Battenberg_files/battenberg_wgs_gc_correction_1000g_v3_noNA/1000_genomes_GC_corr_chr_"
REPLICCORRECTPREFIX = "/cluster/projects/kridelgroup/RAP_ANALYSIS/Battenberg_files/battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_"
IMPUTE_EXE = "impute2"

PLATFORM_GAMMA = 1
PHASING_GAMMA = 1
SEGMENTATION_GAMMA = 10
SEGMENTATIIN_KMIN = 3
PHASING_KMIN = 1
CLONALITY_DIST_METRIC = 0
ASCAT_DIST_METRIC = 1
MIN_PLOIDY = 1.6
MAX_PLOIDY = 4.8
MIN_RHO = 0.1
MIN_GOODNESS_OF_FIT = 0.63
BALANCED_THRESHOLD = 0.51
MIN_NORMAL_DEPTH = 10
MIN_BASE_QUAL = 20
MIN_MAP_QUAL = 35
CALC_SEG_BAF_OPTION = 3

# WGS specific static
ALLELECOUNTER = "alleleCounter"
PROBLEMLOCI = "/cluster/projects/kridelgroup/RAP_ANALYSIS/Battenberg_files/probloci_270415.txt.gz"

# Change to work directory and load the chromosome information
setwd(RUN_DIR)

battenberg(tumourname=TUMOURNAME,
           normalname=NORMALNAME,
           tumour_data_file=TUMOURBAM,
           normal_data_file=NORMALBAM,
           ismale=IS.MALE,
           imputeinfofile=IMPUTEINFOFILE,
           g1000prefix=G1000PREFIX,
           g1000allelesprefix=G1000PREFIX_AC,
           gccorrectprefix=GCCORRECTPREFIX,
           repliccorrectprefix=REPLICCORRECTPREFIX,
           problemloci=PROBLEMLOCI,
           data_type="wgs",
           impute_exe=IMPUTE_EXE,
           allelecounter_exe=ALLELECOUNTER,
           nthreads=NTHREADS,
           platform_gamma=PLATFORM_GAMMA,
           phasing_gamma=PHASING_GAMMA,
           segmentation_gamma=SEGMENTATION_GAMMA,
           segmentation_kmin=SEGMENTATIIN_KMIN,
           phasing_kmin=PHASING_KMIN,
           clonality_dist_metric=CLONALITY_DIST_METRIC,
           ascat_dist_metric=ASCAT_DIST_METRIC,
           min_ploidy=MIN_PLOIDY,
           max_ploidy=MAX_PLOIDY,
           min_rho=MIN_RHO,
           min_goodness=MIN_GOODNESS_OF_FIT,
           uninformative_BAF_threshold=BALANCED_THRESHOLD,
           min_normal_depth=MIN_NORMAL_DEPTH,
           min_base_qual=MIN_BASE_QUAL,
           min_map_qual=MIN_MAP_QUAL,
           calc_seg_baf_option=CALC_SEG_BAF_OPTION,
           skip_allele_counting=SKIP_ALLELECOUNTING,
           skip_preprocessing=SKIP_PREPROCESSING,
           skip_phasing=SKIP_PHASING,
           prior_breakpoints_file=PRIOR_BREAKPOINTS_FILE)
