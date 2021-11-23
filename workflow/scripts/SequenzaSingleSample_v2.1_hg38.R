
#############################
# script to run             #
# sequenza on snps and cnv  #
#files from Varscan         #
############################

## Modification of the script to take in a Cancer type and modify the ploidy prior probabilities

## Sequenza version in /pughlab/bin/r_lib/library/3.1 has "depth.tumor" bug fixed.
#.libPaths(c("/mnt/work1/users/pughlab/bin/r_lib/library/3.1",.libPaths()))

library(sequenza)
library(optparse)
library(GenomicRanges)

print(sessionInfo())

#############
# FUNCTIONS #
#############


BedBedOverlap <- function(bed.frame1,
                          bed.frame2,
                          pad=0) {

  colnames(bed.frame1)[1:3] <- c("Chromsome","Start_Position","End_Position")
  bed.frame1$Start_Position <- as.numeric(as.character(bed.frame1$Start_Position))
  bed.frame1$End_Position <- as.numeric(as.character(bed.frame1$End_Position))
  bed.frame1$Chromsome <- gsub("chr","",bed.frame1$Chromsome)

  colnames(bed.frame2)[1:3] <- c("Chromsome","Start_Position","End_Position")
  bed.frame2$Start_Position <- as.numeric(as.character(bed.frame2$Start_Position))
  bed.frame2$End_Position <- as.numeric(as.character(bed.frame2$End_Position))
  bed.frame2$Chromsome <- gsub("chr","",bed.frame2$Chromsome)

  ## add padding to bed if desired
  bed.frame1$Start_Position <- bed.frame1$Start_Position - pad
  bed.frame1$End_Position <- bed.frame1$End_Position + pad

  bed.frame2$Start_Position <- bed.frame2$Start_Position - pad
  bed.frame2$End_Position <- bed.frame2$End_Position + pad

  bed_ranges1 <- with(bed.frame1, GRanges(Chromsome, IRanges(start=Start_Position, end=End_Position)))
  bed_ranges2 <- with(bed.frame2, GRanges(Chromsome, IRanges(start=Start_Position, end=End_Position)))


  hits = overlapsAny(bed_ranges1,bed_ranges2)

  return(hits)
}

alternative.cp.solutions <- function(cp.table) {
  ci <- get.ci(cp.table)
  p.alt <- which(diff(sign(diff(ci$values.ploidy$y))) == -2) + 1
  get.alt <- function(idx.p, cp.table) {
    idx.c <- which.max(cp.table$lpp[idx.p,])
    c(cellularity = cp.table$cellularity[idx.c],
      ploidy = cp.table$ploidy[idx.p],
      SLPP = cp.table$lpp[idx.p, idx.c])
  }
  res <- lapply(p.alt, FUN = function (x) get.alt(x, cp.table))
  res <- as.data.frame(do.call(rbind, res))
  if (nrow(res) > 0 ){
    res[order(res$SLPP, decreasing = TRUE), ]
  } else {
    data.frame(cellularity = ci$max.cellularity,
               ploidy = ci$max.ploidy,
               SLPP =  cp.table$lpp[which(cp.table$ploidy == ci$max.ploidy),
                                    which(cp.table$cellularity == ci$max.cellularity)])
  }
}

#################
# PARSE OPTIONS #
#################

## For now, test if commands are in original, trailing format, or new opt-parse format
option_list = list(
  make_option(c("-s", "--snp_file"), type="character", default=NULL,
              help="varscan snp file name", metavar="character"),
  make_option(c("-c", "--cnv_file"), type="character", default=NULL,
              help="varscan copy number file name [default= %default]", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL,
              help="output directory [default= %default]", metavar="character"),
  make_option(c("-p", "--purity"), type="double", default=1,
              help="alternate purity [default= %default]", metavar="double"),
  make_option(c("-P", "--ploidy"), type="double", default=2,
              help="alternate purity [default= %default]", metavar="double"),
  make_option(c("-m", "--min.reads.normal"), type="double", default=10,
              help="min reads to make genotype call [default= %default]", metavar="double"),
  make_option(c("-n", "--remove_centromeres"), type="logical", default=TRUE,
              help="remove all snp and copy-number data from centromeric regions",
              metavar="logical"),
  make_option(c("-t", "--cancer_type"), type="character", default="all",
              help="Cancer type to define priors [default= %default]", metavar="character"),
  make_option(c("-b", "--breaks_method"), type="character", default="het",
              help="breaks method for segmentation [default= %default]", metavar="character"),
  make_option(c("-w", "--window"), type="integer", default=1000000,
              help="window size for plotting [default= %default]", metavar="integer"),
  make_option(c("-r", "--ratio_priority"), type="logical", default=FALSE,
              help="Use only depth ratio for segmentation? [default= %default]", metavar="logical"),
  make_option(c("-y", "--remove_Y_vars"), type="logical", default=TRUE,
              help="Remove Y chromosome variants? [default= %default]", metavar="logical"),
  make_option(c("-f", "--min_reads_baf"), type="integer", default=1,
              help="Threshold on the depth of the positions included to calculate the average BAF for segment. Set to extreme value (ex. 1x10^6) to exclude BAF from calculations [default= %default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


snp.file <- opt$snp_file

cnv.file <- opt$cnv_file

outdir <- opt$out_dir

purity <- opt$purity

usr_cancer_type <- opt$cancer_type

ud_ploidy <- opt$ploidy

min.reads.normal <- opt$min.reads.normal

min.reads.baf <- opt$min_reads_baf
#n_var_freq <- opt$normal_var_freq

remove_Y <- opt$remove_Y_vars

print("Running Sequenza with the following options:")
print(opt)

########
# MAIN #
########

## add another output directory for user defined purity
output_udp <- paste(outdir,"/output_udp_",purity,"/",sep="")
dir.create(output_udp, showWarnings = FALSE, recursive=TRUE)

## read the snp file

snp <- read.table(snp.file,sep="\t",header = TRUE)
cnv <- read.table(cnv.file,sep="\t",header = TRUE)

## remove any non-canonical chromosomes
chr_list <- c(paste0("chr",1:22),"chrX","chrY")

snp <- subset(snp,chrom %in% chr_list)
cnv <- subset(cnv,chrom %in% chr_list)

## remove chromosome Y because it seems to cause problems
if(remove_Y){
  snp<-subset(snp,chrom!="chrY")
  cnv<-subset(cnv,chrom!="chrY")
}

## remove centromeric regions
if (opt$remove_centromeres){
 cent <- read.table("/cluster/home/selghamr/workflows/ExomeSeq/resources/centromeres.txt",
            sep = "\t", header = F)

 del_snp <- BedBedOverlap(snp[c(1,2,2)],
                          cent[c(2:4)])

 snp <- snp[!del_snp,]

 del_cnv <- BedBedOverlap(cnv[c(1:3)],
                          cent[c(2:4)])

 cnv <- cnv[!del_cnv,]

}



## bug where "chr" causes problems with "fast" breaks method

if (opt$breaks_method=="fast") {
  snp$chrom <- gsub("chr","",snp$chrom)
  cnv$chrom <- gsub("chr","",cnv$chrom)
}

##### prepare sequenza data file
seqz.data <- VarScan2seqz(varscan.somatic = snp,
                          varscan.copynumber = cnv)

file <- gsub(".snp", "",snp.file)
write.table(seqz.data, file , col.names = TRUE, row.names = FALSE, sep = "\t")

### full vs. het methods, applicable in different scenarios
data  <- sequenza.extract(file, gz = FALSE,
                          breaks.method = opt$breaks_method,
                          window = opt$window,
                          min.reads.normal = min.reads.normal,
                          min.reads.baf = min.reads.baf)


#### TEST ####
#### TEST ####
# data  <- sequenza.extract(file,gz = FALSE)
#### TEST ####
#### TEST ####



### run sequenza analysis

## feed in prior probabilities based on cancer type

load(file="/cluster/projects/pughlab/src/sequenza_wrapper/PANCAN_ASCAT_ploidy_prob.Rdata")

priors <- subset(ploidy_table,cancer_type==toupper(usr_cancer_type))[c("CN","value")]

CP.example <- sequenza.fit(data,
                           priors.table = priors,
                           ratio.priority = opt$ratio_priority)
#### TEST ####
#### TEST ####
#
# CP.example <- sequenza.fit(data,
#                            priors.table = priors)
#### TEST ####
#### TEST ####

## output series of results files (images and text files)
sequenza.results(sequenza.extract = data, cp.table = CP.example,
                 sample.id = basename(file) , out.dir=outdir)

## output results for all original alternative solutions

### Extract all alternative solutions###

print ("Getting alternative purity solutions");

alt <- alternative.cp.solutions(CP.example);

purities <- alt$cellularity;
ploidies <- alt$ploidy;

## output series of results files (images and text files) for the alternative purity solutions
## index starts at 2 because the first solution is already output as the top hit.

if(length(purities)>1){ ## bug fix to not run alt solutions if there are no alt solutions.
  for (i in 2:length(purities)){
    print (paste("Creating new directories and printing solution:", i, sep = " "));
    output_alt <- paste(outdir,"/sol",i,"_",purities[i],"/",sep="")
    dir.create(output_udp, showWarnings = FALSE, recursive=TRUE)
    sequenza.results(sequenza.extract = data,
                     cp.table = CP.example,
                     sample.id = basename(file) ,
                     out.dir=output_alt,
                     cellularity=purities[i],
                     ploidy=ploidies[i])
  }
}

## output results using forced purity parameter
## If you input purity alone, the ploidy result will still be selected from the previous
## fit, unless you re-fit with the user-defined purity as one of the (using the command below)
## options.  I have opted to instead [force the user to instead provide both the purity and ploidy
## values.

#CP.example <- sequenza.fit(data, cellularity=c(purity,0.9,0.95))

sequenza.results(sequenza.extract = data,
                 cp.table = CP.example,
                 sample.id = basename(file) ,
                 out.dir=output_udp,
                 cellularity=purity,
                 ploidy = ud_ploidy)

## Ouput a total Copy-number seg file for viewing in IGV etc.

data.seg <-data$segments[[1]]
for (i in 2:length(data$segments))
{
  data.seg <- rbind(data.seg,data$segments[[i]])
}

colnames(data.seg)<-c("chrom","loc.start","loc.end","Bf","N.BAF","sd.BAF","seg.mean","num.mark","sd.ratio")
data.seg$seg.mean<-log2(data.seg$seg.mean)
data.seg$ID<-rep(basename(file),nrow(data.seg))
data.seg<-data.seg[,c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")]

setwd(outdir)
write.table(data.seg,paste(basename(file),"_Total_CN.seg",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
