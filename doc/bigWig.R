## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  #install devtools if necessary
#  install.packages("devtools")
#  library('devtools')
#  #location of bigWig package and subfolder
#  pkgLoc='andrelmartins/bigWig'
#  subFld='bigWig'
#  devtools::install_github(pkgLoc, subdir=subFld)

## ----eval=FALSE----------------------------------------------------------
#  #install devtools if necessary
#  install.packages("devtools")
#  library('devtools')
#  #Set the working directory to the directory where the source files are located
#  setwd('~/Dir')
#  build()

## ----setup---------------------------------------------------------------
library(bigWig)

## ----eval=FALSE----------------------------------------------------------
#  #directory where data is stored
#  dtDir='/home/directory'
#  # specific bigWig file being used
#  dtFn='GSM3618124_HEK293T_TIR1_Cl4_3hrDMSO_rep1_minus_body_0-mer.bigWig'

## ----echo=FALSE----------------------------------------------------------
#directory where data is stored
dtDir='/home/lutus/_projects/_dev/guertin/bigWig/workflow/_data/GSE126919_RAW/'
# specific bigWig file being used
dtFn='GSM3618124_HEK293T_TIR1_Cl4_3hrDMSO_rep1_minus_body_0-mer.bigWig'

## ----eval=FALSE----------------------------------------------------------
#  load.bigWig(filename, udcDir = NULL)

## ------------------------------------------------------------------------
bw=load.bigWig(paste0(dtDir, dtFn))

## ----eval=FALSE----------------------------------------------------------
#  unload.bigWig(bw)

## ------------------------------------------------------------------------
unload.bigWig(bw)
ls()
remove(bw)
ls()

## ----eval = FALSE--------------------------------------------------------
#  query.bigWig(bw, chrom, start, end, clip = TRUE)

## ----echo=FALSE----------------------------------------------------------
bw=load.bigWig(paste0(dtDir, dtFn))


## ------------------------------------------------------------------------
query.bigWig(bw, chrom='chr1', start=1, end=12000)

## ------------------------------------------------------------------------
bwQ=query.bigWig(bw, chrom='chr1', start=1, end=20000)
bwQ[3]

## ------------------------------------------------------------------------
bwQ[1,]

## ------------------------------------------------------------------------
bwQ[1,2]

## ------------------------------------------------------------------------
bwQ[1,'start']

## ----eval=FALSE----------------------------------------------------------
#  print.bigWig(bw)

## ----echo=FALSE----------------------------------------------------------
cat("bigWig\n")
cat(" version:", bw$version, "\n")
cat(" isCompressed", bw$isCompressed, "\n")
cat(" isSwapped", bw$isSwapped, "\n")
cat(" primaryDataSize:", prettyNum(bw$primaryDataSize, big.mark=','), "\n")
cat(" primaryIndexSize:", prettyNum(bw$primaryIndexSize, big.mark=','), "\n")
cat(" zoomLevels:", bw$zoomLevels, "\n")
cat(" chromCount:", length(bw$chroms), "\n")
for (i in 1:5)
  cat("    ", bw$chroms[i], bw$chromSizes[i], "\n")
cat("    ", "...", "\n")
z=length(bw$chroms)-5
for (i in z:length(bw$chroms))
  cat("    ", bw$chroms[i], bw$chromSizes[i], "\n")

cat(" basesCovered:", prettyNum(bw$basesCovered, big.mark=','), "\n")
cat(" mean: ", bw$mean, "\n")
cat(" min: ", bw$min, "\n")
cat(" max: ", bw$max, "\n")
cat(" std: ", bw$std, "\n")

## ----eval=FALSE----------------------------------------------------------
#  region.bpQuery.bigWig(bw, chrom, start, end,
#                         op = "sum", abs.value = FALSE
#                        bwMap = NULL)
#  region.probeQuery.bigWig(bw, chrom, start, end,
#                        op = "wavg", abs.value = FALSE, gap.value = NA)

## ------------------------------------------------------------------------
query.bigWig(bw, chrom='chr2', start=229990, end=230235)


## ------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='sum')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='sum')

## ------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='max')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='max')

## ------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='min')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='min')

## ------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='avg')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='avg')

## ------------------------------------------------------------------------
230235-229990

## ----echo=FALSE----------------------------------------------------------
#negative bw files
dtDirNeg='/home/lutus/_projects/_dev/guertin/bigWig/workflow/_data/neg/'

# specific bigWig file being used
dtFnPlus='GSM3452725_K562_Nuc_NoRNase_plus.bw'
dtFnMinus='GSM3452725_K562_Nuc_NoRNase_minus.bw'
bw.plus=load.bigWig(paste0(dtDirNeg, dtFnPlus))
bw.minus=load.bigWig(paste0(dtDirNeg, dtFnMinus))

## ------------------------------------------------------------------------
query.bigWig(bw.minus, chrom='chr1', start=10140, end=10190)

## ------------------------------------------------------------------------
region.probeQuery.bigWig(bw.minus,chrom='chr1',start=10140, end=10190, op='avg')
region.probeQuery.bigWig(bw.minus,chrom='chr1',start=10140, end=10190, op='avg', abs.value=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  bed.region.bpQuery.bigWig(bw, bed,
#                            strand = NA, op = "sum", abs.value = FALSE, gap.value = 0,
#                            bwMap = NULL)
#  bed.region.probeQuery.bigWig(bw, bed,
#                            op = "wavg", abs.value = FALSE, gap.value = NA)

## ------------------------------------------------------------------------
bed=data.frame('chr1',10496,10497)
#set column headers
colnames(bed)=c('chrom','start', 'end')

## ------------------------------------------------------------------------
rbind(bed, c('chr2', 10000, 20000))

## ------------------------------------------------------------------------
levels(bed$chrom)=c('chr1', 'chr2')

## ------------------------------------------------------------------------
dim(bed)
attributes(bed)
bed

## ------------------------------------------------------------------------
# note: If you leave out op='', it will default to op='sum'
bed.region.bpQuery.bigWig(bw, bed)

## ------------------------------------------------------------------------
bed=rbind(bed, c('chr2', 10500,10501))

## ------------------------------------------------------------------------
bed.region.bpQuery.bigWig(bw, bed)

## ------------------------------------------------------------------------
bed2=rbind(bed, c('chr1', 13000,14001))
bed.region.bpQuery.bigWig(bw, bed2)

## ----eval=FALSE----------------------------------------------------------
#  step.bpQuery.bigWig(bw, chrom, start, end, step,
#                      strand = NA, op = "sum", abs.value = FALSE, gap.value = 0,
#                      bwMap = NULL, with.attributes = TRUE)
#  
#  step.probeQuery.bigWig(bw, chrom, start, end, step,
#                      op = "wavg", abs.value = FALSE, gap.value = NA,
#                      with.attributes = TRUE)

## ------------------------------------------------------------------------
step.bpQuery.bigWig(bw,chrom='chr1',start=1, end=20001, op='sum', step=1000)

## ------------------------------------------------------------------------
step.bpQuery.bigWig(bw,chrom='chr1',start=1, end=20001, op='sum', step=10000, gap.value=10)

## ----eval=FALSE----------------------------------------------------------
#  bed.step.bpQuery.bigWig(bw, chrom, start, end, step,
#                      strand = NA, op = "sum", abs.value = FALSE, gap.value = 0,
#                      bwMap = NULL, with.attributes = TRUE)
#  
#  bed.step.probeQuery.bigWig(bw, bed, step,
#                      op = "wavg", abs.value = FALSE, gap.value = NA,
#                      with.attributes = TRUE, as.matrix = FALSE)

## ------------------------------------------------------------------------
#Create bed dataframe
bed3 = data.frame('chr1', 15000, 25000)
colnames(bed3)=c('chrom', 'start', 'end')
bed3=rbind(bed3, c("chr1", 30000, 35000))
bed.step.bpQuery.bigWig(bw, bed3, step=1000, op='avg', with.attributes=FALSE)

## ------------------------------------------------------------------------
bed.step.probeQuery.bigWig(bw, bed3, step=1000, op='avg', with.attributes=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  bed6=data.frame('chr1',1,100000,'','','+')
#  colnames=c('chrom', 'start', 'end', 'name', 'score', 'strand')

## ----eval=FALSE----------------------------------------------------------
#  bed6.region.bpQuery.bigWig(bw.plus, bw.minus, bed6,
#                             op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL)
#  
#  bed6.region.probeQuery.bigWig(bw.plus, bw.minus, bed6, step,
#                            op = "wavg", abs.value = FALSE, gap.value = NA,
#                            with.attributes = TRUE, as.matrix = FALSE,
#                            follow.strand = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  dtDir = '/home/directory'
#  dtFnPlus='GSM3452725_K562_Nuc_NoRNase_plus.bw'
#  dtFnMinus='GSM3452725_K562_Nuc_NoRNase_minus.bw'
#  bw.plus=load.bigWig(paste0(dtDirNeg, dtFnPlus))
#  bw.minus=load.bigWig(paste0(dtDirNeg, dtFnMinus))

## ------------------------------------------------------------------------
query.bigWig(bw.minus, chrom='chr1', start=25000, end=50000)
query.bigWig(bw.plus, chrom='chr1', start=25000, end=50000)

## ------------------------------------------------------------------------
bed6=data.frame('chr1',25000,50000,'','','+')
colnames(bed6)=c('chrom', 'start', 'end', 'name', 'score', 'strand')

## ------------------------------------------------------------------------
bed6.region.probeQuery.bigWig(bw.plus, bw.minus,
                  bed6, op='wavg', abs.value = FALSE, gap.value=0)

## ------------------------------------------------------------------------
levels(bed6$strand)=c('+', '-')
bed6=rbind(bed6, c('chr1', 25000, 50000, '', '', '-'))
bed6.region.probeQuery.bigWig(bw.plus, bw.minus, bed6, op='sum', abs.value = FALSE, gap.value=0)

## ------------------------------------------------------------------------
bed6.region.probeQuery.bigWig(bw.plus, bw.minus, bed6,
                  op='sum', abs.value = TRUE, gap.value=0)

## ----eval=FALSE----------------------------------------------------------
#  
#  bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step,
#                           op = "sum", abs.value = FALSE, gap.value = 0,
#                           bwMap = NULL, with.attributes = TRUE, as.matrix = FALSE,
#                           follow.strand = FALSE)
#  
#  bed6.step.probeQuery.bigWig(bw.plus, bw.minus, bed6, step,
#                            op = "wavg", abs.value = FALSE, gap.value = NA,
#                            with.attributes = TRUE, as.matrix = FALSE,
#                            follow.strand = FALSE)

## ------------------------------------------------------------------------
bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step=5000,
                         op = "sum", abs.value = FALSE, gap.value = 0,
                         bwMap = NULL, with.attributes = TRUE, as.matrix = FALSE,
                         follow.strand = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step=5000,
#                           op = "sum", abs.value = FALSE, gap.value = 0,
#                           bwMap = NULL, with.attributes = TRUE, as.matrix = TRUE,
#                           follow.strand = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  #find length of chr1
#  
#  #query region start=len-50000 and end=len-25000

