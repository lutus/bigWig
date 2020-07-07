## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tinytex.verbose = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  #install devtools if necessary
#  install.packages("devtools")
#  library('devtools')
#  #location of bigWig package and subfolder
#  pkgLoc='andrelmartins/bigWig'
#  subFld='bigWig'
#  devtools::install_github(pkgLoc, subdir=subFld)

## ----eval=FALSE---------------------------------------------------------------
#  #install devtools if necessary
#  install.packages("devtools")
#  library('devtools')
#  #Set the working directory to the directory where the source files are located
#  setwd('~/Dir')
#  build()

## ----setup--------------------------------------------------------------------
library(bigWig)

## ----eval=FALSE---------------------------------------------------------------
#  #directory where data is stored
#  dtDir='/home/directory'
#  # specific bigWig file being used
#  dtFn='GSM3618124_HEK293T_TIR1_Cl4_3hrDMSO_rep1_minus_body_0-mer.bigWig'

## ----echo=FALSE---------------------------------------------------------------
#directory where data is stored
dtDir='/home/lutus/_projects/_dev/guertin/bigWig/workflow/_data/GSE126919_RAW/'
# specific bigWig file being used
dtFn='GSM3618124_HEK293T_TIR1_Cl4_3hrDMSO_rep1_minus_body_0-mer.bigWig'

## ----eval=FALSE---------------------------------------------------------------
#  load.bigWig(filename, udcDir = NULL)

## -----------------------------------------------------------------------------
#load bigWig into variable bw
bw=load.bigWig(paste0(dtDir, dtFn))

## -----------------------------------------------------------------------------
# list all attributes
attributes(bw)

#access individual attribute
bw$basesCovered

## ----eval=FALSE---------------------------------------------------------------
#  unload.bigWig(bw)

## -----------------------------------------------------------------------------
#destroy C object
unload.bigWig(bw)
ls()
#remove variable in R
remove(bw)
ls()

## ----eval = FALSE-------------------------------------------------------------
#  query.bigWig(bw, chrom, start, end, clip = TRUE)

## ----echo=FALSE---------------------------------------------------------------
bw=load.bigWig(paste0(dtDir, dtFn))


## -----------------------------------------------------------------------------
query.bigWig(bw, chrom='chr1', start=1, end=12000)

## -----------------------------------------------------------------------------
bwQ=query.bigWig(bw, chrom='chr1', start=1, end=20000)
bwQ[3]

## -----------------------------------------------------------------------------
bwQ[1,]

## -----------------------------------------------------------------------------
bwQ[1,2]

## -----------------------------------------------------------------------------
bwQ[1,'start']

## ----eval=FALSE---------------------------------------------------------------
#  print.bigWig(bw)

## ----echo=FALSE---------------------------------------------------------------
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

## ----echo=FALSE---------------------------------------------------------------
floc='/home/lutus/_projects/_dev/guertin/bigWig/workflow/_data/testBed.bed'


## -----------------------------------------------------------------------------
bed=read.table(floc, header=FALSE, sep='\t', stringsAsFactors=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  bed=read.table(floc, header=FALSE, sep='\t', stringsAsFactors=FALSE, skip=1)

## -----------------------------------------------------------------------------
chrom=c('chr1')
start=c(300)
end=c(310)
bedT=data.frame(chrom, start, end, stringsAsFactors=FALSE)
bedT

## -----------------------------------------------------------------------------
bedT=rbind(bedT,c('chr1',400,402))
bedT=transform(bedT, start=as.numeric(start), end=as.numeric(end))
bedT

## -----------------------------------------------------------------------------
chrom=c('chr1', 'chr1','chr1', 'chr1')
start=c(300,400,500,600)
end=c(310,410,510,610)
name=c('na','na','na','na')
score=c(1,1,1,1)
strand=c('+','-','+','-')
bed6=data.frame(chrom,start,end,name,score,strand,stringsAsFactors=FALSE)
bed6

## ----eval=FALSE---------------------------------------------------------------
#  center.bed(bed, upstreamWindow, downstreamWindow)
#  
#  fiveprime.bed(bed, upstreamWindow, downstreamWindow)
#  
#  threeprime.bed(bed, upstreamWindow, downstreamWindow)

## -----------------------------------------------------------------------------
bed

## -----------------------------------------------------------------------------
center.bed(bed, upstreamWindow = 0, downstreamWindow = 0)

## -----------------------------------------------------------------------------
center.bed(bed, upstreamWindow = 5, downstreamWindow = 5)

## -----------------------------------------------------------------------------
center.bed(bed, upstreamWindow = -1, downstreamWindow = 4)

## -----------------------------------------------------------------------------
fiveprime.bed(bed, upstreamWindow = 0, downstreamWindow = 0)
threeprime.bed(bed, upstreamWindow = 0, downstreamWindow = 0)

## -----------------------------------------------------------------------------
fiveprime.bed(bed, upstreamWindow = 1, downstreamWindow = 5)
threeprime.bed(bed, upstreamWindow = 1, downstreamWindow = 5)

# negative value
fiveprime.bed(bed, upstreamWindow = -1, downstreamWindow = 5)
threeprime.bed(bed, upstreamWindow = -1, downstreamWindow = 5)

## -----------------------------------------------------------------------------
fiveprime.bed(bed6, upstreamWindow=4, downstreamWindow=2)
threeprime.bed(bed6, upstreamWindow=4, downstreamWindow=2)

## ----eval=FALSE---------------------------------------------------------------
#  downstream.bed(bed, downstreamWindow)
#  
#  upstream.bed(bed, upstreamWindow)
#  

## -----------------------------------------------------------------------------
downstream.bed(bed, downstreamWindow = 5)
upstream.bed(bed, upstreamWindow = 5)

## -----------------------------------------------------------------------------
downstream.bed(bed6,5)
upstream.bed(bed6,5)

## ----eval=FALSE---------------------------------------------------------------
#  foreach.bed(bed, func, envir = parent.frame())

## -----------------------------------------------------------------------------
sizes.bed <- function(bed) {
  N = dim(bed)[1]
  sizes = vector(mode="integer", length=N)

  foreach.bed(bed, function(i, chrom, start, end, strand) {
    sizes[i] <<- end - start
  })

  return(sizes)
}
sizes.bed(bed)

## ----eval=FALSE---------------------------------------------------------------
#  func <- function(i, chrom, start, end, strand) {
#    sizes[i] <<- end - start
#  }
#  
#  sizes.bed <- function(bed) {
#    N = dim(bed)[1]
#    sizes = vector(mode="integer", length=N)
#  
#    foreach.bed(bed, func)
#  
#    return(sizes)
#  }
#  sizes.bed(bed)

## ----eval=FALSE---------------------------------------------------------------
#  region.bpQuery.bigWig(bw, chrom, start, end, strand = NA
#                         op = "sum", abs.value = FALSE,
#                         bwMap = NULL, gap.value = 0)
#  
#  region.probeQuery.bigWig(bw, chrom, start, end,
#                         op = "wavg", abs.value = FALSE,
#                         gap.value = NA)

## -----------------------------------------------------------------------------
query.bigWig(bw, chrom='chr2', start=229990, end=230235)


## -----------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='sum')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='sum')

## -----------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='max')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='max')

## -----------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='min')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='min')

## -----------------------------------------------------------------------------
region.bpQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='avg')
region.probeQuery.bigWig(bw,chrom='chr2',start=229990, end=230235, op='avg')

## -----------------------------------------------------------------------------
230235-229990

## ----echo=FALSE---------------------------------------------------------------
#negative bw files
dtDirNeg='/home/lutus/_projects/_dev/guertin/bigWig/workflow/_data/neg/'

# specific bigWig file being used
dtFnPlus='GSM3452725_K562_Nuc_NoRNase_plus.bw'
dtFnMinus='GSM3452725_K562_Nuc_NoRNase_minus.bw'
bw.plus=load.bigWig(paste0(dtDirNeg, dtFnPlus))
bw.minus=load.bigWig(paste0(dtDirNeg, dtFnMinus))

## -----------------------------------------------------------------------------
query.bigWig(bw.minus, chrom='chr1', start=10140, end=10190)

## -----------------------------------------------------------------------------
region.probeQuery.bigWig(bw.minus,chrom='chr1',start=10140, end=10190, op='avg')
region.probeQuery.bigWig(bw.minus,chrom='chr1',start=10140, end=10190, op='avg', abs.value=TRUE)

## -----------------------------------------------------------------------------
query.bigWig(bw, chrom='chr2', start=229990, end=230235)

## -----------------------------------------------------------------------------
query.bigWig(bw, chrom='chr2', start=229993, end=230001)

## -----------------------------------------------------------------------------
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='min')
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='max')
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='sum')
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='avg')

## -----------------------------------------------------------------------------
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='min', gap.value=1)
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='max', gap.value=1)
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='sum', gap.value=1)
region.bpQuery.bigWig(bw, chrom='chr2', start=229993, end=230001, op='avg', gap.value=1)

## ----eval=FALSE---------------------------------------------------------------
#  bed.region.bpQuery.bigWig(bw, bed, strand = NA,
#                            op = "sum", abs.value = FALSE,
#                            gap.value = 0, bwMap = NULL)
#  bed.region.probeQuery.bigWig(bw, bed, op = "wavg",
#                            abs.value = FALSE, gap.value = NA)

## -----------------------------------------------------------------------------
bed=data.frame('chr1',10496,10497)
#set column headers
colnames(bed)=c('chrom','start', 'end')

## -----------------------------------------------------------------------------
rbind(bed, c('chr2', 10000, 20000))

## -----------------------------------------------------------------------------
levels(bed$chrom)=c('chr1', 'chr2')

## -----------------------------------------------------------------------------
dim(bed)
attributes(bed)
bed

## -----------------------------------------------------------------------------
# note: If you leave out op='', it will default to op='sum'
bed.region.bpQuery.bigWig(bw, bed)

## -----------------------------------------------------------------------------
bed=rbind(bed, c('chr2', 10500,10501))

## -----------------------------------------------------------------------------
bed.region.bpQuery.bigWig(bw, bed)

## -----------------------------------------------------------------------------
bed2=rbind(bed, c('chr1', 13000,14001))
bed.region.bpQuery.bigWig(bw, bed2)

## -----------------------------------------------------------------------------
query.bigWig(bw, chrom='chr2',start=229990, end=229992)

## -----------------------------------------------------------------------------
query.bigWig(bw, chrom='chr2',start=229993, end=230001)

## -----------------------------------------------------------------------------
bedWgap =data.frame('chr2', 229990, 229992)
bedWgap=rbind(bedWgap,c('chr2', 229993, 230001))
colnames(bedWgap)=c('chrom', 'start', 'end')

## -----------------------------------------------------------------------------
bed.region.bpQuery.bigWig(bw, bedWgap, op='avg', gap.value=270)

## -----------------------------------------------------------------------------
bed.region.bpQuery.bigWig(bw, bedWgap, op='avg', gap.value=NA)
bed.region.bpQuery.bigWig(bw, bedWgap, op='avg', gap.value=0)
bed.region.probeQuery.bigWig(bw, bedWgap, op='avg', gap.value=NA)
bed.region.probeQuery.bigWig(bw, bedWgap, op='avg', gap.value=0)

## ----eval=FALSE---------------------------------------------------------------
#  step.bpQuery.bigWig(bw, chrom, start, end, step,
#                      strand = NA, op = "sum", abs.value = FALSE, gap.value = 0,
#                      bwMap = NULL, with.attributes = TRUE)
#  
#  step.probeQuery.bigWig(bw, chrom, start, end, step,
#                      op = "wavg", abs.value = FALSE, gap.value = NA,
#                      with.attributes = TRUE)

## -----------------------------------------------------------------------------
step.bpQuery.bigWig(bw,chrom='chr1',start=1, end=20001, op='sum', step=1000)

## -----------------------------------------------------------------------------
#gap.value=0
step.bpQuery.bigWig(bw,chrom='chr1',start=1, end=20001, op='sum', step=10000,
                    gap.value=0, with.attributes=FALSE)

#gap.value=10
step.bpQuery.bigWig(bw,chrom='chr1',start=1, end=20001, op='sum', step=10000,
                    gap.value=10, with.attributes=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  bed.step.bpQuery.bigWig(bw, chrom, start, end, step,
#                      strand = NA, op = "sum", abs.value = FALSE, gap.value = 0,
#                      bwMap = NULL, with.attributes = TRUE)
#  
#  bed.step.probeQuery.bigWig(bw, bed, step,
#                      op = "wavg", abs.value = FALSE, gap.value = NA,
#                      with.attributes = TRUE, as.matrix = FALSE)

## -----------------------------------------------------------------------------
#Create bed dataframe
bed3 = data.frame('chr1', 15000, 25000)
colnames(bed3)=c('chrom', 'start', 'end')
bed3=rbind(bed3, c("chr1", 30000, 35000))
bed.step.bpQuery.bigWig(bw, bed3, step=1000, op='avg', with.attributes=FALSE)

## -----------------------------------------------------------------------------
bed.step.probeQuery.bigWig(bw, bed3, step=1000, op='avg', with.attributes=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  bed6=data.frame('chr1',1,100000,'','','+')
#  colnames(bed6)=c('chrom', 'start', 'end', 'name', 'score', 'strand')

## ----eval=FALSE---------------------------------------------------------------
#  bed6.region.bpQuery.bigWig(bw.plus, bw.minus, bed6,
#                             op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL)
#  
#  bed6.region.probeQuery.bigWig(bw.plus, bw.minus, bed6, step,
#                            op = "wavg", abs.value = FALSE, gap.value = NA,
#                            with.attributes = TRUE, as.matrix = FALSE,
#                            follow.strand = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  dtDir = '/home/directory'
#  dtFnPlus='GSM3452725_K562_Nuc_NoRNase_plus.bw'
#  dtFnMinus='GSM3452725_K562_Nuc_NoRNase_minus.bw'
#  bw.plus=load.bigWig(paste0(dtDirNeg, dtFnPlus))
#  bw.minus=load.bigWig(paste0(dtDirNeg, dtFnMinus))

## -----------------------------------------------------------------------------
query.bigWig(bw.minus, chrom='chr1', start=25000, end=50000)
query.bigWig(bw.plus, chrom='chr1', start=25000, end=50000)

## -----------------------------------------------------------------------------
bed6=data.frame('chr1',25000,50000,'','','+')
colnames(bed6)=c('chrom', 'start', 'end', 'name', 'score', 'strand')

## -----------------------------------------------------------------------------
bed6.region.probeQuery.bigWig(bw.plus, bw.minus,
                  bed6, op='wavg', abs.value = FALSE, gap.value=0)

## -----------------------------------------------------------------------------
levels(bed6$strand)=c('+', '-')
bed6=rbind(bed6, c('chr1', 25000, 50000, '', '', '-'))
bed6.region.probeQuery.bigWig(bw.plus, bw.minus, bed6, op='sum', abs.value = FALSE, gap.value=0)

## -----------------------------------------------------------------------------
bed6.region.probeQuery.bigWig(bw.plus, bw.minus, bed6,
                  op='sum', abs.value = TRUE, gap.value=0)

## ----eval=FALSE---------------------------------------------------------------
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

## -----------------------------------------------------------------------------
bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step=5000,
                         op = "sum", abs.value = FALSE, gap.value = 0,
                         bwMap = NULL, with.attributes = FALSE, as.matrix = FALSE,
                         follow.strand = FALSE)

## -----------------------------------------------------------------------------
bed6=data.frame('chr1', 1, 100001, 'a', 'c', '+')
colnames(bed6)=c('chrom', 'start', 'end', 'name', 'score', 'strand')
bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step=5000,
                         op = "sum", abs.value = FALSE, gap.value = 0,
                         bwMap = NULL, with.attributes = TRUE, as.matrix = TRUE,
                         follow.strand = FALSE)

## ----eval=TRUE----------------------------------------------------------------
#follow.strand = FALSE
bed6=data.frame('chr1', 1, 100001, 'a', 'c', '-')
colnames(bed6)=c('chrom', 'start', 'end', 'name', 'score', 'strand')
bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step =5000,
                         op = "sum", abs.value = FALSE, gap.value = 0,
                         bwMap = NULL, with.attributes = FALSE, as.matrix = TRUE,
                         follow.strand = FALSE)


## ----eval=TRUE----------------------------------------------------------------
#follow.strand = TRUE
bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step =5000,
                         op = "sum", abs.value = FALSE, gap.value = 0,
                         bwMap = NULL, with.attributes = FALSE, as.matrix = TRUE,
                         follow.strand = TRUE)


## ----eval=FALSE---------------------------------------------------------------
#  quantiles.metaprofile(mat, quantiles = c(0.875, 0.5, 0.125))
#  
#  subsampled.quantiles.metaprofile(mat, quantiles = c(0.875, 0.5, 0.125), fraction = 0.10,
#                                   n.samples = 1000)
#  
#  confinterval.metaprofile(mat, alpha = 0.05)
#  
#  bootstrapped.confinterval.metaprofile(mat, alpha = 0.05, n.samples = 300)
#  
#  metaprofile.bigWig(bed, bw.plus, bw.minus = NULL, step = 1, name = "Signal",
#                     matrix.op = NULL, profile.op = subsampled.quantiles.metaprofile, ...)

## -----------------------------------------------------------------------------
bed6=data.frame('chr1', 1, 100001, 'a', 'c', '+')
colnames(bed6)=c('chrom', 'start', 'end', 'name', 'score', 'strand')
bed6=rbind(bed6, c('chr1', 200001, 300001, 'a', 'c', '+'))
bed6=transform(bed6, start=as.numeric(start), end=as.numeric(end))
bed6
mat=bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step=50000,
                         op = "sum", abs.value = FALSE, gap.value = 0,
                         bwMap = NULL, with.attributes = TRUE, as.matrix = TRUE,
                         follow.strand = FALSE)
mat

## -----------------------------------------------------------------------------
quantiles.metaprofile(mat, quantiles = c(0.875, 0.5, 0.125))

## -----------------------------------------------------------------------------
subsampled.quantiles.metaprofile(mat, quantiles = c(0.875, 0.5, 0.125), fraction = 0.90,
                                 n.samples = 5000)

## -----------------------------------------------------------------------------
confinterval.metaprofile(mat, alpha = 0.05)

## -----------------------------------------------------------------------------
bootstrapped.confinterval.metaprofile(mat, alpha = 0.05, n.samples = 300)

## -----------------------------------------------------------------------------
metaprofile.bigWig(bed6, bw.plus, bw.minus = bw.minus, step = 50000, name = "Signal",
                  matrix.op = NULL, profile.op = quantiles.metaprofile)

## -----------------------------------------------------------------------------
metaprofile.bigWig(bed6, bw.plus, bw.minus = bw.minus, step = 50000, name = "Signal",
                  matrix.op = NULL, profile.op = bootstrapped.confinterval.metaprofile, alpha=0.05, n.samples=1000)

## -----------------------------------------------------------------------------
# Original mat
mat

rpkm.scale(mat, step=50000, libSize=1000000)

## -----------------------------------------------------------------------------
densityToOne.scale(mat, na.on.zero = TRUE)

## ----echo=FALSE---------------------------------------------------------------
mat1=matrix(c(0,4,0,9), nrow=2, ncol=2)

## -----------------------------------------------------------------------------
#Original Matrix
mat1

densityToOne.scale(mat1, na.on.zero = TRUE)
densityToOne.scale(mat1, na.on.zero = FALSE)

## -----------------------------------------------------------------------------
#Original Matrix
mat
maxToOne.scale(mat)
mat1
maxToOne.scale(mat1)

## -----------------------------------------------------------------------------
mat
maxToOne.scale(mat)

## -----------------------------------------------------------------------------
mat1
maxToOne.scale(mat1)

## ----echo=FALSE---------------------------------------------------------------
mat2=matrix(c(2,4,2,9), nrow=2, ncol=2)

## -----------------------------------------------------------------------------
mat2
zeroToOne.scale(mat2)

## -----------------------------------------------------------------------------
metaprofile.bigWig(bed6, bw.plus, bw.minus = bw.minus, step = 50000, name = "Signal",
                  matrix.op = zeroToOne.scale,
                  profile.op = bootstrapped.confinterval.metaprofile)

## ----eval=FALSE---------------------------------------------------------------
#  
#  plot.metaprofile(x, minus.profile = NULL, X0 = x$X0,
#        draw.error = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"),
#        ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = x$name)

## -----------------------------------------------------------------------------
x=metaprofile.bigWig(bed6, bw.plus, bw.minus = bw.minus, step = 50000, name = "Signal",
                  matrix.op = NULL, profile.op = confinterval.metaprofile)

plot.metaprofile(x, minus.profile = NULL, X0 = x$X0,
      draw.error = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"),
      ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = x$name)

## -----------------------------------------------------------------------------
x=metaprofile.bigWig(bed6, bw.plus, bw.minus = bw.minus, step = 50000, name = "Signal",
                  matrix.op = NULL, profile.op = confinterval.metaprofile)
xr=metaprofile.bigWig(bed6, bw.minus, bw.minus = bw.plus, step = 50000, name = "Signal",
                  matrix.op = NULL, profile.op = confinterval.metaprofile)
plot.metaprofile(x, minus.profile = xr, X0 = x$X0,
       draw.error = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"),
       ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = x$name)

## -----------------------------------------------------------------------------
plot.metaprofile(x, minus.profile = xr, X0 = 25000,
      draw.error = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"),
      ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = x$name)

## -----------------------------------------------------------------------------
plot.metaprofile(x, minus.profile = xr, X0 = 25000,
      draw.error = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"),
      ylim = c(-500,500), xlim = c(-1000, 50000), xlab = "Distance (bp)", ylab = x$name)

## -----------------------------------------------------------------------------
plot.metaprofile(x, minus.profile = xr, X0 = 25000,
                draw.error = FALSE, col = c("red", "blue", "lightgrey", "lightgrey"),
                ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = x$name)

## -----------------------------------------------------------------------------
colors()[1:25]

## -----------------------------------------------------------------------------
plot.metaprofile(x, minus.profile = xr, X0 = 25000,
                draw.error = TRUE, col = c("green", "yellow", "purple", "grey"),
                ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = x$name)

