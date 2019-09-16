## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(bigWig)

## ------------------------------------------------------------------------
#directory where data is stored
dtDir='/home/directory'
# specific bigWig file being used
dtFn='GSM3618124_HEK293T_TIR1_Cl4_3hrDMSO_rep1_minus_body_0-mer.bigWig'

## ----echo=FALSE----------------------------------------------------------
#directory where data is stored
dtDir='/home/lutus/_projects/_dev/guertin/bigWig/workflow/_data/GSE126919_RAW/'
# specific bigWig file being used
dtFn='GSM3618124_HEK293T_TIR1_Cl4_3hrDMSO_rep1_minus_body_0-mer.bigWig'

## ------------------------------------------------------------------------
bw=load.bigWig(paste0(dtDir, dtFn))

## ------------------------------------------------------------------------
unload.bigWig(bw)

## ----echo=FALSE----------------------------------------------------------
bw=load.bigWig(paste0(dtDir, dtFn))


## ------------------------------------------------------------------------
query.bigWig(bw, chrom='chr1', start=1, end=20000)

## ------------------------------------------------------------------------
bwQ=query.bigWig(bw, chrom='chr1', start=1, end=20000)
bwQ[3]

