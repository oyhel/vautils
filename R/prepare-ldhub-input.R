#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'

prepare_ldhub_input(df, source){
  if(source == "BOLT" | source == "BOLT-LMM" | source == "BOLTLMM" | source == "bolt"){

  }

}







library(data.table)
test <- fread('/home/oyvind/tmp/z_bmi0-stats-bgen.gz', header = T, stringsAsFactors = F, data.table = F)


snplist <- read.table('/home/oyvind/gdrive/Plasma/Downloads/tempkvitt/w_hm3.noMHC.snplist', header = T, stringsAsFactors = F)


#!/usr/bin/Rscript

##################################################################################
#
#   Script for formatting snptest output to LD hub format
#
##################################################################################

#-------------------------------------Arguments----------------------------------#
args <- commandArgs(TRUE)
input.file <- args[1]
output.file <- args[2]
ldhub.snplist <- args[3]

#-------------------------------------Libraries----------------------------------#
library(data.table)
library(dplyr)

#-------------------------------------Script-------------------------------------#

df <- fread(input.file, header = T, stringsAsFactors = F)
snplist <- fread(ldhub.snplist, header=T, stringsAsFactors = F)

# get the number of samples by adding all genotypes and extracting max==samples==N
df$Nacum <- df$cohort_1_AA + df$cohort_1_AB + df$cohort_1_BB
df$N <- max(df$Nacum)

# get Zscore
df$Zscore <- df$frequentist_add_beta_1 / df$frequentist_add_se_1

# Extract columns
df.out.full <- dplyr::select(df, rsid, alleleB, alleleA, Zscore, N, frequentist_add_pvalue)

# rename columns for LD hub
names(df.out.full) <- c('snpid','A1','A2','Zscore','N','P-value')

# extract only SNPs found in ldhub snplist
df.out <- subset(df.out.full, snpid %in% snplist$SNP)

# Write output file
write.table(x = df.out, file = output.file, col.names = T, row.names = F, quote = F, sep='\t')
