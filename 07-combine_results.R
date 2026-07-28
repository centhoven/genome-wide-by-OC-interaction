# set working directory and path
setwd("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results")
path<-"C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results"

# chromosomes and SNP chunk size
chroms <- 1:23
chunk_size <- 2000
max_snps <- 40000

# read all RDS files
dfs <- lapply(chroms, function(chr) {
  chr_path <- file.path(path, paste0("Chr", chr))
  
  files <- list.files(
    chr_path,
    pattern = paste0("cox_results_chr", chr, "_"),
    full.names = TRUE
  )
  
  do.call(cbind, lapply(files, function(f) as.data.frame(readRDS(f))))
})

# combine all chromosomes
df <- do.call(cbind, dfs)

# transpose and clean
df_t <- as.data.frame(t(df))
rownames(df_t) <- NULL

# rename columns
colnames(df_t) <- c("snp_name", "beta_snp_m1", "se_snp_m1", "ro_se_snp_m1", "p_snp_m1", "beta_OC_m1", 
                    "se_OC_m1", "ro_se_OC_m1", "p_OC_m1", "beta_int_m1", "se_int_m1", "ro_se_int_m1", "p_int_m1", 
                    "beta_snp_m2", "se_snp_m2", "ro_se_snp_m2", "p_snp_m2", "beta_OC_m2", "se_OC_m2", "ro_se_OC_m2", "p_OC_m2")


# convert numeric columns
df_t[, -1] <- lapply(df_t[, -1], as.numeric)

sum_stats <- df_t

# Save file
#saveRDS(sum_stats, file = file.path(path, "sum_stats.rds"))
# Open file with summary statistics
sum_stats <- readRDS(file = file.path(path, "sum_stats.rds"))

# Change variable names
tmp <- do.call(rbind, strsplit(sum_stats$snp_name, "[:_]"))

sum_stats$CHR_hg38  <- tmp[, 2]
sum_stats$BP_hg38   <- as.integer(tmp[, 3])
sum_stats$A1        <- tmp[, 4] # effect allele
sum_stats$A2        <- tmp[, 5] # reference allele

# Obtain 10 most significant snps GWAS
hits_gwas <- sum_stats[order(sum_stats$p_snp_m2),]
hits_gwas <- hits_gwas[1:10,c("snp_name", "beta_snp_m2", "se_snp_m2", "p_snp_m2", "beta_int_m1", "se_int_m1", "p_int_m1")]

# Save top 10 as csv file
write.table(
  hits_gwas,
  file = "10snps_hits_gwas.csv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Obtain 10 most significant snps GWDIS
hits_gwdis <- sum_stats[order(sum_stats$p_int_m1),]
hits_gwdis <- hits_gwdis[1:10,c("snp_name", "beta_int_m1", "se_int_m1", "p_int_m1", "beta_snp_m2", "se_snp_m2", "p_snp_m2")]

# Save top 10 as csv file
write.table(
  hits_gwdis,
  file = "10snps_hits_gwdis.csv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Change coördinates from hg38 to Hg19

# Install rtracklayer package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

library(rtracklayer)
library(GenomicRanges)

# Download from UCSC
# download.file(
#   "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
#   destfile = "hg38ToHg19.over.chain.gz",
#   mode = "wb"
# )

chain <- import.chain("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results\\hg38ToHg19.over.chain\\hg38ToHg19.over.chain")

# Convert summary stats to GRanges
gr38 <- GRanges(
  seqnames = sum_stats$CHR_hg38,
  ranges   = IRanges(start = sum_stats$BP_hg38, end = sum_stats$BP_hg38)
)

# Liftover
gr19 <- liftOver(gr38, chain)

# Keep only successfully mapped SNPs
mapped <- elementNROWS(gr19) == 1
sumstats_hg19 <- sum_stats[mapped, ]

dim(sum_stats) # 518643 snps
dim(sumstats_hg19) # 516646 snps --> n=1,997 variants without unique mapping
str(sumstats_hg19)

# Change and add variable names
sumstats_hg19$BP <- start(unlist(gr19[mapped]))
sumstats_hg19$CHR <- sub("chr", "", as.character(seqnames(unlist(gr19[mapped]))))
sumstats_hg19[sumstats_hg19$CHR=="Y",] # # 1 snp (DRAGEN:chrX:2134621:C:T_T) on chr Y, so we additionally need to remove this one 
sumstats_hg19 <- sumstats_hg19[!sumstats_hg19$CHR=="Y",]
#sumstats_hg19$CHR <- as.numeric(ifelse(sumstats_hg19$CHR=='X', 23, sumstats_hg19$CHR))

# Save as txt file
write.table(
  sumstats_hg19,
  file = "sumstats_hg19.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

sumstats_hg19 <- read.table("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results\\sumstats_hg19.txt", header = T)
str(sumstats_hg19)

##################################################
# Make two files, one for GWAS and one for GWDIS #
##################################################

############
# For GWAS #
sumstats_gwas <- sumstats_hg19[, c(
  "CHR",
  "BP",
  "A1",
  "A2",
  "beta_snp_m2",
  "se_snp_m2",
  "p_snp_m2")]

names(sumstats_gwas) <- c("CHR", "BP", "A1", "A2", "BETA", "SE", "P")
str(sumstats_gwas)

# Save as txt file
write.table(
  sumstats_gwas,
  file = "sumstats_gwas_FUMA.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)

#############
# For GWDIS #
sumstats_gwdis <- sumstats_hg19[, c(
  "CHR",
  "BP",
  "A1",
  "A2",
  "beta_int_m1",
  "se_int_m1",
  "p_int_m1")]

names(sumstats_gwdis) <- c("CHR", "BP", "A1", "A2", "BETA", "SE", "P")
str(sumstats_gwdis)

# Save as txt file
write.table(
  sumstats_gwdis,
  file = "sumstats_gwdis_FUMA.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)

#####################################################################
# Add rsIDs and N for heritability calculation using ldsc in python #
#####################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Download Plink files (1000G_Phase3_plinkfiles.tgz) from: https://zenodo.org/records/10515792
# https://doi.org/10.5281/zenodo.7768713 (I used version 4)

# Inlezen bim files
bim_files <- list.files("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\1000G_EUR_Phase3_plink", 
                        pattern = "^1000G\\.EUR\\.QC\\.[0-9]+\\.bim$", full.names = TRUE)

bim_list <- lapply(bim_files, function(f) {
  # standaard plink .bim formaat, geen header: CHR SNP CM BP A1 A2
  dt <- fread(f, header = FALSE,
              col.names = c("CHR", "SNP", "CM", "BP", "A1_ref", "A2_ref"))
  cat("  ", basename(f), ": ", nrow(dt), " SNPs\n", sep = "")
  dt[, .(CHR, BP, SNP)]
})

bim_all <- rbindlist(bim_list)
nrow(bim_all) # 9997231

# Inlezen gwas sumstats file
sumstats_gwas <- fread("sumstats_gwas_FUMA.txt") # 516645 SNPs
head(sumstats_gwas)

# Verwijder CHR X uit sumstats
sumstats_gwas <- sumstats_gwas[sumstats_gwas$CHR != "X", ] # 500364 SNPs

# Zorg dat CHR in beide data.tables integer is
sumstats_gwas[, CHR := as.integer(CHR)]
bim_all[, CHR := as.integer(CHR)]

# Match rs ID met BP
setkey(bim_all, CHR, BP)
setkey(sumstats_gwas, CHR, BP)

merged <- merge(sumstats_gwas, bim_all, by = c("CHR", "BP"))
cat("  ", nrow(merged), " van de ", nrow(sumstats_gwas), " SNPs gematcht (",
    round(100 * nrow(merged) / nrow(sumstats_gwas), 1), "%)\n", sep = "")
# 475779 van de 500364 SNPs gematcht (95.1%)

# Add N cases and N controls
merged$N_CAS <- 2067
merged$N_CON <- 200176

# Save file
fwrite(merged, "C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\sumstats_gwas_ldsc.txt", 
       sep = "\t")

# Save SNPlist
writeLines(as.character(merged$SNP), "C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\my_snps.txt")

# Save SNP+allele list
merged <- fread("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\sumstats_gwas_ldsc.txt", 
       sep = "\t")

custom_merge <- merged[, c("SNP", "A1", "A2")]
write.table(custom_merge, "C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\custom_merge_alleles.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#######################################################################
# Add rsIDs and N for heritability calculation using PIGEON in python #
#######################################################################
library(data.table)

# Download Plink files (1000G_Phase3_plinkfiles.tgz) from: https://zenodo.org/records/10515792
# https://doi.org/10.5281/zenodo.7768713 (I used version 4)

# Inlezen bim files
bim_files <- list.files("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\1000G_EUR_Phase3_plink", 
                        pattern = "^1000G\\.EUR\\.QC\\.[0-9]+\\.bim$", full.names = TRUE)

bim_list <- lapply(bim_files, function(f) {
  # standaard plink .bim formaat, geen header: CHR SNP CM BP A1 A2
  dt <- fread(f, header = FALSE,
              col.names = c("CHR", "SNP", "CM", "BP", "A1_ref", "A2_ref"))
  cat("  ", basename(f), ": ", nrow(dt), " SNPs\n", sep = "")
  dt[, .(CHR, BP, SNP)]
})

bim_all <- rbindlist(bim_list)
nrow(bim_all) # 9997231

# Inlezen gwas sumstats file
sumstats_gwdis <- fread("sumstats_gwdis_FUMA.txt") # 516645 SNPs
head(sumstats_gwdis)

# Verwijder CHR X uit sumstats
sumstats_gwdis <- sumstats_gwdis[sumstats_gwdis$CHR != "X", ] # 500364 SNPs

# Zorg dat CHR in beide data.tables integer is
sumstats_gwdis[, CHR := as.integer(CHR)]
bim_all[, CHR := as.integer(CHR)]

# Match rs ID met BP
setkey(bim_all, CHR, BP)
setkey(sumstats_gwdis, CHR, BP)

merged <- merge(sumstats_gwdis, bim_all, by = c("CHR", "BP"))
cat("  ", nrow(merged), " van de ", nrow(sumstats_gwdis), " SNPs gematcht (",
    round(100 * nrow(merged) / nrow(sumstats_gwdis), 1), "%)\n", sep = "")
# 475779 van de 500364 SNPs gematcht (95.1%)

# Add N cases and N controls
merged$N_CAS <- 2067
merged$N_CON <- 200176

# Save file
fwrite(merged, "C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\sumstats_gwdis_ldsc.txt", 
       sep = "\t")

# Save SNP+allele list
merged <- fread("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\sumstats_gwdis_ldsc.txt", 
                sep = "\t")

custom_merge <- merged[, c("SNP", "A1", "A2")]
write.table(custom_merge, "C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\ldsc\\gwdis_merge_alleles.txt", sep = "\t", row.names = FALSE, quote = FALSE)

