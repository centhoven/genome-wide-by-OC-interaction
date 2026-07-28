#!/usr/bin/env Rscript

install.packages("BEDMatrix")
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(BEDMatrix)
  library(parallel)
})

# ---------------------------
# User inputs (via envvars)
# ---------------------------
bedfile   <- Sys.getenv("BEDFILE")
phenofile <- Sys.getenv("PHENOFILE")
snp_start <- as.integer(Sys.getenv("SNP_START"))
snp_end   <- as.integer(Sys.getenv("SNP_END"))

if (is.na(snp_start) || is.na(snp_end)) stop("Set SNP_START and SNP_END env vars")

# ---------------------------
# Load phenotype + fam
# ---------------------------
ph <- fread(phenofile)
fam <- fread(sub("\\.bed$", ".fam", bedfile), header = FALSE)

# Ensure eid types match
ph[, eid := as.character(eid)]
fam[, V2 := as.character(V2)]

# Create row index that maps phenotype rows to genotype rows
idx_rows <- match(ph$eid, fam$V2)
if (any(is.na(idx_rows))) {
  warning("Some phenotype eids not found in fam; their genotypes will be NA")}

# ---------------------------
# Load genotypes for chunk (columns)
# ---------------------------
geno <- BEDMatrix(bedfile)  
snp_end <- min(snp_end, ncol(geno)) # in case snp_end has a larger number than the number of snps in the chr file
snp_idx <- snp_start:snp_end

# Read the chunk of SNP columns, and only rows corresponding to ph order
# This returns a matrix with rows in same order as fam; we subset to idx_rows
# Note: BEDMatrix supports geno[rows,cols] efficient mapping
geno_chunk <- as.matrix(geno[idx_rows, snp_idx, drop = FALSE])

# Convert column names and prepare
snp_names <- colnames(geno)[snp_idx]
colnames(geno_chunk) <- snp_names

# Precompute numeric event fields (avoid repeated as.numeric calls)
ph[, tstart_num := as.numeric(tstart)]
ph[, tstop_num  := as.numeric(tstop)]
ph[, event_num  := as.numeric(event)]

# Output filename 
chr <- sub(".*/([^_/]+).*", "\\1", bedfile)
outfile <- sprintf("cox_results_%s_%d_%d.RDS", chr, snp_start, snp_end)

# Helper to safely pull coefficient
pull_coef <- function(summ, name, colname) {
  tryCatch({
    # some summaries use rownames like "snp" or "factor(OCstart)1"
    summ$coefficients[name, colname]
  }, error = function(e) NA)
}

beta_snp_m1  <- NA
se_snp_m1    <- NA
ro_se_snp_m1 <- NA
p_snp_m1     <- NA
beta_OC_m1   <- NA
se_OC_m1     <- NA
ro_se_OC_m1  <- NA
p_OC_m1      <- NA
beta_int_m1  <- NA
se_int_m1    <- NA
ro_int_snp_m1<- NA
p_int_m1     <- NA

beta_snp_m2  <- NA
se_snp_m2    <- NA
ro_se_snp_m2 <- NA
p_snp_m2     <- NA
beta_OC_m2   <- NA
se_OC_m2     <- NA
ro_se_OC_m2  <- NA
p_OC_m2      <- NA

# Loop SNPs using mclapply()
cox_func <- function(snp_no){
  snp_name <- snp_names[snp_no]
  ph$snp <- geno_chunk[, snp_no]   # fast assignment; copies only one numeric vector
  
  # Fit models (wrapped in tryCatch so single-snp failure won't stop the loop)
  s1 <- s2 <- NULL
  tryCatch({  
    fit1 <- coxph(Surv(tstart_num, tstop_num, event_num) ~
                    factor(OCstart) * snp +
                    yob + deprivation + menar + sex + acne + ovarydys + endom + dysmen +
                    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + array + cluster(eid),
                  data = ph)
    s1 <- summary(fit1)
    
    beta_snp_m1 <- pull_coef(s1, "snp", "coef")
    se_snp_m1   <- pull_coef(s1, "snp", "se(coef)")
    ro_se_snp_m1<- pull_coef(s1, "snp", "robust se")
    p_snp_m1    <- pull_coef(s1, "snp", "Pr(>|z|)")
    beta_OC_m1  <- pull_coef(s1, "factor(OCstart)1", "coef")
    se_OC_m1    <- pull_coef(s1, "factor(OCstart)1", "se(coef)")
    ro_se_OC_m1 <- pull_coef(s1, "factor(OCstart)1", "robust se")
    p_OC_m1     <- pull_coef(s1, "factor(OCstart)1", "Pr(>|z|)")
    beta_int_m1 <- pull_coef(s1, "factor(OCstart)1:snp", "coef")
    se_int_m1   <- pull_coef(s1, "factor(OCstart)1:snp", "se(coef)")
    ro_se_int_m1<- pull_coef(s1, "factor(OCstart)1:snp", "robust se")
    p_int_m1    <- pull_coef(s1, "factor(OCstart)1:snp", "Pr(>|z|)")
    
  }, error=function(err){
    
    #error message try catch
    message("Unknown error in ", snp_name)
  })
  
  tryCatch({
    fit2 <- coxph(Surv(tstart_num, tstop_num, event_num) ~
                    factor(OCstart) + snp +
                    yob + deprivation + menar + sex + acne + ovarydys + endom + dysmen +
                    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + array + cluster(eid),
                  data = ph)
    s2 <- summary(fit2)
    
    beta_snp_m2 <- pull_coef(s2, "snp", "coef")
    se_snp_m2   <- pull_coef(s2, "snp", "se(coef)")
    ro_se_snp_m2<- pull_coef(s2, "snp", "robust se")
    p_snp_m2    <- pull_coef(s2, "snp", "Pr(>|z|)")
    beta_OC_m2  <- pull_coef(s2, "factor(OCstart)1", "coef")
    se_OC_m2    <- pull_coef(s2, "factor(OCstart)1", "se(coef)")
    ro_se_OC_m2 <- pull_coef(s2, "factor(OCstart)1", "robust se")
    p_OC_m2     <- pull_coef(s2, "factor(OCstart)1", "Pr(>|z|)")
    
  }, error=function(err){
    
    #error message try catch
    message("Unknown error in ", snp_name)
    
  })
  
  if ((snp_no %% 100) == 0) cat(sprintf("Processed %d SNPs (last: %s)\n", snp_no, snp_name))
  
  save <- c(snp_name, beta_snp_m1, se_snp_m1, ro_se_snp_m1, p_snp_m1, beta_OC_m1, se_OC_m1, ro_se_OC_m1,  p_OC_m1, beta_int_m1, se_int_m1, ro_se_int_m1, p_int_m1, 
            beta_snp_m2, se_snp_m2, ro_se_snp_m2, p_snp_m2, beta_OC_m2, se_OC_m2, ro_se_OC_m2, p_OC_m2)
}

cox_results <- mclapply(1:ncol(geno_chunk), FUN=cox_func,  mc.cores=8)

saveRDS(cox_results, file=outfile)

# optional upload
system(sprintf("dx upload %s", cox_results))

cat("Done chunk", snp_start, "-", snp_end, "\n")
