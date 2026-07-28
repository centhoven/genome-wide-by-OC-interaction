# Create QQ, Manhattan and calculate lambda

# Read in files
sumstats_gwas <- read.table("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results\\sumstats_gwas_FUMA.txt", header = T)
str(sumstats_gwas)

sumstats_gwdis <- read.table("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results\\sumstats_gwdis_FUMA.txt", header = T)
str(sumstats_gwdis)

#install.packages("qqman")
library(qqman)

# Create QQ plots

png("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results\\GWAScombined_qq.png", width = 8, height = 4, units = "in", res = 300)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))  # 1 row, 2 cols


qq(sumstats_gwdis$P)
mtext("A", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)

qq(sumstats_gwas$P)
mtext("B", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)


dev.off()

# Lambda for gwas
chisq <- qchisq(1-sumstats_gwas$P,1)
lambda_gwas <- median(chisq)/qchisq(0.5,1)

# Lambda for gwdis
chisq <- qchisq(1-sumstats_gwdis$P,1)
lambda_gwdis <- median(chisq)/qchisq(0.5,1)


# CHR must be numeric to create manhattan plot
# GWAS
sumstats_gwas$CHR <- as.numeric(ifelse(sumstats_gwas$CHR == "X", 23, sumstats_gwas$CHR))
sumstats_gwas$SNP <- 0 # to just add this column so that the function will work
str(sumstats_gwas)

manhattan_gwas <- manhattan(sumstats_gwas)

#GWDIS
sumstats_gwdis$CHR <- as.numeric(ifelse(sumstats_gwdis$CHR == "X", 23, sumstats_gwdis$CHR))
sumstats_gwdis$SNP <- 0 # to just add this column so that the function will work
str(sumstats_gwdis)

manhattan_gwdis <- manhattan(sumstats_gwdis)

png("C:\\Users\\clair\\OneDrive - Erasmus MC\\Vancouver\\Project\\Analyses\\GWAS_results\\GWAS_resultscombined_manhattan.png", width = 10, height = 10, units = "in", res = 300)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))  # 2 rows, 1 col

manhattan(sumstats_gwdis)
mtext("A", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)

manhattan(sumstats_gwas)
mtext("B", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)

dev.off()
