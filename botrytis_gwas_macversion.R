#11/7/2016

#Purpose: to run BigRR on Botrytis SNP and expression data, in order to do a
# genome with association study between Botrytis SNPs and the gene expression

##########
# Set Up

setwd("~/Dropbox/UCD/rotations/Kliebenstein/botrytis_gwas")
library(bigRR)

data.path <- "~/Box Sync/Botrytis_DKrotation/data/"

##########
# Data

snps.all <- read.csv(paste(data.path, "binSNP_bigRR_MAF20hp/binSNP_bigRR_MAF20hp.csv", sep=""), 
                    header=T, sep=",", row.names=1, check.names = FALSE)
#make the chromosome position a facotr
snps.all$POS <- as.factor(as.character(snps.all$POS))
expression.y <- read.table(paste(data.path, "result.lsm.csv", sep=""), header=T, sep=",", row.names = 1)
#make the Isolate a factor
expression.y$Isolate <- as.character(expression.y$Isolate)

####################

#there's a discrepency in the data - there are a few more isolates in the expression data than
# there are in the SNP data set

#library(xlsx)

#supp <- read.xlsx("data/FinalSupplemental-3.xlsx", sheetIndex = 1)
supp <- read.table(paste(data.path, "FinalSupplemental-3.txt", sep=""), sep="\t", header=T)
supp$Name <- as.character(supp$Name)

snps.iso <- colnames(snps.all)[-c(1:3)]

iso.index.keep <- as.character(supp[supp$Name %in% snps.iso,]$Isolate.Index)

#dropped:
setdiff(supp$Isolate.Index, iso.index.keep)

#subset the expression data into the host genotypes
exp.coi.1 <- expression.y[expression.y$HostGenotype == "coi.1",]
exp.coi.1 <- exp.coi.1[exp.coi.1$Isolate %in% iso.index.keep,]
exp.coi.1 <- exp.coi.1[,-c(1:2)]

exp.col.0 <- expression.y[expression.y$HostGenotype == "col.0",]
exp.col.0 <- exp.col.0[exp.col.0$Isolate %in% iso.index.keep,]
exp.col.0 <- exp.col.0[,-c(1:2)]

exp.npr.1 <- expression.y[expression.y$HostGenotype == "npr.1",]
exp.npr.1 <- exp.npr.1[exp.npr.1$Isolate %in% iso.index.keep,]
exp.npr.1 <- exp.npr.1[,-c(1:2)]

#subset the expression.y data frame to only have the expression data
exp.y <- expression.y[,3:9269]

iso.names.keep <- supp[supp$Isolate.Index %in% iso.index.keep,]$Name

#subset the SNP data
snps.Z <- snps.all[,colnames(snps.all) %in% iso.names.keep]
snps.all <- cbind(snps.all[,c(1:3)], snps.Z)

#subset the snp.Z to only have SNP information
#snp.Z <- snps.all[,4:92]

###################
# Test of the Variance
# for each of the genes

exp.variance <- data.frame()

for (j in 1:dim(exp.coi.1)[2]) {
  gene.var <- cbind(colnames(exp.coi.1)[j], var(exp.coi.1[,j]))
  exp.variance <- rbind(exp.variance, gene.var)
}
colnames(exp.variance) <- c("gene", "exp.var")
exp.variance$exp.var <- as.numeric(as.character(exp.variance$exp.var))

####################

# Running BigRR

# the assumption associated with the ridge regression is that markers have little to no effect

#for each of the phenotype conditions (expression data for a gene), 
# perform bigRR using the SNP data

#make the snp.Z a matrix
snp.Z.mat <- as.matrix(snp.Z)
dim(snp.Z.mat)

###################

#group the chromosomes by number
bot.chroms <- levels(snps.all$X.CHROM)

x.chrom.numbers <- data.frame()

for (i in 1:length(bot.chroms)) {
  (count <- nrow(snps.all[snps.all$X.CHROM == bot.chroms[i],]))
  row <- cbind(bot.chroms[i], count)
  x.chrom.numbers <- rbind(x.chrom.numbers, row)
}

chrom.group <- as.factor(c(rep(1, 24531+1301), rep(2, 10562+2956+8732+3332+6303+1792), 
                 rep(3, 1643+9674+10979), rep(4, 5056+670), rep(5, 20103+16518),
                 rep(6, 19776+1506+1452+3649), rep(7, 7027+3368+9430), 
                 rep(8, 2032+9428+7324), rep(9, 998+21661), rep(10, 499+20650),
                 rep(11, 2995+2365), rep(12, 7693+6126), rep(13, 8032+7087+2233),
                 rep(14, 6315+11221+1750), rep(15, 765+3280+218+1538+3892),
                 rep(16, 12309+7082+690+6748+1343+8115+983+2886+425+3470+1894+1078)))

snps.all$CHROM.GROUP <- chrom.group

###################
# FUNCTIONS

bigRR_function <- function(resp.variable, param.design.mat, shrink.param.design) {
  bigRR.result <- bigRR(y = resp.variable, X=param.design.mat, Z=shrink.param.design, 
        family = gaussian(link="identity"),
        tol.err = 1e-06,
        tol.conv = 1e-08,
        GPU = FALSE)
  return(bigRR.result)
}

bigRR_update_function <- function(blup.result, shrink.param.design) {
  hem.result <- bigRR_update(blup.results, shrink.param.design)
  return(hem.result)
}

###################
#Coi.1

snps.chrom1 <- snps.all[snps.all$CHROM.GROUP == 1,]
snps.chrom1 <- snps.chrom1[,-96]

coi.output.BLUP <- snps.chrom1[,c(1,2)]
coi.output.HEM <- snps.chrom1[,c(1,2)]

snps.chrom1.mat <- as.matrix(snps.chrom1[,-c(1:3)])
dim(snps.chrom1.mat)

for (j in 1:dim(exp.coi.1)[2]) {
  #j is the expression data for each gene across all tested conditions
  print(j)
  print(colnames(exp.coi.1)[j])
  gene <- colnames(exp.coi.1)[j]
      
  #make the design matrix
  j.X <- matrix(1, dim(exp.coi.1)[1], 1)
      
  #fit the ridge regression
  BLUP.result <- try(bigRR_function(exp.coi.1[,j], j.X, t(snps.chrom1.mat)), TRUE)
  if (class(BLUP.result) == "try-error") {
    print(paste(colnames(exp.coi.1)[j], "at j =", j, "encountered an error and was skipped"))
    print("ERROR MESSAGE:")
    print(BLUP.result)
    next
  }
  #BLUP.result <- bigRR(y = exp.coi.1[,j], X=j.X, Z=t(snps.chrom1.mat), 
  #                     family = gaussian(link="identity"), #poisson for RNA expression data
  #                     tol.err = 1e-06,
  #                     tol.conv = 1e-08,
  #                     GPU = FALSE)
      #updates the obtained bigRR object into a new object with heteroscedastic assumptions
      #meaning that the variability of a variable is unequal acress the range of values of the variable
      # that is being used to predict it. (here, the variability of expression 
      # data for a gene would be unequal across the SNP values)
  HEM.result <- try(bigRR_update(BLUP.result, t(snps.chrom1.mat)), TRUE)
      
  #capture results
  coi.output.BLUP <- cbind(coi.output.BLUP, BLUP.result$u)
  colnames(coi.output.BLUP) <- c(colnames(coi.output.BLUP[1:(j-1+2)]), gene)
  coi.output.HEM <- cbind(coi.output.HEM, HEM.result$u)
  colnames(coi.output.HEM) <- c(colnames(coi.output.HEM[1:(j-1+2)]), gene)
}
  
# There are instances where convergence for the ridge regression does not occur.

###############
# Ploting
row.names(coi.output.BLUP) <- coi.output.BLUP[,1] 
coi.output.BLUP <- coi.output.BLUP[,-1]

plot(abs(coi.output.BLUP), cex = .6, col = 'slateblue')

