#11/7/2016

#Purpose: to run BigRR on Botrytis SNP and expression data, in order to do a
# genome with association study between Botrytis SNPs and the gene expression

##########
# Set Up

#CHANGE THIS ALWAYS

date <- "12052016"

setwd("~/Dropbox/UCD/rotations/Kliebenstein/botrytis_gwas")
sink(paste(paste("bot_bigRR_500SNPS", date, sep="_"), ".txt", sep=""))
library(bigRR)

data.path <- "~/Box Sync/Botrytis_DKrotation/data/"

##########
# Data

snps.all <- read.csv(paste(data.path, "binSNP_bigRR_MAF20hp/binSNP_bigRR_MAF20hp.csv", sep=""), 
                     header=T, sep=",", row.names=1, check.names = FALSE)
#make the chromosome position a facotr
snps.all$POS <- as.factor(as.character(snps.all$POS))
expression.y <- read.table(paste(data.path, "result.lsm.csv", sep=""), header=T, sep=",", row.names = 1)

#subset the snp.Z to only have SNP information
#snp.Z <- snps.all[,4:96]

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
(setdiff(supp$Isolate.Index, iso.index.keep))

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

#subset the sample genotypes from the isolate index
iso.names.keep <- supp[supp$Isolate.Index %in% iso.index.keep,]$Name
(setdiff(supp$Name, iso.names.keep))

#subset the SNP data
snp.Z <- snps.all[,colnames(snps.all) %in% iso.names.keep]
(dim(snp.Z))
snps.all <- cbind(snps.all[,c(1:3)], snp.Z)

###################
# Test of the Variance
# for each of the genes

#exp.variance <- data.frame(gene = as.character(), gene.var = as.numeric())

#for (j in 1:dim(exp.coi.1)[2]) {
#  gene.var <- cbind(colnames(exp.coi.1)[j], var(exp.coi.1[,j]))
#  exp.variance <- rbind(exp.variance, gene.var)
#}
#colnames(exp.variance) <- c("gene", "exp.var")
#exp.variance$exp.var <- as.numeric(as.character(exp.variance$exp.var))
#write.table(exp.variance, file="expression_variance.txt", sep=",")

###################

#group the chromosomes by number
#bot.chroms <- levels(snps.all$X.CHROM)

#x.chrom.numbers <- data.frame()

#for (i in 1:length(bot.chroms)) {
#  (count <- nrow(snps.all[snps.all$X.CHROM == bot.chroms[i],]))
#  row <- cbind(bot.chroms[i], count)
#  x.chrom.numbers <- rbind(x.chrom.numbers, row)
#}

#chrom.group <- as.factor(c(rep(1, 24531+1301), rep(2, 10562+2956+8732+3332+6303+1792), 
#                           rep(3, 1643+9674+10979), rep(4, 5056+670), rep(5, 20103+16518),
#                           rep(6, 19776+1506+1452+3649), rep(7, 7027+3368+9430), 
#                           rep(8, 2032+9428+7324), rep(9, 998+21661), rep(10, 499+20650),
#                           rep(11, 2995+2365), rep(12, 7693+6126), rep(13, 8032+7087+2233),
#                           rep(14, 6315+11221+1750), rep(15, 765+3280+218+1538+3892),
#                           rep(16, 12309+7082+690+6748+1343+8115+983+2886+425+3470+1894+1078)))

#snps.all$CHROM.GROUP <- chrom.group

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

####################

#Random Sample of genes

#random.traits <- sample(3:length(exp.coi.1), 500, replace=FALSE)
#add one and two
#random.traits <- c(1,2,random.traits)
#order it?
#random.traits <- random.traits[order(random.traits)]
#print(random.traits)

random.traits <- c(1,2,4269,5891,5016,7628,8213,7061,4977,8042,2093,5864,6050,6114,1273,6777,670,2644,4174,4280,1253,1257,6652,2310,4083,1202,2050,
                   8355,8758,3775,4216,356,8521,2034,4201,6080,3426,2218,5960,2181,7831,3316,9196,3416,1387,5000,3941,2132,7925,119,5842,9233,
                   3924,6295,888,6954,3596,5292,3931,309,4305,2249,8715,466,2965,55,530,629,6351,777,4190,2804,3854,4611,6782,490,17,
                   4333,2665,554,4035,4956,6684,7796,3647,6493,1446,7859,5925,4809,4661,3302,4993,2040,8040,960,2696,6974,971,6368,920,8686,
                   1402,3670,758,6027,4402,152,8840,1585,4217,5138,7157,1900,1237,2583,3774,1916,6537,5239,3306,5002,7170,6556,6916,8163,2409,
                   4236,2012,1290,9056,5890,3871,7386,4592,4562,4368,7896,2602,4505,1616,7897,1169,1091,7480,1903,1479,7874,6821,8576,8392,1456,
                   7992,6520,4055,3001,2687,4267,1279,8699,4488,9120,3716,2009,5118,4004,2485,4737,6693,1391,6938,1277,7815,3818,8984,2419,5260,
                   2015,2568,2201,4265,1256,471,4633,5013,4328,3045,3174,5972,2365,3578,5876,760,941,4670,6532,799,9256,7395,2559,306,231,
                   4023,4158,6133,396,583,5870,8328,5868,1703,1373,2775,4534,5115,2421,3332,4750,1359,2112,5962,926,2715,2707,6275,5105,2776,
                   9124,1530,5843,4179,4753,6706,5423,5383,8289,9209,3946,233,6206,8367,6978,4370,3081,1749,1401,3791,8248,5441,3465,582,6384,
                   6884,5738,5387,6985,7981,6500,7551,1959,2661,8656,2810,469,1439,3389,3231,5092,7449,6408,8671,5969,7409,2131,303,4658,7145,
                   5084,7643,9130,6576,7765,1765,7189,7696,4099,6892,4140,6702,713,3156,4336,7193,1995,4887,3736,4245,5034,8357,2314,6956,269,
                   19,5803,3002,3229,614,1404,6470,897,2948,8423,7312,5992,6090,4814,4556,980,4172,7710,3213,8421,5181,5788,1198,2422,6703,
                   7171,3101,5832,2043,8970,3327,9243,8664,8330,352,5826,3248,2929,4844,532,2638,1428,3039,2571,5791,448,5246,742,2232,6768,
                   6374,4790,5093,1113,6719,6084,822,162,1972,8818,3773,8323,7403,8162,2657,6853,1343,2848,6941,6867,6043,3482,8952,4495,8713,
                   6468,8296,4634,733,9219,8748,7007,373,1426,4803,2704,1127,2301,3136,1476,3228,1800,5562,7631,1158,9018,26,6049,4097,7243,
                   2019,4722,5173,6855,1756,3947,6964,7674,4533,7952,6438,973,1416,4127,3686,5408,5617,2125,4110,5518,8413,8395,220,2126,6856,
                   813,4472,998,7318,7667,4071,1427,28,6966,5856,5643,3653,6579,1424,9087,2258,1229,4609,2580,8761,2439,5069,5403,1443,4535,
                   8566,3441,2914,7592,1793,6813,936,5280,2205,4863,8574,5950,7200,7227,6773,7833,4914,2387,5381,3683,6353,7127,5150,8004,4939,
                   7779,7296,4465,3816,8230,5278,2647,1623,288,8286,4732,8063,5430,308,689,3637,398,6220,4434,6013,6506,6812,7019,2876,5409)

exp.500 <- expression.y[,random.traits]
print("Length of Colnames of exp.500 data frame")
print(length(colnames(exp.500)))

#check the order of this....
#expected dimension....

exp.coi1.500 <- exp.500[exp.500$HostGenotype == "coi.1",]
exp.coi1.500 <- exp.coi1.500[exp.coi1.500$Isolate %in% iso.index.keep,]
exp.coi1.500 <- exp.coi1.500[,-c(1:2)]

exp.col0.500 <- exp.500[exp.500$HostGenotype == "col.0",]
exp.col0.500 <- exp.col0.500[exp.col0.500$Isolate %in% iso.index.keep,]
exp.col0.500 <- exp.col0.500[,-c(1:2)]

exp.npr1.500 <- exp.500[exp.500$HostGenotype == "npr.1",]
exp.npr1.500 <- exp.npr1.500[exp.npr1.500$Isolate %in% iso.index.keep,]
exp.npr1.500 <- exp.npr1.500[,-c(1:2)]

####################

#Z-score expression

#Coi
z.score.list <- list()

for (j in 1:length(exp.coi1.500)) {
  gene <- colnames(exp.coi1.500)[j]
  mu <- mean(exp.coi1.500[,j])
  stddev <- sd(x = exp.coi1.500[,j])
  
  z.vector <- as.matrix((exp.coi1.500[,j] - mu)/stddev)
  z.score.list[[j]] <- z.vector
  colnames(z.score.list[[j]]) <- colnames(exp.coi1.500)[j]
}

z.coi1.500 <- as.data.frame(do.call(cbind, z.score.list))

#Col
z.score.list <- list()

for (j in 1:length(exp.col0.500)) {
  gene <- colnames(exp.col0.500)[j]
  mu <- mean(exp.col0.500[,j])
  stddev <- sd(x = exp.col0.500[,j])
  
  z.vector <- as.matrix((exp.col0.500[,j] - mu)/stddev)
  z.score.list[[j]] <- z.vector
  colnames(z.score.list[[j]]) <- colnames(exp.col0.500)[j]
}

z.col0.500 <- as.data.frame(do.call(cbind, z.score.list))

#npr
z.score.list <- list()

for (j in 1:length(exp.npr1.500)) {
  gene <- colnames(exp.npr1.500)[j]
  mu <- mean(exp.npr1.500[,j])
  stddev <- sd(x = exp.npr1.500[,j])
  
  z.vector <- as.matrix((exp.npr1.500[,j] - mu)/stddev)
  z.score.list[[j]] <- z.vector
  colnames(z.score.list[[j]]) <- colnames(exp.npr1.500)[j]
}

z.npr1.500 <- as.data.frame(do.call(cbind, z.score.list))


####################

# Running BigRR

#for each of the phenotype conditions (expression data for a gene), 
# perform bigRR using the SNP data

#make data strucutures to store bigRR output

#make the snp.Z a matrix
snp.Z.mat <- as.matrix(snp.Z)
dim(snp.Z.mat)

###################
#Coi.1

#outdir <- "~/Botrytis"

strt <- 1

outfile.blup <- paste(paste("coi1_BLUP_results", date, sep="_"), ".csv", sep="")
outfile.hem <- paste(paste("coi1_HEM_results", date, sep="_"), ".csv", sep="")

write.table(t(snps.all[,c(1,2)]), file=outfile.blup, sep=",",
            quote=FALSE, row.names = TRUE, col.names = FALSE, append =FALSE)
write.table(t(snps.all[,c(1,2)]), file=outfile.hem, sep=",", 
            quote = FALSE, row.names =TRUE, col.names=FALSE, append=FALSE)

#do gwas analysis in chunks, to not fill the R memory environment
for (p in 1:10) {
  print(paste("This is outer iteration:", p, sep=" "))
  #iterate through the trait data in chunks of 500
  # the reason for this is, the first time we ran it, it maxed out the memory after 705
  # iterations. I want to stay well below that.
  #This is not going to get all of the Botrytis genes. There are 9267. We'll have to add
  # the 267 genes later because I don't know how to to iteratively set an odd number in this chunk fashion.
  stp <- p*50
  
  print(stp)
  #create numberic vector for columns for the chunk
  chunk <- c(strt:stp)
  print("Here's your chunk:")
  print(chunk)
  
  #subset the larger expression data set by that chunk
  exp.chunk <- exp.coi1.500[,chunk]
  print("Chunk dimensions")
  print(dim(exp.chunk))
  
  #create lists to store the output
  chunk.list.BLUP <- list()
  chunk.list.HEM <- list()
  
  #loop through expression subset by the genes and run bigRR
  for (j in 1:dim(exp.chunk)[2]) {
    #j is the expression data for each gene across all tested conditions
    print(j)
    print(colnames(exp.chunk)[j])
    gene <- colnames(exp.chunk)[j]
    print(paste("The trait here is:", gene, sep=" "))
    #print the variance
    print(paste("Expression variance:", var(exp.chunk[,j])), sep=" ")
    
    #make the design matrix
    j.X <- matrix(1, dim(exp.chunk)[1], 1)
    
    #fit the ridge regression
    BLUP.result <- try(bigRR_function(exp.chunk[,j], j.X, t(snp.Z.mat)), TRUE)
    
    if (class(BLUP.result) == "try-error") {
      print(paste(colnames(exp.chunk)[j], "at j =", j, "encountered an error in bigRR() and was skipped"))
      print("ERROR MESSAGE:")
      print(BLUP.result)
      
      #add to chunk list
      chunk.list.BLUP[[gene]] <- as.matrix(rep(9999, nrow(snp.Z.mat)))
      colnames(chunk.list.BLUP[[gene]]) <- paste(gene, "FAIL", sep="_")
      chunk.list.HEM[[gene]] <- as.matrix(rep(9999, nrow(snp.Z.mat)))
      colnames(chunk.list.HEM[[gene]]) <- paste(gene, "FAIL", sep="_")
      
      #go to next iteration
      next
    }
    
    #BLUP.result <- bigRR(y = exp.chunk[,j], X=j.X, Z=t(snps.chrom1.mat), 
    #                     family = gaussian(link="identity"), #poisson for RNA expression data
    #                     tol.err = 1e-06,
    #                     tol.conv = 1e-08,
    #                     GPU = FALSE)
    #updates the obtained bigRR object into a new object with heteroscedastic assumptions
    #meaning that the variability of a variable is unequal acress the range of values of the variable
    # that is being used to predict it. (here, the variability of expression 
    # data for a gene would be unequal across the SNP values)
    HEM.result <- try(bigRR_update(BLUP.result, t(snp.Z.mat)), TRUE)
    if (class(HEM.result) == "try-error") {
      print(paste(colnames(exp.chunk)[j], "at j =", j, "encountered an error in bigRR_update() and was skipped"))
      print("ERROR MESSAGE:")
      print(HEM.result)
      
      #add to chunk list
      #add the BLUP result, which was not an error at this point
      chunk.list.BLUP[[gene]] <- BLUP.result$u
      colnames(chunk.list.BLUP[[gene]]) <- gene
      #add the filler for the gene
      chunk.list.HEM[[gene]] <- as.matrix(rep(9999,nrow(snp.Z.mat)))
      colnames(chunk.list.HEM[[gene]]) <- paste(gene, "FAIL", sep="_")
      
      #go to next iteration
      next
    }
    
    #if it makes it through bigRR() and bigRR_update() without any error, store the data
    chunk.list.BLUP[[gene]] <- BLUP.result$u
    colnames(chunk.list.BLUP[[gene]]) <- gene
    chunk.list.HEM[[gene]] <- HEM.result$u
    colnames(chunk.list.HEM[[gene]]) <- gene
  }
  
  #After the loop, form list into data frame
  BLUP.results.df <- do.call(cbind, chunk.list.BLUP)
  dim(BLUP.results.df)
  
  HEM.results.df <- do.call(cbind, chunk.list.HEM)
  dim(BLUP.results.df)
  
  #and append to a file
  write.table(t(BLUP.results.df), file=outfile.blup, sep=",",
              quote=FALSE, row.names = TRUE, col.names = FALSE, append =TRUE)
  write.table(t(HEM.results.df), file=outfile.hem, sep=",", 
              quote=FALSE, row.names = TRUE, col.names = FALSE, append =TRUE)
  
  #update start
  strt <- stp+1
  print(strt)
}

#checking execution time
proc.time() -timing

###################
#Col0

#outdir <- "~/Botrytis"

strt <- 1

outfile.blup <- paste(paste("col0_BLUP_results", date, sep="_"), ".csv", sep="")
outfile.hem <- paste(paste("col0_HEM_results", date, sep="_"), ".csv", sep="")

write.table(t(snps.all[,c(1,2)]), file=outfile.blup, sep=",",
            quote=FALSE, row.names = TRUE, col.names = FALSE, append =FALSE)
write.table(t(snps.all[,c(1,2)]), file=outfile.hem, sep=",", 
            quote = FALSE, row.names =TRUE, col.names=FALSE, append=FALSE)

#do gwas analysis in chunks, to not fill the R memory environment
for (p in 1:10) {
  print(paste("This is outer iteration:", p, sep=" "))
  #iterate through the trait data in chunks of 500
  # the reason for this is, the first time we ran it, it maxed out the memory after 705
  # iterations. I want to stay well below that.
  #This is not going to get all of the Botrytis genes. There are 9267. We'll have to add
  # the 267 genes later because I don't know how to to iteratively set an odd number in this chunk fashion.
  stp <- p*50
  
  print(stp)
  #create numberic vector for columns for the chunk
  chunk <- c(strt:stp)
  print("Here's your chunk:")
  print(chunk)
  
  #subset the larger expression data set by that chunk
  exp.chunk <- exp.col0.500[,chunk]
  print("Chunk dimensions")
  print(dim(exp.chunk))
  
  #create lists to store the output
  chunk.list.BLUP <- list()
  chunk.list.HEM <- list()
  
  #loop through expression subset by the genes and run bigRR
  for (j in 1:dim(exp.chunk)[2]) {
    #j is the expression data for each gene across all tested conditions
    print(j)
    print(colnames(exp.chunk)[j])
    gene <- colnames(exp.chunk)[j]
    print(paste("The trait here is:", gene, sep=" "))
    #print the variance
    print(paste("Expression variance:", var(exp.chunk[,j])), sep=" ")
    
    #make the design matrix
    j.X <- matrix(1, dim(exp.chunk)[1], 1)
    
    #fit the ridge regression
    BLUP.result <- try(bigRR_function(exp.chunk[,j], j.X, t(snp.Z.mat)), TRUE)
    
    if (class(BLUP.result) == "try-error") {
      print(paste(colnames(exp.chunk)[j], "at j =", j, "encountered an error in bigRR() and was skipped"))
      print("ERROR MESSAGE:")
      print(BLUP.result)
      
      #add to chunk list
      chunk.list.BLUP[[gene]] <- as.matrix(rep(9999, nrow(snp.Z.mat)))
      colnames(chunk.list.BLUP[[gene]]) <- paste(gene, "FAIL", sep="_")
      chunk.list.HEM[[gene]] <- as.matrix(rep(9999, nrow(snp.Z.mat)))
      colnames(chunk.list.HEM[[gene]]) <- paste(gene, "FAIL", sep="_")
      
      #go to next iteration
      next
    }
    
    #BLUP.result <- bigRR(y = exp.chunk[,j], X=j.X, Z=t(snps.chrom1.mat), 
    #                     family = gaussian(link="identity"), #poisson for RNA expression data
    #                     tol.err = 1e-06,
    #                     tol.conv = 1e-08,
    #                     GPU = FALSE)
    #updates the obtained bigRR object into a new object with heteroscedastic assumptions
    #meaning that the variability of a variable is unequal acress the range of values of the variable
    # that is being used to predict it. (here, the variability of expression 
    # data for a gene would be unequal across the SNP values)
    HEM.result <- try(bigRR_update(BLUP.result, t(snp.Z.mat)), TRUE)
    if (class(HEM.result) == "try-error") {
      print(paste(colnames(exp.chunk)[j], "at j =", j, "encountered an error in bigRR_update() and was skipped"))
      print("ERROR MESSAGE:")
      print(HEM.result)
      
      #add to chunk list
      #add the BLUP result, which was not an error at this point
      chunk.list.BLUP[[gene]] <- BLUP.result$u
      colnames(chunk.list.BLUP[[gene]]) <- gene
      #add the filler for the gene
      chunk.list.HEM[[gene]] <- as.matrix(rep(9999,nrow(snp.Z.mat)))
      colnames(chunk.list.HEM[[gene]]) <- paste(gene, "FAIL", sep="_")
      
      #go to next iteration
      next
    }
    
    #if it makes it through bigRR() and bigRR_update() without any error, store the data
    chunk.list.BLUP[[gene]] <- BLUP.result$u
    colnames(chunk.list.BLUP[[gene]]) <- gene
    chunk.list.HEM[[gene]] <- HEM.result$u
    colnames(chunk.list.HEM[[gene]]) <- gene
  }
  
  #After the loop, form list into data frame
  BLUP.results.df <- do.call(cbind, chunk.list.BLUP)
  dim(BLUP.results.df)
  
  HEM.results.df <- do.call(cbind, chunk.list.HEM)
  dim(BLUP.results.df)
  
  #and append to a file
  write.table(t(BLUP.results.df), file=outfile.blup, sep=",",
              quote=FALSE, row.names = TRUE, col.names = FALSE, append =TRUE)
  write.table(t(HEM.results.df), file=outfile.hem, sep=",", 
              quote=FALSE, row.names = TRUE, col.names = FALSE, append =TRUE)
  
  #update start
  strt <- stp+1
  print(strt)
}

#checking execution time
proc.time() -timing

###################
#npr1

#outdir <- "~/Botrytis"

strt <- 1

outfile.blup <- paste(paste("npr1_BLUP_results", date, sep="_"), ".csv", sep="")
outfile.hem <- paste(paste("npr1_HEM_results", date, sep="_"), ".csv", sep="")

write.table(t(snps.all[,c(1,2)]), file=outfile.blup, sep=",",
            quote=FALSE, row.names = TRUE, col.names = FALSE, append =FALSE)
write.table(t(snps.all[,c(1,2)]), file=outfile.hem, sep=",", 
            quote = FALSE, row.names =TRUE, col.names=FALSE, append=FALSE)

#do gwas analysis in chunks, to not fill the R memory environment
for (p in 1:10) {
  print(paste("This is outer iteration:", p, sep=" "))
  #iterate through the trait data in chunks of 500
  # the reason for this is, the first time we ran it, it maxed out the memory after 705
  # iterations. I want to stay well below that.
  #This is not going to get all of the Botrytis genes. There are 9267. We'll have to add
  # the 267 genes later because I don't know how to to iteratively set an odd number in this chunk fashion.
  stp <- p*50
  
  print(stp)
  #create numberic vector for columns for the chunk
  chunk <- c(strt:stp)
  print("Here's your chunk:")
  print(chunk)
  
  #subset the larger expression data set by that chunk
  exp.chunk <- exp.npr1.500[,chunk]
  print("Chunk dimensions")
  print(dim(exp.chunk))
  
  #create lists to store the output
  chunk.list.BLUP <- list()
  chunk.list.HEM <- list()
  
  #loop through expression subset by the genes and run bigRR
  for (j in 1:dim(exp.chunk)[2]) {
    #j is the expression data for each gene across all tested conditions
    print(j)
    print(colnames(exp.chunk)[j])
    gene <- colnames(exp.chunk)[j]
    print(paste("The trait here is:", gene, sep=" "))
    #print the variance
    print(paste("Expression variance:", var(exp.chunk[,j])), sep=" ")
    
    #make the design matrix
    j.X <- matrix(1, dim(exp.chunk)[1], 1)
    
    #fit the ridge regression
    BLUP.result <- try(bigRR_function(exp.chunk[,j], j.X, t(snp.Z.mat)), TRUE)
    
    if (class(BLUP.result) == "try-error") {
      print(paste(colnames(exp.chunk)[j], "at j =", j, "encountered an error in bigRR() and was skipped"))
      print("ERROR MESSAGE:")
      print(BLUP.result)
      
      #add to chunk list
      chunk.list.BLUP[[gene]] <- as.matrix(rep(9999, nrow(snp.Z.mat)))
      colnames(chunk.list.BLUP[[gene]]) <- paste(gene, "FAIL", sep="_")
      chunk.list.HEM[[gene]] <- as.matrix(rep(9999, nrow(snp.Z.mat)))
      colnames(chunk.list.HEM[[gene]]) <- paste(gene, "FAIL", sep="_")
      
      #go to next iteration
      next
    }
    
    #BLUP.result <- bigRR(y = exp.chunk[,j], X=j.X, Z=t(snps.chrom1.mat), 
    #                     family = gaussian(link="identity"), #poisson for RNA expression data
    #                     tol.err = 1e-06,
    #                     tol.conv = 1e-08,
    #                     GPU = FALSE)
    #updates the obtained bigRR object into a new object with heteroscedastic assumptions
    #meaning that the variability of a variable is unequal acress the range of values of the variable
    # that is being used to predict it. (here, the variability of expression 
    # data for a gene would be unequal across the SNP values)
    HEM.result <- try(bigRR_update(BLUP.result, t(snp.Z.mat)), TRUE)
    if (class(HEM.result) == "try-error") {
      print(paste(colnames(exp.chunk)[j], "at j =", j, "encountered an error in bigRR_update() and was skipped"))
      print("ERROR MESSAGE:")
      print(HEM.result)
      
      #add to chunk list
      #add the BLUP result, which was not an error at this point
      chunk.list.BLUP[[gene]] <- BLUP.result$u
      colnames(chunk.list.BLUP[[gene]]) <- gene
      #add the filler for the gene
      chunk.list.HEM[[gene]] <- as.matrix(rep(9999, nrow(snp.Z.mat)))
      colnames(chunk.list.HEM[[gene]]) <- paste(gene, "FAIL", sep="_")
      
      #go to next iteration
      next
    }
    
    #if it makes it through bigRR() and bigRR_update() without any error, store the data
    chunk.list.BLUP[[gene]] <- BLUP.result$u
    colnames(chunk.list.BLUP[[gene]]) <- gene
    chunk.list.HEM[[gene]] <- HEM.result$u
    colnames(chunk.list.HEM[[gene]]) <- gene
  }
  
  #After the loop, form list into data frame
  BLUP.results.df <- do.call(cbind, chunk.list.BLUP)
  dim(BLUP.results.df)
  
  HEM.results.df <- do.call(cbind, chunk.list.HEM)
  dim(BLUP.results.df)
  
  #and append to a file
  write.table(t(BLUP.results.df), file=outfile.blup, sep=",",
              quote=FALSE, row.names = TRUE, col.names = FALSE, append =TRUE)
  write.table(t(HEM.results.df), file=outfile.hem, sep=",", 
              quote=FALSE, row.names = TRUE, col.names = FALSE, append =TRUE)
  
  #update start
  strt <- stp+1
  print(strt)
}

#checking execution time
proc.time() -timing


sink()