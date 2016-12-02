#11/7/2016

#Purpose: to run BigRR on Botrytis SNP and expression data, in order to do a
# genome with association study between Botrytis SNPs and the gene expression

##########
# Set Up

setwd("~/Dropbox/UCD/rotations/Kliebenstein/botrytis_gwas")
sink("bot_bigRR_allSNPS_11292016.txt")
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

#random.genes <- sample(3:length(exp.coi.1), 500, replace=FALSE)
#add one and two
#random.genes <- c(1,2,random.genes)
#order it?
#random.genes <- random.genes[order(random.genes)]
#print(random.genes)

random.genes <-  c(1,2,7,16,36,48,52,84,96,113,115,130,156,171,236,263,277,279,298,319,321,
                   325,333,354,375,376,389,417,439,446,459,463,470,476,494,501,516,542,596,618,633,636,
                   737,757,764,779,785,801,837,869,885,928,932,936,939,1006,1030,1056,1063,1087,1089,1094,1111,
                   1119,1133,1142,1143,1178,1180,1188,1191,1215,1234,1270,1271,1278,1288,1330,1337,1351,1382,1390,1401,1436,
                   1490,1492,1523,1529,1559,1610,1623,1668,1687,1689,1720,1732,1733,1749,1755,1766,1767,1768,1796,1799,1839,
                   1861,1875,1882,1891,1907,1929,1944,1961,1973,2051,2067,2072,2076,2095,2096,2100,2102,2138,2176,2177,2204,
                   2265,2290,2291,2315,2321,2336,2356,2412,2453,2462,2471,2490,2492,2503,2553,2561,2571,2580,2587,2603,2604,
                   2648,2704,2706,2710,2711,2737,2764,2780,2852,2873,2888,2907,2917,2965,2966,2992,3008,3040,3055,3101,3123,
                   3133,3136,3137,3142,3146,3180,3193,3204,3273,3276,3284,3287,3308,3331,3333,3341,3343,3348,3353,3371,3381,
                   3383,3412,3442,3449,3458,3503,3506,3526,3532,3571,3608,3627,3676,3680,3682,3755,3763,3772,3792,3848,3883,
                   3885,3901,3921,3971,3998,4007,4010,4028,4038,4059,4083,4108,4121,4135,4137,4145,4168,4213,4214,4272,4275,
                   4299,4316,4319,4322,4340,4366,4412,4426,4450,4478,4495,4506,4551,4556,4566,4590,4601,4650,4705,4745,4771,
                   4798,4809,4822,4823,4828,4829,4834,4845,4846,4853,4891,4901,4911,4972,4985,4988,4993,5069,5070,5101,5110,
                   5140,5153,5162,5187,5204,5224,5247,5275,5280,5307,5313,5337,5340,5353,5359,5365,5386,5429,5435,5441,5448,
                   5455,5459,5461,5463,5464,5476,5477,5488,5489,5494,5501,5510,5523,5545,5556,5557,5560,5571,5582,5620,5632,
                   5644,5668,5707,5745,5795,5807,5839,5899,5909,5917,5939,5969,5970,5972,5976,5977,5983,5986,5994,6034,6059,
                   6068,6072,6109,6137,6168,6177,6195,6212,6219,6224,6227,6249,6253,6256,6259,6288,6303,6313,6314,6373,6392,
                   6481,6492,6513,6584,6631,6677,6701,6708,6721,6737,6755,6759,6778,6781,6795,6798,6816,6944,6950,6960,6966,
                   6972,6978,6992,7000,7001,7029,7043,7046,7050,7072,7094,7100,7101,7146,7193,7209,7219,7238,7298,7323,7350,
                   7361,7371,7382,7402,7412,7414,7443,7460,7462,7512,7513,7515,7544,7599,7606,7638,7644,7669,7683,7720,7722,
                   7757,7758,7766,7767,7775,7794,7810,7825,7835,7841,7891,7898,7914,7938,7973,7999,8020,8026,8071,8076,8109,
                   8121,8138,8147,8161,8202,8212,8215,8261,8283,8304,8313,8345,8347,8356,8370,8372,8383,8401,8410,8424,8436,
                   8476,8478,8487,8500,8523,8531,8568,8573,8605,8635,8665,8688,8729,8737,8746,8749,8775,8784,8787,8790,8805,
                   8833,8849,8850,8858,8863,8918,8937,8942,8943,8944,8964,8969,8977,9066,9076,9079,9164,9165,9175)

exp.500 <- expression.y[,random.genes]
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

# Running BigRR

#for each of the phenotype conditions (expression data for a gene), 
# perform bigRR using the SNP data

#make data strucutures to store bigRR output

#make the snp.Z a matrix
snp.Z.mat <- as.matrix(snp.Z)
dim(snp.Z.mat)

###################
#Coi.1
  
#snps.chrom1 <- snps.all[snps.all$CHROM.GROUP == 1,]
#snps.chrom1 <- snps.chrom1[,-96]

#coi.output.BLUP <- snps.all[,c(1,2)]
#coi.output.HEM <- snps.all[,c(1,2)]

#snps.chrom1.mat <- as.matrix(snps.chrom1[,-c(1:3)])
#dim(snps.chrom1.mat)

#outdir <- "~/Botrytis"

strt <- 1

date <- "11292016"

outfile.blup <- paste(paste("coi_BLUP_results", date, sep="_"), ".csv", sep="")
outfile.hem <- paste(paste("coi_HEM_results", date, sep="_"), ".csv", sep="")

write.table(t(snps.all[,c(1,2)]), file=outfile.blup, sep=",",
            quote=FALSE, row.names = TRUE, col.names = FALSE, append =FALSE)
write.table(t(snps.all[,c(1,2)]), file=outfile.hem, sep=",", 
            quote = FALSE, row.names =TRUE, col.names=FALSE, append=FALSE)

#do gwas analysis in chunks, to not fill the R memory environment
#loop through this twise
for (p in 1:10) {
  print(paste("This is outer iteration:", p, sep=" "))
  #iterate through the trait data in chunks of 100
  # the reason for this is, the first time we ran it, it maxed out the memory after 705
  # iterations. I want to stay well below that.
 
  stp <- p*50

  #create numberic vector for columns for the chunk
  chunk <- c(strt:stp)
  print("Here's your chunk:")
  print(chunk)
  
  #subset the larger expression data set by that chunk
  exp.chunk <- exp.coi1.500[,chunk]
  print("Chunk dimensions")
  print(dim(exp.chunk))
  
  print("Names of the genes in this chunk:")
  print(colnames(exp.chunk))
 
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
      chunk.list.HEM[[gene]] <- as.matrix(rep(9999, nrow(coi.output.HEM)))
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
  #print(strt)
}

#list.test.blup <- read.table(outfile.blup, sep=",")

#list.test.hem <- read.table(outfile.hem)

sink()
