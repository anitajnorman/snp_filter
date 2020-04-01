#read in file from Stacks output file
sumstats <- read.delim('batch_1.sumstats.tsv', skip=1,stringsAsFactors = FALSE, sep='\t')
seqs <- read.delim('batch_1.catalog.tags.tsv', skip=1,stringsAsFactors = FALSE, header=FALSE,sep='\t')

head(sumstats)
seqs[1,]
length(seqs[10])

#limiting snps to those that have one snp per read
countSNPs <- table(sumstats$Locus.ID)
singleSNPs <- which(countSNPs == 1)

runningSNPs1 <- sumstats[singleSNPs,]

#limiting snps to those that fall in the middle of the read
#assuming reads are 130, limiting by 40 nucleotides on each side
midseq <- which(runningSNPs1$Col %in% 40:90)

runningSNPs2 <- runningSNPs1[midseq,]

#limiting snps to equal allele freqs
allelefreq <- which(runningSNPs2$P <= 0.70)

runningSNPs3 <- runningSNPs2[allelefreq,]

#limiting SNPs to min number of inds it was detected in
minN <- which(runningSNPs3$N >= 10)

runningSNPs4 <- runningSNPs3[minN,]

write.table(runningSNPs4,'./snpStats.tsv', sep='\t',quote=FALSE, row.names = FALSE)
