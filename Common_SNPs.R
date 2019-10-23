
### library vcfR (not well due to one chromsome was avaiable)
library(vcfR)
vcf<-read.vcfR("SNP.vcf")
dna <- ape::read.dna("dna/Brapa_sequence_v3.0.fasta",format = "fasta")
gff <- read.table("gff/modified.gff")

chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff) #Creating chromR objects
plot(chrom)
chrom <- masker(chrom) ## filter variation
write.vcf(chrom, file = "filtered.vcf",mask = TRUE) ## error


### statistics gene distribution

library(GenomicFeatures)
library(VariantAnnotation)
tr<-makeTxDbFromGFF("gff/modified.gff")
genes <-genes(tr)
filter_vcf<-readVcf("3SNP_filtered_by_quality.vcf")
snp_gr<-rowRanges(filter_vcf)

fo<-findOverlaps(genes,snp_gr)
unique(queryHits(fo)) #16807 genes with at least one snp, 12336 genes with at least two snps
summary(as.numeric(table(queryHits(fo))))


#### statistics coverge of resequencing with cds
cds_gr<-unlist(cdsBy(tr,"gene"))
library(bamsignals)
bam1_count<-bamCount("resequence_result/Cole_R428-02-R01_good.bam",cds_gr)
bam2_count<-bamCount("resequence_result/Cole_R428-02-R02_good.bam",cds_gr)
cds_gr$reads_count1 <-bam1_count
cds_gr$reads_count2 <-bam2_count
cds_gr$gene_name <-names(cds_gr)
names(cds_gr)<-NULL

library(dplyr)
cds_gr<-as.data.frame(cds_gr)
a<-cds_gr%>%group_by(gene_name)%>%summarise_at(c("reads_count1","reads_count2"),mean)%>%filter(reads_count1>1&reads_count2>1)
#44749/46238
par(mar = c(5,5,4,2))
plot(cex = 1, 
     pch = 16,
     x=log(b$reads_count1+1,2),
     y=log(b$reads_count2+1,2),
     col = ifelse(b$reads_count1>1&b$reads_count2>1,rgb(0, 0, 255, 40, maxColorValue=255),rgb(255, 200, 0, 40, maxColorValue=255)),
     xlab="R01 reads number (log2)",
     ylab="R02 reads number (log2)",
     main="Reads covering cds region of Rapa",
     ylim=c(0,14),
     xlim=c(0,14),
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.5)

legend(0, 14, c("Covered (n = 44749)", "Not covered (n = 1489)"),
       col=c(rgb(0, 0, 255, 200, maxColorValue=255),rgb(255, 200, 0, 200, maxColorValue=255)), 
       pch=16)
