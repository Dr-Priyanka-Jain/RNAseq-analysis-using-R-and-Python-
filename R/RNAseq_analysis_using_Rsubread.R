\\if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")

library(Rsubread)
library(limma)
library(edgeR)
buildindex(basename="chr1_mm10",reference="chr1.fa")

#########Alignment for multiple files
reads1 <- list.files(path=".",pattern = "*.fastq.gz$")
#reads2 <- list.files(path=".",pattern = "_R2_.*.fq$" )
#all.equal(length(reads1),length(reads2))

align(index="chr1_mm10", readfile1=reads1,input_format="FASTQ",output_format="BAM")

bamfiles <-list.files(path=".",pattern = "*.BAM$")

#Feature Count
fc <- featureCounts(files=bamfiles,annot.inbuilt="mm10")
names(fc)
fc$stat
head(fc$annotation)
write.csv(fc$counts, file = "path/to/your/directory")

write.csv(fc$annotation, file = "path/to/your/directory")

write.csv(fc$targets, file = "path/to/your/directory")

sampleInfo <- read.table("sample_info.csv",header=TRUE,sep=",",row.names=1)

#Differential Gene Expression Analysis
dgeFull <-DGEList(counts=fc$counts, gene=fc$annotation[,c("GeneID","Length")],group=sampleInfo$condition)

dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
head(dgeFull$counts)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
head(dgeFull$counts)
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors
normCounts <- cpm(dgeFull)
pseudoNormCounts <- log2(normCounts + 1)
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull
write.csv(dgeFull, file = "C:/RNAseq using_Rsubread/3219673/dgeFull.csv")
dgeTest <- exactTest(dgeFull)
dgeTest
write.csv(dgeTest, file = "C:/RNAseq using_Rsubread/3219673/dgeTest.csv")
hist(dgeTest$table[,"PValue"], breaks=50)
hist(dgeTestFilt$table[,"PValue"], breaks=50)
https://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html



