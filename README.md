# RNAseq-analysis-using-R-and-Python

**Bulk RNA-seq analysis using python (pyrpipe) and R (Rsubread)**

**RNA SEQUENCING**

RNA sequencing (RNA-Seq) uses the capabilities of high-throughput sequencing methods to provide insight into the transcriptome of a cell. 
RNA-Seq provides far higher coverage and greater resolution of the dynamic nature of the transcriptome. 

**PREREQUISITES FOR R**

 Windows OS (at least 8GB RAM) with working Power shell
 
 Install R(version "4.4")
 
 After installing R run the following command to install Rsubread:
 
 if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
 
 BiocManager::install("Rsubread")

**LOADING R PACKAGES**

Rsubread provides functions for read alignment and feature counting. It is particularly useful for handling large RNA-seq datasets efficiently.

Limma stands for Linear Models for Microarray Data. It is widely used for the analysis of gene expression data.

edgeR stands for Empirical Analysis of Digital Gene Expression Data in R. It is designed for differential expression analysis of RNA-seq count data. It uses statistical methods to model the read counts and test for differential expression.

library(Rsubread)

library(limma)

library(edgeR)

**LISTING FASTQ FILES FOR ANALYSIS**

FASTQ files are a common format used to store raw sequence data from high-throughput sequencing experiments. Each FASTQ file contains sequences of nucleotides along with their corresponding quality scores. 
FastQ files contains raw sequencing reads. Each file represents reads from a specific sample. 
This step is important because it allows us to identify all the FASTQ files in the directory, ensuring that we process all available samples.

reads1 <- list.files(path=".", pattern="*.fastq.gz$")

**BUILDING INDEX FOR REFERENCE GENOME**

Indexing a reference genome is the process of generating a set of data structures that facilitate rapid access to specific locations within the genomic sequence. 
This process involves mapping the sequence of nucleotides in the genome into a format that can be efficiently searched, enabling bioinformatics tools to quickly locate and align sequencing reads to the reference genome. Indexing is crucial for improving the performance and speed of downstream analyses, such as alignment and variant detection.

buildindex(basename="chr1_mm10",reference="chr1.fa")

**ALIGNMENT**

The next step is to align the RNA-seq reads to the reference genome using the “align” function from the Rsubread package.
The input format is indicated as “FASTQ” to show the input files are in FASTQ format.
The output format is specifies as “BAM” to indicate that the output should be in BAM format.
The RNA seq reads are aligned with indexed reference genome.

align(index="chr1_mm10", readfile1=reads1, input_format="FASTQ", output_format="BAM")

**OUTPUT FILE**

A BAM file (*.bam) is the compressed binary version of a SAM file that is used to represent aligned sequences up to 128 Mb.
BAM files store aligned sequence data, which includes information on where each read maps to the reference genome.
BAM files contain a header section and an alignment section:
Header—Contains information about the entire file, such as sample name, sample length, and alignment method. Alignments in the alignments section are associated with specific information in the header section.
Alignments—Contains read name, read sequence, read quality, alignment information, and custom tags. The read name includes the chromosome, start coordinate, alignment quality, and the match descriptor string.

The alignments section includes the following information for each or read pair:

RG: Read group, which indicates the number of reads for a specific sample.

BC: Barcode tag, which indicates the demultiplexed sample ID associated with the read.

SM: Single-end alignment quality.

AS: Paired-end alignment quality.

NM: Edit distance tag, which records the Levenshtein distance between the read and the reference.

XN: Amplicon name tag, which records the amplicon tile ID associated with the read.

bamfiles <-list.files(path=".",pattern = "*.BAM$")

**FEATURE COUNTS**

FeatureCounts is a function under RSubread used in bioinformatics for counting the number of reads (from RNA sequencing) that map to genomic features such as genes, exons, or genomic regions. It is part of the Subread package.
The next step is to count the number of reads that map to each gene using the “featureCounts” function.
This step is crucial because it provides the raw data for differential expression analysis. 
Feature counting quantifies the number of reads that map to each gene, providing the raw data for subsequent analysis of gene expression.

fc <- featureCounts(files=bamfiles,annot.inbuilt="mm10")

write.csv(fc$counts, file = "C:/Users/SUPER/Documents/fc_data.csv")

write.csv(fc$annotation, file = "C:/RNAsequsing_Rsubread/3219673/FC_annotation.csv")

write.csv(fc$targets, file = "C:/RNAseq_using_Rsubread/3219673/FC_targets.csv")

names(fc)

fc$stat

head(fc$annotation)

**LOADING SAMPLE INFORMATION FROM CSV FILE**

Reading sample information provides the metadata needed for differential expression analysis, linking read counts to experimental conditions and ensuring accurate statistical analysis.
A sample file needs to be created with the information given in the image and save the file with .csv extension.
Sample Info is used to describe the experimental condition associated with each sample. Conditions being control and treatment.

Control : Untreated or baseline state
Treatment : Manipulated for experiment

sampleInfo <- read.table("sample_info.csv", header=TRUE, sep=",", row.names=1)



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

boxplot(pseudoNormCounts, col="gray", las=3)
plotMDS(pseudoNormCounts)

dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull
dgeTest <- exactTest(dgeFull)
dgeTest

write.csv(dgeTest, file = "C:/RNAseq using_Rsubread/3219673/dgeTest.csv")

hist(dgeTest$table[,"PValue"], breaks=50)

hist(dgeTestFilt$table[,"PValue"], breaks=50)


