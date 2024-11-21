# RNAseq analysis using R

## Bulk RNA-seq analysis using R (Rsubread)

### RNA Sequencing

RNA sequencing (RNA-Seq) uses the capabilities of high-throughput sequencing methods to provide insight into the transcriptome of a cell. 
RNA-Seq provides far higher coverage and greater resolution of the dynamic nature of the transcriptome. 

### Prequisites For R

Windows OS (at least 8GB RAM) with working command line
 
Install R(version "4.4") - https://cran.r-project.org/bin/windows/base/
 
Install RStudio - https://posit.co/download/rstudio-desktop/
 
After installing R run the following command to install Rsubread:
``` 
if (!require("BiocManager", quietly = TRUE))
 
install.packages("BiocManager")
 
BiocManager::install("Rsubread")

BiocManager::install("edgeR")

BiocManager::install("limma")
```
The input files for the RNAseq analysis are to be  downlaoded from link : https://figshare.com/s/f5d63d8c265a05618137 OR from the R folder

- SRR1552444.fastq.gz

- SRR1552445.fastq.gz

- SRR1552454.fastq.gz

- SRR1552455.fastq.gz

**Download the reference genome from this link : 

https://figshare.com/s/f5d63d8c265a05618137 
The following files are te refrence files named as chr1_mm10 and the index file named: chr1_mm10.files, chr1_mm10.00.b.tab and chr1_mm10.00.b.array.

Download all the files given under the R folder.

Set the working directory in R to where you have downloaded all the files:

`setwd("C:/path/to/your/directory")`

### Step 1 : Loading R Packages

Rsubread provides functions for read alignment and feature counting. It is particularly useful for handling large RNA-seq datasets efficiently.

Limma is an R/Bioconductor software package that provides an integrated solution for analysing data from gene expression experiments. It contains rich features for handling complex experimental designs and for information borrowing to overcome the problem of small sample sizes.

edgeR is a Bioconductor package for differential expression analysis of digital gene expression data.
```
library(Rsubread)

library(limma)

library(edgeR)
```

### Step 2 : QC of the raw reads
FastQc to check the quality of raw reads

Download from the link below: 

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


### Step 3 : Listing FASTQ Files for Analysis

FASTQ files are a common format used to store raw sequence data obtained from high-throughput sequencing experiments. Each FASTQ file contains sequences of nucleotides along with their corresponding quality scores. 
FastQ files contains raw sequencing reads. Each file represents reads from a specific sample. 
The format is widely used because it efficiently represents sequencing data and is compatible with many bioinformatics tools for downstream analysis.
This step is important because it allows us to identify all the FASTQ files in the directory, ensuring that we process all available samples.
```
reads1 <- list.files(path=".", pattern="*.fastq.gz$")
```
### Step 4 : Building Index For Reference Genome

Indexing a reference genome is the process of generating a set of data structures that facilitate rapid access to specific locations within the genomic sequence. 
Indexing is crucial for improving the performance and speed of downstream analyses, such as alignment and variant detection.

It is akin to creating a table of contents for a book, where each chapter or section corresponds to a specific region in the genome. This process organizes the sequence data into a structure that allows tools to quickly locate and retrieve specific segments, much like flipping to a chapter by its page number.
```
buildindex(basename="chr1_mm10",reference="chr1.fa")
```
### Step 5 : Alignment

The next step is to align the RNA-seq reads which is the FASTQ files to the reference genome using the “align” function from the Rsubread package.
The input format is indicated as “FASTQ” to show the input files are in FASTQ format.
This process involves mapping the sequence of nucleotides in the genome into a format that can be efficiently searched, enabling bioinformatics tools to quickly locate and align sequencing reads to the reference genome.

The output format is specifies as “BAM” to indicate that the output should be in BAM format.
The RNA seq reads are aligned with indexed reference genome.
```
align(index="chr1_mm10", readfile1=reads1, input_format="FASTQ", output_format="BAM")
```
### Step 6 : BAM File

A BAM (Binary Alignment Map) file (*.bam) is the compressed binary version of a SAM(Sequence Alignment Map) file that is used to represent aligned sequences.
BAM files store aligned sequence data, which includes information on where each read maps to the reference genome.
BAM files contain a header section and an alignment section:
Header—Contains information about the entire file, such as sample name, sample length, and alignment method. Alignments in the alignments section are associated with specific information in the header section.
Alignments—Contains read name, read sequence, read quality, alignment information, and custom tags. The read name includes the chromosome, start coordinate, alignment quality, and the match descriptor string.

The alignments section includes the following information for each or read pair:

- RG: Read group, which indicates the number of reads for a specific sample.

- BC: Barcode tag, which indicates the demultiplexed sample ID associated with the read.

- SM: Single-end alignment quality.

- AS: Paired-end alignment quality.

- NM: Edit distance tag, which records the Levenshtein distance between the read and the reference.

- XN: Amplicon name tag, which records the amplicon tile ID associated with the read.
```
bamfiles <-list.files(path=".",pattern = "*.BAM$")
```
### Step 7 : Feature Counts

FeatureCounts is a function under RSubread used in bioinformatics for counting the number of reads (from RNA sequencing) that map to genomic features such as genes, exons, or genomic regions. It is part of the Subread package.

This step is to count the number of reads that map to each gene using the “featureCounts” function.

This step is crucial because it provides the raw data for differential expression analysis. 
```
fc <- featureCounts(files=bamfiles,annot.inbuilt="mm10")

names(fc)

fc$stat

head(fc$annotation)

write.csv(fc$counts, file = "path/to/your/directory")

write.csv(fc$annotation, file = "path/to/your/directory")

write.csv(fc$targets, file = "path/to/your/directory")
```

### Step 8 : Loading Sample Information from the CSV File

Reading sample information provides the metadata needed for differential expression analysis, linking read counts to experimental conditions and ensuring accurate statistical analysis.
A sample file needs to be created with the information given in the image and save the file with .csv extension.
Sample Info is used to describe the experimental condition associated with each sample. Conditions being control and treatment.

Control : Untreated or baseline state

Treatment : Manipulated for experiment

File	                               Condition
SRR1552444.fastq.gz.subread.BAM	       V

SRR1552445.fastq.gz.subread.BAM	       V

SRR1552454.fastq.gz.subread.BAM"	      L

SRR1552455.fastq.gz.subread.BAM	       L
```
sampleInfo <- read.table("sample_info.csv", header=TRUE, sep=",", row.names=1)
```
### STEP 9:  Differential Gene Expression
Differential expression is the process of determining the differences in gene expression levels between different biological conditions. It identifies which set of genes are expressed at different levels under varying experimental conditions, such as treatments, time points, or disease states.

edgeR stores data in a simple list-based data object called a DGEList. This type of object is
easy to use because it can be manipulated like any list in R. 

dgeFull is a variable here in which we are saving DGEList.
```
dgeFull <-DGEList(counts=fc$counts, gene=fc$annotation[,c("GeneID","Length")],group=sampleInfo$condition)
```
The code below filters out genes that have zero counts across all samples.
```
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],group=dgeFull$samples$group)

head(dgeFull$counts)
```
The normLibSizes function normalizes the library sizes in such a way to minimize the log-fold
changes between the samples for most genes. The default method for computing these scale
factors uses a trimmed mean of M-values (TMM) between each pair of samples.

TMM stands for Trimmed Mean of M-values. It is a normalization method that adjusts for compositional differences between samples. TMM aims to make the majority of genes have similar expression levels across samples by trimming extreme values and calculating a scaling factor for each sample.

We call the product of the original library size and the scaling factor the **effective library size**,
i.e., the normalized library size. The effective library size replaces the original library size in
all downstream analyses.
```
dgeFull <- calcNormFactors(dgeFull, method="TMM")

dgeFull$samples

head(dgeFull$counts)

eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors

normCounts <- cpm(dgeFull)
```
The pseudo-counts represent the equivalent counts would have been
observed had the library sizes all been equal, assuming the fitted model. The pseudo-counts
are computed for a specific purpose, and their computation depends on the experimental
design as well as the library sizes. 
Log transformation reduces the influence of highly expressed genes, making patterns in lower-expressed genes easier to identify.
```
pseudoNormCounts <- log2(normCounts + 1)
```
The function plotMDS draws a multi-dimensional scaling plot of the RNA samples in which
distances correspond to leading log-fold-changes between each pair of RNA samples. The
leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold changes between each pair of samples.

It generates an MDS plot, where the distance between points (samples) represents the similarity of their gene expression profiles

```
plotMDS(pseudoNormCounts)
```
**estimateCommonDisp Estimate Common Negative Binomial Dispersion by Conditional Maximum Likelihood**

The estimateCommonDisp function in the edgeR package is used to estimate a common dispersion parameter for a set of counts following a negative binomial distribution. This helps in accurately modeling the data and improving the reliability of downstream analyses, such as identifying differentially expressed genes.

It estimates a common dispersion parameter for all genes in the dataset. Dispersion reflects the variability of gene expression across replicates and is critical for modeling RNA-Seq data using a negative binomial distribution.
```
dgeFull <- estimateCommonDisp(dgeFull)
```
**EstimateTagwiseDisp Estimate Empirical Bayes Tagwise Dispersion Values**

The estimateTagwiseDisp function refines the RNA-seq data analysis by providing gene-specific dispersion estimates using the empirical Bayes method. This process enhances the accuracy of differential expression analysis by accounting for the unique variability of each gene.

It computes tagwise dispersion values for each gene using an empirical Bayes approach. Unlike the common dispersion, tagwise dispersion accounts for gene-specific variability.

```
dgeFull <- estimateTagwiseDisp(dgeFull)

dgeFull
```
The exact test is a statistical method used in RNA-seq data analysis to identify differentially expressed genes between experimental groups. This test compares the read counts for each gene between groups, taking into account the estimated dispersion.By performing the exact test, one can determine which genes show statistically significant differences in expression between conditions, providing insights into the underlying biological processes and responses.

This test compares read counts for each gene while accounting for the estimated dispersion.

```
dgeTest <- exactTest(dgeFull)

dgeTest

write.csv(dgeTest, file = "path/to/your/directory")
```
The P-value histogram helps evaluate the analysis quality. A peak near zero indicates a strong signal (many genes are differentially expressed), 
while a uniform distribution suggests no significant differences.

```
hist(dgeTest$table[,"PValue"], breaks=50)

hist(dgeTestFilt$table[,"PValue"], breaks=50)
```
### Contact

For any questions or issues, please contact Priyanka at priyankathareja10@gmail.com.




