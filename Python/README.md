
# Pyrpipe Setup and Usage Guide

## Introduction

This guide will walk you through setting up a Conda environment for using **Pyrpipe** and related bioinformatics tools. Pyrpipe is a Python package designed for the reproducible analysis of RNA-Seq data. It integrates various popular RNA-Seq analysis tools into a streamlined workflow.

## Prerequisites

Ensure you have Conda installed on your system. If you don't have Conda, you can install it by following the instructions on the Conda website.

## Step 1: Configure Conda Channels

Before installing Pyrpipe, ensure that Conda channels are added in the correct order:

bash

`conda config --add channels defaults`
`conda config --add channels bioconda`
`conda config --add channels conda-forge`

This will ensure that the necessary packages are sourced from the correct channels.

## Step 2: Create a New Conda Environment

Next, create a new Conda environment named `pyrpipe` with Python 3.8:

bash

`conda create -n pyrpipe python=3.8`
`conda activate pyrpipe` 

## Step 3: Install Pyrpipe and Required Tools

With the environment activated, install Pyrpipe along with the required bioinformatics tools:

bash

`conda install -c bioconda pyrpipe star=2.7.7a sra-tools=2.10.9 stringtie=2.1.4 trim-galore=0.6.6 orfipy=0.0.3 salmon=1.4.0` 

This command installs the following tools:

-   **Pyrpipe**: Main Python package for RNA-Seq analysis.
-   **STAR**: RNA-Seq read aligner.
-   **SRA-Tools**: Tools for working with sequence data from NCBI's Sequence Read Archive.
-   **StringTie**: Transcriptome assembly and quantification.
-   **Trim Galore**: A wrapper around Cutadapt and FastQC for quality control and trimming.
-   **Orfipy**: Tool for finding open reading frames (ORFs).
-   **Salmon**: Quantification of transcript abundance.

## Step 4: Download Test Sample Files

To get started with Pyrpipe, you'll need a sample dataset. Run the following commands in your bash terminal to download the necessary files:

bash

`wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz`
`gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz`
`wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz`
`gunzip Arabidopsis_thaliana.TAIR10.46.gtf.gz` 

These commands download the reference genome and GTF file for _Arabidopsis thaliana_ from the Ensembl Plants database.

## Step 5: Simple RNA-Seq Processing Pipeline

RNA-Seq processing with Pyrpipe is straightforward. The following Python script provides a basic example of using Pyrpipe on publicly available RNA-Seq data:

python


```from pyrpipe import sra, qc, mapping, assembly

#define some variables
run_id = 'SRR976159'
working_dir = 'example_output'
gen = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'
ann = 'Arabidopsis_thaliana.TAIR10.46.gtf'
star_index = 'star_index/athaliana'

#initialize objects
#creates a STAR object to use with threads
star = mapping.Star(index=star_index, genome=gen, threads=4)

#use Trim Galore for trimming
trim_galore = qc.Trimgalore()

#Stringtie for assembly
stringtie = assembly.Stringtie(guide=ann)

#create SRA object; this will download FASTQ if it doesn't exist
srr_object = sra.SRA(run_id, directory=working_dir)

#create a pipeline using the objects
srr_object.trim(trim_galore).align(star).assemble(stringtie)

#The assembled transcripts are in srr_object.gtf
print('Final result', srr_object.gtf)
```


### Explanation of the Code

1.  **Imports the required Pyrpipe modules.**
2.  **Lines 3 to 7** define variables for reference files, the output directory, and the STAR index. The output directory will be used to store the downloaded RNA-Seq data and will be the default directory for all the results.
3.  **Creates a STAR object**: It takes the index and genome as parameters. It will automatically verify the index, and if an index is not found, it will use the genome to build one and save it to the index path provided.
4.  **Creates a Trimgalore object**.
5.  **Creates a Stringtie object**.
6.  **Creates an SRA object**: This represents the RNA-Seq data. If the raw data is not available on disk, it auto-downloads it via `fasterq-dump`.
7.  **Pipeline creation**:
    -   `trim()`: Takes a QC type object and performs trimming via `qc.perform_qc` method. The trimmed FASTQ files are updated in the SRA object.
    -   `align()`: Takes a mapping type object and performs alignment via `mapping.perform_alignment` method. The resulting BAM file is stored in `SRA.bam_path`.
    -   `assemble()`: Takes an assembly type object and performs assembly via `mapping.perform_assembly` method. The resulting GTF file is stored in `SRA.gtf`.

This pipeline downloads FASTQ files from NCBI-SRA, uses Trim Galore for trimming, STAR for alignment to the reference genome, and Stringtie for assembly.

## Conclusion

You are now ready to use Pyrpipe for your RNA-Seq analysis. For more information on how to use Pyrpipe, please refer to the [official documentation](https://pyrpipe.readthedocs.io/).

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

## Contact

For any questions or issues, please contact Aryaman at aryamansajwan2002@gmail.com.
