# Pyrpipe Setup and Usage Guide

## Introduction

This guide will walk you through setting up a Conda environment for using **Pyrpipe** and related bioinformatics tools. Pyrpipe is a Python package designed for the reproducible analysis of RNA-Seq data. It integrates various popular RNA-Seq analysis tools into a streamlined workflow.

## Prerequisites

Ensure you have Conda installed on your system. If you don't have Conda, you can install it by following the instructions on the Conda website.

## Step 1: Configure Conda Channels

Before installing Pyrpipe, ensure that Conda channels are added in the correct order:

bash
``` conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge ```

`conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge` 

This will ensure that the necessary packages are sourced from the correct channels.

## Step 2: Create a New Conda Environment

Next, create a new Conda environment named `pyrpipe` with Python 3.8:

bash

`conda create -n pyrpipe python=3.8
conda activate pyrpipe`

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

`wget ftp://ftp.ensemblgenomes.org/pub/release46/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/release46/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.46.gtf.gz` 

These commands download the reference genome and GTF file for _Arabidopsis thaliana_ from the Ensembl Plants database.

## Conclusion

You are now ready to use Pyrpipe for your RNA-Seq analysis. For more information on how to use Pyrpipe, please refer to the [official documentation](https://pyrpipe.readthedocs.io/).

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

## Contact

For any questions or issues, please contact Aryaman at aryamansajwan2002@gmail.com.

----------
