RNA-Seq processing is as easy as creating required objects and executing required functions. The following python
script provides a basic example of using pyrpipe on publicly available RNA-Seq data


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

