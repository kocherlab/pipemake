.. filetypes:

###################
Pipemake File Types
###################

Pipemake operates using a combination of two file types: Snakemake files and pipeline configuration files.

*************************
Snakemake files (Modules)
*************************

Snakemake files are used to define Pipemake **Modules**. In general, **Modules** follow the same structure and nomenclature as typical Snakemake files. However, Pipemake **Modules** are focused on being reusable. This is achieved by following a few key principles:
* Limiting a **Module** to a collection of rules used to perform a particular task (align reads, call variants, annotate a genome, etc.)
* Consistent input and output usage, which may be either hard-coded or generated from a set of configurable terms
* Using singularity containers to ensure a consistent software environment

By following these principles, Pipemake **Modules** may be easily used in multiple pipelines. For example, a basic snakemake file that aligns reads from a sample to a reference genome might have the following rules:

.. code-block:: python

    rule index_reference:
        input:
            "myref.fasta"
        output:
            "myref.fasta.bwt"
        singularity:
		    "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa index {input}"

    rule align_reads:
        input:
            reads="Sample1.R1.fastq.gz",
            ref="myref.fasta",
            index="myref.fasta.bwt"
        output:
            "Sample1.bam"
        singularity:
		    "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa mem -t 8 {input.ref} {input.reads} | samtools view -bS - > {output}"

To convert this into a Pipemake **Module** with configurable files, we would only need to modify the input and output files as follows:

.. code-block:: python

    rule index_reference:
        input:
            os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa")
        output:
            os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa.bwt")
        singularity:
		    "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa index {input}"

    rule align_reads:
        input:
            reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
            ref=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa"),
            index=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa.bwt")
        output:
            os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.bam")
        singularity:
		    "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa mem -t 8 {input.ref} {input.reads} | samtools view -bS - > {output}"

In this example, we have replaced the unique filenames with the following set of configurable terms: 
* `config['species']`
* `config['assembly_version']`
* `config['paths']['assembly_dir']`
* `config['paths']['rnaseq_fastq_dir']`
* `config['paths']['rnaseq_aligned_bam_dir']`

By consistently using these configurable term (or standardized filenames), it is possible to easily connect multiple **Modules** together to form pipelines. For example, the output of the `align_reads` rule may be used as the input for another **Module** to count reads:

.. code-block:: python

    rule count_reads:
        input:
            bam=os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.bam"),
            gtf=os.path.join(config['paths']['annotation_dir'], f"{config['species']}_{config['annotation_version']}.gtf")
        output:
            os.path.join(config['paths']['rnaseq_count_dir'], "{sample}.counts")
        singularity:
            "docker://quay.io/biocontainers/subread:1.6.4--py36pl5.22.0_0"
        shell:
            "featureCounts -a {input.gtf} -o {output} {input.bam}"

In the above example, the `count_reads` rule uses BAM files stored within the `config['paths']['rnaseq_aligned_bam_dir']` directory. As the name implies, the directory contains aligned RNAseq reads in BAM format and thus we may use this path whenever we need to gain access to them. By consistently using the same terms, such as `config['paths']['rnaseq_aligned_bam_dir']`, we are easily able to connect **Modules** together to form pipelines.

.. note::

    Pipemake is designed to detect configurable terms and will ensure the terms are properly assigned in the configuration file. Configurable terms may also be grouped together in the configuration file. For example, the filepath terms `config['paths']['assembly_dir']`, `config['paths']['rnaseq_fastq_dir']`, and `config['paths']['rnaseq_aligned_bam_dir']` will be stored together within `config['paths']`. Grouping related terms together allows for a more organized configuration file, but is not required.

.. attention::

    While the usage of configurable terms is not required, it is highly recommended.

****************************
Pipeline configuration files
****************************

Pipemake uses YAML-formatted files to define **Pipelines**. These files are used to define the following aspects of a pipeline:
* The **Pipeline** name and description
* Input files
* Configurable terms
* **Pipeline** settings and/or parameters
* Steps needed to standardize the input files for the **Pipeline**
* And lastly, the **Modules** used within the **Pipeline**

The following is an example of a **Pipeline** configuration file:

.. code-block:: bash

    pipeline: rnaseq-counts-star
    parser:
    help: Count RNAseq reads within a genome assembly using STAR and featureCounts
    groups:
        input_parser:
        type: mutually_exclusive
        args:
            required: True
    args:
        params:
        rnaseq-wildcard:
            help: "Wildcard statement to represent RNAseq FASTQs"
            type: str
            group: input_parser
        rnaseq-table:
            help: "Table with sample and FASTQs filenames"
            type: str
            group: input_parser
            action: confirmFile
        rnaseq-copy-method:
            help: "Socifies if RNAseq FASTQs should be copied or symbolically linked."
            choices:
            - 'symbolic_link'
            - 'copy'
            default: 'symbolic_link'
        rnaseq-standardized-wildcard:
            help: "Standardized wildcard statement used to store RNAseq FASTQs"
            type: str
            default: 
            str: "{sample}_{read}.fq.gz"
        assembly-fasta:
            help: "Assembly fasta"
            type: str
            required: True
            action: confirmFile
        assembly-gtf:
            help: "Assembly GTF"
            type: str
            required: True
            action: confirmFile
        read-len:
            help: "Read Length"
            type: int
            required: True
        assembly-version:
            help: "Assembly Version"
            type: str
            default:
            str: "v"
            suffix:
                - function: jobRandomString
        species:
            help: "Species name"
            type: str
            default:
            str: "Sp"
            suffix:
                - function: jobRandomString
        work-dir:
            help: "Assign the working directory for snakemake"
            type: str
            default:
            str: "RNAseqCounts"
            suffix:
                - function: jobTimeStamp
                - function: jobRandomString
        snakemake-job-prefix:
            help: "Assign the snakemake job prefix"
            type: str
            default:
            str: "countSTAR"
            suffix:
                - function: jobTimeStamp
                - function: jobRandomString
        paths:
        assembly-dir:
            help: "Directory to store assembly"
            type: str
            default: "Assembly"
        index-dir:
            help: "Directory to store indices"
            type: str
            default: "Indices"
        rnaseq-fastq-dir:
            help: "Directory to store the FASTQs files"
            type: str
            default: "RNAseq/FASTQs"
        rnaseq-splice-aligned-dir:
            help: "Directory to store BAM files"
            type: str
            default: "RNAseq/SpliceJunctions/Aligned"
        rnaseq-bam-dir:
            help: "Directory to store BAM files"
            type: str
            default: "RNAseq/BAMs"
        rnaseq-aligned-bam-dir:
            help: "Directory to store sorted BAM files"
            type: str
            default: "RNAseq/BAMs/Aligned"
        rnaseq-sorted-bam-dir:
            help: "Directory to store sorted BAM files"
            type: str
            default: "RNAseq/BAMs/Sorted"
        rnaseq-count-dir:
            help: "Directory to store RNAseq counts"
            type: str
            default: "RNAseq/Counts" 
    setup:
    rnaseq_input:
        wildcard-method:
        input:
            args:
            - "work-dir"
            - "rnaseq-wildcard"
            - "rnaseq-standardized-wildcard"
            - "rnaseq-fastq-dir"
        standardize:
            method: "wildcard-str"
            args:
            wildcard_str: "{rnaseq-wildcard}"
            standardized_filename: "{rnaseq-standardized-wildcard}"
            out_dir: "{rnaseq-fastq-dir}"
            work_dir: '{work-dir}'
            copy_method: '{rnaseq_copy_method}'
            gzipped: True
        samples:
            method: "wildcard-str"
            args:
            wildcard_str: "{rnaseq-wildcard}"
            sample_wildcard: 'sample'

        table-method:
        input:
            args:
            - "work-dir"
            - "rnaseq-table"
            - "rnaseq-standardized-wildcard"
            - "rnaseq-fastq-dir"
        standardize:
            method: "table-file"
            args:
            table_filename: "{rnaseq-table}"
            standardized_filename: "{rnaseq-standardized-wildcard}"
            out_dir: "{rnaseq-fastq-dir}"
            work_dir: '{work-dir}'
            copy_method: '{rnaseq_copy_method}'
            gzipped: True
        samples:
            method: "table-file"
            args:
            table_filename: "{rnaseq-table}"
    
    assembly_input:
        file-method:
        input:
            args:
            - "work-dir"
            - "assembly-fasta"
            - "assembly-dir"
        standardize:
            method: "file-str"
            args:
            input_filename: "{assembly-fasta}"
            standardized_filename: "{species}_{assembly_version}.fa"
            out_dir: "{assembly-dir}"
            work_dir: '{work-dir}'
            gzipped: False
    
    gtf_input:
        file-method:
        input:
            args:
            - "work-dir"
            - "assembly-gtf"
            - "assembly-dir"
        standardize:
            method: "file-str"
            args:
            input_filename: "{assembly-gtf}"
            standardized_filename: "{species}_{assembly_version}.gtf"
            out_dir: "{assembly-dir}"
            work_dir: '{work-dir}'
            gzipped: False

    snakefiles:
    - rna_seq_2pass_star
    - rna_seq_sort
    - rna_seq_feature_counts
