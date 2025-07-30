.. _usage:

#####
Usage
#####

Calling pipemake without any arguments or with the help flags will display the currently available pipelines (see :ref:`pipelines`):

.. code-block:: bash

    pipemake

Once a pipeline has been selected, you may view the pipeline's arguments by simply adding the pipeline name to the command. 

For example, to view the arguments for the `trim-fastqs` pipeline, you would use the following command:

.. code-block:: bash

    pipemake trim-fastqs

This would display the following information:

.. code-block:: text

    error: one of the arguments --fastq-wildcard --fastq-table --paired-end-sra --single-end-sra is required

    usage: pipemake trim-fastqs (--fastq-wildcard FASTQ_WILDCARD | --fastq-table FASTQ_TABLE | --paired-end-sra PAIRED_END_SRA | --single-end-sra SINGLE_END_SRA) 
                                [--fastq-copy-method {symbolic_link,copy}] 
                                [--min-length MIN_LENGTH] 
                                [--cut-front] 
                                [--cut-tail] 
                                [--cut-right] 
                                [--workflow-dir WORKFLOW_DIR]
                                [--scale-threads SCALE_THREADS] 
                                [--scale-mem SCALE_MEM] 
                                [--resource-yml] 
                                [--singularity-dir SINGULARITY_DIR] 
                                [--no-overwrite]
                                [-h]

    trim-fastqs required arguments:
    --fastq-wildcard FASTQ_WILDCARD
                            Wildcard statement to represent FASTQs (supported wildcards: samples, reads)
    --fastq-table FASTQ_TABLE
                            Table with sample and FASTQs filenames (supported wildcards: samples, reads)
    --paired-end-sra PAIRED_END_SRA
                            Table with SRA accession numbers for paired-end reads (supported wildcards: samples)
    --single-end-sra SINGLE_END_SRA
                            Table with SRA accession numbers for single-end reads (supported wildcards: samples)

    trim-fastqs optional arguments:
    --fastq-copy-method {symbolic_link,copy}
                            Specifies if FASTQs should be copied or symbolically linked. (default: symbolic_link)
    --min-length MIN_LENGTH
                            Minimum length of reads to keep (default: 50)
    --cut-front           Trim from the front of the reads
    --cut-tail            Trim from the tail of the reads
    --cut-right           Trim the right side of the reads, starting from the 5' end
    --workflow-dir WORKFLOW_DIR
                            Assign the workflow directory. If not provided a default will be used.
    --scale-threads SCALE_THREADS
                            Scale the threads for each task
    --scale-mem SCALE_MEM
                            Scale the memory (RAM) for each task
    --resource-yml        Create a seperate resource yaml
    --singularity-dir SINGULARITY_DIR
                            Assign different directory of singularity images
    --no-overwrite        Do not overwrite existing files
    -h, --help            show this help message and exit

You may notice that the arguments are divided into two categories: required and optional.

* **required** is the category that holds arguments required by the pipeline to operate
* **optional** is the category that holds optional arguments that modify the pipeline's behavior.

For the `trim-fastqs` pipeline to operate, the user must provided either: `--fastq-wildcard`, `--fastq-table`, `--paired-end-sra`, or `--single-end-sra` argument.

* `--fastq-wildcard` is used to specify a wildcard statement that represents the FASTQ files
* `--fastq-table` is used to specify a table with sample and FASTQ filenames
* `--paired-end-sra` is used to specify a table with SRA accession numbers for paired-end reads
* `--single-end-sra` is used to specify a table with SRA accession numbers for single-end reads

For example, let's examine the files within the example directory `example/trim-fastqs`:

.. code-block:: bash

    test1_R1.fq.gz
    test1_R2.fq.gz
    test2_R1.fq.gz

As you can see, the directory contains:

* A set of paired-end FASTQ files: `test1_R1.fq.gz`, `test1_R2.fq.gz`
* A single-end FASTQ file: `test2_R1.fq.gz`.

Since these files share a similar naming convention, we can use the `--fastq-wildcard` argument to assign the files. To do this, we must use the supported wildcards as given in the help message: `samples` and `reads`.
If you wanted to perform the `trim-fastqs` pipeline on these files, you may use the following command:

.. code-block:: bash

    pipemake trim-fastqs --fastq-wildcard example/trim-fastqs/{samples}_{reads}.fq.gz --workflow-dir TrimTest

This would generate a snakemake workflow called **TrimTest**.

.. note::
    
    A warning will be displayed if:

    * The input files have inconsistent wildcard usage (such as using paired-end alongside single-end files).
    * No singularity path was specified. Which may result in redundant singularity containers being downloaded.

The workflow directory contains all neccessary files to execute the `trim-fastqs` pipeline on the provided FASTQ samples: **test1** and **test2**. 

The workflow could then be executed using the following command:

.. code-block:: bash
    
    cd TrimTest
    snakemake --use-singularity --cores 4

