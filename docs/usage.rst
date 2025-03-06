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

    usage: pipemake trim-fastqs (--fastq-wildcard FASTQ_WILDCARD | --fastq-table FASTQ_TABLE)
                                 [--fastq-copy-method {symbolic_link,copy}] 
                                 [--fastq-standardized-wildcard FASTQ_STANDARDIZED_WILDCARD] 
                                 [--min-length MIN_LENGTH] 
                                 [--untrimmed-fastq-dir UNTRIMMED_FASTQ_DIR] 
                                 [--trimmed-fastq-dir TRIMMED_FASTQ_DIR]
                                 [--workflow-prefix WORKFLOW_PREFIX]
                                 [--work-dir WORK_DIR]
                                 [--scale-threads SCALE_THREADS]
                                 [--scale-mem SCALE_MEM]
                                 [--resource-yml]
                                 [--singularity-dir SINGULARITY_DIR] 
                                 [-h]

    trim-fastqs required arguments:
    --fastq-wildcard FASTQ_WILDCARD
                            Wildcard statement to represent FASTQs
    --fastq-table FASTQ_TABLE
                            Table with sample and FASTQs filenames

    trim-fastqs paths arguments:
    --untrimmed-fastq-dir UNTRIMMED_FASTQ_DIR
                            Directory to store untrimmed FASTQ files (default: FASTQ/Untrimmed)
    --trimmed-fastq-dir TRIMMED_FASTQ_DIR
                            Directory to store trimmed FASTQ files (default: FASTQ/Trimmed)

    trim-fastqs optional arguments:
    --fastq-copy-method {symbolic_link,copy}
                            Specifies if FASTQs should be copied or symbolically linked.
    --fastq-standardized-wildcard FASTQ_STANDARDIZED_WILDCARD
                            Standardized wildcard statement used to store FASTQs (default: {samples}_{reads}.fq.gz)
    --min-length MIN_LENGTH
                            Minimum length of reads to keep (default: 50)
    --workflow-prefix WORKFLOW_PREFIX
                            Assign the workflow prefix
    --work-dir WORK_DIR   Assign the working directory for snakemake. If not provided, the current directory will be used.
    --scale-threads SCALE_THREADS
                            Scale the threads for each task
    --scale-mem SCALE_MEM
                            Scale the memory (RAM) for each task
    --resource-yml        Create a seperate resource yaml
    --singularity-dir SINGULARITY_DIR
                            Assign different directory of singularity images
    -h, --help            show this help message and exit

You may notice that the arguments are divided into three categories: required, paths, and optional.

* **required** is the category that holds arguments required by the pipeline to operate
* **paths** is a pipeline-defined category that holds path-defining arguments for convenience. Other pipeline-defined categories may exist depending on the pipeline.
* **optional** is the category that holds optional arguments that modify the pipeline's behavior.

For the `trim-fastqs` pipeline to operate, the user must provided either the `--fastq-wildcard` or `--fastq-table` argument.

* `--fastq-wildcard` is used to specify a wildcard statement that represents the FASTQ files
* `--fastq-table` is used to specify a table with sample and FASTQ filenames

For example, let's examine the files within the example directory `example/trim-fastqs`:

.. code-block:: bash

    test1_R1.fq.gz
    test1_R2.fq.gz
    test2_R1.fq.gz

As you can see, the directory contains:

* A set of paired-end FASTQ files: `test1_R1.fq.gz`, `test1_R2.fq.gz`
* A single-end FASTQ file: `test2_R1.fq.gz`.

Since these files share a similar naming convention, we can use the `--fastq-wildcard` argument to assign the files. To do this, we must use the same wildcards as `--fastq-standardized-wildcard`, `samples` and `reads` (see default).

If you wanted to perform the `trim-fastqs` pipeline on these files, you may  use the following command:

.. code-block:: bash

    pipemake trim-fastqs --fastq-wildcard example/trim-fastqs/{samples}_{reads}.fq.gz --workflow-prefix TrimTest

This would generate a snakemake workflow called **TrimTest** that includes the snakemake file **TrimTest.smk**, the configuaration file **TrimTest.yaml**, and the workflow directory **TrimTest**.

.. note::
    
    A warning will be displayed if the input files have inconsistent wildcard usage (such as using paired-end alongside single-end files).

The workflow includes all neccessary files to execute the `trim-fastqs` pipeline on the provided FASTQ samples: **test1** and **test2**. 

The workflow could then be executed using the following command:

.. code-block:: bash

    snakemake -s TrimTest.smk --use-singularity --cores 4

