.. _usage:

#####
Usage
#####

The following sections provide a brief overview of how to use pipemake.

***********
Basic Usage
***********


Calling pipemake without any arguments or with the help flags will display the currently available pipelines (see :ref:`pipelines`):

.. code-block:: bash

    pipemake

Once a pipeline has been selected, you may view the pipeline's arguments by simply adding the pipeline name to the command. 

For example, to view the arguments for the `fastq-filter` pipeline, you would use the following command:

.. code-block:: bash

    pipemake fastq-filter

This would display the following information:

.. code-block:: text

    error: one of the arguments --fastq-wildcard --fastq-table is required

    usage: pipemake fastq-filter (--fastq-wildcard FASTQ_WILDCARD | --fastq-table FASTQ_TABLE)
                                [--fastq-copy-method {symbolic_link,copy}]
                                [--fastq-standardized-wildcard FASTQ_STANDARDIZED_WILDCARD]
                                [--min-length MIN_LENGTH] 
                                [--work-dir WORK_DIR]
                                [--snakemake-job-prefix SNAKEMAKE_JOB_PREFIX]
                                [--unfiltered-fastq-dir UNFILTERED_FASTQ_DIR]
                                [--filtered-fastq-dir FILTERED_FASTQ_DIR]
                                [--scale-threads SCALE_THREADS]
                                [--scale-mem SCALE_MEM] 
                                [--resource-yml]
                                [--singularity-dir SINGULARITY_DIR] [-h]

    fastq-filter required arguments:
    --fastq-wildcard FASTQ_WILDCARD
                            Wildcard statement to represent FASTQs
    --fastq-table FASTQ_TABLE
                            Table with sample and FASTQs filenames

    fastq-filter paths arguments:
    --unfiltered-fastq-dir UNFILTERED_FASTQ_DIR
                            Directory to store unfiltered FASTQ files
    --filtered-fastq-dir FILTERED_FASTQ_DIR
                            Directory to store filtered FASTQ files

    fastq-filter optional arguments:
    --fastq-copy-method {symbolic_link,copy}
                            Socifies if FASTQs should be copied or symbolically linked.
    --fastq-standardized-wildcard FASTQ_STANDARDIZED_WILDCARD
                            Standardized wildcard statement used to store FASTQs
    --min-length MIN_LENGTH
                            Minimum length of reads to keep
    --work-dir WORK_DIR   Assign the working directory for snakemake
    ---workflow-prefix WORKFLOW_PREFIX
                            Assign the worflow prefix
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
* **paths** is a pipeline-defined category that hold path-defining arguments for convenience. Other pipeline-defined categories may exist depending on the pipeline.
* **optional** is the category that holds optional arguments that modify the pipeline's behavior.

For the `fastq-filter` pipeline to operate, the user must provided either the `--fastq-wildcard` or `--fastq-table` argument.

* `--fastq-wildcard` is used to specify a wildcard statement that represents the FASTQ files
* `--fastq-table` is used to specify a table with sample and FASTQ filenames

For example, lets say you have the following FASTQ files in your directory:

.. code-block:: bash

    A1_R1.fq.gz
    B4_R1.fq.gz
    C12_R1.fq.gz
    D04_R1.fq.gz

If you wanted to perform the `fastq-filter` pipeline on these files, you could use the following command:

.. code-block:: bash

    pipemake fastq-filter --fastq-wildcard "{sample}_R1.fq.gz" ---workflow-prefix "FilterTest"

This would generate a snakemake workflow called **FilterTest** that includes the snakemake file **FilterTest.smk**, the configuaration file **FilterTest.yaml**, and the workflow directory **FilterTest**.

The workflow includes all neccessary files to execute the `fastq-filter` pipeline on the provided FASTQ samples: **A1**, **B4**, **C12**, and **D04**.

The workflow could then be executed using the following command:

.. code-block:: bash

    snakemake -s FilterTest.smk --use-singularity --cores 4

