.. _usage:

#####
Usage
#####

Calling pipemake without any arguments or with the help flags will display the currently available pipelines (see :ref:`pipelines`):

.. code-block:: bash

    pipemake

Once a pipeline has been selected, you may view the pipeline's arguments by simply adding the pipeline name to the command. 

For example, to view the arguments for the `fastq-filter` pipeline, you would use the following command:

.. code-block:: bash

    pipemake fastq-filter

This would display the following information:

.. code-block:: text

    usage: pipemake fastq-filter (--fastq-wildcard FASTQ_WILDCARD | --fastq-table FASTQ_TABLE)
                                 [--fastq-copy-method {symbolic_link,copy}] 
                                 [--fastq-standardized-wildcard FASTQ_STANDARDIZED_WILDCARD] 
                                 [--min-length MIN_LENGTH] 
                                 [--unfiltered-fastq-dir UNFILTERED_FASTQ_DIR] 
                                 [--filtered-fastq-dir FILTERED_FASTQ_DIR]
                                 [--workflow-prefix WORKFLOW_PREFIX]
                                 [--work-dir WORK_DIR]
                                 [--scale-threads SCALE_THREADS]
                                 [--scale-mem SCALE_MEM]
                                 [--resource-yml]
                                 [--singularity-dir SINGULARITY_DIR] 
                                 [-h]

    fastq-filter required arguments:
    --fastq-wildcard FASTQ_WILDCARD
                            Wildcard statement to represent FASTQs
    --fastq-table FASTQ_TABLE
                            Table with sample and FASTQs filenames

    fastq-filter paths arguments:
    --unfiltered-fastq-dir UNFILTERED_FASTQ_DIR
                            Directory to store unfiltered FASTQ files (default: FASTQ/Unfiltered)
    --filtered-fastq-dir FILTERED_FASTQ_DIR
                            Directory to store filtered FASTQ files (default: FASTQ/Filtered)

    fastq-filter optional arguments:
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

For the `fastq-filter` pipeline to operate, the user must provided either the `--fastq-wildcard` or `--fastq-table` argument.

* `--fastq-wildcard` is used to specify a wildcard statement that represents the FASTQ files
* `--fastq-table` is used to specify a table with sample and FASTQ filenames

For example, let's examine the files within the example directory `example/fastq-filter`:

.. code-block:: bash

    test1_R1.fq.gz
    test1_R2.fq.gz
    test2_R1.fq.gz

As you can see, the directory contains:

* A set of paired-end FASTQ files: `test1_R1.fq.gz`, `test1_R2.fq.gz`
* A single-end FASTQ file: `test2_R1.fq.gz`.

Since these files share a similar naming convention, we can use the `--fastq-wildcard` argument to assign the files. To do this, we must use the same wildcards as `--fastq-standardized-wildcard`, `samples` and `reads` (see default).

If you wanted to perform the `fastq-filter` pipeline on these files, you may  use the following command:

.. code-block:: bash

    pipemake fastq-filter --fastq-wildcard example/fastq-filter/{samples}_{reads}.fq.gz --workflow-prefix FilterTest

This would generate a snakemake workflow called **FilterTest** that includes the snakemake file **FilterTest.smk**, the configuaration file **FilterTest.yaml**, and the workflow directory **FilterTest**.

.. note::
    
    A warning will be displayed if the input files have inconsistent wildcard usage (such as using paired-end alongside single-end files).

The workflow includes all neccessary files to execute the `fastq-filter` pipeline on the provided FASTQ samples: **test1** and **test2**. 

The workflow could then be executed using the following command:

.. code-block:: bash

    snakemake -s FilterTest.smk --use-singularity --cores 4

