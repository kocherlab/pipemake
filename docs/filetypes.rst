.. _filetypes:

###################
Pipeline File Types
###################

***********
Table Files
***********

While many pipemake pipeline may be run using a wildcard argument (e.g. `--fastq-wildcard`) to assign input files, it's also possible to use a table file to assign input files. 

Table files are a convenient way to assign input files to a pipeline, especially when the files are not easily assigned using a wildcard argument due to inconsistent naming conventions or file locations.

Table files are tab-delimited files that typically contain two coluumn types: **sample** and **file**.

* The **sample** column is used to define the sample name. Please note that sample names should be unique and only a single **sample** column is allowed in a table file.
* **file** columns that contains the file paths for the samples. Multiple **file** columns are allowed in a table file.

Columns should be named based on the relevant standardized wildcard argument. These names can be viewed by running the pipeline with the `--help` argument. See the `--fastq-table` argument below.

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

In this example, a table file would contain the following columns: **samples** and **reads**. Let's look at an example table file:

.. literalinclude:: _static/table_file.tsv

As you can see, the table file actually contains three columns: **samples**, **reads:R1**, and **reads:R2**. File columns may use `:` to assign a file to a specific wildcard. 

For example, the file `lane01_test1_read_1.fq.gz` is assigned as the `reads = R1` for `samples = test1`. When these values are standardized using the template `{samples}_{reads}.fq.gz`, the file will be stored as `test1_R1.fq.gz` in the appropriate directory.

***********
Model Files
***********

Model files are used to define the population models. For more details on the model files, please refer to the `Popgen Pipeline Platform <https://ppp.readthedocs.io/en/latest/PPP_pages/model.html>`_.

In addition to the JSON files used by the Popgen Pipeline Platform, pipemake supports Model files in the yaml format. Below is an example of a model file in the yaml format:

.. code-block:: bash

    - name: 2Pop
      pops:
        Verus:
          inds:
          - "Pan_troglodytes_verus-9668_Bosco"
          - 'Pan_troglodytes_verus-9730_Donald'
          - 'Pan_troglodytes_verus-A956_Jimmie'
          - 'Pan_troglodytes_verus-Clint'
          - 'Pan_troglodytes_verus-X00100_Koby'
        Troglodytes:
          inds:
          - 'Pan_troglodytes_troglodytes-A957_Vaillant'
          - 'Pan_troglodytes_troglodytes-A958_Doris'
          - 'Pan_troglodytes_troglodytes-A959_Julie'
          - 'Pan_troglodytes_troglodytes-A960_Clar'
          
    - name: 3Pop
        pops:
          Verus:
            inds:
            - "Pan_troglodytes_verus-9668_Bosco"
            - 'Pan_troglodytes_verus-9730_Donald'
            - 'Pan_troglodytes_verus-A956_Jimmie'
            - 'Pan_troglodytes_verus-Clint'
            - 'Pan_troglodytes_verus-X00100_Koby'
          Troglodytes:
            inds:
            - 'Pan_troglodytes_troglodytes-A957_Vaillant'
            - 'Pan_troglodytes_troglodytes-A958_Doris'
            - 'Pan_troglodytes_troglodytes-A959_Julie'
            - 'Pan_troglodytes_troglodytes-A960_Clar'
          Schweinfurthii:
            inds:
            - 'Pan_troglodytes_schweinfurthii-100037_Vincent'
            - 'Pan_troglodytes_schweinfurthii-100040_Andromeda'
            - 'Pan_troglodytes_schweinfurthii-9729_Harriet'
            - 'Pan_troglodytes_schweinfurthii-A910_Bwambale'
            - 'Pan_troglodytes_schweinfurthii-A911_Kidongo'
            - 'Pan_troglodytes_schweinfurthii-A912_Nakuu'