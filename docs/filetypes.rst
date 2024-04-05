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
* Command-line arguments (input files, configurable terms, pipeline parameters, etc.)
* Steps needed to standardize the input files for the **Pipeline**
* And lastly, the **Modules** used within the **Pipeline**

The following is an example of a **Pipeline** configuration file:

.. code-block:: bash

    pipeline: rnaseq-counts-star
      parser:
        help: Count RNAseq reads within a genome assembly using STAR and featureCounts
        arg-groups:
          rnaseq-input:
            group:
              type: mutually_exclusive
              required: True
            args:
              rnaseq-wildcard:
                help: "Wildcard statement to represent RNAseq FASTQs"
                type: str
              rnaseq-table:
                help: "Table with sample and FASTQs filenames"
                type: str
                action: confirmFile
          basic:
            args:
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
            args:
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
                copy_method: '{rnaseq-copy-method}'
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
                copy_method: '{rnaseq-copy-method}'
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

****************************
Pipeline configuration guide
****************************

A pipeline configuration file begins with the `pipeline` keyword, which is used to define the name of the pipeline. As this name is used to identify a pipeline within `pipemake`, it must be unique. The configuration file then consists of the following required sections: `parser`, `setup`, and `snakefiles`.

.. code-block:: bash

    pipeline: rnaseq-counts-star
      parser:
        ...
      setup:
        ...
      snakefiles:
        ...

parser:
#######

The parser section is used to create the command-line interface for a pipeline. It is divided into the following sub-sections: `help` and `arg-groups`.

help:
*****

The help sub-section is used to define the description of the pipeline, which is displayed when `Pipemake` is run with the `--help` flag.

.. code-block:: bash

    pipeline: rnaseq-counts-star
      parser:
        help: Count RNAseq reads within a genome assembly using STAR and featureCounts

arg-groups:
***********

The `arg-groups` sub-section is used by `pipemake` to define command-line argument groups. The `basic` group is reserved by `pipemake`, arguments within this group will be automatically grouped within `required` or `optional` based on their `required` keyword. Users may place all arguments within the `basic` group or create additional groups as desired. Additional `arg-groups` may be defined as needed.

.. code-block:: bash

    pipeline: rnaseq-counts-star
      parser:
        help: Count RNAseq reads within a genome assembly using STAR and featureCounts
        arg-groups:
          rnaseq-input:
            group:
              type: mutually_exclusive
              required: True
            args:
              rnaseq-wildcard:
                help: "Wildcard statement to represent RNAseq FASTQs"
                type: str
              rnaseq-table:
                help: "Table with sample and FASTQs filenames"
                type: str
                action: confirmFile
          basic:
            args:
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
              assembly-version:
                help: "Assembly Version"
                type: str
                default:
                  str: "v"
                  suffix:
                    - function: jobRandomString
          paths:
            args:
              assembly-dir:
                help: "Directory to store assembly"
                type: str
                default: "Assembly"

group:
======

By default, `arg-groups` are defined as `argparse` argument groups. However, by including the `group` keyword, it becomes possibile to alter this behavior. For example, the `group` keyword within the `rnaseq-input` argument group includes two keywords: 

* `type`: Defines the type of argument group. Currently supported types are `argument_group` and `mutually_exclusive`.

  * `argument_group`: Standard `argparse` argument group. Groups arguments together within the help message (default)
  * `mutually_exclusive`: Only one of the arguments within this group may be used at a time

* `required`: Where applicable, defines if the argument group is required.

Argument groups may be used to organize related arguments within the pipeline help message, for example grouping all path arguments together in `paths`. Mutually exclusive groups are recommended when two or more arguments cannot be used together, such as when a pipeline may accept either a wildcard statement or a table of input files.

args:
=====

Each `arg-groups` also includes a list of `args` that define the command-line arguments. Each argument must have the following keywords:

* `help`: A description of the argument
* `type`: The type of the argument

And the following optional keywords are also supported:

* `required`: If the argument is required (default is False)
* `action`: An action to perform on the argument (see below for supported actions)
* `default`: The default value of the argument (see below for additional options)

.. note::

    Arguments are parsed using `argparse <https://docs.python.org/3/library/argparse.html>`_ and therefore support may be added to allow all of the same options as `argparse`.

action:
-------

At present, `pipemake` supports the following actions:

* `confirmFile`: Require the given string to be a file. If the file does not exist, an error will be raised.
* `confirmDir`: Require the given string to be a directory. If the directory does not exist, an error will be raised.

.. note::

    Additional actions may be added in the future, or updates to `pipemake` to allow for custom actions.

default:
--------

The `default` keyword may be used to define the default value of an argument. In general, the default value may share the same type as the `type` keyword. However, it's also possible to define more complex default values.

.. code-block:: bash

    pipeline: rnaseq-counts-star
      parser:
        help: Count RNAseq reads within a genome assembly using STAR and featureCounts
        arg-groups:
          basic:
            args:
              assembly-version:
                help: "Assembly Version"
                type: str
                default:
                  str: "v"
                  suffix:
                    - function: jobRandomString

In the above example, the `assembly-version` argument has a default value of `v` followed by a random string. This is achieved by using the `suffix` keyword. The `suffix` keyword allows for a list of values to be concatenated to the default value. These values may be either strings or one of the following functions: `jobRandomString` or `jobTimeStamp`.

setup:
######

The `setup` section is used to define the steps needed to standardize the input files for the pipeline. Within the `setup` section, each sub-section is used to define a separate input file(s).

.. code-block:: bash

    pipeline: rnaseq-counts-star
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
                copy_method: '{rnaseq-copy-method}'
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
                copy_method: '{rnaseq-copy-method}'
                gzipped: True
            samples:
              method: "table-file"
              args:
                table_filename: "{rnaseq-table}"

In the above example, the `setup` section includes a sub-section called `rnaseq_input` to standardize RNAseq input files. `rnaseq_input` includes two methods to standardize input files: `wildcard-method` and `table-method`. 

Standardization methods are defined by the following required keywords:

* `input`: Keywords related to the input files to be standardized
  
  * `args`: Contains the command-line arguments needed to standardize the input file(s)

* `standardize`: Keywords that define the standardized method
  
    * `method`: The method used to standardize the input files. Currently supported: `wildcard-str`, `table-file`, and `file-str`.

      * `wildcard-str`: Standardize input file(s) using a wildcard statement
      * `table-file`: Standardize input files within a table file
      * `file-str`: Standardize a single file string

    * `args`: The standardization arguments, which may include the following keywords:

      * `wildcard_str`: The command-line argument of the wildcard statement (only usable with the `wildcard-str` method)
      * `table_filename`: The table filename command-line argument (only usable with the `table-file` method)
      * `input_filename`: The input filename command-line argument (only usable with the `file-str` method)
      * `standardized_filename`: The standardized filename(s). This may be a string with or without a wildcard statements. Should result in a filename(s) specified in a `snakefile` rule
      * `copy_method`: The method used to copy (`copy`) or symbolically link (`symbolic_link`) the input file(s)
      * `gzipped`: If the input file(s) are gzipped (`True`, `False`) or keep the gzipped status of the input file(s) (`None`)
      * `out_dir`: The output directory
      * `work_dir`: The working directory

Standardization methods may also include the following optional keyword:

* `samples`: Keywords that define the samples (only usable with the `wildcard-str` and `table-file` methods)
  
    * `method`: The method used to define the samples. Currently supported: `wildcard-str` and `table-file`.

      * `wildcard-str`: Define samples using a wildcard statement
      * `table-file`: Define samples using a table file

    * `args`: The sample arguments, which may include the following keywords:

      * `wildcard_str`: The command-line argument of the wildcard statement (only usable with the `wildcard-str` method)
      * `sample_wildcard`: The wildcard element used to define the samples (only usable with the `wildcard-str` method)
      * `table_filename`: The table filename command-line argument. Samples are defined using the `sample` column (only usable with the `table-file` method)

snakefiles:
###########

.. code-block:: bash

    pipeline: rnaseq-counts-star
      snakefiles:
        - rna_seq_2pass_star
        - rna_seq_sort
        - rna_seq_feature_counts

The `snakefiles` section is used to define the **Modules** used within the pipeline. **Modules** are defined by the name of the Snakemake file and must be included within the `snakefiles` list.