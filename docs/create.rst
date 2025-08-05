.. filetypes:

##################
Creating Pipelines
##################

pipemake pipelines operates using a combination of two file types: Snakemake files and pipeline configuration files.

*************************
Snakemake files (Modules)
*************************

Snakemake files are used to define pipemake **Modules**. In general, **Modules** follow the same structure and nomenclature as typical Snakemake files. However, pipemake **Modules** are focused on being reusable. This is achieved by following a few key principles:

* Limiting a **Module** to a collection of rules used to perform a particular task (align reads, call variants, annotate a genome, etc.)
* Consistent input and output usage 
* Defining configurable terms (e.g. samples, wildcards, parameters, etc.) from the config
* Using singularity containers to ensure a consistent software environment

By following these principles, pipemake **Modules** may be easily used in multiple pipelines. For example, the following are examples of a **Module** that aligns RNAseq reads to a reference genome using BWA and outputs the results in BAM format.

In this first example, the **Module** only requires `samples` to be defined in the configuration file.

.. code-block::

    rule all:
        input:
            expand("reSEQ/BAM/Aligned/{sample}.bam", sample=config['samples'])

    rule index_reference:
        input:
            f"Assembly/assembly.fa"
        output:
            f"Assembly/assembly.fa.bwt"
        singularity:
        "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa index {input}"

    rule align_reads:
        input:
            reads="reSEQ/FASTQ/{sample}_R1.fastq.gz",
            ref=f"Assembly/assembly.fa"
            index=f"Assembly/assembly.fa.bwt"
        output:
            "reSEQ/BAM/Aligned/{sample}.bam"
        singularity:
        "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa mem -t 8 {input.ref} {input.reads} | samtools view -bS - > {output}"

In this next example, the **Module** uses additional configurable terms to define the `species` and `assembly_version`. While additional configurable do require additional command-line arguments, they allow for greater flexibility and easier reporting.

.. code-block::

    rule all:
        input:
            expand("reSEQ/BAM/Aligned/{sample}.bam", sample=config['samples'])

    rule index_reference:
        input:
            f"Assembly/{config['species']}_{config['assembly_version']}.fa"
        output:
            f"Assembly/{config['species']}_{config['assembly_version']}.fa.bwt"
        singularity:
		    "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa index {input}"

    rule align_reads:
        input:
            reads="reSEQ/FASTQ/{sample}_R1.fastq.gz",
            ref=f"Assembly/{config['species']}_{config['assembly_version']}.fa"
            index=f"Assembly/{config['species']}_{config['assembly_version']}.fa.bwt"
        output:
            "reSEQ/BAM/Aligned/{sample}.bam"
        singularity:
		    "docker://quay.io/biocontainers/bwa:0.7.8"
        shell:
            "bwa mem -t 8 {input.ref} {input.reads} | samtools view -bS - > {output}"

By consistently using configurable terms (or standardized filenames), it is possible to easily connect multiple **Modules** together to form pipelines. For example, the output of the `align_reads` rule may be used as the input for another **Module** to sort the BAM files:

.. code-block::

    rule all:
        input:
            expand("reSEQ/BAM/Sorted/{sample}.bam", sample=config['samples'])

    rule sort_reads:
        input:
            "reSEQ/BAM/Aligned/{sample}.bam",
        output:
            "reSEQ/BAM/Sorted/{sample}.bam",
        singularity:
            "docker://quay.io/biocontainers/samtools:1.9"
        shell:
            "samtools sort {input} -o {output}"

.. note::

    pipemake is designed to detect configurable terms and will ensure the terms are properly assigned in the configuration file. Configurable terms may also be grouped together, `config['assembly']['species']` and `config['assembly']['assembly_version']`, if desited.

****************************
Pipeline configuration files
****************************

pipemake uses YAML-formatted files to define **Pipelines**. These files are used to define the following aspects of a pipeline:

* The **Pipeline** name and version
* Command-line arguments (description, input files, configurable terms, pipeline parameters, etc.)
* Steps needed to standardize the input files for the **Pipeline**
* And lastly, the **Modules** and **Links** requried for the **Pipeline**

The following is an example of a **Pipeline** configuration file:

.. code-block::

    pipeline: rnaseq-counts-star
    version: 1.0
    parser:
      help: Count RNAseq reads within a genome assembly using STAR and featureCounts
      arg-groups:
        basic:
          mutually-exclusive-groups:
            input-parser:
              required: True
          wildcards-args:
            rnaseq-standardized-wildcard: "{sample}_{read}.fq.gz"
          args:
            rnaseq-wildcard:
              help: "Wildcard statement to represent RNAseq FASTQs"
              type: str
              mutually-exclusive: 'input-parser'
              wildcards: rnaseq-standardized-wildcard
            rnaseq-table:
              help: "Table with sample and FASTQs filenames"
              type: str
              action: confirmFile
              mutually-exclusive: 'input-parser'
              wildcards: rnaseq-standardized-wildcard
            rnaseq-copy-method:
              help: "Specifies if RNAseq FASTQs should be copied or symbolically linked."
              choices:
                - 'symbolic_link'
                - 'copy'
              default: 'symbolic_link'
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
    setup:
      rnaseq_input:
        methods:
          wildcard-str: "{rnaseq-wildcard}"
          table-file: "{rnaseq-table}"
        args:
          standardized_filename: "RNAseq/FASTQ/{rnaseq-standardized-wildcard}"
          copy_method: '{rnaseq-copy-method}'
          gzipped: True
          sample_keywords:
            'samples'

      assembly_input:
        methods:
          file-str: "{assembly-fasta}"
        args:
          standardized_filename: "Assembly/{species}_{assembly_version}.fa"
          gzipped: False
      
      gtf_input:
        methods:
          file-str: "{assembly-gtf}"
        args:
          standardized_filename: "Assembly/{species}_{assembly_version}.gtf"
          gzipped: False
    
    snakemake:
      modules:
        - fastq_trim_fastp
        - rna_seq_2pass_star
        - rna_seq_sort
        - rna_seq_feature_counts
      links:
        - input: fastp_single_end
          output: star_single_end_p1
        - input: fastp_pair_end
          output: star_pair_end_p1
          file_mappings:
          - input: r1_reads
            output: r1_reads
          - input: r2_reads
            output: r2_reads


****************************
Pipeline configuration guide
****************************

A pipeline configuration file begins with the `pipeline` keyword, which is used to define the name of the pipeline. As this name is used to identify a pipeline within pipemake, it must be unique. Next is the `version` keyword, which is used to define the version of the pipeline and is included to track changes to the pipeline over time. 

The configuration file then consists of the following required sections: `parser`, `setup`, and `snakemake`.

.. code-block::

    pipeline: rnaseq-counts-star
    version: 1.0
    parser:
      ...
    setup:
      ...
    snakemake:
      ...

parser:
#######

The parser section is used to create the command-line interface for a pipeline. It is divided into the following sub-sections: `help` and `arg-groups`.

help:
*****

The help sub-section is used to define the description of the pipeline, which is displayed when pipemake is run with the `--help` flag.

.. code-block::

    parser:
      help: Count RNAseq reads within a genome assembly using STAR and featureCounts

arg-groups:
***********

The `arg-groups` sub-section is used by pipemake to define command-line argument groups. The `basic` group is reserved by pipemake, arguments within this group will be automatically grouped within `required` or `optional` based on their `required` keyword. Users may place all arguments within the `basic` group or create additional groups as desired. Additional `arg-groups` may be defined as needed to organize related arguments, such as parameters for a particular software package.

.. code-block::

    parser:
      arg-groups:
        basic:
          mutually-exclusive-groups:
            input-parser:
              required: True
          wildcards-args:
            rnaseq-standardized-wildcard: "{sample}_{read}.fq.gz"
          args:
            rnaseq-wildcard:
              help: "Wildcard statement to represent RNAseq FASTQs"
              type: str
              mutually-exclusive: 'input-parser'
              wildcards: rnaseq-standardized-wildcard
            rnaseq-table:
              help: "Table with sample and FASTQs filenames"
              type: str
              action: confirmFile
              mutually-exclusive: 'input-parser'
              wildcards: rnaseq-standardized-wildcard
            rnaseq-copy-method:
              help: "Specifies if RNAseq FASTQs should be copied or symbolically linked."
              choices:
                - 'symbolic_link'
                - 'copy'
              default: 'symbolic_link'
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

mutually-exclusive-groups:
==========================

Each `arg-groups` may use the `mutually-exclusive-groups` keyword to define mutually exclusive arguments to ensure that only one of the arguments within a group may be used at a time. This is useful when a pipeline accepts different types of input, such as a wildcard statement or a table of input files. To create a `mutually-exclusive-group`, a user is only required to name the group.

.. code-block::

    parser:
      arg-groups:
        basic:
          mutually-exclusive-groups:
            input-parser:
              required: True

In this example, pipemake will create a single `mutually-exclusive-group` called `input-parser`. Currently, `mutually-exclusive-groups` supports the following keywords:

* `required`: Defines if the `mutually-exclusive-group` is required (default is `False`)

.. note::

    Please note that if a `mutually-exclusive-group` is placed within the `basic` group the `required` keyword will be used to place the arguments within `required` or `optional`.

.. attention::

    At present, pipemake requires that the name of `mutually-exclusive-groups` to be unique among all `arg-groups`.


wildcards-args:
===============

The `wildcards-args` keyword is used to define a wildcard statement that may then be used by multiple arguments within the `arg-groups`. This is useful when a pipeline supports multiple input methods that should be standardized to a common naming convention.

.. code-block::

    parser:
      arg-groups:
        basic:
          wildcards-args:
            rnaseq-standardized-wildcard: "{sample}_{read}.fq.gz"

In this example, we defined a wildcard statement called `rnaseq-standardized-wildcard`. This wildcard statement is used to standardize the naming of RNAseq FASTQ files. The wildcard statement is defined as `"{sample}_{read}.fq.gz"`, where `{sample}` represents the sample name and `{read}` represents the read type (e.g. R1, R2).


args:
=====

Each `arg-groups` requires a list of `args` that define the command-line arguments. Each argument must have the following keywords:

* `help`: A description of the argument
* `type`: The type of the argument

And the following optional keywords are also supported:

* `required`: If the argument is required (default is `False``)
* `choices`: A list of choices for the argument
* `mutually-exclusive`: Adds the argument to the specified `mutually-exclusive-group`
* `wildcards`: Requires the argument to use the wildcards given in the specified `wildcards-args`
* `action`: An action to perform on the argument (see below for supported actions)
* `default`: The default value of the argument (see below for additional options)

.. note::

    Arguments are parsed using `argparse <https://docs.python.org/3/library/argparse.html>`_ and therefore support may be added to allow all of the same options as `argparse`.

action:
-------

Beyond built-in actions, `pipemake` also supports the following actions:

* `confirmFile`: Require the given string to be a file. If the file does not exist, an error will be raised.
* `confirmDir`: Require the given string to be a directory. If the directory does not exist, an error will be raised.

.. note::

    Additional actions may be added in the future, or updates to pipemake to allow for custom actions.

default:
--------

The `default` keyword may be used to define the default value of an argument. In general, the default value may share the same type as the `type` keyword. However, it's also possible to define more complex default values.

.. code-block::

    parser:
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

The `setup` section is used to define arguments that require input standardization. Each standardization argument includes the following keywords: `methods`, `args`, and optionally `snakefiles`.

.. code-block::

    setup:
      rnaseq_input:
        methods:
          wildcard-str: "{rnaseq-wildcard}"
          table-file: "{rnaseq-table}"
        args:
          standardized_filename: "RNAseq/FASTQ/{rnaseq-standardized-wildcard}"
          copy_method: '{rnaseq-copy-method}'
          gzipped: True
          sample_keywords:
            'samples'

In the above example, the `setup` section includes a single argument called `rnaseq_input`, which includes two methods to standardize input files: `wildcard-str` and `table-file`. 

Standardization methods are defined by the following required keywords:

* `methods`: The supported methods and their associated command-line argument

    * `wildcard-str`: Standardize input file(s) using the specified wildcard statement - `"{rnaseq-wildcard}"`
    * `table-file`: Standardize input files within the specified table file - `"{rnaseq-table}"`
    * `file-str`: Standardize the specified file - `"{assembly-fasta}"`
    * `dir-str`: Standardize the specified directory - e.g. `"{index-dir}"`
  
* `args`: Contains the command-line arguments needed for the methods standardize the input file(s)

    * `standardized_filename` or `standardized_directory`: The standardized filename(s) or directory. This may be a string with or without a wildcard statements. Should result in a filename(s) specified in a `snakemake` module
    * `copy_method`: The method used to copy (`copy`) or symbolically link (`symbolic_link`) the input file(s)
    * `gzipped`: If the standardized file(s) should be gzipped (`True`, `False`) or keep the gzipped status of the input file(s) (`None`)
    * `sample_keywords`: A list of keywords that should be treated as samples (optional)

* `snakefiles`: An optional list of Snakemake modules if the current standardization method is used

snakemake:
##########

The `snakemake` section is used to define the Snakemake modules used within the pipeline. This section includes the following sub-sections: `modules` and `links`.

.. code-block::

    snakemake:
      modules:
        - fastq_trim_fastp
        - rna_seq_2pass_star
        - rna_seq_sort
        - rna_seq_feature_counts
      links:
        - input: fastp_single_end
          output: star_single_end_p1
        - input: fastp_pair_end
          output: star_pair_end_p1
          file_mappings:
          - input: r1_reads
            output: r1_reads
          - input: r2_reads
            output: r2_reads


modules:
--------

The `modules` sub-section is used to define the Snakemake modules used within the pipeline. Each module is defined by the name of the Snakemake file (e.g. `rna_seq_2pass_star.smk`), inclusion of the `.smk` extension is optional.

links:
------

The `links` sub-section is used to define the links between Snakemake rules. You think of links as a way to connect two rules that normally would not be connected. This is useful when rules have inconsistent filenames. Let's examine a single link:

.. code-block::

    snakemake:
      links:
        - input: fastp_pair_end
          output: star_pair_end_p1
          file_mappings:
          - input: r1_reads
            output: r1_reads
          - input: r2_reads
            output: r2_reads

In this link, the `input` keyword indicates the input rule for the link (`fastp_pair_end`), whereas `output` keyword indicates the output rule for the link (`star_pair_end_p1`). This rule would then connect the output of the `fastp_pair_end` rule to the input of the `star_pair_end_p1` rule. The `file_mappings` keyword is used to define the mapping of input and output files between the two rules. This is useful when the keywords used by the rules differ.